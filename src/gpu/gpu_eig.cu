// CUDA-C batched generalized symmetric-definite eigensolver, callable from the
// gfortran-built xtb via iso_c_binding. This is plain CUDA C compiled with nvcc
// into an object linked into xtb -- it sidesteps nvfortran (which cannot build
// xtb's dependency tree) while still giving real GPU diagonalization.
//
//   solves  H_k C_k = S_k C_k diag(W_k),  k = 1..nbatch   (itype=1, upper)
//
//   H : in  = symmetric Hamiltonian blocks (n*n column-major, nbatch of them)
//       out = eigenvectors C_k (same layout) -- overwritten in place
//   S : symmetric-positive-definite overlap blocks (consumed on device)
//   W : out = eigenvalues, ascending, n per system
//
// Per-system cusolverDnDsygvd with device buffers + workspace + handle reused
// across the batch (data staged per system). Returns 0 on success.
#include <cusolverDn.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <mutex>
#include <cstdio>
#include <cstdlib>
#include <chrono>

namespace {

struct SolverContext {
   cusolverDnHandle_t handle = nullptr;
   cublasHandle_t blas = nullptr;
   double *dA = nullptr;
   double *dB = nullptr;
   double *dW = nullptr;
   double *dWork = nullptr;
   double *dScaled = nullptr;
   double *dP = nullptr;
   double *dH0 = nullptr;
   double *dShift = nullptr;
   double *dQsh = nullptr;
   int *dMatlist = nullptr;
   int *dAo2sh = nullptr;
   int *dInfo = nullptr;
   int capacity = 0;
   int matrix_list_capacity = 0;
   int shell_capacity = 0;
   int lwork = 0;
   unsigned long long solves = 0;

   ~SolverContext()
   {
      cudaFree(dA);
      cudaFree(dB);
      cudaFree(dW);
      cudaFree(dWork);
      cudaFree(dScaled);
      cudaFree(dP);
      cudaFree(dH0);
      cudaFree(dShift);
      cudaFree(dQsh);
      cudaFree(dMatlist);
      cudaFree(dAo2sh);
      cudaFree(dInfo);
      if (blas) cublasDestroy(blas);
      if (handle) cusolverDnDestroy(handle);
   }
};

SolverContext ctx;
std::mutex ctx_mutex;

int cuda_ok(cudaError_t status, int code)
{
   return status == cudaSuccess ? 0 : code;
}

int ensure_capacity(int n)
{
   if (!ctx.handle &&
       cusolverDnCreate(&ctx.handle) != CUSOLVER_STATUS_SUCCESS) return -1;
   if (!ctx.blas && cublasCreate(&ctx.blas) != CUBLAS_STATUS_SUCCESS) return -14;
   if (n <= ctx.capacity) return 0;

   cudaFree(ctx.dA); ctx.dA = nullptr;
   cudaFree(ctx.dB); ctx.dB = nullptr;
   cudaFree(ctx.dW); ctx.dW = nullptr;
   cudaFree(ctx.dWork); ctx.dWork = nullptr;
   cudaFree(ctx.dScaled); ctx.dScaled = nullptr;
   cudaFree(ctx.dP); ctx.dP = nullptr;
   cudaFree(ctx.dAo2sh); ctx.dAo2sh = nullptr;
   cudaFree(ctx.dInfo); ctx.dInfo = nullptr;

   const size_t nn = (size_t)n * (size_t)n;
   if (cuda_ok(cudaMalloc((void**)&ctx.dA, sizeof(double) * nn), -2)) return -2;
   if (cuda_ok(cudaMalloc((void**)&ctx.dB, sizeof(double) * nn), -3)) return -3;
   if (cuda_ok(cudaMalloc((void**)&ctx.dW, sizeof(double) * n), -4)) return -4;
   if (cuda_ok(cudaMalloc((void**)&ctx.dInfo, sizeof(int)), -5)) return -5;
   if (cuda_ok(cudaMalloc((void**)&ctx.dScaled, sizeof(double) * nn), -15)) return -15;
   if (cuda_ok(cudaMalloc((void**)&ctx.dP, sizeof(double) * nn), -16)) return -16;
   if (cuda_ok(cudaMalloc((void**)&ctx.dAo2sh, sizeof(int) * n), -39)) return -39;

   if (cusolverDnDsygvd_bufferSize(
         ctx.handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR,
         CUBLAS_FILL_MODE_UPPER, n, ctx.dA, n, ctx.dB, n, ctx.dW,
         &ctx.lwork) != CUSOLVER_STATUS_SUCCESS) return -6;
   if (cuda_ok(cudaMalloc((void**)&ctx.dWork,
                          sizeof(double) * ctx.lwork), -7)) return -7;
   ctx.capacity = n;
   return 0;
}

int ensure_scf_capacity(int nmat, int nshell)
{
   if (nmat > ctx.matrix_list_capacity) {
      cudaFree(ctx.dH0); ctx.dH0 = nullptr;
      cudaFree(ctx.dMatlist); ctx.dMatlist = nullptr;
      if (cudaMalloc((void**)&ctx.dH0, sizeof(double) * nmat) != cudaSuccess)
         return -22;
      if (cudaMalloc((void**)&ctx.dMatlist, sizeof(int) * 2 * nmat) != cudaSuccess)
         return -23;
      ctx.matrix_list_capacity = nmat;
   }
   if (nshell > ctx.shell_capacity) {
      cudaFree(ctx.dShift); ctx.dShift = nullptr;
      cudaFree(ctx.dQsh); ctx.dQsh = nullptr;
      if (cudaMalloc((void**)&ctx.dShift, sizeof(double) * nshell) != cudaSuccess)
         return -24;
      if (cudaMalloc((void**)&ctx.dQsh, sizeof(double) * nshell) != cudaSuccess)
         return -25;
      ctx.shell_capacity = nshell;
   }
   return 0;
}

__global__ void scale_columns(int n, const double *C, const double *f,
                              double *scaled)
{
   const size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
   const size_t nn = (size_t)n * n;
   if (idx < nn) {
      const int col = (int)(idx / n);
      scaled[idx] = C[idx] * f[col];
   }
}

__global__ void build_isotropic_h1_kernel(
      int n, int nmat, const int *matlist, const double *H0,
      const double *S, const double *shift, const int *ao2sh,
      double autoev, double *H)
{
   const int m = blockIdx.x * blockDim.x + threadIdx.x;
   if (m >= nmat) return;
   const int i = matlist[2 * m] - 1;
   const int j = matlist[2 * m + 1] - 1;
   const int packed = j + (i + 1) * i / 2;
   const int ish = ao2sh[i] - 1;
   const int jsh = ao2sh[j] - 1;
   const double h1 = -0.5 * S[j + (size_t)i * n] * autoev *
                     (shift[ish] + shift[jsh]);
   const double value = H0[packed] + h1;
   H[j + (size_t)i * n] = value;
   H[i + (size_t)j * n] = value;
}

__global__ void mulliken_shell_kernel(
      int n, const int *ao2sh, const double *S, const double *P,
      double *qsh)
{
   const int i = blockIdx.x * blockDim.x + threadIdx.x;
   if (i >= n) return;
   const int ish = ao2sh[i] - 1;
   atomicAdd(&qsh[ish], P[i + (size_t)i * n] * S[i + (size_t)i * n]);
   for (int j = 0; j < i; ++j) {
      const double ps = P[j + (size_t)i * n] * S[j + (size_t)i * n];
      atomicAdd(&qsh[ish], ps);
      atomicAdd(&qsh[ao2sh[j] - 1], ps);
   }
}

} // namespace

extern "C" int gpu_sygvd_batch(int n, int nbatch,
                               double* H, const double* S, double* W)
{
   if (n <= 0 || nbatch <= 0) return 0;
   std::lock_guard<std::mutex> lock(ctx_mutex);

   int rc = ensure_capacity(n);
   if (rc != 0) return rc;
   const size_t nn = (size_t)n * (size_t)n;

   for (int k = 0; k < nbatch; ++k) {
      if (cudaMemcpy(ctx.dA, H + (size_t)k * nn, sizeof(double) * nn,
                     cudaMemcpyHostToDevice) != cudaSuccess) return -8;
      if (cudaMemcpy(ctx.dB, S + (size_t)k * nn, sizeof(double) * nn,
                     cudaMemcpyHostToDevice) != cudaSuccess) return -9;
      cusolverStatus_t st = cusolverDnDsygvd(ctx.handle, CUSOLVER_EIG_TYPE_1,
            CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,
            n, ctx.dA, n, ctx.dB, n, ctx.dW, ctx.dWork, ctx.lwork, ctx.dInfo);
      if (st != CUSOLVER_STATUS_SUCCESS) return (int)st;

      int info = 0;
      if (cudaMemcpy(&info, ctx.dInfo, sizeof(int),
                     cudaMemcpyDeviceToHost) != cudaSuccess) return -10;
      if (info != 0) return 1000 + info;

      // Copy eigenvectors (in dA) and eigenvalues (in dW) back to the host.
      if (cudaMemcpy(H + (size_t)k * nn, ctx.dA, sizeof(double) * nn,
                     cudaMemcpyDeviceToHost) != cudaSuccess) return -11;
      if (cudaMemcpy(W + (size_t)k * n, ctx.dW, sizeof(double) * n,
                     cudaMemcpyDeviceToHost) != cudaSuccess) return -12;
      ++ctx.solves;
   }
   return cudaDeviceSynchronize() == cudaSuccess ? 0 : -13;
}

extern "C" unsigned long long gpu_sygvd_solve_count()
{
   std::lock_guard<std::mutex> lock(ctx_mutex);
   return ctx.solves;
}

extern "C" int gpu_density_matrix(int n, const double *C, const double *f,
                                  double *P)
{
   if (n <= 0) return 0;
   std::lock_guard<std::mutex> lock(ctx_mutex);
   int rc = ensure_capacity(n);
   if (rc != 0) return rc;

   const size_t nn = (size_t)n * n;
   if (cudaMemcpy(ctx.dA, C, sizeof(double) * nn,
                  cudaMemcpyHostToDevice) != cudaSuccess) return -17;
   if (cudaMemcpy(ctx.dW, f, sizeof(double) * n,
                  cudaMemcpyHostToDevice) != cudaSuccess) return -18;

   const int threads = 256;
   const int blocks = (int)((nn + threads - 1) / threads);
   scale_columns<<<blocks, threads>>>(n, ctx.dA, ctx.dW, ctx.dScaled);
   if (cudaGetLastError() != cudaSuccess) return -19;

   const double alpha = 1.0, beta = 0.0;
   if (cublasDgemm(ctx.blas, CUBLAS_OP_N, CUBLAS_OP_T, n, n, n,
                   &alpha, ctx.dScaled, n, ctx.dA, n,
                   &beta, ctx.dP, n) != CUBLAS_STATUS_SUCCESS) return -20;
   if (cudaMemcpy(P, ctx.dP, sizeof(double) * nn,
                  cudaMemcpyDeviceToHost) != cudaSuccess) return -21;
   return 0;
}

extern "C" int gpu_build_isotropic_h1(
      int n, int nmat, int nshell, const int *matlist, const double *H0,
      const double *S, const double *shift, const int *ao2sh,
      double autoev, double *H)
{
   if (n <= 0 || nmat <= 0) return 0;
   std::lock_guard<std::mutex> lock(ctx_mutex);
   int rc = ensure_capacity(n);
   if (rc != 0) return rc;
   rc = ensure_scf_capacity(nmat, nshell);
   if (rc != 0) return rc;
   const size_t nn = (size_t)n * n;
   if (cudaMemcpy(ctx.dMatlist, matlist, sizeof(int) * 2 * nmat,
                  cudaMemcpyHostToDevice) != cudaSuccess) return -26;
   if (cudaMemcpy(ctx.dH0, H0, sizeof(double) * nmat,
                  cudaMemcpyHostToDevice) != cudaSuccess) return -27;
   if (cudaMemcpy(ctx.dB, S, sizeof(double) * nn,
                  cudaMemcpyHostToDevice) != cudaSuccess) return -28;
   if (cudaMemcpy(ctx.dShift, shift, sizeof(double) * nshell,
                  cudaMemcpyHostToDevice) != cudaSuccess) return -29;
   if (cudaMemcpy(ctx.dAo2sh, ao2sh, sizeof(int) * n,
                  cudaMemcpyHostToDevice) != cudaSuccess) return -30;

   const int threads = 256;
   build_isotropic_h1_kernel<<<(nmat + threads - 1) / threads, threads>>>(
      n, nmat, ctx.dMatlist, ctx.dH0, ctx.dB, ctx.dShift, ctx.dAo2sh,
      autoev, ctx.dA);
   if (cudaGetLastError() != cudaSuccess) return -31;
   if (cudaMemcpy(H, ctx.dA, sizeof(double) * nn,
                  cudaMemcpyDeviceToHost) != cudaSuccess) return -32;
   return 0;
}

extern "C" int gpu_mulliken_shell(
      int n, int nshell, const int *ao2sh, const double *S,
      const double *P, double *qsh)
{
   if (n <= 0) return 0;
   std::lock_guard<std::mutex> lock(ctx_mutex);
   int rc = ensure_capacity(n);
   if (rc != 0) return rc;
   rc = ensure_scf_capacity(0, nshell);
   if (rc != 0) return rc;
   const size_t nn = (size_t)n * n;
   if (cudaMemcpy(ctx.dAo2sh, ao2sh, sizeof(int) * n,
                  cudaMemcpyHostToDevice) != cudaSuccess) return -33;
   if (cudaMemcpy(ctx.dB, S, sizeof(double) * nn,
                  cudaMemcpyHostToDevice) != cudaSuccess) return -34;
   if (cudaMemcpy(ctx.dP, P, sizeof(double) * nn,
                  cudaMemcpyHostToDevice) != cudaSuccess) return -35;
   if (cudaMemset(ctx.dQsh, 0, sizeof(double) * nshell) != cudaSuccess) return -36;

   const int threads = 256;
   mulliken_shell_kernel<<<(n + threads - 1) / threads, threads>>>(
      n, ctx.dAo2sh, ctx.dB, ctx.dP, ctx.dQsh);
   if (cudaGetLastError() != cudaSuccess) return -37;
   if (cudaMemcpy(qsh, ctx.dQsh, sizeof(double) * nshell,
                  cudaMemcpyDeviceToHost) != cudaSuccess) return -38;
   return 0;
}

// ============================================================================
//  Resident SCF session: keep S, H0, matlist, ao2sh resident on the device for
//  the whole SCF, and the Hamiltonian/eigenvectors (dH) + density (dP) resident
//  between solve and finish, so per iteration ONLY small vectors cross PCIe
//  (shift/focc in; emo/P/qsh out). Eliminates the ~8 n^2 host<->device copies
//  per iteration of the per-call path. P is returned because the host electro()
//  energy needs it.
// ============================================================================
namespace {
struct ScfSession {
   int n = 0, nshell = 0, nmat = 0, lwork = 0;
   bool active = false, gfn1 = false;
   cusolverDnHandle_t handle = nullptr;
   cublasHandle_t blas = nullptr;
   double *dS = nullptr, *dH0 = nullptr, *dH = nullptr, *dB = nullptr;
   double *dW = nullptr, *dP = nullptr, *dScaled = nullptr;
   double *dShift = nullptr, *dFocc = nullptr, *dQsh = nullptr, *dWork = nullptr;
   int *dMatlist = nullptr, *dAo2sh = nullptr, *dInfo = nullptr;

   void free() {
      cudaFree(dS); cudaFree(dH0); cudaFree(dH); cudaFree(dB); cudaFree(dW);
      cudaFree(dP); cudaFree(dScaled); cudaFree(dShift); cudaFree(dFocc);
      cudaFree(dQsh); cudaFree(dWork); cudaFree(dMatlist); cudaFree(dAo2sh);
      cudaFree(dInfo);
      dS=dH0=dH=dB=dW=dP=dScaled=dShift=dFocc=dQsh=dWork=nullptr;
      dMatlist=dAo2sh=nullptr; dInfo=nullptr;
      if (blas) { cublasDestroy(blas); blas=nullptr; }
      if (handle) { cusolverDnDestroy(handle); handle=nullptr; }
      n=nshell=nmat=lwork=0; active=false; gfn1=false;
   }
};
ScfSession sess;
} // namespace

// Open a session. H0/matlist/ao2sh are only needed for GFN1 (on-device H1 build);
// pass gfn1=0 and they may be null (GFN2 builds H on the host and passes it in).
extern "C" int gpu_scf_open(int n, int nshell, int nmat, int gfn1,
                            const double *H0, const double *S,
                            const int *matlist, const int *ao2sh)
{
   if (n <= 0) return -50;
   std::lock_guard<std::mutex> lock(ctx_mutex);
   sess.free();
   sess.n = n; sess.nshell = nshell; sess.nmat = nmat; sess.gfn1 = (gfn1 != 0);
   const size_t nn = (size_t)n * n;
   if (cusolverDnCreate(&sess.handle) != CUSOLVER_STATUS_SUCCESS) return -51;
   if (cublasCreate(&sess.blas) != CUBLAS_STATUS_SUCCESS) return -52;
   if (cudaMalloc((void**)&sess.dS, sizeof(double)*nn) != cudaSuccess) return -53;
   if (cudaMalloc((void**)&sess.dH, sizeof(double)*nn) != cudaSuccess) return -54;
   if (cudaMalloc((void**)&sess.dB, sizeof(double)*nn) != cudaSuccess) return -55;
   if (cudaMalloc((void**)&sess.dP, sizeof(double)*nn) != cudaSuccess) return -56;
   if (cudaMalloc((void**)&sess.dScaled, sizeof(double)*nn) != cudaSuccess) return -57;
   if (cudaMalloc((void**)&sess.dW, sizeof(double)*n) != cudaSuccess) return -58;
   if (cudaMalloc((void**)&sess.dFocc, sizeof(double)*n) != cudaSuccess) return -59;
   if (cudaMalloc((void**)&sess.dAo2sh, sizeof(int)*n) != cudaSuccess) return -60;
   if (cudaMalloc((void**)&sess.dInfo, sizeof(int)) != cudaSuccess) return -61;
   if (nshell > 0) {
      if (cudaMalloc((void**)&sess.dShift, sizeof(double)*nshell) != cudaSuccess) return -62;
      if (cudaMalloc((void**)&sess.dQsh, sizeof(double)*nshell) != cudaSuccess) return -63;
   }
   if (sess.gfn1 && nmat > 0) {
      if (cudaMalloc((void**)&sess.dH0, sizeof(double)*nmat) != cudaSuccess) return -64;
      if (cudaMalloc((void**)&sess.dMatlist, sizeof(int)*2*nmat) != cudaSuccess) return -65;
      if (cudaMemcpy(sess.dH0, H0, sizeof(double)*nmat, cudaMemcpyHostToDevice) != cudaSuccess) return -66;
      if (cudaMemcpy(sess.dMatlist, matlist, sizeof(int)*2*nmat, cudaMemcpyHostToDevice) != cudaSuccess) return -67;
   }
   // resident constants
   if (cudaMemcpy(sess.dS, S, sizeof(double)*nn, cudaMemcpyHostToDevice) != cudaSuccess) return -68;
   if (ao2sh && cudaMemcpy(sess.dAo2sh, ao2sh, sizeof(int)*n, cudaMemcpyHostToDevice) != cudaSuccess) return -69;
   if (cusolverDnDsygvd_bufferSize(sess.handle, CUSOLVER_EIG_TYPE_1,
         CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, sess.dH, n,
         sess.dB, n, sess.dW, &sess.lwork) != CUSOLVER_STATUS_SUCCESS) return -70;
   if (cudaMalloc((void**)&sess.dWork, sizeof(double)*sess.lwork) != cudaSuccess) return -71;
   sess.active = true;
   return 0;
}

// One SCF diagonalization. GFN1: build H on-device from resident H0/S/shift.
// GFN2: H_in is the host-built Hamiltonian (copied in once). Eigenvectors stay
// resident in dH; only eigenvalues (emo) come back.
extern "C" int gpu_scf_solve(int n, int nshell, const double *H_in,
                             const double *shift, double autoev, double *emo)
{
   if (!sess.active || n != sess.n) return -72;
   std::lock_guard<std::mutex> lock(ctx_mutex);
   const size_t nn = (size_t)n * n;
   const int threads = 256;
   if (sess.gfn1 && shift) {
      if (cudaMemcpy(sess.dShift, shift, sizeof(double)*nshell, cudaMemcpyHostToDevice) != cudaSuccess) return -73;
      build_isotropic_h1_kernel<<<(sess.nmat + threads - 1)/threads, threads>>>(
         n, sess.nmat, sess.dMatlist, sess.dH0, sess.dS, sess.dShift, sess.dAo2sh,
         autoev, sess.dH);
      if (cudaGetLastError() != cudaSuccess) return -74;
   } else if (H_in) {
      if (cudaMemcpy(sess.dH, H_in, sizeof(double)*nn, cudaMemcpyHostToDevice) != cudaSuccess) return -75;
   } else {
      return -76;
   }
   // sygvd destroys B -> refresh from resident dS (device-to-device, no PCIe)
   if (cudaMemcpy(sess.dB, sess.dS, sizeof(double)*nn, cudaMemcpyDeviceToDevice) != cudaSuccess) return -77;
   cusolverStatus_t st = cusolverDnDsygvd(sess.handle, CUSOLVER_EIG_TYPE_1,
      CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, sess.dH, n,
      sess.dB, n, sess.dW, sess.dWork, sess.lwork, sess.dInfo);
   if (st != CUSOLVER_STATUS_SUCCESS) return (int)st;
   int info = 0;
   if (cudaMemcpy(&info, sess.dInfo, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess) return -78;
   if (info != 0) return 2000 + info;
   if (cudaMemcpy(emo, sess.dW, sizeof(double)*n, cudaMemcpyDeviceToHost) != cudaSuccess) return -79;
   return 0;
}

// Density (from resident eigenvectors dH + occupations) and shell-Mulliken
// charges. Returns P (host needs it for electro) and qsh.
extern "C" int gpu_scf_finish(int n, int nshell, const double *focc,
                              double *P, double *qsh)
{
   if (!sess.active || n != sess.n) return -80;
   std::lock_guard<std::mutex> lock(ctx_mutex);
   const size_t nn = (size_t)n * n;
   const int threads = 256;
   if (cudaMemcpy(sess.dFocc, focc, sizeof(double)*n, cudaMemcpyHostToDevice) != cudaSuccess) return -81;
   scale_columns<<<(int)((nn + threads - 1)/threads), threads>>>(n, sess.dH, sess.dFocc, sess.dScaled);
   if (cudaGetLastError() != cudaSuccess) return -82;
   const double alpha = 1.0, beta = 0.0;
   if (cublasDgemm(sess.blas, CUBLAS_OP_N, CUBLAS_OP_T, n, n, n, &alpha,
                   sess.dScaled, n, sess.dH, n, &beta, sess.dP, n) != CUBLAS_STATUS_SUCCESS) return -83;
   if (cudaMemcpy(P, sess.dP, sizeof(double)*nn, cudaMemcpyDeviceToHost) != cudaSuccess) return -84;
   if (nshell > 0) {
      if (cudaMemset(sess.dQsh, 0, sizeof(double)*nshell) != cudaSuccess) return -85;
      mulliken_shell_kernel<<<(n + threads - 1)/threads, threads>>>(
         n, sess.dAo2sh, sess.dS, sess.dP, sess.dQsh);
      if (cudaGetLastError() != cudaSuccess) return -86;
      if (cudaMemcpy(qsh, sess.dQsh, sizeof(double)*nshell, cudaMemcpyDeviceToHost) != cudaSuccess) return -87;
   }
   return cudaDeviceSynchronize() == cudaSuccess ? 0 : -88;
}

// Copy the resident eigenvectors (from the last solve) back to the host. Called
// once after SCF convergence so the host C is available for the gradient's
// energy-weighted density and any orbital-based properties.
extern "C" int gpu_scf_get_vectors(int n, double *C)
{
   if (!sess.active || n != sess.n) return -90;
   std::lock_guard<std::mutex> lock(ctx_mutex);
   const size_t nn = (size_t)n * n;
   if (cudaMemcpy(C, sess.dH, sizeof(double)*nn, cudaMemcpyDeviceToHost) != cudaSuccess) return -91;
   return 0;
}

extern "C" int gpu_scf_close()
{
   std::lock_guard<std::mutex> lock(ctx_mutex);
   sess.free();
   return 0;
}

// ============================================================================
//  Analytical gradient (build_dSDQH0): overlap/dipole/quadrupole integral-
//  derivative contraction with the density and energy-weighted density, ported
//  faithfully from src/intgrad.f90 + src/xtb/hamiltonian.F90 build_dSDQH0.
//  One CUDA thread per atom pair (iat>jat); forces accumulate via atomicAdd.
//  GFN1/2 basis is s/p/d (l<=2); the gradient raises one index to f (l<=3).
// ============================================================================
namespace gradk {

__device__ __constant__ double EVTOAU = 1.0 / 27.21138505;  // overwritten at launch
// itt, llao, llao2 small tables
__device__ __forceinline__ int itt(int l){ const int t[4]={0,1,4,10}; return t[l]; }
__device__ __forceinline__ int llao(int l){ const int t[4]={1,3,6,10}; return t[l]; }
__device__ __forceinline__ int llao2(int l){ const int t[4]={1,3,5,7}; return t[l]; }

// Cartesian powers lx,ly,lz for AO components 1..10 (s,p,d), 0-based columns 0..9
__device__ __forceinline__ void lxyz_of(int col, int &lx, int &ly, int &lz){
   const int LX[10]={0, 1,0,0, 2,0,0,1,1,0};
   const int LY[10]={0, 0,1,0, 0,2,0,1,0,1};
   const int LZ[10]={0, 0,0,1, 0,0,2,0,1,1};
   lx=LX[col]; ly=LY[col]; lz=LZ[col];
}

// olapp: partial 1D overlap (intgrad.f90:95)
__device__ __forceinline__ double olapp(int l, double gama){
   if (l & 1) return 0.0;
   const double dftr[8]={1.,1.,3.,15.,105.,945.,10395.,135135.};
   int lh=l/2; double gm=0.5/gama, s=1.0;
   for(int i=0;i<lh;i++) s*=gm;
   return s*dftr[lh];
}

// horizontal_shift (intgrad.f90:357); c[] is 0-based by power
__device__ __forceinline__ void hshift(double ae, int l, double *c){
   if (l==1){ c[0]+=ae*c[1]; }
   else if (l==2){ c[0]+=ae*ae*c[2]; c[1]+=2.0*ae*c[2]; }
   else if (l==3){ c[0]+=ae*ae*ae*c[3]; c[1]+=3.0*ae*ae*c[3]; c[2]+=3.0*ae*c[3]; }
}

// form_product (intgrad.f90:382), 0-based; la,lb<=3 for GFN gradient
__device__ __forceinline__ void form_product(const double*a,const double*b,int la,int lb,double*d){
   if (la>=3||lb>=3){
      d[0]=a[0]*b[0]; d[1]=a[0]*b[1]+a[1]*b[0]; d[2]=a[0]*b[2]+a[2]*b[0]; d[3]=a[0]*b[3]+a[3]*b[0];
      if (la==0||lb==0) return;
      d[2]+=a[1]*b[1]; d[3]+=a[1]*b[2]+a[2]*b[1]; d[4]=a[1]*b[3]+a[3]*b[1];
      if (la<=1||lb<=1) return;
      d[4]+=a[2]*b[2]; d[5]=a[2]*b[3]+a[3]*b[2];
      if (la<=2||lb<=2) return;
      d[6]=a[3]*b[3];
      return;
   }
   if (la>=2||lb>=2){
      d[0]=a[0]*b[0]; d[1]=a[0]*b[1]+a[1]*b[0]; d[2]=a[0]*b[2]+a[2]*b[0];
      if (la==0||lb==0) return;
      d[2]+=a[1]*b[1]; d[3]=a[1]*b[2]+a[2]*b[1];
      if (la<=1||lb<=1) return;
      d[4]=a[2]*b[2];
      return;
   }
   d[0]=a[0]*b[0];
   if (la==0&&lb==0) return;
   d[1]=a[0]*b[1]+a[1]*b[0];
   if (la==0||lb==0) return;
   d[2]=a[1]*b[1];
}

// multipole_grad_3d (intgrad.f90:607) -> s3d[10], ds3d[3][10]
__device__ void multipole_grad_3d(const double*ri,const double*rj,const double*rc,
      const double*rp,double ai,double aj,const int*li,const int*lj,
      const double*s1d,double*s3d,double ds3d[3][10]){
   double val[3][3]={{0}}, gra[3][3]={{0}};
   for(int k=0;k<3;k++){
      double vv[8]={0}, gg[8]={0}, vi[7]={0}, vj[7]={0}, gi[7]={0};
      double rpc=rp[k]-rc[k];
      vi[li[k]]=1.0; vj[lj[k]]=1.0; gi[li[k]+1]=2.0*ai;
      if (li[k]>0) gi[li[k]-1]=-(double)li[k];
      hshift(rp[k]-ri[k], li[k]-1, gi);
      hshift(rp[k]-ri[k], li[k]+1, gi);
      hshift(rp[k]-ri[k], li[k],   vi);
      hshift(rp[k]-rj[k], lj[k],   vj);
      form_product(vi,vj,li[k],lj[k],vv);
      form_product(gi,vj,li[k]+1,lj[k],gg);
      for(int l=0;l<=li[k]+lj[k]+1;l++){
         double a0=s1d[l], a1=s1d[l+1]+rpc*s1d[l], a2=s1d[l+2]+2*rpc*s1d[l+1]+rpc*rpc*s1d[l];
         val[k][0]+=a0*vv[l]; val[k][1]+=a1*vv[l]; val[k][2]+=a2*vv[l];
         gra[k][0]+=a0*gg[l]; gra[k][1]+=a1*gg[l]; gra[k][2]+=a2*gg[l];
      }
   }
   // s3d (1..10 -> 0..9)
   s3d[0]=val[0][0]*val[1][0]*val[2][0];
   s3d[1]=val[0][1]*val[1][0]*val[2][0];
   s3d[2]=val[0][0]*val[1][1]*val[2][0];
   s3d[3]=val[0][0]*val[1][0]*val[2][1];
   s3d[4]=val[0][2]*val[1][0]*val[2][0];
   s3d[5]=val[0][0]*val[1][2]*val[2][0];
   s3d[6]=val[0][0]*val[1][0]*val[2][2];
   s3d[7]=val[0][1]*val[1][1]*val[2][0];
   s3d[8]=val[0][1]*val[1][0]*val[2][1];
   s3d[9]=val[0][0]*val[1][1]*val[2][1];
   ds3d[0][0]=gra[0][0]*val[1][0]*val[2][0]; ds3d[1][0]=val[0][0]*gra[1][0]*val[2][0]; ds3d[2][0]=val[0][0]*val[1][0]*gra[2][0];
   ds3d[0][1]=gra[0][1]*val[1][0]*val[2][0]; ds3d[1][1]=val[0][1]*gra[1][0]*val[2][0]; ds3d[2][1]=val[0][1]*val[1][0]*gra[2][0];
   ds3d[0][2]=gra[0][0]*val[1][1]*val[2][0]; ds3d[1][2]=val[0][0]*gra[1][1]*val[2][0]; ds3d[2][2]=val[0][0]*val[1][1]*gra[2][0];
   ds3d[0][3]=gra[0][0]*val[1][0]*val[2][1]; ds3d[1][3]=val[0][0]*gra[1][0]*val[2][1]; ds3d[2][3]=val[0][0]*val[1][0]*gra[2][1];
   ds3d[0][4]=gra[0][2]*val[1][0]*val[2][0]; ds3d[1][4]=val[0][2]*gra[1][0]*val[2][0]; ds3d[2][4]=val[0][2]*val[1][0]*gra[2][0];
   ds3d[0][5]=gra[0][0]*val[1][2]*val[2][0]; ds3d[1][5]=val[0][0]*gra[1][2]*val[2][0]; ds3d[2][5]=val[0][0]*val[1][2]*gra[2][0];
   ds3d[0][6]=gra[0][0]*val[1][0]*val[2][2]; ds3d[1][6]=val[0][0]*gra[1][0]*val[2][2]; ds3d[2][6]=val[0][0]*val[1][0]*gra[2][2];
   ds3d[0][7]=gra[0][1]*val[1][1]*val[2][0]; ds3d[1][7]=val[0][1]*gra[1][1]*val[2][0]; ds3d[2][7]=val[0][1]*val[1][1]*gra[2][0];
   ds3d[0][8]=gra[0][1]*val[1][0]*val[2][1]; ds3d[1][8]=val[0][1]*gra[1][0]*val[2][1]; ds3d[2][8]=val[0][1]*val[1][0]*gra[2][1];
   ds3d[0][9]=gra[0][0]*val[1][1]*val[2][1]; ds3d[1][9]=val[0][0]*gra[1][1]*val[2][1]; ds3d[2][9]=val[0][0]*val[1][1]*gra[2][1];
}

// shiftintg (intgrad.f90:703): g[3][19], s[10], r[3]  (indices 0-based: 11..19 -> 10..18)
__device__ void shiftintg(double g[3][19], const double*s, const double*r){
   for(int x=0;x<3;x++){
      g[x][10]=g[x][1]-r[0]*g[x][0];
      g[x][11]=g[x][2]-r[1]*g[x][0];
      g[x][12]=g[x][3]-r[2]*g[x][0];
      g[x][13]=g[x][4]-2*r[0]*g[x][1]+r[0]*r[0]*g[x][0];
      g[x][14]=g[x][5]-2*r[1]*g[x][2]+r[1]*r[1]*g[x][0];
      g[x][15]=g[x][6]-2*r[2]*g[x][3]+r[2]*r[2]*g[x][0];
      g[x][16]=g[x][7]-r[0]*g[x][2]-r[1]*g[x][1]+r[0]*r[1]*g[x][0];
      g[x][17]=g[x][8]-r[0]*g[x][3]-r[2]*g[x][1]+r[0]*r[2]*g[x][0];
      g[x][18]=g[x][9]-r[1]*g[x][3]-r[2]*g[x][2]+r[1]*r[2]*g[x][0];
   }
   g[0][10]-=s[0]; g[1][11]-=s[0]; g[2][12]-=s[0];
   g[0][13]+=-2*s[1]+2*r[0]*s[0];
   g[1][14]+=-2*s[2]+2*r[1]*s[0];
   g[2][15]+=-2*s[3]+2*r[2]*s[0];
   g[0][16]+=-s[2]+r[1]*s[0];
   g[1][16]+=-s[1]+r[0]*s[0];
   g[0][17]+=-s[3]+r[2]*s[0];
   g[2][17]+=-s[1]+r[0]*s[0];
   g[1][18]+=-s[3]+r[2]*s[0];
   g[2][18]+=-s[2]+r[1]*s[0];
}

// trafo matrix (CAO->SAO), trafo[m][jj]
__device__ __forceinline__ double trafo(int m,int jj){
   const double s5=0.4472135954999579, h3=0.8660254037844386; // sqrt(1/5), 0.5*sqrt(3)
   // columns: 0:[s5,s5,s5,0,0,0] 1:[h3,-h3,0,0,0,0] 2:[.5,.5,-1,0,0,0] 3:e4 4:e5 5:e6
   if (jj==0) return (m<3)? s5 : 0.0;
   if (jj==1) return (m==0)? h3 : (m==1? -h3 : 0.0);
   if (jj==2) return (m==0||m==1)? 0.5 : (m==2? -1.0 : 0.0);
   return (m==jj+0 && m>=3 && (m-3)==(jj-3))? 1.0 : (m==jj? 1.0:0.0);
}

// dtrf2 (intgrad.f90:144) in place on s[6][6]; li,lj in {0,1,2}
__device__ void dtrf2(double s[6][6], int li, int lj){
   if (li<2 && lj<2) return;
   double s2[6][6];
   if (li==0){
      for(int jj=0;jj<6;jj++){ double v=0; for(int m=0;m<6;m++) v+=trafo(m,jj)*s[m][0]; s2[jj][0]=v; }
      for(int r=0;r<5;r++) s[r][0]=s2[r+1][0];
      return;
   }
   if (li==1){
      for(int ii=0;ii<3;ii++){ for(int jj=0;jj<6;jj++){ double v=0; for(int m=0;m<6;m++) v+=trafo(m,jj)*s[m][ii]; s2[jj][ii]=v; }
         for(int r=0;r<5;r++) s[r][ii]=s2[r+1][ii]; }
      return;
   }
   if (lj==0){
      for(int jj=0;jj<6;jj++){ double v=0; for(int m=0;m<6;m++) v+=trafo(m,jj)*s[0][m]; s2[0][jj]=v; }
      for(int c=0;c<5;c++) s[0][c]=s2[0][c+1];
      return;
   }
   if (lj==1){
      for(int ii=0;ii<3;ii++){ for(int jj=0;jj<6;jj++){ double v=0; for(int m=0;m<6;m++) v+=trafo(m,jj)*s[ii][m]; s2[ii][jj]=v; }
         for(int c=0;c<5;c++) s[ii][c]=s2[ii][c+1]; }
      return;
   }
   // d-d: dum = trafo^T * s ; s2 = dum * trafo
   double dum[6][6];
   for(int a=0;a<6;a++) for(int b=0;b<6;b++){ double v=0; for(int m=0;m<6;m++) v+=trafo(m,a)*s[m][b]; dum[a][b]=v; }
   for(int a=0;a<6;a++) for(int b=0;b<6;b++){ double v=0; for(int m=0;m<6;m++) v+=dum[a][m]*trafo(m,b); s2[a][b]=v; }
   for(int a=0;a<5;a++) for(int b=0;b<5;b++) s[a][b]=s2[a+1][b+1];
}

// dshellPoly (grad_core.f90:40)
__device__ void dshellPoly(double iPoly,double jPoly,double iRad,double jRad,
      double rab2,const double*x1,const double*x2,double &rf,double*dxyz){
   double dx=x1[0]-x2[0], dy=x1[1]-x2[1], dz=x1[2]-x2[2];
   double rab=sqrt(rab2), r=iRad+jRad, rr=rab/r;
   double k1=iPoly*0.01, k2=jPoly*0.01;
   double t14=sqrt(rr), t15=k1*t14, t17=1.0/rab2, t23=k2*t14;
   rf=(1.0+t15)*(1.0+k2*t14);
   double pref=(t15*(1.0+t23)+(1.0+t15)*k2*t14)*0.5*t17;
   dxyz[0]=pref*dx; dxyz[1]=pref*dy; dxyz[2]=pref*dz;
}

} // namespace gradk

// flattened data for the kernel (column-major Fortran arrays passed straight in)
struct GradData {
   int nat, nao, maxsh, nelem, ntrans;   // maxsh = hData shell leading dim
   int ldSE, ldcao, ldsp;                // selfEnergy/dSEdcn/ves, caoshell/saoshell, shellPoly
   const int *nShell, *at, *angShell, *valShell, *caoshell, *saoshell, *nprim, *primcount;
   const double *xyz, *trans, *slaterExp, *shellPoly, *atomicRad, *en, *enScale,
                *kScale, *pairParam, *selfEnergy, *dSEdcn, *alp, *cont,
                *P, *Pew, *ves, *vs, *vd, *vq;
   double *gout, *sigout, *dhdcn;   // accumulated outputs (seeded with prior contributions)
   double enScale4, kDiff, wExp, intcut, evtoau;
};

// index helpers for column-major Fortran arrays
#define SH2(a,ish,iat,ld) ((a)[(ish)+(size_t)(ld)*(iat)])     // [ld, *]
__device__ __forceinline__ double h0scal_dev(const GradData d,int il,int jl,int izp,int jzp,
      bool vi,bool vj){
   // il,jl 1-based (ishtyp+1); izp,jzp 1-based atomic numbers; enScale/kScale [4,4]
   if (vi&&vj){
      double den=(d.en[izp-1]-d.en[jzp-1]); den*=den;
      double ens=d.enScale[(jl-1)+4*(il-1)];
      double enpoly=1.0+ens*den*(1.0+d.enScale4*den);
      return d.kScale[(jl-1)+4*(il-1)]*enpoly*d.pairParam[(izp-1)+(size_t)d.nelem*(jzp-1)];
   }
   if (!vi&&!vj) return d.kDiff;
   if (!vi&&vj)  return 0.5*(d.kScale[(jl-1)+4*(jl-1)]+d.kDiff);
   return 0.5*(d.kScale[(il-1)+4*(il-1)]+d.kDiff);  // !vj && vi
}

__global__ void grad_pair_kernel(GradData d)
{
   int iat=blockIdx.x*blockDim.x+threadIdx.x + 1;   // 1-based
   int jat=blockIdx.y*blockDim.y+threadIdx.y + 1;
   if (iat>d.nat || jat>=iat) return;
   const double evtoau=d.evtoau;
   int izp=d.at[iat-1], jzp=d.at[jat-1];
   double ri[3]={d.xyz[0+3*(iat-1)],d.xyz[1+3*(iat-1)],d.xyz[2+3*(iat-1)]};
   double gx[3]={0,0,0}, sig[9]={0,0,0,0,0,0,0,0,0}, dhi=0.0, dhj=0.0;
   int nshi=d.nShell[izp-1], nshj=d.nShell[jzp-1];
   for(int ish=1; ish<=nshi; ish++){
      int ishtyp=SH2(d.angShell,ish-1,izp-1,d.maxsh);
      int icao=SH2(d.caoshell,ish-1,iat-1,d.ldcao);
      int naoi=gradk::llao(ishtyp), iptyp=gradk::itt(ishtyp);
      double hii=SH2(d.selfEnergy,ish-1,iat-1,d.ldSE);
      double zi=SH2(d.slaterExp,ish-1,izp-1,d.maxsh);
      bool vali=SH2(d.valShell,ish-1,izp-1,d.maxsh)!=0;
      for(int jsh=1; jsh<=nshj; jsh++){
         int jshtyp=SH2(d.angShell,jsh-1,jzp-1,d.maxsh);
         int jcao=SH2(d.caoshell,jsh-1,jat-1,d.ldcao);
         int naoj=gradk::llao(jshtyp), jptyp=gradk::itt(jshtyp);
         int il=ishtyp+1, jl=jshtyp+1;
         double hjj=SH2(d.selfEnergy,jsh-1,jat-1,d.ldSE);
         double zj=SH2(d.slaterExp,jsh-1,jzp-1,d.maxsh);
         bool valj=SH2(d.valShell,jsh-1,jzp-1,d.maxsh)!=0;
         double zetaij=pow(2.0*sqrt(zi*zj)/(zi+zj), d.wExp);
         double km=h0scal_dev(d,il,jl,izp,jzp,vali,valj);
         double hav=0.5*km*(hii+hjj)*zetaij*evtoau;
         for(int itr=0; itr<d.ntrans; itr++){
            double rj[3]={d.xyz[0+3*(jat-1)]+d.trans[0+3*itr],
                          d.xyz[1+3*(jat-1)]+d.trans[1+3*itr],
                          d.xyz[2+3*(jat-1)]+d.trans[2+3*itr]};
            double rij[3]={ri[0]-rj[0],ri[1]-rj[1],ri[2]-rj[2]};
            double rij2=rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
            if (rij2>1600.0) continue;
            double shpoly, dshpoly[3];
            gradk::dshellPoly(SH2(d.shellPoly,il-1,izp-1,d.ldsp),SH2(d.shellPoly,jl-1,jzp-1,d.ldsp),
               d.atomicRad[izp-1],d.atomicRad[jzp-1],rij2,ri,rj,shpoly,dshpoly);
            // ---- integrals: get_grad_multiint ----
            double sdq[10][6][6]; double sdqg[6][6][3][19];
            // only the [naoj][naoi] (Cartesian) block is ever touched
            for(int a=0;a<10;a++)for(int b=0;b<naoj;b++)for(int c=0;c<naoi;c++) sdq[a][b][c]=0.0;
            for(int b=0;b<naoj;b++)for(int c=0;c<naoi;c++)for(int x=0;x<3;x++)for(int k=0;k<19;k++) sdqg[b][c][x][k]=0.0;
            int npi=d.nprim[icao], npj=d.nprim[jcao];
            for(int ip=0; ip<npi; ip++){
               double alpi=d.alp[ip+d.primcount[icao]];
               for(int jp=0; jp<npj; jp++){
                  double alpj=d.alp[jp+d.primcount[jcao]];
                  double est=alpi*alpj*rij2;
                  if (est > d.intcut*(alpi+alpj)) continue;
                  double ab=1.0/(alpi+alpj); est*=ab;
                  double kab=exp(-est)*pow(sqrt(M_PI)*sqrt(ab),3.0);
                  double rp[3]={(alpi*ri[0]+alpj*rj[0])*ab,(alpi*ri[1]+alpj*rj[1])*ab,(alpi*ri[2]+alpj*rj[2])*ab};
                  double t[12]; for(int k=0;k<=ishtyp+jshtyp+3;k++) t[k]=gradk::olapp(k,alpi+alpj);
                  for(int mli=1;mli<=naoi;mli++){
                     double ci=d.cont[ip+d.primcount[icao+mli-1]];
                     int liv[3]; gradk::lxyz_of(iptyp+mli-1,liv[0],liv[1],liv[2]);
                     for(int mlj=1;mlj<=naoj;mlj++){
                        double cc=kab*d.cont[jp+d.primcount[jcao+mlj-1]]*ci;
                        int ljv[3]; gradk::lxyz_of(jptyp+mlj-1,ljv[0],ljv[1],ljv[2]);
                        double saw[10], sawg[3][10];
                        gradk::multipole_grad_3d(ri,rj,rj,rp,alpi,alpj,liv,ljv,t,saw,sawg);
                        for(int k=0;k<10;k++) sdq[k][mlj-1][mli-1]+=saw[k]*cc;
                        for(int x=0;x<3;x++)for(int k=0;k<10;k++) sdqg[mlj-1][mli-1][x][k]+=sawg[x][k]*cc;
                     }
                  }
               }
            }
            // shiftintg per (mlj,mli): fills sdqg components 11..19 (0-based 10..18)
            for(int mli=0;mli<naoi;mli++)for(int mlj=0;mlj<naoj;mlj++){
               double s10[10]; for(int k=0;k<10;k++) s10[k]=sdq[k][mlj][mli];
               gradk::shiftintg(sdqg[mlj][mli], s10, rij);
            }
            // dtrf2 (CAO->SAO) on the overlap block and on each sdqg(ixyz,k,:,:).
            // CPU layout sdq(1,jj,ii) -> here sdq[0][jj][ii]; build s[i=ii][j=jj].
            // dtrf2 (CAO->SAO) only does work when a d-shell is involved; for
            // s/p it is a no-op, so skip the 57 build-transform-store passes
            // entirely (bit-exact, big saving for C/N/O/H). dtrf2 expects
            // s[mlj][mli], which is how sdq/sdqg are stored (no transpose).
            if (ishtyp>=2 || jshtyp>=2){
               { double s[6][6];
                 for(int i=0;i<6;i++)for(int j=0;j<6;j++) s[i][j]=sdq[0][i][j];
                 gradk::dtrf2(s,ishtyp,jshtyp);
                 for(int i=0;i<6;i++)for(int j=0;j<6;j++) sdq[0][i][j]=s[i][j];
               }
               for(int k=0;k<19;k++)for(int x=0;x<3;x++){
                  double s[6][6];
                  for(int i=0;i<6;i++)for(int j=0;j<6;j++) s[i][j]=sdqg[i][j][x][k];
                  gradk::dtrf2(s,ishtyp,jshtyp);
                  for(int i=0;i<6;i++)for(int j=0;j<6;j++) sdqg[i][j][x][k]=s[i][j];
               }
            }
            // ---- contraction (build_dSDQH0 lines 509-546) ----
            int li2=gradk::llao2(ishtyp), lj2=gradk::llao2(jshtyp);
            double vesi=d.ves[(ish-1)+(size_t)d.ldSE*(iat-1)];
            double vesj=d.ves[(jsh-1)+(size_t)d.ldSE*(jat-1)];
            double dSi=SH2(d.dSEdcn,ish-1,iat-1,d.ldSE);
            double dSj=SH2(d.dSEdcn,jsh-1,jat-1,d.ldSE);
            double gxyz[3]={0,0,0};
            for(int ii=1;ii<=li2;ii++){
               int iao=ii+SH2(d.saoshell,ish-1,iat-1,d.ldcao);  // 1-based SAO
               for(int jj=1;jj<=lj2;jj++){
                  int jao=jj+SH2(d.saoshell,jsh-1,jat-1,d.ldcao);
                  double Pij=d.P[(jao-1)+(size_t)d.nao*(iao-1)];
                  double HPij=hav*shpoly*Pij;
                  double ov=sdq[0][jj-1][ii-1];
                  double Pewji=d.Pew[(jao-1)+(size_t)d.nao*(iao-1)];
                  for(int x=0;x<3;x++) gxyz[x]+=2*HPij*ov*dshpoly[x]/shpoly;
                  for(int x=0;x<3;x++){
                     double sg1=sdqg[jj-1][ii-1][x][0];
                     double stmp=sg1*(2*HPij-2*Pewji - Pij*(vesi+vesj) + Pij*(d.vs[iat-1]+d.vs[jat-1]));
                     double dsum=0,qsum=0;
                     for(int m=0;m<3;m++) dsum+=sdqg[jj-1][ii-1][x][10+m]*d.vd[m+3*(iat-1)]
                                              + sdqg[jj-1][ii-1][x][1+m]*d.vd[m+3*(jat-1)];
                     for(int m=0;m<6;m++) qsum+=sdqg[jj-1][ii-1][x][13+m]*d.vq[m+6*(iat-1)]
                                              + sdqg[jj-1][ii-1][x][4+m]*d.vq[m+6*(jat-1)];
                     gxyz[x]+=stmp+Pij*dsum+Pij*qsum;
                  }
                  double HP2=km*zetaij*shpoly*Pij*ov*evtoau;
                  dhi+=HP2*dSi; dhj+=HP2*dSj;
               }
            }
            for(int x=0;x<3;x++) gx[x]+=gxyz[x];
            // sigma(row,col) += gxyz(col)*rij(row)  (column-major)
            for(int col=0;col<3;col++)for(int row=0;row<3;row++) sig[row+3*col]+=gxyz[col]*rij[row];
         } // itr
      } // jsh
   } // ish
   // accumulate into the shared outputs once per pair
   for(int x=0;x<3;x++){ atomicAdd(&d.gout[x+3*(iat-1)],  gx[x]);
                         atomicAdd(&d.gout[x+3*(jat-1)], -gx[x]); }
   for(int k=0;k<9;k++)  atomicAdd(&d.sigout[k], sig[k]);
   atomicAdd(&d.dhdcn[iat-1], dhi);
   atomicAdd(&d.dhdcn[jat-1], dhj);
}

// diagonal dE/dCN term (build_dSDQH0 lines 557-569): one thread per atom
__global__ void grad_diag_kernel(GradData d)
{
   int iat=blockIdx.x*blockDim.x+threadIdx.x + 1;
   if (iat>d.nat) return;
   int izp=d.at[iat-1];
   double acc=0.0;
   for(int ish=1; ish<=d.nShell[izp-1]; ish++){
      int ishtyp=SH2(d.angShell,ish-1,izp-1,d.maxsh);
      double dS=SH2(d.dSEdcn,ish-1,iat-1,d.ldSE);
      int sao=SH2(d.saoshell,ish-1,iat-1,d.ldcao);
      for(int iao=1; iao<=gradk::llao2(ishtyp); iao++){
         int i=iao+sao;                       // 1-based SAO
         double Pii=d.P[(i-1)+(size_t)d.nao*(i-1)];
         acc+=Pii*dS*d.evtoau;
      }
   }
   atomicAdd(&d.dhdcn[iat-1], acc);
}

// ---- device-side malloc/copy helpers ----
template<typename T> static T* devcopy(const T* host, size_t n){
   T* dev=nullptr; if (n==0) return dev;
   if (cudaMalloc((void**)&dev, n*sizeof(T))!=cudaSuccess) return nullptr;
   cudaMemcpy(dev, host, n*sizeof(T), cudaMemcpyHostToDevice); return dev;
}

extern "C" int gpu_build_dsdqh0(
   int nat,int nao,int nbf,int maxsh,int nelem,int nprimtot,int ntrans,
   int ldSE,int ldcao,int ldsp,
   const int* nShell,const int* at,const double* xyz,const double* trans,
   const int* angShell,const int* valShell,const double* slaterExp,
   const double* shellPoly,const double* atomicRad,const double* en,
   const double* enScale,double enScale4,double kDiff,double wExp,
   const double* kScale,const double* pairParam,
   const double* selfEnergy,const double* dSEdcn,
   const int* caoshell,const int* saoshell,const int* nprim,const int* primcount,
   const double* alp,const double* cont,double intcut,double evtoau,
   const double* P,const double* Pew,const double* ves,
   const double* vs,const double* vd,const double* vq,
   double* g,double* sigma,double* dhdcn)
{
   std::lock_guard<std::mutex> lock(ctx_mutex);
   using clk=std::chrono::high_resolution_clock;
   const bool tprof = getenv("XTB_GRAD_TRACE")!=nullptr;
   auto t0=clk::now();
   GradData d; memset(&d,0,sizeof(d));
   d.nat=nat; d.nao=nao; d.maxsh=maxsh; d.nelem=nelem; d.ntrans=ntrans;
   d.ldSE=ldSE; d.ldcao=ldcao; d.ldsp=ldsp;
   d.enScale4=enScale4; d.kDiff=kDiff; d.wExp=wExp; d.intcut=intcut; d.evtoau=evtoau;
   d.nShell=devcopy(nShell,nelem); d.at=devcopy(at,nat);
   d.xyz=devcopy(xyz,(size_t)3*nat); d.trans=devcopy(trans,(size_t)3*ntrans);
   d.angShell=devcopy(angShell,(size_t)maxsh*nelem); d.valShell=devcopy(valShell,(size_t)maxsh*nelem);
   d.slaterExp=devcopy(slaterExp,(size_t)maxsh*nelem); d.shellPoly=devcopy(shellPoly,(size_t)ldsp*nelem);
   d.atomicRad=devcopy(atomicRad,nelem); d.en=devcopy(en,nelem);
   d.enScale=devcopy(enScale,16); d.kScale=devcopy(kScale,16);
   d.pairParam=devcopy(pairParam,(size_t)nelem*nelem);
   d.selfEnergy=devcopy(selfEnergy,(size_t)ldSE*nat); d.dSEdcn=devcopy(dSEdcn,(size_t)ldSE*nat);
   d.caoshell=devcopy(caoshell,(size_t)ldcao*nat); d.saoshell=devcopy(saoshell,(size_t)ldcao*nat);
   d.nprim=devcopy(nprim,nbf); d.primcount=devcopy(primcount,nbf);
   d.alp=devcopy(alp,nprimtot); d.cont=devcopy(cont,nprimtot);
   d.P=devcopy(P,(size_t)nao*nao); d.Pew=devcopy(Pew,(size_t)nao*nao);
   d.ves=devcopy(ves,(size_t)ldSE*nat);
   d.vs=devcopy(vs,nat); d.vd=devcopy(vd,(size_t)3*nat); d.vq=devcopy(vq,(size_t)6*nat);
   // accumulators seeded with the incoming (non-GPU) contributions
   d.gout=devcopy(g,(size_t)3*nat); d.sigout=devcopy(sigma,9); d.dhdcn=devcopy(dhdcn,nat);
   if (!d.nShell||!d.at||!d.xyz||!d.angShell||!d.alp||!d.cont||!d.P||!d.Pew||!d.gout||!d.dhdcn)
      { if(tprof) fprintf(stderr,"[gpu-grad] alloc failed -> CPU fallback\n"); return -1; }
   auto t1=clk::now();

   dim3 blk(16,16), grd((nat+15)/16,(nat+15)/16);
   grad_pair_kernel<<<grd,blk>>>(d);
   int tb=128; grad_diag_kernel<<<(nat+tb-1)/tb,tb>>>(d);
   cudaError_t kerr=cudaDeviceSynchronize();
   if (kerr!=cudaSuccess){ if(tprof) fprintf(stderr,"[gpu-grad] kernel error: %s -> CPU fallback\n",cudaGetErrorString(kerr)); return -2; }
   auto t2=clk::now();

   cudaMemcpy(g,    d.gout,  sizeof(double)*3*nat, cudaMemcpyDeviceToHost);
   cudaMemcpy(sigma,d.sigout,sizeof(double)*9,     cudaMemcpyDeviceToHost);
   cudaMemcpy(dhdcn,d.dhdcn, sizeof(double)*nat,   cudaMemcpyDeviceToHost);

   void* ptrs[]={(void*)d.nShell,(void*)d.at,(void*)d.xyz,(void*)d.trans,(void*)d.angShell,
      (void*)d.valShell,(void*)d.slaterExp,(void*)d.shellPoly,(void*)d.atomicRad,(void*)d.en,
      (void*)d.enScale,(void*)d.kScale,(void*)d.pairParam,(void*)d.selfEnergy,(void*)d.dSEdcn,
      (void*)d.caoshell,(void*)d.saoshell,(void*)d.nprim,(void*)d.primcount,(void*)d.alp,
      (void*)d.cont,(void*)d.P,(void*)d.Pew,(void*)d.ves,(void*)d.vs,(void*)d.vd,(void*)d.vq,
      (void*)d.gout,(void*)d.sigout,(void*)d.dhdcn};
   for(void* p : ptrs) cudaFree(p);
   if (tprof){
      auto t3=clk::now();
      auto ms=[](clk::time_point a,clk::time_point b){
         return std::chrono::duration<double,std::milli>(b-a).count(); };
      fprintf(stderr,"[gpu-grad] nat=%d nao=%d | upload+malloc %.1fms | kernel %.1fms | download+free %.1fms\n",
              nat,nao,ms(t0,t1),ms(t1,t2),ms(t2,t3));
   }
   return 0;
}

// ============================================================================
//  GFN2 AES (anisotropic electrostatics) on GPU
// ============================================================================
namespace aesk {
// 1-based (l1,l2) -> packed lower-triangle index 1..6 (matches Fortran lin)
__device__ __forceinline__ int lin(int a,int b){
   return (a>=b) ? b + a*(a-1)/2 : a + b*(b-1)/2;
}
}

// setvsdq: AES potentials vs/vd/vq from q/dipm/qp via gab3/gab5 + CT correction.
// One thread per atom i; each writes its own vs(i)/vd(:,i)/vq(:,i) (no atomics).
// Faithful port of aespot.F90:setvsdq.
__global__ void aes_setvsdq_kernel(int nat, const int* at, const double* xyz,
      const double* q, const double* dipm, const double* qp,
      const double* gab3, const double* gab5,
      const double* dipKernel, const double* quadKernel,
      double* vs, double* vd, double* vq)
{
   int i = blockIdx.x*blockDim.x+threadIdx.x;
   if (i>=nat) return;
   double ra[3]={xyz[0+3*i],xyz[1+3*i],xyz[2+3*i]};
   double stmp=0, dtmp[3]={0,0,0}, qtmp[6]={0,0,0,0,0,0};
   for(int j=0;j<nat;j++){
      double g3=gab3[j+(size_t)nat*i], g5=gab5[j+(size_t)nat*i];
      double dra[3]={ra[0]-xyz[0+3*j],ra[1]-xyz[1+3*j],ra[2]-xyz[2+3*j]};
      double dum5a=0,r2a=0,r2ab=0,t1a=0,t2a=0,t3a=0;
      for(int l1=1;l1<=3;l1++){
         double ral1=ra[l1-1], dral1=dra[l1-1];
         r2a += ral1*ral1; r2ab += dral1*dral1; t1a += ral1*dral1;
         t2a += dipm[(l1-1)+3*j]*dral1; t3a += ral1*dipm[(l1-1)+3*j];
         for(int l2=1;l2<=3;l2++){
            int ll=aesk::lin(l1,l2);
            dum5a += -qp[(ll-1)+6*j]*dral1*dra[l2-1]
                     -1.5*q[j]*dral1*dra[l2-1]*ral1*ra[l2-1];
            if (l2>=l1) continue;
            qtmp[(l1+l2+1)-1] += -3.0*q[j]*g5*dra[l2-1]*dral1;
         }
         qtmp[l1-1] += -1.5*q[j]*g5*dral1*dral1;
      }
      double dum3a = -t1a*q[j]-t2a;
      dum5a += t3a*r2ab - 3.0*t1a*t2a + 0.5*q[j]*r2a*r2ab;
      stmp += dum5a*g5 + dum3a*g3;
      for(int l1=1;l1<=3;l1++){
         double dral1=dra[l1-1];
         double d3 = dral1*q[j];
         double d5 = 3.0*dral1*t2a - r2ab*dipm[(l1-1)+3*j] - q[j]*r2ab*ra[l1-1]
                   + 3.0*q[j]*dral1*t1a;
         dtmp[l1-1] += d3*g3 + d5*g5;
         qtmp[l1-1] += 0.5*r2ab*q[j]*g5;
      }
   }
   vs[i]=stmp; for(int k=0;k<3;k++) vd[k+3*i]=dtmp[k];
   for(int k=0;k<6;k++) vq[k+6*i]=qtmp[k];
   // CT correction (per atom)
   double qs1=dipKernel[at[i]-1]*2.0, qs2=quadKernel[at[i]-1]*6.0, ct3=0, ct2=0;
   for(int l1=1;l1<=3;l1++){
      ct3 += ra[l1-1]*dipm[(l1-1)+3*i]*qs1;
      vd[(l1-1)+3*i] -= qs1*dipm[(l1-1)+3*i];
      for(int l2=1;l2<l1;l2++){
         int ll=aesk::lin(l1,l2);
         vq[(l1+l2+1)-1+6*i] -= qp[(ll-1)+6*i]*qs2;
         ct3 -= ra[l1-1]*ra[l2-1]*qp[(ll-1)+6*i]*qs2;
         vd[(l1-1)+3*i] += ra[l2-1]*qp[(ll-1)+6*i]*qs2;
         vd[(l2-1)+3*i] += ra[l1-1]*qp[(ll-1)+6*i]*qs2;
      }
      int lld=aesk::lin(l1,l1);
      vq[(l1-1)+6*i] -= qp[(lld-1)+6*i]*qs2*0.5;
      ct3 -= ra[l1-1]*ra[l1-1]*qp[(lld-1)+6*i]*qs2*0.5;
      vd[(l1-1)+3*i] += ra[l1-1]*qp[(lld-1)+6*i]*qs2;
      ct2 += qp[(lld-1)+6*i];
   }
   vs[i] += ct3;
   ct2 *= quadKernel[at[i]-1];
   for(int l1=1;l1<=3;l1++){
      vq[(l1-1)+6*i] += ct2;
      vd[(l1-1)+3*i] -= 2.0*ra[l1-1]*ct2;
      vs[i] += ct2*ra[l1-1]*ra[l1-1];
   }
}

extern "C" int gpu_setvsdq(int nat, int nelem,
      const int* at, const double* xyz, const double* q, const double* dipm,
      const double* qp, const double* gab3, const double* gab5,
      const double* dipKernel, const double* quadKernel,
      double* vs, double* vd, double* vq)
{
   std::lock_guard<std::mutex> lock(ctx_mutex);
   int *d_at=devcopy(at,nat);
   double *d_xyz=devcopy(xyz,(size_t)3*nat), *d_q=devcopy(q,nat);
   double *d_dipm=devcopy(dipm,(size_t)3*nat), *d_qp=devcopy(qp,(size_t)6*nat);
   double *d_g3=devcopy(gab3,(size_t)nat*nat), *d_g5=devcopy(gab5,(size_t)nat*nat);
   double *d_dk=devcopy(dipKernel,nelem), *d_qk=devcopy(quadKernel,nelem);
   double *d_vs=nullptr,*d_vd=nullptr,*d_vq=nullptr;
   cudaMalloc((void**)&d_vs,sizeof(double)*nat);
   cudaMalloc((void**)&d_vd,sizeof(double)*3*nat);
   cudaMalloc((void**)&d_vq,sizeof(double)*6*nat);
   int rc=0;
   if(!d_at||!d_xyz||!d_q||!d_dipm||!d_qp||!d_g3||!d_g5||!d_dk||!d_qk||!d_vs||!d_vd||!d_vq){
      rc=-1;
   } else {
      int tb=128;
      aes_setvsdq_kernel<<<(nat+tb-1)/tb,tb>>>(nat,d_at,d_xyz,d_q,d_dipm,d_qp,
            d_g3,d_g5,d_dk,d_qk,d_vs,d_vd,d_vq);
      if(cudaDeviceSynchronize()!=cudaSuccess) rc=-2;
      else{
         cudaMemcpy(vs,d_vs,sizeof(double)*nat,cudaMemcpyDeviceToHost);
         cudaMemcpy(vd,d_vd,sizeof(double)*3*nat,cudaMemcpyDeviceToHost);
         cudaMemcpy(vq,d_vq,sizeof(double)*6*nat,cudaMemcpyDeviceToHost);
      }
   }
   cudaFree(d_at);cudaFree(d_xyz);cudaFree(d_q);cudaFree(d_dipm);cudaFree(d_qp);
   cudaFree(d_g3);cudaFree(d_g5);cudaFree(d_dk);cudaFree(d_qk);
   cudaFree(d_vs);cudaFree(d_vd);cudaFree(d_vq);
   return rc;
}
