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
