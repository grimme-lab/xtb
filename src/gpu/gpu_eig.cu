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
#include <cuda_runtime.h>

extern "C" int gpu_sygvd_batch(int n, int nbatch,
                               double* H, const double* S, double* W)
{
   if (n <= 0 || nbatch <= 0) return 0;

   cusolverDnHandle_t handle = nullptr;
   if (cusolverDnCreate(&handle) != CUSOLVER_STATUS_SUCCESS) return -1;

   const size_t nn = (size_t)n * (size_t)n;
   double *dA = nullptr, *dB = nullptr, *dW = nullptr, *dWork = nullptr;
   int *dInfo = nullptr, lwork = 0, rc = 0;

   cudaMalloc((void**)&dA, sizeof(double) * nn);
   cudaMalloc((void**)&dB, sizeof(double) * nn);
   cudaMalloc((void**)&dW, sizeof(double) * n);
   cudaMalloc((void**)&dInfo, sizeof(int));

   cusolverDnDsygvd_bufferSize(handle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR,
                               CUBLAS_FILL_MODE_UPPER, n, dA, n, dB, n, dW, &lwork);
   cudaMalloc((void**)&dWork, sizeof(double) * lwork);

   for (int k = 0; k < nbatch; ++k) {
      cudaMemcpy(dA, H + (size_t)k * nn, sizeof(double) * nn, cudaMemcpyHostToDevice);
      cudaMemcpy(dB, S + (size_t)k * nn, sizeof(double) * nn, cudaMemcpyHostToDevice);
      cusolverStatus_t st = cusolverDnDsygvd(handle, CUSOLVER_EIG_TYPE_1,
            CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,
            n, dA, n, dB, n, dW, dWork, lwork, dInfo);
      // copy eigenvectors (in dA) and eigenvalues (in dW) back to the host
      cudaMemcpy(H + (size_t)k * nn, dA, sizeof(double) * nn, cudaMemcpyDeviceToHost);
      cudaMemcpy(W + (size_t)k * n,  dW, sizeof(double) * n,  cudaMemcpyDeviceToHost);
      if (st != CUSOLVER_STATUS_SUCCESS) rc = (int)st;
   }
   cudaDeviceSynchronize();

   cudaFree(dA); cudaFree(dB); cudaFree(dW); cudaFree(dWork); cudaFree(dInfo);
   cusolverDnDestroy(handle);
   return rc;
}
