// CUDA-C batched generalized symmetric-definite eigensolver, callable from
// gfortran via iso_c_binding. Sidesteps nvfortran (which can't build xtb's
// dependency tree) -- this is plain CUDA C compiled with nvcc into an object
// the normal gfortran xtb can link.
//
//   solves  H_k C_k = S_k C_k diag(W_k),  k = 1..nbatch   (itype=1, upper)
//
// H, S : nbatch contiguous n*n column-major blocks (symmetric; upper used)
// W    : nbatch contiguous length-n blocks (out, ascending eigenvalues)
//
// First version loops cusolverDnDsygvd with device buffers + workspace + handle
// reused across the batch (data staged per system). This is the honest baseline
// to measure GPU-vs-CPU for many small systems before optimizing to a fully
// batched (syevjBatched) pipeline.
#include <cusolverDn.h>
#include <cuda_runtime.h>
#include <cstdio>

extern "C" int gpu_sygvd_batch(int n, int nbatch,
                               const double* H, const double* S, double* W)
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
      cudaMemcpy(W + (size_t)k * n, dW, sizeof(double) * n, cudaMemcpyDeviceToHost);
      if (st != CUSOLVER_STATUS_SUCCESS) rc = (int)st;
   }
   cudaDeviceSynchronize();

   cudaFree(dA); cudaFree(dB); cudaFree(dW); cudaFree(dWork); cudaFree(dInfo);
   cusolverDnDestroy(handle);
   return rc;
}
