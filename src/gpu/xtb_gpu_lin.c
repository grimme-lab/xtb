/* This file is part of xtb.
 *
 * xtb is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * xtb is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with xtb.  If not, see <https://www.gnu.org/licenses/>.
 */

/* cuSOLVER-backed dense eigensolvers for xtb, driven from xtb_gpu_ctrl.
 *
 * Two entry points mirror the LAPACK routines they replace:
 *   xtb_gpu_sygvd  <- dsygvd(itype=1, jobz='v', uplo='u')   SCC
 *   xtb_gpu_syevd  <- dsyevd(jobz='v', uplo='u')            ANC Hessian
 * B (the overlap) is only read, so the caller can keep one copy for a whole
 * SCF.  A non-zero return means "solve on the host instead" with A, B and w
 * untouched, except for a positive return, which means A was partly
 * overwritten by a failed device-to-host copy and must not be reused.
 *
 * The interesting part is the *precision ramp* in xtb_gpu_sygvd.  On a
 * consumer/workstation Ada card FP64 runs at ~1/32 of FP32, so a plain FP64
 * port is barely faster than 32 CPU cores.  What makes the GPU win is running
 * the early SCC iterations in FP32 and switching to FP64 once the density is
 * close.  The switch criterion is max|dH_ii| -- the change in the *input*
 * Hamiltonian diagonal between consecutive solves.  Two reasons for that
 * choice over the SCF energy change:
 *   - it is known *before* the solve, so the current solve can already use the
 *     right precision;
 *   - it is a property of the Hamiltonian, not of the solver, so it is not
 *     polluted by the FP32 noise we are trying to ramp out of.
 * The switch is sticky: once in FP64 the SCF never goes back.
 *
 * Every SCF ramps the same way, the first one and every restart (geometry
 * optimisation calls xtb_gpu_ramp_reset() at the start of each).  A
 * warm-started SCF switches after a few solves, so most of its solves are
 * FP64 anyway; starting it in FP64 as well bought no accuracy and cost one
 * FP64 solve per optimisation cycle.  What guarantees the accuracy instead is
 * xtb_gpu_accept_density(): the SCC asks it before it stops, and if the last
 * solve was FP32 the answer is no, the ramp switches to FP64 and the SCC does
 * one more iteration.  The converged density, energy and gradient therefore
 * always come from an FP64 solve.
 */

/* clock_gettime/CLOCK_MONOTONIC are POSIX, not ISO C; meson compiles with a
 * strict -std= that would otherwise hide them. */
#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <cuda_runtime.h>
#include <cusolverDn.h>

/* ------------------------------------------------------------------ state */

static int gpu_dev = 0;
static int gpu_verbose = 0;
static int dev_open = 0;  /* device chosen and checked (xtb_gpu_open) */
static int cus_ready = 0; /* context and cuSOLVER handle exist (ensure_ready) */
static cusolverDnHandle_t cus;
static char dev_name[128] = "unknown";

/* device / host scratch, grown on demand and reused across solves */
static int cap_n = 0, cap_lwork = 0, cap_lwork32 = 0;
static double *dA, *dB, *dW, *dwork;
static float *dA32, *dB32, *dW32, *dwork32;
static int *dinfo;
static float *hA32, *hB32; /* host staging for the FP32 cast */
static int have32 = 0;

/* precision ramp */
static long solve_idx = 0;
static int scf_id = 0;
static int switched = 0;
static int last_fp32 = 0; /* the latest solve of this SCF was FP32 */
static double *h_prev_diag = NULL;
static int h_prev_n = 0;

/* stats */
static double t_h2d, t_solve, t_d2h, t_cast, t_init;
static long n_fp32, n_fp64, n_syevd;
static double t_syevd;

static double now(void) {
   struct timespec t;
   clock_gettime(CLOCK_MONOTONIC, &t);
   return t.tv_sec + 1e-9 * t.tv_nsec;
}

#define CK(x)                                                                  \
   do {                                                                        \
      cudaError_t e_ = (x);                                                    \
      if (e_ != cudaSuccess) {                                                 \
         fprintf(stderr, "[xtb-gpu] CUDA error: %s (%s:%d)\n",                 \
                 cudaGetErrorString(e_), __FILE__, __LINE__);                  \
         return -1;                                                            \
      }                                                                        \
   } while (0)

/* ------------------------------------------------------------- public API */

int xtb_gpu_compiled(void) { return 1; }

void xtb_gpu_name(char *buf, int buflen) {
   if (buflen <= 0) return;
   strncpy(buf, dev_name, (size_t)buflen - 1);
   buf[buflen - 1] = '\0';
}

/* Largest visible device, by SM count then by memory.  Picking this wrong is
 * expensive rather than merely suboptimal: on a host with a compute card plus
 * a small display adapter, the display adapter is slower than the CPU for
 * every path we offload, so a bad default turns --gpu into a slowdown.  CUDA's
 * own default ordering is "fastest first", which is a heuristic on clock and
 * capability and does not reliably put the compute card first. */
static int pick_device(int ndev) {
   int best = -1, best_sm = -1;
   size_t best_mem = 0;
   for (int i = 0; i < ndev; i++) {
      struct cudaDeviceProp p;
      if (cudaGetDeviceProperties(&p, i) != cudaSuccess) continue;
      if (p.major < 6) continue; /* pre-Pascal: not worth the transfers */
      if (p.multiProcessorCount > best_sm ||
          (p.multiProcessorCount == best_sm && p.totalGlobalMem > best_mem)) {
         best = i;
         best_sm = p.multiProcessorCount;
         best_mem = p.totalGlobalMem;
      }
   }
   return best;
}

int xtb_gpu_open(int dev, int verbose) {
   struct cudaDeviceProp p;
   int ndev = 0;
   double t = now();

   gpu_verbose = verbose;

   if (cudaGetDeviceCount(&ndev) != cudaSuccess || ndev <= 0) {
      fprintf(stderr, "[xtb-gpu] no CUDA device found\n");
      return -1;
   }
   if (dev < 0) {
      dev = pick_device(ndev);
      if (dev < 0) {
         fprintf(stderr, "[xtb-gpu] no device with compute capability 6.0 or "
                         "newer\n");
         return -1;
      }
   }
   gpu_dev = dev;
   if (dev >= ndev) {
      fprintf(stderr, "[xtb-gpu] device %d out of range (%d present)\n", dev, ndev);
      return -1;
   }
   if (cudaGetDeviceProperties(&p, dev) == cudaSuccess)
      snprintf(dev_name, sizeof(dev_name), "%s, sm%d%d, %.1f GB", p.name,
               p.major, p.minor, p.totalGlobalMem / 1073741824.0);
   dev_open = 1;
   t_init += now() - t;
   return 0;
}

/* The context, the cuSOLVER handle and the info word are created by the first
 * solve that actually goes to the device, not by xtb_gpu_open: creating the
 * context and tearing it down at exit costs ~0.17 s, which a run whose
 * matrices all stay below the size threshold would pay for nothing.  A
 * failure here is sticky and sends every later solve to the host. */
static int ensure_ready(void) {
   static int failed = 0;
   double t;

   if (cus_ready) return 0;
   if (!dev_open || failed) return -1;
   t = now();
   failed = 1;
   CK(cudaSetDevice(gpu_dev));
   if (cusolverDnCreate(&cus) != CUSOLVER_STATUS_SUCCESS) {
      fprintf(stderr, "[xtb-gpu] cusolverDnCreate failed, solving on the host\n");
      return -1;
   }
   CK(cudaMalloc((void **)&dinfo, sizeof(int)));
   failed = 0;
   cus_ready = 1;
   t_init += now() - t;
   return 0;
}

/* The device the solvers run on, for the banner; -1 before a successful open. */
int xtb_gpu_device(void) { return dev_open ? gpu_dev : -1; }

void xtb_gpu_ramp_reset(void) {
   solve_idx = 0;
   scf_id++;
   switched = 0;
   last_fp32 = 0;
}

/* May the SCC stop on the density of the latest solve?  Not if that solve was
 * FP32: then switch to FP64 for the rest of this SCF and let the SCC iterate
 * once more. */
int xtb_gpu_accept_density(void) {
   if (!last_fp32) return 1;
   switched = 1;
   last_fp32 = 0;
   if (gpu_verbose)
      fprintf(stderr, "[xtb-gpu] scf=%d converged on an FP32 solve, one more in FP64\n", scf_id);
   return 0;
}

void xtb_gpu_report(void) {
   if (!dev_open) return;
   if (!cus_ready) {
      fprintf(stderr, "[xtb-gpu] device %d (%s) not used: no matrix reached the "
                      "size threshold\n", gpu_dev, dev_name);
      return;
   }
   fprintf(stderr,
           "[xtb-gpu] device %d (%s)\n"
           "[xtb-gpu] sygvd: fp32=%ld fp64=%ld | h2d=%.3f solve=%.3f d2h=%.3f "
           "cast=%.3f init=%.3f\n",
           gpu_dev, dev_name, n_fp32, n_fp64, t_h2d, t_solve, t_d2h, t_cast,
           t_init);
   if (n_syevd)
      fprintf(stderr, "[xtb-gpu] syevd: calls=%ld t=%.3f\n", n_syevd, t_syevd);
}

/* -------------------------------------------------------------- scratch */

/* Release the FP32 buffers together with their workspace. */
static void free32(void) {
   cudaFree(dA32); cudaFree(dB32); cudaFree(dW32); cudaFree(dwork32);
   free(hA32); free(hB32);
   dA32 = dB32 = dW32 = dwork32 = NULL; hA32 = hB32 = NULL;
   have32 = 0;
   cap_lwork32 = 0;
}

static int ensure_n(int n, int fp32) {
   if (n > cap_n) {
      /* Free everything sized for the old dimension, the workspaces included:
       * their sizes restart from zero, so they would never be freed later.
       * (An optimization grows from the SCC dimension to the 3N of the ANC
       * Hessian.) */
      cudaFree(dA); cudaFree(dB); cudaFree(dW); cudaFree(dwork);
      dA = dB = dW = dwork = NULL;
      cap_n = cap_lwork = 0;
      free32();
      size_t nn = (size_t)n * n;
      CK(cudaMalloc((void **)&dA, nn * sizeof(double)));
      CK(cudaMalloc((void **)&dB, nn * sizeof(double)));
      CK(cudaMalloc((void **)&dW, (size_t)n * sizeof(double)));
      cap_n = n;
   }
   if (fp32 && !have32) {
      size_t nn = (size_t)cap_n * cap_n;
      if (cudaMalloc((void **)&dA32, nn * sizeof(float)) != cudaSuccess ||
          cudaMalloc((void **)&dB32, nn * sizeof(float)) != cudaSuccess ||
          cudaMalloc((void **)&dW32, (size_t)cap_n * sizeof(float)) != cudaSuccess ||
          !(hA32 = (float *)malloc(nn * sizeof(float))) ||
          !(hB32 = (float *)malloc(nn * sizeof(float)))) {
         free32();
         return -1;
      }
      have32 = 1;
   }
   return 0;
}

/* ---------------------------------------------------------------- solvers */

static int solve_fp64(int n, double *A, const double *B, double *w) {
   double t;
   int lw = 0, info = 0;
   size_t nn = (size_t)n * n * sizeof(double);

   if (ensure_n(n, 0)) return -1;
   if (cusolverDnDsygvd_bufferSize(cus, CUSOLVER_EIG_TYPE_1,
                                   CUSOLVER_EIG_MODE_VECTOR,
                                   CUBLAS_FILL_MODE_UPPER, n, dA, n, dB, n, dW,
                                   &lw) != CUSOLVER_STATUS_SUCCESS)
      return -1;
   if (lw > cap_lwork) {
      cudaFree(dwork);
      dwork = NULL;
      cap_lwork = 0;
      CK(cudaMalloc((void **)&dwork, (size_t)lw * sizeof(double)));
      cap_lwork = lw;
   }

   t = now();
   CK(cudaMemcpy(dA, A, nn, cudaMemcpyHostToDevice));
   CK(cudaMemcpy(dB, B, nn, cudaMemcpyHostToDevice));
   CK(cudaDeviceSynchronize());
   t_h2d += now() - t;

   t = now();
   if (cusolverDnDsygvd(cus, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR,
                        CUBLAS_FILL_MODE_UPPER, n, dA, n, dB, n, dW, dwork,
                        cap_lwork, dinfo) != CUSOLVER_STATUS_SUCCESS)
      return -1;
   CK(cudaDeviceSynchronize());
   t_solve += now() - t;

   CK(cudaMemcpy(&info, dinfo, sizeof(int), cudaMemcpyDeviceToHost));
   if (info != 0) return -1;

   t = now();
   CK(cudaMemcpy(w, dW, (size_t)n * sizeof(double), cudaMemcpyDeviceToHost));
   if (cudaMemcpy(A, dA, nn, cudaMemcpyDeviceToHost) != cudaSuccess)
      return 1; /* A partly overwritten: the caller must not fall back on it */
   t_d2h += now() - t;
   return 0;
}

static int solve_fp32(int n, double *A, const double *B, double *w) {
   double t;
   int lw = 0, info = 0, i;
   size_t nn = (size_t)n * n, k;
   float *w32;

   if (ensure_n(n, 1)) return -1;

   t = now();
#pragma omp parallel for schedule(static)
   for (k = 0; k < nn; k++) {
      hA32[k] = (float)A[k];
      hB32[k] = (float)B[k];
   }
   t_cast += now() - t;

   if (cusolverDnSsygvd_bufferSize(cus, CUSOLVER_EIG_TYPE_1,
                                   CUSOLVER_EIG_MODE_VECTOR,
                                   CUBLAS_FILL_MODE_UPPER, n, dA32, n, dB32, n,
                                   dW32, &lw) != CUSOLVER_STATUS_SUCCESS)
      return -1;
   if (lw > cap_lwork32) {
      cudaFree(dwork32);
      dwork32 = NULL;
      cap_lwork32 = 0;
      CK(cudaMalloc((void **)&dwork32, (size_t)lw * sizeof(float)));
      cap_lwork32 = lw;
   }

   t = now();
   CK(cudaMemcpy(dA32, hA32, nn * sizeof(float), cudaMemcpyHostToDevice));
   CK(cudaMemcpy(dB32, hB32, nn * sizeof(float), cudaMemcpyHostToDevice));
   CK(cudaDeviceSynchronize());
   t_h2d += now() - t;

   t = now();
   if (cusolverDnSsygvd(cus, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR,
                        CUBLAS_FILL_MODE_UPPER, n, dA32, n, dB32, n, dW32,
                        dwork32, cap_lwork32, dinfo) != CUSOLVER_STATUS_SUCCESS)
      return -1;
   CK(cudaDeviceSynchronize());
   t_solve += now() - t;

   CK(cudaMemcpy(&info, dinfo, sizeof(int), cudaMemcpyDeviceToHost));
   if (info != 0) return -1;

   w32 = (float *)malloc((size_t)n * sizeof(float));
   if (!w32) return -1;
   t = now();
   CK(cudaMemcpy(hA32, dA32, nn * sizeof(float), cudaMemcpyDeviceToHost));
   CK(cudaMemcpy(w32, dW32, (size_t)n * sizeof(float), cudaMemcpyDeviceToHost));
   t_d2h += now() - t;

   t = now();
#pragma omp parallel for schedule(static)
   for (k = 0; k < nn; k++) A[k] = (double)hA32[k];
   for (i = 0; i < n; i++) w[i] = (double)w32[i];
   t_cast += now() - t;
   free(w32);
   return 0;
}

/* --------------------------------------------------------------- dispatch */

int xtb_gpu_sygvd(int n, double *A, const double *B, double *w, double fp32_tol) {
   double dh = -1.0;
   int use32, rc, i;

   last_fp32 = 0; /* a failed or host solve is FP64 */
   if (n <= 0 || ensure_ready()) return -1;
   solve_idx++;

   /* max|dH_ii| against the previous solve.  A still holds the *input*
    * Hamiltonian here -- read it before cuSOLVER overwrites it. */
   if (fp32_tol > 0.0) {
      if (h_prev_n != n) {
         free(h_prev_diag);
         h_prev_diag = (double *)malloc((size_t)n * sizeof(double));
         h_prev_n = h_prev_diag ? n : 0;
         solve_idx = 1; /* no usable history for this dimension */
      }
      if (h_prev_diag) {
         if (solve_idx > 1) {
            dh = 0.0;
            for (i = 0; i < n; i++) {
               double d = fabs(A[(size_t)i * n + i] - h_prev_diag[i]);
               if (d > dh) dh = d;
            }
         }
         for (i = 0; i < n; i++) h_prev_diag[i] = A[(size_t)i * n + i];
      }
   }

   if (fp32_tol > 0.0) {
      if (!switched && dh >= 0.0 && dh < fp32_tol) switched = 1; /* sticky */
      use32 = !switched;
   } else {
      use32 = 0;
   }

   if (gpu_verbose)
      fprintf(stderr, "[xtb-gpu] scf=%d solve=%ld n=%d max|dH_ii|=%.3e -> %s\n",
              scf_id, solve_idx, n, dh, use32 ? "fp32" : "fp64");

   rc = use32 ? solve_fp32(n, A, B, w) : solve_fp64(n, A, B, w);
   if (rc < 0 && use32) {
      /* The FP32 solve failed before touching A: S is positive definite in
       * FP64 but not in FP32, or the FP32 buffers do not fit.  FP64 on the
       * device for the rest of this SCF, instead of a host solve each time. */
      if (gpu_verbose)
         fprintf(stderr, "[xtb-gpu] scf=%d FP32 solve failed, FP64 from here\n", scf_id);
      switched = 1;
      use32 = 0;
      rc = solve_fp64(n, A, B, w);
   }
   if (rc != 0) return rc;
   if (use32) n_fp32++; else n_fp64++;
   last_fp32 = use32;
   return 0;
}

int xtb_gpu_syevd(int n, double *A, double *w) {
   double t = now();
   int lw = 0, info = 1;
   size_t nn = (size_t)n * n * sizeof(double);

   if (n <= 0 || ensure_ready()) return -1;
   if (ensure_n(n, 0)) return -1;
   if (cusolverDnDsyevd_bufferSize(cus, CUSOLVER_EIG_MODE_VECTOR,
                                   CUBLAS_FILL_MODE_UPPER, n, dA, n, dW,
                                   &lw) != CUSOLVER_STATUS_SUCCESS)
      return -1;
   if (lw > cap_lwork) {
      cudaFree(dwork);
      dwork = NULL;
      cap_lwork = 0;
      CK(cudaMalloc((void **)&dwork, (size_t)lw * sizeof(double)));
      cap_lwork = lw;
   }
   CK(cudaMemcpy(dA, A, nn, cudaMemcpyHostToDevice));
   if (cusolverDnDsyevd(cus, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,
                        n, dA, n, dW, dwork, cap_lwork,
                        dinfo) != CUSOLVER_STATUS_SUCCESS)
      return -1;
   CK(cudaMemcpy(&info, dinfo, sizeof(int), cudaMemcpyDeviceToHost));
   if (info != 0) return -1;
   CK(cudaMemcpy(w, dW, (size_t)n * sizeof(double), cudaMemcpyDeviceToHost));
   if (cudaMemcpy(A, dA, nn, cudaMemcpyDeviceToHost) != cudaSuccess)
      return 1; /* A partly overwritten: the caller must not fall back on it */
   t_syevd += now() - t;
   n_syevd++;
   return 0;
}
