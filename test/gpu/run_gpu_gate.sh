#!/usr/bin/env bash
# GPU numerical gate for the batched eigensolver.
#
# Validates that cuSolver (cusolverDnDsygvd, the GPU backend in
# src/gpu/batched_eig.F90) reproduces LAPACK's eigenvalues on REAL GFN0 H/S
# matrices, including the bucket-padding scheme xtb --gpu-batch uses. This is the
# GPU gate from GPU_ACCELERATION_PLAN.md (energy/eigenvalue parity <= 1e-6).
#
# Why standalone: nvfortran (HPC SDK) currently ICEs on xtb's modern-Fortran
# dependency tree (toml-f deferred-length strings etc.), so the *full* xtb cannot
# yet be built with nvfortran. This test exercises the exact cuSolver path in
# isolation, which is the numerically meaningful part of the GPU port.
#
# Prereqs:
#   - NVIDIA HPC SDK (nvfortran) installed; set NVHPC below to its versioned dir.
#   - A CPU xtb build at ../../build (for dumping real matrices); XTBPATH set.
#   - A CUDA-capable GPU (developed/tested on an RTX 3050, compute capability 86).
set -e

NVHPC=${NVHPC:-/opt/nvidia/hpc_sdk/Linux_x86_64/26.3}
GPU_ARCH=${GPU_ARCH:-86}
XTB=${XTB:-../../build/xtb}
MOLS=${MOLS:-/tmp/divmols}
HSFILE=${HSFILE:-/tmp/gfn0_hs.txt}

export PATH="$NVHPC/compilers/bin:/usr/local/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="$NVHPC/compilers/lib:$NVHPC/math_libs/lib64:$LD_LIBRARY_PATH"
cd "$(dirname "$0")"

# 1. Dump real GFN0 (H,S) from the CPU build (skip if the file already exists).
if [ ! -f "$HSFILE" ]; then
  echo "==> dumping real GFN0 H/S via the CPU build"
  XTB_DUMP_HS="$HSFILE" "$XTB" --gfn 0 --gpu-batch "$MOLS"/*.xyz >/dev/null 2>&1 || true
fi
[ -f "$HSFILE" ] || { echo "ERROR: no $HSFILE (run a CPU --gpu-batch with XTB_DUMP_HS set)"; exit 1; }

# 2. Compile the gate with nvfortran + cuSolver/cuBLAS + OpenACC.
echo "==> compiling cusolver_gate (nvfortran -acc -gpu=cc$GPU_ARCH -cudalib=cusolver,cublas)"
nvfortran -acc -gpu=cc"$GPU_ARCH" -Mallocatable=03 -cudalib=cusolver,cublas \
  -o cusolver_gate cusolver_gate.f90 -llapack -lblas

# 3. Run on the GPU.
echo "==> running on GPU"
XTB_DUMP_HS="$HSFILE" ./cusolver_gate
