#!/usr/bin/env bash
set -euo pipefail
export PATH="$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin"
export FC=gfortran CC=gcc
cd /mnt/e/Prasanna/xTB/xtb || exit 1
echo "=== build (incremental) ==="
ok=0; for i in 1 2 3 4; do ninja -C build-gpushim xtb && { ok=1; break; }; done
[ "$ok" = 1 ] || { echo BUILD_FAILED; exit 1; }
echo "BUILD OK"

NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
CUDA_VER=13.1
MATH=$NVHPC/math_libs/$CUDA_VER/targets/x86_64-linux/lib
MATH_REDIST=$NVHPC/REDIST/math_libs/$CUDA_VER/targets/x86_64-linux/lib
CUDART=$NVHPC/REDIST/cuda/$CUDA_VER/targets/x86_64-linux/lib
NVJIT=$NVHPC/REDIST/cuda/$CUDA_VER/targets/x86_64-linux/lib
export LD_LIBRARY_PATH="$MATH_REDIST:$MATH:$CUDART:$NVJIT:$NVHPC/compilers/lib:${LD_LIBRARY_PATH:-}"
export XTBPATH=/mnt/e/Prasanna/xTB/xtb
export OMP_NUM_THREADS=8
XTB=/mnt/e/Prasanna/xTB/xtb/build-gpushim/xtb
cd /mnt/e/Prasanna/xTB/win/test

echo "=== linked CUDA runtime deps ==="
ldd "$XTB" | grep -E 'cuda|cusolver|cublas|nvJit|not found' || true

echo "=== correctness + determinism: 3 parallel GPU runs (OMP=8) ==="
echo "expected CPU: ammonia -4.575742782498  methane -4.359317912047  water -4.366769919234"
for run in 1 2 3; do
  echo "--- run $run ---"
  "$XTB" --gfn 0 --gpu ammonia.xyz methane.xyz water.xyz | tee "/tmp/xtb_parallel_gpu_run_${run}.out"
  grep -aE "ammonia|methane|water|processed" "/tmp/xtb_parallel_gpu_run_${run}.out" | grep -avE "building" || true
done

echo "=== CPU fallback same path (OMP=8) ==="
XTB_BATCH_CPU=1 "$XTB" --gfn 0 --gpu ammonia.xyz methane.xyz water.xyz | tee /tmp/xtb_parallel_cpu_fallback.out
