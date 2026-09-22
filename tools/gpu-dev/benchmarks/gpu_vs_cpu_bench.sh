#!/usr/bin/env bash
# Fair GPU-vs-CPU benchmark on the SAME energy path (build-gpushim binary):
# only the batched diagonalization differs (GPU cuSolver vs CPU LAPACK).
NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
MATH=$NVHPC/math_libs/13.1/targets/x86_64-linux/lib
CUDART=$(dirname "$(find $NVHPC -name 'libcudart.so*' 2>/dev/null | head -1)")
NVJIT=$(dirname "$(find $NVHPC -name 'libnvJitLink.so*' 2>/dev/null | head -1)")
export LD_LIBRARY_PATH="$MATH:$CUDART:$NVJIT:$NVHPC/compilers/lib:${LD_LIBRARY_PATH:-}"
export XTBPATH=/mnt/e/Prasanna/xTB/xtb
XTB=/mnt/e/Prasanna/xTB/xtb/build-gpushim/xtb
SRC="/mnt/e/Prasanna/Research/CO2 Capture/Data/test"

WORK=/tmp/gcbench; rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
i=0
for f in "$SRC"/*.xyz; do
  cp "$f" "$WORK/" 2>/dev/null; i=$((i+1)); [ "$i" -ge 300 ] && break
done
N=$(ls *.xyz | wc -l)
echo "N=$N (single-process batched GFN0 energy path; only diag differs)"

t0=$(date +%s.%N); XTB_BATCH_CPU=1 "$XTB" --gfn 0 --gpu *.xyz >/dev/null 2>&1; t1=$(date +%s.%N)
echo "CPU $(awk "BEGIN{printf \"%.2f\", $t1-$t0}")"

t0=$(date +%s.%N); "$XTB" --gfn 0 --gpu *.xyz >/dev/null 2>&1; t1=$(date +%s.%N)
echo "GPU $(awk "BEGIN{printf \"%.2f\", $t1-$t0}")"
echo DONE
