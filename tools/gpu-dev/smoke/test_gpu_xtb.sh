#!/usr/bin/env bash
NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
MATH=$NVHPC/math_libs/13.1/targets/x86_64-linux/lib
CUDART=$(dirname "$(find $NVHPC -name 'libcudart.so*' 2>/dev/null | head -1)")
NVJIT=$(dirname "$(find $NVHPC -name 'libnvJitLink.so*' 2>/dev/null | head -1)")
export LD_LIBRARY_PATH="$MATH:$CUDART:$NVJIT:$NVHPC/compilers/lib:${LD_LIBRARY_PATH:-}"
export XTBPATH=/mnt/e/Prasanna/xTB/xtb
XTB=/mnt/e/Prasanna/xTB/xtb/build-gpushim/xtb
cd /mnt/e/Prasanna/xTB/win/test

echo "=== GPU path: xtb --gfn 0 --gpu *.xyz ==="
"$XTB" --gfn 0 --gpu *.xyz 2>&1 | grep -aE "batched GFN0 path|diagonalized|#|water|methane|ammonia|processed|throughput" | head -20

echo
echo "expected (CPU reference): ammonia -4.575742782498  methane -4.359317912047  water -4.366769919234"
