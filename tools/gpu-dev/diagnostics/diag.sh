#!/usr/bin/env bash
NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
V=13.1
MATH=$NVHPC/math_libs/$V/targets/x86_64-linux
cd /mnt/e/Prasanna/xTB/xtb
echo "=== direct nvcc compile of gpu_eig.cu ==="
"$NVHPC/cuda/bin/nvcc" -O2 -Xcompiler -fPIC -c src/gpu/gpu_eig.cu \
  -I"$MATH/include" -o /tmp/gpu_eig_test.o 2>&1 | head -40
echo "nvcc_exit=${PIPESTATUS[0]}"
