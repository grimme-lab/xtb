#!/usr/bin/env bash
set -e
NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3; V=13.1
export LD_LIBRARY_PATH="$NVHPC/math_libs/$V/targets/x86_64-linux/lib:$NVHPC/REDIST/cuda/$V/targets/x86_64-linux/lib:$NVHPC/compilers/lib:${LD_LIBRARY_PATH:-}"
export PATH="$NVHPC/compilers/bin:$NVHPC/cuda/$V/bin:$PATH"
cd /mnt/e/Prasanna/xTB/xtb
echo "=== ninja build (serial xtb target) ==="
ninja -C build-gpushim -j1 xtb 2>&1 | tail -30
echo "=== build exit: ${PIPESTATUS[0]} ==="
ls -la build-gpushim/xtb
