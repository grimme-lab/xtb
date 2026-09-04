#!/usr/bin/env bash
set -e
export NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
export PATH="$NVHPC/compilers/bin:$HOME/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
export FC=nvfortran
export CC=nvc
cd /mnt/e/Prasanna/xTB/xtb

# Fresh GPU build dir (RTX 3050 = compute capability 8.6).
rm -rf build-gpu
meson setup build-gpu \
  --buildtype release \
  -Doptimization=1 \
  -Ddefault_library=static \
  -Dgpu=true \
  -Dgpu_arch=86 \
  -Dcusolver=true \
  2>&1 | tail -25
