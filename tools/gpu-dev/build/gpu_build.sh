#!/usr/bin/env bash
export NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
export PATH="$NVHPC/compilers/bin:$HOME/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
cd /mnt/e/Prasanna/xTB/xtb
ninja -C build-gpu -j2 xtb
echo "NINJA_EXIT=$?"
