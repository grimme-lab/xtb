#!/usr/bin/env bash
# Convenience launcher for the GPU-accelerated xtb build (run inside WSL2).
#
#   ./xtb-gpu.sh mol.xyz --gfn 2 --opt
#
# By default the GPU runs the analytical gradient + SCF diagonalization (the big
# wins). GFN2 AES stays on the multi-threaded CPU because it is faster there at
# pocket scale. To force the (bit-exact) GPU AES path on very large systems:
#
#   XTB_GPU_AES=1 ./xtb-gpu.sh mol.xyz --gfn 2 --opt
#
set -euo pipefail
NVHPC=${NVHPC:-/opt/nvidia/hpc_sdk/Linux_x86_64/26.3}
V=${CUDA_VER:-13.1}
ROOT=$(cd "$(dirname "$0")" && pwd)

export LD_LIBRARY_PATH="/usr/lib/wsl/lib:$NVHPC/math_libs/$V/targets/x86_64-linux/lib:$NVHPC/REDIST/cuda/$V/targets/x86_64-linux/lib:$NVHPC/compilers/lib:${LD_LIBRARY_PATH:-}"
export XTBPATH="${XTBPATH:-$ROOT}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"

exec "$ROOT/build-gpushim/xtb" "$@" --gpu
