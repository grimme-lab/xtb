#!/usr/bin/env bash
# Run the end-to-end GFN1/GFN2 CUDA single-point, gradient and --opt gate.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
NVHPC=${NVHPC:-/opt/nvidia/hpc_sdk/Linux_x86_64/26.3}
CUDA_VER=${CUDA_VER:-13.1}
CPU_XTB=${CPU_XTB:-"$ROOT/build/xtb"}
GPU_XTB=${GPU_XTB:-"$ROOT/build-gpushim/xtb"}
FIXTURE=${FIXTURE:-"$ROOT/test/gpu/water.xyz"}

export LD_LIBRARY_PATH="$NVHPC/math_libs/$CUDA_VER/targets/x86_64-linux/lib:$NVHPC/REDIST/cuda/$CUDA_VER/targets/x86_64-linux/lib:${LD_LIBRARY_PATH:-}"

# GPU AES (GFN2 setvsdq/mmompop/buildH1/aniso) is opt-in; force it on so the gate
# validates the GPU AES path bit-exactly even though it is off by default.
export XTB_GPU_AES=${XTB_GPU_AES:-1}

python3 "$ROOT/test/gpu/gfn12_gpu_gate.py" \
  --cpu "$CPU_XTB" \
  --gpu "$GPU_XTB" \
  --fixture "$FIXTURE"
