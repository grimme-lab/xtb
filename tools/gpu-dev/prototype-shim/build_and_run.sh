#!/usr/bin/env bash
set -e
NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
CUDA=$NVHPC/cuda
MATH=$NVHPC/math_libs/13.1/targets/x86_64-linux
CUDART_DIR=$(dirname "$(find "$NVHPC" -name 'libcudart.so' 2>/dev/null | head -1)")
cd "$(dirname "$0")"

echo "=== nvcc: compile CUDA-C shim ==="
"$CUDA/bin/nvcc" -O2 -c gpu_eig.cu -I"$MATH/include" -o gpu_eig.o
echo "  -> gpu_eig.o"

echo "=== gfortran: compile + link (CUDA libs into a gfortran program) ==="
gfortran -O2 gpu_bench.f90 gpu_eig.o -o gpu_bench \
  -L"$MATH/lib" -lcusolver \
  -L"$CUDART_DIR" -lcudart \
  -llapack -lblas -lstdc++ -lm
echo "  -> gpu_bench"

export LD_LIBRARY_PATH="$MATH/lib:$CUDART_DIR:$NVHPC/compilers/lib:${LD_LIBRARY_PATH:-}"

echo
echo "######## n=48,  nbatch=200  (typical small molecule, screening batch) ########"
./gpu_bench 48 200
echo
echo "######## n=120, nbatch=200  (medium) ########"
./gpu_bench 120 200
echo
echo "######## n=400, nbatch=50   (larger single systems) ########"
./gpu_bench 400 50
