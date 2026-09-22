#!/usr/bin/env bash
NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
echo "=== nvcc ==="
ls "$NVHPC"/cuda/bin/nvcc 2>/dev/null && "$NVHPC"/cuda/bin/nvcc --version 2>/dev/null | tail -2 || echo "no nvcc in HPC SDK cuda"
echo "=== cuda include (cusolverDn.h) ==="
find "$NVHPC"/cuda -name cusolverDn.h 2>/dev/null | head -1
find "$NVHPC"/math_libs -name cusolverDn.h 2>/dev/null | head -1
echo "=== cuda runtime + cusolver .so for linking ==="
find "$NVHPC" -name 'libcudart.so' 2>/dev/null | head -1
find "$NVHPC" -name 'libcusolver.so' 2>/dev/null | head -1
echo "=== gfortran + lapack ==="
which gfortran; ls /usr/lib/x86_64-linux-gnu/liblapack.so* 2>/dev/null | head -1
