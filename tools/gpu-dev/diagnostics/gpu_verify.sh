#!/usr/bin/env bash
# Clean PATH (avoid the Windows interop PATH which has spaces/parens).
export NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
export PATH="$NVHPC/compilers/bin:$HOME/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

echo "=== compilers ==="
nvfortran --version | head -2
nvc --version | head -2
echo
echo "=== meson / ninja ==="
meson --version 2>/dev/null && echo "meson ok" || echo "meson MISSING"
ninja --version 2>/dev/null && echo "ninja ok" || echo "ninja MISSING"
echo
echo "=== cuSolver / cuBLAS shared libs in SDK ==="
find "$NVHPC" -name 'libcusolver*.so*' -o -name 'libcublas*.so*' 2>/dev/null | head
echo
echo "=== cusolverDn Fortran module (cusolverdn.mod) ==="
find "$NVHPC" -iname 'cusolverdn*.mod' 2>/dev/null | head
echo
echo "=== math_libs / cuda dirs ==="
ls -d "$NVHPC"/math_libs* "$NVHPC"/cuda* 2>/dev/null
