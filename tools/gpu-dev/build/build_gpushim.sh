#!/usr/bin/env bash
export PATH="$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin"
export FC=gfortran CC=gcc
cd /mnt/e/Prasanna/xTB/xtb || exit 1

echo "=== configure (gpu_shim, static, -O2, gfortran) ==="
rm -rf build-gpushim
meson setup build-gpushim --buildtype release -Doptimization=2 \
  -Ddefault_library=static -Dgpu_shim=true 2>&1 | tail -15

echo "=== build xtb (resume-loop for Fortran module race) ==="
ok=0
for i in $(seq 1 8); do
  echo "--- attempt $i ---"
  if ninja -C build-gpushim xtb; then ok=1; echo "BUILD OK ($i)"; break; fi
done
echo "----"
[ "$ok" = 1 ] && ls -l build-gpushim/xtb && echo BINARY_OK || echo BINARY_MISSING
