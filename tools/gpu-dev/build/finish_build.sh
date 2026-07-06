#!/usr/bin/env bash
export PATH="$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin"
cd /mnt/e/Prasanna/xTB/xtb || exit 1
# Resume the build until it converges (works around the meson/gfortran Fortran
# module-ordering race on a fresh parallel build: each pass builds more .mod).
for i in $(seq 1 8); do
  echo "===== build attempt $i ====="
  if ninja -C build -j4 xtb; then
    echo "BUILD OK on attempt $i"
    break
  fi
done
echo "----"
ls -l build/xtb 2>/dev/null && echo "BINARY_OK" || echo "BINARY_MISSING"
