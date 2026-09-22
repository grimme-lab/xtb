#!/usr/bin/env bash
export PATH="$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin"
cd /mnt/e/Prasanna/xTB/xtb || exit 1
echo "=== clean configure (release, -O2, gfortran) ==="
meson setup build --buildtype release -Doptimization=2 2>&1 | tail -8
echo "=== build (ninja -j4) ==="
ninja -C build -j4
echo "NINJA_EXIT=$?"
ls -l build/xtb 2>/dev/null
