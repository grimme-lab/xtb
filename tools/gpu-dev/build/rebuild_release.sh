#!/usr/bin/env bash
export PATH="$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin"
cd /mnt/e/Prasanna/xTB/xtb
echo "=== bump optimization to -O2 (release) ==="
meson configure build -Doptimization=2 2>&1 | grep -i optimization | tail -2
echo "=== full rebuild (-j4) ==="
ninja -C build -j4 xtb
echo "NINJA_EXIT=$?"
ls -l build/xtb
