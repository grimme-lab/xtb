#!/usr/bin/env bash
export PATH="$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin"
export FC=gfortran CC=gcc
cd /mnt/e/Prasanna/xTB/xtb || exit 1

echo "=== rebuild CPU reference (build/) ==="
cok=0; for i in 1 2 3 4 5; do ninja -C build xtb && { cok=1; break; }; done
[ "$cok" = 1 ] && echo "CPU_OK" || echo "CPU_FAILED"

echo "=== rebuild GPU shim (build-gpushim/) ==="
gok=0; for i in 1 2 3 4 5; do ninja -C build-gpushim xtb && { gok=1; break; }; done
[ "$gok" = 1 ] && echo "GPU_OK" || echo "GPU_FAILED"

ls -l build/xtb build-gpushim/xtb 2>/dev/null
echo "BUILD_DONE cpu=$cok gpu=$gok"
