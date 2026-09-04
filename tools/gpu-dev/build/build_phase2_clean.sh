#!/usr/bin/env bash
export PATH="$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin"
export FC=gfortran CC=gcc
cd /mnt/e/Prasanna/xTB/xtb || exit 1

echo "=== fresh configure CPU (build/) ==="
rm -rf build
meson setup build --buildtype release -Doptimization=2 2>&1 | tail -3
echo "=== build CPU ==="
cok=0; for i in 1 2 3 4 5 6; do ninja -C build xtb && { cok=1; break; }; done
[ "$cok" = 1 ] && echo "CPU_OK" || echo "CPU_FAILED"

echo "=== fresh configure GPU (build-gpushim/, sm_86) ==="
rm -rf build-gpushim
meson setup build-gpushim --buildtype release -Doptimization=2 \
  -Ddefault_library=static -Dgpu_shim=true -Dgpu_arch=86 2>&1 | tail -3
echo "=== build GPU ==="
gok=0; for i in 1 2 3 4 5 6; do ninja -C build-gpushim xtb && { gok=1; break; }; done
[ "$gok" = 1 ] && echo "GPU_OK" || echo "GPU_FAILED"

ls -l build/xtb build-gpushim/xtb 2>/dev/null
echo "BUILD_DONE cpu=$cok gpu=$gok"
