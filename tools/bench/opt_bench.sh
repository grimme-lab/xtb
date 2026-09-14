#!/usr/bin/env bash
# GFN2 --opt CPU vs GPU across pocket-relevant sizes. Bounded to 8 geometry
# cycles so we measure time-per-cycle (the metric that governs opt speed).
NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3; V=13.1
export LD_LIBRARY_PATH="$NVHPC/math_libs/$V/targets/x86_64-linux/lib:$NVHPC/REDIST/cuda/$V/targets/x86_64-linux/lib:$NVHPC/compilers/lib:${LD_LIBRARY_PATH:-}"
export XTBPATH=/mnt/e/Prasanna/xTB/xtb
CPU=/mnt/e/Prasanna/xTB/xtb/build/xtb
GPU=/mnt/e/Prasanna/xTB/xtb/build-gpushim/xtb
WORK=/tmp/optb; rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
printf '$opt\n   maxcycle=8\n$end\n' > opt.inp

cp /mnt/e/Prasanna/xTB/xtb/assets/inputs/xyz/taxol.xyz taxol.xyz
mkwater() { python3 - "$1" <<'PY'
import sys
N=int(sys.argv[1]); a=3.2; A=[]
for i in range(N):
 for j in range(N):
  for k in range(N):
   x,y,z=i*a,j*a,k*a
   A+=[("O",x,y,z),("H",x+0.757,y+0.586,z),("H",x-0.757,y+0.586,z)]
print(len(A)); print("w")
for s,x,y,z in A: print(f"{s} {x:.4f} {y:.4f} {z:.4f}")
PY
}
mkwater 4 > w64.xyz   # 192 atoms
mkwater 5 > w125.xyz  # 375 atoms

bench() {
  local file="$1" label="$2" bin="$3"; shift 3
  rm -f xtbopt.xyz xtbrestart
  local t0 t1 cyc tot
  t0=$(date +%s.%N)
  OMP_NUM_THREADS=8 "$bin" "$file" --gfn 2 --opt loose --input opt.inp "$@" > o.txt 2>&1
  t1=$(date +%s.%N)
  cyc=$(grep -ac "CYCLE" o.txt)
  tot=$(awk "BEGIN{print $t1-$t0}")
  local per="?"; [ "${cyc:-0}" -gt 0 ] && per=$(awk "BEGIN{printf \"%.2f\", $tot/$cyc}")
  printf "  %-4s %7.1fs  %2s cyc  %6ss/cyc\n" "$label" "$tot" "${cyc:-0}" "$per"
}
for fa in "taxol.xyz:113" "w64.xyz:192" "w125.xyz:375"; do
  f="${fa%%:*}"; n="${fa##*:}"
  echo "=== $f ($n atoms) GFN2 --opt loose, 8 cycles ==="
  bench "$f" CPU "$CPU"
  bench "$f" GPU "$GPU" --gpu
done
