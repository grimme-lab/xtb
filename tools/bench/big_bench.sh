#!/usr/bin/env bash
# Generate a large water cluster and benchmark CPU vs GPU for GFN1 and GFN2.
NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3; V=13.1
export LD_LIBRARY_PATH="$NVHPC/math_libs/$V/targets/x86_64-linux/lib:$NVHPC/REDIST/cuda/$V/targets/x86_64-linux/lib:$NVHPC/compilers/lib:${LD_LIBRARY_PATH:-}"
export XTBPATH=/mnt/e/Prasanna/xTB/xtb
CPU=/mnt/e/Prasanna/xTB/xtb/build/xtb
GPU=/mnt/e/Prasanna/xTB/xtb/build-gpushim/xtb
WORK=/tmp/bigbench; rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"

N=${1:-6}   # grid NxNxN waters
python3 - "$N" > cluster.xyz <<'PY'
import sys
N=int(sys.argv[1]); a=3.1
atoms=[]
for i in range(N):
 for j in range(N):
  for k in range(N):
   ox,oy,oz=i*a,j*a,k*a
   atoms.append(("O",ox,oy,oz))
   atoms.append(("H",ox+0.757,oy+0.586,oz))
   atoms.append(("H",ox-0.757,oy+0.586,oz))
print(len(atoms)); print("water cluster")
for s,x,y,z in atoms: print(f"{s} {x:.4f} {y:.4f} {z:.4f}")
PY
NAT=$(head -1 cluster.xyz)
echo "cluster: $NAT atoms  (~$((NAT*2)) basis fns est.)"
run() {
  local label="$1" bin="$2" gfn="$3" omp="$4"; shift 4
  local t0 t1 e
  t0=$(date +%s.%N)
  OMP_NUM_THREADS="$omp" "$bin" cluster.xyz --gfn "$gfn" --sp "$@" > out.txt 2>&1
  local rc=$?
  t1=$(date +%s.%N)
  e=$(grep -a "TOTAL ENERGY" out.txt | grep -aoE "\-[0-9]+\.[0-9]+" | head -1)
  printf "%-14s gfn%s OMP=%-2s  %8.2f s   E=%s  rc=%s\n" "$label" "$gfn" "$omp" "$(awk "BEGIN{print $t1-$t0}")" "$e" "$rc"
}
for G in 1 2; do
  echo "--- GFN$G ---"
  run "CPU"  "$CPU" "$G" 8
  run "GPU"  "$GPU" "$G" 8 --gpu
done
