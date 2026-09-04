#!/usr/bin/env bash
NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3; V=13.1
export LD_LIBRARY_PATH="$NVHPC/math_libs/$V/targets/x86_64-linux/lib:$NVHPC/REDIST/cuda/$V/targets/x86_64-linux/lib:$NVHPC/compilers/lib:${LD_LIBRARY_PATH:-}"
export XTBPATH=/mnt/e/Prasanna/xTB/xtb
CPU=/mnt/e/Prasanna/xTB/xtb/build/xtb
GPU=/mnt/e/Prasanna/xTB/xtb/build-gpushim/xtb
WORK=/tmp/taxbench; rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
cp /mnt/e/Prasanna/xTB/xtb/assets/inputs/xyz/taxol.xyz .
echo "taxol: $(head -1 taxol.xyz | tr -d ' \r') atoms, GFN2 single point"
run() {
  local label="$1" bin="$2" omp="$3"; shift 3
  local t0 t1 e
  t0=$(date +%s.%N)
  OMP_NUM_THREADS="$omp" "$bin" taxol.xyz --gfn 2 --sp "$@" > out.txt 2>&1
  t1=$(date +%s.%N)
  e=$(grep -a "TOTAL ENERGY" out.txt | grep -aoE "\-[0-9]+\.[0-9]+" | head -1)
  printf "%-16s OMP=%-2s  %7.2f s   E=%s\n" "$label" "$omp" "$(awk "BEGIN{print $t1-$t0}")" "$e"
}
echo "--- single thread ---"
run "CPU"  "$CPU" 1
run "GPU --gpu" "$GPU" 1 --gpu
echo "--- 8 threads ---"
run "CPU"  "$CPU" 8
run "GPU --gpu" "$GPU" 8 --gpu
