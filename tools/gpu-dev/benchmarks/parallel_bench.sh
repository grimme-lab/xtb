#!/usr/bin/env bash
# Benchmark the parallel-build GPU path vs CPU-same-path vs serial, on 300 real
# molecules, and verify GPU==CPU energies (parity) on the full set.
NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3; V=13.1
export LD_LIBRARY_PATH="$NVHPC/REDIST/math_libs/$V/targets/x86_64-linux/lib:$NVHPC/math_libs/$V/targets/x86_64-linux/lib:$NVHPC/REDIST/cuda/$V/targets/x86_64-linux/lib:$NVHPC/compilers/lib:${LD_LIBRARY_PATH:-}"
export XTBPATH=/mnt/e/Prasanna/xTB/xtb
XTB=/mnt/e/Prasanna/xTB/xtb/build-gpushim/xtb
SRC="/mnt/e/Prasanna/Research/CO2 Capture/Data/test"

WORK=/tmp/pbench; rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
i=0; for f in "$SRC"/*.xyz; do cp "$f" .; i=$((i+1)); [ "$i" -ge 300 ] && break; done
N=$(ls *.xyz | wc -l); echo "N=$N"

run() { # $1 label  $2 omp  $3 extra-env
  local t0 t1
  t0=$(date +%s.%N)
  env $3 OMP_NUM_THREADS=$2 "$XTB" --gfn 0 --gpu *.xyz > "out_$1.txt" 2>&1
  t1=$(date +%s.%N)
  printf "%-22s %6.2f s\n" "$1" "$(awk "BEGIN{print $t1-$t0}")"
}

echo "=== timings ==="
run "GPU_serial(OMP1)"   1  ""
run "GPU_parallel(OMP8)" 8  ""
run "CPU_parallel(OMP8)" 8  "XTB_BATCH_CPU=1"

echo "=== parity: GPU(OMP8) vs CPU(OMP8) energy tables ==="
ecol() { grep -aE '\.xyz ' "$1" | awk '{print $2, $(NF-1)}' | sort; }
if diff <(ecol out_GPU_parallel\(OMP8\).txt) <(ecol out_CPU_parallel\(OMP8\).txt) >/dev/null; then
  echo "PARITY OK: all $N GPU energies == CPU energies"
else
  echo "PARITY MISMATCH:"; diff <(ecol out_GPU_parallel\(OMP8\).txt) <(ecol out_CPU_parallel\(OMP8\).txt) | head
fi
