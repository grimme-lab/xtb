#!/bin/bash
# Throughput benchmark for the --gpu-batch driver.
#
# Times the traditional one-process-per-molecule loop (A) against --gpu-batch
# (B) using the SAME xtb binary. On the CPU build this isolates the batch
# driver's amortization of process startup, parameter-file load and OpenMP team
# spin-up (the per-molecule diagonalization is identical). On a GPU build it also
# captures the eigensolver speedup once the batched solve is wired in.
#
# Usage:
#   XTB=build/xtb SRC=/path/to/xyz_dir [GFN=0] [COPIES=1] benchmark/gpu_batch_bench.sh
#
#   XTB     path to the xtb binary           (default: build/xtb)
#   SRC     directory containing *.xyz        (default: ./test)
#   GFN     GFN method                        (default: 0)
#   COPIES  replicate the set this many times (default: 1, to grow the workload)
set -u

XTB=${XTB:-build/xtb}
SRC=${SRC:-./test}
GFN=${GFN:-0}
COPIES=${COPIES:-1}

[ -x "$XTB" ] || { echo "error: xtb binary not found/executable: $XTB" >&2; exit 1; }
[ -d "$SRC" ] || { echo "error: structure dir not found: $SRC" >&2; exit 1; }

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
i=0
while [ "$i" -lt "$COPIES" ]; do
  i=$((i+1))
  for f in "$SRC"/*.xyz; do
    [ -e "$f" ] || continue
    cp "$f" "$WORK/${i}_$(basename "$f")"
  done
done
cd "$WORK" || exit 1

N=$(ls ./*.xyz 2>/dev/null | wc -l)
[ "$N" -gt 0 ] || { echo "error: no .xyz files under $SRC" >&2; exit 1; }
echo "structures: $N  (GFN$GFN)"

echo "--- A: traditional (one xtb process per molecule) ---"
A0=$(date +%s.%N)
for f in ./*.xyz; do "$XTB" --gfn "$GFN" "$f" >/dev/null 2>&1; done
A1=$(date +%s.%N)
A=$(awk "BEGIN{print $A1-$A0}")

echo "--- B: --gpu-batch (one process, all molecules) ---"
B0=$(date +%s.%N)
"$XTB" --gfn "$GFN" --gpu-batch ./*.xyz >/dev/null 2>&1
B1=$(date +%s.%N)
B=$(awk "BEGIN{print $B1-$B0}")

echo "================ RESULT ================"
awk "BEGIN{
  printf \"structures   : %d\n\", $N;
  printf \"A per-process: %8.2f s   (%6.2f mol/s)\n\", $A, $N/$A;
  printf \"B --gpu-batch: %8.2f s   (%6.2f mol/s)\n\", $B, $N/$B;
  printf \"speedup A/B  : %6.2fx\n\", $A/$B;
}"
echo "========================================"
