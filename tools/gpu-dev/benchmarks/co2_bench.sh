#!/usr/bin/env bash
# Time benchmark on the real CO2-capture test set across execution strategies.
# Method: GFN0 single point (fast; the parallel speedup ratios transfer to
# GFN2/--opt). Runs in a scratch dir so the user's own results/ is untouched.
SRC="/mnt/e/Prasanna/Research/CO2 Capture/Data/test"
HELP=/mnt/e/Prasanna/xTB/win/xtbfolder_run.sh
XTB=/mnt/e/Prasanna/xTB/xtb/build/xtb
export XTBPATH=/mnt/e/Prasanna/xTB/xtb
WORK=/tmp/co2bench
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"

N=$(ls "$SRC"/*.xyz 2>/dev/null | wc -l)
echo "N=$N"

wall() { date +%s.%N; }
elapsed() { awk "BEGIN{printf \"%.2f\", $2-$1}"; }

# parallel sweep via xtbfolder helper (XTB_JOBS override)
for j in 1 4 8 16; do
  rm -rf "$WORK/results"
  t0=$(wall); XTB_JOBS=$j bash "$HELP" "$SRC" --gfn 0 >/dev/null 2>&1; t1=$(wall)
  echo "JOBS=$j $(elapsed $t0 $t1)"
done

# single-process batch driver (CPU batched-eigensolver path)
t0=$(wall); "$XTB" --gfn 0 --gpu-batch "$SRC"/*.xyz >/dev/null 2>&1; t1=$(wall)
echo "GPUBATCH $(elapsed $t0 $t1)"

echo "DONE"
