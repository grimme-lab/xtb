#!/usr/bin/env bash
cd /tmp
echo "=== energy lines, each GPU run ==="
for r in 1 2 3; do
  echo "--- run $r ---"
  grep -aE 'ammonia|methane|water' "xtb_parallel_gpu_run_${r}.out" | grep -avi building
done
echo
echo "=== determinism: unique energy-line sets across the 3 runs ==="
{ for r in 1 2 3; do grep -aE 'ammonia|methane|water' "xtb_parallel_gpu_run_${r}.out" | grep -avi building; done ; } | sort | uniq -c
echo
echo "=== full-precision: diff run1 vs run2, run1 vs run3 (whole files minus progress bar) ==="
clean() { grep -av $'\r' "$1" | grep -avi 'building\|wall time\|throughput'; }
if diff <(clean xtb_parallel_gpu_run_1.out) <(clean xtb_parallel_gpu_run_2.out) >/dev/null && \
   diff <(clean xtb_parallel_gpu_run_1.out) <(clean xtb_parallel_gpu_run_3.out) >/dev/null; then
  echo "ALL THREE GPU RUNS IDENTICAL (deterministic; no race observed)"
else
  echo "RUNS DIFFER -- possible race:"
  diff <(clean xtb_parallel_gpu_run_1.out) <(clean xtb_parallel_gpu_run_2.out) | head
fi
