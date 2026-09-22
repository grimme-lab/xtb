#!/usr/bin/env bash
cd /mnt/e/Prasanna/xTB/win/test
rm -rf results
XTB_QUEUE_TRACE=1 bash /mnt/e/Prasanna/xTB/win/xtbx_run.sh . --gfn 2 --gpu > /tmp/queue.out 2>&1
echo "=== dispatch line ==="
grep -aE "GPU dynamic queue|GPU devices:|queue: (launch|complete)|done:" /tmp/queue.out
echo "=== per-compound: GPU cuSolver marker present in each xtb.out? ==="
for d in results/*/; do
  n=$(basename "$d")
  m=$(grep -ac "GPU cuSolver SCC enabled" "$d/xtb.out" 2>/dev/null)
  echo "  $n : cuSolver-marker=$m"
done
rm -rf results
