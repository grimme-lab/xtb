#!/usr/bin/env bash
H=/mnt/e/Prasanna/xTB/win/xtbx_run.sh
cd /mnt/e/Prasanna/xTB/win/test

echo "### single molecule GFN2 --gpu --opt  -> GPU build ###"
bash "$H" water.xyz --gfn 2 --gpu --opt 2>&1 | grep -aE "GPU cuSolver SCC enabled|GEOMETRY OPTIMIZATION CONVERGED|TOTAL ENERGY" | head -3

echo
echo "### single molecule GFN2 (no --gpu) -> CPU build ###"
bash "$H" water.xyz --gfn 2 --sp 2>&1 | grep -aE "GPU cuSolver SCC enabled|TOTAL ENERGY" | head -2
echo "(no 'GPU cuSolver' line above = CPU, correct)"

echo
echo "### folder GFN2 --gpu -> note + parallel CPU ###"
rm -rf results
bash "$H" . --gfn 2 --gpu 2>&1 | grep -aE "note:|parallel-CPU|processing|done:" | head -4
rm -rf results

echo
echo "### folder GFN0 --gpu -> GPU batch ###"
bash "$H" . --gfn 0 --gpu 2>&1 | grep -aE "GPU cuSolver|processed" | head -2
