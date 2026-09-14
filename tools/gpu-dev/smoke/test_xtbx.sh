#!/usr/bin/env bash
H=/mnt/e/Prasanna/xTB/win/xtbx_run.sh
cd /mnt/e/Prasanna/xTB/win/test
rm -rf results

echo "############ MODE 1: single molecule (xtbx water.xyz --gfn 0) ############"
bash "$H" water.xyz --gfn 0 2>&1 | grep -aE "TOTAL ENERGY|normal termination" | head -2

echo
echo "############ MODE 2: folder, CPU parallel (xtbx . --gfn 0) ############"
bash "$H" . --gfn 0 2>&1 | grep -aE "processing|done:"
echo "summary.csv:"; cat results/summary.csv 2>/dev/null
rm -rf results

echo
echo "############ MODE 3: folder, GPU screen (xtbx . --gfn 0 --gpu) ############"
bash "$H" . --gfn 0 --gpu 2>&1 | grep -aE "GPU cuSolver|diagonalized|ammonia|methane|water|processed" | head -8
