#!/usr/bin/env bash
export XTBPATH=/mnt/e/Prasanna/xTB/xtb
T=/mnt/e/Prasanna/xTB/xtb/build/test/unit/tester
for t in gfn2 gfn1; do
  echo "=== $t (direct, no sanitizer wrapper) ==="
  t0=$(date +%s)
  "$T" "$t" > "/tmp/${t}_test.log" 2>&1
  rc=$?
  t1=$(date +%s)
  echo "exit=$rc  time=$((t1 - t0))s"
  tail -4 "/tmp/${t}_test.log"
  echo
done
