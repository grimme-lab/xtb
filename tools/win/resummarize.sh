#!/usr/bin/env bash
# Rebuild summary.csv from existing per-compound xtb.out files (no re-run).
# Usage: resummarize.sh "<results-dir>"
R="$1"
case "$R" in [A-Za-z]:/*) R="$(wslpath -a "$R" 2>/dev/null || printf '%s' "$R")";; esac
[ -d "$R" ] || { echo "ERROR: not a folder: $R"; exit 1; }
{ echo "structure,energy_Eh,gap_eV,status"
  for d in "$R"/*/; do
    [ -f "$d/xtb.out" ] || continue
    name="$(basename "$d")"
    if grep -qa "normal termination" "$d/xtb.out"; then
      e=$(grep -a "TOTAL ENERGY"  "$d/xtb.out" | tail -1 | awk '{print $(NF-2)}')
      g=$(grep -a "HOMO-LUMO GAP" "$d/xtb.out" | tail -1 | awk '{print $(NF-2)}')
      echo "$name,$e,$g,ok"
    else
      echo "$name,,,FAILED"
    fi
  done | sort
} > "$R/summary.csv"
echo "rebuilt: $R/summary.csv  ($(($(wc -l < "$R/summary.csv")-1)) compounds)"
