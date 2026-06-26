#!/usr/bin/env bash
# Scan the system and write the recommended parallel-job count to xtbg.conf.
# Heuristic: use the physical core count (best for compute-bound xtb), but cap
# it so each parallel job has ~0.7 GB RAM headroom.
WINDIR="$(cd "$(dirname "$0")" && pwd)"
CONF="$WINDIR/xtbg.conf"

logical=$(nproc 2>/dev/null || echo 1)
phys=$(lscpu -p=Core 2>/dev/null | grep -v '^#' | sort -u | wc -l)
{ [ "$phys" -ge 1 ]; } 2>/dev/null || phys=$logical
ram=$(free -g | awk '/^Mem:/{print $2}')
{ [ "$ram" -ge 1 ]; } 2>/dev/null || ram=4

# RAM-bound cap at ~0.7 GB/job  (ram / 0.7  ==  ram * 10 / 7)
membound=$(( ram * 10 / 7 ))
[ "$membound" -lt 1 ] && membound=1

jobs=$phys
[ "$jobs" -gt "$membound" ] && jobs=$membound
[ "$jobs" -lt 1 ] && jobs=1

printf 'JOBS=%d\n' "$jobs" > "$CONF"

echo "==================== xtb parallel setup ===================="
echo " logical CPUs   : $logical"
echo " physical cores : $phys"
echo " RAM            : ${ram} GB  (-> at most $membound jobs by memory)"
echo " ----------------------------------------------------------"
echo " parallel jobs  : $jobs        (written to xtbg.conf)"
echo "============================================================"
echo "xtbfolder will now use $jobs compounds in parallel."
echo "Override anytime:  set XTB_JOBS=<n>   (or edit $CONF)"
