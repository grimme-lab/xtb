#!/usr/bin/env bash
# Unified xtb front-end. One command, dispatches to the best engine:
#   xtbx mol.xyz [opts]               -> single molecule; small -> CPU, large ->
#                                        GPU automatically (with CPU->GPU fallback)
#   xtbx mol.xyz --gpu [--opt ...]    -> force the GPU (any GFN0/1/2 / task)
#   xtbx <folder|files> [opts]        -> parallel per-compound (results/<name>/ + CSV)
#   xtbx <folder|files> --gpu         -> dynamic per-GPU queue
#   xtbx <folder|files> --gpu --gpu-batch
#                                     -> legacy all-at-once GFN0 batch
#
# When does the GPU pay off?  The GPU wins big on LARGE systems (measured: a
# 648-atom cluster runs ~10x (GFN2) to ~32x (GFN1) faster than 8 CPU threads)
# and on HIGH THROUGHPUT (>~100 small molecules keep it saturated). For 1-2
# small molecules the CPU is faster (GPU launch overhead dominates). xtbx applies
# this automatically: a single large molecule is routed to the GPU even without
# --gpu, and any failed CPU run is retried on the GPU. Tune the size cutoff with
# XTB_GPU_AUTO_ATOMS (default 350); set it to 0 to disable size-based auto-GPU.
#
# A path argument is anything that exists as a file/dir (or a glob). Everything
# else (e.g. --gfn 0, --opt, --chrg -1) is passed straight through to xtb.
ROOT=/mnt/e/Prasanna/xTB/xtb
export XTBPATH="$ROOT"
CPU_XTB="$ROOT/build/xtb"
GPU_XTB="$ROOT/build-gpushim/xtb"
CONF="$(dirname "$0")/xtbg.conf"

# Set up the CUDA runtime env (libs + default thread count) for the GPU build.
setup_gpu_env() {
  local NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3 V=13.1 j
  export LD_LIBRARY_PATH="$NVHPC/REDIST/math_libs/$V/targets/x86_64-linux/lib:$NVHPC/math_libs/$V/targets/x86_64-linux/lib:$NVHPC/REDIST/cuda/$V/targets/x86_64-linux/lib:$NVHPC/compilers/lib:${LD_LIBRARY_PATH:-}"
  if [ -z "${OMP_NUM_THREADS:-}" ]; then
    j=$(awk -F= '/^JOBS=/{gsub(/[^0-9]/,"",$2);print $2}' "$CONF" 2>/dev/null)
    case "$j" in ''|*[!0-9]*) j=8;; esac
    export OMP_NUM_THREADS="$j"
  fi
}

# Is a usable GPU build + runtime device present?
gpu_available() { [ -x "$GPU_XTB" ] && nvidia-smi -L >/dev/null 2>&1; }

# Atom count of an .xyz (first line). Returns "" for non-xyz / unknown so the
# caller falls back to the safe (CPU) default.
count_atoms() {
  local n
  n=$(head -1 "$1" 2>/dev/null | tr -d '[:space:]')
  case "$n" in ''|*[!0-9]*) n="";; esac
  printf '%s' "$n"
}

# GFN-FF pre-relaxation: a cheap O(N) force-field optimization to a sane geometry
# BEFORE the expensive GFN2 optimization, cutting the number of GFN2 cycles (and
# giving a fast bulk relax for protein-sized inputs). Runs in a temp dir to avoid
# clobbering the GFN2 output files. Echoes the relaxed .xyz path on success, or
# the original input on failure (so the caller is never worse off). stderr only.
prerelax_geometry() {
  local input tag pre tmp
  input="$(realpath "$1" 2>/dev/null || printf '%s' "$1")"
  tag="$(basename "${input%.*}")"
  pre="$(pwd)/${tag}_gfnff.xyz"
  tmp="$(mktemp -d /tmp/xtbx-prerelax.XXXXXX)" || { printf '%s' "$input"; return; }
  echo "xtbx: GFN-FF pre-relaxation (fast O(N) force field) before GFN2 ..." >&2
  # NOTE: this build's GFN-FF xtbopt.xyz writer corrupts the element-symbol
  # column (coordinates are correct). GFN-FF preserves atom order, so rebuild a
  # clean geometry from the ORIGINAL symbols + the optimized coordinates.
  if ( cd "$tmp" && "$CPU_XTB" "$input" --gfnff --opt > prerelax.out 2>&1 ) \
     && [ -f "$tmp/xtbopt.xyz" ] \
     && python3 - "$input" "$tmp/xtbopt.xyz" "$pre" <<'PY'
import sys, re
orig, opt, out = sys.argv[1:4]
# Element symbols: from the ORIGINAL input (clean). split on '\n' only.
Lo = open(orig, encoding="latin-1").read().split("\n")
n = int(Lo[0].split()[0])
syms = [Lo[2+i].split()[0] for i in range(n)]
# Coordinates: from the GFN-FF output. Its element-symbol column is corrupt with
# arbitrary bytes (incl. newlines), so DON'T rely on line structure -- skip the
# clean count+comment lines, then pull every long-decimal float in order. The
# garbage bytes never form a 14-digit decimal, so this recovers exactly 3N coords.
raw = open(opt, encoding="latin-1").read().split("\n", 2)
body = raw[2] if len(raw) > 2 else ""
nums = re.findall(r'[-+]?\d+\.\d{6,}', body)
if len(nums) != 3*n: sys.exit(1)
with open(out, "w") as o:
    o.write("%d\nxtbx GFN-FF prerelaxed\n" % n)
    for i in range(n):
        o.write("%s %s %s %s\n" % (syms[i], nums[3*i], nums[3*i+1], nums[3*i+2]))
PY
  then
    cp -f "$tmp/prerelax.out" "$(pwd)/${tag}_gfnff.out" 2>/dev/null
    rm -rf "$tmp"
    echo "xtbx: pre-relaxed geometry -> $pre" >&2
    printf '%s' "$pre"
  else
    rm -rf "$tmp"
    echo "xtbx: GFN-FF pre-relax failed; continuing from the original geometry" >&2
    printf '%s' "$input"
  fi
}

# ------------------------------- pretty terminal -----------------------------
# ANSI styling, an animated/continuous progress bar, and a final energy table.
# Colours are only emitted to a real terminal so redirected output stays clean.
if [ -t 1 ]; then
  C_RST=$'\033[0m'; C_B=$'\033[1m'; C_DIM=$'\033[2m'
  C_GRN=$'\033[32m'; C_RED=$'\033[31m'; C_CYN=$'\033[36m'; C_YEL=$'\033[33m'
  CLR=$'\033[K'; HIDE=$'\033[?25l'; SHOW=$'\033[?25h'
else
  C_RST=; C_B=; C_DIM=; C_GRN=; C_RED=; C_CYN=; C_YEL=; CLR=; HIDE=; SHOW=
fi
SPIN='⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏'; _spin=0
_BLK=' ▏▎▍▌▋▊▉█'   # 1/8-cell fill steps for a smooth (sub-character) bar

# smartbar DONE TOTAL T0 [LABEL] : one in-place animated bar with %/count/ETA.
smartbar() {
  local done=$1 total=$2 t0=$3 label="${4:-}" width=30
  local pct eighths full rem bar="" i frame now el extra=""
  [ "$total" -le 0 ] && total=1
  [ "$done" -gt "$total" ] && done=$total
  pct=$(( done * 100 / total ))
  eighths=$(( done * width * 8 / total ))
  full=$(( eighths / 8 )); rem=$(( eighths % 8 ))
  for ((i=0; i<full; i++)); do bar+="█"; done
  if [ "$full" -lt "$width" ]; then
    [ "$rem" -gt 0 ] && { bar+="${_BLK:$rem:1}"; i=$((full+1)); } || i=$full
    for (( ; i<width; i++)); do bar+=" "; done
  fi
  frame="${SPIN:$_spin:1}"; _spin=$(( (_spin+1) % 10 ))
  now=$(date +%s); el=$(( now - t0 )); [ "$el" -lt 0 ] && el=0
  if [ "$done" -gt 0 ] && [ "$el" -gt 0 ]; then
    local eta=$(( (total-done) * el / done ))
    extra=$(printf '%d/s  ETA %d:%02d' "$(( done / el ))" "$(( eta/60 ))" "$(( eta%60 ))")
  fi
  printf '\r%s%s%s %s%s%s %s%3d%%%s %s%d/%d%s  %s%s %s%s%s' \
    "$C_CYN" "$frame" "$C_RST" \
    "$C_GRN" "$bar" "$C_RST" \
    "$C_B" "$pct" "$C_RST" \
    "$C_B" "$done" "$total" "$C_RST" \
    "$C_DIM" "$label" "$extra" "$C_RST" "$CLR" >&2
}

# print_energy_table CSV : aligned per-compound energy table + summary footer.
print_energy_table() {
  local csv="$1"
  [ -f "$csv" ] || return 0
  awk -F, -v b="$C_B" -v r="$C_RST" -v d="$C_DIM" -v g="$C_GRN" -v rd="$C_RED" -v cy="$C_CYN" '
    NR==1 { next }
    {
      name=$1; gsub(/^"|"$/,"",name); n=split(name,a,/[\/\\]/); name=a[n];
      sub(/\.(xyz|sdf|coord|mol)$/,"",name);
      e=$2; gp=$3; st=$NF;
      gsub(/[ \t\r]/,"",e); gsub(/[ \t\r]/,"",gp); gsub(/[ \t\r]/,"",st);
      rows++; idx[rows]=name; en[rows]=e; ga[rows]=gp; sta[rows]=st;
      if (st=="ok") { okc++; ev=e+0; if(!have||ev<lo){lo=ev; loname=name; have=1} } else failc++;
    }
    END {
      sep=""; for(i=0;i<78;i++) sep=sep"-";
      printf "\n%s  #   %-26s %18s %10s   status%s\n", b, "compound", "energy / Eh", "gap / eV", r;
      printf "%s%s%s\n", d, sep, r;
      for(i=1;i<=rows;i++){
        edisp = (en[i]=="") ? "        --" : sprintf("%.6f", en[i]+0);
        gdisp = (ga[i]=="") ? "    --"     : sprintf("%.3f", ga[i]+0);
        col = (sta[i]=="ok") ? g : rd;
        printf "%4d   %-26.26s %18s %10s   %s%s%s\n", i, idx[i], edisp, gdisp, col, sta[i], r;
      }
      printf "%s%s%s\n", d, sep, r;
      printf "  %s%d ok%s", g, okc+0, r;
      printf ", %s%d failed%s", (failc? rd:d), failc+0, r;
      if(have) printf "    %slowest%s %s%s%s @ %s%.6f Eh%s", cy, r, b, loname, r, b, lo, r;
      printf "\n";
    }
  ' "$csv"
}

# ---- classify arguments into paths vs options ----
declare -a copts batch
# gfn defaults to 2 (xtb's default method) when --gfn is not given, so dispatch
# routing matches the actual calculation. We do NOT inject --gfn into the command;
# xtb itself applies the GFN2 default.
folder=""; has_gpu=0; has_gpu_batch=0; has_heavy=0; has_task=0; has_prerelax=0; gfn="2"; prev=""
for a in "$@"; do
  ta="$a"
  case "$ta" in [A-Za-z]:/*) ta="$(wslpath -a "$ta" 2>/dev/null || printf '%s' "$ta")";; esac
  if [ -d "$ta" ]; then
    folder="$ta"
  elif [ -f "$ta" ]; then
    batch+=("$ta")
  elif [[ "$ta" == *'*'* || "$ta" == *'?'* ]]; then
    for f in $ta; do [ -f "$f" ] && batch+=("$f"); done
  else
    case "$a" in
      --gpu) has_gpu=1;;                                   # dispatch flag (not for CPU xtb)
      --gpu-batch) has_gpu=1; has_gpu_batch=1;;            # explicit legacy all-at-once path
      --prerelax) has_prerelax=1;;                         # GFN-FF pre-opt before GFN2 (xtbx only)
      --opt|--ohess|--hess|--grad|--md|--omd|--metadyn|--modef|--esp|--stm) has_heavy=1; has_task=1; copts+=("$a");;
      --sp|--vip|--vea|--vipea|--vfukui|--vomega) has_task=1; copts+=("$a");;
      *) copts+=("$a");;
    esac
    [ "$prev" = "--gfn" ] && gfn="$a"
  fi
  prev="$a"
done

# expand a folder into its input structures (skip xtb's own outputs)
if [ -n "$folder" ]; then
  shopt -s nullglob
  for f in "$folder"/*.xyz "$folder"/*.sdf; do
    case "$(basename "$f")" in xtbopt.xyz|xtbhess.xyz|xtblast.xyz) continue;; esac
    batch+=("$f")
  done
fi
N=${#batch[@]}

# For a folder / multi-file run with no explicit task, optimize each compound by
# default (so every input geometry yields an optimized structure in its folder).
# Single files keep xtb's single-point default; --gpu-batch is an explicit
# single-point fast-screen and is left alone.
if [ "$has_task" = 0 ] && [ "$has_gpu_batch" = 0 ] && { [ -n "$folder" ] || [ "$N" -gt 1 ]; }; then
  copts+=("--opt"); has_heavy=1; has_task=1
fi

# High-throughput GFN1/GFN2 single-point pool. Each worker is one persistent
# native --gpu-batch process, so CUDA handles/workspaces and xTB parameters are
# reused across many molecules instead of being recreated for every input.
run_gpu_batch_pool() {
  local workers devices_count i idx gpu pool worker_dir rc
  local -a devices worker_pid
  workers=${XTB_GPU_JOBS:-${XTB_JOBS:-8}}
  case "$workers" in ''|*[!0-9]*) workers=8;; esac
  [ "$workers" -gt "$N" ] && workers=$N
  [ "$workers" -lt 1 ] && workers=1

  if [ -n "${XTB_GPU_DEVICES:-}" ]; then
    IFS=', ' read -r -a devices <<< "$XTB_GPU_DEVICES"
  elif [ -n "${CUDA_VISIBLE_DEVICES:-}" ] && [ "$CUDA_VISIBLE_DEVICES" != "-1" ]; then
    IFS=',' read -r -a devices <<< "$CUDA_VISIBLE_DEVICES"
  else
    while IFS= read -r gpu; do [ -n "$gpu" ] && devices+=("$gpu"); done < <(
      nvidia-smi --query-gpu=index --format=csv,noheader 2>/dev/null | sed 's/[[:space:]]//g'
    )
  fi
  [ "${#devices[@]}" -eq 0 ] && devices=(0)
  devices_count=${#devices[@]}

  results="$(pwd)/results"
  mkdir -p "$results/_gpu_workers"
  pool=$(mktemp -d /tmp/xtbx-gpu-pool.XXXXXX) || return 1

  # Round-robin assignment balances mixed-size folders while each worker keeps
  # one CUDA context alive for its complete share of the queue.
  for ((i=0; i<workers; i++)); do
    mkdir -p "$pool/$i"
    : > "$pool/$i/files"
  done
  for ((idx=0; idx<N; idx++)); do
    printf '%s\0' "$(realpath "${batch[$idx]}")" >> "$pool/$((idx % workers))/files"
  done

  echo "processing $N compounds -> $results   (persistent GPU batch pool, $workers workers)"
  echo "GPU devices: ${devices[*]}   (tuned persistent workers; one CUDA context per worker)"

  run_pool_worker() {
    local slot="$1" device="$2" dir="$pool/$1"
    local -a files
    mapfile -d '' -t files < "$dir/files"
    CUDA_VISIBLE_DEVICES="$device" XTB_BATCH_CSV="$dir/results.csv" \
      "$GPU_XTB" --gpu --gpu-batch "${copts[@]}" "${files[@]}" \
      > "$results/_gpu_workers/worker_${slot}.log" 2>&1
    rc=$?
    if [ ! -s "$dir/results.csv" ]; then
      {
        echo "structure,energy_Eh,gap_eV,status"
        for f in "${files[@]}"; do printf '"%s",,,FAILED\n' "$f"; done
      } > "$dir/results.csv"
    fi
    return "$rc"
  }

  for ((i=0; i<workers; i++)); do
    gpu="${devices[$((i % devices_count))]}"
    run_pool_worker "$i" "$gpu" &
    worker_pid[$i]=$!
  done

  # Each worker draws its own \r progress bar into its log; sum the latest
  # "done/total" across all workers to show ONE live aggregate bar on screen.
  pool_done() {
    local f seg n tot=0
    for ((f=0; f<workers; f++)); do
      [ -f "$results/_gpu_workers/worker_${f}.log" ] || continue
      seg=$(tr '\r' '\n' < "$results/_gpu_workers/worker_${f}.log" 2>/dev/null \
            | grep -oE '[0-9]+/[0-9]+' | tail -1)
      n=${seg%/*}; case "$n" in ''|*[!0-9]*) n=0;; esac
      tot=$((tot+n))
    done
    printf '%s' "$tot"
  }
  local t0; t0=$(date +%s)
  printf '%s' "$HIDE" >&2
  smartbar 0 "$N" "$t0" "compounds"
  while :; do
    local alive=0
    for ((i=0; i<workers; i++)); do
      [ -n "${worker_pid[$i]:-}" ] && kill -0 "${worker_pid[$i]}" 2>/dev/null && alive=$((alive+1))
    done
    smartbar "$(pool_done)" "$N" "$t0" "compounds"
    [ "$alive" -eq 0 ] && break
    sleep 0.25
  done
  rc=0
  for ((i=0; i<workers; i++)); do wait "${worker_pid[$i]}" || rc=1; done
  smartbar "$N" "$N" "$t0" "compounds"; printf '%s\n' "$SHOW" >&2

  {
    echo "structure,energy_Eh,gap_eV,status"
    for ((i=0; i<workers; i++)); do
      tail -n +2 "$pool/$i/results.csv"
    done | sort
  } > "$results/summary.csv"

  print_energy_table "$results/summary.csv"
  echo "  ${C_DIM}-> $results/summary.csv${C_RST}"
  case "$pool" in
    /tmp/xtbx-gpu-pool.*) rm -rf -- "$pool";;
    *) echo "WARNING: refusing to remove unexpected pool path: $pool" >&2;;
  esac
  return "$rc"
}

# ---- dispatch ----
# (0) no input -> pass through (e.g. --version, --help)
if [ "$N" -eq 0 ]; then exec "$CPU_XTB" "${copts[@]}"; fi

# (1) single molecule -> full output in the current folder.
#       --gpu              : force the GPU (any GFN0/1/2 and --sp/--grad/--opt).
#       (no flag), large   : route to the GPU automatically -- it is far faster
#                            on large systems and CPU would be impractically slow.
#       (no flag), small   : run on the CPU (faster for small systems), but if
#                            that run fails, retry on the GPU so we still get a
#                            result on hard/large molecules ("dynamic switch").
if [ "$N" -eq 1 ] && [ -z "$folder" ]; then
  MOL="${batch[0]}"
  # Optional GFN-FF pre-relaxation before a GFN2/GFN1 optimization.
  if [ "$has_prerelax" = 1 ] && printf ' %s ' "${copts[*]}" | grep -q -- ' --opt'; then
    MOL="$(prerelax_geometry "$MOL")"
  fi

  if [ "$has_gpu" = 1 ]; then
    setup_gpu_env
    exec "$GPU_XTB" --gpu "$MOL" "${copts[@]}"
  fi

  thresh=${XTB_GPU_AUTO_ATOMS:-350}
  natoms=$(count_atoms "$MOL")
  case "$thresh" in ''|*[!0-9]*) thresh=350;; esac
  if [ "$thresh" -gt 0 ] && [ -n "$natoms" ] && [ "$natoms" -ge "$thresh" ] && gpu_available; then
    echo "xtbx: large system ($natoms atoms >= $thresh) -> GPU (much faster here; CPU would be slow)" >&2
    setup_gpu_env
    exec "$GPU_XTB" --gpu "$MOL" "${copts[@]}"
  fi

  # Small system: CPU first. On failure, fall back to the GPU.
  "$CPU_XTB" "$MOL" "${copts[@]}"; rc=$?
  if [ "$rc" -ne 0 ] && gpu_available; then
    echo "xtbx: CPU run failed (exit $rc) -> retrying on GPU" >&2
    setup_gpu_env
    exec "$GPU_XTB" --gpu "$MOL" "${copts[@]}"
  fi
  exit "$rc"
fi

# (2) GFN0 has a purpose-built cross-molecule path: build H/S in parallel,
#     bucket by AO size, and diagonalize the buckets with one persistent CUDA
#     context. Always use it for multi-file GPU single points; the generic
#     process queue is several times slower for this non-SCF method.
if [ "$has_gpu" = 1 ] && [ "$gfn" = "0" ] && [ "$has_heavy" = 0 ]; then
  setup_gpu_env
  results="$(pwd)/results"
  mkdir -p "$results"
  export XTB_BATCH_CSV="$results/summary.csv"
  exec "$GPU_XTB" --gpu --gpu-batch "${copts[@]}" "${batch[@]}"
fi

# (3) Explicit fast-screening batch (--gpu-batch): submit the whole list to ONE
#     persistent xtb process (CUDA handles/workspaces reused across molecules),
#     then print the energy table. This is the fastest single-point screen but
#     writes only summary.csv -- no per-compound folders / geometries. The
#     default --gpu folder path (below) gives per-compound folders instead.
if [ "$has_gpu" = 1 ] && [ "$has_gpu_batch" = 1 ] && [ "$has_heavy" = 0 ]; then
  setup_gpu_env
  results="$(pwd)/results"
  mkdir -p "$results"
  export XTB_BATCH_CSV="${XTB_BATCH_CSV:-$results/summary.csv}"
  "$GPU_XTB" --gpu --gpu-batch "${copts[@]}" "${batch[@]}"; rc=$?
  print_energy_table "$results/summary.csv"
  echo "  ${C_DIM}-> $results/summary.csv${C_RST}"
  exit "$rc"
fi

# (5) folder / many molecules, per-compound -> results/<name>/ + summary.csv.
#     Choose the engine:
#       --gpu                             -> one pinned worker slot per GPU;
#                                            refill the exact slot as soon as
#                                            that GPU's process finishes.
#       otherwise                          -> parallel-CPU engine.
#     Either way the pool is completion-driven, so a freed slot is filled at once.
ENGINE=("$CPU_XTB"); engine_label="parallel-CPU"
CPU_FALLBACK_GPU=0
if [ "$has_gpu" = 1 ]; then
  setup_gpu_env
  ENGINE=("$GPU_XTB" --gpu); engine_label="GPU dynamic queue"
elif gpu_available; then
  # CPU engine, but a GPU is present: arm the per-compound CPU->GPU fallback and
  # set up the CUDA runtime libs so the fallback can launch.
  setup_gpu_env
  CPU_FALLBACK_GPU=1
  if [ "$N" -ge 100 ]; then
    echo "xtbx: $N compounds on CPU. For high throughput add --gpu (GPU queue)." >&2
  fi
fi

# Build GPU worker slots. CUDA_VISIBLE_DEVICES can contain physical indexes or
# UUIDs; each child sees only its assigned device as CUDA device 0. By default
# keep eight jobs in flight per GPU (measured optimum on the RTX 3050) so CPU
# setup and GPU solves overlap. Slots
# are distributed round-robin across devices and refilled as they complete.
# XTB_GPU_JOBS/XTB_JOBS can override the total slot count.
declare -a gpu_devices slot_gpu slot_pid
if [ "$has_gpu" = 1 ]; then
  if [ -n "${XTB_GPU_DEVICES:-}" ]; then
    IFS=', ' read -r -a gpu_devices <<< "$XTB_GPU_DEVICES"
  elif [ -n "${CUDA_VISIBLE_DEVICES:-}" ] && [ "$CUDA_VISIBLE_DEVICES" != "-1" ]; then
    IFS=',' read -r -a gpu_devices <<< "$CUDA_VISIBLE_DEVICES"
  else
    while IFS= read -r dev; do
      [ -n "$dev" ] && gpu_devices+=("$dev")
    done < <(nvidia-smi --query-gpu=index --format=csv,noheader 2>/dev/null | sed 's/[[:space:]]//g')
  fi
  [ "${#gpu_devices[@]}" -eq 0 ] && gpu_devices=(0)
fi

# queue width (worker slots / jobs in flight)
JOBS=${XTB_JOBS:-}
if [ -z "$JOBS" ]; then
  if [ "$has_gpu" = 1 ]; then
    JOBS=${XTB_GPU_JOBS:-$((8 * ${#gpu_devices[@]}))}
  else
    JOBS=$(awk -F= '/^JOBS=/{gsub(/[^0-9]/,"",$2);print $2}' "$CONF" 2>/dev/null)
  fi
fi
case "$JOBS" in ''|*[!0-9]*) JOBS=$(nproc 2>/dev/null || echo 4);; esac
[ "$JOBS" -gt "$N" ] && JOBS=$N; [ "$JOBS" -lt 1 ] && JOBS=1
export OMP_NUM_THREADS=${XTB_OMP:-1}   # 1 thread/process; the queue provides the parallelism

if [ "$has_gpu" = 1 ]; then
  for ((slot=0; slot<JOBS; slot++)); do
    slot_gpu[$slot]="${gpu_devices[$((slot % ${#gpu_devices[@]}))]}"
  done
fi

results="$(pwd)/results"; mkdir -p "$results"
task_label="single point"; printf '%s\n' "${copts[*]}" | grep -q -- '--opt' && task_label="optimize"
printf '%s┌─ xtbx · %d compounds · %s · %s · %d in flight%s\n' \
  "$C_CYN$C_B" "$N" "$task_label" "$engine_label" "$JOBS" "$C_RST"
if [ "$has_gpu" = 1 ]; then
  printf '%s└─ GPU %s · per-compound output -> %s%s\n' "$C_DIM" "${gpu_devices[*]}" "$results" "$C_RST"
else
  printf '%s└─ per-compound output -> %s%s\n' "$C_DIM" "$results" "$C_RST"
fi

run_one() {
  local f="$1" gpu="${2:-}" input name outdir
  input="$(realpath "$f" 2>/dev/null || readlink -f "$f" 2>/dev/null)"
  if [ -z "$input" ] || [ ! -f "$input" ]; then
    name="$(basename "${f%.*}")"; outdir="$results/$name"
    mkdir -p "$outdir"
    printf '%s,,,FAILED\n' "$name" > "$outdir/.csvline"
    printf 'ERROR: input file not found: %s\n' "$f" > "$outdir/xtb.out"
    return 1
  fi
  name="$(basename "${f%.*}")"; outdir="$results/$name"
  mkdir -p "$outdir"; cp -f "$input" "$outdir/" 2>/dev/null
  if [ -n "$gpu" ]; then
    ( cd "$outdir" && CUDA_VISIBLE_DEVICES="$gpu" "${ENGINE[@]}" "$input" "${copts[@]}" > xtb.out 2>&1 )
  else
    ( cd "$outdir" && "${ENGINE[@]}" "$input" "${copts[@]}" > xtb.out 2>&1 )
  fi
  # dynamic switch: a compound that fails on the CPU engine is retried on the
  # GPU (large/hard molecules) so the folder run still yields a result for it.
  if ! grep -qa "normal termination" "$outdir/xtb.out" 2>/dev/null \
     && [ "$has_gpu" = 0 ] && [ "${CPU_FALLBACK_GPU:-0}" = 1 ]; then
    printf '\n[xtbx] CPU run did not converge -> retrying on GPU\n' >> "$outdir/xtb.out"
    ( cd "$outdir" && "$GPU_XTB" --gpu "$input" "${copts[@]}" >> xtb.out 2>&1 )
  fi
  if grep -qa "normal termination" "$outdir/xtb.out" 2>/dev/null; then
    local e g
    e=$(grep -a "TOTAL ENERGY"  "$outdir/xtb.out" | tail -1 | awk '{print $(NF-2)}')
    g=$(grep -a "HOMO-LUMO GAP" "$outdir/xtb.out" | tail -1 | awk '{print $(NF-2)}')
    printf '%s,%s,%s,ok\n' "$name" "$e" "$g" > "$outdir/.csvline"
  else
    printf '%s,,,FAILED\n' "$name" > "$outdir/.csvline"
  fi
}
launch_slot() {
  local slot="$1" f="$2" gpu="" name
  [ "$has_gpu" = 1 ] && gpu="${slot_gpu[$slot]}"
  name="$(basename "$f")"
  [ "${XTB_QUEUE_TRACE:-0}" = 1 ] && printf '\nqueue: launch slot=%d gpu=%s job=%s\n' "$slot" "${gpu:-cpu}" "$name"
  run_one "$f" "$gpu" &
  slot_pid[$slot]=$!
}

reap_slot() {
  local finished="" slot pid
  if help wait 2>/dev/null | grep -q -- '-p var'; then
    wait -n -p finished 2>/dev/null || true
    for ((slot=0; slot<JOBS; slot++)); do
      [ "${slot_pid[$slot]:-}" = "$finished" ] && { REAPED_SLOT=$slot; return; }
    done
  fi

  # Portable fallback for shells without `wait -n -p`, or if a shell did not
  # return the completed PID. Poll only active children and reap the first done.
  while :; do
    for ((slot=0; slot<JOBS; slot++)); do
      pid="${slot_pid[$slot]:-}"
      if [ -n "$pid" ] && ! kill -0 "$pid" 2>/dev/null; then
        wait "$pid" 2>/dev/null || true
        REAPED_SLOT=$slot
        return
      fi
    done
    sleep 0.05
  done
}

# A background ticker animates the bar continuously (every 0.25 s) from a shared
# completion count, so it keeps moving even while a long single-compound job
# runs. The main loop only updates the count; the ticker owns the drawing.
t0=$(date +%s)
tickf=$(mktemp /tmp/xtbx-tick.XXXXXX); echo 0 > "$tickf"; : > "$tickf.run"
printf '%s' "$HIDE" >&2
( while [ -e "$tickf.run" ]; do
    smartbar "$(cat "$tickf" 2>/dev/null || echo 0)" "$N" "$t0" "compounds"
    sleep 0.25
  done ) &
ticker_pid=$!

# Seed one job per worker slot. Thereafter a completion immediately feeds the
# next queued molecule to that same slot (and therefore the same pinned GPU).
next=0; running=0; done=0
for ((slot=0; slot<JOBS && next<N; slot++)); do
  launch_slot "$slot" "${batch[$next]}"
  next=$((next+1)); running=$((running+1))
done
while [ "$running" -gt 0 ]; do
  REAPED_SLOT=-1
  reap_slot
  slot=$REAPED_SLOT
  [ "$slot" -lt 0 ] && { echo; echo "ERROR: dynamic queue could not identify a completed worker"; exit 1; }
  [ "${XTB_QUEUE_TRACE:-0}" = 1 ] && printf '\nqueue: complete slot=%d gpu=%s\n' "$slot" "${slot_gpu[$slot]:-cpu}"
  slot_pid[$slot]=""
  running=$((running-1)); done=$((done+1)); echo "$done" > "$tickf"
  if [ "$next" -lt "$N" ]; then
    launch_slot "$slot" "${batch[$next]}"
    next=$((next+1)); running=$((running+1))
  fi
done
rm -f "$tickf.run"; kill "$ticker_pid" 2>/dev/null; wait "$ticker_pid" 2>/dev/null
smartbar "$N" "$N" "$t0" "compounds"; printf '%s\n' "$SHOW" >&2
rm -f "$tickf"

{ echo "structure,energy_Eh,gap_eV,status"
  for d in "$results"/*/; do [ -f "$d/.csvline" ] && cat "$d/.csvline"; done | sort; } > "$results/summary.csv"
print_energy_table "$results/summary.csv"
echo "  ${C_DIM}-> $results/summary.csv   (per-compound folders + xtbopt.xyz in $results/<name>/)${C_RST}"
