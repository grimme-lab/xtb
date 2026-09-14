# Windows / WSL launchers for xtbx

A snapshot of the one-command `xtbx` front-end used to drive this GPU build of
xtb from Windows (the calculation runs in WSL2; the GPU is reached through WSL).

> **Note:** these scripts contain machine-specific absolute paths
> (`/mnt/e/Prasanna/xTB/...`, the NVHPC SDK location, etc.). They are checked in
> as a backup of the working setup, not as a portable install. Adjust the paths
> at the top of `xtbx_run.sh` and the WSL path in `xtbx.cmd` for another machine.

## Files
- `xtbx.cmd` — Windows entry point. Forwards args (and `XTB_*` env vars) into WSL.
- `xtbx_run.sh` — the dispatcher. Picks the engine and output layout:
  - `xtbx mol.xyz` — one molecule; small → CPU, large (≥ `XTB_GPU_AUTO_ATOMS`,
    default 350 atoms) → GPU automatically, with a CPU→GPU fallback on failure.
  - `xtbx mol.xyz --gpu` — force the GPU (any GFN0/1/2, `--sp`/`--grad`/`--opt`).
  - `xtbx <folder> [--gpu]` — per-compound dynamic queue; **optimizes each
    compound by default** and writes `results/<name>/` (with `xtbopt.xyz`),
    a live animated progress bar, a final energy table, and `results/summary.csv`.
  - `xtbx <folder> --gpu --gpu-batch` — fastest single-point screen (one
    persistent process, `summary.csv` + table, no per-compound folders).
- `xtb-setup.cmd` / `xtb-setup.sh` — one-time setup: add the launcher folder to
  the Windows user PATH and auto-detect cores/RAM into `xtbg.conf`.
- `detect_cores.sh`, `resummarize.sh` — helpers (core detection; rebuild
  `summary.csv` from existing `xtb.out` files).

## Useful environment variables
- `XTB_GPU_AUTO_ATOMS` — single-molecule auto-GPU size cutoff (default 350; 0 disables).
- `XTB_JOBS` / `XTB_GPU_JOBS` — queue width (jobs in flight).
- `XTB_OMP` — OpenMP threads per queued process (default 1).
- `XTB_GPU_DEVICES` — comma-separated GPU indices to use.
- `XTB_QUEUE_TRACE=1` — log queue launch/complete events.
