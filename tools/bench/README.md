# GPU benchmark / characterization scripts

Reproducible performance characterization for the GPU build (`build-gpushim`,
the cuSolver SCF shim) vs the CPU build (`build`). Used to establish the
CPU↔GPU crossover and the memory ceiling.

> Paths inside are machine-specific (`/mnt/e/Prasanna/xTB/...`, the NVHPC SDK
> location). Adjust the header of each script for another machine. They run in
> WSL and assume `build/xtb` (CPU) and `build-gpushim/xtb` (GPU) exist.

## Scripts
- `build_gpu.sh` — rebuild `build-gpushim/xtb` (serial `ninja` to avoid the
  Fortran module-ordering race), with the NVHPC CUDA libs on `LD_LIBRARY_PATH`.
- `taxol_bench.sh` — taxol (113 atoms) GFN2 single point, CPU vs GPU at OMP 1/8.
  Small-molecule case: CPU wins (GPU launch overhead dominates).
- `big_bench.sh N` — NxNxN water cluster (N=6 → 648 atoms) GFN1 & GFN2 single
  point, CPU vs GPU. Large-system case: GPU wins ~10–32×.
- `opt_bench.sh` — GFN2 `--opt` time-per-cycle across pocket sizes
  (taxol + water clusters), bounded to 8 cycles. Establishes the optimization
  crossover (~350 atoms).

## Key results (RTX 3050, 6 GB)
| system | n (basis) | result |
|---|---|---|
| taxol 113 atoms, GFN2 sp | ~300 | CPU faster (GPU launch-bound) |
| 648-atom cluster, GFN2 sp | ~1296 | **GPU 9.5×** |
| 648-atom cluster, GFN1 sp | ~1296 | **GPU 32×** |
| GFN2 `--opt` | — | crossover ~350 atoms; GPU wins above |

- VRAM: ~173 bytes/n² measured → dense-GFN2 ceiling ≈ **n ≈ 5,500 basis
  functions** (~1,500–2,800 atoms by composition) on a 6 GB card.
- In a *warm* optimization cycle: SCF ≈ 47%, analytical gradient ≈ 42%
  (the gradient is the next offload target; see GPU_PHASE2_BATCHED_SCF.md).

The numerical correctness gate is `test/gpu/gfn12_gpu_gate.py`
(run via `test/gpu/run_gfn12_gpu_gate.sh`).
