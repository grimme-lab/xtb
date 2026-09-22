# GPU development helpers

This directory contains machine-specific scripts used while developing and
validating the xTB CUDA backend. They are not part of the normal xTB command
surface or release build.

Canonical GPU code and tests live in the main repository:

- `src/gpu/` — CUDA and Fortran GPU implementation
- `test/gpu/cusolver_gate.f90` — cuSolver versus LAPACK numerical gate
- `test/gpu/run_gpu_gate.sh` — builds and runs the numerical gate
- `test/gpu/run_gfn12_gpu_gate.sh` — GFN1/GFN2 GPU validation
- `benchmark/gpu_batch_bench.sh` — maintained GPU batch benchmark

The helper folders here are retained for development history:

- `build/` — local WSL build recipes
- `diagnostics/` — compiler, CUDA, and NVIDIA HPC SDK probes
- `benchmarks/` — experimental or machine-specific benchmarks
- `smoke/` — ad hoc integration checks for local wrappers
- `prototype-shim/` — the original standalone CUDA-C proof of concept

Most helpers assume the local checkout is `/mnt/e/Prasanna/xTB/xtb` and the
NVIDIA HPC SDK is installed below `/opt/nvidia/hpc_sdk`. Prefer the canonical
repository tests above for current validation.
