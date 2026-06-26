# CUDA-C analytical-gradient kernel — implementation plan

Goal: move the GFN1/GFN2 analytical-gradient hotspot onto the GPU via a CUDA-C
kernel in the shim (`src/gpu/gpu_eig.cu`), linked through `iso_c_binding` like
the SCF solver — the reliable nvcc+gfortran toolchain, no OpenACC.

## Why
Profiling a *warm* GFN2 optimization cycle (375 atoms, 6-iter SCF): SCF ≈ 47%,
**analytical gradient ≈ 42%**. The gradient's dominant cost is `build_dSDQH0`
(overlap/dipole/quadrupole integral-derivative contraction). The SCF
diagonalization is already on the GPU (cuSolver resident session); the gradient
is the next lever for fast pocket/protein optimization.

## Reference
`src/xtb/hamiltonian_gpu.f90 :: build_dSDQH0_gpu` is already a GPU-structured
(OpenACC) version of the algorithm — parallel over atom pairs, with the integral
derivative math inline (`sdq`, `sdqg`, `dtrf`, multipole derivatives). It does
NOT run under gfortran (nvfortran/managed-memory clauses), but it is the exact
algorithm to translate to CUDA-C. The CPU `build_dSDQH0` (`hamiltonian.F90`) is
the bit-exact oracle.

## Interface (C, called from scf_module via xtb_gpu_runtime)
```
int gpu_build_dSDQH0(int nat, int nao, int nbf,
                     const int* nShell, const int* at, const double* xyz,
                     /* basis tables */ const int* caoshell, const int* saoshell,
                     const int* nprim, const int* primcount,
                     const double* alp, const double* cont,
                     /* Hamiltonian data (hData fields, flattened) */ ...,
                     const double* selfEnergy, const double* dSEdcn,
                     double intcut,
                     const double* P, const double* Pew,            /* densities */
                     const double* vs, const double* vd, const double* vq, /* AES pot */
                     double* g /*3*nat*/, double* sigma /*9*/, double* dhdcn /*nat*/);
```
Returns 0 on success; non-zero ⇒ caller uses the CPU `build_dSDQH0` (fallback).

## Parallelization
- One CUDA block per atom pair (iat>jat); shells/primitives looped inside.
- Per-pair scratch (`sdq`, `sdqg`, the 6×6 shell blocks) in registers/shared mem.
- Force accumulation into `g(3,nat)` and `sigma(9)` via `atomicAdd` (matches the
  OpenACC `!$acc atomic`); `dhdcn` likewise.
- `dtrf` (cartesian→spherical) and `h0scal` ported as `__device__` helpers.

## Incremental milestones (each gated bit-exact before the next)
1. **Scaffold** — C entry + `xtb_gpu_runtime` wrapper + `scf_module` routing
   behind a runtime flag, kernel stubbed to return non-zero ⇒ CPU fallback.
   Verifies plumbing without changing results. ✅ low risk.
2. **Overlap-derivative term only** on GPU (the `Pew`/`P`·dS Pulay force),
   everything else CPU; validate gradient ≤1e-8 vs CPU on water + taxol.
3. **+ H0/self-energy derivative + shellPoly** terms; re-validate.
4. **+ dipole/quadrupole integral derivatives** (GFN2 AES `vs/vd/vq` coupling);
   re-validate GFN2.
5. **Full kernel**; gate (grad ≤1e-6, opt Δx≈0) + opt-cycle benchmark.

## Validation
`test/gpu/gfn12_gpu_gate.py` (extend with an explicit `--grad` max-abs check at
1e-6). Bit-exactness is unrealistic (atomic-add reordering); target ≤1e-8 per
term early, ≤1e-6 overall — within the existing gate tolerance.

## Risk / honesty
This is a ~250-line numerical kernel with subtle Gaussian-integral-derivative
math. It must be built and validated **one milestone at a time** (each needs a
GPU rebuild + gate run); shipping it as a single unvalidated blob is not
acceptable. Until milestone 5 passes, the CPU gradient stays the default and the
kernel is opt-in behind a flag.
