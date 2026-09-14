# Phase 2 — Batched GFN1/GFN2 SCF on GPU (design)

Goal: a **general** GPU path for the self-consistent methods (GFN1, GFN2), not
tied to any one dataset — process many molecules' SCF iterations together so the
GPU stays saturated. This is the hard part of the GPU plan; this document is the
concrete, code-grounded design + an incremental, *validatable* path.

## Why GFN1/2 are different from GFN0

GFN0 (`xtb_peeq`) is one build + one diagonalization → already batched on GPU
(`--gpu`). GFN1/2 are **self-consistent**: the SCF loop in
`scc` (`src/scc_core.f90:255`, `scc_iterator: do iter = 1, thisiter` at line 434)
repeats, per molecule, until convergence:

```
addShift(q,qsh -> potential)                      ! charges -> potential
buildIsotropicH1 / buildIsoAnisotropicH1          ! charge-dependent H   (line 458/461)
solver%fact_solve(H, S_factorized, emo)           ! diagonalize          (line 475)
dmat(focc, C -> P)                                ! density              (line 513)
mpopsh(S,P -> qsh)                                ! Mulliken charges     (line 516)
Broyden mixing of qsh (+ GFN2 multipoles)         ! charge update
convergence test (econverged .and. qconverged)
```

Different molecules need different iteration counts, and the solver is
`fact_solve` (caches the overlap Cholesky once, reuses it every iteration) — not
the GFN0 `solve`. So batching means a **per-molecule SCF state machine** that
steps all molecules in lockstep with an `active`/converged mask.

## Algorithm: lockstep batched SCF with an active mask

For one size-bucket of `B` molecules (padded to order `n`):

```
init each molecule: H0, S, S-factor, initial qsh, Broyden state, energies
active(1:B) = .true.
for iter = 1 .. maxiter while any(active):
    # all on the batch, skipping inactive molecules
    batched build  H1_k from qsh_k            (buildIsotropicH1 / AES)   [GPU/OMP]
    batched solve  H1_k C_k = S_k C_k eps_k                              [GPU cuSolver, have it]
    batched dmat   P_k = C_k focc_k C_k^T      (cublasDgemmStridedBatched)[GPU]
    batched mpopsh qsh_k                                                  [GPU/OMP]
    per molecule:  Broyden mix, test convergence, set active_k
finalize each: assemble GFN1/2 energy terms, gap, properties
```

GPU wins because every heavy step (diag, density GEMM, H build) runs once for the
whole batch per iteration, and converged molecules drop out (mask), keeping the
batch dense.

## State the batch must carry (per molecule)

From `scc`'s arguments/locals: `qsh, q, dipm, qp` (GFN2 multipoles), `P`, `H/H0`,
`S` + `S_factorized`, `emo`, `focc`, Broyden history (`broyden` module / `nbr`),
the running energy terms, and `active/converged` flags. Wrap in a
`TScfBatchState` (array-of-struct or struct-of-arrays).

## Refactor required (single source of truth)

`scc` is monolithic (SCF loop inside). Split it into:
- `scc_init`   — everything before the loop (setup, S factorization, Broyden init)
- `scc_step`   — **one** SCF iteration (the body above), returns convergence flags
- `scc_final`  — energy assembly + properties after convergence

Then today's `scc` becomes `scc_init; do ...; scc_step; ...; scc_final` —
behaviourally identical, so the regression suite validates the refactor before
any batching. The batched driver calls `scc_init` ×B, then loops `scc_step`
across the active set (with batched kernels), then `scc_final` ×B.

## Increments (each independently validated vs CPU, gate: energy ≤ 1e-6 Eh)

1. **Refactor `scc` → init/step/final**, no behaviour change. Gate: full `meson
   test` unchanged; GFN1/GFN2 energies bit-identical. (Touches core SCF — do
   first, in isolation.)
2. **Batched GFN1 driver, CPU backend.** Drive `scc_step` over a bucket in
   lockstep with an `active` mask, CPU diag (`TBatchedEigensolver`). Gate:
   batched GFN1 energies == per-molecule `scc` for a diverse set.
3. **GPU diag + batched `dmat`.** Route step-2's per-bucket diagonalization
   through the cuSolver shim; batch the density GEMM
   (`cublasDgemmStridedBatched`). Gate: unchanged energies; benchmark vs CPU.
4. **Batched H1 build (GFN1 isotropic).** Make `buildIsotropicH1` gang over
   molecules (or OMP across the batch). Gate + benchmark.
5. **GFN2 extras.** Anisotropic/multipole electrostatics (AES: `mmompop`,
   `aniso_electro` in `aespot.F90`) batched, self-consistent D4 inside the loop.
   Gate vs CPU GFN2.
6. **Gradients (Phase 3 overlap)** — separate, later.

## Honest scope + where it pays off

- Effort: weeks (this is the plan's Phase 2, 6–10 wk estimate). Increment 1 alone
  is a careful core-SCF refactor.
- GPU wins for **large systems** (n ≳ 120, where each diag/GEMM dominates) and
  for **big batches of medium molecules**. For tiny molecules (n ≈ 16–50) the
  per-iteration non-diag work and launch overhead dominate — parallel-CPU
  processes (current `xtbx`) will stay competitive or better there.
- Risk: increment 1 touches `scc`, used by every GFN1/2 run. Mitigation: keep the
  thin wrapper behaviourally identical and gate on the regression suite.

## Status

**In progress (2026-06-24).**

- Increment 1A is implemented: the existing `scc` entry point now delegates to
  three behavior-preserving internal seams, `scc_init`, `scc_step` (exactly one
  SCF iteration), and `scc_final`. The public call contract and numerical order
  are unchanged.
- Increment 1B is implemented: all data that must survive between SCF
  iterations now lives in an explicit `TScfBatchState`. This includes the
  factorized overlap, Broyden matrices/history, GFN2 potential intermediates,
  damping/convergence values, iteration counters, and an `active` flag. The
  legacy `scc` wrapper owns one state; the batched driver can next own an array
  of states.
- GFN1 and GFN2 taxol single-points were compared with an independently linked
  pre-refactor binary. All printed SCF iteration rows, total energies, and
  HOMO-LUMO gaps were identical:
  - GFN1: `-195.911950900256 Eh`, gap `2.982933771696 eV`
  - GFN2: `-186.500449124213 Eh`, gap `2.501305953618 eV`
- The new `xtb_gpu_runtime` leaf module is wired into both Meson/CMake source
  lists and is now the single owner of the `gpu_use` flag and CUDA eigensolver
  interface used by the batch driver. This removes the duplicate runtime state
  before the core SCC code starts depending on it.
- The production GFN1/GFN2 SCC loop now routes generalized diagonalization
  through cuSolver and density construction through cuBLAS when `--gpu` is set.
  CUDA allocations, solver handles, and workspaces persist across SCF and
  optimization cycles.
- `--gpu` no longer implies `--gpu-batch`, so it works with the normal
  single-point, `--grad`, and `--opt` call paths. Combining `--gpu` with
  `--gpu-batch` enables the same CUDA SCC path for every GFN1/GFN2 structure in
  the multi-file driver.
- A hardware gate on the RTX 3050 passes for both methods:
  - single-point energy: CPU/GPU difference `0.0 Eh` at printed precision
  - gradient: max difference `9.992e-16 Eh/a0` (GFN1),
    `1.313e-16 Eh/a0` (GFN2)
  - loose optimization: both converge in 5 steps; optimized energy difference
    `0.0 Eh`; maximum coordinate difference `0.0 A` (GFN1) and `8e-14 A`
    (GFN2)
  The reproducible gate is `test/gpu/gfn12_gpu_gate.py`.

The state refactor was also checked at `OMP_NUM_THREADS=1` against the
pre-refactor binary: all 204 printed GFN1 and 92 printed GFN2 numerical rows
were identical.

Remaining before the original batched Phase-2 architecture is complete:
lift the step procedure behind a batch-callable system context, compact active
states into true cross-molecule SCF launches, and move GFN2 AES potentials off
the CPU. The current `--opt` path is numerically validated and uses
GPU diagonalization+density at every SCF cycle, while analytical force assembly
is still CPU-side.

## Increment 2 — resident GPU SCF session (performance)

The per-call GPU offload above was correct but *slower* than the CPU for small
molecules: each SCF iteration copied H/S/C/P host<->device through several
separate kernels (~8 n^2 transfers/iter), and the constant overlap S was
re-uploaded ~3x per iteration. cuSolver's per-call generalized solve also
re-factorized S every time.

A resident **SCF session** (`gpu_scf_open/solve/finish/get_vectors/close`,
src/gpu/gpu_eig.cu + xtb_gpu_runtime) now keeps S — and, for GFN1, H0 / matlist
/ ao2sh — resident on the device for the whole SCF. The Hamiltonian and density
stay resident between solve and finish, so per iteration only small vectors
cross PCIe: `shift`/`focc` in, `emo`/`qsh` out, plus `P` out (the host
`electro()` energy needs it). GFN1 rebuilds H1 on-device each cycle (zero n^2
in); GFN2 still uploads the host-built anisotropic H once per cycle. The
converged eigenvectors are fetched back once after convergence for the
gradient's energy-weighted density. `scc_core.f90` opens the session in
`scc_init`, routes solve/finish in `scc_step`, and closes in `scc_final`, with a
full CPU fallback whenever the session does not open or a solve fails.

The gate stays bit-exact (GFN1/GFN2 |dE_sp| = 0, max|dG| ~1e-16, opt dx = 0).

Benchmarks (RTX 3050, GPU uses 8 CPU helper threads):

| system | method | CPU | GPU | speedup |
|---|---|---|---|---|
| taxol, 113 atoms | GFN2 sp | 0.74 s | 3.6 s | 0.2x (CPU wins) |
| water cluster, 648 atoms | GFN2 sp | 95.2 s | 10.1 s | **9.5x** |
| water cluster, 648 atoms | GFN1 sp | 433.8 s | 13.6 s | **32x** |

Crossover is near a few hundred atoms: below it CPU wins (GPU launch overhead
dominates), above it the GPU wins by ~10-32x. The `xtbx` front-end encodes this
— a single large molecule is auto-routed to the GPU, small ones run on the CPU
with a CPU->GPU fallback on failure, and folders advise `--gpu` for throughput.

## Increment 3 — analytical gradient / AES offload (investigation, shelved)

Profiling a *warm* GFN2 optimization cycle (375-atom system, 6-iter SCF) showed
the analytical gradient is ~42% of each cycle (SCF ~47%), so the gradient is
worth offloading for optimization. The dominant gradient cost is
`build_dSDQH0` (overlap/dipole/quadrupole integral-derivative contraction).

A pre-existing OpenACC port exists (`#ifdef XTB_GPU` `!$acc` regions in
`hamiltonian_gpu.f90`, `aespot.F90`, `dftd4.F90`, `intgrad.f90`, `repulsion.F90`,
`property.F90`) and is selected by `scf_module.F90`. The `gpu_acc` meson option
was added to drive it (`-DXTB_GPU -fopenacc -foffload=nvptx-none`), and the
gfortran NVPTX offload toolchain (`gcc-11-offload-nvptx`) was verified working
on the RTX 3050. **However the OpenACC code does not build/run under gfortran**:
it was authored for nvfortran with managed/unified memory. Concrete blockers:
- nested `vector` parallelism in `aespot.F90` (fixed: inner loops -> `seq`);
- data clauses listing a derived type **and** its components together
  (`hData, hData%...`, `dispm, dispm%...`) — gfortran error "mixed component and
  non-component accesses"; needs splitting into manual deep-copy (parent shallow,
  then components) across every enter/exit pair;
- `!$acc routine` with gang/worker/vector in PURE procedures (dftd4);
- module parameters in copy clauses (dftd4 `zeff`);
- `default(present)` regions in `intgrad.f90`/`repulsion.F90` (not even gated by
  `XTB_GPU`) that assume managed memory and would fault without explicit
  enter/exit data.

Verdict: enabling it is a real OpenACC port (rewrite all derived-type data
management for gfortran + validate gcc-11 nvptx deep-copy at runtime), not a free
reuse. **Shelved.** Production GPU builds use `gpu_shim` (cuSolver SCF) only;
`build-gpushim` stays the validated, fast build. The realistic fallback for the
gradient is a CUDA-C kernel for `build_dSDQH0` (like the SCF shim) — a sizeable
effort, deferred pending a decision. The `gpu_acc` scaffolding is left in (off by
default, marked experimental) for whoever completes the port.
