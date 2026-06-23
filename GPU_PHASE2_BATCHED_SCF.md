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
Plan only. Recommended start: **Increment 1** (the `scc` init/step/final split),
validated against the regression suite before any GPU work.
