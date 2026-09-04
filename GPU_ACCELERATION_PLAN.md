# xTB GPU Acceleration — Implementation Plan

**Goal:** A new version of xTB that (a) runs faster per molecule on large systems and (b) processes many molecules at high throughput, using NVIDIA GPUs natively.

**Targets chosen:** GFN0 / GFN1 / GFN2-xTB · NVIDIA-only (CUDA / cuSolver / cuBLAS, NVIDIA HPC SDK `nvfortran`) · throughput batching first, then large single systems.

---

## 1. Current state of GPU support in this tree

A **partial OpenACC port from NVIDIA (2020)** already exists but is bit-rotted and single-molecule only:

| Component | File | What it does |
|---|---|---|
| H0 / overlap / multipole integral build | `xtb/src/xtb/hamiltonian_gpu.f90` | OpenACC `build_SDQH0_gpu`, `build_dSDQH0_gpu` |
| Generalized eigensolver (per SCF iter) | `xtb/src/mctc/lapack/eigensolve.F90` | routes `dsygvd` → `cusolverDnDsygvd` |
| AES anisotropic electrostatics | `xtb/src/aespot.F90` | `mmompop_gpu`, `aniso_electro_gpu` |
| D4 dispersion ATM term | `xtb/src/disp/dftd4.F90` | `atm_gradient_latp_gpu` |
| Integral driver | `xtb/src/intgrad.f90` | `sdqint_gpu`, `build_sdq_ints_gpu` |
| Build wiring (meson only) | `xtb/meson/meson.build:85-96`, `xtb/meson_options.txt:71-75` | `-Dgpu=true -Dgpu_arch=NN` |

**Problems with the existing port:**
- Deprecated compiler flags (`-ta=tesla:cc70` → must become `-acc -gpu=ccNN`).
- Wired into **meson only**, not CMake (`xtb/src/xtb/CMakeLists.txt` always compiles the file but no `-DXTB_GPU`).
- **One molecule at a time.** For xTB's typical 20–300-atom molecules, a single dense GPU eigensolve loses to CPU LAPACK (PCIe transfer + kernel-launch overhead, GPU <5% occupied).
- No CI, no validated numerics, likely won't compile as-is.

**Conclusion:** Reuse the existing kernels as a *starting point for the large-molecule path*, but the throughput goal requires a new **batched** architecture that does not exist yet.

---

## 2. Where the time actually goes (cost model)

SCF iteration loop — `xtb/src/scc_core.f90:434-523`. Per iteration, for `nao` basis functions:

1. Build charge-dependent H — `buildIsoAnisotropicH1` / `buildIsotropicH1` — **O(nao²)**
2. **`solver%fact_solve` — generalized symmetric eigensolve — O(nao³) ← dominant**
3. Density matrix `dmat` — GEMM — **O(nao³)**
4. Mulliken pop `mpopsh`, energy `electro`, Broyden mixing — O(nao²) and smaller

Typical: 10–20 SCF iterations. So per molecule ≈ `(2 × O(nao³)) × ~15`.

**The two regimes:**
- **Large molecule** (`nao` ≳ 2000): each eigensolve/GEMM is big → GPU offload of the *single* operation wins. This is what the 2020 port targets.
- **Many small molecules** (`nao` ~ 50–1000): each operation is tiny; the only way to use the GPU is to **run many molecules concurrently** (batched kernels + streams). This is the throughput win and is the bigger lever for screening workloads.

GFN0 is special: **no SCF loop** — one EEQ charge solve + one H build + one diagonalization. Embarrassingly batchable; the highest-throughput target and the best Phase-1 proving ground.

---

## 3. Target architecture

```
                         ┌─────────────────────────────────────────┐
   N molecules  ───────► │  Host scheduler (new)                    │
   (files / API)         │  • parse + basis setup per molecule      │
                         │  • bin by nao into size buckets          │
                         │  • pad-to-bucket, pack into batched tensors
                         └───────────────┬─────────────────────────┘
                                         │  per bucket
              ┌──────────────────────────┼───────────────────────────┐
              ▼                          ▼                            ▼
     Batched H/S build           Batched eigensolver          Batched density
     (port hamiltonian_gpu  →   (cuSolver syevjBatched /      (cublasDgemmBatched /
      batched, one mol/gang)     potrfBatched+trsm+syevj,      strided-batched)
                                  or MAGMA batched, or
                                  custom batched Jacobi)
              └──────────────────────────┬───────────────────────────┘
                                         ▼
                       Per-mol SCF state machine on host
              (Broyden mixing per mol; "active mask" — molecules drop
               out of the batch as they converge; repack remaining)
```

Two code paths share the same kernels:
- **`--gpu` (latency):** one molecule, batch size 1, large-N kernels (modernized 2020 port).
- **`--gpu-batch` (throughput):** many molecules, the scheduler above.

### Key design decisions
- **Mixed precision:** keep SCF charges/eigensolve in fp64 for numerical parity (xTB validates to µHartree). Optionally explore fp32 H/S build with fp64 accumulation later — *not* in v1.
- **Batched generalized → standard transform:** `HC = SCε` → Cholesky `S = LLᵀ` (`cusolverDnDpotrfBatched`), form `L⁻¹ H L⁻ᵀ` (`cublasDtrsmBatched` ×2), solve standard symmetric eig (`cusolverDnDsyevjBatched` — batched Jacobi, ideal for many small matrices), back-transform vectors. S is constant across SCF iters, so **factor S once per molecule**, reuse every iteration (the code already does this via `S_factorized` / `fact_solve`).
- **Size binning + padding:** bucket molecules by `nao` (e.g. ≤64, ≤128, ≤256, …). Pad each matrix to bucket size; set padded diagonal of H to a large value and S to identity so padded eigenvalues sort above the occupied window and never mix in. `log()` the padding waste so it's visible.
- **Convergence divergence:** molecules in a batch finish in different iteration counts. Maintain an `active` mask; compact the batch (drop converged molecules, optionally refill from a queue) to keep the GPU full.
- **Streams:** one CUDA stream per bucket so independent buckets/copies overlap with compute.

---

## 4. Phased implementation plan

### Phase 0 — Resurrect & validate the existing port (1–2 weeks)
Get a known-good single-molecule GPU baseline before building anything new.
- Modernize flags in `xtb/meson/meson.build:85-96`: replace `-ta=tesla:ccNN` with `-acc -gpu=ccNN`; verify `-Mcudalib=cusolver,cublas`.
- Add the same `gpu`/`cusolver` options to the **CMake** build (`xtb/CMakeLists.txt`, `xtb/src/xtb/CMakeLists.txt`) so both build systems work; gate `hamiltonian_gpu.f90` on `XTB_GPU`.
- Build with NVIDIA HPC SDK `nvfortran`; fix compile errors (e.g. `eigensolve.F90:95` `real(dp) :: dummy(:)` is malformed — must be an allocatable/assumed-shape dummy).
- **Numerical gate:** run the existing test suite (`ninja -C build test`) on GPU vs CPU; energies must match to ≤1e-6 Eh, gradients ≤1e-6 Eh/a₀. This is the regression bar for everything after.
- Benchmark single-molecule speedup vs molecule size to find the crossover where GPU beats CPU.

**Deliverable:** a compiling, test-passing `xtb --gpu` on one molecule (latency path foundation).

### Phase 1 — GFN0 batched single-point (3–5 weeks)
GFN0 has no SCF loop → simplest path to real throughput; validates the whole batched pipeline.
- New module `xtb/src/gpu/batch_scheduler.f90`: parse a multi-molecule input (dir of files, or list), build basis per molecule, bin by `nao`, pack batched arrays.
- Batched EEQ solve (linear system per molecule) — `cusolverDnDpotrfBatched` + batched solve, or batched LU.
- Batched H0 build: refactor `build_SDQH0_gpu` (`hamiltonian_gpu.f90`) so the outer OpenACC `gang` loop is **over molecules in the batch** (currently over atom pairs of one molecule).
- Batched eigensolver module `xtb/src/gpu/batched_eig.f90` (potrf→trsm→syevj as in §3).
- Batched property extraction (energies, gradients, Mulliken).
- **CLI:** `xtb --gpu-batch <inputs...>`; emit per-molecule results + a throughput summary.

**Deliverable:** GFN0 screening of thousands of molecules with 1 GPU vs N CPU cores benchmarked. Numerics gated against CPU GFN0.

### Phase 2 — GFN1/GFN2 batched SCF (6–10 weeks)
The hard part: a batched, per-molecule SCF state machine.
- Host SCF driver over a batch: per-molecule Broyden state (`xtb/src/broyden.f90`), `active` mask, batch compaction on convergence.
- Per iteration on GPU: batched H build (`buildIsoAnisotropicH1`/`buildIsotropicH1` ported), batched `fact_solve` (reuse S factorization across iters), batched `dmat` (`cublasDgemmStridedBatched`), batched `mpopsh`.
- **GFN2 extras:** AES multipole electrostatics — port/​batch `mmompop_gpu` + `aniso_electro_gpu` (`aespot.F90`) and self-consistent D4 (`dftd4.F90`); these run inside the SCF loop.
- GFN1 first (isotropic electrostatics, halogen term in `xtb/src/xtb/halogen.f90`), then GFN2 (adds anisotropic/multipole).

**Deliverable:** `--gpu-batch` for GFN1 and GFN2 single-points, numerics-gated.

### Phase 3 — Gradients, opt & MD (4–6 weeks)
Most xTB usage is geometry optimization / MD, not single points.
- Batched gradient: port/batch `build_dSDQH0_gpu` and the D4/AES gradient kernels.
- Batched ANC/L-BFGS optimizer (`xtb/src/lbfgs_anc/`, `xtb/src/relaxation_engine.f90`) — independent optimizers running in lockstep over the batch with an `active`/converged mask.
- Keep geometry + density resident on GPU across opt steps to avoid PCIe round-trips.

**Deliverable:** batched `--opt` (and optionally MD) on GPU.

### Phase 4 — Productization (2–4 weeks)
- **Batch C API** in `xtb/src/api/` (extend `interface.f90`/`calculator.f90`) accepting arrays of molecules → arrays of results.
- Python bindings for the batch API (mirrors the existing `xtb-python` C-API consumer).
- Multi-GPU: shard buckets across devices.
- Docs + benchmark suite + CI job on a GPU runner.

---

## 5. Build & toolchain

- **Compiler:** NVIDIA HPC SDK (`nvfortran`/`nvc`), v24+ recommended. The CPU-only build stays on gfortran/ifx.
- **Libraries:** cuSolver, cuBLAS (and optionally **MAGMA** for batched `dsygvd` if cuSolver's batched coverage is insufficient for generalized problems of varying size).
- **Flags (modernized):** `-acc -gpu=ccNN -Minfo=accel -DXTB_GPU -cudalib=cusolver,cublas -DUSE_CUSOLVER -DUSE_CUBLAS`.
- Add a CMake `XTB_GPU` option mirroring the meson one so both build systems are first-class.
- Keep `XTB_GPU` / `USE_CUSOLVER` macros so a single source tree builds CPU-only or GPU.

---

## 6. Validation & benchmarking (non-negotiable)

- **Numerical parity gate** in CI: GPU vs CPU on the existing test suite — energy ≤1e-6 Eh, gradient ≤1e-6 Eh/a₀, dipoles/charges within tolerance. Run on every PR.
- **Throughput benchmark:** a fixed set (e.g. 10k drug-like molecules from a public set) — molecules/sec and Wh/molecule, GPU-batch vs CPU (1 core and full node).
- **Latency benchmark:** single large molecules across sizes to locate and report the GPU/CPU crossover.
- Track GPU occupancy / padding waste (`log()` it) to guide bucket boundaries.

---

## 7. Risks & mitigations

| Risk | Mitigation |
|---|---|
| Batched **generalized** eig coverage in cuSolver is thin | Transform to standard (potrf+trsm+syevjBatched); fall back to MAGMA batched if needed |
| Convergence divergence wastes GPU as batch empties | `active` mask + batch compaction + refill from a queue |
| Padding waste for mixed sizes | Size bucketing; tune bucket edges from real distributions |
| fp64 throughput limited on consumer GPUs (RTX) | Target data-center GPUs (A100/H100) for fp64; investigate mixed precision later |
| Numerical drift vs CPU breaking validated results | Hard CI parity gate from Phase 0 onward; keep fp64 in SCF |
| `nvfortran`-only narrows the dev/CI surface | Keep CPU build on gfortran/ifx; GPU build is additive behind macros |
| Existing 2020 port doesn't compile | Phase 0 explicitly budgets for fixing it before building new code |

---

## 8. Effort & sequencing

| Phase | Scope | Est. |
|---|---|---|
| 0 | Resurrect + validate existing port (latency baseline) | 1–2 wk |
| 1 | GFN0 batched single-point (throughput pipeline) | 3–5 wk |
| 2 | GFN1/GFN2 batched SCF | 6–10 wk |
| 3 | Batched gradients / opt / MD | 4–6 wk |
| 4 | Batch C API + Python + multi-GPU + CI | 2–4 wk |

**Recommended first move:** Phase 0 → Phase 1. GFN0 batching is the fastest route to a demonstrable "runs faster, accepts more molecules" result and de-risks the batched eigensolver before the harder SCF batching.
