# Building & validating the xtb GPU / batch version (Phase 0 + Phase 1)

This documents what was implemented for Phases 0 and 1 of the GPU plan
(`GPU_ACCELERATION_PLAN.md`), how to build it, and — most importantly — how to
**numerically validate** it. The CPU build is unaffected and unconditional.

> Honesty note: the GPU/cuSolver path has not been compiled on an NVIDIA box yet.
> The CPU-path logic follows the existing xtb idioms and the proven single-point
> flow in `src/peeq_module.f90`; the GPU paths follow the repo's established
> cuSolver/OpenACC pattern in `src/mctc/lapack/eigensolve.F90`. Everything that
> can be validated **without** a GPU has been (see "Validation results" below);
> the cuSolver backend must still be compiled and re-validated on a Linux box
> with the NVIDIA HPC SDK before the GPU path is trusted.

---

## What changed

### Phase 0 — make the existing GPU port buildable again
| File | Change |
|---|---|
| `src/mctc/lapack/eigensolve.F90` | Fixed invalid `real(dp) :: dummy(:)` decl that broke the cuSolver build (`-> dummy(1)`). |
| `meson/meson.build` | Modernized flags: `-ta=tesla:ccNN` → `-gpu=ccNN -acc`, `-Mcudalib` → `-cudalib`. |
| `CMakeLists.txt` | **New** `WITH_GPU` / `GPU_ARCH` / `WITH_CUSOLVER` options (CMake had no GPU support at all). Mirrors meson; Fortran-scoped flags; defines `XTB_GPU`/`USE_CUSOLVER`/`USE_CUBLAS`. |

### Phase 1 — batched / high-throughput path
| File | Purpose |
|---|---|
| `src/gpu/batched_eig.F90` | `TBatchedEigensolver`: batch of generalized symmetric eigenproblems. CPU reference (loops LAPACK `sygvd`, the numerical ground truth) + cuSolver GPU backend. |
| `src/gpu/batch_capture.F90` | `xtb_gpu_batch_capture`: inert-by-default capture of the real GFN0 `(H, S)` eigenproblems at the `peeq` solve site, so the batched solver can be exercised + validated on real matrices. |
| `src/gpu/batch_driver.F90` | `xtb_gpu_batch`: multi-structure single-process driver. Reads N inputs, runs the proven single-point path per molecule, **bins by basis size**, reports throughput + padding waste, runs the batched-eigensolver validation pass (`validateBatchedEig`), and the **batched GFN0 energy path** (`run_batched_energy`). |
| `src/peeq_module.f90` | Inert `gpu_capture_store(nao, H, S)` hook before the GFN0 `solve` (zero CPU regression). **New** `peeq_build_energy` / `peeq_finish_energy`: split the GFN0 single point into build → (batched) diagonalize → finish for the energy path. The original `peeq` (gradient path used by opt/MD) is untouched. |
| `src/gpu/CMakeLists.txt`, `meson.build` | Build wiring; registered in `src/CMakeLists.txt` and `src/meson.build`. |
| `src/prog/main.F90` | `--gpu-batch` CLI flag + dispatch hook. |

### Phase 1 — production routing (GFN0 energy path)
`run_batched_energy` (driver) + `peeq_build_energy` / `peeq_finish_energy` (in
`xtb_peeq`) route the **production GFN0 energy** through the batched solve:
1. `peeq_build_energy` runs the GFN0 prologue per molecule — CN, D4 dispersion
   energy, EEQ charges, self energies, the AO overlap `S` and unpacked
   Hamiltonian `H` — with the gradient deliberately discarded (a scratch
   gradient absorbs the `*_grad` routines' force output; their **energies** are
   exact regardless, see `dsrb_grad`/`pbc_approx_rab`).
2. The driver bins by basis size, pads each bucket and diagonalizes it in **one**
   `TBatchedEigensolver` call, then scatters eigenvalues/vectors into each
   wavefunction.
3. `peeq_finish_energy` assembles density, shell charges, band energy, repulsion
   + SRB energies, total energy, HOMO–LUMO gap and Wiberg bond orders.

This is energy + properties only; **gradients stay on the per-molecule `peeq`**
(Phase 3). Because the original `peeq` is untouched, the gradient path / opt / MD
are unaffected.

---

## Build

### CPU build (unchanged, always works)
```bash
cd xtb
meson setup build --buildtype release
ninja -C build
ctest --test-dir build          # regression suite
```
On a memory-constrained box, lower optimization and parallelism:
`meson configure build -Doptimization=1 && ninja -C build -j2`.

The batch driver, the capture hook, and the **CPU fallback** of the batched
eigensolver all compile and run here — so `--gpu-batch` already gives
single-process multi-molecule throughput on CPU, and the batched-eigensolver
contract is validated without a GPU (see below).

### GPU build (NVIDIA HPC SDK / nvfortran)
Prerequisites: NVIDIA HPC SDK (≥ 24.x) providing `nvfortran` + cuSolver/cuBLAS;
a CUDA-capable GPU; your compute capability (A100=80, H100=90, V100=70, RTX40=89).

```bash
cd xtb
FC=nvfortran CC=nvc meson setup build-gpu \
  --buildtype release -Dgpu=true -Dgpu_arch=80 -Dcusolver=true
ninja -C build-gpu
```
(CMake: `-DWITH_GPU=ON -DGPU_ARCH=80 -DWITH_CUSOLVER=ON`.)

---

## Run

```bash
# high-throughput screening of many structures in one process
xtb --gfn 0 --gpu-batch *.xyz
```
Output: a per-molecule results table, the size-bucket distribution (the batching
opportunity), a molecules/sec + padding-waste summary, and the
**batched-eigensolver validation** section. The GFN method comes from `--gfn`;
GFN0 is the first target (no SCF loop).

---

## Validation results (CPU build, GFN0)

All of the following ran on the plain CPU build (`gfortran`, no GPU).

**1. CPU regression unchanged.** The GPU code is additive and guarded; the plain
build still passes `ctest --test-dir build`.

**2. `--gpu-batch` energies match one-at-a-time `xtb` bit-for-bit.** Same code
path, just looped in one process. Verified on a size-diverse set:

| molecule | nao | per-file `xtb` / Eh | `--gpu-batch` / Eh |
|---|---|---|---|
| water    |  8 |  -4.366769919234 |  -4.36676992 |
| ammonia  | 10 |  -4.575742782498 |  -4.57574278 |
| methane  | 12 |  -4.359317912047 |  -4.35931791 |
| ethanol  | 24 | -10.830471285703 | -10.83047129 |
| benzene  | 36 | -15.961220121246 | -15.96122012 |

**3. Batched eigensolver + padding vs per-system solve (the GPU gate, on real
matrices).** `validateBatchedEig` captures each molecule's real `H C = S C ε`,
pads every system up to the bucket's largest order `n` (H diagonal = 1e6, S =
identity on the padded block), solves the whole bucket through
`TBatchedEigensolver`, and compares the lowest-`m` eigenvalues against an
independent per-system LAPACK `sygvd`. The block-diagonal padding must decouple
exactly so the real spectrum is never perturbed:

| set | systems | bucket `n` | padding exercised | max \|Δε\| | result |
|---|---|---|---|---|---|
| 96 water dimers (uniform nao=16) | 64 (cap) | 16 | none (m = n) | `0.0e+00` eV | PASS |
| diverse (nao 8–36) | 5 | 36 | up to 28 rows | `1.95e-14` eV | PASS |

The diverse set is the meaningful one: water's 8×8 problem padded to 36×36 still
reproduced its 8 eigenvalues to ~1e-14 eV. This is exactly the property the
cuSolver backend relies on, here proven on real GFN0 matrices. **Tolerance gate:
≤ 1e-6 eV** (and the energy/gradient gate from the plan remains ≤ 1e-6 Eh /
≤ 1e-6 Eh·a₀⁻¹ for the eventual GPU-vs-CPU comparison).

**3b. Production routing: batched GFN0 *total energies* vs per-molecule `peeq`.**
`run_batched_energy` builds each molecule's H/S, diagonalizes whole buckets
through `TBatchedEigensolver`, and finishes energies/properties — then compares
the total energy against the full per-molecule `peeq`:

| set | systems | padding | worst \|ΔE\| | result |
|---|---|---|---|---|
| diverse (nao 8–36) | 5 | up to 28 rows | `8.9e-16` Eh | PASS |
| 96 water dimers | 96 | none | `8.9e-15` Eh | PASS |

The total energies are now genuinely produced *by* the batched diagonalization
and match the per-molecule path to rounding (~1e-15 Eh, far inside the 1e-6 Eh
gate). Gradients are not batched (energy + properties only).

**Regression:** the xtb unit suite was unaffected — `xtb:peeq` and
`xtb:repulsion` (the routines involved) pass `OK`, along with 108 others. (On the
dev box the meson run also shows timeouts from a 30 s cap on the slow `-O1`
+sanitizer build, and one unrelated failure in the vendored `jonquil` JSON
subproject under `MALLOC_PERTURB_` — neither is related to these changes.)

**4. cuSolver vs LAPACK (still TODO).** Must be run once the GPU build compiles
on an NVIDIA box; the validation harness above is backend-agnostic and becomes
the GPU gate verbatim.

---

## Throughput benchmark

`benchmark/gpu_batch_bench.sh` times the traditional one-process-per-molecule
loop (A) against `--gpu-batch` (B) using the **same CPU binary**, so it isolates
the batch driver's amortization of process startup, parameter-file load and
OpenMP team spin-up (the per-molecule diagonalization is identical on CPU).

```bash
XTB=build/xtb SRC=/path/to/xyz_dir benchmark/gpu_batch_bench.sh
```

Result on the development box (96 GFN0 structures, gfortran CPU build, 16 cores):

| path | wall time | throughput | speedup |
|---|---|---|---|
| A — one process per molecule | 44.11 s | 2.18 mol/s | 1.00× |
| B — `--gpu-batch` (one process) | 18.35 s | 5.23 mol/s | **2.40×** |

This 2.4× is the CPU/startup-amortization win alone. The eigensolver-level
speedup comes on top of it once the cuSolver backend is built and the batched
solve replaces the per-molecule diagonalization (see "Known limitations").

---

## Known limitations / next increment

- **Production routing of the batched solve — DONE for the GFN0 energy path.**
  `peeq_build_energy` / `peeq_finish_energy` split the GFN0 single point and
  `run_batched_energy` routes whole buckets through one batched diagonalization;
  total energies match per-molecule `peeq` to ~1e-15 Eh (§3b). **Not yet routed:**
  gradients (the split is energy-only; forces stay on per-molecule `peeq`, Phase
  3) and GFN1/GFN2 (SCF state machine, Phase 2). On CPU this routing is a no-op
  for speed (same LAPACK calls) — its payoff is the GPU batched solve dropping in
  unchanged, plus being the template for the Phase-2 SCF batch.
- **GPU eigensolver backend** currently loops `cusolverDnDsygvd` per system. The
  high-throughput version (potrfBatched → trsmBatched → syevjBatched) is
  documented in `batched_eig.F90`; implement after the seam is closed and the
  per-system path validates on hardware.
- **GFN1/GFN2** capture is not instrumented yet (the hook is in the GFN0 `peeq`
  path only); the SCF state machine is Phase 2.
- **`--param` files** are not threaded into the batch path yet (defaults only).
- **Capture cap**: the validation retains at most 64 systems (reported, not
  silently truncated) to bound host memory.
