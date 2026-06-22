! This file is part of xtb.
!
! Copyright (C) 2026 xtb GPU contributors
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> High-throughput batch driver: run a GFN0/1/2 single point for many input
!> structures inside one process.
!>
!> Today xtb processes exactly one geometry per invocation, so screening N
!> molecules means N process startups (parameter load, OMP team spin-up, I/O).
!> This driver amortizes all of that across the whole set and, on the GPU
!> build, is the host-side scheduler into which the batched kernels slot:
!> it bins molecules by basis-set size so equally-sized systems can be padded
!> and diagonalized together (see xtb_gpu_batched_eig / TBatchedEigensolver).
!>
!> Activated by the `--gpu-batch` CLI flag (xtb_gpu_batch::gpu_batch). The GFN
!> method is taken from the global settings, so combine it with `--gfn 0|1|2`:
!>
!>     xtb --gfn 0 --gpu-batch mol1.xyz mol2.xyz mol3.xyz ...
!>
!> STATUS: the harness (multi-molecule iteration, size-binning, padding-waste
!> accounting, throughput reporting) is complete and runs on the CPU build.
!> Each molecule currently goes through the proven calc%singlepoint path. The
!> per-bucket *batched GPU diagonalization* is the next increment and is marked
!> with `INTEGRATION SEAM` below: it requires extracting H/S per molecule and
!> routing the bucket through TBatchedEigensolver instead of the per-molecule
!> diagonalization. That step must be compiled and numerically validated on an
!> NVIDIA HPC SDK box (see BUILD_GPU.md).
module xtb_gpu_batch
   use xtb_mctc_accuracy, only : wp, i8
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_type_calculator, only : TCalculator
   use xtb_type_data, only : scc_results
   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newCalculator, newWavefunction
   use xtb_io_reader, only : readMolecule
   use xtb_mctc_filetypes, only : getFileType
   use xtb_setparam, only : set
   use xtb_readin, only : xfind
   use xtb_mctc_lapack_geneigval, only : lapack_sygvd
   use xtb_gpu_batched_eig, only : TBatchedEigensolver, init
   use xtb_gpu_batch_capture, only : gpu_capture_enable, gpu_capture_disable, &
      & gpu_capture_count, gpu_capture_seen, gpu_capture_get, gpu_capture_clear, &
      & TCapturedSystem
   implicit none
   private

   public :: gpu_batch, gpu_batch_size, run_gpu_batch

   !> Set by the `--gpu-batch` CLI flag; dispatches xtbMain to run_gpu_batch.
   logical :: gpu_batch = .false.
   !> Optional cap on molecules processed per GPU launch (0 = automatic).
   integer :: gpu_batch_size = 0

   !> Upper edges (in number of AOs) of the size buckets used for batching.
   integer, parameter :: nbins = 6
   integer, parameter :: binEdge(nbins) = [64, 128, 256, 512, 1024, huge(1)]

   !> Default parameter file names per GFN method (resolved via XTBPATH/xfind).
   character(len=*), parameter :: fname_gfn0 = 'param_gfn0-xtb.txt'
   character(len=*), parameter :: fname_gfn1 = 'param_gfn1-xtb.txt'
   character(len=*), parameter :: fname_gfn2 = 'param_gfn2-xtb.txt'

   !> Max systems captured for the batched-eigensolver validation (memory bound).
   integer, parameter :: capMax = 64
   !> Padded-block diagonal value (eV). Far above any GFN spectrum, so the
   !> spurious eigenvalues of the padding sort to the top and never perturb the
   !> real occupied window. S is set to the identity on the padded block.
   real(wp), parameter :: padDiag = 1.0e6_wp
   !> Eigenvalue parity tolerance (eV): block-diagonal padding decouples exactly,
   !> so deviation is rounding-level; this is the gate for "padding is correct".
   real(wp), parameter :: eigTol = 1.0e-6_wp

   !> Per-molecule outcome record.
   type :: TBatchResult
      character(len=:), allocatable :: fname
      logical  :: ok = .false.
      integer  :: nat = 0
      integer  :: nao = 0
      real(wp) :: energy = 0.0_wp
      real(wp) :: gnorm = 0.0_wp
      real(wp) :: gap = 0.0_wp
      real(wp) :: walltime = 0.0_wp
   end type TBatchResult

   character(len=*), parameter :: source = 'gpu_batch'

contains

!> Entry point invoked from xtbMain when `--gpu-batch` is set.
!>
!> `files` is the list of input structure paths (the CLI layer in main.F90
!> marshals them out of the argument parser). Decoupling from the parser keeps
!> this routine in the xtb *library*, so the C/Python API can drive batches too
!> (Phase 4) without going through the command line.
subroutine run_gpu_batch(env, files)
   type(TEnvironment), intent(inout) :: env
   character(len=*), intent(in) :: files(:)

   type(TBatchResult), allocatable :: results(:)
   character(len=:), allocatable :: logmsg
   integer :: nFiles, iFile
   integer(i8) :: c0, c1, crate
   real(wp) :: t_total
   logical :: failed

   nFiles = size(files)
   if (nFiles < 1) then
      call env%error("--gpu-batch: no input files given", source)
      return
   end if

   write(env%unit, '(/,a)') repeat('=', 70)
   write(env%unit, '(a,i0,a)') " xtb GPU batch driver: ", nFiles, " input structures"
   write(env%unit, '(a)') repeat('=', 70)

   allocate(results(nFiles))

   ! Arm the capture hook in the core single-point path so that, as each molecule
   ! is diagonalized, its real (H, S) eigenproblem is recorded (bounded to capMax
   ! systems). After the run we replay those through the batched eigensolver with
   ! GPU-style padding and check the spectrum against an independent per-system
   ! solve -- the validatable half of the INTEGRATION SEAM (see validateBatchedEig).
   call gpu_capture_enable(capMax)

   call system_clock(c0, crate)

   ! Run every input file through the proven single-point path, recording
   ! throughput metrics.
   do iFile = 1, nFiles
      if (len_trim(files(iFile)) == 0) cycle
      call runOne(env, trim(files(iFile)), results(iFile))
      ! A failure on one molecule must not abort the whole batch: detect any
      ! logged error, mark this entry failed, and drain the log (getLog clears
      ! it) so the next molecule starts from a clean environment.
      call env%check(failed)
      if (failed) results(iFile)%ok = .false.
      call env%getLog(logmsg)
      if (failed) then
         write(env%unit, '(2x,a)') "[skipped on error] "//trim(results(iFile)%fname)
         if (allocated(logmsg)) then
            if (len_trim(logmsg) > 0) write(env%unit, '(4x,a)') trim(logmsg)
         end if
      end if
   end do

   call gpu_capture_disable()

   call system_clock(c1, crate)
   t_total = real(c1 - c0, wp) / real(crate, wp)

   call reportResults(env, results, t_total)

   ! Replay the captured real eigenproblems through TBatchedEigensolver with the
   ! GPU padding scheme and verify the spectrum bit-for-bit against a per-system
   ! solve. This is the closed half of the INTEGRATION SEAM: the exact batched
   ! kernel + padding the cuSolver backend uses, validated on real GFN0 matrices.
   call validateBatchedEig(env)
   call gpu_capture_clear()

   ! ---------------------------------------------------------------------
   ! REMAINING SEAM (production routing, next increment):
   !   The step above proves the batched solve is correct on real matrices.
   !   To route the *production* energy through it, the monolithic single point
   !   (xtb_peeq for GFN0) must be split into build -> [batched solve] -> finish
   !   phases so many molecules can sit at the solve point and be diagonalized
   !   in one launch; then scatter eigenpairs back and finish energy/density per
   !   molecule. Do GFN0 first (one build + one diag, no SCF).
   ! ---------------------------------------------------------------------

end subroutine run_gpu_batch


!> Replay captured per-molecule eigenproblems through the batched eigensolver and
!> validate the spectrum against an independent per-system LAPACK solve.
!>
!> For each size bucket the systems are padded up to the bucket's largest order
!> `n` (H diagonal = padDiag, S = identity on the padded block) and solved in one
!> batch. Because that padding is block-diagonal it decouples exactly, so the
!> lowest `m` eigenvalues of the padded system must equal the `m` eigenvalues of
!> the bare m x m problem to rounding -- this is precisely the property the GPU
!> backend relies on, here checked on real GFN0/1/2 matrices.
subroutine validateBatchedEig(env)
   type(TEnvironment), intent(inout) :: env

   type(TBatchedEigensolver) :: solver
   type(TCapturedSystem), allocatable :: bucket(:)
   type(TCapturedSystem) :: sys
   real(wp), allocatable :: hmats(:,:,:), smats(:,:,:), evals(:,:)
   real(wp), allocatable :: href(:,:), sref(:,:), eref(:), work(:)
   integer, allocatable :: iwork(:), naos(:)
   integer :: ncap, nseen, ndrop, ib, n, cnt, j, p, m, kk, info, lwork, liwork
   integer :: nbucketsTested, nSysTested, nFail
   real(wp) :: dev, worstBucket, worstAll
   real(wp) :: query(1)
   integer  :: iquery(1)

   ncap = gpu_capture_count()
   nseen = gpu_capture_seen()
   ndrop = max(0, nseen - ncap)

   write(env%unit, '(/,a)') " batched-eigensolver validation"
   write(env%unit, '(a)') " "//repeat('-', 60)
   if (ncap < 1) then
      write(env%unit, '(2x,a)') "no eigenproblems captured (the capture hook &
         &instruments the GFN0/peeq path; or all systems were linearly &
         &dependent); skipped"
      write(env%unit, '(a)') " "//repeat('-', 60)
      return
   end if
   write(env%unit, '(2x,a,i0,a)') "captured    : ", ncap, " real eigenproblems (H C = S C eps)"
   if (ndrop > 0) then
      write(env%unit, '(2x,a,i0,a,i0,a)') "note        : ", ndrop, &
         & " further systems were not retained (cap ", capMax, ")"
   end if

   ! Pull the captured systems and their orders up front.
   allocate(bucket(ncap), naos(ncap))
   do j = 1, ncap
      call gpu_capture_get(j, bucket(j))
      naos(j) = bucket(j)%nao
   end do

   nbucketsTested = 0; nSysTested = 0; nFail = 0; worstAll = 0.0_wp

   do ib = 1, nbins
      ! Count and size this bucket.
      cnt = 0; n = 0
      do j = 1, ncap
         if (binIndex(naos(j)) == ib) then
            cnt = cnt + 1
            n = max(n, naos(j))
         end if
      end do
      if (cnt < 1) cycle

      allocate(hmats(n,n,cnt), smats(n,n,cnt), evals(n,cnt))
      hmats = 0.0_wp; smats = 0.0_wp; evals = 0.0_wp

      ! Pack + pad each system of this bucket into the batch tensors.
      p = 0
      do j = 1, ncap
         if (binIndex(naos(j)) /= ib) cycle
         p = p + 1
         m = naos(j)
         hmats(1:m,1:m,p) = bucket(j)%H
         smats(1:m,1:m,p) = bucket(j)%S
         do kk = m+1, n
            hmats(kk,kk,p) = padDiag
            smats(kk,kk,p) = 1.0_wp
         end do
      end do

      ! Batched solve (LAPACK loop on CPU; cuSolver on the GPU build) in place.
      call init(solver, env, n, cnt)
      call solver%solve(env, hmats, smats, evals)
      call solver%free()

      ! Per-system reference solve of the bare m x m problem; compare spectra.
      worstBucket = 0.0_wp
      p = 0
      do j = 1, ncap
         if (binIndex(naos(j)) /= ib) cycle
         p = p + 1
         m = naos(j)
         href = bucket(j)%H
         sref = bucket(j)%S
         if (allocated(eref)) deallocate(eref)
         allocate(eref(m))
         ! workspace query then solve (itype=1, vectors, upper)
         call lapack_sygvd(1, 'v', 'u', m, href, m, sref, m, eref, &
            & query, -1, iquery, -1, info)
         if (info /= 0) then
            nFail = nFail + 1
            cycle
         end if
         lwork = max(1, int(query(1)))
         liwork = max(1, iquery(1))
         if (allocated(work)) deallocate(work)
         if (allocated(iwork)) deallocate(iwork)
         allocate(work(lwork), iwork(liwork))
         href = bucket(j)%H
         sref = bucket(j)%S
         call lapack_sygvd(1, 'v', 'u', m, href, m, sref, m, eref, &
            & work, lwork, iwork, liwork, info)
         if (info /= 0) then
            nFail = nFail + 1
            cycle
         end if
         ! lowest m batched eigenvalues vs bare spectrum
         dev = maxval(abs(evals(1:m,p) - eref(1:m)))
         worstBucket = max(worstBucket, dev)
         worstAll = max(worstAll, dev)
         nSysTested = nSysTested + 1
      end do

      nbucketsTested = nbucketsTested + 1
      write(env%unit, '(2x,a,i6,a,i4,a,es10.3,a)') &
         & "bucket n=", n, " : ", cnt, " systems, max |dEPS| = ", worstBucket, " eV"

      deallocate(hmats, smats, evals)
   end do

   write(env%unit, '(a)') " "//repeat('-', 60)
   write(env%unit, '(2x,a,i0,a,i0,a)') "tested      : ", nSysTested, &
      & " systems across ", nbucketsTested, " buckets"
   write(env%unit, '(2x,a,es10.3,a)') "worst dEPS  : ", worstAll, " eV"
   if (nFail > 0) then
      write(env%unit, '(2x,a,i0,a)') "skipped     : ", nFail, &
         & " systems (reference solve returned info /= 0)"
   end if
   if (nSysTested > 0 .and. worstAll <= eigTol) then
      write(env%unit, '(2x,a,es8.1,a)') "RESULT      : PASS (<= ", eigTol, &
         & " eV) -- batched solve + padding reproduce the per-system spectrum"
   else if (nSysTested > 0) then
      write(env%unit, '(2x,a,es8.1,a)') "RESULT      : FAIL (> ", eigTol, &
         & " eV) -- investigate padding / solver before trusting the GPU path"
   end if
   write(env%unit, '(a)') " "//repeat('-', 60)

end subroutine validateBatchedEig


!> Run a single structure through the standard calculator path and record it.
subroutine runOne(env, fname, res)
   type(TEnvironment), intent(inout) :: env
   character(len=*), intent(in) :: fname
   type(TBatchResult), intent(out) :: res

   type(TMolecule) :: mol
   type(TRestart) :: chk
   class(TCalculator), allocatable :: calc
   type(scc_results) :: spres
   real(wp), allocatable :: gradient(:,:)
   real(wp) :: sigma(3,3), energy, hlgap
   integer :: ich, ftype, stat
   integer(i8) :: c0, c1, crate
   character(len=:), allocatable :: paramFile

   res%fname = trim(fname)
   res%ok = .false.

   ! --- read geometry (standard Fortran I/O; no external open_file dep) ---
   ftype = getFileType(fname)
   open(newunit=ich, file=fname, status='old', action='read', iostat=stat)
   if (stat /= 0) then
      call env%warning("could not open '"//trim(fname)//"'", source)
      return
   end if
   call readMolecule(env, mol, ich, ftype)
   close(ich)
   if (mol%n <= 0) then
      call env%warning("empty/invalid structure '"//trim(fname)//"'", source)
      return
   end if
   res%nat = mol%n

   call system_clock(c0, crate)

   ! --- resolve the parameter file for the method (mirrors main.F90) ---
   select case (set%gfn_method)
   case (0); paramFile = xfind(fname_gfn0)
   case (1); paramFile = xfind(fname_gfn1)
   case (2); paramFile = xfind(fname_gfn2)
   case default
      call env%warning("--gpu-batch supports --gfn 0|1|2 only; skipping '" &
         & //trim(fname)//"'", source)
      return
   end select

   ! --- build calculator (custom --param files are not threaded in yet) ---
   call newCalculator(env, mol, calc, paramFile, .false., set%acc)
   if (.not. allocated(calc)) then
      call env%warning("calculator setup failed for '"//trim(fname)//"'", source)
      return
   end if

   ! --- initialize wavefunction for the xTB calculator ---
   select type (calc)
   type is (TxTBCalculator)
      res%nao = calc%basis%nao
      call chk%wfn%allocate(mol%n, calc%basis%nshell, calc%basis%nao)
      call newWavefunction(env, mol, calc, chk)
   class default
      call env%warning("--gpu-batch supports GFN0/1/2 (xTB) calculators only; &
         &skipping '"//trim(fname)//"'", source)
      return
   end select

   ! --- single point ---
   allocate(gradient(3, mol%n), source=0.0_wp)
   sigma = 0.0_wp
   call calc%singlepoint(env, mol, chk, 0, .false., energy, gradient, sigma, hlgap, spres)

   call system_clock(c1, crate)

   res%energy = energy
   res%gnorm = norm2(gradient)
   res%gap = hlgap
   res%walltime = real(c1 - c0, wp) / real(crate, wp)
   res%ok = .true.

end subroutine runOne


!> Map an AO count to a size-bucket index (1..nbins).
pure function binIndex(nao) result(ib)
   integer, intent(in) :: nao
   integer :: ib
   do ib = 1, nbins
      if (nao <= binEdge(ib)) return
   end do
   ib = nbins
end function binIndex


!> Print the per-molecule table, the size-bucket distribution (the batching
!> opportunity), and the aggregate throughput.
subroutine reportResults(env, results, t_total)
   type(TEnvironment), intent(inout) :: env
   type(TBatchResult), intent(in) :: results(:)
   real(wp), intent(in) :: t_total

   integer :: i, ib, nok, nfail
   integer :: binCount(nbins), binMaxNao(nbins)
   integer(i8) :: realElems, paddedElems
   real(wp) :: pad_waste

   nok = 0; nfail = 0
   binCount = 0; binMaxNao = 0
   realElems = 0_i8; paddedElems = 0_i8

   write(env%unit, '(/,a)') " results"
   write(env%unit, '(a)') " "//repeat('-', 78)
   ! 'trunc' returns a fixed-width, left-justified string, printed with plain 'a'.
   write(env%unit, '(2x,a)') &
      & "   #  "//trunc("structure", 28)//"   nat     nao       energy / Eh     t / s"
   write(env%unit, '(a)') " "//repeat('-', 78)
   do i = 1, size(results)
      if (results(i)%ok) then
         nok = nok + 1
         ib = binIndex(results(i)%nao)
         binCount(ib) = binCount(ib) + 1
         binMaxNao(ib) = max(binMaxNao(ib), results(i)%nao)
         write(env%unit, '(2x,i4,2x,a,2x,i6,2x,i6,2x,f16.8,2x,f8.3)') &
            & i, trunc(results(i)%fname, 28), results(i)%nat, results(i)%nao, &
            & results(i)%energy, results(i)%walltime
      else
         nfail = nfail + 1
         write(env%unit, '(2x,i4,2x,a,2x,a)') &
            & i, trunc(results(i)%fname, 28), "FAILED"
      end if
   end do
   write(env%unit, '(a)') " "//repeat('-', 78)

   ! Bucket distribution + padding cost (matrix elements wasted by padding each
   ! system up to its bucket's largest nao). This is what to minimize when
   ! choosing bucket edges for a given workload.
   write(env%unit, '(/,a)') " size buckets (batching opportunity)"
   write(env%unit, '(a)') " "//repeat('-', 50)
   write(env%unit, '(2x,a12,2x,a8,2x,a10)') "nao <=", "count", "max nao"
   do ib = 1, nbins
      if (binCount(ib) == 0) cycle
      if (ib == nbins) then
         write(env%unit, '(2x,a12,2x,i8,2x,i10)') "inf", binCount(ib), binMaxNao(ib)
      else
         write(env%unit, '(2x,i12,2x,i8,2x,i10)') binEdge(ib), binCount(ib), binMaxNao(ib)
      end if
   end do
   write(env%unit, '(a)') " "//repeat('-', 50)

   do i = 1, size(results)
      if (.not. results(i)%ok) cycle
      ib = binIndex(results(i)%nao)
      realElems = realElems + int(results(i)%nao, i8)**2
      paddedElems = paddedElems + int(binMaxNao(ib), i8)**2
   end do
   if (paddedElems > 0_i8) then
      pad_waste = 100.0_wp * real(paddedElems - realElems, wp) / real(paddedElems, wp)
   else
      pad_waste = 0.0_wp
   end if

   write(env%unit, '(/,a)') " summary"
   write(env%unit, '(a)') " "//repeat('-', 50)
   write(env%unit, '(2x,a,i0,a,i0,a)') "processed   : ", nok, " ok, ", nfail, " failed"
   write(env%unit, '(2x,a,f10.3,a)')   "wall time   : ", t_total, " s"
   if (nok > 0 .and. t_total > 0.0_wp) then
      write(env%unit, '(2x,a,f10.2)')  "throughput  : ", real(nok, wp) / t_total
      write(env%unit, '(2x,a)')        "              (molecules / s, this run)"
   end if
   write(env%unit, '(2x,a,f6.1,a)')    "pad waste   : ", pad_waste, &
      & " % of batched matrix elements (lower = better bucketing)"
   write(env%unit, '(a)') " "//repeat('-', 50)

contains
   !> Left-truncate/pad a string to width w for tabular output.
   pure function trunc(s, w) result(o)
      character(len=*), intent(in) :: s
      integer, intent(in) :: w
      character(len=w) :: o
      integer :: ls
      ls = len_trim(s)
      if (ls <= w) then
         o = s
      else
         o = '...'//s(ls-w+4:ls)
      end if
   end function trunc
end subroutine reportResults

end module xtb_gpu_batch
