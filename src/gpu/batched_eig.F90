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

!> Batched solver for generalized symmetric-definite eigenproblems
!>
!>     H_k C_k = S_k C_k diag(eps_k),   k = 1 .. nbatch
!>
!> where every system in the batch shares the same leading dimension `n`
!> (the caller pads each molecule's Hamiltonian/overlap up to the bucket
!> size; see xtb_gpu_batch). Solving many small systems together is what
!> keeps the GPU saturated for high-throughput screening of small molecules.
!>
!> Two backends are provided behind the USE_CUSOLVER macro:
!>
!>   * CPU reference (default): loops LAPACK ?sygvd over the batch. This is
!>     bit-for-bit the same routine the rest of xtb uses, so it is the
!>     numerical ground truth the GPU path is validated against, and it lets
!>     the batch driver run end-to-end without an NVIDIA toolchain.
!>
!>   * GPU (USE_CUSOLVER): offloads each system's diagonalization to
!>     cuSolverDn (cusolverDnDsygvd), modeled on the existing single-system
!>     path in src/mctc/lapack/eigensolve.F90. Data is staged once for the
!>     whole batch to amortize PCIe transfers.
!>
!> NOTE (performance roadmap): the highest-throughput backend replaces the
!> per-system cusolverDnDsygvd loop with a fully batched pipeline --
!> cusolverDnDpotrfBatched (Cholesky of S) -> cublasDtrsmBatched (reduce to
!> standard form) -> cusolverDnDsyevjBatched (batched Jacobi) -> trsm
!> (back-transform). That is the documented next optimization; the loop below
!> is the correct-by-construction first step that reuses the proven kernel.
module xtb_gpu_batched_eig
   use xtb_mctc_accuracy, only : wp, dp
   use xtb_type_environment, only : TEnvironment
   use xtb_mctc_lapack_geneigval, only : lapack_sygvd
#ifdef USE_CUSOLVER
   use xtb_mctc_global, only : cusolverDnH
   use cusolverDn
#endif
   implicit none
   private

   public :: TBatchedEigensolver, init

   !> Reusable workspace + metadata for a batch of identically-sized eigenproblems
   type :: TBatchedEigensolver
      !> Padded leading dimension, uniform across the batch
      integer :: n = 0
      !> Number of systems in the batch
      integer :: nbatch = 0
      !> Real workspace (LAPACK or cuSolver, depending on backend)
      real(dp), allocatable :: dwork(:)
#ifdef USE_CUSOLVER
      !> cuSolver workspace length returned by the bufferSize query
      integer :: lwork = 0
#else
      !> Integer workspace for LAPACK divide-and-conquer
      integer, allocatable :: iwork(:)
#endif
   contains
      procedure :: solve => solveBatch
      procedure :: free  => freeBatch
   end type TBatchedEigensolver

   interface init
      module procedure :: initBatchedEigensolver
   end interface init

   character(len=*), parameter :: source = 'gpu_batched_eig'

contains

!> Allocate workspace for a batch of `nbatch` eigenproblems of order `n`.
subroutine initBatchedEigensolver(self, env, n, nbatch)
   class(TBatchedEigensolver), intent(out) :: self
   type(TEnvironment), intent(inout) :: env
   integer, intent(in) :: n
   integer, intent(in) :: nbatch
#ifdef USE_CUSOLVER
   integer :: istat, lwork
   ! Placeholder; cusolver only needs valid sizes for the bufferSize query.
   real(dp) :: dummy(1)
#endif

   self%n = n
   self%nbatch = nbatch
   if (n <= 0 .or. nbatch <= 0) return

#ifdef USE_CUSOLVER
   istat = cusolverDnDsygvd_bufferSize(cusolverDnH, CUSOLVER_EIG_TYPE_1, &
      & CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, dummy, n, &
      & dummy, n, dummy, lwork)
   if (istat /= 0) then
      call env%error("cusolverDnDsygvd_bufferSize failed", source)
      return
   end if
   self%lwork = lwork
   allocate(self%dwork(lwork))
#else
   ! LAPACK divide-and-conquer workspace for ?sygvd with jobz='v'
   allocate(self%dwork(1 + 6*n + 2*n**2))
   allocate(self%iwork(3 + 5*n))
#endif

end subroutine initBatchedEigensolver


!> Solve the whole batch in place.
!>
!>   hmats(n,n,nbatch)  in : symmetric Hamiltonian (upper triangle used)
!>                     out : eigenvectors C_k (column k-block)
!>   smats(n,n,nbatch)  in : symmetric-positive-definite overlap
!>                     out : overwritten (Cholesky factor / cuSolver scratch)
!>   evals(n,nbatch)   out : eigenvalues, ascending, per system
!>
!> For padded systems the caller must set the padded diagonal of H to a large
!> value and S to the identity there, so the spurious eigenvalues sort to the
!> top of the spectrum and never enter the occupied window.
subroutine solveBatch(self, env, hmats, smats, evals)
   class(TBatchedEigensolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(wp), intent(inout), contiguous :: hmats(:,:,:)
   real(wp), intent(inout), contiguous :: smats(:,:,:)
   real(wp), intent(out),   contiguous :: evals(:,:)

   integer :: k, n, nbatch, info
#ifdef USE_CUSOLVER
   integer :: istat
#else
   integer :: ldwork, liwork
#endif

   n = self%n
   nbatch = self%nbatch
   if (n <= 0 .or. nbatch <= 0) return

#ifdef USE_CUSOLVER
   ! Stage the entire batch on the device once, then issue one diagonalization
   ! per system reusing the shared workspace. (See module note for the fully
   ! batched cuSolver pipeline that supersedes this loop.)
   !$acc enter data copyin(hmats, smats) create(evals, self%dwork, info)
   do k = 1, nbatch
      !$acc host_data use_device(hmats, smats, evals, self%dwork, info)
      istat = cusolverDnDsygvd(cusolverDnH, CUSOLVER_EIG_TYPE_1, &
         & CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, &
         & hmats(1,1,k), n, smats(1,1,k), n, evals(1,k), &
         & self%dwork, self%lwork, info)
      !$acc end host_data
      if (istat /= 0) then
         call env%error("cusolverDnDsygvd failed in batch", source)
         exit
      end if
   end do
   !$acc wait
   !$acc exit data copyout(hmats, evals) delete(smats, self%dwork, info)
#else
   ldwork = size(self%dwork)
   liwork = size(self%iwork)
   do k = 1, nbatch
      call lapack_sygvd(1, 'v', 'u', n, hmats(:,:,k), n, smats(:,:,k), n, &
         & evals(:,k), self%dwork, ldwork, self%iwork, liwork, info)
      if (info /= 0) then
         call env%error("LAPACK sygvd failed in batch", source)
         return
      end if
   end do
#endif

end subroutine solveBatch


subroutine freeBatch(self)
   class(TBatchedEigensolver), intent(inout) :: self
   if (allocated(self%dwork)) deallocate(self%dwork)
#ifndef USE_CUSOLVER
   if (allocated(self%iwork)) deallocate(self%iwork)
#endif
   self%n = 0
   self%nbatch = 0
end subroutine freeBatch

end module xtb_gpu_batched_eig
