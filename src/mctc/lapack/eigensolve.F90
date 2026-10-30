! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
! Copyright (C) 2020, NVIDIA CORPORATION. All rights reserved.
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

!> Wrapper for eigensolver routines
module xtb_mctc_lapack_eigensolve
   use xtb_mctc_accuracy, only : sp, dp
   use xtb_mctc_blas_level3, only : blas_trsm, blas_trmm, blas_gemm
   use xtb_mctc_lapack_geneigval, only : lapack_sygvd
   use xtb_mctc_lapack_stdeigval, only : lapack_syevd, lapack_sytrd, &
      & lapack_stedc, lapack_ormtr
   use xtb_mctc_lapack_gst, only : lapack_sygst
   use xtb_mctc_lapack_trf, only : mctc_potrf
   use xtb_mctc_lapack_tri, only : mctc_trtri
   use xtb_type_environment, only : TEnvironment
#ifdef USE_CUSOLVER
   use xtb_mctc_global
   use cusolverDn
#endif
   implicit none
   private

   public :: TEigenSolver, init


   !> Column-panel width for upper-triangle congruence A = R^T H R (R = U^-1).
   !> Because standard BLAS has no symmetric R^T H R routine, two full-matrix
   !> blas_trmm calls cost 2*N^3 FLOPs; evaluating only rows 1:jj per column
   !> panel j:jj exploits R^T being lower triangular to compute just the upper
   !> triangle of A in 1*N^3 FLOPs (matching dsygst while using gemm + trmm).
   integer, parameter :: nb = 96


   type :: TEigenSolver
      private
      integer :: n
      integer, allocatable :: iwork(:)
      real(sp), allocatable :: swork(:)
      real(sp), allocatable :: sbmat(:, :)
      real(sp), allocatable :: szmat(:, :)
      real(sp), allocatable :: stau(:)
      real(sp), allocatable :: sesub(:)
      real(sp), allocatable :: swork_trd(:)
      real(dp), allocatable :: dwork(:)
      real(dp), allocatable :: dbmat(:, :)
      real(dp), allocatable :: dzmat(:, :)
      real(dp), allocatable :: dtau(:)
      real(dp), allocatable :: desub(:)
      real(dp), allocatable :: dwork_trd(:)
#ifdef USE_CUSOLVER
      integer :: lwork
#endif
   contains
      generic :: solve => sgen_solve, dgen_solve
      procedure :: sgen_solve => mctc_ssygvd
      procedure :: dgen_solve => mctc_dsygvd
      generic :: fact_solve => sfact_solve, dfact_solve
      procedure :: sfact_solve => mctc_ssygvd_factorized
      procedure :: dfact_solve => mctc_dsygvd_factorized
      generic :: fact_solve_subspace => sfact_solve_subspace, dfact_solve_subspace
      procedure :: sfact_solve_subspace => mctc_ssygvd_fact_subspace
      procedure :: dfact_solve_subspace => mctc_dsygvd_fact_subspace
      generic :: fact_backtransform => sfact_backtransform, dfact_backtransform
      procedure :: sfact_backtransform => mctc_ssygvd_fact_backtransform
      procedure :: dfact_backtransform => mctc_dsygvd_fact_backtransform
   end type TEigenSolver


   interface init
      module procedure :: initSEigenSolver
      module procedure :: initDEigenSolver
   end interface init


contains


subroutine initSEigenSolver(self, env, bmat)
   character(len=*), parameter :: source = 'mctc_lapack_sygvd'
   class(TEigenSolver), intent(out) :: self
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: bmat(:, :)
   integer :: info, lwork_trd, lda
   real(sp) :: work_trd_q(1), work_orm_q(1)
   logical :: exitRun

   self%n = size(bmat, 1)

   allocate(self%swork(1 + 6*self%n + 2*self%n**2))
   allocate(self%iwork(3 + 5*self%n))
   allocate(self%szmat(self%n, self%n), source=0.0_sp)
   allocate(self%stau(self%n), self%sesub(self%n))

   self%sbmat = bmat
   lda = max(1, self%n)
   call lapack_sytrd('U', self%n, self%szmat, lda, self%swork, self%sesub, &
      & self%stau, work_trd_q, -1, info)
   call lapack_ormtr('L', 'U', 'N', self%n, self%n, self%sbmat, lda, self%stau, &
      & self%szmat, lda, work_orm_q, -1, info)
   lwork_trd = max(1, self%n, int(work_trd_q(1)), int(work_orm_q(1)))
   allocate(self%swork_trd(lwork_trd))

   ! Check for Cholesky factorisation
   call mctc_potrf(env, self%sbmat)
   call env%check(exitRun)
   if (exitRun) return
   call mctc_trtri(env, self%sbmat)

end subroutine initSEigenSolver


subroutine initDEigenSolver(self, env, bmat)
   character(len=*), parameter :: source = 'mctc_lapack_sygvd'
   class(TEigenSolver), intent(out) :: self
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: bmat(:, :)
#ifdef USE_CUSOLVER
   integer :: istat, lwork
   ! dummy is only a dummy argument used to query the workspace size needed
   ! for cuSolverDnDsygvd -- it is okay to pass an empty array to cuSolverDnDsygvd_bufferSize
   real(dp) :: dummy(1)
#endif
   integer :: info, ldwork_cpu, lwork_trd, lda
   real(dp) :: work_trd_q(1), work_orm_q(1)
   logical :: exitRun

   self%n = size(bmat, 1)
   ldwork_cpu = 1 + 6*self%n + 2*self%n**2

#ifdef USE_CUSOLVER
   istat = cusolverDnDsygvd_bufferSize(cusolverDnH, CUSOLVER_EIG_TYPE_1, &
     CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, self%n, dummy,    &
     self%n, dummy, self%n, dummy, lwork)
   if (istat /= 0) then
      call env%error("failed to get dygvd buffer size", source)
   end if

   self%lwork = lwork
   allocate(self%dwork(max(lwork, ldwork_cpu)))
#else
   allocate(self%dwork(ldwork_cpu))
#endif
   allocate(self%iwork(3 + 5*self%n))
   allocate(self%dzmat(self%n, self%n), source=0.0_dp)
   allocate(self%dtau(self%n), self%desub(self%n))

   self%dbmat = bmat
   lda = max(1, self%n)
   call lapack_sytrd('U', self%n, self%dzmat, lda, self%dwork, self%desub, &
      & self%dtau, work_trd_q, -1, info)
   call lapack_ormtr('L', 'U', 'N', self%n, self%n, self%dbmat, lda, self%dtau, &
      & self%dzmat, lda, work_orm_q, -1, info)
   lwork_trd = max(1, self%n, int(work_trd_q(1)), int(work_orm_q(1)))
   allocate(self%dwork_trd(lwork_trd))

   ! Check for Cholesky factorisation
   call mctc_potrf(env, self%dbmat)
   call env%check(exitRun)
   if (exitRun) return
   call mctc_trtri(env, self%dbmat)

end subroutine initDEigenSolver


subroutine mctc_ssygvd(self, env, amat, bmat, eval)
   character(len=*), parameter :: source = 'mctc_lapack_sygvd'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(out) :: eval(:)
   integer :: info, lswork, liwork

   self%sbmat(:, :) = bmat

   lswork = size(self%swork)
   liwork = size(self%iwork)
   call lapack_sygvd(1, 'v', 'u', self%n, amat, self%n, self%sbmat, self%n, eval, &
      & self%swork, lswork, self%iwork, liwork, info)

   if (info /= 0) then
      call env%error("Failed to solve eigenvalue problem", source)
   end if

end subroutine mctc_ssygvd


subroutine mctc_dsygvd(self, env, amat, bmat, eval)
   character(len=*), parameter :: source = 'mctc_lapack_sygvd'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(out) :: eval(:)
   integer :: info, ldwork, liwork
#ifdef USE_CUSOLVER
   integer :: istat
#endif

   self%dbmat(:, :) = bmat

#ifdef USE_CUSOLVER
   !$acc enter data copyin(amat, self%dbmat, eval, info) create(self%dwork)

   !$acc host_data use_device(amat, self%dbmat, eval, self%dwork, info)
   istat = cusolverDnDsygvd(cusolverDnH, CUSOLVER_EIG_TYPE_1, &
     CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, self%n, amat, self%n, &
     self%dbmat, self%n, eval, self%dwork, self%lwork, info)
   !$acc end host_data

   !$acc exit data copyout(amat, self%dbmat, eval, info) delete(self%dwork)

   if (istat /= 0) then
      call env%error("cuSovlerDnDsygvd failed", source)
   end if
#else
   ldwork = size(self%dwork)
   liwork = size(self%iwork)
   call lapack_sygvd(1, 'v', 'u', self%n, amat, self%n, self%dbmat, self%n, eval, &
      & self%dwork, ldwork, self%iwork, liwork, info)
#endif

   if (info /= 0) then
      call env%error("Failed to solve eigenvalue problem", source)
   end if

end subroutine mctc_dsygvd


subroutine mctc_ssygvd_factorized(self, env, amat, bmat_factorized, eval)
   character(len=*), parameter :: source = 'mctc_lapack_ssygvd_factorized'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(in) :: bmat_factorized(:, :)
   real(sp), intent(out) :: eval(:)
   integer :: info, lswork, liwork

   lswork = size(self%swork)
   liwork = size(self%iwork)

   CALL lapack_sygst( 1, 'u', self%n, amat, self%n, bmat_factorized, self%n, info )

   if (info /= 0) then
      call env%error("Failed to reduce eigenvalue problem", source)
      return
   end if

   CALL lapack_syevd( 'v', 'u', self%n, amat, self%n, eval, self%swork, lswork, self%iwork, liwork, info )

   if (info /= 0) then
      call env%error("Failed to compute eigenvalues and eigenvectors", source)
      return
   end if

   CALL blas_trsm( 'l', 'u', 'n', 'n', self%n, self%n, 1.0_sp, bmat_factorized, self%n, amat, self%n )

end subroutine mctc_ssygvd_factorized


subroutine mctc_dsygvd_factorized(self, env, amat, bmat_factorized, eval)
   character(len=*), parameter :: source = 'mctc_lapack_dsygvd_factorized'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(in) :: bmat_factorized(:, :)
   real(dp), intent(out) :: eval(:)
   integer :: info, ldwork, liwork

   ldwork = size(self%dwork)
   liwork = size(self%iwork)

   CALL lapack_sygst( 1, 'u', self%n, amat, self%n, bmat_factorized, self%n, info )

   if (info /= 0) then
      call env%error("Failed to reduce eigenvalue problem", source)
      return
   end if

   CALL lapack_syevd( 'v', 'u', self%n, amat, self%n, eval, self%dwork, ldwork, self%iwork, liwork, info )

   if (info /= 0) then
      call env%error("Failed to compute eigenvalues and eigenvectors", source)
      return
   end if

   CALL blas_trsm( 'l', 'u', 'n', 'n', self%n, self%n, 1.0_dp, bmat_factorized, self%n, amat, self%n )

end subroutine mctc_dsygvd_factorized


subroutine mctc_ssygvd_fact_subspace(self, env, amat, eval, cmat, m_sub)
   character(len=*), parameter :: source = 'mctc_lapack_ssygvd_fact_subspace'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(out) :: eval(:)
   real(sp), intent(inout), optional :: cmat(:, :)
   integer, intent(in), optional :: m_sub
   integer :: j, jj, k, info, lswork, liwork, lswork_trd

   lswork_trd = size(self%swork_trd)
   lswork = size(self%swork)
   liwork = size(self%iwork)

   do j = 1, self%n
      do jj = 1, j - 1
         amat(j, jj) = amat(jj, j)
      end do
   end do
   do j = 1, self%n, nb
      jj = min(self%n, j + nb - 1)
      k = jj - j + 1
      self%szmat(1:jj, j:jj) = amat(1:jj, j:jj)
      call blas_trmm('R', 'U', 'N', 'N', jj, k, 1.0_sp, self%sbmat(j:jj, j:jj), k, &
         & self%szmat(:, j:jj), self%n)
      if (j > 1) then
         call blas_gemm('N', 'N', jj, k, j - 1, 1.0_sp, amat, self%n, &
            & self%sbmat(:, j:jj), self%n, 1.0_sp, self%szmat(:, j:jj), self%n)
      end if
      call blas_trmm('L', 'U', 'T', 'N', jj, k, 1.0_sp, self%sbmat, self%n, &
         & self%szmat(:, j:jj), self%n)
   end do
   do j = 1, self%n
      amat(1:j, j) = self%szmat(1:j, j)
   end do

   call lapack_sytrd('U', self%n, amat, self%n, eval, self%sesub, self%stau, &
      & self%swork_trd, lswork_trd, info)
   if (info == 0) then
      call lapack_stedc('I', self%n, eval, self%sesub, self%szmat, self%n, &
         & self%swork, lswork, self%iwork, liwork, info)
   end if
   if (info /= 0) then
      call env%error("Failed to compute eigenvalues and eigenvectors", source)
      return
   end if

   if (present(m_sub) .and. present(cmat)) then
      if (m_sub > 0) then
         call self%fact_backtransform(env, amat, cmat, 1, m_sub)
      end if
   end if

end subroutine mctc_ssygvd_fact_subspace


subroutine mctc_dsygvd_fact_subspace(self, env, amat, eval, cmat, m_sub)
   character(len=*), parameter :: source = 'mctc_lapack_dsygvd_fact_subspace'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(out) :: eval(:)
   real(dp), intent(inout), optional :: cmat(:, :)
   integer, intent(in), optional :: m_sub
   integer :: j, jj, k, info, ldwork, liwork, ldwork_trd

   ldwork_trd = size(self%dwork_trd)
   ldwork = size(self%dwork)
   liwork = size(self%iwork)

   do j = 1, self%n
      do jj = 1, j - 1
         amat(j, jj) = amat(jj, j)
      end do
   end do
   do j = 1, self%n, nb
      jj = min(self%n, j + nb - 1)
      k = jj - j + 1
      self%dzmat(1:jj, j:jj) = amat(1:jj, j:jj)
      call blas_trmm('R', 'U', 'N', 'N', jj, k, 1.0_dp, self%dbmat(j:jj, j:jj), k, &
         & self%dzmat(:, j:jj), self%n)
      if (j > 1) then
         call blas_gemm('N', 'N', jj, k, j - 1, 1.0_dp, amat, self%n, &
            & self%dbmat(:, j:jj), self%n, 1.0_dp, self%dzmat(:, j:jj), self%n)
      end if
      call blas_trmm('L', 'U', 'T', 'N', jj, k, 1.0_dp, self%dbmat, self%n, &
         & self%dzmat(:, j:jj), self%n)
   end do
   do j = 1, self%n
      amat(1:j, j) = self%dzmat(1:j, j)
   end do

   call lapack_sytrd('U', self%n, amat, self%n, eval, self%desub, self%dtau, &
      & self%dwork_trd, ldwork_trd, info)
   if (info == 0) then
      call lapack_stedc('I', self%n, eval, self%desub, self%dzmat, self%n, &
         & self%dwork, ldwork, self%iwork, liwork, info)
   end if
   if (info /= 0) then
      call env%error("Failed to compute eigenvalues and eigenvectors", source)
      return
   end if

   if (present(m_sub) .and. present(cmat)) then
      if (m_sub > 0) then
         call self%fact_backtransform(env, amat, cmat, 1, m_sub)
      end if
   end if

end subroutine mctc_dsygvd_fact_subspace


subroutine mctc_ssygvd_fact_backtransform(self, env, amat, cmat, ilo, ihi)
   character(len=*), parameter :: source = 'mctc_lapack_ssygvd_fact_backtransform'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: cmat(:, :)
   integer, intent(in) :: ilo, ihi
   integer :: m, info, lswork_trd

   m = ihi - ilo + 1
   if (m <= 0) return
   lswork_trd = size(self%swork_trd)

   call lapack_ormtr('L', 'U', 'N', self%n, m, amat, self%n, self%stau, &
      & self%szmat(:, ilo:ihi), self%n, self%swork_trd, lswork_trd, info)
   if (info /= 0) then
      call env%error("Failed to back-transform eigenvectors", source)
      return
   end if
   call blas_trmm('L', 'U', 'N', 'N', self%n, m, 1.0_sp, self%sbmat, self%n, &
      & self%szmat(:, ilo:ihi), self%n)
   cmat(:, ilo:ihi) = self%szmat(:, ilo:ihi)

end subroutine mctc_ssygvd_fact_backtransform


subroutine mctc_dsygvd_fact_backtransform(self, env, amat, cmat, ilo, ihi)
   character(len=*), parameter :: source = 'mctc_lapack_dsygvd_fact_backtransform'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: cmat(:, :)
   integer, intent(in) :: ilo, ihi
   integer :: m, info, ldwork_trd

   m = ihi - ilo + 1
   if (m <= 0) return
   ldwork_trd = size(self%dwork_trd)

   call lapack_ormtr('L', 'U', 'N', self%n, m, amat, self%n, self%dtau, &
      & self%dzmat(:, ilo:ihi), self%n, self%dwork_trd, ldwork_trd, info)
   if (info /= 0) then
      call env%error("Failed to back-transform eigenvectors", source)
      return
   end if
   call blas_trmm('L', 'U', 'N', 'N', self%n, m, 1.0_dp, self%dbmat, self%n, &
      & self%dzmat(:, ilo:ihi), self%n)
   cmat(:, ilo:ihi) = self%dzmat(:, ilo:ihi)

end subroutine mctc_dsygvd_fact_backtransform

end module xtb_mctc_lapack_eigensolve
