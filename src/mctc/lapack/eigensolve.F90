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
   use xtb_mctc_blas_level3, only : blas_trsm
   use xtb_mctc_lapack_geneigval, only : lapack_sygvd
   use xtb_mctc_lapack_stdeigval, only : lapack_syevd
   use xtb_mctc_lapack_gst, only : lapack_sygst
   use xtb_mctc_lapack_trf, only : mctc_potrf
   use xtb_type_environment, only : TEnvironment
#ifdef USE_CUSOLVER
   use xtb_mctc_global
   use cusolverDn
#endif
   implicit none
   private

   public :: TEigenSolver, init


   type :: TEigenSolver
      private
      integer :: n
      integer, allocatable :: iwork(:)
      real(sp), allocatable :: swork(:)
      real(sp), allocatable :: sbmat(:, :)
      real(dp), allocatable :: dwork(:)
      real(dp), allocatable :: dbmat(:, :)
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

   self%n = size(bmat, 1)

   allocate(self%swork(1 + 6*self%n + 2*self%n**2))
   allocate(self%iwork(3 + 5*self%n))

   self%sbmat = bmat
   ! Check for Cholesky factorisation
   call mctc_potrf(env, self%sbmat)

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
   real(dp) :: dummy(:) 
#endif

   self%n = size(bmat, 1)

#ifdef USE_CUSOLVER
   istat = cusolverDnDsygvd_bufferSize(cusolverDnH, CUSOLVER_EIG_TYPE_1, &
     CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, self%n, dummy,    &
     self%n, dummy, self%n, dummy, lwork)
   if (istat /= 0) then
      call env%error("failed to get dygvd buffer size", source)
   end if

   self%lwork = lwork
   allocate(self%dwork(lwork))
#else
   allocate(self%dwork(1 + 6*self%n + 2*self%n**2))
   allocate(self%iwork(3 + 5*self%n))
#endif

   self%dbmat = bmat
   ! Check for Cholesky factorisation
   call mctc_potrf(env, self%dbmat)

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
   use xtb_tracying
   character(len=*), parameter :: source = 'mctc_lapack_dsygvd_factorized'
   class(TEigenSolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(in) :: bmat_factorized(:, :)
   real(dp), intent(out) :: eval(:)
   integer :: info, ldwork, liwork

   type(xtb_zone) :: ctx

   ldwork = size(self%dwork)
   liwork = size(self%iwork)

   call ctx%start("src/mctc/lapack/eigensolve.F90", source, __LINE__, zone_name="sygst", color=TracyColors%Pink)

   CALL lapack_sygst( 1, 'u', self%n, amat, self%n, bmat_factorized, self%n, info )
   if (info /= 0) then
      call env%error("Failed to reduce eigenvalue problem", source)
      return
   end if

   call ctx%end()
   call ctx%start("src/mctc/lapack/eigensolve.F90", source, __LINE__, zone_name="syevd", color=TracyColors%Purple)

   CALL lapack_syevd( 'v', 'u', self%n, amat, self%n, eval, self%dwork, ldwork, self%iwork, liwork, info )
   if (info /= 0) then
      call env%error("Failed to compute eigenvalues and eigenvectors", source)
      return
   end if

   call ctx%end()
   call ctx%start("src/mctc/lapack/eigensolve.F90", source, __LINE__, zone_name="trsm", color=TracyColors%Purple)

   CALL blas_trsm( 'l', 'u', 'n', 'n', self%n, self%n, 1.0_dp, bmat_factorized, self%n, amat, self%n )

   call ctx%end()

end subroutine mctc_dsygvd_factorized

end module xtb_mctc_lapack_eigensolve
