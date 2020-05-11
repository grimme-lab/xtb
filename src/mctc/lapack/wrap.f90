! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> LAPACK linear equation solvers.
module xtb_mctc_lapack_wrap
   use xtb_mctc_accuracy, only : sp, dp
   use xtb_mctc_lapack_gst, only : mctc_sygst, mctc_spgst
   use xtb_mctc_lapack_trf, only : mctc_getrf, mctc_sytrf, mctc_sptrf, &
      & mctc_potrf, mctc_pptrf
   use xtb_mctc_lapack_tri, only : mctc_getri, mctc_sytri, mctc_sptri, &
      & mctc_potri, mctc_pptri
   use xtb_mctc_lapack_trs, only : mctc_getrs, mctc_sytrs, mctc_sptrs, &
      & mctc_potrs, mctc_pptrs
   use xtb_type_environment, only : TEnvironment
   implicit none
   private

   public :: mctc_sygst, mctc_spgst
   public :: mctc_getrs, mctc_sytrs, mctc_sptrs, mctc_potrs, mctc_pptrs
   public :: mctc_getrf, mctc_sytrf, mctc_sptrf, mctc_potrf, mctc_pptrf
   public :: mctc_getri, mctc_sytri, mctc_sptri, mctc_potri, mctc_pptri


   interface mctc_getrs
      module procedure :: mctc_sgetrs1
      module procedure :: mctc_sgetrs3
      module procedure :: mctc_dgetrs1
      module procedure :: mctc_dgetrs3
   end interface mctc_getrs

   interface mctc_sytrs
      module procedure :: mctc_ssytrs1
      module procedure :: mctc_ssytrs3
      module procedure :: mctc_dsytrs1
      module procedure :: mctc_dsytrs3
   end interface mctc_sytrs

   interface mctc_sptrs
      module procedure :: mctc_ssptrs1
      module procedure :: mctc_ssptrs3
      module procedure :: mctc_dsptrs1
      module procedure :: mctc_dsptrs3
   end interface mctc_sptrs

   interface mctc_potrs
      module procedure :: mctc_spotrs1
      module procedure :: mctc_spotrs3
      module procedure :: mctc_dpotrs1
      module procedure :: mctc_dpotrs3
   end interface mctc_potrs

   interface mctc_pptrs
      module procedure :: mctc_spptrs1
      module procedure :: mctc_spptrs3
      module procedure :: mctc_dpptrs1
      module procedure :: mctc_dpptrs3
   end interface mctc_pptrs


contains


subroutine mctc_sgetrs1(env, amat, bvec, ipiv, trans)
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout), target :: bvec(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call mctc_getrs(env, amat, bptr, ipiv, trans)
end subroutine mctc_sgetrs1


subroutine mctc_sgetrs3(env, amat, bmat, ipiv, trans)
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout), contiguous, target :: bmat(:, :, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call mctc_getrs(env, amat, bptr, ipiv, trans)
end subroutine mctc_sgetrs3


subroutine mctc_dgetrs1(env, amat, bvec, ipiv, trans)
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout), target :: bvec(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call mctc_getrs(env, amat, bptr, ipiv, trans)
end subroutine mctc_dgetrs1


subroutine mctc_dgetrs3(env, amat, bmat, ipiv, trans)
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout), contiguous, target :: bmat(:, :, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call mctc_getrs(env, amat, bptr, ipiv, trans)
end subroutine mctc_dgetrs3


subroutine mctc_ssytrs1(env, amat, bvec, ipiv, uplo)
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout), target :: bvec(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call mctc_sytrs(env, amat, bptr, ipiv, uplo)
end subroutine mctc_ssytrs1


subroutine mctc_ssytrs3(env, amat, bmat, ipiv, uplo)
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout), contiguous, target :: bmat(:, :, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call mctc_sytrs(env, amat, bptr, ipiv, uplo)
end subroutine mctc_ssytrs3


subroutine mctc_dsytrs1(env, amat, bvec, ipiv, uplo)
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout), target :: bvec(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call mctc_sytrs(env, amat, bptr, ipiv, uplo)
end subroutine mctc_dsytrs1


subroutine mctc_dsytrs3(env, amat, bmat, ipiv, uplo)
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout), contiguous, target :: bmat(:, :, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call mctc_sytrs(env, amat, bptr, ipiv, uplo)
end subroutine mctc_dsytrs3


subroutine mctc_ssptrs1(env, amat, bvec, ipiv, uplo)
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:)
   real(sp), intent(inout), target :: bvec(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call mctc_sptrs(env, amat, bptr, ipiv, uplo)
end subroutine mctc_ssptrs1


subroutine mctc_ssptrs3(env, amat, bmat, ipiv, uplo)
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:)
   real(sp), intent(inout), contiguous, target :: bmat(:, :, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call mctc_sptrs(env, amat, bptr, ipiv, uplo)
end subroutine mctc_ssptrs3


subroutine mctc_dsptrs1(env, amat, bvec, ipiv, uplo)
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:)
   real(dp), intent(inout), target :: bvec(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call mctc_sptrs(env, amat, bptr, ipiv, uplo)
end subroutine mctc_dsptrs1


subroutine mctc_dsptrs3(env, amat, bmat, ipiv, uplo)
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:)
   real(dp), intent(inout), contiguous, target :: bmat(:, :, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call mctc_sptrs(env, amat, bptr, ipiv, uplo)
end subroutine mctc_dsptrs3


subroutine mctc_spotrs1(env, amat, bvec, uplo)
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout), target :: bvec(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call mctc_potrs(env, amat, bptr, uplo)
end subroutine mctc_spotrs1


subroutine mctc_spotrs3(env, amat, bmat, uplo)
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout), contiguous, target :: bmat(:, :, :)
   character(len=1), intent(in), optional :: uplo
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call mctc_potrs(env, amat, bptr, uplo)
end subroutine mctc_spotrs3


subroutine mctc_dpotrs1(env, amat, bvec, uplo)
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout), target :: bvec(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call mctc_potrs(env, amat, bptr, uplo)
end subroutine mctc_dpotrs1


subroutine mctc_dpotrs3(env, amat, bmat, uplo)
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout), contiguous, target :: bmat(:, :, :)
   character(len=1), intent(in), optional :: uplo
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call mctc_potrs(env, amat, bptr, uplo)
end subroutine mctc_dpotrs3


subroutine mctc_spptrs1(env, amat, bvec, uplo)
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:)
   real(sp), intent(inout), target :: bvec(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call mctc_pptrs(env, amat, bptr, uplo)
end subroutine mctc_spptrs1


subroutine mctc_spptrs3(env, amat, bmat, uplo)
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:)
   real(sp), intent(inout), contiguous, target :: bmat(:, :, :)
   character(len=1), intent(in), optional :: uplo
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call mctc_pptrs(env, amat, bptr, uplo)
end subroutine mctc_spptrs3


subroutine mctc_dpptrs1(env, amat, bvec, uplo)
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:)
   real(dp), intent(inout), target :: bvec(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call mctc_pptrs(env, amat, bptr, uplo)
end subroutine mctc_dpptrs1


subroutine mctc_dpptrs3(env, amat, bmat, uplo)
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:)
   real(dp), intent(inout), contiguous, target :: bmat(:, :, :)
   character(len=1), intent(in), optional :: uplo
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call mctc_pptrs(env, amat, bptr, uplo)
end subroutine mctc_dpptrs3


end module xtb_mctc_lapack_wrap
