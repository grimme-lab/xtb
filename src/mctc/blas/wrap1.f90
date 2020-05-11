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

!> Wrappers for BLAS
module xtb_mctc_blas_wrap1
   use xtb_mctc_accuracy, only : sp, dp
   use xtb_mctc_blas_level1, only : mctc_asum, mctc_axpy, mctc_copy, mctc_dot, &
      & mctc_nrm2, mctc_rot, mctc_scal, mctc_swap, mctc_iamax
   implicit none
   private

   public :: mctc_asum, mctc_axpy, mctc_copy, mctc_dot, mctc_nrm2, mctc_rot
   public :: mctc_scal, mctc_swap, mctc_iamax

   !> Takes the sum of the absolute values.
   interface mctc_asum
      module procedure :: mctc_sasum2
      module procedure :: mctc_dasum2
   end interface mctc_asum

   !> Constant times a vector plus a vector.
   interface mctc_axpy
      module procedure :: mctc_saxpy12
      module procedure :: mctc_saxpy21
      module procedure :: mctc_saxpy22
      module procedure :: mctc_daxpy12
      module procedure :: mctc_daxpy21
      module procedure :: mctc_daxpy22
   end interface mctc_axpy

   !> Copies a vector, x, to a vector, y.
   interface mctc_copy
      module procedure :: mctc_scopy12
      module procedure :: mctc_scopy21
      module procedure :: mctc_scopy22
      module procedure :: mctc_dcopy12
      module procedure :: mctc_dcopy21
      module procedure :: mctc_dcopy22
   end interface mctc_copy

   !> Forms the dot product of two vectors.
   interface mctc_dot
      module procedure :: mctc_sdot12
      module procedure :: mctc_sdot21
      module procedure :: mctc_sdot22
      module procedure :: mctc_ddot12
      module procedure :: mctc_ddot21
      module procedure :: mctc_ddot22
   end interface mctc_dot

   !> Returns the euclidean norm of a vector via the function
   !  name, so that
   !
   !    NRM2 := sqrt( x'*x )
   interface mctc_nrm2
      module procedure :: mctc_snrm22
      module procedure :: mctc_dnrm22
   end interface mctc_nrm2

   !> Applies a plane rotation.
   interface mctc_rot
      module procedure :: mctc_srot22
      module procedure :: mctc_drot22
   end interface mctc_rot

   !> Scales a vector by a constant.
   interface mctc_scal
      module procedure :: mctc_sscal2
      module procedure :: mctc_dscal2
   end interface mctc_scal

   !> Interchanges two vectors.
   interface mctc_swap
      module procedure :: mctc_sswap22
      module procedure :: mctc_dswap22
   end interface mctc_swap

   !> Finds the index of the first element having maximum absolute value.
   interface mctc_iamax
      module procedure :: mctc_isamax2
      module procedure :: mctc_idamax2
   end interface mctc_iamax


contains


function mctc_sasum2(xvec) result(asum)
   real(sp) :: asum
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   asum = mctc_asum(xptr)
end function mctc_sasum2


function mctc_dasum2(xvec) result(asum)
   real(dp) :: asum
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   asum = mctc_asum(xptr)
end function mctc_dasum2


subroutine mctc_saxpy12(xvec, yvec, alpha)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout), contiguous, target :: yvec(:, :)
   real(sp), intent(in), optional :: alpha
   real(sp), pointer :: yptr(:)
   yptr(1:size(yvec)) => yvec
   call mctc_axpy(xvec, yptr, alpha)
end subroutine mctc_saxpy12


subroutine mctc_saxpy21(xvec, yvec, alpha)
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   call mctc_axpy(xptr, yvec, alpha)
end subroutine mctc_saxpy21


subroutine mctc_saxpy22(xvec, yvec, alpha)
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(inout), contiguous, target :: yvec(:, :)
   real(sp), intent(in), optional :: alpha
   real(sp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   call mctc_axpy(xptr, yptr, alpha)
end subroutine mctc_saxpy22


subroutine mctc_daxpy12(xvec, yvec, alpha)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), intent(in), optional :: alpha
   real(dp), pointer :: yptr(:)
   yptr(1:size(yvec)) => yvec
   call mctc_axpy(xvec, yptr, alpha)
end subroutine mctc_daxpy12


subroutine mctc_daxpy21(xvec, yvec, alpha)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   call mctc_axpy(xptr, yvec, alpha)
end subroutine mctc_daxpy21


subroutine mctc_daxpy22(xvec, yvec, alpha)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), intent(in), optional :: alpha
   real(dp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   call mctc_axpy(xptr, yptr, alpha)
end subroutine mctc_daxpy22


subroutine mctc_scopy12(xvec, yvec)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout), contiguous, target :: yvec(:, :)
   real(sp), pointer :: yptr(:)
   yptr(1:size(yvec)) => yvec
   call mctc_copy(xvec, yptr)
end subroutine mctc_scopy12


subroutine mctc_scopy21(xvec, yvec)
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(inout) :: yvec(:)
   real(sp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   call mctc_copy(xptr, yvec)
end subroutine mctc_scopy21


subroutine mctc_scopy22(xvec, yvec)
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(inout), contiguous, target :: yvec(:, :)
   real(sp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   call mctc_copy(xptr, yptr)
end subroutine mctc_scopy22


subroutine mctc_dcopy12(xvec, yvec)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), pointer :: yptr(:)
   yptr(1:size(yvec)) => yvec
   call mctc_copy(xvec, yptr)
end subroutine mctc_dcopy12


subroutine mctc_dcopy21(xvec, yvec)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(inout) :: yvec(:)
   real(dp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   call mctc_copy(xptr, yvec)
end subroutine mctc_dcopy21


subroutine mctc_dcopy22(xvec, yvec)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   call mctc_copy(xptr, yptr)
end subroutine mctc_dcopy22


function mctc_sdot12(xvec, yvec) result(dot)
   real(sp) :: dot
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(in), contiguous, target :: yvec(:, :)
   real(sp), pointer :: yptr(:)
   yptr(1:size(yvec)) => yvec
   dot = mctc_dot(xvec, yptr)
end function mctc_sdot12


function mctc_sdot21(xvec, yvec) result(dot)
   real(sp) :: dot
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(in) :: yvec(:)
   real(sp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   dot = mctc_dot(xptr, yvec)
end function mctc_sdot21


function mctc_sdot22(xvec, yvec) result(dot)
   real(sp) :: dot
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(in), contiguous, target :: yvec(:, :)
   real(sp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   dot = mctc_dot(xptr, yptr)
end function mctc_sdot22


function mctc_ddot12(xvec, yvec) result(dot)
   real(dp) :: dot
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(in), contiguous, target :: yvec(:, :)
   real(dp), pointer :: yptr(:)
   yptr(1:size(yvec)) => yvec
   dot = mctc_dot(xvec, yptr)
end function mctc_ddot12


function mctc_ddot21(xvec, yvec) result(dot)
   real(dp) :: dot
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(in) :: yvec(:)
   real(dp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   dot = mctc_dot(xptr, yvec)
end function mctc_ddot21


function mctc_ddot22(xvec, yvec) result(dot)
   real(dp) :: dot
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(in), contiguous, target :: yvec(:, :)
   real(dp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   dot = mctc_dot(xptr, yptr)
end function mctc_ddot22


function mctc_snrm22(xvec) result(nrm2)
   real(sp) :: nrm2
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   nrm2 = mctc_nrm2(xptr)
end function mctc_snrm22


function mctc_dnrm22(xvec) result(nrm2)
   real(dp) :: nrm2
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   nrm2 = mctc_nrm2(xptr)
end function mctc_dnrm22


subroutine mctc_srot22(xvec, yvec, c, s)
   real(sp), intent(in) :: c
   real(sp), intent(in) :: s
   real(sp), intent(inout), contiguous, target :: xvec(:, :)
   real(sp), intent(inout), contiguous, target :: yvec(:, :)
   real(sp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   call mctc_rot(xptr, yptr, c, s)
end subroutine mctc_srot22


subroutine mctc_drot22(xvec, yvec, c, s)
   real(dp), intent(in) :: c
   real(dp), intent(in) :: s
   real(dp), intent(inout), contiguous, target :: xvec(:, :)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   call mctc_rot(xptr, yptr, c, s)
end subroutine mctc_drot22


subroutine mctc_sscal2(xvec, alpha)
   real(sp), intent(in) :: alpha
   real(sp), intent(inout), contiguous, target :: xvec(:, :)
   real(sp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   call mctc_scal(xptr, alpha)
end subroutine mctc_sscal2


subroutine mctc_dscal2(xvec, alpha)
   real(dp), intent(in) :: alpha
   real(dp), intent(inout), contiguous, target :: xvec(:, :)
   real(dp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   call mctc_scal(xptr, alpha)
end subroutine mctc_dscal2


subroutine mctc_sswap22(xvec, yvec)
   real(sp), intent(inout), contiguous, target :: xvec(:, :)
   real(sp), intent(inout), contiguous, target :: yvec(:, :)
   real(sp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   call mctc_swap(xptr, yptr)
end subroutine mctc_sswap22


subroutine mctc_dswap22(xvec, yvec)
   real(dp), intent(inout), contiguous, target :: xvec(:, :)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   call mctc_swap(xptr, yptr)
end subroutine mctc_dswap22


function mctc_isamax2(xvec) result(iamax)
   integer :: iamax
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   iamax = mctc_iamax(xptr)
end function mctc_isamax2


function mctc_idamax2(xvec) result(iamax)
   integer :: iamax
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   iamax = mctc_iamax(xptr)
end function mctc_idamax2


end module xtb_mctc_blas_wrap1
