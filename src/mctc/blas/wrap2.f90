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
module xtb_mctc_blas_wrap2
   use xtb_mctc_accuracy, only : sp, dp
   use xtb_mctc_blas_level2, only : mctc_gemv, mctc_ger, mctc_spmv, mctc_spr, &
      & mctc_spr2, mctc_symv, mctc_syr, mctc_syr2
   implicit none
   private

   public :: mctc_gemv, mctc_ger, mctc_spmv, mctc_spr, mctc_spr2, mctc_symv
   public :: mctc_syr, mctc_syr2


   !> Performs one of the matrix-vector operations
   !
   !    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !
   ! where alpha and beta are scalars, x and y are vectors and A is an
   ! m by n matrix.
   interface mctc_gemv
      module procedure :: mctc_sgemv312
      module procedure :: mctc_sgemv321
      module procedure :: mctc_dgemv312
      module procedure :: mctc_dgemv321
   end interface mctc_gemv


contains


subroutine mctc_sgemv312(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout), contiguous, target :: yvec(:, :)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: aptr(:, :), yptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      yptr(1:size(yvec, 1) * size(yvec, 2)) => yvec
   end if
   call mctc_gemv(aptr, xvec, yptr, alpha, beta, tra)
end subroutine mctc_sgemv312


subroutine mctc_sgemv321(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call mctc_gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine mctc_sgemv321


subroutine mctc_dgemv312(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), yptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      yptr(1:size(yvec, 1) * size(yvec, 2)) => yvec
   end if
   call mctc_gemv(aptr, xvec, yptr, alpha, beta, tra)
end subroutine mctc_dgemv312


subroutine mctc_dgemv321(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call mctc_gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine mctc_dgemv321


end module xtb_mctc_blas_wrap2
