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
module xtb_mctc_blas_wrap3
   use xtb_mctc_accuracy, only : sp, dp
   use xtb_mctc_blas_level3, only : mctc_gemm, mctc_symm, mctc_syrk, mctc_syr2k, &
      & mctc_trmm, mctc_trsm
   implicit none
   private

   public :: mctc_gemm, mctc_symm, mctc_syrk, mctc_syr2k, mctc_trmm, mctc_trsm


   !> Performs one of the matrix-matrix operations
   !
   !    C := alpha*op( A )*op( B ) + beta*C,
   !
   ! where  op( X ) is one of
   !
   !    op( X ) = X   or   op( X ) = X**T,
   !
   ! alpha and beta are scalars, and A, B and C are matrices, with op( A )
   ! an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
   interface mctc_gemm
      module procedure :: mctc_sgemm323
      module procedure :: mctc_sgemm233
      module procedure :: mctc_sgemm332
      module procedure :: mctc_dgemm323
      module procedure :: mctc_dgemm233
      module procedure :: mctc_dgemm332
   end interface mctc_gemm


contains


subroutine mctc_sgemm323(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   real(sp), pointer :: aptr(:, :), cptr(:, :)
   character(len=1) :: tra, trb
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   end if
   cptr(1:size(cmat, 1)*size(cmat, 2), 1:size(cmat, 3)) => cmat
   call mctc_gemm(aptr, bmat, cptr, tra, transb, alpha, beta)
end subroutine mctc_sgemm323


subroutine mctc_sgemm233(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in), contiguous, target :: bmat(:, :, :)
   real(sp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   real(sp), pointer :: bptr(:, :), cptr(:, :)
   character(len=1) :: trb
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   end if
   call mctc_gemm(amat, bptr, cptr, transa, trb, alpha, beta)
end subroutine mctc_sgemm233


subroutine mctc_sgemm332(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in), contiguous, target :: bmat(:, :, :)
   real(sp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   real(sp), pointer :: aptr(:, :), bptr(:, :)
   character(len=1) :: tra, trb
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   end if
   call mctc_gemm(aptr, bptr, cmat, tra, trb, alpha, beta)
end subroutine mctc_sgemm332


subroutine mctc_dgemm323(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   real(dp), pointer :: aptr(:, :), cptr(:, :)
   character(len=1) :: tra
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   end if
   cptr(1:size(cmat, 1)*size(cmat, 2), 1:size(cmat, 3)) => cmat
   call mctc_gemm(aptr, bmat, cptr, tra, transb, alpha, beta)
end subroutine mctc_dgemm323


subroutine mctc_dgemm233(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in), contiguous, target :: bmat(:, :, :)
   real(dp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   real(dp), pointer :: bptr(:, :), cptr(:, :)
   character(len=1) :: trb
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   end if
   cptr(1:size(cmat, 1), 1:size(cmat, 2)*size(cmat, 3)) => cmat
   call mctc_gemm(amat, bptr, cptr, transa, trb, alpha, beta)
end subroutine mctc_dgemm233


subroutine mctc_dgemm332(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in), contiguous, target :: bmat(:, :, :)
   real(dp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   real(dp), pointer :: aptr(:, :), bptr(:, :)
   character(len=1) :: tra, trb
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   end if
   call mctc_gemm(aptr, bptr, cmat, tra, trb, alpha, beta)
end subroutine mctc_dgemm332


end module xtb_mctc_blas_wrap3
