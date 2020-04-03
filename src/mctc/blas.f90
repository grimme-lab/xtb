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

!> Interfaces to BLAS
module xtb_mctc_blas
   use xtb_mctc_accuracy, only : sp, dp
   implicit none
   private

   !> Level 3
   public :: symm, gemm
   !> Level 2
   public :: symv, gemv
   !> Level 1
   public :: axpy, scal, copy, dot


   interface symm
      pure subroutine ssymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
         import :: sp
         real(sp),intent(in)    :: a(lda,*)
         real(sp),intent(in)    :: b(ldb,*)
         real(sp),intent(inout) :: c(ldc,*)
         character,intent(in)   :: side
         character,intent(in)   :: uplo
         real(sp),intent(in)    :: alpha
         real(sp),intent(in)    :: beta
         integer, intent(in)    :: m
         integer, intent(in)    :: n
         integer, intent(in)    :: lda
         integer, intent(in)    :: ldb
         integer, intent(in)    :: ldc
      end subroutine ssymm
      pure subroutine dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
         import :: dp
         real(dp),intent(in)    :: a(lda,*)
         real(dp),intent(in)    :: b(ldb,*)
         real(dp),intent(inout) :: c(ldc,*)
         character,intent(in)   :: side
         character,intent(in)   :: uplo
         real(dp),intent(in)    :: alpha
         real(dp),intent(in)    :: beta
         integer, intent(in)    :: m
         integer, intent(in)    :: n
         integer, intent(in)    :: lda
         integer, intent(in)    :: ldb
         integer, intent(in)    :: ldc
      end subroutine dsymm
   end interface symm

   interface gemm
      pure subroutine sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta, &
            &                     c,ldc)
         import :: sp
         real(sp),intent(in)    :: a(lda,*)
         real(sp),intent(in)    :: b(ldb,*)
         real(sp),intent(inout) :: c(ldc,*)
         character,intent(in)   :: transa
         character,intent(in)   :: transb
         real(sp),intent(in)    :: alpha
         real(sp),intent(in)    :: beta
         integer, intent(in)    :: m
         integer, intent(in)    :: n
         integer, intent(in)    :: k
         integer, intent(in)    :: lda
         integer, intent(in)    :: ldb
         integer, intent(in)    :: ldc
      end subroutine sgemm
      pure subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta, &
            &                     c,ldc)
         import :: dp
         real(dp),intent(in)    :: a(lda,*)
         real(dp),intent(in)    :: b(ldb,*)
         real(dp),intent(inout) :: c(ldc,*)
         character,intent(in)   :: transa
         character,intent(in)   :: transb
         real(dp),intent(in)    :: alpha
         real(dp),intent(in)    :: beta
         integer, intent(in)    :: m
         integer, intent(in)    :: n
         integer, intent(in)    :: k
         integer, intent(in)    :: lda
         integer, intent(in)    :: ldb
         integer, intent(in)    :: ldc
      end subroutine dgemm
   end interface gemm

   interface symv
      pure subroutine ssymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
         import :: sp
         real(sp),intent(in)    :: a(lda,*)
         real(sp),intent(in)    :: x(*)
         real(sp),intent(inout) :: y(*)
         character,intent(in)   :: uplo
         real(sp),intent(in)    :: alpha
         real(sp),intent(in)    :: beta
         integer, intent(in)    :: incx
         integer, intent(in)    :: incy
         integer, intent(in)    :: n
         integer, intent(in)    :: lda
      end subroutine ssymv
      pure subroutine dsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
         import :: dp
         real(dp),intent(in)    :: a(lda,*)
         real(dp),intent(in)    :: x(*)
         real(dp),intent(inout) :: y(*)
         character,intent(in)   :: uplo
         real(dp),intent(in)    :: alpha
         real(dp),intent(in)    :: beta
         integer, intent(in)    :: incx
         integer, intent(in)    :: incy
         integer, intent(in)    :: n
         integer, intent(in)    :: lda
      end subroutine dsymv
   end interface symv

   interface gemv
      pure subroutine sgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
         import sp
         real(sp),intent(in)    :: a(lda,*)
         real(sp),intent(in)    :: x(*)
         real(sp),intent(inout) :: y(*)
         real(sp),intent(in)    :: alpha
         real(sp),intent(in)    :: beta
         character,intent(in)   :: trans
         integer, intent(in)    :: incx
         integer, intent(in)    :: incy
         integer, intent(in)    :: m
         integer, intent(in)    :: n
         integer, intent(in)    :: lda
      end subroutine sgemv
      pure subroutine dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
         import dp
         real(dp),intent(in)    :: a(lda,*)
         real(dp),intent(in)    :: x(*)
         real(dp),intent(inout) :: y(*)
         real(dp),intent(in)    :: alpha
         real(dp),intent(in)    :: beta
         character,intent(in)   :: trans
         integer, intent(in)    :: incx
         integer, intent(in)    :: incy
         integer, intent(in)    :: m
         integer, intent(in)    :: n
         integer, intent(in)    :: lda
      end subroutine dgemv
   end interface gemv

   interface axpy
      pure subroutine saxpy(n,a,x,incx,y,incy)
         import sp
         real(sp),intent(in) :: x(*)
         real(sp),intent(inout) :: y(*)
         real(sp),intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine saxpy
      pure subroutine daxpy(n,a,x,incx,y,incy)
         import dp
         real(dp),intent(in) :: x(*)
         real(dp),intent(inout) :: y(*)
         real(dp),intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine daxpy
   end interface axpy

   interface scal
      pure subroutine sscal(n,a,x,incx)
         import sp
         real(sp),intent(inout) :: x(*)
         real(sp),intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine sscal
      pure subroutine dscal(n,a,x,incx)
         import dp
         real(dp),intent(inout) :: x(*)
         real(dp),intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine dscal
   end interface scal

   interface copy
      pure subroutine scopy(n,x,incx,y,incy)
         import sp
         real(sp),intent(in) :: x(*)
         real(sp),intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine scopy
      pure subroutine dcopy(n,x,incx,y,incy)
         import dp
         real(dp),intent(in) :: x(*)
         real(dp),intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine dcopy
   end interface copy

   interface dot
      pure function sdot(n,x,incx,y,incy)
         import sp
         real(sp) :: sdot
         real(sp),intent(in) :: x(*)
         real(sp),intent(in) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end function sdot
      pure function ddot(n,x,incx,y,incy)
         import dp
         real(dp) :: ddot
         real(dp),intent(in) :: x(*)
         real(dp),intent(in) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end function ddot
   end interface dot


end module xtb_mctc_blas
