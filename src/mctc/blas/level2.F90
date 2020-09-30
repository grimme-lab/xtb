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
module xtb_mctc_blas_level2
   use xtb_mctc_accuracy, only : sp, dp
#ifdef USE_CUBLAS
   use cublas, only : cublasSgemv, cublasDgemv, cublasSger, cublasDger, &
      & cublasSspmv, cublasDspmv, cublasSspr, cublasDspr, &
      & cublasSspr2, cublasDspr2, cublasSsymv, cublasDsymv, &
      & cublasSsyr, cublasDsyr, cublasSsyr2, cublasDsyr2
#endif
   implicit none
   private

   public :: mctc_gemv, mctc_ger, mctc_spmv, mctc_spr, mctc_spr2, mctc_symv
   public :: mctc_syr, mctc_syr2

   public :: blas_gbmv, blas_gemv, blas_ger, blas_gerc, blas_geru
   public :: blas_hbmv, blas_hemv, blas_her, blas_her2, blas_hpmv, blas_hpr, blas_hpr2
   public :: blas_sbmv, blas_spmv, blas_spr, blas_spr2, blas_symv, blas_syr, blas_syr2
   public :: blas_tbmv, blas_tbsv, blas_tpmv, blas_tpsv, blas_trmv, blas_trsv


   !> Performs one of the matrix-vector operations
   !
   !    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !
   ! where alpha and beta are scalars, x and y are vectors and A is an
   ! m by n matrix.
   interface mctc_gemv
      module procedure :: mctc_sgemv
      module procedure :: mctc_dgemv
   end interface mctc_gemv

   !> Performs the rank 1 operation
   !
   !    A := alpha*x*y**T + A,
   !
   ! where alpha is a scalar, x is an m element vector, y is an n element
   ! vector and A is an m by n matrix.
   interface mctc_ger
      module procedure :: mctc_sger
      module procedure :: mctc_dger
   end interface mctc_ger

   !> Performs the matrix-vector operation
   !
   !    y := alpha*A*x + beta*y,
   !
   ! where alpha and beta are scalars, x and y are n element vectors and
   ! A is an n by n symmetric matrix, supplied in packed form.
   interface mctc_spmv
      module procedure :: mctc_sspmv
      module procedure :: mctc_dspmv
   end interface mctc_spmv

   !> Performs the symmetric rank 1 operation
   !
   !    A := alpha*x*x**T + A,
   !
   ! where alpha is a real scalar, x is an n element vector and A is an
   ! n by n symmetric matrix, supplied in packed form.
   interface mctc_spr
      module procedure :: mctc_sspr
      module procedure :: mctc_dspr
   end interface mctc_spr

   !> Performs the symmetric rank 2 operation
   !
   !    A := alpha*x*y**T + alpha*y*x**T + A,
   !
   ! where alpha is a scalar, x and y are n element vectors and A is an
   ! n by n symmetric matrix, supplied in packed form.
   interface mctc_spr2
      module procedure :: mctc_sspr2
      module procedure :: mctc_dspr2
   end interface mctc_spr2

   !> Performs the matrix-vector  operation
   !
   !    y := alpha*A*x + beta*y,
   !
   ! where alpha and beta are scalars, x and y are n element vectors and
   ! A is an n by n symmetric matrix.
   interface mctc_symv
      module procedure :: mctc_ssymv
      module procedure :: mctc_dsymv
   end interface mctc_symv

   !> Performs the symmetric rank 1 operation
   !
   !    A := alpha*x*x**T + A,
   !
   ! where alpha is a real scalar, x is an n element vector and A is an
   ! n by n symmetric matrix.
   interface mctc_syr
      module procedure :: mctc_ssyr
      module procedure :: mctc_dsyr
   end interface mctc_syr

   !> Performs the symmetric rank 2 operation
   !
   !    A := alpha*x*y**T + alpha*y*x**T + A,
   !
   ! where alpha is a scalar, x and y are n element vectors and A is an n
   ! by n symmetric matrix.
   interface mctc_syr2
      module procedure :: mctc_ssyr2
      module procedure :: mctc_dsyr2
   end interface mctc_syr2


   !> Performs one of the matrix-vector operations
   !
   !    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !
   ! where alpha and beta are scalars, x and y are vectors and A is an
   ! m by n band matrix, with kl sub-diagonals and ku super-diagonals.
   interface blas_gbmv
      pure subroutine sgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, &
            & y, incy)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         integer, intent(in) :: kl
         integer, intent(in) :: m
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: ku
         integer, intent(in) :: lda
      end subroutine sgbmv
      pure subroutine dgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, &
            & y, incy)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: kl
         integer, intent(in) :: m
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: ku
         integer, intent(in) :: lda
      end subroutine dgbmv
      pure subroutine cgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, &
            & y, incy)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(inout) :: y(*)
         integer, intent(in) :: kl
         integer, intent(in) :: m
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: ku
         integer, intent(in) :: lda
      end subroutine cgbmv
      pure subroutine zgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, &
            & y, incy)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(inout) :: y(*)
         integer, intent(in) :: kl
         integer, intent(in) :: m
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: ku
         integer, intent(in) :: lda
      end subroutine zgbmv
   end interface blas_gbmv

   !> Performs one of the matrix-vector operations
   !
   !    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !
   ! where alpha and beta are scalars, x and y are vectors and A is an
   ! m by n matrix.
   interface blas_gemv
      pure subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine sgemv
      pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dgemv
      pure subroutine cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(inout) :: y(*)
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine cgemv
      pure subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(inout) :: y(*)
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine zgemv
      pure subroutine scgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(inout) :: y(*)
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine scgemv
      pure subroutine dzgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(inout) :: y(*)
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dzgemv
   end interface blas_gemv

   !> Performs the rank 1 operation
   !
   !    A := alpha*x*y**T + A,
   !
   ! where alpha is a scalar, x is an m element vector, y is an n element
   ! vector and A is an m by n matrix.
   interface blas_ger
      pure subroutine sger(m, n, alpha, x, incx, y, incy, a, lda)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(in) :: y(*)
         real(sp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine sger
      pure subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(in) :: y(*)
         real(dp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dger
   end interface blas_ger

   !> Performs the rank 1 operation
   !
   !    A := alpha*x*y**H + A,
   !
   ! where alpha is a scalar, x is an m element vector, y is an n element
   ! vector and A is an m by n matrix.
   interface blas_gerc
      pure subroutine cgerc(m, n, alpha, x, incx, y, incy, a, lda)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(in) :: y(*)
         complex(sp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine cgerc
      pure subroutine zgerc(m, n, alpha, x, incx, y, incy, a, lda)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(in) :: y(*)
         complex(dp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine zgerc
   end interface blas_gerc

   !> Performs the rank 1 operation
   !
   !    A := alpha*x*y**T + A,
   !
   ! where alpha is a scalar, x is an m element vector, y is an n element
   ! vector and A is an m by n matrix.
   interface blas_geru
      pure subroutine cgeru(m, n, alpha, x, incx, y, incy, a, lda)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(in) :: y(*)
         complex(sp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine cgeru
      pure subroutine zgeru(m, n, alpha, x, incx, y, incy, a, lda)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(in) :: y(*)
         complex(dp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine zgeru
   end interface blas_geru

   !> Performs the matrix-vector  operation
   !
   !    y := alpha*A*x + beta*y,
   !
   ! where alpha and beta are scalars, x and y are n element vectors and
   ! A is an n by n hermitian band matrix, with k super-diagonals.
   interface blas_hbmv
      pure subroutine chbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
      end subroutine chbmv
      pure subroutine zhbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
      end subroutine zhbmv
   end interface blas_hbmv

   !> Performs the matrix-vector  operation
   !
   !    y := alpha*A*x + beta*y,
   !
   ! where alpha and beta are scalars, x and y are n element vectors and
   ! A is an n by n hermitian matrix.
   interface blas_hemv
      pure subroutine chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine chemv
      pure subroutine zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine zhemv
   end interface blas_hemv

   !> Performs the hermitian rank 1 operation
   !
   !    A := alpha*x*x**H + A,
   !
   ! where alpha is a real scalar, x is an n element vector and A is an
   ! n by n hermitian matrix.
   interface blas_her
      pure subroutine cher(uplo, n, alpha, x, incx, a, lda)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         complex(sp), intent(in) :: x(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine cher
      pure subroutine zher(uplo, n, alpha, x, incx, a, lda)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         complex(dp), intent(in) :: x(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine zher
   end interface blas_her

   !> Performs the hermitian rank 2 operation
   !
   !    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
   !
   ! where alpha is a scalar, x and y are n element vectors and A is an n
   ! by n hermitian matrix.
   interface blas_her2
      pure subroutine cher2(uplo, n, alpha, x, incx, y, incy, a, lda)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(in) :: y(*)
         character(len=1), intent(in) :: uplo
         complex(sp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine cher2
      pure subroutine zher2(uplo, n, alpha, x, incx, y, incy, a, lda)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(in) :: y(*)
         character(len=1), intent(in) :: uplo
         complex(dp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine zher2
   end interface blas_her2

   !> Performs the matrix-vector operation
   !
   !    y := alpha*A*x + beta*y,
   !
   ! where alpha and beta are scalars, x and y are n element vectors and
   ! A is an n by n hermitian matrix, supplied in packed form.
   interface blas_hpmv
      pure subroutine chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
         import :: sp
         complex(sp), intent(in) :: ap(*)
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine chpmv
      pure subroutine zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
         import :: dp
         complex(dp), intent(in) :: ap(*)
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine zhpmv
   end interface blas_hpmv

   !> Performs the hermitian rank 1 operation
   !
   !    A := alpha*x*x**H + A,
   !
   ! where alpha is a real scalar, x is an n element vector and A is an
   ! n by n hermitian matrix, supplied in packed form.
   interface blas_hpr
      pure subroutine chpr(uplo, n, alpha, x, incx, ap)
         import :: sp
         complex(sp), intent(inout) :: ap(*)
         complex(sp), intent(in) :: x(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine chpr
      pure subroutine zhpr(uplo, n, alpha, x, incx, ap)
         import :: dp
         complex(dp), intent(inout) :: ap(*)
         complex(dp), intent(in) :: x(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine zhpr
   end interface blas_hpr

   !> Performs the hermitian rank 2 operation
   !
   !    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
   !
   ! where alpha is a scalar, x and y are n element vectors and A is an
   ! n by n hermitian matrix, supplied in packed form.
   interface blas_hpr2
      pure subroutine chpr2(uplo, n, alpha, x, incx, y, incy, ap)
         import :: sp
         complex(sp), intent(inout) :: ap(*)
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(in) :: y(*)
         character(len=1), intent(in) :: uplo
         complex(sp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine chpr2
      pure subroutine zhpr2(uplo, n, alpha, x, incx, y, incy, ap)
         import :: dp
         complex(dp), intent(inout) :: ap(*)
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(in) :: y(*)
         character(len=1), intent(in) :: uplo
         complex(dp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine zhpr2
   end interface blas_hpr2

   !> Performs the matrix-vector  operation
   !
   !    y := alpha*A*x + beta*y,
   !
   ! where alpha and beta are scalars, x and y are n element vectors and
   ! A is an n by n symmetric band matrix, with k super-diagonals.
   interface blas_sbmv
      pure subroutine ssbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
      end subroutine ssbmv
      pure subroutine dsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
      end subroutine dsbmv
   end interface blas_sbmv

   !> Performs the matrix-vector operation
   !
   !    y := alpha*A*x + beta*y,
   !
   ! where alpha and beta are scalars, x and y are n element vectors and
   ! A is an n by n symmetric matrix, supplied in packed form.
   interface blas_spmv
      pure subroutine sspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
         import :: sp
         real(sp), intent(in) :: ap(*)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine sspmv
      pure subroutine dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
         import :: dp
         real(dp), intent(in) :: ap(*)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine dspmv
   end interface blas_spmv

   !> Performs the symmetric rank 1 operation
   !
   !    A := alpha*x*x**T + A,
   !
   ! where alpha is a real scalar, x is an n element vector and A is an
   ! n by n symmetric matrix, supplied in packed form.
   interface blas_spr
      pure subroutine sspr(uplo, n, alpha, x, incx, ap)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         real(sp), intent(in) :: x(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine sspr
      pure subroutine dspr(uplo, n, alpha, x, incx, ap)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         real(dp), intent(in) :: x(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine dspr
   end interface blas_spr

   !> Performs the symmetric rank 2 operation
   !
   !    A := alpha*x*y**T + alpha*y*x**T + A,
   !
   ! where alpha is a scalar, x and y are n element vectors and A is an
   ! n by n symmetric matrix, supplied in packed form.
   interface blas_spr2
      pure subroutine sspr2(uplo, n, alpha, x, incx, y, incy, ap)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         real(sp), intent(in) :: x(*)
         real(sp), intent(in) :: y(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine sspr2
      pure subroutine dspr2(uplo, n, alpha, x, incx, y, incy, ap)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         real(dp), intent(in) :: x(*)
         real(dp), intent(in) :: y(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine dspr2
   end interface blas_spr2

   !> Performs the matrix-vector  operation
   !
   !    y := alpha*A*x + beta*y,
   !
   ! where alpha and beta are scalars, x and y are n element vectors and
   ! A is an n by n symmetric matrix.
   interface blas_symv
      pure subroutine ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine ssymv
      pure subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dsymv
   end interface blas_symv

   !> Performs the symmetric rank 1 operation
   !
   !    A := alpha*x*x**T + A,
   !
   ! where alpha is a real scalar, x is an n element vector and A is an
   ! n by n symmetric matrix.
   interface blas_syr
      pure subroutine ssyr(uplo, n, alpha, x, incx, a, lda)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine ssyr
      pure subroutine dsyr(uplo, n, alpha, x, incx, a, lda)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dsyr
   end interface blas_syr

   !> Performs the symmetric rank 2 operation
   !
   !    A := alpha*x*y**T + alpha*y*x**T + A,
   !
   ! where alpha is a scalar, x and y are n element vectors and A is an n
   ! by n symmetric matrix.
   interface blas_syr2
      pure subroutine ssyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(in) :: y(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine ssyr2
      pure subroutine dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(in) :: y(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(in) :: alpha
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dsyr2
   end interface blas_syr2

   !> Performs one of the matrix-vector operations
   !
   !    x := A*x,   or   x := A**T*x,
   !
   ! where x is an n element vector and  A is an n by n unit, or non-unit,
   ! upper or lower triangular band matrix, with ( k + 1 ) diagonals.
   interface blas_tbmv
      pure subroutine stbmv(uplo, trans, diag, n, k, a, lda, x, incx)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
      end subroutine stbmv
      pure subroutine dtbmv(uplo, trans, diag, n, k, a, lda, x, incx)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
      end subroutine dtbmv
      pure subroutine ctbmv(uplo, trans, diag, n, k, a, lda, x, incx)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
      end subroutine ctbmv
      pure subroutine ztbmv(uplo, trans, diag, n, k, a, lda, x, incx)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
      end subroutine ztbmv
   end interface blas_tbmv

   !> Solves one of the systems of equations
   !
   !    A*x = b,   or   A**T*x = b,
   !
   ! where b and x are n element vectors and A is an n by n unit, or
   ! non-unit, upper or lower triangular band matrix, with ( k + 1 )
   ! diagonals.
   !
   ! No test for singularity or near-singularity is included in this
   ! routine. Such tests must be performed before calling this routine.
   interface blas_tbsv
      pure subroutine stbsv(uplo, trans, diag, n, k, a, lda, x, incx)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
      end subroutine stbsv
      pure subroutine dtbsv(uplo, trans, diag, n, k, a, lda, x, incx)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
      end subroutine dtbsv
      pure subroutine ctbsv(uplo, trans, diag, n, k, a, lda, x, incx)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
      end subroutine ctbsv
      pure subroutine ztbsv(uplo, trans, diag, n, k, a, lda, x, incx)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
      end subroutine ztbsv
   end interface blas_tbsv

   !> Performs one of the matrix-vector operations
   !
   !    x := A*x,   or   x := A**T*x,
   !
   ! where x is an n element vector and  A is an n by n unit, or non-unit,
   ! upper or lower triangular matrix, supplied in packed form.
   interface blas_tpmv
      pure subroutine stpmv(uplo, trans, diag, n, ap, x, incx)
         import :: sp
         real(sp), intent(in) :: ap(*)
         real(sp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine stpmv
      pure subroutine dtpmv(uplo, trans, diag, n, ap, x, incx)
         import :: dp
         real(dp), intent(in) :: ap(*)
         real(dp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine dtpmv
      pure subroutine ctpmv(uplo, trans, diag, n, ap, x, incx)
         import :: sp
         complex(sp), intent(in) :: ap(*)
         complex(sp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine ctpmv
      pure subroutine ztpmv(uplo, trans, diag, n, ap, x, incx)
         import :: dp
         complex(dp), intent(in) :: ap(*)
         complex(dp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine ztpmv
   end interface blas_tpmv

   !> Solves one of the systems of equations
   !
   !    A*x = b,   or   A**T*x = b,
   !
   ! where b and x are n element vectors and A is an n by n unit, or
   ! non-unit, upper or lower triangular matrix, supplied in packed form.
   !
   ! No test for singularity or near-singularity is included in this
   ! routine. Such tests must be performed before calling this routine.
   interface blas_tpsv
      pure subroutine stpsv(uplo, trans, diag, n, ap, x, incx)
         import :: sp
         real(sp), intent(in) :: ap(*)
         real(sp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine stpsv
      pure subroutine dtpsv(uplo, trans, diag, n, ap, x, incx)
         import :: dp
         real(dp), intent(in) :: ap(*)
         real(dp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine dtpsv
      pure subroutine ctpsv(uplo, trans, diag, n, ap, x, incx)
         import :: sp
         complex(sp), intent(in) :: ap(*)
         complex(sp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine ctpsv
      pure subroutine ztpsv(uplo, trans, diag, n, ap, x, incx)
         import :: dp
         complex(dp), intent(in) :: ap(*)
         complex(dp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine ztpsv
   end interface blas_tpsv

   !> Performs one of the matrix-vector operations
   !
   !     x := A*x,   or   x := A**T*x,
   !
   !  where x is an n element vector and  A is an n by n unit, or non-unit,
   !  upper or lower triangular matrix.
   interface blas_trmv
      pure subroutine strmv(uplo, trans, diag, n, a, lda, x, incx)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine strmv
      pure subroutine dtrmv(uplo, trans, diag, n, a, lda, x, incx)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dtrmv
      pure subroutine ctrmv(uplo, trans, diag, n, a, lda, x, incx)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine ctrmv
      pure subroutine ztrmv(uplo, trans, diag, n, a, lda, x, incx)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine ztrmv
   end interface blas_trmv

   !> Solves one of the systems of equations
   !
   !     A*x = b,   or   A**T*x = b,
   !
   !  where b and x are n element vectors and A is an n by n unit, or
   !  non-unit, upper or lower triangular matrix.
   !
   !  No test for singularity or near-singularity is included in this
   !  routine. Such tests must be performed before calling this routine.
   interface blas_trsv
      pure subroutine strsv(uplo, trans, diag, n, a, lda, x, incx)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine strsv
      pure subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dtrsv
      pure subroutine ctrsv(uplo, trans, diag, n, a, lda, x, incx)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine ctrsv
      pure subroutine ztrsv(uplo, trans, diag, n, a, lda, x, incx)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(inout) :: x(*)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         character(len=1), intent(in) :: diag
         integer, intent(in) :: incx
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine ztrsv
   end interface blas_trsv


contains


pure subroutine mctc_sgemv(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp) :: a, b
   character(len=1) :: tra
   integer :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec, amat)
   !$acc host_data use_device(xvec, yvec, amat)
   call cublasSgemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
   !$acc end host_data
   !$acc exit data copyout(yvec) delete(amat, xvec)
#else
   call blas_gemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
#endif
end subroutine mctc_sgemv


pure subroutine mctc_dgemv(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp) :: a, b
   character(len=1) :: tra
   integer :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec, amat)
   !$acc host_data use_device(xvec, yvec, amat)
   call cublasDgemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
   !$acc end host_data
   !$acc exit data copyout(yvec) delete(amat, xvec)
#else
   call blas_gemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
#endif
end subroutine mctc_dgemv


pure subroutine mctc_sger(amat, xvec, yvec, alpha)
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(in) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp) :: a
   integer :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec, amat)
   !$acc host_data use_device(xvec, yvec, amat)
   call cublasSger(m, n, a, xvec, incx, yvec, incy, amat, lda)
   !$acc end host_data
   !$acc exit data copyout(amat) delete(xvec, yvec)
#else
   call blas_ger(m, n, a, xvec, incx, yvec, incy, amat, lda)
#endif
end subroutine mctc_sger


pure subroutine mctc_dger(amat, xvec, yvec, alpha)
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(in) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp) :: a
   integer :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec, amat)
   !$acc host_data use_device(xvec, yvec, amat)
   call cublasDger(m, n, a, xvec, incx, yvec, incy, amat, lda)
   !$acc end host_data
   !$acc exit data copyout(amat) delete(xvec, yvec)
#else
   call blas_ger(m, n, a, xvec, incx, yvec, incy, amat, lda)
#endif
end subroutine mctc_dger


pure subroutine mctc_sspmv(amat, xvec, yvec, uplo, alpha, beta)
   character(len=1), intent(in), optional :: uplo
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   real(sp), intent(in) :: amat(:)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   character(len=1) :: ula
   real(sp) :: a, b
   integer :: incx, incy, n
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   incy = 1
   n = size(xvec)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec, amat)
   !$acc host_data use_device(xvec, yvec, amat)
   call cublasSspmv(ula, n, a, amat, xvec, incx, b, yvec, incy)
   !$acc end host_data
   !$acc exit data copyout(yvec) delete(amat, xvec)
#else
   call blas_spmv(ula, n, a, amat, xvec, incx, b, yvec, incy)
#endif
end subroutine mctc_sspmv


pure subroutine mctc_dspmv(amat, xvec, yvec, uplo, alpha, beta)
   character(len=1), intent(in), optional :: uplo
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   real(dp), intent(in) :: amat(:)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   character(len=1) :: ula
   real(dp) :: a, b
   integer :: incx, incy, n
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   incy = 1
   n = size(xvec)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec, amat)
   !$acc host_data use_device(xvec, yvec, amat)
   call cublasDspmv(ula, n, a, amat, xvec, incx, b, yvec, incy)
   !$acc end host_data
   !$acc exit data copyout(yvec) delete(amat, xvec)
#else
   call blas_spmv(ula, n, a, amat, xvec, incx, b, yvec, incy)
#endif
end subroutine mctc_dspmv


pure subroutine mctc_sspr(amat, xvec, uplo, alpha)
   real(sp), intent(inout) :: amat(:)
   real(sp), intent(in) :: xvec(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), intent(in), optional :: alpha
   character(len=1) :: ula
   real(sp) :: a
   integer :: incx, n
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   n = size(xvec)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, amat)
   !$acc host_data use_device(xvec, amat)
   call cublasSspr(ula, n, a, xvec, incx, amat)
   !$acc end host_data
   !$acc exit data copyout(amat) delete(xvec)
#else
   call blas_spr(ula, n, a, xvec, incx, amat)
#endif
end subroutine mctc_sspr


pure subroutine mctc_dspr(amat, xvec, uplo, alpha)
   real(dp), intent(inout) :: amat(:)
   real(dp), intent(in) :: xvec(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), intent(in), optional :: alpha
   character(len=1) :: ula
   real(dp) :: a
   integer :: incx, n
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   n = size(xvec)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, amat)
   !$acc host_data use_device(xvec, amat)
   call cublasDspr(ula, n, a, xvec, incx, amat)
   !$acc end host_data
   !$acc exit data copyout(amat) delete(xvec)
#else
   call blas_spr(ula, n, a, xvec, incx, amat)
#endif
end subroutine mctc_dspr


pure subroutine mctc_sspr2(amat, xvec, yvec, uplo, alpha)
   real(sp), intent(inout) :: amat(:)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(in) :: yvec(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), intent(in), optional :: alpha
   character(len=1) :: ula
   real(sp) :: a
   integer :: incx, incy, n
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   incy = 1
   n = size(xvec)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec, amat)
   !$acc host_data use_device(xvec, yvec, amat)
   call cublasSspr2(ula, n, a, xvec, incx, yvec, incy, amat)
   !$acc end host_data
   !$acc exit data copyout(amat) delete(xvec, yvec)
#else
   call blas_spr2(ula, n, a, xvec, incx, yvec, incy, amat)
#endif
end subroutine mctc_sspr2


pure subroutine mctc_dspr2(amat, xvec, yvec, uplo, alpha)
   real(dp), intent(inout) :: amat(:)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(in) :: yvec(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), intent(in), optional :: alpha
   character(len=1) :: ula
   real(dp) :: a
   integer :: incx, incy, n
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   incy = 1
   n = size(xvec)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec, amat)
   !$acc host_data use_device(xvec, yvec, amat)
   call cublasDspr2(ula, n, a, xvec, incx, yvec, incy, amat)
   !$acc end host_data
   !$acc exit data copyout(amat) delete(xvec, yvec)
#else
   call blas_spr2(ula, n, a, xvec, incx, yvec, incy, amat)
#endif
end subroutine mctc_dspr2


pure subroutine mctc_ssymv(amat, xvec, yvec, uplo, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1) :: ula
   real(sp) :: a, b
   integer :: incx, incy, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec, amat)
   !$acc host_data use_device(xvec, yvec, amat)
   call cublasSsymv(ula, n, a, amat, lda, xvec, incx, b, yvec, incy)
   !$acc end host_data
   !$acc exit data copyout(yvec) delete(amat, xvec)
#else
   call blas_symv(ula, n, a, amat, lda, xvec, incx, b, yvec, incy)
#endif
end subroutine mctc_ssymv


pure subroutine mctc_dsymv(amat, xvec, yvec, uplo, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1) :: ula
   real(dp) :: a, b
   integer :: incx, incy, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec, amat)
   !$acc host_data use_device(xvec, yvec, amat)
   call cublasDsymv(ula, n, a, amat, lda, xvec, incx, b, yvec, incy)
   !$acc end host_data
   !$acc exit data copyout(yvec) delete(amat, xvec)
#else
   call blas_symv(ula, n, a, amat, lda, xvec, incx, b, yvec, incy)
#endif
end subroutine mctc_dsymv


pure subroutine mctc_ssyr(amat, xvec, uplo, alpha)
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(in) :: xvec(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), intent(in), optional :: alpha
   character(len=1) :: ula
   real(sp) :: a
   integer :: incx, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, amat)
   !$acc host_data use_device(xvec, amat)
   call cublasSsyr(ula, n, a, xvec, incx, amat, lda)
   !$acc end host_data
   !$acc exit data copyout(amat) delete(xvec)
#else
   call blas_syr(ula, n, a, xvec, incx, amat, lda)
#endif
end subroutine mctc_ssyr


pure subroutine mctc_dsyr(amat, xvec, uplo, alpha)
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(in) :: xvec(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), intent(in), optional :: alpha
   character(len=1) :: ula
   real(dp) :: a
   integer :: incx, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, amat)
   !$acc host_data use_device(xvec, amat)
   call cublasDsyr(ula, n, a, xvec, incx, amat, lda)
   !$acc end host_data
   !$acc exit data copyout(amat) delete(xvec)
#else
   call blas_syr(ula, n, a, xvec, incx, amat, lda)
#endif
end subroutine mctc_dsyr


pure subroutine mctc_ssyr2(amat, xvec, yvec, uplo, alpha)
   character(len=1), intent(in), optional :: uplo
   real(sp), intent(in), optional :: alpha
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(in) :: yvec(:)
   character(len=1) :: ula
   real(sp) :: a
   integer :: incx, incy, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec, amat)
   !$acc host_data use_device(xvec, yvec, amat)
   call cublasSsyr2(ula, n, a, xvec, incx, yvec, incy, amat, lda)
   !$acc end host_data
   !$acc exit data copyout(amat) delete(xvec, yvec)
#else
   call blas_syr2(ula, n, a, xvec, incx, yvec, incy, amat, lda)
#endif
end subroutine mctc_ssyr2


pure subroutine mctc_dsyr2(amat, xvec, yvec, uplo, alpha)
   character(len=1), intent(in), optional :: uplo
   real(dp), intent(in), optional :: alpha
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(in) :: yvec(:)
   character(len=1) :: ula
   real(dp) :: a
   integer :: incx, incy, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   incx = 1
   incy = 1
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec, amat)
   !$acc host_data use_device(xvec, yvec, amat)
   call cublasDsyr2(ula, n, a, xvec, incx, yvec, incy, amat, lda)
   !$acc end host_data
   !$acc exit data copyout(amat) delete(xvec, yvec)
#else
   call blas_syr2(ula, n, a, xvec, incx, yvec, incy, amat, lda)
#endif
end subroutine mctc_dsyr2


end module xtb_mctc_blas_level2
