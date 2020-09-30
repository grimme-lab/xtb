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
module xtb_mctc_blas_level3
   use xtb_mctc_accuracy, only : sp, dp
#ifdef USE_CUBLAS
   use cublas, only : cublasSgemm, cublasDgemm, cublasSsymm, cublasDsymm, &
      & cublasSsyrk, cublasDsyrk, cublasSsyr2k, cublasDsyr2k, &
      & cublasStrmm, cublasDtrmm, cublasStrsm, cublasDtrsm
#endif
   implicit none
   private

   public :: mctc_gemm, mctc_symm, mctc_syrk, mctc_syr2k, mctc_trmm, mctc_trsm

   public :: blas_gemm, blas_hemm, blas_herk, blas_her2k
   public :: blas_symm, blas_syrk, blas_syr2k, blas_trmm, blas_trsm


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
      module procedure :: mctc_sgemm
      module procedure :: mctc_dgemm
   end interface mctc_gemm

   !> Performs one of the matrix-matrix operations
   !
   !    C := alpha*A*B + beta*C,
   !
   ! or
   !
   !    C := alpha*B*A + beta*C,
   !
   ! where alpha and beta are scalars,  A is a symmetric matrix and  B and
   ! C are  m by n matrices.
   interface mctc_symm
      module procedure :: mctc_ssymm
      module procedure :: mctc_dsymm
   end interface mctc_symm

   !> Performs one of the symmetric rank k operations
   !
   !    C := alpha*A*A**T + beta*C,
   !
   ! or
   !
   !    C := alpha*A**T*A + beta*C,
   !
   ! where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
   ! and  A  is an  n by k  matrix in the first case and a  k by n  matrix
   ! in the second case.
   interface mctc_syrk
      module procedure :: mctc_ssyrk
      module procedure :: mctc_dsyrk
   end interface mctc_syrk

   !> Performs one of the symmetric rank 2k operations
   !
   !    C := alpha*A*B**T + alpha*B*A**T + beta*C,
   !
   ! or
   !
   !    C := alpha*A**T*B + alpha*B**T*A + beta*C,
   !
   ! where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
   ! and  A and B  are  n by k  matrices  in the  first  case  and  k by n
   ! matrices in the second case.
   interface mctc_syr2k
      module procedure :: mctc_ssyr2k
      module procedure :: mctc_dsyr2k
   end interface mctc_syr2k

   !> Performs one of the matrix-matrix operations
   !
   !    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
   !
   ! where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
   ! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
   !
   !    op( A ) = A   or   op( A ) = A**T.
   interface mctc_trmm
      module procedure :: mctc_strmm
      module procedure :: mctc_dtrmm
   end interface mctc_trmm

   !> Solves one of the matrix equations
   !
   !    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
   !
   ! where alpha is a scalar, X and B are m by n matrices, A is a unit, or
   ! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
   !
   !    op( A ) = A   or   op( A ) = A**T.
   !
   ! The matrix X is overwritten on B.
   interface mctc_trsm
      module procedure :: mctc_strsm
      module procedure :: mctc_dtrsm
   end interface mctc_trsm


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
   interface blas_gemm
      pure subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            & beta, c, ldc)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: b(ldb, *)
         real(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine sgemm
      pure subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            & beta, c, ldc)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: b(ldb, *)
         real(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine dgemm
      pure subroutine cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            & beta, c, ldc)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(in) :: b(ldb, *)
         complex(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine cgemm
      pure subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            & beta, c, ldc)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(in) :: b(ldb, *)
         complex(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine zgemm
      pure subroutine scgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            & beta, c, ldc)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         complex(sp), intent(in) :: b(ldb, *)
         complex(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine scgemm
      pure subroutine dzgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            & beta, c, ldc)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         complex(dp), intent(in) :: b(ldb, *)
         complex(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine dzgemm
   end interface blas_gemm

   !> Performs one of the matrix-matrix operations
   !
   !    C := alpha*A*B + beta*C,
   !
   ! or
   !
   !    C := alpha*B*A + beta*C,
   !
   ! where alpha and beta are scalars, A is an hermitian matrix and  B and
   ! C are m by n matrices.
   interface blas_hemm
      pure subroutine chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(in) :: b(ldb, *)
         complex(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine chemm
      pure subroutine zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(in) :: b(ldb, *)
         complex(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine zhemm
   end interface blas_hemm

   !> Performs one of the hermitian rank k operations
   !
   !    C := alpha*A*A**H + beta*C,
   !
   ! or
   !
   !    C := alpha*A**H*A + beta*C,
   !
   ! where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
   ! matrix and  A  is an  n by k  matrix in the  first case and a  k by n
   ! matrix in the second case.
   interface blas_herk
      pure subroutine cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldc
      end subroutine cherk
      pure subroutine zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldc
      end subroutine zherk
   end interface blas_herk

   !> Performs one of the hermitian rank 2k operations
   !
   !    C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
   !
   ! or
   !
   !    C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
   !
   ! where  alpha and beta  are scalars with  beta  real,  C is an  n by n
   ! hermitian matrix and  A and B  are  n by k matrices in the first case
   ! and  k by n  matrices in the second case.
   interface blas_her2k
      pure subroutine cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, &
            & ldc)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(in) :: b(ldb, *)
         complex(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         complex(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine cher2k
      pure subroutine zher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, &
            & ldc)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(in) :: b(ldb, *)
         complex(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         complex(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine zher2k
   end interface blas_her2k

   !> Performs one of the matrix-matrix operations
   !
   !    C := alpha*A*B + beta*C,
   !
   ! or
   !
   !    C := alpha*B*A + beta*C,
   !
   ! where alpha and beta are scalars,  A is a symmetric matrix and  B and
   ! C are  m by n matrices.
   interface blas_symm
      pure subroutine ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: b(ldb, *)
         real(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine ssymm
      pure subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: b(ldb, *)
         real(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine dsymm
      pure subroutine csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(in) :: b(ldb, *)
         complex(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine csymm
      pure subroutine zsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(in) :: b(ldb, *)
         complex(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine zsymm
   end interface blas_symm

   !> Performs one of the symmetric rank k operations
   !
   !    C := alpha*A*A**T + beta*C,
   !
   ! or
   !
   !    C := alpha*A**T*A + beta*C,
   !
   ! where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
   ! and  A  is an  n by k  matrix in the first case and a  k by n  matrix
   ! in the second case.
   interface blas_syrk
      pure subroutine ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldc
      end subroutine ssyrk
      pure subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldc
      end subroutine dsyrk
      pure subroutine csyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldc
      end subroutine csyrk
      pure subroutine zsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldc
      end subroutine zsyrk
   end interface blas_syrk

   !> Performs one of the symmetric rank 2k operations
   !
   !    C := alpha*A*B**T + alpha*B*A**T + beta*C,
   !
   ! or
   !
   !    C := alpha*A**T*B + alpha*B**T*A + beta*C,
   !
   ! where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
   ! and  A and B  are  n by k  matrices  in the  first  case  and  k by n
   ! matrices in the second case.
   interface blas_syr2k
      pure subroutine ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, &
            & ldc)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: b(ldb, *)
         real(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine ssyr2k
      pure subroutine dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, &
            & ldc)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: b(ldb, *)
         real(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine dsyr2k
      pure subroutine csyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, &
            & ldc)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(in) :: b(ldb, *)
         complex(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         complex(sp), intent(in) :: alpha
         complex(sp), intent(in) :: beta
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine csyr2k
      pure subroutine zsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, &
            & ldc)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(in) :: b(ldb, *)
         complex(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
         integer, intent(in) :: n
         integer, intent(in) :: k
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         integer, intent(in) :: ldc
      end subroutine zsyr2k
   end interface blas_syr2k

   !> Performs one of the matrix-matrix operations
   !
   !    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
   !
   ! where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
   ! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
   !
   !    op( A ) = A   or   op( A ) = A**T.
   interface blas_trmm
      pure subroutine strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: diag
         real(sp), intent(in) :: alpha
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine strmm
      pure subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: diag
         real(dp), intent(in) :: alpha
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine dtrmm
      pure subroutine ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: diag
         complex(sp), intent(in) :: alpha
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine ctrmm
      pure subroutine ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: diag
         complex(dp), intent(in) :: alpha
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine ztrmm
   end interface blas_trmm

   !> Solves one of the matrix equations
   !
   !    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
   !
   ! where alpha is a scalar, X and B are m by n matrices, A is a unit, or
   ! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
   !
   !    op( A ) = A   or   op( A ) = A**T.
   !
   ! The matrix X is overwritten on B.
   interface blas_trsm
      pure subroutine strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: diag
         real(sp), intent(in) :: alpha
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine strsm
      pure subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: diag
         real(dp), intent(in) :: alpha
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine dtrsm
      pure subroutine ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: diag
         complex(sp), intent(in) :: alpha
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine ctrsm
      pure subroutine ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: diag
         complex(dp), intent(in) :: alpha
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine ztrsm
   end interface blas_trsm


contains


pure subroutine mctc_sgemm(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1) :: tra, trb
   real(sp) :: a, b
   integer :: m, n, k, lda, ldb, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_sp
   end if
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
   if ((tra.eq.'n'.or.tra.eq.'N')) then
      k = size(amat, 2)
   else
      k = size(amat, 1)
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   m = size(cmat, 1)
   n = size(cmat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(amat, bmat, cmat)
   !$acc host_data use_device(amat, bmat, cmat)
   call cublasSgemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
   !$acc end host_data
   !$acc exit data copyout(cmat) delete(amat, bmat)
#else
   call blas_gemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
#endif
end subroutine mctc_sgemm


pure subroutine mctc_dgemm(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1) :: tra, trb
   real(dp) :: a, b
   integer :: m, n, k, lda, ldb, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_dp
   end if
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
   if ((tra.eq.'n'.or.tra.eq.'N')) then
      k = size(amat, 2)
   else
      k = size(amat, 1)
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   m = size(cmat, 1)
   n = size(cmat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(amat, bmat, cmat)
   !$acc host_data use_device(amat, bmat, cmat)
   call cublasDgemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
   !$acc end host_data
   !$acc exit data copyout(cmat) delete(amat, bmat)
#else
   call blas_gemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
#endif
end subroutine mctc_dgemm


pure subroutine mctc_ssymm(amat, bmat, cmat, side, uplo, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: side
   character(len=1), intent(in), optional :: uplo
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1) :: lra, ula
   real(sp) :: a, b
   integer :: m, n, lda, ldb, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_sp
   end if
   if (present(side)) then
      lra = side
   else
      lra = 'l'
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   m = size(cmat, 1)
   n = size(cmat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(amat, bmat, cmat)
   !$acc host_data use_device(amat, bmat, cmat)
   call cublasSsymm(lra, ula, m, n, a, amat, lda, bmat, ldb, b, cmat, ldc)
   !$acc end host_data
   !$acc exit data copyout(cmat) delete(amat, bmat)
#else
   call blas_symm(lra, ula, m, n, a, amat, lda, bmat, ldb, b, cmat, ldc)
#endif
end subroutine mctc_ssymm


pure subroutine mctc_dsymm(amat, bmat, cmat, side, uplo, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: side
   character(len=1), intent(in), optional :: uplo
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1) :: lra, ula
   real(dp) :: a, b
   integer :: m, n, lda, ldb, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_dp
   end if
   if (present(side)) then
      lra = side
   else
      lra = 'l'
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   m = size(cmat, 1)
   n = size(cmat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(amat, bmat, cmat)
   !$acc host_data use_device(amat, bmat, cmat)
   call cublasDsymm(lra, ula, m, n, a, amat, lda, bmat, ldb, b, cmat, ldc)
   !$acc end host_data
   !$acc exit data copyout(cmat) delete(amat, bmat)
#else
   call blas_symm(lra, ula, m, n, a, amat, lda, bmat, ldb, b, cmat, ldc)
#endif
end subroutine mctc_dsymm


pure subroutine mctc_ssyrk(amat, cmat, uplo, trans, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: uplo
   character(len=1), intent(in), optional :: trans
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1) :: ula, tra
   real(sp) :: a, b
   integer :: n, k, lda, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_sp
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   if ((tra.eq.'n'.or.tra.eq.'N')) then
      k = size(amat, 2)
   else
      k = size(amat, 1)
   end if
   lda = max(1, size(amat, 1))
   ldc = max(1, size(cmat, 1))
   n = size(cmat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(amat, cmat)
   !$acc host_data use_device(amat, cmat)
   call cublasSsyrk(ula, tra, n, k, a, amat, lda, b, cmat, ldc)
   !$acc end host_data
   !$acc exit data copyout(cmat) delete(amat)
#else
   call blas_syrk(ula, tra, n, k, a, amat, lda, b, cmat, ldc)
#endif
end subroutine mctc_ssyrk


pure subroutine mctc_dsyrk(amat, cmat, uplo, trans, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: uplo
   character(len=1), intent(in), optional :: trans
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1) :: ula, tra
   real(dp) :: a, b
   integer :: n, k, lda, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_dp
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   if ((tra.eq.'n'.or.tra.eq.'N')) then
      k = size(amat, 2)
   else
      k = size(amat, 1)
   end if
   lda = max(1, size(amat, 1))
   ldc = max(1, size(cmat, 1))
   n = size(cmat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(amat, cmat)
   !$acc host_data use_device(amat, cmat)
   call cublasDsyrk(ula, tra, n, k, a, amat, lda, b, cmat, ldc)
   !$acc end host_data
   !$acc exit data copyout(cmat) delete(amat)
#else
   call blas_syrk(ula, tra, n, k, a, amat, lda, b, cmat, ldc)
#endif
end subroutine mctc_dsyrk


pure subroutine mctc_ssyr2k(amat, bmat, cmat, uplo, trans, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: uplo
   character(len=1), intent(in), optional :: trans
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1) :: ula, tra
   real(sp) :: a, b
   integer :: n, k, lda, ldb, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_sp
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   if ((tra.eq.'n'.or.tra.eq.'N')) then
      k = size(amat, 2)
   else
      k = size(amat, 1)
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   n = size(cmat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(amat, bmat, cmat)
   !$acc host_data use_device(amat, bmat, cmat)
   call cublasSsyr2k(ula, tra, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
   !$acc end host_data
   !$acc exit data copyout(cmat) delete(amat, bmat)
#else
   call blas_syr2k(ula, tra, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
#endif
end subroutine mctc_ssyr2k


pure subroutine mctc_dsyr2k(amat, bmat, cmat, uplo, trans, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: uplo
   character(len=1), intent(in), optional :: trans
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1) :: ula, tra
   real(dp) :: a, b
   integer :: n, k, lda, ldb, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_dp
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   if ((tra.eq.'n'.or.tra.eq.'N')) then
      k = size(amat, 2)
   else
      k = size(amat, 1)
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   n = size(cmat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(amat, bmat, cmat)
   !$acc host_data use_device(amat, bmat, cmat)
   call cublasDsyr2k(ula, tra, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
   !$acc end host_data
   !$acc exit data copyout(cmat) delete(amat, bmat)
#else
   call blas_syr2k(ula, tra, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
#endif
end subroutine mctc_dsyr2k


pure subroutine mctc_strsm(amat, bmat, side, uplo, transa, diag, alpha)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   character(len=1), intent(in), optional :: side
   character(len=1), intent(in), optional :: uplo
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: diag
   real(sp), intent(in), optional :: alpha
   character(len=1) :: lra, ula, tra, dia
   real(sp) :: a
   integer :: m, n, lda, ldb
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(diag)) then
      dia = diag
   else
      dia = 'n'
   end if
   if (present(side)) then
      lra = side
   else
      lra = 'l'
   end if
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   m = size(bmat, 1)
   n = size(bmat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(amat, bmat)
   !$acc host_data use_device(amat, bmat)
   call cublasStrsm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
   !$acc end host_data
   !$acc exit data copyout(bmat) delete(amat)
#else
   call blas_trsm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
#endif
end subroutine mctc_strsm


pure subroutine mctc_dtrsm(amat, bmat, side, uplo, transa, diag, alpha)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   character(len=1), intent(in), optional :: side
   character(len=1), intent(in), optional :: uplo
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: diag
   real(dp), intent(in), optional :: alpha
   character(len=1) :: lra, ula, tra, dia
   real(dp) :: a
   integer :: m, n, lda, ldb
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(diag)) then
      dia = diag
   else
      dia = 'n'
   end if
   if (present(side)) then
      lra = side
   else
      lra = 'l'
   end if
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   m = size(bmat, 1)
   n = size(bmat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(amat, bmat)
   !$acc host_data use_device(amat, bmat)
   call cublasDtrsm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
   !$acc end host_data
   !$acc exit data copyout(bmat) delete(amat)
#else
   call blas_trsm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
#endif
end subroutine mctc_dtrsm


pure subroutine mctc_strmm(amat, bmat, side, uplo, transa, diag, alpha)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   character(len=1), intent(in), optional :: side
   character(len=1), intent(in), optional :: uplo
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: diag
   real(sp), intent(in), optional :: alpha
   character(len=1) :: lra, ula, tra, dia
   real(sp) :: a
   integer :: m, n, lda, ldb
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(diag)) then
      dia = diag
   else
      dia = 'n'
   end if
   if (present(side)) then
      lra = side
   else
      lra = 'l'
   end if
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   m = size(bmat, 1)
   n = size(bmat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(amat, bmat)
   !$acc host_data use_device(amat, bmat)
   call cublasStrmm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
   !$acc end host_data
   !$acc exit data copyout(bmat) delete(amat)
#else
   call blas_trmm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
#endif
end subroutine mctc_strmm


pure subroutine mctc_dtrmm(amat, bmat, side, uplo, transa, diag, alpha)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   character(len=1), intent(in), optional :: side
   character(len=1), intent(in), optional :: uplo
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: diag
   real(dp), intent(in), optional :: alpha
   character(len=1) :: lra, ula, tra, dia
   real(dp) :: a
   integer :: m, n, lda, ldb
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(diag)) then
      dia = diag
   else
      dia = 'n'
   end if
   if (present(side)) then
      lra = side
   else
      lra = 'l'
   end if
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   m = size(bmat, 1)
   n = size(bmat, 2)
#ifdef USE_CUBLAS
   !$acc enter data copyin(amat, bmat)
   !$acc host_data use_device(amat, bmat)
   call cublasDtrmm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
   !$acc end host_data
   !$acc exit data copyout(bmat) delete(amat)
#else
   call blas_trmm(lra, ula, tra, dia, m, n, a, amat, lda, bmat, ldb)
#endif
end subroutine mctc_dtrmm


end module xtb_mctc_blas_level3
