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
module xtb_mctc_lapack_trs
   use xtb_mctc_accuracy, only : sp, dp
   use xtb_type_environment, only : TEnvironment
   implicit none
   private

   public :: mctc_getrs, mctc_sytrs, mctc_sptrs, mctc_potrs, mctc_pptrs

   public :: lapack_getrs, lapack_sytrs, lapack_sptrs, lapack_potrs, lapack_pptrs


   interface mctc_getrs
      module procedure :: mctc_sgetrs
      module procedure :: mctc_dgetrs
   end interface mctc_getrs

   interface mctc_sytrs
      module procedure :: mctc_ssytrs
      module procedure :: mctc_dsytrs
   end interface mctc_sytrs

   interface mctc_sptrs
      module procedure :: mctc_ssptrs
      module procedure :: mctc_dsptrs
   end interface mctc_sptrs

   interface mctc_potrs
      module procedure :: mctc_spotrs
      module procedure :: mctc_dpotrs
   end interface mctc_potrs

   interface mctc_pptrs
      module procedure :: mctc_spptrs
      module procedure :: mctc_dpptrs
   end interface mctc_pptrs


   interface lapack_getrs
      pure subroutine sgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         real(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: trans
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine sgetrs
      pure subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         real(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: trans
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine dgetrs
      pure subroutine cgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         complex(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: trans
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine cgetrs
      pure subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         complex(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: trans
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine zgetrs
   end interface lapack_getrs

   interface lapack_sytrs
      pure subroutine ssytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine ssytrs
      pure subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine dsytrs
      pure subroutine csytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(inout) :: b(ldb, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine csytrs
      pure subroutine zsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(inout) :: b(ldb, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine zsytrs
   end interface lapack_sytrs

   interface lapack_sptrs
      pure subroutine ssptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info)
         import :: sp
         real(sp), intent(in) :: ap(*)
         real(sp), intent(inout) :: b(ldb, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: ldb
      end subroutine ssptrs
      pure subroutine dsptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info)
         import :: dp
         real(dp), intent(in) :: ap(*)
         real(dp), intent(inout) :: b(ldb, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: ldb
      end subroutine dsptrs
      pure subroutine csptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info)
         import :: sp
         complex(sp), intent(in) :: ap(*)
         complex(sp), intent(inout) :: b(ldb, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: ldb
      end subroutine csptrs
      pure subroutine zsptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info)
         import :: dp
         complex(dp), intent(in) :: ap(*)
         complex(dp), intent(inout) :: b(ldb, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: ldb
      end subroutine zsptrs
   end interface lapack_sptrs

   interface lapack_potrs
      pure subroutine spotrs(uplo, n, nrhs, a, lda, b, ldb, info)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine spotrs
      pure subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine dpotrs
      pure subroutine cpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
         import :: sp
         complex(sp), intent(in) :: a(lda, *)
         complex(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine cpotrs
      pure subroutine zpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
         import :: dp
         complex(dp), intent(in) :: a(lda, *)
         complex(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine zpotrs
   end interface lapack_potrs

   interface lapack_pptrs
      pure subroutine spptrs(uplo, n, nrhs, ap, b, ldb, info)
         import :: sp
         real(sp), intent(in) :: ap(*)
         real(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: ldb
      end subroutine spptrs
      pure subroutine dpptrs(uplo, n, nrhs, ap, b, ldb, info)
         import :: dp
         real(dp), intent(in) :: ap(*)
         real(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: ldb
      end subroutine dpptrs
      pure subroutine cpptrs(uplo, n, nrhs, ap, b, ldb, info)
         import :: sp
         complex(sp), intent(in) :: ap(*)
         complex(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: ldb
      end subroutine cpptrs
      pure subroutine zpptrs(uplo, n, nrhs, ap, b, ldb, info)
         import :: dp
         complex(dp), intent(in) :: ap(*)
         complex(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: ldb
      end subroutine zpptrs
   end interface lapack_pptrs


contains


subroutine mctc_sgetrs(env, amat, bmat, ipiv, trans)
   character(len=*), parameter :: source = 'mctc_lapack_getrs'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: trans
   character(len=1) :: tra
   integer :: info, n, nrhs, lda, ldb
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_getrs(tra, n, nrhs, amat, lda, ipiv, bmat, ldb, info)
   if (info /= 0) then
      call env%error("Solving linear system failed", source)
   end if
end subroutine mctc_sgetrs


subroutine mctc_dgetrs(env, amat, bmat, ipiv, trans)
   character(len=*), parameter :: source = 'mctc_lapack_getrs'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: trans
   character(len=1) :: tra
   integer :: info, n, nrhs, lda, ldb
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_getrs(tra, n, nrhs, amat, lda, ipiv, bmat, ldb, info)
   if (info /= 0) then
      call env%error("Solving linear system failed", source)
   end if
end subroutine mctc_dgetrs


subroutine mctc_ssytrs(env, amat, bmat, ipiv, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sytrs'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, nrhs, lda, ldb
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_sytrs(ula, n, nrhs, amat, lda, ipiv, bmat, ldb, info)
   if (info /= 0) then
      call env%error("Solving linear system failed", source)
   end if
end subroutine mctc_ssytrs


subroutine mctc_dsytrs(env, amat, bmat, ipiv, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sytrs'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, nrhs, lda, ldb
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_sytrs(ula, n, nrhs, amat, lda, ipiv, bmat, ldb, info)
   if (info /= 0) then
      call env%error("Solving linear system failed", source)
   end if
end subroutine mctc_dsytrs


subroutine mctc_ssptrs(env, amat, bmat, ipiv, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sptrs'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:)
   real(sp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, nrhs, ldb, nn
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   ldb = max(1, size(bmat, 1))
   nn = size(amat)
   if (nn <= 0) then
      n = nn
   else
      n = int((-1+sqrt(1+8*real(nn, sp))))/2
   end if
   nrhs = size(bmat, 2)
   call lapack_sptrs(ula, n, nrhs, amat, ipiv, bmat, ldb, info)
   if (info /= 0) then
      call env%error("Solving linear system failed", source)
   end if
end subroutine mctc_ssptrs


subroutine mctc_dsptrs(env, amat, bmat, ipiv, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sptrs'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:)
   real(dp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, nrhs, ldb, nn
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   ldb = max(1, size(bmat, 1))
   nn = size(amat)
   if (nn <= 0) then
      n = nn
   else
      n = int((-1+sqrt(1+8*real(nn, dp))))/2
   end if
   nrhs = size(bmat, 2)
   call lapack_sptrs(ula, n, nrhs, amat, ipiv, bmat, ldb, info)
   if (info /= 0) then
      call env%error("Solving linear system failed", source)
   end if
end subroutine mctc_dsptrs


subroutine mctc_spotrs(env, amat, bmat, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_potrs'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, nrhs, lda, ldb
   intrinsic max, present, size
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_potrs(ula, n, nrhs, amat, lda, bmat, ldb, info)
   if (info /= 0) then
      call env%error("Solving linear system failed", source)
   end if
end subroutine mctc_spotrs


subroutine mctc_dpotrs(env, amat, bmat, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_potrs'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, nrhs, lda, ldb
   intrinsic max, present, size
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_potrs(ula, n, nrhs, amat, lda, bmat, ldb, info)
   if (info /= 0) then
      call env%error("Solving linear system failed", source)
   end if
end subroutine mctc_dpotrs


subroutine mctc_spptrs(env, amat, bmat, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_pptrs'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(in) :: amat(:)
   real(sp), intent(inout) :: bmat(:, :)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, nrhs, ldb, nn
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   ldb = max(1, size(bmat, 1))
   nn = size(amat)
   if (nn <= 0) then
      n = nn
   else
      n = int((-1+sqrt(1+8*real(nn, sp))))/2
   end if
   nrhs = size(bmat, 2)
   call lapack_pptrs(ula, n, nrhs, amat, bmat, ldb, info)
   if (info /= 0) then
      call env%error("Solving linear system failed", source)
   end if
end subroutine mctc_spptrs


subroutine mctc_dpptrs(env, amat, bmat, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_pptrs'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: amat(:)
   real(dp), intent(inout) :: bmat(:, :)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, nrhs, ldb, nn
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   ldb = max(1, size(bmat, 1))
   nn = size(amat)
   if (nn <= 0) then
      n = nn
   else
      n = int((-1+sqrt(1+8*real(nn, dp))))/2
   end if
   nrhs = size(bmat, 2)
   call lapack_pptrs(ula, n, nrhs, amat, bmat, ldb, info)
   if (info /= 0) then
      call env%error("Solving linear system failed", source)
   end if
end subroutine mctc_dpptrs


end module xtb_mctc_lapack_trs
