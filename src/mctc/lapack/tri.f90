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
module xtb_mctc_lapack_tri
   use xtb_mctc_accuracy, only : sp, dp
   implicit none
   private

   public :: mctc_getri, mctc_sytri, mctc_sptri, mctc_potri, mctc_pptri, mctc_trtri


   interface mctc_getri
      module procedure :: mctc_sgetri
      module procedure :: mctc_dgetri
   end interface mctc_getri

   interface mctc_sytri
      module procedure :: mctc_ssytri
      module procedure :: mctc_dsytri
   end interface mctc_sytri

   interface mctc_sptri
      module procedure :: mctc_ssptri
      module procedure :: mctc_dsptri
   end interface mctc_sptri

   interface mctc_potri
      module procedure :: mctc_spotri
      module procedure :: mctc_dpotri
   end interface mctc_potri

   interface mctc_pptri
      module procedure :: mctc_spptri
      module procedure :: mctc_dpptri
   end interface mctc_pptri

   interface mctc_trtri
      module procedure :: mctc_strtri
      module procedure :: mctc_dtrtri
   end interface mctc_trtri


   interface lapack_getri
      pure subroutine sgetri(n, a, lda, ipiv, work, lwork, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine sgetri
      pure subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine dgetri
      pure subroutine cgetri(n, a, lda, ipiv, work, lwork, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         complex(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine cgetri
      pure subroutine zgetri(n, a, lda, ipiv, work, lwork, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         complex(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine zgetri
   end interface lapack_getri

   interface lapack_sytri
      pure subroutine ssytri(uplo, n, a, lda, ipiv, work, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(in) :: work(*)
      end subroutine ssytri
      pure subroutine dsytri(uplo, n, a, lda, ipiv, work, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(in) :: work(*)
      end subroutine dsytri
      pure subroutine csytri(uplo, n, a, lda, ipiv, work, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         complex(sp), intent(in) :: work(*)
      end subroutine csytri
      pure subroutine zsytri(uplo, n, a, lda, ipiv, work, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         complex(dp), intent(in) :: work(*)
      end subroutine zsytri
   end interface lapack_sytri

   interface lapack_sptri
      pure subroutine ssptri(uplo, n, ap, ipiv, work, info)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         real(sp), intent(in) :: work(*)
      end subroutine ssptri
      pure subroutine dsptri(uplo, n, ap, ipiv, work, info)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         real(dp), intent(in) :: work(*)
      end subroutine dsptri
      pure subroutine csptri(uplo, n, ap, ipiv, work, info)
         import :: sp
         complex(sp), intent(inout) :: ap(*)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         complex(sp), intent(in) :: work(*)
      end subroutine csptri
      pure subroutine zsptri(uplo, n, ap, ipiv, work, info)
         import :: dp
         complex(dp), intent(inout) :: ap(*)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         complex(dp), intent(in) :: work(*)
      end subroutine zsptri
   end interface lapack_sptri

   interface lapack_trtri
      pure subroutine strtri(uplo, diag, n, a, lda, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: diag
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine strtri
      pure subroutine dtrtri(uplo, diag, n, a, lda, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: diag
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dtrtri
      pure subroutine ctrtri(uplo, diag, n, a, lda, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: diag
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine ctrtri
      pure subroutine ztrtri(uplo, diag, n, a, lda, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: diag
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine ztrtri
   end interface lapack_trtri

   interface lapack_potri
      pure subroutine spotri(uplo, n, a, lda, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine spotri
      pure subroutine dpotri(uplo, n, a, lda, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dpotri
      pure subroutine cpotri(uplo, n, a, lda, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine cpotri
      pure subroutine zpotri(uplo, n, a, lda, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine zpotri
   end interface lapack_potri

   interface lapack_pptri
      pure subroutine spptri(uplo, n, ap, info)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine spptri
      pure subroutine dpptri(uplo, n, ap, info)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine dpptri
      pure subroutine cpptri(uplo, n, ap, info)
         import :: sp
         complex(sp), intent(inout) :: ap(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine cpptri
      pure subroutine zpptri(uplo, n, ap, info)
         import :: dp
         complex(dp), intent(inout) :: ap(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine zpptri
   end interface lapack_pptri


contains


pure subroutine mctc_sgetri(amat, ipiv, info)
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(in) :: ipiv(:)
   integer, intent(out), optional :: info
   integer :: info_, n, lda, lwork, stat_alloc
   real(sp), allocatable :: work(:)
   real(sp) :: test(1)
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   lwork = -1
   call lapack_getri(n, amat, lda, ipiv, test, lwork, info_)
   if (info_ == 0) then
      lwork = max(1, nint(test(1)))
      allocate(work(lwork), stat=stat_alloc)
      if (stat_alloc == 0) then
         call lapack_getri(n, amat, lda, ipiv, work, lwork, info_)
      else
         info_ = -1000
      end if
   end if
   if (present(info)) info = info_
end subroutine mctc_sgetri


pure subroutine mctc_dgetri(amat, ipiv, info)
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(in) :: ipiv(:)
   integer, intent(out), optional :: info
   integer :: info_, n, lda, lwork, stat_alloc
   real(dp), allocatable :: work(:)
   real(dp) :: test(1)
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   lwork = -1
   call lapack_getri(n, amat, lda, ipiv, test, lwork, info_)
   if (info_ == 0) then
      lwork = max(1, nint(test(1)))
      allocate(work(lwork), stat=stat_alloc)
      if (stat_alloc == 0) then
         call lapack_getri(n, amat, lda, ipiv, work, lwork, info_)
      else
         info_ = -1000
      end if
   end if
   if (present(info)) info = info_
end subroutine mctc_dgetri


pure subroutine mctc_ssytri(amat, ipiv, uplo, info)
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: info_, n, lda, stat_alloc
   real(sp), allocatable :: work(:)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   allocate(work(max(1, n)), stat=stat_alloc)
   if (stat_alloc == 0) then
      call lapack_sytri(ula, n, amat, lda, ipiv, work, info_)
   else
      info_ = -1000
   end if
   if (present(info)) info = info_
end subroutine mctc_ssytri


pure subroutine mctc_dsytri(amat, ipiv, uplo, info)
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: info_, n, lda, stat_alloc
   real(dp), allocatable :: work(:)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   allocate(work(max(1, n)), stat=stat_alloc)
   if (stat_alloc == 0) then
      call lapack_sytri(ula, n, amat, lda, ipiv, work, info_)
   else
      info_ = -1000
   end if
   if (present(info)) info = info_
end subroutine mctc_dsytri


pure subroutine mctc_ssptri(amat, ipiv, uplo, info)
   real(sp), intent(inout) :: amat(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: info_, n, stat_alloc, nn
   real(sp), allocatable :: work(:)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   nn = size(amat)
   if (nn <= 0) then
      n = nn
   else
      n = int((-1+sqrt(1+8*real(nn, sp))))/2
   end if
   allocate(work(max(1, n)), stat=stat_alloc)
   if (stat_alloc == 0) then
      call lapack_sptri(ula, n, amat, ipiv, work, info_)
   else
      info_ = -1000
   end if
   if (present(info)) info = info_
end subroutine mctc_ssptri


pure subroutine mctc_dsptri(amat, ipiv, uplo, info)
   real(dp), intent(inout) :: amat(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: info_, n, stat_alloc, nn
   real(dp), allocatable :: work(:)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   nn = size(amat)
   if (nn <= 0) then
      n = nn
   else
      n = int((-1+sqrt(1+8*real(nn, dp))))/2
   end if
   allocate(work(max(1, n)), stat=stat_alloc)
   if (stat_alloc == 0) then
      call lapack_sptri(ula, n, amat, ipiv, work, info_)
   else
      info_ = -1000
   end if
   if (present(info)) info = info_
end subroutine mctc_dsptri


pure subroutine mctc_spotri(amat, uplo, info)
   real(sp), intent(inout) :: amat(:, :)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: info_, n, lda
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_potri(ula, n, amat, lda, info_)
   if (present(info)) info = info_
end subroutine mctc_spotri


pure subroutine mctc_dpotri(amat, uplo, info)
   real(dp), intent(inout) :: amat(:, :)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: info_, n, lda
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_potri(ula, n, amat, lda, info_)
   if (present(info)) info = info_
end subroutine mctc_dpotri


pure subroutine mctc_spptri(amat, uplo, info)
   real(sp), intent(inout) :: amat(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: info_, n, nn
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   nn = size(amat)
   if (nn <= 0) then
      n = nn
   else
      n = int((-1+sqrt(1+8*real(nn, sp))))/2
   end if
   call lapack_pptri(ula, n, amat, info_)
   if (present(info)) info = info_
end subroutine mctc_spptri


pure subroutine mctc_dpptri(amat, uplo, info)
   real(dp), intent(inout) :: amat(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: info_, n, nn
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   nn = size(amat)
   if (nn <= 0) then
      n = nn
   else
      n = int((-1+sqrt(1+8*real(nn, dp))))/2
   end if
   call lapack_pptri(ula, n, amat, info_)
   if (present(info)) info = info_
end subroutine mctc_dpptri


pure subroutine mctc_strtri(amat, uplo, diag, info)
   real(sp), intent(inout) :: amat(:, :)
   character(len=1), intent(in), optional :: uplo
   character(len=1), intent(in), optional :: diag
   integer, intent(out), optional :: info
   character(len=1) :: ula, dia
   integer :: info_, n, lda
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   if (present(diag)) then
      dia = diag
   else
      dia = 'n'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_trtri(ula, dia, n, amat, lda, info_)
   if (present(info)) info = info_
end subroutine mctc_strtri


pure subroutine mctc_dtrtri(amat, uplo, diag, info)
   real(dp), intent(inout) :: amat(:, :)
   character(len=1), intent(in), optional :: uplo
   character(len=1), intent(in), optional :: diag
   integer, intent(out), optional :: info
   character(len=1) :: ula, dia
   integer :: info_, n, lda
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   if (present(diag)) then
      dia = diag
   else
      dia = 'n'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_trtri(ula, dia, n, amat, lda, info_)
   if (present(info)) info = info_
end subroutine mctc_dtrtri


end module xtb_mctc_lapack_tri
