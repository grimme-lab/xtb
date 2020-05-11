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
   use xtb_type_environment, only : TEnvironment
   implicit none
   private

   public :: mctc_getri, mctc_sytri, mctc_sptri, mctc_potri, mctc_pptri

   public :: lapack_getri, lapack_sytri, lapack_sptri, lapack_potri, lapack_pptri


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


subroutine mctc_sgetri(env, amat, ipiv)
   character(len=*), parameter :: source = 'mctc_lapack_getri'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(in) :: ipiv(:)
   integer :: info, n, lda, lwork, stat_alloc, stat_dealloc
   real(sp), allocatable :: work(:)
   real(sp) :: test(1)
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   lwork = -1
   call lapack_getri(n, amat, lda, ipiv, test, lwork, info)
   if (info == 0) then
      lwork = nint(test(1))
      allocate(work(lwork), stat=stat_alloc)
      if (stat_alloc==0) then
         call lapack_getri(n, amat, lda, ipiv, work, lwork, info)
      else
         info = -1000
      end if
      deallocate(work, stat=stat_dealloc)
   end if
   if (info /= 0) then
      call env%error("Inverting matrix failed", source)
   end if
end subroutine mctc_sgetri


subroutine mctc_dgetri(env, amat, ipiv)
   character(len=*), parameter :: source = 'mctc_lapack_getri'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(in) :: ipiv(:)
   integer :: info, n, lda, lwork, stat_alloc, stat_dealloc
   real(dp), allocatable :: work(:)
   real(dp) :: test(1)
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   lwork = -1
   call lapack_getri(n, amat, lda, ipiv, test, lwork, info)
   if (info == 0) then
      lwork = nint(test(1))
      allocate(work(lwork), stat=stat_alloc)
      if (stat_alloc==0) then
         call lapack_getri(n, amat, lda, ipiv, work, lwork, info)
      else
         info = -1000
      end if
      deallocate(work, stat=stat_dealloc)
   end if
   if (info /= 0) then
      call env%error("Inverting matrix failed", source)
   end if
end subroutine mctc_dgetri


subroutine mctc_ssytri(env, amat, ipiv, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sytri'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, lda, stat_alloc, stat_dealloc
   real(sp), allocatable :: work(:)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   allocate(work(n), stat=stat_alloc)
   if (stat_alloc==0) then
      call lapack_sytri(ula, n, amat, lda, ipiv, work, info)
   else
      info = -1000
   end if
   deallocate(work, stat=stat_dealloc)
   if (info /= 0) then
      call env%error("Inverting matrix failed", source)
   end if
end subroutine mctc_ssytri


subroutine mctc_dsytri(env, amat, ipiv, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sytri'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, lda, stat_alloc, stat_dealloc
   real(dp), allocatable :: work(:)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   allocate(work(n), stat=stat_alloc)
   if (stat_alloc==0) then
      call lapack_sytri(ula, n, amat, lda, ipiv, work, info)
   else
      info = -1000
   end if
   deallocate(work, stat=stat_dealloc)
   if (info /= 0) then
      call env%error("Inverting matrix failed", source)
   end if
end subroutine mctc_dsytri


subroutine mctc_ssptri(env, amat, ipiv, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sptri'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, stat_alloc, stat_dealloc, nn
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
   stat_alloc = 0
   allocate(work(n), stat=stat_alloc)
   if (stat_alloc==0) then
      call lapack_sptri(ula, n, amat, ipiv, work, info)
   else
      info = -1000
   end if
   deallocate(work, stat=stat_dealloc)
   if (info /= 0) then
      call env%error("Inverting matrix failed", source)
   end if
end subroutine mctc_ssptri


subroutine mctc_dsptri(env, amat, ipiv, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sptri'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, stat_alloc, stat_dealloc, nn
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
   stat_alloc = 0
   allocate(work(n), stat=stat_alloc)
   if (stat_alloc==0) then
      call lapack_sptri(ula, n, amat, ipiv, work, info)
   else
      info = -1000
   end if
   deallocate(work, stat=stat_dealloc)
   if (info /= 0) then
      call env%error("Inverting matrix failed", source)
   end if
end subroutine mctc_dsptri


subroutine mctc_spotri(env, amat, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_potri'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, lda
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_potri(ula, n, amat, lda, info)
   if (info /= 0) then
      call env%error("Inverting matrix failed", source)
   end if
end subroutine mctc_spotri


subroutine mctc_dpotri(env, amat, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_potri'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, lda
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_potri(ula, n, amat, lda, info)
   if (info /= 0) then
      call env%error("Inverting matrix failed", source)
   end if
end subroutine mctc_dpotri


subroutine mctc_spptri(env, amat, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_pptri'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, nn
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
   call lapack_pptri(ula, n, amat, info)
   if (info /= 0) then
      call env%error("Inverting matrix failed", source)
   end if
end subroutine mctc_spptri


subroutine mctc_dpptri(env, amat, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_pptri'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, nn
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
   call lapack_pptri(ula, n, amat, info)
   if (info /= 0) then
      call env%error("Inverting matrix failed", source)
   end if
end subroutine mctc_dpptri


end module xtb_mctc_lapack_tri
