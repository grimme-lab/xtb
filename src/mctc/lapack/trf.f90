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
module xtb_mctc_lapack_trf
   use xtb_mctc_accuracy, only : sp, dp
   use xtb_type_environment, only : TEnvironment
   implicit none
   private

   public :: mctc_getrf, mctc_sytrf, mctc_sptrf, mctc_potrf, mctc_pptrf

   public :: lapack_getrf, lapack_sytrf, lapack_sptrf, lapack_potrf, lapack_pptrf


   interface mctc_getrf
      module procedure :: mctc_sgetrf
      module procedure :: mctc_dgetrf
   end interface mctc_getrf

   interface mctc_sytrf
      module procedure :: mctc_ssytrf
      module procedure :: mctc_dsytrf
   end interface mctc_sytrf

   interface mctc_sptrf
      module procedure :: mctc_ssptrf
      module procedure :: mctc_dsptrf
   end interface mctc_sptrf

   interface mctc_potrf
      module procedure :: mctc_spotrf
      module procedure :: mctc_dpotrf
   end interface mctc_potrf

   interface mctc_pptrf
      module procedure :: mctc_spptrf
      module procedure :: mctc_dpptrf
   end interface mctc_pptrf


   interface lapack_getrf
      pure subroutine sgetrf(m, n, a, lda, ipiv, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine sgetrf
      pure subroutine dgetrf(m, n, a, lda, ipiv, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dgetrf
      pure subroutine cgetrf(m, n, a, lda, ipiv, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine cgetrf
      pure subroutine zgetrf(m, n, a, lda, ipiv, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine zgetrf
   end interface lapack_getrf

   interface lapack_sytrf
      pure subroutine ssytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine ssytrf
      pure subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine dsytrf
      pure subroutine csytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         complex(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine csytrf
      pure subroutine zsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         complex(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine zsytrf
   end interface lapack_sytrf

   interface lapack_sptrf
      pure subroutine ssptrf(uplo, n, ap, ipiv, info)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine ssptrf
      pure subroutine dsptrf(uplo, n, ap, ipiv, info)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine dsptrf
      pure subroutine csptrf(uplo, n, ap, ipiv, info)
         import :: sp
         complex(sp), intent(inout) :: ap(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine csptrf
      pure subroutine zsptrf(uplo, n, ap, ipiv, info)
         import :: dp
         complex(dp), intent(inout) :: ap(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine zsptrf
   end interface lapack_sptrf

   interface lapack_potrf
      pure subroutine spotrf(uplo, n, a, lda, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine spotrf
      pure subroutine dpotrf(uplo, n, a, lda, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine dpotrf
      pure subroutine cpotrf(uplo, n, a, lda, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine cpotrf
      pure subroutine zpotrf(uplo, n, a, lda, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
      end subroutine zpotrf
   end interface lapack_potrf

   interface lapack_pptrf
      pure subroutine spptrf(uplo, n, ap, info)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine spptrf
      pure subroutine dpptrf(uplo, n, ap, info)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine dpptrf
      pure subroutine cpptrf(uplo, n, ap, info)
         import :: sp
         complex(sp), intent(inout) :: ap(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine cpptrf
      pure subroutine zpptrf(uplo, n, ap, info)
         import :: dp
         complex(dp), intent(inout) :: ap(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine zpptrf
   end interface lapack_pptrf


contains


subroutine mctc_sgetrf(env, amat, ipiv)
   character(len=*), parameter :: source = 'mctc_lapack_getrf'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(out) :: ipiv(:)
   integer :: info, m, n, lda
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call lapack_getrf(m, n, amat, lda, ipiv, info)
   if (info /= 0) then
      call env%error("Factorisation of matrix failed", source)
   end if
end subroutine mctc_sgetrf


subroutine mctc_dgetrf(env, amat, ipiv)
   character(len=*), parameter :: source = 'mctc_lapack_getrf'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(out) :: ipiv(:)
   integer :: info, m, n, lda
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call lapack_getrf(m, n, amat, lda, ipiv, info)
   if (info /= 0) then
      call env%error("Factorisation of matrix failed", source)
   end if
end subroutine mctc_dgetrf


subroutine mctc_ssytrf(env, amat, ipiv, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sytrf'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(out) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, lda, lwork, stat_alloc, stat_dealloc
   real(sp), allocatable :: work(:)
   real(sp) :: test(1)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   lwork = -1
   call lapack_sytrf(ula, n, amat, lda, ipiv, test, lwork, info)
   if (info == 0) then
      lwork = nint(test(1))
      if (stat_alloc==0) then
         allocate(work(lwork), stat=stat_alloc)
      end if
      if (stat_alloc==0) then
         call lapack_sytrf(ula, n, amat, lda, ipiv, work, lwork, info)
      else
         info = -1000
      end if
      deallocate(work, stat=stat_dealloc)
   end if
   if (info /= 0) then
      call env%error("Factorisation of matrix failed", source)
   end if
end subroutine mctc_ssytrf


subroutine mctc_dsytrf(env, amat, ipiv, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sytrf'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(out) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, lda, lwork, stat_alloc, stat_dealloc
   real(dp), allocatable :: work(:)
   real(dp) :: test(1)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   lwork = -1
   call lapack_sytrf(ula, n, amat, lda, ipiv, test, lwork, info)
   if (info == 0) then
      lwork = nint(test(1))
      if (stat_alloc==0) then
         allocate(work(lwork), stat=stat_alloc)
      end if
      if (stat_alloc==0) then
         call lapack_sytrf(ula, n, amat, lda, ipiv, work, lwork, info)
      else
         info = -1000
      end if
      deallocate(work, stat=stat_dealloc)
   end if
   if (info /= 0) then
      call env%error("Factorisation of matrix failed", source)
   end if
end subroutine mctc_dsytrf


subroutine mctc_ssptrf(env, amat, ipiv, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sptrf'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:)
   integer, intent(out) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, stat_alloc, stat_dealloc, nn
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
   if (stat_alloc==0) then
      call lapack_sptrf(ula, n, amat, ipiv, info)
   else
      info = -1000
   end if
   if (info /= 0) then
      call env%error("Factorisation of matrix failed", source)
   end if
end subroutine mctc_ssptrf


subroutine mctc_dsptrf(env, amat, ipiv, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sptrf'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:)
   integer, intent(out) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: info, n, stat_alloc, stat_dealloc, nn
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
   if (stat_alloc==0) then
      call lapack_sptrf(ula, n, amat, ipiv, info)
   else
      info = -1000
   end if
   if (info /= 0) then
      call env%error("Factorisation of matrix failed", source)
   end if
end subroutine mctc_dsptrf


subroutine mctc_spotrf(env, amat, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_potrf'
   character(len=1), intent(in), optional :: uplo
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   character(len=1) :: ula
   integer :: info, n, lda
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_potrf(ula, n, amat, lda, info)
   if (info /= 0) then
      call env%error("Factorisation of matrix failed", source)
   end if
end subroutine mctc_spotrf


subroutine mctc_dpotrf(env, amat, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_potrf'
   character(len=1), intent(in), optional :: uplo
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   character(len=1) :: ula
   integer :: info, n, lda
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_potrf(ula, n, amat, lda, info)
   if (info /= 0) then
      call env%error("Factorisation of matrix failed", source)
   end if
end subroutine mctc_dpotrf


subroutine mctc_spptrf(env, amat, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_pptrf'
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
   call lapack_pptrf(ula, n, amat, info)
   if (info /= 0) then
      call env%error("Factorisation of matrix failed", source)
   end if
end subroutine mctc_spptrf


subroutine mctc_dpptrf(env, amat, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_pptrf'
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
   call lapack_pptrf(ula, n, amat, info)
   if (info /= 0) then
      call env%error("Factorisation of matrix failed", source)
   end if
end subroutine mctc_dpptrf


end module xtb_mctc_lapack_trf
