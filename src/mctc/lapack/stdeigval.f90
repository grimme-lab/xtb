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

!> LAPACK eigenproblem solvers
module xtb_mctc_lapack_stdeigval
   use xtb_mctc_accuracy, only : sp, dp, i8
   implicit none
   private

   public :: mctc_syev, mctc_syevd, mctc_syevx, mctc_spevx
   public :: lapack_syevd, lapack_sytrd, lapack_stedc, lapack_ormtr


   interface mctc_syev
      module procedure :: mctc_ssyev
      module procedure :: mctc_dsyev
   end interface mctc_syev

   interface mctc_syevd
      module procedure :: mctc_ssyevd
      module procedure :: mctc_dsyevd
   end interface mctc_syevd

   interface mctc_syevx
      module procedure :: mctc_ssyevx
      module procedure :: mctc_dsyevx
   end interface mctc_syevx

   interface mctc_spevx
      module procedure :: mctc_sspevx
      module procedure :: mctc_dspevx
   end interface mctc_spevx

   interface lapack_sytrd
      pure subroutine ssytrd(uplo, n, a, lda, d, e, tau, work, lwork, info)
         import :: sp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: d(*)
         real(sp), intent(out) :: e(*)
         real(sp), intent(out) :: tau(*)
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(out) :: info
      end subroutine ssytrd
      pure subroutine dsytrd(uplo, n, a, lda, d, e, tau, work, lwork, info)
         import :: dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: d(*)
         real(dp), intent(out) :: e(*)
         real(dp), intent(out) :: tau(*)
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(out) :: info
      end subroutine dsytrd
   end interface lapack_sytrd

   interface lapack_stedc
      pure subroutine sstedc(compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info)
         import :: sp
         character(len=1), intent(in) :: compz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(sp), intent(inout) :: d(*)
         real(sp), intent(inout) :: e(*)
         real(sp), intent(inout) :: z(ldz, *)
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
         integer, intent(out) :: info
      end subroutine sstedc
      pure subroutine dstedc(compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info)
         import :: dp
         character(len=1), intent(in) :: compz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(dp), intent(inout) :: d(*)
         real(dp), intent(inout) :: e(*)
         real(dp), intent(inout) :: z(ldz, *)
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
         integer, intent(out) :: info
      end subroutine dstedc
   end interface lapack_stedc

   interface lapack_ormtr
      pure subroutine sormtr(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info)
         import :: sp
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldc
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: tau(*)
         real(sp), intent(inout) :: c(ldc, *)
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(out) :: info
      end subroutine sormtr
      pure subroutine dormtr(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info)
         import :: dp
         character(len=1), intent(in) :: side
         character(len=1), intent(in) :: uplo
         character(len=1), intent(in) :: trans
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldc
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: tau(*)
         real(dp), intent(inout) :: c(ldc, *)
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(out) :: info
      end subroutine dormtr
   end interface lapack_ormtr


   interface lapack_syev
      pure subroutine ssyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine ssyev
      pure subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine dsyev
   end interface lapack_syev

   interface lapack_heev
      pure subroutine cheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         complex(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(sp), intent(in) :: rwork(*)
      end subroutine cheev
      pure subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         complex(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(dp), intent(in) :: rwork(*)
      end subroutine zheev
   end interface lapack_heev

   interface lapack_syevd
      pure subroutine ssyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, &
            & liwork, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine ssyevd
      pure subroutine dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, &
            & liwork, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine dsyevd
   end interface lapack_syevd

   interface lapack_heevd
      pure subroutine cheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, &
            & lrwork, iwork, liwork, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         complex(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(sp), intent(inout) :: rwork(*)
         integer, intent(in) :: lrwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine cheevd
      pure subroutine zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, &
            & lrwork, iwork, liwork, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         complex(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(dp), intent(inout) :: rwork(*)
         integer, intent(in) :: lrwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine zheevd
   end interface lapack_heevd

   interface lapack_syevx
      pure subroutine ssyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
            & abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(out) :: z(ldz, *)
         real(sp), intent(in) :: vl
         real(sp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: ifail(*)
         real(sp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldz
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(out) :: iwork(*)
      end subroutine ssyevx
      pure subroutine dsyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
            & abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(out) :: z(ldz, *)
         real(dp), intent(in) :: vl
         real(dp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: ifail(*)
         real(dp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldz
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(out) :: iwork(*)
      end subroutine dsyevx
   end interface lapack_syevx

   interface lapack_heevx
      pure subroutine cheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
            & abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         complex(sp), intent(out) :: z(ldz, *)
         real(sp), intent(in) :: vl
         real(sp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: ifail(*)
         real(sp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldz
         complex(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(sp), intent(in) :: rwork(*)
         integer, intent(in) :: iwork(*)
      end subroutine cheevx
      pure subroutine zheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
            & abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         complex(dp), intent(out) :: z(ldz, *)
         real(dp), intent(in) :: vl
         real(dp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: ifail(*)
         real(dp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldz
         complex(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(dp), intent(in) :: rwork(*)
         integer, intent(in) :: iwork(*)
      end subroutine zheevx
   end interface lapack_heevx

   interface lapack_syevr
      pure subroutine ssyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
            & abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(out) :: z(ldz, *)
         real(sp), intent(in) :: vl
         real(sp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: isuppz(*)
         real(sp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldz
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine ssyevr
      pure subroutine dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
            & abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(out) :: z(ldz, *)
         real(dp), intent(in) :: vl
         real(dp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: isuppz(*)
         real(dp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldz
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine dsyevr
   end interface lapack_syevr

   interface lapack_heevr
      pure subroutine cheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
            & abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, &
            & liwork, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         complex(sp), intent(out) :: z(ldz, *)
         real(sp), intent(in) :: vl
         real(sp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: isuppz(*)
         real(sp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldz
         complex(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(sp), intent(inout) :: rwork(*)
         integer, intent(in) :: lrwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine cheevr
      pure subroutine zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
            & abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, &
            & liwork, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         complex(dp), intent(out) :: z(ldz, *)
         real(dp), intent(in) :: vl
         real(dp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: isuppz(*)
         real(dp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldz
         complex(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(dp), intent(inout) :: rwork(*)
         integer, intent(in) :: lrwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine zheevr
   end interface lapack_heevr

   interface lapack_spev
      pure subroutine sspev(jobz, uplo, n, ap, w, z, ldz, work, info)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(out) :: z(ldz, *)
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(sp), intent(in) :: work(*)
      end subroutine sspev
      pure subroutine dspev(jobz, uplo, n, ap, w, z, ldz, work, info)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(out) :: z(ldz, *)
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(dp), intent(in) :: work(*)
      end subroutine dspev
   end interface lapack_spev

   interface lapack_hpev
      pure subroutine chpev(jobz, uplo, n, ap, w, z, ldz, work, rwork, info)
         import :: sp
         complex(sp), intent(inout) :: ap(*)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         complex(sp), intent(out) :: z(ldz, *)
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         complex(sp), intent(in) :: work(*)
         real(sp), intent(in) :: rwork(*)
      end subroutine chpev
      pure subroutine zhpev(jobz, uplo, n, ap, w, z, ldz, work, rwork, info)
         import :: dp
         complex(dp), intent(inout) :: ap(*)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         complex(dp), intent(out) :: z(ldz, *)
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         complex(dp), intent(in) :: work(*)
         real(dp), intent(in) :: rwork(*)
      end subroutine zhpev
   end interface lapack_hpev

   interface lapack_spevd
      pure subroutine sspevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, &
            & liwork, info)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(out) :: z(ldz, *)
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine sspevd
      pure subroutine dspevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, &
            & liwork, info)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(out) :: z(ldz, *)
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine dspevd
   end interface lapack_spevd

   interface lapack_hpevd
      pure subroutine chpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, &
            & lrwork, iwork, liwork, info)
         import :: sp
         complex(sp), intent(inout) :: ap(*)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         complex(sp), intent(out) :: z(ldz, *)
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         complex(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(sp), intent(inout) :: rwork(*)
         integer, intent(in) :: lrwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine chpevd
      pure subroutine zhpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, &
            & lrwork, iwork, liwork, info)
         import :: dp
         complex(dp), intent(inout) :: ap(*)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         complex(dp), intent(out) :: z(ldz, *)
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         complex(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(dp), intent(inout) :: rwork(*)
         integer, intent(in) :: lrwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine zhpevd
   end interface lapack_hpevd

   interface lapack_spevx
      pure subroutine sspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, &
            & m, w, z, ldz, work, iwork, ifail, info)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(out) :: z(ldz, *)
         real(sp), intent(in) :: vl
         real(sp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: ifail(*)
         real(sp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(sp), intent(in) :: work(*)
         integer, intent(in) :: iwork(*)
      end subroutine sspevx
      pure subroutine dspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, &
            & m, w, z, ldz, work, iwork, ifail, info)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(out) :: z(ldz, *)
         real(dp), intent(in) :: vl
         real(dp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: ifail(*)
         real(dp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(dp), intent(in) :: work(*)
         integer, intent(in) :: iwork(*)
      end subroutine dspevx
   end interface lapack_spevx

   interface lapack_hpevx
      pure subroutine chpevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, &
            & m, w, z, ldz, work, rwork, iwork, ifail, info)
         import :: sp
         complex(sp), intent(inout) :: ap(*)
         real(sp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         complex(sp), intent(out) :: z(ldz, *)
         real(sp), intent(in) :: vl
         real(sp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: ifail(*)
         real(sp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         complex(sp), intent(in) :: work(*)
         real(sp), intent(in) :: rwork(*)
         integer, intent(in) :: iwork(*)
      end subroutine chpevx
      pure subroutine zhpevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, &
            & m, w, z, ldz, work, rwork, iwork, ifail, info)
         import :: dp
         complex(dp), intent(inout) :: ap(*)
         real(dp), intent(out) :: w(*)
         character(len=1), intent(in) :: uplo
         complex(dp), intent(out) :: z(ldz, *)
         real(dp), intent(in) :: vl
         real(dp), intent(in) :: vu
         integer, intent(in) :: il
         integer, intent(in) :: iu
         integer, intent(out) :: m
         integer, intent(out) :: ifail(*)
         real(dp), intent(in) :: abstol
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: range
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         complex(dp), intent(in) :: work(*)
         real(dp), intent(in) :: rwork(*)
         integer, intent(in) :: iwork(*)
      end subroutine zhpevx
   end interface lapack_hpevx


contains


pure subroutine mctc_ssyev(amat, eval, jobz, uplo, info)
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(out) :: eval(:)
   character(len=1), intent(in), optional :: jobz, uplo
   integer, intent(out), optional :: info

   character(len=1) :: job, ula
   integer :: info_, lda, lwork, n, stat
   real(sp) :: query(1)
   real(sp), allocatable :: work(:)

   job = 'V'
   if (present(jobz)) job = jobz
   ula = 'U'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   info_ = 0
   if (size(amat, 1) /= n) then
      info_ = -5
   else if (size(eval) < n) then
      info_ = -6
   else if (n > 0) then
      call lapack_syev(job, ula, n, amat, lda, eval, query, -1, info_)
      if (info_ == 0) then
         lwork = max(1, nint(query(1)))
         allocate(work(lwork), stat=stat)
         if (stat == 0) then
            call lapack_syev(job, ula, n, amat, lda, eval, work, lwork, info_)
         else
            info_ = -1000
         end if
      end if
   end if
   if (present(info)) info = info_
end subroutine mctc_ssyev


pure subroutine mctc_dsyev(amat, eval, jobz, uplo, info)
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(out) :: eval(:)
   character(len=1), intent(in), optional :: jobz, uplo
   integer, intent(out), optional :: info

   character(len=1) :: job, ula
   integer :: info_, lda, lwork, n, stat
   real(dp) :: query(1)
   real(dp), allocatable :: work(:)

   job = 'V'
   if (present(jobz)) job = jobz
   ula = 'U'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   info_ = 0
   if (size(amat, 1) /= n) then
      info_ = -5
   else if (size(eval) < n) then
      info_ = -6
   else if (n > 0) then
      call lapack_syev(job, ula, n, amat, lda, eval, query, -1, info_)
      if (info_ == 0) then
         lwork = max(1, nint(query(1)))
         allocate(work(lwork), stat=stat)
         if (stat == 0) then
            call lapack_syev(job, ula, n, amat, lda, eval, work, lwork, info_)
         else
            info_ = -1000
         end if
      end if
   end if
   if (present(info)) info = info_
end subroutine mctc_dsyev


pure subroutine mctc_ssyevd(amat, eval, jobz, uplo, info)
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(out) :: eval(:)
   character(len=1), intent(in), optional :: jobz, uplo
   integer, intent(out), optional :: info

   character(len=1) :: job, ula
   integer :: info_, lda, lwork, liwork, n, stat, iquery(1)
   real(sp) :: query(1)
   real(sp), allocatable :: work(:)
   integer, allocatable :: iwork(:)

   job = 'V'
   if (present(jobz)) job = jobz
   ula = 'U'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   info_ = 0
   if (size(amat, 1) /= n) then
      info_ = -5
   else if (size(eval) < n) then
      info_ = -6
   else if (n > 0) then
      call lapack_syevd(job, ula, n, amat, lda, eval, query, -1, iquery, -1, info_)
      if (info_ == 0) then
         lwork = max(1, nint(query(1)))
         liwork = max(1, iquery(1))
         allocate(work(lwork), iwork(liwork), stat=stat)
         if (stat == 0) then
            call lapack_syevd(job, ula, n, amat, lda, eval, work, lwork, iwork, liwork, info_)
         else
            info_ = -1000
         end if
      end if
   end if
   if (present(info)) info = info_
end subroutine mctc_ssyevd


pure subroutine mctc_dsyevd(amat, eval, jobz, uplo, info)
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(out) :: eval(:)
   character(len=1), intent(in), optional :: jobz, uplo
   integer, intent(out), optional :: info

   character(len=1) :: job, ula
   integer :: info_, lda, lwork, liwork, n, stat, iquery(1)
   real(dp) :: query(1)
   real(dp), allocatable :: work(:)
   integer, allocatable :: iwork(:)

   job = 'V'
   if (present(jobz)) job = jobz
   ula = 'U'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   info_ = 0
   if (size(amat, 1) /= n) then
      info_ = -5
   else if (size(eval) < n) then
      info_ = -6
   else if (n > 0) then
      call lapack_syevd(job, ula, n, amat, lda, eval, query, -1, iquery, -1, info_)
      if (info_ == 0) then
         lwork = max(1, nint(query(1)))
         liwork = max(1, iquery(1))
         allocate(work(lwork), iwork(liwork), stat=stat)
         if (stat == 0) then
            call lapack_syevd(job, ula, n, amat, lda, eval, work, lwork, iwork, liwork, info_)
         else
            info_ = -1000
         end if
      end if
   end if
   if (present(info)) info = info_
end subroutine mctc_dsyevd


pure subroutine mctc_ssyevx(amat, eval, uplo, z, vl, vu, il, iu, abstol, m, info, n)
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(out) :: eval(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), intent(out), optional :: z(:, :)
   real(sp), intent(in), optional :: vl, vu, abstol
   integer, intent(in), optional :: il, iu, n
   integer, intent(out), optional :: m, info

   character(len=1) :: rng, ula
   integer :: nn, lda, ldz, ill, iuu, mm, lwork, info_, stat, idum1(1), idum2(1)
   real(sp) :: vll, vuu, tol, query(1), zdum(1, 1)
   real(sp), allocatable :: work(:)
   integer, allocatable :: iwork(:), ifail(:)

   ula = 'U'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   nn = size(amat, 2)
   if (present(n)) nn = n
   vll = 0.0_sp
   if (present(vl)) vll = vl
   vuu = 0.0_sp
   if (present(vu)) vuu = vu
   ill = 1
   if (present(il)) ill = il
   iuu = nn
   if (present(iu)) iuu = iu
   tol = 0.0_sp
   if (present(abstol)) tol = abstol
   rng = merge('I', merge('V', 'A', present(vl) .or. present(vu)), present(il) .or. present(iu))

   mm = 0
   info_ = 0
   if ((present(il) .or. present(iu)) .and. (present(vl) .or. present(vu))) then
      info_ = -2
   else if (nn < 0) then
      info_ = -4
   else if (size(amat, 1) < nn .or. size(amat, 2) < nn &
         & .or. (.not. present(n) .and. size(amat, 1) /= nn)) then
      info_ = -6
   else if (rng == 'I' .and. (ill < 1 .or. ill > max(1, nn))) then
      info_ = -9
   else if (rng == 'I' .and. (iuu < min(nn, ill) .or. iuu > nn)) then
      info_ = -10
   else if (size(eval) < nn) then
      info_ = -13
   else if (present(z)) then
      if (size(z, 1) < nn) then
         info_ = -15
      else if (nn > 0 .and. size(z, 2) < merge(iuu - ill + 1, nn, rng == 'I')) then
         info_ = -14
      end if
   end if

   if (info_ == 0 .and. nn > 0) then
      if (present(z)) then
         ldz = max(1, size(z, 1))
         call lapack_syevx('V', rng, ula, nn, amat, lda, vll, vuu, ill, iuu, tol, &
            & mm, eval, z, ldz, query, -1, idum1, idum2, info_)
      else
         call lapack_syevx('N', rng, ula, nn, amat, lda, vll, vuu, ill, iuu, tol, &
            & mm, eval, zdum, 1, query, -1, idum1, idum2, info_)
      end if
      if (info_ == 0) then
         lwork = max(1, 8*nn, nint(query(1)))
         allocate(work(lwork), iwork(max(1, 5*nn)), ifail(max(1, nn)), stat=stat)
         if (stat /= 0) then
            info_ = -1000
         else if (present(z)) then
            call lapack_syevx('V', rng, ula, nn, amat, lda, vll, vuu, ill, iuu, tol, &
               & mm, eval, z, ldz, work, lwork, iwork, ifail, info_)
         else
            call lapack_syevx('N', rng, ula, nn, amat, lda, vll, vuu, ill, iuu, tol, &
               & mm, eval, zdum, 1, work, lwork, iwork, ifail, info_)
         end if
      end if
   end if
   if (present(m)) m = mm
   if (present(info)) info = info_
end subroutine mctc_ssyevx


pure subroutine mctc_dsyevx(amat, eval, uplo, z, vl, vu, il, iu, abstol, m, info, n)
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(out) :: eval(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), intent(out), optional :: z(:, :)
   real(dp), intent(in), optional :: vl, vu, abstol
   integer, intent(in), optional :: il, iu, n
   integer, intent(out), optional :: m, info

   character(len=1) :: rng, ula
   integer :: nn, lda, ldz, ill, iuu, mm, lwork, info_, stat, idum1(1), idum2(1)
   real(dp) :: vll, vuu, tol, query(1), zdum(1, 1)
   real(dp), allocatable :: work(:)
   integer, allocatable :: iwork(:), ifail(:)

   ula = 'U'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   nn = size(amat, 2)
   if (present(n)) nn = n
   vll = 0.0_dp
   if (present(vl)) vll = vl
   vuu = 0.0_dp
   if (present(vu)) vuu = vu
   ill = 1
   if (present(il)) ill = il
   iuu = nn
   if (present(iu)) iuu = iu
   tol = 0.0_dp
   if (present(abstol)) tol = abstol
   rng = merge('I', merge('V', 'A', present(vl) .or. present(vu)), present(il) .or. present(iu))

   mm = 0
   info_ = 0
   if ((present(il) .or. present(iu)) .and. (present(vl) .or. present(vu))) then
      info_ = -2
   else if (nn < 0) then
      info_ = -4
   else if (size(amat, 1) < nn .or. size(amat, 2) < nn &
         & .or. (.not. present(n) .and. size(amat, 1) /= nn)) then
      info_ = -6
   else if (rng == 'I' .and. (ill < 1 .or. ill > max(1, nn))) then
      info_ = -9
   else if (rng == 'I' .and. (iuu < min(nn, ill) .or. iuu > nn)) then
      info_ = -10
   else if (size(eval) < nn) then
      info_ = -13
   else if (present(z)) then
      if (size(z, 1) < nn) then
         info_ = -15
      else if (nn > 0 .and. size(z, 2) < merge(iuu - ill + 1, nn, rng == 'I')) then
         info_ = -14
      end if
   end if

   if (info_ == 0 .and. nn > 0) then
      if (present(z)) then
         ldz = max(1, size(z, 1))
         call lapack_syevx('V', rng, ula, nn, amat, lda, vll, vuu, ill, iuu, tol, &
            & mm, eval, z, ldz, query, -1, idum1, idum2, info_)
      else
         call lapack_syevx('N', rng, ula, nn, amat, lda, vll, vuu, ill, iuu, tol, &
            & mm, eval, zdum, 1, query, -1, idum1, idum2, info_)
      end if
      if (info_ == 0) then
         lwork = max(1, 8*nn, nint(query(1)))
         allocate(work(lwork), iwork(max(1, 5*nn)), ifail(max(1, nn)), stat=stat)
         if (stat /= 0) then
            info_ = -1000
         else if (present(z)) then
            call lapack_syevx('V', rng, ula, nn, amat, lda, vll, vuu, ill, iuu, tol, &
               & mm, eval, z, ldz, work, lwork, iwork, ifail, info_)
         else
            call lapack_syevx('N', rng, ula, nn, amat, lda, vll, vuu, ill, iuu, tol, &
               & mm, eval, zdum, 1, work, lwork, iwork, ifail, info_)
         end if
      end if
   end if
   if (present(m)) m = mm
   if (present(info)) info = info_
end subroutine mctc_dsyevx


pure subroutine mctc_sspevx(amat, eval, uplo, z, vl, vu, il, iu, abstol, m, info, n)
   real(sp), intent(inout) :: amat(:)
   real(sp), intent(out) :: eval(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), intent(out), optional :: z(:, :)
   real(sp), intent(in), optional :: vl, vu, abstol
   integer, intent(in), optional :: il, iu, n
   integer, intent(out), optional :: m, info

   character(len=1) :: rng, ula
   integer :: nn, ill, iuu, mm, info_, stat
   integer(i8) :: np
   real(sp) :: vll, vuu, tol, zdum(1, 1)
   real(sp), allocatable :: work(:)
   integer, allocatable :: iwork(:), ifail(:)

   ula = 'U'
   if (present(uplo)) ula = uplo
   if (present(n)) then
      nn = n
   else
      nn = nint((-1.0_dp + sqrt(1.0_dp + 8.0_dp*real(size(amat), dp)))*0.5_dp)
   end if
   vll = 0.0_sp
   if (present(vl)) vll = vl
   vuu = 0.0_sp
   if (present(vu)) vuu = vu
   ill = 1
   if (present(il)) ill = il
   iuu = nn
   if (present(iu)) iuu = iu
   tol = 0.0_sp
   if (present(abstol)) tol = abstol
   rng = merge('I', merge('V', 'A', present(vl) .or. present(vu)), present(il) .or. present(iu))

   mm = 0
   info_ = 0
   if ((present(il) .or. present(iu)) .and. (present(vl) .or. present(vu))) then
      info_ = -2
   else if (nn < 0) then
      info_ = -4
   else
      np = int(nn, i8) * int(nn + 1, i8) / 2_i8
      if (size(amat, kind=i8) < np &
            & .or. (.not. present(n) .and. size(amat, kind=i8) /= np)) then
         info_ = -5
      else if (rng == 'I' .and. (ill < 1 .or. ill > max(1, nn))) then
         info_ = -8
      else if (rng == 'I' .and. (iuu < min(nn, ill) .or. iuu > nn)) then
         info_ = -9
      else if (size(eval) < nn) then
         info_ = -12
      else if (present(z)) then
         if (size(z, 1) < nn) then
            info_ = -14
         else if (nn > 0 .and. size(z, 2) < merge(iuu - ill + 1, nn, rng == 'I')) then
            info_ = -13
         end if
      end if
   end if

   if (info_ == 0 .and. nn > 0) then
      allocate(work(max(1, 8*nn)), iwork(max(1, 5*nn)), ifail(max(1, nn)), stat=stat)
      if (stat /= 0) then
         info_ = -1000
      else if (present(z)) then
         call lapack_spevx('V', rng, ula, nn, amat, vll, vuu, ill, iuu, tol, &
            & mm, eval, z, max(1, size(z, 1)), work, iwork, ifail, info_)
      else
         call lapack_spevx('N', rng, ula, nn, amat, vll, vuu, ill, iuu, tol, &
            & mm, eval, zdum, 1, work, iwork, ifail, info_)
      end if
   end if
   if (present(m)) m = mm
   if (present(info)) info = info_
end subroutine mctc_sspevx


pure subroutine mctc_dspevx(amat, eval, uplo, z, vl, vu, il, iu, abstol, m, info, n)
   real(dp), intent(inout) :: amat(:)
   real(dp), intent(out) :: eval(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), intent(out), optional :: z(:, :)
   real(dp), intent(in), optional :: vl, vu, abstol
   integer, intent(in), optional :: il, iu, n
   integer, intent(out), optional :: m, info

   character(len=1) :: rng, ula
   integer :: nn, ill, iuu, mm, info_, stat
   integer(i8) :: np
   real(dp) :: vll, vuu, tol, zdum(1, 1)
   real(dp), allocatable :: work(:)
   integer, allocatable :: iwork(:), ifail(:)

   ula = 'U'
   if (present(uplo)) ula = uplo
   if (present(n)) then
      nn = n
   else
      nn = nint((-1.0_dp + sqrt(1.0_dp + 8.0_dp*real(size(amat), dp)))*0.5_dp)
   end if
   vll = 0.0_dp
   if (present(vl)) vll = vl
   vuu = 0.0_dp
   if (present(vu)) vuu = vu
   ill = 1
   if (present(il)) ill = il
   iuu = nn
   if (present(iu)) iuu = iu
   tol = 0.0_dp
   if (present(abstol)) tol = abstol
   rng = merge('I', merge('V', 'A', present(vl) .or. present(vu)), present(il) .or. present(iu))

   mm = 0
   info_ = 0
   if ((present(il) .or. present(iu)) .and. (present(vl) .or. present(vu))) then
      info_ = -2
   else if (nn < 0) then
      info_ = -4
   else
      np = int(nn, i8) * int(nn + 1, i8) / 2_i8
      if (size(amat, kind=i8) < np &
            & .or. (.not. present(n) .and. size(amat, kind=i8) /= np)) then
         info_ = -5
      else if (rng == 'I' .and. (ill < 1 .or. ill > max(1, nn))) then
         info_ = -8
      else if (rng == 'I' .and. (iuu < min(nn, ill) .or. iuu > nn)) then
         info_ = -9
      else if (size(eval) < nn) then
         info_ = -12
      else if (present(z)) then
         if (size(z, 1) < nn) then
            info_ = -14
         else if (nn > 0 .and. size(z, 2) < merge(iuu - ill + 1, nn, rng == 'I')) then
            info_ = -13
         end if
      end if
   end if

   if (info_ == 0 .and. nn > 0) then
      allocate(work(max(1, 8*nn)), iwork(max(1, 5*nn)), ifail(max(1, nn)), stat=stat)
      if (stat /= 0) then
         info_ = -1000
      else if (present(z)) then
         call lapack_spevx('V', rng, ula, nn, amat, vll, vuu, ill, iuu, tol, &
            & mm, eval, z, max(1, size(z, 1)), work, iwork, ifail, info_)
      else
         call lapack_spevx('N', rng, ula, nn, amat, vll, vuu, ill, iuu, tol, &
            & mm, eval, zdum, 1, work, iwork, ifail, info_)
      end if
   end if
   if (present(m)) m = mm
   if (present(info)) info = info_
end subroutine mctc_dspevx


end module xtb_mctc_lapack_stdeigval
