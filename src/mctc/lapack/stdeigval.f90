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
   use xtb_mctc_accuracy, only : sp, dp
   implicit none
   private

   public :: lapack_syev, lapack_syevd, lapack_syevx, lapack_syevr, lapack_spev
   public :: lapack_spevd, lapack_spevx
   public :: lapack_heev, lapack_heevd, lapack_heevx, lapack_heevr, lapack_hpev
   public :: lapack_hpevd, lapack_hpevx


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
         integer, intent(in) :: iwork(*)
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
         integer, intent(in) :: iwork(*)
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


end module xtb_mctc_lapack_stdeigval
