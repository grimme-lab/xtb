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
module xtb_mctc_lapack_geneigval
   use xtb_mctc_accuracy, only : sp, dp
   implicit none
   private

   public :: lapack_sygv, lapack_sygvd, lapack_sygvx, lapack_spgv, lapack_spgvd
   public :: lapack_spgvx
   public :: lapack_hegv, lapack_hegvd, lapack_hegvx, lapack_hpgv, lapack_hpgvd
   public :: lapack_hpgvx


   interface lapack_sygv
      pure subroutine ssygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
            & info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine ssygv
      pure subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
            & info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         real(dp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine dsygv
   end interface lapack_sygv

   interface lapack_hegv
      pure subroutine chegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
            & rwork, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         complex(sp), intent(inout) :: b(ldb, *)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         complex(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(sp), intent(in) :: rwork(*)
      end subroutine chegv
      pure subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
            & rwork, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         complex(dp), intent(inout) :: b(ldb, *)
         real(dp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         complex(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(dp), intent(in) :: rwork(*)
      end subroutine zhegv
   end interface lapack_hegv

   interface lapack_sygvd
      pure subroutine ssygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
            & iwork, liwork, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine ssygvd
      pure subroutine dsygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
            & iwork, liwork, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         real(dp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine dsygvd
   end interface lapack_sygvd

   interface lapack_hegvd
      pure subroutine chegvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
            & rwork, lrwork, iwork, liwork, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         complex(sp), intent(inout) :: b(ldb, *)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         complex(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(sp), intent(inout) :: rwork(*)
         integer, intent(in) :: lrwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine chegvd
      pure subroutine zhegvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
            & rwork, lrwork, iwork, liwork, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         complex(dp), intent(inout) :: b(ldb, *)
         real(dp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         complex(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(dp), intent(inout) :: rwork(*)
         integer, intent(in) :: lrwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine zhegvd
   end interface lapack_hegvd

   interface lapack_sygvx
      pure subroutine ssygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, &
            & il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
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
         integer, intent(in) :: ldb
         integer, intent(in) :: ldz
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(in) :: iwork(*)
      end subroutine ssygvx
      pure subroutine dsygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, &
            & il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         real(dp), intent(out) :: w(*)
         integer, intent(in) :: itype
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
         integer, intent(in) :: ldb
         integer, intent(in) :: ldz
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(in) :: iwork(*)
      end subroutine dsygvx
   end interface lapack_sygvx

   interface lapack_hegvx
      pure subroutine chegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, &
            & il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         complex(sp), intent(inout) :: b(ldb, *)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
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
         integer, intent(in) :: ldb
         integer, intent(in) :: ldz
         complex(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(sp), intent(in) :: rwork(*)
         integer, intent(in) :: iwork(*)
      end subroutine chegvx
      pure subroutine zhegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, &
            & il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         complex(dp), intent(inout) :: b(ldb, *)
         real(dp), intent(out) :: w(*)
         integer, intent(in) :: itype
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
         integer, intent(in) :: ldb
         integer, intent(in) :: ldz
         complex(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         real(dp), intent(in) :: rwork(*)
         integer, intent(in) :: iwork(*)
      end subroutine zhegvx
   end interface lapack_hegvx

   interface lapack_spgv
      pure subroutine sspgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         real(sp), intent(inout) :: bp(*)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         real(sp), intent(out) :: z(ldz, *)
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(sp), intent(in) :: work(*)
      end subroutine sspgv
      pure subroutine dspgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         real(dp), intent(inout) :: bp(*)
         real(dp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         real(dp), intent(out) :: z(ldz, *)
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(dp), intent(in) :: work(*)
      end subroutine dspgv
   end interface lapack_spgv

   interface lapack_hpgv
      pure subroutine chpgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, rwork, &
            & info)
         import :: sp
         complex(sp), intent(inout) :: ap(*)
         complex(sp), intent(inout) :: bp(*)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         complex(sp), intent(out) :: z(ldz, *)
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         complex(sp), intent(in) :: work(*)
         real(sp), intent(in) :: rwork(*)
      end subroutine chpgv
      pure subroutine zhpgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, rwork, &
            & info)
         import :: dp
         complex(dp), intent(inout) :: ap(*)
         complex(dp), intent(inout) :: bp(*)
         real(dp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         complex(dp), intent(out) :: z(ldz, *)
         integer, intent(out) :: info
         character(len=1), intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         complex(dp), intent(in) :: work(*)
         real(dp), intent(in) :: rwork(*)
      end subroutine zhpgv
   end interface lapack_hpgv

   interface lapack_spgvd
      pure subroutine sspgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, &
            & iwork, liwork, info)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         real(sp), intent(inout) :: bp(*)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
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
      end subroutine sspgvd
      pure subroutine dspgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, &
            & iwork, liwork, info)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         real(dp), intent(inout) :: bp(*)
         real(dp), intent(out) :: w(*)
         integer, intent(in) :: itype
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
      end subroutine dspgvd
   end interface lapack_spgvd

   interface lapack_hpgvd
      pure subroutine chpgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, &
            & rwork, lrwork, iwork, liwork, info)
         import :: sp
         complex(sp), intent(inout) :: ap(*)
         complex(sp), intent(inout) :: bp(*)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
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
      end subroutine chpgvd
      pure subroutine zhpgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, &
            & rwork, lrwork, iwork, liwork, info)
         import :: dp
         complex(dp), intent(inout) :: ap(*)
         complex(dp), intent(inout) :: bp(*)
         real(dp), intent(out) :: w(*)
         integer, intent(in) :: itype
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
      end subroutine zhpgvd
   end interface lapack_hpgvd

   interface lapack_spgvx
      pure subroutine sspgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, &
            & abstol, m, w, z, ldz, work, iwork, ifail, info)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         real(sp), intent(inout) :: bp(*)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
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
      end subroutine sspgvx
      pure subroutine dspgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, &
            & abstol, m, w, z, ldz, work, iwork, ifail, info)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         real(dp), intent(inout) :: bp(*)
         real(dp), intent(out) :: w(*)
         integer, intent(in) :: itype
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
      end subroutine dspgvx
   end interface lapack_spgvx

   interface lapack_hpgvx
      pure subroutine chpgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, &
            & abstol, m, w, z, ldz, work, rwork, iwork, ifail, info)
         import :: sp
         complex(sp), intent(inout) :: ap(*)
         complex(sp), intent(inout) :: bp(*)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
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
      end subroutine chpgvx
      pure subroutine zhpgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, &
            & abstol, m, w, z, ldz, work, rwork, iwork, ifail, info)
         import :: dp
         complex(dp), intent(inout) :: ap(*)
         complex(dp), intent(inout) :: bp(*)
         real(dp), intent(out) :: w(*)
         integer, intent(in) :: itype
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
      end subroutine zhpgvx
   end interface lapack_hpgvx


contains


end module xtb_mctc_lapack_geneigval
