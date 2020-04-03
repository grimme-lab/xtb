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

!> Interfaces to LAPACK
module xtb_mctc_lapack
   use xtb_mctc_accuracy, only : sp, dp
   implicit none
   private

   public :: syev, spev, syevd, spevd, syevx, spevx, sygvd
   public :: sysv, sytrf, sytri


   interface syev
      pure subroutine ssyev(jobz,uplo,n,a,lda,w,work,lwork,info)
         import :: sp
         integer, intent(in)    :: lda
         real(sp),intent(inout) :: a(lda,*)
         real(sp),intent(out)   :: w(*)
         character,intent(in)   :: jobz
         character,intent(in)   :: uplo
         integer, intent(out)   :: info
         integer, intent(in)    :: n
         real(sp),intent(inout) :: work(*)
         integer, intent(in)    :: lwork
      end subroutine ssyev
      pure subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
         import :: dp
         integer, intent(in)    :: lda
         real(dp),intent(inout) :: a(lda,*)
         real(dp),intent(out)   :: w(*)
         character,intent(in)   :: jobz
         character,intent(in)   :: uplo
         integer, intent(out)   :: info
         integer, intent(in)    :: n
         real(dp),intent(inout) :: work(*)
         integer, intent(in)    :: lwork
      end subroutine dsyev
   end interface syev

   interface spev
      pure subroutine sspev(jobz,uplo,n,ap,w,z,ldz,work,info)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         real(sp), intent(out) :: w(*)
         character, intent(in) :: uplo
         real(sp), intent(out) :: z(ldz,*)
         integer, intent(out) :: info
         character, intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(sp), intent(in) :: work(*)
      end subroutine sspev
      pure subroutine dspev(jobz,uplo,n,ap,w,z,ldz,work,info)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         real(dp), intent(out) :: w(*)
         character, intent(in) :: uplo
         real(dp), intent(out) :: z(ldz,*)
         integer, intent(out) :: info
         character, intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(dp), intent(in) :: work(*)
      end subroutine dspev
   end interface spev

   interface syevd
      pure subroutine ssyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info)
         import :: sp
         real(sp), intent(inout) :: a(lda,*)
         real(sp), intent(out) :: w(*)
         character, intent(in) :: jobz
         character, intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine ssyevd
      pure subroutine dsyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info)
         import :: dp
         real(dp), intent(inout) :: a(lda,*)
         real(dp), intent(out) :: w(*)
         character, intent(in) :: jobz
         character, intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine dsyevd
   end interface syevd

   interface spevd
      pure subroutine sspevd(jobz,uplo,n,ap,w,z,ldz,work,lwork,iwork,liwork,&
            &                                                             info)
         import sp
         real(sp), intent(inout) :: ap(*)
         real(sp), intent(out) :: w(*)
         character, intent(in) :: uplo
         real(sp), intent(out) :: z(ldz,*)
         integer, intent(out) :: info
         character, intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine sspevd
      pure subroutine dspevd(jobz,uplo,n,ap,w,z,ldz,work,lwork,iwork,liwork,&
            &                                                             info)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         real(dp), intent(out) :: w(*)
         character, intent(in) :: uplo
         real(dp), intent(out) :: z(ldz,*)
         integer, intent(out) :: info
         character, intent(in) :: jobz
         integer, intent(in) :: n
         integer, intent(in) :: ldz
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine dspevd
   end interface spevd

   interface syevx
      pure subroutine ssyevx(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,&
            &                                z,ldz,work,lwork,iwork,ifail,info)
         import sp
         real(sp),intent(inout) :: a(lda,*)
         real(sp),intent(out)   :: w(*)
         character,intent(in)   :: uplo
         real(sp),intent(out)   :: z(ldz,*)
         real(sp),intent(in)    :: vl
         real(sp),intent(in)    :: vu
         integer, intent(in)    :: il
         integer, intent(in)    :: iu
         integer, intent(out)   :: m
         integer, intent(out)   :: ifail(*)
         real(sp),intent(in)    :: abstol
         integer, intent(out)   :: info
         character,intent(in)   :: jobz
         character,intent(in)   :: range
         integer, intent(in)    :: n
         integer, intent(in)    :: lda
         integer, intent(in)    :: ldz
         real(sp),intent(inout) :: work(*)
         integer, intent(in)    :: lwork
         integer, intent(in)    :: iwork(*)
      end subroutine ssyevx
      pure subroutine dsyevx(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,&
            &                                z,ldz,work,lwork,iwork,ifail,info)
         import :: dp
         real(dp),intent(inout) :: a(lda,*)
         real(dp),intent(out)   :: w(*)
         character,intent(in)   :: uplo
         real(dp),intent(out)   :: z(ldz,*)
         real(dp),intent(in)    :: vl
         real(dp),intent(in)    :: vu
         integer, intent(in)    :: il
         integer, intent(in)    :: iu
         integer, intent(out)   :: m
         integer, intent(out)   :: ifail(*)
         real(dp),intent(in)    :: abstol
         integer, intent(out)   :: info
         character,intent(in)   :: jobz
         character,intent(in)   :: range
         integer, intent(in)    :: n
         integer, intent(in)    :: lda
         integer, intent(in)    :: ldz
         real(dp),intent(inout) :: work(*)
         integer, intent(in)    :: lwork
         integer, intent(in)    :: iwork(*)
      end subroutine dsyevx
   end interface syevx

   interface spevx
      pure subroutine sspevx(jobz,range,uplo,n,ap,vl,vu,il,iu,abstol,m,w,z, &
            &                                        ldz,work,iwork,ifail,info)
         import :: sp
         real(sp),intent(inout) :: ap(*)
         real(sp),intent(out)   :: w(*)
         character,intent(in)   :: uplo
         real(sp),intent(out)   :: z(ldz,*)
         real(sp),intent(in)    :: vl
         real(sp),intent(in)    :: vu
         integer, intent(in)    :: il
         integer, intent(in)    :: iu
         integer, intent(out)   :: m
         integer, intent(out)   :: ifail(*)
         real(sp),intent(in)    :: abstol
         integer, intent(out)   :: info
         character,intent(in)   :: jobz
         character,intent(in)   :: range
         integer, intent(in)    :: n
         integer, intent(in)    :: ldz
         real(sp),intent(in)    :: work(*)
         integer, intent(in)    :: iwork(*)
      end subroutine sspevx
      pure subroutine dspevx(jobz,range,uplo,n,ap,vl,vu,il,iu,abstol,m,w,z, &
            &                                        ldz,work,iwork,ifail,info)
         import :: dp
         real(dp),intent(inout) :: ap(*)
         real(dp),intent(out)   :: w(*)
         character,intent(in)   :: uplo
         real(dp),intent(out)   :: z(ldz,*)
         real(dp),intent(in)    :: vl
         real(dp),intent(in)    :: vu
         integer, intent(in)    :: il
         integer, intent(in)    :: iu
         integer, intent(out)   :: m
         integer, intent(out)   :: ifail(*)
         real(dp),intent(in)    :: abstol
         integer, intent(out)   :: info
         character,intent(in)   :: jobz
         character,intent(in)   :: range
         integer, intent(in)    :: n
         integer, intent(in)    :: ldz
         real(dp),intent(in)    :: work(*)
         integer, intent(in)    :: iwork(*)
      end subroutine dspevx
   end interface spevx

   interface sygvd
      pure subroutine ssygvd(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,&
            &                      iwork,liwork,info)
         import :: sp
         real(sp),intent(inout) :: a(lda,*)
         real(sp),intent(inout) :: b(ldb,*)
         real(sp),intent(out)   :: w(*)
         integer, intent(in)    :: itype
         character,intent(in)   :: jobz
         character,intent(in)   :: uplo
         integer, intent(out)   :: info
         integer, intent(in)    :: n
         integer, intent(in)    :: lda
         integer, intent(in)    :: ldb
         real(sp),intent(inout) :: work(*)
         integer, intent(in)    :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in)    :: liwork
      end subroutine ssygvd
      pure subroutine dsygvd(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,&
            &                      iwork,liwork,info)
         import :: dp
         real(dp),intent(inout) :: a(lda,*)
         real(dp),intent(inout) :: b(ldb,*)
         real(dp),intent(out)   :: w(*)
         integer, intent(in)    :: itype
         character,intent(in)   :: jobz
         character,intent(in)   :: uplo
         integer, intent(out)   :: info
         integer, intent(in)    :: n
         integer, intent(in)    :: lda
         integer, intent(in)    :: ldb
         real(dp),intent(inout) :: work(*)
         integer, intent(in)    :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in)    :: liwork
      end subroutine dsygvd
   end interface sygvd


   interface sysv
      pure subroutine ssysv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
         import :: sp
         real(sp), intent(inout) :: a(lda,*)
         real(sp), intent(inout) :: b(ldb,*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine ssysv
      pure subroutine dsysv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
         import :: dp
         real(dp), intent(inout) :: a(lda,*)
         real(dp), intent(inout) :: b(ldb,*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine dsysv
   end interface sysv

   interface sytrf
      pure subroutine ssytrf(uplo,n,a,lda,ipiv,work,lwork,info)
         import :: sp
         real(sp), intent(inout) :: a(lda,*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine ssytrf
      pure subroutine dsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
         import dp
         real(dp), intent(inout) :: a(lda,*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine dsytrf
   end interface sytrf

   interface sytri
      pure subroutine ssytri(uplo,n,a,lda,ipiv,work,info)
         import :: sp
         real(sp), intent(inout) :: a(lda,*)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(in) :: work(*)
      end subroutine ssytri
      pure subroutine dsytri(uplo,n,a,lda,ipiv,work,info)
         import dp
         real(dp), intent(inout) :: a(lda,*)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(in) :: work(*)
      end subroutine dsytri
   end interface sytri


end module xtb_mctc_lapack
