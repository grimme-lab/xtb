! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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

module mctc_la
   implicit none

interface syev
   pure subroutine ssyev(jobz,uplo,n,a,lda,w,work,lwork,info)
   use iso_fortran_env, only : wp => real32
   integer, intent(in)    :: lda
   real(wp),intent(inout) :: a(lda,*)
   real(wp),intent(out)   :: w(*)
   character,intent(in)   :: jobz
   character,intent(in)   :: uplo
   integer, intent(out)   :: info
   integer, intent(in)    :: n
   real(wp),intent(inout) :: work(*)
   integer, intent(in)    :: lwork
   end subroutine ssyev
   pure subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
   use iso_fortran_env, only : wp => real64
   integer, intent(in)    :: lda
   real(wp),intent(inout) :: a(lda,*)
   real(wp),intent(out)   :: w(*)
   character,intent(in)   :: jobz
   character,intent(in)   :: uplo
   integer, intent(out)   :: info
   integer, intent(in)    :: n
   real(wp),intent(inout) :: work(*)
   integer, intent(in)    :: lwork
   end subroutine dsyev
end interface syev

interface spev
  pure subroutine sspev(jobz,uplo,n,ap,w,z,ldz,work,info)
   use iso_fortran_env, only : wp => real32
    real(wp), intent(inout) :: ap(*)
    real(wp), intent(out) :: w(*)
    character, intent(in) :: uplo
    real(wp), intent(out) :: z(ldz,*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    integer, intent(in) :: n
    integer, intent(in) :: ldz
    real(wp), intent(in) :: work(*)
  end subroutine sspev
  pure subroutine dspev(jobz,uplo,n,ap,w,z,ldz,work,info)
   use iso_fortran_env, only : wp => real64
    real(wp), intent(inout) :: ap(*)
    real(wp), intent(out) :: w(*)
    character, intent(in) :: uplo
    real(wp), intent(out) :: z(ldz,*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    integer, intent(in) :: n
    integer, intent(in) :: ldz
    real(wp), intent(in) :: work(*)
  end subroutine dspev
end interface spev

interface syevd
  pure subroutine ssyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,   &
     &                                                             info)
   use iso_fortran_env, only : wp => real32
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(out) :: w(*)
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    real(wp), intent(inout) :: work(*)
    integer, intent(in) :: lwork
    integer, intent(inout) :: iwork(*)
    integer, intent(in) :: liwork
  end subroutine ssyevd
  pure subroutine dsyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,   &
     &                                                             info)
   use iso_fortran_env, only : wp => real64
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(out) :: w(*)
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    real(wp), intent(inout) :: work(*)
    integer, intent(in) :: lwork
    integer, intent(inout) :: iwork(*)
    integer, intent(in) :: liwork
  end subroutine dsyevd
end interface syevd

interface spevd
  pure subroutine sspevd(jobz,uplo,n,ap,w,z,ldz,work,lwork,iwork,liwork,&
     &                                                             info)
   use iso_fortran_env, only : wp => real32
    real(wp), intent(inout) :: ap(*)
    real(wp), intent(out) :: w(*)
    character, intent(in) :: uplo
    real(wp), intent(out) :: z(ldz,*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    integer, intent(in) :: n
    integer, intent(in) :: ldz
    real(wp), intent(inout) :: work(*)
    integer, intent(in) :: lwork
    integer, intent(inout) :: iwork(*)
    integer, intent(in) :: liwork
  end subroutine sspevd
  pure subroutine dspevd(jobz,uplo,n,ap,w,z,ldz,work,lwork,iwork,liwork,&
     &                                                             info)
   use iso_fortran_env, only : wp => real64
    real(wp), intent(inout) :: ap(*)
    real(wp), intent(out) :: w(*)
    character, intent(in) :: uplo
    real(wp), intent(out) :: z(ldz,*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    integer, intent(in) :: n
    integer, intent(in) :: ldz
    real(wp), intent(inout) :: work(*)
    integer, intent(in) :: lwork
    integer, intent(inout) :: iwork(*)
    integer, intent(in) :: liwork
  end subroutine dspevd
end interface spevd

interface syevx
  pure subroutine ssyevx(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,&
     &                                z,ldz,work,lwork,iwork,ifail,info)
  use iso_fortran_env, only : wp => real32
  real(wp),intent(inout) :: a(lda,*)
  real(wp),intent(out)   :: w(*)
  character,intent(in)   :: uplo
  real(wp),intent(out)   :: z(ldz,*)
  real(wp),intent(in)    :: vl
  real(wp),intent(in)    :: vu
  integer, intent(in)    :: il
  integer, intent(in)    :: iu
  integer, intent(out)   :: m
  integer, intent(out)   :: ifail(*)
  real(wp),intent(in)    :: abstol
  integer, intent(out)   :: info
  character,intent(in)   :: jobz
  character,intent(in)   :: range
  integer, intent(in)    :: n
  integer, intent(in)    :: lda
  integer, intent(in)    :: ldz
  real(wp),intent(inout) :: work(*)
  integer, intent(in)    :: lwork
  integer, intent(in)    :: iwork(*)
  end subroutine ssyevx
  pure subroutine dsyevx(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,&
     &                                z,ldz,work,lwork,iwork,ifail,info)
  use iso_fortran_env, only : wp => real64
  real(wp),intent(inout) :: a(lda,*)
  real(wp),intent(out)   :: w(*)
  character,intent(in)   :: uplo
  real(wp),intent(out)   :: z(ldz,*)
  real(wp),intent(in)    :: vl
  real(wp),intent(in)    :: vu
  integer, intent(in)    :: il
  integer, intent(in)    :: iu
  integer, intent(out)   :: m
  integer, intent(out)   :: ifail(*)
  real(wp),intent(in)    :: abstol
  integer, intent(out)   :: info
  character,intent(in)   :: jobz
  character,intent(in)   :: range
  integer, intent(in)    :: n
  integer, intent(in)    :: lda
  integer, intent(in)    :: ldz
  real(wp),intent(inout) :: work(*)
  integer, intent(in)    :: lwork
  integer, intent(in)    :: iwork(*)
  end subroutine dsyevx
end interface syevx

interface spevx
  pure subroutine sspevx(jobz,range,uplo,n,ap,vl,vu,il,iu,abstol,m,w,z, &
     &                                        ldz,work,iwork,ifail,info)
  use iso_fortran_env, only : wp => real32
  real(wp),intent(inout) :: ap(*)
  real(wp),intent(out)   :: w(*)
  character,intent(in)   :: uplo
  real(wp),intent(out)   :: z(ldz,*)
  real(wp),intent(in)    :: vl
  real(wp),intent(in)    :: vu
  integer, intent(in)    :: il
  integer, intent(in)    :: iu
  integer, intent(out)   :: m
  integer, intent(out)   :: ifail(*)
  real(wp),intent(in)    :: abstol
  integer, intent(out)   :: info
  character,intent(in)   :: jobz
  character,intent(in)   :: range
  integer, intent(in)    :: n
  integer, intent(in)    :: ldz
  real(wp),intent(in)    :: work(*)
  integer, intent(in)    :: iwork(*)
  end subroutine sspevx
  pure subroutine dspevx(jobz,range,uplo,n,ap,vl,vu,il,iu,abstol,m,w,z, &
     &                                        ldz,work,iwork,ifail,info)
  use iso_fortran_env, only : wp => real64
  real(wp),intent(inout) :: ap(*)
  real(wp),intent(out)   :: w(*)
  character,intent(in)   :: uplo
  real(wp),intent(out)   :: z(ldz,*)
  real(wp),intent(in)    :: vl
  real(wp),intent(in)    :: vu
  integer, intent(in)    :: il
  integer, intent(in)    :: iu
  integer, intent(out)   :: m
  integer, intent(out)   :: ifail(*)
  real(wp),intent(in)    :: abstol
  integer, intent(out)   :: info
  character,intent(in)   :: jobz
  character,intent(in)   :: range
  integer, intent(in)    :: n
  integer, intent(in)    :: ldz
  real(wp),intent(in)    :: work(*)
  integer, intent(in)    :: iwork(*)
  end subroutine dspevx
end interface spevx

interface sygvd
   pure subroutine ssygvd(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,&
   &                      iwork,liwork,info)
   use iso_fortran_env, only : wp => real32
   real(wp),intent(inout) :: a(lda,*)
   real(wp),intent(inout) :: b(ldb,*)
   real(wp),intent(out)   :: w(*)
   integer, intent(in)    :: itype
   character,intent(in)   :: jobz
   character,intent(in)   :: uplo
   integer, intent(out)   :: info
   integer, intent(in)    :: n
   integer, intent(in)    :: lda
   integer, intent(in)    :: ldb
   real(wp),intent(inout) :: work(*)
   integer, intent(in)    :: lwork
   integer, intent(inout) :: iwork(*)
   integer, intent(in)    :: liwork
   end subroutine ssygvd
   pure subroutine dsygvd(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,&
   &                      iwork,liwork,info)
   use iso_fortran_env, only : wp => real64
   real(wp),intent(inout) :: a(lda,*)
   real(wp),intent(inout) :: b(ldb,*)
   real(wp),intent(out)   :: w(*)
   integer, intent(in)    :: itype
   character,intent(in)   :: jobz
   character,intent(in)   :: uplo
   integer, intent(out)   :: info
   integer, intent(in)    :: n
   integer, intent(in)    :: lda
   integer, intent(in)    :: ldb
   real(wp),intent(inout) :: work(*)
   integer, intent(in)    :: lwork
   integer, intent(inout) :: iwork(*)
   integer, intent(in)    :: liwork
   end subroutine dsygvd
end interface sygvd

interface symm
   pure subroutine ssymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
   use iso_fortran_env, only : wp => real32
   real(wp),intent(in)    :: a(lda,*)
   real(wp),intent(in)    :: b(ldb,*)
   real(wp),intent(inout) :: c(ldc,*)
   character,intent(in)   :: side
   character,intent(in)   :: uplo
   real(wp),intent(in)    :: alpha
   real(wp),intent(in)    :: beta
   integer, intent(in)    :: m
   integer, intent(in)    :: n
   integer, intent(in)    :: lda
   integer, intent(in)    :: ldb
   integer, intent(in)    :: ldc
   end subroutine ssymm
   pure subroutine dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
   use iso_fortran_env, only : wp => real64
   real(wp),intent(in)    :: a(lda,*)
   real(wp),intent(in)    :: b(ldb,*)
   real(wp),intent(inout) :: c(ldc,*)
   character,intent(in)   :: side
   character,intent(in)   :: uplo
   real(wp),intent(in)    :: alpha
   real(wp),intent(in)    :: beta
   integer, intent(in)    :: m
   integer, intent(in)    :: n
   integer, intent(in)    :: lda
   integer, intent(in)    :: ldb
   integer, intent(in)    :: ldc
   end subroutine dsymm
end interface symm

interface gemm
   pure subroutine sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta, &
   &                     c,ldc)
   use iso_fortran_env, only : wp => real32
   real(wp),intent(in)    :: a(lda,*)
   real(wp),intent(in)    :: b(ldb,*)
   real(wp),intent(inout) :: c(ldc,*)
   character,intent(in)   :: transa
   character,intent(in)   :: transb
   real(wp),intent(in)    :: alpha
   real(wp),intent(in)    :: beta
   integer, intent(in)    :: m
   integer, intent(in)    :: n
   integer, intent(in)    :: k
   integer, intent(in)    :: lda
   integer, intent(in)    :: ldb
   integer, intent(in)    :: ldc
   end subroutine sgemm
   pure subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta, &
   &                     c,ldc)
   use iso_fortran_env, only : wp => real64
   real(wp),intent(in)    :: a(lda,*)
   real(wp),intent(in)    :: b(ldb,*)
   real(wp),intent(inout) :: c(ldc,*)
   character,intent(in)   :: transa
   character,intent(in)   :: transb
   real(wp),intent(in)    :: alpha
   real(wp),intent(in)    :: beta
   integer, intent(in)    :: m
   integer, intent(in)    :: n
   integer, intent(in)    :: k
   integer, intent(in)    :: lda
   integer, intent(in)    :: ldb
   integer, intent(in)    :: ldc
   end subroutine dgemm
end interface gemm

interface symv
   pure subroutine ssymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
   use iso_fortran_env, only : wp => real32
   real(wp),intent(in)    :: a(lda,*)
   real(wp),intent(in)    :: x(*)
   real(wp),intent(inout) :: y(*)
   character,intent(in)   :: uplo
   real(wp),intent(in)    :: alpha
   real(wp),intent(in)    :: beta
   integer, intent(in)    :: incx
   integer, intent(in)    :: incy
   integer, intent(in)    :: n
   integer, intent(in)    :: lda
   end subroutine ssymv
   pure subroutine dsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
   use iso_fortran_env, only : wp => real64
   real(wp),intent(in)    :: a(lda,*)
   real(wp),intent(in)    :: x(*)
   real(wp),intent(inout) :: y(*)
   character,intent(in)   :: uplo
   real(wp),intent(in)    :: alpha
   real(wp),intent(in)    :: beta
   integer, intent(in)    :: incx
   integer, intent(in)    :: incy
   integer, intent(in)    :: n
   integer, intent(in)    :: lda
   end subroutine dsymv
end interface symv

interface gemv
   pure subroutine sgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
   use iso_fortran_env, only : wp => real32
   real(wp),intent(in)    :: a(lda,*)
   real(wp),intent(in)    :: x(*)
   real(wp),intent(inout) :: y(*)
   real(wp),intent(in)    :: alpha
   real(wp),intent(in)    :: beta
   character,intent(in)   :: trans
   integer, intent(in)    :: incx
   integer, intent(in)    :: incy
   integer, intent(in)    :: m
   integer, intent(in)    :: n
   integer, intent(in)    :: lda
   end subroutine sgemv
   pure subroutine dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
   use iso_fortran_env, only : wp => real64
   real(wp),intent(in)    :: a(lda,*)
   real(wp),intent(in)    :: x(*)
   real(wp),intent(inout) :: y(*)
   real(wp),intent(in)    :: alpha
   real(wp),intent(in)    :: beta
   character,intent(in)   :: trans
   integer, intent(in)    :: incx
   integer, intent(in)    :: incy
   integer, intent(in)    :: m
   integer, intent(in)    :: n
   integer, intent(in)    :: lda
   end subroutine dgemv
end interface gemv

interface axpy
   pure subroutine saxpy(n,a,x,incx,y,incy)
   use iso_fortran_env, only : wp => real32
   real(wp),intent(in) :: x(*)
   real(wp),intent(inout) :: y(*)
   real(wp),intent(in) :: a
   integer, intent(in) :: incx
   integer, intent(in) :: incy
   integer, intent(in) :: n
   end subroutine saxpy
   pure subroutine daxpy(n,a,x,incx,y,incy)
   use iso_fortran_env, only : wp => real64
   real(wp),intent(in) :: x(*)
   real(wp),intent(inout) :: y(*)
   real(wp),intent(in) :: a
   integer, intent(in) :: incx
   integer, intent(in) :: incy
   integer, intent(in) :: n
   end subroutine daxpy
end interface axpy

interface scal
   pure subroutine sscal(n,a,x,incx)
   use iso_fortran_env, only : wp => real32
   real(wp),intent(inout) :: x(*)
   real(wp),intent(in) :: a
   integer, intent(in) :: incx
   integer, intent(in) :: n
   end subroutine sscal
   pure subroutine dscal(n,a,x,incx)
   use iso_fortran_env, only : wp => real64
   real(wp),intent(inout) :: x(*)
   real(wp),intent(in) :: a
   integer, intent(in) :: incx
   integer, intent(in) :: n
   end subroutine dscal
end interface scal

interface copy
  pure subroutine scopy(n,x,incx,y,incy)
   use iso_fortran_env, only : wp => real32
    real(wp),intent(in) :: x(*)
    real(wp),intent(inout) :: y(*)
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    integer, intent(in) :: n
  end subroutine scopy
  pure subroutine dcopy(n,x,incx,y,incy)
   use iso_fortran_env, only : wp => real64
    real(wp),intent(in) :: x(*)
    real(wp),intent(inout) :: y(*)
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    integer, intent(in) :: n
  end subroutine dcopy
end interface copy

interface dot
   pure function sdot(n,x,incx,y,incy)
   use iso_fortran_env, only : wp => real32
   real(wp) :: sdot
   real(wp),intent(in) :: x(*)
   real(wp),intent(in) :: y(*)
   integer, intent(in) :: incx
   integer, intent(in) :: incy
   integer, intent(in) :: n
   end function sdot
   pure function ddot(n,x,incx,y,incy)
   use iso_fortran_env, only : wp => real64
   real(wp) :: ddot
   real(wp),intent(in) :: x(*)
   real(wp),intent(in) :: y(*)
   integer, intent(in) :: incx
   integer, intent(in) :: incy
   integer, intent(in) :: n
   end function ddot
end interface dot

interface htosq
   module procedure dhtosq
end interface htosq

interface syprj
   module procedure dsyprj
end interface syprj

interface blckmgs
   module procedure dblckmgs
end interface blckmgs

contains

subroutine dhtosq(n,a,b)
   use iso_fortran_env, only : wp => real64
   implicit none
! ---------------------------------------------------------------------
!      expand trigonal matrix b to full matrix a
!      a and b may refer to the same address
! ---------------------------------------------------------------------
   integer, intent(in)  :: n
   real(wp),intent(out) :: a(n,n)
   real(wp),intent(in)  :: b(n*(n+1)/2)
   integer :: i,j,ioff

   do i=n,1,-1
      ioff=i*(i-1)/2
      do j=i,1,-1
         a(j,i)=b(ioff+j)
      enddo
   enddo

   do i=1,n
      do j=1,i-1
         a(i,j)=a(j,i)
      enddo
   enddo

end subroutine dhtosq

subroutine dsyprj(nbdim,m,bmat,n,asym)
!----------------------------------------------------------------------
! Purpose:
! Performs projection of a real symmetric matrix asym packed in an
! upper triangular form:
! asym = (1 - bmat*bmat')*asym*(1 - bmat*bmat')
! where (1 - bmat*bmat') is a projector constructed from bmat
! The matrix asym is overwritten on output; bmat is left unchanged.
! Operation count is proportional to n*n*m
!
! Input/Output:
! nbdim - First dimension of bmat as declared in the calling routine
! m     - Number of columns in the matrix bmat
! bmat  - n*m matrix used to build projector
! n     - Dimension of the matrix asym
! asym  - Symmetric matrix in packed upper triangular form
!----------------------------------------------------------------------
   use iso_fortran_env, only : wp => real64
  implicit none

! Input/Output:
  integer, intent(in)    :: nbdim,m,n
  real(wp),intent(in)    :: bmat(nbdim,m)
  real(wp),intent(inout) :: asym(n*(n+1)/2)

! Local:
  integer  :: i,j,ij
  real(wp),allocatable :: scrb(:,:)
  real(wp),allocatable :: scra(:,:)
!----------------------------------------------------------------------
  allocate( scrb(n,m), scra(n,n) )

! Expand trigonal matrix asym to full matrix on scra
  call htosq(n,scra,asym)

! Calculate scrb = asym*bmat
  call symm('l','u',n,m,1.0_wp,scra,n,bmat,nbdim,0.0_wp,scrb,n)
  
! Calculate scra = scrb*bmat'
  call gemm('n','t',n,n,m,1.0_wp,scrb,n,bmat,nbdim,0.0_wp,scra,n)

! Calculate asym = asym - scra - scra'
  do i=1,n
     do j=1,i
        ij = i*(i-1)/2 + j
        asym(ij) = asym(ij) - scra(i,j) - scra(j,i)
     end do
  end do
     
! Calculate scrb' = scra'*bmat
  call gemm('t','n',n,m,n,1.0_wp,scra,n,bmat,nbdim,0.0_wp,scrb,n)

! Calculate scra = bmat*scrb'
  call gemm('n','t',n,n,m,1.0_wp,bmat,nbdim,scrb,n,0.0_wp,scra,n)
  
! Calculate asym = asym + scra
  do i=1,n
     do j=1,i
        ij = i*(i-1)/2 + j
        asym(ij) = asym(ij) + scra(i,j)
     end do
  end do

  deallocate( scra, scrb )
     
end subroutine dsyprj

subroutine dblckmgs(m,n,ndim,darray)
!----------------------------------------------------------------------
! Purpose:
! Subroutine performs modified Gramm-Schmidt orthonormalization of a 
! real matrix. Orthonormalization is done in-place, so the darray is 
! overwritten on exit. Linearly dependent vectors are set to zero.
!
! Input:
! m      - Number of rows in the matrix darray
! n      - Number of columns in the matrix darray
! ndim   - First array dimension as declared in the calling routine
! darray - Array to be orthonormalized
!
! Output:
! darray - Orthonormalized array
!----------------------------------------------------------------------
   use iso_fortran_env, only : wp => real64
  implicit none

! Input/Output:
  integer, intent(in) :: m,n,ndim
  real(wp), dimension(ndim,n), intent(out) :: darray
!----------------------------------------------------------------------  
! Local:
  integer :: ii,jj,kk,ll,ibsize,nblcks,istrt,jstrt,iend,ncol,ierr
  real(wp) :: tmp
  real(wp), dimension(:,:), allocatable :: smat
!---------------------------------------------------------------------- 
! Local parameters
  real(wp) thr
! use BLAS directly without MCTC interfaces, for dirty hacks
  external dgemm ! I hope you know what you are doing...
!----------------------------------------------------------------------  

  thr = epsilon(1.0_wp) !Threshold for zero vectors
!----------------------------------------------------------------------  
! Block size optimized for Athlon 1200 MHz with 2.0GB memory for 
! matrices up to 5000x5000
  ibsize = 60
!----------------------------------------------------------------------  

!-----------------------------------------------------------------
! Allocate overlap matrix
!-----------------------------------------------------------------
  allocate(smat(ibsize,ibsize),stat=ierr)
  if(ierr /= 0)  call raise('E','Memory allocation error in blckmgs',1)

!-----------------------------------------------------------------
! Calculate the number of blocks
!-----------------------------------------------------------------
  nblcks = (n+ibsize-1)/ibsize
  ibsize = min(n,ibsize)

!-----------------------------------------------------------------
! Orthogonalize the first block using modified schmidt
!-----------------------------------------------------------------
  do ii=1,ibsize

     tmp = dot(m,darray(:,ii),1,darray(:,ii),1)

! Linear dependence
     if(tmp < thr) then
        darray(1:m,ii) = 0.0_wp
        cycle
     end if

     tmp = 1.0_wp/sqrt(tmp)
     call scal(m,tmp,darray(:,ii),1)

     do jj=ii+1,ibsize
        tmp = dot(m,darray(:,ii),1,darray(:,jj),1)
        call axpy(m,-tmp,darray(:,ii),1,darray(:,jj),1)
     end do

  end do

!-----------------------------------------------------------------
! Loop over remaining blocks
!-----------------------------------------------------------------
  do ii=1,nblcks-1

! Initial and final column and number of columns in the block ii+1
     istrt = ii*ibsize+1
     iend  = (ii+1)*ibsize
     iend  = min(n,iend)
     ncol  = iend - istrt + 1

! Orthogonalize the block ii+1 against the previous ones
     do jj=1,ii

! Initial index of the block jj
        jstrt = (jj-1)*ibsize+1

        call dgemm('t','n',ibsize,ncol,m,1.0_wp,darray(1,jstrt),ndim,darray(1,istrt),ndim,0.0_wp,smat,ibsize)
        call dgemm('n','n',m,ncol,ibsize,-1.0_wp,darray(1,jstrt),ndim,smat,ibsize,1.0_wp,darray(1,istrt),ndim)

     end do

! Othogonalize vectors on the block ii+1 among themself using modified schmidt
     do kk=istrt,iend

        tmp = dot(m,darray(:,kk),1,darray(:,kk),1)

! Linear dependence
        if(tmp < thr) then
           darray(1:m,kk) = 0.0_wp
           cycle
        end if
     
        tmp = 1.0_wp/sqrt(tmp)
        call scal(m,tmp,darray(:,kk),1)

        do ll=kk+1,iend
           tmp = dot(m,darray(:,kk),1,darray(:,ll),1)
           call axpy(m,-tmp,darray(:,kk),1,darray(:,ll),1)
        end do

     end do

  end do

! Clean up
  deallocate(smat,stat=ierr)
  if(ierr /= 0) call raise('E','Memory deallocation error in blckmgs',1)

end subroutine dblckmgs

function ssyluinv(Amat,m) result(info)
   use iso_fortran_env, only : wp => real32
   implicit none
   integer, intent(in) :: m
   real(wp), intent(inout) :: Amat(m,m)
   integer, allocatable :: ipiv(:)
   real(wp),allocatable :: temp(:)
   real(wp),allocatable :: work(:)
   integer  :: lwork
   integer  :: info
   real(wp) :: test(1)
   integer  :: i,j

   ! assume work space query, set best value to test after first dsytrf call
   call ssytrf('L',m,Amat,m,ipiv,test,-1,info)
   lwork=int(test(1))
   allocate( work(lwork), source = 0.0_wp )

   ! Bunch-Kaufman factorization A = L*D*L**T
   call ssytrf('L',m,Amat,m,ipiv,work,lwork,info)
   if(info > 0) return

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Amat matrix
   ! Amat matrix is overwritten with lower triangular part of A⁻¹
   call ssytri('L',m,Amat,m,ipiv,work,info)
   if (info > 0) return

   ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
   do i = 1, m
      do j = i+1, m
         Amat(i,j)=Amat(j,i)
      enddo
   enddo

end function ssyluinv

function dsyluinv(Amat,m) result(info)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in) :: m
   real(wp), intent(inout) :: Amat(m,m)
   integer, allocatable :: ipiv(:)
   real(wp),allocatable :: temp(:)
   real(wp),allocatable :: work(:)
   integer  :: lwork
   integer  :: info
   real(wp) :: test(1)
   integer  :: i,j

   allocate( ipiv(m), source = 0 )

   ! assume work space query, set best value to test after first dsytrf call
   call dsytrf('L',m,Amat,m,ipiv,test,-1,info)
   lwork=int(test(1))
   allocate( work(lwork), source = 0.0_wp )

   ! Bunch-Kaufman factorization A = L*D*L**T
   call dsytrf('L',m,Amat,m,ipiv,work,lwork,info)
   if(info > 0) return

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Amat matrix
   ! Amat matrix is overwritten with lower triangular part of A⁻¹
   call dsytri('L',m,Amat,m,ipiv,work,info)
   if (info > 0) return

   ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
   do i = 1, m
      do j = i+1, m
         Amat(i,j)=Amat(j,i)
      enddo
   enddo

end function dsyluinv

end module mctc_la
