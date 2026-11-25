! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

module xtb_mctc_la
   use xtb_mctc_accuracy, only : sp, dp
   use xtb_mctc_lapack
   use xtb_mctc_blas
   use xtb_mctc_blas_wrap3, only : contract => mctc_gemm
   use xtb_mctc_blas_wrap2, only : contract => mctc_gemv
   implicit none
   public ! Forward lapack/blas module

   private :: sp, dp




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
   implicit none
! ---------------------------------------------------------------------
!      expand trigonal matrix b to full matrix a
!      a and b may refer to the same address
! ---------------------------------------------------------------------
   integer, intent(in)  :: n
   real(dp),intent(out) :: a(n,n)
   real(dp),intent(in)  :: b(n*(n+1)/2)
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
  implicit none

! Input/Output:
  integer, intent(in)    :: nbdim,m,n
  real(dp),intent(in)    :: bmat(nbdim,m)
  real(dp),intent(inout) :: asym(n*(n+1)/2)

! Local:
  integer  :: i,j,ij
  real(dp),allocatable :: scrb(:,:)
  real(dp),allocatable :: scra(:,:)
!----------------------------------------------------------------------
  allocate( scrb(n,m), scra(n,n) )

! Expand trigonal matrix asym to full matrix on scra
  call htosq(n,scra,asym)

! Calculate scrb = asym*bmat
  call mctc_symm(scra, bmat, scrb, side='l', uplo='u', m=n, n=m)
  
! Calculate scra = scrb*bmat'
  call mctc_gemm(scrb, bmat, scra, transb='t', m=n, n=n, k=m)

! Calculate asym = asym - scra - scra'
  do i=1,n
     do j=1,i
        ij = i*(i-1)/2 + j
        asym(ij) = asym(ij) - scra(i,j) - scra(j,i)
     end do
  end do
     
! Calculate scrb' = scra'*bmat
  call mctc_gemm(scra, bmat, scrb, transa='t', m=n, n=m, k=n)

! Calculate scra = bmat*scrb'
  call mctc_gemm(bmat, scrb, scra, transb='t', m=n, n=n, k=m)
  
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
  implicit none

! Input/Output:
  integer, intent(in) :: m,n,ndim
  real(dp), dimension(ndim,n), intent(out) :: darray
!----------------------------------------------------------------------  
! Local:
  integer :: ii,jj,kk,ll,ibsize,nblcks,istrt,jstrt,iend,ncol,ierr
  real(dp) :: tmp
  real(dp), dimension(:,:), allocatable :: smat
!---------------------------------------------------------------------- 
! Local parameters
  real(dp) thr
!----------------------------------------------------------------------  

  thr = epsilon(1.0_dp) !Threshold for zero vectors
!----------------------------------------------------------------------  
! Block size optimized for Athlon 1200 MHz with 2.0GB memory for 
! matrices up to 5000x5000
  ibsize = 60
!----------------------------------------------------------------------  

!-----------------------------------------------------------------
! Allocate overlap matrix
!-----------------------------------------------------------------
  allocate(smat(ibsize,ibsize),stat=ierr)
  if(ierr /= 0)  call raise('E','Memory allocation error in blckmgs')

!-----------------------------------------------------------------
! Calculate the number of blocks
!-----------------------------------------------------------------
  nblcks = (n+ibsize-1)/ibsize
  ibsize = min(n,ibsize)

!-----------------------------------------------------------------
! Orthogonalize the first block using modified schmidt
!-----------------------------------------------------------------
  do ii=1,ibsize

     tmp = mctc_dot(darray(:,ii), darray(:,ii), n=m)

! Linear dependence
     if(tmp < thr) then
        darray(1:m,ii) = 0.0_dp
        cycle
     end if

     tmp = 1.0_dp/sqrt(tmp)
     call mctc_scal(darray(:,ii), tmp, n=m)

     do jj=ii+1,ibsize
        tmp = mctc_dot(darray(:,ii), darray(:,jj), n=m)
        call mctc_axpy(darray(:,ii), darray(:,jj), alpha=-tmp, n=m)
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

        call mctc_gemm(darray(:,jstrt:), darray(:,istrt:), smat, &
           & transa='t', m=ibsize, n=ncol, k=m)
        call mctc_gemm(darray(:,jstrt:), smat, darray(:,istrt:), &
           & alpha=-1.0_dp, beta=1.0_dp, m=m, n=ncol, k=ibsize)

     end do

! Othogonalize vectors on the block ii+1 among themself using modified schmidt
     do kk=istrt,iend

        tmp = mctc_dot(darray(:,kk), darray(:,kk), n=m)

! Linear dependence
        if(tmp < thr) then
           darray(1:m,kk) = 0.0_dp
           cycle
        end if
     
        tmp = 1.0_dp/sqrt(tmp)
        call mctc_scal(darray(:,kk), tmp, n=m)

        do ll=kk+1,iend
           tmp = mctc_dot(darray(:,kk), darray(:,ll), n=m)
           call mctc_axpy(darray(:,kk), darray(:,ll), alpha=-tmp, n=m)
        end do

     end do

  end do

! Clean up
  deallocate(smat,stat=ierr)
  if(ierr /= 0) call raise('E','Memory deallocation error in blckmgs')

end subroutine dblckmgs

function ssyluinv(Amat,m) result(info)
   implicit none
   integer, intent(in) :: m
   real(sp), intent(inout) :: Amat(m,m)
   integer, allocatable :: ipiv(:)
   integer  :: info
   integer  :: i,j

   allocate( ipiv(m), source = 0 )

   ! Bunch-Kaufman factorization A = L*D*L**T
   call mctc_sytrf(Amat, ipiv, uplo='L', info=info)
   if (info /= 0) return

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Amat matrix
   ! Amat matrix is overwritten with lower triangular part of A⁻¹
   call mctc_sytri(Amat, ipiv, uplo='L', info=info)
   if (info /= 0) return

   ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
   do i = 1, m
      do j = i+1, m
         Amat(i,j)=Amat(j,i)
      enddo
   enddo

end function ssyluinv

function dsyluinv(Amat,m) result(info)
   implicit none
   integer, intent(in) :: m
   real(dp), intent(inout) :: Amat(m,m)
   integer, allocatable :: ipiv(:)
   integer  :: info
   integer  :: i,j

   allocate( ipiv(m), source = 0 )

   ! Bunch-Kaufman factorization A = L*D*L**T
   call mctc_sytrf(Amat, ipiv, uplo='L', info=info)
   if (info /= 0) return

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Amat matrix
   ! Amat matrix is overwritten with lower triangular part of A⁻¹
   call mctc_sytri(Amat, ipiv, uplo='L', info=info)
   if (info /= 0) return

   ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
   do i = 1, m
      do j = i+1, m
         Amat(i,j)=Amat(j,i)
      enddo
   enddo

end function dsyluinv


end module xtb_mctc_la
