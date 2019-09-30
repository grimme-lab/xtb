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

! Modified Boryden mixer according to
! D.D. Johnson, PRB 38, 12807 (1988)
! Input Variables:
! iter = current iteration
! maxiter = max number of iterations
! n = number of atoms (charges)
! q = current iteration INPUT charges (before diag)
! qlast = last iteration INPUT charges (before diag)
! dqlast = charge difference of prev. iteration (q_in - q_out)
! dq = current iteration charge difference (q_in - q_out)
! df = previous |dF>
! u = previous |u>
! a = a matrix
! Note: The modified Broyden method works on the INPUT charges of the scf iteration, and
! the difference produced by the SCC step. Especially the INPUT charge of the previous
! iteration is needed, which is usually not stored!
       
      subroutine broyden(n,q,qlast,dq,dqlast,iter,maxiter,
     &                   alpha,omega,df,u,a)
         use iso_fortran_env, wp => real64
      implicit none
      
      integer i,j
      integer n,iter,maxiter,it1
      real(wp) inv
      
!     charge arrays      
      real(wp) q(n),qlast(n),dq(n),dqlast(n)
      

!     |dF>, |u> and a from prev. iterations      
      real(wp) df(maxiter,n),u(maxiter,n),a(maxiter,maxiter)
      real(wp) omega(maxiter)
!     on the fly variables needed for the mixing
      real(wp), allocatable :: beta(:,:), c(:), gamma(:)
!      real(wp), allocatable :: beta1(:,:)

!     mixing parameters

      real(wp) alpha
      real(wp) omega0
      real(wp) minw,maxw,wfac
      real(wp) tmp(n)




      it1 = iter - 1



!     set parameters 

!     alpha = 0.25d0
      omega0 = 0.01d0
      minw = 1.0d0
      maxw = 100000.0d0
      wfac = 0.01d0
!     wfac = 0.05d0

!     if case for first iteration: simple damping
      if (iter == 1 ) then
         dqlast = dq
         qlast = q
         q = q + alpha * dq
         return
      endif

!     allocate
      allocate(beta(it1,it1),c(it1),gamma(it1))
!     allocate(beta1(it1,it1))



!     create omega (weight) for the current iteration
      omega(it1) = sqrt(dot_product(dq,dq))

      if (omega(it1) .gt. (wfac / maxw)) then 
         omega(it1) = wfac / omega(it1)
      else
         omega(it1) = maxw
      endif
      if (omega(it1) .lt. minw) then
         omega(it1) = minw
      endif




!     Build dF(iter-1)
      df(it1,:)=dq-dqlast
      inv=1.0d0 / sqrt(dot_product(df(it1,:),df(it1,:)))
      df(it1,:)=inv*df(it1,:)



!     Next: build a, beta, c, gamma
      do i=1,it1
         a(i,it1) = dot_product(df(i,:),df(it1,:))
         a(it1,i) = a(i,it1)
         c(i) = omega(i) * dot_product(df(i,:),dq) 
      enddo



!     Build beta from a and omega      
      do i = 1, it1
         beta(:it1, i) = omega(:it1) * omega(i) * a(:it1,i)
         beta(i, i) = beta(i, i) + omega0*omega0
      enddo

!     build beta^-1
      call matinv(beta,it1) 

!     build gamma
      gamma = matmul(c,beta)

!     Build |u>      
      u(it1,:) = alpha * df(it1,:) + inv * (q-qlast) !!!

!     save charges and deltas
      dqlast = dq
      qlast = q

!     calculate new charges      
      q = q + alpha * dq

      do i=1, it1
         q = q - omega(i) * gamma(i) * u(i,:)
      enddo

      deallocate(beta,c,gamma)
      
      end subroutine broyden
     
   
   
   
      subroutine matinv(a, nrow)
         use iso_fortran_env, wp => real64
      integer nRow
      real(wp) a(nrow,nrow)
      
      integer info
      integer, allocatable :: ipiv(:)
      real(wp), allocatable :: work(:)
      
      
      allocate(ipiv(nrow),work(nrow))
!     LU decomoposition of a general matrix      
      call dgetrf(nrow,nrow,a,nrow,ipiv,info)
      if (info == 0) then
!     generate inverse of a matrix given its LU decomposition      
         call dgetri(nrow,a,nrow,ipiv,work,nrow,info)
      endif
      deallocate(ipiv,work)
      if(info .ne. 0)then
         write(*,*)"Error in Broyden matrix inversion!"
         write(*,*)"Error code",info
         stop
      endif
 
   
      end subroutine matinv
      
