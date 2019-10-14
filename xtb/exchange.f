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

! compute open-shell exchange correction for high-spin case
! exchange integrals are computed ala sTDA but using Mulliken
! transition densities for convenience

      subroutine exch(nat,at,nao,nopen,ihomo,xyz,focc,s,cmo,xint,aoat)
      use setparam
      use aoparam
      implicit none          
      integer nao,nat,at(nat),nopen,ihomo
      real*8  cmo(nao,nao),focc(nao),xyz(3,nat),s(nao,nao)
      real*8  xint
      integer, intent(in) :: aoat(nao)

      real*8, allocatable ::q12(:),qint(:)
      real*8, allocatable ::gamk(:,:)

      integer mo,n,i,j,k,idum,olist(nao),o1,o2,ii,jj,oo1,oo2
      real*8 xk,rabx,alphak,ddot

      alphak=2.0 ! sTDA-xTB value

      allocate(gamk(nat,nat),q12(nat),qint(nat))

      do i=1,nat
         ii=at(i)
         do j=1,i      
            jj=at(j)
            xk  =0.50d0*(gam(ii)+gam(jj)) 
            rabx=sqrt((xyz(1,i)-xyz(1,j))**2
     .               +(xyz(2,i)-xyz(2,j))**2
     .               +(xyz(3,i)-xyz(3,j))**2)
            gamk(j,i)=1./(rabx**alphak+1./xk**alphak)**(1.0d0/alphak)
            gamk(i,j)=gamk(j,i)
         enddo
      enddo

! active MOs      
      j=0
      do i=ihomo-nopen+1,ihomo
         j=j+1          
         olist(j)=i
      enddo

! sum over all open-shell pairs the exchange integrals
      xint=0
      do o1=1,nopen-1
         oo1=olist(o1)
         do o2=o1+1,nopen
            oo2=olist(o2)
! Mulliken transition density for MOs oo1, oo2
            call pop12(nat,nao,aoat,oo1,oo2,S,cmo,q12)
            call dsymv('l',nat,1.0d0,gamk,nat,q12,1,0.0,qint,1)
            xk=ddot(nat,q12,1,qint,1)
            xint=xint+xk
         enddo
      enddo

! empirical scaling obtained by adjusting to a few PAH S0-T1 gaps     
! values between 0.25 and 0.35 seem to be reasonable for the standard
! GFN-xTB parameters for T1 states, HOMO-LUMO S1 excitation energies in dye12
! set are reproduced to an MAD of 0.2 eV with a factor of -1.4     

      xint=xint*ex_open

      end

      subroutine pop12(n,nao,aoat,mo1,mo2,S,C,q)
      implicit none
      real*8  C(nao,nao)
      real*8  S(nao,nao)
      real*8  occ(nao)
      real*8  q(n)
      integer nao,n,aoat(nao),mo1,mo2
   
      integer i,j,ii,jj
      real*8  cc

      q=0
      do i=1,nao 
         ii=aoat(i)
         do j=1,i-1
            jj=aoat(j)
            cc=S(j,i)*C(i,mo1)*C(j,mo2)
            q(jj)=q(jj)+cc
            q(ii)=q(ii)+cc
         enddo
         cc=S(i,i)*C(i,mo1)*C(i,mo2)
         q(ii)=q(ii)+cc
      enddo

      end
