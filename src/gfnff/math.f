! This file is part of xtb.
!
! Copyright (C) 2019-2020 Stefan Grimme
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


      subroutine vsub(a,b,c,n)
      implicit none
      integer n
      real*8 a(n),b(n),c(n)
      integer i

      do i=1,n
         c(i) = a(i) - b(i)
      end do

      return
      end subroutine vsub

c  .....................................................................

      real*8 function vlen(a)
      implicit none !double precision (a-h,o-z)
      real*8 a(3)
      real*8 tot

      tot = a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
      vlen = 0.0d0
      if (tot.gt.0.0d0) vlen = dsqrt(tot)

      return
      end function vlen

c  .....................................................................

      real*8 function valijklff(natoms,xyz,i,j,k,l)

      implicit none

      external
     ,     vecnorm,valijk

      integer
     ,     ic,i,j,k,l,natoms

      real*8
     ,     xyz(3,natoms),
     ,     eps,ra(3),rb(3),rc(3),na(3),nb(3),
     ,     rab,rbc,thab,thbc,valijk,
     ,     vecnorm,nan,nbn,rcn,snanb,deter,pi

      parameter (eps=1.0d-14)
      data pi/3.1415926535897932384626433832795029d0/

c ... get torsion coordinate
      do ic=1,3
        ra(ic)=xyz(ic,j)-xyz(ic,i)
        rb(ic)=xyz(ic,k)-xyz(ic,j)
        rc(ic)=xyz(ic,l)-xyz(ic,k)
      end do

c ... determinante of rb,ra,rc
      deter= ra(1)*(rb(2)*rc(3)-rb(3)*rc(2))
     ,      -ra(2)*(rb(1)*rc(3)-rb(3)*rc(1))
     ,      +ra(3)*(rb(1)*rc(2)-rb(2)*rc(1))

      thab=valijk(natoms,xyz,i,k,j)
      thbc=valijk(natoms,xyz,j,l,k)
      call crossprod(ra,rb,na)
      call crossprod(rb,rc,nb)
      nan=vecnorm(na,3,1)
      nbn=vecnorm(nb,3,1)

      snanb=0.0d0
      do ic=1,3
        snanb=snanb+na(ic)*nb(ic)
      end do
      if (abs(abs(snanb)-1.d0).lt.eps) then
        snanb=sign(1.d0,snanb)
      end if

      valijklff=acos(snanb)

      end function valijklff
