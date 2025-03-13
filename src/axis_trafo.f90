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

module xtb_axis
   use xtb_mctc_accuracy, only : wp
contains
  subroutine axis(numat,nat,xyz,aa,bb,cc)
    use xtb_splitparam
    implicit integer (i-n)
    implicit double precision (a-h,o-z)
    dimension xyz(3,numat)
    integer nat(numat)

    PARAMETER (BOHR=0.52917726)
    dimension t(6), rot(3), xyzmom(3), eig(3), evec(3,3)
    dimension x(numat),y(numat),z(numat),coord(3,numat)
    data t /6*0.d0/
    !************************************************************************
    !*     const1 =  10**40/(n*a*a)
    !*               n = avergadro's number
    !*               a = cm in an angstrom
    !*               10**40 is to allow units to be 10**(-40)gram-cm**2
    !*
    !************************************************************************
    const1 = 1.66053d0
    !************************************************************************
    !*
    !*     const2 = conversion factor from angstrom-amu to cm**(-1)
    !*
    !*            = (planck's constant*n*10**16)/(8*pi*pi*c)
    !*            = 6.62618*10**(-27)[erg-sec]*6.02205*10**23*10**16/
    !*              (8*(3.1415926535)**2*2.997925*10**10[cm/sec])
    !*
    !***********************************************************************
    const2=16.8576522d0

    sumw=1.d-20
    sumwx=0.d0
    sumwy=0.d0
    sumwz=0.d0

    coord(1:3,1:numat)=xyz(1:3,1:numat)*bohr

    do i=1,numat
       sumw=sumw+atmass(i)
       sumwx=sumwx+atmass(i)*coord(1,i)
       sumwy=sumwy+atmass(i)*coord(2,i)
       sumwz=sumwz+atmass(i)*coord(3,i)
    enddo

    sumwx=sumwx/sumw
    sumwy=sumwy/sumw
    sumwz=sumwz/sumw
    f=1.0d0/bohr
    do  i=1,numat
       x(i)=coord(1,i)-sumwx
       y(i)=coord(2,i)-sumwy
       z(i)=coord(3,i)-sumwz
    enddo

    !************************************************************************
    !*    matrix for moments of inertia is of form
    !*
    !*           |   y**2+z**2                         |
    !*           |    -y*x       z**2+x**2             | -i =0
    !*           |    -z*x        -z*y       x**2+y**2 |
    !*
    !************************************************************************
    do i=1,6
       t(i)=dble(i)*1.0d-10
    enddo
    do i=1,numat
       t(1)=t(1)+atmass(i)*(y(i)**2+z(i)**2)
       t(2)=t(2)-atmass(i)*x(i)*y(i)
       t(3)=t(3)+atmass(i)*(z(i)**2+x(i)**2)
       t(4)=t(4)-atmass(i)*z(i)*x(i)
       t(5)=t(5)-atmass(i)*y(i)*z(i)
       t(6)=t(6)+atmass(i)*(x(i)**2+y(i)**2)
    enddo
    call rsp(t,3,3,eig,evec)
    !     xsum=evec(1,1)*(evec(2,2)*evec(3,3)-evec(3,2)*evec(2,3)) +
    !    1     evec(1,2)*(evec(2,3)*evec(3,1)-evec(2,1)*evec(3,3)) +
    !    2     evec(1,3)*(evec(2,1)*evec(3,2)-evec(2,2)*evec(3,1))
    do i=1,3
       if(eig(i).lt.3.d-4) then
          eig(i)=0.d0
          rot(i)=0.d0
       else
          rot(i)=2.9979245d+4*const2/eig(i)
       endif
       xyzmom(i)=eig(i)*const1
    enddo

    aa=rot(3)/2.9979245d+4
    bb=rot(2)/2.9979245d+4
    cc=rot(1)/2.9979245d+4
    avmom=1.d-47*(xyzmom(1)+xyzmom(2)+xyzmom(3))/3.

    return
  end subroutine axis



subroutine axis2(n,xyz,aa,bb,cc,avmom,sumw)
   
   use xtb_splitparam, only : atmass
   use xtb_mctc_convert, only : autoaa
   
   implicit double precision (a-h,o-z)
   
   !> number of atoms 
   integer, intent(in) :: n
   
   !> cartesian coordinates
   real(wp), intent(in) :: xyz(:,:)

   !> rotational constants
   real(wp), intent(out) :: aa,bb,cc

   !> average moment of inertia
   real(wp), intent(out) :: avmom

   !> sum of atomic masses
   real(wp), intent(out) :: sumw
   
   !> const1 =  10**40/(n*a*a)
   !>            n = avergadro's number
   !>            a = cm in an angstrom
   !>            10**40 is to allow units to be 10**(-40)gram-cm**2
   real(wp), parameter :: const1 = 1.66053_wp

   !> const2 = conversion factor from angstrom-amu to cm**(-1)
   !>        = (planck's constant*n*10**16)/(8*pi*pi*c)
   !>        = 6.62618*10**(-27)[erg-sec]*6.02205*10**23*10**16/
   !>         (8*(3.1415926535)**2*2.997925*10**10[cm/sec])
   real(wp), parameter :: const2 = 16.8576522_wp

   !> temporary variables
   real(wp) :: xsum,eps,ams
   
   !> weighted sum of coordinates
   real(wp) :: sumwx,sumwy,sumwz

   !> loop variables
   integer :: i

   !> temporary arrays
   real(wp) :: t(6), rot(3), xyzmom(3), eig(3), evec(3,3)
   real(wp) :: x(n),y(n),z(n)

   !> Cartesian coordinates in Bohr
   real(wp) :: coord(3,n)  
   
   t(:)  = 0.0_wp 
   sumw  = 0.0_wp
   sumwx = 0.0_wp
   sumwy = 0.0_wp
   sumwz = 0.0_wp

   ! convert Bohr to Angstrom ! 
   coord(:,:) = xyz(:,:) * autoaa
   
   ! cma !
   do i=1,n
      sumw=sumw+atmass(i)
      sumwx=sumwx+atmass(i)*coord(1,i)
      sumwy=sumwy+atmass(i)*coord(2,i)
      sumwz=sumwz+atmass(i)*coord(3,i)
   enddo

   sumwx=sumwx/sumw
   sumwy=sumwy/sumw
   sumwz=sumwz/sumw

   do i=1,n
      x(i)=coord(1,i)-sumwx
      y(i)=coord(2,i)-sumwy
      z(i)=coord(3,i)-sumwz
   enddo

   !************************************************************************
   !*    matrix for moments of inertia is of form
   !*
   !*           |   y**2+z**2                         |
   !*           |    -y*x       z**2+x**2             | -i =0
   !*           |    -z*x        -z*y       x**2+y**2 |
   !*
   !************************************************************************
   do i=1,6
      t(i)=dble(i)*1.0d-10
   enddo

   do i=1,n
      t(1)=t(1)+atmass(i)*(y(i)**2+z(i)**2)
      t(2)=t(2)-atmass(i)*x(i)*y(i)
      t(3)=t(3)+atmass(i)*(z(i)**2+x(i)**2)
      t(4)=t(4)-atmass(i)*z(i)*x(i)
      t(5)=t(5)-atmass(i)*y(i)*z(i)
      t(6)=t(6)+atmass(i)*(x(i)**2+y(i)**2)
   enddo

   call rsp(t,3,3,eig,evec)
   
   do i=1,3

      if(eig(i).lt.3.d-4) then
         eig(i)=0.d0
         rot(i)=0.d0
      else
         rot(i)=2.9979245d+4*const2/eig(i)
      endif

      xyzmom(i)=eig(i)*const1

   enddo

   aa=rot(3)/2.9979245d+4
   bb=rot(2)/2.9979245d+4
   cc=rot(1)/2.9979245d+4

   avmom=1.d-47*(xyzmom(1)+xyzmom(2)+xyzmom(3))/3.

end subroutine axis2

  subroutine axisvec(numat,nat,xyz,aa,bb,cc,evec)
    use xtb_splitparam
    implicit integer (i-n)
    implicit double precision (a-h,o-z)
    dimension xyz(3,numat)
    integer nat(numat)

    PARAMETER (BOHR=0.52917726)
    dimension t(6), rot(3), xyzmom(3), eig(3), evec(3,3)
    dimension x(numat),y(numat),z(numat),coord(3,numat)
    data t /6*0.d0/
    !************************************************************************
    !*     const1 =  10**40/(n*a*a)
    !*               n = avergadro's number
    !*               a = cm in an angstrom
    !*               10**40 is to allow units to be 10**(-40)gram-cm**2
    !*
    !************************************************************************
    const1 = 1.66053d0
    !************************************************************************
    !*
    !*     const2 = conversion factor from angstrom-amu to cm**(-1)
    !*
    !*            = (planck's constant*n*10**16)/(8*pi*pi*c)
    !*            = 6.62618*10**(-27)[erg-sec]*6.02205*10**23*10**16/
    !*              (8*(3.1415926535)**2*2.997925*10**10[cm/sec])
    !*
    !************************************************************************
    const2=16.8576522d0

    sumw=1.d-20
    sumwx=0.d0
    sumwy=0.d0
    sumwz=0.d0

    coord(1:3,1:numat)=xyz(1:3,1:numat)*bohr

    do i=1,numat
       sumw=sumw+atmass(i)
       sumwx=sumwx+atmass(i)*coord(1,i)
       sumwy=sumwy+atmass(i)*coord(2,i)
       sumwz=sumwz+atmass(i)*coord(3,i)
    enddo

    sumwx=sumwx/sumw
    sumwy=sumwy/sumw
    sumwz=sumwz/sumw
    f=1.0d0/bohr
    do i=1,numat
       x(i)=coord(1,i)-sumwx
       y(i)=coord(2,i)-sumwy
       z(i)=coord(3,i)-sumwz
    enddo

    !************************************************************************
    !*    matrix for moments of inertia is of form
    !*
    !*           |   y**2+z**2                         |
    !*           |    -y*x       z**2+x**2             | -i =0
    !*           |    -z*x        -z*y       x**2+y**2 |
    !*
    !************************************************************************
    do i=1,6
       t(i)=dble(i)*1.0d-10
    enddo
    do i=1,numat
       t(1)=t(1)+atmass(i)*(y(i)**2+z(i)**2)
       t(2)=t(2)-atmass(i)*x(i)*y(i)
       t(3)=t(3)+atmass(i)*(z(i)**2+x(i)**2)
       t(4)=t(4)-atmass(i)*z(i)*x(i)
       t(5)=t(5)-atmass(i)*y(i)*z(i)
       t(6)=t(6)+atmass(i)*(x(i)**2+y(i)**2)
    enddo
    call rsp(t,3,3,eig,evec)
    !     xsum=evec(1,1)*(evec(2,2)*evec(3,3)-evec(3,2)*evec(2,3)) +
    !    1     evec(1,2)*(evec(2,3)*evec(3,1)-evec(2,1)*evec(3,3)) +
    !    2     evec(1,3)*(evec(2,1)*evec(3,2)-evec(2,2)*evec(3,1))
    do i=1,3
       if(eig(i).lt.3.d-4) then
          eig(i)=0.d0
          rot(i)=0.d0
       else
          rot(i)=2.9979245d+4*const2/eig(i)
       endif
       xyzmom(i)=eig(i)*const1
    enddo

    aa=rot(3)/2.9979245d+4
    bb=rot(2)/2.9979245d+4
    cc=rot(1)/2.9979245d+4
    avmom=1.d-47*(xyzmom(1)+xyzmom(2)+xyzmom(3))/3.

    return
  end subroutine axisvec

!*****************************************************************
! transform to CMA for output on trj xyz file
! input coords remain unchanged
!*****************************************************************

   subroutine axis3(mode,numat,nat,coord,coordout,eig)
      use xtb_splitparam
      implicit none
      integer, intent(in)  :: numat,nat(numat),mode
      real(wp),intent(in)  :: coord(3,numat)
      real(wp),intent(out) :: coordout(3,numat),eig(3)

      real(wp) :: sumw,sumwx,sumwy,sumwz,xsum,eps,ams
      real(wp) :: t(6), evec(3,3)
      real(wp) :: x(numat),y(numat),z(numat),coordtmp(3,numat)
      real(wp) :: coord1(3,numat)
      integer  :: i,j,k
      t = (/(0.0_wp,i=1,6)/)

      sumw=1.e-20_wp
      sumwx=0.0_wp
      sumwy=0.0_wp
      sumwz=0.0_wp
      do i=1,numat
         if(mode.eq.0)then
            ams = atmass(i)
         else
            ams = 1._wp/atmass(i)
         endif
         sumw=sumw+ams
         sumwx=sumwx+ams*coord(1,i)
         sumwy=sumwy+ams*coord(2,i)
         sumwz=sumwz+ams*coord(3,i)
      enddo

      eps=1.e-3_wp
      sumwx=sumwx/sumw
      sumwy=sumwy/sumw
      sumwz=sumwz/sumw
      do i=1,numat
         x(i)=coord(1,i)-sumwx
         y(i)=coord(2,i)-sumwy
         z(i)=coord(3,i)-sumwz
         coordtmp(1,i)=x(i)
         coordtmp(2,i)=y(i)
         coordtmp(3,i)=z(i)
      enddo

      do i=1,6
         t(i)=real(i,wp)*1.0e-10_wp
      enddo

      do i=1,numat
         t(1)=t(1)+atmass(i)*(y(i)**2+z(i)**2)+eps
         t(2)=t(2)-atmass(i)*x(i)*y(i)
         t(3)=t(3)+atmass(i)*(z(i)**2+x(i)**2)+eps
         t(4)=t(4)-atmass(i)*z(i)*x(i)
         t(5)=t(5)-atmass(i)*y(i)*z(i)
         t(6)=t(6)+atmass(i)*(x(i)**2+y(i)**2)+eps
      enddo

      call rsp(t,3,3,eig,evec)

!   now to orient the molecule so the chirality is preserved
      xsum=evec(1,1)*(evec(2,2)*evec(3,3)-evec(3,2)*evec(2,3)) + &
     &     evec(1,2)*(evec(2,3)*evec(3,1)-evec(2,1)*evec(3,3)) + &
     &     evec(1,3)*(evec(2,1)*evec(3,2)-evec(2,2)*evec(3,1))
      if( xsum .lt. 0) then
         do j=1,3
   80       evec(j,1)=-evec(j,1)
         enddo
      endif
!     call prmat(6,evec,3,3,'Rmat')

      do i=1,numat
         do j=1,3
            xsum=0.d0
            do k=1,3
               xsum=xsum+coordtmp(k,i)*evec(k,j)
            enddo
            coord1(j,i)=xsum
         enddo
      enddo

      do i=1,numat
         coordout(1,i)=coord1(1,i)
         coordout(2,i)=coord1(2,i)
         coordout(3,i)=coord1(3,i)
      enddo

      return
   end subroutine axis3

end module xtb_axis
