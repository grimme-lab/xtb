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
! You should have received a copy of the GNU Lesser General Public Licen
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! but including H in OH

subroutine heavyrmsd(n,nall,k,l,xyz,at,rmsdval)
   use xtb_lsrmsd
   implicit none
   integer n,at(n),j,nall,k,l,nn
   real*8 xyz(3,n,nall),rmsdval

   logical oh,ohbonded

   real*8 ,allocatable :: xyz1(:,:), xyz2(:,:)
   !Dummys:
   real*8  g(3,3), U(3,3), x_center(3), y_center(3)
   real*8, parameter ::bohr=0.52917726d0
   integer i

   nn=0
   do j=1,n
      oh=ohbonded(n,j,xyz(1,1,k),at)
      if(at(j).gt.2.or.oh) nn=nn+1
   enddo
   allocate(xyz1(3,nn),xyz2(3,nn))

   i=0
   do j=1,n
      oh=ohbonded(n,j,xyz(1,1,k),at)
      if(at(j).gt.2.or.oh) then
         i=i+1
         xyz1(1:3,i)=xyz(1:3,j,k)*bohr
         xyz2(1:3,i)=xyz(1:3,j,l)*bohr
      endif
   enddo

   call rmsd(i,xyz1,xyz2,0,U,x_center,y_center,rmsdval,.false.,g)

   deallocate(xyz1,xyz2)
end subroutine

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine heavyrmsdfile(k,l,rmsdval)
   use xtb_lsrmsd
   implicit none
   integer k,l
   real*8 rmsdval
   character*80 atmp1,atmp2,atmp
   !Dummys:
   real*8  g(3,3), U(3,3), x_center(3), y_center(3)
   real*8, parameter ::bohr=0.52917726d0
   integer n,idum
   real*8 ,allocatable :: xyz1(:,:), xyz2(:,:)
   integer,allocatable :: at(:)

   !call getname1(k,atmp1)
   write(atmp1,'(''scoord.'',i0)')k
   !call getname1(l,atmp2)
   write(atmp2,'(''scoord.'',i0)')l
   call rdatomnumber(atmp1,n)
   allocate(xyz1(3,n),xyz2(3,n),at(n))

   !print*,trim(atmp1)
   call rdcoord(trim(atmp1),n,xyz1,at)
   !print*,trim(atmp2)
   call rdcoord(trim(atmp2),n,xyz2,at)
   xyz1 = xyz1 * bohr
   xyz2 = xyz2 * bohr

   call rmsd(n,xyz1,xyz2,0,U,x_center,y_center,rmsdval,.false.,g)

   deallocate(xyz1,xyz2)
end subroutine

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

logical function ohbonded(n,m,xyz,at)
   implicit none
   integer n,at(n),m,i
   real*8 xyz(3,n),r

   ohbonded=.false.
   if(at(m).ne.1) return

   do i=1,n
      if(i.eq.m.or.at(i).ne.8) cycle
      r  =sqrt((xyz(1,i)-xyz(1,m))**2&
         &           +(xyz(2,i)-xyz(2,m))**2&
         &           +(xyz(3,i)-xyz(3,m))**2)
      if(r*0.52917726d0.lt.1.1) ohbonded=.true.
   enddo

end function
