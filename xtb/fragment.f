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

      subroutine mrec(molcount,xyz,cn,bond,nat,at,molvec)
! molcount: number of total fragments (increased during search)
! xyz: overall Cart. coordinates
! nat: overall number of atoms
! at: atomic number array
! molvec: assignment vector of atom to fragment

      use iso_fortran_env, wp => real64
      implicit none
      real(wp),intent(in)    :: xyz(3,nat),cn(nat)
      integer, intent(in)    :: nat,at(nat)
      integer, intent(inout) :: molvec(nat),molcount
      real(wp),intent(in)    :: bond(nat,nat)
      logical, allocatable   :: taken(:)
      integer :: i
      allocate( taken(nat) )
      molvec=0
      molcount=1
      taken=.false.
      do i=1,nat
       if(.not.taken(i)) then
         molvec(i)=molcount
         taken(i)=.true.
         call neighbours(i,xyz,cn,at,taken,nat,bond,molvec,molcount)
         molcount=molcount+1
      endif
      enddo
      molcount=molcount-1
      end subroutine mrec

      recursive subroutine neighbours(i,xyz,cn,iat,taken,nat,bond,
     &                                molvec,molcnt)
      use iso_fortran_env, wp => real64
      implicit none
      real(wp) xyz(3,nat),cn(nat)
      real(wp)  bond(nat,nat)
      integer i,nat, molcnt,molvec(nat),j,iat(nat),icn,k
      logical taken(nat)
      
      icn=nint(cn(i))
      do k=1,icn
         j=maxloc(bond(:,i),1)
         bond(j,i)=0
         if (i .eq. j) cycle
         if (.not.taken(j)) then
            molvec(j)=molcnt
            taken(j)=.true.
            call neighbours(j,xyz,cn,iat,taken,nat,bond,molvec,molcnt)
         endif
      enddo
      end subroutine neighbours

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine cutcov(n,at,xyz,cn,wb,bond)
      use iso_fortran_env, wp => real64
      implicit none
      integer i,n,j,at(n),k,nb(20,n),iring,c(8,n),s(n)
      real(wp) bond(n,n)
      real(wp) wb(n,n),cn(n),xyz(3,n)
      logical samering

      bond = 0

      call neighborh(n,at,xyz,nb)
      do i=1,n
         call getring(n,nb,i,c(1,i),s(i))
      enddo

      do i=1,n
         do j=1,n
            if(wb(j,i).lt.0.5)  cycle                        
            if(wb(j,i).gt.1.3)               bond(j,i)=1
            if(cn(i).lt.1.2.or.cn(j).lt.1.2) bond(j,i)=1
            if(samering(n,i,j,c,s))          bond(j,i)=1
         enddo
      enddo

      end subroutine cutcov

