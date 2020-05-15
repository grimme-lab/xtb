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
module xtb_gfnff_mrec
      private
      public :: mrecgff
contains

      subroutine mrecgff(nat,nb,molcount,molvec)
      ! molcount: number of total fragments (increased during search)
      ! nat: overall number of atoms
      ! molvec: assignment vector of atom to fragment
      implicit none
      integer nat,molvec(nat),i,j,molcount,nb(20,nat)
      real*8, allocatable:: bond(:,:)
      logical,allocatable:: taken(:)

      allocate(taken(nat),bond(nat,nat))
      bond=0
      do i=1,nat
         do j=1,nb(20,i)
            bond(i,nb(j,i))=1
            bond(nb(j,i),i)=1
         enddo
      enddo
      molvec=0
      molcount=1
      taken=.false.
      do i=1,nat
       if(.not.taken(i)) then
         molvec(i)=molcount
         taken(i)=.true.
         call mrecgff2(nb,i,taken,nat,bond,molvec,molcount)
         molcount=molcount+1
      endif
      enddo
      molcount=molcount-1
      end subroutine mrecgff

      recursive subroutine mrecgff2(nb,i,taken,nat,bond,molvec,molcnt)
      implicit none
      integer i,nat, molcnt,molvec(nat),j,icn,k,nb(20,nat)
      real*8 bond(nat,nat)
      logical taken(nat)

      icn=nb(20,i)
      do k=1,icn
         j=maxloc(bond(:,i),1)
         bond(j,i)=0
         if (i .eq. j) cycle
         if (.not.taken(j)) then
            molvec(j)=molcnt
            taken(j)=.true.
            call mrecgff2(nb,j,taken,nat,bond,molvec,molcnt)
         endif
      enddo
      end subroutine mrecgff2
end module xtb_gfnff_mrec
