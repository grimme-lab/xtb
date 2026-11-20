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
      public :: mrecgff, mrecgffPBC, canonical_charge_fragment
contains

  !> Select a unique largest fragment without using its label or atom indices.
  !> Equal sizes are ranked lexicographically by element counts, high Z first.
  !> Return zero for an unresolved tie; charge localization then needs input.
  integer function canonical_charge_fragment(nat,at,nfrag,fraglist) result(selected)
    integer,intent(in) :: nat,nfrag,at(nat),fraglist(nat)
    integer,allocatable :: counts(:,:)
    integer :: i,f,z,comparison,best,maxz
    logical :: tied

    selected=0
    if(nat<=0.or.nfrag<=0) return
    if(any(fraglist<1).or.any(fraglist>nfrag)) return
    maxz=maxval(at)
    allocate(counts(0:maxz,nfrag),source=0)
    do i=1,nat
      f=fraglist(i)
      counts(0,f)=counts(0,f)+1
      counts(at(i),f)=counts(at(i),f)+1
    enddo
    best=1
    tied=.false.
    do f=2,nfrag
      comparison=counts(0,f)-counts(0,best)
      if(comparison==0) then
        do z=maxz,1,-1
          comparison=counts(z,f)-counts(z,best)
          if(comparison/=0) exit
        enddo
      endif
      if(comparison>0) then
        best=f
        tied=.false.
      else if(comparison==0) then
        tied=.true.
      endif
    enddo
    if(.not.tied) selected=best
  end function canonical_charge_fragment



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


      subroutine mrecgffPBC(nat,numctr,numnb,nb,molcount,molvec)
      ! molcount: number of total fragments (increased during search)
      ! nat: overall number of atoms
      ! molvec: assignment vector of atom to fragment
      implicit none
      integer, intent(in) :: nat, numctr, numnb, nb(numnb,nat,numctr)
      integer, intent(inout) ::molvec(nat),molcount
      integer :: i,j, iTr
      real*8, allocatable:: bond(:,:,:)
      logical,allocatable:: taken(:)

      allocate(taken(nat),bond(nat,nat,numctr))
      ! create array with 1's as marks for bond partners
      bond=0
      do i=1,nat
        do iTr=1, numctr  
          do j=1,nb(numnb,i,iTr)
            bond(nb(j,i,iTr), i, iTr)=1
          enddo
        enddo
      enddo
     
      ! double-check
      if (int(sum(bond)).ne.sum(nb(numnb,:,:))) then
        write(*,*) 
        write(*,*) 'Warning: Check mrec.f90',int(sum(bond)),sum(nb(numnb,:,:))
        write(*,*)
      endif

      ! recursive search for fragments 
      molvec=0 
      molcount=1
      taken=.false.
      do i=1,nat
          if(.not.taken(i)) then
            molvec(i)=molcount
            taken(i)=.true.
            call mrecgff2PBC(numctr,numnb,nat,nb,i,taken,bond,molvec,molcount)
            molcount=molcount+1
          endif
      enddo
      molcount=molcount-1
      end subroutine mrecgffPBC

      recursive subroutine mrecgff2PBC(numctr,numnb,nat,nb,i,taken,bond,molvec,molcnt)
      implicit none
      integer, intent(in) :: numctr, numnb, nat, nb(numnb,nat,numctr),i, molcnt
      integer j,icn,k,iTr, j_iTr(2)
      real*8, intent(inout) :: bond(nat,nat,numctr)
      integer, intent(inout) :: molvec(nat)
      logical, intent(inout) :: taken(nat)

      icn=sum(nb(numnb,i,:))
      do k=1,icn
         j_iTr=maxloc(bond(:,i,:))  ! get positions of the 1's one after another
         j  = j_iTr(1)
         iTr= j_iTr(2)
         bond(j,i,iTr)=0            ! set to zero to get to next 1-entry
         if (i .eq. j.and.iTr.eq.1) cycle
         if (.not.taken(j)) then
            molvec(j)=molcnt
            taken(j)=.true.
            call mrecgff2PBC(numctr,numnb,nat,nb,j,taken,bond,molvec,molcnt)
         endif
      enddo
      end subroutine mrecgff2PBC
end module xtb_gfnff_mrec
