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

!> generate a Wigner--Seitz cell from a given structure
subroutine generate_wsc(mol,wsc,rep)
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use tbdef_wsc
   implicit none
   !> molecular structure informtion
   type(tb_molecule),intent(inout) :: mol
   !> Wigner--Seitz cell data type (might be contained in mol)
   type(tb_wsc),     intent(inout) :: wsc
   !> images of the unit cell to consider
   integer,          intent(in)    :: rep(3)
! ------------------------------------------------------------------------
!  Variables
! ------------------------------------------------------------------------
   integer  :: ii,jj,ich
   integer  :: aa,bb,cc
   integer  :: c,wc
   integer  :: minpos
   integer  :: nminpos
   real(wp) :: t(3)
   real(wp) :: mindist
   real(wp) :: nmindist
   !> overall WSC tolerance to consider atoms as WSC-images
   real(wp),parameter :: tol = 0.01_wp
   real(wp),allocatable,dimension(:,:,:) :: txyz
   real(wp),allocatable,dimension(:)     :: dist
   logical, allocatable,dimension(:)     :: trans

! ------------------------------------------------------------------------
!  allocate space for the WSC first
   call wsc%allocate(mol%n,rep,mol%lattice)

! ------------------------------------------------------------------------
! initialize
! ------------------------------------------------------------------------
   allocate( txyz(3,wsc%cells,mol%n) ); txyz = 0.0_wp
   allocate( dist(wsc%cells) );     dist = 0.0_wp
   allocate( trans(wsc%cells) );    trans = .true.
! ------------------------------------------------------------------------
! Create the Wigner-Seitz Cell (WSC)
! ------------------------------------------------------------------------
   wsc%at  = 0
   wsc%itbl= 0
!$omp parallel default(none) &
!$omp private(ii,jj,wc,c,dist,trans,t) &
!$omp shared(mol,wsc,rep,txyz) &
!$omp shared(mindist,minpos,nmindist,nminpos)
!$omp do schedule(dynamic)
   ! Each WSC of one atom consists of n atoms
   do ii=1,mol%n
      do jj=1,mol%n
         !if (ii.eq.jj) cycle
         ! find according neighbours
         c=0
         dist = 0.0_wp
         do aa=-rep(1),rep(1),1
            do bb=-rep(2),rep(2),1
               do cc=-rep(3),rep(3),1
                  if ((aa.eq.0 .and. bb.eq.0 .and. cc.eq.0).and.ii.eq.jj) cycle
                  t = [aa,bb,cc]
                  c=c+1
                  txyz(:,c,ii) = mol%xyz(:,jj) + matmul(mol%lattice,t)
                  dist(c)=norm2(mol%xyz(:,ii)-txyz(:,c,ii))
               end do
            end do
         end do
         ! get first image with same dist
         ! find minimum in dist-array and assign it to minpos = minimum position
         trans=.true.
         minpos=minloc(dist(:c),dim=1)
         ! minimum distance saved in mindist
         mindist=dist(minpos)
         trans(minpos)=.false.
         wc=1
         wsc%xyz(:,wc,jj,ii)=txyz(:,minpos,ii)
         ! get other images with same distance
         find_images : do
            nminpos=minloc(dist(:c),dim=1,mask=trans(:c))
            nmindist=dist(nminpos)
            if(abs(mindist-nmindist).lt.tol)then
               trans(nminpos)=.false.
               wc=wc+1
               wsc%xyz(:,wc,jj,ii)=txyz(:,nminpos,ii)
            else
               wsc%w(jj,ii)=1.0_wp/real(wc,wp)
               wsc%itbl(jj,ii)=wc
               wsc%at(jj,ii)=jj
               exit find_images
            end if
         end do find_images
      end do
   end do
!$omp enddo
!$omp endparallel

end subroutine generate_wsc

