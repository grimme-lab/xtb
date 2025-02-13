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

!> generate a Wigner--Seitz cell from a given structure
subroutine generate_wsc(mol,wsc)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule
   use xtb_type_wsc
   implicit none
   !> molecular structure informtion
   type(TMolecule),intent(inout) :: mol
   !> Wigner--Seitz cell data type (might be contained in mol)
   type(tb_wsc),     intent(inout) :: wsc
! ------------------------------------------------------------------------
!  Variables
! ------------------------------------------------------------------------
   integer  :: rep(3)
   integer  :: ii,jj
   integer  :: aa,bb,cc
   integer  :: c,wc
   integer  :: minpos
   integer  :: nminpos
   real(wp) :: t(3),rw(3)
   real(wp) :: mindist
   real(wp) :: nmindist
   !> overall WSC tolerance to consider atoms as WSC-images
   real(wp),parameter :: tol = 0.01_wp
   integer, allocatable,dimension(:,:)   :: lattr
   real(wp),allocatable,dimension(:)     :: dist
   logical, allocatable,dimension(:)     :: trans

   where(mol%pbc)
      rep = 1
   elsewhere
      rep = 0
   endwhere

! ------------------------------------------------------------------------
!  allocate space for the WSC first
   call wsc%allocate(mol%n,rep)

! ------------------------------------------------------------------------
! Create the Wigner-Seitz Cell (WSC)
! ------------------------------------------------------------------------
   wsc%at  = 0
   wsc%itbl= 0

   !$omp parallel default(none) &
   !$omp private(ii,jj,aa,bb,cc,wc,c,dist,trans,t,lattr,rw) &
   !$omp shared(mol,wsc,rep) &
   !$omp private(mindist,minpos,nmindist,nminpos)

   allocate( lattr(3,wsc%cells),      source = 0 )
   allocate( dist(wsc%cells),         source = 0.0_wp )
   allocate( trans(wsc%cells),        source = .true. )

   ! Each WSC of one atom consists of n atoms
   !$omp do collapse(2) schedule(dynamic,32)
   do ii=1,mol%n
      do jj=1,mol%n
         !if (ii.eq.jj) cycle
         ! find according neighbours
         c=0
         dist = 0.0_wp
         lattr = 0
         do aa=-rep(1),rep(1),1
            do bb=-rep(2),rep(2),1
               do cc=-rep(3),rep(3),1
                  if ((aa.eq.0 .and. bb.eq.0 .and. cc.eq.0).and.ii.eq.jj) cycle
                  t = [aa,bb,cc]
                  c=c+1
                  lattr(:,c) = [aa,bb,cc]
                  rw = mol%xyz(:,jj) + matmul(mol%lattice,t)
                  dist(c)=sqrt(sum((mol%xyz(:,ii)-rw)**2))
               end do
            end do
         end do
         ! sanity check; otherwise code below crashes sometimes
         if (c .eq. 0) cycle
         ! get first image with same dist
         ! find minimum in dist-array and assign it to minpos = minimum position
         trans=.true.
         minpos=minloc(dist(:c),dim=1)
         ! minimum distance saved in mindist
         mindist=dist(minpos)
         trans(minpos)=.false.
         wc=1
         wsc%lattr(:,wc,jj,ii)=lattr(:,minpos)
         ! get other images with same distance
         if (c > 1) then
            find_images : do
               nminpos=minloc(dist(:c),dim=1,mask=trans(:c))
               nmindist=dist(nminpos)
               if(abs(mindist-nmindist).lt.tol)then
                  trans(nminpos)=.false.
                  wc=wc+1
                  wsc%lattr(:,wc,jj,ii)=lattr(:,nminpos)
               else
                  wsc%w(jj,ii)=1.0_wp/real(wc,wp)
                  wsc%itbl(jj,ii)=wc
                  wsc%at(jj,ii)=jj
                  exit find_images
               end if
            end do find_images
         else
            wsc%w(jj,ii) = 1.0_wp
            wsc%itbl(jj,ii) = 1
            wsc%at(jj,ii)=jj
         endif
      end do
   end do
   !$omp end do
   deallocate(lattr, dist, trans)
   !$omp end parallel

end subroutine generate_wsc
