! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> Wigner-Seitz cell based on a neighbour list
module xtb_type_wignerseitzcell
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_search, only : bisectSearch
   use xtb_mctc_sort, only : indexHeapSort
   use xtb_mctc_thresholds, only : tolSameDist2, minNeighDist2
   use xtb_type_environment, only : TEnvironment
   use xtb_type_neighbourlist, only : TNeighbourList, initialNumberOfNeighbours, &
      & resizeNeigh, resizeImage
   implicit none
   private

   public :: TWignerSeitzCell, init


   !> Wigner-Seitz cell based on the infrastructure for the neighbour list
   type, extends(TNeighbourList) :: TWignerSeitzCell
   contains

      !> Generate Wigner-Seitz cell
      procedure :: generate

   end type TWignerSeitzCell


   !> Initialize Wigner-Seitz cell
   interface init
      module procedure :: initWignerSeitzCell
   end interface init


   !> Tolerance for including images in the Wigner-Seitz cell
   real(wp), parameter :: wsTol2 = 0.01_wp**2


contains


!> Initialize Wigner-Seitz cell
subroutine initWignerSeitzCell(self, nAtom)

   !> Instance of the Wigner-Seitz cell
   type(TWignerSeitzCell), intent(out) :: self

   !> Number of atoms in the system
   integer, intent(in) :: nAtom

   allocate(self%neighs(nAtom))
   allocate(self%image(nAtom))
   allocate(self%coords(3, nAtom))
   allocate(self%trans(nAtom))

   allocate(self%iNeigh(0:initialNumberOfNeighbours, nAtom))
   allocate(self%dist2(0:initialNumberOfNeighbours, nAtom))
   allocate(self%weight(0:initialNumberOfNeighbours, nAtom))

   self%cutoff = -1.0_wp

end subroutine initWignerSeitzCell


!> Generate Wigner-Seitz cell
subroutine generate(self, env, coords0, cutoff, latPoint, symmetric)

   !> Source for generating errors
   character(len=*), parameter :: source = 'type_wignerseitzcell_generate'

   !> Instance of the Wigner-Seitz cell
   class(TWignerSeitzCell), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Folded coordinates of the atoms inside the unit cell
   real(wp), intent(in) :: coords0(:, :)

   !> Realspace cutoff for the generation of neighbours
   real(wp), intent(in) :: cutoff

   !> Lattice points for all supercells within the cutoff
   real(wp), intent(in) :: latPoint(:, :)

   !> Generate a symmetric map
   logical, intent(in) :: symmetric

   !> Number of atoms in the central cell
   integer :: nAtom

   !> Number of translations to consider
   integer :: nTrans

   !> Capacity of the atom related arrays
   integer :: mImage

   !> Capacity of the neighbour related arrays
   integer :: mNeigh

   !> Index of the current image
   integer :: nImage

   !> Quadratic cutoff
   real(wp) :: cutoff2

   integer :: iTr, iAt, iAtOld, jAt, jAtEnd, iWs, nWs
   real(wp) :: ri(3), dist2, trans(3), wsWeight
   integer, allocatable :: wsTrans(:), indx(:)
   real(wp), allocatable :: wsDist2(:)

   !> Get dimensions for the generator loop
   nAtom = size(self%neighs, dim=1)
   nTrans = size(latPoint, dim=2)
   nImage = nAtom*nTrans

   !> Get initial capacities
   mImage = size(self%coords, dim=2)
   mNeigh = ubound(self%iNeigh, dim=1)

   if (nImage > mImage) then
      call resizeImage(nImage, self%image, self%trans, self%coords)
   end if

   !> Set the value of the cutoff
   self%cutoff = cutoff
   cutoff2 = cutoff**2

   !> Initialize neighbour list with distances from every atom to itself
   self%neighs(:) = 0
   self%iNeigh(:, :) = 0
   self%dist2(:, :) = 0.0_wp
   self%weight(:, :) = 0.0_wp
   do iAt = 1, nAtom
      self%iNeigh(0, iAt) = iAt
      self%weight(0, iAt) = 1.0_wp
      self%coords(:, iAt) = coords0(:, iAt)
      self%image(iAt) = iAt
      self%trans(iAt) = 1
   end do

   !> Allocate scratch space for Wigner-Seitz neighbours
   allocate(wsDist2(nTrans))
   allocate(wsTrans(nTrans))
   allocate(indx(nTrans))

   !> Start with images in the central cell already generated
   nImage = nAtom
   lpIAt: do iAt = 1, nAtom
      if (symmetric) then
         jAtEnd = nAtom
      else
         jAtEnd = iAt
      end if
      lpJAt: do jAt = 1, jAtEnd

         nWs = 0
         lpTrans: do iTr = 1, nTrans
            !> Generate position of image
            ri(:) = coords0(:, iAt) + latPoint(:, iTr)
            dist2 = sum((ri - coords0(:, jAt))**2)
            if (dist2 > cutoff2) then
               cycle lpTrans
            end if

            if (dist2 < minNeighDist2) then
               !> Small distances between the atom and itself are okay
               if (iTr == 1 .and. iAt == jAt) then
                  cycle lpTrans
               else
                  call env%error("Very close distance found", source)
               end if
            end if

            nWs = nWs + 1
            wsDist2(nWs) = dist2
            wsTrans(nWs) = iTr
         end do lpTrans

         if (nWs == 0) cycle lpJAt

         !> Sort the neighbours in the Wigner-Seitz cell
         call indexHeapSort(indx(1:nWs), wsDist2(1:nWs), tolSameDist2)
         wsDist2(1:nWs) = wsDist2(indx(1:nWs))
         wsTrans(1:nWs) = wsTrans(indx(1:nWs))

         !> Drop all but the closest neighbours
         call bisectSearch(nWs, wsDist2(1:nWs), wsDist2(1)+wsTol2, tolSameDist2)
         wsWeight = 1.0_wp/real(nWs, wp)

         !> Make sure the capacity of the arrays is large enough
         if (nImage + nWs > mImage) then
            mImage = max(2*mImage, mImage + 2*nWs)
            call resizeImage(mImage, self%image, self%trans, self%coords)
         end if
         if (self%neighs(jAt) + nWs > mNeigh) then
            mNeigh = max(2*mNeigh, mNeigh + 2*nWs)
            call resizeNeigh(mNeigh, self%iNeigh, self%dist2, self%weight)
         end if

         !> Store the Wigner-Seitz cell
         do iWs = 1, nWs
            iTr = wsTrans(iWs)
            if (iTr /= 1) then
               nImage = nImage + 1
               self%coords(:, nImage) = coords0(:, iAt) + latPoint(:, wsTrans(iWs))
               self%image(nImage) = iAt
               self%trans(nImage) = iTr
               self%iNeigh(self%neighs(jAt)+iWs, jAt) = nImage
            else
               self%iNeigh(self%neighs(jAt)+iWs, jAt) = iAt
            end if
            self%dist2(self%neighs(jAt)+iWs, jAt) = wsDist2(iWs)
            self%weight(self%neighs(jAt)+iWs, jAt) = wsWeight
         end do
         self%neighs(jAt) = self%neighs(jAt) + nWs

      end do lpJAt
   end do lpIAt

   call resizeImage(nImage, self%image, self%trans, self%coords)

   call self%sort

end subroutine generate


end module xtb_type_wignerseitzcell
