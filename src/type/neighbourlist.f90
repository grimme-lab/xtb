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

!> Static neighbour list implementation
!>
!> Based on the implementation used in DFTB+ (see module dftbp_periodic there)
module xtb_type_neighbourlist
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_search, only : bisectSearch
   use xtb_mctc_sort, only : indexHeapSort
   use xtb_mctc_thresholds, only : tolSameDist2, minNeighDist2
   use xtb_type_environment, only : TEnvironment
   implicit none
   private

   public :: TNeighbourList, init
   public :: initialNumberOfNeighbours, resizeNeigh, resizeImage


   !> Static neighbour list
   type :: TNeighbourList

      !> Maximum number of neighbours for each atom
      integer, allocatable :: neighs(:)

      !> Mapping from neighbours to atoms in the central cell
      integer, allocatable :: image(:)

      !> Index of the translation vector generating the neighbour
      integer, allocatable :: trans(:)

      !> Coordinates of the atoms and their images
      real(wp), allocatable :: coords(:, :)

      !> Neighbours of each atom
      integer, allocatable :: iNeigh(:, :)

      !> Quadratic distances to all neighbours
      real(wp), allocatable :: dist2(:, :)

      !> Weights of the neighbour contributions
      real(wp), allocatable :: weight(:, :)

      !> Maximum cutoff supported by this neighbour list
      real(wp) :: cutoff

   contains

      !> Generate neighbour list
      procedure :: generate

      !> Update neighbour list
      procedure :: update

      !> Get number of neighbours
      generic :: getNeighs => getNeighsAll, getNeighsOne

      !> Get number of all neighbours
      procedure, private :: getNeighsAll

      !> Get number of neighbours for one atom
      procedure, private :: getNeighsOne

      !> Sort the neighbour list
      procedure :: sort

   end type TNeighbourList


   !> Initialize neighbour list
   interface init
      module procedure :: initNeighbourList
   end interface init


   !> Initial number of neighbours
   integer, parameter :: initialNumberOfNeighbours = 40


contains


!> Initialize neighbour list
subroutine initNeighbourList(self, nAtom)

   !> Instance of the neighbour list
   type(TNeighbourList), intent(out) :: self

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

end subroutine initNeighbourList


!> Generate neighbour list
subroutine generate(self, env, coords0, cutoff, latPoint, symmetric)

   !> Instance of the neighbour list
   class(TNeighbourList), intent(inout) :: self

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

   call generateSequential(self, env, coords0, cutoff, latPoint, symmetric)

   call self%sort

end subroutine generate


!> Generate neighbour list in serial
subroutine generateSequential(self, env, coords0, cutoff, latPoint, symmetric)

   !> Source for generating errors
   character(len=*), parameter :: source = 'type_neighbourlist_generateSerial'

   !> Instance of the neighbour list
   class(TNeighbourList), intent(inout) :: self

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

   integer :: iTr, iAt, iAtOld, jAt, jAtEnd
   real(wp) :: ri(3), dist2, trans(3)

   !> Get dimensions for the generator loop
   nAtom = size(self%neighs, dim=1)
   nTrans = size(latPoint, dim=2)

   !> Get initial capacities
   mImage = size(self%coords, dim=2)
   mNeigh = ubound(self%iNeigh, dim=1)

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
   end do

   nImage = 0
   lpTrans: do iTr = 1, nTrans
      trans(:) = latPoint(:, iTr)
      iAtOld = 0
      lpIAt: do iAt = 1, nAtom
         !> Generate position of image
         ri(:) = coords0(:, iAt) + trans(:)
         if (symmetric) then
            jAtEnd = nAtom
         else
            jAtEnd = iAt
         end if
         lpJAt: do jAt = 1, jAtEnd
            !> Calculate distance between atoms in cell and the image
            dist2 = sum((coords0(:, jAt) - ri(:))**2)
            if (dist2 > cutoff2) then
               cycle lpJAt
            end if

            !> New atom in range found, add image to list
            if (iAt /= iAtOld) then
               nImage = nImage + 1
               if (nImage > mImage) then
                  mImage = 2*mImage
                  call resizeImage(mImage, self%image, self%trans, self%coords)
               end if
               self%coords(:, nImage) = ri(:)
               self%image(nImage) = iAt
               self%trans(nImage) = iTr
               iAtOld = iAt
            end if

            !> Check for particular close distance
            if (dist2 < minNeighDist2) then
               !> Small distances between the atom and itself are okay
               if (iTr == 1 .and. iAt == jAt) then
                  cycle lpJAt
               else
                  call env%error("Very close distance found", source)
               end if
            end if

            self%neighs(jAt) = self%neighs(jAt) + 1
            if (self%neighs(jAt) > mNeigh) then
               mNeigh = 2*mNeigh
               call resizeNeigh(mNeigh, self%iNeigh, self%dist2, self%weight)
            end if
            self%iNeigh(self%neighs(jAt), jAt) = nImage
            self%dist2(self%neighs(jAt), jAt) = dist2
            self%weight(self%neighs(jAt), jAt) = 1.0_wp

         end do lpJAt
      end do lpIAt
   end do lpTrans

   call resizeImage(nImage, self%image, self%trans, self%coords)

end subroutine generateSequential


!> Update neighbour list with new coordinates
subroutine update(self, coords0, cellVec)

   !> Initialized instance on exit
   class(TNeighbourList), intent(inout) :: self

   !> Absolute coordinates of the shifted super cells
   real(wp), intent(in) :: cellVec(:, :)

   !> Coordinates of atoms (folded into unit cell, if periodic)
   real(wp), intent(in) :: coords0(:, :)

   !> Number of atoms in the system
   integer :: nAtom

   !> Maximum number of atom without reallocation
   integer :: nImage

   integer :: iAt, iAtf, jAt, iTr, iNeigh
   real(wp) :: vec(3)

   nAtom = size(self%neighs, dim=1)
   nImage = size(self%coords, dim=2)

   !> Update the coordinates of all images
   do iAt = 1, nImage
      iAtf = self%image(iAt)
      iTr = self%trans(iAt)
      self%coords(:, iAt) = coords0(:, iAtf) + cellVec(:, iTr)
   end do

   !> Update neighbour list for all atoms
   do iAt = 1, nAtom
      do iNeigh = 1, self%neighs(iAt)
         jAt = self%iNeigh(iNeigh, iAt)
         vec(:) = self%coords(:, iAt) - self%coords(:, jAt)
         self%dist2(iNeigh, iAt) = sum(vec**2)
      end do
   end do

   !> Sort neighbour lists
   call self%sort

end subroutine update


!> Sort neighbours for all atom by distance
subroutine sort(self)

   !> Instance of the neighbour list
   class(TNeighbourList), intent(inout) :: self

   integer, allocatable :: indx(:)
   integer :: nAtom, mNeigh, iAt, nn

   nAtom = size(self%neighs, dim=1)
   mNeigh = ubound(self%iNeigh, dim=1)

   allocate(indx(mNeigh))

   do iAt = 1, nAtom
      nn = self%neighs(iAt)
      call indexHeapSort(indx(1:nn), self%dist2(1:nn, iAt), tolSameDist2)
      self%iNeigh(1:nn, iAt) = self%iNeigh(indx(:nn), iAt)
      self%dist2(1:nn, iAt) = self%dist2(indx(:nn), iAt)
      self%weight(1:nn, iAt) = self%weight(indx(:nn), iAt)
   end do

end subroutine sort


!> Get all neighbours for a given cutoff
pure subroutine getNeighsAll(self, neighs, cutoff)

   !> Instance of the neighbour list
   class(TNeighbourList), intent(in) :: self

   !> Realspace cutoff radius
   real(wp), intent(in), optional :: cutoff

   !> Number of neighbour.
   integer, intent(out) :: neighs(:)

   integer :: iatom
   real(wp) :: cutoff0

   if (present(cutoff)) then
      cutoff0 = min(cutoff, self%cutoff)
   else
      cutoff0 = self%cutoff
   end if

   do concurrent(iatom = 1:size(neighs))
      call self%getNeighs(neighs(iatom), iatom, cutoff0)
   end do

end subroutine getNeighsAll


!> Get neighbour for a given cutoff
elemental subroutine getNeighsOne(self, neighs, iatom, cutoff)

   !> Instance of the neighbour list
   class(TNeighbourList), intent(in) :: self

   !> Atom of interest.
   integer, intent(in) :: iatom

   !> Realspace cutoff radius
   real(wp), intent(in), optional :: cutoff

   !> Number of neighbour.
   integer, intent(out) :: neighs

   real(wp) :: cutoff2

   if (present(cutoff)) then
      cutoff2 = min(cutoff, self%cutoff)**2
   else
      cutoff2 = self%cutoff**2
   endif

   call bisectSearch(neighs, self%dist2(1:self%neighs(iatom), iatom), &
      & cutoff2, tolSameDist2)

end subroutine getNeighsOne


!> Reallocate arrays which depends on the maximal nr. of all atoms.
subroutine resizeImage(mImage, image, trans, coords)

   !> Maximum number of new atoms
   integer, intent(in) :: mImage

   !> Array mapping images of atoms to originals in the central cell
   integer, allocatable, intent(inout) :: image(:)

   !> Index of unit cell containing atom
   integer, allocatable, intent(inout) :: trans(:)

   !> Coordinates of all atoms (actual and image)
   real(wp), allocatable, intent(inout) :: coords(:, :)

   integer :: nImage
   integer, allocatable :: iTmp(:)
   real(wp), allocatable :: rTmp(:, :)

   nImage = min(size(image), mImage)

   call move_alloc(image, iTmp)
   allocate(image(mImage))
   image(:) = 0
   image(:nImage) = iTmp(:nImage)

   iTmp(:) = trans(:)
   deallocate(trans)
   allocate(trans(mImage))
   trans(:nImage) = iTmp(:nImage)

   call move_alloc(coords, rTmp)
   allocate(coords(3, mImage))
   coords(:, :nImage) = rTmp(:, :nImage)

end subroutine resizeImage


!> Reallocate array which depends on the maximal nr. of neighbours.
subroutine resizeNeigh(mNeigh, iNeigh, dist2, weight)

   !> Maximum number of new atoms
   integer, intent(in) :: mNeigh

   !> List of neighbours
   integer, allocatable, intent(inout) :: iNeigh(:, :)

   !> Square of distances between atoms
   real(wp), allocatable, intent(inout) :: dist2(:,:)

   !> Weights for all neighbours
   real(wp), allocatable, intent(inout) :: weight(:,:)

   integer :: nNeigh, nAtom
   integer, allocatable :: iTmp(:,:)
   real(wp), allocatable :: rTmp(:,:)

   nNeigh = ubound(iNeigh, dim=1)
   nAtom = size(iNeigh, dim=2)

   call move_alloc(iNeigh, iTmp)
   allocate(iNeigh(0:mNeigh, nAtom))
   iNeigh(:,:) = 0
   iNeigh(:nNeigh, :nAtom) = iTmp

   call move_alloc(dist2, rTmp)
   allocate(dist2(0:mNeigh, nAtom))
   dist2(:,:) = 0.0_wp
   dist2(:nNeigh, :nAtom) = rTmp

   rTmp(:,:) = weight
   deallocate(weight)
   allocate(weight(0:mNeigh, nAtom))
   weight(:,:) = 0.0_wp
   weight(:nNeigh, :nAtom) = rTmp

end subroutine resizeNeigh


end module xtb_type_neighbourlist
