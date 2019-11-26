! This file is part of xtb.
! Based on the DFTB+ implementation: https://github.com/dftbplus/dftbplus
!
! Copyright (C) 2006-2019 DFTB+ developers group
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

!> Implements neighbour list generator with static neighbour list.
!
!  The implementation is adapted from DFTB+ and is a merge of the dynamic
!  neighbour list (here the generator) and the neighbour list from the
!  periodic module. Originally the dynamic neighbour list iterator was also
!  implemented but it turned out to not parallize well together with OpenMP
!  and was discarded in favor of a static neighbourlist.
!
!  Other changes include mainly naming style and indentation width, as well
!  as the perference of type bound procedures over public procedures for
!  better namespacing.
module tbdef_neighbourlist
   use iso_fortran_env, only: wp => real64
   use tbdef_latp_generator
   implicit none
   private

   public :: tb_neighbourlist
   public :: tb_neighlist_generator

   integer, parameter :: iter_chunk_size = 512

   !> Tolerance for atomic distances. Atoms closer than that are regarded to sit
   !  on the same positions. (Dummy atoms)
   real(wp), parameter :: tolSameDist = 1.0e-5_wp
   !> Tolerance for atomic square distances
   real(wp), parameter :: tolSameDist2 = tolSameDist**2

   !> Minimal distance between neighbours. (Neighbour distances smaller than that
   !  are meaningless because the parametrisation usually do not cover this
   !  region.)
   real(wp), parameter :: minNeighDist = 1.0e-2_wp
   !> Minimal square distance between neighbours
   real(wp), parameter :: minNeighDist2 = minNeighDist**2


   !> Dynamic neighbour list
   type :: tb_neighlist_generator
      private
      !> Cutoff for neighbour generation
      real(wp) :: cutoff
      !> Nr. of atoms
      integer :: nAtom
      !> Coordinates of atoms (folded into unit cell, if periodic)
      real(wp), allocatable :: coords0(:, :)
      !> Whether system is periodic
      logical :: periodic
      !> Lattice vectors, if system is periodic
      real(wp), allocatable :: lattice(:, :)
      !> Inverse lattice vectors, if system is periodic
      real(wp), allocatable :: inv_lat(:, :)
      !> Returns cell translation vectors in relative coordinates.
      real(wp), allocatable :: latp(:, :)
   contains
      procedure :: new => neighgen_new
      procedure :: set_lattice => neighgen_set_lattice
      procedure :: set_coords => neighgen_set_coords
      procedure :: set_cutoff => neighgen_set_cutoff
      procedure, private :: get_latp => neighgen_get_latp
   end type tb_neighlist_generator


   !> Dynamic neighbour list
   type :: tb_neighbourlist
      !> Number of atoms
      integer :: nAtom
      !> Cutoff for neighbour generation
      real(wp) :: cutoff
      !> index of neighbour atoms
      integer, allocatable :: ineigh(:,:)
      !> nr. of neighbours
      integer, allocatable :: neighs(:)
      !> temporary array for neighbour distances
      real(wp), allocatable :: dists2(:,:)
      !> Coordinates of the objects interacting with the atoms in the central cell
      real(wp), allocatable :: coords(:,:)
      !> Returns for all atoms the index of an atom in the central cell which
      !  the atom is mapped on to.
      integer, allocatable :: image(:)
      !> Returns the index of the translating superlattice vector for each atom.
      integer, allocatable :: iCellVec(:)
   contains
      generic :: new => new_default, new_from_generator
      procedure, private :: new_default => neighlist_new_default
      procedure, private :: new_from_generator => neighlist_new_from_generator
      procedure :: update => neighlist_update
      generic :: get_neighs => get_neighs_all, get_neighs_one
      procedure, private :: get_neighs_all => neighlist_get_neighs_all
      procedure, private :: get_neighs_one => neighlist_get_neighs_one
      procedure, private :: sort => neighlist_sort
      procedure, private :: resize_neigh => neighlist_resize_neigh
      procedure, private :: resize_coord => neighlist_resize_coord
   end type tb_neighbourlist


contains

!> Initializes a dynamic neighbour list.
subroutine neighlist_new_default(self, nAtom, cutoff)
   !> Initialized instance on exit.
   class(tb_neighbourlist), intent(out) :: self
   !> Number of atoms
   integer, intent(in) :: nAtom
   !> Cutoff up to which the neighbours should be generated
   real(wp), intent(in) :: cutoff

   self%nAtom = nAtom
   self%cutoff = cutoff

   allocate(self%neighs(nAtom), source=0)
   allocate(self%ineigh(0:40, nAtom), source=0)
   allocate(self%dists2(0:40, nAtom), source=0.0_wp)

   allocate(self%coords(3, nAtom), source=0.0_wp)
   allocate(self%image(nAtom), source=0)
   allocate(self%iCellVec(nAtom), source=0)

end subroutine neighlist_new_default

subroutine neighlist_new_from_generator(self, neighgen)
   !> Initialized instance on exit.
   class(tb_neighbourlist), intent(out) :: self
   !> Neighbourlist generator to initialize and fill neighbourlist
   class(tb_neighlist_generator), intent(in) :: neighgen

   integer :: nAll

   call self%new(neighgen%nAtom, neighgen%cutoff)

   call self%update(neighgen, nAll)

end subroutine neighlist_new_from_generator

!> Initializes a dynamic neighbour list.
subroutine neighgen_new(self, cutoff, nAtom, coords, lattice, periodic)
   !> Initialized instance on exit.
   class(tb_neighlist_generator), intent(out) :: self
   !> Cutoff up to which the neighbours should be generated
   real(wp), intent(in) :: cutoff
   !> Nr. of atoms in the system.
   integer, intent(in) :: nAtom
   !> New coordinates.
   real(wp), intent(in), optional :: coords(:, :)
   !> New lattice vectors.
   real(wp), intent(in), optional :: lattice(:, :)
   !> Whether the system is periodic
   logical, intent(in) :: periodic

   self%cutoff = cutoff
   self%nAtom = nAtom
   self%periodic = periodic
   allocate(self%coords0(3, self%nAtom))
   if (present(coords)) then
      call self%set_coords(coords)
   endif

   if (self%periodic) then
      allocate(self%lattice(3, 3))
      allocate(self%inv_lat(3, 3))
      if (present(lattice)) then
         call self%set_lattice(lattice)
      end if
   end if

   call self%get_latp

end subroutine neighgen_new

!> Calculates the translation vectors for cells, which could contain atoms
!  interacting with any of the atoms in the central cell.
!  This subroutine uses a simple guess to get the necessary translation vectors.
!  This results in a set of vectors wich could for very asymmetric cells a large
!  amount bigger than the real necessary one.
subroutine neighgen_get_latp(self, without_origin)
   !> Initialized instance on exit.
   class(tb_neighlist_generator), intent(inout) :: self
   !> whether to exclude the (0,0,0) point
   logical,  intent(in), optional :: without_origin

   type(tb_latp_generator) :: latpgen
   real(wp), allocatable :: latp(:, :)
   logical :: with_origin

   if (self%periodic) then
      call latpgen%new(self%lattice, self%inv_lat, self%cutoff, &
         &             pos_extension=1, neg_extension=1, exclude_origin=.true.)
      call latpgen%get_all_points(latp)

      with_origin = .true.
      if (present(without_origin)) with_origin = .not.without_origin
      if (with_origin) then
         allocate(self%latp(3, size(latp, dim=2) + 1), source=0.0_wp)
         self%latp(:, 2:) = latp
      else
         call move_alloc(latp, self%latp)
      endif
   else
      allocate(self%latp(3, 1), source=0.0_wp)
   endif

end subroutine neighgen_get_latp


!> Updates the lattice vectors.
subroutine neighgen_set_lattice(self, lattice)
   use pbc_tools
   !> Instance.
   class(tb_neighlist_generator), intent(inout) :: self
   !> New lattice vectors.
   real(wp), intent(in) :: lattice(:, :)

   if (self%periodic) then
      self%lattice = lattice
      self%inv_lat = mat_inv_3x3(lattice)
   endif

end subroutine neighgen_set_lattice

!> Updates the coordiantes.
subroutine neighgen_set_coords(self, coords)
   !> Instance.
   class(tb_neighlist_generator), intent(inout) :: self
   !> New coordinates.
   real(wp), intent(in) :: coords(:, :)

   self%coords0 = coords

end subroutine neighgen_set_coords

!> Updates the cutoff radius.
subroutine neighgen_set_cutoff(self, cutoff)
   use pbc_tools
   !> Instance.
   class(tb_neighlist_generator), intent(inout) :: self
   !> New lattice vectors.
   real(wp), intent(in) :: cutoff

   if (cutoff > 0.0_wp) then
      self%cutoff = cutoff
   end if

end subroutine neighgen_set_cutoff

!> Updates the neighbour list according a given geometry.
!  The neighbourlist for the given cutoff is calculated. Arrays are resized if
!  necessary. The neighbour list determination is a simple N^2 algorithm,
!  calculating the distance between the possible atom pairs.
subroutine neighlist_update(self, neighgen, nAllAtom, symmetric)
   !> Neighbourlist.
   class(tb_neighlist_generator), intent(in) :: neighgen
   !> Neighbourlist.
   class(tb_neighbourlist), intent(inout) :: self
   !> Returns the nr. of all objects (including those in the translated cells.)
   integer, intent(out) :: nAllAtom
   !> Optional, whether the map should be symmetric (dftb default = .false.)
   logical, intent(in), optional :: symmetric

   !> Nr. of atoms in the system
   integer :: nAtom
   !> Max. nr. of atom without reallocation
   integer :: mAtom
   !> Max. nr. of neighbours without reallocation
   integer :: maxNeighbour
   !> Nr. of cell translation vectors
   integer :: nCellVec
   !> Square of the diatomic interaction cutoffs
   real(wp) :: cutoff2
   !> Absolute coordinates of the shifted supercells which could have interacting
   !  atoms with the central cell.
   real(wp), allocatable :: rCellVec(:,:)

   real(wp) :: dist2
   real(wp) :: rCell(3), rr(3)
   integer :: ii, iAtom1, oldIAtom1, iAtom2
   integer :: nn1, iAtom2End
   logical :: symm
   integer, allocatable :: indx(:)


   allocate(rCellVec(3, size(neighgen%latp, dim=2)))
   do ii = 1, size(rCellVec, dim=2)
      rCellVec(:, ii) = matmul(neighgen%lattice, neighgen%latp(:, ii))
   end do

   nAtom = size(self%neighs, dim=1)
   mAtom = size(self%coords, dim=2)
   maxNeighbour = ubound(self%ineigh, dim=1)
   nCellVec = size(rCellVec, dim=2)

   !:ASSERT(nAtom <= mAtom)
   !:ASSERT(allocated(coords))
   !:ASSERT(size(coords, dim=1) == 3)
   !:ASSERT(allocated(image))
   !:ASSERT(size(image) == mAtom)
   !:ASSERT(allocated(iCellVec))
   !:ASSERT(size(iCellVec) == mAtom)
   !:ASSERT(size(self%ineigh, dim=2) == nAtom)
   !:ASSERT((size(coord0, dim=1) == 3) .and. size(coord0, dim=2) >= nAtom)
   !:ASSERT((size(rCellVec, dim=1) == 3))
   !:ASSERT(neighgen%cutoff >= 0.0_wp)

   symm = .false.
   if (present(symmetric)) then
      symm = symmetric
   end if
   self%cutoff = neighgen%cutoff
   cutoff2 = neighgen%cutoff**2
   nAllAtom = 0

   ! Clean arrays.
   !  (Every atom is the 0th neighbour of itself with zero distance square.)
   self%neighs = 0
   self%ineigh = 0
   do ii = 1, nAtom
      self%ineigh(0, ii) = ii
   end do
   self%dists2 = 0.0_wp

   ! Loop over all possible neighbours for all atoms in the central cell.
   ! Only those neighbours are considered which map on atom with a higher
   ! or equal index in the central cell.
   ! Outer two loops: all atoms in all cells.
   ! Inner loop: all atoms in the central cell.
   do ii = 1, nCellVec
      rCell = rCellVec(:, ii)
      oldIAtom1 = 0
      do iAtom1 = 1, nAtom
         rr = neighgen%coords0(:, iAtom1) + rCell
         if (symm) then
            iAtom2End = nAtom
         else
            iAtom2End = iAtom1
         end if
         do iAtom2 = 1, iAtom2End
            !  If distance greater than cutoff -> skip
            dist2 = sum((neighgen%coords0(:, iAtom2) - rr)**2)
            if (dist2 > cutoff2) then
               cycle
            end if
            ! New interacting atom -> append
            ! We need that before checking for interaction with dummy atom or
            ! with itself to make sure that atoms in the central cell are
            ! appended  exactly in the same order as found in the coord0 array.
            if (iAtom1 /= oldIAtom1) then
               nAllAtom = nAllAtom + 1
               if (nAllAtom > mAtom) then
                  mAtom = 2*mAtom
                  call self%resize_coord(mAtom)
               end if
               self%coords(:, nAllAtom) = rr(:)
               self%image(nAllAtom) = iAtom1
               self%iCellVec(nAllAtom) = ii
               oldIAtom1 = iAtom1
            end if

            if (dist2 < minNeighDist2) then
               if (ii == 1 .and. iAtom1 == iAtom2) then
                  ! We calculated the distance between the same atom in the unit cell
                  cycle
               end if
            end if

            self%neighs(iAtom2) = self%neighs(iAtom2) + 1
            if (self%neighs(iAtom2) > maxNeighbour) then
               maxNeighbour = 2*maxNeighbour
               call self%resize_neigh(maxNeighbour)
            end if
            self%ineigh(self%neighs(iAtom2), iAtom2) = nAllAtom
            self%dists2(self%neighs(iAtom2), iAtom2) = dist2

         end do
      end do
   end do

   call self%sort

   call self%resize_coord(nAllAtom)

end subroutine neighlist_update

!> Sort neighbours according to distance.
subroutine neighlist_sort(self)
   !> Neighbourlist.
   class(tb_neighbourlist), intent(inout) :: self
   integer :: iAtom1, nn1
   integer, allocatable :: indx(:)
   ! Sort neighbours for all atom by distance
   allocate(indx(ubound(self%ineigh, dim=1)))
   do concurrent(iAtom1 = 1:size(self%ineigh, dim=2))
      nn1 = self%neighs(iAtom1)
      call index_heap_sort(indx(1:nn1), self%dists2(1:nn1, iAtom1), tolSameDist2)
      self%ineigh(1:nn1, iAtom1) = self%ineigh(indx(:nn1), iAtom1)
      self%dists2(1:nn1, iAtom1) = self%dists2(indx(:nn1), iAtom1)
   end do
end subroutine neighlist_sort

pure subroutine neighlist_get_neighs_all(self, neighs, cutoff)
   !> Neighbourlist.
   class(tb_neighbourlist), intent(in) :: self
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
   endif

   do concurrent(iatom = 1:self%nAtom)
      call self%get_neighs(neighs(iatom), iatom, cutoff0)
   enddo

end subroutine neighlist_get_neighs_all

elemental pure subroutine neighlist_get_neighs_one(self, neighs, iatom, cutoff)
   !> Neighbourlist.
   class(tb_neighbourlist), intent(in) :: self
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

   call bisect_search(neighs, self%dists2(1:self%neighs(iatom), iatom), &
      &               cutoff2, tolSameDist2)

end subroutine neighlist_get_neighs_one

!> Real case heap sort returning an index.
!  Based on Numerical Recipes Software 1986-92
pure subroutine index_heap_sort(indx, array, tolerance)
   !> Indexing array on return
   integer, intent(out) :: indx(:)
   !> Array of values to be sorted
   real(wp), intent(in) :: array(:)
   !> Tolerance for equality of two elements
   real(wp), intent(in), optional :: tolerance

   integer :: n, ir, ij, il, ii
   integer :: indxTmp
   real(wp) :: arrayTmp, tol

   !:ASSERT(size(array)==size(indx))

   if (present(tolerance)) then
      tol = tolerance
   else
      tol = epsilon(0.0_wp)
   end if

   do ii = 1, size(indx)
      indx(ii) = ii
   end do
   n = size(array)
   if (n <= 1) return
   il = n/2 + 1
   ir = n
   do
      if (il > 1) then
         il = il - 1
         indxTmp = indx(il)
         arrayTmp = array(indxTmp)
      else
         indxTmp = indx(ir)
         arrayTmp = array(indxTmp)
         indx(ir) = indx(1)
         ir = ir - 1
         if (ir < 1) then
            indx(1) = indxTmp
            return
         end if
      end if
      ii = il
      ij = 2 * il
      do while (ij <= ir)
         if (ij < ir) then
            if (array(indx(ij)) < array(indx(ij+1)) - tol) then
               ij = ij + 1
            end if
         end if
         if(arrayTmp < array(indx(ij)) - tol) then
            indx(ii) = indx(ij)
            ii = ij
            ij = 2*ij
         else
            ij = ir + 1
         end if
      end do
      indx(ii)=indxTmp
   end do

end subroutine index_heap_sort

!> real case for bisection search
pure subroutine bisect_search(j, xx, x, tol)
   !> located element such that xx(j) < x < xx(j+1)
   integer, intent(out) :: j
   !> array of values in monotonic order to search through
   real(wp), intent(in) :: xx(:)
   !> value to locate j for
   real(wp), intent(in) :: x
   !> Tolerance for equality comparision
   real(wp), intent(in), optional :: tol

   integer :: n
   integer :: jlower, jupper, jcurr
   real(wp) :: rTol
   logical :: ascending

   n = size(xx)
   if (n == 0) then
      j = 0
      return
   end if

   if (present(tol)) then
      rTol = tol
   else
      rTol = epsilon(0.0_wp)
   end if

   if (x < xx(1) - rTol) then
      j = 0
   else if (abs(x - xx(1)) <= rTol) then
      j = 1
   else if (abs(x - xx(n)) <= rTol) then
      j = n - 1
   else if (x > xx(n) + rTol) then
      j = n
   else
      ascending = (xx(n) >=  xx(1))
      jlower = 0
      jcurr = n+1
      do while ((jcurr-jlower) > 1)
         jupper = (jcurr+jlower)/2
         if (ascending .eqv. (x >= xx(jupper) + rTol)) then
            jlower = jupper
         else
            jcurr = jupper
         end if
      end do
      j = jlower
   end if
end subroutine bisect_search

!> Reallocate array which depends on the maximal nr. of neighbours.
subroutine neighlist_resize_neigh(self, mNewNeighbour)
   !> Neighbourlist.
   class(tb_neighbourlist), intent(inout) :: self
   !> maximum number of new atoms
   integer, intent(in) :: mNewNeighbour

   integer :: mNeighbour, mAtom
   integer, allocatable :: tmpIntR2(:,:)
   real(wp), allocatable :: tmpRealR2(:,:)

   mNeighbour = ubound(self%ineigh, dim=1)
   mAtom = size(self%ineigh, dim=2)

   !:ASSERT(mNewNeighbour > 0 .and. mNewNeighbour > mNeighbour)
   !:ASSERT(all(shape(self%dists2) == shape(self%ineigh)))

   call move_alloc(self%ineigh, tmpIntR2)
   allocate(self%ineigh(0:mNewNeighbour, mAtom), source=0)
   self%ineigh(:mNeighbour, :mAtom) = tmpIntR2

   call move_alloc(self%dists2, tmpRealR2)
   allocate(self%dists2(0:mNewNeighbour, mAtom), source=0.0_wp)
   self%dists2(:mNeighbour, :mAtom) = tmpRealR2

end subroutine neighlist_resize_neigh

!> Reallocate arrays which depends on the maximal nr. of all atoms.
subroutine neighlist_resize_coord(self, mNewAtom)
   !> Neighbourlist.
   class(tb_neighbourlist), intent(inout) :: self
   !> maximum number of new atoms
   integer, intent(in) :: mNewAtom

   integer :: mAtom
   integer, allocatable :: itmp(:)
   real(wp), allocatable :: atmp(:, :)

   mAtom = size(self%image)

   !:ASSERT(size(self%iCellVec) == mAtom)
   !:ASSERT(all(shape(self%coords) == (/ 3, mAtom /)))
   !:ASSERT((mNewAtom > 0) .and. (mNewAtom > mAtom))
   !:ASSERT((mNewAtom > 0))
   mAtom = min(mAtom, mNewAtom)

   call move_alloc(self%image, itmp)
   allocate(self%image(mNewAtom), source=0)
   self%image(:mAtom) = itmp(:mAtom)

   itmp = self%iCellVec
   deallocate(self%iCellVec)
   allocate(self%iCellVec(mNewAtom), source=0)
   self%iCellVec(:mAtom) = itmp(:mAtom)

   call move_alloc(self%coords, atmp)
   allocate(self%coords(3, mNewAtom), source=0.0_wp)
   self%coords(:, :mAtom) = atmp(:, :mAtom)

end subroutine neighlist_resize_coord

end module tbdef_neighbourlist
