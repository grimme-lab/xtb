! This file is part of xtb.
!
! Copyright (C) 2026 Leopold M. Seidler
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

!> Redundant internal coordinates of a system.
!>
!> The coordinate values are evaluated once, at the reference geometry handed
!> to init, and stored in q. They are a snapshot and go stale when the atoms
!> move; build the set again from the geometry of interest.
module xtb_internals_redundant
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_math, only : crossProd
   use xtb_basic_geo, only : bangl
   use xtb_internals_graph, only : graph_type
   implicit none
   private

   public :: redundant_type, init

   !> Bond coordinate kind, one coordinate per edge
   integer, parameter :: coord_bond = 1
   !> Angle coordinate kind, one coordinate per pair of neighbours of a centre
   integer, parameter :: coord_angle = 2
   !> Dihedral coordinate kind, one coordinate per pair of outer bonds of an edge
   integer, parameter :: coord_dihedral = 3

   !> Redundant internal coordinates of a system.
   type :: redundant_type

      !> Number of atoms in the system
      integer :: n = 0

      !> Number of coordinates, the length of q and the row count of B
      integer :: nint = 0

      !> Number of bonds, angles and dihedrals, nint = nbond + nangle + ndihedral
      integer :: nbond = 0, nangle = 0, ndihedral = 0

      !> Atoms of each coordinate, coord(:, ic) = [i, j, k, l] with zero for
      !> atoms the coordinate kind does not use: a bond is [i, j, 0, 0], an
      !> angle [i, j, k, 0] with j the centre, a dihedral [i, j, k, l]
      integer, allocatable :: coord(:, :)

      !> Coordinate values at the reference geometry, q(self%nint): lengths in
      !> bohr, angles and dihedrals in radians
      real(wp), allocatable :: q(:)

   end type redundant_type

   !> Build the redundant coordinate set of a molecular graph
   interface init
      module procedure :: new_redundant
   end interface init

contains

!> Build the coordinate definitions from the molecular graph and evaluate them
!> at the reference geometry.
pure subroutine new_redundant(self, graph, xyz)
   !> Coordinate set to fill, its allocatable components are allocated here
   type(redundant_type), intent(out) :: self
   !> Molecular graph of the system
   type(graph_type), intent(in) :: graph
   !> Cartesian reference coordinates, dimension (3, graph%n), in the length
   !> unit the coordinate values are reported in
   real(wp), intent(in) :: xyz(:, :)

   integer :: i, j, k, a, b, ia, la, ndeg

   self%n = graph%n

   ! count the coordinates of each kind
   self%nbond = 0
   do i = 1, graph%n
      do k = 1, graph%neighs(i)
         if (graph%i_neigh(k, i) > i) self%nbond = self%nbond + 1
      end do
   end do

   self%nangle = 0
   self%ndihedral = 0
   do i = 1, graph%n
      ndeg = graph%neighs(i)
      self%nangle = self%nangle + ndeg*(ndeg - 1)/2
      do k = 1, ndeg
         ! every edge is enumerated once, from its smaller atom
         j = graph%i_neigh(k, i)
         if (j < i) cycle
         ! a dihedral needs an outer bond at both ends of the edge
         do a = 1, ndeg
            ia = graph%i_neigh(a, i)
            if (ia == j) cycle
            do b = 1, graph%neighs(j)
               la = graph%i_neigh(b, j)
               if (la == i .or. la == ia) cycle
               self%ndihedral = self%ndihedral + 1
            end do
         end do
      end do
   end do

   self%nint = self%nbond + self%nangle + self%ndihedral
   allocate(self%coord(4, self%nint), source = 0)

   ! fill bonds, then angles, then dihedrals
   ia = 0
   do i = 1, graph%n
      do k = 1, graph%neighs(i)
         j = graph%i_neigh(k, i)
         if (j < i) cycle
         ia = ia + 1
         self%coord(1, ia) = i
         self%coord(2, ia) = j
      end do
   end do

   do i = 1, graph%n
      ndeg = graph%neighs(i)
      do a = 1, ndeg - 1
         do b = a + 1, ndeg
            ia = ia + 1
            self%coord(1, ia) = graph%i_neigh(a, i)
            self%coord(2, ia) = i
            self%coord(3, ia) = graph%i_neigh(b, i)
         end do
      end do
   end do

   do i = 1, graph%n
      do k = 1, graph%neighs(i)
         j = graph%i_neigh(k, i)
         if (j < i) cycle
         do a = 1, graph%neighs(i)
            if (graph%i_neigh(a, i) == j) cycle
            do b = 1, graph%neighs(j)
               la = graph%i_neigh(b, j)
               if (la == i .or. la == graph%i_neigh(a, i)) cycle
               ia = ia + 1
               self%coord(1, ia) = graph%i_neigh(a, i)
               self%coord(2, ia) = i
               self%coord(3, ia) = j
               self%coord(4, ia) = la
            end do
         end do
      end do
   end do

   ! every coordinate of the set was written
   if (ia /= self%nint) error stop "redundant: coordinate count mismatch"

   ! evaluate the coordinate values at the reference geometry
   allocate(self%q(self%nint), source = 0.0_wp)
   do ia = 1, self%nbond
      self%q(ia) = dist(xyz, self%coord(1, ia), self%coord(2, ia))
   end do
   do ia = self%nbond + 1, self%nbond + self%nangle
      call bangl(xyz, self%coord(1, ia), self%coord(2, ia), self%coord(3, ia), self%q(ia))
   end do
   do ia = self%nbond + self%nangle + 1, self%nint
      self%q(ia) = dihedral_value(xyz(:, self%coord(1, ia)), xyz(:, self%coord(2, ia)), &
         & xyz(:, self%coord(3, ia)), xyz(:, self%coord(4, ia)))
   end do

end subroutine new_redundant

!> Distance between atoms i and j
pure function dist(xyz, i, j) result(r)
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atom indices of the pair
   integer, intent(in) :: i, j

   real(wp) :: r
   r = sqrt((xyz(1, i) - xyz(1, j))**2 + (xyz(2, i) - xyz(2, j))**2 + (xyz(3, i) - xyz(3, j))**2)
end function dist

!> Signed dihedral i-j-k-l in (-pi, pi]
pure function dihedral_value(xi, xj, xk, xl) result(tau)
   !> Cartesian coordinates of the four dihedral atoms i-j-k-l
   real(wp), intent(in) :: xi(3), xj(3), xk(3), xl(3)

   real(wp) :: tau, b1(3), b2(3), b3(3), n1(3), n2(3), m(3)

   b1 = xj - xi; b2 = xk - xj; b3 = xl - xk
   n1 = crossProd(b1, b2)
   n2 = crossProd(b2, b3)
   m = crossProd(n1, b2/norm2(b2))
   tau = atan2(dot_product(m, n2), dot_product(n1, n2))
end function dihedral_value

end module xtb_internals_redundant
