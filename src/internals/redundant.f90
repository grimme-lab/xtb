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
   use xtb_mctc_constants, only : pi
   use xtb_mctc_math, only : crossProd
   use xtb_basic_geo, only : bangl
   use xtb_bmatrix, only : linbend_frame
   use xtb_internals_graph, only : graph_type
   use xtb_internals_type, only : internal_coords_set_type, coord_bond, &
      & coord_angle, coord_dihedral, coord_linbend
   implicit none
   private

   public :: redundant_type, init

   !> A centre is linear when its reference angle is within five degrees of pi
   real(wp), parameter :: lin_tol = 5.0_wp*pi/180.0_wp

   !> Redundant internal coordinates of a system. The coordinate definitions
   !> are the inherited ones, the counts and the reference values are the
   !> redundant set's own state.
   type, extends(internal_coords_set_type) :: redundant_type

      !> Number of atoms in the system
      integer :: n = 0

      !> Number of bonds
      integer :: nbond = 0

      !> Number of angular coordinate rows; a linear angle contributes two
      integer :: nangle = 0

      !> Number of dihedrals
      integer :: ndihedral = 0

      !> Coordinate values at the reference geometry, q(self%ncoords): lengths
      !> in bohr, regular angles and dihedrals in radians, and linear-bend
      !> projections dimensionless
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

   integer :: i, j, k, a, b, ia, la, ndeg, component
   real(wp) :: angle, b1(3), b2(3), b3(3), n1(3), n2(3), m(3)
   real(wp) :: v1(3), v2(3), bend_frame(3, 2)

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
      do a = 1, ndeg - 1
         do b = a + 1, ndeg
            call bangl(xyz, graph%i_neigh(a, i), i, graph%i_neigh(b, i), angle)
            if (angle > pi - lin_tol) then
               self%nangle = self%nangle + 2
            else
               self%nangle = self%nangle + 1
            end if
         end do
      end do
      do k = 1, ndeg
         ! every edge is enumerated once, from its smaller atom
         j = graph%i_neigh(k, i)
         if (j < i) cycle
         ! a dihedral needs a valid outer bond angle at both ends of the edge
         do a = 1, ndeg
            ia = graph%i_neigh(a, i)
            if (ia == j) cycle
            do b = 1, graph%neighs(j)
               la = graph%i_neigh(b, j)
               if (la == i .or. la == ia) cycle
               if (.not. valid_torsion(xyz, ia, i, j, la)) cycle
               self%ndihedral = self%ndihedral + 1
            end do
         end do
      end do
   end do

   self%ncoords = self%nbond + self%nangle + self%ndihedral
   allocate(self%kind(self%ncoords), self%atoms(4, self%ncoords), source = 0)
   allocate(self%frame(3, self%ncoords), self%q(self%ncoords), source = 0.0_wp)

   ! fill bonds, then angular rows, then dihedrals
   ia = 0
   do i = 1, graph%n
      do k = 1, graph%neighs(i)
         j = graph%i_neigh(k, i)
         if (j < i) cycle
         ia = ia + 1
         self%kind(ia) = coord_bond
         self%atoms(1, ia) = i
         self%atoms(2, ia) = j
         self%q(ia) = norm2(xyz(:, i) - xyz(:, j))
      end do
   end do

   do i = 1, graph%n
      ndeg = graph%neighs(i)
      do a = 1, ndeg - 1
         do b = a + 1, ndeg
            call bangl(xyz, graph%i_neigh(a, i), i, graph%i_neigh(b, i), angle)
            if (angle > pi - lin_tol) then
               v1 = xyz(:, graph%i_neigh(a, i)) - xyz(:, i)
               v1 = v1/norm2(v1)
               v2 = xyz(:, graph%i_neigh(b, i)) - xyz(:, i)
               v2 = v2/norm2(v2)
               call linbend_frame(v1, bend_frame(:, 1), bend_frame(:, 2))
               do component = 1, 2
                  ia = ia + 1
                  self%kind(ia) = coord_linbend
                  self%atoms(1, ia) = graph%i_neigh(a, i)
                  self%atoms(2, ia) = i
                  self%atoms(3, ia) = graph%i_neigh(b, i)
                  self%frame(:, ia) = bend_frame(:, component)
                  self%q(ia) = dot_product(self%frame(:, ia), v1 + v2)
               end do
            else
               ia = ia + 1
               self%kind(ia) = coord_angle
               self%atoms(1, ia) = graph%i_neigh(a, i)
               self%atoms(2, ia) = i
               self%atoms(3, ia) = graph%i_neigh(b, i)
               self%q(ia) = angle
            end if
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
               if (.not. valid_torsion(xyz, graph%i_neigh(a, i), i, j, la)) cycle
               ia = ia + 1
               self%kind(ia) = coord_dihedral
               self%atoms(1, ia) = graph%i_neigh(a, i)
               self%atoms(2, ia) = i
               self%atoms(3, ia) = j
               self%atoms(4, ia) = la
               b1 = xyz(:, i) - xyz(:, self%atoms(1, ia))
               b2 = xyz(:, j) - xyz(:, i)
               b3 = xyz(:, la) - xyz(:, j)
               n1 = crossProd(b1, b2)
               n2 = crossProd(b2, b3)
               m = crossProd(n1, b2/norm2(b2))
               self%q(ia) = atan2(dot_product(m, n2), dot_product(n1, n2))
            end do
         end do
      end do
   end do

   ! every coordinate of the set was written
   if (ia /= self%ncoords) error stop "redundant: coordinate count mismatch"

end subroutine new_redundant

!> Whether both interior angles define a nonsingular dihedral.
pure logical function valid_torsion(xyz, i, j, k, l)
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Candidate dihedral atoms
   integer, intent(in) :: i, j, k, l

   real(wp) :: angle1, angle2

   call bangl(xyz, i, j, k, angle1)
   call bangl(xyz, j, k, l, angle2)
   valid_torsion = angle1 >= lin_tol .and. angle1 <= pi - lin_tol &
      & .and. angle2 >= lin_tol .and. angle2 <= pi - lin_tol
end function valid_torsion


end module xtb_internals_redundant
