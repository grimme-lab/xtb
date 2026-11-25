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
!
!> Internal coordinates of a system.
!>
!> A coordinate set is a plain definition of which coordinates exist and
!> which atoms they involve; the values of the coordinates are not part of
!> it. Concrete sets extend this type and add their own state.
module xtb_internals_type
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_math, only : crossProd
   implicit none
   private

   public :: internal_coords_set_type, coord_bond, coord_angle, coord_dihedral, &
      & coord_linbend, angle_value, dihedral_value

   !> Bond coordinate kind, one coordinate per edge
   integer, parameter :: coord_bond = 1
   !> Ordinary-angle coordinate kind, one coordinate per nonlinear neighbour pair
   integer, parameter :: coord_angle = 2
   !> Dihedral coordinate kind, one coordinate per pair of outer bonds of an edge
   integer, parameter :: coord_dihedral = 3
   !> Linear-bend coordinate kind, two fixed-frame coordinates per linear angle
   integer, parameter :: coord_linbend = 4

   !> Definition of a set of internal coordinates
   type :: internal_coords_set_type

      !> Number of coordinates in the set
      integer :: ncoords = 0

      !> Kind of each coordinate, kind(self%ncoords), one of the coord_* kinds
      integer, allocatable :: kind(:)

      !> Atoms of each coordinate, allocated as atoms(4, self%ncoords).
      !> atoms(:, ic) = [i, j, 0, 0] for a bond, [i, j, k, 0] for an angle
      !> or linear bend with j the centre, and [i, j, k, l] for a dihedral,
      !> with zero for the atoms the coordinate kind does not use
      integer, allocatable :: atoms(:, :)

      !> Fixed unit direction of each linear bend, allocated as
      !> frame(3, self%ncoords), zero for all other coordinate kinds
      real(wp), allocatable :: frame(:, :)

   end type internal_coords_set_type


contains


!> Valence angle i-j-k with centre j, in radians in [0, pi]
pure function angle_value(xyz, i, j, k) result(val)
   !> Cartesian coordinates, dimension (3, nat)
   real(wp), intent(in) :: xyz(:, :)
   !> Atom indices, centre at j
   integer, intent(in) :: i, j, k
   !> Valence angle in radians
   real(wp) :: val

   real(wp) :: v1(3), v2(3)

   v1 = xyz(:, i) - xyz(:, j)
   v2 = xyz(:, k) - xyz(:, j)
   val = atan2(norm2(crossProd(v1, v2)), dot_product(v1, v2))
end function angle_value


!> IUPAC (Klyne-Prelog) dihedral angle i-j-k-l, in radians in [-pi, pi]
pure function dihedral_value(xyz, i, j, k, l) result(val)
   !> Cartesian coordinates, dimension (3, nat)
   real(wp), intent(in) :: xyz(:, :)
   !> Atom indices along the torsion
   integer, intent(in) :: i, j, k, l
   !> Dihedral angle in radians
   real(wp) :: val

   real(wp) :: b1(3), b2(3), b3(3), n1(3), n2(3), m(3)

   b1 = xyz(:, j) - xyz(:, i)
   b2 = xyz(:, k) - xyz(:, j)
   b3 = xyz(:, l) - xyz(:, k)
   n1 = crossProd(b1, b2)
   n2 = crossProd(b2, b3)
   m = crossProd(b2 / norm2(b2), n1)
   val = atan2(dot_product(m, n2), dot_product(n1, n2))
end function dihedral_value

end module xtb_internals_type
