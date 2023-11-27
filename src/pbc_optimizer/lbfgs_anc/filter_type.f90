! This file is part of ancopt.
! SPDX-Identifier: LGPL-3.0-or-later
!
! ancopt is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! ancopt is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with ancopt.  If not, see <https://www.gnu.org/licenses/>.

module xtb_pbc_optimizer_filter_type
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule, only : TMolecule
!   use mctc_io, only : structure_type
   implicit none
   private

   public :: filter_type

   type, abstract :: filter_type
      logical :: optcell = .false.
   contains
      procedure(get_dimension), deferred :: get_dimension
      procedure(transform_structure), deferred :: transform_structure
      procedure(transform_derivative), deferred :: transform_derivative
   end type filter_type

   abstract interface
      pure function get_dimension(self) result(n)
         import :: filter_type
         class(filter_type), intent(in) :: self
         integer :: n
      end function get_dimension

      subroutine transform_structure(self, mol, displacement)
         import :: filter_type, TMolecule, wp
         class(filter_type), intent(in) :: self
         type(TMolecule), intent(inout) :: mol
         real(wp), intent(in) :: displacement(:)
      end subroutine transform_structure

      subroutine transform_derivative(self, mol, energy, gradient, sigma, val, deriv)
         import :: filter_type, TMolecule, wp
         class(filter_type), intent(in) :: self
         type(TMolecule), intent(in) :: mol
         real(wp), intent(in) :: energy
         real(wp), intent(in) :: gradient(:, :)
         real(wp), intent(in) :: sigma(:, :)
         real(wp), intent(out) :: val
         real(wp), intent(out) :: deriv(:)
      end subroutine transform_derivative
   end interface

end module xtb_pbc_optimizer_filter_type
