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

module xtb_pbc_optimizer_filter_cart
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule, only : TMolecule
   use xtb_mctc_math, only : matInv3x3
   use xtb_pbc_optimizer_filter_type, only : filter_type
   implicit none
   private

   public :: cartesian_filter, new_cartesian_filter

   type, extends(filter_type) :: cartesian_filter
      integer :: nvar
   contains
      procedure :: get_dimension
      procedure :: transform_structure
      procedure :: transform_derivative
   end type cartesian_filter

contains


subroutine new_cartesian_filter(self, mol, optcell)
   !> Instance of the transformation filter
   type(cartesian_filter), intent(out) :: self
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol
   !> Relax lattice parameters
   logical, intent(in) :: optcell

   self%optcell = optcell
   self%nvar = size(mol%xyz)
   if (self%optcell) self%nvar = self%nvar + size(mol%lattice)
end subroutine new_cartesian_filter


pure function get_dimension(self) result(n)
   !> Instance of the transformation filter
   class(cartesian_filter), intent(in) :: self
   integer :: n

   n = self%nvar
end function get_dimension


subroutine transform_structure(self, mol, displacement)
   !> Instance of the transformation filter
   class(cartesian_filter), intent(in) :: self
   !> Molecular structure data
   type(TMolecule), intent(inout) :: mol
   real(wp), intent(in) :: displacement(:)

   mol%xyz(:, :) = mol%xyz + reshape(displacement(:size(mol%xyz)), shape(mol%xyz))
   if (self%optcell) then
      mol%lattice(:, :) = mol%lattice &
         & + reshape(displacement(size(mol%xyz)+1:), shape(mol%lattice))
   end if
end subroutine transform_structure


subroutine transform_derivative(self, mol, energy, gradient, sigma, val, deriv)
   !> Instance of the transformation filter
   class(cartesian_filter), intent(in) :: self
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol
   real(wp), intent(in) :: energy
   real(wp), intent(in) :: gradient(:, :)
   real(wp), intent(in) :: sigma(:, :)
   real(wp), intent(out) :: val
   real(wp), intent(out) :: deriv(:)

   real(wp) :: inv_lat(3, 3), lat_grad(3, 3)

   val = energy
   deriv(:size(gradient)) = reshape(gradient, [size(gradient)])
   if (self%optcell) then
      inv_lat(:, :) = transpose(matInv3x3(mol%lattice))
      lat_grad(:, :) = matmul(sigma, inv_lat) &
         & - matmul(gradient, matmul(transpose(mol%xyz), inv_lat))
      deriv(size(gradient)+1:) = reshape(lat_grad, [size(lat_grad)])
   end if
end subroutine transform_derivative


end module xtb_pbc_optimizer_filter_cart
