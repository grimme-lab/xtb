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

module xtb_pbc_optimizer_filter
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule, only : TMolecule
   use xtb_pbc_optimizer_filter_cart, only : cartesian_filter, new_cartesian_filter
   !use xtb_pbc_optimizer_filter_anc, only : anc_input
   !use ancopt_filter_cartesian, only : cartesian_filter, new_cartesian_filter
   use xtb_pbc_optimizer_filter_type, only : filter_type
   implicit none
   private

   public :: filter_type, filter_input, new_filter

   type :: filter_input
      !type(anc_input), allocatable :: anc
      logical :: optcell = .false.
   end type filter_input

contains


subroutine new_filter(self, mol, input)
   !> Instance of the transformation filter
   class(filter_type), allocatable, intent(out) :: self
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol
   !> Input to generate transformation filter from
   type(filter_input), intent(in) :: input

   call new_cartesian_filter_(self, mol, input%optcell)

end subroutine new_filter


subroutine new_cartesian_filter_(self, mol, optcell)
   !> Instance of the transformation filter
   class(filter_type), allocatable, intent(out) :: self
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol
   !> Relax lattice parameters
   logical, intent(in) :: optcell

   type(cartesian_filter), allocatable :: new

   allocate(new)
   call new_cartesian_filter(new, mol, optcell .and. any(mol%pbc))
   call move_alloc(new, self)
end subroutine new_cartesian_filter_


end module xtb_pbc_optimizer_filter
