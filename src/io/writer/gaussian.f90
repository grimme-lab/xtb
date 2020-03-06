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

module xtb_io_writer_gaussian
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule, only : TMolecule, len
   implicit none
   private

   public :: writeMoleculeGaussianExternal


contains


subroutine writeMoleculeGaussianExternal(mol, unit)
   type(TMolecule), intent(in) :: mol
   integer, intent(in) :: unit
   integer :: iat

   write(unit, '(4i10)') len(mol), 1, nint(mol%chrg), mol%uhf
   do iat = 1, len(mol)
      write(unit, '(i10,4f20.12)') mol%at(iat), mol%xyz(:, iat), 0.0_wp
   end do

end subroutine writeMoleculeGaussianExternal


end module xtb_io_writer_gaussian
