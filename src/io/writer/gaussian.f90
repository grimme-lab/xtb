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

   public :: writeResultsGaussianExternal


contains


subroutine writeResultsGaussianExternal(unit, energy, dipole, gradient)
   integer, intent(in) :: unit
   real(wp), intent(in) :: energy
   real(wp), intent(in) :: dipole(:)
   real(wp), intent(in) :: gradient(:, :)

   write(unit, '(4D20.12)') energy, dipole
   write(unit, '(3D20.12)') gradient

end subroutine writeResultsGaussianExternal


end module xtb_io_writer_gaussian
