! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> Implementation of orca output formats
module xtb_io_writer_orca
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule
   implicit none
   private

   public :: writeResultsOrca


contains


subroutine writeResultsOrca(unit, mol, energy, gradient)
   type(TMolecule), intent(in) :: mol
   integer, intent(in) :: unit
   real(wp), intent(in) :: energy
   real(wp), intent(in) :: gradient(:, :)
   integer :: i

   write(unit, '(a)') "#", "# Number of atoms", "#"
   write(unit, '(i10)') mol%n
   write(unit, '(a)') "#", "# The current total energy in Eh", "#"
   write(unit, '(f20.12)') energy
   write(unit, '(a)') "#", "# The current gradient in Eh/bohr", "#"
   write(unit, '(1x,f20.12)') gradient
   write(unit, '(a)') "#", "# The atomic numbers and current coordinates in Bohr", "#"
   do i = 1, mol%n
      write(unit, '(1x,i3,1x,3(1x,f12.7))') mol%at(i), mol%xyz(:, i)
   end do

end subroutine writeResultsOrca


end module xtb_io_writer_orca
