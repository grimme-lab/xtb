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

module xtb_io_writer_turbomole
   use xtb_type_molecule, only : TMolecule, len
   implicit none
   private

   public :: writeMoleculeCoord


contains


subroutine writeMoleculeCoord(mol, unit)
   class(TMolecule), intent(in) :: mol
   integer, intent(in) :: unit
   integer :: i

   write(unit,'(a)') "$coord"
   do i = 1, len(mol)
      write(unit,'(3f20.14,6x,a)') mol%xyz(:,i), trim(mol%sym(i))
   enddo
   write(unit,'(a,1x,i0)') "$periodic", mol%npbc
   if (mol%npbc > 0) then
      write(unit,'(a)') "$lattice bohr"
      write(unit,'(3f20.14)') mol%lattice
   endif
   write(unit,'(a)') "$end"

end subroutine writeMoleculeCoord


end module xtb_io_writer_turbomole
