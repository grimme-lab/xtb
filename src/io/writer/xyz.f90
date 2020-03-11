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

module xtb_io_writer_xyz
   use xtb_type_molecule, only : TMolecule, len
   implicit none
   private

   public :: writeMoleculeXYZ

contains


subroutine writeMoleculeXYZ(mol, unit, comment_line)
   use xtb_mctc_convert
   class(TMolecule), intent(in) :: mol
   integer, intent(in) :: unit
   character(len=*), intent(in) :: comment_line
   integer :: i

   write(unit, '(i0)') len(mol)
   write(unit, '(a)') comment_line
   do i = 1, len(mol)
      write(unit, '(a4,2x,3f20.14)') mol%sym(i), mol%xyz(:,i)*autoaa
   enddo

end subroutine writeMoleculeXYZ


end module xtb_io_writer_xyz
