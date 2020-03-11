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

module xtb_io_writer_vasp
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autoaa
   use xtb_type_molecule, only : TMolecule, len
   implicit none
   private

   public :: writeMoleculeVasp


contains


subroutine writeMoleculeVasp(mol, unit, comment_line)
   class(TMolecule), intent(in) :: mol
   integer, intent(in) :: unit
   character(len=*), intent(in) :: comment_line
   integer :: i,j,iat
   integer,allocatable :: kinds(:)

   allocate(kinds(mol%n), source = 1)

   ! use vasp 5.x format
   write(unit,'(a)') comment_line

   ! scaling factor for lattice parameters is always one
   write(unit,'(f20.14)') mol%vasp%scale
   ! write the lattice parameters
   do i = 1, 3
      write(unit,'(3f20.14)') mol%lattice(:,i)*autoaa/mol%vasp%scale
   enddo

   j = 0
   iat = 0
   do i = 1, len(mol)
      if (iat.eq.mol%at(i)) then
         kinds(j) = kinds(j)+1
      else
         j = j+1
         iat = mol%at(i)
         write(unit,'(1x,a)',advance='no') mol%sym(i)
      endif
   enddo
   write(unit,'(a)')

   ! write the count of the consequtive atom types
   do i = 1, j
      write(unit,'(1x,i0)',advance='no') kinds(i)
   enddo
   write(unit,'(a)')
   deallocate(kinds)

   if (mol%vasp%selective) write(unit,'("Selective")')

   ! we write cartesian coordinates
   if (mol%vasp%cartesian) then
      write(unit,'("Cartesian")')

      ! now write the cartesian coordinates
      do i = 1, len(mol)
         write(unit,'(3f20.14)') mol%xyz(:,i)*autoaa/mol%vasp%scale
      enddo
   else
      write(unit,'("Direct")')

      ! now write the fractional coordinates
      do i = 1, len(mol)
         write(unit,'(3f20.14)') mol%abc(:,i)
      enddo
   endif

end subroutine writeMoleculeVasp


end module xtb_io_writer_vasp
