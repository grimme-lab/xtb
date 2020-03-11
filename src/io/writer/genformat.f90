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

module xtb_io_writer_genformat
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autoaa
   use xtb_mctc_symbols, only : toSymbol
   use xtb_type_molecule, only : TMolecule, len
   implicit none
   private

   public :: writeMoleculeGenFormat


contains


subroutine writeMoleculeGenFormat(mol, unit, comment_line)
   class(TMolecule), intent(in) :: mol
   integer, intent(in) :: unit
   character(len=*), intent(in) :: comment_line
   integer :: i,j,iat
   real(wp), parameter :: zero3(3) = 0.0_wp
   integer, allocatable :: species(:)
   character(len=4) :: sym

   write(unit, '(i0,1x)', advance='no') len(mol)
   if (mol%npbc == 0) then
      write(unit, '("C")') ! cluster
   else
      if (mol%vasp%cartesian) then
         write(unit, '("S")') ! supercell
      else
         write(unit, '("F")') ! fractional
      endif
   endif
   write(unit, '("#",1x,a)') comment_line

   allocate(species(len(mol)), source=0)
   j = 0
   do iat = minval(mol%at), maxval(mol%at)
      if (any(mol%at == iat)) then
         j = j+1
         sym = toSymbol(iat)
         write(unit,'(1x,a)',advance='no') trim(sym)
         where(mol%at == iat)
            species = j
         endwhere
      endif
   enddo
   write(unit,'(a)')

   if (mol%npbc == 0 .or. mol%vasp%cartesian) then
      ! now write the cartesian coordinates
      do i = 1, len(mol)
         write(unit,'(2i5,3es24.15)') i, species(i), mol%xyz(:,i)*autoaa
      enddo
   else
      ! now write the fractional coordinates
      do i = 1, len(mol)
         write(unit,'(2i5,3es24.15)') i, species(i), mol%abc(:,i)
      enddo
   endif
   deallocate(species)

   if (mol%npbc > 0) then
      ! scaling factor for lattice parameters is always one
      write(unit,'(3f20.14)') zero3
      ! write the lattice parameters
      do i = 1, 3
         write(unit,'(3f20.14)') mol%lattice(:,i)*autoaa
      enddo
   endif

end subroutine writeMoleculeGenFormat


end module xtb_io_writer_genformat
