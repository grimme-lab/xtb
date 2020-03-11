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

module xtb_io_writer_pdb
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autoaa
   use xtb_type_molecule, only : TMolecule, len
   implicit none
   private

   public :: writeMoleculePDB


contains


subroutine writeMoleculePDB(mol, unit, number)
   type(TMolecule), intent(in) :: mol
   integer, intent(in) :: unit
   integer, intent(in), optional :: number
   character(len=6) :: w1
   character(len=2) :: a_charge
   character(len=1) :: last_chain
   logical :: last_het
   integer :: offset
   integer :: iatom, jatom
   real(wp) :: xyz(3)
! ATOM   2461  HA3 GLY A 153     -10.977  -7.661   2.011  1.00  0.00           H
! TER    2462      GLY A 153
! a6----i5---xa4--aa3-xai4--axxxf8.3----f8.3----f8.3----f6.2--f6.2--xxxxxxa4--a2a2
! HETATM 2463  CHA HEM A 154       9.596 -13.100  10.368  1.00  0.00           C
   character(len=*), parameter :: pdb_format = &
      &  '(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)'


   offset = 0
   last_chain = mol%pdb(1)%chains
   last_het = mol%pdb(1)%het
   if (present(number)) write(unit, '("MODEL ",4x,i4)') number
   do iatom = 1, len(mol)

      ! handle the terminator
      if (mol%pdb(iatom)%het .neqv. last_het) then
         write(unit, '("TER   ",i5,6x,a3,1x,a1,i4)') iatom + offset, &
            &  mol%pdb(iatom-1)%residue, last_chain, mol%pdb(iatom)%residue_number
         last_het = .not.last_het
         offset = offset+1
      else if (mol%pdb(iatom)%chains /= last_chain) then
         write(unit, '("TER   ",i5,6x,a3,1x,a1,i4)') iatom + offset, &
            &  mol%pdb(iatom-1)%residue, last_chain, mol%pdb(iatom)%residue_number
         offset = offset+1
      endif

      jatom = iatom + offset
      if (mol%pdb(iatom)%het) then
         w1 = 'HETATM'
      else
         w1 = 'ATOM  '
      endif


      xyz = mol%xyz(:,iatom) * autoaa
      if (mol%pdb(iatom)%charge < 0) then
         write(a_charge, '(i1,"-")') abs(mol%pdb(iatom)%charge)
      else if (mol%pdb(iatom)%charge > 0) then
         write(a_charge, '(i1,"+")') abs(mol%pdb(iatom)%charge)
      else
         a_charge = '  '
      endif

      write(unit, pdb_format) &
         &  w1, jatom, mol%pdb(iatom)%name, mol%pdb(iatom)%loc, &
         &  mol%pdb(iatom)%residue, mol%pdb(iatom)%chains, mol%pdb(iatom)%residue_number, &
         &  mol%pdb(iatom)%code, xyz, 1.0_wp, 0.0_wp, mol%pdb(iatom)%segid, &
         &  mol%sym(iatom), a_charge
   enddo
   if (present(number)) then
      write(unit, '("ENDMDL")')
   else
      write(unit, '("END")')
   endif

end subroutine writeMoleculePDB


end module xtb_io_writer_pdb
