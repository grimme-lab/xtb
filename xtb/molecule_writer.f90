! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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

submodule(tbdef_molecule) molecule_writer
   use tbdef_molecule
   implicit none

contains

module subroutine write_molecule_generic(self, unit, format, energy, gnorm, number)
   use tbmod_file_utils
   include 'xtb_version.fh'
   class(tb_molecule), intent(in) :: self
   integer, intent(in) :: unit
   integer, intent(in), optional :: format
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gnorm
   integer, intent(in), optional :: number
   character(len=:), allocatable :: comment_line
   character(len=20) :: energy_line
   character(len=20) :: gnorm_line
   integer :: ftype
   if (present(format)) then
      ftype = format
   else
      ftype = self%ftype
   endif

   comment_line = ''
   if (present(energy)) then
      write(energy_line, '(f20.12)') energy
      comment_line = comment_line // " energy: " // trim(adjustl(energy_line))
   endif
   if (present(gnorm)) then
      write(gnorm_line, '(f20.12)') gnorm
      comment_line = comment_line // " gnorm: " // trim(adjustl(gnorm_line))
   endif
   comment_line = comment_line // " xtb: " // version

   select case(ftype)
   case(p_ftype%xyz)
      call write_xyz(self, unit, trim(comment_line))
   case(p_ftype%tmol)
      call write_tmol(self, unit)
   case(p_ftype%molfile)
      call write_molfile(self, unit, trim(comment_line))
   case(p_ftype%sdf)
      call write_sdf(self, unit, energy, gnorm)
   case(p_ftype%vasp)
      call write_vasp(self, unit, trim(comment_line))
   case(p_ftype%pdb)
      if (present(number)) then
         call write_pdb(self, unit, number)
      else
         call write_pdb(self, unit)
      endif
   end select

end subroutine write_molecule_generic

subroutine write_xyz(mol, unit, comment_line)
   use mctc_econv
   class(tb_molecule), intent(in) :: mol
   integer, intent(in) :: unit
   character(len=*), intent(in) :: comment_line
   integer :: i

   write(unit, '(i0)') mol%n
   write(unit, '(a)') comment_line
   do i = 1, mol%n
      write(unit, '(a2,6x,3f20.14)') mol%sym(i), mol%xyz(:,i)*autoaa
   enddo

end subroutine write_xyz

subroutine write_tmol(mol, unit)
   class(tb_molecule), intent(in) :: mol
   integer, intent(in) :: unit
   integer :: i

   write(unit,'(a)') "$coord"
   do i = 1, mol%n
      write(unit,'(3f20.14,6x,a2)') mol%xyz(:,i), mol%sym(i)
   enddo
   write(unit,'(a,1x,i0)') "$periodic", mol%npbc
   if (mol%npbc > 0) then
      write(unit,'(a)') "$lattice bohr"
      write(unit,'(3f20.14)') mol%lattice
   endif
   write(unit,'(a)') "$end"

end subroutine write_tmol

subroutine write_sdf(mol, unit, energy, gnorm)
   use tbdef_buffer
   include 'xtb_version.fh'
   class(tb_molecule), intent(in) :: mol
   integer, intent(in) :: unit
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gnorm
   type(tb_buffer) :: sd_values
   character(len=:), allocatable :: line
   character(len=*), parameter :: sd_format = &
      & '("> <",a,">",/,f20.12,/)'

   call write_molfile(mol, unit, "xtb: "//version)

   sd_values = mol%info
   call sd_values%reset
   do while(sd_values%next())
      call sd_values%getline(line)
      write(unit, '(a)') line
   enddo

   if (present(energy)) then
      write(unit, sd_format) "total energy / Eh", energy
   endif

   if (present(gnorm)) then
      write(unit, sd_format) "gradient norm / Eh/a0", gnorm
   endif

   write(unit, '("$$$$")')

end subroutine write_sdf

subroutine write_molfile(mol, unit, comment_line)
   use mctc_econv
   class(tb_molecule), intent(in) :: mol
   integer, intent(in) :: unit
   character(len=*), intent(in) :: comment_line
   integer, parameter :: list4(4) = 0
   integer :: iatom, ibond, iatoms(2), list12(12)
   logical :: has_sdf_data
   integer, parameter :: charge_to_ccc(-3:3) = [7, 6, 5, 0, 3, 2, 1]
   character(len=8)  :: date
   character(len=10) :: time

   call date_and_time(date, time)

   write(unit, '(a)') mol%name
   write(unit, '(2x,"xtb",5x,3a2,a4,"3D")') &
      &  date(5:6), date(7:8), date(3:4), time(:4)
   write(unit, '(a)') comment_line
   write(unit, '(3i3,3x,2i3,12x,i3,1x,a5)') &
      &  len(mol), len(mol%bonds), 0, 0, 0, 999, 'V2000'

   has_sdf_data = allocated(mol%sdf)

   do iatom = 1, mol%n
      if (has_sdf_data) then
         list12 = [mol%sdf(iatom)%isotope, charge_to_ccc(mol%sdf(iatom)%charge), &
            &      0, 0, 0, mol%sdf(iatom)%valence, 0, 0, 0, 0, 0, 0]
      else
         list12 = 0
      endif
      write(unit, '(3f10.4,1x,a3,i2,11i3)') &
         & mol%xyz(:, iatom)*autoaa, mol%sym(iatom), list12
   enddo

   do ibond = 1, len(mol%bonds)
      call mol%bonds%get_item(ibond, iatoms)
      write(unit, '(7i3)') &
         & iatoms, 1, list4
   enddo

   if (nint(mol%chrg) /= 0) then
      write(unit, '(a,*(1x,i3))') "M  CHG", 1, 1, nint(mol%chrg)
   endif

   write(unit, '(a)') "M  END"

end subroutine write_molfile

subroutine write_vasp(mol, unit, comment_line)
   use mctc_econv
   class(tb_molecule), intent(in) :: mol
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

end subroutine write_vasp

subroutine write_pdb(mol, unit, number)
   use mctc_econv
   type(tb_molecule), intent(in) :: mol
   integer, intent(in) :: unit
   integer, intent(in), optional :: number
   character(len=6) :: w1
   character(len=2) :: a_charge
   character(len=1) :: last_chain
   logical :: last_het
   integer :: offset
   integer :: iatom, jatom, iresidue
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
            &  mol%pdb(iatom-1)%residue, last_chain, iresidue
         last_het = .not.last_het
         offset = offset+1
      else if (mol%pdb(iatom)%chains /= last_chain) then
         write(unit, '("TER   ",i5,6x,a3,1x,a1,i4)') iatom + offset, &
            &  mol%pdb(iatom-1)%residue, last_chain, iresidue
         offset = offset+1
      endif

      jatom = iatom + offset
      if (mol%pdb(iatom)%het) then
         w1 = 'HETATM'
      else
         w1 = 'ATOM  '
      endif

      iresidue = mol%frag%list(iatom)

      xyz = mol%xyz(:,iatom) * autoaa
      if (mol%pdb(iatom)%charge /= 0) then
         write(a_charge, '(i2)') mol%pdb(iatom)%charge
      else
         a_charge = '  '
      endif

      write(unit, pdb_format) &
         &  w1, jatom, mol%pdb(iatom)%name, mol%pdb(iatom)%loc, &
         &  mol%pdb(iatom)%residue, mol%pdb(iatom)%chains, iresidue, &
         &  mol%pdb(iatom)%code, xyz, 1.0_wp, 0.0_wp, mol%pdb(iatom)%segid, &
         &  mol%sym(iatom), a_charge
   enddo
   if (present(number)) then
      write(unit, '("ENDMDL")')
   else
      write(unit, '("END")')
   endif

end subroutine write_pdb

end submodule molecule_writer
