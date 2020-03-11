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

module xtb_io_reader_genformat
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   use xtb_mctc_strings
   use xtb_mctc_symbols, only : toNumber, toSymbol
   use xtb_mctc_systools
   use xtb_pbc_tools
   use xtb_type_molecule
   use xtb_type_vendordata, only : vasp_info
   implicit none
   private

   public :: readMoleculeGenFormat


contains


subroutine readMoleculeGenFormat(mol, unit, status, iomsg)
   logical, parameter :: debug = .false.
   type(TMolecule),intent(out) :: mol
   integer,intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg

   character(len=:), allocatable :: line
   integer :: natoms, nspecies, iatom, dummy, isp, ilat, error
   logical :: cartesian
   real(wp) :: coord(3)
   integer, allocatable :: species(:)
   character(len=1) :: variant
   character(len=4), allocatable :: symbols(:)

   status = .false.

   call next_line(unit, line, error)
   read(line, *, iostat=error) natoms, variant
   if (error /= 0 .or. natoms < 1) then
      iomsg = 'could not read number of atoms'
      return
   endif

   call mol%allocate(natoms)
   allocate(symbols(len(mol)), source='    ')

   select case(variant)
   case('c', 'C')
      cartesian = .true.
   case('s', 'S')
      cartesian = .true.
      mol%npbc = 3
      mol%pbc = .true.
   case('f', 'F')
      cartesian = .false.
      mol%npbc = 3
      mol%pbc = .true.
   case default
      iomsg = 'invalid input version'
      return
   endselect

   call next_line(unit, line, error)
   call parse(line, ' ', symbols, nspecies)
   allocate(species(nspecies), source=0)
   species = toNumber(symbols(:nspecies))
   if (any(species == 0)) then
      iomsg = 'unknown atom type present'
      return
   endif

   do iatom = 1, len(mol)
      call next_line(unit, line, error)
      read(line, *, iostat=error) dummy, isp, coord
      if (error /= 0) then
         iomsg = 'could not read coordinates from file'
         return
      endif
      mol%at(iatom) = species(isp)
      mol%sym(iatom) = symbols(isp)
      if (cartesian) then
         mol%xyz(:, iatom) = coord * aatoau
      else
         mol%abc(:, iatom) = coord
      endif
   enddo

   if (mol%npbc > 0) then
      call next_line(unit, line, error)
      if (error /= 0) then
         iomsg = 'missing lattice information'
         return
      endif
      do ilat = 1, 3
         call next_line(unit, line, error)
         read(line, *, iostat=error) coord
         if (error /= 0) then
            iomsg = 'could not read lattice from file'
            return
         endif
         mol%lattice(:, ilat) = coord * aatoau
      enddo
      if (cartesian) then
         call xyz_to_abc(len(mol), mol%lattice, mol%xyz, mol%abc, mol%pbc)
      else
         call abc_to_xyz(len(mol), mol%lattice, mol%abc, mol%xyz)
      endif
   endif

   mol%vasp = vasp_info(cartesian=cartesian)

   status = .true.

contains

subroutine next_line(unit, line, error)
   integer,intent(in) :: unit
   character(len=:), allocatable, intent(out) :: line
   integer, intent(out) :: error
   integer :: ihash

   error = 0
   do while(error == 0)
      call getline(unit, line, error)
      ihash = index(line, '#')
      if (ihash > 0) line = line(:ihash-1)
      if (len_trim(line) > 0) exit
   enddo
   line = trim(adjustl(line))
end subroutine next_line

end subroutine readMoleculeGenFormat


end module xtb_io_reader_genformat
