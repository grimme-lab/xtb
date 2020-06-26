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

module xtb_io_reader_ctfile
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   use xtb_mctc_symbols, only : toNumber, toSymbol, symbolLength
   use xtb_mctc_systools
   use xtb_type_molecule
   use xtb_type_vendordata, only : sdf_data
   implicit none
   private

   public :: readMoleculeSDF, readMoleculeMolfile


contains


subroutine readMoleculeSDF(mol, unit, status, iomsg)
   class(TMolecule), intent(out) :: mol
   integer, intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   character(len=:), allocatable :: line
   integer :: error

   call readMoleculeMolfile(mol, unit, status, iomsg)
   if (.not.status) return
   status = .false.

   error = 0
   do while(error == 0)
      call getline(unit, line, error)
      if (index(line, '$$$$') == 1) exit
      call mol%info%push_back(line)
   enddo
   if (error /= 0) then
      iomsg = "failed while reading SDF key-value pairs"
      return
   endif

   status = .true.

end subroutine readMoleculeSDF


subroutine readMoleculeMolfile(mol, unit, status, iomsg)
   class(TMolecule), intent(out) :: mol
   integer, intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   character(len=:), allocatable :: line
   character(len=:), allocatable :: name
   integer :: i, iatom, jatom, ibond, btype, atomtype
   integer :: error, length, charge(2,15)
   integer :: number_of_atoms, number_of_bonds, number_of_atom_lists, &
      &       chiral_flag, number_of_stext_entries, i999
   integer :: list4(4), list12(12)
   real(wp) :: x, y, z
   character(len=2) :: sdf_dim
   character(len=3) :: symbol
   character(len=5) :: v2000
   integer, parameter :: ccc_to_charge(0:7) = [0, +3, +2, +1, 0, -1, -2, -3]
   logical :: two_dim

   !> Element symbols
   character(len=symbolLength),allocatable :: sym(:)

   !> SDF specific information about atom types
   type(sdf_data), allocatable :: sdf(:)

   !> Cartesian coordinates in bohr
   real(wp),allocatable :: xyz(:,:)

   status = .false.
   two_dim = .false.

   call getline(unit, name, error)
   call getline(unit, line, error)
   read(line, '(20x,a2)', iostat=error) sdf_dim
   if (error == 0) then
      two_dim = sdf_dim == '2D' .or. sdf_dim == '2d'
   endif
   call getline(unit, line, error)
   call getline(unit, line, error)
   read(line, '(3i3,3x,2i3,12x,i3,1x,a5)', iostat=error) &
      & number_of_atoms, number_of_bonds, number_of_atom_lists, &
      & chiral_flag, number_of_stext_entries, i999, v2000
   if (error /= 0) then
      iomsg = "could not read header of molfile"
      return
   endif

   allocate(sdf(number_of_atoms), source=sdf_data())
   allocate(xyz(3, number_of_atoms))
   allocate(sym(number_of_atoms))

   do iatom = 1, number_of_atoms
      call getline(unit, line, error)
      read(line, '(3f10.4,1x,a3,i2,11i3)', iostat=error) &
         & x, y, z, symbol, list12
      if (error /= 0) then
         iomsg = "could not coordinates from connection table"
         return
      endif
      atomtype = toNumber(symbol)
      if (atomtype == 0) then
         iomsg = "unknown atom type '"//trim(symbol)//"' in connection table"
         return
      end if
      xyz(:, iatom) = [x, y, z] * aatoau
      sym(iatom) = trim(symbol)
      sdf(iatom)%isotope = list12(1)
      sdf(iatom)%charge = ccc_to_charge(list12(2)) ! drop doublet radical
      sdf(iatom)%hydrogens = list12(4)
      sdf(iatom)%valence = list12(6)
   enddo

   call init(mol, sym, xyz)
   call move_alloc(sdf, mol%sdf)
   if (len(name) > 0) mol%name = name
   if (two_dim) mol%struc%two_dimensional = .true.

   call mol%bonds%allocate(size=number_of_bonds, order=3)
   do ibond = 1, number_of_bonds
      call getline(unit, line, error)
      read(line, '(7i3)', iostat=error) &
         & iatom, jatom, btype, list4
      if (error /= 0) then
         iomsg = "could not topology from connection table"
         return
      endif
      call mol%bonds%push_back([iatom, jatom, btype])
   enddo

   do while(error == 0)
      call getline(unit, line, error)
      if (index(line, 'M  END') == 1) exit
      if (index(line, 'M  CHG') == 1) then
         read(line(7:9), *) length
         read(line(10:), '(*(1x,i3,1x,i3))') (charge(:, i), i=1, length)
         do i = 1, length
            mol%sdf(charge(1, i))%charge = charge(2, i)
         enddo
      endif
   enddo
   if (error /= 0) then
      iomsg = "could not read connection table"
      return
   endif

   mol%chrg = sum(mol%sdf%charge)

   if (any(mol%sdf%hydrogens > 1)) then
      iomsg = "explicit hydrogen atoms are required"
      return
   endif

   status = .true.

end subroutine readMoleculeMolfile


end module xtb_io_reader_ctfile
