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
   use xtb_mctc_symbols, only : toNumber, symbolLength
   use xtb_mctc_systools
   use xtb_pbc_tools
   use xtb_type_molecule
   use xtb_type_reader
   use xtb_type_vendordata, only : vasp_info
   implicit none
   private

   public :: readMoleculeGenFormat
   public :: readHessianDFTBPlus


contains


subroutine readMoleculeGenFormat(mol, unit, status, iomsg)
   logical, parameter :: debug = .false.
   type(TMolecule),intent(out) :: mol
   integer,intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg

   character(len=:), allocatable :: line
   integer :: natoms, nspecies, iatom, dummy, isp, ilat, error
   logical :: cartesian, periodic
   real(wp) :: coord(3), lattice(3, 3)
   character(len=1) :: variant
   character(len=symbolLength), allocatable :: species(:), sym(:)
   real(wp), allocatable :: xyz(:, :), abc(:, :)

   status = .false.

   call next_line(unit, line, error)
   read(line, *, iostat=error) natoms, variant
   if (error /= 0 .or. natoms < 1) then
      iomsg = 'could not read number of atoms'
      return
   endif

   allocate(species(natoms))
   allocate(sym(natoms))
   allocate(xyz(3, natoms))
   allocate(abc(3, natoms))

   select case(variant)
   case('c', 'C')
      cartesian = .true.
      periodic = .false.
   case('s', 'S')
      cartesian = .true.
      periodic = .true.
   case('f', 'F')
      cartesian = .false.
      periodic = .true.
   case default
      iomsg = 'invalid input version'
      return
   endselect

   call next_line(unit, line, error)
   call parse(line, ' ', species, nspecies)
   if (any(toNumber(species(:nspecies)) == 0)) then
      iomsg = 'unknown atom type present'
      return
   endif

   do iatom = 1, natoms
      call next_line(unit, line, error)
      read(line, *, iostat=error) dummy, isp, coord
      if (error /= 0) then
         iomsg = 'could not read coordinates from file'
         return
      endif
      sym(iatom) = species(isp)
      if (cartesian) then
         xyz(:, iatom) = coord * aatoau
      else
         abc(:, iatom) = coord
      endif
   enddo

   if (periodic) then
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
         lattice(:, ilat) = coord * aatoau
      enddo
      if (.not.cartesian) then
         call abc_to_xyz(natoms, lattice, abc, xyz)
      endif
      call init(mol, sym, xyz, lattice=lattice)
   else
      call init(mol, sym, xyz)
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


subroutine readHessianDFTBPlus(hessian, reader, mol, status, iomsg)
   real(wp), intent(out) :: hessian(:, :)
   type(TMolecule), intent(in) :: mol
   type(TReader), intent(inout) :: reader
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   character(len=:), allocatable :: line
   character(len=32) :: buffer
   integer :: error, ndim, ii, jj, jbatch, iline

   status = .false.
   iline = 0

   if (error /= 0) then
      iomsg = "Could not find $hessian data group"
      return
   end if

   ndim = 3*len(mol)

   read(reader%unit, *, iostat=error) ((hessian(jj, ii), jj=1, ndim), ii=1, ndim)

   if (error /= 0) then
      if (is_iostat_end(error)) then
         iomsg = "Unexpected end of file while reading hessian"
      else
         write(buffer, '(i0)') iline
         iomsg = "Failed to read hessian in line "//trim(buffer)
      end if
      return
   end if

   status = .true.

end subroutine readHessianDFTBPlus


end module xtb_io_reader_genformat
