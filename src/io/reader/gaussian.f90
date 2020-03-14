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

module xtb_io_reader_gaussian
   use xtb_mctc_symbols, only : toNumber, toSymbol
   use xtb_type_molecule
   implicit none
   private

   public :: readMoleculeGaussianExternal


contains


subroutine readMoleculeGaussianExternal(mol, unit, status, iomsg)
   use xtb_mctc_accuracy, only : wp
   logical, parameter :: debug = .false.
   class(TMolecule), intent(out) :: mol
   integer, intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   integer :: err
   integer :: n, mode, chrg, spin, iat, ii
   integer, allocatable :: at(:)
   real(wp), allocatable :: xyz(:,:)
   real(wp) :: coord(3)
   real(wp) :: q
   real(wp) :: conv

   status = .false.

   read(unit, '(4i10)', iostat=err) n, mode, chrg, spin
   if (err.ne.0) then
      iomsg = "Could not read number of atoms, check format!"
      return
   endif

   if (n <= 0) then
      iomsg = "Found no atoms, cannot work without atoms!"
      return
   end if

   allocate(xyz(3, n))
   allocate(at(n))

   ii = 0
   do while (ii < n)
      read(unit, '(i10, 4f20.12)', iostat=err) iat, coord, q
      if (is_iostat_end(err)) exit
      if (err.ne.0) then
         iomsg = "Could not read geometry from Gaussian file"
         return
      endif
      if (iat > 0) then
         ii = ii+1
         at(ii) = iat
         xyz(:, ii) = coord
      else
         iomsg = "Invalid atomic number"
         return
      end if
   end do

   call init(mol, at, xyz, real(chrg, wp), spin)

   if (ii /= n) then
      iomsg = "Atom number missmatch in Gaussian file"
      return
   endif

   status = .true.

end subroutine readMoleculeGaussianExternal


end module xtb_io_reader_gaussian
