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
   integer :: n, mode, chrg, spin, iat
   real(wp) :: xyz(3)
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

   call mol%allocate(n)
   mol%npbc = 0
   mol%pbc = .false.
   mol%chrg = chrg
   mol%uhf = spin

   n = 0
   do while (n < mol%n)
      read(unit, '(i10, 4f20.12)', iostat=err) iat, xyz, q
      if (is_iostat_end(err)) exit
      if (err.ne.0) then
         iomsg = "Could not read geometry from Gaussian file"
         return
      endif
      if (iat > 0) then
         n = n+1
         mol%at(n) = iat
         mol%sym(n) = toSymbol(iat)
         mol%xyz(:, n) = xyz
      else
         iomsg = "Invalid atomic number"
         return
      end if
   end do

   if (n /= mol%n) then
      iomsg = "Atom number missmatch in Gaussian file"
      return
   endif

   status = .true.

end subroutine readMoleculeGaussianExternal


end module xtb_io_reader_gaussian
