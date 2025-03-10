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
   use xtb_mctc_symbols, only : toNumber
   use xtb_mctc_systools
   use xtb_pbc_tools
   use xtb_type_molecule
   use xtb_type_reader
   use xtb_type_vendordata, only : vasp_info
   implicit none
   private

   public :: readHessianDFTBPlus


contains


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
