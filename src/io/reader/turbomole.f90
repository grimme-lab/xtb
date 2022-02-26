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

module xtb_io_reader_turbomole
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_mctc_convert
   use xtb_mctc_resize
   use xtb_mctc_symbols, only : toNumber, symbolLength
   use xtb_pbc_tools
   use xtb_readin, getline => strip_line
   use xtb_type_molecule
   use xtb_type_reader
   use xtb_type_vendordata, only : turbo_info
   implicit none
   private

   public :: readHessianTurbomole


   logical, parameter :: debug = .false.


contains


subroutine readHessianTurbomole(hessian, reader, mol, status, iomsg)
   real(wp), intent(out) :: hessian(:, :)
   type(TMolecule), intent(in) :: mol
   type(TReader), intent(inout) :: reader
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   character(len=:), allocatable :: line
   character(len=32) :: buffer
   integer :: error, ndim, ii, jj, jbatch, iline

   status = .false.
   iline = 1

   call reader%read(line, error)
   do while(error == 0)
      if (index(line, '$hessian') == 1) exit
      call reader%read(line, error)
      iline = iline + 1
   end do

   if (error /= 0) then
      iomsg = "Could not find $hessian data group"
      return
   end if

   ndim = 3*len(mol)

   rdlp: do ii = 1, ndim
      do jj = 1, ndim, 5
         call reader%read(line, error)
         iline = iline + 1
         if (error /= 0) exit rdlp
         jbatch = min(jj+4, ndim)
         read(line, '(5x, 5f15.10)', iostat=error) hessian(jj:jbatch, ii)
         if (error /= 0) exit rdlp
      end do
   end do rdlp

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

end subroutine readHessianTurbomole


end module xtb_io_reader_turbomole
