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

module xtb_io_reader_orca
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   use xtb_mctc_symbols, only : toNumber, symbolLength
   use xtb_type_molecule
   use xtb_pbc_tools
   use xtb_type_reader
   use xtb_readin, only : getline => strip_line
   implicit none
   private

   public :: readHessianOrca


   logical, parameter :: debug = .false.


contains


subroutine readHessianOrca(hessian, reader, mol, status, iomsg)

   real(wp), intent(out) :: hessian(:, :)
   type(TMolecule), intent(in) :: mol
   type(TReader), intent(inout) :: reader
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   character(len=:), allocatable :: line
   character(len=32) :: buffer
   integer :: nline
   real(wp),allocatable :: values(:)
   integer, allocatable :: list(:)
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

   allocate(list(5))
   allocate(values(5))
   call reader%read(line, error)
   iline = iline + 1
   read(line, *, iostat=error) ndim

   if (error /= 0) then
      iomsg = "Could not find dimension of ORCA Hessian"
      return
   end if

   rdlp: do 
      call reader%read(line, error)
      iline = iline + 1
      if (error /= 0) exit rdlp
      call linelen(line,nline)
      list=spread(0,1,nline)
      values=spread(0.0_wp,1,nline)
      if (nline.eq.0) exit rdlp
      read(line, *,iostat=error) list
      if (error /= 0) exit rdlp
      do ii=1, ndim
         call reader%read(line, error)
         iline = iline + 1
         if (error /= 0) exit rdlp
         read(line, *, iostat=error) jj,values
         if (error /= 0) exit rdlp
         hessian(list+1,jj+1) = values
      end do
      if (nline.lt.5) exit rdlp
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

end subroutine readHessianOrca


subroutine linelen(line,nline)
character(len=*), intent(in) :: line
integer, intent(out) :: nline
integer :: istart, iend

istart = 0
iend = 0
nline = 0
do
   istart = verify(line(iend+1:), " ") - 1 + iend
   if (istart < iend) istart = len(line)
   iend = scan(line(istart+1:), " ") - 1 + istart
   if (iend < istart) iend = len(line)
   if (iend == len(line)) exit
   nline = nline +1
end do

end subroutine linelen

end module xtb_io_reader_orca
