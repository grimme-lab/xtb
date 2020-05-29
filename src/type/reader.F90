! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

#ifdef __GNUC__
#if GCC_VERSION < 70500
#define HAS_GCC_BUG_84412
#endif
#endif

!> Wrapper to read from Fortran units
module xtb_type_reader
   implicit none
   private

   public :: TReader


   !> Wrapper for reading from Fortran units
   type :: TReader

      !> Connected unit
      integer :: unit = -1

   contains

      !> Connect reader
      generic :: open => openUnit, openFile

      !> Connect reader to unit
      procedure, private :: openUnit

      !> Connect reader to file
      procedure, private :: openFile

      !> Read from unit
      generic :: read => readLine, readInteger

      !> Read entire record from line
      procedure, private :: readLine

      !> Read an integer from a record
      procedure, private :: readInteger

      !> Close unit in reader
      generic :: close => closeUnit

      !> Close connected unit
      ! Should be private, but must be declared as public due to PGI bug #28452
      procedure :: closeUnit

   end type TReader


contains


!> Connect reader to unit
subroutine openUnit(self, unit)

   !> Instance of reader
   class(TReader), intent(inout) :: self

   !> IO unit connect to
   integer, intent(in) :: unit

#ifdef HAS_GCC_BUG_84412
   self%unit = unit
#else
   logical :: opened

   inquire(opened=opened, unit=unit)
   if (opened) then
      self%unit = unit
   end if
#endif

end subroutine openUnit


!> Connect reader to unit
subroutine openFile(self, file)

   !> Instance of reader
   class(TReader), intent(inout) :: self

   !> IO unit connect t
   character(len=*), intent(in) :: file

   logical :: exist

   inquire(exist=exist, file=file)
   if (exist) then
      open(file=file, newunit=self%unit, status='old')
   end if

end subroutine openFile


!> Read from unit
subroutine readLine(self, line, iostat)

   !> Instance of reader
   class(TReader), intent(inout) :: self

   !> Character buffer to read into
   character(len=:), allocatable, intent(out) :: line

   !> Status of read
   integer, intent(out) :: iostat

   integer, parameter :: bufferSize = 512
   character(len=bufferSize) :: buffer
   integer :: size

   line = ''
   do
      read(self%unit, '(a)', advance='no', iostat=iostat, size=size) buffer
      if (iostat > 0) then
         exit
      end if
      line = line // buffer(:size)
      if (iostat < 0) then
         if (is_iostat_eor(iostat)) then
            iostat = 0
         end if
         exit
      end if
   end do

end subroutine readLine


subroutine readInteger(self, val, iostat)

   !> Instance of reader
   class(TReader), intent(inout) :: self

   !> Value of integer read
   integer, intent(out) :: val

   !> Status of read
   integer, intent(out) :: iostat

   read(self%unit, *, iostat=iostat) val

end subroutine readInteger


!> Close unit in reader
subroutine closeUnit(self)

   !> Instance of reader
   class(TReader), intent(inout) :: self

#ifdef HAS_GCC_BUG_84412
   close(self%unit)
#else
   logical :: opened

   inquire(opened=opened, unit=self%unit)
   if (opened) then
      close(self%unit)
   end if
#endif

end subroutine closeUnit


end module xtb_type_reader
