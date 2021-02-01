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

!> Implementation of the command line argument parser
module xtb_prog_argparser
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_chartools, only : toLowercase
   implicit none
   private

   public :: TArgParser, init


   !> Command line argument
   type :: TArgument

      !> Has the argument been processed
      logical :: unused

      !> Corresponds to an existing file
      logical :: isFile

      !> Looks like a command line flag
      logical :: isFlag

      !> Raw value of the argument
      character(len=:), allocatable :: raw

   end type TArgument


   !> Command line argument parser
   type :: TArgParser
      private

      !> All command line arguments
      type(TArgument), allocatable :: arg(:)

      !> Start of command line flags
      integer :: flagStart

      !> End of command line flags
      integer :: flagEnd

      !> Position of the flag iterator
      integer :: flagPos

   contains

      !> Basic search of argument list
      procedure, private :: findArg

      !> Basic search of argument list, also markes argument as read
      procedure, private :: getArg

      procedure :: countFlags

      procedure :: countFiles

      procedure :: nextFlag

      procedure :: nextFile

      procedure :: nextArg

      procedure :: reset

   end type TArgParser


   !> Initialize command line argument parser
   interface init
      module procedure :: initArgParser
      module procedure :: initArgument
   end interface init


contains


!> Initialized command line argument
subroutine initArgument(self, iArg)

   !> Command argument
   type(TArgument), intent(out) :: self

   !> Position of argument
   integer, intent(in) :: iArg

   integer :: length
   character(len=*), parameter :: numbers = &
      '0123456789.eE+-'
   character(len=*), parameter :: flagchars = &
      'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-'

   self%unused = .true.

   call get_command_argument(iArg, length=length)
   allocate(character(len=length) :: self%raw)
   call get_command_argument(iArg, self%raw)

   inquire(file=self%raw, exist=self%isFile)
   self%isFlag = index(self%raw, '-') == 1 &  ! starts with a dash
      & .and. verify(self%raw, flagchars) == 0 &  ! looks like a flag
      & .and. verify(self%raw, numbers) > 0  ! but not like a number

end subroutine initArgument


!> Initialize command line argument parser
subroutine initArgParser(self)

   !> Instance of the parser
   type(TArgParser), intent(out) :: self

   integer :: nArg, iArg

   nArg = command_argument_count()
   allocate(self%arg(nArg))

   do iArg = 1, nArg
      call init(self%arg(iArg), iArg)
   end do

   self%flagStart = 1
   iArg = self%getArg('--')
   if (iArg > 0) then
      self%flagEnd = iArg
   else
      self%flagEnd = nArg
   end if

   self%flagPos = 0

end subroutine initArgParser


!> Reset the command line argument parser
subroutine reset(self)

   !> Instance of the parser
   class(TArgParser), intent(inout) :: self

   self%arg%unused = .true.
   self%flagPos = 0

end subroutine reset


!> Count number of unprocessed flags
function countFlags(self) result(nFlags)

   !> Instance of the parser
   class(TArgParser), intent(in) :: self

   !> Number of unprocessed flags
   integer :: nFlags

   nFlags = count(self%arg%unused .and. self%arg%isFlag)

end function countFlags


!> Count number of unprocessed files
function countFiles(self) result(nFiles)

   !> Instance of the parser
   class(TArgParser), intent(in) :: self

   !> Number of unprocessed files
   integer :: nFiles

   nFiles = count(self%arg%unused .and. self%arg%isFile)

end function countFiles


!> Return file name
subroutine nextFlag(self, flag)

   !> Instance of the parser
   class(TArgParser), intent(inout) :: self

   character(len=:), allocatable, intent(out) :: flag

   integer :: iArg

   do iArg = self%flagStart, self%flagEnd
      if (self%arg(iArg)%unused .and. self%arg(iArg)%isFlag) then
         flag = self%arg(iArg)%raw
         self%flagPos = iArg
         self%arg(iArg)%unused = .false.
         exit
      end if
   end do

end subroutine nextFlag


!> Return file name
subroutine nextArg(self, arg)

   !> Instance of the parser
   class(TArgParser), intent(inout) :: self

   character(len=:), allocatable, intent(out) :: arg

   integer :: iArg

   iArg = self%flagPos + 1
   if (iArg <= self%flagEnd) then
      if (self%arg(iArg)%unused .and. .not.self%arg(iArg)%isFlag) then
         self%flagPos = iArg
         arg = self%arg(iArg)%raw
         self%arg(iArg)%unused = .false.
      end if
   end if

end subroutine nextArg


!> Return file name
subroutine nextFile(self, file)

   !> Instance of the parser
   class(TArgParser), intent(inout) :: self

   character(len=:), allocatable, intent(out) :: file

   integer :: iArg

   do iArg = 1, size(self%arg)
      if (self%arg(iArg)%unused .and. self%arg(iArg)%isFile) then
         file = self%arg(iArg)%raw
         self%arg(iArg)%unused = .false.
         exit
      end if
   end do

end subroutine nextFile


!> Basic search of argument list
function findArg(self, arg, iStart, iEnd) result(position)

   !> Instance of the parser
   class(TArgParser), intent(in) :: self

   !> Argument
   character(len=*), intent(in) :: arg

   !> Start position
   integer, intent(in), optional :: iStart

   !> End position
   integer, intent(in), optional :: iEnd

   !> Position of argument
   integer :: position

   integer :: iArg, jStart, jEnd

   if (present(iStart)) then
      jStart = max(1, iStart)
   else
      jStart = 1
   end if

   if (present(iEnd)) then
      jEnd = min(iEnd, size(self%arg))
   else
      jEnd = size(self%arg)
   end if

   position = 0
   do iArg = jStart, jEnd
      if (self%arg(iArg)%unused) then
         if (toLowercase(arg) == toLowercase(self%arg(iArg)%raw)) then
            position = iArg
            exit
         end if
      end if
   end do

end function findArg


!> Basic search of argument list, also markes argument as read
function getArg(self, arg, iStart, iEnd) result(position)

   !> Instance of the parser
   class(TArgParser), intent(inout) :: self

   !> Argument
   character(len=*), intent(in) :: arg

   !> Start position
   integer, intent(in), optional :: iStart

   !> End position
   integer, intent(in), optional :: iEnd

   !> Position of argument
   integer :: position

   position = self%findArg(arg, iStart, iEnd)
   if (position > 0) then
      self%arg(position)%unused = .false.
   end if

end function getArg


end module xtb_prog_argparser
