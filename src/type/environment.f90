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

!> Calculation environment
module xtb_type_environment
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_systools, only : rdvar, rdarg
   implicit none
   private

   public :: TEnvironment, init


   !> Wrapper for characters
   type :: TMessage
      logical :: error
      character(len=:), allocatable :: message
   end type TMessage


   !> Calculation environment
   type :: TEnvironment

      !> Unit for IO
      integer :: unit

      !> Number of logged messages
      integer :: nLog

      !> Message log
      type(TMessage), allocatable :: log(:)

      !> Actual executable name
      character(len=:),allocatable :: whoami

      !> Name of hosting machine
      character(len=:),allocatable :: hostname

      !> User's home directory
      character(len=:),allocatable :: home

      !> Systems path variable
      character(len=:),allocatable :: path

      !> Calculations path variable
      character(len=:),allocatable :: xtbpath

      !> Calculations home directory
      character(len=:),allocatable :: xtbhome

   contains

      !> Add a warning to the message log
      procedure :: warning

      !> Add an error to the message log
      procedure :: error

      !> Check status of environment
      procedure :: check

      !> Check status of environment
      procedure :: checkpoint

   end type TEnvironment


   !> Constructor for the calculation environment
   interface init
      module procedure :: initEnvironment
   end interface


   !> Overloaded resize interface
   interface resize
      module procedure :: resizeMessage
   end interface resize


   !> Initial size of the message log
   integer, parameter :: initialSize = 20


contains


!> Construct the calculation environment
subroutine initEnvironment(self)

   !> Calculation environment
   type(TEnvironment), intent(out) :: self

   integer :: err

   self%unit = stdout
   self%nLog = 0
   allocate(self%log(initialSize))

   call rdarg(0, self%whoami, err)
   call rdvar('HOSTNAME', self%hostname, err)
   call rdvar('HOME', self%home, err)
   call rdvar('PATH', self%path, err)
   call rdvar('XTBHOME', self%xtbhome, err)
   if (err /= 0 .or. len(self%xtbhome) <= 0) then
      self%xtbhome = self%home
   end if
   call rdvar('XTBPATH', self%xtbpath, err)
   if (err /= 0 .or. len(self%xtbpath) <= 0) then
      self%xtbpath = self%xtbhome
   end if

end subroutine initEnvironment


subroutine check(self, terminate)

   !> Calculation environment
   class(TEnvironment), intent(in) :: self

   !> Recommendation for terminating run
   logical, intent(out) :: terminate

   terminate = any(self%log(:self%nLog)%error)

end subroutine check


!> Check status of calculation environment
subroutine checkpoint(self, message, terminate)

   !> Calculation environment
   class(TEnvironment), intent(inout) :: self

   !> Message in case of error
   character(len=*), intent(in) :: message

   !> Recommendation for terminating run
   logical, intent(out) :: terminate

   integer :: iLog

   call self%check(terminate)

   if (self%nLog > 0) then
      if (terminate) then
         write(self%unit, '(/, "[ERROR]", 1x, a)') message
      else
         write(self%unit, '(/, "[WARNING]", 1x, a)') message
      end if
      do iLog = self%nLog, 1, -1
         write(self%unit, '("-", i0, "-", 1x, a)') iLog, self%log(iLog)%message
      end do
   end if

end subroutine checkpoint


!> Create and push back a new error to the message log
subroutine error(self, message, source)

   !> Calculation environment
   class(TEnvironment), intent(inout) :: self

   !> Error message
   character(len=*), intent(in) :: message

   !> Source of the error message
   character(len=*), intent(in), optional :: source

   if (self%nLog >= size(self%log)) then
      call resize(self%log)
   end if

   self%nLog = self%nLog + 1
   if (present(source)) then
      self%log(self%nLog) = TMessage(.true., source // ": " // message)
   else
      self%log(self%nLog) = TMessage(.true., message)
   end if

end subroutine error


!> Create and push back a new error to the message log
subroutine warning(self, message, source)

   !> Calculation environment
   class(TEnvironment), intent(inout) :: self

   !> Error message
   character(len=*), intent(in) :: message

   !> Source of the error message
   character(len=*), intent(in), optional :: source

   if (self%nLog >= size(self%log)) then
      call resize(self%log)
   end if

   self%nLog = self%nLog + 1
   if (present(source)) then
      self%log(self%nLog) = TMessage(.false., source // ": " // message)
   else
      self%log(self%nLog) = TMessage(.false., message)
   end if

end subroutine warning


!> Reallocate message list
subroutine resizeMessage(var, n)

   !> Message list
   type(TMessage), allocatable, intent(inout) :: var(:)

   !> Desired size
   integer, intent(in), optional :: n

   type(TMessage), allocatable :: tmp(:)
   integer :: length, currentLength

   currentLength = size(var)
   if (currentLength > 0) then
      if (present(n)) then
         if (n <= currentLength) return
         length = n
      else
         length = currentLength + currentLength/2 + 1
      endif
      allocate(tmp(length))
      tmp(:currentLength) = var(:currentLength)
      deallocate(var)
      call move_alloc(tmp, var)
   else
      if (present(n)) then
         length = n
      else
         length = 64
      endif
      allocate(var(length))
   endif

end subroutine resizeMessage


end module xtb_type_environment
