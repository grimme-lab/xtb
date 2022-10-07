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

!> 
module xtb_type_iohandler
   use xtb_mctc_filetypes, only : generateFileMetaInfo
   implicit none
   private

   public :: TIOHandler, init


   type :: TFileStatusEnum
      integer :: readin = 1
      integer :: append = 2
      integer :: written = 3
      integer :: replaced = 4
      integer :: deleted = 5
      integer :: created = 6
   end type

   type(TFileStatusEnum), parameter :: fileStatus = TFileStatusEnum()


   type :: TFileHandle
      character(len=:),allocatable :: name
      integer :: status
      integer :: unit
      logical :: open
   end type TFileHandle

   interface TFileHandle
      module procedure :: initFileHandle
   end interface TFileHandle


   type :: TIOHandler

      character(len=:), allocatable :: namespace

      integer :: count
      type(TFileHandle), allocatable :: log(:)

   contains

      procedure :: readFile
      procedure :: writeFile
      procedure :: readBinary
      procedure :: writeBinary
      procedure :: touchFile
      procedure :: closeFile
      procedure :: deleteFile
      procedure :: list

      procedure, private :: getName
      procedure, private :: pushBack
      procedure, private :: findUnit

   end type TIOHandler


   interface init
      module procedure :: initIOHandler
   end interface init


contains


subroutine initIOHandler(self, namespace)
   type(TIOHandler), intent(out) :: self
   character(len=*), intent(in), optional :: namespace

   if (present(namespace)) then
      self%namespace = namespace
   end if
   self%count = 0
   allocate(self%log(20))

end subroutine initIOHandler


elemental function initFileHandle(name, status, unit, open) result(self)
   character(len=*), intent(in) :: name
   integer, intent(in) :: status
   integer, intent(in) :: unit
   logical, intent(in) :: open
   type(TFileHandle) :: self

   self%name = name
   self%status = status
   self%unit = unit
   self%open = open

end function initFileHandle


subroutine getName(self, file, filename)
   class(TIOHandler), intent(inout) :: self
   character(len=*), intent(in) :: file
   character(len=:), allocatable, intent(out) :: filename

   if (allocated(self%namespace)) then
      if (index(file, '/') == 0) then
         if (file(1:1) == '.') then
            filename = '.'//self%namespace//file
         else
            filename = self%namespace//'.'//file
         end if
      else
         filename = file
      end if
   else
      filename = file
   end if

end subroutine getName


subroutine pushBack(self, fileHandle)
   class(TIOHandler), intent(inout) :: self
   type(TFileHandle), intent(in) :: fileHandle
   type(TFileHandle), allocatable :: tmp(:)
   integer :: n

   self%count = self%count + 1
   if (self%count > size(self%log)) then
      n = size(self%log)
      allocate(tmp(n))
      tmp=self%log
      deallocate(self%log)
      allocate(self%log(n + n/2 + 1))
      self%log(1:n) = tmp
      deallocate(tmp)
   end if
   self%log(self%count) = fileHandle

end subroutine pushBack


subroutine findUnit(self, unit, pos)
   class(TIOHandler), intent(in) :: self
   integer, intent(in) :: unit
   integer, intent(out) :: pos
   integer :: i

   pos = 0
   do i = 1, self%count
      if (self%log(i)%open .and. self%log(i)%unit == unit) then
         pos = i
         exit
      end if
   end do

end subroutine findUnit


subroutine list(self, unit)
   class(TIOHandler), intent(in) :: self
   integer, intent(in) :: unit
   character(len=8) :: status
   integer :: i

   if (self%count > 0) then
      write(unit,'(1x,"unit",2x,"open",3x,"action",5x,"filename")')
   end if
   do i = 1, self%count
      select case(self%log(i)%status)
      case(fileStatus%readin);   status = "read    "
      case(fileStatus%written);  status = "written "
      case(fileStatus%deleted);  status = "deleted "
      case(fileStatus%replaced); status = "replaced"
      case(fileStatus%created);  status = "created "
      case default;              status = "unknown "
      end select
      if (self%log(i)%open) then
         write(unit,'(i5,1x,a5,1x,":",1x,a8,3x,a)') &
            & abs(self%log(i)%unit), "true", status, self%log(i)%name
      else
         write(unit,'(i5,1x,a5,1x,":",1x,a8,3x,a)') &
            & abs(self%log(i)%unit), "false", status, self%log(i)%name
      end if
   end do

end subroutine list


subroutine readFile(self, unit, file, iostat)
   class(TIOHandler), intent(inout) :: self
   integer, intent(out) :: unit
   character(len=*), intent(in) :: file
   integer, intent(out), optional :: iostat
   character(len=:), allocatable :: name
   logical :: exist
   integer :: error

   unit = -1
   error = 0

   call self%getName(file, name)
   !$omp critical(io)
   inquire(file=name, exist=exist)

   if (exist) then
      open(file=name, newunit=unit, status='old', action='read', iostat=error)
      if (error == 0) then
         call self%pushBack(TFileHandle(name, fileStatus%readin, unit, .true.))
      else
         unit = -1
      end if
   else
      inquire(file=file, exist=exist)
      if (exist) then
         open(file=file, newunit=unit, status='old', action='read', iostat=error)
         if (error == 0) then
            call self%pushBack(TFileHandle(file, fileStatus%readin, unit, .true.))
         else
            unit = -1
         end if
      else
         error = 1
      end if
   end if
   !$omp end critical(io)

   if (present(iostat)) then
      iostat = error
   end if

end subroutine readFile


subroutine writeFile(self, unit, file, iostat)
   class(TIOHandler), intent(inout) :: self
   integer, intent(out) :: unit
   character(len=*), intent(in) :: file
   integer, intent(out), optional :: iostat
   character(len=:), allocatable :: name
   logical :: exist
   integer :: error

   unit = -1
   error = 0

   call self%getName(file, name)
   !$omp critical(io)
   inquire(file=name, exist=exist)

   open(file=name, newunit=unit, action='write', iostat=error)
   if (error == 0) then
      if (exist) then
         call self%pushBack(TFileHandle(name, fileStatus%replaced, unit, .true.))
      else
         call self%pushBack(TFileHandle(name, fileStatus%written, unit, .true.))
      end if
   else
      unit = -1
   end if
   !$omp end critical(io)

   if (present(iostat)) then
      iostat = error
   end if

end subroutine writeFile


subroutine readBinary(self, unit, file, iostat)
   class(TIOHandler), intent(inout) :: self
   integer, intent(out) :: unit
   character(len=*), intent(in) :: file
   integer, intent(out), optional :: iostat
   character(len=:), allocatable :: name
   logical :: exist
   integer :: error

   unit = -1
   error = 0

   call self%getName(file, name)
   !$omp critical(io)
   inquire(file=name, exist=exist)

   if (exist) then
      open(file=name, newunit=unit, status='old', action='read', &
         & form='unformatted', iostat=error)
      if (error == 0) then
         call self%pushBack(TFileHandle(name, fileStatus%readin, unit, .true.))
      else
         unit = -1
      end if
   else
      inquire(file=file, exist=exist)
      if (exist) then
         open(file=file, newunit=unit, status='old', action='read', &
            & form='unformatted', iostat=error)
         if (error == 0) then
            call self%pushBack(TFileHandle(file, fileStatus%readin, unit, .true.))
         else
            unit = -1
         end if
      else
         error = 1
      end if
   end if
   !$omp end critical(io)

   if (present(iostat)) then
      iostat = error
   end if

end subroutine readBinary


subroutine writeBinary(self, unit, file, iostat)
   class(TIOHandler), intent(inout) :: self
   integer, intent(out) :: unit
   character(len=*), intent(in) :: file
   integer, intent(out), optional :: iostat
   character(len=:), allocatable :: name
   logical :: exist
   integer :: error

   unit = -1
   error = 0

   call self%getName(file, name)
   !$omp critical(io)
   inquire(file=name, exist=exist)

   open(file=name, newunit=unit, action='write', form='unformatted', iostat=error)
   if (error == 0) then
      if (exist) then
         call self%pushBack(TFileHandle(name, fileStatus%replaced, unit, .true.))
      else
         call self%pushBack(TFileHandle(name, fileStatus%written, unit, .true.))
      end if
   else
      unit = -1
   end if
   !$omp end critical(io)

   if (present(iostat)) then
      iostat = error
   end if

end subroutine writeBinary


subroutine touchFile(self, file, iostat)
   class(TIOHandler), intent(inout) :: self
   character(len=*), intent(in) :: file
   integer, intent(out), optional :: iostat
   character(len=:), allocatable :: name
   integer :: unit, error
   logical :: exist

   error = 0

   call self%getName(file, name)
   !$omp critical(io)
   inquire(file=name, exist=exist)
   if (exist) then
      error = 1
   else
      open(file=name, newunit=unit, status='new', iostat=error)
      if (error == 0) then
         call self%pushBack(TFileHandle(name, fileStatus%created, unit, .false.))
         close(unit, iostat=error)
      end if
   end if
   !$omp end critical(io)

   if (present(iostat)) then
      iostat = error
   end if

end subroutine touchFile


subroutine closeFile(self, unit, iostat, remove)
   class(TIOHandler), intent(inout) :: self
   integer, intent(in) :: unit
   integer, intent(out), optional :: iostat
   logical, intent(in), optional :: remove
   integer :: pos, error
   logical :: delete

   if (present(remove)) then
      delete = remove
   else
      delete = .false.
   end if

   error = 0

   call self%findUnit(unit, pos)

   !$omp critical(io)
   if (pos > 0) then
      if (delete) then
         close(unit, iostat=error, status='delete')
      else
         close(unit, iostat=error)
      end if
      if (error == 0) then
         self%log(pos)%open = .false.
         if (delete) then
            self%log(pos)%status = fileStatus%deleted
         end if
      end if
   else
      error = 1
   end if
   !$omp end critical(io)

   if (present(iostat)) then
      iostat = error
   end if

end subroutine closeFile


subroutine deleteFile(self, file, iostat)
   class(TIOHandler), intent(inout) :: self
   character(len=*), intent(in) :: file
   integer, intent(out), optional :: iostat
   character(len=:), allocatable :: name
   integer :: unit
   logical :: exist
   integer :: error

   unit = -1
   error = 0

   call self%getName(file, name)
   !$omp critical(io)
   inquire(file=name, exist=exist)

   if (exist) then
      open(file=name, newunit=unit, status='old', iostat=error)
      if (error == 0) then
         call self%pushBack(TFileHandle(name, fileStatus%deleted, unit, .false.))
         close(unit, status='delete', iostat=error)
      end if
   else
      error = 1
   end if
   !$omp end critical(io)

   if (present(iostat)) then
      iostat = error
   end if

end subroutine deleteFile


end module xtb_type_iohandler
