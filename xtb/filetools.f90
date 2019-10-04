! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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

subroutine print_filelist(iunit)
   use mctc_filetools
   use readin
   implicit none
   integer,intent(in) :: iunit
   character(len=8) :: status
   integer :: i

   write(iunit,'(1x,"unit",2x,"open",3x,"action",5x,"filename")')
   do i = 1, nfiles
      select case(filelist(i)%status)
      case('r','R'); status = "read    "
      case('a','A'); status = "appended"
      case('w','W'); status = "written "
      case('s','S'); status = "replaced"
      case('d','D'); status = "deleted "
      case('t','T'); status = "touched "
      case default;  status = "        "
      end select
      write(iunit,'(i5,1x,a5,1x,":",1x,a8,3x,a)') &
         abs(filelist(i)%unit), bool2string(filelist(i)%open), status, &
         filelist(i)%name
   enddo
end subroutine print_filelist

subroutine delete_file(name)
   use mctc_filetools
   use setparam
   implicit none
   character(len=*),intent(in)  :: name
   character(len=:),allocatable :: file
   logical :: exist,opened
   integer :: iunit,err

   file = get_namespace(name)

   !$omp critical (io)
   inquire(file=file, exist=exist, opened=opened, number=iunit)
   if (opened) then
      close(iunit,status='delete')
      call pop_file(iunit,'d')
   else if (exist) then
      inquire(file=file, opened=opened, number=iunit)
      open(newunit=iunit, file=file, iostat=err, status='old')
      if (err.eq.0) then
         close(iunit,status='delete')
         call push_file(iunit,file,'d')
      endif
   endif
   !$omp end critical (io)
end subroutine delete_file

subroutine touch_file(name)
   use mctc_filetools
   use setparam
   implicit none
   character(len=*),intent(in)  :: name
   character(len=:),allocatable :: file
   logical :: exist
   integer :: iunit,err

   file = get_namespace(name)

   !$omp critical (io)
   inquire(file=file, exist=exist)
   if (.not.exist) then
      open(newunit=iunit, file=file, iostat=err, status='new')
      if (err.eq.0) then
         close(iunit)
         call push_file(iunit,file,'t')
      endif
   endif
   !$omp end critical (io)
end subroutine touch_file

subroutine open_binary(iunit,name,status)
   use mctc_filetools
   use setparam
   implicit none
   character(len=*),intent(in)  :: name
   character(len=1),intent(in)  :: status
   character(len=:),allocatable :: file
   integer,intent(out) :: iunit
   character(len=1)    :: cstatus
   integer :: err
   logical :: exist

   iunit = -1 ! failed to open file
   cstatus = status

   file = get_namespace(name)

   ! acquire mutex for IO
   !$omp critical (io)
   select case(status)
   case default
      err = -1
   case('a','A')
      inquire(file=file, exist=exist)
      if (exist) then
         open(newunit=iunit, file=file, iostat=err, action='write', status='old',&
            form='unformatted')
         if (err.ne.0) iunit = -1
      else
         open(newunit=iunit, file=file, iostat=err, action='write', status='new',&
            form='unformatted')
         cstatus = 'w'
         if (err.ne.0) iunit = -1
      endif
   case('r','R')
      inquire(file=file, exist=exist)
      if (.not.exist) then
         file = name
         inquire(file=file, exist=exist)
      endif
      if (exist) then
         open(newunit=iunit, file=file, iostat=err, action='read', status='old',&
            form='unformatted')
         if (err.ne.0) iunit = -1
      endif
   case('w','W')
      inquire(file=file, exist=exist)
      if (exist) cstatus = 's'
      open(newunit=iunit, file=file, iostat=err, action='write',&
         form='unformatted')
      if (err.ne.0) iunit = -1
   end select

   if (iunit.ne.-1) then
      call push_file(iunit,file,cstatus)
   endif
   !$omp end critical (io)

end subroutine open_binary

subroutine open_file(iunit,name,status)
   use mctc_filetools
   use setparam
   implicit none
   character(len=*),intent(in)  :: name
   character(len=1),intent(in)  :: status
   character(len=:),allocatable :: file
   integer,intent(out) :: iunit
   character(len=1)    :: cstatus
   integer :: err
   logical :: exist

   iunit = -1 ! failed to open file
   cstatus = status

   file = get_namespace(name)

   ! acquire mutex for IO
   !$omp critical (io)
   select case(status)
   case default
      err = -1
   case('a','A')
      inquire(file=file, exist=exist)
      if (exist) then
         open(newunit=iunit, file=file, iostat=err, action='write', status='old')
         if (err.ne.0) iunit = -1
      else
         open(newunit=iunit, file=file, iostat=err, action='write', status='new')
         cstatus = 'w'
         if (err.ne.0) iunit = -1
      endif
   case('r','R')
      inquire(file=file, exist=exist)
      if (.not.exist) then
         file = name
         inquire(file=file, exist=exist)
      endif
      if (exist) then
         open(newunit=iunit, file=file, iostat=err, action='read', status='old')
         if (err.ne.0) iunit = -1
      endif
   case('w','W')
      inquire(file=file, exist=exist)
      if (exist) cstatus = 's'
      open(newunit=iunit, file=file, iostat=err, action='write')
      if (err.ne.0) iunit = -1
   end select

   if (iunit.ne.-1) then
      call push_file(iunit,file,cstatus)
   endif
   !$omp end critical (io)

end subroutine open_file

subroutine close_file(unit)
   use mctc_filetools
   implicit none
   integer,intent(in) :: unit
   logical :: opened
   !$omp critical (io)
   inquire(unit=unit,opened=opened)
   if (opened) then
      close(unit)
      call pop_file(unit)
   endif
   !$omp end critical (io)
end subroutine close_file

subroutine remove_file(unit)
   use mctc_filetools
   implicit none
   integer,intent(in) :: unit
   logical :: opened
   !$omp critical (io)
   inquire(unit=unit,opened=opened)
   if (opened) then
      close(unit,status='delete')
      call pop_file(unit,'d')
   endif
   !$omp end critical (io)
end subroutine remove_file
