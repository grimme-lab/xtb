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

module mctc_systools

   character,parameter :: space = ' '
   character,parameter :: colon = ':'
   character,parameter :: slash = '/'

contains

subroutine getline(unit,line,iostat)
   use iso_fortran_env, only : iostat_eor
   integer,intent(in) :: unit
   character(len=:),allocatable,intent(out) :: line
   integer,intent(out),optional :: iostat

   integer,parameter  :: buffersize=144
   character(len=buffersize) :: buffer
   integer :: size
   integer :: err

   line = ''
   do
      read(unit,'(a)',advance='no',iostat=err,size=size)  &
      &    buffer
      if (err.gt.0) then
         if (present(iostat)) iostat=err
         return ! an error occurred
      endif
      line = line // buffer(:size)
      if (err.lt.0) then
         if (err.eq.iostat_eor) err = 0
         if (present(iostat)) iostat=err
         return
      endif
   enddo

end subroutine getline

subroutine rdpath(path,arg,fname,ex)
   implicit none
   character(len=*),intent(in)  :: arg
   character(len=*),intent(in)  :: path
   character(len=:),allocatable,intent(out) :: fname
   logical,intent(out),optional :: ex

!* temporary variables
   character(len=:),allocatable :: scratch1
   character(len=:),allocatable :: scratch2
   character(len=:),allocatable :: fpath
   logical :: exist
   integer :: i

   scratch1 = path
   do
      i = index(scratch1,colon)
      if (i.eq.0) then
         scratch2 = scratch1
      else
         scratch2 = scratch1(:i-1)
         scratch1 = scratch1(i+1:)
      endif
      fpath = scratch2//slash//arg
      inquire(file=fpath,exist=exist)
      if (exist) exit
!     print*,fpath,space,scratch1,exist
      if (i.eq.0) exit
   enddo

!  print*,fpath,exist

   if (exist) fname = fpath
   if (present(ex)) ex = exist
   
end subroutine rdpath

subroutine rdarg(i,arg,iostat)
   integer,intent(in) :: i
   character(len=:),allocatable,intent(out) :: arg
   integer,intent(out),optional :: iostat
   integer :: l,err
   if (allocated(arg)) deallocate(arg)
   call get_command_argument(i,length=l,status=err)
   if (err.ne.0) then
      if (present(iostat)) then
         iostat = err
         return
      else
         call raise('E','Command argument corrupted',1)
      endif
   endif
   allocate( character(len=l) :: arg, stat=err )
   if (err.ne.0) then
      if (present(iostat)) then
         iostat = err
         return
      else
         call raise('E','could not be allocated',1)
      endif
   endif
   call get_command_argument(i,arg,status=err)
   if (err.ne.0) then
      if (present(iostat)) then
         iostat = err
         return
      else
         call raise('E','Command argument corrupted',1)
      endif
   endif
   if (present(iostat)) iostat=0
end subroutine rdarg

subroutine rdvar(name,var,iostat)
   character(len=*),intent(in) :: name
   character(len=:),allocatable,intent(out) :: var
   integer,intent(out),optional :: iostat
   integer :: l,err
   if (allocated(var)) deallocate(var)
   call get_environment_variable(name,length=l,status=err)
   if (err.ne.0) then
      if (present(iostat)) then
         iostat = err
         return
      else
         call raise('E','System variable unassigned',1)
      endif
   endif
   allocate( character(len=l) :: var, stat=err )
   if (err.ne.0) then
      if (present(iostat)) then
         iostat = err
         return
      else
         call raise('E','could not be allocated',1)
      endif
   endif
   call get_environment_variable(name,var,status=err)
   if (err.ne.0) then
      if (present(iostat)) then
         iostat = err
         return
      else
         call raise('E','System variable corrupted',1)
      endif
   endif
   if (present(iostat)) iostat=0
end subroutine rdvar

end module mctc_systools
