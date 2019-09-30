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

subroutine raise(mode,message,code)
   use iso_fortran_env
   use mctc_global
   character,       intent(in) :: mode
   character(len=*),intent(in) :: message
   integer,         intent(in) :: code
   select case(mode)
   case('S','s') ! save to message buffer
      if (msgid == maxmsg) then ! check for overflow and flush if necessary
         call flush_warning
      endif
      if (strict) then
         call error
      else
         msgid = msgid + 1
         errorbuffer(msgid) % msg = message
         errorbuffer(msgid) % len = len(message)
      endif
   case('F','f') ! flush message buffer
      if (msgid.gt.0) then
         call flush_warning
      endif
   case('W','w') ! print warning directly
      if (strict) then
         call error
      else
         call warning
      endif
   case('E','e')
      call error
   end select
contains
subroutine flush_warning
   write(output_unit,'(a)')
   write(output_unit,'(72("#"))')
   write(output_unit,'("# WARNING!",1x,a)') message
   do i = 1, msgid
      write(output_unit,'("#  -",1x,a)')  &
         &  errorbuffer(i) % msg
   enddo
   write(output_unit,'(72("#"))')
   write(output_unit,'(a)')
   ! after we are done with printing, clear the buffer
   call init_errorbuffer
end subroutine flush_warning
subroutine warning
   write(output_unit,'("#WARNING!",1x,a)') message
end subroutine warning
subroutine error
   write(output_unit,'("#ERROR!",1x,a)')   message
   call terminate(code)
end subroutine error
end subroutine raise

subroutine terminate(signal)
   use iso_fortran_env, only : error_unit
   use mctc_global, only : name
   integer,intent(in) :: signal
   integer,parameter  :: p_exit_success = 0
   integer,parameter  :: p_exit_failure = 1
   integer,parameter  :: p_exit_external = -1
   if (.not.allocated(name)) name = 'program'
   select case(signal)
   case(p_exit_success)
      write(error_unit,'(  "normal termination of",1x,a)') name
      stop
   case(p_exit_external)
      write(error_unit,'("external termination of",1x,a)') name
      error stop
   case default
      write(error_unit,'("abnormal termination of",1x,a)') name
      error stop
   end select
end subroutine terminate

subroutine wsigint
   use iso_fortran_env, only : error_unit
   write(error_unit,'("recieved SIGINT, terminating...")')
   call terminate(-1)
end subroutine wsigint

subroutine wsigterm
   use iso_fortran_env, only : error_unit
   write(error_unit,'("recieved SIGTERM, terminating...")')
   call terminate(-1)
end subroutine wsigterm
