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

subroutine raise(mode,message)
   use xtb_mctc_global
   implicit none
   character, intent(in) :: mode
   character(len=*), intent(in) :: message

   select case(mode)
   case('S','s') ! save to message buffer
      if (strict) then
         call persistentEnv%error(message)
         call persistentEnv%terminate("Global environment terminated")
      else
         call persistentEnv%warning(message)
      endif
   case('F','f') ! flush message buffer
      call persistentEnv%checkpoint(message)
   case('W','w') ! print warning directly
      if (strict) then
         call persistentEnv%error(message)
         call persistentEnv%terminate("Global environment terminated")
      else
         call persistentEnv%warning(message)
      endif
   case('E','e')
      call persistentEnv%error(message)
      call persistentEnv%terminate("Global environment terminated")
   end select

end subroutine raise

subroutine terminate(signal)
   use xtb_mctc_io, only : stderr
   use xtb_mctc_global, only : name
   integer,intent(in) :: signal
   integer,parameter  :: p_exit_success = 0
   integer,parameter  :: p_exit_failure = 1
   integer,parameter  :: p_exit_external = -1
   if (.not.allocated(name)) name = 'program'
   select case(signal)
   case(p_exit_success)
      write(stderr,'(  "normal termination of",1x,a)') name
      stop
   case(p_exit_external)
      write(stderr,'("external termination of",1x,a)') name
      error stop
   case default
      write(stderr,'("abnormal termination of",1x,a)') name
      error stop
   end select
end subroutine terminate

subroutine xtb_wsigint(signal) bind(c)
   use, intrinsic :: iso_c_binding, only : c_int
   use xtb_mctc_io, only : stderr
   implicit none
   integer(c_int), value :: signal
   write(stderr,'("recieved SIGINT, terminating...")')
   call terminate(-1)
end subroutine xtb_wsigint

subroutine xtb_wsigterm(signal) bind(c)
   use, intrinsic :: iso_c_binding, only : c_int
   use xtb_mctc_io, only : stderr
   implicit none
   integer(c_int), value :: signal
   write(stderr,'("recieved SIGTERM, terminating...")')
   call terminate(-1)
end subroutine xtb_wsigterm
