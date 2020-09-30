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

module xtb_mctc_global
   use xtb_type_environment, only : TEnvironment
#ifdef USE_CUSOLVER
   use cublas
   use cusolverDn
#endif
   character(len=:),allocatable :: name  !< name of the currently running program
   character(len=:),allocatable,target :: msgbuffer !< error message buffer
   integer :: msgid !< number of generated errors
   integer :: maxmsg = 100 !< size of error buffer
   logical :: is_release = .false.
   ! set strict to turn all warnings into hard errors
   logical :: strict = .false.
   ! mute the error messages and generate warnings instead
   logical :: mute = .false.
   ! sanity status of the mctc environment, beware if it goes insane
   logical :: good = .true.

   !> Global environment
   type(TEnvironment), allocatable :: persistentEnv
#ifdef USE_CUSOLVER
   type(cuSolverDnHandle) :: cusolverDnH
#endif

   type :: errormsg
      character(len=:),allocatable :: msg
      integer :: len
   endtype errormsg

   type(errormsg),allocatable :: errorbuffer(:)

contains

subroutine init_errorbuffer
   implicit none
   if (allocated(errorbuffer)) deallocate ( errorbuffer )
   allocate ( errorbuffer(maxmsg) )
   msgid = 0
   msgbuffer = ' '
end subroutine init_errorbuffer

end module xtb_mctc_global

