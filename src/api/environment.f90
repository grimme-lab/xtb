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

!> API for dealing with the calculation environment
module xtb_api_environment
   use, intrinsic :: iso_c_binding
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io
   use xtb_api_utils
   use xtb_type_environment
   use xtb_xtb_calculator
   implicit none
   private

   public :: VEnvironment
   public :: newEnvironment_api, delEnvironment_api, checkEnvironment_api, &
      & showEnvironment_api, setOutput_api, releaseOutput_api, setVerbosity_api, &
      & getError_api


   !> Void pointer to calculation environment
   type :: VEnvironment
      type(TEnvironment) :: ptr
      integer :: verbosity
   end type VEnvironment


contains


!> Create new xtb calculation environment object
function newEnvironment_api() result(venv) &
      & bind(C, name="xtb_newEnvironment")
   !DEC$ ATTRIBUTES DLLEXPORT :: newEnvironment_api
   type(VEnvironment), pointer :: env
   type(c_ptr) :: venv

   call checkGlobalEnv

   allocate(env)
   call init(env%ptr)
   env%verbosity = 1
   venv = c_loc(env)

end function newEnvironment_api


!> Delete a xtb calculation environment object
subroutine delEnvironment_api(venv) &
      & bind(C, name="xtb_delEnvironment")
   !DEC$ ATTRIBUTES DLLEXPORT :: delEnvironment_api
   type(c_ptr), intent(inout) :: venv
   type(VEnvironment), pointer :: env

   call releaseOutput_api(venv)

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv
      deallocate(env)
      venv = c_null_ptr
   end if

end subroutine delEnvironment_api


!> Check current status of calculation environment
function checkEnvironment_api(venv) result(status) &
      & bind(C, name="xtb_checkEnvironment")
   !DEC$ ATTRIBUTES DLLEXPORT :: checkEnvironment_api
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   integer(c_int) :: status
   logical :: exitRun

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      call env%ptr%check(exitRun)
      if (exitRun) then
         status = 1
      else
         status = 0
      end if
   end if

end function checkEnvironment_api


!> Show and empty error stack
subroutine showEnvironment_api(venv, charptr) &
      & bind(C, name="xtb_showEnvironment")
   !DEC$ ATTRIBUTES DLLEXPORT :: showEnvironment_api
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   character(kind=c_char), intent(in), optional :: charptr(*)
   character(len=:, kind=c_char), allocatable :: message

   if (present(charptr)) then
      call c_f_character(charptr, message)
   else
      message = ''
   end if

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      call env%ptr%show(message)
   end if

end subroutine showEnvironment_api


!> Bind output from this environment
subroutine setOutput_api(venv, charptr) &
      & bind(C, name="xtb_setOutput")
   !DEC$ ATTRIBUTES DLLEXPORT :: setOutput_api
   character(len=*), parameter :: source = 'xtb_api_setOutput'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   character(kind=c_char), intent(in) :: charptr(*)
   character(len=:, kind=c_char), allocatable :: filename
   integer :: stat
   character(len=512) :: message

   call c_f_character(charptr, filename)

   call releaseOutput_api(venv)

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (filename == c_char_'-' .or. filename == c_char_'STDOUT') then
         env%ptr%unit = stdout
      else if (filename == c_char_'STDERR') then
         env%ptr%unit = stderr
      else
         open(file=filename, newunit=env%ptr%unit, iostat=stat, iomsg=message)
         if (stat /= 0) then
            call env%ptr%error(trim(message), source)
         end if
      end if
   end if

end subroutine setOutput_api


!> Release output unit from this environment
subroutine releaseOutput_api(venv) &
      & bind(C, name="xtb_releaseOutput")
   !DEC$ ATTRIBUTES DLLEXPORT :: releaseOutput_api
   character(len=*), parameter :: source = 'xtb_api_releaseOutput'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   integer :: stat
   character(len=512) :: message

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.any(env%ptr%unit /= [-1, stdout, stderr])) then
         close(unit=env%ptr%unit, iostat=stat, iomsg=message)
         if (stat /= 0) then
            call env%ptr%error(trim(message), source)
         end if
      end if
   end if

end subroutine releaseOutput_api


!> Set verbosity of calculation output
subroutine setVerbosity_api(venv, verbosity) &
      & bind(C, name="xtb_setVerbosity")
   !DEC$ ATTRIBUTES DLLEXPORT :: setVerbosity_api
   character(len=*), parameter :: source = 'xtb_api_setVerbosity'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   integer, value :: verbosity

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.any(verbosity /= [0_c_int, 1_c_int, 2_c_int])) then
         call env%ptr%warning("Verbosity specified is out of range", source)
      end if

      env%verbosity = verbosity
   end if

end subroutine setVerbosity_api


!> Set verbosity of calculation output
subroutine getError_api(venv, charptr, buffersize) &
      & bind(C, name="xtb_getError")
   !DEC$ ATTRIBUTES DLLEXPORT :: setVerbosity_api, getError_api
   character(len=*), parameter :: source = 'xtb_api_getError'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   character(kind=c_char), intent(inout) :: charptr(*)
   character(len=:, kind=c_char), allocatable :: buffer
   integer(c_int), intent(in), optional :: buffersize
   logical :: exitRun
   integer :: maxLength

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (present(buffersize)) then
         maxLength = buffersize
      else
         maxLength = huge(maxlength) - 2
      end if

      call env%ptr%check(exitRun)
      if (exitRun) then
         call env%ptr%getLog(buffer)
         if (len(buffer) > 0) then
            call f_c_character(buffer, charptr, maxLength)
         end if
      end if
   end if

end subroutine getError_api


end module xtb_api_environment
