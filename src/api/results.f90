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

!> API for dealing with the calculation results
module xtb_api_results
   use, intrinsic :: iso_c_binding
   use xtb_mctc_accuracy, only : wp
   use xtb_api_environment
   use xtb_api_utils
   use xtb_type_wavefunction
   implicit none
   private

   public :: VResults
   public :: newResults_api, delResults_api, getEnergy_api, getGradient_api, &
      & getVirial_api, getDipole_api, getCharges_api, getBondOrders_api


   !> Void pointer to wavefunction and result objects
   type :: VResults
      type(TWavefunction), allocatable :: wfn
      real(wp), allocatable :: energy
      real(wp), allocatable :: gradient(:, :)
      real(wp), allocatable :: sigma(:, :)
      real(wp), allocatable :: dipole(:)
      real(wp), allocatable :: egap
   end type VResults


contains


!> Create new singlepoint results object
function newResults_api() result(vres) &
      & bind(C, name="xtb_newResults")
   type(VResults), pointer :: res
   type(c_ptr) :: vres

   call checkGlobalEnv

   allocate(res)
   vres = c_loc(res)

end function newResults_api


!> Delete singlepoint results object
subroutine delResults_api(vres) &
      & bind(C, name="xtb_delResults")
   type(c_ptr), intent(inout) :: vres
   type(VResults), pointer :: res

   call checkGlobalEnv

   if (c_associated(vres)) then
      call c_f_pointer(vres, res)
      deallocate(res)
      vres = c_null_ptr
   end if

end subroutine delResults_api


!> Query singlepoint results object for energy
subroutine getEnergy_api(venv, vres, dptr) &
      & bind(C, name="xtb_getEnergy")
   character(len=*), parameter :: source = "xtb_api_getEnergy"
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vres
   type(VResults), pointer :: res
   real(c_double), intent(inout) :: dptr

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vres)) then
         call env%ptr%error("Results object is not allocated", source)
         return
      end if
      call c_f_pointer(vres, res)

      if (.not.allocated(res%energy)) then
         call env%ptr%error("Energy is not available in results", source)
         return
      end if

      dptr = res%energy

   end if

end subroutine getEnergy_api


!> Query singlepoint results object for gradient
subroutine getGradient_api(venv, vres, dptr) &
      & bind(C, name="xtb_getGradient")
   character(len=*), parameter :: source = "xtb_api_getGradient"
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vres
   type(VResults), pointer :: res
   real(c_double), intent(inout) :: dptr(3, *)

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vres)) then
         call env%ptr%error("Results object is not allocated", source)
         return
      end if
      call c_f_pointer(vres, res)

      if (.not.allocated(res%gradient)) then
         call env%ptr%error("Gradient is not available in results", source)
         return
      end if

      dptr(1:3, 1:size(res%gradient, 2)) = res%gradient

   end if

end subroutine getGradient_api


!> Query singlepoint results object for virial
subroutine getVirial_api(venv, vres, dptr) &
      & bind(C, name="xtb_getVirial")
   character(len=*), parameter :: source = "xtb_api_getVirial"
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vres
   type(VResults), pointer :: res
   real(c_double), intent(inout) :: dptr(3, 3)

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vres)) then
         call env%ptr%error("Results object is not allocated", source)
         return
      end if
      call c_f_pointer(vres, res)

      if (.not.allocated(res%sigma)) then
         call env%ptr%error("Virial is not available in results", source)
         return
      end if

      dptr(1:3, 1:3) = res%sigma

   end if

end subroutine getVirial_api


!> Query singlepoint results object for dipole moment
subroutine getDipole_api(venv, vres, dptr) &
      & bind(C, name="xtb_getDipole")
   character(len=*), parameter :: source = "xtb_api_getDipole"
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vres
   type(VResults), pointer :: res
   real(c_double), intent(inout) :: dptr(3)

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vres)) then
         call env%ptr%error("Results object is not allocated", source)
         return
      end if
      call c_f_pointer(vres, res)

      if (.not.allocated(res%dipole)) then
         call env%ptr%error("Dipole moment is not available in results", source)
         return
      end if

      dptr(1:3) = res%dipole

   end if

end subroutine getDipole_api


!> Query singlepoint results object for partial charges
subroutine getCharges_api(venv, vres, dptr) &
      & bind(C, name="xtb_getCharges")
   character(len=*), parameter :: source = "xtb_api_getCharges"
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vres
   type(VResults), pointer :: res
   real(c_double), intent(inout) :: dptr(*)

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vres)) then
         call env%ptr%error("Results object is not allocated", source)
         return
      end if
      call c_f_pointer(vres, res)

      if (.not.allocated(res%wfn)) then
         call env%ptr%error("Partial charges are not available in results", source)
         return
      end if

      dptr(1:size(res%wfn%q)) = res%wfn%q

   end if

end subroutine getCharges_api


!> Query singlepoint results object for bond orders
subroutine getBondOrders_api(venv, vres, dptr) &
      & bind(C, name="xtb_getBondOrders")
   character(len=*), parameter :: source = "xtb_api_getBondOrders"
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vres
   type(VResults), pointer :: res
   real(c_double), intent(inout) :: dptr(*)

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vres)) then
         call env%ptr%error("Results object is not allocated", source)
         return
      end if
      call c_f_pointer(vres, res)

      if (.not.allocated(res%wfn)) then
         call env%ptr%error("Bond orders are not available in results", source)
         return
      end if

      dptr(1:size(res%wfn%wbo)) = reshape(res%wfn%wbo, [size(res%wfn%wbo)])

   end if

end subroutine getBondOrders_api


end module xtb_api_results
