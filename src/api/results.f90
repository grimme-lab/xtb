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
   use xtb_mctc_convert, only : evtoau
   use xtb_api_environment
   use xtb_api_utils
   use xtb_type_restart
   implicit none
   private

   public :: VResults
   public :: newResults_api, delResults_api, getEnergy_api, getGradient_api, &
      & getPCGradient_api, getVirial_api, getDipole_api, getCharges_api, &
      & getBondOrders_api, getNao_api, getOrbitalEigenvalues_api, &
      & getOrbitalOccupations_api, getOrbitalCoefficients_api, copyResults_api


   !> Void pointer to wavefunction and result objects
   type :: VResults
      type(TRestart), allocatable :: chk
      real(wp), allocatable :: energy
      real(wp), allocatable :: gradient(:, :)
      real(wp), allocatable :: pcgradient(:, :)
      real(wp), allocatable :: sigma(:, :)
      real(wp), allocatable :: dipole(:)
      real(wp), allocatable :: egap
   end type VResults


contains


!> Create new singlepoint results object
function newResults_api() result(vres) &
      & bind(C, name="xtb_newResults")
   !DEC$ ATTRIBUTES DLLEXPORT :: newResults_api
   type(VResults), pointer :: res
   type(c_ptr) :: vres

   call checkGlobalEnv

   allocate(res)
   vres = c_loc(res)

end function newResults_api


!> Delete singlepoint results object
subroutine delResults_api(vres) &
      & bind(C, name="xtb_delResults")
   !DEC$ ATTRIBUTES DLLEXPORT :: delResults_api
   type(c_ptr), intent(inout) :: vres
   type(VResults), pointer :: res

   call checkGlobalEnv

   if (c_associated(vres)) then
      call c_f_pointer(vres, res)
      deallocate(res)
      vres = c_null_ptr
   end if

end subroutine delResults_api


!> Create copy from a singlepoint results object
function copyResults_api(vold) result(vres) &
      & bind(C, name="xtb_copyResults")
   !DEC$ ATTRIBUTES DLLEXPORT :: copyResults_api
   type(VResults), pointer :: res
   type(c_ptr) :: vres
   type(VResults), pointer :: old
   type(c_ptr), value :: vold

   vres = c_null_ptr

   if (c_associated(vold)) then
      call c_f_pointer(vold, old)
      call checkGlobalEnv

      allocate(res, source=old)
      vres = c_loc(res)
   end if

end function copyResults_api


!> Query singlepoint results object for energy
subroutine getEnergy_api(venv, vres, dptr) &
      & bind(C, name="xtb_getEnergy")
   !DEC$ ATTRIBUTES DLLEXPORT :: getEnergy_api
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
   !DEC$ ATTRIBUTES DLLEXPORT :: getGradient_api
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


!> Query singlepoint results object for gradients on external charges
subroutine getPCGradient_api(venv, vres, dptr) &
      & bind(C, name="xtb_getPCGradient")
   !DEC$ ATTRIBUTES DLLEXPORT :: getPCGradient_api
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

      if (.not.allocated(res%pcgradient)) then
         call env%ptr%error("External charge gradients are not available in results", source)
         return
      end if

      dptr(1:3, 1:size(res%pcgradient, 2)) = res%pcgradient

   end if

end subroutine getPCGradient_api


!> Query singlepoint results object for virial
subroutine getVirial_api(venv, vres, dptr) &
      & bind(C, name="xtb_getVirial")
   !DEC$ ATTRIBUTES DLLEXPORT :: getVirial_api
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
   !DEC$ ATTRIBUTES DLLEXPORT :: getDipole_api
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
   !DEC$ ATTRIBUTES DLLEXPORT :: getCharges_api
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

      if (.not.allocated(res%chk)) then
         call env%ptr%error("Partial charges are not available in results", source)
         return
      end if

      dptr(1:size(res%chk%wfn%q)) = res%chk%wfn%q

   end if

end subroutine getCharges_api


!> Query singlepoint results object for bond orders
subroutine getBondOrders_api(venv, vres, dptr) &
      & bind(C, name="xtb_getBondOrders")
   !DEC$ ATTRIBUTES DLLEXPORT :: getBondOrders_api
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

      if (.not.allocated(res%chk)) then
         call env%ptr%error("Bond orders are not available in results", source)
         return
      end if

      dptr(1:size(res%chk%wfn%wbo)) = reshape(res%chk%wfn%wbo, [size(res%chk%wfn%wbo)])

   end if

end subroutine getBondOrders_api


!> Query singlepoint results object for bond orders
subroutine getNao_api(venv, vres, iptr) &
      & bind(C, name="xtb_getNao")
   !DEC$ ATTRIBUTES DLLEXPORT :: getNao_api
   character(len=*), parameter :: source = "xtb_api_getNao"
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vres
   type(VResults), pointer :: res
   integer(c_int), intent(inout) :: iptr

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vres)) then
         call env%ptr%error("Results object is not allocated", source)
         return
      end if
      call c_f_pointer(vres, res)

      if (allocated(res%chk)) then
         iptr = res%chk%wfn%nao
      else
         iptr = 0
      end if

   end if

end subroutine getNao_api


!> Query singlepoint results object for orbital energies
subroutine getOrbitalEigenvalues_api(venv, vres, dptr) &
      & bind(C, name="xtb_getOrbitalEigenvalues")
   !DEC$ ATTRIBUTES DLLEXPORT :: getOrbitalEigenvalues_api
   character(len=*), parameter :: source = "xtb_api_getOrbitalEigenvalues"
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

      if (.not.allocated(res%chk)) then
         call env%ptr%error("Orbital eigenvalues are not available in results", &
            & source)
         return
      end if

      dptr(1:size(res%chk%wfn%emo)) = res%chk%wfn%emo * evtoau

   end if

end subroutine getOrbitalEigenvalues_api


!> Query singlepoint results object for occupation numbers
subroutine getOrbitalOccupations_api(venv, vres, dptr) &
      & bind(C, name="xtb_getOrbitalOccupations")
   !DEC$ ATTRIBUTES DLLEXPORT :: getOrbitalOccupations_api 
   character(len=*), parameter :: source = "xtb_api_getOrbitalOccupations"
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

      if (.not.allocated(res%chk)) then
         call env%ptr%error("Occupation numbers are not available in results", &
            & source)
         return
      end if

      dptr(1:size(res%chk%wfn%focc)) = res%chk%wfn%focc

   end if

end subroutine getOrbitalOccupations_api


!> Query singlepoint results object for orbital coefficients
subroutine getOrbitalCoefficients_api(venv, vres, dptr) &
      & bind(C, name="xtb_getOrbitalCoefficients")
   !DEC$ ATTRIBUTES DLLEXPORT :: getOrbitalCoefficients_api 
   character(len=*), parameter :: source = "xtb_api_getOrbitalCoefficients"
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

      if (.not.allocated(res%chk)) then
         call env%ptr%error("Orbital coefficients are not available in results", &
            & source)
         return
      end if

      dptr(1:size(res%chk%wfn%C)) = reshape(res%chk%wfn%C, [size(res%chk%wfn%C)])

   end if

end subroutine getOrbitalCoefficients_api


end module xtb_api_results
