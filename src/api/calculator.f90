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

!> API for dealing with the single point calculator
module xtb_api_calculator
   use, intrinsic :: iso_c_binding
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io
   use xtb_mctc_systools, only : rdpath
   use xtb_api_environment
   use xtb_api_molecule
   use xtb_api_utils
   use xtb_gfnff_calculator
   use xtb_type_pcem
   use xtb_main_setup
   use xtb_solv_kernel
   use xtb_solv_input
   use xtb_solv_state
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_calculator
   use xtb_xtb_calculator
   implicit none
   private

   public :: VCalculator
   public :: newCalculator_api, delCalculator_api
   public :: loadGFNFF_api, loadGFN0xTB_api, loadGFN1xTB_api, loadGFN2xTB_api
   public :: setSolvent_api, releaseSolvent_api
   public :: setExternalCharges_api, releaseExternalCharges_api
   public :: setAccuracy_api, setElectronicTemp_api, setMaxIter_api


   !> Void pointer to single point calculator
   type :: VCalculator
      class(TCalculator), allocatable :: ptr
   end type VCalculator


contains


function newCalculator_api() result(vcalc) &
      & bind(C, name="xtb_newCalculator")
   !DEC$ ATTRIBUTES DLLEXPORT :: newCalculator_api
   type(VCalculator), pointer :: calc
   type(c_ptr) :: vcalc

   call checkGlobalEnv

   allocate(calc)
   vcalc = c_loc(calc)

end function newCalculator_api


subroutine delCalculator_api(vcalc) &
      & bind(C, name="xtb_delCalculator")
   !DEC$ ATTRIBUTES DLLEXPORT :: delCalculator_api
   type(c_ptr), intent(inout) :: vcalc
   type(VCalculator), pointer :: calc

   call checkGlobalEnv

   if (c_associated(vcalc)) then
      call c_f_pointer(vcalc, calc)
      deallocate(calc)
      vcalc = c_null_ptr
   end if

end subroutine delCalculator_api


subroutine loadGFNFF_api(venv, vmol, vcalc, charptr) &
      & bind(C, name="xtb_loadGFNFF")
   !DEC$ ATTRIBUTES DLLEXPORT :: loadGFNFF_api
   character(len=*), parameter :: source = 'xtb_api_loadGFNFF'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vmol
   type(VMolecule), pointer :: mol
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc
   character(kind=c_char), intent(in), optional :: charptr(*)
   character(len=:, kind=c_char), allocatable :: dummy, filename
   character(len=*), parameter :: pFilename = '.param_gfnff.xtb'
   type(TGFFCalculator), allocatable :: gff
   logical :: exist, exitRun

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vmol)) then
         call env%ptr%error("Molecular structure data is not allocated", source)
         return
      end if
      call c_f_pointer(vmol, mol)

      if (.not.c_associated(vcalc)) then
         call env%ptr%error("Singlepoint calculator is not allocated", source)
         return
      end if
      call c_f_pointer(vcalc, calc)

      if (present(charptr)) then
         call c_f_character(charptr, dummy)
      else
         dummy = pFilename
      end if
      inquire(file=dummy, exist=exist)
      if (.not.exist) then
         call rdpath(env%ptr%xtbpath, dummy, filename, exist)
         if (.not.exist) then
            filename = dummy
         end if
      else
         filename = dummy
      end if

      allocate(gff)
      call newGFFCalculator(env%ptr, mol%ptr, gff, filename, .false.)

      call env%ptr%check(exitRun)
      if (exitRun) then
         call env%ptr%error("Could not construct GFN-FF calculator", source)
         return
      end if

      call move_alloc(gff, calc%ptr)

   end if

end subroutine loadGFNFF_api


subroutine loadGFN0xTB_api(venv, vmol, vcalc, charptr) &
      & bind(C, name="xtb_loadGFN0xTB")
   !DEC$ ATTRIBUTES DLLEXPORT :: loadGFN0xTB_api
   character(len=*), parameter :: source = 'xtb_api_loadGFN0xTB'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vmol
   type(VMolecule), pointer :: mol
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc
   character(kind=c_char), intent(in), optional :: charptr(*)
   character(len=:, kind=c_char), allocatable :: dummy, filename
   character(len=*), parameter :: pFilename = 'param_gfn0-xtb.txt'
   type(TxTBCalculator), allocatable :: xtb
   logical :: exist, exitRun

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vmol)) then
         call env%ptr%error("Molecular structure data is not allocated", source)
         return
      end if
      call c_f_pointer(vmol, mol)

      if (.not.c_associated(vcalc)) then
         call env%ptr%error("Singlepoint calculator is not allocated", source)
         return
      end if
      call c_f_pointer(vcalc, calc)

      if (present(charptr)) then
         call c_f_character(charptr, dummy)
      else
         dummy = pFilename
      end if
      inquire(file=dummy, exist=exist)
      if (.not.exist) then
         call rdpath(env%ptr%xtbpath, dummy, filename, exist)
         if (.not.exist) then
            filename = dummy
         end if
      else
         filename = dummy
      end if

      allocate(xtb)
      call newXTBCalculator(env%ptr, mol%ptr, xtb, filename)

      call env%ptr%check(exitRun)
      if (exitRun) then
         call env%ptr%error("Could not construct GFN0-xTB calculator", source)
         return
      end if

      call move_alloc(xtb, calc%ptr)

   end if

end subroutine loadGFN0xTB_api


subroutine loadGFN1xTB_api(venv, vmol, vcalc, charptr) &
      & bind(C, name="xtb_loadGFN1xTB")
   !DEC$ ATTRIBUTES DLLEXPORT :: loadGFN1xTB_api
   character(len=*), parameter :: source = 'xtb_api_loadGFN1xTB'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vmol
   type(VMolecule), pointer :: mol
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc
   character(kind=c_char), intent(in), optional :: charptr(*)
   character(len=:, kind=c_char), allocatable :: dummy, filename
   character(len=*), parameter :: pFilename = 'param_gfn1-xtb.txt'
   type(TxTBCalculator), allocatable :: xtb
   logical :: exist, exitRun

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vmol)) then
         call env%ptr%error("Molecular structure data is not allocated", source)
         return
      end if
      call c_f_pointer(vmol, mol)

      if (.not.c_associated(vcalc)) then
         call env%ptr%error("Singlepoint calculator is not allocated", source)
         return
      end if
      call c_f_pointer(vcalc, calc)

      if (present(charptr)) then
         call c_f_character(charptr, dummy)
      else
         dummy = pFilename
      end if
      inquire(file=dummy, exist=exist)
      if (.not.exist) then
         call rdpath(env%ptr%xtbpath, dummy, filename, exist)
         if (.not.exist) then
            filename = dummy
         end if
      else
         filename = dummy
      end if

      allocate(xtb)
      call newXTBCalculator(env%ptr, mol%ptr, xtb, filename)

      call env%ptr%check(exitRun)
      if (exitRun) then
         call env%ptr%error("Could not construct GFN1-xTB calculator", source)
         return
      end if

      call move_alloc(xtb, calc%ptr)

   end if

end subroutine loadGFN1xTB_api


subroutine loadGFN2xTB_api(venv, vmol, vcalc, charptr) &
      & bind(C, name="xtb_loadGFN2xTB")
   !DEC$ ATTRIBUTES DLLEXPORT :: loadGFN2xTB_api
   character(len=*), parameter :: source = 'xtb_api_loadGFN2xTB'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vmol
   type(VMolecule), pointer :: mol
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc
   character(kind=c_char), intent(in), optional :: charptr(*)
   character(len=:, kind=c_char), allocatable :: dummy, filename
   character(len=*), parameter :: pFilename = 'param_gfn2-xtb.txt'
   type(TxTBCalculator), allocatable :: xtb
   logical :: exist, exitRun

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vmol)) then
         call env%ptr%error("Molecular structure data is not allocated", source)
         return
      end if
      call c_f_pointer(vmol, mol)

      if (.not.c_associated(vcalc)) then
         call env%ptr%error("Singlepoint calculator is not allocated", source)
         return
      end if
      call c_f_pointer(vcalc, calc)

      if (present(charptr)) then
         call c_f_character(charptr, dummy)
      else
         dummy = pFilename
      end if
      inquire(file=dummy, exist=exist)
      if (.not.exist) then
         call rdpath(env%ptr%xtbpath, dummy, filename, exist)
         if (.not.exist) then
            filename = dummy
         end if
      else
         filename = dummy
      end if

      allocate(xtb)
      call newXTBCalculator(env%ptr, mol%ptr, xtb, filename)

      call env%ptr%check(exitRun)
      if (exitRun) then
         call env%ptr%error("Could not construct GFN2-xTB calculator", source)
         return
      end if

      call move_alloc(xtb, calc%ptr)

   end if

end subroutine loadGFN2xTB_api


!> Add a solvation model to calculator (requires loaded parametrisation)
subroutine setSolvent_api(venv, vcalc, charptr, state, temperature, grid) &
      & bind(C, name="xtb_setSolvent")
   !DEC$ ATTRIBUTES DLLEXPORT :: setSolvent_api
   character(len=*), parameter :: source = 'xtb_api_setSolvent'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc
   character(kind=c_char), intent(in) :: charptr(*)
   integer(c_int), intent(in), optional :: state
   real(c_double), intent(in), optional :: temperature
   integer(c_int), intent(in), optional :: grid
   character(len=:), allocatable :: solvent
   type(TSolvInput) :: input
   integer :: gsolvstate, nang
   real(wp) :: temp
   logical :: exitRun

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vcalc)) then
         call env%ptr%error("Singlepoint calculator is not allocated", source)
         return
      end if
      call c_f_pointer(vcalc, calc)

      if (.not.allocated(calc%ptr)) then
         call env%ptr%error("No calculator loaded to add solvation model", source)
         return
      end if

      if (present(state)) then
         gsolvstate = state
      else
         gsolvstate = solutionState%gsolv
      end if

      if (present(temperature)) then
         temp = temperature
      else
         temp = 298.15_wp
      end if

      if (present(grid)) then
         nang = grid
      else
         nang = 230
      end if

      call c_f_character(charptr, solvent)

      ! PGI 20.5 cannot use default constructor with deferred-length characters:
      ! input = TSolvInput(solvent=solvent, temperature=temp, state=gsolvstate, &
      !    & nang=nang)
      input%solvent = solvent
      input%temperature = temp
      input%state = gsolvstate
      input%nang = nang
      input%alpb = .false.
      input%kernel = gbKernel%still
      call addSolvationModel(env%ptr, calc%ptr, input)

      call env%ptr%check(exitRun)
      if (exitRun) then
         call env%ptr%error("Could not add solvation model for '"//solvent//"'", &
            & source)
         return
      end if

   end if

end subroutine setSolvent_api


!> Unset the solvation model
subroutine releaseSolvent_api(venv, vcalc) &
      & bind(C, name="xtb_releaseSolvent")
   !DEC$ ATTRIBUTES DLLEXPORT :: releaseSolvent_api
   character(len=*), parameter :: source = 'xtb_api_setRelease'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vcalc)) then
         call env%ptr%error("Singlepoint calculator is not allocated", source)
         return
      end if
      call c_f_pointer(vcalc, calc)

      if (allocated(calc%ptr)) then
         calc%ptr%lSolv = .false.
         if (allocated(calc%ptr%solvation)) then
            deallocate(calc%ptr%solvation)
         end if
      end if

   end if

end subroutine releaseSolvent_api


!> Add a external charge potential to calculator (only supported in GFN1/2-xTB)
subroutine setExternalCharges_api(venv, vcalc, npc, numbers, charges, positions) &
      & bind(C, name="xtb_setExternalCharges")
   !DEC$ ATTRIBUTES DLLEXPORT :: setExternalCharges_api
   character(len=*), parameter :: source = 'xtb_api_setExternalCharges'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc
   integer(c_int), intent(in) :: npc
   integer(c_int), intent(in) :: numbers(*)
   real(c_double), intent(in) :: charges(*)
   real(c_double), intent(in) :: positions(3, *)
   integer :: ii

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (npc <= 0) then
         call env%ptr%error("Negative number of point charges provided", source)
         return
      end if

      if (.not.c_associated(vcalc)) then
         call env%ptr%error("Singlepoint calculator is not allocated", source)
         return
      end if
      call c_f_pointer(vcalc, calc)

      if (.not.allocated(calc%ptr)) then
         call env%ptr%error("No calculator loaded to add external potential", &
            & source)
         return
      end if

      select type(xtb => calc%ptr)
      class default
         call env%ptr%error("Calculator does not support external potentials", &
            & source)
         return

      type is(TxTBCalculator)
         if (xtb%xtbData%level == 0) then
            call env%ptr%error("GFN0-xTB does not support external potentials", &
               & source)
            return
         end if

         call xtb%pcem%allocate(npc)
         do ii = 1, npc
            xtb%pcem%xyz(:, ii) = positions(:, ii)
            xtb%pcem%gam(ii) = xtb%xtbData%coulomb%chemicalHardness(numbers(ii))
            xtb%pcem%q(ii) = charges(ii)
         end do

      end select

   end if

end subroutine setExternalCharges_api


!> Unset the external charge potential
subroutine releaseExternalCharges_api(venv, vcalc) &
      & bind(C, name="xtb_releaseExternalCharges")
   !DEC$ ATTRIBUTES DLLEXPORT :: releaseExternalCharges_api
   character(len=*), parameter :: source = 'xtb_api_releaseExternalCharges'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vcalc)) then
         call env%ptr%error("Singlepoint calculator is not allocated", source)
         return
      end if
      call c_f_pointer(vcalc, calc)

      ! There are only limited cases when we need to perform any action.
      ! Deconstructors should behave nicely, therefore no errors on wrong input.
      if (allocated(calc%ptr)) then
         select type(xtb => calc%ptr)
         type is(TxTBCalculator)
            call xtb%pcem%deallocate
         end select
      end if

   end if

end subroutine releaseExternalCharges_api


subroutine setAccuracy_api(venv, vcalc, accuracy) &
      & bind(C, name="xtb_setAccuracy")
   !DEC$ ATTRIBUTES DLLEXPORT :: setAccuracy_api
   character(len=*), parameter :: source = 'xtb_api_setAccuracy'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc
   real(c_double), value, intent(in) :: accuracy

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vcalc)) then
         call env%ptr%error("Singlepoint calculator is not allocated", source)
         return
      end if
      call c_f_pointer(vcalc, calc)

      if (.not.allocated(calc%ptr)) then
         call env%ptr%error("Setting accuracy not possible, no calculator loaded", &
            & source)
         return
      end if

      if (accuracy < 1.e-4_c_double) then
         call env%ptr%warning("We cannot provide this level of accuracy, "//&
            & "resetted accuracy to 0.0001", source)
         calc%ptr%accuracy = 1.e-4_wp
      else if (accuracy > 1.e+3_c_double) then
         call env%ptr%warning("We cannot provide this level of accuracy, "//&
            & "resetted accuracy to 1000", source)
         calc%ptr%accuracy = 1.e+3_wp
      else
         calc%ptr%accuracy = accuracy
      end if

   end if

end subroutine setAccuracy_api


subroutine setMaxIter_api(venv, vcalc, maxiter) &
      & bind(C, name="xtb_setMaxIter")
   !DEC$ ATTRIBUTES DLLEXPORT :: setMaxIter_api
   character(len=*), parameter :: source = 'xtb_api_setMaxIter'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc
   integer(c_int), value, intent(in) :: maxiter

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vcalc)) then
         call env%ptr%error("Singlepoint calculator is not allocated", source)
         return
      end if
      call c_f_pointer(vcalc, calc)

      if (.not.allocated(calc%ptr)) then
         call env%ptr%error("Setting accuracy not possible, no calculator loaded", &
            & source)
         return
      end if

      select type(xtb => calc%ptr)
      class default
         call env%ptr%warning("Cannot set iterations for non-iterative method", &
            & source)
      type is(TxTBCalculator)
         xtb%maxiter = max(maxiter, 1)
      end select

   end if

end subroutine setMaxIter_api


subroutine setElectronicTemp_api(venv, vcalc, temperature) &
      & bind(C, name="xtb_setElectronicTemp")
   !DEC$ ATTRIBUTES DLLEXPORT :: setElectronicTemp_api
   character(len=*), parameter :: source = 'xtb_api_setElectronicTemp'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc
   real(c_double), value, intent(in) :: temperature

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vcalc)) then
         call env%ptr%error("Singlepoint calculator is not allocated", source)
         return
      end if
      call c_f_pointer(vcalc, calc)

      if (.not.allocated(calc%ptr)) then
         call env%ptr%error("Setting accuracy not possible, no calculator loaded", &
            & source)
         return
      end if

      select type(xtb => calc%ptr)
      class default
         call env%ptr%warning("Calculator does not support electronic temperature", &
            & source)
      type is(TxTBCalculator)
         xtb%etemp = max(temperature, 1.0e-6_c_double)
      end select

   end if

end subroutine setElectronicTemp_api


end module xtb_api_calculator
