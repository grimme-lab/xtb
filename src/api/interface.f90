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

!> Actual calculation interface
module xtb_api_interface
   use, intrinsic :: iso_c_binding
   use xtb_mctc_accuracy, only : wp
   use xtb_api_calculator
   use xtb_api_environment
   use xtb_api_molecule
   use xtb_api_results
   use xtb_api_utils
   use xtb_gfnff_calculator, only : TGFFCalculator
   use xtb_scc_core, only : iniqshell
   use xtb_type_data, only : scc_results
   use xtb_xtb_calculator, only : TxTBCalculator, newWavefunction
   implicit none
   private

   public :: singlepoint_api


contains


subroutine singlepoint_api(venv, vmol, vcalc, vres) &
      & bind(C, name="xtb_singlepoint")
   !DEC$ ATTRIBUTES DLLEXPORT :: singlepoint_api
   character(len=*), parameter :: source = 'xtb_api_singlepoint'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vmol
   type(VMolecule), pointer :: mol
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc
   type(c_ptr), value :: vres
   type(VResults), pointer :: res
   type(scc_results) :: spRes

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

      if (.not.allocated(calc%ptr)) then
         call env%ptr%error("No calculator loaded for single point", &
            & source)
         return
      end if

      if (.not.c_associated(vres)) then
         call env%ptr%error("Calculation results are not allocated", source)
         return
      end if
      call c_f_pointer(vres, res)

      ! check cache, automatically invalidate missmatched data
      if (allocated(res%chk)) then
         select type(xtb => calc%ptr)
         type is(TxTBCalculator)
            if (res%chk%wfn%n /= mol%ptr%n .or. res%chk%wfn%n /= xtb%basis%n .or. &
               & res%chk%wfn%nao /= xtb%basis%nao .or. &
               & res%chk%wfn%nshell /= xtb%basis%nshell) then
               deallocate(res%chk)
            end if
         end select
      end if

      if (.not.allocated(res%chk)) then
         allocate(res%chk)
         ! in case of a new wavefunction cache we have to perform an initial guess
         select type(xtb => calc%ptr)
         type is(TxTBCalculator)
            call newWavefunction(env%ptr, mol%ptr, xtb, res%chk)
         end select
      end if

      if (.not.allocated(res%energy)) then
         allocate(res%energy)
      end if

      if (.not.allocated(res%egap)) then
         allocate(res%egap)
      end if

      if (allocated(res%pcgradient)) then
         deallocate(res%pcgradient)
      end if

      if (allocated(res%gradient)) then
         if (any(shape(res%gradient) /= [3, mol%ptr%n])) then
            call env%ptr%warning("Shape missmatch in gradient, reallocating", source)
            deallocate(res%gradient)
         end if
      end if
      if (.not.allocated(res%gradient)) then
         allocate(res%gradient(3, mol%ptr%n))
      end if

      if (allocated(res%sigma)) then
         if (any(shape(res%sigma) /= [3, 3])) then
            call env%ptr%warning("Shape missmatch in virial, reallocating", source)
            deallocate(res%sigma)
         end if
      end if
      if (.not.allocated(res%sigma)) then
         allocate(res%sigma(3, 3))
      end if

      ! singlepoint calculation
      call calc%ptr%singlepoint(env%ptr, mol%ptr, res%chk, env%verbosity, .true., &
         & res%energy, res%gradient, res%sigma, res%egap, spRes)

      ! invalidate cache for properties not produced in GFN-FF
      select type(gfnff => calc%ptr)
      type is(TGFFCalculator)
         deallocate(res%chk)
         deallocate(res%egap)
         deallocate(res%sigma)
      end select

      ! check if external charge gradients have been calculated
      if (allocated(spRes%pcem%grd) .and. spRes%pcem%n > 0) then
         res%pcgradient = spRes%pcem%grd
      end if

      res%dipole = spRes%dipole

   end if

end subroutine singlepoint_api


end module xtb_api_interface
