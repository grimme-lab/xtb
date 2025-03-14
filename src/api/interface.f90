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

subroutine cpcmx_calc_api(venv, vmol, vcalc, vres) &
      & bind(C, name="xtb_cpcmx_calc")
   !DEC$ ATTRIBUTES DLLEXPORT :: cpcmx_calc_api
   use xtb_solv_cpx, only: TCpcmx
   use xtb_type_calculator, only: TCalculator

   character(len=*), parameter :: source = 'xtb_api_cpcmx_calc'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vmol
   type(VMolecule), pointer :: mol
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc
   type(c_ptr), value :: vres
   type(VResults), pointer :: res
   type(scc_results) :: spRes

   type(TCpcmx) :: cpx
   type(VCalculator) :: cpx_calc
   real(c_double) :: energy_gas

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

      ! Fail early if CPCM-X solvation model is not set
      if (allocated(calc%ptr%solvation)) then
         if (.not. allocated(calc%ptr%solvation%cpxsolvent)) then
            call env%ptr%error("CPCM-X solvent not set", source)
            return
         end if
      else
         call env%ptr%error("No solvation input given", source)
         return
      end if

      ! Fail early if not using xTB
      select type(xtb => calc%ptr)
        type is (TGFFCalculator)
           call env%ptr%error("CPCM-X is not possible with a force field.", source)
           return
      end select

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

      ! Initial COSMO singlepoint calculation
      call calc%ptr%singlepoint(env%ptr, mol%ptr, res%chk, env%verbosity, .true., &
         & res%energy, res%gradient, res%sigma, res%egap, spRes)

      ! CPCM-X calculation
      call cpx%setup(env%ptr, calc%ptr%solvation%cpxsolvent)
      cpx_calc = calc
      deallocate(cpx_calc%ptr%solvation)

      energy_gas = res%energy
      call generic_header(env%ptr%unit, "CPCM-X post-SCF solvation evaluation", 49, 10)
      call cpx_calc%ptr%singlepoint(env%ptr, mol%ptr, res%chk, env%verbosity, .false., &
         & energy_gas, res%gradient, res%sigma, res%egap, spRes)


      call cpx%calc_solv(env%ptr, calc%ptr%solvation%cpxsolvent, energy_gas, &
             &           0.4_wp, 298.15_wp, 500, 0.0001_wp, spRes%e_total)
      call cpx%print(.true.)

      res%energy = spRes%e_total
      res%solvation_energy = res%energy - energy_gas
      res%dipole = spRes%dipole

      ! Zero out the gradient and sigma (not yet implemented for CPCM-X)
      res%gradient = 0.0_wp
      res%sigma = 0.0_wp

   endif

end subroutine cpcmx_calc_api


subroutine hessian_api(venv, vmol, vcalc, vres, c_hess, &
                     & c_step, c_list, c_dipgrad, c_polgrad) &
      & bind(C, name="xtb_hessian")
   !DEC$ ATTRIBUTES DLLEXPORT :: hessian_api
   character(len=*), parameter :: source = 'xtb_api_hessian'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vmol
   type(VMolecule), pointer :: mol
   type(c_ptr), value :: vcalc
   type(VCalculator), pointer :: calc
   type(c_ptr), value :: vres
   type(VResults), pointer :: res

   !> Array to add Hessian to
   real(c_double), intent(inout) :: c_hess(*)
   real(wp), allocatable :: hess(:, :)
   !> List of atoms to displace
   integer(c_int), intent(in), optional :: c_list(:)
   integer, allocatable :: list(:)
   !> Step size for numerical differentiation
   real(c_double), intent(in), optional :: c_step
   real(wp) :: step
   !> Array to add dipole gradient to
   real(c_double), intent(inout), optional :: c_dipgrad(*)
   real(wp), allocatable :: dipgrad(:, :)
   !> Array to add polarizability gradient to
   real(c_double), intent(inout), optional :: c_polgrad(*)
   real(wp), allocatable :: polgrad(:, :)

   integer :: natom, natsq, i, j
   logical :: has_polgrad, has_dipgrad

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vmol)) then
         call env%ptr%error("Molecular structure data is not allocated", source)
         return
      end if
      call c_f_pointer(vmol, mol)
      natom = mol%ptr%n
      natsq = natom * natom

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

      hess = reshape(c_hess(:9*natsq), &
                    &(/3*natom, 3*natom/))
      ! Need to initialize, as the subroutine increments the values
      hess = 0.0_wp

      if (.not.present(c_step)) then
         step = 0.005_wp
      else
         step = c_step
      end if

      if (.not.present(c_list)) then
         list = [(i, i=1, natom)]
      else
         list = c_list
      end if

      ! Dipole gradient is required by the hessian method,
      ! so we have to allocate it
      has_dipgrad = present(c_dipgrad)
      if (.not. has_dipgrad) then
         allocate(dipgrad(3, 3*natom))
      else
         dipgrad = reshape(c_dipgrad(:9*natom), &
                     &(/3, 3*natom/))
      end if

      has_polgrad = present(c_polgrad)
      if (has_polgrad) then
         polgrad = reshape(c_polgrad(:18*natom), &
                         &(/6, 3*natom/))
      end if

      ! hessian calculation
      if (has_polgrad) then
         call calc%ptr%hessian(env%ptr, mol%ptr, res%chk, list, step, &
                            &  hess, dipgrad, polgrad)
      else
         call calc%ptr%hessian(env%ptr, mol%ptr, res%chk, list, step, &
                           &   hess, dipgrad)
      end if

      ! Symmetrize the hessian
      do i = 1, 3*natom
         do j = i+1, 3*natom
            hess(i, j) = 0.5_wp * (hess(i, j) + hess(j, i))
            hess(j, i) = hess(i, j)
         end do
      end do

      ! copy back the results
      c_hess(:9*natsq) = reshape(hess, (/9*natsq/))
      if (has_dipgrad) then
         c_dipgrad(:9*natom) = reshape(dipgrad, (/9*natom/))
      end if
      if (has_polgrad) then
         c_polgrad(:18*natom) = reshape(polgrad, (/18*natom/))
      end if
   end if

end subroutine hessian_api

end module xtb_api_interface
