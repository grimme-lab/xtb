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
submodule(xtb_calculators) gfn1_calc_implementation
   implicit none
contains
! ========================================================================
!> GFN1-xTB calculation
module subroutine gfn1_calculation &
      (iunit,env,opt,mol,pcem,wfn,hl_gap,energy,gradient)
   use xtb_mctc_accuracy, only : wp

   use xtb_mctc_systools

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_wavefunction
   use xtb_type_basisset
   use xtb_type_param
   use xtb_type_data
   use xtb_type_pcem

   use xtb_setparam, only : gfn_method, ngrida

   use xtb_basis
   use xtb_eeq
   use xtb_chargemodel
   use xtb_disp_ncoord
   use xtb_scc_core
   use xtb_scf
   use xtb_solv_gbobc
   use xtb_embedding
   use xtb_restart
   use xtb_readparam
   use xtb_paramset

   use xtb_main_setup
   use xtb_xtb_calculator
   use xtb_xtb_data
   use xtb_xtb_gfn1

   implicit none

   character(len=*), parameter :: source = 'calculator_gfn1'

   integer, intent(in) :: iunit

   type(TMolecule),    intent(inout) :: mol
   type(scc_options),    intent(in)    :: opt
   type(TEnvironment), intent(inout)    :: env
   type(tb_pcem),        intent(inout) :: pcem
   type(TWavefunction),intent(inout) :: wfn

   real(wp), intent(out) :: energy
   real(wp), intent(out) :: hl_gap
   real(wp), intent(out) :: gradient(3,mol%n)

   integer, parameter    :: wsc_rep(3) = [1,1,1] ! FIXME

   type(TxTBCalculator) :: calc
   type(scc_results)     :: res
   type(chrg_parameter)  :: chrgeq

   real(wp), allocatable :: cn(:)

   character(len=*),parameter :: outfmt = &
      '(9x,"::",1x,a,f24.12,1x,a,1x,"::")'
   character(len=*), parameter   :: p_fnv_gfn1 = 'param_gfn1-xtb.txt'
   character(len=:), allocatable :: fnv
   type(TxTBParameter) :: globpar
   integer  :: ipar
   logical  :: exist
   logical :: exitRun

   logical  :: okbas,diff

   gfn_method = 1
   call init_pcem

   ! ====================================================================
   !  STEP 1: prepare geometry input
   ! ====================================================================
   ! we assume that the user provides a resonable molecule input
   ! -> all atoms are inside the unit cell, all data is set and consistent

   wfn%nel = nint(sum(mol%z) - mol%chrg)
   wfn%nopen = mol%uhf
   ! at this point, don't complain about odd multiplicities for even electron
   ! systems and just fix it silently, the API is supposed catch this
   if (mod(wfn%nopen,2) == 0.and.mod(wfn%nel,2) /= 0) wfn%nopen = 1
   if (mod(wfn%nopen,2) /= 0.and.mod(wfn%nel,2) == 0) wfn%nopen = 0

   ! give an optional summary on the geometry used
   if (opt%prlevel > 2) then
      call main_geometry(iunit,mol)
   endif

   ! ====================================================================
   !  STEP 2: get the parametrisation
   ! ====================================================================
   ! we could require our user to perform this step, but if we want
   ! to be sure about getting the correct parameters, we should do it here
   call rdpath(env%xtbpath, p_fnv_gfn1, fnv, exist)
   if (exist) then
      call newXTBCalculator(env, mol, calc, fnv)
   else
      call newXTBCalculator(env, mol, calc, p_fnv_gfn1)
   end if

   if (opt%prlevel > 1) then
      call gfn1_header(iunit)
   endif

   lgbsa = len_trim(opt%solvent).gt.0 .and. opt%solvent.ne."none"
   if (lgbsa) then
      call init_gbsa(iunit,trim(opt%solvent),0,opt%etemp,gfn_method,ngrida, &
         & opt%prlevel > 0)
   endif

   ! ====================================================================
   !  STEP 4: setup the initial wavefunction
   ! ====================================================================

   call wfn%allocate(mol%n,calc%basis%nshell,calc%basis%nao)

   if (mol%npbc == 0) then
      ! do a EEQ guess
      allocate( cn(mol%n), source = 0.0_wp )
      call new_charge_model_2019(chrgeq,mol%n,mol%at)
      call ncoord_erf(mol%n,mol%at,mol%xyz,cn)
      call eeq_chrgeq(mol,env,chrgeq,cn,wfn%q)
      deallocate(cn)
      call env%check(exitRun)
      if (exitRun) then
         call env%error("EEQ quess failed", source)
      end if
   else
      wfn%q(:) = 0.0_wp
   end if

   call iniqshell(calc%xtbData,mol%n,mol%at,mol%z,calc%basis%nshell,wfn%q,wfn%qsh,gfn_method)

   if (opt%restart) &
      call readRestart(env,wfn,'xtbrestart',mol%n,mol%at,gfn_method,exist,.false.)

   ! ====================================================================
   !  STEP 5: do the calculation
   ! ====================================================================
   call scf(env,mol,wfn,calc%basis,pcem,calc%xtbData,hl_gap, &
      &     opt%etemp,opt%maxiter,opt%prlevel,.false.,opt%grad,opt%acc, &
      &     energy,gradient,res)

   call env%check(exitRun)
   if (exitRun) then
      call env%error("SCF calculation terminated", source)
   end if

   if (opt%restart) then
      call writeRestart(env,wfn,'xtbrestart',gfn_method)
   endif

   if (opt%prlevel > 0) then
      write(iunit,'(9x,53(":"))')
      write(iunit,outfmt) "total energy      ", res%e_total,"Eh  "
      write(iunit,outfmt) "gradient norm     ", res%gnorm,  "Eh/α"
      write(iunit,outfmt) "HOMO-LUMO gap     ", res%hl_gap, "eV  "
      write(iunit,'(9x,53(":"))')
   endif

end subroutine gfn1_calculation
end submodule gfn1_calc_implementation
