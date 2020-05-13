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
submodule(xtb_calculators) gfn0_calc_implementation
   implicit none
contains
! ========================================================================
!> periodic GFN0-xTB (PEEQ) calculation
module subroutine gfn0_calculation &
      (iunit,env,opt,mol,hl_gap,energy,gradient,stress,lattice_gradient)
   use xtb_mctc_accuracy, only : wp

   use xtb_mctc_systools

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_wavefunction
   use xtb_type_basisset
   use xtb_type_param
   use xtb_type_data

   use xtb_setparam, only : gfn_method, ngrida

   use xtb_pbc_tools
   use xtb_basis
   use xtb_peeq
   use xtb_solv_gbobc
   use xtb_readparam
   use xtb_paramset

   use xtb_main_setup
   use xtb_xtb_calculator
   use xtb_xtb_data
   use xtb_xtb_gfn0

   implicit none

   character(len=*), parameter :: source = 'calculator_gfn0'

   integer, intent(in) :: iunit

   type(TMolecule),    intent(inout) :: mol
   type(peeq_options),   intent(in)    :: opt
   type(TEnvironment), intent(inout)    :: env

   real(wp), intent(out) :: energy
   real(wp), intent(out) :: hl_gap
   real(wp), intent(out) :: gradient(3,mol%n)
   real(wp), intent(out) :: stress(3,3)
   real(wp), intent(out) :: lattice_gradient(3,3)
   real(wp)              :: sigma(3,3)
   real(wp)              :: inv_lat(3,3)

   integer, parameter    :: wsc_rep(3) = [1,1,1] ! FIXME

   type(TWavefunction) :: wfn
   type(TxTBCalculator) :: calc
   type(scc_results)     :: res

   character(len=*),parameter :: outfmt = &
      '(9x,"::",1x,a,f24.12,1x,a,1x,"::")'
   character(len=*), parameter   :: p_fnv_gfn0 = 'param_gfn0-xtb.txt'
   character(len=:), allocatable :: fnv
   type(TxTBParameter) :: globpar
   integer  :: ipar,i
   logical  :: exist
   logical :: exitRun

   logical  :: okbas,diff

   gfn_method = 0
   sigma = 0.0_wp

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

   if (mol%npbc > 0) then
      ! get Wigner-Seitz cell if necessary
      ! -> this modifies the molecule, since the WSC is bound to a molecule
      call generate_wsc(mol,mol%wsc,wsc_rep)
   endif

   ! give an optional summary on the geometry used
   if (opt%prlevel > 2) then
      call main_geometry(iunit,mol)
   endif

   ! ====================================================================
   !  STEP 2: get the parametrisation
   ! ====================================================================
   ! we could require our user to perform this step, but if we want
   ! to be sure about getting the correct parameters, we should do it here
   ! let's check if we can find the parameter file
   call rdpath(env%xtbpath, p_fnv_gfn0, fnv, exist)
   if (exist) then
      call newXTBCalculator(env, mol, calc, fnv)
   else
      call newXTBCalculator(env, mol, calc, p_fnv_gfn0)
   end if

   if (opt%prlevel > 1) then
      call gfn0_header(iunit)
   endif

   lgbsa = len_trim(opt%solvent).gt.0 .and. opt%solvent.ne."none" &
      &    .and. mol%npbc == 0 ! GBSA is not yet periodic
   if (lgbsa) then
      call init_gbsa(iunit,trim(opt%solvent),0,opt%etemp,gfn_method,ngrida)
   endif

   ! ====================================================================
   !  STEP 4: setup the initial wavefunction
   ! ====================================================================

   call wfn%allocate(mol%n,calc%basis%nshell,calc%basis%nao)
   ! do a SAD guess since we are not need any of this information later
   wfn%q = mol%chrg / real(mol%n,wp)

   ! ====================================================================
   !  STEP 5: do the calculation
   ! ====================================================================

   call peeq(env,mol,wfn,calc%basis,calc%xtbData,hl_gap,opt%etemp,opt%prlevel, &
      & opt%grad,opt%ccm,opt%acc,energy,gradient,sigma,res)

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Single point calculation terminated", source)
      return
   end if

   if (mol%npbc > 0) then
      inv_lat = mat_inv_3x3(mol%lattice)
      call sigma_to_latgrad(sigma,inv_lat,lattice_gradient)
      stress = sigma/mol%volume
   endif

   if (opt%prlevel > 0) then
      write(iunit,'(9x,53(":"))')
      write(iunit,outfmt) "total energy      ", res%e_total,"Eh  "
      write(iunit,outfmt) "gradient norm     ", res%gnorm,  "Eh/α"
      if (mol%npbc > 0) &
      write(iunit,outfmt) "gradlatt norm     ", norm2(lattice_gradient),  "Eh/α"
      write(iunit,outfmt) "HOMO-LUMO gap     ", res%hl_gap, "eV  "
      write(iunit,'(9x,53(":"))')
   endif

end subroutine gfn0_calculation
end submodule gfn0_calc_implementation
