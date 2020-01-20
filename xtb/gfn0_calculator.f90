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
submodule(tb_calculators) gfn0_calc_implementation
   implicit none
contains
! ========================================================================
!> periodic GFN0-xTB (PEEQ) calculation
module subroutine gfn0_calculation &
      (iunit,env,opt,mol,hl_gap,energy,gradient,stress,lattice_gradient)
   use iso_fortran_env, wp => real64

   use mctc_systools

   use tbdef_options
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data

   use setparam, only : gfn_method, ngrida
   use aoparam,  only : use_parameterset

   use pbc_tools
   use xbasis
   use peeq_module
   use gbobc

   implicit none

   integer, intent(in) :: iunit

   type(tb_molecule),    intent(inout) :: mol
   type(peeq_options),   intent(in)    :: opt
   type(tb_environment), intent(in)    :: env

   real(wp), intent(out) :: energy
   real(wp), intent(out) :: hl_gap
   real(wp), intent(out) :: gradient(3,mol%n)
   real(wp), intent(out) :: stress(3,3)
   real(wp), intent(out) :: lattice_gradient(3,3)
   real(wp)              :: sigma(3,3)
   real(wp)              :: inv_lat(3,3)

   integer, parameter    :: wsc_rep(3) = [1,1,1] ! FIXME

   type(tb_wavefunction) :: wfn
   type(tb_basisset)     :: basis
   type(scc_parameter)   :: param
   type(scc_results)     :: res

   character(len=*),parameter :: outfmt = &
      '(9x,"::",1x,a,f24.12,1x,a,1x,"::")'
   character(len=*), parameter   :: p_fnv_gfn0 = '.param_gfn0.xtb'
   character(len=:), allocatable :: fnv
   real(wp) :: globpar(25)
   integer  :: ipar,i
   logical  :: exist

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

   ! we will try an internal parameter file first to avoid IO
   call use_parameterset(p_fnv_gfn0,globpar,exist)
   ! no luck, we have to fire up some IO to get our parameters
   if (.not.exist) then
      ! let's check if we can find the parameter file
      call rdpath(env%xtbpath,p_fnv_gfn0,fnv,exist)
      ! maybe the user provides a local parameter file, this was always
      ! an option in `xtb', so we will give it a try
      if (.not.exist) fnv = p_fnv_gfn0
      call open_file(ipar,fnv,'r')
      if (ipar.eq.-1) then
         ! at this point there is no chance to recover from this error
         ! THEREFORE, we have to kill the program
         call raise('E',"Parameter file '"//fnv//"' not found!",1)
         return
      endif
      call read_gfn_param(ipar,globpar,.true.)
      call close_file(ipar)
   endif
   call set_gfn0_parameter(param,globpar)
   if (opt%prlevel > 1) then
      call gfn0_header(iunit)
      call gfn0_prparam(iunit,mol%n,mol%at,param)
   endif

   lgbsa = len_trim(opt%solvent).gt.0 .and. opt%solvent.ne."none" &
      &    .and. mol%npbc == 0 ! GBSA is not yet periodic
   if (lgbsa) then
      call init_gbsa(iunit,trim(opt%solvent),0,opt%etemp,gfn_method,ngrida)
   endif

   ! ====================================================================
   !  STEP 3: expand our Slater basis set in contracted Gaussians
   ! ====================================================================

   call xbasis0(mol%n,mol%at,basis)
   call xbasis_gfn0(mol%n,mol%at,basis,okbas,diff)

   ! ====================================================================
   !  STEP 4: setup the initial wavefunction
   ! ====================================================================

   call wfn%allocate(mol%n,basis%nshell,basis%nao)
   ! do a SAD guess since we are not need any of this information later
   wfn%q = mol%chrg / real(mol%n,wp)

   ! ====================================================================
   !  STEP 5: do the calculation
   ! ====================================================================

   call peeq(iunit,mol,wfn,basis,param,hl_gap,opt%etemp,opt%prlevel,opt%grad, &
      &      opt%ccm,opt%acc,energy,gradient,sigma,res)

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
