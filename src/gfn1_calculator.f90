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
submodule(tb_calculators) gfn1_calc_implementation
   implicit none
contains
! ========================================================================
!> GFN1-xTB calculation
module subroutine gfn1_calculation &
      (iunit,env,err,opt,mol,pcem,wfn,hl_gap,energy,gradient)
   use xtb_mctc_accuracy, only : wp

   use mctc_systools
   use mctc_logging

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_wavefunction
   use xtb_type_basisset
   use xtb_type_param
   use xtb_type_data
   use xtb_type_pcem

   use setparam, only : gfn_method, ngrida
   use aoparam,  only : use_parameterset

   use xbasis
   use eeq_model
   use ncoord
   use scc_core
   use scf_module
   use gbobc
   use embedding

   implicit none

   integer, intent(in) :: iunit

   type(TMolecule),    intent(inout) :: mol
   type(scc_options),    intent(in)    :: opt
   type(tb_environment), intent(in)    :: env
   type(mctc_error), allocatable, intent(inout) :: err
   type(tb_pcem),        intent(inout) :: pcem
   type(TWavefunction),intent(inout) :: wfn

   real(wp), intent(out) :: energy
   real(wp), intent(out) :: hl_gap
   real(wp), intent(out) :: gradient(3,mol%n)

   integer, parameter    :: wsc_rep(3) = [1,1,1] ! FIXME

   type(TBasisset)     :: basis
   type(scc_parameter)   :: param
   type(scc_results)     :: res
   type(chrg_parameter)  :: chrgeq

   real(wp), allocatable :: cn(:)

   character(len=*),parameter :: outfmt = &
      '(9x,"::",1x,a,f24.12,1x,a,1x,"::")'
   character(len=*), parameter   :: p_fnv_gfn1 = '.param_gfn.xtb'
   character(len=:), allocatable :: fnv
   real(wp) :: globpar(25)
   integer  :: ipar
   logical  :: exist

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

   ! we will try an internal parameter file first to avoid IO
   call use_parameterset(p_fnv_gfn1,globpar,exist)
   ! no luck, we have to fire up some IO to get our parameters
   if (.not.exist) then
      ! let's check if we can find the parameter file
      call rdpath(env%xtbpath,p_fnv_gfn1,fnv,exist)
      ! maybe the user provides a local parameter file, this was always
      ! an option in `xtb', so we will give it a try
      if (.not.exist) fnv = p_fnv_gfn1
      call open_file(ipar,fnv,'r')
      if (ipar.eq.-1) then
         ! at this point there is no chance to recover from this error
         err = mctc_error("Parameter file '"//fnv//"' not found")
         return
      endif
      call read_gfn_param(ipar,globpar,.true.)
      call close_file(ipar)
   endif
   call set_gfn1_parameter(param,globpar)
   if (opt%prlevel > 1) then
      call gfn1_header(iunit)
      call gfn1_prparam(iunit,mol%n,mol%at,param)
   endif

   lgbsa = len_trim(opt%solvent).gt.0 .and. opt%solvent.ne."none"
   if (lgbsa) then
      call init_gbsa(iunit,trim(opt%solvent),0,opt%etemp,gfn_method,ngrida)
   endif

   ! ====================================================================
   !  STEP 3: expand our Slater basis set in contracted Gaussians
   ! ====================================================================

   call xbasis0(mol%n,mol%at,basis)
   call xbasis_gfn1(mol%n,mol%at,basis,okbas,diff)

   ! ====================================================================
   !  STEP 4: setup the initial wavefunction
   ! ====================================================================

   call wfn%allocate(mol%n,basis%nshell,basis%nao)

   ! do a EEQ guess
   allocate( cn(mol%n), source = 0.0_wp )
   call new_charge_model_2019(chrgeq,mol%n,mol%at)
   call ncoord_erf(mol%n,mol%at,mol%xyz,cn)
   call eeq_chrgeq(mol,err,chrgeq,cn,wfn%q)
   deallocate(cn)
   if (allocated(err)) return

   call iniqshell(mol%n,mol%at,mol%z,basis%nshell,wfn%q,wfn%qsh,gfn_method)

   ! ====================================================================
   !  STEP 5: do the calculation
   ! ====================================================================
   call scf(iunit,err,mol,wfn,basis,param,pcem,hl_gap, &
      &     opt%etemp,opt%maxiter,opt%prlevel,.false.,opt%grad,opt%acc, &
      &     energy,gradient,res)
   if (allocated(err)) return

   if (opt%prlevel > 0) then
      write(iunit,'(9x,53(":"))')
      write(iunit,outfmt) "total energy      ", res%e_total,"Eh  "
      write(iunit,outfmt) "gradient norm     ", res%gnorm,  "Eh/α"
      write(iunit,outfmt) "HOMO-LUMO gap     ", res%hl_gap, "eV  "
      write(iunit,'(9x,53(":"))')
   endif

end subroutine gfn1_calculation
end submodule gfn1_calc_implementation
