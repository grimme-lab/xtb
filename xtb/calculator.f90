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
module tb_calculators
   implicit none
contains

! ========================================================================
!> periodic GFN0-xTB (PEEQ) calculation
subroutine gfn0_calculation &
      (iunit,env,opt,mol,gfn,hl_gap,energy,gradient,stress,lattice_gradient)
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
   type(gfn_parameter),  intent(in)    :: gfn
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
      call print_pbcsum(iunit,mol)
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
   call set_gfn0_parameter(param,globpar,mol%n,mol%at)
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
   call xbasis_cao2sao(mol%n,mol%at,basis)

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

! ========================================================================
!> GFN2-xTB calculation
subroutine gfn2_calculation &
      (iunit,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)
   use iso_fortran_env, wp => real64

   use mctc_systools

   use tbdef_options
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data
   use tbdef_pcem

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

   type(tb_molecule),    intent(inout) :: mol
   type(tb_wavefunction),intent(inout) :: wfn
   type(gfn_parameter),  intent(in)    :: gfn
   type(scc_options),    intent(in)    :: opt
   type(tb_environment), intent(in)    :: env
   type(tb_pcem),        intent(inout) :: pcem

   real(wp), intent(out) :: energy
   real(wp), intent(out) :: hl_gap
   real(wp), intent(out) :: gradient(3,mol%n)

   integer, parameter    :: wsc_rep(3) = [1,1,1] ! FIXME

   type(tb_basisset)     :: basis
   type(scc_parameter)   :: param
   type(scc_results)     :: res
   type(chrg_parameter)  :: chrgeq

   real(wp), allocatable :: cn(:)

   character(len=*),parameter :: outfmt = &
      '(9x,"::",1x,a,f24.12,1x,a,1x,"::")'
   character(len=*), parameter   :: p_fnv_gfn2 = '.param_gfn2.xtb'
   character(len=:), allocatable :: fnv
   real(wp) :: globpar(25)
   integer  :: ipar
   logical  :: exist

   logical  :: okbas

   gfn_method = 2
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
   call use_parameterset(p_fnv_gfn2,globpar,exist)
   ! no luck, we have to fire up some IO to get our parameters
   if (.not.exist) then
      ! let's check if we can find the parameter file
      call rdpath(env%xtbpath,p_fnv_gfn2,fnv,exist)
      ! maybe the user provides a local parameter file, this was always
      ! an option in `xtb', so we will give it a try
      if (.not.exist) fnv = p_fnv_gfn2
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
   call set_gfn2_parameter(param,globpar,mol%n,mol%at)
   if (opt%prlevel > 1) then
      call gfn2_header(iunit)
      call gfn2_prparam(iunit,mol%n,mol%at,param)
   endif

   lgbsa = len_trim(opt%solvent).gt.0 .and. opt%solvent.ne."none"
   if (lgbsa) then
      call init_gbsa(iunit,trim(opt%solvent),0,opt%etemp,gfn_method,ngrida)
   endif

   ! ====================================================================
   !  STEP 3: expand our Slater basis set in contracted Gaussians
   ! ====================================================================

   call xbasis0(mol%n,mol%at,basis)
   call xbasis_gfn2(mol%n,mol%at,basis,okbas)
   call xbasis_cao2sao(mol%n,mol%at,basis)

   ! ====================================================================
   !  STEP 4: setup the initial wavefunction
   ! ====================================================================

   call wfn%allocate(mol%n,basis%nshell,basis%nao)

   ! do an EEQ guess
   allocate( cn(mol%n), source = 0.0_wp )
   call new_charge_model_2019(chrgeq,mol%n,mol%at)
   call ncoord_erf(mol%n,mol%at,mol%xyz,cn)
   call eeq_chrgeq(mol,chrgeq,cn,wfn%q)
   deallocate(cn)

   call iniqshell(mol%n,mol%at,mol%z,basis%nshell,wfn%q,wfn%qsh,gfn_method)

   if (opt%restart) &
      call read_restart(wfn,'xtbrestart',mol%n,mol%at,gfn_method,exist,.false.)

   ! ====================================================================
   !  STEP 5: do the calculation
   ! ====================================================================
   call scf(iunit,mol,wfn,basis,param,pcem,hl_gap, &
      &     opt%etemp,opt%maxiter,opt%prlevel,.false.,opt%grad,opt%acc, &
      &     energy,gradient,res)

   if (opt%restart) then
      call write_restart(wfn,'xtbrestart',gfn_method)
   endif 

   if (opt%prlevel > 0) then
      write(iunit,'(9x,53(":"))')
      write(iunit,outfmt) "total energy      ", res%e_total,"Eh  "
      write(iunit,outfmt) "gradient norm     ", res%gnorm,  "Eh/α"
      write(iunit,outfmt) "HOMO-LUMO gap     ", res%hl_gap, "eV  "
      write(iunit,'(9x,53(":"))')
   endif

end subroutine gfn2_calculation

! ========================================================================
!> GFN1-xTB calculation
subroutine gfn1_calculation &
      (iunit,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)
   use iso_fortran_env, wp => real64

   use mctc_systools

   use tbdef_options
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data
   use tbdef_pcem

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

   type(tb_molecule),    intent(inout) :: mol
   type(gfn_parameter),  intent(in)    :: gfn
   type(scc_options),    intent(in)    :: opt
   type(tb_environment), intent(in)    :: env
   type(tb_pcem),        intent(inout) :: pcem
   type(tb_wavefunction),intent(inout) :: wfn

   real(wp), intent(out) :: energy
   real(wp), intent(out) :: hl_gap
   real(wp), intent(out) :: gradient(3,mol%n)

   integer, parameter    :: wsc_rep(3) = [1,1,1] ! FIXME

   type(tb_basisset)     :: basis
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
         ! THEREFORE, we have to kill the program
         call raise('E',"Parameter file '"//fnv//"' not found!",1)
         return
      endif
      call read_gfn_param(ipar,globpar,.true.)
      call close_file(ipar)
   endif
   call set_gfn1_parameter(param,globpar,mol%n,mol%at)
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
   call xbasis_cao2sao(mol%n,mol%at,basis)

   ! ====================================================================
   !  STEP 4: setup the initial wavefunction
   ! ====================================================================

   call wfn%allocate(mol%n,basis%nshell,basis%nao)

   ! do a EEQ guess
   allocate( cn(mol%n), source = 0.0_wp )
   call new_charge_model_2019(chrgeq,mol%n,mol%at)
   call ncoord_erf(mol%n,mol%at,mol%xyz,cn)
   call eeq_chrgeq(mol,chrgeq,cn,wfn%q)
   deallocate(cn)

   call iniqshell(mol%n,mol%at,mol%z,basis%nshell,wfn%q,wfn%qsh,gfn_method)

   ! ====================================================================
   !  STEP 5: do the calculation
   ! ====================================================================
   call scf(iunit,mol,wfn,basis,param,pcem,hl_gap, &
      &     opt%etemp,opt%maxiter,opt%prlevel,.false.,opt%grad,opt%acc, &
      &     energy,gradient,res)

   if (opt%prlevel > 0) then
      write(iunit,'(9x,53(":"))')
      write(iunit,outfmt) "total energy      ", res%e_total,"Eh  "
      write(iunit,outfmt) "gradient norm     ", res%gnorm,  "Eh/α"
      write(iunit,outfmt) "HOMO-LUMO gap     ", res%hl_gap, "eV  "
      write(iunit,'(9x,53(":"))')
   endif

end subroutine gfn1_calculation

!> interface to the DFT-D4 module
subroutine d4_calculation(iunit,opt,mol,dparam,energy,gradient)
   use iso_fortran_env, wp => real64
!$ use omp_lib

! ------------------------------------------------------------------------
!  class definitions
! ------------------------------------------------------------------------
   use tbdef_options
   use tbdef_param
   use tbdef_molecule

! ------------------------------------------------------------------------
!  interfaces
! ------------------------------------------------------------------------
   use ncoord
   use eeq_model
   use dftd4
   !use dfuncpar

   implicit none

   !> output unit, usually bound to STDOUT
   !  not used if opt%silent is set (maybe for fatal errors)
   integer, intent(in) :: iunit

! ------------------------------------------------------------------------
!  class declarations
! ------------------------------------------------------------------------
   !> molecule structure information
   type(tb_molecule), intent(inout) :: mol
   !> calculation options for DFT-D4 method
   type(dftd_options), intent(in) :: opt
   !> damping parameter for the DFT-D4 method
   type(dftd_parameter),intent(in) :: dparam
   type(chrg_parameter)            :: chrgeq

! ------------------------------------------------------------------------
!  output variables
! ------------------------------------------------------------------------
   !> dispersion energy, always referenced
   real(wp),intent(out) :: energy
   !> nuclear gradient, only references if opt%lgradient is set
   real(wp),intent(out) :: gradient(3,mol%n)

! ------------------------------------------------------------------------
!  local variables
! ------------------------------------------------------------------------
   integer  :: ndim                      ! matrix dimension
   integer  :: i,j,k,l,ii,jj
   integer  :: err
   real(wp) :: memory
   real(wp) :: etmp,etwo,emany,er,el,es
   real(wp) :: molpol,molc6,molc8        ! molecular Polarizibility
   real(wp) :: sigma(3,3)
   real(wp),allocatable :: q(:)          ! partial charges
   real(wp),allocatable :: dqdr(:,:,:)   ! partial charges
   real(wp),allocatable :: dqdL(:,:,:)   ! partial charges
   real(wp),allocatable :: covcn(:)      ! covalent coordination number
   real(wp),allocatable :: dcovcndr(:,:,:) ! covalent coordination number derivative
   real(wp),allocatable :: dcovcndL(:,:,:) ! covalent coordination number derivative
   real(wp),allocatable :: cn(:)         ! erf coordination number
   real(wp),allocatable :: dcndr(:,:,:)  ! erf coordination number derivative
   real(wp),allocatable :: dcndL(:,:,:)  ! erf coordination number derivative
   real(wp),allocatable :: gweights(:)   ! gaussian weights
   real(wp),allocatable :: refc6(:,:)    ! reference C6 coeffients
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: aw(:,:)
   real(wp),allocatable :: ges(:,:)
   real(wp),allocatable :: gr(:,:)
   real(wp),allocatable :: gl(:,:)
   real(wp),parameter   :: step = 1.0e-5_wp, step2 = 0.5_wp/step

! ------------------------------------------------------------------------
!  Output: Initialization and Parameter setup
! ------------------------------------------------------------------------
   call d4init(mol%n,mol%at,opt%g_a,opt%g_c,p_refq_goedecker,ndim)
   memory = ((mol%n)+(mol%n)+(ndim)+(ndim*ndim) &
            +(mol%n*mol%n)+(23*mol%n)+(mol%n)+(3*mol%n) &
            +(3*mol%n*mol%n)+(3*mol%n*(mol%n+1)) &
            +(3*mol%n*mol%n)) * wp / (1024.0_wp**2)

   call generic_header(iunit,'Calculation Setup',49,10)
   write(iunit,'(3x,a," : ",i6)')   'number of atoms     ', mol%n
   write(iunit,'(3x,a," : ",i6)')   'charge              ', nint(mol%chrg)
   write(iunit,'(3x,a," : ",a)')    'non-additivity corr.', lmbd2string(opt%lmbd)
   write(iunit,'(3x,a," : ",a)')    'charge model        ', refq2string(p_refq_goedecker)
!$ write(iunit,'(3x,a," : ",i6)')   'omp threads         ',omp_get_num_threads()
   write(iunit,'(3x,a," : ",f6.1,1x,a)') &
      "memory needed (est.)",memory,"Mb"
   if (opt%wf  /= 6.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f6.1)')   'weighting factor    ', opt%wf
   if (opt%g_a /= 3.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f6.1)')   'q-scale height      ', opt%g_a
   if (opt%g_c /= 2.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f6.1)')   'q-scale steepness   ', opt%g_c

   write(iunit,'(a)')

   !if (.not.opt%silent) &
   !call prd4ref(mol)

   allocate( q(mol%n),covcn(mol%n),gweights(ndim),refc6(ndim,ndim),&
             c6ab(mol%n,mol%n),aw(23,mol%n),cn(mol%n), &
             dcndr(3,mol%n,mol%n),dqdr(3,mol%n,mol%n+1), &
             dcovcndr(3,mol%n,mol%n), stat = err )
   if (err /= 0) then
      call raise('E','Memory allocation failed',1)
      return
   endif

   call get_d4_cn(mol,covcn,dcovcndr)
   call d4(mol%n,ndim,mol%at,opt%wf,opt%g_a,opt%g_c,covcn,gweights,refc6)

   call get_erf_cn(mol,cn,dcndr)

! ------------------------------------------------------------------------
!  get partial charges
! ------------------------------------------------------------------------
   !if (opt%verbose) &
   !call eeq_header
   call new_charge_model_2019(chrgeq,mol%n,mol%at)
   call eeq_chrgeq(mol,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
                   .false.,.false.,.true.)
   !if (opt%verbose) &
   !call print_chrgeq(iunit,chrgeq,mol,q,cn)

! ------------------------------------------------------------------------
!  calculate properties
! ------------------------------------------------------------------------
dispersion_properties: if (opt%lmolpol) then
   call generic_header(iunit,'Molecular Properties',49,10)
   call mdisp(mol%n,ndim,mol%at,q,mol%xyz,opt%g_a,opt%g_c,gweights,refc6, &
      &       molc6,molc8,molpol,aw,c6ab)
   call prmolc6(molc6,molc8,molpol,mol%n,mol%at, &
      &         covcn=covcn,q=q,c6ab=c6ab,alpha=aw(1,:))
endif dispersion_properties

if (.not.opt%silent.and.(opt%lenergy.or.opt%lgradient)) then
   call generic_header(iunit,'Damping Parameters',49,10)
   write(iunit,'(3x,a," : ",f10.4)')   's6                  ', dparam%s6
   write(iunit,'(3x,a," : ",f10.4)')   's8                  ', dparam%s8
   if (dparam%s10 /= 0.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f10.4)')   's10                 ', dparam%s10
   if (dparam%s9  /= 1.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f10.4)')   's9                  ', dparam%s9
   write(iunit,'(3x,a," : ",f10.4)')   'a1                  ', dparam%a1
   write(iunit,'(3x,a," : ",f10.4)')   'a2                  ', dparam%a2
   write(iunit,'(a)')
endif

! ------------------------------------------------------------------------
!  calculate energy
! ------------------------------------------------------------------------
dispersion_energy: if (opt%lenergy) then
   call edisp(mol%n,ndim,mol%at,q,mol%xyz,dparam,opt%g_a,opt%g_c, &
              gweights,refc6,opt%lmbd,energy,etwo=etwo,emany=emany)
endif dispersion_energy

! ------------------------------------------------------------------------
!  calculate gradient
! ------------------------------------------------------------------------
dispersion_gradient: if (opt%lgradient) then
   call dispgrad(mol%n,ndim,mol%at,q,mol%xyz,dparam,opt%wf,opt%g_a,opt%g_c, &
   &             refc6,opt%lmbd,gradient,etmp,dqdr)
   if (.not.opt%lenergy) energy = etmp

endif dispersion_gradient

end subroutine d4_calculation

!> interface to the DFT-D4 module for periodic boundary conditions
subroutine d4_pbc_calculation(iunit,opt,mol,dparam,energy,gradient,latgrad)
   use iso_fortran_env, wp => real64
!$ use omp_lib

! ------------------------------------------------------------------------
!  class definitions
! ------------------------------------------------------------------------
   use tbdef_options
   use tbdef_param
   use tbdef_molecule

! ------------------------------------------------------------------------
!  interfaces
! ------------------------------------------------------------------------
   use ncoord
   use eeq_model
   use dftd4
   use pbc_tools
   use pbc, only : get_realspace_cutoff

   implicit none

   !> output unit, usually bound to STDOUT
   !  not used if opt%silent is set (maybe for fatal errors)
   integer, intent(in) :: iunit

! ------------------------------------------------------------------------
!  class declarations
! ------------------------------------------------------------------------
   !> molecule structure information
   type(tb_molecule), intent(in) :: mol
   !> calculation options for DFT-D4 method
   type(dftd_options), intent(in) :: opt
   !> damping parameter for the DFT-D4 method
   type(dftd_parameter),intent(in) :: dparam
   type(chrg_parameter)            :: chrgeq

! ------------------------------------------------------------------------
!  output variables
! ------------------------------------------------------------------------
   !> dispersion energy, always referenced
   real(wp),intent(out) :: energy
   !> nuclear gradient, only references if opt%lgradient is set
   real(wp),intent(out) :: gradient(3,mol%n)
   !> nuclear hessian, only references if opt%lhessian is set
   real(wp),intent(out) :: latgrad(3,3)

! ------------------------------------------------------------------------
!  local variables
! ------------------------------------------------------------------------
   integer  :: ndim                      ! matrix dimension
   integer  :: i,j,k,l,ii,jj
   integer  :: err
   integer  :: rep_cn(3),rep_vdw(3),rep_mbd(3)
   real(wp) :: memory
   real(wp) :: etmp,etwo,emany,er,el,es
   real(wp) :: molpol,molc6,molc8        ! molecular Polarizibility
   real(wp) :: sigma(3,3),inv_lat(3,3)
   real(wp),allocatable :: q(:)          ! partial charges
   real(wp),allocatable :: dqdr(:,:,:)   ! partial charges
   real(wp),allocatable :: dqdL(:,:,:)   ! partial charges
   real(wp),allocatable :: covcn(:)      ! covalent coordination number
   real(wp),allocatable :: dcovcndr(:,:,:) ! covalent coordination number derivative
   real(wp),allocatable :: dcovcndL(:,:,:) ! covalent coordination number derivative
   real(wp),allocatable :: cn(:)         ! erf coordination number
   real(wp),allocatable :: dcndr(:,:,:)  ! erf coordination number derivative
   real(wp),allocatable :: dcndL(:,:,:)  ! erf coordination number derivative
   real(wp),allocatable :: gweights(:)   ! gaussian weights
   real(wp),allocatable :: refc6(:,:)    ! reference C6 coeffients
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: aw(:,:)
   real(wp),allocatable :: gtmp(:,:)
   real(wp),allocatable :: stmp(:,:)

   real(wp),parameter   :: rthr_cn  =  900.0_wp
   real(wp),parameter   :: rthr_mbd = 1600.0_wp
   real(wp),parameter   :: rthr_vdw = 4000.0_wp

! ------------------------------------------------------------------------
!  Output: Initialization and Parameter setup
! ------------------------------------------------------------------------
   call d4init(mol%n,mol%at,opt%g_a,opt%g_c,p_refq_goedecker,ndim)
   memory = ((mol%n)+(mol%n)+(ndim)+(ndim*ndim) &
            +(mol%n*mol%n)+(23*mol%n)+(mol%n)+(3*mol%n) &
            +(3*mol%n*mol%n)+(3*mol%n*(mol%n+1)) &
            +(3*3*mol%n)+(3*3*(mol%n+1))+(3*3*mol%n) &
            +(3*mol%n*mol%n)) * wp / (1024.0_wp**2)

   call get_realspace_cutoff(mol%lattice,rthr_cn, rep_cn)
   call get_realspace_cutoff(mol%lattice,rthr_vdw,rep_vdw)
   call get_realspace_cutoff(mol%lattice,rthr_mbd,rep_mbd)

   call generic_header(iunit,'Calculation Setup',49,10)
   write(iunit,'(3x,a," : ",i6)')   'number of atoms     ', mol%n
   write(iunit,'(3x,a," : ",i6)')   'charge              ', nint(mol%chrg)
   write(iunit,'(3x,a," : ",a)')    'non-additivity corr.', lmbd2string(opt%lmbd)
   write(iunit,'(3x,a," : ",a)')    'charge model        ', refq2string(p_refq_goedecker)
!$ write(iunit,'(3x,a," : ",i6)')   'omp threads         ',omp_get_num_threads()
   write(iunit,'(3x,a," : ",f6.1,1x,a)') &
      "memory needed (est.)",memory,"Mb"
   if (opt%wf  /= 6.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f6.1)')   'weighting factor    ', opt%wf
   if (opt%g_a /= 3.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f6.1)')   'q-scale height      ', opt%g_a
   if (opt%g_c /= 2.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f6.1)')   'q-scale steepness   ', opt%g_c
   write(iunit,'(3x,a," : ",i0,2(" × ",i0))') 'supercell for CN    ',2*rep_cn +1
   write(iunit,'(3x,a," : ",i0,2(" × ",i0))') 'supercell for disp. ',2*rep_vdw+1
   write(iunit,'(3x,a," : ",i0,2(" × ",i0))') 'supercell for ATM   ',2*rep_mbd+1

   write(iunit,'(a)')

   !if (.not.opt%silent) &
   !call prd4ref(mol)

   allocate( q(mol%n),covcn(mol%n),gweights(ndim),refc6(ndim,ndim),&
             c6ab(mol%n,mol%n),aw(23,mol%n),cn(mol%n), &
             dcndr(3,mol%n,mol%n),dcndL(3,3,mol%n), &
             dqdr(3,mol%n,mol%n+1),dqdL(3,3,mol%n+1), &
             dcovcndr(3,mol%n,mol%n), dcovcndL(3,3,mol%n), &
             source = 0.0_wp, stat = err )
   if (err /= 0) then
      call raise('E','Memory allocation failed',1)
      return
   endif
   sigma = 0.0_wp

   call get_d4_cn(mol,covcn,dcovcndr,dcovcndL,thr=rthr_cn)
   call pbc_d4(mol%n,ndim,mol%at,opt%wf,opt%g_a,opt%g_c,covcn,gweights,refc6)

   call get_erf_cn(mol,cn,dcndr,dcndL,thr=rthr_cn)

! ------------------------------------------------------------------------
!  get partial charges
! ------------------------------------------------------------------------
   !if (opt%verbose) &
   !call eeq_header
   call new_charge_model_2019(chrgeq,mol%n,mol%at)
   call eeq_chrgeq(mol,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,es,gtmp,stmp, &
                   .false.,.false.,.true.)
   !if (opt%verbose) &
   !call print_chrgeq(iunit,chrgeq,mol,q,cn)

! ------------------------------------------------------------------------
!  calculate properties
! ------------------------------------------------------------------------
dispersion_properties: if (opt%lmolpol) then
   call generic_header(iunit,'Molecular Properties',49,10)
   call mdisp(mol%n,ndim,mol%at,q,mol%xyz,opt%g_a,opt%g_c,gweights,refc6, &
      &       molc6,molc8,molpol,aw,c6ab)
   call prmolc6(molc6,molc8,molpol,mol%n,mol%at, &
      &         covcn=covcn,q=q,c6ab=c6ab,alpha=aw(1,:))
endif dispersion_properties

if (.not.opt%silent.and.(opt%lenergy.or.opt%lgradient)) then
   call generic_header(iunit,'Damping Parameters',49,10)
   write(iunit,'(3x,a," : ",f10.4)')   's6                  ', dparam%s6
   write(iunit,'(3x,a," : ",f10.4)')   's8                  ', dparam%s8
   if (dparam%s10 /= 0.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f10.4)')   's10                 ', dparam%s10
   if (dparam%s9  /= 1.0_wp .or. opt%verbose) &
   write(iunit,'(3x,a," : ",f10.4)')   's9                  ', dparam%s9
   write(iunit,'(3x,a," : ",f10.4)')   'a1                  ', dparam%a1
   write(iunit,'(3x,a," : ",f10.4)')   'a2                  ', dparam%a2
   write(iunit,'(a)')
endif

! ------------------------------------------------------------------------
!  calculate energy
! ------------------------------------------------------------------------
dispersion_energy: if (opt%lenergy) then
   call edisp_3d(mol,ndim,q,rep_vdw,rep_mbd,rthr_vdw,rthr_mbd,dparam, &
      &          opt%g_a,opt%g_c,gweights,refc6,opt%lmbd,energy, &
      &          etwo=etwo,embd=emany)
endif dispersion_energy

! ------------------------------------------------------------------------
!  calculate gradient
! ------------------------------------------------------------------------
dispersion_gradient: if (opt%lgradient) then
   call dispgrad_3d(mol,ndim,q,covcn,dcovcndr,dcovcndL,rep_vdw,rep_mbd, &
      &             rthr_vdw,rthr_mbd,dparam,opt%wf,opt%g_a,opt%g_c, &
      &             refc6,opt%lmbd,gradient,sigma,etmp,dqdr,dqdL)
   if (.not.opt%lenergy) energy = etmp
   inv_lat = mat_inv_3x3(mol%lattice)
   call sigma_to_latgrad(sigma,inv_lat,latgrad)
endif dispersion_gradient

end subroutine d4_pbc_calculation
end module tb_calculators
