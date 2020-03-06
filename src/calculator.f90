! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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
module xtb_calculators
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment
   implicit none
   private

   public :: gfn0_calculation, gfn1_calculation, gfn2_calculation
   public :: d4_calculation, d4_pbc_calculation

   interface
      module subroutine gfn0_calculation &
            (iunit,env,err,opt,mol,hl_gap,energy,gradient,stress,lattice_gradient)
         use xtb_mctc_logging
         use xtb_type_options
         use xtb_type_molecule
         use xtb_type_wavefunction
         use xtb_type_basisset
         use xtb_type_param
         use xtb_type_data
         integer, intent(in) :: iunit
         type(TMolecule),    intent(inout) :: mol
         type(peeq_options),   intent(in)    :: opt
         type(TEnvironment), intent(in)    :: env
         type(mctc_error), allocatable, intent(inout) :: err
         real(wp), intent(out) :: energy
         real(wp), intent(out) :: hl_gap
         real(wp), intent(out) :: gradient(3,mol%n)
         real(wp), intent(out) :: stress(3,3)
         real(wp), intent(out) :: lattice_gradient(3,3)
      end subroutine gfn0_calculation
      module subroutine gfn1_calculation &
            (iunit,env,err,opt,mol,pcem,wfn,hl_gap,energy,gradient)
         use xtb_mctc_logging
         use xtb_type_options
         use xtb_type_molecule
         use xtb_type_wavefunction
         use xtb_type_basisset
         use xtb_type_param
         use xtb_type_data
         use xtb_type_pcem
         integer, intent(in) :: iunit
         type(TMolecule),    intent(inout) :: mol
         type(scc_options),    intent(in)    :: opt
         type(TEnvironment), intent(in)    :: env
         type(mctc_error), allocatable, intent(inout) :: err
         type(tb_pcem),        intent(inout) :: pcem
         type(TWavefunction),intent(inout) :: wfn
         real(wp), intent(out) :: energy
         real(wp), intent(out) :: hl_gap
         real(wp), intent(out) :: gradient(3,mol%n)
      end subroutine gfn1_calculation
      module subroutine gfn2_calculation &
            (iunit,env,err,opt,mol,pcem,wfn,hl_gap,energy,gradient)
         use xtb_mctc_logging
         use xtb_type_options
         use xtb_type_molecule
         use xtb_type_wavefunction
         use xtb_type_basisset
         use xtb_type_param
         use xtb_type_data
         use xtb_type_pcem
         integer, intent(in) :: iunit
         type(TMolecule),    intent(inout) :: mol
         type(TWavefunction),intent(inout) :: wfn
         type(scc_options),    intent(in)    :: opt
         type(TEnvironment), intent(in)    :: env
         type(mctc_error), allocatable, intent(inout) :: err
         type(tb_pcem),        intent(inout) :: pcem
         real(wp), intent(out) :: energy
         real(wp), intent(out) :: hl_gap
         real(wp), intent(out) :: gradient(3,mol%n)
      end subroutine gfn2_calculation
   end interface
contains



!> interface to the DFT-D4 module
subroutine d4_calculation(iunit,opt,mol,dparam,energy,gradient)
   use xtb_mctc_accuracy, only : wp
!$ use omp_lib

! ------------------------------------------------------------------------
!  class definitions
! ------------------------------------------------------------------------
   use xtb_mctc_logging
   use xtb_type_options
   use xtb_type_param
   use xtb_type_molecule

! ------------------------------------------------------------------------
!  interfaces
! ------------------------------------------------------------------------
   use xtb_disp_ncoord
   use xtb_eeq
   use xtb_disp_dftd4
   !use dfuncpar

   implicit none

   !> output unit, usually bound to STDOUT
   !  not used if opt%silent is set (maybe for fatal errors)
   integer, intent(in) :: iunit

! ------------------------------------------------------------------------
!  class declarations
! ------------------------------------------------------------------------
   !> molecule structure information
   type(TMolecule), intent(inout) :: mol
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
   integer  :: stat
   real(wp) :: memory
   real(wp) :: etmp,etwo,emany,er,el,es
   real(wp) :: molpol,molc6,molc8        ! molecular Polarizibility
   real(wp) :: sigma(3,3)
   type(mctc_error), allocatable :: err
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
   call d4init(opt%g_a,opt%g_c,p_refq_goedecker)
   call d4dim(mol%n,mol%at,ndim)
   memory = ((mol%n)+(mol%n)+(ndim)+(ndim*ndim) &
            +(mol%n*mol%n)+(23*mol%n)+(mol%n)+(3*mol%n) &
            +(3*mol%n*mol%n)+(3*mol%n*(mol%n+1)) &
            +(3*mol%n*mol%n)) * wp / (1024.0_wp**2)

   call generic_header(iunit,'Calculation Setup',49,10)
   write(iunit,'(3x,a," : ",i6)')   'number of atoms     ', mol%n
   write(iunit,'(3x,a," : ",i6)')   'charge              ', nint(mol%chrg)
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
             dcovcndr(3,mol%n,mol%n), stat = stat )
   if (stat /= 0) then
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
   call eeq_chrgeq(mol,err,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
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
   use xtb_mctc_accuracy, only : wp
!$ use omp_lib

! ------------------------------------------------------------------------
!  class definitions
! ------------------------------------------------------------------------
   use xtb_mctc_logging
   use xtb_type_options
   use xtb_type_param
   use xtb_type_molecule

! ------------------------------------------------------------------------
!  interfaces
! ------------------------------------------------------------------------
   use xtb_disp_ncoord
   use xtb_eeq
   use xtb_disp_dftd4
   use xtb_pbc_tools
   use xtb_pbc, only : get_realspace_cutoff

   implicit none

   !> output unit, usually bound to STDOUT
   !  not used if opt%silent is set (maybe for fatal errors)
   integer, intent(in) :: iunit

! ------------------------------------------------------------------------
!  class declarations
! ------------------------------------------------------------------------
   !> molecule structure information
   type(TMolecule), intent(in) :: mol
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
   type(mctc_error), allocatable :: err
   integer  :: ndim                      ! matrix dimension
   integer  :: i,j,k,l,ii,jj
   integer  :: stat
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
   call d4init(opt%g_a,opt%g_c,p_refq_goedecker)
   call d4dim(mol%n,mol%at,ndim)
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
             source = 0.0_wp, stat = stat )
   if (stat /= 0) then
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
   call eeq_chrgeq(mol,err,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,es,gtmp,stmp, &
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
end module xtb_calculators
