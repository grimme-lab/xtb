subroutine test_peeq_sp
   use iso_fortran_env, wp => real64, istdout => output_unit
   use assertion

   use mctc_systools

   use tbdef_options
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data

   use setparam, only : gfn_method
   use aoparam,  only : use_parameterset

   use pbc_tools
   use xbasis
   use peeq_module

   implicit none

   real(wp),parameter :: thr = 1.0e-7_wp
   ! SiO2
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,14,8]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[4.51379029_wp,    8.29324256_wp,   -0.60000000_wp, &
      & 0.56666667_wp,    0.56666667_wp,   -0.60000000_wp, &
      & 0.96666667_wp,    0.96666667_wp,    1.20000001_wp],&
      & shape(xyz))
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[6.04712363_wp,    0.00000000_wp,    0.00000000_wp,    &
      & 0.00000000_wp,    9.82657590_wp,    0.00000000_wp,    &
      & 0.00000000_wp,    0.00000000_wp,   13.22808290_wp],   &
      & shape(lattice))
   real(wp),parameter :: et = 300.0_wp
   integer, parameter :: prlevel = 2
   logical, parameter :: lgrad = .true.
   logical, parameter :: ccm = .true.
   logical, parameter :: restart = .false.
   real(wp),parameter :: acc = 1.0_wp

   type(tb_molecule)     :: mol
   type(gfn_parameter)   :: gfn
   type(tb_environment)  :: env
   type(tb_wavefunction) :: wfn
   type(tb_basisset)     :: basis
   type(scc_parameter)   :: param
   type(scc_results)     :: res

   real(wp)              :: energy
   real(wp)              :: hl_gap
   real(wp), allocatable :: gradient(:,:)
   real(wp)              :: sigma(3,3)

   integer, parameter    :: wsc_rep(3) = [1,1,1] ! FIXME

   character(len=*), parameter   :: p_fnv_gfn0 = '.param_gfn0.xtb'
   character(len=:), allocatable :: fnv
   real(wp) :: globpar(25)
   integer  :: ipar
   logical  :: exist

   logical  :: okbas,diff

   gfn_method = 0
   call env%setup

   call mol%allocate(nat)
   mol%at   = at
   mol%xyz  = xyz
   mol%npbc = 3
   mol%pbc  = .true.
   mol%lattice = lattice
   call mol%update
   call mol%set_nuclear_charge
   call generate_wsc(mol,mol%wsc,wsc_rep)

   allocate( gradient(3,mol%n) )
   energy = 0.0_wp
   gradient = 0.0_wp
   sigma = 0.0_wp

   ! give an optional summary on the geometry used
   call print_pbcsum(istdout,mol)

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
         call raise('E',"Parameter file '"//fnv//"' not found!")
      endif
      call read_gfn_param(ipar,globpar,.true.)
      call close_file(ipar)
   endif
   call set_gfn0_parameter(param,globpar,mol%n,mol%at)
   call gfn0_prparam(istdout,mol%n,mol%at,param)

   call xbasis0(mol%n,mol%at,basis)
   call xbasis_gfn0(mol%n,mol%at,basis,okbas,diff)
   call xbasis_cao2sao(mol%n,mol%at,basis)

   call wfn%allocate(mol%n,basis%nshell,basis%nao)
   wfn%nel = idint(sum(mol%z)) - mol%chrg
   wfn%nopen = 0
   wfn%q = mol%chrg / real(mol%n,wp)

   call mctc_mute

   call peeq(output_unit,mol,wfn,basis,param,hl_gap,et,prlevel,lgrad,.true.,acc, &
      &      energy,gradient,sigma,res)

   call assert_close(energy,-7.3591722560663_wp,thr)
   call assert_close(hl_gap, 2.0722850435118_wp,1.0e-4_wp)
   call assert_close(norm2(gradient),4.6530280011052_wp,thr)
   call assert_close(gradient(1,1),0.0204805697292_wp,thr)
   call assert_close(gradient(3,2),3.1001489574039_wp,thr)

   call assert_close(sigma(1,1),-0.46780039586900_wp,thr)
   call assert_close(sigma(2,1),-0.34304914894764_wp,thr)
   call assert_close(sigma(2,3), -1.3504530498122_wp,thr)

   call assert_close(res%e_elec,-8.323348293826_wp,thr)
   call assert_close(res%e_es,  -0.106932060267_wp,thr)
   call assert_close(res%e_disp,-0.005606475405_wp,thr)
   call assert_close(res%e_rep,  1.076714573432_wp,thr)
   call assert_close(res%e_xb,   0.000000000000_wp,thr)

   ! reset for reevulatuation without CCM
   energy = 0.0_wp
   gradient = 0.0_wp
   sigma = 0.0_wp

   call peeq(output_unit,mol,wfn,basis,param,hl_gap,et,prlevel,lgrad,.false.,acc, &
      &      energy,gradient,sigma,res)

   call assert_close(energy,-7.3529949176686_wp,thr)
   call assert_close(hl_gap, 2.1721883949504_wp,1.0e-4_wp)
   call assert_close(norm2(gradient),4.6425135475316_wp,thr)
   call assert_close(gradient(1,1),0.0109650455537_wp,thr)
   call assert_close(gradient(3,2),3.0958544267552_wp,thr)

   call assert_close(sigma(1,1),-0.49443191637757_wp,thr)
   call assert_close(sigma(2,1),-0.32463867696413_wp,thr)
   call assert_close(sigma(2,3), -1.3397725889616_wp,thr)

   call assert_close(res%e_elec,-8.317170955429_wp,thr)
   call assert_close(res%e_es,  -0.106932060267_wp,thr)
   call assert_close(res%e_disp,-0.005606475405_wp,thr)
   call assert_close(res%e_rep,  1.076714573432_wp,thr)
   call assert_close(res%e_xb,   0.000000000000_wp,thr)

   call terminate(afail)

end subroutine test_peeq_sp

subroutine test_peeq_api
   use iso_fortran_env, wp => real64, istdout => output_unit
   use assertion

   use tbdef_options
   use tbdef_molecule
   use tbdef_param

   use pbc_tools

   use tb_calculators

   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
   ! CaF2
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [9,9,20]
   real(wp),parameter :: abc(3,nat) = reshape(&
      &[0.25_wp, 0.25_wp, 0.25_wp, &
      & 0.75_wp, 0.75_wp, 0.75_wp, &
      & 0.00_wp, 0.00_wp, 0.00_wp], shape(abc))
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[5.9598811567890_wp,      2.1071361905157_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      6.3214085715472_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      0.0000000000000_wp,      7.2993338808807_wp],   &
      & shape(lattice))
   type(peeq_options),parameter :: opt = peeq_options( &
      &  prlevel = 2, ccm = .true., acc = 1.0_wp, etemp = 300.0_wp, grad = .true. )

   type(tb_molecule)    :: mol
   type(tb_environment) :: env
   type(gfn_parameter)  :: gfn

   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp) :: gradlatt(3,3)
   real(wp) :: stress(3,3)
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call env%setup

   call mol%allocate(nat)
   mol%at   = at
   mol%abc  = abc
   mol%npbc = 3
   mol%pbc  = .true.
   mol%lattice = lattice
   mol%volume = dlat_to_dvol(lattice)
   call dlat_to_cell(lattice,mol%cellpar)
   call dlat_to_rlat(lattice,mol%rec_lat)
   call coord_trafo(nat,lattice,abc,mol%xyz)
   call mol%set_nuclear_charge
   call mol%update

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp
   gradlatt = 0.0_wp

   call mctc_mute

   call gfn0_calculation &
      (istdout,env,opt,mol,gfn,hl_gap,energy,gradient,stress,gradlatt)

   call assert_close(hl_gap, 4.8620892163953_wp,thr)
   call assert_close(energy,-8.4930019025474_wp,thr)
   call assert_close(norm2(gradient),0.00000000000000E+00_wp,thr)
   call assert_close(norm2(gradlatt),0.54550838330373E-02_wp,thr)

   call terminate(afail)

end subroutine test_peeq_api
