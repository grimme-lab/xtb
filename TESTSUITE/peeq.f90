subroutine test_peeq_sp
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use assertion

   use xtb_mctc_systools

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_wavefunction
   use xtb_type_basisset
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment

   use xtb_setparam, only : gfn_method
   use xtb_aoparam,  only : use_parameterset

   use xtb_pbc_tools
   use xtb_basis
   use xtb_peeq
   use xtb_readparam

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

   type(TMolecule)     :: mol
   type(TEnvironment)  :: env
   type(TWavefunction) :: wfn
   type(TBasisset)     :: basis
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
   call init(env)

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
   call print_pbcsum(stdout,mol)

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
         call env%error("Parameter file '"//fnv//"' not found!")
         call terminate(1)
         return
      endif
      call readParam(env,ipar,globpar,.true.)
      call close_file(ipar)
   endif
   call set_gfn0_parameter(param,globpar,mol%n,mol%at)
   call gfn0_prparam(stdout,mol%n,mol%at,param)

   call xbasis0(mol%n,mol%at,basis)
   call xbasis_gfn0(mol%n,mol%at,basis,okbas,diff)

   call wfn%allocate(mol%n,basis%nshell,basis%nao)
   wfn%nel = idint(sum(mol%z)) - mol%chrg
   wfn%nopen = 0
   wfn%q = mol%chrg / real(mol%n,wp)

   call mctc_mute

   call peeq(env,mol,wfn,basis,param,hl_gap,et,prlevel,lgrad,.true.,acc, &
      &      energy,gradient,sigma,res)

   call assert_close(energy,-7.3576550429483_wp,thr)
   call assert_close(hl_gap, 2.0722850435118_wp,1.0e-4_wp)
   call assert_close(norm2(gradient),4.6530106590558_wp,thr)
   call assert_close(gradient(1,1),0.0204787285706_wp,thr)
   call assert_close(gradient(3,2),3.1001342337520_wp,thr)

   call assert_close(sigma(1,1),-0.46949660136776_wp,thr)
   call assert_close(sigma(2,1),-0.34304476307923_wp,thr)
   call assert_close(sigma(2,3), -1.3504527773901_wp,thr)

   call assert_close(res%e_elec,-8.323348293826_wp,thr)
   call assert_close(res%e_es,  -0.106932060267_wp,thr)
   call assert_close(res%e_disp,-0.004089262287_wp,thr)
   call assert_close(res%e_rep,  1.076714573432_wp,thr)
   call assert_close(res%e_xb,   0.000000000000_wp,thr)

   ! reset for reevaluatuation without CCM
   energy = 0.0_wp
   gradient = 0.0_wp
   sigma = 0.0_wp

   call peeq(env,mol,wfn,basis,param,hl_gap,et,prlevel,lgrad,.false.,acc, &
      &      energy,gradient,sigma,res)

   call assert_close(energy,-7.3514777045762_wp,thr)
   call assert_close(hl_gap, 2.1721883949504_wp,1.0e-4_wp)
   call assert_close(norm2(gradient),4.6424961940410_wp,thr)
   call assert_close(gradient(1,1),0.0109632042166_wp,thr)
   call assert_close(gradient(3,2),3.0958397030189_wp,thr)

   call assert_close(sigma(1,1),-0.49612812087251_wp,thr)
   call assert_close(sigma(2,1),-0.32463429156960_wp,thr)
   call assert_close(sigma(2,3), -1.3397723168730_wp,thr)

   call assert_close(res%e_elec,-8.317170955429_wp,thr)
   call assert_close(res%e_es,  -0.106932060267_wp,thr)
   call assert_close(res%e_disp,-0.004089262287_wp,thr)
   call assert_close(res%e_rep,  1.076714573432_wp,thr)
   call assert_close(res%e_xb,   0.000000000000_wp,thr)

   call terminate(afail)

end subroutine test_peeq_sp

subroutine test_peeq_api
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use assertion

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_environment

   use xtb_pbc_tools

   use xtb_calculators

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

   type(TMolecule)    :: mol
   type(TEnvironment) :: env

   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp) :: gradlatt(3,3)
   real(wp) :: stress(3,3)
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call init(env)

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
      (stdout,env,opt,mol,hl_gap,energy,gradient,stress,gradlatt)

   call assert_close(hl_gap, 4.8620892163953_wp,thr)
   call assert_close(energy,-8.4898922181241_wp,thr)
   call assert_close(norm2(gradient),0.00000000000000E+00_wp,thr)
   call assert_close(norm2(gradlatt),0.45059748320564E-02_wp,thr)

   call terminate(afail)

end subroutine test_peeq_api

subroutine test_peeq_api_srb
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use assertion

   use xtb_mctc_convert
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_environment

   use xtb_pbc_tools

   use xtb_calculators

   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
   ! CaF2
   integer, parameter :: nat = 32
   integer, parameter :: at(nat) = [spread(8, 1, 4), spread(6, 1, 12),  spread(1, 1, 16)]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[4.5853168464880421_wp,  4.9392326929575878_wp,  4.1894081210748118_wp,  &
      & 5.8862267152491423_wp,  1.1425258978245871_wp,  6.5015058204768126_wp,  &
      & 1.4284279616220412_wp,  4.6017511285540875_wp,  3.1465884436348119_wp,  &
      & 2.3323404704521411_wp,  1.4471801154820869_wp,  0.7121932185858125_wp,  &
      & 4.7333543155561415_wp,  2.7747291872305868_wp,  3.1951352976178122_wp,  &
      & 5.6617754419101418_wp,  2.2133191164485870_wp,  2.1235404838618126_wp,  &
      & 5.3107598618381422_wp,  2.6902988056185868_wp,  0.7384466319968125_wp,  &
      & 4.4947071071761426_wp,  3.9530790692635867_wp,  0.6801776747598124_wp,  &
      & 4.8005171923760424_wp,  4.9185102874975870_wp,  1.8186363449528122_wp,  &
      & 4.6951362687070421_wp,  4.2781752812835867_wp,  3.1816411821728123_wp,  &
      & 1.3838574419160412_wp,  5.2817805008910863_wp,  4.1482702947948136_wp,  &
      & 1.0268974195990415_wp,  4.7234752637800881_wp,  5.4989995400388123_wp,  &
      & 2.0852659694760409_wp,  5.0956317453800875_wp,  6.5351699846458127_wp,  &
      & 2.3344644666691412_wp, -0.0736561690909131_wp,  6.4245628001158135_wp,  &
      & 2.4894017448231409_wp,  0.6213510313930869_wp,  5.0967297417158131_wp,  &
      & 1.5745272273791413_wp,  0.1243470825760870_wp,  3.9731040773988129_wp,  &
      & 5.8221065925130420_wp,  5.3013563342055878_wp,  1.7264876737078123_wp,  &
      & 3.4487807319551416_wp,  3.6355832152975864_wp,  0.7429568016758125_wp,  &
      & 4.8499393376520423_wp,  3.4713855169305874_wp,  6.4691872586348129_wp,  &
      & 0.2495364434351412_wp,  2.4795455690160870_wp,  2.1043557230378123_wp,  &
      & 5.6691068338331423_wp,  1.1234174220755870_wp,  2.1414388326468128_wp,  &
      & 3.7072009289431418_wp,  2.4357632918535872_wp,  3.0094700999208119_wp,  &
      & 4.1414520030430415_wp,  5.7877262477775879_wp,  1.7803680119358125_wp,  &
      & 5.0142851411171421_wp,  2.4165926460955873_wp,  4.1857610486448129_wp,  &
      & 3.0280930003030413_wp,  4.6201081184690871_wp,  6.2533190952188136_wp,  &
      & 0.5863628696651412_wp,  0.5757236365910867_wp,  4.1021714214668128_wp,  &
      & 2.3776130524831411_wp,  1.6969724987740866_wp,  5.2327688986668139_wp,  &
      & 1.9486148363011413_wp,  0.4390675147070869_wp,  2.9999022491838123_wp,  &
      & 3.5312997625581413_wp,  0.4467415528495868_wp,  4.8114121395028135_wp,  &
      & 6.5089895990100421_wp,  5.2100409408535882_wp,  6.0066553789008132_wp,  &
      & 0.9001165013630412_wp,  3.6420787128610868_wp,  5.4413106648508132_wp,  &
      & 1.6012116650460413_wp,  5.6845471271780879_wp,  0.7675566847298124_wp], &
      & shape(xyz)) * aatoau
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[6.4411018522600_wp,    0.0492571261505_wp,    0.2192046129910_wp,  &
      & 0.0462076831749_wp,    6.6435057067500_wp,    0.1670513770770_wp,  &
      & 0.2262248220170_wp,   -0.9573234940220_wp,    6.7608039126200_wp], &
      & shape(lattice)) * aatoau
   type(peeq_options),parameter :: opt = peeq_options( &
      &  prlevel = 2, ccm = .true., acc = 1.0_wp, etemp = 300.0_wp, grad = .true. )

   type(TMolecule)    :: mol
   type(TEnvironment) :: env

   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp) :: gradlatt(3,3)
   real(wp) :: stress(3,3)
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call init(env)

   call mol%allocate(nat)
   mol%at   = at
   mol%xyz  = xyz
   mol%npbc = 3
   mol%pbc  = .true.
   mol%lattice = lattice
   call mol%set_nuclear_charge
   call mol%update

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp
   gradlatt = 0.0_wp

   call mctc_mute

   call gfn0_calculation &
      (stdout,env,opt,mol,hl_gap,energy,gradient,stress,gradlatt)

   call assert_close(hl_gap, 3.3073195156202_wp,thr)
   call assert_close(energy,-47.310985782789_wp,thr)
   call assert_close(norm2(gradient),0.35661228276887E-01_wp,thr)
   call assert_close(norm2(gradlatt),0.19204614066516E-01_wp,thr)

   call assert_close(gradient(1, 3), 0.71669969919753E-03_wp,thr)
   call assert_close(gradient(2,11), 0.86141903688108E-02_wp,thr)
   call assert_close(gradient(1, 6),-0.23307620462240E-02_wp,thr)
   call assert_close(gradient(3, 5), 0.90733350522692E-03_wp,thr)
   call assert_close(gradient(1, 8), 0.30461459040208E-02_wp,thr)

   call assert_close(gradlatt(1,3),0.49893356364686E-03_wp,thr)
   call assert_close(gradlatt(2,2),0.94665150568032E-02_wp,thr)
   call assert_close(gradlatt(1,2),0.85097891190246E-03_wp,thr)

   call terminate(afail)

end subroutine test_peeq_api_srb
