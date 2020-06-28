subroutine test_gfn2_scc
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use assertion

   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_wavefunction
   use xtb_type_basisset
   use xtb_type_param
   use xtb_type_data
   use xtb_type_pcem
   use xtb_type_environment
   use xtb_type_solvation
   use xtb_solv_gbsa

   use xtb_setparam
   use xtb_basis
   use xtb_scf
   use xtb_scc_core
   use xtb_paramset

   use xtb_xtb_data
   use xtb_xtb_gfn2

   implicit none
   real(wp),parameter :: thr = 1.0e-7_wp
   real(wp),parameter :: thr2 = 1.0e-5_wp
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.00000000000000_wp,  0.00000000000000_wp, -0.73578586109551_wp, &
      & 1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp, &
      &-1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp  &
      & ],shape(xyz))
   real(wp),parameter :: et = 300.0_wp
   integer, parameter :: maxiter = 50
   integer, parameter :: prlevel = 2
   logical, parameter :: lgrad = .true.
   logical, parameter :: restart = .false.
   real(wp),parameter :: acc = 1.0_wp

   type(TEnvironment) :: env
   type(TMolecule)     :: mol
   type(scc_results)     :: res
   type(TBasisset)     :: basis
   type(TWavefunction) :: wfn
   type(tb_pcem)         :: pcem
   type(TxTBData) :: xtbData
   class(TSolvation), allocatable :: solvation

   real(wp) :: etot,egap
   real(wp), allocatable :: g(:,:)

   type(TxTBParameter) :: globpar
   logical  :: okpar,okbas
   logical :: exitRun

   gfn_method = 2
   call init(env)

   call init(mol, at, xyz)

   wfn%nel = idint(sum(mol%z))
   wfn%nopen = 0

   allocate( g(3,mol%n), source = 0.0_wp )

   call use_parameterset('param_gfn2-xtb.txt',globpar,xtbData,okpar)
   call assert(okpar)

   call newBasisset(xtbData,mol%n,mol%at,basis,okbas)
   call assert(okbas)

   call assert_eq(basis%nshell,4)
   call assert_eq(basis%nbf,   6)
   call assert_eq(basis%nao,   6)

   call wfn%allocate(mol%n,basis%nshell,basis%nao)
   wfn%q = mol%chrg/real(mol%n,wp)

   call iniqshell(xtbData,mol%n,mol%at,mol%z,basis%nshell,wfn%q,wfn%qsh,gfn_method)

   g = 0.0_wp

   call scf(env,mol,wfn,basis,pcem,xtbData,solvation, &
      &   egap,et,maxiter,prlevel,restart,lgrad,acc,etot,g,res)

   call env%check(exitRun)
   call assert(.not.exitRun)

   call assert(res%converged)

   call assert_close(res%e_total,-5.070451355118_wp,thr)
   call assert_close(res%gnorm,   0.006457420125_wp,thr)
   ! value in electron volt
   call assert_close(res%hl_gap, 14.450372368833_wp,1.0e-4_wp)
   call assert_close(res%e_elec, -5.104362813671_wp,thr)
   ! es and aes are not welldefined at this accuracy level
   call assert_close(res%e_es,    0.031408028914_wp,thr*10)
   call assert_close(res%e_aes,   0.000563300591_wp,thr*10)
   call assert_close(res%e_axc,  -0.000813611055_wp,thr)
   call assert_close(res%e_disp, -0.000141250480_wp,thr)
   call assert_close(res%e_rep,   0.033911458523_wp,thr)

   call assert_close(wfn%q(2),0.28158903353422_wp,thr2)
   call assert_close(wfn%q(3),wfn%q(2),thr2)
   call assert_close(wfn%qsh(1),0.25886041477578_wp,thr2)
   call assert_close(wfn%dipm(2,2), 0.00000000000000E+00_wp,thr2)
   call assert_close(wfn%dipm(1,3),-0.55238805565739E-01_wp,thr2)
   call assert_close(wfn%qp(6,1),0.41207363331933E-01_wp,thr2)
   call assert_close(wfn%qp(1,3),0.95084219046109E-01_wp,thr2)

   call assert_eq(wfn%ihomo,4)
   call assert_eq(wfn%ihomoa,wfn%ihomob)
   call assert_close(wfn%emo(wfn%ihomo),-12.166283951806_wp,thr2)
   call assert_close(wfn%focca(wfn%ihomo),wfn%foccb(wfn%ihomo),thr2)
   call assert_close(wfn%focc(wfn%ihomo),2.0_wp,thr2)

   call mol%deallocate
   call wfn%deallocate
   call basis%deallocate

   call terminate(afail)
end subroutine test_gfn2_scc

subroutine test_gfn2_api
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use assertion

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction

   implicit none

   real(wp),parameter :: thr = 1.0e-9_wp
   integer, parameter :: nat = 7
   integer, parameter :: at(nat) = [6,6,6,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[0.00000000000000_wp, 0.00000000000000_wp,-1.79755622305860_wp, &
      & 0.00000000000000_wp, 0.00000000000000_wp, 0.95338756106749_wp, &
      & 0.00000000000000_wp, 0.00000000000000_wp, 3.22281255790261_wp, &
      &-0.96412815539807_wp,-1.66991895015711_wp,-2.53624948351102_wp, &
      &-0.96412815539807_wp, 1.66991895015711_wp,-2.53624948351102_wp, &
      & 1.92825631079613_wp, 0.00000000000000_wp,-2.53624948351102_wp, &
      & 0.00000000000000_wp, 0.00000000000000_wp, 5.23010455462158_wp], shape(xyz))
   type(scc_options),parameter :: opt = scc_options( &
      &  prlevel = 2, maxiter = 30, acc = 1.0_wp, etemp = 300.0_wp, grad = .true.,&
      &  solvent = "none")

   type(TMolecule)    :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   real(wp) :: energy, sigma(3, 3)
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call init(env)

   call init(mol, at, xyz)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=2)
   call newWavefunction(env, mol, calc, chk)

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)
   call assert_close(hl_gap, 7.0005867526665_wp,thr)
   call assert_close(energy,-8.3824793849585_wp,thr)
   call assert_close(norm2(gradient),0.11544410028854E-01_wp,thr)

   call assert_close(gradient(2,3), 0.00000000000000E+00_wp,thr)
   call assert_close(gradient(3,1),-0.74649908147372E-03_wp,thr)
   call assert_close(gradient(1,4),-0.28433755158510E-03_wp,thr)
   call assert_close(gradient(3,7), 0.99750545315944E-02_wp,thr)

   call terminate(afail)

end subroutine test_gfn2_api

subroutine test_gfn2gbsa_api
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use assertion

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment
   use xtb_solv_input

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction, addSolvationModel

   implicit none

   real(wp),parameter :: thr = 5.0e-7_wp
   integer, parameter :: nat = 11
   integer, parameter :: at(nat) = [7,1,1,1,15,6,6,6,7,7,7]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[4.84112686339230_wp,-1.15671648374837_wp, 0.00000000000000_wp, &
      & 5.87217277985151_wp,-0.70323701391366_wp, 1.54511189584035_wp, &
      & 5.87217277985151_wp,-0.70323701391366_wp,-1.54511189584035_wp, &
      & 4.61871932563459_wp,-3.05696083228965_wp, 0.00000000000000_wp, &
      & 0.23751589819013_wp, 1.28722286021158_wp, 0.00000000000000_wp, &
      & 0.00161949529588_wp,-0.89075264835377_wp, 2.57562301980395_wp, &
      &-3.06580664031769_wp, 2.16487781673832_wp, 0.00000000000000_wp, &
      & 0.00161949529588_wp,-0.89075264835377_wp,-2.57562301980395_wp, &
      &-5.06585309511506_wp, 3.00827014874972_wp, 0.00000000000000_wp, &
      & 0.00161949529588_wp,-2.15043086045241_wp, 4.34183098346579_wp, &
      & 0.00161949529588_wp,-2.15043086045241_wp,-4.34183098346579_wp], shape(xyz))
   type(scc_options),parameter :: opt = scc_options( &
      &  prlevel = 2, maxiter = 30, acc = 1.0_wp, etemp = 300.0_wp, grad = .true.,&
      &  solvent = "h2o")

   type(TMolecule)    :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   real(wp) :: energy, sigma(3,3)
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call init(env)

   call init(mol, at, xyz)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=2)
   call newWavefunction(env, mol, calc, chk)
   call addSolvationModel(env, calc, TSolvInput(solvent=opt%solvent))

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call assert_close(hl_gap, 3.408607724814_wp,1e-5_wp)
   call assert_close(energy,-22.002501380096_wp,thr)
   call assert_close(norm2(gradient),0.19441977481008E-01_wp,thr)
   call assert_close(gradient(1,1), .9038674439127e-02_wp,thr)
   call assert_close(gradient(3,2),-.1394693523214e-02_wp,thr)
   call assert_close(gradient(3,11),-gradient(3,10),thr)
   call assert_close(gradient(1,8),0.22890674680144E-02_wp,thr)

   call terminate(afail)

end subroutine test_gfn2gbsa_api

subroutine test_gfn2salt_api
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_convert, only : aatoau

   use assertion

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment
   use xtb_solv_input

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction, addSolvationModel

   implicit none

   real(wp),parameter :: thr = 5.0e-7_wp
   integer, parameter :: nat = 8
   integer, parameter :: at(nat) = [7,7,9,1,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[0.00000000000000_wp, 1.30633555031959_wp, 0.00000000000000_wp,  &
      &-1.26402080275201_wp,-4.06871453465823_wp, 0.00000000000000_wp,  &
      & 1.42937939828404_wp, 3.57212474255400_wp, 0.00000000000000_wp,  &
      &-2.92500314615638_wp,-3.12058224200551_wp, 0.00000000000000_wp,  &
      & 0.71133449039308_wp, 0.35812010971703_wp, 1.51289962417918_wp,  &
      &-1.25697212427635_wp,-5.20406577474818_wp,-1.53711646457211_wp,  &
      & 0.71133449039308_wp, 0.35812010971703_wp,-1.51289962417918_wp,  &
      &-1.25697212427635_wp,-5.20406577474818_wp, 1.53711646457211_wp], shape(xyz))
   type(scc_options),parameter :: opt = scc_options( &
      &  prlevel = 2, maxiter = 30, acc = 1.0_wp, etemp = 300.0_wp, grad = .true.,&
      &  solvent = 'ch2cl2' )

   type(TMolecule)    :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   real(wp) :: energy, sigma(3,3)
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call init(env)

   call init(mol, at, xyz)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=2)
   call newWavefunction(env, mol, calc, chk)
   call addSolvationModel(env, calc, TSolvInput(solvent=opt%solvent, &
      & ionRad=1.0_wp*aatoau, ionStrength=1.0e-3_wp))

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call assert_close(hl_gap, 6.895830675032_wp,5e-5_wp)
   call assert_close(energy,-13.027106170796_wp,thr)
   call assert_close(norm2(gradient),0.014655633557_wp,thr)

   call assert_close(gradient(1,1),-0.6696141112967e-02_wp,thr)
   call assert_close(gradient(2,4),-0.7608884457863e-03_wp,thr)
   call assert_close(gradient(1,5), 0.2359793849677e-02_wp,thr)
   call assert_close(gradient(3,7), 0.6882248903260e-02_wp,thr)

   call terminate(afail)

end subroutine test_gfn2salt_api

subroutine test_gfn2_pcem_api
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use assertion

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_pcem
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction, addSolvationModel

   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 12, nat2 = nat/2
   integer, parameter :: at(nat) = [8,1,1, 8,1,1, 8,1,1, 8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape([&
      &-2.75237178376284_wp, 2.43247309226225_wp,-0.01392519847964_wp, &
      &-0.93157260886974_wp, 2.79621404458590_wp,-0.01863384029005_wp, &
      &-3.43820531288547_wp, 3.30583608421060_wp, 1.42134539425148_wp, &
      &-2.43247309226225_wp,-2.75237178376284_wp, 0.01392519847964_wp, &
      &-2.79621404458590_wp,-0.93157260886974_wp, 0.01863384029005_wp, &
      &-3.30583608421060_wp,-3.43820531288547_wp,-1.42134539425148_wp, &
      & 2.75237178376284_wp,-2.43247309226225_wp,-0.01392519847964_wp, &
      & 0.93157260886974_wp,-2.79621404458590_wp,-0.01863384029005_wp, &
      & 3.43820531288547_wp,-3.30583608421060_wp, 1.42134539425148_wp, &
      & 2.43247309226225_wp, 2.75237178376284_wp, 0.01392519847964_wp, &
      & 2.79621404458590_wp, 0.93157260886974_wp, 0.01863384029005_wp, &
      & 3.30583608421060_wp, 3.43820531288547_wp,-1.42134539425148_wp], shape(xyz))
   real(wp),parameter :: q(nat2) = [&
      &-0.69645733_wp,       0.36031084_wp,       0.33614649_wp, &
      &-0.69645733_wp,       0.36031084_wp,       0.33614649_wp]
   real(wp),parameter :: gam(nat2) = [&
      & 0.451896_wp,       0.405771_wp,       0.405771_wp, &
      & 0.451896_wp,       0.405771_wp,       0.405771_wp]
   type(scc_options),parameter :: opt = scc_options( &
      &  prlevel = 2, maxiter = 30, acc = 1.0_wp, etemp = 300.0_wp, grad = .true. )

   type(TMolecule)    :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   real(wp) :: energy, sigma(3,3)
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call init(env)

   call init(mol, at, xyz)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=2)
   call newWavefunction(env, mol, calc, chk)

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call assert_close(hl_gap, 12.391144584178_wp,thr)
   call assert_close(energy,-20.323978512117_wp,thr)
   call assert_close(norm2(gradient),0.78119239557115E-02_wp,thr)

   call assert_close(gradient(1,5),-0.22192122053513E-02_wp,thr)
   call assert_close(gradient(2,2), 0.22192122053512E-02_wp,thr)
   call assert_close(gradient(1,4), 0.95621597761913E-03_wp,thr)
   call assert_close(gradient(3,6),-0.11904153838296E-02_wp,thr)

   ! reset
   call mol%deallocate
   deallocate(gradient)

   call init(mol, at(:nat2), xyz(:, :nat2))
   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=2)

   call calc%pcem%allocate(nat2)
   calc%pcem%xyz = xyz(:,nat2+1:)
   ! gam from xtb_aoparam is now filled with GFN2-xTB hardnesses
   calc%pcem%gam = gam
   calc%pcem%q   = q
   calc%pcem%grd = 0.0_wp

   call newWavefunction(env, mol, calc, chk)
   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call assert_close(hl_gap, 12.718203165741_wp,thr)
   call assert_close(energy,-10.160927754235_wp,thr)
   call assert_close(norm2(gradient),0.21549557655285E-01_wp,thr)

   call assert_close(gradient(1,5),-0.20326749264006E-02_wp,thr)
   call assert_close(gradient(2,2), 0.33125459724368E-02_wp,thr)
   call assert_close(gradient(1,4),-0.11929405659085E-02_wp,thr)
   call assert_close(gradient(3,6),-0.16607682747438E-02_wp,thr)

   call assert_close(norm2(calc%pcem%grd),0.10043976337709E-01_wp,thr)
   call assert_close(calc%pcem%grd(1,5),-0.24831958012524E-03_wp,thr)
   call assert_close(calc%pcem%grd(2,2), 0.14208444722973E-02_wp,thr)
   call assert_close(calc%pcem%grd(1,4), 0.37466852704082E-02_wp,thr)
   call assert_close(calc%pcem%grd(3,6), 0.65161344334732E-03_wp,thr)

   ! reset
   energy = 0.0_wp
   gradient = 0.0_wp
   calc%pcem%grd = 0.0_wp
   calc%pcem%gam = 999.0_wp ! point charges

   call newWavefunction(env, mol, calc, chk)
   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call assert_close(hl_gap, 13.024345612330_wp,thr)
   call assert_close(energy,-10.168788268962_wp,thr)
   call assert_close(norm2(gradient),0.18113624400926E-01_wp,thr)

   call assert_close(gradient(1,5),-0.50646289499094E-03_wp,thr)
   call assert_close(gradient(2,2), 0.25468320656932E-02_wp,thr)
   call assert_close(gradient(1,4),-0.74927880093291E-02_wp,thr)
   call assert_close(gradient(3,6),-0.13248514654811E-02_wp,thr)

   call assert_close(norm2(calc%pcem%grd),0.18721791896294E-01_wp,thr)
   call assert_close(calc%pcem%grd(1,5),-0.21573703225712E-02_wp,thr)
   call assert_close(calc%pcem%grd(2,2), 0.25653662154150E-02_wp,thr)
   call assert_close(calc%pcem%grd(1,4), 0.12124342986218E-01_wp,thr)
   call assert_close(calc%pcem%grd(3,6), 0.12140575062433E-02_wp,thr)

   call terminate(afail)

end subroutine test_gfn2_pcem_api


subroutine test_gfn2_mindless_basic
   use assertion
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction

   implicit none

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)

   character(len=*), parameter :: mindless(10) = [&
      & "mindless01", "mindless02", "mindless03", "mindless04", "mindless05", &
      & "mindless06", "mindless07", "mindless08", "mindless09", "mindless10"]
   real(wp), parameter :: ref_energies(10) = &
      &[-30.348902328339_wp, -24.069929673504_wp, -23.718708414840_wp, &
      & -22.751497890244_wp, -27.735987617741_wp, -18.559531652263_wp, &
      & -33.423246445911_wp, -29.982410709899_wp, -20.549458882257_wp, &
      & -25.648093031152_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.06081938165963_wp, 0.05480379303387_wp, 0.04179957306749_wp, &
      & 0.06932568214054_wp, 0.05133081215262_wp, 0.05897000784363_wp, &
      & 0.04600837003153_wp, 0.06088718433105_wp, 0.05314451272641_wp, &
      & 0.04588800056975_wp]
   real(wp), parameter :: ref_hlgaps(10) = &
      &[2.2611412120964_wp, 1.2455369848896_wp, 0.1490345448030_wp, &
      & 1.2694246332861_wp, 2.2362028272626_wp, 2.6053844166160_wp, &
      & 2.6312636988491_wp, 0.8502228061534_wp, 3.2960732039022_wp, &
      & 0.4185654147579_wp]

   call init(env)
   do iMol = 1, 10
      if (afail > 0) exit

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newXTBCalculator(env, mol, calc, 'param_gfn2-xtb.txt', 2)
      call newWavefunction(env, mol, calc, chk)

      call env%check(exitRun)
      call assert(.not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call assert(.not.exitRun)
      if (exitRun) exit

      call assert_close(energy, ref_energies(iMol), thr)
      call assert_close(norm2(gradient), ref_gnorms(iMol), thr)
      call assert_close(hl_gap, ref_hlgaps(iMol), thr)

   end do

   call terminate(afail)

end subroutine test_gfn2_mindless_basic


subroutine test_gfn2_mindless_solvation
   use assertion
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction, addSolvationModel
   use xtb_solv_input, only : TSolvInput

   implicit none

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)

   character(len=*), parameter :: mindless(10) = [&
      & "mindless01", "mindless02", "mindless03", "mindless04", "mindless05", &
      & "mindless06", "mindless07", "mindless08", "mindless09", "mindless10"]
   character(len=*), parameter :: solvents(10) = [character(len=20) ::&
      & "h2o", "chcl3", "thf", "acetonitrile", "toluene", &
      & "ch2cl2", "ether", "methanol", "cs2", "dmso"]
   real(wp), parameter :: ref_energies(10) = &
      &[-30.384994095620_wp, -24.092790321708_wp, -23.760000587724_wp, &
      & -22.759033676135_wp, -27.755756860587_wp, -18.578906535415_wp, &
      & -33.448544467530_wp, -29.993777657049_wp, -20.567481283491_wp, &
      & -25.665458003353_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.064382736361151_wp, 0.055824257363022_wp, 0.040696787567357_wp, &
      & 0.073197307180171_wp, 0.049852029310233_wp, 0.053051084845551_wp, &
      & 0.046463198340453_wp, 0.057936724779825_wp, 0.055081716352177_wp, &
      & 0.045339090185313_wp]
   real(wp), parameter :: ref_hlgaps(10) = &
      &[2.5109978045214_wp, 1.3559515037838_wp, 0.1719065524538_wp, &
      & 1.2428812469179_wp, 2.1428443722148_wp, 2.2548300583012_wp, &
      & 2.6935601499212_wp, 0.8310005969097_wp, 3.3270919380043_wp, &
      & 0.3859217429034_wp]

   call init(env)
   do iMol = 1, 10
      if (afail > 0) exit

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newXTBCalculator(env, mol, calc, 'param_gfn2-xtb.txt', 2)
      call newWavefunction(env, mol, calc, chk)
      call addSolvationModel(env, calc, TSolvInput(solvent=trim(solvents(iMol)), &
         & alpb=mod(iMol, 2)==0))

      call env%check(exitRun)
      call assert(.not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call assert(.not.exitRun)
      if (exitRun) exit

      call assert_close(energy, ref_energies(iMol), thr)
      call assert_close(norm2(gradient), ref_gnorms(iMol), thr)
      call assert_close(hl_gap, ref_hlgaps(iMol), thr)

   end do

   call terminate(afail)

end subroutine test_gfn2_mindless_solvation


subroutine test_gfn2_dmetal
   use assertion
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction, addSolvationModel
   use xtb_solv_input, only : TSolvInput

   implicit none

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)

   real(wp), parameter :: ref_energies(3) = &
      &[-34.066853522474_wp, -34.078580189077_wp, -34.078529449408_wp]
   real(wp), parameter :: ref_gnorms(3) = &
      &[0.26985488354859_wp, 0.26591831790033_wp, 0.26613440125774_wp]
   real(wp), parameter :: ref_hlgaps(3) = &
      &[ 2.9239308006972_wp,  2.9163891345974_wp,  2.9168589130125_wp]

   call init(env)
   do iMol = 1, 3
      if (afail > 0) exit

      call getMolecule(mol, 'feco5')

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newXTBCalculator(env, mol, calc, 'param_gfn2-xtb.txt', 2)
      call newWavefunction(env, mol, calc, chk)
      if (iMol > 1) then
         call addSolvationModel(env, calc, TSolvInput(solvent='ch2cl2', &
            & alpb=iMol==3))
      end if

      call env%check(exitRun)
      call assert(.not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 1, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call assert(.not.exitRun)
      if (exitRun) exit

      call assert_close(energy, ref_energies(iMol), thr)
      call assert_close(norm2(gradient), ref_gnorms(iMol), thr)
      call assert_close(hl_gap, ref_hlgaps(iMol), thr)

   end do

   call terminate(afail)

end subroutine test_gfn2_dmetal
