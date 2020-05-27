subroutine test_gfn1_scc
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
   use xtb_type_solvent

   use xtb_setparam
   use xtb_basis
   use xtb_scf
   use xtb_scc_core
   use xtb_paramset
   use xtb_xtb_data
   use xtb_xtb_gfn1

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
   type(TSolvent), allocatable :: gbsa

   real(wp) :: etot,egap
   real(wp), allocatable :: g(:,:)

   type(TxTBParameter) :: globpar
   logical  :: okpar,okbas,diff,exitRun

   gfn_method = 1
   call init(env)

   call init(mol, at, xyz)

   wfn%nel = idint(sum(mol%z))
   wfn%nopen = 0

   allocate( g(3,mol%n), source = 0.0_wp )

   call use_parameterset('param_gfn1-xtb.txt',globpar,xtbData,okpar)
   call assert(okpar)

   call newBasisset(xtbData,mol%n,mol%at,basis,okbas)
   call assert(okbas)

   call assert_eq(basis%nshell,6)
   call assert_eq(basis%nbf,   8)
   call assert_eq(basis%nao,   8)

   call wfn%allocate(mol%n,basis%nshell,basis%nao)
   wfn%q = mol%chrg/real(mol%n,wp)

   call iniqshell(xtbData,mol%n,mol%at,mol%z,basis%nshell,wfn%q,wfn%qsh,gfn_method)

   g = 0.0_wp

   call scf(env,mol,wfn,basis,pcem,xtbData,gbsa, &
      &   egap,et,maxiter,prlevel,restart,lgrad,acc,etot,g,res)

   call env%check(exitRun)
   call assert(.not.exitRun)

   call assert(res%converged)

   call assert_close(res%e_total,-5.768694907835_wp,thr)
   call assert_close(res%gnorm,   0.006990472552_wp,thr*10)
   ! value in electron volt
   call assert_close(res%hl_gap,  9.314571212134_wp,1.0e-4_wp)
   call assert_close(res%e_elec, -5.805476374506_wp,thr)
   call assert_close(res%e_es,    0.050114150158_wp,thr)
   call assert_close(res%e_disp, -0.000137581088_wp,thr)
   call assert_close(res%e_rep,   0.036919047759_wp,thr)

   call assert_close(wfn%q(2),0.33430923478162_wp,thr2)
   call assert_close(wfn%q(3),wfn%q(2),thr2)
   call assert_close(wfn%qsh(1),0.30697460552546_wp,thr2)
   call assert_close(wfn%qsh(4),-0.31974031743775E-01_wp,thr2)

   call assert_eq(wfn%ihomo,4)
   call assert_eq(wfn%ihomoa,wfn%ihomob)
   call assert_close(wfn%emo(wfn%ihomo),-13.612017816475_wp,thr2)
   call assert_close(wfn%focca(wfn%ihomo),wfn%foccb(wfn%ihomo),thr2)
   call assert_close(wfn%focc(wfn%ihomo),2.0_wp,thr2)

   call mol%deallocate
   call wfn%deallocate
   call basis%deallocate

   call terminate(afail)
end subroutine test_gfn1_scc

subroutine test_gfn1_api
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use assertion

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_data
   use xtb_type_wavefunction
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction

   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 7
   integer, parameter :: at(nat) = [6,6,6,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[0.00000000000000_wp, 0.00000000000000_wp, 0.00000000000000_wp, &
      & 0.00000000000000_wp, 0.00000000000000_wp,-2.45883087071663_wp, &
      & 0.00000000000000_wp, 0.00000000000000_wp, 2.45883087071663_wp, &
      &-1.23873299308709_wp,-1.23873299308709_wp,-3.52402313872811_wp, &
      & 1.23873299308709_wp, 1.23873299308709_wp,-3.52402313872811_wp, &
      &-1.23873299308709_wp, 1.23873299308709_wp, 3.52402313872811_wp, &
      & 1.23873299308709_wp,-1.23873299308709_wp, 3.52402313872811_wp], shape(xyz))
   type(scc_options),parameter :: opt = scc_options( &
      &  prlevel = 2, maxiter = 30, acc = 1.0_wp, etemp = 300.0_wp, grad = .true. )

   type(TMolecule)    :: mol
   type(TEnvironment) :: env
   type(TWavefunction):: wfn
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

   call newXTBCalculator(env, mol, calc, method=1)
   call newWavefunction(env, mol, calc, wfn)

   call calc%singlepoint(env, mol, wfn, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call assert_close(hl_gap, 5.6067613073468_wp,thr)
   call assert_close(energy,-8.4156335928615_wp,thr)
   call assert_close(norm2(gradient),0.95790240549503E-02_wp,thr)

   call assert_close(gradient(1,5), 0.18116310445596E-02_wp,thr)
   call assert_close(gradient(2,2), 0.00000000000000E+00_wp,thr)
   call assert_close(gradient(1,4),-0.18116310445594E-02_wp,thr)
   call assert_close(gradient(3,6),-0.76256813454808E-03_wp,thr)

   call terminate(afail)

end subroutine test_gfn1_api

subroutine test_gfn1gbsa_api
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use assertion

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_data
   use xtb_type_wavefunction
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction, addSolvationModel

   implicit none

   real(wp),parameter :: thr = 1.0e-7_wp
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
   type(TEnvironment) :: env
   type(TWavefunction):: wfn
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

   call newXTBCalculator(env, mol, calc, method=1)
   call newWavefunction(env, mol, calc, wfn)
   call addSolvationModel(env, calc, opt%solvent)

   call calc%singlepoint(env, mol, wfn, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call assert_close(hl_gap, 6.641641300724_wp,1e-4_wp)
   call assert_close(energy,-14.215790820910_wp,thr)
   call assert_close(norm2(gradient),0.14758139345468E-01_wp,thr)
   call assert_close(gradient(2,3),0.1002863160985e-01_wp,thr)
   call assert_close(gradient(3,5),-gradient(3,7),thr)
   call assert_close(gradient(1,7),-0.6983782950712e-03_wp,thr)
   call assert_close(gradient(3,8),0.9313074280892e-03_wp,thr)

   call terminate(afail)

end subroutine test_gfn1gbsa_api

subroutine test_gfn1_pcem_api
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use assertion

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_data
   use xtb_type_wavefunction
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction

   implicit none

   real(wp),parameter :: thr = 1.0e-9_wp
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
   type(scc_options),parameter :: opt = scc_options( &
      &  prlevel = 2, maxiter = 30, acc = 1.0_wp, etemp = 300.0_wp, grad = .true. )
   real(wp),parameter :: gam(nat2) = [&
      & 0.583349_wp,0.470099_wp,0.470099_wp, 0.583349_wp,0.470099_wp,0.470099_wp]

   type(TMolecule)    :: mol
   type(TEnvironment) :: env
   type(TWavefunction):: wfn
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

   call newXTBCalculator(env, mol, calc, method=1)
   call newWavefunction(env, mol, calc, wfn)

   call calc%singlepoint(env, mol, wfn, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call assert_close(hl_gap, 9.0155275960407_wp,thr*10)
   call assert_close(energy,-23.113490914998_wp,thr)
   call assert_close(norm2(gradient),0.11143014174684E-01_wp,thr)

   call assert_close(gradient(1,5),-0.17083259496397E-02_wp,thr)
   call assert_close(gradient(2,2), 0.17083259496398E-02_wp,thr)
   call assert_close(gradient(1,4), 0.27992413917387E-02_wp,thr)
   call assert_close(gradient(3,6),-0.10966149569550E-02_wp,thr)

   ! reset
   call mol%deallocate
   energy = 0.0_wp
   gradient = 0.0_wp

   call init(mol, at(:nat2), xyz(:, :nat2))

   call newXTBCalculator(env, mol, calc, 'param_gfn1-xtb.txt', 1)

   call calc%pcem%allocate(nat2)
   calc%pcem%xyz = xyz(:,nat2+1:)
   ! gam from xtb_aoparam is now filled with GFN1-xTB hardnesses
   calc%pcem%gam = gam
   calc%pcem%q   = q
   calc%pcem%grd = 0.0_wp

   call newWavefunction(env, mol, calc, wfn)
   call calc%singlepoint(env, mol, wfn, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call assert_close(hl_gap, 8.7253450652232_wp,thr)
   call assert_close(energy,-11.559896105984_wp,thr)
   call assert_close(norm2(gradient),0.24211484942219E-01_wp,thr)

   call assert_close(gradient(1,5),-0.16368872701817E-02_wp,thr)
   call assert_close(gradient(2,2), 0.34511834874966E-02_wp,thr)
   call assert_close(gradient(1,4),-0.27846344196073E-02_wp,thr)
   call assert_close(gradient(3,6),-0.12174494093948E-02_wp,thr)

   call assert_close(norm2(calc%pcem%grd),0.12965281862178E-01_wp,thr)
   call assert_close(calc%pcem%grd(1,5),-0.65532598920592E-03_wp,thr)
   call assert_close(calc%pcem%grd(2,2), 0.17820246510446E-02_wp,thr)
   call assert_close(calc%pcem%grd(1,4), 0.60888638785130E-02_wp,thr)
   call assert_close(calc%pcem%grd(3,6), 0.86753094430381E-03_wp,thr)

   ! reset
   energy = 0.0_wp
   gradient = 0.0_wp
   calc%pcem%grd = 0.0_wp
   calc%pcem%gam = 999.0_wp ! point charges

   call newWavefunction(env, mol, calc, wfn)
   call calc%singlepoint(env, mol, wfn, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call assert_close(hl_gap, 8.9183046283326_wp,thr)
   call assert_close(energy,-11.565012263827_wp,thr)
   call assert_close(norm2(gradient),0.23134284179991E-01_wp,thr)

   call assert_close(gradient(1,5),-0.63000570838230E-03_wp,thr)
   call assert_close(gradient(2,2), 0.28673054867063E-02_wp,thr)
   call assert_close(gradient(1,4),-0.75488974649673E-02_wp,thr)
   call assert_close(gradient(3,6),-0.12128428341685E-02_wp,thr)

   call assert_close(norm2(calc%pcem%grd),0.18251544072073E-01_wp,thr)
   call assert_close(calc%pcem%grd(1,5),-0.16079631752423E-02_wp,thr)
   call assert_close(calc%pcem%grd(2,2), 0.23749979339001E-02_wp,thr)
   call assert_close(calc%pcem%grd(1,4), 0.10140193991067E-01_wp,thr)
   call assert_close(calc%pcem%grd(3,6), 0.11833638475792E-02_wp,thr)

   call terminate(afail)

end subroutine test_gfn1_pcem_api

subroutine test_gfn1_xb
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use assertion

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_data
   use xtb_type_wavefunction
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction

   implicit none

   real(wp),parameter :: thr = 1.0e-9_wp
   integer, parameter :: nat = 6
   integer, parameter :: at(nat) = [35,35,8,6,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape([&
      &-1.785333747_wp,    -3.126082999_wp,     0.000000000_wp, &
      & 0.000000000_wp,     0.816042264_wp,     0.000000000_wp, &
      & 2.658286999_wp,     5.297075806_wp,     0.000000000_wp, &
      & 4.885971586_wp,     4.861161373_wp,     0.000000000_wp, &
      & 5.615509753_wp,     2.908222159_wp,     0.000000000_wp, &
      & 6.289076126_wp,     6.399636435_wp,     0.000000000_wp], shape(xyz))
   type(scc_options),parameter :: opt = scc_options( &
      &  prlevel = 2, maxiter = 30, acc = 1.0_wp, etemp = 300.0_wp, grad = .true. )

   type(TMolecule)    :: mol
   type(TEnvironment) :: env
   type(TWavefunction):: wfn
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

   call newXTBCalculator(env, mol, calc, 'param_gfn1-xtb.txt', 1)
   call newWavefunction(env, mol, calc, wfn)

   call calc%singlepoint(env, mol, wfn, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call assert_close(hl_gap, 2.4991941159068_wp,thr)
   call assert_close(energy,-15.606233877972_wp,thr)
   call assert_close(norm2(gradient),0.23014320345408E-01_wp,thr)

   call assert_close(gradient(1,5),-0.39000047348209E-02_wp,thr)
   call assert_close(gradient(2,2),-0.49294645520340E-02_wp,thr)
   call assert_close(gradient(1,4), 0.17228152301357E-01_wp,thr)
   call assert_close(gradient(3,6), 0.00000000000000E+00_wp,thr)


   call terminate(afail)

end subroutine test_gfn1_xb


subroutine test_gfn1_pbc3d
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use assertion

   use xtb_mctc_convert
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data
   use xtb_type_environment
   use xtb_type_wavefunction

   use xtb_pbc_tools

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction

   implicit none

   real(wp),parameter :: thr = 1.0e-9_wp
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
   type(scc_options),parameter :: opt = scc_options( &
      &  prlevel = 2, maxiter = 30, acc = 1.0_wp, etemp = 300.0_wp, grad = .true. )

   type(TMolecule)    :: mol
   type(TEnvironment) :: env
   type(tb_pcem)      :: pcem
   type(TWavefunction):: wfn
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   real(wp) :: energy, sigma(3,3)
   real(wp) :: hl_gap
   real(wp) :: gradlatt(3,3)
   real(wp) :: stress(3,3)
   real(wp),allocatable :: gradient(:,:)
   real(wp),allocatable :: xyz(:,:)

   ! setup the environment variables
   call init(env)

   allocate(xyz(3, nat))
   call coord_trafo(nat,lattice,abc,xyz)
   call init(mol, at, xyz, lattice=lattice)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp
   gradlatt = 0.0_wp

   call mctc_mute

   call newXTBCalculator(env, mol, calc, 'param_gfn1-xtb.txt', 1)
   call newWavefunction(env, mol, calc, wfn)

   call calc%singlepoint(env, mol, wfn, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call assert_close(hl_gap, 7.5302549721955_wp,thr)
   call assert_close(energy,-11.069452578476_wp,thr)
   call assert_close(norm2(gradient), 0.28766497266274E-01_wp,thr)

   call assert_close(gradient(1,1),-0.18972454825631E-03_wp,thr)
   call assert_close(gradient(1,2),-0.66834750631031E-02_wp,thr)
   call assert_close(gradient(1,3), 0.68731996113594E-02_wp,thr)
   call assert_close(gradient(2,1),-0.26831102888425E-03_wp,thr)
   call assert_close(gradient(3,3), 0.16835831947879E-01_wp,thr)

   call terminate(afail)

end subroutine test_gfn1_pbc3d
