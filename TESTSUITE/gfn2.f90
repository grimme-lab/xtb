subroutine test_gfn2_scc
   use iso_fortran_env, wp => real64

   use assertion

   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data
   use tbdef_pcem

   use setparam
   use aoparam
   use xbasis
   use scf_module
   use scc_core

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

   type(tb_molecule)     :: mol
   type(scc_results)     :: res
   type(tb_basisset)     :: basis
   type(tb_wavefunction) :: wfn
   type(scc_parameter)   :: param
   type(tb_pcem)         :: pcem

   real(wp) :: etot,egap
   real(wp), allocatable :: g(:,:)

   real(wp) :: globpar(25)
   logical  :: okpar,okbas

   gfn_method = 2

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp
   call mol%set_nuclear_charge
   call mol%update

   wfn%nel = idint(sum(mol%z))
   wfn%nopen = 0

   allocate( g(3,mol%n), source = 0.0_wp )
 
   call use_parameterset('.param_gfn2.xtb',globpar,okpar)
   call assert(okpar)

   call set_gfn2_parameter(param,globpar,mol%n,mol%at)

   call xbasis0(mol%n,mol%at,basis)
   call xbasis_gfn2(mol%n,mol%at,basis,okbas)
   call assert(okbas)

   call xbasis_cao2sao(mol%n,mol%at,basis)

   call assert_eq(basis%nshell,4)
   call assert_eq(basis%nbf,   6)
   call assert_eq(basis%nao,   6)

   call wfn%allocate(mol%n,basis%nshell,basis%nao)
   wfn%q = mol%chrg/real(mol%n,wp)

   call iniqshell(mol%n,mol%at,mol%z,basis%nshell,wfn%q,wfn%qsh,gfn_method)

   g = 0.0_wp

   call scf(output_unit,mol,wfn,basis,param,pcem, &
      &   egap,et,maxiter,prlevel,restart,lgrad,acc,etot,g,res)

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
   use iso_fortran_env, wp => real64, istdout => output_unit

   use assertion

   use tbdef_options
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_param
   use tbdef_pcem

   use tb_calculators

   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
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

   type(tb_molecule)    :: mol
   type(tb_wavefunction):: wfn
   type(tb_environment) :: env
   type(gfn_parameter)  :: gfn
   type(tb_pcem)        :: pcem

   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call env%setup

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   call mol%set_nuclear_charge
   call mol%update

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call gfn2_calculation &
      (istdout,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

   call assert_close(hl_gap, 7.0005867526665_wp,thr)
   call assert_close(energy,-8.3824793818504_wp,thr)
   call assert_close(norm2(gradient),0.11544410028854E-01_wp,thr)

   call assert_close(gradient(2,3), 0.00000000000000E+00_wp,thr)
   call assert_close(gradient(3,1),-0.74649908147372E-03_wp,thr)
   call assert_close(gradient(1,4),-0.28433755158510E-03_wp,thr)
   call assert_close(gradient(3,7), 0.99750545315944E-02_wp,thr)

   call terminate(afail)

end subroutine test_gfn2_api

subroutine test_gfn2gbsa_api
   use iso_fortran_env, wp => real64, istdout => output_unit

   use assertion

   use tbdef_options
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_param
   use tbdef_pcem

   use tb_calculators

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

   type(tb_molecule)    :: mol
   type(tb_wavefunction):: wfn
   type(tb_environment) :: env
   type(gfn_parameter)  :: gfn
   type(tb_pcem)        :: pcem

   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call env%setup

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   call mol%set_nuclear_charge
   call mol%update

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call gfn2_calculation &
      (istdout,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

   call assert_close(hl_gap, 3.408607724814_wp,1e-5_wp)
   call assert_close(energy,-22.002501380096_wp,thr)
   call assert_close(norm2(gradient),0.019344118754_wp,thr)
   call assert_close(gradient(1,1), .9038674439127e-02_wp,thr)
   call assert_close(gradient(3,2),-.1394693523214e-02_wp,thr)
   call assert_close(gradient(3,11),-gradient(3,10),thr)
   call assert_close(gradient(1,8),.1383384688560e-02_wp,thr)

   call terminate(afail)

end subroutine test_gfn2gbsa_api

subroutine test_gfn2salt_api
   use iso_fortran_env, wp => real64, istdout => output_unit

   use assertion

   use tbdef_options
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_param
   use tbdef_pcem

   use gbobc

   use tb_calculators

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

   type(tb_molecule)    :: mol
   type(tb_wavefunction):: wfn
   type(tb_environment) :: env
   type(gfn_parameter)  :: gfn
   type(tb_pcem)        :: pcem

   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call env%setup

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   call mol%set_nuclear_charge
   call mol%update

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   lsalt = .true.
   ion_rad = 1.0_wp
   ionst = 1.0e-3_wp

   call gfn2_calculation &
      (istdout,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

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
   use iso_fortran_env, wp => real64, istdout => output_unit

   use assertion

   use tbdef_options
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_param
   use tbdef_pcem

   use aoparam

   use tb_calculators

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
   type(scc_options),parameter :: opt = scc_options( &
      &  prlevel = 2, maxiter = 30, acc = 1.0_wp, etemp = 300.0_wp, grad = .true. )

   type(tb_molecule)    :: mol
   type(tb_wavefunction):: wfn
   type(tb_environment) :: env
   type(gfn_parameter)  :: gfn
   type(tb_pcem)        :: pcem

   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call env%setup

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   call mol%set_nuclear_charge
   call mol%update

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call gfn2_calculation &
      (istdout,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

   call assert_close(hl_gap, 12.391144583778_wp,thr)
   call assert_close(energy,-20.323978513218_wp,thr)
   call assert_close(norm2(gradient),0.78119239557115E-02_wp,thr)

   call assert_close(gradient(1,5),-0.22192122053513E-02_wp,thr)
   call assert_close(gradient(2,2), 0.22192122053512E-02_wp,thr)
   call assert_close(gradient(1,4), 0.95621597761913E-03_wp,thr)
   call assert_close(gradient(3,6),-0.11904153838296E-02_wp,thr)

   ! reset
   call mol%deallocate
   energy = 0.0_wp
   gradient = 0.0_wp

   call mol%allocate(nat2)
   mol%at  = at(:nat2)
   mol%xyz = xyz(:,:nat2)
   call mol%set_nuclear_charge
   call mol%update

   call pcem%allocate(nat2)
   pcem%xyz = xyz(:,nat2+1:)
   ! gam from aoparam is now filled with GFN2-xTB hardnesses
   pcem%gam = gam(at(nat2+1:))
   pcem%q   = q
   pcem%grd = 0.0_wp

   call gfn2_calculation &
      (istdout,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

   call assert_close(hl_gap, 12.718203165741_wp,thr)
   call assert_close(energy,-10.160927752124_wp,thr)
   call assert_close(norm2(gradient),0.21549557655285E-01_wp,thr)

   call assert_close(gradient(1,5),-0.20326749264006E-02_wp,thr)
   call assert_close(gradient(2,2), 0.33125459724368E-02_wp,thr)
   call assert_close(gradient(1,4),-0.11929405659085E-02_wp,thr)
   call assert_close(gradient(3,6),-0.16607682747438E-02_wp,thr)

   call assert_close(norm2(pcem%grd),0.10043976337709E-01_wp,thr)
   call assert_close(pcem%grd(1,5),-0.24831958012524E-03_wp,thr)
   call assert_close(pcem%grd(2,2), 0.14208444722973E-02_wp,thr)
   call assert_close(pcem%grd(1,4), 0.37466852704082E-02_wp,thr)
   call assert_close(pcem%grd(3,6), 0.65161344334732E-03_wp,thr)

   ! reset
   energy = 0.0_wp
   gradient = 0.0_wp
   pcem%grd = 0.0_wp
   pcem%gam = 999.0_wp ! point charges

   call gfn2_calculation &
      (istdout,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

   call assert_close(hl_gap, 13.024345612330_wp,thr)
   call assert_close(energy,-10.168788269555_wp,thr)
   call assert_close(norm2(gradient),0.18113624400926E-01_wp,thr)

   call assert_close(gradient(1,5),-0.50646289499094E-03_wp,thr)
   call assert_close(gradient(2,2), 0.25468320656932E-02_wp,thr)
   call assert_close(gradient(1,4),-0.74927880093291E-02_wp,thr)
   call assert_close(gradient(3,6),-0.13248514654811E-02_wp,thr)

   call assert_close(norm2(pcem%grd),0.18721791896294E-01_wp,thr)
   call assert_close(pcem%grd(1,5),-0.21573703225712E-02_wp,thr)
   call assert_close(pcem%grd(2,2), 0.25653662154150E-02_wp,thr)
   call assert_close(pcem%grd(1,4), 0.12124342986218E-01_wp,thr)
   call assert_close(pcem%grd(3,6), 0.12140575062433E-02_wp,thr)

   call terminate(afail)

end subroutine test_gfn2_pcem_api
