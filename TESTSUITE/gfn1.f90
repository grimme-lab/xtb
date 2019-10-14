subroutine test_gfn1_scc
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
   logical  :: okpar,okbas,diff

   gfn_method = 1

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp
   call mol%set_nuclear_charge
   call mol%update

   wfn%nel = idint(sum(mol%z))
   wfn%nopen = 0

   allocate( g(3,mol%n), source = 0.0_wp )

   call use_parameterset('.param_gfn.xtb',globpar,okpar)
   call assert(okpar)

   call set_gfn1_parameter(param,globpar,mol%n,mol%at)

   call xbasis0(mol%n,mol%at,basis)
   call xbasis_gfn1(mol%n,mol%at,basis,okbas,diff)
   call assert(okbas)

   call xbasis_cao2sao(mol%n,mol%at,basis)

   call assert_eq(basis%nshell,6)
   call assert_eq(basis%nbf,   8)
   call assert_eq(basis%nao,   8)

   call wfn%allocate(mol%n,basis%nshell,basis%nao)
   wfn%q = mol%chrg/real(mol%n,wp)

   call iniqshell(mol%n,mol%at,mol%z,basis%nshell,wfn%q,wfn%qsh,gfn_method)

   g = 0.0_wp

   call scf(output_unit,mol,wfn,basis,param,pcem, &
      &   egap,et,maxiter,prlevel,restart,lgrad,acc,etot,g,res)

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
   use iso_fortran_env, wp => real64, istdout => output_unit

   use assertion

   use tbdef_options
   use tbdef_molecule
   use tbdef_param
   use tbdef_pcem
   use tbdef_wavefunction

   use tb_calculators

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

   type(tb_molecule)    :: mol
   type(tb_environment) :: env
   type(gfn_parameter)  :: gfn
   type(tb_pcem)        :: pcem
   type(tb_wavefunction):: wfn

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

   call gfn1_calculation &
      (istdout,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

   call assert_close(hl_gap, 5.6067613075402_wp,thr)
   call assert_close(energy,-8.4156335932985_wp,thr)
   call assert_close(norm2(gradient),0.95790240549503E-02_wp,thr)

   call assert_close(gradient(1,5), 0.18116310445596E-02_wp,thr)
   call assert_close(gradient(2,2), 0.00000000000000E+00_wp,thr)
   call assert_close(gradient(1,4),-0.18116310445594E-02_wp,thr)
   call assert_close(gradient(3,6),-0.76256813454808E-03_wp,thr)

   call terminate(afail)

end subroutine test_gfn1_api

subroutine test_gfn1gbsa_api
   use iso_fortran_env, wp => real64, istdout => output_unit

   use assertion

   use tbdef_options
   use tbdef_molecule
   use tbdef_param
   use tbdef_pcem
   use tbdef_wavefunction

   use tb_calculators

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

   type(tb_molecule)    :: mol
   type(tb_environment) :: env
   type(gfn_parameter)  :: gfn
   type(tb_pcem)        :: pcem
   type(tb_wavefunction):: wfn

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

   call gfn1_calculation &
      (istdout,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

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
   use iso_fortran_env, wp => real64, istdout => output_unit

   use assertion

   use tbdef_options
   use tbdef_molecule
   use tbdef_param
   use tbdef_pcem
   use tbdef_wavefunction

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
   type(tb_environment) :: env
   type(gfn_parameter)  :: gfn
   type(tb_pcem)        :: pcem
   type(tb_wavefunction):: wfn

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

   call gfn1_calculation &
      (istdout,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

   call assert_close(hl_gap, 9.0155275960407_wp,thr*10)
   call assert_close(energy,-23.113490916186_wp,thr)
   call assert_close(norm2(gradient),0.11143014174684E-01_wp,thr)

   call assert_close(gradient(1,5),-0.17083259496397E-02_wp,thr)
   call assert_close(gradient(2,2), 0.17083259496398E-02_wp,thr)
   call assert_close(gradient(1,4), 0.27992413917387E-02_wp,thr)
   call assert_close(gradient(3,6),-0.10966149569550E-02_wp,thr)

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
   ! gam from aoparam is now filled with GFN1-xTB hardnesses
   pcem%gam = gam(at(nat2+1:))
   pcem%q   = q
   pcem%grd = 0.0_wp

   call gfn1_calculation &
      (istdout,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

   call assert_close(hl_gap, 8.7253450666347_wp,thr)
   call assert_close(energy,-11.559896105984_wp,thr)
   call assert_close(norm2(gradient),0.24211484942219E-01_wp,thr)

   call assert_close(gradient(1,5),-0.16368872701817E-02_wp,thr)
   call assert_close(gradient(2,2), 0.34511834874966E-02_wp,thr)
   call assert_close(gradient(1,4),-0.27846344196073E-02_wp,thr)
   call assert_close(gradient(3,6),-0.12174494093948E-02_wp,thr)

   call assert_close(norm2(pcem%grd),0.12965281862178E-01_wp,thr)
   call assert_close(pcem%grd(1,5),-0.65532598920592E-03_wp,thr)
   call assert_close(pcem%grd(2,2), 0.17820246510446E-02_wp,thr)
   call assert_close(pcem%grd(1,4), 0.60888638785130E-02_wp,thr)
   call assert_close(pcem%grd(3,6), 0.86753094430381E-03_wp,thr)

   ! reset
   energy = 0.0_wp
   gradient = 0.0_wp
   pcem%grd = 0.0_wp
   pcem%gam = 999.0_wp ! point charges

   call gfn1_calculation &
      (istdout,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

   call assert_close(hl_gap, 8.9183046297437_wp,thr)
   call assert_close(energy,-11.565012263827_wp,thr)
   call assert_close(norm2(gradient),0.23134284179991E-01_wp,thr)

   call assert_close(gradient(1,5),-0.63000570838230E-03_wp,thr)
   call assert_close(gradient(2,2), 0.28673054867063E-02_wp,thr)
   call assert_close(gradient(1,4),-0.75488974649673E-02_wp,thr)
   call assert_close(gradient(3,6),-0.12128428341685E-02_wp,thr)

   call assert_close(norm2(pcem%grd),0.18251544072073E-01_wp,thr)
   call assert_close(pcem%grd(1,5),-0.16079631752423E-02_wp,thr)
   call assert_close(pcem%grd(2,2), 0.23749979339001E-02_wp,thr)
   call assert_close(pcem%grd(1,4), 0.10140193991067E-01_wp,thr)
   call assert_close(pcem%grd(3,6), 0.11833638475792E-02_wp,thr)

   call terminate(afail)

end subroutine test_gfn1_pcem_api
