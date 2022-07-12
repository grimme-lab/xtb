! This file is part of xtb.
! SPDX-Identifier: LGPL-3.0-or-later
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

module test_gfn1
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check, test_failed
   implicit none
   private

   public :: collect_gfn1

contains

!> Collect all exported unit tests
subroutine collect_gfn1(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("scc", test_gfn1_scc), &
      new_unittest("api", test_gfn1_api), &
      new_unittest("gbsa", test_gfn1gbsa_api), &
      new_unittest("pcem", test_gfn1_pcem_api), &
      new_unittest("xb", test_gfn1_xb), &
      new_unittest("pbc3d", test_gfn1_pbc3d), &
      new_unittest("mindless-basic", test_gfn1_mindless_basic), &
      new_unittest("mindless-solvation", test_gfn1_mindless_solvation), &
      new_unittest("ipea-indole", test_ipea_indole), &
      new_unittest("mindless-cosmo", test_gfn1_mindless_cosmo) &
      ]

end subroutine collect_gfn1


subroutine test_gfn1_scc(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

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
   use xtb_xtb_gfn1

   type(error_type), allocatable, intent(out) :: error
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
   logical  :: okpar,okbas,diff,exitRun

   set%gfn_method = 1
   call init(env)

   call init(mol, at, xyz)

   wfn%nel = idint(sum(mol%z))
   wfn%nopen = 0

   allocate( g(3,mol%n), source = 0.0_wp )

   call use_parameterset('param_gfn1-xtb.txt',globpar,xtbData,okpar)
   call check_(error, okpar)

   call newBasisset(xtbData,mol%n,mol%at,basis,okbas)
   call check_(error, okbas)

   call check_(error, basis%nshell,6)
   call check_(error, basis%nbf,   8)
   call check_(error, basis%nao,   8)

   call wfn%allocate(mol%n,basis%nshell,basis%nao)
   wfn%q = mol%chrg/real(mol%n,wp)

   call iniqshell(xtbData,mol%n,mol%at,mol%z,basis%nshell,wfn%q,wfn%qsh,set%gfn_method)

   g = 0.0_wp

   call scf(env,mol,wfn,basis,pcem,xtbData,solvation, &
      &   egap,et,maxiter,prlevel,restart,lgrad,acc,etot,g,res)

   call env%check(exitRun)
   call check_(error, .not.exitRun)

   call check_(error, res%converged)

   call check_(error, res%e_total,-5.768694907835_wp, thr=thr)
   call check_(error, res%gnorm,   0.006990472552_wp, thr=thr*10)
   ! value in electron volt
   call check_(error, res%hl_gap,  9.314571212134_wp, thr=1.0e-4_wp)
   call check_(error, res%e_elec, -5.805476374506_wp, thr=thr)
   call check_(error, res%e_es,    0.050114150158_wp, thr=thr)
   call check_(error, res%e_disp, -0.000137581088_wp, thr=thr)
   call check_(error, res%e_rep,   0.036919047759_wp, thr=thr)

   call check_(error, wfn%q(2),0.33430923478162_wp, thr=thr2)
   call check_(error, wfn%q(3),wfn%q(2), thr=thr2)
   call check_(error, wfn%qsh(1),0.30697460552546_wp, thr=thr2)
   call check_(error, wfn%qsh(4),-0.31974031743775E-01_wp, thr=thr2)

   call check_(error, wfn%ihomo,4)
   call check_(error, wfn%ihomoa,wfn%ihomob)
   call check_(error, wfn%emo(wfn%ihomo),-13.612017816475_wp, thr=thr2)
   call check_(error, wfn%focca(wfn%ihomo),wfn%foccb(wfn%ihomo), thr=thr2)
   call check_(error, wfn%focc(wfn%ihomo),2.0_wp, thr=thr2)

   call mol%deallocate
   call wfn%deallocate
   call basis%deallocate

end subroutine test_gfn1_scc

subroutine test_gfn1_api(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_data
   use xtb_type_restart
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   type(error_type), allocatable, intent(out) :: error

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
   type(TRestart) :: wfn
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

   call check_(error, hl_gap, 5.6067613073468_wp, thr=thr)
   call check_(error, energy,-8.4156335928615_wp, thr=thr)
   call check_(error, norm2(gradient),0.95790240549503E-02_wp, thr=thr)

   call check_(error, gradient(1,5), 0.18116310445596E-02_wp, thr=thr)
   call check_(error, gradient(2,2), 0.00000000000000E+00_wp, thr=thr)
   call check_(error, gradient(1,4),-0.18116310445594E-02_wp, thr=thr)
   call check_(error, gradient(3,6),-0.76256813454808E-03_wp, thr=thr)

end subroutine test_gfn1_api

subroutine test_gfn1gbsa_api(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_data
   use xtb_type_restart
   use xtb_type_environment
   use xtb_solv_input

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_kernel, only : gbKernel

   type(error_type), allocatable, intent(out) :: error

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
   type(TRestart):: wfn
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
   call addSolvationModel(env, calc, TSolvInput(solvent=opt%solvent, alpb=.false., kernel=gbKernel%still))

   call calc%singlepoint(env, mol, wfn, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call check_(error, hl_gap, 6.641641300724_wp, thr=1e-4_wp)
   call check_(error, energy,-14.215790820910_wp, thr=thr)
   call check_(error, norm2(gradient),0.14758139345468E-01_wp, thr=thr)
   call check_(error, gradient(2,3),0.1002863160985e-01_wp, thr=thr)
   call check_(error, gradient(3,5),-gradient(3,7), thr=thr)
   call check_(error, gradient(1,7),-0.6983782950712e-03_wp, thr=thr)
   call check_(error, gradient(3,8),0.9313074280892e-03_wp, thr=thr)

end subroutine test_gfn1gbsa_api

subroutine test_gfn1_pcem_api(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_data
   use xtb_type_restart
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   type(error_type), allocatable, intent(out) :: error

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
   type(TRestart) :: wfn
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

   call check_(error, hl_gap, 9.0155275960407_wp, thr=thr*10)
   call check_(error, energy,-23.113490914998_wp, thr=thr)
   call check_(error, norm2(gradient),0.11143014174684E-01_wp, thr=thr)

   call check_(error, gradient(1,5),-0.17083259496397E-02_wp, thr=thr)
   call check_(error, gradient(2,2), 0.17083259496398E-02_wp, thr=thr)
   call check_(error, gradient(1,4), 0.27992413917387E-02_wp, thr=thr)
   call check_(error, gradient(3,6),-0.10966149569550E-02_wp, thr=thr)

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

   call check_(error, hl_gap, 8.7253450652232_wp, thr=thr)
   call check_(error, energy,-11.559896105984_wp, thr=thr)
   call check_(error, norm2(gradient),0.24211484942219E-01_wp, thr=thr)

   call check_(error, gradient(1,5),-0.16368872701817E-02_wp, thr=thr)
   call check_(error, gradient(2,2), 0.34511834874966E-02_wp, thr=thr)
   call check_(error, gradient(1,4),-0.27846344196073E-02_wp, thr=thr)
   call check_(error, gradient(3,6),-0.12174494093948E-02_wp, thr=thr)

   call check_(error, norm2(calc%pcem%grd),0.12965281862178E-01_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,5),-0.65532598920592E-03_wp, thr=thr)
   call check_(error, calc%pcem%grd(2,2), 0.17820246510446E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,4), 0.60888638785130E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(3,6), 0.86753094430381E-03_wp, thr=thr)

   ! reset
   energy = 0.0_wp
   gradient = 0.0_wp
   calc%pcem%grd = 0.0_wp
   calc%pcem%gam = 999.0_wp ! point charges

   call newWavefunction(env, mol, calc, wfn)
   call calc%singlepoint(env, mol, wfn, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call check_(error, hl_gap, 8.9183046283326_wp, thr=thr)
   call check_(error, energy,-11.565012263827_wp, thr=thr)
   call check_(error, norm2(gradient),0.23134284179991E-01_wp, thr=thr)

   call check_(error, gradient(1,5),-0.63000570838230E-03_wp, thr=thr)
   call check_(error, gradient(2,2), 0.28673054867063E-02_wp, thr=thr)
   call check_(error, gradient(1,4),-0.75488974649673E-02_wp, thr=thr)
   call check_(error, gradient(3,6),-0.12128428341685E-02_wp, thr=thr)

   call check_(error, norm2(calc%pcem%grd),0.18251544072073E-01_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,5),-0.16079631752423E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(2,2), 0.23749979339001E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,4), 0.10140193991067E-01_wp, thr=thr)
   call check_(error, calc%pcem%grd(3,6), 0.11833638475792E-02_wp, thr=thr)

end subroutine test_gfn1_pcem_api

subroutine test_gfn1_xb(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_data
   use xtb_type_restart
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   type(error_type), allocatable, intent(out) :: error

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
   type(TRestart) :: wfn
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

   call check_(error, hl_gap, 2.4991941159068_wp, thr=thr)
   call check_(error, energy,-15.606233877972_wp, thr=thr)
   call check_(error, norm2(gradient),0.23014320345408E-01_wp, thr=thr)

   call check_(error, gradient(1,5),-0.39000047348209E-02_wp, thr=thr)
   call check_(error, gradient(2,2),-0.49294645520340E-02_wp, thr=thr)
   call check_(error, gradient(1,4), 0.17228152301357E-01_wp, thr=thr)
   call check_(error, gradient(3,6), 0.00000000000000E+00_wp, thr=thr)

end subroutine test_gfn1_xb


subroutine test_gfn1_pbc3d(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_mctc_convert
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data
   use xtb_type_environment
   use xtb_type_restart

   use xtb_pbc_tools

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   type(error_type), allocatable, intent(out) :: error

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
   type(TRestart) :: wfn
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

   call check_(error, hl_gap, 7.5302549612743_wp, thr=thr)
   call check_(error, energy,-11.069452578476_wp, thr=thr)
   call check_(error, norm2(gradient), 0.28766497266274E-01_wp, thr=thr)

   call check_(error, gradient(1,1),-0.18972454825631E-03_wp, thr=thr)
   call check_(error, gradient(1,2),-0.66834750631031E-02_wp, thr=thr)
   call check_(error, gradient(1,3), 0.68731996113594E-02_wp, thr=thr)
   call check_(error, gradient(2,1),-0.26831102888425E-03_wp, thr=thr)
   call check_(error, gradient(3,3), 0.16835831947879E-01_wp, thr=thr)

end subroutine test_gfn1_pbc3d


subroutine test_gfn1_mindless_basic(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   type(error_type), allocatable, intent(out) :: error

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
      &[-33.040345103560_wp, -26.858637887392_wp, -25.818131312372_wp, &
      & -24.896623131166_wp, -29.038403257541_wp, -20.596982571512_wp, &
      & -35.754243805311_wp, -33.088749099589_wp, -21.414820157390_wp, &
      & -26.650774557476_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.049510354228196_wp, 0.062164738458454_wp, 0.050970505747783_wp, &
      & 0.089329244570940_wp, 0.049188889391502_wp, 0.059353809155596_wp, &
      & 0.063932193893054_wp, 0.059154356097925_wp, 0.043307936447560_wp, &
      & 0.043009212776427_wp]
   real(wp), parameter :: ref_hlgaps(10) = &
      &[3.4192962150687_wp, 1.6702597448498_wp, 2.4899572196060_wp, &
      & 1.6706491472352_wp, 2.4613273978074_wp, 3.9778768185046_wp, &
      & 2.5370236462819_wp, 1.3590769445407_wp, 4.0816020166740_wp, &
      & 1.3640382002752_wp]

   call init(env)
   do iMol = 1, 10

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newXTBCalculator(env, mol, calc, 'param_gfn1-xtb.txt', 1)
      call newWavefunction(env, mol, calc, chk)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)
      call check_(error, hl_gap, ref_hlgaps(iMol), thr=thr)

   end do

end subroutine test_gfn1_mindless_basic


subroutine test_gfn1_mindless_solvation(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel

   type(error_type), allocatable, intent(out) :: error

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
      &[-33.075972501756_wp, -26.890550325491_wp, -25.875403448921_wp, &
      & -24.905567384699_wp, -29.056610710932_wp, -20.626564707704_wp, &
      & -35.787363117115_wp, -33.103345507512_wp, -21.425126566837_wp, &
      & -26.683086381488_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.058132809354921_wp, 0.061771913881691_wp, 0.056320324680277_wp, &
      & 0.094108546452594_wp, 0.049107097494380_wp, 0.056792563143050_wp, &
      & 0.063361086949237_wp, 0.059227312312587_wp, 0.041889795619126_wp, &
      & 0.045018190836542_wp]
   real(wp), parameter :: ref_hlgaps(10) = &
      &[3.5484557435287_wp, 1.6696997912942_wp, 2.5331364339285_wp, &
      & 1.6610531980301_wp, 2.3956600668633_wp, 3.7354634449786_wp, &
      & 2.6533136247624_wp, 1.3546253158640_wp, 4.0956047417239_wp, &
      & 1.3095666818230_wp]


   call init(env)
   do iMol = 1, 10

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newXTBCalculator(env, mol, calc, 'param_gfn1-xtb.txt', 1)
      call newWavefunction(env, mol, calc, chk)
      call addSolvationModel(env, calc, TSolvInput(solvent=trim(solvents(iMol)), &
         & alpb=mod(iMol, 2)==0, kernel=gbKernel%still))

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)
      call check_(error, hl_gap, ref_hlgaps(iMol), thr=thr)

   end do

end subroutine test_gfn1_mindless_solvation


subroutine test_ipea_indole(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule
   use xtb_mctc_systools, only : rdpath

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp

   integer, parameter :: nat = 15
   character(len=*), parameter :: sym(nat) = [character(len=4)::&
      & "C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "C", "C", "N", "H", &
      & "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & -5.35892196512242_wp,  1.73312414112064_wp, -0.57717892384120_wp, &
      & -7.65881264141790_wp,  0.49826397815934_wp, -0.98802419062021_wp, &
      & -3.18596416481952_wp,  0.23723616573036_wp, -0.37947581998019_wp, &
      & -3.24547162745174_wp, -2.39135340069686_wp, -0.57504353378802_wp, &
      & -5.58372373567590_wp, -3.59229298714979_wp, -0.98874228639030_wp, &
      & -7.77510636165883_wp, -2.12452413033580_wp, -1.19171772498436_wp, &
      & -9.38333853055452_wp,  1.60452830650258_wp, -1.15099413591722_wp, &
      & -5.68917043033699_wp, -5.63501089146626_wp, -1.14933117729174_wp, &
      & -9.58383842846746_wp, -3.04470716793929_wp, -1.51155380153258_wp, &
      & -5.27586852013397_wp,  3.77546410029495_wp, -0.41965138859121_wp, &
      & -0.72015557111823_wp, -3.26425550101682_wp, -0.28287304164551_wp, &
      &  0.80122480411036_wp, -1.18587847253807_wp,  0.07759213768431_wp, &
      & -0.70611490908719_wp,  0.91783977772213_wp,  0.01598707951313_wp, &
      & -0.09943736690084_wp,  2.71570592169382_wp,  0.22824107134708_wp, &
      & -0.09051786154600_wp, -5.20876436014445_wp, -0.33266731412436_wp],&
      & shape(xyz))
   real(wp), parameter :: charge = -1.0_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun, exist
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)
   character(len=:), allocatable :: paramfile

   call init(env)
   call init(mol, sym, xyz, chrg=charge)

   allocate(gradient(3, len(mol)))

   call rdpath(env%xtbpath, 'param_ipea-xtb.txt', paramfile, exist)
   if (.not.exist) paramfile = 'param_ipea-xtb.txt'

   call newXTBCalculator(env, mol, calc, paramfile, 1)
   call env%check(exitRun)
   call check_(error, .not.exitRun)
   if (.not.exitRun) then

      call newWavefunction(env, mol, calc, chk)

   end if

   call env%check(exitRun)
   call check_(error, .not.exitRun)
   if (.not.exitRun) then

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

   end if

   call env%check(exitRun)
   call check_(error, .not.exitRun)
   if (.not.exitRun) then

      call check_(error, energy, -26.590861716652_wp, thr=thr)
      call check_(error, norm2(gradient), 0.84641833840045E-01_wp, thr=thr)
      call check_(error, hl_gap, 2.5955961749533_wp, thr=thr)

   end if

end subroutine test_ipea_indole


subroutine test_gfn1_mindless_cosmo(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput

   type(error_type), allocatable, intent(out) :: error

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
      &[-33.086826283006_wp,-26.892421337119_wp,-25.862366080135_wp, &
      & -24.929812068243_wp,-29.068535381601_wp,-20.631250827831_wp, &
      & -35.780608195840_wp,-33.119209493948_wp,-21.444833411574_wp, &
      & -26.674508354162_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.056910263744_wp, 0.059732400683_wp, 0.052197041661_wp, &
      & 0.097844384419_wp, 0.047961348424_wp, 0.057230114090_wp, &
      & 0.060889237549_wp, 0.068600226343_wp, 0.042102755940_wp, &
      & 0.042610278777_wp]
   real(wp), parameter :: ref_hlgaps(10) = &
      &[3.560961680882_wp, 1.701178515210_wp, 2.575575698307_wp, &
      & 1.642444939544_wp, 2.389954524574_wp, 3.624912973388_wp, &
      & 2.770751688023_wp, 1.357191094207_wp, 4.102623151933_wp, &
      & 1.391461153653_wp]

   call init(env)
   do iMol = 1,10

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newXTBCalculator(env, mol, calc, 'param_gfn1-xtb.txt', 1)
      call newWavefunction(env, mol, calc, chk)
      call addSolvationModel(env, calc, TSolvInput(solvent=trim(solvents(iMol)), &
         & cosmo=.true., nang=170))

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)
      call check_(error, hl_gap, ref_hlgaps(iMol), thr=thr)

   end do

end subroutine test_gfn1_mindless_cosmo

end module test_gfn1
