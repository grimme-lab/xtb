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
      new_unittest("heavy", test_gfn1_heavy), &
      new_unittest("gbsa", test_gfn1gbsa_api), &
      new_unittest("pcem", test_gfn1_pcem_api), &
      new_unittest("xb", test_gfn1_xb), &
      new_unittest("pbc3d", test_gfn1_pbc3d), &
      new_unittest("mindless-basic", test_gfn1_mindless_basic), &
      new_unittest("mindless-solvation", test_gfn1_mindless_solvation), &
      new_unittest("ipea-indole", test_ipea_indole), &
      new_unittest("mindless-cosmo", test_gfn1_mindless_cosmo), &
      new_unittest("restart", test_gfn1_restart) &
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
   call check_(error, res%e_es,    0.050115468904_wp, thr=thr)
   call check_(error, res%e_disp, -0.000137581088_wp, thr=thr)
   call check_(error, res%e_rep,   0.036919047759_wp, thr=thr)

   call check_(error, wfn%q(2),0.33430923478162_wp, thr=thr2)
   call check_(error, wfn%q(3),wfn%q(2), thr=thr2)
   call check_(error, wfn%qsh(1),0.30697460552546_wp, thr=thr2)
   call check_(error, wfn%qsh(4),-0.31974031743775E-01_wp, thr=thr2)

   call check_(error, wfn%ihomo,4)
   call check_(error, wfn%ihomoa,wfn%ihomob)
   call check_(error, wfn%emo(wfn%ihomo),-13.612098312207_wp, thr=thr2)
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

   call check_(error, hl_gap, 5.6067618151110_wp, thr=thr)
   call check_(error, energy,-8.4156335928615_wp, thr=thr)
   call check_(error, norm2(gradient),0.95790087998601E-02_wp, thr=thr)

   call check_(error, gradient(1,5), 0.18116554909835E-02_wp, thr=thr)
   call check_(error, gradient(2,2), 0.00000000000000E+00_wp, thr=thr)
   call check_(error, gradient(1,4),-0.18116554909835E-02_wp, thr=thr)
   call check_(error, gradient(3,6),-0.76253144701942E-03_wp, thr=thr)

end subroutine test_gfn1_api

subroutine test_gfn1_heavy(error)
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

   ! Elements with 5sp valence shells (I, Sn) cover the STO-6G expansion
   ! of the 5s and 5p functions, which is not exercised by any other test
   real(wp),parameter :: thr = 1.0e-7_wp
   integer, parameter :: nat1 = 5
   integer, parameter :: at1(nat1) = [6, 53, 1, 1, 1]
   real(wp),parameter :: xyz1(3,nat1) = reshape(&
      &[0.00000000000000_wp, 0.00000000000000_wp, 0.00000000000000_wp, &
      & 0.00000000000000_wp, 0.00000000000000_wp,-4.05000000000000_wp, &
      & 1.94000000000000_wp, 0.00000000000000_wp, 0.69000000000000_wp, &
      &-0.97000000000000_wp, 1.68009920000000_wp, 0.69000000000000_wp, &
      &-0.97000000000000_wp,-1.68009920000000_wp, 0.69000000000000_wp], shape(xyz1))
   integer, parameter :: nat2 = 5
   integer, parameter :: at2(nat2) = [50, 1, 1, 1, 1]
   real(wp),parameter :: xyz2(3,nat2) = reshape(&
      &[0.00000000000000_wp, 0.00000000000000_wp, 0.00000000000000_wp, &
      & 1.85000000000000_wp, 1.85000000000000_wp, 1.85000000000000_wp, &
      &-1.85000000000000_wp,-1.85000000000000_wp, 1.85000000000000_wp, &
      &-1.85000000000000_wp, 1.85000000000000_wp,-1.85000000000000_wp, &
      & 1.85000000000000_wp,-1.85000000000000_wp,-1.85000000000000_wp], shape(xyz2))
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

   ! methyl iodide
   call init(mol, at1, xyz1)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=1)
   call newWavefunction(env, mol, calc, wfn)

   call calc%singlepoint(env, mol, wfn, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call check_(error, hl_gap, 4.5821368916337_wp, thr=thr)
   call check_(error, energy,-7.6689568639736_wp, thr=thr)
   call check_(error, norm2(gradient), 0.10777047947064E-01_wp, thr=thr)

   call check_(error, gradient(3,1),-0.72530697669606E-02_wp, thr=thr)
   call check_(error, gradient(3,2), 0.57697004960942E-02_wp, thr=thr)
   call check_(error, gradient(1,3), 0.31347786926901E-02_wp, thr=thr)

   deallocate(gradient)

   ! stannane
   call init(mol, at2, xyz2)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=1)
   call newWavefunction(env, mol, calc, wfn)

   call calc%singlepoint(env, mol, wfn, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call check_(error, hl_gap, 8.8297303007391_wp, thr=thr)
   call check_(error, energy,-4.6287486160793_wp, thr=thr)
   call check_(error, norm2(gradient), 0.25503852666812E-01_wp, thr=thr)

   call check_(error, gradient(1,2), 0.73623281012783E-02_wp, thr=thr)
   call check_(error, gradient(2,3),-0.73623281012783E-02_wp, thr=thr)

end subroutine test_gfn1_heavy

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

   call check_(error, hl_gap, 6.641469982355_wp, thr=1e-4_wp)
   call check_(error, energy,-14.215790820910_wp, thr=thr)
   call check_(error, norm2(gradient),0.14758421349075E-01_wp, thr=thr)
   call check_(error, gradient(2,3),0.1002882653346e-01_wp, thr=thr)
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

   call check_(error, hl_gap, 9.0155262664573_wp, thr=thr*10)
   call check_(error, energy,-23.113490914998_wp, thr=thr)
   call check_(error, norm2(gradient),0.11143043525454E-01_wp, thr=thr)

   call check_(error, gradient(1,5),-0.17083320570111E-02_wp, thr=thr)
   call check_(error, gradient(2,2), 0.17083320570110E-02_wp, thr=thr)
   call check_(error, gradient(1,4), 0.27992723747300E-02_wp, thr=thr)
   call check_(error, gradient(3,6),-0.10966273400897E-02_wp, thr=thr)

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

   call check_(error, hl_gap, 8.7251771874211_wp, thr=thr)
   call check_(error, energy,-11.559896105984_wp, thr=thr)
   call check_(error, norm2(gradient),0.24211924303314E-01_wp, thr=thr)

   call check_(error, gradient(1,5),-0.16369208859684E-02_wp, thr=thr)
   call check_(error, gradient(2,2), 0.34511615007829E-02_wp, thr=thr)
   call check_(error, gradient(1,4),-0.27846663969147E-02_wp, thr=thr)
   call check_(error, gradient(3,6),-0.12172970039721E-02_wp, thr=thr)

   call check_(error, norm2(calc%pcem%grd),0.12965272381631E-01_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,5),-0.65535083622257E-03_wp, thr=thr)
   call check_(error, calc%pcem%grd(2,2), 0.17820485074561E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,4), 0.60889168360438E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(3,6), 0.86753420850359E-03_wp, thr=thr)

   ! reset
   energy = 0.0_wp
   gradient = 0.0_wp
   calc%pcem%grd = 0.0_wp
   calc%pcem%gam = 999.0_wp ! point charges

   call newWavefunction(env, mol, calc, wfn)
   call calc%singlepoint(env, mol, wfn, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call check_(error, hl_gap, 8.9181165007295_wp, thr=thr)
   call check_(error, energy,-11.565012263827_wp, thr=thr)
   call check_(error, norm2(gradient),0.23134746641164E-01_wp, thr=thr)

   call check_(error, gradient(1,5),-0.63003810023097E-03_wp, thr=thr)
   call check_(error, gradient(2,2), 0.28672774098732E-02_wp, thr=thr)
   call check_(error, gradient(1,4),-0.75489084958916E-02_wp, thr=thr)
   call check_(error, gradient(3,6),-0.12126926086201E-02_wp, thr=thr)

   call check_(error, norm2(calc%pcem%grd),0.18251542582434E-01_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,5),-0.16079978457523E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(2,2), 0.23750272542147E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,4), 0.10140274089826E-01_wp, thr=thr)
   call check_(error, calc%pcem%grd(3,6), 0.11833684312000E-02_wp, thr=thr)

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

   call check_(error, hl_gap, 2.4992555928479_wp, thr=thr)
   call check_(error, energy,-15.606233893431_wp, thr=thr)
   call check_(error, norm2(gradient),0.23014500667889E-01_wp, thr=thr)

   call check_(error, gradient(1,5),-0.39000062479073E-02_wp, thr=thr)
   call check_(error, gradient(2,2),-0.49295777191606E-02_wp, thr=thr)
   call check_(error, gradient(1,4), 0.17228260050308E-01_wp, thr=thr)
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

   call check_(error, hl_gap, 7.5302581099357_wp, thr=thr)
   call check_(error, energy,-11.069452578476_wp, thr=thr)
   call check_(error, norm2(gradient), 0.28766518636475E-01_wp, thr=thr)

   call check_(error, gradient(1,1),-0.18972454825631E-03_wp, thr=thr)
   call check_(error, gradient(1,2),-0.66834799702025E-02_wp, thr=thr)
   call check_(error, gradient(1,3), 0.68732047705008E-02_wp, thr=thr)
   call check_(error, gradient(2,1),-0.26831102888425E-03_wp, thr=thr)
   call check_(error, gradient(3,3), 0.16835844585390E-01_wp, thr=thr)

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
      &[0.049511012828239_wp, 0.062164535741469_wp, 0.050970379218436_wp, &
      & 0.089329119085910_wp, 0.049187948039353_wp, 0.059354443034588_wp, &
      & 0.063931958282716_wp, 0.059153850571371_wp, 0.043308153196369_wp, &
      & 0.043009186658671_wp]
   real(wp), parameter :: ref_hlgaps(10) = &
      &[3.4192922278275_wp, 1.6702720653143_wp, 2.4899353609074_wp, &
      & 1.6706542694313_wp, 2.4613663260650_wp, 3.9778145408256_wp, &
      & 2.5369739726974_wp, 1.3590705098817_wp, 4.0815997449053_wp, &
      & 1.3640393430408_wp]

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
      &[0.058132319775704_wp, 0.061771791928070_wp, 0.056320658647760_wp, &
      & 0.094108321785749_wp, 0.049107732913765_wp, 0.056792270927399_wp, &
      & 0.063360697329010_wp, 0.059227387236917_wp, 0.041889911457433_wp, &
      & 0.045018593389814_wp]
   real(wp), parameter :: ref_hlgaps(10) = &
      &[3.5484491978907_wp, 1.6697112965583_wp, 2.5331239090835_wp, &
      & 1.6610557534739_wp, 2.3957225956189_wp, 3.7354677402279_wp, &
      & 2.6533561727140_wp, 1.3546217990187_wp, 4.0956192322644_wp, &
      & 1.3095350888232_wp]


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
      call check_(error, norm2(gradient), 0.84641977700707E-01_wp, thr=thr)
      call check_(error, hl_gap, 2.5956236874585_wp, thr=thr)

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
      &[-33.086826328982_wp,-26.892421146445_wp,-25.862366002296_wp, &
      & -24.929812238139_wp,-29.068535337284_wp,-20.631251045425_wp, &
      & -35.780608221292_wp,-33.119209598199_wp,-21.444833460634_wp, &
      & -26.674507760379_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.056910443984_wp, 0.059732368118_wp, 0.052196930249_wp, &
      & 0.097844394536_wp, 0.047961212224_wp, 0.057229874514_wp, &
      & 0.060889251332_wp, 0.068599771521_wp, 0.042102585891_wp, &
      & 0.042610069723_wp]
   real(wp), parameter :: ref_hlgaps(10) = &
      &[3.560961563706_wp, 1.701181755113_wp, 2.575588152573_wp, &
      & 1.642442603588_wp, 2.389958024194_wp, 3.624942110407_wp, &
      & 2.770747631357_wp, 1.357190976011_wp, 4.102622313077_wp, &
      & 1.391454491037_wp]

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


subroutine test_gfn1_restart(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_wavefunction, only : TWavefunction
   use xtb_restart, only : readRestart, writeRestart

   type(error_type), allocatable, intent(out) :: error
   character(len=*), parameter :: fname = '.xtb_issue_1146_restart'

   type(TEnvironment) :: env
   type(TWavefunction) :: wfn, wfn_bad
   logical :: success

   call init(env)
   call wfn%allocate(2, 2, 2)
   wfn%nel = 2
   wfn%nopen = 0
   wfn%qsh = [0.1_wp, 0.2_wp]

   call writeRestart(env, wfn, fname, 1)

   call readRestart(env, wfn, fname, 2, [1, 1], 1, success, .false.)
   call check_(error, success)

   call readRestart(env, wfn, fname, 2, [1, 1], 2, success, .false.)
   call check_(error, success)

   wfn%nel = 3
   call readRestart(env, wfn, fname, 2, [1, 1], 1, success, .false.)
   call check_(error, success)

   wfn%nel = 2
   wfn%nopen = 1
   call readRestart(env, wfn, fname, 2, [1, 1], 1, success, .false.)
   call check_(error, success)

   call wfn_bad%allocate(3, 3, 3)
   wfn_bad%nel = 2
   wfn_bad%nopen = 0
   call readRestart(env, wfn_bad, fname, 3, [1, 1, 1], 1, success, .false.)
   call check_(error, .not.success)

   call delete_file(fname)
   call wfn_bad%deallocate
   call wfn%deallocate

end subroutine test_gfn1_restart

end module test_gfn1
