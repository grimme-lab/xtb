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

module test_gfn2
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check, test_failed
   implicit none
   private

   public :: collect_gfn2

contains

!> Collect all exported unit tests
subroutine collect_gfn2(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("scc", test_gfn2_scc), &
      new_unittest("api", test_gfn2_api), &
      new_unittest("rotation", test_gfn2_rotation), &
      new_unittest("translation", test_gfn2_translation), &
      new_unittest("gbsa", test_gfn2gbsa_api), &
      new_unittest("salt", test_gfn2salt_api), &
      new_unittest("pcem", test_gfn2_pcem_api), &
      new_unittest("pcem_io", test_gfn2_pcem_io), &
      new_unittest("mindless-basic", test_gfn2_mindless_basic), &
      new_unittest("mindless-solvation", test_gfn2_mindless_solvation), &
      new_unittest("dmetal", test_gfn2_dmetal), &
      new_unittest("mindless-cosmo", test_gfn2_mindless_cosmo), &
      new_unittest("wbo", test_gfn2_wbo), &
      new_unittest("lmo-molden", test_gfn2_lmo_molden) &
      ]

end subroutine collect_gfn2


subroutine test_gfn2_scc(error)
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
   use xtb_xtb_gfn2

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
   logical  :: okpar,okbas
   logical :: exitRun

   set%gfn_method = 2
   call init(env)

   call init(mol, at, xyz)

   wfn%nel = idint(sum(mol%z))
   wfn%nopen = 0

   allocate( g(3,mol%n), source = 0.0_wp )

   call use_parameterset('param_gfn2-xtb.txt',globpar,xtbData,okpar)
   call check_(error, okpar)

   call newBasisset(xtbData,mol%n,mol%at,basis,okbas)
   call check_(error, okbas)

   call check_(error, basis%nshell,4)
   call check_(error, basis%nbf,   6)
   call check_(error, basis%nao,   6)

   call wfn%allocate(mol%n,basis%nshell,basis%nao)
   wfn%q = mol%chrg/real(mol%n,wp)

   call iniqshell(xtbData,mol%n,mol%at,mol%z,basis%nshell,wfn%q,wfn%qsh,set%gfn_method)

   g = 0.0_wp

   call scf(env,mol,wfn,basis,pcem,xtbData,solvation, &
      &   egap,et,maxiter,prlevel,restart,lgrad,acc,etot,g,res)

   call env%check(exitRun)
   call check_(error, .not.exitRun)

   call check_(error, res%converged)

   call check_(error, res%e_total,-5.070451355118_wp, thr=thr)
   call check_(error, res%gnorm,   0.006457865418_wp, thr=thr)
   ! value in electron volt
   call check_(error, res%hl_gap, 14.450372368833_wp, thr=1.0e-4_wp)
   call check_(error, res%e_elec, -5.104362813671_wp, thr=thr)
   ! es and aes are not welldefined at this accuracy level
   call check_(error, res%e_es,    0.031408028914_wp, thr=thr*10)
   call check_(error, res%e_aes,   0.000563300591_wp, thr=thr*10)
   call check_(error, res%e_axc,  -0.000813611055_wp, thr=thr)
   call check_(error, res%e_disp, -0.000141250480_wp, thr=thr)
   call check_(error, res%e_rep,   0.033911458523_wp, thr=thr)

   call check_(error, wfn%q(2),0.28158903353422_wp, thr=thr2)
   call check_(error, wfn%q(3),wfn%q(2), thr=thr2)
   call check_(error, wfn%qsh(1),0.25886041477578_wp, thr=thr2)
   call check_(error, wfn%dipm(2,2), 0.00000000000000E+00_wp, thr=thr2)
   call check_(error, wfn%dipm(1,3),-0.55238805565739E-01_wp, thr=thr2)
   call check_(error, wfn%qp(6,1),0.41207363331933E-01_wp, thr=thr2)
   call check_(error, wfn%qp(1,3),0.95084219046109E-01_wp, thr=thr2)

   call check_(error, wfn%ihomo,4)
   call check_(error, wfn%ihomoa,wfn%ihomob)
   call check_(error, wfn%emo(wfn%ihomo),-12.166283951806_wp, thr=thr2)
   call check_(error, wfn%focca(wfn%ihomo),wfn%foccb(wfn%ihomo), thr=thr2)
   call check_(error, wfn%focc(wfn%ihomo),2.0_wp, thr=thr2)

   call mol%deallocate
   call wfn%deallocate
   call basis%deallocate

end subroutine test_gfn2_scc

subroutine test_gfn2_api(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   type(error_type), allocatable, intent(out) :: error

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
   call check_(error, hl_gap, 7.0005867673492_wp, thr=thr)
   call check_(error, energy,-8.3824793849585_wp, thr=thr)
   call check_(error, norm2(gradient),0.11544364742340E-01_wp, thr=thr)

   call check_(error, gradient(2,3), 0.00000000000000E+00_wp, thr=thr)
   call check_(error, gradient(3,1),-0.74644788414501E-03_wp, thr=thr)
   call check_(error, gradient(1,4),-0.28436357035188E-03_wp, thr=thr)
   call check_(error, gradient(3,7), 0.99750262249210E-02_wp, thr=thr)

end subroutine test_gfn2_api


!> The SCC mixing vector has to transform orthogonally under rotations of the
!> molecule, otherwise the Broyden iterations depend on the orientation and
!> the SCC converges to different charges for rotated copies of one molecule.
!> This requires the off-diagonal quadrupole components to be weighted by
!> sqrt(2), which is checked here with default SCC accuracy.
subroutine test_gfn2_rotation(error)
   use xtb_mctc_accuracy, only : wp

   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_data
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr_energy = 1.0e-10_wp
   real(wp), parameter :: thr_charge = 1.0e-8_wp
   real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
   integer, parameter :: nat = 5
   integer, parameter :: at(nat) = [6, 8, 8, 1, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp, 0.75589045018313_wp, 0.00000000000000_wp, &
      & 2.17318504427650_wp, 1.54957542287542_wp, 0.18897261254578_wp, &
      &-1.98421243173072_wp, 2.26767135054939_wp,-0.18897261254578_wp, &
      & 0.05669178376373_wp,-1.32280828782048_wp, 0.09448630627289_wp, &
      &-3.59047963836987_wp, 1.51178090036626_wp, 0.00000000000000_wp], &
      & shape(xyz))
   !> Rotation axis, generic to produce all components of the quadrupole
   real(wp), parameter :: axis(3) = [1.0_wp, 2.0_wp, 3.0_wp] / sqrt(14.0_wp)
   real(wp), parameter :: angle = 0.25_wp * pi

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   integer :: irot, i
   real(wp) :: rot(3, 3), cross(3, 3), xyz_rot(3, nat)
   real(wp) :: energy(2), charges(nat, 2)
   real(wp) :: sigma(3, 3), hl_gap
   real(wp), allocatable :: gradient(:, :)

   call init(env)

   ! Rotation matrix from the Rodrigues formula
   cross = reshape([ &
      &  0.0_wp,    axis(3), -axis(2), &
      & -axis(3),   0.0_wp,   axis(1), &
      &  axis(2),  -axis(1),  0.0_wp], [3, 3])
   rot = (1.0_wp - cos(angle)) * spread(axis, 2, 3) * spread(axis, 1, 3) &
      & + sin(angle) * cross
   do i = 1, 3
      rot(i, i) = rot(i, i) + cos(angle)
   end do
   xyz_rot = matmul(rot, xyz)

   allocate(gradient(3, nat))

   do irot = 1, 2
      if (irot == 1) then
         call init(mol, at, xyz)
      else
         call init(mol, at, xyz_rot)
      end if
      energy(irot) = 0.0_wp
      gradient = 0.0_wp

      call newXTBCalculator(env, mol, calc, method=2)
      call newWavefunction(env, mol, calc, chk)

      call calc%singlepoint(env, mol, chk, 0, .false., energy(irot), gradient, &
         & sigma, hl_gap, res)
      call check_(error, res%converged)
      if (allocated(error)) return

      charges(:, irot) = chk%wfn%q
      call mol%deallocate
   end do

   call check_(error, maxval(abs(charges(:, 2) - charges(:, 1))), 0.0_wp, &
      & thr=thr_charge)
   if (allocated(error)) return

   call check_(error, energy(2), energy(1), thr=thr_energy)

end subroutine test_gfn2_rotation

subroutine test_gfn2_translation(error)
   use xtb_mctc_accuracy, only : wp

   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_data
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   type(error_type), allocatable, intent(out) :: error

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
   real(wp),parameter :: shift(3) = [50.0_wp, -30.0_wp, 20.0_wp]

   type(TMolecule) :: mol, moll
   type(TRestart) :: chk, chkl
   type(TEnvironment) :: env
   type(scc_results) :: res, resl
   type(TxTBCalculator) :: calc

   integer :: iat
   real(wp) :: energy, energyl, sigma(3, 3), hl_gap
   real(wp), allocatable :: gradient(:, :), gradientl(:, :), xyzl(:, :)

   call init(env)
   allocate(gradient(3, nat), gradientl(3, nat), xyzl(3, nat), source=0.0_wp)

   call init(mol, at, xyz)
   call newXTBCalculator(env, mol, calc, method=2)
   call newWavefunction(env, mol, calc, chk)
   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   ! the multipole integrals and potentials must not depend on the molecular position
   do iat = 1, nat
      xyzl(:, iat) = xyz(:, iat) + shift
   end do
   call init(moll, at, xyzl)
   call newWavefunction(env, moll, calc, chkl)
   call calc%singlepoint(env, moll, chkl, 2, .false., energyl, gradientl, sigma, &
      & hl_gap, resl)

   call check_(error, energyl, energy, thr=thr)
   call check_(error, maxval(abs(gradientl - gradient)), 0.0_wp, thr=thr)
   call check_(error, maxval(abs(resl%dipole - res%dipole)), 0.0_wp, thr=thr)

end subroutine test_gfn2_translation

subroutine test_gfn2gbsa_api(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment
   use xtb_solv_input

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_kernel, only : gbKernel

   type(error_type), allocatable, intent(out) :: error

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
   call addSolvationModel(env, calc, TSolvInput(solvent=opt%solvent, alpb=.false., kernel=gbKernel%still))

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call check_(error, hl_gap, 3.408569859125_wp, thr=1e-5_wp)
   call check_(error, energy,-22.002501380096_wp, thr=thr)
   call check_(error, norm2(gradient),0.19367559715184E-01_wp, thr=thr)
   call check_(error, gradient(1,1), .9038674439127e-02_wp, thr=thr)
   call check_(error, gradient(3,2),-.1394693523214e-02_wp, thr=thr)
   call check_(error, gradient(3,11),-gradient(3,10), thr=thr)
   call check_(error, gradient(1,8),0.17878350605942E-02_wp, thr=thr)

end subroutine test_gfn2gbsa_api

subroutine test_gfn2salt_api(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_convert, only : aatoau

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment
   use xtb_solv_input

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_kernel, only : gbKernel

   type(error_type), allocatable, intent(out) :: error

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
      & ionRad=1.0_wp*aatoau, ionStrength=1.0e-3_wp, alpb=.false., kernel=gbKernel%still))

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call check_(error, hl_gap, 6.895830675032_wp, thr=5e-5_wp)
   call check_(error, energy,-13.027106170796_wp, thr=thr)
   call check_(error, norm2(gradient),0.014657211478_wp, thr=thr)

   call check_(error, gradient(1,1),-0.6694856649275e-02_wp, thr=thr)
   call check_(error, gradient(2,4),-0.7608884457863e-03_wp, thr=thr)
   call check_(error, gradient(1,5), 0.2358927421507e-02_wp, thr=thr)
   call check_(error, gradient(3,7), 0.6883271590152e-02_wp, thr=thr)

end subroutine test_gfn2salt_api

subroutine test_gfn2_pcem_api(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_pcem
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction
   use xtb_main_setup, only : addSolvationModel

   type(error_type), allocatable, intent(out) :: error

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

   call check_(error, hl_gap, 12.391067835855_wp, thr=thr)
   call check_(error, energy,-20.323978512016_wp, thr=thr)
   call check_(error, norm2(gradient),0.78076070662211E-02_wp, thr=thr)

   call check_(error, gradient(1,5),-0.22189420554446E-02_wp, thr=thr)
   call check_(error, gradient(2,2), 0.22189420554445E-02_wp, thr=thr)
   call check_(error, gradient(1,4), 0.95540323278656E-03_wp, thr=thr)
   call check_(error, gradient(3,6),-0.11891232885042E-02_wp, thr=thr)

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

   call check_(error, hl_gap, 12.718455231981_wp, thr=thr)
   call check_(error, energy,-10.160927754776_wp, thr=thr)
   call check_(error, norm2(gradient),0.21551799327098E-01_wp, thr=thr)

   call check_(error, gradient(1,5),-0.20334032406558E-02_wp, thr=thr)
   call check_(error, gradient(2,2), 0.33126833827722E-02_wp, thr=thr)
   call check_(error, gradient(1,4),-0.11907187586151E-02_wp, thr=thr)
   call check_(error, gradient(3,6),-0.16630879952139E-02_wp, thr=thr)

   call check_(error, norm2(calc%pcem%grd),0.10043973324176E-01_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,5),-0.24831204136192E-03_wp, thr=thr)
   call check_(error, calc%pcem%grd(2,2), 0.14208406102449E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,4), 0.37466631500699E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(3,6), 0.65161176362747E-03_wp, thr=thr)

   ! reset
   energy = 0.0_wp
   gradient = 0.0_wp
   calc%pcem%grd = 0.0_wp
   calc%pcem%gam = 999.0_wp ! point charges

   call newWavefunction(env, mol, calc, chk)
   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call check_(error, hl_gap, 13.024343620925_wp, thr=thr)
   call check_(error, energy,-10.168788268962_wp, thr=thr)
   call check_(error, norm2(gradient),0.18113593320713E-01_wp, thr=thr)

   call check_(error, gradient(1,5),-0.50646432321609E-03_wp, thr=thr)
   call check_(error, gradient(2,2), 0.25468331741091E-02_wp, thr=thr)
   call check_(error, gradient(1,4),-0.74927895629019E-02_wp, thr=thr)
   call check_(error, gradient(3,6),-0.13248434928006E-02_wp, thr=thr)

   call check_(error, norm2(calc%pcem%grd),0.18721775625220E-01_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,5),-0.21573659448600E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(2,2), 0.25653642490271E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,4), 0.12124325222776E-01_wp, thr=thr)
   call check_(error, calc%pcem%grd(3,6), 0.12140563197835E-02_wp, thr=thr)

end subroutine test_gfn2_pcem_api

subroutine test_gfn2_pcem_io(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_pcem
   use xtb_type_environment
   use xtb_embedding, only: read_pcem

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction
   use xtb_main_setup, only : addSolvationModel

   type(error_type), allocatable, intent(out) :: error

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

   integer :: tmp_unit

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

   call check_(error, hl_gap, 12.391067835855_wp, thr=thr)
   call check_(error, energy,-20.323978512016_wp, thr=thr)
   call check_(error, norm2(gradient),0.78076070662222E-02_wp, thr=thr)

   call check_(error, gradient(1,5),-0.22189420554448E-02_wp, thr=thr)
   call check_(error, gradient(2,2), 0.22189420554443E-02_wp, thr=thr)
   call check_(error, gradient(1,4), 0.95540323278606E-03_wp, thr=thr)
   call check_(error, gradient(3,6),-0.11891232885043E-02_wp, thr=thr)

   ! reset
   call mol%deallocate
   deallocate(gradient)

   call init(mol, at(:nat2), xyz(:, :nat2))
   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=2)

   open(newunit=tmp_unit, Status="Scratch")
      write(tmp_unit,'(a)') &
      "6",&
      "-0.69645733    2.75237178376284    -2.43247309226225    -0.01392519847964  O", &
      " 0.36031084    0.93157260886974    -2.79621404458590    -0.01863384029005  H", &
      " 0.33614649    3.43820531288547    -3.30583608421060     1.42134539425148  H", &
      "-0.69645733    2.43247309226225     2.75237178376284     0.01392519847964  O", &
      " 0.36031084    2.79621404458590     0.93157260886974     0.01863384029005  H", &
      " 0.33614649    3.30583608421060     3.43820531288547    -1.42134539425148  H"
   Call read_pcem(tmp_unit,env,calc%pcem,calc%xtbData%coulomb)
   close(tmp_unit,Status="Delete")

   calc%pcem%grd = 0.0_wp

   call newWavefunction(env, mol, calc, chk)
   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call check_(error, hl_gap, 12.718455231981_wp, thr=thr)
   call check_(error, energy,-10.160927754776_wp, thr=thr)
   call check_(error, norm2(gradient),0.21551799327098E-01_wp, thr=thr)

   call check_(error, gradient(1,5),-0.20334032406558E-02_wp, thr=thr)
   call check_(error, gradient(2,2), 0.33126833827722E-02_wp, thr=thr)
   call check_(error, gradient(1,4),-0.11907187586151E-02_wp, thr=thr)
   call check_(error, gradient(3,6),-0.16630879952139E-02_wp, thr=thr)

   call check_(error, norm2(calc%pcem%grd),0.10043973324176E-01_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,5),-0.24831204136192E-03_wp, thr=thr)
   call check_(error, calc%pcem%grd(2,2), 0.14208406102449E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,4), 0.37466631500699E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(3,6), 0.65161176362747E-03_wp, thr=thr)

   ! reset
   energy = 0.0_wp
   gradient = 0.0_wp
   calc%pcem%grd = 0.0_wp
   calc%pcem%gam = 999.0_wp ! point charges

   call newWavefunction(env, mol, calc, chk)
   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call check_(error, hl_gap, 13.024343620925_wp, thr=thr)
   call check_(error, energy,-10.168788268962_wp, thr=thr)
   call check_(error, norm2(gradient),0.18113593320713E-01_wp, thr=thr)

   call check_(error, gradient(1,5),-0.50646432321609E-03_wp, thr=thr)
   call check_(error, gradient(2,2), 0.25468331741091E-02_wp, thr=thr)
   call check_(error, gradient(1,4),-0.74927895629019E-02_wp, thr=thr)
   call check_(error, gradient(3,6),-0.13248434928006E-02_wp, thr=thr)

   call check_(error, norm2(calc%pcem%grd),0.18721775625220E-01_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,5),-0.21573659448600E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(2,2), 0.25653642490271E-02_wp, thr=thr)
   call check_(error, calc%pcem%grd(1,4), 0.12124325222776E-01_wp, thr=thr)
   call check_(error, calc%pcem%grd(3,6), 0.12140563197835E-02_wp, thr=thr)

end subroutine test_gfn2_pcem_io

subroutine test_gfn2_mindless_basic(error)
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
      &[-30.348902328339_wp, -24.069929673504_wp, -23.718708414840_wp, &
      & -22.751497890244_wp, -27.735987617741_wp, -18.559531652263_wp, &
      & -33.423246445911_wp, -29.982410709899_wp, -20.549458882257_wp, &
      & -25.648093031152_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.06081860391196_wp, 0.05480374997753_wp, 0.04180003046237_wp, &
      & 0.06932469610481_wp, 0.05133165781937_wp, 0.05897028984233_wp, &
      & 0.04600837003153_wp, 0.06088707462221_wp, 0.05314386793712_wp, &
      & 0.04588827433979_wp]
   real(wp), parameter :: ref_hlgaps(10) = &
      &[2.2611050499278_wp, 1.2455353260482_wp, 0.1489849070154_wp, &
      & 1.2694337641912_wp, 2.2361895850260_wp, 2.6052851670765_wp, &
      & 2.6312109717862_wp, 0.8502447036074_wp, 3.2960764806338_wp, &
      & 0.4185678259178_wp]

   call init(env)
   do iMol = 1, 10

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newXTBCalculator(env, mol, calc, 'param_gfn2-xtb.txt', 2)
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

end subroutine test_gfn2_mindless_basic


subroutine test_gfn2_mindless_solvation(error)
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
      &[-30.384994095620_wp, -24.092790321708_wp, -23.760000587724_wp, &
      & -22.759033676135_wp, -27.755756860587_wp, -18.578906535415_wp, &
      & -33.448544467530_wp, -29.993777657049_wp, -20.567481283491_wp, &
      & -25.665458003353_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.064382549110824_wp, 0.055824640225377_wp, 0.040696600402333_wp, &
      & 0.073196929823830_wp, 0.049853680430954_wp, 0.053051001578728_wp, &
      & 0.046461139146085_wp, 0.057937529580509_wp, 0.055081761535570_wp, &
      & 0.045339325291251_wp]
   real(wp), parameter :: ref_hlgaps(10) = &
      &[2.5109870734258_wp, 1.3559660535296_wp, 0.1719324902715_wp, &
      & 1.2428846924083_wp, 2.1427985709246_wp, 2.2547542623367_wp, &
      & 2.6934479162794_wp, 0.83095433836574_wp, 3.3270810805007_wp, &
      & 0.3858722147956_wp]

   call init(env)
   do iMol = 1, 10

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newXTBCalculator(env, mol, calc, 'param_gfn2-xtb.txt', 2)
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

end subroutine test_gfn2_mindless_solvation


subroutine test_gfn2_dmetal(error)
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

   real(wp), parameter :: ref_energies(3) = &
      &[-34.066853522474_wp, -34.078580189077_wp, -34.078529449408_wp]
   real(wp), parameter :: ref_gnorms(3) = &
      &[0.26985508088903_wp, 0.26593471388897_wp, 0.26615193356412_wp]
   real(wp), parameter :: ref_hlgaps(3) = &
      &[ 2.9239302552159_wp,  2.9165099394194_wp,  2.9169886443340_wp]

   call init(env)
   do iMol = 1, 3

      call getMolecule(mol, 'feco5')

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newXTBCalculator(env, mol, calc, 'param_gfn2-xtb.txt', 2)
      call newWavefunction(env, mol, calc, chk)
      if (iMol > 1) then
         call addSolvationModel(env, calc, TSolvInput(solvent='ch2cl2', &
            & alpb=iMol==3, kernel=gbKernel%still))
      end if

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 1, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)
      call check_(error, hl_gap, ref_hlgaps(iMol), thr=thr)

   end do

end subroutine test_gfn2_dmetal


subroutine test_gfn2_mindless_cosmo(error)
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
      &[-30.398257705621_wp,-24.125868388467_wp,-23.787176871118_wp, &
      & -22.779539008049_wp,-27.783054816813_wp,-18.610486787528_wp, &
      & -33.456895964521_wp,-30.008556295313_wp,-20.588957661926_wp, &
      & -25.684671213225_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.070118274688_wp, 0.060272025636_wp, 0.048063907692_wp, &
      & 0.076295318322_wp, 0.050041390331_wp, 0.051542844593_wp, &
      & 0.046067454533_wp, 0.060435250964_wp, 0.058023794880_wp, &
      & 0.050440859372_wp]
   real(wp), parameter :: ref_hlgaps(10) = &
      &[3.173017030723_wp, 1.553005205039_wp, 0.324800762527_wp, &
      & 1.233104389193_wp, 1.918285342697_wp, 2.067687622302_wp, &
      & 2.757568550511_wp, 0.846322891292_wp, 3.336532492782_wp, &
      & 0.315091448813_wp]

   call init(env)
   do iMol = 1, 10

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newXTBCalculator(env, mol, calc, 'param_gfn2-xtb.txt', 2)
      call newWavefunction(env, mol, calc, chk)
      call addSolvationModel(env, calc, TSolvInput(solvent=trim(solvents(iMol)), &
         & cosmo=.true., nang=74))

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
end subroutine test_gfn2_mindless_cosmo

subroutine test_gfn2_wbo(error)
   use xtb_mctc_accuracy, only : wp
   use mctc_io_convert, only : aatoau
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

   energy=0.0_wp
   call init(env)
   
   !-----!
   ! RHF !
   !-----!
   
   ! initialization of mol str !
   call getMolecule(mol, 'co_cnx6')
   mol%xyz = mol%xyz*aatoau

   ! allocate GFN2 calculator !
   if (allocated(gradient)) deallocate(gradient)
   allocate(gradient(3, len(mol)))
   call newXTBCalculator(env, mol, calc, 'param_gfn2-xtb.txt', 2)
   call newWavefunction(env, mol, calc, chk)

   call env%check(exitRun)
   call check_(error, .not.exitRun)

   ! actual SP !
   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)
   
   call env%check(exitRun)
   call check_(error, .not.exitRun)

   ! check  scc !
   call check_(error, energy, -34.191439852632_wp, thr=thr)
   call check_(error, norm2(gradient), 0.041767312316_wp, thr=thr)
   call check_(error, hl_gap, 3.769053698090_wp, thr=thr)

   ! check wbo !
   call check_(error, chk%wfn%wbo(1,2), 0.889766823730745_wp,thr=thr) 
   call check_(error, chk%wfn%wbo(3,4), 1.728709690632781E-002_wp,thr=thr) 
   call check_(error, chk%wfn%wbo(2,10), 2.73379950999972_wp,thr=thr) 

   !-----!
   ! UHF !
   !-----!
   
   ! initialization of mol str !
   call getMolecule(mol, 'fe_cnx6')
   mol%xyz = mol%xyz*aatoau
   
   ! allocate GFN2 calculator !
   if (allocated(gradient)) deallocate(gradient)
   allocate(gradient(3, len(mol)))
   call newXTBCalculator(env, mol, calc, 'param_gfn2-xtb.txt', 2)
   call newWavefunction(env, mol, calc, chk)

   call env%check(exitRun)
   call check_(error, .not.exitRun)
  
   ! actual SP !
   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)
 
   call env%check(exitRun)
   call check_(error, .not.exitRun)

   ! check  scc !
   call check_(error, energy, -33.314958144107_wp, thr=thr)
   call check_(error, norm2(gradient), 0.018945278169_wp, thr=thr)
   call check_(error, hl_gap, 1.805933456291_wp, thr=thr)


   ! check wbo !
   call check_(error, chk%wfn%wbo(1,2), 0.587096200577129_wp,thr=thr) 
   call check_(error, chk%wfn%wbo(3,4), 3.798528974276730E-002_wp,thr=thr) 
   call check_(error, chk%wfn%wbo(2,10), 2.81892815139406_wp,thr=thr) 

end subroutine test_gfn2_wbo

!> Check LMO center file and optional the molden output
subroutine test_gfn2_lmo_molden(error)
   use xtb_mctc_accuracy, only : wp
   use mctc_io_convert, only : aatoau
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart
   use xtb_setparam, only : set
   use xtb_solv_model, only : TSolvModel

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction
   use xtb_propertyoutput, only : main_property

   type(error_type), allocatable, intent(out) :: error

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   type(scc_results) :: res
   type(TSolvModel), allocatable :: solvModel

   real(wp) :: energy, hl_gap, sigma(3,3)
   real(wp), allocatable :: gradient(:,:)
   integer :: unit1, unit2, stat1, stat2
   character(len=1024) :: line1, line2
   logical :: exitRun
   logical :: exist_lmoinfo, exist_lmocent, exist_molden, exist_molden_lmo
   logical :: pr_lmo_save, pr_molden_save
   logical :: differ
   logical :: in_mo1, in_mo2
   
   energy = 0.0_wp
   call init(env)
   call getMolecule(mol, 'co_cnx6')
   mol%xyz = mol%xyz*aatoau
   allocate(gradient(3, mol%n))

   call newXTBCalculator(env, mol, calc, 'param_gfn2-xtb.txt', 2)
   call newWavefunction(env, mol, calc, chk)

   call env%check(exitRun)
   call check_(error, .not.exitRun)
   if (allocated(error)) return

   ! Mirror the command line with molden and lmo options
   pr_lmo_save = set%pr_lmo
   pr_molden_save = set%pr_molden_input
   set%pr_lmo = .true.
   set%pr_molden_input = .true.

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call env%check(exitRun)
   call check_(error, .not.exitRun)
   if (allocated(error)) return

   ! The canonical Molden file is written by the property printout
   call main_property(6, env, mol, chk%wfn, calc%basis, calc%xtbData, res, &
      & solvModel, 1.0_wp)

   set%pr_lmo = pr_lmo_save
   set%pr_molden_input = pr_molden_save

   ! Check LMO classification results
   call check_(error, allocated(res%iff_results))
   if (allocated(error)) return
   call check_(error, res%iff_results%nlmo, 33)
   if (allocated(error)) return

   ! Check LMO center files
   inquire(file='xtblmoinfo', exist=exist_lmoinfo)
   inquire(file='lmocent.coord', exist=exist_lmocent)
   call check_(error, exist_lmoinfo)
   if (allocated(error)) return
   call check_(error, exist_lmocent)
   if (allocated(error)) return

   ! Check that both canonical and localized-MO Molden files are written
   inquire(file='molden.input', exist=exist_molden)
   inquire(file='molden-lmo.input', exist=exist_molden_lmo)
   call check_(error, exist_molden)
   if (allocated(error)) return
   call check_(error, exist_molden_lmo)
   if (allocated(error)) return

   ! Check that localized and canonical coefficients are different
   differ = .false.
   in_mo1 = .false.
   in_mo2 = .false.
   open(newunit=unit1, file='molden.input', status='old', action='read')
   open(newunit=unit2, file='molden-lmo.input', status='old', action='read')
   do
      ! Next coefficient line of the canonical file
      do
         read(unit1, '(a)', iostat=stat1) line1
         if (stat1 /= 0) exit
         if (index(line1, '[MO]') == 1) then
            in_mo1 = .true.
            cycle
         end if
         if (in_mo1 .and. .not.any(line1(1:4) == ['Sym=', 'Ene=', 'Spin', 'Occu'])) exit
      end do
      ! Next coefficient line of the localized file
      do
         read(unit2, '(a)', iostat=stat2) line2
         if (stat2 /= 0) exit
         if (index(line2, '[MO]') == 1) then
            in_mo2 = .true.
            cycle
         end if
         if (in_mo2 .and. .not.any(line2(1:4) == ['Sym=', 'Ene=', 'Spin', 'Occu'])) exit
      end do
      if (stat1 /= 0 .or. stat2 /= 0) then
         if (stat1 /= stat2) differ = .true.
         exit
      end if
      if (line1 /= line2) then
         differ = .true.
         exit
      end if
   end do
   close(unit1)
   close(unit2)
   call check_(error, differ)
   if (allocated(error)) return

   call delete_file('xtblmoinfo')
   call delete_file('lmocent.coord')
   call delete_file('molden.input')
   call delete_file('molden-lmo.input')
   call delete_file('charges')

end subroutine test_gfn2_lmo_molden

end module test_gfn2
