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

module test_tblite
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check, &
      & test_failed, skip_test
   use xtb_mctc_accuracy, only : wp
   use xtb_features, only : get_xtb_feature
   use xtb_tblite_calculator, only : TTBLiteCalculator, TTBLiteInput, newTBLiteCalculator, &
      & newTBLiteWavefunction
   implicit none
   private

   public :: collect_tblite

   real(wp),parameter :: thr = 5.0e-8_wp, thrg = 1.0e-6_wp

contains

!> Collect all exported unit tests
subroutine collect_tblite(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gfn1", test_gfn1), &
      new_unittest("gfn2", test_gfn2), &
      new_unittest("gfn1-mindless", test_gfn2_mindless_basic), &
      new_unittest("gfn2-mindless", test_gfn2_mindless_basic) &
      ]

end subroutine collect_tblite

subroutine test_gfn2(error)
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment

   type(error_type), allocatable, intent(out) :: error

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
   type(TTBLiteCalculator) :: calc

   real(wp) :: energy, sigma(3, 3)
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   if (.not.get_xtb_feature('tblite')) then
      call skip_test(error, "xtb not compiled with tblite support")
      return
   end if

   ! setup the environment variables
   call init(env)

   call init(mol, at, xyz)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newTBLiteCalculator(env, mol, calc, TTBLiteInput(method="gfn2"))
   call newTBLiteWavefunction(env, mol, calc, chk)

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)
   call check_(error, energy,-8.3824793849585_wp, thr=thr)
   call check_(error, norm2(gradient),0.11544410028854E-01_wp, thr=thrg)
end subroutine test_gfn2


subroutine test_gfn1(error)

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_data
   use xtb_type_restart
   use xtb_type_environment

   type(error_type), allocatable, intent(out) :: error

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
   type(TRestart) :: chk
   type(scc_results) :: res
   type(TTBLiteCalculator) :: calc

   real(wp) :: energy, sigma(3, 3)
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   if (.not.get_xtb_feature('tblite')) then
      call skip_test(error, "xtb not compiled with tblite support")
      return
   end if

   ! setup the environment variables
   call init(env)

   call init(mol, at, xyz)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newTBLiteCalculator(env, mol, calc, TTBLiteInput(method="gfn1"))
   call newTBLiteWavefunction(env, mol, calc, chk)

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call check_(error, energy,-8.4156335928615_wp, thr=thr)
   call check_(error, norm2(gradient),0.95790240549503E-02_wp, thr=thrg)

end subroutine test_gfn1


subroutine test_gfn2_mindless_basic(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   type(error_type), allocatable, intent(out) :: error

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TTBLiteCalculator) :: calc
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

   if (.not.get_xtb_feature('tblite')) then
      call skip_test(error, "xtb not compiled with tblite support")
      return
   end if

   call init(env)
   do iMol = 2, 10, 2

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newTBLiteCalculator(env, mol, calc, TTBLiteInput(method="gfn2", accuracy=0.01_wp))
      call newTBLiteWavefunction(env, mol, calc, chk)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thrg)

   end do

end subroutine test_gfn2_mindless_basic


subroutine test_gfn1_mindless_basic(error)
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   type(error_type), allocatable, intent(out) :: error

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TTBLiteCalculator) :: calc
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

   if (.not.get_xtb_feature('tblite')) then
      call skip_test(error, "xtb not compiled with tblite support")
      return
   end if

   call init(env)
   do iMol = 1, 10

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newTBLiteCalculator(env, mol, calc, TTBLiteInput(method="gfn1", accuracy=0.01_wp))
      call newTBLiteWavefunction(env, mol, calc, chk)

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

end module test_tblite
