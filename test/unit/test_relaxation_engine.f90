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

module test_relaxation_engine
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check
   implicit none
   private

   public :: collect_relaxation_engine

contains

subroutine collect_relaxation_engine(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gfnff-fix-lancopt", test_gfnff_fix_lancopt) &
      ]

end subroutine collect_relaxation_engine


subroutine test_gfnff_fix_lancopt(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init
   use xtb_type_restart, only : TRestart
   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_relaxation_engine, only : l_ancopt
   use xtb_setparam, only : set, p_ext_gfnff, p_olev_normal
   use xtb_fixparam, only : init_fix, clear_fix, fixset
   use xtb_splitparam, only : init_split, atmass
   use mctcpar_atomic_masses, only : atomic_mass
   use xtb_mctc_convert, only : autoamu

   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 510
   integer, parameter :: nfix = 500
   real(wp), parameter :: spacing = 2.8_wp
   real(wp), parameter :: zigzag = 0.5_wp
   real(wp), parameter :: frag_sep = 200.0_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc

   integer, allocatable :: at(:)
   real(wp), allocatable :: xyz(:,:)
   real(wp), allocatable :: gradient(:, :)
   real(wp) :: energy, egap, sigma(3, 3)

   logical :: fail, exitRun
   integer :: saved_mode_extrun, saved_micro_opt

   integer :: i

   ! Regression test for: GFN-FF + L-ANC + $fix + (N > 500) used to abort in
   ! ANC generation with e.g. "k=30 nvar=1530" / "ANC generation failed".
   ! Core invariant: degrees of freedom must account for fixed atoms,
   !   nvar = 3*N - 3*nfix - 3,
   ! even when using the fragmented Hessian path in L-ANC.

   saved_mode_extrun = set%mode_extrun
   saved_micro_opt = set%optset%micro_opt
   set%mode_extrun = p_ext_gfnff
   set%optset%micro_opt = 1 ! keep runtime low (1 LBFGS micro step)

   allocate(at(nat), source=6)
   allocate(xyz(3, nat), source=0.0_wp)

   do i = 1, nat/2
      xyz(1, i) = spacing*real(i-1, wp)
      xyz(2, i) = merge(zigzag, -zigzag, mod(i, 2) == 0)
   end do
   do i = nat/2 + 1, nat
      xyz(1, i) = spacing*real(i-1-nat/2, wp)
      xyz(2, i) = frag_sep + merge(zigzag, -zigzag, mod(i, 2) == 0)
   end do

   call init(env)
   call init(mol, at, xyz)

   call init_split(mol%n)
   atmass = atomic_mass(mol%at) * autoamu

   call delete_file('charges')
   call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
   call env%check(exitRun)
   call check_(error, .not.exitRun)
   if (allocated(error)) goto 100

   call init_fix(mol%n)
   fixset%n = nfix
   do i = 1, nfix
      fixset%atoms(i) = i
   end do

   allocate(gradient(3, nat), source=0.0_wp)
   energy = 0.0_wp
   egap = 0.0_wp
   sigma = 0.0_wp

   call l_ancopt(env, -1, mol, chk, calc, p_olev_normal, 1, energy, egap, gradient, sigma, 0, fail)
   call check_(error, .not.fail)

100 continue
   call clear_fix
   set%mode_extrun = saved_mode_extrun
   set%optset%micro_opt = saved_micro_opt
   call mol%deallocate

end subroutine test_gfnff_fix_lancopt


end module test_relaxation_engine
