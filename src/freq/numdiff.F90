! This file is part of xtb.
!
! Copyright (C) 2021 Sebastian Ehlert
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

!> Implementation of numerical second derivatives for single point calculators
module xtb_freq_numdiff
   use xtb_mctc_accuracy, only : wp
   use xtb_type_calculator, only : TCalculator
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   implicit none
   private

   public :: numdiff2


   !> Evaluate hessian by finite difference
   interface numdiff2
      module procedure :: numdiff2_all
      module procedure :: numdiff2_list
   end interface numdiff2


contains


!> Evaluate hessian by finite difference for all atoms
subroutine numdiff2_all(env, mol0, chk0, calc, step, hessian, dipgrad, parallelize)
   character(len=*), parameter :: source = "hessian_numdiff_numdiff2"
   !> Computation environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol0
   !> Restart data
   type(TRestart), intent(in) :: chk0
   !> Single point calculator
   class(TCalculator), intent(inout) :: calc
   !> Step size for numerical differentiation
   real(wp), intent(in) :: step
   !> Array to add Hessian to
   real(wp), intent(inout) :: hessian(:, :)
   !> Array to add dipole gradient to
   real(wp), intent(inout) :: dipgrad(:, :)
   !> Parallelize over displacements
   logical, intent(in) :: parallelize

   integer :: iat, jat, ic, jc, ii, jj
   type(TMolecule), allocatable :: mol
   type(TRestart), allocatable :: chk
   real(wp) :: er, el, dr(3), dl(3), sr(3, 3), sl(3, 3), egap, step2
   real(wp) :: t0, t1, w0, w1
   real(wp), allocatable :: gr(:, :), gl(:, :)
   type(scc_results) :: rr, rl

   call timing(t0, w0)
   step2 = 0.5_wp / step
   allocate(gr(3, mol0%n), gl(3, mol0%n))

! This OpenMP statements lead to invalid LLVM-IR with NVHPC (20.7 to 20.11)
#ifndef __PGIC__
   !$omp parallel do if(parallelize) schedule(runtime) collapse(2) &
   !$omp private(jat, jc, jj, ii, er, el, gr, gl, sr, sl, rr, rl, dr, dl, &
   !$omp& mol, chk, t1, w1)
#endif
   do iat = 1, mol0%n
      do ic = 1, 3
         ii = 3*(iat - 1) + ic
         er = 0.0_wp
         el = 0.0_wp
         gr = 0.0_wp
         gl = 0.0_wp

         mol = mol0
         mol%xyz(ic, iat) = mol0%xyz(ic, iat) + step
         chk = chk0
         call calc%singlepoint(env, mol, chk, -1, .true., er, gr, sr, egap, rr)
         dr = rr%dipole

         mol = mol0
         mol%xyz(ic, iat) = mol0%xyz(ic, iat) - step
         chk = chk0
         call calc%singlepoint(env, mol, chk, -1, .true., el, gl, sl, egap, rl)
         dl = rl%dipole

         dipgrad(:, ii) = (dr - dl) * step2
         do jat = 1, mol0%n
            do jc = 1, 3
               jj = 3*(jat - 1) + jc
               hessian(jj, ii) = hessian(jj, ii) &
                  & + (gr(jc, jat) - gl(jc, jat)) * step2
            end do
         end do

         if (iat == 3 .and. ic == 3) then
            !$omp critical(xtb_numdiff2)
            call timing(t1, w1)
            write(*,'("estimated CPU  time",F10.2," min")') &
               & 0.3333333_wp*mol0%n*(t1-t0)/60.0_wp
            write(*,'("estimated wall time",F10.2," min")') &
               & 0.3333333_wp*mol0%n*(w1-w0)/60.0_wp
            !$omp end critical(xtb_numdiff2)
         endif

      end do
   end do

end subroutine numdiff2_all


!> Evaluate hessian by finite difference for a list of atoms
subroutine numdiff2_list(env, mol0, chk0, calc, list, step, hessian, dipgrad, parallelize)
   character(len=*), parameter :: source = "hessian_numdiff_numdiff2"
   !> Computation environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol0
   !> Restart data
   type(TRestart), intent(in) :: chk0
   !> Single point calculator
   class(TCalculator), intent(inout) :: calc
   !> List of atoms to displace
   integer, intent(in) :: list(:)
   !> Step size for numerical differentiation
   real(wp), intent(in) :: step
   !> Array to add Hessian to
   real(wp), intent(inout) :: hessian(:, :)
   !> Array to add dipole gradient to
   real(wp), intent(inout) :: dipgrad(:, :)
   !> Parallelize over displacements
   logical, intent(in) :: parallelize

   integer :: iat, jat, kat, ic, jc, ii, jj
   type(TMolecule), allocatable :: mol
   type(TRestart), allocatable :: chk
   real(wp) :: er, el, dr(3), dl(3), sr(3, 3), sl(3, 3), egap, step2
   real(wp) :: t0, t1, w0, w1
   real(wp), allocatable :: gr(:, :), gl(:, :)
   type(scc_results) :: rr, rl

   call timing(t0, w0)
   step2 = 0.5_wp / step
   allocate(gr(3, mol0%n), gl(3, mol0%n))

! This OpenMP statements lead to invalid LLVM-IRwith NVHPC (20.7 to 20.11)
#ifndef __PGIC__
   !$omp parallel do if(parallelize) schedule(runtime) collapse(2) &
   !$omp private(jat, jc, jj, ii, er, el, gr, gl, sr, sl, rr, rl, dr, dl, &
   !$omp& mol, chk, t1, w1)
#endif
   do kat = 1, size(list)
      do ic = 1, 3
         iat = list(kat)
         ii = 3*(iat - 1) + ic
         er = 0.0_wp
         el = 0.0_wp
         gr = 0.0_wp
         gl = 0.0_wp

         mol = mol0
         mol%xyz(ic, iat) = mol0%xyz(ic, iat) + step
         chk = chk0
         call calc%singlepoint(env, mol, chk, -1, .true., er, gr, sr, egap, rr)
         dr = rr%dipole

         mol = mol0
         mol%xyz(ic, iat) = mol0%xyz(ic, iat) - step
         chk = chk0
         call calc%singlepoint(env, mol, chk, -1, .true., el, gl, sl, egap, rl)
         dl = rl%dipole

         dipgrad(:, ii) = (dr - dl) * step2
         do jat = 1, mol0%n
            do jc = 1, 3
               jj = 3*(jat - 1) + jc
               hessian(jj, ii) = hessian(jj, ii) &
                  & + (gr(jc, jat) - gl(jc, jat)) * step2
            end do
         end do

         if (kat == 3 .and. ic == 3) then
            !$omp critical(xtb_numdiff2)
            call timing(t1, w1)
            write(*,'("estimated CPU  time",F10.2," min")') &
               & 0.3333333_wp*size(list)*(t1-t0)/60.0_wp
            write(*,'("estimated wall time",F10.2," min")') &
               & 0.3333333_wp*size(list)*(w1-w0)/60.0_wp
            !$omp end critical(xtb_numdiff2)
         endif

      end do
   end do

end subroutine numdiff2_list


end module xtb_freq_numdiff
