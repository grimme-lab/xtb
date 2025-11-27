! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

!> abstract calculator that hides implementation details from calling codes
module xtb_type_calculator
   use xtb_mctc_accuracy, only : wp
   use xtb_solv_model, only : TSolvModel
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   implicit none

   public :: TCalculator
   private


   !> Base calculator
   type, abstract :: TCalculator

      real(wp) :: accuracy
      logical :: lSolv = .false.
      type(TSolvModel), allocatable :: solvation
      logical :: threadsafe = .true.

   contains

      !> Perform single point calculation
      procedure(singlepoint), deferred :: singlepoint

      !> Perform hessian calculation
      procedure :: hessian

      !> Write informative printout
      procedure(writeInfo), deferred :: writeInfo

   end type TCalculator


   abstract interface
      subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
            & energy, gradient, sigma, hlgap, results)
         import :: TCalculator, TEnvironment, TMolecule, TRestart, wp
         import :: scc_results

         !> Calculator instance
         class(TCalculator), intent(inout) :: self

         !> Computational environment
         type(TEnvironment), intent(inout) :: env

         !> Molecular structure data
         type(TMolecule), intent(inout) :: mol

         !> Wavefunction data
         type(TRestart), intent(inout) :: chk

         !> Print level for IO
         integer, intent(in) :: printlevel

         !> Restart from previous results
         logical, intent(in) :: restart

         !> Total energy
         real(wp), intent(out) :: energy

         !> Molecular gradient
         real(wp), intent(out) :: gradient(:, :)

         !> Strain derivatives
         real(wp), intent(out) :: sigma(:, :)

         !> HOMO-LUMO gap
         real(wp), intent(out) :: hlgap

         !> Detailed results
         type(scc_results), intent(out) :: results

      end subroutine singlepoint


      subroutine writeInfo(self, unit, mol)
         import :: TCalculator, TMolecule

         !> Calculator instance
         class(TCalculator), intent(in) :: self

         !> Unit for I/O
         integer, intent(in) :: unit

         !> Molecular structure data
         type(TMolecule), intent(in) :: mol

      end subroutine writeInfo
   end interface


contains


!> Evaluate hessian by finite difference for all atoms
subroutine hessian(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad)
   character(len=*), parameter :: source = "hessian_numdiff_numdiff2"
   !> Single point calculator
   class(TCalculator), intent(inout) :: self
   !> Computation environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol0
   !> Restart data
   type(TRestart), intent(in) :: chk0
   !> List of atoms to displace
   integer, intent(in) :: list(:)
   !> Step size for numerical differentiation
   real(wp), intent(in) :: step
   !> Array to add Hessian to
   real(wp), intent(inout) :: hess(:, :)
   !> Array to add dipole gradient to
   real(wp), intent(inout) :: dipgrad(:, :)
   !> Array to add polarizability gradient to
   real(wp), intent(inout), optional :: polgrad(:, :)

   integer :: iat, jat, kat, ic, jc, ii, jj
   real(wp) :: er, el, dr(3), dl(3), sr(3, 3), sl(3, 3), egap, step2
   real(wp) :: alphal(3, 3), alphar(3, 3)
   real(wp) :: t0, t1, w0, w1
   real(wp), allocatable :: gr(:, :), gl(:, :)

   call timing(t0, w0)
   step2 = 0.5_wp / step

   !$omp parallel if(self%threadsafe) default(none) &
   !$omp shared(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad, step2, t0, w0) &
   !$omp private(kat, iat, jat, jc, jj, ii, er, el, egap, gr, gl, sr, sl, dr, dl, alphar, alphal, &
   !$omp& t1, w1)

   allocate(gr(3, mol0%n), gl(3, mol0%n))

   !$omp do collapse(2) schedule(runtime)
   do kat = 1, size(list)
      do ic = 1, 3

         iat = list(kat)
         ii = 3*(iat - 1) + ic
         er = 0.0_wp
         el = 0.0_wp
         gr = 0.0_wp
         gl = 0.0_wp

         call hessian_point(self, env, mol0, chk0, iat, ic, +step, er, gr, sr, egap, dr, alphar)
         call hessian_point(self, env, mol0, chk0, iat, ic, -step, el, gl, sl, egap, dl, alphal)

         if (present(polgrad)) then
            polgrad(1, ii) = (alphar(1, 1) - alphal(1, 1)) * step2
            polgrad(2, ii) = (alphar(1, 2) - alphal(1, 2)) * step2
            polgrad(3, ii) = (alphar(2, 2) - alphal(2, 2)) * step2
            polgrad(4, ii) = (alphar(1, 3) - alphal(1, 3)) * step2
            polgrad(5, ii) = (alphar(2, 3) - alphal(2, 3)) * step2
            polgrad(6, ii) = (alphar(3, 3) - alphal(3, 3)) * step2
         endif

         dipgrad(:, ii) = (dr - dl) * step2

         do jat = 1, mol0%n
            do jc = 1, 3
               jj = 3*(jat - 1) + jc
               hess(jj, ii) = hess(jj, ii) &
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
   !$omp end parallel
end subroutine hessian

subroutine hessian_point(self, env, mol0, chk0, iat, ic, step, energy, gradient, sigma, egap, dipole, alpha)
   class(TCalculator), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol0
   type(TRestart), intent(in) :: chk0
   integer, intent(in) :: ic, iat
   real(wp), intent(in) :: step
   real(wp), intent(out) :: energy
   real(wp), intent(out) :: gradient(:, :)
   real(wp), intent(out) :: sigma(3, 3)
   real(wp), intent(out) :: egap
   real(wp), intent(out) :: dipole(3)
   real(wp), intent(out) :: alpha(3, 3)

   ! internal variables
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(scc_results) :: res

   call mol%copy(mol0)
   mol%xyz(ic, iat) = mol0%xyz(ic, iat) + step
   call chk%copy(chk0)
   call self%singlepoint(env, mol, chk, -1, .true., energy, gradient, sigma, egap, res)

   dipole = res%dipole
   alpha(:, :) = res%alpha

end subroutine hessian_point


!> Evaluate hessian using O1NumHess algorithm
!> Implementation according to Wang et al. (https://doi.org/10.48550/arXiv.2508.07544)
subroutine odlrhessian(self, env, mol0, chk0, list, step, hess)
   character(len=*), parameter :: source = "hessian_odlr"
   !> Single point calculator
   class(TCalculator), intent(inout) :: self
   !> Computation environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol0
   !> Restart data
   type(TRestart), intent(in) :: chk0
   !> List of atoms to displace
   integer, intent(in) :: list(:)
   !> Step size for numerical differentiation
   real(wp), intent(in) :: step
   !> Array to add Hessian to
   real(wp), intent(inout) :: hess(:, :)

   type(scc_results) :: res
   real(wp), allocatable :: distmat(:, :), h0(:, :), displdir(:, :), g(:, :), g0(:, :)
   real(wp) :: energy, sigma(3, 3), egap, dist, barycenter(3), inertia(3), ax(3, 3), cross(3)
   real(wp) :: identity3(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   logical :: linear
   integer :: N, i, j, k, Ntr
   
   ! ========== INITIALIZATION ==========
   N = 3*mol0%n
   allocate(distmat(N, N), h0(N, N), displdir(N, N), g(N, N))

   ! hessian initial guess
   do i = 1, N
      ! unity as initial guess for hessian
      ! TODO: do we have a better guess? Swart model hessian
      h0(i, i) = 1.0_wp
   end do

   ! calculate unperturbed gradient
   call self%singlepoint(env, mol0, chk0, -1, .true., energy, g0, sigma, egap, res)

   ! setup distmat
   do i = 0, mol0%n - 1
      do j = i, mol0%n - 1
         dist = mol0%dist(i, j) ! substract vdw radii
         distmat(3*i+1:3*i+3+1, 3*j+1:3*j+3+1) = dist
         distmat(3*j+1:3*j+3+1, 3*i+1:3*i+3+1) = dist
      end do
   end do

   ! set up initial displdir with trans, rot, and totally symmetric vib mode
   ! get trans, rot displacements
   Imat0 = 0.0_wp
   do i = 1, mol0%n
      vec = mol0%xyz(:, i) - barycenter(:)
      Imat0 = Imat0 + matmul(vec, transpose(vec))
   end do
   Imat = Imat0 * identity3
   do i = 1, 3
      do j = 1, 3
         do k = 1, mol0%n
            Imat(i, j) = (mol0%xyz(i, k) - barycenter(i)) * (mol0%xyz(j, k) - barycenter(j))
         end do
      end do
   end do

   ! TODO: compute interial moment and inertial tensor

   ! translational displacements
   do i = 0, N - 1
      displdir(3*i+1, 1) = 1.0_wp/sqrt(mol0%n)
      displdir(3*i+1+1, 2) = 1.0_wp/sqrt(mol0%n)
      displdir(3*i+2+1, 3) = 1.0_wp/sqrt(mol0%n)
   end do
   
   ! rotational displacements
   Ntr = 3
   do i = 1, 3
      if (inertia(i) < 1e-4_wp) cycle
      do j = 0, mol0%n - 1
         crossprod(ax(:, i), mol0%xyz(:, j) - barycenter(:), cross)
         displdir(3*j+1:3*j+3+1, Ntr) = cross
      end do
      displdir(:, Ntr) = displdir(:, Ntr) / norm2(displdir(:, Ntr))
      Ntr = Ntr + 1
      ! NOTE: is this really correct to normalize before adding more Ntr?
   end do

   ! totally symmetric vibrational displacement
   do i = 0, mol0%n - 1
      displdir(3*i+1:3*i+3+1, Ntr) = mol0%xyz(:, i) - barycenter(:) ! FIXME: Ntr index
      displdir(3*i+1:3*i+3+1, Ntr) = displdir(3*i+1:3*i+3+1, Ntr) / norm2(displdir(3*i+1:3*i+3+1, Ntr))
   end do

   ! generate remaining displdirs based on distmat and dmax
   ! TODO: compute neighbor list
   ! TODO: populate displdir
   
   ! ========== GRADIENT DERIVATIVES ==========
   ! get gradient derivatives for trans, rot and symmetric vib mode
   ! calculate remaining gradient derivatives

   ! ========== FINAL HESSIAN ==========
   ! construct hessian from local hessian and odlr correction

end subroutine odlrhessian

end module xtb_type_calculator
