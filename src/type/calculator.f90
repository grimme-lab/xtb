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
   use xtb_mctc_convert, only : autoaa
   use xtb_mctc_math, only : crossProd
   use xtb_mctc_accuracy, only : wp
   use xtb_solv_model, only : TSolvModel
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_o1numhess, only : adj_list, gen_local_hessian, &
   & lr_loop, gen_displdir, get_neighbor_list
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

      !> Perform ODLR approximated numerical hessian
      procedure :: odlrhessian

      procedure :: get_gradient_derivs

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

   !$omp parallel if (self%threadsafe) default(none) &
   !$omp shared(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad, step2, t0, w0) &
   !$omp private(kat, iat, jat, jc, jj, ii, er, el, egap, gr, gl, sr, sl, dr, dl, alphar, alphal, &
   !$omp& t1, w1)

   allocate(gr(3, mol0%n), gl(3, mol0%n))

   !$omp do collapse(2) schedule(runtime)
   do kat = 1, size(list)
      do ic = 1, 3

         iat = list(kat)
         ii = 3 * (iat - 1) + ic
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
         end if

         dipgrad(:, ii) = (dr - dl) * step2

         do jat = 1, mol0%n
            do jc = 1, 3
               jj = 3 * (jat - 1) + jc
               hess(jj, ii) = hess(jj, ii) &
                  & + (gr(jc, jat) - gl(jc, jat)) * step2
            end do
         end do

         if (kat == 3 .and. ic == 3) then
            !$omp critical(xtb_numdiff2)
            call timing(t1, w1)
            write(*, '("estimated CPU  time",F10.2," min")') &
               & 0.3333333_wp * size(list) * (t1 - t0) / 60.0_wp
            write(*, '("estimated wall time",F10.2," min")') &
               & 0.3333333_wp * size(list) * (w1 - w0) / 60.0_wp
            !$omp end critical(xtb_numdiff2)
         end if

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

!> Implementation according to Wang et al. (https://doi.org/10.48550/arXiv.2508.07544)
subroutine odlrhessian(self, env, mol0, chk0, list, step, displdir0, g, hess)
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
   !> Displacement directions
   real(wp), intent(in) :: displdir0(:, :)
   !> Gradients
   real(wp), allocatable, intent(inout) :: g(:, :)
   !> Array to add Hessian to
   real(wp), intent(inout) :: hess(:, :)
   !> Array for displacement directions
   real(wp), allocatable :: displdir(:, :)
   
   ! UFF vdw radii - could be replaced with any other vdw radii i guess
   real(wp), parameter :: vdw_radii(1:103) = [ &
      2.886_wp, 2.362_wp, 2.451_wp, 2.745_wp, 4.083_wp, 3.851_wp, 3.66_wp, 3.5_wp, 3.364_wp, &
      3.243_wp, 2.983_wp, 3.021_wp, 4.499_wp, 4.295_wp, 4.147_wp, 4.035_wp, 3.947_wp, 3.868_wp, &
      3.812_wp, 3.399_wp, 3.295_wp, 3.175_wp, 3.144_wp, 3.023_wp, 2.961_wp, 2.912_wp, 2.872_wp, &
      2.834_wp, 3.495_wp, 2.763_wp, 4.383_wp, 4.28_wp, 4.23_wp, 4.205_wp, 4.189_wp, 4.141_wp, &
      4.114_wp, 3.641_wp, 3.345_wp, 3.124_wp, 3.165_wp, 3.052_wp, 2.998_wp, 2.963_wp, 2.929_wp, &
      2.899_wp, 3.148_wp, 2.848_wp, 4.463_wp, 4.392_wp, 4.42_wp, 4.47_wp, 4.5_wp, 4.404_wp, &
      4.517_wp, 3.703_wp, 3.522_wp, 3.556_wp, 3.606_wp, 3.575_wp, 3.547_wp, 3.52_wp, 3.493_wp, &
      3.368_wp, 3.451_wp, 3.428_wp, 3.409_wp, 3.391_wp, 3.374_wp, 3.355_wp, 3.64_wp, 3.141_wp, &
      3.17_wp, 3.069_wp, 2.954_wp, 3.12_wp, 2.84_wp, 2.754_wp, 3.293_wp, 2.705_wp, 4.347_wp, &
      4.297_wp, 4.37_wp, 4.709_wp, 4.75_wp, 4.765_wp, 4.9_wp, 3.677_wp, 3.478_wp, 3.396_wp, &
      3.424_wp, 3.395_wp, 3.424_wp, 3.424_wp, 3.381_wp, 3.326_wp, 3.339_wp, 3.313_wp, 3.299_wp, &
      3.286_wp, 3.274_wp, 3.248_wp, 3.236_wp &
   ] / 2.0_wp / autoaa
   real(wp), parameter :: dmax = 1.0_wp, eps = 1.0e-8_wp, eps2 = 1.0e-15_wp

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(scc_results) :: res
   type(adj_list), allocatable :: neighborlist(:)
   real(wp), allocatable :: distmat(:, :), h0(:, :), h0v(:), tmp_grad(:, :), g0(:), x(:), xyz(:, :), gr(:, :), gl(:, :), gtmp(:, :)
   real(wp) :: energy, sigma(3, 3), egap, dist, barycenter(3), inertia(3), ax(3, 3), cross(3), Imat0, query(1), displmax
   real(wp) :: identity3(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1],[3, 3]), final_err
   logical :: linear
   integer, allocatable :: nbcounts(:)
   integer :: N, i, j, k, Ntr, info, lwork, ndispl_final, max_nb, ginit
   
   ! ========== INITIALIZATION ==========
   ! NOTE: maybe this needs to go to numhess?
   N = 3 * mol0%n
   ginit = size(g, 2)

   call mol%copy(mol0)
   call chk%copy(chk0)

   ! hessian initial guess
   ! allocate(h0v(N*(N+1)/2))
   ! call ddvopt(mol0%xyz, mol0%n, h0v, mol0%at, 20.0_wp)
   ! h0 = unpack_sym(h0v, mask, N) ! TODO: mask??
   allocate(h0(N, N))
   h0 = 0.0_wp
   do i = 1, N
      h0(i, i) = 1.0_wp
   end do

   ! calculate unperturbed gradient
   write(env%unit, '(A)') "Calculating unperturbed gradient"
   call self%singlepoint(env, mol, chk, -1, .true., energy, tmp_grad, sigma, egap, res)

   ! gradients need to be flattened since hessian is also "flat"
   g0 = reshape(tmp_grad,[N])

   ! setup distmat
   write(env%unit, '(A)') "Distmat setup"
   allocate(distmat(N, N))
   ! do i = 1, mol0%n
   !    do j = i, mol0%n
   !       ! effective distmat
   !       dist = mol0%dist(i, j) - vdw_radii(mol0%at(i)) - vdw_radii(mol0%at(j))
   !       distmat(3 * i - 2:3 * i, 3 * j - 2:3 * j) = dist
   !       distmat(3 * j - 2:3 * j, 3 * i - 2:3 * i) = dist
   !    end do
   ! end do
   do i = 1, mol0%n
      do j = 1, mol0%n
         distmat(i, j) = abs(i - j)
      end do
   end do

   ! allocate(displdir(N, N))
   ! ! set up initial displdir with trans, rot, and totally symmetric vib mode first
   ! ! translational displacements
   ! do i = 1, N
   !    displdir(3 * i - 2, 1) = 1.0_wp / sqrt(real(mol0%n, wp))
   !    displdir(3 * i - 1, 2) = 1.0_wp / sqrt(real(mol0%n, wp))
   !    displdir(3 * i, 3) = 1.0_wp / sqrt(real(mol0%n, wp))
   ! end do
   !
   ! ! calculate inertial moment and axes
   ! barycenter = sum(mol0%xyz, dim=2) / real(mol0%n, wp)
   ! Imat0 = 0.0_wp
   ! do i = 1, mol0%n
   !    vec = mol0%xyz(:, i) - barycenter(:)
   !    Imat0 = Imat0 + matmul(vec, transpose(vec))
   ! end do
   ! ax = Imat0 * identity3
   ! do i = 1, 3
   !    do j = 1, 3
   !       do k = 1, mol0%n
   !          ax(i, j) = ax(i, j) - (mol0%xyz(i, k) - barycenter(i)) * (mol0%xyz(j, k) - barycenter(j))
   !       end do
   !    end do
   ! end do
   !
   ! lwork = -1
   ! call dsyev('V', 'U', 3, ax, 3, inertia, query, lwork, info)
   ! lwork = int(query(1))
   ! allocate(aux(lwork))
   ! call dsyev('V', 'U', 3, ax, 3, inertia, aux, lwork, info)
   !
   ! ! rotational displacements
   ! Ntr = 3 ! number of translational and rotational degrees of freedom
   ! do i = 1, 3
   !    if (inertia(i) < 1e-4_wp) cycle
   !    do j = 1, mol0%n
   !       crossProd(ax(:, i), mol0%xyz(:, j) - barycenter(:), cross)
   !       displdir(3 * j - 2:3 * j, Ntr + 1) = cross
   !    end do
   !    displdir(:, Ntr + 1) = displdir(:, Ntr + 1) / norm2(displdir(:, Ntr + 1))
   !    Ntr = Ntr + 1
   ! end do
   !
   ! ! totally symmetric vibrational displacement
   ! do i = 1, mol0%n
   !    displdir(3 * i - 2:3 * i, Ntr + 1) = mol0%xyz(:, i) - barycenter(:)
   !    displdir(3 * i - 2:3 * i, Ntr + 1) = displdir(3 * i - 2:3 * i, Ntr + 1) / norm2(displdir(3 * i - 2:3 * i, Ntr + 1))
   ! end do
   !
   ! TODO: get gradient derivs along rotational displacements - in numhess

   ! generate remaining displdirs based on distmat and dmax
   ! compute neighbor list
   write(env%unit, '(A)') "Getting neighbor list"
   call get_neighbor_list(distmat, dmax, neighborlist)
   allocate(nbcounts(N))
   max_nb = 0
   do i = 1, N
      nbcounts(i) = size(neighborlist(i)%neighbors)
      if (nbcounts(i) > max_nb) max_nb = nbcounts(i)
   end do
   
   ! TODO: orthonormalize displdir?
   ! populate displdir
   write(env%unit, '(A)') "Generating displacements"
   Ntr = 0
   call gen_displdir(N, Ntr, h0, max_nb, neighborlist, nbcounts, eps, eps2, displdir, ndispl_final)


   ! allocate g with correct size
   allocate(gtmp(N, ndispl_final))
   gtmp(:, 1:ginit) = g
   call move_alloc(gtmp, g)
   g(:, ginit+1:ndispl_final) = 0.0_wp
   
   ! ========== GRADIENT DERIVATIVES ==========
   write(env%unit, '(A)') "Calculating gradient derivatives"
   call get_gradient_derivs(self, env, step, Ntr, ndispl_final, displdir, mol0, chk0, g0, g)

   ! ========== FINAL HESSIAN ==========
   ! construct hessian from local hessian and odlr correction
   ! compute local hessian
   write(env%unit, '(A)') "Computing local Hessian"
   call gen_local_hessian(distmat, displdir, g, dmax, hess)

   ! compute low rank correction
   write(env%unit, '(A)') "Computing low rank correction"
   call lr_loop(ndispl_final, g, hess, displdir, final_err)

end subroutine odlrhessian

subroutine get_gradient_derivs(self, env, step, Ntr, ndispl_final, displdir, mol0, chk0, g0, g)
   class(TCalculator), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(wp), intent(in) :: step
   integer, intent(in) :: Ntr, ndispl_final
   real(wp), intent(in) :: displdir(:, :)
   type(TMolecule), intent(in) :: mol0
   type(TRestart), intent(in) :: chk0
   real(wp), intent(in) :: g0(:)
   real(wp), intent(inout) :: g(:, :)

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(scc_results) :: res
   integer :: i, N, j
   real(wp) :: displmax, sigma(3, 3), energy, egap
   real(wp), allocatable :: tmp_grad(:, :)
   real(wp), allocatable :: x(:)

   allocate(tmp_grad(3, mol0%n))
   tmp_grad = 0.0_wp
   N = 3 * mol0%n
   do i = Ntr + 1, ndispl_final
      displmax = maxval(abs(displdir(:, i)))
      ! TODO: what about double sided stuff?
      call mol%copy(mol0)
      call chk%copy(chk0)
      x = reshape(mol0%xyz,[N])
      x = x + step * displdir(:, i) / displmax
      mol%xyz = reshape(x,[3, mol0%n])
      call self%singlepoint(env, mol, chk, -1, .false., energy, tmp_grad, sigma, egap, res)
      g(:, i) = reshape(tmp_grad,[N])
      g(:, i) = (g(:, i) - g0(:)) / step * displmax
   end do
end subroutine get_gradient_derivs

end module xtb_type_calculator
