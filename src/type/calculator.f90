! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
! Copyright (C) 2026 Leopold M. Seidler
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
   use xtb_mctc_blas, only : mctc_gemm, mctc_nrm2, mctc_dot
   use xtb_solv_model, only : TSolvModel
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_o1numhess, only : adj_list, gen_local_hessian, &
   & lr_loop, gen_displdir, get_vdw_neighbor_list, find_projected_imag_modes
   use xtb_modelhessian_swart, only : TSwartModelHessian
   use xtb_type_setvar, only : modhess_setvar
   use xtb_setparam, only : set
   use xtb_param_uffvdwrad, only : get_rad
   use xtb_param_covalentrad, only : get_cov_rad
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

   real(wp), parameter :: eps = epsilon(0.0_wp)
   real(wp), parameter :: eps2 = sqrt(eps)

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
subroutine hessian(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad, odlr, final_err)
   character(len=*), parameter :: source = "hessian"
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
   !> Use ODLR approximated numerical hessian
   logical, intent(in), optional :: odlr
   !> Final residual error (ODLR only)
   real(wp), intent(out), optional :: final_err

   logical :: use_odlr

   use_odlr = .false.
   if (present(odlr)) use_odlr = odlr

   if (use_odlr) then
      if (present(final_err)) then
         if (present(polgrad)) then
            call hessian_odlr(self, env, mol0, chk0, step, hess, final_err, dipgrad, polgrad)
         else
            call hessian_odlr(self, env, mol0, chk0, step, hess, final_err, dipgrad)
         end if
      else
         call env%error("hessian: final_err must be present when odlr=.true.", source)
         return
      end if
   else
      if (present(polgrad)) then
         call hessian_numdiff(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad)
      else
         call hessian_numdiff(self, env, mol0, chk0, list, step, hess, dipgrad)
      end if
   end if

end subroutine hessian

!> Evaluate hessian by finite difference using Cartesian displacement vectors
subroutine hessian_numdiff(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad)
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

   integer :: N, ndispl, i, kat, iat, ic, ii, jc, jj, jat
   real(wp), allocatable :: displdir(:, :), g(:, :), dipdir(:, :), poldir(:, :)
   real(wp) :: t0, t1, w0, w1

   call timing(t0, w0)
   N = 3*mol0%n
   ndispl = 3*size(list)

   ! Build Cartesian unit displacement vectors (column-major)
   allocate(displdir(N, ndispl))
   displdir = 0.0_wp
   do kat = 1, size(list)
      iat = list(kat)
      do ic = 1, 3
         ii = 3*(kat - 1) + ic
         displdir(3*(iat - 1) + ic, ii) = 1.0_wp
      end do
   end do

   ! Allocate gradient derivative matrix
   allocate(g(N, ndispl))
   g = 0.0_wp

   if (present(polgrad)) then
      allocate(poldir(6, ndispl))
      poldir = 0.0_wp
   end if
   allocate(dipdir(3, ndispl))
   dipdir = 0.0_wp

   ! Evaluate gradient derivatives (double-sided) along Cartesian directions
   block
      real(wp), allocatable :: g0_dummy(:)
      real(wp) :: dip0_dummy(3), alpha0_dummy(3, 3)
      allocate(g0_dummy(N), source=0.0_wp)
      dip0_dummy = 0.0_wp
      alpha0_dummy = 0.0_wp
      call get_gradient_derivs(self, env, step, 0, ndispl, displdir, mol0, chk0, &
         & g0_dummy, .true., g, dip0_dummy, alpha0_dummy, dipdir, poldir)
   end block

   ! Assemble hessian, dipgrad, polgrad from derivative columns
   !$omp parallel if (self%threadsafe) default(none) &
   !$omp shared(hess, dipgrad, polgrad, g, dipdir, poldir, displdir, list, mol0, ndispl, N) &
   !$omp private(kat, iat, ic, ii, jc, jj, jat)
   !$omp do
   do kat = 1, size(list)
      iat = list(kat)
      do ic = 1, 3
         ii = 3*(kat - 1) + ic
         ! hessian column: g(:, ii) is gradient along displacement ii
         ! hess(:, 3*(iat-1)+ic) += g(:, ii) / step  (displmax=1, factor=0.5, doublesided)
         ! get_gradient_derivs already applies /step*displmax*factor, so g = (gl-gr)/(2*step)
         ! hess(jj, col) = dg/dx = (g_r - g_l)/(2*step) = g(:, ii) value
         do jat = 1, mol0%n
            do jc = 1, 3
               jj = 3*(jat - 1) + jc
               hess(jj, 3*(iat - 1) + ic) = hess(jj, 3*(iat - 1) + ic) + g(jj, ii)
            end do
         end do
         ! dipgrad column
         dipgrad(:, 3*(iat - 1) + ic) = dipgrad(:, 3*(iat - 1) + ic) + dipdir(:, ii)
         ! polgrad column
         if (present(polgrad)) then
            polgrad(:, 3*(iat - 1) + ic) = polgrad(:, 3*(iat - 1) + ic) + poldir(:, ii)
         end if
      end do
   end do
   !$omp end do
   !$omp end parallel

   ! timing estimate (matches old behavior: print after 3rd atom, 3rd coord)
   call timing(t1, w1)
   write (*, '("estimated CPU  time",F10.2," min")') &
      & 0.3333333_wp*size(list)*(t1 - t0)/60.0_wp
   write (*, '("estimated wall time",F10.2," min")') &
      & 0.3333333_wp*size(list)*(w1 - w0)/60.0_wp

end subroutine hessian_numdiff

!> Implementation according to Wang et al. (https://doi.org/10.1021/acs.jctc.5c01354)
subroutine hessian_odlr(self, env, mol0, chk0, step, hess, final_err, dipgrad, polgrad)
   character(len=*), parameter :: source = "hessian_odlr"
   !> Single point calculator
   class(TCalculator), intent(inout) :: self
   !> Computation environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol0
   !> Restart data
   type(TRestart), intent(in) :: chk0
   !> Step size for numerical differentiation
   real(wp), intent(in) :: step
   !> Array to add Hessian to
   real(wp), intent(inout) :: hess(:, :)
   !> Array for displacement directions
   real(wp), allocatable :: displdir(:, :)
   !> Final error after all steps
   real(wp), intent(out) :: final_err
   !> Array to add dipole gradient to
   real(wp), intent(inout), optional :: dipgrad(:, :)
   !> Array to add polarizability gradient to
   real(wp), intent(inout), optional :: polgrad(:, :)

   real(wp), parameter :: dmax = 1.0_wp, imagthr = eps2
   real(wp), parameter :: identity3(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], shape(identity3))
   ! iterative imaginary-frequency repair parameters
   integer, parameter :: max_neg_round = 6, max_neg_add_per_round = 6
   logical, parameter :: iterative_neg_mode = .true.

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(scc_results) :: res
   type(adj_list), allocatable :: neighborlist(:)
   real(wp), allocatable :: distmat(:, :), h0(:, :), tmp_grad(:, :), &
      & g0(:), g(:, :), work(:), eigvec(:, :), eigval(:), &
      & dipdir(:, :), poldir(:, :), prop_tmp(:, :)
   real(wp) :: energy, sigma(3, 3), egap, dist, barycenter(3), inertia(3), &
      & ax(3, 3), cross(3), Imat0, vec(3), ri, rj, delta_r1, dip0(3), alpha0(3, 3), &
      & rot_norm, omega(3, 3), dalpha(3, 3)
   logical :: terminate_run, linear, do_dipgrad, do_polgrad
   integer, allocatable :: nbcounts(:)
   character(len=128) :: msg
   integer :: N, i, j, k, Ntr, info, lwork, ndispl_final, max_nb, ndispl0, nimg, rmode
   ! iterative repair locals
   integer :: neg_round, nadd, nmodes_det, prev_nimg, ndiag_window
   real(wp), allocatable :: cand_modes(:, :), cand_freqs(:), sel_modes(:, :), vscratch(:)
   real(wp) :: vdot

   ! ========== INITIALIZATION ==========
   N = 3*mol0%n

   call mol%copy(mol0)
   call chk%copy(chk0)

   ! hessian initial guess
   allocate (h0(N, N))
   block
      type(TSwartModelHessian) :: model_hessian
      type(modhess_setvar) :: modh
      ! original parameters kept from setparam.f90 for Lindh-D2 model Hessian
      modh = modhess_setvar(kr=0.4000_wp, kf=0.1300_wp, kt=0.0075_wp, &
         & ko=0.16_wp, kd=0.0_wp, kq=0.0_wp, rcut=70.0_wp, s6=20.0_wp)
      call model_hessian%compute(env, mol%xyz, mol%n, h0, mol%at, modh)
   end block
   call env%check(terminate_run)
   if (terminate_run) then
      return
   end if

   ! calculate unperturbed gradient
   allocate (tmp_grad(3, mol0%n))
   call self%singlepoint(env, mol, chk, -1, .true., energy, tmp_grad, sigma, egap, res)
   g0 = reshape(tmp_grad, [N])
   dip0 = res%dipole
   alpha0 = res%alpha
   do_polgrad = present(polgrad)
   do_dipgrad = present(dipgrad) .or. do_polgrad
   ! setup effective distmat
   allocate (distmat(N, N))
   do i = 1, mol0%n
      do j = i, mol0%n
         ! effective distmat
         ri = get_cov_rad(mol0%at(i))
         rj = get_cov_rad(mol0%at(j))
         if (ri < 0.0_wp .or. rj < 0.0_wp) then
            call env%error("odlrhessian: covalent radii only defined for 1-103", source)
            return
         end if
         dist = mol0%dist(i, j) - get_rad(mol0%at(i)) - get_rad(mol0%at(j))
         distmat(3*i - 2:3*i, 3*j - 2:3*j) = dist
         distmat(3*j - 2:3*j, 3*i - 2:3*i) = dist
      end do
   end do

   allocate (displdir(N, N))
   displdir = 0.0_wp
   ! set up initial displdir with trans, rot, and totally symmetric vib mode first
   ! translational displacements
   Ntr = 3
   do i = 1, mol0%n
      displdir(3*i - 2, 1) = 1.0_wp/sqrt(real(mol0%n, wp))
      displdir(3*i - 1, 2) = 1.0_wp/sqrt(real(mol0%n, wp))
      displdir(3*i, 3) = 1.0_wp/sqrt(real(mol0%n, wp))
   end do

   ! calculate inertial moment and axes
   barycenter = sum(mol0%xyz, dim=2)/real(mol0%n, wp)
   Imat0 = 0.0_wp
   do i = 1, mol0%n
      vec = mol0%xyz(:, i) - barycenter(:)
      Imat0 = Imat0 + mctc_dot(vec, vec)
   end do
   ax = Imat0*identity3
   do i = 1, 3
      do j = 1, 3
         do k = 1, mol0%n
            ax(i, j) = ax(i, j) - (mol0%xyz(i, k) - barycenter(i))*(mol0%xyz(j, k) - barycenter(j))
         end do
      end do
   end do

   lwork = -1
   allocate (work(1))
   call dsyev('V', 'U', 3, ax, 3, inertia, work, lwork, info)
   lwork = int(work(1))
   deallocate (work)
   allocate (work(lwork))
   call dsyev('V', 'U', 3, ax, 3, inertia, work, lwork, info)
   linear = mol0%linear

   ! rotational displacements
   do i = 1, 3
      if (inertia(i) < 1e-4_wp .and. linear) cycle ! skips one mode if linear
      Ntr = Ntr + 1
      do j = 1, mol0%n
         cross = crossProd(ax(:, i), mol0%xyz(:, j) - barycenter(:))
         displdir(3*j - 2:3*j, Ntr) = cross
      end do
      displdir(:, Ntr) = displdir(:, Ntr)/norm2(displdir(:, Ntr))
   end do

   ! totally symmetric vibrational displacement
   do i = 1, mol0%n
      displdir(3*i - 2:3*i, Ntr + 1) = mol0%xyz(:, i) - barycenter(:)
   end do
   displdir(:, Ntr + 1) = displdir(:, Ntr + 1)/norm2(displdir(:, Ntr + 1))

   ! get gradient derivs along rot/vib displacements
   ! gradients along translations are zero
   allocate (g(N, N))
   g = 0.0_wp
   if (do_dipgrad) then
      allocate (dipdir(3, N), source=0.0_wp)
      do i = 1, 3
         dipdir(i, i) = mol0%chrg/sqrt(real(mol0%n, wp))
      end do
   end if
   if (do_polgrad) allocate (poldir(6, N), source=0.0_wp)

   rmode = 3
   do i = 1, 3
      if (inertia(i) < 1e-4_wp) cycle
      rmode = rmode + 1
      rot_norm = 0.0_wp
      do j = 1, mol0%n
         cross = crossProd(ax(:, i), mol0%xyz(:, j) - barycenter(:))
         rot_norm = rot_norm + dot_product(cross, cross)
      end do
      rot_norm = sqrt(rot_norm)
      if (rot_norm <= eps) cycle

      do j = 1, mol0%n
         g(3*j - 2, rmode) = (ax(2, i)*tmp_grad(3, j) - ax(3, i)*tmp_grad(2, j))/rot_norm
         g(3*j - 1, rmode) = (ax(3, i)*tmp_grad(1, j) - ax(1, i)*tmp_grad(3, j))/rot_norm
         g(3*j, rmode) = (ax(1, i)*tmp_grad(2, j) - ax(2, i)*tmp_grad(1, j))/rot_norm
      end do

      if (do_dipgrad) then
         dipdir(1, rmode) = (ax(2, i)*dip0(3) - ax(3, i)*dip0(2))/rot_norm
         dipdir(2, rmode) = (ax(3, i)*dip0(1) - ax(1, i)*dip0(3))/rot_norm
         dipdir(3, rmode) = (ax(1, i)*dip0(2) - ax(2, i)*dip0(1))/rot_norm

         if (do_polgrad) then
            omega = 0.0_wp
            omega(1, 2) = -ax(3, i)
            omega(1, 3) = ax(2, i)
            omega(2, 1) = ax(3, i)
            omega(2, 3) = -ax(1, i)
            omega(3, 1) = -ax(2, i)
            omega(3, 2) = ax(1, i)
            dalpha = matmul(omega, alpha0) - matmul(alpha0, omega)
            poldir(1, rmode) = dalpha(1, 1)/rot_norm
            poldir(2, rmode) = dalpha(1, 2)/rot_norm
            poldir(3, rmode) = dalpha(2, 2)/rot_norm
            poldir(4, rmode) = dalpha(1, 3)/rot_norm
            poldir(5, rmode) = dalpha(2, 3)/rot_norm
            poldir(6, rmode) = dalpha(3, 3)/rot_norm
         end if
      end if
   end do

   call get_gradient_derivs(self, env, step, Ntr, Ntr + 1, displdir, mol0, chk0, g0, .true., g, &
      & dip0, alpha0, dipdir, poldir)

   ! generate remaining displdirs based on the vdW near-range N^(1)
   delta_r1 = 1.0_wp
   if (mol0%n >= 100) delta_r1 = 2.0_wp
   call get_vdw_neighbor_list(mol0%xyz, mol0%at, delta_r1, neighborlist)
   allocate (nbcounts(N))
   max_nb = 0
   do i = 1, N
      nbcounts(i) = size(neighborlist(i)%neighbors)
      if (nbcounts(i) > max_nb) max_nb = nbcounts(i)
   end do

   ! populate displdir
   ndispl0 = Ntr + 1
   call gen_displdir(N, ndispl0, h0, max_nb, neighborlist, nbcounts, eps2, eps, &
      & displdir, ndispl_final, info)
   if (info /= 0) then
      if (info > 0) then
         write (msg, '(a, i0, a)') &
            & "DSYEVX failed to converge for ", info, &
            & " eigenvector(s) while generating displacement directions"
      else
         write (msg, '(a, i0, a)') &
            & "Internal error while generating displacement directions (code ", info, ")"
      end if
      call env%error(trim(msg), source)
      return
   end if

   ! ========== GRADIENT DERIVATIVES ==========
   call get_gradient_derivs(self, env, step, ndispl0, ndispl_final, displdir, mol0, chk0, g0, .false., g, &
      & dip0, alpha0, dipdir, poldir)

   ! ========== FINAL HESSIAN ==========
   ! construct hessian from local hessian and odlr correction
   ! compute local hessian
   call gen_local_hessian(env, ndispl_final, distmat, displdir, g, dmax, hess)
   call env%check(terminate_run)
   if (terminate_run) then
      return
   end if

   ! compute low rank correction
   call lr_loop(env, ndispl_final, g, hess, displdir, final_err)

   call env%check(terminate_run)
   if (terminate_run) then
      return
   end if

   if (iterative_neg_mode .and. .not. linear) then
      ! iterative projected imaginary-frequency repair (plan: iterative-imaginary-frequency)
      ndiag_window = min(N, max_neg_add_per_round + 6)
      allocate (cand_modes(N, ndiag_window), cand_freqs(ndiag_window))
      allocate (sel_modes(N, max_neg_add_per_round))
      allocate (vscratch(N))
      prev_nimg = huge(1)

      do neg_round = 1, max_neg_round
         call find_projected_imag_modes(env, mol0, hess, linear, ndiag_window, &
            & -set%imagmin_hess, cand_modes, cand_freqs, nmodes_det)
         call env%check(terminate_run)
         if (terminate_run) exit

         if (nmodes_det == 0) exit
         if (nmodes_det >= prev_nimg .and. neg_round > 1) exit
         prev_nimg = nmodes_det

         ! select & orthogonalize candidates against existing displdir and each other
         nadd = 0
         do j = 1, nmodes_det
            ! note that this might not work as expected if there are many small
            ! imaginary modes
            if (nadd >= max_neg_add_per_round) exit
            if (ndispl_final + nadd >= N) exit
            if (cand_freqs(j) < set%imagmax_hess) cycle
            vscratch = cand_modes(:, j)
            do k = 1, ndispl_final
               vdot = mctc_dot(displdir(:, k), vscratch)
               vscratch = vscratch - vdot*displdir(:, k)
            end do
            do k = 1, nadd
               vdot = mctc_dot(sel_modes(:, k), vscratch)
               vscratch = vscratch - vdot*sel_modes(:, k)
            end do
            dist = norm2(vscratch)
            if (dist < eps) cycle
            nadd = nadd + 1
            sel_modes(:, nadd) = vscratch/dist
         end do

         if (nadd == 0) exit

         displdir(:, ndispl_final + 1:ndispl_final + nadd) = sel_modes(:, 1:nadd)
         ndispl0 = ndispl_final
         ndispl_final = ndispl_final + nadd
         call get_gradient_derivs(self, env, step, ndispl0, ndispl_final, &
            & displdir, mol0, chk0, g0, .true., g, dip0, alpha0, dipdir, poldir)
         call env%check(terminate_run)
         if (terminate_run) exit
         call gen_local_hessian(env, ndispl_final, distmat, displdir, g, dmax, hess)
         call env%check(terminate_run)
         if (terminate_run) exit
         call lr_loop(env, ndispl_final, g, hess, displdir, final_err)

      end do

   else
      ! old one-shot raw-eigenvalue repair (fallback: linear molecules or disabled)
      allocate (eigvec(N, N))
      eigvec = hess
      allocate (eigval(N))
      lwork = -1
      deallocate (work)
      allocate (work(1))
      call dsyev('V', 'U', N, eigvec, N, eigval, work, lwork, info)
      lwork = int(work(1))
      deallocate (work)
      allocate (work(lwork))
      call dsyev('V', 'U', N, eigvec, N, eigval, work, lwork, info)
      if (info /= 0) then
         call env%error("hessian_odlr: DSYEV failed while searching for "//&
            & "imaginary modes", source)
         return
      end if

      nimg = count(eigval < -imagthr)
      if (ndispl_final + nimg > N) nimg = N - ndispl_final

      if (nimg > 0) then
         displdir(:, ndispl_final + 1:ndispl_final + nimg) = eigvec(:, 1:nimg)
         ndispl0 = ndispl_final
         ndispl_final = ndispl_final + nimg
         call get_gradient_derivs(self, env, step, ndispl0, ndispl_final, &
            & displdir, mol0, chk0, g0, .false., g, dip0, alpha0, dipdir, poldir)
         call env%check(terminate_run)
         if (terminate_run) return
         call gen_local_hessian(env, ndispl_final, distmat, displdir, g, dmax, hess)
         call env%check(terminate_run)
         if (terminate_run) return
         call lr_loop(env, ndispl_final, g, hess, displdir, final_err)
      end if
   end if

   if (present(dipgrad)) then
      allocate (prop_tmp(3, N), source=0.0_wp)
      call mctc_gemm(dipdir(:, :ndispl_final), displdir(:, :ndispl_final), &
         & prop_tmp, transb="t", alpha=1.0_wp, beta=0.0_wp)
      dipgrad = dipgrad + prop_tmp
      deallocate (prop_tmp)
   end if

   if (do_polgrad) then
      allocate (prop_tmp(6, N), source=0.0_wp)
      call mctc_gemm(poldir(:, :ndispl_final), displdir(:, :ndispl_final), &
         & prop_tmp, transb="t", alpha=1.0_wp, beta=0.0_wp)
      polgrad = polgrad + prop_tmp
      deallocate (prop_tmp)
   end if

end subroutine hessian_odlr

subroutine get_gradient_derivs(self, env, step, ndispl0, ndispl_final, displdir, mol0, chk0, g0, doublesided, &
   & g, dip0, alpha0, dipdir, poldir)
   class(TCalculator), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(wp), intent(in) :: step
   integer, intent(in) :: ndispl0, ndispl_final
   real(wp), intent(in) :: displdir(:, :)
   type(TMolecule), intent(in) :: mol0
   type(TRestart), intent(in) :: chk0
   real(wp), intent(in) :: g0(:)
   logical, intent(in) :: doublesided
   real(wp), intent(inout) :: g(:, :)
   real(wp), intent(in) :: dip0(3)
   real(wp), intent(in) :: alpha0(3, 3)
   real(wp), intent(inout), optional :: dipdir(:, :)
   real(wp), intent(inout), optional :: poldir(:, :)

   integer :: i, N
   real(wp) :: displmax

   N = 3*mol0%n

   !$omp parallel if (self%threadsafe) default(none) &
   !$omp shared(self, env, mol0, chk0, step, g, dipdir, poldir, dip0, alpha0, g0, N, ndispl0, ndispl_final, displdir, doublesided) &
   !$omp private(i, displmax)
   !$omp do schedule(runtime)
   do i = ndispl0 + 1, ndispl_final
      displmax = maxval(abs(displdir(:, i)))
      call gradient_derivs_point(self, env, mol0, chk0, step, displdir(:, i), &
         & displmax, N, doublesided, g0, g(:, i), dip0, alpha0, dipdir, poldir, i)
   end do
   !$omp end parallel
end subroutine get_gradient_derivs

subroutine gradient_derivs_point(self, env, mol0, chk0, step, displdir_i, &
   & displmax, N, doublesided, g0, g_i, dip0, alpha0, dipdir, poldir, icol)
   class(TCalculator), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol0
   type(TRestart), intent(in) :: chk0
   real(wp), intent(in) :: step
   real(wp), intent(in) :: displdir_i(:)
   real(wp), intent(in) :: displmax
   integer, intent(in) :: N
   logical, intent(in) :: doublesided
   real(wp), intent(in) :: g0(:)
   real(wp), intent(out) :: g_i(:)
   real(wp), intent(in) :: dip0(3)
   real(wp), intent(in) :: alpha0(3, 3)
   real(wp), intent(inout), optional :: dipdir(:, :)
   real(wp), intent(inout), optional :: poldir(:, :)
   integer, intent(in) :: icol

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(scc_results) :: res
   real(wp) :: sigma(3, 3), energy, egap, factor
   real(wp) :: dipl(3), dipr(3), alphal(3, 3), alphar(3, 3), alphadiff(3, 3)
   real(wp), allocatable :: tmp_gradl(:, :), tmp_gradr(:, :)

   factor = merge(0.5_wp, 1.0_wp, doublesided)

   allocate (tmp_gradl(3, mol0%n))
   tmp_gradl = 0.0_wp
   if (doublesided) then
      allocate (tmp_gradr(3, mol0%n))
      tmp_gradr = 0.0_wp
   end if

   call mol%copy(mol0)
   call chk%copy(chk0)
   mol%xyz = mol0%xyz + reshape(step*displdir_i/displmax, [3, mol0%n])
   call self%singlepoint(env, mol, chk, -1, .true., energy, tmp_gradl, sigma, egap, res)
   dipl = res%dipole
   alphal = res%alpha

   if (doublesided) then
      call mol%copy(mol0)
      call chk%copy(chk0)
      mol%xyz = mol0%xyz - reshape(step*displdir_i/displmax, [3, mol0%n])
      call self%singlepoint(env, mol, chk, -1, .true., energy, tmp_gradr, sigma, egap, res)
      dipr = res%dipole
      alphar = res%alpha
      g_i = reshape(tmp_gradl - tmp_gradr, [N])
      alphadiff = alphal - alphar
   else
      g_i = reshape(tmp_gradl, [N]) - g0
      alphadiff = alphal - alpha0
   end if
   g_i = g_i/step*displmax*factor

   if (present(dipdir)) then
      if (doublesided) then
         dipdir(:, icol) = (dipl - dipr)/step*displmax*factor
      else
         dipdir(:, icol) = (dipl - dip0)/step*displmax*factor
      end if
   end if

   if (present(poldir)) then
      poldir(:, icol) = [alphadiff(1, 1), alphadiff(1, 2), alphadiff(2, 2), &
         & alphadiff(1, 3), alphadiff(2, 3), alphadiff(3, 3)]/step*displmax*factor
   end if
end subroutine gradient_derivs_point

end module xtb_type_calculator
