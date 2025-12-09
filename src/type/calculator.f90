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

      !> Write informative printout
      procedure(writeInfo), deferred :: writeInfo

   end type TCalculator

   type :: adj_list
      integer, allocatable :: neighbors(:)
   end type adj_list

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
subroutine odlrhessian(self, env, mol0, chk0, list, step, displdir, g, hess)
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
   real(wp), intent(inout) :: displdir(:, :)
   !> Gradients
   real(wp), intent(inout) :: g(:, :)
   !> Array to add Hessian to
   real(wp), intent(inout) :: hess(:, :)
   
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
   real(wp), allocatable :: distmat(:, :), h0(:, :), h0v(:), tmp_grad(:, :), g0(:), x(:), xyz(:, :), gr(:, :), gl(:, :)
   real(wp) :: energy, sigma(3, 3), egap, dist, barycenter(3), inertia(3), ax(3, 3), cross(3), Imat0, query(1), displmax
   real(wp) :: identity3(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1],[3, 3]), final_err
   logical :: linear
   integer, allocatable :: nbcounts(:)
   integer :: N, i, j, k, Ntr, info, lwork, ndispl_final, max_nb
   
   ! ========== INITIALIZATION ==========
   ! NOTE: maybe this needs to go to numhess?
   N = 3 * mol0%n

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
   call gen_displdir(N, Ntr, h0, displdir, max_nb, neighborlist, nbcounts, eps, eps2, displdir, ndispl_final)


   g = 0.0_wp ! TODO: this should not be done in general since gradient derivs might be input
   
   ! ========== GRADIENT DERIVATIVES ==========
   write(env%unit, '(A)') "Calculating gradient derivatives"
   do i = Ntr + 1, ndispl_final
      displmax = maxval(abs(displdir(:, i)))
      ! TODO: what about double sided stuff?
      call mol%copy(mol0)
      call chk%copy(chk0)
      x = reshape(mol0%xyz,[N])
      x = x + step * displdir(:, i) / displmax
      mol%xyz = reshape(x,[3, mol0%n])
      call self%singlepoint(env, mol, chk, -1, .true., energy, tmp_grad, sigma, egap, res)
      g(:, i) = reshape(tmp_grad,[N])
      g(:, i) = (g(:, i) - g0(:)) / step * displmax
   end do

   ! ========== FINAL HESSIAN ==========
   ! construct hessian from local hessian and odlr correction
   ! compute local hessian
   write(env%unit, '(A)') "Computing local Hessian"
   call gen_local_hessian(distmat, displdir, g, dmax, hess)

   ! compute low rank correction
   write(env%unit, '(A)') "Computing low rank correction"
   call lr_loop(ndispl_final, g, hess, displdir, final_err)

end subroutine odlrhessian

!> Main routine to recover local Hessian
subroutine gen_local_hessian(distmat, displdir, g, dmax, hess_out)
   !> Distance matrix between atoms
   real(wp), intent(in) :: distmat(:, :)
   !> Displacement directions
   real(wp), intent(in) :: displdir(:, :)
   !> Gradient derivatives
   real(wp), intent(in) :: g(:, :)
   !> Maximum distance threshold
   real(wp), intent(in) :: dmax
   !> Output Hessian matrix
   real(wp), intent(out) :: hess_out(:, :)

   real(wp), parameter :: lam = 1.0e-2_wp, bet = 1.5_wp, ddmax = 5.0_wp
   
   ! Local work arrays
   integer, allocatable :: idx_map(:, :)
   real(wp), allocatable :: W2(:, :), rhs(:, :), rhsv(:), A(:, :), D(:, :)
   logical, allocatable :: mask(:, :), mask_ut(:, :)
   integer :: i, j, k, l, ndim, ndispl, N, idx_ij, idx_kl, info
   real(wp) :: Dij_dot_Dkl

   ndispl = size(displdir, 1)
   N = size(distmat, 1)

   ! Calculate Regularization Term W2
   allocate(W2(N, N))
   W2 = lam * max(0.0_wp, distmat(:, :) - dmax)**(2.0_wp * bet)

   ! Calculate rhs
   allocate(rhs(N, N))
   call dgemm('N', 'T', N, N, ndispl, 1.0_wp, g, N, displdir, N, 0.0_wp, rhs, N)
   rhs = 0.5_wp * (rhs + transpose(rhs))

   ! Masks and Packing
   allocate(mask(N, N), mask_ut(N, N))
   mask = (distmat < (dmax + ddmax))
   mask_ut = .false.
   do j = 2, N
      do i = 1, j - 1
         mask_ut(i, j) = mask(i, j)
      end do
   end do
   
   ! RHS Vector (b in Ax=b)
   rhsv = pack(rhs, mask_ut)
   ndim = size(rhsv)

   !  Build index mapping from (i,j) to packed index
   idx_map = 0
   k = 0
   do j = 2, n
      do i = 1, j - 1
            if (mask_ut(i, j)) then
               k = k + 1
               idx_map(i, j) = k
               idx_map(j, i) = k  ! Symmetric access
            end if
      end do
   end do

   allocate(D(N, N))
   call dgemm('N', 'T', N, N, ndispl, 1.0_wp, displdir, N, displdir, N, 0.0_wp, D, N)

   ! Compute A
   allocate(A(ndim, ndim))

   A = 0.0_wp
      
   do j = 2, n
      do i = 1, j - 1
         if (.not. mask_ut(i, j)) cycle
         idx_ij = idx_map(i, j)
         
         do l = 2, n
            do k = 1, l - 1
               if (.not. mask_ut(k, l)) cycle
               idx_kl = idx_map(k, l)
               
               !  Only compute upper triangle of A (it's symmetric)
               if (idx_kl < idx_ij) cycle
               
               !  A[(i,j),(k,l)] = D(i,k)*D(j,l) + D(i,l)*D(j,k)
               Dij_dot_Dkl = D(i,k)*D(j,l) + D(i,l)*D(j,k)
               
               A(idx_ij, idx_kl) = Dij_dot_Dkl
               
               !  Add regularization on diagonal
               if (idx_ij == idx_kl) then
                     A(idx_ij, idx_kl) = A(idx_ij, idx_kl) + W2(i, j)
               end if
               
               ! Fill lower triangle
               if (idx_kl > idx_ij) then
                     A(idx_kl, idx_ij) = A(idx_ij, idx_kl)
               end if
            end do
         end do
      end do
   end do

   ! Solve
   call dposv('U', ndim, 1, A, ndim, rhsv, ndim, info)

   ! Recover Hessian from vector
   hess_out = unpack_sym(rhsv, mask_ut, N)

   ! Fill lower triangle
   do j = 2, N
      do i = 1, j - 1
         hess_out(i, j) = hess_out(j, i)
      end do
   end do

   ! Free memory
   deallocate(idx_map, mask, mask_ut, A, rhsv, rhs)

end subroutine gen_local_hessian

!> Corrects Hessian hnum using a symmetric, low-rank update
!> so that g approx hnum * displdir
subroutine lr_loop(ndispl, g_in, hess_out, d_in, final_err)
   integer, intent(in) :: ndispl
   real(wp), intent(in) :: g_in(:, :)    ! Input Gradients
   real(wp), intent(inout) :: hess_out(:, :) ! Hessian to correct
   real(wp), intent(in) :: d_in(:, :)    ! Displacement directions
   real(wp), intent(out) :: final_err ! Final residual error

   real(wp), parameter :: mingrad_LR = 1.0e-3_wp
   real(wp), parameter :: thresh_LR = 1.0e-8_wp
   integer, parameter :: maxiter_LR = 100

   ! Local variables
   real(wp), allocatable :: g_local(:, :), displdir_local(:, :)
   real(wp), allocatable :: resid(:, :), hcorr(:, :)
   real(wp) :: dampfac, err0, err, norm_g, norm_gi, val
   integer :: i, j, it, N
   
   ! BLAS Helper
   real(wp), external :: dnrm2

   N = size(g_in, 1)

   allocate(g_local(N, N), displdir_local(N, N))
   allocate(resid(N, N), hcorr(N, N))

   ! 1. Weighting Step
   g_local = 0.0_wp
   displdir_local = 0.0_wp
   
   do i = 1, ndispl
      ! Python: norm_gi = max(norm(g[:,i]), mingrad)
      norm_gi = max(dnrm2(n, g_in(:, i), 1), mingrad_LR)
      
      ! Scale columns
      g_local(:, i) = (mingrad_LR / norm_gi) * g_in(:, i)
      displdir_local(:, i) = (mingrad_LR / norm_gi) * d_in(:, i)
   end do

   dampfac = 1.0_wp
   err0 = huge(1.0_wp)
   
   ! Calculate Frobenius norm of entire matrix G (treated as vector of size n*ndispl)
   norm_g = dnrm2(N * ndispl, g_local, 1)

   ! 2. Iterative Correction Loop
   loop_lr: do it = 1, maxiter_LR
      
      resid = g_local ! Copy g to resid first
      
      ! Call DGEMM: C = alpha*A*B + beta*C
      ! resid = (-1.0) * hnum * d + (1.0) * resid
      call dgemm('N', 'N', N, ndispl, N, -1.0_wp, hess_out, N, displdir_local, N, 1.0_wp, resid, N)
      
      err = dnrm2(N * N, resid, 1)
      
      if (err < thresh_LR) then
            ! Converged successfully
            exit loop_lr
            
      else if (abs(err - err0) < thresh_LR * err0) then
            print *, 'Warning: Gradients cannot be reproduced by symmetric Hessian (Stagnation).'
            exit loop_lr
            
      else if (err > err0 .and. err > norm_g) then
            ! Divergence detected
            dampfac = dampfac * 0.5_wp
            ! (Optional: Print warning if verbose)
            ! print *, 'Damping factor reduced to', dampfac
      end if
      
      ! (Optional: Print iteration status)
      ! print *, 'Iter', it, 'Error', err

      ! hcorr = resid * d^T
      ! DGEMM: C = alpha*A*B^T + beta*C
      call dgemm('N', 'T', N, N, ndispl, 1.0_wp, resid, N, displdir_local, N, 0.0_wp, hcorr, N)
      
      ! Symmetrize hcorr and Apply Update to hnum
      ! hnum = hnum + dampfac * 0.5 * (hcorr + hcorr^T)
      do j = 1, N
            do i = 1, N
               ! Average off-diagonals
               val = 0.5_wp * (hcorr(i, j) + hcorr(j, i))
               hess_out(i, j) = hess_out(i, j) + (dampfac * val)
            end do
      end do
      
      err0 = err
      
   end do loop_lr

   final_err = err

end subroutine lr_loop

!> Helper to unpack vector to Symmetric Matrix
function unpack_sym(v, mask, n) result(H)
   real(wp), intent(in) :: v(:)
   logical, intent(in) :: mask(n, n)
   integer, intent(in) :: n
   real(wp) :: H(n, n)
   integer :: i, j
   
   H = 0.0_wp
   ! Unpack lower triangle based on mask
   H = unpack(v, mask, field=0.0_wp)
   
   ! Symmetrize (Copy lower to upper)
   do j = 1, n - 1
      do i = j + 1, n
            H(j, i) = H(i, j)
      end do
   end do
end function unpack_sym

subroutine get_neighbor_list(distmat, dmax, nblist)
   real(wp), intent(in) :: distmat(:, :)
   real(wp), intent(in) :: dmax
   type(adj_list), allocatable, intent(out) :: nblist(:)

   integer, allocatable :: labels(:)
   real(wp), allocatable :: comp_dist(:, :)
   integer, allocatable :: mst_matrix(:, :)
   real(wp) :: d, min_d
   real(wp), parameter :: eps = 1.0e-8_wp
   integer :: i, j, ncomp, N

   N = size(distmat, 1)
   allocate(nblist(N))

   ! 1. Calculate Initial Neighbors
   do i = 1, N - 1
      do j = i + 1, N
            if (d < dmax) then
               call add_neighbor(nblist(i), j)
               call add_neighbor(nblist(j), i)
            end if
      end do
   end do

   ! 2. Identify Connected Components (DFS)
   allocate(labels(N))
   labels = 0
   ncomp = 0
   
   do i = 1, N
      if (labels(i) == 0) then
            ncomp = ncomp + 1
            call dfs_label(i, N, nblist, labels, ncomp)
      end if
   end do

   if (ncomp == 1) return ! Graph is already connected

   ! 3. Distance between components
   !    For every pair of components, find the minimum distance
   allocate(comp_dist(ncomp, ncomp))
   comp_dist = huge(1.0_wp)

   do i = 1, N - 1
      do j = i + 1, N
            if (labels(i) /= labels(j)) then
               d = distmat(i, j)
               if (d < comp_dist(labels(i), labels(j))) then
                  comp_dist(labels(i), labels(j)) = d
                  comp_dist(labels(j), labels(i)) = d
               end if
            end if
      end do
   end do

   ! 4. Minimum Spanning Tree (Prim's Algorithm) on Components
   !    Returns symmetric matrix: 1 if connected in MST, 0 otherwise
   call prim_mst(ncomp, comp_dist, mst_matrix)

   ! 5. Stitching: Add necessary links closer than MST distance + eps
   do i = 1, N - 1
      do j = i + 1, N
            if (labels(i) /= labels(j)) then
               ! If these components are connected in the MST
               if (mst_matrix(labels(i), labels(j)) == 1) then
                  ! Get the min distance required to bridge them
                  min_d = comp_dist(labels(i), labels(j))
                  
                  ! If this pair provides that bridge (handling degeneracy)
                  if (distmat(i, j) <= min_d + eps) then
                        call add_neighbor_unique(nblist(i), j)
                        call add_neighbor_unique(nblist(j), i)
                  end if
               end if
            end if
      end do
   end do
end subroutine get_neighbor_list

! --- Helper: Add to dynamic array ---
subroutine add_neighbor(list, val)
   type(adj_list), intent(inout) :: list
   integer, intent(in) :: val
   integer, allocatable :: tmp(:)
   integer :: sz

   if (.not. allocated(list%neighbors)) then
      allocate(list%neighbors(1))
      list%neighbors(1) = val
   else
      sz = size(list%neighbors)
      allocate(tmp(sz + 1))
      tmp(1:sz) = list%neighbors
      tmp(sz + 1) = val
      call move_alloc(tmp, list%neighbors)
   end if
end subroutine add_neighbor

! --- Helper: Add only if not present ---
subroutine add_neighbor_unique(list, val)
   type(adj_list), intent(inout) :: list
   integer, intent(in) :: val
   integer :: k
   
   if (allocated(list%neighbors)) then
      do k = 1, size(list%neighbors)
            if (list%neighbors(k) == val) return
      end do
   end if
   call add_neighbor(list, val)
end subroutine add_neighbor_unique

! --- Helper: Recursive DFS for labeling ---
recursive subroutine dfs_label(u, n, nblist, labels, comp_id)
   integer, intent(in) :: u, n, comp_id
   type(adj_list), intent(in) :: nblist(:)
   integer, intent(inout) :: labels(:)
   integer :: k, v

   labels(u) = comp_id
   if (.not. allocated(nblist(u)%neighbors)) return

   do k = 1, size(nblist(u)%neighbors)
      v = nblist(u)%neighbors(k)
      if (labels(v) == 0) then
            call dfs_label(v, n, nblist, labels, comp_id)
      end if
   end do
end subroutine dfs_label

! --- Helper: Prim's Algorithm for MST ---
subroutine prim_mst(nc, dists, adj_mst)
   integer, intent(in) :: nc
   real(8), intent(in) :: dists(nc, nc)
   integer, allocatable, intent(out) :: adj_mst(:, :)
   
   real(8) :: min_val, key(nc)
   integer :: parent(nc)
   logical :: mst_set(nc)
   integer :: i, count, u, v

   allocate(adj_mst(nc, nc))
   adj_mst = 0
   
   key = huge(1.0d0)
   parent = 0
   mst_set = .false.
   
   key(1) = 0.0d0
   parent(1) = -1

   do count = 1, nc - 1
      ! Pick minimum key vertex not yet including in MST
      min_val = huge(1.0d0)
      u = -1
      do i = 1, nc
            if (.not. mst_set(i) .and. key(i) < min_val) then
               min_val = key(i)
               u = i
            end if
      end do
      
      if (u == -1) exit
      mst_set(u) = .true.

      ! Update adjacent vertices
      do v = 1, nc
            if (dists(u, v) > 0.0d0 .and. .not. mst_set(v) .and. dists(u, v) < key(v)) then
               parent(v) = u
               key(v) = dists(u, v)
            end if
      end do
   end do

   ! Convert parent array to adjacency matrix
   do i = 2, nc
      if (parent(i) /= -1) then
            adj_mst(i, parent(i)) = 1
            adj_mst(parent(i), i) = 1
      end if
   end do
end subroutine prim_mst

subroutine gen_displdir(n, ndispl0, h0, displdir0, max_nb, nblist, nbcounts, &
                        eps, eps2, displdir, ndispl_final)
   integer, intent(in) :: n, ndispl0, max_nb
   real(wp), intent(in) :: h0(n,n)
   real(wp), intent(in) :: displdir0(n, ndispl0)
   type(adj_list), intent(in) :: nblist(:)
   integer, intent(in) :: nbcounts(n)           ! Actual number of neighbors per atom
   real(wp), intent(in) :: eps, eps2
   
   real(wp), intent(out) :: displdir(n,n)
   integer, intent(out) :: ndispl_final

   ! Local variables
   integer :: i, j, k, p, q, nnb, info, n_curr, idx, local_max_ind
   integer :: nb_idx(max_nb)
   real(wp) :: ev(n), coverage(n), locev(max_nb), submat(max_nb, max_nb)
   real(wp) :: projmat(max_nb, max_nb), eye(max_nb, max_nb)
   real(wp) :: vec_subset(max_nb, n), U(max_nb, max_nb), VT(max_nb, max_nb)
   real(wp) :: S(max_nb), loceigs(max_nb)
   real(wp) :: ev1, ev2, norm1, norm2, v_norm, d_dot
   integer, allocatable :: iwork(:)
   real(wp), allocatable :: work(:)
   real(wp) :: norm_locev_max
   logical :: early_break

   ! Initialize
   displdir = 0.0_wp
   displdir(:, 1:ndispl0) = displdir0
   
   ! Workspace for LAPACK (allocate generously)
   allocate(work(10*max_nb + 10*n)) 
   allocate(iwork(8*max_nb))

   early_break = .true.
   ndispl_final = ndispl0

   ! --- Outer Loop: Generate new directions ---
   ! Corresponds to Python: for i in range(N-Ndispl0)
   do i = 1, n - ndispl0
      
      n_curr = ndispl0 + i - 1 ! Number of existing vectors
      ev = 0.0_wp
      coverage = 0.0_wp

      ! --- Inner Loop: Iterate over atoms/DoFs ---
      do j = 1, n
            nnb = nbcounts(j)
            nb_idx(1:nnb) = nblist(j)%neighbors(1:nnb)

            ! Skip if subspace saturated (heuristic from Python code)
            if (nnb <= n_curr) cycle

            ! 1. Extract submatrix H0 (submat)
            do p = 1, nnb
               do q = 1, nnb
                  submat(q, p) = h0(nb_idx(q), nb_idx(p))
               end do
            end do

            ! 2. Local Projection (orth replacement)
            ! Form matrix A = displdir[neighbors, 0:n_curr]
            do p = 1, n_curr
               do q = 1, nnb
                  vec_subset(q, p) = displdir(nb_idx(q), p)
               end do
            end do

            ! Perform SVD on vec_subset to find basis: A = U * S * VT
            ! We define rank based on S > eps
            if (n_curr > 0) then
               call dgesdd('S', nnb, n_curr, vec_subset, max_nb, S, U, max_nb, &
                           VT, max_nb, work, -1, iwork, info) ! Query size
               call dgesdd('S', nnb, n_curr, vec_subset, max_nb, S, U, max_nb, &
                           VT, max_nb, work, int(work(1)), iwork, info)
            else
               U = 0.0_wp ! No existing vectors
               S = 0.0_wp
            end if

            ! Construct Projector: P = I - sum(u*u.T) for significant singular values
            eye = 0.0_wp
            do p = 1, nnb; eye(p,p) = 1.0_wp; end do
            
            if (n_curr > 0) then
               do k = 1, min(nnb, n_curr)
                  if (S(k) > 1.0d-13) then ! Numerical threshold for rank
                        do p = 1, nnb
                           do q = 1, nnb
                              eye(q, p) = eye(q, p) - U(q, k) * U(p, k)
                           end do
                        end do
                  end if
               end do
            end if
            
            ! submat = P * submat * P.T
            ! Step A: temp = P * submat
            call dgemm('N', 'N', nnb, nnb, nnb, 1.0_wp, eye, max_nb, submat, max_nb, 0.0_wp, projmat, max_nb)
            ! Step B: submat = temp * P.T
            call dgemm('N', 'T', nnb, nnb, nnb, 1.0_wp, projmat, max_nb, eye, max_nb, 0.0_wp, submat, max_nb)
            
            ! Symmetrize
            do p = 1, nnb
               do q = 1, nnb
                     submat(q,p) = 0.5_wp * (submat(q,p) + submat(p,q))
               end do
            end do

            ! 3. Diagonalization (Eigen decomposition)
            ! dsyev: computes eigenvalues and eigenvectors
            call dsyev('V', 'U', nnb, submat, max_nb, loceigs, work, 10*max_nb, info)
            
            ! Eigenvectors stored in submat columns. Last one is max due to ascending sort in dsyev.
            locev(1:nnb) = submat(1:nnb, nnb)

            ! 4. Patching / Phase fixing
            local_max_ind = 1
            norm_locev_max = -1.0_wp
            
            do k = 1, nnb
               idx = nb_idx(k)
               
               ! Calc candidates
               ev1 = (coverage(idx) * ev(idx) + locev(k)) / (coverage(idx) + 1.0_wp)
               ev2 = (coverage(idx) * ev(idx) - locev(k)) / (coverage(idx) + 1.0_wp)
               
               ! We calculate norms squared for decision
               norm1 = 0.0_wp
               norm2 = 0.0_wp
               do p = 1, nnb
                  idx = nb_idx(p)
                  if (p == k) then ! This index helps us find max element for later
                        if (abs(locev(p)) > norm_locev_max) then
                           norm_locev_max = abs(locev(p))
                           local_max_ind = p
                        end if
                  end if
                  
                  ! Recompute tentative slices for norm
                  norm1 = norm1 + ((coverage(idx)*ev(idx) + locev(p))/(coverage(idx)+1.0_wp))**2
                  norm2 = norm2 + ((coverage(idx)*ev(idx) - locev(p))/(coverage(idx)+1.0_wp))**2
               end do
               norm1 = sqrt(norm1)
               norm2 = sqrt(norm2)
               
               ! Apply decision to ALL neighbors at once (loop break/structure needed)
               exit ! We have the norms for the patch, stop k loop, update all p
            end do
            
            ! Apply update
            if (norm1 > norm2 + eps) then
               do p = 1, nnb
                     idx = nb_idx(p)
                     ev(idx) = (coverage(idx)*ev(idx) + locev(p))/(coverage(idx)+1.0_wp)
               end do
            else if (norm1 < norm2 - eps) then
               do p = 1, nnb
                     idx = nb_idx(p)
                     ev(idx) = (coverage(idx)*ev(idx) - locev(p))/(coverage(idx)+1.0_wp)
               end do
            else
               ! Deterministic sign fix based on max element
               if (locev(local_max_ind) > 0.0_wp) then
                  do p = 1, nnb
                        idx = nb_idx(p)
                        ev(idx) = (coverage(idx)*ev(idx) + locev(p))/(coverage(idx)+1.0_wp)
                  end do
               else
                  do p = 1, nnb
                        idx = nb_idx(p)
                        ev(idx) = (coverage(idx)*ev(idx) - locev(p))/(coverage(idx)+1.0_wp)
                  end do
               end if
            end if

            ! Update coverage
            do p = 1, nnb
               coverage(nb_idx(p)) = coverage(nb_idx(p)) + 1.0_wp
            end do
      end do ! End J loop

      ! --- Gram-Schmidt Orthogonalization ---
      ! Project out previous columns from global ev
      do k = 1, n_curr
            ! d = dot(ev, displdir(:,k))
            d_dot = dot_product(ev, displdir(:, k))
            ! ev = ev - d * displdir(:,k)
            ev = ev - d_dot * displdir(:, k)
      end do

      ! --- Check Norm ---
      v_norm = sqrt(dot_product(ev, ev))
      
      if (v_norm < eps2) then
            early_break = .true.
            exit ! Break out of i loop
      else
            early_break = .false.
      end if

      ! Normalize and store
      ev = ev / v_norm
      displdir(:, n_curr + 1) = ev

      ndispl_final = n_curr + 1
   end do

end subroutine gen_displdir

end module xtb_type_calculator
