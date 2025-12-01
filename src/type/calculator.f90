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
   use xtb_mctc_math, only: crossProd
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
subroutine odlrhessian(self, env, mol0, chk0, list, step, hess) ! TODO: this needs to return displdir and g for more manipulations
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

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(scc_results) :: res
   real(wp), allocatable :: distmat(:, :), h0(:, :), displdir(:, :), g(:, :), tmp_grad(:, :), g0(:), x(:), xyz(:, :), gr(:, :), gl(:, :)
   real(wp) :: energy, sigma(3, 3), egap, dist, barycenter(3), inertia(3), ax(3, 3), cross(3), Imat0, query(1), displmax
   real(wp) :: identity3(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   logical :: linear
   integer :: N, i, j, k, Ntr, info, lwork
   
   ! ========== INITIALIZATION ==========
   ! NOTE: maybe this needs to go to numhess?
   N = 3*mol0%n
   allocate(distmat(N, N), h0(N, N), displdir(N, N), g(N, N))

   ! hessian initial guess
   do i = 1, N
      ! unity as initial guess for hessian
      ! TODO: do we have a better guess? Swart model hessian
      h0(i, i) = 1.0_wp
   end do

   ! calculate unperturbed gradient
   call self%singlepoint(env, mol0, chk0, -1, .true., energy, tmp_grad, sigma, egap, res)
   ! gradients need to be flattened since hessian is also "flat"
   g0 = reshape(tmp_grad, [N])

   ! setup distmat
   do i = 1, mol0%n
      do j = i, mol0%n
         dist = mol0%dist(i, j) ! substract vdw radii
         distmat(3*i-2:3*i, 3*j-2:3*j) = dist
         distmat(3*j-2:3*j, 3*i-2:3*i) = dist
      end do
   end do

   ! set up initial displdir with trans, rot, and totally symmetric vib mode first
   ! translational displacements
   do i = 1, N
      displdir(3*i-2, 1) = 1.0_wp/sqrt(real(mol0%n, wp))
      displdir(3*i-1, 2) = 1.0_wp/sqrt(real(mol0%n, wp))
      displdir(3*i, 3) = 1.0_wp/sqrt(real(mol0%n, wp))
   end do
   
   ! calculate inertial moment and axes
   barycenter = sum(mol0%xyz, dim=2) / real(mol0%n, wp)
   Imat0 = 0.0_wp
   do i = 1, mol0%n
      vec = mol0%xyz(:, i) - barycenter(:)
      Imat0 = Imat0 + matmul(vec, transpose(vec))
   end do
   ax = Imat0 * identity3
   do i = 1, 3
      do j = 1, 3
         do k = 1, mol0%n
            ax(i, j) = ax(i, j) - (mol0%xyz(i, k) - barycenter(i)) * (mol0%xyz(j, k) - barycenter(j))
         end do
      end do
   end do

   lwork = -1
   call dsyev('V', 'U', 3, ax, 3, inertia, query, lwork, info)
   lwork = int(query(1))
   allocate(aux(lwork))
   call dsyev('V', 'U', 3, ax, 3, inertia, aux, lwork, info)

   ! rotational displacements
   Ntr = 3 ! number of translational and rotational degrees of freedom
   do i = 1, 3
      if (inertia(i) < 1e-4_wp) cycle
      do j = 1, mol0%n
         crossProd(ax(:, i), mol0%xyz(:, j) - barycenter(:), cross)
         displdir(3*j-2:3*j, Ntr+1) = cross
      end do
      displdir(:, Ntr+1) = displdir(:, Ntr+1) / norm2(displdir(:, Ntr+1))
      Ntr = Ntr + 1
   end do

   ! totally symmetric vibrational displacement
   do i = 1, mol0%n
      displdir(3*i-2:3*i, Ntr+1) = mol0%xyz(:, i) - barycenter(:)
      displdir(3*i-2:3*i, Ntr+1) = displdir(3*i-2:3*i, Ntr+1) / norm2(displdir(3*i-2:3*i, Ntr+1))
   end do

   ! generate remaining displdirs based on distmat and dmax
   ! TODO: compute neighbor list
   ! TODO: populate displdir
   ! NOTE: orthonormalize displdir?

   g = 0.0_wp
   ! TODO: get gradients along rotational displacements
   
   ! ========== GRADIENT DERIVATIVES ==========
   do i = Ntr + 1, N
      displmax = maxval(abs(displdir(:, i)))
      ! TODO: what about double sided stuff?
      call mol%copy(mol0)
      call chk%copy(chk0)
      x = reshape(mol0%xyz, [N])
      x = x + step*displdir(:, i)/displmax
      mol%xyz = reshape(x, [3, mol0%n])
      call self%singlepoint(env, mol, chk, -1, .true., energy, tmp_grad, sigma, egap, res)
      g(:, i) = reshape(tmp_grad, [N])
      g(:, i) = (g(:, i) - g0(:))/step*displmax
   end do

   ! ========== FINAL HESSIAN ==========
   ! construct hessian from local hessian and odlr correction
   ! compute local hessian
   call gen_local_hessian(distmat, displdir, g, dmax, hess)

   ! compute low rank correction
   call lr_loop(g, hess, displdir, final_err)

   contains
   
   !> Main routine to recover local Hessian using Conjugate Gradient
   subroutine gen_local_hessian(distmat_local, displdir_local, g_local, dmax, hess_out)
      !> Distance matrix between atoms
      real(wp), intent(in) :: distmat_local(:, :)
      !> Displacement directions
      real(wp), intent(in) :: displdir_local(:, :)
      !> Gradient derivatives
      real(wp), intent(in) :: g_local(:, :)
      !> Maximum distance threshold
      real(wp), intent(in) :: dmax
      !> Output Hessian matrix
      real(wp), intent(out) :: hess_out(:, :)

      real(wp), parameter :: lam = 1.0e-2_wp, bet = 1.5_wp, ddmax = 5.0_wp
      
      ! Local work arrays
      real(wp), allocatable :: rhsv(:), hnumv(:)
      real(wp), allocatable :: W2(:, :), RHS(:, :), D_DT(:, :)
      logical, allocatable :: mask(:,:), mask_tril(:,:)
      integer :: l, m, ndim, ndispl = size(displdir_local, 1)
      
      ! 1. Calculate Regularization Term W2
      allocate(W2(N, N))
      do m = 1, N
         do l = 1, N
               W2(l,m) = lam * (max(0.0_wp, distmat_local(l,m) - dmax))**(2.0_wp * bet)
         end do
      end do

      ! 2. Calculate RHS = (g * displdir^T)
      allocate(RHS(N, N))
      ! Call BLAS DGEMM: RHS = 1.0 * g * displdir^T
      call dgemm('N', 'T', N, N, ndispl, 1.0_wp, g_local, N, displdir_local, N, 0.0_wp, RHS, N)
      
      ! Force Symmetry
      do m = 1, N
         do l = 1, N
               RHS(l,m) = 0.5_wp * (RHS(l,m) + RHS(m,l))
         end do
      end do

      ! 3. Precompute D * D^T for the operator
      allocate(D_DT(N, N))
      call dgemm('N', 'T', N, N, ndispl, 1.0_wp, displdir_local, N, displdir_local, N, 0.0_wp, D_DT, N)

      ! 4. Masks and Packing
      allocate(mask(N,N), mask_tril(N,N))
      mask = (distmat_local < (dmax + ddmax))
      forall(l=1:N, m=1:N) mask_tril(l,m) = mask(l,m) .and. (l >= m)
      
      ! RHS Vector (b in Ax=b)
      rhsv = pack(RHS, mask_tril)
      ndim = size(rhsv)
      
      allocate(hnumv(ndim))
      hnumv = rhsv ! Initial guess = RHS

      ! 5. Call Conjugate Gradient Solver
      call cg_solver(ndim, rhsv, hnumv, matvec_wrapper)

      ! 6. Recover Hessian from vector
      hess_out = unpack_sym(hnumv, mask_tril, N)

   end subroutine gen_local_hessian

   !> Corrects Hessian hnum using a symmetric, low-rank update
   !> so that g approx hnum * displdir
   subroutine lr_loop(g_in, hess_out, d_in, final_err)
      real(wp), intent(in) :: g_in(N, N)    ! Input Gradients
      real(wp), intent(inout) :: hess_out(N, N) ! Hessian to correct
      real(wp), intent(in) :: d_in(N, N)    ! Displacement directions
      real(wp), intent(out) :: final_err ! Final residual error

      real(wp), parameter :: mingrad_LR = 1.0e-3_wp
      real(wp), parameter :: thresh_LR = 1.0_e-8_wp
      integer, parameter :: maxiter_LR = 100

      ! Local variables
      real(wp), allocatable :: g_local(:, :), displdir_local(:, :)
      real(wp), allocatable :: resid(:, :), hcorr(:, :)
      real(wp) :: dampfac, err0, err, norm_g, norm_gi
      integer :: i, j, it
      
      ! BLAS Helper
      real(wp), external :: dnrm2

      allocate(g_local(N, N), displdir_local(N, N))
      allocate(resid(N, N), hcorr(N,N))

      ! 1. Weighting Step
      g_local = 0.0_wp
      displdir_local = 0.0_wp
      
      do i = 1, k
         ! Python: norm_gi = max(norm(g[:,i]), mingrad)
         norm_gi = max(dnrm2(n, g_in(:, i), 1), mingrad_LR)
         
         ! Scale columns
         g_local(:, i) = (mingrad_LR / norm_gi) * g_in(:, i)
         displdir_local(:, i) = (mingrad_LR / norm_gi) * d_in(:, i)
      end do

      dampfac = 1.0_wp
      err0 = huge(1.0_wp)
      
      ! Calculate Frobenius norm of entire matrix G (treated as vector of size n*k)
      norm_g = dnrm2(N*k, g_local, 1)

      ! 2. Iterative Correction Loop
      loop_lr: do it = 1, maxiter_LR
         
         resid = g_local ! Copy g to resid first
         
         ! Call DGEMM: C = alpha*A*B + beta*C
         ! resid = (-1.0) * hnum * d + (1.0) * resid
         call dgemm('N', 'N', N, k, N, -1.0_wp, hess_out, N, displdir_local, N, 1.0_wp, resid, N)
         
         err = dnrm2(N*N, resid, 1)
         
         if (err < thresh_LR) then
               ! Converged successfully
               exit loop_lr
               
         elseif (abs(err - err0) < thresh_LR * err0) then
               print *, 'Warning: Gradients cannot be reproduced by symmetric Hessian (Stagnation).'
               exit loop_lr
               
         elseif (err > err0 .and. err > norm_g) then
               ! Divergence detected
               dampfac = dampfac * 0.5_wp
               ! (Optional: Print warning if verbose)
               ! print *, 'Damping factor reduced to', dampfac
         end if
         
         ! (Optional: Print iteration status)
         ! print *, 'Iter', it, 'Error', err

         ! hcorr = resid * d^T
         ! DGEMM: C = alpha*A*B^T + beta*C
         call dgemm('N', 'T', N, N, k, 1.0_wp, resid, N, displdir_local, N, 0.0_wp, hcorr, N)
         
         ! Symmetrize hcorr and Apply Update to hnum
         ! hnum = hnum + dampfac * 0.5 * (hcorr + hcorr^T)
         do j = 1, N
               do i = 1, N
                  ! Average off-diagonals
                  real(wp) :: val
                  val = 0.5_wp * (hcorr(i, j) + hcorr(j, i))
                  hess_out(i, j) = hess_out(i, j) + (dampfac * val)
               end do
         end do
         
         err0 = err
         
      end do loop_lr

      final_err = err

   end subroutine lr_loop

end subroutine odlrhessian

!> Assumes A is Symmetric Positive Definite
subroutine cg_solver(n, b, x, matvec)
   integer, intent(in) :: n
   real(wp), intent(in) :: b(n)       ! RHS
   real(wp), intent(inout) :: x(n)    ! Solution (Input: Initial Guess)
   ! Procedure pointer for the matrix-vector product
   interface
      subroutine matvec(v_in, v_out)
            import :: wp
            real(wp), intent(in) :: v_in(:)
            real(wp), intent(out) :: v_out(:)
      end subroutine matvec
   end interface

   ! Work arrays
   real(wp), allocatable :: r(:), p(:), Ap(:)
   real(wp) :: alpha, beta, rsold, rsnew
   integer :: k
   real(wp), parameter :: tol = 1.0e-14_wp
   integer, parameter :: max_iter = 1000
   
   ! BLAS functions
   real(wp), external :: ddot, daxpy, dscal
   
   allocate(r(n), p(n), Ap(n))
   
   ! 1. r = b - A * x
   call matvec(x, Ap)
   r = b - Ap
   
   ! 2. p = r
   p = r
   
   ! 3. rsold = r' * r
   rsold = ddot(n, r, 1, r, 1)
   
   ! Check for immediate convergence
   if (sqrt(rsold) < tol) return

   ! 4. Iteration loop
   do k = 1, max_iter
      
      ! Ap = A * p
      call matvec(p, Ap)
      
      ! alpha = rsold / (p' * Ap)
      alpha = rsold / ddot(n, p, 1, Ap, 1)
      
      ! x = x + alpha * p
      call daxpy(n, alpha, p, 1, x, 1)
      
      ! r = r - alpha * Ap
      call daxpy(n, -alpha, Ap, 1, r, 1)
      
      ! rsnew = r' * r
      rsnew = ddot(n, r, 1, r, 1)
      
      if (sqrt(rsnew) < tol) exit
      
      ! beta = rsnew / rsold
      beta = rsnew / rsold
      
      ! p = r + beta * p
      ! We do p = beta*p + r by scaling p first then adding r
      call dscal(n, beta, p, 1)
      call daxpy(n, 1.0_wp, r, 1, p, 1)
      
      rsold = rsnew
   end do
   
   if (k > max_iter) then
      print *, "Warning: CG did not converge within max iterations"
   end if

end subroutine cg_solver

!> Helper to unpack vector to Symmetric Matrix
function unpack_sym(v, mask, n) result(H)
   real(wp), intent(in) :: v(:)
   logical, intent(in) :: mask(n,n)
   integer, intent(in) :: n
   real(wp) :: H(n,n)
   integer :: i, j
   
   H = 0.0_wp
   ! Unpack lower triangle based on mask
   H = unpack(v, mask, field=0.0_wp)
   
   ! Symmetrize (Copy lower to upper)
   do j = 1, n - 1
      do i = j+1, n
            H(j,i) = H(i,j)
      end do
   end do
end function unpack_sym

end module xtb_type_calculator
