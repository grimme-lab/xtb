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
         crossprod(ax(:, i), mol0%xyz(:, j) - barycenter(:), cross)
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
   ! TODO: conjugate gradient solver

   ! compute low rank correction

   contains
   
   !> Main routine to recover local Hessian using Conjugate Gradient
   subroutine gen_local_hessian(n, ndispl, distmat, displdir, g, dmax, ddmax, lam, bet, hess_out)
      real(dp), allocatable :: rhsv(:), hnumv(:)
      logical, allocatable :: mask(:,:), mask_tril(:,:)
      integer :: i, j, ndim
      
      ! 1. Calculate Regularization Term W2
      allocate(W2(n,n))
      do j = 1, n
         do i = 1, n
               W2(i,j) = lam * (max(0.0_dp, distmat(i,j) - dmax))**(2.0_dp * bet)
         end do
      end do

      ! 2. Calculate RHS = (g * displdir.T). Symmetrized.
      allocate(RHS(n,n))
      ! Call BLAS DGEMM: RHS = 1.0 * g * displdir^T
      call dgemm('N', 'T', n, n, ndispl, 1.0_dp, g, n, displdir, n, 0.0_dp, RHS, n)
      
      ! Force Symmetry
      do j = 1, n
         do i = 1, n
               RHS(i,j) = 0.5_dp * (RHS(i,j) + RHS(j,i))
         end do
      end do

      ! 3. Precompute D * D^T for the operator
      allocate(D_DT(n,n))
      call dgemm('N', 'T', n, n, ndispl, 1.0_dp, displdir, n, displdir, n, 0.0_dp, D_DT, n)

      ! 4. Masks and Packing
      allocate(mask(n,n), mask_tril(n,n))
      mask = (distmat < (dmax + ddmax))
      forall(i=1:n, j=1:n) mask_tril(i,j) = mask(i,j) .and. (i >= j)
      
      ! RHS Vector (b in Ax=b)
      rhsv = pack(RHS, mask_tril)
      ndim = size(rhsv)
      
      allocate(hnumv(ndim))
      hnumv = rhsv ! Initial guess = RHS (common heuristic)

      ! 5. Call Conjugate Gradient Solver
      call cg_solver(ndim, rhsv, hnumv, matvec_wrapper)

      ! 6. Recover Hessian from vector
      hess_out = unpack_sym(hnumv, mask_tril, n)

   contains

   !> Assumes A is Symmetric Positive Definite
   subroutine cg_solver(n, b, x, matvec)
      integer, intent(in) :: n
      real(dp), intent(in) :: b(n)       ! RHS
      real(dp), intent(inout) :: x(n)    ! Solution (Input: Initial Guess)
      ! Procedure pointer for the matrix-vector product
      interface
         subroutine matvec(v_in, v_out)
               import :: dp
               real(dp), intent(in) :: v_in(:)
               real(dp), intent(out) :: v_out(:)
         end subroutine matvec
      end interface

      ! Work arrays
      real(dp), allocatable :: r(:), p(:), Ap(:)
      real(dp) :: alpha, beta, rsold, rsnew
      integer :: k
      real(dp), parameter :: tol = 1.0e-14_dp
      integer, parameter :: max_iter = 1000
      
      ! BLAS functions
      real(dp), external :: ddot
      
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
         call daxpy(n, 1.0_dp, r, 1, p, 1)
         
         rsold = rsnew
      end do
      
      if (k > max_iter) then
         print *, "Warning: CG did not converge within max iterations"
      end if

   end subroutine cg_solver

   !> Helper to unpack vector to Symmetric Matrix
   function unpack_sym(v, mask, n) result(H)
      real(dp), intent(in) :: v(:)
      logical, intent(in) :: mask(n,n)
      integer, intent(in) :: n
      real(dp) :: H(n,n)
      integer :: i, j
      
      H = 0.0_dp
      ! Unpack lower triangle based on mask
      H = unpack(v, mask, field=0.0_dp)
      
      ! Symmetrize (Copy lower to upper)
      do j = 1, n
         do i = j+1, n
               H(j,i) = H(i,j)
         end do
      end do
   end function unpack_sym

   !> Simple BLAS wrapper for daxpy (y = alpha*x + y)
   !> Note: Modern Fortran compilers link BLAS automatically or via -lblas
   subroutine daxpy(n, alpha, x, incx, y, incy)
      use iso_c_binding
      integer, intent(in) :: n, incx, incy
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: x(*)
      real(dp), intent(inout) :: y(*)
      ! This interface is usually provided by the external library, 
      ! explicitly defined here just in case specific linking is needed.
      ! In practice, remove this body and use 'external'.
   end subroutine daxpy

   !> BLAS scaler (x = alpha*x)
   subroutine dscal(n, alpha, x, incx)
      integer, intent(in) :: n, incx
      real(dp), intent(in) :: alpha
      real(dp), intent(inout) :: x(*)
   end subroutine dscal

   end subroutine gen_local_hessian

end subroutine odlrhessian

end module xtb_type_calculator
