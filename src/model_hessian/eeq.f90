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

!> EEQ model Hessian contribution
module xtb_modelhessian_eeq
   use xtb_bmatrix, only : bmat_accum_pairblock_packed
   use xtb_mctc_blas_wrap3, only : mctc_gemm
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_mctc_lapack_trf, only : lapack_sytrf
   use xtb_mctc_lapack_tri, only : lapack_sytri
   use xtb_type_environment, only : TEnvironment
   use xtb_type_param, only : chrg_parameter
   implicit none(type, external)

   interface
      subroutine dsysv(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
         import :: wp
         implicit none(type, external)
         character(len = 1), intent(in) :: uplo
         integer, intent(in) :: n, nrhs, lda, ldb, lwork
         ! allow(C071)
         real(wp), intent(inout) :: a(lda, *)
         ! allow(C071)
         integer, intent(out) :: ipiv(*)
         ! allow(C071)
         real(wp), intent(inout) :: b(ldb, *)
         ! allow(C071)
         real(wp), intent(inout) :: work(*)
         integer, intent(out) :: info
      end subroutine dsysv
   end interface

   public :: add_eeq_hessian
   private

contains

!> Add the EEQ contribution to an existing packed Hessian
subroutine add_eeq_hessian(env, n, at, xyz, chrg, chrgeq, kq, hess)
   !> Calculation environment
   type(TEnvironment), intent(inout) :: env
   !> Number of atoms
   integer, intent(in) :: n
   !> Atomic numbers
   integer, intent(in) :: at(n)
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Total molecular charge
   real(wp), intent(in) :: chrg
   !> Charge-model parameters
   type(chrg_parameter), intent(in) :: chrgeq
   !> EEQ Hessian scaling factor
   real(wp), intent(in) :: kq
   !> Packed Hessian updated in place
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)

   character(len=*), parameter :: source = "model_hessian_eeq"

   !  √π
   real(wp), parameter :: sqrtpi = sqrt(pi)
   !  √(2/π)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)

   ! Charge model
   integer :: m ! dimension of the Lagrangian
   real(wp), allocatable :: Amat(:, :)
   real(wp), allocatable :: Xvec(:)
   real(wp), allocatable :: Ainv(:, :)
   real(wp), allocatable :: dAmat(:, :, :)
   real(wp), allocatable :: dqdr(:, :, :)

   ! Local variables
   integer :: i, j, k
   real(wp) :: r, rij(3), r2
   real(wp) :: gamij, gamij2
   real(wp) :: arg, arg2, dtmp
   real(wp) :: expterm, erfterm
   real(wp) :: rxr(3, 3)

   ! Scratch variables
   real(wp), allocatable :: alpha(:)
   real(wp), allocatable :: xtmp(:)
   real(wp), allocatable :: atmp(:, :)

   ! LAPACK work variables
   integer, allocatable :: ipiv(:)
   real(wp), allocatable :: work(:)
   integer :: lwork
   integer :: info
   real(wp) :: test(1)

   ! Initialization
   m = n + 1
   allocate( ipiv(m), source = 0 )
   allocate( Amat(m, m), Xvec(m), alpha(n), dqdr(3, n, m), source = 0.0_wp )

   ! Set up the A matrix and X vector
   !  αi -> alpha(i), ENi -> xi(i), κi -> kappa(i), Jii -> gam(i)
   !  γij = 1/√(αi+αj)
   !  Xi  = -ENi + κi·√CNi
   !  Aii = Jii + 2/√π·γii
   !  Aij = erf(γij·Rij)/Rij = 2/√π·F0(γ²ij·R²ij)
   !  prepare some arrays
   !$omp parallel default(none) &
   !$omp shared(n,at,chrgeq) &
   !$omp private(i) &
   !$omp shared(Xvec,alpha)
   !$omp do schedule(dynamic)
   do i = 1, n
      Xvec(i) = -chrgeq%en(i)
      alpha(i) = chrgeq%alpha(i)**2
   end do
   !$omp enddo
   !$omp endparallel

   !$omp parallel default(none) &
   !$omp shared(n,at,xyz,chrgeq,alpha) &
   !$omp private(i,j,r,gamij) &
   !$omp shared(Amat)
   !$omp do schedule(dynamic)
   ! prepare A matrix
   do i = 1, n
      ! EN of atom i
      do j = 1, i - 1
         r = sqrt(sum((xyz(:, j) - xyz(:, i))**2))
         gamij = 1.0_wp / sqrt(alpha(i) + alpha(j))
         Amat(j, i) = erf(gamij*r) / r
         Amat(i, j) = Amat(j, i)
      end do
      Amat(i, i) = chrgeq%gam(i) + sqrt2pi / sqrt(alpha(i))
   end do
   !$omp enddo
   !$omp endparallel

   ! Solve the linear equations to obtain partial charges
   Amat(m, 1:m) = 1.0_wp
   Amat(1:m, m) = 1.0_wp
   Amat(m, m  ) = 0.0_wp
   Xvec(m) = chrg
   ! generate temporary copy
   allocate( Atmp(m, m), source = Amat )
   allocate( Xtmp(m), source = Xvec )

   ! assume work space query, set best value to test after first dsysv call
   call dsysv("u", m, 1, Atmp, m, ipiv, Xtmp, m, test, -1, info)
   if (info /= 0) then
      call env%error("DSYSV workspace query failed for the EEQ matrix", source)
      return
   end if
   lwork = max(1, int(test(1)))
   allocate( work(lwork), source = 0.0_wp )

   call dsysv("u", m, 1, Atmp, m, ipiv, Xtmp, m, work, lwork, info)
   if (info /= 0) then
      call env%error("DSYSV failed while solving EEQ charges", source)
      return
   end if

   if (abs(sum(Xtmp(:n)) - chrg) > 1.0e-6_wp) then
      call env%error("EEQ charge constraint was not satisfied", source)
      return
   end if
   !print'(3f20.14)',Xtmp

   ! Calculate isotropic electrostatic (IES) energy
   !  E = ∑i (ENi - κi·√CNi)·qi + ∑i (Jii + 2/√π·γii)·q²i
   !      + ½ ∑i ∑j,j≠i qi·qj·2/√π·F0(γ²ij·R²ij)
   !    = q·(½A·q - X)
   !   work(:m) = Xvec
   !   call dsymv('u',m,0.5_wp,Amat,m,Xtmp,1,-1.0_wp,work,1)
   !   es = dot_product(Xtmp,work(:m))
   !   energy = es + energy

   ! Calculate molecular gradient of the IES energy
   !  dE/dRj -> g(:,j), ∂Xi/∂Rj -> -dcn(:,i,j), ½∂Aij/∂Rj -> dAmat(:,j,i)
   !  dE/dR = (½∂A/∂R·q - ∂X/∂R)·q
   !  ∂Aij/∂Rj = ∂Aij/∂Ri
   allocate( dAmat(3, n, m), source = 0.0_wp )
   !$omp parallel default(none) &
   !$omp shared(n,xyz,alpha,Amat,Xtmp) &
   !$omp private(i,j,rij,r2,gamij,arg,dtmp) &
   !$omp reduction(+:dAmat)
   !$omp do schedule(dynamic)
   do i = 1, n
      do j = 1, i - 1
         rij = xyz(:, i) - xyz(:, j)
         r2 = sum(rij**2)
         gamij = 1.0_wp / sqrt(alpha(i) + alpha(j))
         arg = gamij**2 * r2
         dtmp = 2.0_wp * gamij * exp(-arg) / (sqrtpi*r2) - Amat(j, i) / r2
         dAmat(:, i, i) = +dtmp * rij * Xtmp(j) + dAmat(:, i, i)
         dAmat(:, j, j) = -dtmp * rij * Xtmp(i) + dAmat(:, j, j)
         dAmat(:, i, j) = +dtmp * rij * Xtmp(i)
         dAmat(:, j, i) = -dtmp * rij * Xtmp(j)
      end do
   end do
   !$omp enddo
   !$omp endparallel

   ! Invert the A matrix using a Bunch-Kaufman factorization
!  A⁻¹ = (L·D·L^T)⁻¹ = L^T·D⁻¹·L
   allocate( Ainv(m, m), source = Amat )

   ! assume work space query, set best value to test after first dsytrf call
   call lapack_sytrf("L", m, Ainv, m, ipiv, test, -1, info)
   if (info /= 0) then
      call env%error("DSYTRF workspace query failed for the EEQ matrix", source)
      return
   end if
   if (int(test(1)) > lwork) then
      deallocate(work)
      lwork = int(test(1))
      allocate( work(lwork), source = 0.0_wp )
   end if

   ! Bunch-Kaufman factorization A = L*D*L**T
   call lapack_sytrf("L", m, Ainv, m, ipiv, work, lwork, info)
   if (info /= 0) then
      call env%error("DSYTRF failed while factorizing the EEQ matrix", source)
      return
   end if

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Ainv matrix
   ! Ainv matrix is overwritten with lower triangular part of A⁻¹
   call lapack_sytri("L", m, Ainv, m, ipiv, work, info)
   if (info /= 0) then
      call env%error("DSYTRI failed while inverting the EEQ matrix", source)
      return
   end if

   ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
   do i = 1, m
      do j = i + 1, m
         Ainv(i, j) = Ainv(j, i)
      end do
   end do

   ! Calculate gradient of the partial charge with respect to nuclear coordinates
   !call dsymm('r','l',3*n,m,-1.0_wp,Ainv,m,dAmat,3*n,1.0_wp,dqdr,3*n)
   call mctc_gemm(dAmat, Ainv, dqdr, alpha = -1.0_wp, beta = 1.0_wp)
   !print'(/,"analytical gradient")'
   !print'(3f20.14)',dqdr(:,:,:n)

   ! Molecular Hessian calculation
   do i = 1, n
      do j = 1, i - 1
         rij = xyz(:, j) - xyz(:, i)
         r2 = sum(rij**2)
         r = sqrt(r2)
         gamij = 1.0_wp / sqrt(alpha(i) + alpha(j))
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         erfterm = Xtmp(i) * Xtmp(j) * erf(arg) / r
         expterm = Xtmp(i) * Xtmp(j) * 2 * gamij * exp(-arg2) / sqrtpi
         ! ∂²(qAq)/(∂Ri∂Rj):
         ! ∂²(qAq)/(∂Xi∂Xi) = (1-3X²ij/R²ij-2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         !                  - (R²ij-3X²ij) erf[γij·Rij]/R⁵ij
         ! ∂²(qAq)/(∂Xi∂Xj) = (R²ij-3X²ij) erf[γij·Rij]/R⁵ij
         !                  - (1-3X²ij/R²ij-2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         ! ∂²(qAq)/(∂Xi∂Yi) = 3X²ij erf[γij·Rij]/R⁵ij
         !                  - (3X²ij/R²ij+2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         ! ∂²(qAq)/(∂Xi∂Yj) = (3X²ij/R²ij+2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         !                  - 3X²ij erf[γij·Rij]/R⁵ij
         rxr(1, 1) = erfterm * ( 3*rij(1)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(1)**2/r2**2 + 2*gamij2*rij(1)**2/r2 - 1/r2 )
         rxr(2, 2) = erfterm * ( 3*rij(2)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(2)**2/r2**2 + 2*gamij2*rij(2)**2/r2 - 1/r2 )
         rxr(3, 3) = erfterm * ( 3*rij(3)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(3)**2/r2**2 + 2*gamij2*rij(3)**2/r2 - 1/r2 )
         rxr(2, 1) = erfterm * 3 * rij(2) * rij(1) / r2**2 &
                  - expterm * ( 3*rij(2)*rij(1)/r2**2 + 2*gamij2*rij(2)*rij(1)/r2 )
         rxr(3, 1) = erfterm * 3 * rij(3) * rij(1) / r2**2 &
                  - expterm * ( 3*rij(3)*rij(1)/r2**2 + 2*gamij2*rij(3)*rij(1)/r2 )
         rxr(3, 2) = erfterm * 3 * rij(3) * rij(2) / r2**2 &
                  - expterm * ( 3*rij(3)*rij(2)/r2**2 + 2*gamij2*rij(3)*rij(2)/r2 )

         do k = 1, m
            rxr(1, 1) = rxr(1, 1) + 0.5_wp * dqdr(1, i, k) * dAmat(1, j, k) &
                                + 0.5_wp * dqdr(1, j, k) * dAmat(1, i, k)
            rxr(2, 1) = rxr(2, 1) + 0.5_wp * dqdr(2, i, k) * dAmat(1, j, k) &
                                + 0.5_wp * dqdr(2, j, k) * dAmat(1, i, k)
            rxr(3, 1) = rxr(3, 1) + 0.5_wp * dqdr(3, i, k) * dAmat(1, j, k) &
                                + 0.5_wp * dqdr(3, j, k) * dAmat(1, i, k)
            rxr(2, 2) = rxr(2, 2) + 0.5_wp * dqdr(2, i, k) * dAmat(2, j, k) &
                                + 0.5_wp * dqdr(2, j, k) * dAmat(2, i, k)
            rxr(3, 2) = rxr(3, 2) + 0.5_wp * dqdr(3, i, k) * dAmat(2, j, k) &
                                + 0.5_wp * dqdr(3, j, k) * dAmat(2, i, k)
            rxr(3, 3) = rxr(3, 3) + 0.5_wp * dqdr(3, i, k) * dAmat(3, j, k) &
                                + 0.5_wp * dqdr(3, j, k) * dAmat(3, i, k)
         end do
         ! symmetrize
         rxr(1, 2) = rxr(2, 1)
         rxr(1, 3) = rxr(3, 1)
         rxr(2, 3) = rxr(3, 2)

         call bmat_accum_pairblock_packed(n, hess, i, j, -kq*rxr)
      end do
   end do

   ! ∂²(qA)/(∂Ri∂q)·∂q/∂Rj
   ! hessian = hessian + reshape(matmul(reshape(dqdr,(/3*n,m/)),&
   !    transpose(reshape(dAmat,(/3*n,m/)))),(/3,n,3,n/))
   !call dgemm('n','t',3*n,m,3*n,+1.0_wp,dqdr,3*n,dAmat,3*n,1.0_wp,hessian,3*n)
   !call dgemm('n','t',3*n,m,3*n,+1.0_wp,dAmat,3*n,dqdr,3*n,1.0_wp,hessian,3*n)

end subroutine add_eeq_hessian

end module xtb_modelhessian_eeq
