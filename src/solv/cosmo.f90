! This file is part of xtb.
!
! Copyright (C) 2020 Sebastian Ehlert
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

!> COSMO implementation based on ddCOSMO backend
module xtb_solv_cosmo
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas, only : mctc_dot
   use xtb_mctc_constants, only : fourpi
   use xtb_param_vdwradd3, only : getVanDerWaalsRadD3
   use xtb_solv_ddcosmo_core
   use xtb_solv_ddcosmo_solver
   use xtb_type_environment, only : TEnvironment
   use xtb_type_solvation, only : TSolvation
   implicit none
   private

   public :: TCosmo, init


   !> The COSMO solvation model
   type, extends(TSolvation) :: TCosmo

      !> Number of atoms
      integer :: nat

      !> Actual COSMO calculator
      type(TddCosmo) :: ddCosmo

      !> electrostatic potential phi(ncav)
      real(wp), allocatable :: phi(:)

      !> psi vector psi(nylm, n)
      real(wp), allocatable :: psi(:, :)

      !> ddcosmo solution sigma(nylm, n)
      real(wp), allocatable :: sigma(:, :)

      !> ddcosmo adjoint solution s(nylm, n)
      real(wp), allocatable :: s(:, :)

      !> dielectric constant
      real(wp) :: dielectricConst

      !> Van-der-Waal radii
      real(wp), allocatable :: rvdw(:)

      !> Number of grid point
      integer :: nAng

   contains

      !> Update coordinates and internal state
      procedure :: update

      !> Add potential shift
      procedure :: addShift

      !> Calculate solvation energy
      procedure :: getEnergy

      !> Calculate derivatives of solvation energy
      procedure :: addGradient

   end type TCosmo


   !> Initialize data straucture
   interface init
      module procedure :: initCosmo
   end interface init


contains


!> Initialize data straucture
subroutine initCosmo(self, env, num, dielectricConst, nAng, radScale)

   !> Error source
   character(len=*), parameter :: source = 'solv_cosmo_initCosmo'

   !> Instance of the solvation model
   type(TCosmo), intent(out) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Atomic numbers
   integer, intent(in) :: num(:)

   !> Dielectric constant
   real(wp), intent(in) :: dielectricConst

   !> Number of grid point
   integer, intent(in) :: nAng

   !> Scale van-der-Waals radii
   real(wp), intent(in) :: radScale

   self%nat = size(num)
   allocate(self%rvdw(self%nat))
   self%rvdw(:) = radScale*getVanDerWaalsRadD3(num)
   self%dielectricConst = dielectricConst
   self%nAng = nAng

end subroutine initCosmo


subroutine update(self, env, num, xyz)

   !> Instance of the solvation model
   class(TCosmo), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Atomic numbers
   integer, intent(in) :: num(:)

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   if (allocated(self%phi)) deallocate(self%phi)
   if (allocated(self%psi)) deallocate(self%psi)
   if (allocated(self%sigma)) deallocate(self%sigma)
   if (allocated(self%s)) deallocate(self%s)

   self%ddCosmo%TddControl = TddControl(iprint=0, lmax=6, ngrid=self%nAng, &
      & iconv=8, igrad=1, eps=self%dielectricConst, eta=0.2_wp)

   call ddinit(self%ddCosmo, env, self%nat, xyz, self%rvdw)

   allocate(self%phi(self%ddCosmo%ncav), self%psi(self%ddCosmo%nylm, self%nat))
   allocate(self%s(self%ddCosmo%nylm, self%nat))

end subroutine update


!> Add potential shift
subroutine addShift(self, env, qat, qsh, atomicShift, shellShift)

   !> Error source
   character(len=*), parameter :: source = 'solv_cosmo_addShift'

   !> Instance of the solvation model
   class(TCosmo), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)

   !> Atomic potential shift
   real(wp), intent(inout) :: atomicShift(:)

   !> Shell-resolved potential shift
   real(wp), intent(inout) :: shellShift(:)

   real(wp) :: xx(1), keps
   logical :: restart

   restart = allocated(self%sigma)
   if (.not.allocated(self%sigma)) then
      allocate(self%sigma(self%ddCosmo%nylm, self%nat))
   end if

   call mkrhs(self%nat, qat, self%ddCosmo%csph, self%ddCosmo%ncav, &
      & self%ddCosmo%ccav, self%phi, self%ddCosmo%nylm, self%psi)

   call cosmo(self%ddCosmo, env, .true., self%phi, xx, self%sigma, restart)

   keps = 0.5_wp * ((self%ddCosmo%eps - 1.0_wp)/self%ddCosmo%eps) * sqrt(fourpi)
   atomicShift(:) = atomicShift + keps * self%sigma(1, :)

end subroutine addShift


!> Calculate solvation energy
subroutine getEnergy(self, env, qat, qsh, energy)

   !> Error source
   character(len=*), parameter :: source = 'solv_cosmo_getEnergy'

   !> Instance of the solvation model
   class(TCosmo), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)

   !> Total solvation energy
   real(wp), intent(out) :: energy

   energy = 0.5_wp * ((self%ddCosmo%eps - 1.0_wp)/self%ddCosmo%eps) &
      & * mctc_dot(self%sigma, self%psi)

end subroutine getEnergy


!> Calculate derivatives of solvation energy
subroutine addGradient(self, env, num, xyz, qat, qsh, gradient)

   !> Error source
   character(len=*), parameter :: source = 'solv_cosmo_addGradient'

   !> Instance of the solvation model
   class(TCosmo), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Atomic numbers
   integer, intent(in) :: num(:)

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)

   !> Molecular gradient
   real(wp), intent(inout) :: gradient(:, :)

   integer :: ii, isph, ig
   real(wp) :: xx(1), esolv
   real(wp), allocatable :: fx(:, :), zeta(:), ef(:, :)

   allocate(fx(3, self%nat), zeta(self%ddCosmo%ncav), &
      & ef(3, max(self%nat, self%ddCosmo%ncav)))
   call cosmoStar(self%ddCosmo, env, self%psi, self%s)

   ! now call the routine that computes the ddcosmo specific contributions
   ! to the forces.
   call forces(self%ddCosmo, self%nat, self%phi, self%sigma, self%s, fx)

   ! form the "zeta" intermediate
   call ddmkzeta(self%ddCosmo, self%s, zeta)

   ! 1. solute's electric field at the cav points times zeta:
   !    compute the electric field
   call efld(self%nat, qat, self%ddCosmo%csph, self%ddCosmo%ncav, &
      & self%ddCosmo%ccav, ef)

   ! contract it with the zeta intermediate
   ii = 0
   do isph = 1, self%ddCosmo%nsph
      do ig = 1, self%ddCosmo%ngrid
         if (self%ddCosmo%ui(ig, isph).gt.0.0_wp) then
            ii = ii + 1
            fx(:, isph) = fx(:, isph) - zeta(ii)*ef(:, ii)
         end if
      end do
   end do

   ! 2. "zeta's" electric field at the nuclei times the charges.
   !    compute the "electric field"
   call efld(self%ddCosmo%ncav, zeta, self%ddCosmo%ccav, self%nat, &
      & self%ddCosmo%csph, ef)

   ! contract it with the solute's charges.
   do isph = 1, self%ddCosmo%nsph
      fx(:, isph) = fx(:, isph) - ef(:, isph)*qat(isph)
   end do

   gradient(:, :) = gradient(:, :) - fx

end subroutine addGradient


!> Routine to compute the potential and psi vector
subroutine mkrhs(n, charge, xyz, ncav, ccav, phi, nylm, psi)
   integer, intent(in) :: n, ncav, nylm
   real(wp), intent(in) :: xyz(3, n), charge(n)
   real(wp), intent(in) :: ccav(3, ncav)
   real(wp), intent(inout) :: phi(ncav)
   real(wp), intent(inout) :: psi(nylm, n)

   integer :: isph, ic, j
   real(wp) :: v
   real(wp) :: vec(3), d2, d, pi, fac
   real(wp), parameter :: zero=0.0_wp, one=1.0_wp, four=4.0_wp

   fac = sqrt(fourpi)
   phi = 0.0_wp
   psi = 0.0_wp

   !$omp parallel do default(shared) private(ic, v, j, vec, d2, d)
   do ic = 1, ncav
      v  = 0.0_wp
      do j = 1, n
         vec(:) = ccav(:, ic) - xyz(:, j)
         d2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         d = sqrt(d2)
         v = v + charge(j)/d
      end do
      phi(ic) = v
   end do

   ! psi vector:
   do isph = 1, n
      psi(1, isph) = fac*charge(isph)
   end do

end subroutine mkrhs


!> Wrapper for the linear solvers for COSMO equation
!    L sigma = G
!  This routine performs the following operations :
!   - allocates memory for the linear solvers
!   - if star is false and cart is true, assembles the right-hand side for the COSMO
!     equations.
!   - computes a guess for the solution (using the inverse diagonal);
!   - calls the iterative solver;
subroutine cosmo(ddCosmo, env, cart, phi, glm, sigma, restart)

   !> Error source
   character(len=*), parameter :: source = 'solv_cosmo_cosmo'
   type(TEnvironment), intent(inout) :: env
   type(TddCosmo), intent(in) :: ddCosmo

   !> true:  the right-hand side for the COSMO has to be assembled
   !         inside this routine and the unscaled potential at the
   !         external points of the cavity is provided in phi.
   !  false: the right-hand side for the COSMO equations is provided
   !         in glm.
   logical, intent(in) :: cart

   !> Contains the potential at the external cavity points if cart is true.
   !  phi is not referenced in any other case.
   real(wp), intent(in) :: phi(ddCosmo%ncav)

   !> Contains the right-hand side for the COSMO equations if cart is false.
   !  glm is not referenced in any other case
   real(wp), intent(in) :: glm(ddCosmo%nylm, ddCosmo%nsph)

   !> The solution to the COSMO (adjoint) equations
   real(wp), intent(inout) :: sigma(ddCosmo%nylm, ddCosmo%nsph)

   !> Initial guess is provided on sigma
   logical, intent(in) :: restart

   integer :: isph, istatus, n_iter, info, c1, c2, cr
   real(wp) :: tol, r_norm
   logical :: ok

   real(wp), allocatable  :: g(:, :), rhs(:, :), work(:, :)

   ! parameters for the solver and matvec routine
   tol     = 10.0_wp**(-ddCosmo%iconv)
   n_iter  = 200

   ! DIRECT COSMO EQUATION L X = g

   ! allocate workspace for rhs
   allocate(rhs(ddCosmo%nylm, ddCosmo%nsph), stat=istatus)
   if (istatus .ne. 0) then
      write(*, *) ' cosmo: [2] failed allocation'
   endif

   ! 1. RHS
   ! assemble rhs
   if (cart) then

      ! allocate workspace for weighted potential
      allocate(g(ddCosmo%ngrid, ddCosmo%nsph) , stat=istatus)
      if (istatus .ne. 0) then
         write(*, *) ' cosmo: [3] failed allocation'
      endif

      ! weight the potential...
      call wghpot(ddCosmo, phi, g)

      ! ... and compute its multipolar expansion
      do isph = 1, ddCosmo%nsph
         call intrhs(ddCosmo, isph, g(:, isph), rhs(:, isph))
      enddo

      ! deallocate workspace
      deallocate(g , stat=istatus)
      if (istatus.ne.0) then
         write(*, *) 'cosmo: [1] failed deallocation'
      endif

   else
      ! no need to manipulate rhs
      rhs = glm
   end if

   ! 2. INITIAL GUESS
   if (.not.restart) then
      do isph = 1, ddCosmo%nsph
         sigma(:, isph) = ddCosmo%facl(:)*rhs(:, isph)
      end do
   end if

   ! 3. SOLVER CALL
   ! Jacobi method :
   ! L X = (diag + offdiag) X = g   ==>    X = diag^-1 (g - offdiag X_guess)
   ! action of  diag^-1 :  ldm1x
   ! action of  offdiag :  lx
   call jacobi_diis(ddCosmo, env, ddCosmo%nsph*ddCosmo%nylm, ddCosmo%iprint, &
      & ndiis, 4, tol, rhs, sigma, n_iter, ok, lx, ldm1x, hnorm)

   ! check solution
   if (.not.ok) then
      call env%error('direct ddCOSMO did not converge!', source)
      return
   endif

end subroutine cosmo


!> Wrapper for the linear solvers for adjoint COSMO equation
!    L^* sigma = Psi
!  This routine performs the following operations :
!   - allocates memory for the linear solvers
!   - computes a guess for the solution (using the inverse diagonal);
!   - calls the iterative solver;
subroutine cosmoStar(ddCosmo, env, psi, sigma)
   !> Error source
   character(len=*), parameter :: source = 'solv_cosmo_cosmoStar'
   type(TEnvironment), intent(inout) :: env
   type(TddCosmo), intent(in) :: ddCosmo

   !> The psi vector. it is used as a right-hand side
   real(wp), intent(in) :: psi(ddCosmo%nylm, ddCosmo%nsph)

   !> The solution to the COSMO (adjoint) equations
   real(wp), intent(inout) :: sigma(ddCosmo%nylm, ddCosmo%nsph)

   integer :: isph, istatus, n_iter, info, c1, c2, cr
   real(wp) :: tol, r_norm
   logical :: ok

   real(wp), allocatable  :: g(:, :), rhs(:, :), work(:, :)

   ! parameters for the solver and matvec routine
   tol     = 10.0_wp**(-ddCosmo%iconv-3)
   n_iter  = 200

   ! 1. INITIAL GUESS
   do isph = 1, ddCosmo%nsph
      sigma(:, isph) = ddCosmo%facl(:)*psi(:, isph)
   enddo

   ! 2. SOLVER CALL
   ! Jacobi method : see above
   call jacobi_diis(ddCosmo, env, ddCosmo%nsph*ddCosmo%nylm, ddCosmo%iprint, &
      & ndiis, 4, tol, psi, sigma, n_iter, ok, lstarx, ldm1x, hnorm)

   ! check solution
   if (.not.ok) then
      call env%error('adjoint ddCOSMO did not converge!', source)
      return
   endif

end subroutine cosmoStar


!> Sample driver for the calculation of the ddCOSMO forces.
subroutine forces(ddCosmo, n, phi, sigma, s, fx)
   type(TddCosmo), intent(in) :: ddCosmo
   integer, intent(in) :: n
   real(wp), intent(in) :: phi(ddCosmo%ncav)
   real(wp), intent(in) :: sigma(ddCosmo%nylm, ddCosmo%nsph)
   real(wp), intent(in) :: s(ddCosmo%nylm, ddCosmo%nsph)
   real(wp), intent(inout) :: fx(3, n)

   integer :: isph, ig, ii, c1, c2, cr
   real(wp) :: fep

   real(wp), allocatable :: xi(:, :), phiexp(:, :), zeta(:), ef(:, :)
   real(wp), allocatable :: basloc(:), dbsloc(:, :), vplm(:), vcos(:), vsin(:)

   allocate (xi(ddCosmo%ngrid, ddCosmo%nsph), phiexp(ddCosmo%ngrid, ddCosmo%nsph))
   allocate (basloc(ddCosmo%nylm), dbsloc(3, ddCosmo%nylm), vplm(ddCosmo%nylm), &
      & vcos(ddCosmo%lmax+1), vsin(ddCosmo%lmax+1))

   ! initialize the timer:
   call system_clock(count_rate=cr)
   call system_clock(count=c1)

   ! compute xi:
   !$omp parallel do default(shared) private(isph, ig)
   do isph = 1, ddCosmo%nsph
      do ig = 1, ddCosmo%ngrid
         xi(ig, isph) = dot_product(s(:, isph), ddCosmo%basis(:, ig))
      end do
   end do
   !$omp end parallel do

   if (ddCosmo%iprint.ge.4) call ptcart(ddCosmo, 'xi', ddCosmo%nsph, 0, xi)

   ! expand the potential on a sphere-by-sphere basis (needed for parallelism):

   ii = 0
   phiexp = 0.0_wp
   do isph = 1, ddCosmo%nsph
      do ig = 1, ddCosmo%ngrid
         if (ddCosmo%ui(ig, isph).gt.0.0_wp) then
            ii = ii + 1
            phiexp(ig, isph) = phi(ii)
         end if
      end do
   end do

   fx = 0.0_wp
   do isph = 1, ddCosmo%nsph
      call fdoka(ddCosmo, isph, sigma, xi(:, isph), basloc, dbsloc, vplm, &
         & vcos, vsin, fx(:, isph))
      call fdokb(ddCosmo, isph, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, &
         & fx(:, isph))
      call fdoga(ddCosmo, isph, xi, phiexp, fx(:, isph))
   end do

   2000 format(1x, 'ddCOSMO-only contributions to the forces (atomic units):', /, &
      1x, ' atom', 15x, 'x', 15x, 'y', 15x, 'z')

   if (ddCosmo%iprint.ge.4) then
      write(iout, 2000)
      do isph = 1, ddCosmo%nsph
         write(6, '(1x, i5, 3f16.8)') isph, fx(:, isph)
      end do
   end if

   deallocate (basloc, dbsloc, vplm, vcos, vsin)

   call system_clock(count=c2)
   if (ddCosmo%iprint.gt.0) then
      write(iout, 1010) dble(c2-c1)/dble(cr)
      1010 format(' the computation of the ddCOSMO part of the forces took ', f8.3, ' seconds.')
   end if

   deallocate (xi, phiexp)

   ! scale the forces time the cosmo factor:

   fep = pt5*(ddCosmo%eps-1.0_wp)/ddCosmo%eps
   fx  = fep*fx

end subroutine forces


!> Computes the electric field produced by the sources src (nsrc point charges
!  with coordinates csrc) at the ntrg target points ctrg:
subroutine efld(nsrc, src, csrc, ntrg, ctrg, ef)
   integer, intent(in) :: nsrc, ntrg
   real(wp), intent(in) :: src(nsrc)
   real(wp), intent(in) :: csrc(3, nsrc)
   real(wp), intent(in) :: ctrg(3, ntrg)
   real(wp), intent(inout) :: ef(3, ntrg)

   integer :: i, j
   real(wp) :: vec(3), r2, rr, r3, f, e(3)
   real(wp), parameter :: zero=0.0_wp

   ef(:, :) = 0.0_wp
   !$omp parallel do default(shared) private(j, i, vec, r2, rr, r3, e)
   do j = 1, ntrg
      e(:) = 0.0_wp
      do i = 1, nsrc
         vec(:) = ctrg(:, j) - csrc(:, i)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         rr = sqrt(r2)
         r3 = r2*rr
         f = src(i)/r3
         e(:) = e(:) + f*vec
      end do
      ef(:, j) = e
   end do

end subroutine efld


end module xtb_solv_cosmo
