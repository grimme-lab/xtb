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
   use xtb_mctc_blas, only : mctc_dot, mctc_gemv
   use xtb_mctc_constants, only : fourpi, pi
   use xtb_mctc_convert, only : aatoau, autoaa
   use xtb_mctc_search, only : bisectSearch
   use xtb_param_vdwradd3, only : getVanDerWaalsRadD3
   use xtb_solv_ddcosmo_core
   use xtb_solv_ddcosmo_solver
   use xtb_solv_lebedev, only : gridSize, getAngGrid
   use xtb_solv_sasa, only : compute_numsa
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
      type(TDomainDecomposition) :: ddCosmo

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

      !> Interaction matrix with surface charges jmat(ncav, n)
      real(wp), allocatable :: jmat(:, :)

      !> angular grid
      real(wp), allocatable :: angGrid(:, :)
      real(wp), allocatable :: angWeight(:)

      !> cut-off radius for the SASA NN list
      real(wp) :: srcut

      !> number of neighbors for SASA computation
      integer, allocatable :: nnsas(:)

      !> neighbors of an atom for SASA
      integer, allocatable :: nnlists(:, :)

      !> Atom specific surface data
      real(wp), allocatable :: vdwsa(:)
      real(wp), allocatable :: wrp(:)
      real(wp), allocatable :: trj2(:, :)

      !> Atomic surfaces
      real(wp) :: gsasa
      real(wp), allocatable :: gamsasa(:)
      real(wp), allocatable :: sasa(:)

      !> Molecular Surface gradient
      real(wp), allocatable :: dsdr(:, :)
      real(wp), allocatable :: dsdrt(:, :, :)

   contains

      !> Update coordinates and internal state
      procedure :: update

      !> Add potential shift
      procedure :: addShift

      !> Calculate solvation energy
      procedure :: getEnergy

      !> Calculate derivatives of solvation energy
      procedure :: addGradient

      !> Write cavity information as cosmo file
      procedure :: writeCosmoFile

   end type TCosmo


   !> Initialize data straucture
   interface init
      module procedure :: initCosmo
   end interface init

   integer, parameter :: ndiis = 25

   !> Smoothing dielectric function parameters
   real(wp), parameter :: w = 0.3_wp*aatoau
   real(wp), parameter :: w3 = w*(w*w)
   real(wp), parameter :: ah0 = 0.5_wp
   real(wp), parameter :: ah1 = 3._wp/(4.0_wp*w)
   real(wp), parameter :: ah3 = -1._wp/(4.0_wp*w3)

   !> Surface tension (in au)
   real(wp), parameter :: gammas = 1.0e-5_wp

   !> Salt screening
   real(wp), parameter :: kappaConst = 0.7897e-3_wp

contains


!> Initialize data straucture
subroutine initCosmo(self, env, num, dielectricConst, nAng, radScale, vdwRad, &
      & surfaceTension, probeRad, rOffset)

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

   !> Van-der-Waals Radii
   real(wp), intent(in) :: vdwRad(:)

   !> Surface tension scaling
   real(wp), intent(in) :: surfaceTension(:)

   !> Probe radius of the solvent
   real(wp), intent(in) :: probeRad

   !> Offset for surface integration cutoff
   real(wp), intent(in) :: rOffset

   integer :: igrid, stat, iat, izp
   real(wp) :: r
   real(wp), allocatable :: angGrid(:, :), angWeight(:)
   type(TDomainDecompositionInput) :: input

   self%nat = size(num)
   allocate(self%rvdw(self%nat))
   self%rvdw(:) = radScale*vdwRad(num)
   self%dielectricConst = dielectricConst
   ! choose the lebedev grid with number of points closest to nAng:
   call bisectSearch(igrid, gridSize, nAng/2)
   allocate(angGrid(3, gridSize(igrid)))
   allocate(angWeight(gridSize(igrid)))
   call getAngGrid(igrid, angGrid, angWeight, stat)
   if (stat /= 0) then
      call env%error("Could not initialize angular grid for COSMO model", source)
      return
   end if

   input = TDomainDecompositionInput(lmax=6, conv=1.0e-8_wp, eta=0.2_wp)

   call initDomainDecomposition(self%ddCosmo, input, self%rvdw, angWeight, angGrid)

   allocate(self%nnsas(self%nat))
   allocate(self%nnlists(self%nat, self%nat))
   allocate(self%vdwsa(self%nat))
   allocate(self%trj2(2, self%nat))
   allocate(self%wrp(self%nat))
   allocate(self%gamsasa(self%nat))
   allocate(self%sasa(self%nat))
   allocate(self%dsdr(3, self%nat))
   allocate(self%dsdrt(3, self%nat, self%nat))
   do iat = 1, self%nat
      izp = num(iat)
      self%vdwsa(iat) = vdwRad(izp) + probeRad
      self%trj2(1, iat) = (self%vdwsa(iat)-w)**2
      self%trj2(2, iat) = (self%vdwsa(iat)+w)**2
      r=self%vdwsa(iat)+w
      self%wrp(iat)=(0.25_wp/w+ &
         &            3.0_wp*ah3*(0.2_wp*r*r-0.5_wp*r*self%vdwsa(iat)+ &
         &            self%vdwsa(iat)*self%vdwsa(iat)/3.0_wp))*r*r*r
      r=self%vdwsa(iat)-w
      self%wrp(iat)=self%wrp(iat)-(0.25/w+ &
         &    3.0_wp*ah3*(0.2_wp*r*r-0.5_wp*r*self%vdwsa(iat)+ &
         &            self%vdwsa(iat)*self%vdwsa(iat)/3.0_wp))*r*r*r
      self%gamsasa(iat) = surfaceTension(izp)
   end do
   self%srcut = 2 * (w + maxval(self%vdwsa)) + rOffset

   ! choose the lebedev grid with number of points closest to nAng:
   call bisectSearch(igrid, gridSize, nAng)
   self%nAng = gridSize(igrid)
   allocate(self%angGrid(3, gridSize(igrid)))
   allocate(self%angWeight(gridSize(igrid)))
   call getAngGrid(igrid, self%angGrid, self%angWeight, stat)
   if (stat /= 0) then
      call env%error("Could not initialize angular grid for SASA model", source)
      return
   end if

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

   ! initialize the neighbor list
   call update_nnlist_sasa(self%nat, xyz, self%srcut, self%nnsas, self%nnlists)

   ! compute solvent accessible surface and its derivatives
   call compute_numsa(self%nat, self%nnsas, self%nnlists, xyz, self%vdwsa, &
      & self%wrp, self%trj2, self%angWeight, self%angGrid, &
      & self%sasa, self%dsdrt)

   ! contract surface gradient
   call mctc_gemv(self%dsdrt, self%gamsasa, self%dsdr)
   self%gsasa = mctc_dot(self%sasa, self%gamsasa)

   if (allocated(self%phi)) deallocate(self%phi)
   if (allocated(self%psi)) deallocate(self%psi)
   if (allocated(self%sigma)) deallocate(self%sigma)
   if (allocated(self%s)) deallocate(self%s)

   call ddupdate(self%ddcosmo, xyz)

   allocate(self%phi(self%ddCosmo%ncav), self%psi(self%ddCosmo%nylm, self%nat))
   allocate(self%jmat(self%ddCosmo%ncav, self%nat))

   call getCoulombMatrix(xyz, self%ddCosmo%ccav, self%jmat)

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

   real(wp) :: xx(1, 1), keps
   logical :: restart

   keps = 0.5_wp * (1.0_wp - 1.0_wp/self%dielectricConst)
   restart = allocated(self%sigma)
   if (.not.allocated(self%sigma)) then
      allocate(self%sigma(self%ddCosmo%nylm, self%nat))
   end if

   call getPhi(qat, self%jmat, self%phi)

   call solveCosmoDirect(self%ddCosmo, env, .true., self%phi, xx, self%sigma, restart)

   restart = allocated(self%s)
   if (.not.allocated(self%s)) then
      allocate(self%s(self%ddCosmo%nylm, self%nat))
   end if

   call getPsi(qat, self%psi)

   ! solve adjoint ddCOSMO equation to get full potential contributions
   call solveCosmoAdjoint(self%ddCosmo, env, self%psi, self%s, restart)

   ! we abuse Phi to store the unpacked and scaled value of s
   call getZeta(self%ddCosmo, keps, self%s, self%phi)
   ! and contract with the Coulomb matrix
   call mctc_gemv(self%jmat, self%phi, atomicShift, alpha=-1.0_wp, &
      & beta=1.0_wp, trans='t')

   atomicShift(:) = atomicShift + (keps * sqrt(fourpi)) * self%sigma(1, :)

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

   energy = 0.5_wp * (1.0_wp - 1.0_wp/self%dielectricConst) &
      & * mctc_dot(self%sigma, self%psi) + self%gsasa

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

   integer :: ii, iat, ig
   real(wp) :: xx(1), esolv, keps
   real(wp), allocatable :: fx(:, :), zeta(:), ef(:, :)

   gradient = gradient + self%dsdr

   keps = 0.5_wp * (1.0_wp - 1.0_wp/self%dielectricConst)

   allocate(fx(3, self%nat), zeta(self%ddCosmo%ncav), &
      & ef(3, max(self%nat, self%ddCosmo%ncav)))

   call solveCosmoAdjoint(self%ddCosmo, env, self%psi, self%s, .true., &
      & accuracy=self%ddCosmo%conv*1e-3_wp)

   ! reset Phi
   call getPhi(qat, self%jmat, self%phi)

   ! now call the routine that computes the ddcosmo specific contributions
   ! to the forces.
   call forces(self%ddCosmo, keps, self%phi, self%sigma, self%s, fx)

   ! form the "zeta" intermediate
   call getZeta(self%ddCosmo, keps, self%s, zeta)

   ! 1. solute's electric field at the cav points times zeta:
   !    compute the electric field
   call efld(self%nat, qat, self%ddCosmo%xyz, self%ddCosmo%ncav, &
      & self%ddCosmo%ccav, ef)

   ! contract it with the zeta intermediate
   ii = 0
   do iat = 1, self%ddCosmo%nat
      do ig = 1, self%ddCosmo%ngrid
         if (self%ddCosmo%ui(ig, iat) > 0.0_wp) then
            ii = ii + 1
            fx(:, iat) = fx(:, iat) - zeta(ii)*ef(:, ii)
         end if
      end do
   end do

   ! 2. "zeta's" electric field at the nuclei times the charges.
   !    compute the "electric field"
   call efld(self%ddCosmo%ncav, zeta, self%ddCosmo%ccav, self%nat, &
      & self%ddCosmo%xyz, ef)

   ! contract it with the solute's charges.
   do iat = 1, self%ddCosmo%nat
      fx(:, iat) = fx(:, iat) - ef(:, iat)*qat(iat)
   end do

   gradient(:, :) = gradient(:, :) - fx

end subroutine addGradient


!> Evaluate the Coulomb interactions between the atomic sides (xyz) and the
!> surface elements of the cavity (ccav).
subroutine getCoulombMatrix(xyz, ccav, jmat)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: ccav(:, :)
   real(wp), intent(inout) :: jmat(:, :)

   integer :: ic, j
   real(wp) :: vec(3), d2, d

   jmat(:, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(ccav, xyz, jmat) private(ic, j, vec, d2, d)
   do ic = 1, size(ccav, 2)
      do j = 1, size(xyz, 2)
         vec(:) = ccav(:, ic) - xyz(:, j)
         d2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         d = sqrt(d2)
         jmat(ic, j) = 1.0_wp / d
      end do
   end do

end subroutine getCoulombMatrix


!> Routine to compute the psi vector
subroutine getPsi(charge, psi)
   real(wp), intent(in) :: charge(:)
   real(wp), intent(out) :: psi(:, :)

   integer :: iat
   real(wp) :: fac

   fac = sqrt(fourpi)
   psi(:,:) = 0.0_wp

   do iat = 1, size(charge)
      psi(1, iat) = fac*charge(iat)
   end do

end subroutine getPsi


!> Routine to compute the potential vector
subroutine getPhi(charge, jmat, phi)
   real(wp), intent(in) :: charge(:)
   real(wp), intent(in) :: jmat(:, :)
   real(wp), intent(out) :: phi(:)

   phi(:) = 0.0_wp

   call mctc_gemv(jmat, charge, phi)

end subroutine getPhi


!> Wrapper for the linear solvers for COSMO equation
!    L sigma = G
!  This routine performs the following operations :
!   - allocates memory for the linear solvers
!   - if star is false and cart is true, assembles the right-hand side for the COSMO
!     equations.
!   - computes a guess for the solution (using the inverse diagonal);
!   - calls the iterative solver;
subroutine solveCosmoDirect(ddCosmo, env, cart, phi, glm, sigma, restart)

   !> Error source
   character(len=*), parameter :: source = 'cosmo::solveCosmoDirect'
   type(TDomainDecomposition), intent(in) :: ddCosmo

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> true:  the right-hand side for the COSMO has to be assembled
   !         inside this routine and the unscaled potential at the
   !         external points of the cavity is provided in phi.
   !  false: the right-hand side for the COSMO equations is provided
   !         in glm.
   logical, intent(in) :: cart

   !> Contains the potential at the external cavity points if cart is true.
   !  phi is not referenced in any other case.
   real(wp), intent(in) :: phi(:)

   !> Contains the right-hand side for the COSMO equations if cart is false.
   !  glm is not referenced in any other case
   real(wp), intent(in) :: glm(:, :)

   !> The solution to the COSMO (adjoint) equations
   real(wp), intent(inout) :: sigma(:, :)

   !> Initial guess is provided on sigma
   logical, intent(in) :: restart

   integer :: iat, istatus, n_iter, info, c1, c2, cr
   real(wp) :: tol, r_norm
   logical :: ok

   real(wp), allocatable  :: g(:, :), rhs(:, :), work(:, :)

   ! parameters for the solver and matvec routine
   tol     = ddCosmo%conv
   n_iter  = 200

   ! DIRECT COSMO EQUATION L X = g

   ! allocate workspace for rhs
   allocate(rhs(ddCosmo%nylm, ddCosmo%nat), stat=istatus)
   if (istatus /= 0) then
      write(*, *) ' cosmo: [2] failed allocation'
   endif

   ! 1. RHS
   ! assemble rhs
   if (cart) then

      ! allocate workspace for weighted potential
      allocate(g(ddCosmo%ngrid, ddCosmo%nat) , stat=istatus)
      if (istatus /= 0) then
         write(*, *) ' cosmo: [3] failed allocation'
      endif

      ! weight the potential...
      call wghpot(ddCosmo, phi, g)

      ! ... and compute its multipolar expansion
      do iat = 1, ddCosmo%nat
         call intrhs(ddCosmo, iat, g(:, iat), rhs(:, iat))
      enddo

      ! deallocate workspace
      deallocate(g , stat=istatus)
      if (istatus /= 0) then
         write(*, *) 'cosmo: [1] failed deallocation'
      endif

   else
      ! no need to manipulate rhs
      rhs = glm
   end if

   ! 2. INITIAL GUESS
   if (.not.restart) then
      do iat = 1, ddCosmo%nat
         sigma(:, iat) = ddCosmo%facl(:)*rhs(:, iat)
      end do
   end if

   ! 3. SOLVER CALL
   ! Jacobi method :
   ! L X = (diag + offdiag) X = g   ==>    X = diag^-1 (g - offdiag X_guess)
   ! action of  diag^-1 :  ldm1x
   ! action of  offdiag :  lx
   call jacobi_diis(ddCosmo, ddCosmo%nat*ddCosmo%nylm, ddCosmo%iprint, &
      & ndiis, 4, tol, rhs, sigma, n_iter, ok, lx, ldm1x, hnorm)

   ! check solution
   if (.not.ok) then
      call env%error('direct ddCOSMO did not converge!', source)
      return
   endif

end subroutine solveCosmoDirect


!> Wrapper for the linear solvers for adjoint COSMO equation
!>   L^* sigma = Psi
!> This routine performs the following operations :
!>  - allocates memory for the linear solvers
!>  - computes a guess for the solution (using the inverse diagonal);
!>  - calls the iterative solver;
subroutine solveCosmoAdjoint(ddCosmo, env, psi, sigma, restart, accuracy)
   ! Error source
   character(len=*), parameter :: source = 'cosmo::solveCosmoAdjoint'
   type(TDomainDecomposition), intent(in) :: ddCosmo

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> The psi vector. it is used as a right-hand side
   real(wp), intent(in) :: psi(:, :)

   !> The solution to the COSMO (adjoint) equations
   real(wp), intent(inout) :: sigma(:, :)

   !> Initial guess is provided on sigma
   logical, intent(in) :: restart

   !> Overwrite accuracy
   real(wp), intent(in), optional :: accuracy

   integer :: iat, istatus, n_iter, info, c1, c2, cr
   real(wp) :: tol, r_norm
   logical :: ok

   real(wp), allocatable  :: g(:, :), rhs(:, :), work(:, :)

   ! parameters for the solver and matvec routine
   if (present(accuracy)) then
      tol = accuracy
   else
      tol = ddCosmo%conv
   end if
   n_iter  = 200

   ! 1. INITIAL GUESS
   if (.not.restart) then
      do iat = 1, ddCosmo%nat
         sigma(:, iat) = ddCosmo%facl(:)*psi(:, iat)
      end do
   end if

   ! 2. SOLVER CALL
   ! Jacobi method : see above
   call jacobi_diis(ddCosmo, ddCosmo%nat*ddCosmo%nylm, ddCosmo%iprint, &
      & ndiis, 4, tol, psi, sigma, n_iter, ok, lstarx, ldm1x, hnorm)

   ! check solution
   if (.not.ok) then
      call env%error('adjoint ddCOSMO did not converge!', source)
      return
   endif

end subroutine solveCosmoAdjoint


!> Compute
!
! \zeta(n, i) =
!
!  1/2 f(\eps) sum w_n U_n^i Y_l^m(s_n) [S_i]_l^m
!              l, m
!
subroutine getZeta(ddCosmo, keps, s, zeta)
   type(TDomainDecomposition), intent(in) :: ddCosmo
   real(wp), intent(in) :: keps
   real(wp), intent(in) :: s(:, :) ! [ddCosmo%nylm, ddCosmo%nat]
   real(wp), intent(inout) :: zeta(:) ! [ddCosmo%ncav]

   integer :: its, iat, ii

   ii = 0
   do iat = 1, ddCosmo%nat
      do its = 1, ddCosmo%ngrid
         if (ddCosmo%ui(its, iat) > 0.0_wp) then
            ii = ii + 1
            zeta(ii) = keps * ddCosmo%w(its) * ddCosmo%ui(its, iat) &
               & * dot_product(ddCosmo%basis(:, its), s(:, iat))
         end if
      end do
   end do

end subroutine getZeta


!> Sample driver for the calculation of the ddCOSMO forces.
subroutine forces(ddCosmo, keps, phi, sigma, s, fx)
   type(TDomainDecomposition), intent(in) :: ddCosmo
   real(wp), intent(in) :: keps
   real(wp), intent(in) :: phi(:)
   real(wp), intent(in) :: sigma(:, :)
   real(wp), intent(in) :: s(:, :)
   real(wp), intent(inout) :: fx(:, :)

   integer :: iat, ig, ii, c1, c2, cr
   real(wp) :: fep

   real(wp), allocatable :: xi(:, :), phiexp(:, :), zeta(:), ef(:, :)
   real(wp), allocatable :: basloc(:), dbsloc(:, :), vplm(:), vcos(:), vsin(:)

   allocate (xi(ddCosmo%ngrid, ddCosmo%nat), phiexp(ddCosmo%ngrid, ddCosmo%nat))
   allocate (basloc(ddCosmo%nylm), dbsloc(3, ddCosmo%nylm), vplm(ddCosmo%nylm), &
      & vcos(ddCosmo%lmax+1), vsin(ddCosmo%lmax+1))

   ! compute xi:
   !$omp parallel do default(none) collapse(2) schedule(runtime) &
   !$omp shared(ddCosmo, s, xi) private(iat, ig)
   do iat = 1, ddCosmo%nat
      do ig = 1, ddCosmo%ngrid
         xi(ig, iat) = dot_product(s(:, iat), ddCosmo%basis(:, ig))
      end do
   end do

   ! expand the potential on a sphere-by-sphere basis (needed for parallelism):
   ii = 0
   phiexp(:, :) = 0.0_wp
   do iat = 1, ddCosmo%nat
      do ig = 1, ddCosmo%ngrid
         if (ddCosmo%ui(ig, iat) > 0.0_wp) then
            ii = ii + 1
            phiexp(ig, iat) = phi(ii)
         end if
      end do
   end do

   fx(:, :) = 0.0_wp
   do iat = 1, ddCosmo%nat
      call fdoka(ddCosmo, iat, sigma, xi(:, iat), basloc, dbsloc, vplm, &
         & vcos, vsin, fx(:, iat))
      call fdokb(ddCosmo, iat, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, &
         & fx(:, iat))
      call fdoga(ddCosmo, iat, xi, phiexp, fx(:, iat))
   end do

   deallocate (basloc, dbsloc, vplm, vcos, vsin)
   deallocate (xi, phiexp)

   ! scale the forces time the cosmo factor:
   fx  = keps*fx

end subroutine forces


!> Computes the electric field produced by the sources src (nsrc point charges
!  with coordinates csrc) at the ntrg target points ctrg:
subroutine efld(nsrc, src, csrc, ntrg, ctrg, ef)
   integer, intent(in) :: nsrc, ntrg
   real(wp), intent(in) :: src(:)
   real(wp), intent(in) :: csrc(:, :)
   real(wp), intent(in) :: ctrg(:, :)
   real(wp), intent(inout) :: ef(:, :)

   integer :: i, j
   real(wp) :: vec(3), r2, rr, r3, f
   real(wp), parameter :: zero=0.0_wp

   ef(:, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp reduction(+:ef) shared(ntrg, nsrc, ctrg, csrc, src) &
   !$omp private(j, i, f, vec, r2, rr, r3)
   do j = 1, ntrg
      do i = 1, nsrc
         vec(:) = ctrg(:, j) - csrc(:, i)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         rr = sqrt(r2)
         r3 = r2*rr
         f = src(i)/r3
         ef(:, j) = ef(:, j) + f*vec
      end do
   end do

end subroutine efld


!> setup a pairlist and compute pair distances of all neighbors
!  within thresholds lrcut and srcut.
subroutine update_nnlist_sasa(nat, xyz, srcut, nnsas, nnlists)

   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: srcut
   integer, intent(out) :: nnsas(:)
   integer, intent(out) :: nnlists(:, :)

   integer i1, i2
   real(wp) rcutn2, srcut2
   real(wp) x, y, z, dr2

   srcut2 = srcut*srcut

   nnsas=0
   nnlists=0
   do i1 = 1, nat
      do i2 = 1, i1-1
         x=xyz(1,i1)-xyz(1,i2)
         y=xyz(2,i1)-xyz(2,i2)
         z=xyz(3,i1)-xyz(3,i2)
         dr2=x**2+y**2+z**2
         if(dr2.lt.srcut2) then
            nnsas(i1) = nnsas(i1) + 1
            nnsas(i2) = nnsas(i2) + 1
            nnlists(nnsas(i1),i1)=i2
            nnlists(nnsas(i2),i2)=i1
         endif
      end do
   end do

end subroutine update_nnlist_sasa

!> Write a COSMO file output
subroutine writeCosmoFile(self, unit, num, sym, xyz, qat, energy)

   !> COSMO container
   class(TCosmo), intent(in) :: self

   !> Formatted unit for output
   integer, intent(in) :: unit

   !> Atomic numbers for each atom
   integer, intent(in) :: num(:)

   !> Element symbols for each atom
   character(len=*), intent(in) :: sym(:)

   !> Atomic coordinates
   real(wp), intent(in) :: xyz(:,:)

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Total energy
   real(wp), intent(in) :: energy

   integer :: ii, ig, iat
   real(wp) :: dielEnergy, keps
   real(wp), allocatable :: phi(:), zeta(:), area(:)

   allocate(phi(self%ddCosmo%ncav), zeta(self%ddCosmo%ncav), area(self%ddCosmo%ncav))
   ! Reset potential on the cavity, note that the potential is expected in e/Å
   call getPhi(qat, self%jmat, phi)
   ii = 0
   do iat = 1, self%ddCosmo%nat
      do ig = 1, self%ddCosmo%ngrid
         if (self%ddCosmo%ui(ig, iat) > 0.0_wp) then
            ii = ii + 1
            ! Calculate surface charge per area
            zeta(ii) = self%ddCosmo%w(ig) * self%ddCosmo%ui(ig, iat) &
               & * dot_product(self%ddCosmo%basis(:, ig), self%s(:, iat))
            ! Save surface area in Ångström²
            area(ii) = self%ddCosmo%w(ig) * autoaa**2 * self%rvdw(iat)**2
         end if
      end do
   end do


    ! Dielectric energy is the energy on the dielectric continuum
   keps = 0.5_wp * (1.0_wp - 1.0_wp/self%dielectricConst)
   dielEnergy = keps * dot_product(zeta, phi)

   write(unit, '(a)') &
      & "$info", &
      & "prog.: xtb"

   write(unit, '(a)') &
      & "$cosmo"
   write(unit, '(2x, a:, "=", g0)') &
      & "epsilon", self%dielectricConst

   write(unit, '(a)') &
      & "$cosmo_data"
   write(unit, '(2x, a:, "=", g0)') &
      & "fepsi", keps, &
      & "area", sum(area)

   write(unit, '(a)') &
      & "$coord_rad", &
      & "#atom   x                  y                  z             element  radius [A]"
   do iat = 1, size(xyz, 2)
      write(unit, '(i4, 3(1x, f18.14), 2x, a4, 1x, f9.5)') &
         & iat, xyz(:, iat), sym(iat), self%rvdw(iat)*autoaa
   end do

   write(unit, '(a)') &
      & "$coord_car", &
      & "!BIOSYM archive 3", &
   & "PBC=OFF", &
      & "coordinates from COSMO calculation", &
      & "!DATE"
   do iat = 1, size(xyz, 2)
      write(unit, '(a, i0, t5, 3(1x, f14.9), 1x, "COSM 1", 2(6x, a2), 1x, f6.3)') &
         & trim(sym(iat)), iat, xyz(:, iat)*autoaa, sym(iat), sym(iat), 0.0_wp
   end do
   write(unit, '(a)') &
      & "end", "end"

   write(unit, '(a)') &
      & "$screening_charge"
   write(unit, '(2x, a:, "=", g0)') &
      & "cosmo", sum(zeta), &
      & "correction", 0.0_wp, &
      & "total", sum(zeta)

   write(unit, '(a)') &
      & "$cosmo_energy"
   write(unit, '(2x, a:, "=", f21.10)') &
      & "Total energy [a.u.]            ", energy, &
      & "Total energy + OC corr. [a.u.] ", energy, &
      & "Total energy corrected [a.u.]  ", energy, &
      & "Dielectric energy [a.u.]       ", dielEnergy, &
      & "Diel. energy + OC corr. [a.u.] ", dielEnergy

   write(unit, '(a)') &
      & "$segment_information", &
      & "# n             - segment number", &
      & "# atom          - atom associated with segment n", &
      & "# position      - segment coordinates [a.u.]", &
      & "# charge        - segment charge (corrected)", &
      & "# area          - segment area [A**2]", &
      & "# potential     - solute potential on segment (A length scale)", &
      & "#", &
      & "#  n   atom              position (X, Y, Z)                   charge         area        charge/area     potential", &
      & "#", &
      & "#"

   ii = 0
   do iat = 1, self%ddCosmo%nat
      do ig = 1, self%ddCosmo%ngrid
         if (self%ddCosmo%ui(ig, iat) > 0.0_wp) then
            ii = ii + 1
            write(unit, '(2i5, 7(1x, f14.9))') &
               & ii, iat, self%ddCosmo%ccav(:, ii), &
               & zeta(ii), area(ii), zeta(ii)/area(ii), phi(ii)/autoaa
         end if
      end do
   end do
end subroutine writeCosmoFile

end module xtb_solv_cosmo
