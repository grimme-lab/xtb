!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!      888      888  .d8888b.   .d88888b.   .d8888b.  888b     d888  .d88888b.
!      888      888 d88P  Y88b d88P" "Y88b d88P  Y88b 8888b   d8888 d88P" "Y88b
!      888      888 888    888 888     888 Y88b.      88888b.d88888 888     888
!  .d88888  .d88888 888        888     888  "Y888b.   888Y88888P888 888     888
! d88" 888 d88" 888 888        888     888     "Y88b. 888 Y888P 888 888     888
! 888  888 888  888 888    888 888     888       "888 888  Y8P  888 888     888
! Y88b 888 Y88b 888 Y88b  d88P Y88b. .d88P Y88b  d88P 888   "   888 Y88b. .d88P
!  "Y88888  "Y88888  "Y8888P"   "Y88888P"   "Y8888P"  888       888  "Y88888P"
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  COPYRIGHT (C) 2015 by Filippo Lipparini, Benjamin Stamm, Paolo Gatto        !
!  Eric Cancès, Yvon Maday, Jean-Philip Piquemal, Louis Lagardère and          !
!  Benedetta Mennucci.                                                         !
!                             ALL RIGHT RESERVED.                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! A modular implementation of COSMO using a domain decomposition linear scaling
! strategy.
!
! This code is governed by the LGPL license and abiding by the rules of
! distribution of free software.
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Lesser General Public License for more details.
!
! Users of this code are asked to include the following references in their
! publications:
!
! [1] E. Cancès, Y. Maday, B. Stamm
!     "Domain decomposition for implicit solvation models"
!     J. Chem. Phys. 139, 054111 (2013)
!
! [2] F. Lipparini, B. Stamm, E. Cancès, Y. Maday, B. Mennucci
!     "Fast Domain Decomposition Algorithm for Continuum Solvation Models:
!      Energy and First Derivatives"
!     J. Chem. Theory Comput. 9, 3637–3648 (2013)
!
! Also, include one of the three following reference depending on whether you
! use this code in conjunction with a QM [3], Semiempirical [4] or Classical [5]
! description of the solute:
!
! [3] F. Lipparini, G. Scalmani, L. Lagardère, B. Stamm, E. Cancès, Y. Maday,
!     J.-P. Piquemal, M. J. Frisch, B. Mennucci
!     "Quantum, classical, and hybrid QM/MM calculations in solution: General
!      implementation of the ddCOSMO linear scaling strategy"
!     J. Chem. Phys. 141, 184108 (2014)
!     (for quantum mechanical models)
!
! [4] F. Lipparini, L. Lagardère, G. Scalmani, B. Stamm, E. Cancès, Y. Maday,
!     J.-P. Piquemal, M. J. Frisch, B. Mennucci
!     "Quantum Calculations in Solution for Large to Very Large Molecules:
!      A New Linear Scaling QM/Continuum Approach"
!     J. Phys. Chem. Lett. 5, 953-958 (2014)
!     (for semiempirical models)
!
! [5] F. Lipparini, L. Lagardère, C. Raynaud, B. Stamm, E. Cancès, B. Mennucci
!     M. Schnieders, P. Ren, Y. Maday, J.-P. Piquemal
!     "Polarizable Molecular Dynamics in a Polarizable Continuum Solvent"
!     J. Chem. Theory Comput. 11, 623-634 (2015)
!     (for classical models, including polarizable force fields
!
! The users of this code should also include the appropriate reference to the
! COSMO model. This distribution includes the routines to generate lebedev
! grids by D. Laikov and C. van Wuellen, as publicly available on CCL. If the
! routines are used, the following reference should also be included:
!
! [6] V.I. Lebedev, and D.N. Laikov
!     "A quadrature formula for the sphere of the 131st
!      algebraic order of accuracy"
!     Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!
! Written by Filippo Lipparini, October 2015.
! Adapted by Sebastian Ehlert, June 2020.
module xtb_solv_ddcosmo_core
   implicit none
   private

   public :: TDomainDecomposition, TDomainDecompositionInput, initDomainDecomposition
   public :: hsnorm, calcv, intrhs, prtsph, adjrhs, wghpot, ddupdate, fdoka, fdokb, fdoga


   integer, parameter :: wp = selected_real_kind(15)
   integer, parameter :: ndiis=25, iout=6, nngmax=100
   real(wp), parameter :: zero=0._wp, pt5=0.5_wp, one=1._wp, two=2._wp, four=4._wp
   real(wp), parameter :: se = -1.0_wp
   real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)
   real(wp), parameter :: sq2 = sqrt(2.0_wp)


   !> Quantities contained in control file
   type :: TDomainDecompositionInput

      !> Max angular momentum of spherical harmonics basis
      integer :: lmax

      !> Threshold for iterative solver
      real(wp) :: conv

      !> Regularization parameters
      real(wp) :: eta

   end type TDomainDecompositionInput


   !> Cosmo calculator
   type :: TDomainDecomposition

      !> Printing flag
      integer :: iprint

      !> Max angular momentum of spherical harmonics basis
      integer :: lmax

      !> Desired number of Lebedev integration points
      integer :: ngrid

      !> Threshold for iterative solver
      real(wp) :: conv

      !> 1) compute forces ; 0) do not compute forces
      integer :: igrad

      !> Regularization parameters
      real(wp) :: eta

      !> Number of spheres/atoms
      integer :: nat

      !> Number of integration points on cavity's boundary
      integer :: ncav

      !> Number of basis functions, i.e., spherical harmonics
      integer :: nylm

      !> Workspaces
      logical :: grad
      integer, allocatable :: inl(:), nl(:)
      real(wp), allocatable :: rvdw(:), xyz(:, :), ccav(:, :)
      real(wp), allocatable :: w(:), grid(:, :), basis(:, :)
      real(wp), allocatable :: fact(:), facl(:), facs(:)
      real(wp), allocatable :: fi(:, :), ui(:, :), zi(:, :, :)

   end type TDomainDecomposition


contains


!> Initialize data structure
!
!  allocate the various arrays needed for ddcosmo,
!  assemble the cavity and the various associated geometrical quantities.
subroutine initDomainDecomposition(self, input, rvdw, wang, grid)

   character(len=*), parameter :: source = 'ddcosmo_core::initDomainDecomposition'
   type(TDomainDecomposition), intent(out) :: self
   type(TDomainDecompositionInput), intent(in) :: input
   real(wp), intent(in) :: rvdw(:)
   real(wp), intent(in) :: wang(:), grid(:, :)

   integer :: iat, jat, i, ii, lnl, l, ind, m, igrid, inear, jnear
   integer :: istatus
   real(wp) :: fac, fl, ffl, fnorm, d2, r2, v(3), vv, t, xt, swthr

   real(wp), allocatable :: vcos(:), vsin(:), vplm(:)

   self%nat = size(rvdw)
   self%ngrid = size(wang)
   self%iprint = 0
   self%lmax = input%lmax
   self%conv = input%conv
   self%eta = input%eta
   self%igrad = 1

   ! allocate:
   self%grad   = self%igrad /= 0
   self%nylm   = (self%lmax+1)*(self%lmax+1)
   allocate(self%rvdw(self%nat), self%xyz(3, self%nat), self%w(self%ngrid), &
      & self%grid(3, self%ngrid), self%basis(self%nylm, self%ngrid), &
      & self%inl(self%nat+1), self%nl(self%nat*nngmax), &
      & self%fi(self%ngrid, self%nat), self%ui(self%ngrid, self%nat), &
      & self%fact(max(2*self%lmax+1, 2)), self%facl(self%nylm), &
      & self%facs(self%nylm))
   if (self%grad) allocate(self%zi(3, self%ngrid, self%nat))

   ! precompute various quantities (diagonal blocks of ddcosmo matrix,
   ! normalization factors for spherical harmonics, factorials)
   self%fact(1) = 1.0_wp
   self%fact(2) = 1.0_wp
   do i = 3, 2*self%lmax + 1
      self%fact(i) = real(i-1, wp)*self%fact(i-1)
   end do

   do l = 0, self%lmax
      ind = l*l + l + 1
      fl  = (2.0_wp*real(l, wp) + 1.0_wp)/(4.0_wp*pi)
      ffl = sqrt(fl)
      self%facl(ind-l:ind+l) = fl
      self%facs(ind) = ffl
      do m = 1, l
         fnorm = sq2*ffl*sqrt(self%fact(l-m+1)/self%fact(l+m+1))
         if (mod(m, 2) == 1) fnorm = -fnorm
         self%facs(ind+m) = fnorm
         self%facs(ind-m) = fnorm
      end do
   end do

   ! copy the angular grid:
   self%grid = grid
   ! scaling because the weights are normalised
   self%w = wang * four*pi

   ! set the centers and radii of the spheres:
   self%rvdw(:) = rvdw

   ! build a basis of spherical harmonics at the griwpoints:
   allocate(vplm(self%nylm), vcos(self%lmax+1), vsin(self%lmax+1))
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self) private(i, vplm, vcos, vsin)
   do i = 1, self%ngrid
      call ylmbas(self, self%grid(:, i), self%basis(:, i), vplm, vcos, vsin)
   end do
   deallocate(vplm, vcos, vsin)

   if (self%iprint >= 4) then
      call prtsph(self, 'facs', 1, 0, self%facs)
      call prtsph(self, 'facl', 1, 0, self%facl)
      call prtsph(self, 'basis', self%ngrid, 0, self%basis)
      call ptcart(self, 'grid', 3, 0, self%grid)
      call ptcart(self, 'weights', 1, 0, self%w)
   end if

end subroutine initDomainDecomposition


subroutine ddupdate(self, xyz)

   type(TDomainDecomposition), intent(inout) :: self

   real(wp), intent(in) :: xyz(:, :)

   self%xyz(:, :) = xyz

   call mknnl(self)
   call mkfiui(self)
   call mkccav(self)

end subroutine ddupdate


! build neighbors list (CSR format)
!
!  \\  jat  |
! iat   \\  |  1   2   3   4   5   6
! -----------------------------------
!         1 |      x       x   x
!         2 |  x       x       x   x
!         3 |      x       x       x
!         4 |  x       x       x   x
!         5 |  x   x       x
!         6 |      x   x   x
!
!
!  inl =  1,       4,          8,      11,         15,      18, 21        pointer to 1st neighbor
!
!         |        |           |        |           |        |
!         v        V           V        V           V        V
!
!         1| 2| 3| 4| 5| 6| 7| 8| 9|10|11|12|13|14|15|16|17|18|19|20
!
!  nl  =  2, 4, 5, 1, 3, 5, 6, 2, 4, 6, 1, 3, 5, 6, 1, 2, 4, 2, 3, 4     neighbors list
subroutine mknnl(self)

   type(TDomainDecomposition), intent(inout) :: self

   integer :: ii, lnl, iat, jat
   real(wp) :: v(3), d2, r2

   character(len=*), parameter :: f1000 = "(t3, 'neighbours of sphere ', i6)"
   character(len=*), parameter :: f1010 = "(t5, 12i6)"

   ! index of nl
   ii  = 1
   lnl = 0
   do iat = 1, self%nat
      self%inl(iat) = lnl + 1
      do jat = 1, self%nat
         if (iat /= jat) then
            v(:) = self%xyz(:, iat) - self%xyz(:, jat)
            d2 = v(1)**2 + v(2)**2 + v(3)**2
            r2 = (self%rvdw(iat) + self%rvdw(jat))**2
            if (d2 <= r2) then
               self%nl(ii) = jat
               ii  = ii + 1
               lnl = lnl + 1
            end if
         end if
      end do
   end do
   self%inl(self%nat+1) = lnl+1

   if (self%iprint >= 4) then
      write(iout, *) '   inl:'
      write(iout, '(10i6)') self%inl(1:self%nat+1)
      write(iout, *)
      do iat = 1, self%nat
         write(iout, f1000) iat
         write(iout, f1010) self%nl(self%inl(iat):self%inl(iat+1)-1)
      end do
      write(iout, *)
   end if

end subroutine mknnl


! Define :
!
!   N_i = list of neighbors of i-sphere [ excluding i-sphere ]
!
!            | r_i + rho_i s_n - r_j |
!   t_n^ij = -------------------------
!                      rho_j
!
!   fi(n, i) =    sum    \chi(t_n^ij)
!             j \in N_i
!
! Notice that the derivative of fi(n, i) wrt to r_k is (possibly) nonzero
! when either k = i, or k \in N_j .
subroutine mkfiui(self)

   type(TDomainDecomposition), intent(inout) :: self

   integer :: iat, i, ii, jat
   real(wp) :: v(3), vv, t, xt, swthr, fac


   ! build arrays fi, ui, zi
   self%fi = 0.0_wp
   self%ui = 0.0_wp
   if (self%grad) self%zi = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(self) private(iat, i, ii, jat, v, vv, t, xt, swthr, fac)
   ! loop over spheres
   do iat = 1, self%nat

      ! loop over integration points
      do i = 1, self%ngrid

         ! loop over neighbors of i-sphere
         do ii = self%inl(iat), self%inl(iat+1)-1

            ! neighbor's number
            jat = self%nl(ii)

            ! compute t_n^ij
            v(:) = self%xyz(:, iat) + self%rvdw(iat)*self%grid(:, i) - self%xyz(:, jat)
            vv = v(1)**2 + v(2)**2 + v(3)**2
            vv = sqrt(vv)
            t = vv/self%rvdw(jat)

            ! compute \chi(t_n^ij)
            xt = fsw(t, se, self%eta)

            ! upper bound of switch region
            swthr = 1.0_wp + (se + 1._wp)*self%eta / 2._wp

            ! t_n^ij belongs to switch region
            if (self%grad .and. (t < swthr .and. t > swthr-self%eta)) then
               fac = dfsw(t, se, self%eta) / self%rvdw(jat)

               ! accumulate for zi
               self%zi(:, i, iat) = self%zi(:, i, iat) + fac*v(:)/vv
            end if

            ! accumulate for fi
            self%fi(i, iat) = self%fi(i, iat) + xt
         end do

         ! compute ui
         if (self%fi(i, iat) <= 1.0_wp)  self%ui(i, iat) = 1.0_wp - self%fi(i, iat)
      end do
   end do

   if (self%iprint >= 4) then
      call ptcart(self, 'fi', self%nat, 0, self%fi)
      call ptcart(self, 'ui', self%nat, 0, self%ui)
   end if

end subroutine mkfiui


! build cavity array
subroutine mkccav(self)

   type(TDomainDecomposition), intent(inout) :: self

   integer :: iat, i, ii
   character(len=*), parameter :: f1100 = "(t3, i8, 3f14.6)"

   if (allocated(self%ccav)) deallocate(self%ccav)

   ! initialize number of cavity points
   self%ncav=0

   ! loop over spheres
   do iat = 1, self%nat

      ! loop over integration points
      do i = 1, self%ngrid

         ! positive contribution from integration point
         if (self%ui(i, iat) > 0.0_wp) then

            ! accumulate
            self%ncav = self%ncav + 1

         end if
      end do
   end do

   ! allocate cavity array
   allocate(self%ccav(3, self%ncav))

   ! initialize cavity array index
   ii = 0

   ! loop over spheres
   do iat = 1, self%nat

      ! loop over integration points
      do i = 1, self%ngrid

         ! positive contribution from integration point
         if (self%ui(i, iat) > 0.0_wp) then

            ! advance cavity array index
            ii = ii + 1

            ! store point
            self%ccav(:, ii) = self%xyz(:, iat) + self%rvdw(iat)*self%grid(:, i)

         end if
      end do
   end do

   if (self%iprint >= 4) then
      write(iout, *) '   external cavity points:'
      do ii = 1, self%ncav
         write(iout, f1100) ii, self%ccav(:, ii)
      end do
      write(iout, *)
   end if

end subroutine mkccav


!> Scalar product
pure function sprod(n, u, v)

   integer, intent(in) :: n
   real(wp), intent(in) :: u(n), v(n)
   real(wp) :: sprod

   integer :: i
   real(wp) :: ss

   ! initialize
   ss = 0.0_wp

   do i = 1, n
      ! accumulate
      ss = ss + u(i)*v(i)
   end do

   ! redirect
   sprod = ss

end function sprod


!> Switch function (5th degree polynomial)
elemental function fsw(t, s, eta)

   real(wp) :: fsw
   real(wp), intent(in) :: t
   real(wp), intent(in) :: s
   real(wp), intent(in) :: eta

   real(wp) :: a, b, flow, x
   real(wp), parameter :: f6=6.0_wp, f10=10._wp, f12=12._wp, f15=15._wp

   ! shift :
   ! s =  0   =>   t - eta/2  [ CENTERED ]
   ! s =  1   =>   t - eta    [ EXTERIOR ]
   ! s = -1   =>   t          [ INTERIOR ]

   ! apply shift
   x = t - (s + 1._wp)*eta / 2._wp

   ! lower bound of switch region
   flow = 1.0_wp - eta

   ! define switch function \chi
   if (x >= 1.0_wp) then
      fsw = 0.0_wp
   elseif (x <= flow) then
      fsw = 1.0_wp
   else
      a = f15*eta - f12
      b = f10*eta*eta - f15*eta + f6
      fsw = ((x-1.0_wp)*(x-1.0_wp)*(1.0_wp-x)*(f6*x*x + a*x + b)) / (eta**5)
   end if

end function fsw


!> First derivative of switch function
elemental function dfsw(t, s, eta)

   real(wp), intent(in) :: t
   real(wp), intent(in) :: s
   real(wp), intent(in) :: eta
   real(wp) :: dfsw

   real(wp) :: flow, x
   real(wp), parameter :: f30=30.0_wp

   ! shift :
   ! s =  0   =>   t - eta/2  [ CENTERED ]
   ! s =  1   =>   t - eta    [ EXTERIOR ]
   ! s = -1   =>   t          [ INTERIOR ]
   !
   ! apply shift
   x = t - (s + 1._wp)*eta / 2._wp

   ! lower bound of switch region
   flow = 1.0_wp - eta

   ! define derivative of switch function \chi
   if (x >= 1.0_wp) then
      dfsw = 0.0_wp
   elseif (x <= flow) then
      dfsw = 0.0_wp
   else
      dfsw = f30*(1.0_wp-x)*(x-1.0_wp)*(x-1.0_wp+eta)*(x-1.0_wp+eta)/(eta**5)
   end if

end function dfsw


!> Dump an array (ngrid, ncol) or just a column.
subroutine ptcart(self, label, ncol, icol, x)

   type(TDomainDecomposition), intent(in) :: self
   character(len=*), intent(in) :: label
   integer, intent(in) :: ncol, icol
   real(wp), intent(in) :: x(self%ngrid, ncol)

   integer :: ig, noff, nprt, ic, j
   character(len=*), parameter :: f1000 = "(1x, i5, f14.8)"
   character(len=*), parameter :: f1010 = "(6x, 5i14)"
   character(len=*), parameter :: f1020 = "(1x, i5, 5f14.8)"

   ! print header :
   if (ncol == 1) then
      write (iout, '(3x, a, 1x, "(column ", i4, ")")') label, icol
   else
      write (iout, '(3x, a)') label
   end if

   ! print entries :
   if (ncol == 1) then
      do ig = 1, self%ngrid
         write(iout, f1000) ig, x(ig, 1)
      end do

   else
      noff = mod (ncol, 5)
      nprt = max(ncol - noff, 0)
      do ic = 1, nprt, 5
         write(iout, f1010) (j, j = ic, ic+4)
         do ig = 1, self%ngrid
            write(iout, f1020) ig, x(ig, ic:ic+4)
         end do
      end do
      write (iout, f1010) (j, j = nprt+1, nprt+noff)
      do ig = 1, self%ngrid
         write(iout, f1020) ig, x(ig, nprt+1:nprt+noff)
      end do
   end if

end subroutine ptcart


!> Dump an array (nylm, ncol) or just a column.
subroutine prtsph(self, label, ncol, icol, x)

   type(TDomainDecomposition), intent(in) :: self
   character (len=*), intent(in) :: label
   integer, intent(in) :: ncol, icol
   real(wp), intent(in) :: x(self%nylm, ncol)

   integer :: l, m, ind, noff, nprt, ic, j
   character(len=*), parameter :: f1000 = "(1x, i3, i4, f14.8)"
   character(len=*), parameter :: f1010 = "(8x, 5i14)"
   character(len=*), parameter :: f1020 = "(1x, i3, i4, 5f14.8)"

   ! print header :
   if (ncol == 1) then
      write (iout, '(3x, a, 1x, "(column ", i4, ")")') label, icol
   else
      write (iout, '(3x, a)') label
   end if

   ! print entries :
   if (ncol == 1) then
      do l = 0, self%lmax
         ind = l*l + l + 1
         do m = -l, l
            write(iout, f1000) l, m, x(ind+m, 1)
         end do
      end do

   else
      noff = mod(ncol, 5)
      nprt = max(ncol - noff, 0)
      do ic = 1, nprt, 5
         write(iout, f1010) (j, j = ic, ic+4)
         do l = 0, self%lmax
            ind = l*l + l + 1
            do m = -l, l
               write(iout, f1020) l, m, x(ind+m, ic:ic+4)
            end do
         end do
      end do
      write (iout, f1010) (j, j = nprt+1, nprt+noff)
      do l = 0, self%lmax
         ind = l*l + l + 1
         do m = -l, l
            write(iout, f1020) l, m, x(ind+m, nprt+1:nprt+noff)
         end do
      end do
   end if

end subroutine prtsph


!> Integrate against spherical harmonics
subroutine intrhs(self, iat, x, xlm)

   type(TDomainDecomposition), intent(in) :: self
   integer, intent(in) :: iat
   real(wp), intent(in) :: x(:) ! [self%ngrid]
   real(wp), intent(inout) :: xlm(:) ! [self%nylm]

   integer :: ig

   ! initialize
   xlm = 0.0_wp

   ! accumulate over integration points
   do ig = 1, self%ngrid
      xlm = xlm + self%basis(:, ig)*self%w(ig)*x(ig)
   end do

   ! printing
   if (self%iprint >= 5) then
      call ptcart(self, 'pot', 1, iat, x)
      call prtsph(self, 'vlm', 1, iat, xlm)
   end if

end subroutine intrhs


!> Compute spherical harmonics
pure subroutine ylmbas(self, x, basloc, vplm, vcos, vsin)

   type(TDomainDecomposition), intent(in) :: self
   real(wp), intent(in) :: x(:) ! [3]
   real(wp), intent(out) :: basloc(:), vplm(:) ! [self%nylm]
   real(wp), intent(out) :: vcos(:), vsin(:) ! [self%lmax+1]

   integer :: l, m, ind
   real(wp) :: cthe, sthe, cphi, sphi, plm

   ! get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
   ! coordinates of x.

   ! evaluate cos(theta) ; sin(theta)
   cthe = x(3)
   sthe = sqrt(1.0_wp - cthe*cthe)

   ! evalutate cos(phi) ; sin(phi)
   if (sthe /= 0.0_wp) then
      cphi = x(1)/sthe
      sphi = x(2)/sthe
   else
      cphi = 0.0_wp
      sphi = 0.0_wp
   end if

   ! evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is
   ! pointless if z = 1, as the only non vanishing terms will be the
   ! ones with m=0.
   if(sthe /= 0.0_wp) then
      call trgev(self, cphi, sphi, vcos, vsin)
   else
      vcos = 1.0_wp
      vsin = 0.0_wp
   end if

   ! evaluate the generalized legendre polynomials
   call polleg(self, cthe, sthe, vplm)

   ! now build the spherical harmonics. we will distinguish m=0,
   ! m>0 and m<0:
   do l = 0, self%lmax
      ind = l**2 + l + 1

      ! m = 0
      basloc(ind) = self%facs(ind)*vplm(ind)

      do m = 1, l

         plm = vplm(ind+m)

         ! m > 0
         basloc(ind+m) = self%facs(ind+m)*plm*vcos(m+1)

         ! m < 0
         basloc(ind-m) = self%facs(ind-m)*plm*vsin(m+1)

      end do
   end do

end subroutine ylmbas


!> Compute first derivatives of spherical harmonics
pure subroutine dbasis(self, x, basloc, dbsloc, vplm, vcos, vsin)

   type(TDomainDecomposition), intent(in) :: self
   real(wp), intent(in) :: x(:) ! [3]
   real(wp), intent(inout) :: basloc(:), vplm(:) ! [self%nylm]
   real(wp), intent(inout) :: dbsloc(:, :) ! [3, self%nylm]
   real(wp), intent(inout) :: vcos(:), vsin(:) ! [self%lmax+1]

   integer :: l, m, ind, VC, VS
   real(wp) :: cthe, sthe, cphi, sphi, plm, fln, pp1, pm1, pp
   real(wp) :: et(3), ep(3)

   ! get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
   ! coordinates of x.
   cthe = x(3)
   sthe = sqrt(1.0_wp - cthe*cthe)

   if (sthe /= 0.0_wp) then
      ! not (NORTH or SOUTH pole)
      cphi = x(1)/sthe
      sphi = x(2)/sthe
   else
      ! NORTH or SOUTH pole
      cphi = 1.0_wp
      sphi = 0.0_wp
   end if

   ! evaluate the derivatives of theta and phi:
   et(1) = cthe*cphi
   et(2) = cthe*sphi
   et(3) = -sthe

   if(sthe /= 0.0_wp) then
      ! not (NORTH or SOUTH pole)
      ep(1) = -sphi/sthe
      ep(2) = cphi/sthe
      ep(3) = 0.0_wp
   else
      ! NORTH or SOUTH pole
      ep(1) = 0.0_wp
      ep(2) = 1.0_wp
      ep(3) = 0.0_wp
   end if

   ! evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is
   ! pointless if z = 1, as the only non vanishing terms will be the
   ! ones with m=0.

   if (sthe /= 0.0_wp) then
      ! not (NORTH or SOUTH pole)
      call trgev(self, cphi, sphi, vcos, vsin)
   else
      ! NORTH or SOUTH pole
      vcos(:) = 1.0_wp
      vsin(:) = 0.0_wp
   end if
   VC=0.0_wp
   VS=cthe

   ! evaluate the generalized legendre polynomials.
   call polleg(self, cthe, sthe, vplm)

   ! now build the spherical harmonics. we will distinguish m=0,
   ! m>0 and m<0:

   basloc = 0.0_wp
   dbsloc = 0.0_wp
   do l = 0, self%lmax
      ind = l*l + l + 1
      ! m = 0
      fln = self%facs(ind)
      basloc(ind) = fln*vplm(ind)
      if (l > 0) then
         dbsloc(:, ind) = fln*vplm(ind+1)*et(:)
      else
         dbsloc(:, ind) = 0.0_wp
      end if
      !$omp simd
      do m = 1, l
         fln = self%facs(ind+m)
         plm = fln*vplm(ind+m)
         pp1 = 0.0_wp
         if (m < l) pp1 = -0.5_wp*vplm(ind+m+1)
         pm1 = 0.5_wp*(real(l+m, wp)*real(l-m+1, wp)*vplm(ind+m-1))
         pp  = pp1 + pm1

         ! m > 0
         basloc(ind+m) = plm*vcos(m+1)

         if (sthe /= 0.0_wp) then
            ! not (NORTH or SOUTH pole)
            dbsloc(:, ind+m) = -fln*pp*vcos(m+1)*et(:) - real(m, wp)*plm*vsin(m+1)*ep(:)
         else
            ! NORTH or SOUTH pole
            dbsloc(:, ind+m) = -fln*pp*vcos(m+1)*et(:) - fln*pp*ep(:)*VC
         end if

         ! m < 0
         basloc(ind-m) = plm*vsin(m+1)

         if (sthe /= 0.0_wp) then
            ! not (NORTH or SOUTH pole)
            dbsloc(:, ind-m) = -fln*pp*vsin(m+1)*et(:) + real(m, wp)*plm*vcos(m+1)*ep(:)
         else
            ! NORTH or SOUTH pole
            dbsloc(:, ind-m) = -fln*pp*vsin(m+1)*et(:) - fln*pp*ep(:)*VS
         end if

      end do
   end do

end subroutine dbasis


!> Compute the l, m associated legendre polynomial for -1 <= x <= 1 using the
!  recurrence formula
!
!  (l-m)p(l, m) = x(2l-1)p(l-1, m) - (l+m-1)p(l-2, m)
pure subroutine polleg(self, x, y, plm)

   type(TDomainDecomposition), intent(in) :: self
   real(wp), intent(in) :: x, y
   real(wp), intent(inout) :: plm(:) ! [self%nylm]

   integer :: m, ind, l, ind2
   real(wp) :: fact, pmm, somx2, pmm1, pmmo, pll, fm, fl

   fact  = 1.0_wp
   pmm   = 1.0_wp
   somx2 = y

   do m = 0, self%lmax
      ind      = (m + 1)*(m + 1)
      plm(ind) = pmm
      if(m == self%lmax) return
      fm = real(m, wp)
      pmm1 = x*(2.0_wp*fm + 1.0_wp)*pmm
      ind2 = ind + 2*m + 2
      plm(ind2) = pmm1
      pmmo = pmm
      do l = m+2, self%lmax
         fl = real(l, wp)
         pll   = (x*(2.0_wp*fl - 1.0_wp)*pmm1 - (fl + fm - 1.0_wp)*pmm)/(fl - fm)
         ind = l*l + l + 1
         plm(ind+m) = pll
         pmm  = pmm1
         pmm1 = pll
      end do
      pmm  = -pmmo*fact*somx2
      fact = fact + 2.0_wp
   end do

end subroutine polleg


!> Service routine for computation of spherical harmonics
pure subroutine trgev(self, x, y, cx, sx)

   type(TDomainDecomposition), intent(in) :: self
   real(wp), intent(in) :: x, y
   real(wp), intent(inout) :: cx(:), sx(:) ! [max((self%lmax+1), 2)]

   integer :: m

   cx(1) = 1.0_wp
   sx(1) = 0.0_wp
   cx(2) = x
   sx(2) = y
   do m = 3, self%lmax+1
      cx(m) = 2.0_wp*x*cx(m-1) - cx(m-2)
      sx(m) = 2.0_wp*x*sx(m-1) - sx(m-2)
   end do

end subroutine trgev


!> Compute
!
!  sum   4pi/(2l+1) t^l * Y_l^m(s) * sigma_l^m
!  l, m
!
!  which is need to compute action of COSMO matrix L.
pure function intmlp(self, t, sigma, basloc)

   type(TDomainDecomposition), intent(in) :: self
   real(wp), intent(in) :: t
   real(wp), intent(in) :: sigma(:), basloc(:) ! [self%nylm]
   real(wp) :: intmlp

   integer :: l, ind
   real(wp) :: tt, ss, fac

   ! initialize t^l
   tt = 1.0_wp

   ! initialize
   ss = 0.0_wp

   ! loop over l
   do l = 0, self%lmax

      ind = l*l + l + 1

      ! update factor 4pi / (2l+1) * t^l
      fac = tt / self%facl(ind)

      ! contract over l, m and accumulate
      ss = ss + fac * dot_product(basloc(ind-l:ind+l), sigma(ind-l:ind+l))

      ! update t^l
      tt = tt*t
   end do

   ! redirect
   intmlp = ss

end function intmlp


!> Weigh potential at cavity points by characteristic function "ui"
subroutine wghpot(self, phi, g)

   type(TDomainDecomposition), intent(in) :: self
   real(wp), intent(in)  :: phi(:) ! [self%ncav]
   real(wp), intent(out) :: g(:, :) ! [self%ngrid, self%nat]

   integer :: iat, ig, ic

   !> Initialize
   ic = 0
   g(:, :) = 0._wp

   ! loop over spheres
   do iat = 1, self%nat

      ! loop over points
      do ig = 1, self%ngrid

         ! nonzero contribution from point
         if (self%ui(ig, iat) /= 0.0_wp) then

            ! advance cavity point counter
            ic = ic + 1

            ! weight by (negative) characteristic function
            g(ig, iat) = -self%ui(ig, iat) * phi(ic)

         end if

      end do
   end do

end subroutine wghpot


!> Compute H-norm
pure subroutine hsnorm(self, u, unorm)

   type(TDomainDecomposition), intent(in) :: self
   real(wp), intent(in) :: u(:) !< [nylm]
   real(wp), intent(inout) :: unorm

   integer :: l, m, ind
   real(wp) :: fac

   ! initialize
   unorm = 0.0_wp

   ! loop over l
   do l = 0, self%lmax

      ! first index associated to l
      ind = l*l + l + 1

      ! scaling factor
      fac = 1.0_wp/(1.0_wp + real(l, wp))

      ! loop over m
      do m = -l, l

         ! accumulate
         unorm = unorm + fac*u(ind+m)*u(ind+m)

      end do
   end do

   ! the much neglected square root
   unorm = sqrt(unorm)

end subroutine hsnorm


!> Compute
!
!   v_l^m = v_l^m +
!
!               4 pi           l
!     sum  sum  ---- (t_n^ji)  Y_l^m(s_n^ji) W_n^ji [ \xi_j ]_n
!      j    n   2l+1
!
! which is related to the action of the adjont COSMO matrix L^* in the following
! way. Set
!
!   [ \xi_j ]_n = sum  Y_l^m(s_n) [ s_j ]_l^m
!                 l, m
!
! then
!
!   v_l^m = -   sum    (L^*)_ij s_j
!             j \ne i
!
! The auxiliary quantity [ \xi_j ]_l^m needs to be computed explicitly.
pure subroutine adjrhs(self, iat, xi, vlm, basloc, vplm, vcos, vsin)

   type(TDomainDecomposition), intent(in) :: self
   integer, intent(in) :: iat
   real(wp), intent(in) :: xi(:, :) ! [self%ngrid, self%nat]
   real(wp), intent(inout) :: vlm(:) ! [self%nylm]
   real(wp), intent(inout) :: basloc(:), vplm(:) ! [self%nylm]
   real(wp), intent(inout) :: vcos(:), vsin(:) ! [self%lmax+1]

   integer :: ij, jat, ig, l, ind, m
   real(wp) :: vji(3), vvji, tji, sji(3), xji, oji, fac, ffac, t

   ! loop over neighbors of i-sphere
   do ij = self%inl(iat), self%inl(iat+1)-1

      ! j-sphere is neighbor
      jat = self%nl(ij)

      ! loop over integration points
      do ig = 1, self%ngrid

         ! compute t_n^ji = | r_j + \rho_j s_n - r_i | / \rho_i
         vji  = self%xyz(:, jat) + self%rvdw(jat)*self%grid(:, ig) - self%xyz(:, iat)
         vvji = vji(1)**2 + vji(2)**2 + vji(3)**2
         vvji = sqrt(vvji)
         tji  = vvji/self%rvdw(iat)

         ! point is INSIDE i-sphere (+ transition layer)
         if (tji < (1.0_wp + (se+1.0_wp)/2.0_wp*self%eta)) then

            ! compute s_n^ji
            sji = vji/vvji

            ! compute \chi(t_n^ji)
            xji = fsw(tji, se, self%eta)

            ! compute W_n^ji
            if (self%fi(ig, jat) > 1.0_wp) then
               oji = xji/self%fi(ig, jat)
            else
               oji = xji
            end if

            ! compute Y_l^m(s_n^ji)
            call ylmbas(self, sji, basloc, vplm, vcos, vsin)

            ! initialize (t_n^ji)^l
            t = 1.0_wp

            ! compute w_n * xi(n, j) * W_n^ji
            fac = self%w(ig) * xi(ig, jat) * oji

            ! loop over l
            do l = 0, self%lmax
               ind  = l*l + l + 1

               ! compute 4pi / (2l+1) * (t_n^ji)^l * w_n * xi(n, j) * W_n^ji
               ffac = fac*t/self%facl(ind)

               ! loop over m
               do m = -l, l
                  vlm(ind+m) = vlm(ind+m) + ffac*basloc(ind+m)
               end do

               ! update (t_n^ji)^l
               t = t*tji
            end do

         end if
      end do
   end do

end subroutine adjrhs


subroutine header(self)

   type(TDomainDecomposition), intent(in) :: self

   1000 format(/, &
      ' An implementation of COSMO using a domain decomposition linear scaling strategy.', /)
   1010 format(' Parameters:', /, &
      '   number of grid points:                  ', 8x, i8, /,   &
      '   number of spheres:                      ', 8x, i8, /,   &
      '   lmax for the spherical harmonics basis: ', 8x, i8, /,   &
      '   convergence threshold:                  ', 8x, d8.1, /, &
      '   regularization parameter (eta):         ', 8x, f8.3/)

   if (self%iprint > 0) then

      write(iout, 1000)
      write(iout, 1010) self%ngrid, self%nat, self%lmax, self%conv, self%eta

      if (self%igrad == 1)  write(iout, 1013)
      1013   format(' Compute forces.'//)

   end if

end subroutine header


!> Compute the first part of <S, L^(x)X>
pure subroutine fdoka(self, iat, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, fx)

   type(TDomainDecomposition), intent(in) :: self
   integer, intent(in) :: iat
   real(wp), intent(in) :: sigma(:, :) ! [self%nylm, self%nat]
   real(wp), intent(in) :: xi(:) ! [self%ngrid]
   real(wp), intent(inout) :: basloc(:), vplm(:) ! [self%nylm]
   real(wp), intent(inout) :: dbsloc(:, :) ! [3, self%nylm]
   real(wp), intent(inout) :: vcos(:), vsin(:) ! [self%lmax+1]
   real(wp), intent(inout) :: fx(:) ! [3]

   integer :: ig, ij, jat, l, ind, m
   real(wp) :: vvij, tij, xij, oij, t, fac, fl, f1, f2, f3, beta, tlow, thigh
   real(wp) :: vij(3), sij(3), alp(3), va(3)

   tlow  = 1.0_wp - 0.5_wp*(1.0_wp - se)*self%eta
   thigh = 1.0_wp + 0.5_wp*(1.0_wp + se)*self%eta

   do ig = 1, self%ngrid
      va = 0.0_wp
      do ij = self%inl(iat), self%inl(iat+1) - 1
         jat = self%nl(ij)
         vij = self%xyz(:, iat) + self%rvdw(iat)*self%grid(:, ig) - self%xyz(:, jat)
         vvij = vij(1)**2 + vij(2)**2 + vij(3)**2
         vvij = sqrt(vvij)
         tij = vvij/self%rvdw(jat)

         if (tij >= thigh) cycle

         sij = vij/vvij
         call dbasis(self, sij, basloc, dbsloc, vplm, vcos, vsin)
         alp = 0.0_wp
         t = 1.0_wp
         do l = 1, self%lmax
            ind = l*l + l + 1
            fl = real(l, wp)
            fac = t/self%facl(ind)
            do m = -l, l
               f2 = fac*sigma(ind+m, jat)
               f1 = f2*fl*basloc(ind+m)
               alp(:) = alp(:) + f1*sij(:) + f2*dbsloc(:, ind+m)
            end do
            t = t*tij
         end do
         beta = intmlp(self, tij, sigma(:, jat), basloc)
         xij = fsw(tij, se, self%eta)
         if (self%fi(ig, iat) > 1.0_wp) then
            oij = xij/self%fi(ig, iat)
            f2  = -oij/self%fi(ig, iat)
         else
            oij = xij
            f2  = 0.0_wp
         end if
         f1 = oij/self%rvdw(jat)
         va(:) = va(:) + f1*alp(:) + beta*f2*self%zi(:, ig, iat)
         if (tij > tlow) then
            f3 = beta*dfsw(tij, se, self%eta)/self%rvdw(jat)
            if (self%fi(ig, iat) > 1.0_wp) f3 = f3/self%fi(ig, iat)
            va(:) = va(:) + f3*sij(:)
         end if
      end do
      fx = fx - self%w(ig)*xi(ig)*va(:)
   end do

end subroutine fdoka


!> Compute the the second part of <S, L^(x)X>
pure subroutine fdokb(self, iat, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, fx)

   type(TDomainDecomposition), intent(in) :: self
   integer, intent(in) :: iat
   real(wp), intent(in) :: sigma(:, :) ! [self%nylm, self%nat]
   real(wp), intent(in) :: xi(:, :) ! [self%ngrid, self%nat]
   real(wp), intent(inout) :: basloc(:), vplm(:) ! [self%nylm]
   real(wp), intent(inout) :: dbsloc(:, :) ! [3, self%nylm]
   real(wp), intent(inout) :: vcos(:), vsin(:) ! [self%lmax+1]
   real(wp), intent(inout) :: fx(:) ! [3]

   integer :: ig, ji, jat, l, ind, m, jk, kat
   logical :: proc
   real(wp) :: vvji, tji, xji, oji, t, fac, fl, f1, f2, beta, di, tlow, thigh
   real(wp) :: b, g1, g2, vvjk, tjk, f, xjk
   real(wp) :: vji(3), sji(3), alp(3), vb(3), vjk(3), sjk(3), vc(3)



   tlow  = 1.0_wp - 0.5_wp*(1.0_wp - se)*self%eta
   thigh = 1.0_wp + 0.5_wp*(1.0_wp + se)*self%eta

   do ig = 1, self%ngrid
      vb = 0.0_wp
      vc = 0.0_wp
      do ji = self%inl(iat), self%inl(iat+1) - 1
         jat = self%nl(ji)
         vji = self%xyz(:, jat) + self%rvdw(jat)*self%grid(:, ig) - self%xyz(:, iat)
         vvji = vji(1)**2 + vji(2)**2 + vji(3)**2
         vvji = sqrt(vvji)
         tji = vvji/self%rvdw(iat)

         if (tji > thigh) cycle

         sji = vji/vvji
         call dbasis(self, sji, basloc, dbsloc, vplm, vcos, vsin)

         alp = 0.0_wp
         t   = 1.0_wp
         do l = 1, self%lmax
            ind = l*l + l + 1
            fl = real(l, wp)
            fac = t/self%facl(ind)
            do m = -l, l
               f2 = fac*sigma(ind+m, iat)
               f1 = f2*fl*basloc(ind+m)
               alp = alp + f1*sji + f2*dbsloc(:, ind+m)
            end do
            t = t*tji
         end do
         xji = fsw(tji, se, self%eta)
         if (self%fi(ig, jat) > 1.0_wp) then
            oji = xji/self%fi(ig, jat)
         else
            oji = xji
         end if
         f1 = oji/self%rvdw(iat)
         vb = vb + f1*alp*xi(ig, jat)
         if (tji > tlow) then
            beta = intmlp(self, tji, sigma(:, iat), basloc)
            if (self%fi(ig, jat) > 1.0_wp) then
               di  = 1.0_wp/self%fi(ig, jat)
               fac = di*xji
               proc = .false.
               b    = 0.0_wp
               do jk = self%inl(jat), self%inl(jat+1) - 1
                  kat = self%nl(jk)
                  vjk = self%xyz(:, jat) + self%rvdw(jat)*self%grid(:, ig) - self%xyz(:, kat)
                  vvjk = vjk(1)**2 + vjk(2)**2 + vjk(3)**2
                  vvjk = sqrt(vvjk)
                  tjk = vvjk/self%rvdw(kat)
                  if (kat /= iat) then
                     if (tjk <= thigh) then
                        proc = .true.
                        sjk = vjk/vvjk
                        call ylmbas(self, sjk, basloc, vplm, vcos, vsin)
                        g1 = intmlp(self, tjk, sigma(:, kat), basloc)
                        xjk = fsw(tjk, se, self%eta)
                        b = b + g1*xjk
                     end if
                  end if
               end do
               if (proc) then
                  g1 = di*di*dfsw(tji, se, self%eta)/self%rvdw(iat)
                  g2 = g1*xi(ig, jat)*b
                  vc = vc + g2*sji
               end if
            else
               di = 1.0_wp
               fac = 0.0_wp
            end if
            f2 = (1.0_wp-fac)*di*dfsw(tji, se, self%eta)/self%rvdw(iat)
            vb = vb + f2*xi(ig, jat)*beta*sji
         end if
      end do
      fx = fx + self%w(ig)*(vb - vc)
   end do

end subroutine fdokb


!> Compute the U^(x)\Phi contribution to <S, g^(x)>
pure subroutine fdoga(self, iat, xi, phi, fx)

   type(TDomainDecomposition), intent(in) :: self
   integer, intent(in) :: iat
   real(wp), intent(in) :: xi(:, :), phi(:, :) ! [self%ngrid, self%nat]
   real(wp), intent(inout) :: fx(:) ! [3]

   integer :: ig, ji, jat
   real(wp) :: vvji, tji, fac, swthr
   real(wp) :: alp(3), vji(3), sji(3)

   do ig = 1, self%ngrid
      alp = 0.0_wp
      if (self%ui(ig, iat) > 0.0_wp .and. self%ui(ig, iat) <  1.0_wp) then
         alp = alp + phi(ig, iat)*xi(ig, iat)*self%zi(:, ig, iat)
      end if
      do ji = self%inl(iat), self%inl(iat+1) - 1
         jat = self%nl(ji)
         vji = self%xyz(:, jat) + self%rvdw(jat)*self%grid(:, ig) - self%xyz(:, iat)
         vvji = vji(1)**2 + vji(2)**2 + vji(3)**2
         vvji = sqrt(vvji)
         tji = vvji/self%rvdw(iat)
         swthr = 1.0_wp + (se + 1._wp)*self%eta / 2._wp
         if (tji < swthr .and. tji > swthr-self%eta .and. self%ui(ig, jat) > 0.0_wp) then
            sji = vji/vvji
            fac = - dfsw(tji, se, self%eta)/self%rvdw(iat)
            alp = alp + fac*phi(ig, jat)*xi(ig, jat)*sji
         end if
      end do
      fx = fx - self%w(ig)*alp
   end do

end subroutine fdoga


!> Auxiliary routine for COSMO action
!  compute
!
!   \Phi(n) =
!
!                       4 pi           l
!     sum  W_n^ij  sum  ---- (t_n^ij)  Y_l^m(s_n^ij) [ \sigma_j ]_l^m
!      j           l, m  2l+1
!
! which is related to the action of the COSMO matrix L in the following
! way :
!
!   -   sum    L_ij \sigma_j = sum  w_n Y_l^m(s_n) \Phi(n)
!     j \ne i                   n
!
! This second step is performed by routine "intrhs".
pure subroutine calcv(self, first, iat, pot, sigma, basloc, vplm, vcos, vsin)

   type(TDomainDecomposition), intent(in) :: self
   logical, intent(in) :: first
   integer, intent(in) :: iat
   real(wp), intent(in) :: sigma(:, :) ! [self%nylm, self%nat]
   real(wp), intent(inout) :: pot(:) ! [self%ngrid]
   real(wp), intent(inout) :: basloc(:) ! [self%nylm]
   real(wp), intent(inout) :: vplm(:) ! [self%nylm]
   real(wp), intent(inout) :: vcos(:) ! [self%lmax+1]
   real(wp), intent(inout) :: vsin(:) ! [self%lmax+1]

   integer :: its, ij, jat
   real(wp) :: vij(3), sij(3)
   real(wp) :: vvij2, vvij, tij, xij, oij, stslm, stslm2, stslm3

   ! initialize
   pot(:) = 0.0_wp

   ! if 1st iteration of Jacobi method, then done!
   if (first)  return

   ! loop over grid points
   do its = 1, self%ngrid

      ! contribution from integration point present
      if (self%ui(its, iat) < 1.0_wp) then

         ! loop over neighbors of i-sphere
         do ij = self%inl(iat), self%inl(iat+1)-1

            ! neighbor is j-sphere
            jat = self%nl(ij)

            ! compute t_n^ij = | r_i + \rho_i s_n - r_j | / \rho_j
            vij  = self%xyz(:, iat) + self%rvdw(iat)*self%grid(:, its) - self%xyz(:, jat)
            vvij2 = vij(1)**2 + vij(2)**2 + vij(3)**2
            vvij = sqrt(vvij2)
            tij  = vvij / self%rvdw(jat)

            ! point is INSIDE j-sphere
            if (tij < 1.0_wp) then

               ! compute s_n^ij = (r_i + \rho_i s_n - r_j) / | ... |
               sij = vij / vvij

               ! compute \chi(t_n^ij)
               xij = fsw(tij, se, self%eta)

               ! compute W_n^ij
               if (self%fi(its, iat) > 1.0_wp) then

                  oij = xij / self%fi(its, iat)

               else

                  oij = xij

               end if

               ! compute Y_l^m(s_n^ij)
               call ylmbas(self, sij, basloc, vplm, vcos, vsin)

               ! accumulate over j, l, m
               pot(its) = pot(its) + oij * intmlp(self, tij, sigma(:, jat), basloc)

            end if
         end do
      end if
   end do

end subroutine calcv


!> Compute the xi intermediate
!
! \xi(n, i) =
!
!  sum w_n U_n^i Y_l^m(s_n) [S_i]_l^m
!  l, m
!
pure subroutine ddmkxi(self, s, xi)

   type(TDomainDecomposition), intent(in) :: self
   real(wp), intent(in) :: s(:, :) ! [self%nylm, self%nat]
   real(wp), intent(inout) :: xi(:) ! [self%ncav]

   integer :: its, iat, ii

   ii = 0
   do iat = 1, self%nat
      do its = 1, self%ngrid
         if (self%ui(its, iat) > 0.0_wp) then
            ii = ii + 1
            xi(ii) = self%w(its)*self%ui(its, iat) &
               & * dot_product(self%basis(:, its), s(:, iat))
         end if
      end do
   end do

end subroutine ddmkxi


end module xtb_solv_ddcosmo_core
