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
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas, only : mctc_scal
   use xtb_type_environment, only : TEnvironment
   implicit none


   integer, parameter :: ndiis=25, iout=6, nngmax=100
   real(wp), parameter :: zero=0._wp, pt5=0.5_wp, one=1._wp, two=2._wp, four=4._wp
   real(wp), parameter :: se = -1.0_wp
   real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)
   real(wp), parameter :: sq2 = sqrt(2.0_wp)


   !> Quantities contained in control file
   type :: TddControl

      !> Printing flag
      integer :: iprint

      !> Max angular momentum of spherical harmonics basis
      integer :: lmax

      !> Desired number of Lebedev integration points
      integer :: ngrid

      !> Threshold for iterative solver (10^-iconv)
      integer :: iconv

      !> 1) compute forces ; 0) do not compute forces
      integer :: igrad

      !> Dielectric constant of the solvent
      real(wp) :: eps

      !> Regularization parameters
      real(wp) :: eta

   end type TddControl


   !> Cosmo calculator
   type, extends(TddControl) :: TddCosmo

      !> Number of spheres/atoms
      integer :: nsph

      !> Number of integration points on cavity's boundary
      integer :: ncav

      !> Number of basis functions, i.e., spherical harmonics
      integer :: nylm

      !> Workspaces
      logical :: grad
      integer, allocatable :: inl(:), nl(:)
      real(wp), allocatable :: rsph(:), csph(:, :), ccav(:, :)
      real(wp), allocatable :: w(:), grid(:, :), basis(:, :)
      real(wp), allocatable :: fact(:), facl(:), facs(:)
      real(wp), allocatable :: fi(:, :), ui(:, :), zi(:, :, :)

   end type TddCosmo


contains


!> Initialize data structure
!
!  allocate the various arrays needed for ddcosmo,
!  assemble the cavity and the various associated geometrical quantities.
subroutine ddinit(self, env, n, xyz, rvdw)
   use xtb_solv_lebedev, only : gridSize, getAngGrid

   character(len=*), parameter :: source = 'solv_ddcosmo_ddinit'
   type(TEnvironment), intent(inout) :: env
   type(TddCosmo), intent(inout) :: self
   integer, intent(in) :: n
   real(wp), intent(in) :: xyz(3, n), rvdw(n)

   integer :: isph, jsph, i, ii, lnl, l, ind, m, igrid, inear, jnear
   integer :: istatus
   real(wp) :: fac, fl, ffl, fnorm, d2, r2, v(3), vv, t, xt, swthr

   real(wp), allocatable :: vcos(:), vsin(:), vplm(:)

   character(len=*), parameter :: f1000 = "(t3, 'neighbours of sphere ', i6)"
   character(len=*), parameter :: f1010 = "(t5, 12i6)"
   character(len=*), parameter :: f1100 = "(t3, i8, 3f14.6)"

   self%nsph = n

   ! choose the lebedev grid with number of points closest to ngrid:
   igrid = 0
   inear = 100000
   do i = 1, size(gridSize)
      jnear = iabs(gridSize(i)-self%ngrid)
      if (jnear.lt.inear) then
         inear = jnear
         igrid = i
      end if
   end do
   self%ngrid = gridSize(igrid)

   ! print a nice header:
   call header(self)

   ! allocate:
   self%grad   = self%igrad.ne.0
   self%nylm   = (self%lmax+1)*(self%lmax+1)
   allocate(self%rsph(self%nsph), self%csph(3, self%nsph), self%w(self%ngrid), &
      & self%grid(3, self%ngrid), self%basis(self%nylm, self%ngrid), &
      & self%inl(self%nsph+1), self%nl(self%nsph*nngmax), &
      & self%fi(self%ngrid, self%nsph), self%ui(self%ngrid, self%nsph), &
      & self%fact(max(2*self%lmax+1, 2)), self%facl(self%nylm), &
      & self%facs(self%nylm), stat=istatus)
   if (istatus.ne.0) then
      call env%error('allocation failed', source)
      return
   end if
   if (self%grad) allocate(self%zi(3, self%ngrid, self%nsph), stat=istatus)
   if (istatus.ne.0) then
      call env%error('allocation failed', source)
      return
   end if

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
         if (mod(m, 2).eq.1) fnorm = -fnorm
         self%facs(ind+m) = fnorm
         self%facs(ind-m) = fnorm
      end do
   end do

   ! set the centers and radii of the spheres:
   self%csph(:, :) = xyz
   self%rsph = rvdw

   ! load a lebedev grid:
   call getAngGrid(igrid, self%grid, self%w, istatus)
   if (istatus.ne.0) then
      call env%error('angular grid generation failed', source)
      return
   end if
   ! scaling because the weights are normalised
   call mctc_scal(self%w, four*pi)

   ! build a basis of spherical harmonics at the gridpoints:
   allocate(vplm(self%nylm), vcos(self%lmax+1), vsin(self%lmax+1), stat=istatus)
   if (istatus.ne.0) then
      call env%error('allocation failed', source)
      return
   end if
   !$omp parallel do default(shared) private(i, vplm, vcos, vsin)
   do i = 1, self%ngrid
      call ylmbas(self, self%grid(:, i), self%basis(:, i), vplm, vcos, vsin)
   end do
   !$omp end parallel do
   deallocate(vplm, vcos, vsin, stat=istatus)
   if (istatus.ne.0) then
      call env%error('deallocation failed', source)
      return
   end if

   if (self%iprint.ge.4) then
      call prtsph(self, 'facs', 1, 0, self%facs)
      call prtsph(self, 'facl', 1, 0, self%facl)
      call prtsph(self, 'basis', self%ngrid, 0, self%basis)
      call ptcart(self, 'grid', 3, 0, self%grid)
      call ptcart(self, 'weights', 1, 0, self%w)
   end if

   ! build neighbors list (CSR format)
   !
   !  \\  jsph |
   ! isph  \\  |  1   2   3   4   5   6
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

   ! index of nl
   ii  = 1
   lnl = 0
   do isph = 1, self%nsph
      self%inl(isph) = lnl + 1
      do jsph = 1, self%nsph
         if (isph.ne.jsph) then
            d2 = (self%csph(1, isph) - self%csph(1, jsph))**2 + (self%csph(2, isph) - self%csph(2, jsph))**2 + (self%csph(3, isph) - self%csph(3, jsph))**2
            r2 = (self%rsph(isph) + self%rsph(jsph))**2
            if (d2.le.r2) then
               self%nl(ii) = jsph
               ii  = ii + 1
               lnl = lnl + 1
            end if
         end if
      end do
   end do
   self%inl(self%nsph+1) = lnl+1

   if (self%iprint.ge.4) then
      write(iout, *) '   inl:'
      write(iout, '(10i6)') self%inl(1:self%nsph+1)
      write(iout, *)
      do isph = 1, self%nsph
         write(iout, f1000) isph
         write(iout, f1010) self%nl(self%inl(isph):self%inl(isph+1)-1)
      end do
      write(iout, *)
   end if

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

   ! build arrays fi, ui, zi
   self%fi = 0.0_wp
   self%ui = 0.0_wp
   if (self%grad) self%zi = 0.0_wp

   !$omp parallel do default(shared) &
   !$omp private(isph, i, ii, jsph, v, vv, t, xt, swthr, fac)
   ! loop over spheres
   do isph = 1, self%nsph

      ! loop over integration points
      do i = 1, self%ngrid

         ! loop over neighbors of i-sphere
         do ii = self%inl(isph), self%inl(isph+1)-1

            ! neighbor's number
            jsph = self%nl(ii)

            ! compute t_n^ij
            v(:) = self%csph(:, isph) + self%rsph(isph)*self%grid(:, i) - self%csph(:, jsph)
            vv = v(1)**2 + v(2)**2 + v(3)**2
            vv = sqrt(vv)
            t = vv/self%rsph(jsph)

            ! compute \chi(t_n^ij)
            xt = fsw(t, se, self%eta)

            ! upper bound of switch region
            swthr = 1.0_wp + (se + 1._wp)*self%eta / 2._wp

            ! t_n^ij belongs to switch region
            if (self%grad .and. (t.lt.swthr .and. t.gt.swthr-self%eta)) then
               fac = dfsw(t, se, self%eta) / self%rsph(jsph)

               ! accumulate for zi
               self%zi(:, i, isph) = self%zi(:, i, isph) + fac*v(:)/vv
            end if

            ! accumulate for fi
            self%fi(i, isph) = self%fi(i, isph) + xt
         end do

         ! compute ui
         if (self%fi(i, isph).le.1.0_wp)  self%ui(i, isph) = 1.0_wp - self%fi(i, isph)
      end do
   end do
   !$omp end parallel do

   if (self%iprint.ge.4) then
      call ptcart(self, 'fi', self%nsph, 0, self%fi)
      call ptcart(self, 'ui', self%nsph, 0, self%ui)
   end if

   ! build cavity array
   ! initialize number of cavity points
   self%ncav=0

   ! loop over spheres
   do isph = 1, self%nsph

      ! loop over integration points
      do i = 1, self%ngrid

         ! positive contribution from integration point
         if (self%ui(i, isph).gt.0.0_wp) then

            ! accumulate
            self%ncav = self%ncav + 1

         end if
      end do
   end do

   ! allocate cavity array
   allocate(self%ccav(3, self%ncav), stat=istatus)
   if (istatus .ne. 0) then
      call env%error('allocation failed', source)
      return
   end if


   ! initialize cavity array index
   ii = 0

   ! loop over spheres
   do isph = 1, self%nsph

      ! loop over integration points
      do i = 1, self%ngrid

         ! positive contribution from integration point
         if (self%ui(i, isph).gt.0.0_wp) then

            ! advance cavity array index
            ii = ii + 1

            ! store point
            self%ccav(:, ii) = self%csph(:, isph) + self%rsph(isph)*self%grid(:, i)

         end if
      end do
   end do

   if (self%iprint.ge.4) then
      write(iout, *) '   external cavity points:'
      do ii = 1, self%ncav
         write(iout, f1100) ii, self%ccav(:, ii)
      end do
      write(iout, *)
   end if

end subroutine ddinit


!> Free data structure
subroutine memfree(self, env)
   character(len=*), parameter :: source = 'solv_ddcosmo_memfree'
   type(TEnvironment), intent(inout) :: env
   type(TddCosmo), intent(inout) :: self
   integer :: istatus, istatus0

   ! initialize deallocation flags
   istatus0 = 0
   istatus = 0

   ! deallocate the arrays
   if (allocated(self%rsph))  deallocate(self%rsph, stat=istatus) ; istatus0 = istatus0 + istatus
   if (allocated(self%csph))  deallocate(self%csph, stat=istatus) ; istatus0 = istatus0 + istatus
   if (allocated(self%ccav))  deallocate(self%ccav, stat=istatus) ; istatus0 = istatus0 + istatus
   if (allocated(self%w))  deallocate(self%w, stat=istatus) ; istatus0 = istatus0 + istatus
   if (allocated(self%grid))  deallocate(self%grid, stat=istatus) ; istatus0 = istatus0 + istatus
   if (allocated(self%basis))  deallocate(self%basis, stat=istatus) ; istatus0 = istatus0 + istatus
   if (allocated(self%inl))  deallocate(self%inl, stat=istatus) ; istatus0 = istatus0 + istatus
   if (allocated(self%nl))  deallocate(self%nl, stat=istatus) ; istatus0 = istatus0 + istatus
   if (allocated(self%fact))  deallocate(self%fact, stat=istatus) ; istatus0 = istatus0 + istatus
   if (allocated(self%facl))  deallocate(self%facl, stat=istatus) ; istatus0 = istatus0 + istatus
   if (allocated(self%facs))  deallocate(self%facs, stat=istatus) ; istatus0 = istatus0 + istatus
   if (allocated(self%ui))  deallocate(self%ui, stat=istatus) ; istatus0 = istatus0 + istatus
   if (allocated(self%fi))  deallocate(self%fi, stat=istatus) ; istatus0 = istatus0 + istatus
   if (allocated(self%zi))  deallocate(self%zi, stat=istatus) ; istatus0 = istatus0 + istatus

   if (istatus0.ne.0) then
      call env%error('deallocation failed', source)
      return
   end if

end subroutine memfree


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
   if (x.ge.1.0_wp) then
      fsw = 0.0_wp
   elseif (x.le.flow) then
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
   if (x.ge.1.0_wp) then
      dfsw = 0.0_wp
   elseif (x.le.flow) then
      dfsw = 0.0_wp
   else
      dfsw = f30*(1.0_wp-x)*(x-1.0_wp)*(x-1.0_wp+eta)*(x-1.0_wp+eta)/(eta**5)
   end if

end function dfsw


!> Dump an array (ngrid, ncol) or just a column.
subroutine ptcart(self, label, ncol, icol, x)

   type(TddCosmo), intent(in) :: self
   character(len=*), intent(in) :: label
   integer, intent(in) :: ncol, icol
   real(wp), intent(in) :: x(self%ngrid, ncol)

   integer :: ig, noff, nprt, ic, j
   character(len=*), parameter :: f1000 = "(1x, i5, f14.8)"
   character(len=*), parameter :: f1010 = "(6x, 5i14)"
   character(len=*), parameter :: f1020 = "(1x, i5, 5f14.8)"

   ! print header :
   if (ncol.eq.1) then
      write (iout, '(3x, a, 1x, "(column ", i4")")') label, icol
   else
      write (iout, '(3x, a)') label
   end if

   ! print entries :
   if (ncol.eq.1) then
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

   type(TddCosmo), intent(in) :: self
   character (len=*), intent(in) :: label
   integer, intent(in) :: ncol, icol
   real(wp), intent(in) :: x(self%nylm, ncol)

   integer :: l, m, ind, noff, nprt, ic, j
   character(len=*), parameter :: f1000 = "(1x, i3, i4, f14.8)"
   character(len=*), parameter :: f1010 = "(8x, 5i14)"
   character(len=*), parameter :: f1020 = "(1x, i3, i4, 5f14.8)"

   ! print header :
   if (ncol.eq.1) then
      write (iout, '(3x, a, 1x, "(column ", i4")")') label, icol
   else
      write (iout, '(3x, a)') label
   end if

   ! print entries :
   if (ncol.eq.1) then
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
subroutine intrhs(self, isph, x, xlm)

   type(TddCosmo), intent(in) :: self
   integer, intent(in) :: isph
   real(wp), intent(in) :: x(self%ngrid)
   real(wp), intent(inout) :: xlm(self%nylm)

   integer :: ig

   ! initialize
   xlm = 0.0_wp

   ! accumulate over integration points
   do ig = 1, self%ngrid
      xlm = xlm + self%basis(:, ig)*self%w(ig)*x(ig)
   end do

   ! printing
   if (self%iprint.ge.5) then
      call ptcart(self, 'pot', 1, isph, x)
      call prtsph(self, 'vlm', 1, isph, xlm)
   end if

end subroutine intrhs


!> Compute spherical harmonics
pure subroutine ylmbas(self, x, basloc, vplm, vcos, vsin)

   type(TddCosmo), intent(in) :: self
   real(wp), intent(in) :: x(3)
   real(wp), intent(out) :: basloc(self%nylm), vplm(self%nylm)
   real(wp), intent(out) :: vcos(self%lmax+1), vsin(self%lmax+1)

   integer :: l, m, ind
   real(wp) :: cthe, sthe, cphi, sphi, plm

   ! get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
   ! coordinates of x.

   ! evaluate cos(theta) ; sin(theta)
   cthe = x(3)
   sthe = sqrt(1.0_wp - cthe*cthe)

   ! evalutate cos(phi) ; sin(phi)
   if (sthe.ne.0.0_wp) then
      cphi = x(1)/sthe
      sphi = x(2)/sthe
   else
      cphi = 0.0_wp
      sphi = 0.0_wp
   end if

   ! evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is
   ! pointless if z = 1, as the only non vanishing terms will be the
   ! ones with m=0.
   if(sthe.ne.0.0_wp) then
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

   type(TddCosmo), intent(in) :: self
   real(wp), intent(in) :: x(3)
   real(wp), intent(inout) :: basloc(self%nylm), vplm(self%nylm)
   real(wp), intent(inout) :: dbsloc(3, self%nylm)
   real(wp), intent(inout) :: vcos(self%lmax+1), vsin(self%lmax+1)

   integer :: l, m, ind, VC, VS
   real(wp) :: cthe, sthe, cphi, sphi, plm, fln, pp1, pm1, pp
   real(wp) :: et(3), ep(3)

   ! get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
   ! coordinates of x.
   cthe = x(3)
   sthe = sqrt(1.0_wp - cthe*cthe)

   if (sthe.ne.0.0_wp) then
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

   if(sthe.ne.0.0_wp) then
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

   if (sthe.ne.0.0_wp) then
      ! not (NORTH or SOUTH pole)
      call trgev(self, cphi, sphi, vcos, vsin)
   else
      ! NORTH or SOUTH pole
      vcos = 1.0_wp
      vsin = 0.0_wp
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
      if (l.gt.0) then
         dbsloc(:, ind) = fln*vplm(ind+1)*et(:)
      else
         dbsloc(:, ind) = 0.0_wp
      end if
      !dir$ simd
      do m = 1, l
         fln = self%facs(ind+m)
         plm = fln*vplm(ind+m)
         pp1 = 0.0_wp
         if (m.lt.l) pp1 = -0.5_wp*vplm(ind+m+1)
         pm1 = 0.5_wp*(real(l+m, wp)*real(l-m+1, wp)*vplm(ind+m-1))
         pp  = pp1 + pm1

         ! m > 0
         basloc(ind+m) = plm*vcos(m+1)

         if (sthe.ne.0.0_wp) then
            ! not (NORTH or SOUTH pole)
            dbsloc(:, ind+m) = -fln*pp*vcos(m+1)*et(:) - real(m, wp)*plm*vsin(m+1)*ep(:)
         else
            ! NORTH or SOUTH pole
            dbsloc(:, ind+m) = -fln*pp*vcos(m+1)*et(:) - fln*pp*ep(:)*VC
         end if

         ! m < 0
         basloc(ind-m) = plm*vsin(m+1)

         if (sthe.ne.0.0_wp) then
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

   type(TddCosmo), intent(in) :: self
   real(wp), intent(in) :: x, y
   real(wp), intent(inout) :: plm(self%nylm)

   integer :: m, ind, l, ind2
   real(wp) :: fact, pmm, somx2, pmm1, pmmo, pll, fm, fl

   fact  = 1.0_wp
   pmm   = 1.0_wp
   somx2 = y

   do m = 0, self%lmax
      ind      = (m + 1)*(m + 1)
      plm(ind) = pmm
      if(m.eq.self%lmax) return
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

   type(TddCosmo), intent(in) :: self
   real(wp), intent(in) :: x, y
   real(wp), intent(inout) :: cx(max((self%lmax+1), 2)), sx(max((self%lmax+1), 2))

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

   type(TddCosmo), intent(in) :: self
   real(wp), intent(in) :: t
   real(wp), intent(in) :: sigma(self%nylm), basloc(self%nylm)
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

   type(TddCosmo), intent(in) :: self
   real(wp), intent(in)  :: phi(self%ncav)
   real(wp), intent(out) :: g(self%ngrid, self%nsph)

   integer :: isph, ig, ic

   !> Initialize
   ic = 0
   g(:, :) = 0._wp

   ! loop over spheres
   do isph = 1, self%nsph

      ! loop over points
      do ig = 1, self%ngrid

         ! nonzero contribution from point
         if (self%ui(ig, isph).ne.0.0_wp) then

            ! advance cavity point counter
            ic = ic + 1

            ! weight by (negative) characteristic function
            g(ig, isph) = -self%ui(ig, isph) * phi(ic)

         end if

      end do
   end do

end subroutine wghpot


!> Compute H-norm
pure subroutine hsnorm(self, u, unorm)

   type(TddCosmo), intent(in) :: self
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
pure subroutine adjrhs(self, isph, xi, vlm, basloc, vplm, vcos, vsin)

   type(TddCosmo), intent(in) :: self
   integer, intent(in) :: isph
   real(wp), intent(in) :: xi(self%ngrid, self%nsph)
   real(wp), intent(inout) :: vlm(self%nylm)
   real(wp), intent(inout) :: basloc(self%nylm), vplm(self%nylm)
   real(wp), intent(inout) :: vcos(self%lmax+1), vsin(self%lmax+1)

   integer :: ij, jsph, ig, l, ind, m
   real(wp) :: vji(3), vvji, tji, sji(3), xji, oji, fac, ffac, t

   ! loop over neighbors of i-sphere
   do ij = self%inl(isph), self%inl(isph+1)-1

      ! j-sphere is neighbor
      jsph = self%nl(ij)

      ! loop over integration points
      do ig = 1, self%ngrid

         ! compute t_n^ji = | r_j + \rho_j s_n - r_i | / \rho_i
         vji  = self%csph(:, jsph) + self%rsph(jsph)*self%grid(:, ig) - self%csph(:, isph)
         vvji = vji(1)**2 + vji(2)**2 + vji(3)**2
         vvji = sqrt(vvji)
         tji  = vvji/self%rsph(isph)

         ! point is INSIDE i-sphere (+ transition layer)
         if (tji.lt.(1.0_wp + (se+1.0_wp)/2.0_wp*self%eta)) then

            ! compute s_n^ji
            sji = vji/vvji

            ! compute \chi(t_n^ji)
            xji = fsw(tji, se, self%eta)

            ! compute W_n^ji
            if (self%fi(ig, jsph).gt.1.0_wp) then
               oji = xji/self%fi(ig, jsph)
            else
               oji = xji
            end if

            ! compute Y_l^m(s_n^ji)
            call ylmbas(self, sji, basloc, vplm, vcos, vsin)

            ! initialize (t_n^ji)^l
            t = 1.0_wp

            ! compute w_n * xi(n, j) * W_n^ji
            fac = self%w(ig) * xi(ig, jsph) * oji

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

   type(TddCosmo), intent(in) :: self

   1000 format(/, &
      ' An implementation of COSMO using a domain decomposition linear scaling strategy.', /)
   1010 format(' Parameters:', /, &
      '   number of grid points:                  ', 8x, i8, /,   &
      '   number of spheres:                      ', 8x, i8, /,   &
      '   lmax for the spherical harmonics basis: ', 8x, i8, /,   &
      '   convergence threshold:                  ', 8x, d8.1, /, &
      '   regularization parameter (eta):         ', 8x, f8.3, /, &
      '   dielectric constant:                    ', 8x, f8.4/)

   if (self%iprint.gt.0) then

      write(iout, 1000)
      write(iout, 1010) self%ngrid, self%nsph, self%lmax, 10.0_wp**(-self%iconv), self%eta, self%eps

      if (self%igrad.eq.1)  write(iout, 1013)
      1013   format(' Compute forces.'//)

   end if

end subroutine header


!> Compute the first part of <S, L^(x)X>
pure subroutine fdoka(self, isph, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, fx)

   type(TddCosmo), intent(in) :: self
   integer, intent(in) :: isph
   real(wp), intent(in) :: sigma(self%nylm, self%nsph)
   real(wp), intent(in) :: xi(self%ngrid)
   real(wp), intent(inout) :: basloc(self%nylm), vplm(self%nylm)
   real(wp), intent(inout) :: dbsloc(3, self%nylm)
   real(wp), intent(inout) :: vcos(self%lmax+1), vsin(self%lmax+1)
   real(wp), intent(inout) :: fx(3)

   integer :: ig, ij, jsph, l, ind, m
   real(wp) :: vvij, tij, xij, oij, t, fac, fl, f1, f2, f3, beta, tlow, thigh
   real(wp) :: vij(3), sij(3), alp(3), va(3)

   tlow  = 1.0_wp - 0.5_wp*(1.0_wp - se)*self%eta
   thigh = 1.0_wp + 0.5_wp*(1.0_wp + se)*self%eta

   do ig = 1, self%ngrid
      va = 0.0_wp
      do ij = self%inl(isph), self%inl(isph+1) - 1
         jsph = self%nl(ij)
         vij = self%csph(:, isph) + self%rsph(isph)*self%grid(:, ig) - self%csph(:, jsph)
         vvij = vij(1)**2 + vij(2)**2 + vij(3)**2
         vvij = sqrt(vvij)
         tij = vvij/self%rsph(jsph)

         if (tij.ge.thigh) cycle

         sij = vij/vvij
         call dbasis(self, sij, basloc, dbsloc, vplm, vcos, vsin)
         alp = 0.0_wp
         t = 1.0_wp
         do l = 1, self%lmax
            ind = l*l + l + 1
            fl = real(l, wp)
            fac = t/self%facl(ind)
            do m = -l, l
               f2 = fac*sigma(ind+m, jsph)
               f1 = f2*fl*basloc(ind+m)
               alp(:) = alp(:) + f1*sij(:) + f2*dbsloc(:, ind+m)
            end do
            t = t*tij
         end do
         beta = intmlp(self, tij, sigma(:, jsph), basloc)
         xij = fsw(tij, se, self%eta)
         if (self%fi(ig, isph).gt.1.0_wp) then
            oij = xij/self%fi(ig, isph)
            f2  = -oij/self%fi(ig, isph)
         else
            oij = xij
            f2  = 0.0_wp
         end if
         f1 = oij/self%rsph(jsph)
         va(:) = va(:) + f1*alp(:) + beta*f2*self%zi(:, ig, isph)
         if (tij .gt. tlow) then
            f3 = beta*dfsw(tij, se, self%eta)/self%rsph(jsph)
            if (self%fi(ig, isph).gt.1.0_wp) f3 = f3/self%fi(ig, isph)
            va(:) = va(:) + f3*sij(:)
         end if
      end do
      fx = fx - self%w(ig)*xi(ig)*va(:)
   end do

end subroutine fdoka


!> Compute the the second part of <S, L^(x)X>
pure subroutine fdokb(self, isph, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, fx)

   type(TddCosmo), intent(in) :: self
   integer, intent(in) :: isph
   real(wp), intent(in) :: sigma(self%nylm, self%nsph)
   real(wp), intent(in) :: xi(self%ngrid, self%nsph)
   real(wp), intent(inout) :: basloc(self%nylm), vplm(self%nylm)
   real(wp), intent(inout) :: dbsloc(3, self%nylm)
   real(wp), intent(inout) :: vcos(self%lmax+1), vsin(self%lmax+1)
   real(wp), intent(inout) :: fx(3)

   integer :: ig, ji, jsph, l, ind, m, jk, ksph
   logical :: proc
   real(wp) :: vvji, tji, xji, oji, t, fac, fl, f1, f2, beta, di, tlow, thigh
   real(wp) :: b, g1, g2, vvjk, tjk, f, xjk
   real(wp) :: vji(3), sji(3), alp(3), vb(3), vjk(3), sjk(3), vc(3)



   tlow  = 1.0_wp - 0.5_wp*(1.0_wp - se)*self%eta
   thigh = 1.0_wp + 0.5_wp*(1.0_wp + se)*self%eta

   do ig = 1, self%ngrid
      vb = 0.0_wp
      vc = 0.0_wp
      do ji = self%inl(isph), self%inl(isph+1) - 1
         jsph = self%nl(ji)
         vji = self%csph(:, jsph) + self%rsph(jsph)*self%grid(:, ig) - self%csph(:, isph)
         vvji = vji(1)**2 + vji(2)**2 + vji(3)**2
         vvji = sqrt(vvji)
         tji = vvji/self%rsph(isph)

         if (tji.gt.thigh) cycle

         sji = vji/vvji
         call dbasis(self, sji, basloc, dbsloc, vplm, vcos, vsin)

         alp = 0.0_wp
         t   = 1.0_wp
         do l = 1, self%lmax
            ind = l*l + l + 1
            fl = real(l, wp)
            fac = t/self%facl(ind)
            do m = -l, l
               f2 = fac*sigma(ind+m, isph)
               f1 = f2*fl*basloc(ind+m)
               alp = alp + f1*sji + f2*dbsloc(:, ind+m)
            end do
            t = t*tji
         end do
         xji = fsw(tji, se, self%eta)
         if (self%fi(ig, jsph).gt.1.0_wp) then
            oji = xji/self%fi(ig, jsph)
         else
            oji = xji
         end if
         f1 = oji/self%rsph(isph)
         vb = vb + f1*alp*xi(ig, jsph)
         if (tji .gt. tlow) then
            beta = intmlp(self, tji, sigma(:, isph), basloc)
            if (self%fi(ig, jsph) .gt. 1.0_wp) then
               di  = 1.0_wp/self%fi(ig, jsph)
               fac = di*xji
               proc = .false.
               b    = 0.0_wp
               do jk = self%inl(jsph), self%inl(jsph+1) - 1
                  ksph = self%nl(jk)
                  vjk = self%csph(:, jsph) + self%rsph(jsph)*self%grid(:, ig) - self%csph(:, ksph)
                  vvjk = vjk(1)**2 + vjk(2)**2 + vjk(3)**2
                  vvjk = sqrt(vvjk)
                  tjk = vvjk/self%rsph(ksph)
                  if (ksph.ne.isph) then
                     if (tjk .le. thigh) then
                        proc = .true.
                        sjk = vjk/vvjk
                        call ylmbas(self, sjk, basloc, vplm, vcos, vsin)
                        g1 = intmlp(self, tjk, sigma(:, ksph), basloc)
                        xjk = fsw(tjk, se, self%eta)
                        b = b + g1*xjk
                     end if
                  end if
               end do
               if (proc) then
                  g1 = di*di*dfsw(tji, se, self%eta)/self%rsph(isph)
                  g2 = g1*xi(ig, jsph)*b
                  vc = vc + g2*sji
               end if
            else
               di = 1.0_wp
               fac = 0.0_wp
            end if
            f2 = (1.0_wp-fac)*di*dfsw(tji, se, self%eta)/self%rsph(isph)
            vb = vb + f2*xi(ig, jsph)*beta*sji
         end if
      end do
      fx = fx + self%w(ig)*(vb - vc)
   end do

end subroutine fdokb


!> Compute the U^(x)\Phi contribution to <S, g^(x)>
pure subroutine fdoga(self, isph, xi, phi, fx)

   type(TddCosmo), intent(in) :: self
   integer, intent(in) :: isph
   real(wp), intent(in) :: xi(self%ngrid, self%nsph), phi(self%ngrid, self%nsph)
   real(wp), intent(inout) :: fx(3)

   integer :: ig, ji, jsph
   real(wp) :: vvji, tji, fac, swthr
   real(wp) :: alp(3), vji(3), sji(3)

   do ig = 1, self%ngrid
      alp = 0.0_wp
      if (self%ui(ig, isph) .gt. 0.0_wp .and. self%ui(ig, isph).lt.1.0_wp) then
         alp = alp + phi(ig, isph)*xi(ig, isph)*self%zi(:, ig, isph)
      end if
      do ji = self%inl(isph), self%inl(isph+1) - 1
         jsph = self%nl(ji)
         vji = self%csph(:, jsph) + self%rsph(jsph)*self%grid(:, ig) - self%csph(:, isph)
         vvji = vji(1)**2 + vji(2)**2 + vji(3)**2
         vvji = sqrt(vvji)
         tji = vvji/self%rsph(isph)
         swthr = 1.0_wp + (se + 1._wp)*self%eta / 2._wp
         if (tji.lt.swthr .and. tji.gt.swthr-self%eta .and. self%ui(ig, jsph).gt.0.0_wp) then
            sji = vji/vvji
            fac = - dfsw(tji, se, self%eta)/self%rsph(isph)
            alp = alp + fac*phi(ig, jsph)*xi(ig, jsph)*sji
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
pure subroutine calcv(self, first, isph, pot, sigma, basloc, vplm, vcos, vsin)

   type(TddCosmo), intent(in) :: self
   logical, intent(in) :: first
   integer, intent(in) :: isph
   real(wp), intent(in) :: sigma(self%nylm, self%nsph)
   real(wp), intent(inout) :: pot(self%ngrid)
   real(wp), intent(inout) :: basloc(self%nylm)
   real(wp), intent(inout) :: vplm(self%nylm)
   real(wp), intent(inout) :: vcos(self%lmax+1)
   real(wp), intent(inout) :: vsin(self%lmax+1)

   integer :: its, ij, jsph
   real(wp) :: vij(3), sij(3)
   real(wp) :: vvij2, vvij, tij, xij, oij, stslm, stslm2, stslm3

   ! initialize
   pot(:) = 0.0_wp

   ! if 1st iteration of Jacobi method, then done!
   if (first)  return

   ! loop over grid points
   do its = 1, self%ngrid

      ! contribution from integration point present
      if (self%ui(its, isph).lt.1.0_wp) then

         ! loop over neighbors of i-sphere
         do ij = self%inl(isph), self%inl(isph+1)-1

            ! neighbor is j-sphere
            jsph = self%nl(ij)

            ! compute t_n^ij = | r_i + \rho_i s_n - r_j | / \rho_j
            vij  = self%csph(:, isph) + self%rsph(isph)*self%grid(:, its) - self%csph(:, jsph)
            vvij2 = vij(1)**2 + vij(2)**2 + vij(3)**2
            vvij = sqrt(vvij2)
            tij  = vvij / self%rsph(jsph)

            ! point is INSIDE j-sphere
            if (tij.lt.1.0_wp) then

               ! compute s_n^ij = (r_i + \rho_i s_n - r_j) / | ... |
               sij = vij / vvij

               ! compute \chi(t_n^ij)
               xij = fsw(tij, se, self%eta)

               ! compute W_n^ij
               if (self%fi(its, isph).gt.1.0_wp) then

                  oij = xij / self%fi(its, isph)

               else

                  oij = xij

               end if

               ! compute Y_l^m(s_n^ij)
               call ylmbas(self, sij, basloc, vplm, vcos, vsin)

               ! accumulate over j, l, m
               pot(its) = pot(its) + oij * intmlp(self, tij, sigma(:, jsph), basloc)

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

   type(TddCosmo), intent(in) :: self
   real(wp), intent(in) :: s(self%nylm, self%nsph)
   real(wp), intent(inout) :: xi(self%ncav)

   integer :: its, isph, ii

   ii = 0
   do isph = 1, self%nsph
      do its = 1, self%ngrid
         if (self%ui(its, isph) .gt. 0.0_wp) then
            ii = ii + 1
            xi(ii) = self%w(its)*self%ui(its, isph) &
               & * dot_product(self%basis(:, its), s(:, isph))
         end if
      end do
   end do

end subroutine ddmkxi


!> Compute
!
! \zeta(n, i) =
!
!  1/2 f(\eps) sum w_n U_n^i Y_l^m(s_n) [S_i]_l^m
!              l, m
!
pure subroutine ddmkzeta(self, s, zeta)
   type(TddCosmo), intent(in) :: self
   real(wp), intent(in) :: s(self%nylm, self%nsph)
   real(wp), intent(inout) :: zeta(self%ncav)

   integer :: its, isph, ii

   ii = 0
   do isph = 1, self%nsph
      do its = 1, self%ngrid
         if (self%ui(its, isph) .gt. 0.0_wp) then
            ii = ii + 1
            zeta(ii) = self%w(its) * self%ui(its, isph) &
               & * dot_product(self%basis(:, its), s(:, isph))
         end if
      end do
   end do

   zeta = 0.5_wp*((self%eps-1.0_wp)/self%eps)*zeta

end subroutine ddmkzeta


end module xtb_solv_ddcosmo_core
