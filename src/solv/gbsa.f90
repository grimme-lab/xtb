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

!> Implementation for generalized Born and related solvation models
module xtb_solv_gbsa
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas, only : mctc_dot, mctc_gemv, mctc_symv
   use xtb_mctc_constants, only : fourpi, pi
   use xtb_mctc_convert, only : aatoau
   use xtb_mctc_math, only : matDet3x3
   use xtb_mctc_search, only : bisectSearch
   use xtb_solv_born, only : compute_bornr
   use xtb_solv_kernel, only : gbKernel, addBornMatSaltStill, addBornMatStill, &
      & addBornMatP16, addGradientSaltStill, addGradientStill, addGradientP16, &
      & addBornDerivSaltStill, addBornDerivStill
   use xtb_solv_lebedev, only : gridSize, getAngGrid
   use xtb_solv_sasa, only : compute_numsa
   use xtb_type_environment, only : TEnvironment
   use xtb_type_solvation, only : TSolvation
   implicit none
   private

   public :: TBorn, init
   public :: addGradientHBond
   public :: addHBondDeriv
   public :: getADet, addADetDeriv
   public :: compute_fhb, getDebyeHueckel, update_nnlist_gbsa


   type, extends(TSolvation) :: TBorn

      !> number of atoms
      integer :: nat

      !> atom types
      integer, allocatable :: at(:)

      !> number of pairs
      integer :: ntpair

      !> number of angular grid points
      integer :: nang 

      !> angular grid
      real(wp), allocatable :: angGrid(:, :)
      real(wp), allocatable :: angWeight(:)

      !> van der Waals radii of the particles
      real(wp), allocatable :: vdwr(:)

      !> pair descreening approximation radii
      real(wp), allocatable :: rho(:)

      !> offset van der Waals radii
      real(wp), allocatable :: svdw(:)

      !> cut-off radius for the Born radius NN list
      real(wp) :: lrcut

      !> cut-off radius for the SASA NN list
      real(wp) :: srcut

      !> number of neighbors for Born radii
      integer :: nnrad

      !> number of neighbors for SASA computation
      integer, allocatable :: nnsas(:)

      !> neighbors of an atom for Born radii
      integer, allocatable :: nnlistr(:, :)

      !> neighbors of an atom for SASA
      integer, allocatable :: nnlists(:, :)

      !> all pairs indeces array
      integer, allocatable :: ppind(:, :)

      !> all pairs vector differences and magnitudes array
      real(wp), allocatable :: ddpair(:, :)

      !> Atom specific surface data
      real(wp), allocatable :: vdwsa(:)
      real(wp), allocatable :: wrp(:)
      real(wp), allocatable :: trj2(:, :)

      !> Dielectric data
      real(wp) :: gborn

      !> Born radii
      real(wp), allocatable :: brad(:)

      !> Salt screening
      real(wp), allocatable :: ionscr(:)
      real(wp), allocatable :: discr(:)

      !> Atomic surfaces
      real(wp) :: gsasa
      real(wp) :: sasagam
      real(wp), allocatable :: gamsasa(:)
      real(wp), allocatable :: sasa(:)

      !> Molecular Surface gradient
      real(wp), allocatable :: dsdr(:, :)
      real(wp), allocatable :: dsdrt(:, :, :)

      !> Hydrogen bond contribution
      real(wp), allocatable :: hbmag(:)

      !> Hydrogen bond contribution
      real(wp), allocatable :: hbw(:)

      !> Hydrogen bond gradient
      real(wp), allocatable :: dhbdw(:)

      !> Born radii gradient
      real(wp), allocatable :: brdr(:, :, :)

      !> Shape descriptor
      real(wp) :: aDet

      !> Use salt screening
      logical :: lsalt

      !> Use hydrogen bond correction
      logical :: lhb

      !> Interaction kernel
      integer :: kernel

      !> Scaling factor for Born radii
      real(wp) :: bornScale

      !> Analytical linearized Poisson-Boltzmann constant
      real(wp) :: alpbet

      !> Dielectric constant
      real(wp) :: dielectricConst

      !> Debye screening length
      real(wp) :: kappa

      !> Ion radii
      real(wp) :: ionRad

      !> Dielectric screening
      real(wp) :: keps

      !> Free energy shift
      real(wp) :: gshift

      !> Scratch for potential
      real(wp), allocatable :: shift(:)

      !> Interaction matrix
      real(wp), allocatable :: bornMat(:, :)

   contains

      !> Update coordinates and internal state
      procedure :: update

      !> Get complete interaction matrix
      procedure :: addBornMatrix

      !> Add potential shift
      procedure :: addShift

      !> Calculate solvation energy
      procedure :: getEnergy

      !> Obtain the respective energy contributions
      procedure :: getEnergyParts

      !> Calculate derivatives of solvation energy
      procedure :: addGradient

      !> Get complete interaction matrix
      procedure :: addBornDeriv

   end type TBorn


   !> Initialize data straucture
   interface init
      module procedure :: initBorn
   end interface init

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
subroutine initBorn(self, env, num, vdwRad, dielectricConst, freeEnergyShift, &
      & descreening, bornScale, bornOffset, surfaceTension, probeRad, rCutoff, &
      & rOffset, nAng, hBondStrength, temperature, kernel, alpb, ionSt, ionRad)

   !> Error source
   character(len=*), parameter :: source = 'solv_gbsa_initBorn'

   !> Instance of the solvation model
   type(TBorn), intent(out) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Atomic numbers
   integer, intent(in) :: num(:)

   !> Van-der-Waals Radii
   real(wp), intent(in) :: vdwRad(:)

   !> Dielectric constant of the solvent
   real(wp), intent(in) :: dielectricConst

   !> State specific free energy offset
   real(wp), intent(in) :: freeEnergyShift

   !> Dielectric descreening parameter
   real(wp), intent(in) :: descreening(:)

   !> Scaling factor for Born radii
   real(wp), intent(in) :: bornScale

   !> Offset parameter for Born radii integration
   real(wp), intent(in) :: bornOffset

   !> Surface tension scaling
   real(wp), intent(in) :: surfaceTension(:)

   !> Probe radius of the solvent
   real(wp), intent(in) :: probeRad

   !> Real-space cutoff for Born radii integration
   real(wp), intent(in) :: rCutoff

   !> Offset for surface integration cutoff
   real(wp), intent(in) :: rOffset

   !> Number of angular grid points for integration
   integer, intent(in) :: nAng

   !> Hydrogen bond magnitude
   real(wp), intent(in) :: hBondStrength(:)

   !> Temperature
   real(wp), intent(in) :: temperature

   !> Generalized Born interaction kernel
   integer, intent(in) :: kernel

   !> Use analytical linearized Poisson-Boltzmann model instead of GB
   logical, intent(in) :: alpb

   !> Ion strength
   real(wp), intent(in), optional :: ionSt

   !> Ion radius
   real(wp), intent(in), optional :: ionRad

   integer :: iat, izp, jat, ij, iang, ierr
   real(wp) :: r

   self%kernel = kernel
   self%gshift = freeEnergyShift
   self%dielectricConst = dielectricConst
   if (alpb) then
      self%alpbet = 0.571412_wp / self%dielectricConst
   else
      self%alpbet = 0.0_wp
   end if
   self%keps = (1.0_wp/self%dielectricConst - 1.0_wp) / (1.0_wp + self%alpbet)

   self%nat = size(num)
   self%at = num
   self%ntpair = self%nat*(self%nat-1)/2
   allocate(self%ppind(2, self%ntpair))
   allocate(self%nnsas(self%nat))
   allocate(self%nnlistr(3, self%ntpair))
   allocate(self%nnlists(self%nat, self%nat))
   allocate(self%ddpair(4, self%ntpair))
   ij = 0
   do iat = 1, self%nat
      do jat = 1, iat-1
         ij = ij+1
         self%ppind(1, ij) = iat
         self%ppind(2, ij) = jat
      enddo
   enddo

   if (present(ionSt)) then
      self%lsalt = ionSt > 0.0_wp
   else
      self%lsalt = .false.
   end if
   if (self%lsalt) then
      allocate(self%ionscr(self%nat))
      allocate(self%discr(self%nat))
      if (present(ionRad)) then
         self%ionRad = ionRad
      else
         self%ionRad = 0.0_wp
      end if
      self%kappa = sqrt(self%dielectricConst*temperature*kappaConst/ionSt)*aatoau
      self%kappa = 1.0_wp/self%kappa
   end if

   self%lrcut = rCutoff
   self%bornScale = bornScale
   allocate(self%vdwr(self%nat))
   allocate(self%rho(self%nat))
   allocate(self%svdw(self%nat))
   allocate(self%brad(self%nat))
   allocate(self%brdr(3, self%nat, self%nat))

   do iat = 1, self%nat
      izp = num(iat)
      self%vdwr(iat) = vdwRad(izp)
      self%rho(iat) = self%vdwr(iat) * descreening(izp)
      self%svdw(iat) = self%vdwr(iat) - bornOffset
   end do

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

   call bisectSearch(iang, gridSize, nAng)
   allocate(self%angGrid(3, gridSize(iang)))
   allocate(self%angWeight(gridSize(iang)))
   call getAngGrid(iang, self%angGrid, self%angWeight, ierr)

   allocate(self%hbmag(self%nat))
   do iat = 1, self%nat
      izp = num(iat)
      self%hbmag(iat) = hBondStrength(izp)
   end do
   self%lhb = any(self%hbmag < 0.0_wp)
   if (self%lhb) then
      allocate(self%hbw(self%nat))
      allocate(self%dhbdw(self%nat))
   end if

   allocate(self%shift(self%nat))
   allocate(self%bornMat(self%nat, self%nat))

end subroutine initBorn


subroutine update(self, env, num, xyz)

   !> Instance of the solvation model
   class(TBorn), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Atomic numbers
   integer, intent(in) :: num(:)

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   ! initialize the neighbor list
   call update_nnlist_gbsa(self%nat, self%ntpair, self%ppind, xyz, &
      & self%lrcut, self%srcut, self%nnsas, self%nnlists, self%nnrad, &
      & self%nnlistr, self%ddpair, .false.)

   call compute_bornr(self%nat, self%nnrad, self%nnlistr, self%ddpair, &
      & self%vdwr, self%rho, self%svdw, self%bornScale, self%brad, self%brdr)

   ! compute solvent accessible surface and its derivatives
   call compute_numsa(self%nat, self%nnsas, self%nnlists, xyz, self%vdwsa, &
      & self%wrp, self%trj2, self%angWeight, self%angGrid, &
      & self%sasa, self%dsdrt)

   ! contract surface gradient
   call mctc_gemv(self%dsdrt, self%gamsasa, self%dsdr)
   self%gsasa = mctc_dot(self%sasa, self%gamsasa)

   ! compute the Debye-Hueckel ion exclusion term
   if (self%lsalt) then
      call getDebyeHueckel(self%nat, self%dielectricConst, self%kappa, &
         & self%ionRad, self%brad, self%ionscr, self%discr)
   end if

   ! compute the HB term
   if (self%lhb) then
      call compute_fhb(self%nat, self%hbmag, self%vdwsa, self%sasa, &
         & self%hbw, self%dhbdw)
   end if

   if (self%alpbet > 0.0_wp) then
      call getADet(self%nat, xyz, self%vdwr, self%aDet)
   end if

   self%bornMat(:, :) = 0.0_wp
   call self%addBornMatrix(env, num, xyz, self%bornMat)

end subroutine update


subroutine addBornMatrix(self, env, num, xyz, bornMat)

   !> Instance of the solvation model
   class(TBorn), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Atomic numbers
   integer, intent(in) :: num(:)

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Born interaction matrix
   real(wp), intent(inout) :: bornMat(:, :)

   integer :: iat

   select case(self%kernel)
   case(gbKernel%still)
      if (self%lsalt) then
         call addBornMatSaltStill(self%nat, self%ntpair, self%ppind, self%ddpair, &
            & self%kappa, self%brad, self%ionscr, bornMat)
      else
         call addBornMatStill(self%nat, self%ntpair, self%ppind, self%ddpair, &
            & self%keps, self%brad, bornMat)
      end if
   case(gbKernel%p16)
      call addBornMatP16(self%nat, self%ntpair, self%ppind, self%ddpair, &
         & self%keps, self%brad, bornMat)
   end select

   ! compute the HB term
   if (self%lhb) then
      do iat = 1, self%nat
         bornMat(iat, iat) = bornMat(iat, iat) + 2*self%hbw(iat)
      enddo
   endif

   ! ALPB shape dependent correction for charged systems
   if (self%alpbet > 0.0_wp) then
      bornMat(:self%nat, :self%nat) = bornMat(:self%nat, :self%nat) &
         & + self%keps * self%alpbet / self%aDet
   end if

end subroutine addBornMatrix


!> Add potential shift
subroutine addShift(self, env, qat, qsh, atomicShift, shellShift)

   !> Error source
   character(len=*), parameter :: source = 'solv_gbsa_addShift'

   !> Instance of the solvation model
   class(TBorn), intent(inout) :: self

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

   call mctc_symv(self%bornMat, qat, self%shift)
   atomicShift(:) = atomicShift + self%shift

end subroutine addShift


!> Calculate solvation energy
subroutine getEnergy(self, env, qat, qsh, energy)

   !> Error source
   character(len=*), parameter :: source = 'solv_gbsa_getEnergy'

   !> Instance of the solvation model
   class(TBorn), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)

   !> Total solvation energy
   real(wp), intent(out) :: energy

   call mctc_symv(self%bornMat, qat, self%shift, alpha=0.5_wp)
   energy = mctc_dot(qat, self%shift) + self%gsasa + self%gshift

end subroutine getEnergy


!> Calculate solvation energy
subroutine getEnergyParts(self, env, qat, qsh, gborn, ghb, gsasa, gshift)

   !> Error source
   character(len=*), parameter :: source = 'solv_gbsa_getEnergy'

   !> Instance of the solvation model
   class(TBorn), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)

   !> Born solvation energy
   real(wp), intent(out) :: gborn

   !> Hydrogen bonding energy
   real(wp), intent(out) :: ghb

   !> Non-polar surface energy
   real(wp), intent(out) :: gsasa

   !> State specific shift
   real(wp), intent(out) :: gshift

   integer :: iat

   ghb = 0.0_wp
   if (self%lhb) then
      do iat = 1, self%nat
         ghb = ghb + self%hbw(iat) * qat(iat)**2
      end do
   end if

   call mctc_symv(self%bornMat, qat, self%shift, alpha=0.5_wp)
   gborn = mctc_dot(qat, self%shift) - ghb

   gsasa = self%gsasa
   gshift = self%gshift

end subroutine getEnergyParts


!> Calculate derivatives of solvation energy
subroutine addGradient(self, env, num, xyz, qat, qsh, gradient)

   !> Error source
   character(len=*), parameter :: source = 'solv_gbsa_addGradient'

   !> Instance of the solvation model
   class(TBorn), intent(inout) :: self

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

   real(wp) :: gborn, ghb

   select case(self%kernel)
   case(gbKernel%still)
      if (self%lsalt) then
         call addGradientSaltStill(self%nat, self%ntpair, self%ppind, self%ddpair, &
            & qat, self%kappa, self%brad, self%brdr, self%ionscr, self%discr, &
            & gborn, gradient)
      else
         call addGradientStill(self%nat, self%ntpair, self%ppind, self%ddpair, &
            & qat, self%keps, self%brad, self%brdr, gborn, gradient)
      end if
   case(gbKernel%p16)
      call addGradientP16(self%nat, self%ntpair, self%ppind, self%ddpair, &
         & qat, self%keps, self%brad, self%brdr, gborn, gradient)
   end select

   gradient = gradient + self%dsdr

   if (self%lhb) then
      call addGradientHBond(self%nat, self%at, qat, self%hbw, self%dhbdw, &
         & self%dsdrt, ghb, gradient)
   else
      ghb = 0.0_wp
   endif

   if (self%alpbet > 0.0_wp) then
      gborn = gborn + sum(qat)**2 * self%alpbet / self%aDet * self%kEps
      call addADetDeriv(self%nat, xyz, self%vdwr, self%kEps*self%alpbet, &
         & qat, gradient)
   end if

end subroutine addGradient


subroutine getADet(nAtom, xyz, rad, aDet)

   !> Number of atoms
   integer, intent(in) :: nAtom

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Atomic radii
   real(wp), intent(in) :: rad(:)

   !> Shape descriptor of the structure
   real(wp), intent(out) :: aDet

   integer :: iat
   real(wp) :: r2, rad2, rad3, totRad3, vec(3), center(3), inertia(3, 3)
   real(wp), parameter :: tof = 2.0_wp/5.0_wp, unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   totRad3 = 0.0_wp
   center(:) = 0.0_wp
   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      totRad3 = totRad3 + rad3
      center(:) = center + xyz(:, iat) * rad3
   end do
   center = center / totRad3

   inertia(:, :) = 0.0_wp
   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vec(:) = xyz(:, iat) - center
      r2 = sum(vec**2)
      inertia(:, :) = inertia + rad3 * ((r2 + tof*rad2) * unity &
         & - spread(vec, 1, 3) * spread(vec, 2, 3))
   end do

   aDet = sqrt(matDet3x3(inertia)**(1.0_wp/3.0_wp)/(tof*totRad3))

end subroutine getADet


subroutine addADetDeriv(nAtom, xyz, rad, kEps, qvec, gradient)

   !> Number of atoms
   integer, intent(in) :: nAtom

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Atomic radii
   real(wp), intent(in) :: rad(:)

   real(wp), intent(in) :: kEps
   real(wp), intent(in) :: qvec(:)

   !> Molecular gradient
   real(wp), intent(inout) :: gradient(:, :)

   integer :: iat
   real(wp) :: r2, rad2, rad3, totRad3, vec(3), center(3), inertia(3, 3), aDet
   real(wp) :: aDeriv(3, 3), qtotal
   real(wp), parameter :: tof = 2.0_wp/5.0_wp, unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   qtotal = 0.0_wp
   totRad3 = 0.0_wp
   center(:) = 0.0_wp
   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      totRad3 = totRad3 + rad3
      center(:) = center + xyz(:, iat) * rad3
      qtotal = qtotal + qvec(iat)
   end do
   center = center / totRad3

   inertia(:, :) = 0.0_wp
   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vec(:) = xyz(:, iat) - center
      r2 = sum(vec**2)
      inertia(:, :) = inertia + rad3 * ((r2 + tof*rad2) * unity &
         & - spread(vec, 1, 3) * spread(vec, 2, 3))
   end do
   aDet = sqrt(matDet3x3(inertia)**(1.0_wp/3.0_wp)/(tof*totRad3))

   aDeriv(:, :) = reshape([&
      & inertia(1,1)*(inertia(2,2)+inertia(3,3))-inertia(1,2)**2-inertia(1,3)**2, &
      & inertia(1,2)*inertia(3,3)-inertia(1,3)*inertia(2,3), & ! xy
      & inertia(1,3)*inertia(2,2)-inertia(1,2)*inertia(3,2), & ! xz
      & inertia(1,2)*inertia(3,3)-inertia(1,3)*inertia(2,3), & ! xy
      & inertia(2,2)*(inertia(1,1)+inertia(3,3))-inertia(1,2)**2-inertia(2,3)**2, &
      & inertia(1,1)*inertia(2,3)-inertia(1,2)*inertia(1,3), & ! yz
      & inertia(1,3)*inertia(2,2)-inertia(1,2)*inertia(3,2), & ! xz
      & inertia(1,1)*inertia(2,3)-inertia(1,2)*inertia(1,3), & ! yz
      & inertia(3,3)*(inertia(1,1)+inertia(2,2))-inertia(1,3)**2-inertia(2,3)**2],&
      & shape=[3, 3]) * (250.0_wp / (48.0_wp * totRad3**3 * aDet**5)) &
      & * (-0.5_wp * kEps * qtotal**2 / aDet**2)

   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vec(:) = xyz(:, iat) - center
      gradient(:, iat) = gradient(:, iat) + rad3 * matmul(aderiv, vec)
   end do

end subroutine addADetDeriv


!> compute the Debye-Hueckel ion exclusion term
pure subroutine getDebyeHueckel(nat, dielectricConst, kappa, ionRad, brad, &
      & ionscr, discr)

   integer, intent(in) :: nat
   real(wp), intent(in) :: kappa
   real(wp), intent(in) :: dielectricConst
   real(wp), intent(in) :: brad(:)
   real(wp), intent(in) :: ionRad
   real(wp), intent(out) :: ionscr(:)
   real(wp), intent(out) :: discr(:)

   integer :: i
   real(wp) :: aa, gg

   aa=0.5_wp/dielectricConst
   do i = 1, nat
      gg=kappa*(brad(i)+ionRad)
      ionscr(i)=aa*exp(gg)/(1.0_wp+gg)
      discr(i)=ionscr(i)*kappa*gg/(1.0_wp+gg)
   enddo

end subroutine getDebyeHueckel


!> setup a pairlist and compute pair distances of all neighbors
!  within thresholds lrcut and srcut.
subroutine update_nnlist_gbsa(nat, ntpair, ppind, xyz, lrcut, srcut, &
      & nnsas, nnlists, nnrad, nnlistr, ddpair, parallel)

   integer, intent(in) :: nat
   integer, intent(in) :: ntpair
   integer, intent(in) :: ppind(:, :)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: lrcut
   real(wp), intent(in) :: srcut
   integer, intent(out) :: nnsas(:)
   integer, intent(out) :: nnlists(:, :)
   integer, intent(out) :: nnrad
   integer, intent(out) :: nnlistr(:, :)
   real(wp), intent(out) :: ddpair(:, :)
   logical, intent(in) :: parallel

   if (parallel) then
      call update_nnlist_gbsa_parallel(nat, ntpair, ppind, xyz, &
         & lrcut, srcut, nnsas, nnlists, nnrad, nnlistr, ddpair)
   else
      call update_nnlist_gbsa_sequential(nat, ntpair, ppind, xyz, &
         & lrcut, srcut, nnsas, nnlists, nnrad, nnlistr, ddpair)
   endif

end subroutine update_nnlist_gbsa
!> setup a pairlist and compute pair distances of all neighbors
!  within thresholds lrcut and srcut.
!  OMP parallel version.
subroutine update_nnlist_gbsa_parallel(nat, ntpair, ppind, xyz, lrcut, &
      & srcut, nnsas, nnlists, nnrad, nnlistr, ddpair)
!$ use omp_lib

   integer, intent(in) :: nat
   integer, intent(in) :: ntpair
   integer, intent(in) :: ppind(:, :)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: lrcut
   real(wp), intent(in) :: srcut
   integer, intent(out) :: nnsas(:)
   integer, intent(out) :: nnlists(:, :)
   integer, intent(out) :: nnrad
   integer, intent(out) :: nnlistr(:, :)
   real(wp), intent(out) :: ddpair(:, :)

   integer kk, i1, i2
   real(wp) rcutn2, lrcut2, srcut2
   real(wp) x, y, z, dr2
   integer ip, ip2, thrid, nproc
   integer, allocatable :: npid(:)
   integer, allocatable :: plisttr(:, :, :)
   integer, allocatable :: nntmp(:)
   integer, allocatable :: nnls(:, :)

   nproc=1
!$ nproc=omp_get_max_threads()

   allocate(plisttr(3, ntpair, nproc), nnls(nat, nat))
   allocate(nntmp(nat), npid(nproc))
   npid = 0

   lrcut2 = lrcut*lrcut
   srcut2 = srcut*srcut

   nnsas=0
   nnlists=0
!$omp parallel default(none) &
!$omp&         shared ( xyz,lrcut2,srcut2,ntpair,ppind,nat,nnlists,nnsas,ddpair ) &
!$omp&         private( i1,i2,x,y,z,dr2,ip,ip2,thrid,nntmp,nnls ) &
!$omp&         shared ( plisttr, npid )
   ip=0
   ip2=0
   nntmp=0
   nnls=0
   thrid=1
!$ thrid=omp_get_thread_num() + 1
!$omp do
   do kk=1,ntpair
      i1=ppind(1,kk)
      i2=ppind(2,kk)
      x=xyz(1,i1)-xyz(1,i2)
      y=xyz(2,i1)-xyz(2,i2)
      z=xyz(3,i1)-xyz(3,i2)
      dr2=x**2+y**2+z**2
      ddpair(2,kk)=x
      ddpair(3,kk)=y
      ddpair(4,kk)=z
      ddpair(1,kk)=sqrt(dr2)
      if(dr2.lt.lrcut2) then
         ip = ip + 1
         plisttr(1,ip,thrid)=i1
         plisttr(2,ip,thrid)=i2
         plisttr(3,ip,thrid)=kk
         if(dr2.lt.srcut2) then
            nntmp(i1) = nntmp(i1) + 1
            nntmp(i2) = nntmp(i2) + 1
            nnls(nntmp(i1),i1)=i2
            nnls(nntmp(i2),i2)=i1
         endif
      endif
   enddo
!$omp end do
   npid(thrid)=ip
!$omp critical
   do i1=1,nat
      do i2=1,nntmp(i1)
         nnlists(nnsas(i1)+i2,i1)=nnls(i2,i1)
      enddo
      nnsas(i1)=nnsas(i1)+nntmp(i1)
   enddo
!$omp end critical
!$omp end parallel

   nnrad=0
   do thrid=1,nproc
      do kk = nnrad+1,nnrad+npid(thrid)
         nnlistr(1:3,kk)=plisttr(1:3,kk-nnrad,thrid)
      enddo
      nnrad = nnrad + npid(thrid)
   enddo

   deallocate(nntmp,nnls)

end subroutine update_nnlist_gbsa_parallel
!> setup a pairlist and compute pair distances of all neighbors
!  within thresholds lrcut and srcut
!  Sequential version.
pure subroutine update_nnlist_gbsa_sequential(nat, ntpair, ppind, xyz, lrcut, &
      & srcut, nnsas, nnlists, nnrad, nnlistr, ddpair)

   integer, intent(in) :: nat
   integer, intent(in) :: ntpair
   integer, intent(in) :: ppind(:, :)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: lrcut
   real(wp), intent(in) :: srcut
   integer, intent(out) :: nnsas(:)
   integer, intent(out) :: nnlists(:, :)
   integer, intent(out) :: nnrad
   integer, intent(out) :: nnlistr(:, :)
   real(wp), intent(out) :: ddpair(:, :)

   integer kk, i1, i2
   real(wp) rcutn2, lrcut2, srcut2
   real(wp) x, y, z, dr2
   integer ip, ip2

   lrcut2 = lrcut*lrcut
   srcut2 = srcut*srcut

   nnsas=0
   nnlists=0
   ip=0
   ip2=0
   nnlistr=0
   do kk=1,ntpair
      i1=ppind(1,kk)
      i2=ppind(2,kk)
      x=xyz(1,i1)-xyz(1,i2)
      y=xyz(2,i1)-xyz(2,i2)
      z=xyz(3,i1)-xyz(3,i2)
      dr2=x**2+y**2+z**2
      ddpair(2,kk)=x
      ddpair(3,kk)=y
      ddpair(4,kk)=z
      ddpair(1,kk)=sqrt(dr2)
      if(dr2.lt.lrcut2) then
         ip = ip + 1
         nnlistr(1,ip)=i1
         nnlistr(2,ip)=i2
         nnlistr(3,ip)=kk
         if(dr2.lt.srcut2) then
            nnsas(i1) = nnsas(i1) + 1
            nnsas(i2) = nnsas(i2) + 1
            nnlists(nnsas(i1),i1)=i2
            nnlists(nnsas(i2),i2)=i1
         endif
      endif
   enddo
   nnrad = ip

end subroutine update_nnlist_gbsa_sequential


!> Compute contributions to potential for hydrogen bonding correction
pure subroutine compute_fhb(nat, hbmag, vdwsa, sasa, hbw, dhbdw)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Hydrogen bonding strength
   real(wp), intent(in) :: hbmag(:)

   !> Van-der-Waals radius
   real(wp), intent(in) :: vdwsa(:)

   !> Surface area including surface tension
   real(wp), intent(in) :: sasa(:)

   !> Actual hydrogen bonding contribution
   real(wp), intent(out) :: hbw(:)

   !> Actual hydrogen bonding contribution
   real(wp), intent(out) :: dhbdw(:)

   integer  :: i
   integer  :: iz, nhb
   real(wp) :: hbed, dhbed
   real(wp) :: smaxd, sasad, sasaw
   real(wp) :: sfw, dsfw, w3, w2, w1
   integer  :: j
   real(wp) :: wbh, wah

   hbw(:) = 0.0_wp
   dhbdw(:) = 0.0_wp

   do i = 1, nat
      ! SASA-D for HB
      smaxd = 1.0_wp/(vdwsa(i)*vdwsa(i))
      sasad = sasa(i)*smaxd
      hbw(i) = hbmag(i)*sasad
      dhbdw(i) = hbmag(i)*smaxd
   enddo

end subroutine compute_fhb


!> Compute contributions to energy and gradient for hydrogen bonding correction
pure subroutine addGradientHBond(nat, at, q, hbw, dhbdw, dsdrt, ghb, gradient)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Atomic numbers
   integer, intent(in) :: at(:)

   !> Partial charges
   real(wp), intent(in) :: q(:)

   !> Actual hydrogen bonding contribution
   real(wp), intent(in) :: hbw(:)

   !> Actual hydrogen bonding contribution
   real(wp), intent(in) :: dhbdw(:)

   !> Derivative of the surface area w.r.t. to cartesian coordinates
   real(wp), intent(in) :: dsdrt(:, :, :)

   !> Hydrogen bonding energy
   real(wp), intent(out)   :: ghb

   !> Gradient of hydrogen bonding energy
   real(wp), intent(inout) :: gradient(:, :)

   integer  :: i, j
   real(wp) :: dhbed
   real(wp) :: qq

   ghb=0.0_wp
   do i = 1, nat
      qq = q(i)*q(i)
      ghb = ghb + hbw(i)*qq
   enddo

   do i = 1, nat
      dhbed = dhbdw(i)
      if(abs(dhbed).le.0.0_wp) cycle
      dhbed=dhbed*q(i)*q(i)
      do j = 1, nat
         gradient(:, j) = gradient(:, j) + dsdrt(:, j, i)*dhbed
      enddo
   enddo

end subroutine addGradientHBond


pure subroutine addBornDeriv(self,q,gborn,ghb,dAmatdr,Afac)
   implicit none
   class(TBorn), intent(in) :: self

   real(wp), intent(in)    :: q(self%nat)
   real(wp), intent(inout) :: dAmatdr(3,self%nat,self%nat)
   real(wp), intent(inout) :: Afac(3,self%nat)
   real(wp), intent(out)   :: gborn
   real(wp), intent(out)   :: ghb

   integer :: i,j,k,nnj
   integer :: kk
   real(wp), parameter :: a13=1._wp/3._wp
   real(wp), parameter :: a4=0.25_wp
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp) :: aa,r2,fgb,fgb2,br3
   real(wp) :: qq,dd,expd,dfgb,dfgb2,dfgb3,ap,bp,qfg
   real(wp) :: gg,expa,aii,egb
   real(wp) :: r0vdw,r01,r02,ar02
   real(wp) :: grddbi,grddbj
   real(wp) :: dr(3),r

   egb = 0._wp

   select case(self%kernel)
   case(gbKernel%still)
      if (self%lsalt) then
         call addBornDerivSaltStill(self%nat, self%ntpair, self%ppind, &
            & self%ddpair, q, self%kappa, self%brad, self%brdr, self%ionscr, &
            & self%discr, gborn, dAmatdr, Afac)
      else
         call addBornDerivStill(self%nat, self%ntpair, self%ppind, self%ddpair, &
            & q, self%keps, self%brad, self%brdr, gborn, dAmatdr, Afac)
      endif
   case(gbKernel%p16)
   end select

   if (self%lhb) then
      call addHBondDeriv(self%nat, q, self%hbw, self%dhbdw, self%dsdrt, &
         & ghb, dAmatdr)
   endif

end subroutine addBornDeriv


pure subroutine addHBondDeriv(nat,q,hbw,dhbdw,dsdrt,ghb,dAmatdr)

   integer, intent(in) :: nat
   real(wp), intent(in) :: q(:)
   real(wp), intent(in) :: hbw(:)
   real(wp), intent(in) :: dhbdw(:)
   real(wp), intent(in) :: dsdrt(:, :, :)
   real(wp), intent(out) :: ghb
   real(wp), intent(inout) :: dAmatdr(:, :, :)

   integer  :: i,j
   real(wp) :: dhbed
   real(wp) :: qq

   ghb=0.0_wp
   do i = 1, nat
      qq = q(i)*q(i)
      ghb = ghb + hbw(i)*qq
   enddo

   do i = 1, nat
      dhbed=dhbdw(i)
      if(abs(dhbed).le.0.0_wp) cycle
      dhbed=dhbed*q(i)
      dAmatdr(:,:,i) = dAmatdr(:,:,i) + dsdrt(:,:,i)*dhbed
   enddo

end subroutine addHBondDeriv


end module xtb_solv_gbsa
