! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> TODO
module xtb_xtb_ipea
   use xtb_mctc_accuracy, only : wp
   use xtb_param_atomicrad, only : atomicRad
   use xtb_param_paulingen, only : paulingEN
   use xtb_type_param, only : TxTBParameter
   use xtb_xtb_data
   use xtb_xtb_gfn1, only : initGFN1, setGFN1PairParam, setGFN1ReferenceOcc, &
      & setGFN1NumberOfPrimitives
   use xtb_aoparam
   implicit none
   private

   public :: initIPEA, ipeaGlobals


   interface initIPEA
      module procedure :: initData
      module procedure :: initRepulsion
      module procedure :: initCoulomb
      module procedure :: initHamiltonian
   end interface initIPEA


   type(TxTBParameter), parameter :: ipeaGlobals = TxTBParameter( &
      ks =      2.150000000000000_wp, &
      kp =      2.300000000000000_wp, &
      kd =      2.000000000000000_wp, &
      kf =      .000000000000000_wp, &
      kdiffa =  2.250000000000000_wp, &
      kdiffb =  2.250000000000000_wp, &
      ipeashift = 1.780690000000000_wp, &
      zcnf =    .000000000000000_wp, &
      tscal =   .000000000000000_wp, &
      kcn =     1.000000000000000_wp, &
      fpol =    -.3000000000000000_wp, &
      ken =     -.5000000000000000_wp, &
      lshift =  .000000000000000_wp, &
      lshifta = .000000000000000_wp, &
      split =   .000000000000000_wp, &
      zqf =     .000000000000000_wp, &
      alphaj =  1.000000000000000_wp, &
      kexpo =   1.500000000000000_wp, &
      dispa =   .6300000000000000_wp, &
      dispb =   5.000000000000000_wp, &
      dispc =   2.400000000000000_wp, &
      dispatm = .7000000000000001_wp, &
      xbdamp =  .4400000000000001_wp, &
      xbrad =   1.300000000000000_wp )


   !> Maximum number of elements supported by IPEA-xTB
   integer, parameter :: maxElem = 86


   ! ========================================================================
   ! REPULSION DATA
   !>
   real(wp), parameter :: kExp = 1.5_wp

   !>
   real(wp), parameter :: kExpLight = kExp

   !>
   real(wp), parameter :: rExp = 1.0_wp


   ! ========================================================================
   ! COULOMB DATA
   !>

   !>
   logical, parameter :: thirdOrderShellResolved = .false.


   ! ========================================================================
   ! HAMILTONIAN DATA
   !>



contains


subroutine initData(self)

   !>
   type(TxTBData), intent(out) :: self

   self%nShell = ao_n(:maxElem)

   call initIPEA(self%repulsion)
   call initIPEA(self%coulomb, self%nShell)
   call initIPEA(self%hamiltonian, self%nShell)
   allocate(self%halogen)
   call initGFN1(self%halogen)

end subroutine initData


subroutine initRepulsion(self)

   !>
   type(TRepulsionData), intent(out) :: self

   !call init(self, kExp, kExpLight, rExp, 0.0_wp, repAlpha, repZeff)

end subroutine initRepulsion


subroutine initCoulomb(self, nShell)

   !>
   type(TCoulombData), intent(out) :: self

   !>
   integer, intent(in) :: nShell(:)

   !self%chemicalHardness = gam(:maxElem) ! chemcialHardness
   !self%thirdOrderAtom = gam3(:maxElem) ! thirdOrderAtom
   !self%shellHardness = lpar(0:2, :maxElem)

end subroutine initCoulomb


subroutine initHamiltonian(self, nShell)

   !>
   type(THamiltonianData), intent(out) :: self

   !>
   integer, intent(in) :: nShell(:)

   integer :: mShell, nPrim, lAng
   integer :: iZp, iSh
   logical :: valShell(0:3)

   mShell = maxval(nShell)
   self%angShell = ao_l(:mShell, :maxElem)

   self%electronegativity = en(:maxElem)
   self%atomicRad = atomicRad(:maxElem)
   self%shellPoly = polyr(:, :maxElem)
   self%pairParam = kpair(:maxElem, :maxElem)
   self%selfEnergy = ao_lev(:mShell, :maxElem)
   self%slaterExponent = ao_exp(:mShell, :maxElem)
   self%principalQuantumNumber = ao_pqn(:mShell, :maxElem)

   allocate(self%valenceShell(mShell, maxElem))
   call generateValenceShellData(self%valenceShell, nShell, self%angShell)

   allocate(self%referenceOcc(mShell, maxElem))
   call setGFN1ReferenceOcc(self, nShell)

   allocate(self%numberOfPrimitives(mShell, maxElem))
   call setGFN1NumberOfPrimitives(self, nShell)

end subroutine initHamiltonian


end module xtb_xtb_ipea
