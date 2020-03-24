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

!> TODO: GFN0-xTB parametrisation data
module xtb_xtb_gfn0
   use xtb_mctc_accuracy, only : wp
   use xtb_param_atomicrad, only : atomicRad
   use xtb_xtb_data
   use xtb_xtb_gfn2, only : setGFN2ReferenceOcc
   use xtb_aoparam
   implicit none
   private

   public :: initGFN0
   public :: setGFN0NumberOfPrimitives


   interface initGFN0
      module procedure :: initData
      module procedure :: initRepulsion
      module procedure :: initCoulomb
      module procedure :: initHamiltonian
   end interface initGFN0


   !> Maximum number of elements supported by GFN0-xTB
   integer, parameter :: maxElem = 86

   ! ========================================================================
   ! REPULSION DATA
   !>
   real(wp), parameter :: kExp = 1.5_wp

   !>
   real(wp), parameter :: kExpLight = kExp

   !>
   real(wp), parameter :: rExp = 1.0_wp

   !>
   real(wp), parameter :: repAlpha(1:maxElem) = [&
      & 2.1885472_wp, 2.2714498_wp, 0.6634645_wp, 0.9267640_wp, 1.1164621_wp, &
      & 1.2680750_wp, 1.6211038_wp, 2.1037547_wp, 2.2062651_wp, 1.9166982_wp, &
      & 0.8129781_wp, 0.8408742_wp, 0.8361156_wp, 0.8859465_wp, 1.0684151_wp, &
      & 1.1882871_wp, 1.4429448_wp, 1.1993811_wp, 0.5700050_wp, 0.8345430_wp, &
      & 0.6840185_wp, 0.7915733_wp, 1.0676223_wp, 0.9216746_wp, 1.1151815_wp, &
      & 1.1883881_wp, 1.1895339_wp, 1.2692713_wp, 1.1734165_wp, 1.0018764_wp, &
      & 1.1597304_wp, 1.1708353_wp, 1.2085038_wp, 1.1161800_wp, 1.3193094_wp, &
      & 0.7670615_wp, 0.6171015_wp, 0.8421909_wp, 0.6513468_wp, 0.6906528_wp, &
      & 0.8705783_wp, 0.9711021_wp, 1.0252504_wp, 0.9847071_wp, 1.0559061_wp, &
      & 1.0645317_wp, 0.9139636_wp, 0.9095541_wp, 0.9965441_wp, 1.0676257_wp, &
      & 1.0759855_wp, 0.8659486_wp, 0.9301733_wp, 0.8139884_wp, 0.5842740_wp, &
      & 0.8070627_wp, 0.6961124_wp, 0.7599095_wp, 0.7667071_wp, 0.7735047_wp, &
      & 0.7803023_wp, 0.7870999_wp, 0.7938975_wp, 0.8006951_wp, 0.8074927_wp, &
      & 0.8142903_wp, 0.8210879_wp, 0.8278855_wp, 0.8346831_wp, 0.8414808_wp, &
      & 0.8482784_wp, 0.8803684_wp, 0.9915500_wp, 0.9875716_wp, 1.1535600_wp, &
      & 1.1418384_wp, 1.1434832_wp, 1.1783705_wp, 1.0591477_wp, 0.9794378_wp, &
      & 1.2439938_wp, 1.0437958_wp, 1.1391049_wp, 0.9115474_wp, 0.9157573_wp, &
      & 0.8137168_wp]

   !>
   real(wp), parameter :: repZeff(1:maxElem) = [&
      & 1.2455414_wp,  1.3440060_wp,  1.1710492_wp,  2.9064151_wp,  4.4020866_wp, &
      & 4.3101011_wp,  4.5460146_wp,  4.7850603_wp,  7.3393960_wp,  4.2503997_wp, &
      &10.5220970_wp,  7.7916659_wp, 11.3886282_wp, 13.9495563_wp, 16.7912135_wp, &
      &13.3874290_wp, 13.9700526_wp, 14.4971987_wp, 13.8061512_wp, 13.9719788_wp, &
      &10.9127447_wp, 13.4067871_wp, 16.7322903_wp, 21.8192969_wp, 22.8754319_wp, &
      &25.2196212_wp, 26.9753662_wp, 27.2652026_wp, 26.2195102_wp, 14.3840374_wp, &
      &25.4102208_wp, 43.7565690_wp, 34.9344472_wp, 22.8724870_wp, 34.2378269_wp, &
      &15.1027639_wp, 39.1086736_wp, 32.7340796_wp, 18.6398784_wp, 22.6163764_wp, &
      &27.6545601_wp, 37.8625561_wp, 40.9844265_wp, 30.0686254_wp, 35.5737255_wp, &
      &28.4443233_wp, 25.9740558_wp, 28.8257081_wp, 53.9657064_wp, 88.0203443_wp, &
      &82.7978295_wp, 39.3120212_wp, 49.7072042_wp, 45.1199137_wp, 55.2536842_wp, &
      &50.0381164_wp, 48.0939804_wp, 46.1827790_wp, 46.0844595_wp, 45.9861400_wp, &
      &45.8878205_wp, 45.7895010_wp, 45.6911815_wp, 45.5928620_wp, 45.4945424_wp, &
      &45.3962229_wp, 45.2979034_wp, 45.1995839_wp, 45.1012644_wp, 45.0029449_wp, &
      &44.9046254_wp, 41.1538255_wp, 46.6524574_wp, 53.4995959_wp, 73.8197012_wp, &
      &59.6567627_wp, 50.0720023_wp, 49.4064531_wp, 44.5201114_wp, 39.7677937_wp, &
      &58.8051943_wp,103.0123579_wp, 85.5566053_wp, 70.6036525_wp, 82.8260761_wp, &
      &68.9676875_wp]

   ! ========================================================================
   ! HAMILTONIAN DATA
   !>
   integer, parameter :: valenceShell(3, 1:maxElem) = reshape([&
      & 1, 0, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0, &
      & 1, 1, 0,  1, 1, 0,  1, 1, 1,  1, 1, 0,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 0,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 0,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 0,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 0,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 0,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0, &
      & 1, 1, 1,  1, 1, 1], shape(valenceShell))

   !>
   real(wp), parameter :: referenceOcc(3, 1:maxElem) = reshape([&
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp, &
      & 2.0_wp, 3.0_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp, &
      & 2.0_wp, 6.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 2.0_wp,  2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp, &
      & 2.0_wp, 0.0_wp, 5.0_wp,  2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp, &
      & 2.0_wp, 0.0_wp, 8.0_wp,  2.0_wp, 0.0_wp, 9.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 2.0_wp,  2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp, &
      & 2.0_wp, 0.0_wp, 5.0_wp,  2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp, &
      & 2.0_wp, 0.0_wp, 8.0_wp,  2.0_wp, 0.0_wp, 9.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 2.0_wp, &
      & 2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp,  2.0_wp, 0.0_wp, 5.0_wp, &
      & 2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp,  2.0_wp, 0.0_wp, 8.0_wp, &
      & 2.0_wp, 0.0_wp, 9.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp, &
      & 2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp, &
      & 2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp], shape(referenceOcc))

   !>
   real(wp), parameter :: electronegativity(1:maxElem) = [&
      & 1.92_wp,  3.00_wp,  0.98_wp,  1.57_wp,  2.04_wp, &
      & 2.48_wp,  2.97_wp,  3.44_wp,  3.50_wp,  3.50_wp, &
      & 0.93_wp,  1.31_wp,  1.61_wp,  1.90_wp,  2.19_wp, &
      & 2.58_wp,  3.16_wp,  3.50_wp,  1.45_wp,  1.80_wp, &
      & 1.73_wp,  1.54_wp,  1.63_wp,  1.66_wp,  1.55_wp, &
      & 1.83_wp,  1.88_wp,  1.91_wp,  1.90_wp,  1.65_wp, &
      & 1.81_wp,  2.01_wp,  2.18_wp,  2.55_wp,  2.96_wp, &
      & 3.00_wp,  1.50_wp,  1.50_wp,  1.55_wp,  1.33_wp, &
      & 1.60_wp,  2.16_wp,  1.90_wp,  2.20_wp,  2.28_wp, &
      & 2.20_wp,  1.93_wp,  1.69_wp,  1.78_wp,  1.96_wp, &
      & 2.05_wp,  2.10_wp,  2.66_wp,  2.60_wp,  1.50_wp, &
      & 1.60_wp,  1.50_wp,  1.50_wp,  1.50_wp,  1.50_wp, &
      & 1.50_wp,  1.50_wp,  1.50_wp,  1.50_wp,  1.50_wp, &
      & 1.50_wp,  1.50_wp,  1.50_wp,  1.50_wp,  1.50_wp, &
      & 1.50_wp,  1.30_wp,  1.50_wp,  2.36_wp,  1.90_wp, &
      & 2.20_wp,  2.20_wp,  2.28_wp,  2.54_wp,  2.00_wp, &
      & 1.62_wp,  2.33_wp,  2.02_wp,  2.00_wp,  2.20_wp, &
      & 2.20_wp]
 


contains


subroutine initData(self)

   !> Data instance
   type(TxTBData), intent(out) :: self

   self%nShell = ao_n(:maxElem)

   call initGFN0(self%repulsion)
   call initGFN0(self%coulomb, self%nShell)
   call initGFN0(self%hamiltonian, self%nShell)

end subroutine initData


subroutine initRepulsion(self)

   !> Data instance
   type(TRepulsionData), intent(out) :: self

   call init(self, kExp, kExpLight, rExp, 0.0_wp, repAlpha, repZeff, &
      & electronegativity)
!   self%cutoff = 40.0_wp
!   self%kExp = kExp
!   self%kExpLight = kExpLight
!   self%rExp = rExp
!   self%electronegativity = en(:maxElem)
!   self%alpha = repAlpha(:maxElem)
!   self%zeff = repZeff(:maxElem)

end subroutine initRepulsion


subroutine initCoulomb(self, nShell)

   !> Data instance
   type(TCoulombData), intent(out) :: self

   !>
   integer, intent(in) :: nShell(:)

   self%electronegativity = eeqEN(1:maxElem)
   self%chemicalHardness = gam(1:maxElem)
   self%kCN = eeqkCN(1:maxElem)
   self%chargeWidth = alp0(1:maxElem)

end subroutine initCoulomb


subroutine initHamiltonian(self, nShell)

   !> Data instance
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
   self%kCN = kcnat(:, :maxElem)
   self%selfEnergy = ao_lev(:mShell, :maxElem)
   self%slaterExponent = ao_exp(:mShell, :maxElem)
   self%principalQuantumNumber = ao_pqn(:mShell, :maxElem)
   self%kQShell = kqat(:, :maxElem)
   self%kQAtom = kqat2(:maxElem)

   allocate(self%valenceShell(mShell, maxElem))
   call generateValenceShellData(self%valenceShell, nShell, self%angShell)

   allocate(self%referenceOcc(mShell, maxElem))
   call setGFN2ReferenceOcc(self, nShell)

   allocate(self%numberOfPrimitives(mShell, maxElem))
   call setGFN0NumberOfPrimitives(self, nShell)

end subroutine initHamiltonian


!>
subroutine setGFN0NumberOfPrimitives(self, nShell)

   !> Data instance
   type(THamiltonianData), intent(inout) :: self

   !>
   integer, intent(in) :: nShell(:)

   integer :: nPrim, iZp, iSh

   do iZp = 1, maxElem
      do iSh = 1, nShell(iZp)
         nPrim = 0
         if (iZp <= 2) then
            select case(self%angShell(iSh, iZp))
            case(0)
               if (self%valenceShell(iSh, iZp) /= 0) then
                  nPrim = 3
               else
                  nPrim = 2
               end if
            case(1)
               nPrim = 3
            end select
         else
            select case(self%angShell(iSh, iZp))
            case(0)
               if (self%principalQuantumNumber(iSh, iZp) > 5) then
                  nPrim = 6
               else
                  nPrim = 4
               end if
            case(1)
               if (self%principalQuantumNumber(iSh, iZp) > 5) then
                  nPrim = 6
               else
                  nPrim = 3
               end if
            case(2, 3)
               nPrim = 4
            end select
         end if
         self%numberOfPrimitives(iSh, iZp) = nPrim
      end do
   end do

end subroutine setGFN0NumberOfPrimitives


end module xtb_xtb_gfn0
