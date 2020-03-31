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

!> DFT-D3 Covalent radii for coordination number
module xtb_param_covalentradd3
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_mctc_symbols, only : toNumber
   implicit none
   private

   public :: getCovalentRadD3, covalentRadD3


   !> Get atomic radius for a species
   interface getCovalentRadD3
      module procedure :: getCovalentRadD3Symbol
      module procedure :: getCovalentRadD3Number
   end interface getCovalentRadD3


   !> Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
   !> 188-197), values for metals decreased by 10%.
   real(wp), parameter :: covalentRadD3(1:118) = [ &
      & 0.32_wp,0.46_wp, & ! H,He
      & 1.20_wp,0.94_wp,0.77_wp,0.75_wp,0.71_wp,0.63_wp,0.64_wp,0.67_wp, & ! Li-Ne
      & 1.40_wp,1.25_wp,1.13_wp,1.04_wp,1.10_wp,1.02_wp,0.99_wp,0.96_wp, & ! Na-Ar
      & 1.76_wp,1.54_wp, & ! K,Ca
      &                 1.33_wp,1.22_wp,1.21_wp,1.10_wp,1.07_wp, & ! Sc-
      &                 1.04_wp,1.00_wp,0.99_wp,1.01_wp,1.09_wp, & ! -Zn
      &                 1.12_wp,1.09_wp,1.15_wp,1.10_wp,1.14_wp,1.17_wp, & ! Ga-Kr
      & 1.89_wp,1.67_wp, & ! Rb,Sr
      &                 1.47_wp,1.39_wp,1.32_wp,1.24_wp,1.15_wp, & ! Y-
      &                 1.13_wp,1.13_wp,1.08_wp,1.15_wp,1.23_wp, & ! -Cd
      &                 1.28_wp,1.26_wp,1.26_wp,1.23_wp,1.32_wp,1.31_wp, & ! In-Xe
      & 2.09_wp,1.76_wp, & ! Cs,Ba
      &         1.62_wp,1.47_wp,1.58_wp,1.57_wp,1.56_wp,1.55_wp,1.51_wp, & ! La-Eu
      &         1.52_wp,1.51_wp,1.50_wp,1.49_wp,1.49_wp,1.48_wp,1.53_wp, & ! Gd-Yb
      &                 1.46_wp,1.37_wp,1.31_wp,1.23_wp,1.18_wp, & ! Lu-
      &                 1.16_wp,1.11_wp,1.12_wp,1.13_wp,1.32_wp, & ! -Hg
      &                 1.30_wp,1.30_wp,1.36_wp,1.31_wp,1.38_wp,1.42_wp, & ! Tl-Rn
      & 2.01_wp,1.81_wp, & ! Fr,Ra
      &         1.67_wp,1.58_wp,1.52_wp,1.53_wp,1.54_wp,1.55_wp,1.49_wp, & ! Ac-Am
      &         1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
      &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
      &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
      &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp] & ! Nh-Og
      & * aatoau * 4.0_wp / 3.0_wp


contains


!> Get covalent radius for species with a given symbol
elemental function getCovalentRadD3Symbol(symbol) result(rad)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> atomic radius
   real(wp) :: rad

   rad = getCovalentRadD3(toNumber(symbol))

end function getCovalentRadD3Symbol


!> Get covalent radius for species with a given atomic number
elemental function getCovalentRadD3Number(number) result(rad)

   !> Atomic number
   integer, intent(in) :: number

   !> atomic radius
   real(wp) :: rad

   if (number > 0 .and. number <= size(covalentRadD3, dim=1)) then
      rad = covalentRadD3(number)
   else
      rad = -1.0_wp
   end if

end function getCovalentRadD3Number


end module xtb_param_covalentradd3
