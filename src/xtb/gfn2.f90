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
module xtb_xtb_gfn2
   use xtb_mctc_accuracy, only : wp
   use xtb_xtb_data
   use xtb_aoparam
   implicit none
   private

   public :: initGFN2


   interface initGFN2
      module procedure :: initData
      module procedure :: initRepulsion
   end interface initGFN2

   !> Maximum number of elements supported by GFN2-xTB
   integer, parameter :: gfn2Elem = 86

   !>
   real(wp), parameter :: kExp = 1.5_wp

   !>
   real(wp), parameter :: kExpLight = 1.0_wp

   !>
   real(wp), parameter :: rExp = 1.0_wp

   !>
   real(wp), parameter :: repAlpha(1:gfn2Elem) = [&
      & 2.213717_wp, 3.604670_wp, 0.475307_wp, 0.939696_wp, 1.373856_wp, &
      & 1.247655_wp, 1.682689_wp, 2.165712_wp, 2.421394_wp, 3.318479_wp, &
      & 0.572728_wp, 0.917975_wp, 0.876623_wp, 1.187323_wp, 1.143343_wp, &
      & 1.214553_wp, 1.577144_wp, 0.896198_wp, 0.482206_wp, 0.683051_wp, &
      & 0.574299_wp, 0.723104_wp, 0.928532_wp, 0.966993_wp, 1.071100_wp, &
      & 1.113422_wp, 1.241717_wp, 1.077516_wp, 0.998768_wp, 1.160262_wp, &
      & 1.122923_wp, 1.222349_wp, 1.249372_wp, 1.230284_wp, 1.296174_wp, &
      & 0.908074_wp, 0.574054_wp, 0.697345_wp, 0.706172_wp, 0.681106_wp, &
      & 0.865552_wp, 1.034519_wp, 1.019565_wp, 1.031669_wp, 1.094599_wp, &
      & 1.092745_wp, 0.678344_wp, 0.936236_wp, 1.024007_wp, 1.139959_wp, &
      & 1.122937_wp, 1.000712_wp, 1.017946_wp, 1.012036_wp, 0.585257_wp, &
      & 0.716259_wp, 0.737643_wp, 0.729950_wp, 0.734624_wp, 0.739299_wp, &
      & 0.743973_wp, 0.748648_wp, 0.753322_wp, 0.757996_wp, 0.762671_wp, &
      & 0.767345_wp, 0.772020_wp, 0.776694_wp, 0.781368_wp, 0.786043_wp, &
      & 0.790717_wp, 0.852852_wp, 0.990234_wp, 1.018805_wp, 1.170412_wp, &
      & 1.221937_wp, 1.197148_wp, 1.204081_wp, 0.919210_wp, 1.137360_wp, &
      & 1.399312_wp, 1.179922_wp, 1.130860_wp, 0.957939_wp, 0.963878_wp, &
      & 0.965577_wp]

   !>
   real(wp), parameter :: repZeff(1:gfn2Elem) = [&
      &  1.105388_wp,  1.094283_wp,  1.289367_wp,  4.221216_wp,  7.192431_wp, &
      &  4.231078_wp,  5.242592_wp,  5.784415_wp,  7.021486_wp, 11.041068_wp, &
      &  5.244917_wp, 18.083164_wp, 17.867328_wp, 40.001111_wp, 19.683502_wp, &
      & 14.995090_wp, 17.353134_wp,  7.266606_wp, 10.439482_wp, 14.786701_wp, &
      &  8.004267_wp, 12.036336_wp, 15.677873_wp, 19.517914_wp, 18.760605_wp, &
      & 20.360089_wp, 27.127744_wp, 10.533269_wp,  9.913846_wp, 22.099503_wp, &
      & 31.146750_wp, 42.100144_wp, 39.147587_wp, 27.426779_wp, 32.845361_wp, &
      & 17.363803_wp, 44.338211_wp, 34.365525_wp, 17.326237_wp, 24.263093_wp, &
      & 30.562732_wp, 48.312796_wp, 44.779882_wp, 28.070247_wp, 38.035941_wp, &
      & 28.674700_wp,  6.493286_wp, 26.226628_wp, 63.854240_wp, 80.053438_wp, &
      & 77.057560_wp, 48.614745_wp, 63.319176_wp, 51.188398_wp, 67.249039_wp, &
      & 46.984607_wp, 50.927529_wp, 48.676714_wp, 47.669448_wp, 46.662183_wp, &
      & 45.654917_wp, 44.647651_wp, 43.640385_wp, 42.633120_wp, 41.625854_wp, &
      & 40.618588_wp, 39.611322_wp, 38.604057_wp, 37.596791_wp, 36.589525_wp, &
      & 35.582259_wp, 40.186772_wp, 54.666156_wp, 55.899801_wp, 80.410086_wp, &
      & 62.809871_wp, 56.045639_wp, 53.881425_wp, 14.711475_wp, 51.577544_wp, &
      & 58.801614_wp,102.368258_wp,132.896832_wp, 52.301232_wp, 81.771063_wp, &
      &128.133580_wp]


contains


subroutine initData(self)

   !>
   type(TxTBData), intent(out) :: self

   call initGFN2(self%repulsion)

end subroutine initData


subroutine initRepulsion(self)

   !>
   type(TRepulsionData), intent(out) :: self

   self%cutoff = 40.0_wp
   self%kExp = kExp
   self%kExpLight = kExpLight
   self%rExp = rExp
   self%alpha = rep(1, :gfn2Elem) ! repAlpha
   self%zeff = rep(2, :gfn2Elem) ! repZeff

end subroutine initRepulsion


end module xtb_xtb_gfn2
