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
module xtb_xtb_gfn1
   use xtb_mctc_accuracy, only : wp
   use xtb_xtb_data
   use xtb_aoparam
   implicit none
   private

   public :: initGFN1


   interface initGFN1
      module procedure :: initData
      module procedure :: initRepulsion
   end interface initGFN1


   !> Maximum number of elements supported by GFN1-xTB
   integer, parameter :: gfn1Elem = 86


   !>
   real(wp), parameter :: kExp = 1.5_wp

   !>
   real(wp), parameter :: kExpLight = kExp

   !>
   real(wp), parameter :: rExp = 1.0_wp

   !>
   real(wp), parameter :: repAlpha(1:gfn1Elem) = [&
      & 2.209700_wp, 1.382907_wp, 0.671797_wp, 0.865377_wp, 1.093544_wp, &
      & 1.281954_wp, 1.727773_wp, 2.004253_wp, 2.507078_wp, 3.038727_wp, &
      & 0.704472_wp, 0.862629_wp, 0.929219_wp, 0.948165_wp, 1.067197_wp, &
      & 1.200803_wp, 1.404155_wp, 1.323756_wp, 0.581529_wp, 0.665588_wp, &
      & 0.841357_wp, 0.828638_wp, 1.061627_wp, 0.997051_wp, 1.019783_wp, &
      & 1.137174_wp, 1.188538_wp, 1.399197_wp, 1.199230_wp, 1.145056_wp, &
      & 1.047536_wp, 1.129480_wp, 1.233641_wp, 1.270088_wp, 1.153580_wp, &
      & 1.335287_wp, 0.554032_wp, 0.657904_wp, 0.760144_wp, 0.739520_wp, &
      & 0.895357_wp, 0.944064_wp, 1.028240_wp, 1.066144_wp, 1.131380_wp, &
      & 1.206869_wp, 1.058886_wp, 1.026434_wp, 0.898148_wp, 1.008192_wp, &
      & 0.982673_wp, 0.973410_wp, 0.949181_wp, 1.074785_wp, 0.579919_wp, &
      & 0.606485_wp, 1.311200_wp, 0.839861_wp, 0.847281_wp, 0.854701_wp, &
      & 0.862121_wp, 0.869541_wp, 0.876961_wp, 0.884381_wp, 0.891801_wp, &
      & 0.899221_wp, 0.906641_wp, 0.914061_wp, 0.921481_wp, 0.928901_wp, &
      & 0.936321_wp, 0.853744_wp, 0.971873_wp, 0.992643_wp, 1.132106_wp, &
      & 1.118216_wp, 1.245003_wp, 1.304590_wp, 1.293034_wp, 1.181865_wp, &
      & 0.976397_wp, 0.988859_wp, 1.047194_wp, 1.013118_wp, 0.964652_wp, &
      & 0.998641_wp]

   !>
   real(wp), parameter :: repZeff(1:gfn1Elem) = [&
      &  1.116244_wp,  0.440231_wp,  2.747587_wp,  4.076830_wp,  4.458376_wp, &
      &  4.428763_wp,  5.498808_wp,  5.171786_wp,  6.931741_wp,  9.102523_wp, &
      & 10.591259_wp, 15.238107_wp, 16.283595_wp, 16.898359_wp, 15.249559_wp, &
      & 15.100323_wp, 17.000000_wp, 17.153132_wp, 20.831436_wp, 19.840212_wp, &
      & 18.676202_wp, 17.084130_wp, 22.352532_wp, 22.873486_wp, 24.160655_wp, &
      & 25.983149_wp, 27.169215_wp, 23.396999_wp, 29.000000_wp, 31.185765_wp, &
      & 33.128619_wp, 35.493164_wp, 36.125762_wp, 32.148852_wp, 35.000000_wp, &
      & 36.000000_wp, 39.653032_wp, 38.924904_wp, 39.000000_wp, 36.521516_wp, &
      & 40.803132_wp, 41.939347_wp, 43.000000_wp, 44.492732_wp, 45.241537_wp, &
      & 42.105527_wp, 43.201446_wp, 49.016827_wp, 51.718417_wp, 54.503455_wp, &
      & 50.757213_wp, 49.215262_wp, 53.000000_wp, 52.500985_wp, 65.029838_wp, &
      & 46.532974_wp, 48.337542_wp, 30.638143_wp, 34.130718_wp, 37.623294_wp, &
      & 41.115870_wp, 44.608445_wp, 48.101021_wp, 51.593596_wp, 55.086172_wp, &
      & 58.578748_wp, 62.071323_wp, 65.563899_wp, 69.056474_wp, 72.549050_wp, &
      & 76.041625_wp, 55.222897_wp, 63.743065_wp, 74.000000_wp, 75.000000_wp, &
      & 76.000000_wp, 77.000000_wp, 78.000000_wp, 79.000000_wp, 80.000000_wp, &
      & 81.000000_wp, 79.578302_wp, 83.000000_wp, 84.000000_wp, 85.000000_wp, &
      & 86.000000_wp]


contains


subroutine initData(self)

   !>
   type(TxTBData), intent(out) :: self

   call initGFN1(self%repulsion)

end subroutine initData


subroutine initRepulsion(self)

   !>
   type(TRepulsionData), intent(out) :: self

   self%cutoff = 40.0_wp
   self%kExp = kExp
   self%kExpLight = kExpLight
   self%rExp = rExp
   self%alpha = rep(1, :gfn1Elem) ! repAlpha
   self%zeff = rep(2, :gfn1Elem) ! repZeff

end subroutine initRepulsion


end module xtb_xtb_gfn1
