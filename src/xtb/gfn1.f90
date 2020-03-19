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
   use xtb_param_atomicrad, only : atomicRad
   use xtb_param_paulingen, only : paulingEN
   use xtb_xtb_data
   use xtb_aoparam
   implicit none
   private

   public :: initGFN1


   interface initGFN1
      module procedure :: initData
      module procedure :: initRepulsion
      module procedure :: initCoulomb
      module procedure :: initHamiltonian
      module procedure :: initHalogen
   end interface initGFN1


   !> Maximum number of elements supported by GFN1-xTB
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
   real(wp), parameter :: repZeff(1:maxElem) = [&
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

   ! ========================================================================
   ! COULOMB DATA
   !>
   real(wp), parameter :: chemicalHardness(1:maxElem) = [&
      & 0.470099_wp, 1.441379_wp, 0.205342_wp, 0.274022_wp, 0.340530_wp, &
      & 0.479988_wp, 0.476106_wp, 0.583349_wp, 0.788194_wp, 0.612878_wp, &
      & 0.165908_wp, 0.354151_wp, 0.221658_wp, 0.438331_wp, 0.798319_wp, &
      & 0.643959_wp, 0.519712_wp, 0.529906_wp, 0.114358_wp, 0.134187_wp, &
      & 0.778545_wp, 1.044998_wp, 0.985157_wp, 0.468100_wp, 0.609868_wp, &
      & 0.900000_wp, 0.426680_wp, 0.367019_wp, 0.260192_wp, 0.209459_wp, &
      & 0.193302_wp, 0.800000_wp, 0.732367_wp, 0.714534_wp, 0.732530_wp, &
      & 0.820312_wp, 0.075735_wp, 0.122861_wp, 0.351290_wp, 0.168219_wp, &
      & 0.175875_wp, 0.384677_wp, 0.405474_wp, 0.305394_wp, 0.293973_wp, &
      & 0.280766_wp, 0.472978_wp, 0.130828_wp, 0.132120_wp, 0.480655_wp, &
      & 0.564406_wp, 0.400301_wp, 0.520472_wp, 0.935394_wp, 0.085110_wp, &
      & 0.137819_wp, 0.495969_wp, 0.350000_wp, 0.342306_wp, 0.334612_wp, &
      & 0.326917_wp, 0.319223_wp, 0.311529_wp, 0.303835_wp, 0.296140_wp, &
      & 0.288446_wp, 0.280752_wp, 0.273058_wp, 0.265364_wp, 0.257669_wp, &
      & 0.249975_wp, 0.269977_wp, 0.239696_wp, 0.243663_wp, 0.362512_wp, &
      & 0.354318_wp, 0.290898_wp, 0.370447_wp, 0.496380_wp, 0.334997_wp, &
      & 0.671316_wp, 1.000000_wp, 0.944879_wp, 1.091248_wp, 1.264162_wp, &
      & 0.798170_wp]

   !>
   logical, parameter :: thirdOrderShellResolved = .false.

   !>
   real(wp), parameter :: thirdOrderAtom(1:maxElem) = [&
      & 0.000000_wp, 1.500000_wp, 1.027370_wp, 0.900554_wp, 1.300000_wp, &
      & 1.053856_wp, 0.042507_wp,-0.005102_wp, 1.615037_wp, 1.600000_wp, &
      & 1.200000_wp, 1.100000_wp, 1.200000_wp, 1.500000_wp, 1.500000_wp, &
      & 1.500000_wp, 1.000000_wp, 0.829312_wp, 0.732923_wp, 1.116963_wp, &
      & 1.000000_wp, 0.739203_wp, 0.800000_wp, 0.800000_wp, 0.300000_wp, &
      & 0.500000_wp, 0.300000_wp, 1.000000_wp, 0.237602_wp, 1.400000_wp, &
      & 1.400000_wp, 1.400000_wp, 1.300000_wp, 1.300000_wp,-0.500000_wp, &
      & 1.000000_wp, 1.500000_wp, 1.300000_wp, 1.400000_wp, 0.581478_wp, &
      & 0.280147_wp, 0.041052_wp, 0.500000_wp, 0.001205_wp, 0.622690_wp, &
      & 0.500000_wp,-0.445675_wp, 1.362587_wp, 1.063557_wp,-0.321283_wp, &
      &-0.341503_wp, 0.894388_wp,-0.500000_wp,-0.800000_wp, 1.500000_wp, &
      & 1.500000_wp, 1.500000_wp, 1.200000_wp, 1.200000_wp, 1.200000_wp, &
      & 1.200000_wp, 1.200000_wp, 1.200000_wp, 1.200000_wp, 1.200000_wp, &
      & 1.200000_wp, 1.200000_wp, 1.200000_wp, 1.200000_wp, 1.200000_wp, &
      & 1.200000_wp, 0.847011_wp, 0.064592_wp,-0.014599_wp, 0.300000_wp, &
      &-0.170295_wp, 0.965726_wp, 1.092759_wp, 0.123512_wp,-0.267745_wp, &
      & 0.936157_wp, 1.500000_wp, 0.877488_wp,-0.035874_wp,-0.860502_wp, &
      &-0.838429_wp] * 0.1_wp

   ! ========================================================================
   ! HAMILTONIAN DATA
   !>
   integer, parameter :: nShell(1:maxElem) = [&
      & 2, 1, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, &
      & 2, 2, 2, 3, 3, 3]

   !>
   integer, parameter :: angularMomentum(3, 1:maxElem) = reshape([&
      & 0, 0,-1,  0,-1,-1,  0, 1,-1,  0, 1,-1,  0, 1,-1,  0, 1,-1,  0, 1,-1, &
      & 0, 1,-1,  0, 1,-1,  0, 1, 2,  0, 1,-1,  0, 1,-1,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1,-1,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1,-1,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1,-1,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1,-1,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1,-1,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1,-1,  0, 1,-1,  0, 1,-1,  0, 1,-1,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2], shape(angularMomentum))

   !>
   integer, parameter :: principalQuantumNumber(3, 1:maxElem) = reshape([&
      & 1, 2,-1,  1,-1,-1,  2, 2,-1,  2, 2,-1,  2, 2,-1,  2, 2,-1,  2, 2,-1, &
      & 2, 2,-1,  2, 2,-1,  2, 2, 3,  3, 3,-1,  3, 3,-1,  3, 3, 3,  3, 3, 3, &
      & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  4, 4,-1,  4, 4, 3,  4, 4, 3, &
      & 4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3, &
      & 4, 4, 3,  4, 4,-1,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, &
      & 4, 4, 4,  5, 5,-1,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4, &
      & 5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5,-1,  5, 5, 5, &
      & 5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5,-1,  5, 5, 4, &
      & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, &
      & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, &
      & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, &
      & 6, 6, 5,  6, 6, 5,  6, 6,-1,  6, 6,-1,  6, 6,-1,  6, 6,-1,  6, 6, 5, &
      & 6, 6, 5,  6, 6, 5], shape(principalQuantumNumber))

   !>
   integer, parameter :: numberOfPrimitives(1, 1:maxElem) = reshape([&
      & 4, 3, 0,  4, 0, 0,  6, 6, 0,  6, 6, 0,  6, 6, 0,  6, 6, 0,  6, 6, 0, &
      & 6, 6, 0,  6, 6, 0,  6, 6, 4,  6, 6, 0,  6, 6, 0,  6, 6, 4,  6, 6, 4, &
      & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 0,  6, 6, 4,  6, 6, 4, &
      & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, &
      & 6, 6, 4,  6, 6, 0,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, &
      & 6, 6, 4,  6, 6, 0,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, &
      & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 0,  6, 6, 4, &
      & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 0,  6, 6, 4, &
      & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, &
      & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, &
      & 6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4,  6, 6, 4, &
      & 6, 6, 4,  6, 6, 4,  6, 6, 0,  6, 6, 0,  6, 6, 0,  6, 6, 0,  6, 6, 4, &
      & 6, 6, 4,  6, 6, 4], shape(numberOfPrimitives))

   !>
   integer, parameter :: valenceShell(3, 1:maxElem) = reshape([&
      & 1, 0, 0,  1, 0, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0, &
      & 1, 1, 0,  1, 1, 0,  1, 1, 1,  1, 1, 0,  1, 1, 0,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 0,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 0,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 0,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 0,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 0,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 1, &
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



contains


subroutine initData(self)

   !>
   type(TxTBData), intent(out) :: self

   self%nShell = ao_n(:maxElem)

   call initGFN1(self%repulsion)
   call initGFN1(self%coulomb, self%nShell)
   call initGFN1(self%hamiltonian, self%nShell)
   allocate(self%halogen)
   call initGFN1(self%halogen)

end subroutine initData


subroutine initRepulsion(self)

   !>
   type(TRepulsionData), intent(out) :: self

   self%cutoff = 40.0_wp
   self%kExp = kExp
   self%kExpLight = kExpLight
   self%rExp = rExp
   self%alpha = rep(1, :maxElem) ! repAlpha
   self%zeff = rep(2, :maxElem) ! repZeff

end subroutine initRepulsion


subroutine initCoulomb(self, nShell)

   !>
   type(TCoulombData), intent(out) :: self

   !>
   integer, intent(in) :: nShell(:)

   self%chemicalHardness = gam(:maxElem) ! chemcialHardness
   self%thirdOrderAtom = gam3(:maxElem) ! thirdOrderAtom
   self%shellHardness = lpar(0:2, :maxElem)

end subroutine initCoulomb


subroutine initHamiltonian(self, nShell)

   !>
   type(THamiltonianData), intent(out) :: self

   !>
   integer, intent(in) :: nShell(:)

   integer :: mShell

   mShell = maxval(nShell)
   self%angShell = ao_l(:mShell, :maxElem)

   self%electronegativity = paulingEN(:maxElem)
   self%atomicRad = atomicRad(:maxElem)
   self%shellPoly = polyr(:, :maxElem)
   self%pairParam = kpair(:maxElem, :maxElem)
   self%referenceOcc = referenceOcc(:mShell, :maxElem)
   self%selfEnergy = ao_lev(:mShell, :maxElem)
   self%slaterExponent = ao_exp(:mShell, :maxElem)
   self%valenceShell = valenceShell(:mShell, :maxElem)
   self%principalQuantumNumber = ao_pqn(:mShell, :maxElem)

end subroutine initHamiltonian


subroutine initHalogen(self)

   !>
   type(THalogenData), intent(out) :: self

   self%atomicRad = atomicRad(:maxElem)

end subroutine initHalogen


end module xtb_xtb_gfn1
