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
   use xtb_param_atomicrad, only : atomicRad
   use xtb_param_paulingen, only : paulingEN
   use xtb_xtb_data
   use xtb_aoparam
   implicit none
   private

   public :: initGFN2


   interface initGFN2
      module procedure :: initData
      module procedure :: initRepulsion
      module procedure :: initCoulomb
      module procedure :: initMultipole
      module procedure :: initHamiltonian
   end interface initGFN2

   !> Maximum number of elements supported by GFN2-xTB
   integer, parameter :: maxElem = 86

   ! ========================================================================
   ! REPULSION DATA
   !>
   real(wp), parameter :: kExp = 1.5_wp

   !>
   real(wp), parameter :: kExpLight = 1.0_wp

   !>
   real(wp), parameter :: rExp = 1.0_wp

   !>
   real(wp), parameter :: repAlpha(1:maxElem) = [&
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
   real(wp), parameter :: repZeff(1:maxElem) = [&
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

   ! ========================================================================
   ! COULOMB DATA
   !>
   real(wp), parameter :: chemicalHardness(1:maxElem) = [&
      & 0.405771_wp, 0.642029_wp, 0.245006_wp, 0.684789_wp, 0.513556_wp, &
      & 0.538015_wp, 0.461493_wp, 0.451896_wp, 0.531518_wp, 0.850000_wp, &
      & 0.271056_wp, 0.344822_wp, 0.364801_wp, 0.720000_wp, 0.297739_wp, &
      & 0.339971_wp, 0.248514_wp, 0.502376_wp, 0.247602_wp, 0.320378_wp, &
      & 0.472633_wp, 0.513586_wp, 0.589187_wp, 0.396299_wp, 0.346651_wp, &
      & 0.271594_wp, 0.477760_wp, 0.344970_wp, 0.202969_wp, 0.564152_wp, &
      & 0.432236_wp, 0.802051_wp, 0.571748_wp, 0.235052_wp, 0.261253_wp, &
      & 0.424373_wp, 0.210481_wp, 0.340000_wp, 0.711958_wp, 0.461440_wp, &
      & 0.952957_wp, 0.586134_wp, 0.368054_wp, 0.711205_wp, 0.509183_wp, &
      & 0.273310_wp, 0.263740_wp, 0.392012_wp, 0.461812_wp, 0.900000_wp, &
      & 0.942294_wp, 0.750000_wp, 0.383124_wp, 0.424164_wp, 0.236569_wp, &
      & 0.245937_wp, 0.597716_wp, 0.662889_wp, 0.660710_wp, 0.658531_wp, &
      & 0.656352_wp, 0.654173_wp, 0.651994_wp, 0.649815_wp, 0.647635_wp, &
      & 0.645456_wp, 0.643277_wp, 0.641098_wp, 0.638919_wp, 0.636740_wp, &
      & 0.634561_wp, 0.662597_wp, 0.449812_wp, 0.685426_wp, 0.224623_wp, &
      & 0.364388_wp, 0.548507_wp, 0.353574_wp, 0.438997_wp, 0.457611_wp, &
      & 0.418841_wp, 0.168152_wp, 0.900000_wp, 1.023267_wp, 0.288848_wp, &
      & 0.303400_wp]

   !>
   logical, parameter :: thirdOrderShellResolved = .true.

   !>
   real(wp), parameter :: thirdOrderAtom(1:maxElem) = [&
      & 0.800000_wp, 2.000000_wp, 1.303821_wp, 0.574239_wp, 0.946104_wp, &
      & 1.500000_wp,-0.639780_wp,-0.517134_wp, 1.426212_wp, 0.500000_wp, &
      & 1.798727_wp, 2.349164_wp, 1.400000_wp, 1.936289_wp, 0.711291_wp, &
      &-0.501722_wp, 1.495483_wp,-0.315455_wp, 2.033085_wp, 2.006898_wp, &
      & 0.500000_wp, 1.767268_wp, 0.900000_wp, 0.300000_wp, 0.600000_wp, &
      &-0.500000_wp, 0.300000_wp,-0.200000_wp, 0.500000_wp, 2.312896_wp, &
      & 2.334269_wp,-0.064775_wp, 1.106041_wp, 0.913725_wp, 1.300000_wp, &
      & 0.239815_wp, 2.916203_wp, 1.800000_wp, 0.100000_wp, 0.700000_wp, &
      & 0.500000_wp, 0.919928_wp, 0.600000_wp,-0.500000_wp, 0.300000_wp, &
      & 0.800000_wp, 0.200000_wp, 2.073217_wp, 1.900000_wp,-0.178396_wp, &
      & 1.100000_wp, 0.953683_wp, 1.200000_wp,-0.118925_wp, 2.404185_wp, &
      & 2.069097_wp, 0.012793_wp,-0.100000_wp,-0.100002_wp,-0.100004_wp, &
      &-0.100006_wp,-0.100008_wp,-0.100010_wp,-0.100012_wp,-0.100013_wp, &
      &-0.100015_wp,-0.100017_wp,-0.100019_wp,-0.100021_wp,-0.100023_wp, &
      &-0.100025_wp,-0.100000_wp, 0.200000_wp,-0.200000_wp, 0.800000_wp, &
      & 0.800000_wp,-0.100000_wp, 0.600000_wp, 0.850000_wp,-0.116312_wp, &
      &-0.533933_wp, 0.200000_wp,-0.337508_wp, 1.877978_wp, 1.846485_wp, &
      & 0.097834_wp] * 0.1_wp

   ! ========================================================================
   ! MULTIPOLE DATA
   !>
   real(wp), parameter :: dipDamp = 3.0_wp

   !>
   real(wp), parameter :: quadDamp = 4.0_wp

   !>
   real(wp), parameter :: cnShift = 1.2_wp

   !>
   real(wp), parameter :: cnExp = 4.0_wp

   !>
   real(wp), parameter :: cnRMax = 5.0_wp

   !>
   real(wp), parameter :: dipKernel(1:maxElem) = [&
      & 5.563889_wp,-1.000000_wp,-0.500000_wp,-0.613341_wp,-0.481186_wp, &
      &-0.411674_wp, 3.521273_wp,-4.935670_wp,-8.339183_wp,10.000000_wp, &
      & 0.000000_wp,-0.082005_wp, 2.633341_wp,-0.025750_wp, 2.110225_wp, &
      &-0.151117_wp,-2.536958_wp,-2.077329_wp,-0.103383_wp,-0.236675_wp, &
      &-0.515177_wp,-0.434506_wp,-0.350000_wp, 0.149669_wp,-0.759168_wp, &
      & 0.412929_wp,-0.247938_wp,-1.261887_wp,-0.700000_wp,-0.100000_wp, &
      & 0.267219_wp, 0.108460_wp,-0.201294_wp,-0.288648_wp,-1.088586_wp, &
      &-0.889357_wp,-0.093328_wp,-0.459925_wp,-0.637291_wp,-0.599615_wp, &
      &-0.288729_wp, 0.346327_wp,-0.458416_wp,-0.081922_wp, 0.007016_wp, &
      &-0.310361_wp,-0.800314_wp,-0.105364_wp, 0.951079_wp, 0.085029_wp, &
      &-0.015519_wp,-0.263414_wp,-0.603648_wp,-0.214447_wp,-0.080000_wp, &
      &-0.260000_wp,-0.395198_wp,-0.723806_wp,-0.704819_wp,-0.685832_wp, &
      &-0.666845_wp,-0.647858_wp,-0.628871_wp,-0.609884_wp,-0.590897_wp, &
      &-0.571910_wp,-0.552923_wp,-0.533936_wp,-0.514949_wp,-0.495961_wp, &
      &-0.476974_wp,-0.537685_wp,-0.200343_wp, 0.065886_wp,-0.587636_wp, &
      &-0.510090_wp,-0.673822_wp,-0.423684_wp, 0.393418_wp,-0.250000_wp, &
      & 0.374018_wp, 1.007016_wp,-0.737252_wp,-1.344854_wp,-0.348123_wp, &
      &-0.167597_wp] * 0.01_wp

   !>
   real(wp), parameter :: quadKernel(1:maxElem) = [&
      & 0.027431_wp,-0.337528_wp, 0.020000_wp,-0.058586_wp,-0.058228_wp, &
      & 0.213583_wp, 2.026786_wp,-0.310828_wp,-0.245955_wp,-0.500000_wp, &
      & 0.020000_wp,-0.005516_wp,-0.021887_wp,-0.080000_wp, 0.028679_wp, &
      & 0.442859_wp, 0.122783_wp,-1.083404_wp, 0.025000_wp, 0.010000_wp, &
      &-0.042004_wp, 0.059660_wp, 0.009764_wp, 0.137744_wp, 0.229903_wp, &
      & 0.267734_wp, 0.048237_wp,-0.080000_wp,-0.345631_wp, 0.007658_wp, &
      &-0.003616_wp,-0.003589_wp, 0.014149_wp, 0.085728_wp, 0.216935_wp, &
      &-0.415024_wp, 0.015000_wp, 0.015000_wp, 0.010460_wp,-0.012944_wp, &
      & 0.041491_wp, 0.312549_wp, 0.155242_wp, 0.359228_wp, 0.008570_wp, &
      &-0.040485_wp,-0.020810_wp, 0.012250_wp,-0.002031_wp,-0.008243_wp, &
      &-0.020630_wp,-0.026864_wp, 0.069660_wp,-0.156200_wp, 0.008000_wp, &
      & 0.015000_wp,-0.030000_wp,-0.025000_wp,-0.024615_wp,-0.024231_wp, &
      &-0.023846_wp,-0.023462_wp,-0.023077_wp,-0.022692_wp,-0.022308_wp, &
      &-0.021923_wp,-0.021538_wp,-0.021154_wp,-0.020769_wp,-0.020385_wp, &
      &-0.020000_wp,-0.016478_wp, 0.039599_wp, 1.063309_wp, 0.306870_wp, &
      & 0.759049_wp, 0.322935_wp, 0.098019_wp,-0.020320_wp,-0.032901_wp, &
      &-0.008506_wp,-0.001670_wp, 0.162529_wp, 0.013818_wp, 0.021624_wp, &
      &-0.111556_wp] * 0.01_wp

   !>
   real(wp), parameter :: valenceCN(1:maxElem) = [&
      & 1.0_wp, 1.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 2.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, &
      & 4.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, &
      & 4.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, &
      & 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 4.0_wp, &
      & 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 4.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, &
      & 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp]

   !>
   real(wp), parameter :: multiRad(1:maxElem) = [&
      & 1.4_wp, 3.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.0_wp, 1.9_wp, 1.8_wp, 2.4_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.9_wp, 2.1_wp, 3.1_wp, 2.5_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.9_wp, 4.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp]

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
   real(wp), parameter :: referenceOcc(0:2, 1:maxElem) = reshape([&
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp,  1.0_wp, 3.0_wp, 0.0_wp, &
      & 1.5_wp, 3.5_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp, &
      & 2.0_wp, 6.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  1.5_wp, 2.5_wp, 0.0_wp,  1.5_wp, 3.5_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 2.0_wp,  2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp, &
      & 2.0_wp, 0.0_wp, 5.0_wp,  2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp, &
      & 2.0_wp, 0.0_wp, 8.0_wp,  1.0_wp, 0.0_wp,10.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  1.5_wp, 2.5_wp, 0.0_wp,  1.5_wp, 3.5_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 2.0_wp,  2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp, &
      & 2.0_wp, 0.0_wp, 5.0_wp,  2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp, &
      & 2.0_wp, 0.0_wp, 8.0_wp,  1.0_wp, 0.0_wp,10.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
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
      & 1.0_wp, 0.0_wp,10.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp, &
      & 2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp, &
      & 2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp], shape(referenceOcc))

contains


subroutine initData(self)

   !>
   type(TxTBData), intent(out) :: self

   self%nShell = ao_n(:maxElem)

   call initGFN2(self%repulsion)
   call initGFN2(self%coulomb, self%nShell)
   allocate(self%multipole)
   call initGFN2(self%multipole)
   call initGFN2(self%hamiltonian, self%nShell)

end subroutine initData


subroutine initRepulsion(self)

   !>
   type(TRepulsionData), intent(out) :: self

   self%cutoff = 40.0_wp
   self%kExp = kExp
   self%kExpLight = kExpLight
   self%rExp = rExp
   self%electronegativity = spread(1.0_wp, 1, maxElem)
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


subroutine initMultipole(self)

   !>
   type(TMultipoleData), intent(out) :: self

   self%cnShift = cnShift
   self%cnExp = cnExp
   self%cnRMax = cnRMax
   self%dipDamp = dipDamp
   self%quadDamp = quadDamp
   self%dipKernel = dpolc(:maxElem) ! dipKernel
   self%quadKernel = qpolc(:maxElem) ! quadKernel
   self%valenceCN = valenceCN
   self%multiRad = multiRad

end subroutine initMultipole


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
   self%kCN = kcnat(:, :maxElem)
   self%selfEnergy = ao_lev(:mShell, :maxElem)
   self%slaterExponent = ao_exp(:mShell, :maxElem)
   self%principalQuantumNumber = ao_pqn(:mShell, :maxElem)

   allocate(self%valenceShell(mShell, maxElem))
   call generateValenceShellData(self%valenceShell, nShell, self%angShell)

   allocate(self%referenceOcc(mShell, maxElem))
   self%referenceOcc(:, :) = 0.0_wp
   do iZp = 1, maxElem
      do iSh = 1, nShell(iZp)
         lAng = self%angShell(iSh, iZp)
         if (self%valenceShell(iSh, iZp) /= 0) then
            self%referenceOcc(iSh, iZp) = referenceOcc(lAng, iZp)
         end if
      end do
   end do

   allocate(self%numberOfPrimitives(mShell, maxElem))
   do iZp = 1, maxElem
      do iSh = 1, nShell(iZp)
         nPrim = 0
         if (iZp <= 2) then
            select case(self%angShell(iSh, iZp))
            case(0)
               nPrim = 3
            case(1)
               nPrim = 4
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
                  nPrim = 4
               end if
            case(2)
               nPrim = 3
            case(3)
               nPrim = 4
            end select
         end if
         self%numberOfPrimitives(iSh, iZp) = nPrim
      end do
   end do

end subroutine initHamiltonian


end module xtb_xtb_gfn2
