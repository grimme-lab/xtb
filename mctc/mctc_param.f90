! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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

module mctc_param
   use iso_fortran_env, wp => real64
   use mctc_econv

   implicit none

   public :: rad
   public :: atomic_mass

   private

   integer, parameter :: max_elem = 118

!  covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
!  188-197), values for metals decreased by 10 %
   real(wp),parameter :: rad(max_elem) = [ &
   & 0.32,0.46, & ! H,He
   & 1.20,0.94,0.77,0.75,0.71,0.63,0.64,0.67, & ! Li-Ne
   & 1.40,1.25,1.13,1.04,1.10,1.02,0.99,0.96, & ! Na-Ar
   & 1.76,1.54, & ! K,Ca
   &           1.33,1.22,1.21,1.10,1.07,1.04,1.00,0.99,1.01,1.09, & ! Sc-Zn
   &           1.12,1.09,1.15,1.10,1.14,1.17, & ! Ga-Kr
   & 1.89,1.67, & ! Rb,Sr
   &           1.47,1.39,1.32,1.24,1.15,1.13,1.13,1.08,1.15,1.23, & ! Y-Cd
   &           1.28,1.26,1.26,1.23,1.32,1.31, & ! In-Xe
   & 2.09,1.76, & ! Cs,Ba
   &      1.62,1.47,1.58,1.57,1.56,1.55,1.51, & ! La-Eu
   &      1.52,1.51,1.50,1.49,1.49,1.48,1.53, & ! Gd-Yb
   &           1.46,1.37,1.31,1.23,1.18,1.16,1.11,1.12,1.13,1.32, & ! Lu-Hg
   &           1.30,1.30,1.36,1.31,1.38,1.42, & ! Tl-Rn
   & 2.01,1.81, & ! Fr,Ra
   &      1.67,1.58,1.52,1.53,1.54,1.55,1.49, & ! Ac-Am
   &      1.49,1.51,1.51,1.48,1.50,1.56,1.58, & ! Cm-No
   &           1.45,1.41,1.34,1.29,1.27,1.21,1.16,1.15,1.09,1.22, & ! Lr-Cn
   &           1.36,1.43,1.46,1.58,1.48,1.57 ] * aatoau ! Nh-Og

!  covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
!  188-197), values for metals decreased by 10 %
   real(wp),parameter :: covalent_radius_2009(max_elem) = aatoau * [ &
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
   &      1.67_wp,1.58_wp,1.52_wp,1.53_wp,1.54_wp,1.55_wp,1.49_wp, & ! Ac-Am
   &      1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
   &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
   &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
   &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp ] ! Nh-Og

   real(wp),parameter :: atomic_mass_nist(max_elem) = amutoau * [ &
      &  1.00794075_wp,  4.00260193_wp,  6.94003660_wp,  9.01218307_wp,&
      & 10.81102805_wp, 12.01073590_wp, 14.00670321_wp, 15.99940492_wp,&
      & 18.99840316_wp, 20.18004638_wp, 22.98976928_wp, 24.30505162_wp,&
      & 26.98153853_wp, 28.08549871_wp, 30.97376200_wp, 32.06478741_wp,&
      & 35.45293758_wp, 39.94779856_wp, 39.09830091_wp, 40.07802251_wp,&
      & 44.95590828_wp, 47.86674496_wp, 50.94146504_wp, 51.99613176_wp,&
      & 54.93804391_wp, 55.84514443_wp, 58.93319429_wp, 58.69334711_wp,&
      & 63.54603995_wp, 65.37778253_wp, 69.72306607_wp, 72.62755016_wp,&
      & 74.92159457_wp, 78.95938856_wp, 79.90352778_wp, 83.79800000_wp,&
      & 85.46766360_wp, 87.61664447_wp, 88.90584030_wp, 91.22364160_wp,&
      & 92.90637300_wp, 95.95978854_wp, 97.90721240_wp,101.06494014_wp,&
      &102.90549800_wp,106.41532751_wp,107.86814963_wp,112.41155782_wp,&
      &114.81808663_wp,118.71011259_wp,121.75978367_wp,127.60312648_wp,&
      &126.90447190_wp,131.29276145_wp,132.90545196_wp,137.32689163_wp,&
      &138.90546887_wp,140.11573074_wp,140.90765760_wp,144.24159603_wp,&
      &144.91275590_wp,150.36635571_wp,151.96437813_wp,157.25213065_wp,&
      &158.92535470_wp,162.49947282_wp,164.93032880_wp,167.25908265_wp,&
      &168.93421790_wp,173.05415017_wp,174.96681496_wp,178.48497872_wp,&
      &180.94787564_wp,183.84177755_wp,186.20670455_wp,190.22485963_wp,&
      &192.21605165_wp,195.08445686_wp,196.96656879_wp,200.59916703_wp,&
      &204.38341284_wp,207.21690806_wp,208.98039910_wp,208.98243080_wp,&
      &209.98714790_wp,222.01757820_wp,223.01973600_wp,226.02541030_wp,&
      &227.02775230_wp,232.03805580_wp,231.03588420_wp,238.02891046_wp,&
      &237.04817360_wp,244.06420530_wp,243.06138130_wp,247.07035410_wp,&
      &247.07030730_wp,251.07958860_wp,252.08298000_wp,257.09510610_wp,&
      &258.09843150_wp,259.10103000_wp,262.10961000_wp,267.12179000_wp,&
      &269.12791000_wp,271.13393000_wp,270.13336000_wp,276.14846000_wp,&
      &276.15159000_wp,280.16131000_wp,282.16912000_wp,284.17416000_wp,&
      &284.17873000_wp,289.19042000_wp,288.19274000_wp,293.20449000_wp,&
      &292.20746000_wp,294.21392000_wp]

contains

pure elemental function covalent_radius(iatom) result(rcov)
   implicit none
   integer,intent(in) :: iatom
   real(wp) :: rcov

   if (iatom > 0 .and. iatom <= max_elem) then
      rcov = covalent_radius_2009(iatom)
   else
      rcov = 0.0_wp
   endif
end function covalent_radius

pure elemental function atomic_mass(iatom) result(mass)
   implicit none
   integer,intent(in) :: iatom
   real(wp) :: mass

   if (iatom > 0 .and. iatom <= max_elem) then
      mass = atomic_mass_nist(iatom)
   else
      mass = 0.0_wp
   endif
end function atomic_mass



end module mctc_param
