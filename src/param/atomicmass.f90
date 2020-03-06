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

!> Atomic masses
module xtb_param_atomicmass
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : amutoau
   use xtb_mctc_symbols, only : toNumber
   implicit none
   private

   public :: getAtomicMass


   !> Get atomic mass for a species
   interface getAtomicMass
      module procedure :: getAtomicMassSymbol
      module procedure :: getAtomicMassNumber
   end interface getAtomicMass


   !> Atomic masses from NIST
   real(wp), parameter :: atomicMass(1:118) = amutoau * [ &
      &   1.00794075_wp,   4.00260193_wp,   6.94003660_wp,   9.01218307_wp,&
      &  10.81102805_wp,  12.01073590_wp,  14.00670321_wp,  15.99940492_wp,&
      &  18.99840316_wp,  20.18004638_wp,  22.98976928_wp,  24.30505162_wp,&
      &  26.98153853_wp,  28.08549871_wp,  30.97376200_wp,  32.06478741_wp,&
      &  35.45293758_wp,  39.94779856_wp,  39.09830091_wp,  40.07802251_wp,&
      &  44.95590828_wp,  47.86674496_wp,  50.94146504_wp,  51.99613176_wp,&
      &  54.93804391_wp,  55.84514443_wp,  58.93319429_wp,  58.69334711_wp,&
      &  63.54603995_wp,  65.37778253_wp,  69.72306607_wp,  72.62755016_wp,&
      &  74.92159457_wp,  78.95938856_wp,  79.90352778_wp,  83.79800000_wp,&
      &  85.46766360_wp,  87.61664447_wp,  88.90584030_wp,  91.22364160_wp,&
      &  92.90637300_wp,  95.95978854_wp,  97.90721240_wp, 101.06494014_wp,&
      & 102.90549800_wp, 106.41532751_wp, 107.86814963_wp, 112.41155782_wp,&
      & 114.81808663_wp, 118.71011259_wp, 121.75978367_wp, 127.60312648_wp,&
      & 126.90447190_wp, 131.29276145_wp, 132.90545196_wp, 137.32689163_wp,&
      & 138.90546887_wp, 140.11573074_wp, 140.90765760_wp, 144.24159603_wp,&
      & 144.91275590_wp, 150.36635571_wp, 151.96437813_wp, 157.25213065_wp,&
      & 158.92535470_wp, 162.49947282_wp, 164.93032880_wp, 167.25908265_wp,&
      & 168.93421790_wp, 173.05415017_wp, 174.96681496_wp, 178.48497872_wp,&
      & 180.94787564_wp, 183.84177755_wp, 186.20670455_wp, 190.22485963_wp,&
      & 192.21605165_wp, 195.08445686_wp, 196.96656879_wp, 200.59916703_wp,&
      & 204.38341284_wp, 207.21690806_wp, 208.98039910_wp, 208.98243080_wp,&
      & 209.98714790_wp, 222.01757820_wp, 223.01973600_wp, 226.02541030_wp,&
      & 227.02775230_wp, 232.03805580_wp, 231.03588420_wp, 238.02891046_wp,&
      & 237.04817360_wp, 244.06420530_wp, 243.06138130_wp, 247.07035410_wp,&
      & 247.07030730_wp, 251.07958860_wp, 252.08298000_wp, 257.09510610_wp,&
      & 258.09843150_wp, 259.10103000_wp, 262.10961000_wp, 267.12179000_wp,&
      & 269.12791000_wp, 271.13393000_wp, 270.13336000_wp, 276.14846000_wp,&
      & 276.15159000_wp, 280.16131000_wp, 282.16912000_wp, 284.17416000_wp,&
      & 284.17873000_wp, 289.19042000_wp, 288.19274000_wp, 293.20449000_wp,&
      & 292.20746000_wp, 294.21392000_wp]


contains


!> Get atomic mass for species with a given symbol
elemental function getAtomicMassSymbol(symbol) result(mass)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> Atomic mass
   real(wp) :: mass

   mass = getAtomicMass(toNumber(symbol))

end function getAtomicMassSymbol


!> Get atomic mass for species with a given atomic number
elemental function getAtomicMassNumber(number) result(mass)

   !> Atomic number
   integer, intent(in) :: number

   !> Atomic mass
   real(wp) :: mass

   if (number > 0 .and. number <= size(atomicMass, dim=1)) then
      mass = atomicMass(number)
   else
      mass = -1.0_wp
   end if

end function getAtomicMassNumber


end module xtb_param_atomicmass
