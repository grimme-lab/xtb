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

!> Atomic r4/r2 expectation values
module xtb_param_sqrtzr4r2
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_symbols, only : toNumber
   implicit none
   private

   public :: getSqrtZr4r2, sqrtZr4r2


   !> Get r4/r2 expectation value for a species
   interface getSqrtZr4r2
      module procedure :: getSqrtZr4r2Symbol
      module procedure :: getSqrtZr4r2Number
   end interface getSqrtZr4r2

   integer, private, parameter :: max_elem = 118

   !  PBE0/def2-QZVP atomic values calculated by S. Grimme in Gaussian (2010)
   !  rare gases recalculated by J. Mewes with PBE0/aug-cc-pVQZ in Dirac (2018)
   !  He: 3.4698 -> 3.5544, Ne: 3.1036 -> 3.7943, Ar: 5.6004 -> 5.6638,
   !  Kr: 6.1971 -> 6.2312, Xe: 7.5152 -> 8.8367
   !  not replaced but recalculated (PBE0/cc-pVQZ) were
   !   H: 8.0589 ->10.9359, Li:29.0974 ->39.7226, Be:14.8517 ->17.7460
   !  also new super heavies Cn,Nh,Fl,Lv,Og
   real(wp), parameter :: r4Overr2(max_elem) = [  &
   &  8.0589_wp, 3.4698_wp, & ! H,He
   & 29.0974_wp,14.8517_wp,11.8799_wp, 7.8715_wp, 5.5588_wp, 4.7566_wp, 3.8025_wp, 3.1036_wp, & ! Li-Ne
   & 26.1552_wp,17.2304_wp,17.7210_wp,12.7442_wp, 9.5361_wp, 8.1652_wp, 6.7463_wp, 5.6004_wp, & ! Na-Ar
   & 29.2012_wp,22.3934_wp, & ! K,Ca
   &         19.0598_wp,16.8590_wp,15.4023_wp,12.5589_wp,13.4788_wp, & ! Sc-
   &         12.2309_wp,11.2809_wp,10.5569_wp,10.1428_wp, 9.4907_wp, & ! -Zn
   &                 13.4606_wp,10.8544_wp, 8.9386_wp, 8.1350_wp, 7.1251_wp, 6.1971_wp, & ! Ga-Kr
   & 30.0162_wp,24.4103_wp, & ! Rb,Sr
   &         20.3537_wp,17.4780_wp,13.5528_wp,11.8451_wp,11.0355_wp, & ! Y-
   &         10.1997_wp, 9.5414_wp, 9.0061_wp, 8.6417_wp, 8.9975_wp, & ! -Cd
   &                 14.0834_wp,11.8333_wp,10.0179_wp, 9.3844_wp, 8.4110_wp, 7.5152_wp, & ! In-Xe
   & 32.7622_wp,27.5708_wp, & ! Cs,Ba
   &         23.1671_wp,21.6003_wp,20.9615_wp,20.4562_wp,20.1010_wp,19.7475_wp,19.4828_wp, & ! La-Eu
   &         15.6013_wp,19.2362_wp,17.4717_wp,17.8321_wp,17.4237_wp,17.1954_wp,17.1631_wp, & ! Gd-Yb
   &         14.5716_wp,15.8758_wp,13.8989_wp,12.4834_wp,11.4421_wp, & ! Lu-
   &         10.2671_wp, 8.3549_wp, 7.8496_wp, 7.3278_wp, 7.4820_wp, & ! -Hg
   &                 13.5124_wp,11.6554_wp,10.0959_wp, 9.7340_wp, 8.8584_wp, 8.0125_wp, & ! Tl-Rn
   & 29.8135_wp,26.3157_wp, & ! Fr,Ra
   &         19.1885_wp,15.8542_wp,16.1305_wp,15.6161_wp,15.1226_wp,16.1576_wp, 0.0000_wp, & ! Ac-Am
   &          0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, & ! Cm-No
   &          0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, & ! Lr-
   &          0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 5.4929_wp, & ! -Cn
   &                  6.7286_wp, 6.5144_wp,10.9169_wp,10.3600_wp, 9.4723_wp, 8.6641_wp ] ! Nh-Og
   integer, private :: idum
   real(wp), parameter :: sqrtZr4r2(max_elem) = &
   &  sqrt(0.5_wp*(r4Overr2*[(sqrt(real(idum,wp)),idum=1,max_elem)]))


contains


!> Get r4/r2 expectation value for species with a given symbol
elemental function getSqrtZr4r2Symbol(symbol) result(mass)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> r4/r2 expectation value
   real(wp) :: mass

   mass = getSqrtZr4r2(toNumber(symbol))

end function getSqrtZr4r2Symbol


!> Get r4/r2 expectation value for species with a given atomic number
elemental function getSqrtZr4r2Number(number) result(mass)

   !> Atomic number
   integer, intent(in) :: number

   !> r4/r2 expectation value
   real(wp) :: mass

   if (number > 0 .and. number <= size(sqrtZr4r2, dim=1)) then
      mass = sqrtZr4r2(number)
   else
      mass = -1.0_wp
   end if

end function getSqrtZr4r2Number


end module xtb_param_sqrtzr4r2
