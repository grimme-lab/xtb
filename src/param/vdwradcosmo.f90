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

!> Cosmo van-der-Waals radii

module xtb_param_vdwradcosmo
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_mctc_symbols, only : toNumber
   implicit none
   private

   public :: getVanDerWaalsRadCosmo, vanDerWaalsRadCosmo


   !> Get van-der-Waals Rad for a species
   interface getVanDerWaalsRadCosmo
      module procedure :: getVanDerWaalsRadCosmoSymbol
      module procedure :: getVanDerWaalsRadCosmoNumber
   end interface getVanDerWaalsRadCosmo


   !> Default value for unoptimized van-der-Waals radii
   real(wp), parameter :: cosmostub = 2.223_wp

   !> COSMO optimized van-der-Waals radii
   real(wp), parameter :: vanDerWaalsRadCosmo(94) = aatoau * [ &
       & 1.3000_wp, 1.6380_wp, 1.5700_wp, 1.0530_wp, &   ! h-be
       & 2.0480_wp, 2.0000_wp, 1.8300_wp, 1.7200_wp, &   ! B-O
       & 1.7200_wp, 1.8018_wp, 1.8000_wp, 1.6380_wp, &   ! F-Mg
       & 2.1530_wp, 2.2000_wp, 2.1060_wp, 2.1600_wp, &   ! Al-S
       & 2.0500_wp, 2.2000_wp, 2.2230_wp, cosmostub, &   ! Cl-Ca
       & cosmostub, 2.2930_wp, cosmostub, cosmostub, &   ! Sc-Cr
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Mn-Ni
       & cosmostub, 1.6260_wp, cosmostub, 2.7000_wp, &   ! Cu-Ge
       & 2.3500_wp, 2.2000_wp, 2.1600_wp, 2.3630_wp, &   ! As-Kr
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Rb-Zr
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Nb-Ru
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Rh-Cd
       & 2.2580_wp, 2.5500_wp, 2.4100_wp, 2.4100_wp, &   ! In-Te
       & 2.3200_wp, 2.5270_wp, cosmostub, cosmostub, &   ! I-Ba
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! La-Nd
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Pm-Gd
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Tb-Er
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Tm-Hf
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ta-Os
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ir-Hg
       & cosmostub, 2.3600_wp, 2.4220_wp, 2.3050_wp, &   ! Tl-Po
       & 2.3630_wp, 2.5740_wp, cosmostub, cosmostub, &   ! At-Ra
       & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ac-U
       & cosmostub, cosmostub]                           ! Np-Pu


contains


!> Get van-der-Waals radius for species with a given symbol
elemental function getVanDerWaalsRadCosmoSymbol(symbol) result(rad)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> van-der-Waals radius
   real(wp) :: rad

   rad = getVanDerWaalsRadCosmo(toNumber(symbol))

end function getVanDerWaalsRadCosmoSymbol


!> Get van-der-Waals radius for species with a given atomic number
elemental function getVanDerWaalsRadCosmoNumber(number) result(rad)

   !> Atomic number
   integer, intent(in) :: number

   !> van-der-Waals radius
   real(wp) :: rad

   if (number > 0 .and. number <= size(vanDerWaalsRadCosmo, dim=1)) then
      rad = vanDerWaalsRadCosmo(number)
   else
      rad = -1.0_wp
   end if

end function getVanDerWaalsRadCosmoNumber


end module xtb_param_vdwradcosmo
