! This file is part of xtb.
!
! Copyright (C) 2026-2027 Leopold M. Seidler
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
!
!> Unmodified covalent radii from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
!> 188-197
module xtb_param_covalentrad
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autoaa
   implicit none
   private

   public :: get_cov_rad

   !> covalent radii
   real(wp), parameter :: cov(103) = [&
      & 0.32_wp, 0.46_wp, 1.33_wp, 1.02_wp, 0.85_wp, 0.75_wp, 0.71_wp, 0.63_wp, 0.64_wp, 0.67_wp, &
      & 1.55_wp, 1.39_wp, 1.26_wp, 1.16_wp, 1.11_wp, 1.03_wp, 0.99_wp, 0.96_wp, 1.96_wp, 1.71_wp, 1.48_wp, &
      & 1.36_wp, 1.34_wp, 1.22_wp, 1.19_wp, 1.16_wp, 1.10_wp, 1.11_wp, 1.12_wp, 1.18_wp, 1.24_wp, 1.21_wp, &
      & 1.21_wp, 1.16_wp, 1.14_wp, 1.17_wp, 2.10_wp, 1.85_wp, 1.63_wp, 1.54_wp, 1.47_wp, 1.38_wp, 1.28_wp, &
      & 1.25_wp, 1.25_wp, 1.20_wp, 1.28_wp, 1.36_wp, 1.42_wp, 1.40_wp, 1.40_wp, 1.36_wp, 1.33_wp, 1.31_wp, &
      & 2.32_wp, 1.96_wp, 1.80_wp, 1.63_wp, 1.76_wp, 1.74_wp, 1.73_wp, 1.72_wp, 1.68_wp, 1.69_wp, 1.68_wp, &
      & 1.67_wp, 1.66_wp, 1.65_wp, 1.64_wp, 1.70_wp, 1.62_wp, 1.52_wp, 1.46_wp, 1.37_wp, 1.31_wp, 1.29_wp, &
      & 1.22_wp, 1.23_wp, 1.24_wp, 1.33_wp, 1.44_wp, 1.44_wp, 1.51_wp, 1.45_wp, 1.47_wp, 1.42_wp, 2.23_wp, &
      & 2.01_wp, 1.86_wp, 1.75_wp, 1.69_wp, 1.70_wp, 1.71_wp, 1.72_wp, 1.66_wp, 1.66_wp, 1.68_wp, 1.68_wp, &
      & 1.65_wp, 1.67_wp, 1.73_wp, 1.76_wp, 1.61_wp] / autoaa

   ! NOTE: these should later be moved to mctc-lib
   contains

function get_cov_rad(at) result(rad)
   integer, intent(in) :: at

   real(wp) :: rad

   if (at > 103 .or. at < 1) then
      rad = -1.0_wp
   else
      rad = cov(at)
   end if

end function get_cov_rad
end module xtb_param_covalentrad
