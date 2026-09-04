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
!> UFF vdw radii - could be replaced with any other vdw radii i guess
module xtb_param_uffvdwrad
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autoaa
   implicit none
   private

   public :: get_rad

   real(wp), parameter :: vdw_radii(103) = [ &
      2.886_wp, 2.362_wp, 2.451_wp, 2.745_wp, 4.083_wp, 3.851_wp, 3.66_wp, 3.5_wp, 3.364_wp, &
      3.243_wp, 2.983_wp, 3.021_wp, 4.499_wp, 4.295_wp, 4.147_wp, 4.035_wp, 3.947_wp, 3.868_wp, &
      3.812_wp, 3.399_wp, 3.295_wp, 3.175_wp, 3.144_wp, 3.023_wp, 2.961_wp, 2.912_wp, 2.872_wp, &
      2.834_wp, 3.495_wp, 2.763_wp, 4.383_wp, 4.28_wp, 4.23_wp, 4.205_wp, 4.189_wp, 4.141_wp, &
      4.114_wp, 3.641_wp, 3.345_wp, 3.124_wp, 3.165_wp, 3.052_wp, 2.998_wp, 2.963_wp, 2.929_wp, &
      2.899_wp, 3.148_wp, 2.848_wp, 4.463_wp, 4.392_wp, 4.42_wp, 4.47_wp, 4.5_wp, 4.404_wp, &
      4.517_wp, 3.703_wp, 3.522_wp, 3.556_wp, 3.606_wp, 3.575_wp, 3.547_wp, 3.52_wp, 3.493_wp, &
      3.368_wp, 3.451_wp, 3.428_wp, 3.409_wp, 3.391_wp, 3.374_wp, 3.355_wp, 3.64_wp, 3.141_wp, &
      3.17_wp, 3.069_wp, 2.954_wp, 3.12_wp, 2.84_wp, 2.754_wp, 3.293_wp, 2.705_wp, 4.347_wp, &
      4.297_wp, 4.37_wp, 4.709_wp, 4.75_wp, 4.765_wp, 4.9_wp, 3.677_wp, 3.478_wp, 3.396_wp, &
      3.424_wp, 3.395_wp, 3.424_wp, 3.424_wp, 3.381_wp, 3.326_wp, 3.339_wp, 3.313_wp, 3.299_wp, &
      3.286_wp, 3.274_wp, 3.248_wp, 3.236_wp] * 0.5_wp / autoaa

   contains

function get_rad(at) result(rad)
   integer, intent(in) :: at

   real(wp) :: rad

   if (at > 103 .or. at < 1) then
      rad = -1.0_wp
   else
      rad = vdw_radii(at)
   end if

end function get_rad
end module xtb_param_uffvdwrad
