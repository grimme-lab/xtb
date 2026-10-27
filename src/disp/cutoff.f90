! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

!> Smooth real space cutoff for the dispersion interactions, identical to the
!> switching function used in dftd4, s-dftd3 and tblite
module xtb_disp_cutoff
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: smooth_cutoff


   !> Coefficients of the quintic smoothstep polynomial 10*x^3 - 15*x^4 + 6*x^5
   real(wp), parameter :: smoothstep3 = 10.0_wp
   real(wp), parameter :: smoothstep4 = -15.0_wp
   real(wp), parameter :: smoothstep5 = 6.0_wp
   real(wp), parameter :: smoothstep_deriv = 30.0_wp


contains


!> Smooth polynomial switch for real space cutoffs, switches from one at
!> cutoff - width to zero at cutoff, a width of zero gives a hard cutoff
pure subroutine smooth_cutoff(r, cutoff, width, sw, dswdr)
   !$acc routine seq

   !> Interatomic distance
   real(wp), intent(in) :: r

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Width of smooth cutoff
   real(wp), intent(in) :: width

   !> Switching function value
   real(wp), intent(out) :: sw

   !> Derivative of the switching function with respect to distance
   real(wp), intent(out) :: dswdr

   real(wp) :: inner, effective_width, x

   if (width <= 0.0_wp .or. cutoff <= 0.0_wp) then
      sw = 1.0_wp
      dswdr = 0.0_wp
   else
      effective_width = min(width, cutoff)
      inner = cutoff - effective_width
      if (r <= inner) then
         sw = 1.0_wp
         dswdr = 0.0_wp
      else if (r >= cutoff) then
         sw = 0.0_wp
         dswdr = 0.0_wp
      else
         x = (cutoff - r) / effective_width
         sw = x**3 * (smoothstep3 + x*(smoothstep4 + smoothstep5*x))
         dswdr = -smoothstep_deriv * x**2 * (1.0_wp - x)**2 / effective_width
      end if
   end if

end subroutine smooth_cutoff


end module xtb_disp_cutoff
