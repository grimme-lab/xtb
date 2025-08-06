! This file is part of xtb.
!
! Copyright (C) 2025 Igor S. Gerasimov
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

module random_generators
   use xtb_mctc_accuracy, only: sp, dp

   implicit none

   private

   interface normal_distribution
     module procedure normal_distribution_sp
     module procedure normal_distribution_dp
   end interface normal_distribution

   public normal_distribution

contains

   !>
   !> @brief   Generates random numbers according to the Normal (or Gaussian) random number distribution
   !>
   !> @details Implements Box-Muller transform for converting uniform distribution to normal
   !>
   !> @param[in] sigma   standard deviation
   !> @param[in] mu      mean of distribution
   !> @return            pseudorandom value
   real(sp) function normal_distribution_sp(sigma, mu) result(randval)
      real(sp), intent(in) :: sigma, mu
      real(sp), parameter :: two_pi = 2.0_sp * 4.0_sp * atan(1.0_sp)
      real(sp) :: u(2), mag

      u = 0.0_sp
      do while (u(1) == 0.0_sp)
         call random_number(u)
      end do
      mag = sigma * sqrt(-2.0_sp * log(u(1)))
      randval = mag * cos(two_pi * u(2)) + mu

   end function normal_distribution_sp

   !>
   !> @brief   Generates random numbers according to the Normal (or Gaussian) random number distribution
   !>
   !> @details Implements Box-Muller transform for converting uniform distribution to normal
   !>
   !> @param[in] sigma   standard deviation
   !> @param[in] mu      mean of distribution
   !> @return            pseudorandom value
   real(dp) function normal_distribution_dp(sigma, mu) result(randval)
      real(dp), intent(in) :: sigma, mu
      real(dp), parameter :: two_pi = 2.0_dp * 4.0_dp * atan(1.0_dp)
      real(dp) :: u(2), mag

      u = 0.0_dp
      do while (u(1) == 0.0_dp)
         call random_number(u)
      end do
      mag = sigma * sqrt(-2.0_dp * log(u(1)))
      randval = mag * cos(two_pi * u(2)) + mu

   end function normal_distribution_dp

end module random_generators
