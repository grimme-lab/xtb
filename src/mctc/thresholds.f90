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

!> Numerical thresholds shared between multiple algorithms
module xtb_mctc_thresholds
   use xtb_mctc_accuracy, only : wp
   implicit none
   public


   !> Tolerance for atomic distances
   real(wp), parameter :: tolSameDist = 1.0e-5_wp

   !> Tolerance for atomic square distances
   real(wp), parameter :: tolSameDist2 = tolSameDist**2

   !> Minimal distance between neighbours
   real(wp), parameter :: minNeighDist = 1.0e-2_wp

   !> Minimal square distance between neighbours
   real(wp), parameter :: minNeighDist2 = minNeighDist**2


contains


end module xtb_mctc_thresholds
