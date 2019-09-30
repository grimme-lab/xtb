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

module mctc_constants
   use iso_fortran_env, only : wp => real64
   implicit none
   private
   real(wp),public,parameter :: pi = 3.1415926535897932384626433832795029_wp
   real(wp),public,parameter :: pi4 = 3.1415926535897932384626433832795029_wp*4._wp
!  √π
   real(wp),public,parameter ::  sqrtpi  = sqrt(pi)
!  2×π
   real(wp),public,parameter :: twopi = 2.0_wp * pi
!  4×π
   real(wp),public,parameter :: fourpi = 4.0_wp * pi
!  π/2
   real(wp),public,parameter :: pihalf = 0.5_wp * pi
!  Boltzmann constant in Eh/K
   real(wp),public,parameter :: kB = 3.166808578545117e-06_wp
!  speed of light c in vacuum in a.u.
   real(wp),public,parameter :: lightspeed = 137.0359990740_wp
end module mctc_constants
