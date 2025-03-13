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

module xtb_mctc_constants
   use xtb_mctc_accuracy, only : wp
   implicit none
   private
   real(wp),public,parameter :: pi = 3.1415926535897932384626433832795029_wp
!  √π
   real(wp),public,parameter :: sqrtpi = sqrt(pi)
!  2×π
   real(wp),public,parameter :: twopi = 2.0_wp * pi
!  4×π
   real(wp),public,parameter :: fourpi = 4.0_wp * pi
!  π/2
   real(wp),public,parameter :: pihalf = 0.5_wp * pi
!  2π/3
   real(wp),public,parameter :: twothirdpi = 2.0_wp * pi / 3.0_wp
!  Boltzmann constant in Eh/K
   real(wp),public,parameter :: kB = 3.166808578545117e-06_wp
!  speed of light c in vacuum in a.u.
   real(wp),public,parameter :: lightspeed = 137.0359990740_wp
   ! Boltzmann constant in SI units (codata 2018)
   real(wp),public,parameter :: kB_SI = 1.380649e-23_wp
   ! Speed of light in SI units
   real(wp),public,parameter :: lightspeed_SI = 299792458.0_wp
   ! Planck constant in SI units (codata 2018)
   real(wp),public,parameter :: h_SI = 6.62607015e-34_wp
   !> standard atmosphere (Pa) (codata 2018)
   real(wp),public,parameter :: standard_atmosphere = 101325e0_wp
   !> Bohr radius (m) (codata 2018)
   real(wp),public,parameter :: Bohr_radius = 5.29177210903e-11_wp
   !> molar gas constant (J mol^-1 K^-1) (codata 2018)
   real(wp),public,parameter :: molar_gas_constant = 8.314462618_wp

end module xtb_mctc_constants
