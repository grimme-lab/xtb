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

module mctc_econv
   use iso_fortran_env, only : wp => real64
   use mctc_constants
   implicit none
   private
!  convert bohr (a.u.) to Ångström and back
   real(wp),public,parameter :: autoaa = 0.52917726_wp
   real(wp),public,parameter :: aatoau = 1.0_wp/autoaa
!  convert Hartree to eV and back
   real(wp),public,parameter :: autoev = 27.21138505_wp
   real(wp),public,parameter :: evtoau = 1.0_wp/autoev
!  convert Hartree to kcal/mol and back
   real(wp),public,parameter :: autokcal = 627.50947428_wp
   real(wp),public,parameter :: kcaltoau = 1.0_wp/autokcal
!  convert Hartree to kJ/mol and back
   real(wp),public,parameter :: autokj = 2625.49964038_wp
   real(wp),public,parameter :: kjtoau = 1.0_wp/autokj
!  convert Hartree to reciproce centimeters/wavenumbers and back
   real(wp),public,parameter :: autorcm = 219474.63067_wp
   real(wp),public,parameter :: autowav = autorcm
   real(wp),public,parameter :: rcmtoau = 1.0_wp/autorcm
   real(wp),public,parameter :: wavtoau = 1.0_wp/autowav
!  convert Hartree to nanometers and back
   real(wp),public,parameter :: autonm = 45.56335266_wp
   real(wp),public,parameter :: nmtoau = 1.0_wp/autonm
!  masses
!  amu -> kg :: conversion from atomic mass units to kg
!  me  -> kg :: electron mass (a.u.) in kg
!  amu -> au :: conversion from a.u. to amu
   real(wp),public,parameter :: amutokg = 1.660539040e-27_wp
   real(wp),public,parameter :: kgtoamu = 1.0_wp/amutokg
   real(wp),public,parameter :: metokg  = 9.10938356e-31_wp
   real(wp),public,parameter :: kgtome  = 1.0_wp/metokg
   real(wp),public,parameter :: amutoau = amutokg*kgtome
   real(wp),public,parameter :: autoamu = kgtoamu*metokg
!  femtosectons to atomic time units
   real(wp),public,parameter :: fstoau = 41.3413733365614_wp
   real(wp),public,parameter :: autofs = 1.0_wp/fstoau
!  Coulomb to atomic charge units (electrons)
   real(wp),public,parameter :: autoc = 1.6021766208e-19_wp
   real(wp),public,parameter :: ctoau = 1.0_wp/autoc
!  Debye to atomic units
   real(wp),public,parameter :: autod = autoc * lightspeed * autoaa**2 * fstoau * 1.0e+16_wp
   real(wp),public,parameter :: dtoau = 1.0_wp/autod
end module mctc_econv
