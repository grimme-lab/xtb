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

module mctc_param
   use iso_fortran_env, wp => real64
   use mctcpar_pse
   use mctcpar_electronegativities
   use mctcpar_atomic_masses
   use mctcpar_r4r2
   use mctcpar_covalent_radii
   use mctcpar_chemical_hardnesses

   implicit none
   public
end module mctc_param
