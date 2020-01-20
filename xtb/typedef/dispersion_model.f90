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

module tbdef_dispersion_model
   use iso_fortran_env, only: wp => real64
   implicit none
   public :: tb_dispersion_model
   private

   integer, parameter :: max_elem = 118

   type :: tb_dispersion_model
      real(wp) :: g_a = 0.0_wp
      real(wp) :: g_c = 0.0_wp
      integer, dimension(max_elem) :: atoms = 0
      integer, dimension(max_elem) :: nref = 0
      integer, dimension(7,max_elem) :: ncount = 0
      real(wp),dimension(7,max_elem) :: cn = 0.0_wp
      real(wp),dimension(7,max_elem) :: q = 0.0_wp
      real(wp),dimension(23,7,max_elem) :: alpha = 0.0_wp
      real(wp),dimension(7,7,max_elem,max_elem) :: c6 = 0.0_wp
   end type tb_dispersion_model

end module tbdef_dispersion_model
