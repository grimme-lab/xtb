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

module d3param
   use iso_fortran_env, only : wp => real64
   implicit none
   private :: wp
   public

!  d3 common variable block                               
   real(wp) :: rcov
   real(wp) :: r2r4
   real(wp) :: c6ab
   integer  :: mxc
   dimension :: rcov(94),r2r4(94),c6ab(94,94,5,5,3),mxc(94)


end module d3param
