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

subroutine eval_define(verbose)
   use iso_fortran_env, only : id => output_unit
   use mctc_global, only : msgid, name
   implicit none
   logical,intent(in) :: verbose
   integer :: save_i
   if (verbose) then
      name = 'define'
      call definebanner
   endif
   save_i = msgid
   call raise('F','Your input has following faults:',1)
   call terminate(save_i)
end subroutine eval_define

