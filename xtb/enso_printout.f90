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

subroutine enso_printout(iunit,res)
   use iso_fortran_env, only : wp => real64
   use tbdef_data

   implicit none

   integer,intent(in) :: iunit
   type(freq_results),intent(in) :: res

   write(iunit,'("{")')
   write(iunit,'(3x,"""",a,""":",f21.12,",")') "temperature",res%temp
   write(iunit,'(3x,"""",a,""": """,a,""",")')      "unit","hartree"
   write(iunit,'(3x,"""",a,""":",f21.12,",")') "energy",res%etot
   write(iunit,'(3x,"""",a,""":",f21.12,",")') "free energy",res%gtot+res%etot
   write(iunit,'(3x,"""",a,""":",f21.12,",")') "G(T)",res%gtot
   write(iunit,'(3x,"""",a,""":",i0,",")') "number of imags",res%nimag
   write(iunit,'(3x,"""",a,""": """,a,""",")')      "point group",res%pg
   write(iunit,'(3x,"""",a,""":",f21.12)')     "gradient norm",res%gnorm
   write(iunit,'("}")')

end subroutine enso_printout
