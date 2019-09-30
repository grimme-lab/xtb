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

pure subroutine getsymnum(pgroup,lin,symnum)
   use iso_fortran_env, only : wp => real64
   implicit none
   real(wp),intent(out) :: symnum
   character(len=*),intent(in) :: pgroup
   logical,intent(in) :: lin
   symnum=1
   if(index(pgroup,'ci').ne.0)symnum=1
   if(index(pgroup,'cs').ne.0)symnum=1
   if(index(pgroup,'c2').ne.0)symnum=2
   if(index(pgroup,'c3').ne.0)symnum=3
   if(index(pgroup,'s6').ne.0)symnum=3
   if(index(pgroup,'c4').ne.0)symnum=4
   if(index(pgroup,'c5').ne.0)symnum=5
   if(index(pgroup,'c6').ne.0)symnum=6
   if(index(pgroup,'c7').ne.0)symnum=7
   if(index(pgroup,'c8').ne.0)symnum=8
   if(index(pgroup,'c9').ne.0)symnum=9
   if(index(pgroup,'c10').ne.0)symnum=10
   if(index(pgroup,'c11').ne.0)symnum=11
   if(index(pgroup,'d2').ne.0)symnum=4
   if(index(pgroup,'d3').ne.0)symnum=6
   if(index(pgroup,'d4').ne.0)symnum=8
   if(index(pgroup,'d5').ne.0)symnum=10
   if(index(pgroup,'d6').ne.0)symnum=12
   if(index(pgroup,'d7').ne.0)symnum=14
   if(index(pgroup,'d8').ne.0)symnum=16
   if(index(pgroup,'d9').ne.0)symnum=18
   if(index(pgroup,'d10').ne.0)symnum=20
   if(index(pgroup,'td').ne.0)symnum=12
   if(index(pgroup,'oh').ne.0)symnum=24
   if(index(pgroup,'ih').ne.0)symnum=60
   if(index(pgroup,'c').ne.0.and.lin)symnum=1
   if(index(pgroup,'d').ne.0.and.lin)symnum=2
   if(index(pgroup,'Ci').ne.0)symnum=1
   if(index(pgroup,'Cs').ne.0)symnum=1
   if(index(pgroup,'C2').ne.0)symnum=2
   if(index(pgroup,'C3').ne.0)symnum=3
   if(index(pgroup,'S6').ne.0)symnum=3
   if(index(pgroup,'C4').ne.0)symnum=4
   if(index(pgroup,'C5').ne.0)symnum=5
   if(index(pgroup,'C6').ne.0)symnum=6
   if(index(pgroup,'C7').ne.0)symnum=7
   if(index(pgroup,'C8').ne.0)symnum=8
   if(index(pgroup,'C9').ne.0)symnum=9
   if(index(pgroup,'C10').ne.0)symnum=10
   if(index(pgroup,'C11').ne.0)symnum=11
   if(index(pgroup,'D2').ne.0)symnum=4
   if(index(pgroup,'D3').ne.0)symnum=6
   if(index(pgroup,'D4').ne.0)symnum=8
   if(index(pgroup,'D5').ne.0)symnum=10
   if(index(pgroup,'D6').ne.0)symnum=12
   if(index(pgroup,'D7').ne.0)symnum=14
   if(index(pgroup,'D8').ne.0)symnum=16
   if(index(pgroup,'D9').ne.0)symnum=18
   if(index(pgroup,'D10').ne.0)symnum=20
   if(index(pgroup,'Td').ne.0)symnum=12
   if(index(pgroup,'Oh').ne.0)symnum=24
   if(index(pgroup,'Ih').ne.0)symnum=60
   if(index(pgroup,'C').ne.0.and.lin)symnum=1
   if(index(pgroup,'D').ne.0.and.lin)symnum=2
end subroutine getsymnum
