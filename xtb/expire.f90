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

subroutine expire
   use iso_fortran_env
   implicit none
   intrinsic date_and_time

   integer :: i,upper(3),lower(3),val(8),edate,cdate,adate
   parameter (upper=(/2019, 4,30/)) ! expire date
   parameter (edate=365*upper(1) + 30.5*upper(2) + upper(3))
   parameter (lower=(/2018,11, 1/)) ! build date
   parameter (cdate=365*lower(1) + 30.5*lower(2) + lower(3))

   call date_and_time(values=val)
   adate=365*val(1) + 30.5*val(2) + val(3) ! actual date

   if ((adate.gt.cdate).and.(adate.lt.edate)) then
      write(output_unit,'(3x,i0,1x,a)') edate-adate, &
      &   'days left until this program expires'
   else
      write(output_unit,'(a)') 'code expired'
      call terminate(0)
   endif

end subroutine expire
