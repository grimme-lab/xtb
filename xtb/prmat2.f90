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

subroutine prgemat(mat,d1,d2,name,inunit,instep)
   use iso_fortran_env, only : wp => real64, output_unit
   implicit none
   integer, intent(in) :: d1
   integer, intent(in) :: d2
   real(wp),intent(in) :: mat(d1,d2)
   character(len=*),intent(in),optional :: name
   integer, intent(in),optional :: inunit
   integer, intent(in),optional :: instep
   integer :: i,j,k,l,step,unit
   if (present(inunit)) then
      unit = inunit
   else
      unit = output_unit
   endif
   if (present(instep)) then
      step = instep
   else
      step = 6
   endif
   if(present(name)) write(unit,'(/,''matrix printed:'',x,a)') name
   do i = 1, d2, step
      l = min(i+step-1,d2)
      write(unit,'(/,6x)',advance='no')
      do k = i, l
      write(unit,'(6x,i7,3x)',advance='no') k
      enddo
      write(unit,'(a)')
      do j = 1, d1
         write(unit,'(i6)',advance='no') j
         do k = i, l
         write(unit,'(x,e15.8)',advance='no') mat(j,k)
         enddo
         write(unit,'(a)')
      enddo
   enddo
   return
end subroutine prgemat


subroutine prsymat(mat,d1,name,inunit,instep)
   use iso_fortran_env, only : wp => real64, output_unit
   implicit none
   integer, intent(in) :: d1
   real(wp),intent(in) :: mat(d1*(d1+1))
   character(len=*),intent(in),optional :: name
   integer, intent(in),optional :: inunit
   integer, intent(in),optional :: instep
   integer :: i,j,k,l,step,unit
   integer,external :: lin
   if (present(inunit)) then
      unit = inunit
   else
      unit = output_unit
   endif
   if (present(instep)) then
      step = instep
   else
      step = 6
   endif
   if(present(name)) write(unit,'(/,''matrix printed:'',x,a)') name
   do i = 1, d1, step
      l = min(i+step-1,d1)
      write(unit,'(/,6x)',advance='no')
      do k = i, l
      write(unit,'(6x,i7,3x)',advance='no') k
      enddo
      write(unit,'(a)')
      do j = i, d1
         l = min(i+(step-1),j)
         write(unit,'(i6)',advance='no') j
         do k = i, l
         write(unit,'(x,e15.8)',advance='no') mat(lin(j,k))
         enddo
         write(unit,'(a)')
      enddo
   enddo
   return
end subroutine prsymat
