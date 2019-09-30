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

!cuts the at blanks and tabstops and returns all floats and strings in order of occurence
subroutine cutline(line,floats,strings)
   use iso_fortran_env, only : wp => real64
   implicit none
   real(wp) floats(*),num
   character(len=128) line,str,stmp
   character(len=80) strings(3)
   character(len=1) digit
   integer i,ty,cs,cf

   stmp=''
   cs=1
   cf=1
   strings(:)=''
   do i=1,len(trim(line))
      digit=line(i:i)
      if(digit.ne.' '.and.digit.ne.char(9)) then  !should exclude tabstops and blanks, 9 is ascii code for tab
         stmp=trim(stmp)//trim(digit)
      elseif(stmp.ne.'')then
         call checktype(stmp,num,str,ty)      !get type of string, 0=number, 1=character
         if(ty.eq.0) then
            floats(cf)=num
            cf=cf+1
         elseif(ty.eq.1) then
            strings(cs)=trim(str)
            cs=cs+1
         else
            write(*,*)'Problem in checktype, must abort'
            exit
         endif
         stmp=''
      endif
      if(i.eq.len(trim(line))) then  !special case: end of line
         call checktype(stmp,num,str,ty)
         if(ty.eq.0) then
            floats(cf)=num
            cf=cf+1
         elseif(ty.eq.1) then
            strings(cs)=trim(str)
            cs=cs+1
         else
            write(*,*)'Problem in checktype, must abort'
            exit
         endif
         stmp=''
      endif
   enddo
end subroutine cutline


!this checks the type of the string and returns it cast to real or as string.
subroutine checktype(field,num,str,ty)
   use iso_fortran_env, only : wp => real64
   implicit none
   character(len=*) field,str
   real(wp) num
   integer i,e,ty
   logical is_num

   ty=99
   str=''
   is_num=.false.
   read(field,'(F10.5)',IOSTAT=e)num !cast string on real and get error code; 0 means success.
   if(e.eq.0)is_num=.true.
   if(is_num)then
      if(index(field,'.').ne.0) then  !check for integer/real
         read(field,'(F30.16)')num
         ty=0
      else                       !if integer, add .0 to string; otherwise cast to real does not work
         str=trim(field)//'.0'
         read(str,'(F30.16)')num
         str=''
         ty=0
      endif
   else
      str=field
      ty=1
   endif
end subroutine checktype


