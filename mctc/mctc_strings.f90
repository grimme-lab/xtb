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

module mctc_strings
!use iso_fortran_env, only : kr4 => real32, kr8 => real64, &
!&                           ki4 => int32,  ki8 => int64
implicit none

! Real kinds
!> @brief single precision real
integer, parameter :: kr4 = selected_real_kind(6,37)
!> @brief double precision real
integer, parameter :: kr8 = selected_real_kind(15,307)

! Integer kinds
!> @brief single precision integer
integer, parameter :: ki4 = selected_int_kind(9)
!> @brief double precision integer
integer, parameter :: ki8 = selected_int_kind(18)

! Complex kinds
!> @brief single precision complex
integer, parameter :: kc4 = kr4
!> @brief double precision complex
integer, parameter :: kc8 = kr8

private :: kr4,kr8,ki4,ki8,kc4,kc8

private :: value_dr,value_sr,value_di,value_si
private :: write_dr,write_sr,write_di,write_si
private :: writeq_dr,writeq_sr,writeq_di,writeq_si

!> @brief Generic operator for converting a number string to a
!! number. Calling syntax is 'call value(numstring,number,ios)'
!! where 'numstring' is a number string and 'number' is a
!! real number or an integer (single or double precision).
interface value
   module procedure value_dr
   module procedure value_sr
   module procedure value_di
   module procedure value_si
end interface

!> @brief Generic  interface for writing a number to a string. The
!! number is left justified in the string. The calling syntax
!! is 'call writenum(number,string,format)' where 'number' is
!! a real number or an integer, 'string' is a character string
!! containing the result, and 'format' is the format desired,
!! e.g., 'e15.6' or 'i5'.
interface writenum
   module procedure write_dr
   module procedure write_sr
   module procedure write_di
   module procedure write_si
end interface

!> @brief Generic interface equating a name to a numerical value. The
!! calling syntax is 'call writeq(unit,name,value,format)' where
!! unit is the integer output unit number, 'name' is the variable
!! name, 'value' is the real or integer value of the variable,
!! and 'format' is the format of the value. The result written to
!! the output unit has the form <name> = <value>.
interface writeq
   module procedure writeq_dr
   module procedure writeq_sr
   module procedure writeq_di
   module procedure writeq_si
end interface


!**********************************************************************

contains

!**********************************************************************

pure function capitalize (str)
integer :: il,i
character(len=*),intent(in) :: str
character(len=len(str))     :: capitalize
character(len=26),parameter :: cap = &
&       'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
character(len=26),parameter :: low = &
&       'abcdefghijklmnopqrstuvwxyz'
capitalize = str
il = INDEX(low, str(1:1))
if (il.gt.0) capitalize(1:1) = cap(il:il)
do i = 2, len_trim(str)
   il = INDEX(cap, str(i:i))
   if (il.gt.0) capitalize(i:i) = low(il:il)
enddo
end function capitalize

!**********************************************************************

!> @brief Parses the string 'str' into arguments args(1), ..., args(nargs) based on
!! the delimiters contained in the string 'delims'. Preceding a delimiter in
!! 'str' by a backslash (\) makes this particular instance not a delimiter.
!! The integer output variable nargs contains the number of arguments found.
pure subroutine parse(str,delims,args,nargs)

character(len=*), intent(in) :: str
character(len=*), intent(in) :: delims
character(len=len_trim(str)) :: tmpstr
character(len=*), dimension(:), intent(out) :: args

integer, intent(out) :: nargs
integer :: na,i,lenstr,k

tmpstr=str
call compact(tmpstr)
na=size(args)
do i=1,na
  args(i)=' '
end do
nargs=0
lenstr=len_trim(tmpstr)
if(lenstr==0) return
k=0

do
   if(len_trim(tmpstr) == 0) exit
   nargs=nargs+1
   if (nargs > size(args)) exit
   call split(tmpstr,delims,args(nargs))
   call removebksl(args(nargs))
end do

end subroutine parse

!**********************************************************************

!> @brief Converts multiple spaces and tabs to single spaces; deletes control characters;
!! removes initial spaces.
pure subroutine compact(str)


character(len=*), intent(inout) :: str
character(len=1) :: ch
character(len=len_trim(str)) :: outstr

integer :: lenstr,isp,k,i,ich

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
isp=0
k=0

do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)

  select case(ich)

    case(9,32)     ! space or tab character
      if(isp==0) then
        k=k+1
        outstr(k:k)=' '
      end if
      isp=1

    case(33:)      ! not a space, quote, or control character
      k=k+1
      outstr(k:k)=ch
      isp=0

  end select

end do

str=adjustl(outstr)

end subroutine compact

!**********************************************************************

!> @brief Removes spaces, tabs, and control characters in string str
pure subroutine removesp(str)

character(len=*), intent(inout) :: str
character(len=1) :: ch
character(len=len_trim(str)) :: outstr

integer :: lenstr,k,i,ich

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0

do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)
  select case(ich)
    case(0:32)  ! space, tab, or control character
         cycle
    case(33:)
      k=k+1
      outstr(k:k)=ch
  end select
end do

str=adjustl(outstr)

end subroutine removesp

!**********************************************************************

!> @brief Converts number string to a double precision real number
pure subroutine value_dr(str,rnum,iostat)

character(len=*), intent(in) ::str
real(kr8), intent(out) ::rnum
integer, intent(out), optional :: iostat

integer :: ilen,ipos,ios

ilen=len_trim(str)
ipos=scan(str,'Ee')
if(.not.is_digit(str(ilen:ilen)) .and. ipos/=0) then
   if (present(iostat)) iostat=3
   return
end if
read(str,*,iostat=ios) rnum
if (present(iostat)) iostat = ios

end subroutine value_dr

!**********************************************************************

!> @brief Converts number string to a single precision real number
pure subroutine value_sr(str,rnum,iostat)

character(len=*), intent(in) ::str
real(kr4), intent(out) :: rnum
real(kr8) :: rnumd

integer, intent(out), optional :: iostat
integer :: ios

call value_dr(str,rnumd,ios)
if (present(iostat)) iostat = ios
if( abs(rnumd) > huge(rnum) ) then
   if (present(iostat)) iostat = 15
   return
end if
if( abs(rnumd) < tiny(rnum) ) rnum=0.0_kr4
rnum=rnumd

end subroutine value_sr

!**********************************************************************

!> @brief Converts number string to a double precision integer value
pure subroutine value_di(str,inum,iostat)


character(len=*), intent(in) ::str
integer(ki8), intent(out) :: inum
real(kr8) :: rnum

integer, intent(out), optional :: iostat
integer :: ios

call value_dr(str,rnum,ios)
if (present(iostat)) iostat = ios
if(abs(rnum)>huge(inum)) then
   if (present(iostat)) iostat = 15
   return
end if
inum=nint(rnum,ki8)

end subroutine value_di

!**********************************************************************

!> @brief Converts number string to a single precision integer value
pure subroutine value_si(str,inum,iostat)

character(len=*), intent(in) ::str
integer(ki4), intent(out) :: inum
real(kr8) :: rnum

integer, intent(out), optional :: iostat
integer :: ios

call value_dr(str,rnum,ios)
if (present(iostat)) iostat = ios
if(abs(rnum)>huge(inum)) then
   if (present(iostat)) iostat=15
   return
end if
inum=nint(rnum,ki4)

end subroutine value_si

!**********************************************************************

!> @brief Shifts characters in in the string 'str' n positions (positive values
!! denote a right shift and negative values denote a left shift). Characters
!! that are shifted off the end are lost. Positions opened up by the shift
!! are replaced by spaces.
pure subroutine shiftstr(str,n)

character(len=*), intent(inout) :: str
integer, intent(in) :: n

integer :: lenstr,nabs

lenstr=len(str)
nabs=iabs(n)
if(nabs>=lenstr) then
  str=repeat(' ',lenstr)
  return
end if
if(n<0) str=str(nabs+1:)//repeat(' ',nabs)  ! shift left
if(n>0) str=repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right
return

end subroutine shiftstr

!**********************************************************************

!> @brief Inserts the string 'strins' into the string 'str' at position 'loc'.
!! Characters in 'str' starting at position 'loc' are shifted right to
!! make room for the inserted string. Trailing spaces of 'strins' are
!! removed prior to insertion
pure subroutine insertstr(str,strins,loc)


character(len=*), intent(inout) :: str
character(len=*), intent(in) :: strins
character(len=len(str)) :: tempstr

integer, intent(in) :: loc
integer :: lenstrins

lenstrins=len_trim(strins)
tempstr=str(loc:)
call shiftstr(tempstr,lenstrins)
tempstr(1:lenstrins)=strins(1:lenstrins)
str(loc:)=tempstr
return

end subroutine insertstr

!**********************************************************************

!> @brief Deletes first occurrence of substring 'substr' from string 'str' and
!! shifts characters left to fill hole. Trailing spaces or blanks are
!! not considered part of 'substr'.
pure subroutine delsubstr(str,substr)

character(len=*), intent(inout) :: str
character(len=*), intent(in) :: substr

integer :: lensubstr,ipos

lensubstr=len_trim(substr)
ipos=index(str,substr)
if(ipos==0) return
if(ipos == 1) then
   str=str(lensubstr+1:)
else
   str=str(:ipos-1)//str(ipos+lensubstr:)
end if
return

end subroutine delsubstr

!**********************************************************************

!> @brief Deletes all occurrences of substring 'substr' from string 'str' and
!! shifts characters left to fill holes.
pure subroutine delall(str,substr)


character(len=*), intent(inout) :: str
character(len=*), intent(in) :: substr

integer :: lensubstr
integer :: ipos

lensubstr=len_trim(substr)
do
   ipos=index(str,substr)
   if(ipos == 0) exit
   if(ipos == 1) then
      str=str(lensubstr+1:)
   else
      str=str(:ipos-1)//str(ipos+lensubstr:)
   end if
end do
return

end subroutine delall

!**********************************************************************

!> @brief convert string to upper case
function uppercase(str) result(ucstr)


character (len=*):: str
character (len=len_trim(str)):: ucstr

integer :: ilen,ioffset,iquote,i,iav,iqc

ilen=len_trim(str)
ioffset=iachar('A')-iachar('a')
iquote=0
ucstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('a') .and. iav <= iachar('z')) then
    ucstr(i:i)=achar(iav+ioffset)
  else
    ucstr(i:i)=str(i:i)
  end if
end do
return

end function uppercase

!**********************************************************************

!> @brief convert string to lower case
function lowercase(str) result(lcstr)

character (len=*):: str
character (len=len_trim(str)):: lcstr

integer :: ilen,ioffset,iquote,i,iav,iqc

ilen=len_trim(str)
ioffset=iachar('A')-iachar('a')
iquote=0
lcstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('A') .and. iav <= iachar('Z')) then
    lcstr(i:i)=achar(iav-ioffset)
  else
    lcstr(i:i)=str(i:i)
  end if
end do
return

end function lowercase

!**********************************************************************

!> @brief Reads line from unit=nunitr, ignoring blank lines
!! and deleting comments beginning with an exclamation point(!)
subroutine readline(nunitr,line,ios)

character (len=*):: line

integer :: nunitr,ios,ipos

do
  read(nunitr,'(a)', iostat=ios) line      ! read input line
  if(ios /= 0) return
  line=adjustl(line)
  ipos=index(line,'!')
  if(ipos == 1) cycle
  if(ipos /= 0) line=line(:ipos-1)
  if(len_trim(line) /= 0) exit
end do
return

end subroutine readline

!**********************************************************************

!> @brief Sets imatch to the position in string of the delimiter matching the delimiter
!! in position ipos. Allowable delimiters are (), [], {}, <>.
pure subroutine match(str,ipos,imatch,status)

character(len=*), intent(in) :: str
character :: delim1,delim2,ch

integer, intent(out), optional :: status
integer :: stat

integer, intent(in) :: ipos
integer, intent(out) :: imatch
integer :: lenstr,idelim2,istart,inc,iend,isum,i

lenstr=len_trim(str)
delim1=str(ipos:ipos)
select case(delim1)
   case('(')
      idelim2=iachar(delim1)+1
      istart=ipos+1
      iend=lenstr
      inc=1
   case(')')
      idelim2=iachar(delim1)-1
      istart=ipos-1
      iend=1
      inc=-1
   case('[','{','<')
      idelim2=iachar(delim1)+2
      istart=ipos+1
      iend=lenstr
      inc=1
   case(']','}','>')
      idelim2=iachar(delim1)-2
      istart=ipos-1
      iend=1
      inc=-1
   case default
      stat = 1
      if (present(status)) status = stat
      !write(*,*) delim1,' is not a valid delimiter'
      return
end select
if(istart < 1 .or. istart > lenstr) then
   stat = 2
   if (present(status)) status = stat
   !write(*,*) delim1,' has no matching delimiter'
   return
end if
delim2=achar(idelim2) ! matching delimiter

isum=1
do i=istart,iend,inc
   ch=str(i:i)
   if(ch /= delim1 .and. ch /= delim2) cycle
   if(ch == delim1) isum=isum+1
   if(ch == delim2) isum=isum-1
   if(isum == 0) exit
end do
if(isum /= 0) then
   stat = 3
   if (present(status)) status = stat
   !write(*,*) delim1,' has no matching delimiter'
   return
end if
imatch=i
if (present(status)) status = 0

return

end subroutine match

!**********************************************************************

!> @brief Writes double precision real number rnum to string str using format fmt
pure subroutine write_dr(rnum,str,fmt)

real(kr8), intent(in) :: rnum
character(len=*), intent(out) :: str
character(len=*), intent(in), optional :: fmt
character(len=80) :: formt

if (present(fmt)) then
   formt='('//trim(fmt)//')'
else
   formt='(g0)'
endif
write(str,formt) rnum
str=adjustl(str)

end subroutine write_dr

!***********************************************************************

!> @brief Writes single precision real number rnum to string str using format fmt
pure subroutine write_sr(rnum,str,fmt)

real(kr4), intent(in) :: rnum
character(len=*), intent(out) :: str
character(len=*), intent(in), optional :: fmt
character(len=80) :: formt

if (present(fmt)) then
   formt='('//trim(fmt)//')'
else
   formt='(g0)'
endif
write(str,formt) rnum
str=adjustl(str)

end subroutine write_sr

!***********************************************************************

!> @brief Writes double precision integer inum to string str using format fmt
pure subroutine write_di(inum,str,fmt)

integer(ki8), intent(in) :: inum
character(len=*), intent(out) :: str
character(len=*), intent(in), optional :: fmt
character(len=80) :: formt

if (present(fmt)) then
   formt='('//trim(fmt)//')'
else
   formt='(g0)'
endif
write(str,formt) inum
str=adjustl(str)

end subroutine write_di

!***********************************************************************

!> @brief Writes single precision integer inum to string str using format fmt
pure subroutine write_si(inum,str,fmt)

integer(ki4), intent(in) :: inum
character(len=*), intent(out) :: str
character(len=*), intent(in), optional :: fmt
character(len=80) :: formt

if (present(fmt)) then
   formt='('//trim(fmt)//')'
else
   formt='(g0)'
endif
write(str,formt) inum
str=adjustl(str)

end subroutine write_si

!***********************************************************************


!> @brief Deletes nonsignificant trailing zeroes from number string str. If number
!! string ends in a decimal point, one trailing zero is added.
pure subroutine trimzero(str)

character(len=*), intent(inout) :: str
character :: ch
character(len=10) :: exp

integer :: ipos,i,lstr

ipos=scan(str,'eE')
if(ipos>0) then
   exp=str(ipos:)
   str=str(1:ipos-1)
endif
lstr=len_trim(str)
do i=lstr,1,-1
   ch=str(i:i)
   if(ch=='0') cycle
   if(ch=='.') then
      str=str(1:i)//'0'
      if(ipos>0) str=trim(str)//trim(exp)
      exit
   endif
   str=str(1:i)
   exit
end do
if(ipos>0) str=trim(str)//trim(exp)

end subroutine trimzero

!**********************************************************************

!> @brief Writes a string of the form <name> = value to unit
subroutine writeq_dr(unit,namestr,value,fmt)

real(kr8) :: value
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr

call writenum(value,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_dr

!**********************************************************************

!> @brief Writes a string of the form <name> = value to unit
subroutine writeq_sr(unit,namestr,value,fmt)


real(kr4) :: value
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr

call writenum(value,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_sr

!**********************************************************************

!> @brief Writes a string of the form <name> = ivalue to unit
subroutine writeq_di(unit,namestr,ivalue,fmt)


integer(ki8) :: ivalue
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr
call writenum(ivalue,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_di

!**********************************************************************

!> @brief Writes a string of the form <name> = ivalue to unit
subroutine writeq_si(unit,namestr,ivalue,fmt)


integer(ki4) :: ivalue
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr
call writenum(ivalue,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_si

!**********************************************************************

!> @brief Returns .true. if ch is a letter and .false. otherwise
function is_letter(ch) result(res)


character :: ch
logical :: res

select case(ch)
case('A':'Z','a':'z')
  res=.true.
case default
  res=.false.
end select
return

end function is_letter

!**********************************************************************

!> @brief Returns .true. if ch is a digit (0,1,...,9) and .false. otherwise
pure elemental function is_digit(ch) result(res)

character, intent(in) :: ch
logical :: res

select case(ch)
case('0':'9')
  res=.true.
case default
  res=.false.
end select
return

end function is_digit

!**********************************************************************

!> @brief Routine finds the first instance of a character from 'delims' in the
!! the string 'str'. The characters before the found delimiter are
!! output in 'before'. The characters after the found delimiter are
!! output in 'str'. The optional output character 'sep' contains the
!! found delimiter. A delimiter in 'str' is treated like an ordinary
!! character if it is preceded by a backslash (\). If the backslash
!! character is desired in 'str', then precede it with another backslash.
pure subroutine split(str,delims,before,sep)

character(len=*), intent(inout) :: str
character(len=*), intent(in) :: delims
character(len=*), intent(out) :: before
character, intent(out), optional :: sep
logical :: pres
character :: ch,cha

integer :: lenstr,k,ibsl,i,ipos,iposa

pres=present(sep)
str=adjustl(str)
call compact(str)
lenstr=len_trim(str)
if(lenstr == 0) return        ! string str is empty
k=0
ibsl=0                        ! backslash initially inactive
before=' '
do i=1,lenstr
   ch=str(i:i)
   if(ibsl == 1) then          ! backslash active
      k=k+1
      before(k:k)=ch
      ibsl=0
      cycle
   end if
   if(ch == '\') then          ! backslash with backslash inactive
      k=k+1
      before(k:k)=ch
      ibsl=1
      cycle
   end if
   ipos=index(delims,ch)
   if(ipos == 0) then          ! character is not a delimiter
      k=k+1
      before(k:k)=ch
      cycle
   end if
   if(ch /= ' ') then          ! character is a delimiter that is not a space
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
   cha=str(i+1:i+1)            ! character is a space delimiter
   iposa=index(delims,cha)
   if(iposa > 0) then          ! next character is a delimiter
      str=str(i+2:)
      if(pres) sep=cha
      exit
   else
      str=str(i+1:)
      if(pres) sep=ch
      exit
   end if
end do
if(i >= lenstr) str=''
str=adjustl(str)              ! remove initial spaces
return

end subroutine split

!**********************************************************************

!> @brief Removes backslash (\) characters. Double backslashes (\\) are replaced
!! by a single backslash.
pure subroutine removebksl(str)

character(len=*), intent(inout) :: str
character(len=1) :: ch
character(len=len_trim(str)) :: outstr

integer :: lenstr,k,ibsl,i

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0
ibsl=0                        ! backslash initially inactive

do i=1,lenstr
  ch=str(i:i)
  if(ibsl == 1) then          ! backslash active
   k=k+1
   outstr(k:k)=ch
   ibsl=0
   cycle
  end if
  if(ch == '\') then          ! backslash with backslash inactive
   ibsl=1
   cycle
  end if
  k=k+1
  outstr(k:k)=ch              ! non-backslash with backslash inactive
end do

str=adjustl(outstr)

end subroutine removebksl

!**********************************************************************

end module mctc_strings

