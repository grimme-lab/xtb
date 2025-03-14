! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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


subroutine elem(key1, nat)
   implicit integer (i-n)
   implicit double precision (a-h,o-z)
   character(len=*) :: key1
   character(len=2) :: e
   character(len=2),parameter :: elemnt(118) = (/ &
      & 'h ','he', &
      & 'li','be','b ','c ','n ','o ','f ','ne', &
      & 'na','mg','al','si','p ','s ','cl','ar', &
      & 'k ','ca', &
      & 'sc','ti','v ','cr','mn','fe','co','ni','cu','zn', &
      &           'ga','ge','as','se','br','kr', &
      & 'rb','sr', &
      & 'y ','zr','nb','mo','tc','ru','rh','pd','ag','cd', &
      &           'in','sn','sb','te','i ','xe', &
      & 'cs','ba','la', &
      & 'ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb', &
      & 'lu','hf','ta','w ','re','os','ir','pt','au','hg', &
      &           'tl','pb','bi','po','at','rn', &
      & 'fr','ra','ac', &
      & 'th','pa','u ','np','pu','am','cm','bk','cf','es','fm','md','no', &
      & 'lr','rf','db','sg','bh','hs','mt','ds','rg','cn', &
      &           'nh','fl','mc','lv','ts','og' /)

   nat=0
   e='  '
   do i=1,len(key1)
      if(key1(i:i).ne.' ') L=i
   enddo
   k=1
   do j=1,l
      if (k.gt.2)exit
      n=ichar(key1(J:J))
      if(len_trim(e).ge.1 .and. n.eq.ichar(' '))exit!break if space after elem. symbol
      if(len_trim(e).ge.1 .and. n.eq.9)exit !break if tab after elem. symbol
      if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
         e(k:k)=char(n+ichar('a')-ichar('A'))
         k=k+1
      endif
      if(n.ge.ichar('a') .and. n.le.ichar('z') )then
         e(k:k)=key1(j:j)
         k=k+1
      endif
   enddo

   do i=1,118
      if(e.eq.elemnt(i))then
         nat=i
         return
      endif
   enddo
end subroutine elem



