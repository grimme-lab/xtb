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


subroutine readl(a1,x,n)
   use xtb_mctc_accuracy, only : wp
   implicit integer (i-n)
   implicit real(wp) (a-h,o-z)
   character(*) a1
   dimension x(*)
   i=0
   is=1
   10  i=i+1
   x(i)=readaa(a1,is,ib,ie)
   if(ib.gt.0 .and. ie.gt.0) then
      is=ie
      goto 10
   endif
   n=i-1
   return
end subroutine readl


function readaa(a,istart,iend,iend2)
   use xtb_mctc_accuracy, only : wp
   implicit integer (i-n)
   implicit real(wp) (a-h,o-z)
   real(wp) readaa
   character(*) a

   NINE=ICHAR('9')
   IZERO=ICHAR('0')
   MINUS=ICHAR('-')
   IDOT=ICHAR('.')
   ND=ICHAR('D')
   NE=ICHAR('E')
   IBL=ICHAR(' ')

   iend=0
   iend2=0
   idig=0
   c1=0
   c2=0
   one=1.d0
   x = 1.d0
   nl=len(a)
   do j=istart,nl-1
      n=ichar(a(j:j))
      m=ichar(a(j+1:j+1))
      if(n.le.nine.and.n.ge.izero .or.n.eq.idot)goto 20
      if(n.eq.minus.and.(m.le.nine.and.m.ge.izero &
         & .or. m.eq.idot)) goto 20
   enddo
   readaa=0.d0
   return
   20  continue
   iend=j
   do i=j,nl
      n=ichar(a(i:i))
      if(n.le.nine.and.n.ge.izero) then
         idig=idig+1
         if (idig.gt.10) goto 60
         c1=c1*10+n-izero
      elseif(n.eq.minus.and.i.eq.j) then
         one=-1.d0
      elseif(n.eq.idot) then
         goto 40
      else
         goto 60
      endif
   enddo
   40  continue
   idig=0
   do ii=i+1,nl
      n=ichar(a(ii:ii))
      if(n.le.nine.and.n.ge.izero) then
         idig=idig+1
         if (idig.gt.10) goto 60
         c2=c2*10+n-izero
         x = x /10
      elseif(n.eq.minus.and.ii.eq.i) then
         x=-x
      else
         goto 60
      endif
   enddo
   !
   ! put the pieces together
   !
   60  continue
   readaa= one * ( c1 + c2 * x)
   do j=iend,nl
      n=ichar(a(j:j))
      iend2=j
      if(n.eq.ibl)return
      if(n.eq.nd .or. n.eq.ne) goto 57
   enddo
   return

   57  c1=0.0d0
   one=1.0d0
   do i=j+1,nl
      n=ichar(a(i:i))
      iend2=i
      if(n.eq.ibl)goto 70
      if(n.le.nine.and.n.ge.izero) c1=c1*10.0d0+n-izero
      if(n.eq.minus)one=-1.0d0
   enddo
   61  continue
   70  readaa=readaa*10**(one*c1)
   return
end function readaa



