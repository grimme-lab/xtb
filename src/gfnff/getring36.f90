! This file is part of xtb.
!
! Copyright (C) 2019-2020 Stefan Grimme
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

subroutine getring36(n,at,nbin,a0_in,cout,irout)
      implicit none
      integer cout(10,20),irout(20)  ! output: atomlist, ringsize, # of rings in irout(20)
      integer at(n)
      integer n,nbin(20,n),a0,i,nb(20,n),a0_in
      integer    i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
      integer n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
      integer    a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      integer list(n),m,mm,nn,c(10),cdum(10,600),iring
      integer adum1(0:n),adum2(0:n),kk,j,idum(600),same(600),k
      real*8  w(n),av,sd
      logical chkrng

      if(n.le.2) return

      nn=nbin(20,a0_in)

      cdum=0
      kk=0
      do m=1,nn
!     if(nb(m,a0_in).eq.1) cycle
      nb=nbin
      do i=1,n
         if(nb(20,i).eq.1)nb(20,i)=0
      enddo

      do mm=1,nn
         w(mm)=dble(mm)
         list(mm)=mm
      enddo
      w(m)   =0.0d0
      call ssort(nn,w,list)
      do mm=1,nn
         nb(mm,a0_in)=nbin(list(mm),a0_in)
      enddo

      iring =0
      c     =0

      a0=a0_in
      n0=nb(20,a0)

      do i1=1,n0
         a1=nb(i1,a0)
         if(a1.eq.a0) cycle
         n1=nb(20,a1)
         do i2=1,n1
            a2=nb(i2,a1)
            if(a2.eq.a1) cycle
            n2=nb(20,a2)
            do i3=1,n2
               a3=nb(i3,a2)
               n3=nb(20,a3)
               if(a3.eq.a2) cycle
               c(1)=a1
               c(2)=a2
               c(3)=a3
               if(a3.eq.a0.and.chkrng(n,3,c))then
                iring=3
                kk=kk+1
                cdum(1:iring,kk)=c(1:iring)
                idum(kk)=iring
               endif
               do i4=1,n3
                  a4=nb(i4,a3)
                  n4=nb(20,a4)
                  if(a4.eq.a3) cycle
                  c(4)=a4
                  if(a4.eq.a0.and.chkrng(n,4,c))then
                   iring=4
                   kk=kk+1
                   cdum(1:iring,kk)=c(1:iring)
                   idum(kk)=iring
                  endif
                  do i5=1,n4
                     a5=nb(i5,a4)
                     n5=nb(20,a5)
                     if(a5.eq.a4) cycle
                     c(5)=a5
                     if(a5.eq.a0.and.chkrng(n,5,c))then
                      iring=5
                      kk=kk+1
                      cdum(1:iring,kk)=c(1:iring)
                      idum(kk)=iring
                     endif
                     do i6=1,n5
                        a6=nb(i6,a5)
                        n6=nb(20,a6)
                        if(a6.eq.a5) cycle
                        c(6)=a6
                        if(a6.eq.a0.and.chkrng(n,6,c))then
                         iring=6
                         kk=kk+1
                         cdum(1:iring,kk)=c(1:iring)
                         idum(kk)=iring
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo


 99   continue

      enddo

! compare
      same=0
      do i=1,kk
         do j=i+1,kk
            if(idum(i).ne.idum(j)) cycle ! different ring size
            if(same(j).eq.1      ) cycle ! already double
            adum1=0
            adum2=0
            do m=1,10
               i1=cdum(m,i)
               i2=cdum(m,j)
               adum1(i1)=1
               adum2(i2)=1
            enddo
            if(sum(abs(adum1-adum2)).ne.0) then
                  same(j)=0
            else
                  same(j)=1
            endif
         enddo
      enddo

      m=0
      do i=1,kk
         if(same(i).eq.0) then
            m=m+1
            if(m.gt.20) stop 'too many rings'
            irout(m)=idum(i)     ! number of atoms in ring m
            nn=idum(i)
            cout(1:nn,m)=cdum(1:nn,i)
            i2=0
            do k=1,nn            ! determine if its a hetereo
               i1=at(cdum(k,i))
               i2=i2+i1
            enddo
            av=dble(i2)/dble(nn)
            sd=0
            cout(m,19)=0
            do k=1,nn
               i1=at(cdum(k,i))
               sd=sd+(av-dble(i1))**2
            enddo
            if(sd.gt.1.d-6) cout(m,19)=idint(1000.*sqrt(sd)/dble(nn))
         endif
      enddo
      irout(20)=m  ! number of rings for this atom

      return
      end

      subroutine ssort(n,edum,ind)
      implicit none
      integer n,ii,k,j,m,i,sc1
      real*8 edum(n),pp
      integer ind(n)

      do 140   ii = 2, n
         i = ii - 1
         k = i
         pp= edum(i)
         do 120   j = ii, n
            if (edum(j) .gt. pp) go to 120
            k = j
            pp= edum(j)
  120    continue
         if (k .eq. i) go to 140
         edum(k) = edum(i)
         edum(i) = pp
         sc1=ind(i)
         ind(i)=ind(k)
         ind(k)=sc1
  140 continue

      end
