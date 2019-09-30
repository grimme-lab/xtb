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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c analyze a pot curve with n points for stationary points
c xx are the grid points, yy the energies and yy* the 
c derivatives
c nstat is the number of stationary points found
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine curveanal(n,xx,yy,yy1,yy2,list,ityp,nstat,eps)
      use iso_fortran_env, only : wp => real64
      implicit none 
      integer n
      real(wp) xx(n),yy(n),yy1(n),yy2(n),eps
      integer list(*),ityp(*),nstat

      integer i,j,k,kk,mem,it1,it2,icall
      real(wp) xmin,av,dx,dxmax
      character*15 typ(3)
  
      write(*,*) 'analysis of potential (stationary points)'
      typ( 1)='##saddle       ' 
      typ( 2)='##probably min.'
      typ( 3)='##minimum      '

c     do i=1,n
c        write(142,*) xx(i),yy(i)
c     enddo
c     write(142,*) 
c     do i=1,n
c        write(142,*) xx(i),yy1(i)
c     enddo
c     write(142,*) 
 
      dxmax=xx(n)-xx(1)
     
      av=0
      k =0
      do i=1,n
         if(abs(yy1(i)).gt.1.d-6)then
            av=av+abs(yy1(i))
            k=k+1
         endif
      enddo
      av=av/dble(k)

      eps=0.05*av
      k=0
      icall=0

   1  continue
      write(*,*) 'epsilon (for zero check in Eh/Bohr) : ',eps
      if(eps.ne.eps) stop
      i=1
  10  if(abs(yy1(i)).lt.eps.and.abs(yy2(i)).gt.1.d-9)then
            k=k+1
            xmin=1.d+42
            do j=i,i+50
               if(abs(yy1(j)).lt.xmin)then
                  mem=j
                  xmin=abs(yy1(j))
               endif
            enddo
            list(k)=mem
            ityp(k)=2+int(sign(1.0d0,yy2(list(k))))
            i=j+50  ! skip 
            if(k.gt.1)then   ! check if the the same point is close
               it2=2+int(sign(1.0d0,yy2(list(k-1))))
               dx=abs(xx(list(k))-xx(list(k-1)))
               if(ityp(k).eq.it2.and.dx/dxmax.lt.0.1)
     .         k=k-1 ! invalid because its close
               if(k.gt.100) stop 'too many stationary points'
            endif
      endif
      i=i+1
      if(i.lt.n) goto 10

      if(k.eq.0.and.icall.lt.10) then
         eps=eps*5.
         icall=icall+1
         goto 1
      endif

c     check the borders
      kk=k
      i =1
c     do i=1,k
         if(ityp(i).eq.1)then
            if(yy(list(i)).gt.yy(2))then
               kk=kk+1
               ityp(kk)=2
               do j=1,n
                  if(abs(yy2(j)).gt.eps)then
                     list(kk)=j
                     goto 20
                  endif
               enddo
            endif
         endif 
c     enddo
 20   continue
c     do i=k,1,-1
      i=k
         if(ityp(i).eq.1)then
            if(yy(list(i)).gt.yy(n-1))then
               kk=kk+1
               ityp(kk)=2
               do j=n,1,-1
                  if(abs(yy2(j)).gt.eps)then
                     list(kk)=j
                     goto 30
                  endif
               enddo
            endif
         endif 
c     enddo
 30   continue

      do i=1,kk
         write(*,'(i3,2x,a15,'' at '',f8.3,5x,''Erel /kcal '',f8.2)') 
     .   i,typ(ityp(i)),
     .   xx(list(i)),yy(list(i))*627.509
      enddo

      nstat = kk

      end

