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
! You should have received a copy of the GNU Lesser General Public Licen
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

subroutine geosum(n,ic,xyzin)
   use xtb_mctc_convert, only : autoaa
   use xtb_mctc_symbols, only : toSymbol
   use xtb_param_atomicrad, only : atomicRad
   use xtb_splitparam
   implicit none
   integer n,ic(n)
   real*8 xyzin(3,n)

   real*8 ,allocatable :: x(:),y(:),z(:),rmin(:),rmax(:),dev(:)
   real*8 ,allocatable :: rlist(:,:),av(:)
   integer,allocatable :: bondlist(:,:),nbond(:),nty(:)

   real*8 brad,r,xxx,yyy,zzz,avv
   real*8 xyz(3,3),ang,xyz4(3,4),ra(3),rb(3)
   integer i, j, nn, k, imin, jmin, iat1, iat2
   integer nm,ij,m,iat,jat, nsel
   integer nt,a,b,ia,ib,ja,jb,ii,jj,lin,linij
   character*31  lab,lab2,angsum1 (n*150)
   character*30  lab4,angsum2(n*150)
   character*50  lab3
   character*21  lab21

   ! first nsel atoms are considered in angle/dihedreal xtb_printout
   nsel=min(10,n)

   nn=100*100
   allocate(x(n),y(n),z(n),rmin(nn),rmax(nn),dev(nn),av(nn),&
      &         rlist(n,150),bondlist(n,150),nbond(n),nty(nn))

   x(1:n)=xyzin(1,1:n)
   y(1:n)=xyzin(2,1:n)
   z(1:n)=xyzin(3,1:n)
   rmin=10000.
   k   =0
   do i=1,n
      nbond(i)=0
      do j=1,n
         if(i.eq.j)cycle
         k=k+1
         r=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)&
            &        *0.529177260d0
         ! the max bond dist in angst
         brad=(atomicRad(ic(i))+atomicRad(ic(j)))*autoaa
         brad=brad+0.30*brad
         if(r.le.brad) then
            nbond(i)=nbond(i)+1
            bondlist(i,nbond(i))=j
            rlist   (i,nbond(i))=r
         endif
      enddo
   enddo

   write(*,*)
   write(*,*)'Bond Distances (Angstroems)'
   write(*,*)'---------------------------'
   nt=0
   nty=0
   av=0
   dev=0
   rmin=1000
   rmax=-1
   do i=1,n
      do a=1,nbond(i)
         nt=nt+1
         write(lab21,'(a2,i3,''-'',a2,i3)')&
            &      toSymbol(ic(i)),i,toSymbol(ic(bondlist(i,a))),bondlist(i,a)
         call rmblank(lab21,lab4)
         write(angsum1(nt),'(a,''='',f6.4)')&
            &      trim(lab4),rlist(i,a)
         iat1=ic(i)
         iat2=ic(bondlist(i,a))
         linij=lin(iat1,iat2)
         av(linij)=av(linij)+rlist(i,a)
         if(rlist(i,a).lt.rmin(linij))rmin(linij)=rlist(i,a)
         if(rlist(i,a).gt.rmax(linij))rmax(linij)=rlist(i,a)
         nty(linij)=nty(linij)+1
      enddo
   enddo

   do i=1,n
      do a=1,nbond(i)
         iat1=ic(i)
         iat2=ic(bondlist(i,a))
         linij=lin(iat1,iat2)
         avv=av(linij)/float(nty(linij))
         dev(linij)=dev(linij)+(avv-rlist(i,a))**2
      enddo
   enddo
   write(*,'(6a21)')(angsum1(i),i=1,nt)
   do i=1,100
      do j=1,i
         if(nty(lin(i,j)).gt.0)then
            write(*,'(2A3,'' Rav='',f6.4,'' sigma='',f6.4,&
               &      ''  Rmin='',f6.4,''  Rmax='',f6.4,i6)')&
               &      toSymbol(i),toSymbol(j),&
               &      av(lin(i,j))/nty(lin(i,j)),&
               &      sqrt(dev(lin(i,j))/float(nty(lin(i,j)))),&
               &      rmin(lin(i,j)),rmax(lin(i,j)),int(0.5*nty(lin(i,j))+0.0001)
         endif
      enddo
   enddo

   write(*,*)
   write(*,*)'selected bond angles (degree)'
   write(*,*)'--------------------'
   nt=0
   do i=1,nsel
      m=nbond(i)
      if(m.gt.1)then
         xyz(1,2)=x(i)
         xyz(2,2)=y(i)
         xyz(3,2)=z(i)
         do a=1,m
            ia=bondlist(i,a)
            do b=1,a-1
               ib=bondlist(i,b)
               nt=nt+1
               write(lab,'(a2,i5,''-'',a2,i5,''-'',a2,i5)')&
                  &toSymbol(ic(ia)),ia,toSymbol(ic(i)),i,toSymbol(ic(ib)),ib
               call rmblank(lab,lab2)
               xyz(1,1)=x(ia)
               xyz(2,1)=y(ia)
               xyz(3,1)=z(ia)
               xyz(1,3)=x(ib)
               xyz(2,3)=y(ib)
               xyz(3,3)=z(ib)
               call xbangle(xyz,ang)
               write(angsum1(nt),'(a,''='',f6.2)')trim(lab2),ang
            enddo
         enddo
      endif
   enddo
   write(*,'(4a)')(angsum1(i),i=1,nt)

   write(*,*)
   write(*,*)'selected dihedral angles (degree)'
   write(*,*)'---------------------------------'
   nt=0
   do i=1,nsel
      xyz4(1,2)=x(i)
      xyz4(2,2)=y(i)
      xyz4(3,2)=z(i)
      do j=1,i
         if(i.eq.j)cycle
         xyz4(1,3)=x(j)
         xyz4(2,3)=y(j)
         xyz4(3,3)=z(j)
         k=k+1
         r=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)&
            &        *0.529177260d0
         ! the max bond dist in angst
         brad=(atomicRad(ic(i))+atomicRad(ic(j)))*autoaa
         brad=brad+0.25*brad
         if(r.le.brad) then
            do ii=1,nbond(i)
               ib=bondlist(i,ii)
               do jj=1,nbond(j)
                  jb=bondlist(j,jj)
                  if(ib.ne.jb.and.ib.ne.j.and.jb.ne.i)then
                     nt=nt+1
                     write(lab3,'(a2,i5,''-'',a2,i5,''-'',&
                        &                               a2,i5,''-'',a2,i5)')&
                        &                  toSymbol(ic(ib)),ib,toSymbol(ic(i)),i,&
                        &                  toSymbol(ic(j)),j,  toSymbol(ic(jb)),jb
                     call rmblank(lab3,lab4)
                     xyz4(1,1)=x(ib)
                     xyz4(2,1)=y(ib)
                     xyz4(3,1)=z(ib)
                     xyz4(1,4)=x(jb)
                     xyz4(2,4)=y(jb)
                     xyz4(3,4)=z(jb)
                     call xdihed(xyz4,ang)
                     write(angsum2(nt),'(a,''='',f6.2)')&
                        &                  trim(lab4),ang
                  endif
               enddo
            enddo
         endif
      enddo
   enddo
   write(*,'(4a)')(angsum2(i),i=1,nt)

   if(iatf2.ne.0.and.iatf1.ne.0)then
      write(*,*)
      write(*,*)'CMA Distance (Angstroems)'
      write(*,*)'---------------------------'
      call cmafrag(n,ic,xyzin,ra,rb)
      write(*,'(''R(CMA):'',f8.4)')rcma*0.529177260d0
   endif

end subroutine

!-------------------------------------------------

subroutine rmblank(as,re)
   implicit integer (i-n)
   character(len=*) as
   character(len=*) re

   re=' '
   k=0
   nsp=ichar(' ')
   do i=1,len(as)
      j=ichar(as(i:i))
      if(j.ne.nsp)then
         k=k+1
         re(k:k)=as(i:i)
      endif
   end do
end subroutine rmblank

subroutine xbangle(xyz,angle)
   implicit integer (i-n)
   implicit double precision (a-h,o-z)
   dimension xyz(3,3)
   !********************************************************************
   !
   ! bangle calculates the angle between atoms i,j, and k. the
   !        cartesian coordinates are in xyz.
   !
   !********************************************************************
   i=1
   j=2
   k=3
   d2ij = (xyz(1,i)-xyz(1,j))**2+&
      &       (xyz(2,i)-xyz(2,j))**2+&
      &       (xyz(3,i)-xyz(3,j))**2
   d2jk = (xyz(1,j)-xyz(1,k))**2+&
      &       (xyz(2,j)-xyz(2,k))**2+&
      &       (xyz(3,j)-xyz(3,k))**2
   d2ik = (xyz(1,i)-xyz(1,k))**2+&
      &       (xyz(2,i)-xyz(2,k))**2+&
      &       (xyz(3,i)-xyz(3,k))**2
   xy = sqrt(d2ij*d2jk)
   temp = 0.5d0 * (d2ij+d2jk-d2ik) / xy
   if (temp .gt. 1.0d0) temp=1.0d0
   if (temp .lt. -1.0d0) temp=-1.0d0
   angle = acos( temp )*180.0d0/3.14159265358979d0
   return
end subroutine xbangle
subroutine xdihed(xyz,angle)
   implicit integer (i-n)
   implicit double precision (a-h,o-z)
   dimension xyz(3,4)
   !********************************************************************
   !
   !      dihed calculates the dihedral angle between atoms i, j, k,
   !            and l.  the cartesian coordinates of these atoms
   !            are in array xyz.
   !
   !     dihed is a modified version of a subroutine of the same name
   !           which was written by dr. w. theil in 1973.
   !
   !********************************************************************
   data pi/3.14159265358979d0/
   i=1
   j=2
   k=3
   l=4
   xi1=xyz(1,i)-xyz(1,k)
   xj1=xyz(1,j)-xyz(1,k)
   xl1=xyz(1,l)-xyz(1,k)
   yi1=xyz(2,i)-xyz(2,k)
   yj1=xyz(2,j)-xyz(2,k)
   yl1=xyz(2,l)-xyz(2,k)
   zi1=xyz(3,i)-xyz(3,k)
   zj1=xyz(3,j)-xyz(3,k)
   zl1=xyz(3,l)-xyz(3,k)
   !      rotate around z axis to put kj along y axis
   dist= sqrt(xj1**2+yj1**2+zj1**2)
   cosa=zj1/dist
   if(cosa.gt.1.0d0) cosa=1.0d0
   if(cosa.lt.-1.0d0) cosa=-1.0d0
   ddd=1.0d0-cosa**2
   if(ddd.le.0.0) go to 10
   yxdist=dist* sqrt(ddd)
   if(yxdist.gt.1.0d-9) go to 20
   10 continue
   xi2=xi1
   xl2=xl1
   yi2=yi1
   yl2=yl1
   costh=cosa
   sinth=0.d0
   go to 30
   20 cosph=yj1/yxdist
   sinph=xj1/yxdist
   xi2=xi1*cosph-yi1*sinph
   xj2=xj1*cosph-yj1*sinph
   xl2=xl1*cosph-yl1*sinph
   yi2=xi1*sinph+yi1*cosph
   yj2=xj1*sinph+yj1*cosph
   yl2=xl1*sinph+yl1*cosph
   !      rotate kj around the x axis so kj lies along the z axis
   costh=cosa
   sinth=yj2/dist
   30 continue
   yi3=yi2*costh-zi1*sinth
   yl3=yl2*costh-zl1*sinth
   call xdang(xl2,yl3,xi2,yi3,angle)
   if (angle .lt. 0.) angle=2.0d0*pi+angle
   if (angle .ge. 2.0d0*pi    ) angle=0.d0
   angle=angle*180.0d0/pi
   return
end subroutine xdihed
subroutine xdang(a1,a2,b1,b2,rcos)
   implicit double precision (a-h,o-z)
   !*********************************************************************
   !
   !    dang  determines the angle between the points (a1,a2), (0,0),
   !          and (b1,b2).  the result is put in rcos.
   !
   !*********************************************************************
   data pi/3.14159265358979d0/
   zero=1.0d-6
   if( abs(a1).lt.zero.and. abs(a2).lt.zero) go to 10
   if( abs(b1).lt.zero.and. abs(b2).lt.zero) go to 10
   anorm=1.0d0/ sqrt(a1**2+a2**2)
   bnorm=1.0d0/ sqrt(b1**2+b2**2)
   a1=a1*anorm
   a2=a2*anorm
   b1=b1*bnorm
   b2=b2*bnorm
   sinth=(a1*b2)-(a2*b1)
   costh=a1*b1+a2*b2
   if(costh.gt.1.0d0) costh=1.0d0
   if(costh.lt.-1.0d0) costh=-1.0d0
   rcos= acos(costh)
   if( abs(rcos).lt.4.0d-4) go to 10
   if(sinth.gt.0.d0) rcos=2.0d0*pi-rcos
   rcos=-rcos
   return
   10 rcos=0.0d0
   return
end subroutine xdang
