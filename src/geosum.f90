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

subroutine geosum(n,ic,xyzin)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_internals_type, only : angle_value, dihedral_value
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
   real*8 ang,ra(3),rb(3)
   integer i, j, nn, k, imin, jmin, iat1, iat2
   integer nm,ij,m,iat,jat, nsel
   integer nt,a,b,ia,ib,ja,jb,ii,jj,lin,linij
   character*31  lab,lab2
   character*31, allocatable :: angsum1(:)
   character*30  lab4
   character*30, allocatable :: angsum2(:)
   character*50  lab3
   character*21  lab21

   ! first nsel atoms are considered in angle/dihedreal xtb_printout
   nsel=min(10,n)

   nn=100*100
   allocate(x(n),y(n),z(n),rmin(nn),rmax(nn),dev(nn),av(nn),&
      &         rlist(n,150),bondlist(n,150),nbond(n),nty(nn),&
      &         angsum1(n*150),angsum2(n*150))

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
         do a=1,m
            ia=bondlist(i,a)
            do b=1,a-1
               ib=bondlist(i,b)
               nt=nt+1
               write(lab,'(a2,i5,''-'',a2,i5,''-'',a2,i5)')&
                  &toSymbol(ic(ia)),ia,toSymbol(ic(i)),i,toSymbol(ic(ib)),ib
               call rmblank(lab,lab2)
               ang = angle_value(xyzin, ia, i, ib) * 180.0_wp / pi
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
      do j=1,i
         if(i.eq.j)cycle
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
                     ang = modulo(dihedral_value(xyzin, ib, i, j, jb), 2.0_wp * pi) * 180.0_wp / pi
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

