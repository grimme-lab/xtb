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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c compute localized "normal" modes according to the Pipek-Mezey algo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine locmode(n,n3,at,xyz,vthr,freq,rmass,uin,
     .                   ng,group)
         use iso_fortran_env, wp => real64
      implicit none
      integer n,n3,at(n),ng,group(n)
      real(wp) xyz(3,n),freq(n3),rmass(n3),uin(n3,n3),vthr

      integer nvar,nvsq,r,s,a,ii,ij,ija,k,i,j,highest,ga,n36
      real(wp) dum1,dum2,dd,vtyp(3),tors(n3)
      real(wp),allocatable:: ul(:,:), fl(:), pop(:,:), rl(:), pp(:)
      real(wp),allocatable:: d(:,:),u(:,:),bmat(:,:),rr(:), fl2(:)
      integer,allocatable::ind(:),ind2(:)
      
      write(*,*)
      write(*,*)'Pipek-Mezey localization of normal modes'
      do i=1,n3
         if(freq(i).lt.vthr) highest=i 
      enddo 

      nvar=highest-6
      write(*,*) 'number of localized modes ',nvar
      write(*,*) 'frequency thrshold        ',vthr
      write(*,*) 'number of groups          ',ng    
      write(*,*) 'atom-group assingment :   '
      write(*,'(30i3)') group

      nvsq=nvar*(nvar+1)/2
      allocate ( pop(nvsq,ng),ul(n3,n3),u(n3,nvar),
     .           fl(nvar),rl(nvar),rr(n3),fl2(n3),pp(ng),
     .           d(nvar,nvar),ind(nvar),ind2(n3) )

c skip tr/rot
      do i=7,6+nvar    
         u(1:n3,i-6)=uin(1:n3,i)
      enddo

      pop = 0
c do the popmat
      k=0
      do r=1,nvar
         do s=1,r-1
            k=k+1
            pp=0
            do a=1, n
            ga=group(a)
            do ij=1,3
            ija=ij+(a-1)*3
            pp(ga)=pp(ga)+u(ija,r)*u(ija,s)
            enddo
            enddo
            pop(k,1:ng)=pp(1:ng)*2.0d0
         enddo
         k=k+1
         do a=1, n
         ga=group(a)
         do ij=1,3
         ija=ij+(a-1)*3
         pop(k,ga)=pop(k,ga)+u(ija,r)*u(ija,r)
         enddo
         enddo
      enddo
      write(*,*) 'pop matrix done.'
c do the trafo
      call lopt(.true.,nvar,ng,1.d-6,pop,d)
c new freq, rmass
      do i=1,nvar
         dum1=0.0d0
         dum2=0.0d0
         do k=1,nvar
            dd=d(k,i)*d(k,i)
            dum1=dum1+freq (k+6)*dd           
            dum2=dum2+rmass(k+6)*dd           
         enddo
         fl(i)=dum1
         rl(i)=dum2
      enddo
c the local modes in ul
      CALL DGEMM('N','N',n3,nvar,nvar,1.D0,u,n3,d,nvar,0.D0,ul,n3)

c sort
      do i=1,nvar
         ind(i)=i
      enddo
      call Qsort(fl, 1, nvar, ind) 

      do i=1,nvar
         freq (i+6)=fl(    i )
         rmass(i+6)=rl(ind(i))
         uin(1:n3,i+6)=ul(1:n3,ind(i))
      enddo

      end
