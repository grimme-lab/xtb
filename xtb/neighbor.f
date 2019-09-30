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


      subroutine getring(n,nbin,a0,c,iring)
      implicit none
      integer n,nbin(20,n),a0,iring,i,nb(20,n)
      integer    i1,i2,i3,i4,i5,i6,i7,i8
      integer n0,n1,n2,n3,n4,n5,n6,n7,n8
      integer    a1,a2,a3,a4,a5,a6,a7,a8
      integer c(8)
      logical chk

      nb=nbin
      do i=1,n
         if(nb(20,i).eq.1)nb(20,i)=0
      enddo

      iring =0
      c(1:8)=0

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
               if(a3.eq.a0.and.chk(n,3,c))then
               iring=3
               goto 99
               endif
               do i4=1,n3
                  a4=nb(i4,a3)
                  n4=nb(20,a4)
                  if(a4.eq.a3) cycle
                  c(4)=a4
                  if(a4.eq.a0.and.chk(n,4,c))then
                     iring=4
                     goto 99
                  endif
                  do i5=1,n4
                     a5=nb(i5,a4)
                     n5=nb(20,a5)
                     if(a5.eq.a4) cycle
                     c(5)=a5
                     if(a5.eq.a0.and.chk(n,5,c))then
                        iring=5
                        goto 99
                     endif
                     do i6=1,n5
                        a6=nb(i6,a5)
                        n6=nb(20,a6)
                        if(a6.eq.a5) cycle
                        c(6)=a6
                        if(a6.eq.a0.and.chk(n,6,c))then
                           iring=6
                           goto 99
                        endif
                        do i7=1,n6
                           a7=nb(i7,a6)
                           n7=nb(20,a7)
                           if(a7.eq.a6) cycle
                           c(7)=a7
                           if(a7.eq.a0.and.chk(n,7,c))then
                              iring=7
                              goto 99
                           endif
                           do i8=1,n7
                              a8=nb(i8,a7)
                              n8=nb(20,a8)
                              if(a8.eq.a7) cycle
                              c(8)=a8
                              if(a8.eq.a0.and.chk(n,8,c))then
                                 iring=8
                                 goto 99
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      return

 99   continue
c     write(*,*) 'atom ',a0,'ring found',c(1:iring)
c     call isort(iring,c)
      return

      end subroutine getring

      logical function chk(nn,n,c)
      implicit none
      integer n,idum(nn),nn,c(8),i,j
      chk=.true.
      idum=0
      do i=1,n
         idum(c(i))=idum(c(i))+1
      enddo
      j=0
      do i=1,nn
         if(idum(i).eq.1) j=j+1
      enddo
      if(j.ne.n) chk=.false.
      end function chk


      logical function samering(n,i,j,c,s)
      implicit none
      integer n,i,j,c(8,n),s(n)
      integer m,k

      samering=.false.
      if(s(i).eq.0   ) return
      if(s(i).ne.s(j)) return

      m=0
      do k=1,s(i)
         m=m+c(k,i)-c(k,j)
      enddo
      if(m.eq.0) then
         samering=.true.
         return
      endif
      
      do k=1,s(i)
         if(c(k,i).eq.j)then
            samering=.true.
            return
         endif
      enddo
      do k=1,s(j)
         if(c(k,j).eq.i)then
            samering=.true.
            return
         endif
      enddo

      end function samering

      subroutine neighborh(natoms,at,xyz,nb)
      use iso_fortran_env, wp => real64
      implicit none  
      integer  at(natoms),natoms,nb(20,natoms)
      real(wp) xyz(3,natoms)

      logical  da
      integer  iat,i,j,k,ni,ii,jj,kk,ll
      real(wp) rad(94)
      real(wp) dx,dy,dz,r,damp,xn,rr,rco,r2,f,a1
      data rad/
     . 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865,
     . 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527,
     . 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820,
     . 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730,
     . 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923,
     . 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188,
     . 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349,
     . 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216,
     . 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717,
     . 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967,
     . 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625,
     . 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657,
     . 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833,
     . 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098,
     . 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878,
     . 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790,
     . 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584,
     . 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289,
     . 3.82984466, 3.85504098, 3.88023730, 3.90543362 /

      nb=0

      do i=1,natoms
      f=1.0
      k=0
 100  do iat=1,natoms
         da=.false.
         do j=1,k
            if(nb(j,i).eq.iat)da=.true.
         enddo
         if(iat.ne.i.and.(.not.da))then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r2=dx*dx+dy*dy+dz*dz 
            r=sqrt(r2)
            rco=rad(at(i))+rad(at(iat))
c critical step            
            if(r.lt.f*rco.and.k.lt.19)then
               k=k+1
               nb(k,i)=iat
            endif
         endif
      enddo
      if(k.lt.1.and.f.lt.1.5)then
         f=f*1.1
         goto 100
      endif
      nb(20,i)=k
      enddo

      end subroutine neighborh

      subroutine neighborhi(natoms,at,sqrab,i,nb)
      use iso_fortran_env, wp => real64
      implicit none  
      integer  at(natoms),natoms,nb(20)
      real(wp) sqrab(natoms*(natoms+1)/2)

      logical  da
      integer  iat,i,j,k,ni,ii,jj,kk,ll,lin
      real(wp) rad(94),r,rr,rco,r2,f,a1
      data rad/
     . 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865,
     . 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527,
     . 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820,
     . 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730,
     . 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923,
     . 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188,
     . 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349,
     . 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216,
     . 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717,
     . 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967,
     . 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625,
     . 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657,
     . 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833,
     . 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098,
     . 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878,
     . 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790,
     . 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584,
     . 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289,
     . 3.82984466, 3.85504098, 3.88023730, 3.90543362 /

      nb=0

      f=1.0
      k=0
 100  do iat=1,natoms
         da=.false.
         do j=1,k
            if(nb(j).eq.iat)da=.true.
         enddo
         if(iat.ne.i.and.(.not.da))then
            r2=sqrab(lin(iat,i))   
            r=sqrt(r2)
            rco=rad(at(i))+rad(at(iat))
c critical step            
            if(r.lt.f*rco.and.k.lt.20)then
               k=k+1
               nb(k)=iat
            endif
         endif
      enddo
      if(k.lt.1.and.f.lt.1.5)then
         f=f*1.1
         goto 100
      endif
      nb(20)=k

      end subroutine neighborhi
