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

subroutine Dints(n,nbf,xyz,S1,S2,S3,basis)
   use iso_fortran_env, wp => real64
   use mctc_constants
   use tbdef_basisset
   use intpack, only : opab1,propa
   implicit none
   type(tb_basisset), intent(in) :: basis
   integer, intent(in)  :: n
   integer, intent(in)  :: nbf
   real(wp),intent(out) :: S1(nbf*(nbf+1)/2)
   real(wp),intent(out) :: S2(nbf*(nbf+1)/2)
   real(wp),intent(out) :: S3(nbf*(nbf+1)/2)
   real(wp),intent(in)  :: xyz(3,n)

   integer  :: i,j,k,l
   integer  :: iprimcount,jprimcount
   integer  :: npri,nprj
   integer  :: ii,iii,jj,jjj,ll,m,li,lj,mm,nn
   real(wp) :: xyza(3),xyzb(3),rab,est,ss,sss,lmnfak(84),gama,arg
   real(wp) :: point(3),gm2,ttt(3),tt1,tt2,tt3,intcut

   intcut=20.0_wp

   point=0.0_wp
   s1=0.0_wp
   s2=0.0_wp
   s3=0.0_wp

   k=0
   iprimcount=0
   do i=1,nbf
      ! aufpunkt i
      xyza(1:3)=xyz(1:3,basis%aoat(i))
      ! #prims
      npri=basis%nprim(i)
      jprimcount=0
      do j=1,i
         k=k+1
         nprj=basis%nprim(j)
         ! aufpunkt j
         xyzb(1:3)=xyz(1:3,basis%aoat(j))
         rab=sum((xyza-xyzb)**2)
         if(rab.gt.200) goto 42 ! cut-off gives crambin dipole accurate to 1d-3 Deb
         !prim loop
         tt1=0.0_wp
         tt2=0.0_wp
         tt3=0.0_wp
         do ii=1,npri
            iii=iprimcount+ii
            do jj=1,nprj
               jjj=jprimcount+jj
               gama=1.0_wp/(basis%alp(iii)+basis%alp(jjj))
               est=rab*basis%alp(iii)*basis%alp(jjj)*gama
               ! cutoff
               if(est.lt.intcut)then
                  ttt=0
                  call propa(opab1,xyza,xyzb,point,basis%alp(iii),basis%alp(jjj),&
                     &       basis%lao(i),basis%lao(j),ttt,3)
                  tt1=tt1+ttt(1)*basis%cont(iii)*basis%cont(jjj)
                  tt2=tt2+ttt(2)*basis%cont(iii)*basis%cont(jjj)
                  tt3=tt3+ttt(3)*basis%cont(iii)*basis%cont(jjj)
               endif
            enddo
         enddo
         s1(k)=tt1
         s2(k)=tt2
         s3(k)=tt3
 42      jprimcount=jprimcount+nprj
      enddo
      iprimcount=iprimcount+npri
   enddo

end subroutine Dints

subroutine calc_dipole(n,at,xyz,z,nao,P,dpint,dip,d)
   use iso_fortran_env, wp => real64
   use mctc_econv
   implicit none
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: xyz(3,n)
   real(wp),intent(in) :: z(n)
   integer, intent(in) :: nao
   real(wp),intent(in) :: P(nao,nao)
   real(wp),intent(in) :: dpint(3,nao*(1+nao)/2)

   integer  :: i,j,k
   real(wp),intent(out) :: d(3),dip

   ! core part
   d = 0.0_wp
   do i = 1, n
      d = d + xyz(:,i)*z(i)
   enddo

   ! contraction with P
   k = 0
   do i = 1, nao
      do j = 1, i-1
         k = k+1
         d = d - 2.0_wp*P(j,i)*dpint(:,k)
      enddo
      k = k+1
      d = d - P(i,i)*dpint(:,k)
   enddo

   dip = norm2(d)

end subroutine calc_dipole


