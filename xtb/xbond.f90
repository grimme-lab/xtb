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


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! X..Y bond? (X=halogen but not F, Y=N,O,P,S)

pure elemental function xbond(ati,atj) result(bool)
   integer,intent(in) :: ati,atj
   logical :: bool
   logical :: lx1,lx2,ly1,ly2

   bool=.false.
   lx1=.false.
   lx2=.false.
   ly1=.false.
   ly2=.false.

   if(ati.eq.17.or.ati.eq.35.or.ati.eq.53.or.ati.eq.85) lx1=.true.
   if(atj.eq.17.or.atj.eq.35.or.atj.eq.53.or.atj.eq.85) lx2=.true.
   if(ati.eq. 7.or.ati.eq. 8.or.ati.eq.15.or.ati.eq.16) ly1=.true.
   if(atj.eq. 7.or.atj.eq. 8.or.atj.eq.15.or.atj.eq.16) ly2=.true.

   if(lx1.and.ly2) bool=.true.
   if(lx2.and.ly1) bool=.true.

end function xbond

subroutine xbpot(n,at,xyz,sqrab,xblist,nxb,kk,xbrad,a,exb,g)
   use iso_fortran_env, only : wp => real64
   use aoparam
   use lin_mod, only : lin
   implicit none
   integer n,at(n),nxb
   integer xblist(3,nxb+1)
   real(wp) xyz(3,n)
   real(wp) g  (3,n)
   real(wp) sqrab(n*(n+1)/2)
   real(wp) exb
   real(wp) a,xbrad,kk

   integer  m,k,AA,B,X,ati,atj
   real(wp) cc,r0ax,t13,t14,t16
   real(wp) d2ax,rax,term,aterm,xy,d2bx,d2ab,alp,lj2
   real(wp) er,el,step,t1(3),t2(3)

   alp=6.0
   lj2=0.5*a
   step=0.00001

   exb = 0
   if(nxb.lt.1) return

   ! B-X...A
   do k=1,nxb
      X =xblist(1,k)
      AA=xblist(2,k)
      B=xblist(3,k)
      ati=at(X)
      atj=at(AA)
      cc=cxb(ati)
      ! this sloppy conv. factor has been used in development, keep it
      r0ax=xbrad*(rad(ati)+rad(atj))/0.5291670d0
      d2ax=sqrab(lin(AA,X))
      d2ab=sqrab(lin(AA,B))
      d2bx=sqrab(lin(X, B))
      rax=sqrt(d2ax)
      !        angle part. term = cos angle B-X-A
      XY = SQRT(D2BX*D2AX)
      TERM = (D2BX+D2AX-D2AB) / XY
      aterm = (0.5d0-0.25d0*term)**alp
      t13 = r0ax/rax
      t14 = t13**a
      exb = exb +  aterm*cc*(t14-kk*t13**lj2) / (1.0d0+t14)

      do m=1,3
         xyz(m,X)=xyz(m,X)+step
         t1=xyz(:,AA)-xyz(:,X)
         t2=xyz(:, B)-xyz(:,X)
         d2ax=sum(t1*t1)
         d2bx=sum(t2*t2)
         rax=sqrt(d2ax)
         XY = SQRT(D2BX*D2AX)
         TERM = (D2BX+D2AX-D2AB) / XY
         aterm = (0.5d0-0.25d0*term)**alp
         t13 = r0ax/rax
         t14 = t13**a
         er = aterm*(t14-kk*t13**lj2) / (1.0d0+t14)
         xyz(m,X)=xyz(m,X)-step*2.
         t1=xyz(:,AA)-xyz(:,X)
         t2=xyz(:, B)-xyz(:,X)
         d2ax=sum(t1*t1)
         d2bx=sum(t2*t2)
         rax=sqrt(d2ax)
         XY = SQRT(D2BX*D2AX)
         TERM = (D2BX+D2AX-D2AB) / XY
         aterm = (0.5d0-0.25d0*term)**alp
         t13 = r0ax/rax
         t14 = t13**a
         el = aterm*(t14-kk*t13**lj2) / (1.0d0+t14)
         xyz(m,X)=xyz(m,X)+step
         g  (m,X)=g  (m,X)+cc*(er-el)/(2.*step)
      enddo

      d2bx=sqrab(lin(X, B))
      do m=1,3
         xyz(m,AA)=xyz(m,AA)+step
         t1=xyz(:,AA)-xyz(:,X)
         t2=xyz(:,AA)-xyz(:,B)
         d2ax=sum(t1*t1)
         d2ab=sum(t2*t2)
         rax=sqrt(d2ax)
         XY = SQRT(D2BX*D2AX)
         TERM = (D2BX+D2AX-D2AB) / XY
         aterm = (0.5d0-0.25d0*term)**alp
         t13 = r0ax/rax
         t14 = t13**a
         er = aterm*(t14-kk*t13**lj2) / (1.0d0+t14)
         xyz(m,AA)=xyz(m,AA)-step*2.
         t1=xyz(:,AA)-xyz(:,X)
         t2=xyz(:,AA)-xyz(:,B)
         d2ax=sum(t1*t1)
         d2ab=sum(t2*t2)
         rax=sqrt(d2ax)
         XY = SQRT(D2BX*D2AX)
         TERM = (D2BX+D2AX-D2AB) / XY
         aterm = (0.5d0-0.25d0*term)**alp
         t13 = r0ax/rax
         t14 = t13**a
         el = aterm*(t14-kk*t13**lj2) / (1.0d0+t14)
         xyz(m,AA)=xyz(m,AA)+step
         g  (m,AA)=g  (m,AA)+cc*(er-el)/(2.*step)
      enddo

      d2ax=sqrab(lin(X,AA))
      rax=sqrt(d2ax)
      t13 = r0ax/rax
      t14 = t13**a
      t16 = (t14-kk*t13**lj2) / (1.0d0+t14)
      do m=1,3
         xyz(m,B)=xyz(m,B)+step
         t1=xyz(:, B)-xyz(:,X)
         t2=xyz(:,AA)-xyz(:,B)
         d2bx=sum(t1*t1)
         d2ab=sum(t2*t2)
         XY = SQRT(D2BX*D2AX)
         TERM = (D2BX+D2AX-D2AB) / XY
         aterm = (0.5d0-0.25d0*term)**alp
         er = t16 * aterm
         xyz(m,B)=xyz(m,B)-step*2.
         t1=xyz(:, B)-xyz(:,X)
         t2=xyz(:,AA)-xyz(:,B)
         d2bx=sum(t1*t1)
         d2ab=sum(t2*t2)
         XY = SQRT(D2BX*D2AX)
         TERM = (D2BX+D2AX-D2AB) / XY
         aterm = (0.5d0-0.25d0*term)**alp
         el = t16 * aterm
         xyz(m,B)=xyz(m,B)+step
         g  (m,B)=g  (m,B)+cc*(er-el)/(2.*step)
      enddo

   enddo
end subroutine xbpot

pure elemental function early3d(i) result(bool)
   implicit none
   integer,intent(in) :: i
   logical :: bool
   bool = .false.
   if ((i.ge.21).and.(i.le.24)) bool = .true.
end function early3d

