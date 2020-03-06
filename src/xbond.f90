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


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! X..Y bond? (X=halogen but not F, Y=N,O,P,S) !CB: Cl is formally included, but
! has a parater of zero in GFN1-xTB

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
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_aoparam
   use xtb_lin, only : lin
   implicit none
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   integer, intent(in) :: nxb
   integer, intent(in) :: xblist(3,nxb+1)
   real(wp), intent(in) :: xyz(3,n)
   real(wp), intent(inout) :: g(3,n)
   real(wp), intent(in) :: sqrab(n*(n+1)/2)
   real(wp), intent(inout) :: exb
   real(wp), intent(in) :: a
   real(wp), intent(in) :: xbrad
   real(wp), intent(in) :: kk

   integer :: m,k,AA,B,X,ati,atj
   real(wp) :: cc,r0ax,t13,t14,t16
   real(wp) :: d2ax,rax,term,aterm,xy,d2bx,d2ab,alp,lj2
   real(wp) :: er,el,step,dxa(3),dxb(3),dba(3),dcosterm
   real(wp) :: dtermlj,termlj,prefactor,numerator,denominator,rbx
   alp=6.0_wp
   lj2=0.50_wp*a

   exb = 0.0_wp
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
      r0ax=xbrad*(rad(ati)+rad(atj)) * aatoau
      d2ax=sqrab(lin(AA,X))
      d2ab=sqrab(lin(AA,B))
      d2bx=sqrab(lin(X, B))
      rax=sqrt(d2ax)
      ! angle part. term = cos angle B-X-A
      XY = SQRT(D2BX*D2AX)
      TERM = (D2BX+D2AX-D2AB) / XY
      aterm = (0.5_wp-0.25_wp*term)**alp
      t13 = r0ax/rax
      t14 = t13**a
      exb = exb +  aterm*cc*(t14-kk*t13**lj2) / (1.0_wp+t14)
   enddo

   ! analytic gradient 
   do k=1,nxb
      X =xblist(1,k)
      AA=xblist(2,k)
      B=xblist(3,k)
      ati=at(X)
      atj=at(AA)
      cc=cxb(ati)
      ! this sloppy conv. factor has been used in development, keep it
      r0ax=xbrad*(rad(ati)+rad(atj)) * aatoau

      dxa=xyz(:,AA)-xyz(:,X)   ! acceptor - halogen
      dxb=xyz(:, B)-xyz(:,X)   ! neighbor - halogen 
      dba=xyz(:,AA)-xyz(:, B)  ! acceptor - neighbor

      d2ax=sum(dxa*dxa)
      d2bx=sum(dxb*dxb)
      d2ab=sum(dba*dba)
      rax=sqrt(d2ax)+1.0e-18_wp
      rbx=sqrt(d2bx)+1.0e-18_wp
      
      XY = SQRT(D2BX*D2AX)
      TERM = (D2BX+D2AX-D2AB) / XY
      ! now compute angular damping function
      aterm = (0.5_wp-0.25_wp*term)**alp

      ! set up weighted inverted distance and compute the modified Lennard-Jones potential
      t14 = (r0ax/rax)**lj2 ! (rov/r)^lj2 ; lj2 = 6 in GFN1
      numerator = (t14*t14 - kk*t14)
      denominator = (1.0_wp + t14*t14)
      termLJ= numerator/denominator

      ! ----
      ! LJ derivative
      ! denominator part
      dtermlj=2.0_wp*lj2*numerator*t14*t14/(rax*denominator*denominator)
      ! numerator part
      dtermlj=dtermlj+lj2*t14*(kk - 2.0_wp*t14)/(rax*denominator) 
      ! scale w/ angular damping term
      dtermlj=dtermlj*aterm*cc/rax
      ! gradient for the acceptor 
      g(:,AA)=g(:,AA)+dtermlj*dxa(:)
      ! halogen gradient
      g(:,X)=g(:,X)-dtermlj*dxa(:)
      ! ---- 
      ! cosine term derivative
      prefactor=-0.250_wp*alp*(0.5_wp-0.25_wp*term)**(alp-1.0_wp)
      prefactor=prefactor*cc*termlj
      ! AX part
      dcosterm=2.0_wp/rbx - term/rax
      dcosterm=dcosterm*prefactor/rax
      ! gradient for the acceptor 
      g(:,AA)=g(:,AA)+dcosterm*dxa(:)
      ! halogen gradient
      g(:,X)=g(:,X)-dcosterm*dxa(:)
      ! BX part
      dcosterm=2.0_wp/rax - term/rbx
      dcosterm=dcosterm*prefactor/rbx
      ! gradient for the acceptor 
      g(:,B)=g(:,B)+dcosterm*dxb(:)
      ! halogen gradient
      g(:,X)=g(:,X)-dcosterm*dxb(:)
      ! AB part
      t13=2.0_wp*prefactor/xy
      ! acceptor 
      g(:,AA)=g(:,AA)-t13*dba(:)
      ! neighbor
      g(:,B)=g(:,B)+t13*dba(:)

   enddo
end subroutine xbpot

pure elemental function early3d(i) result(bool)
   implicit none
   integer,intent(in) :: i
   logical :: bool
   bool = .false.
   if ((i.ge.21).and.(i.le.24)) bool = .true.
end function early3d

