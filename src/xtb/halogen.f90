! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> TODO
module xtb_xtb_halogen
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_lin, only : lin
   use xtb_xtb_data
   implicit none
   private

   public :: xbpot


contains


subroutine xbpot(halData,n,at,xyz,xblist,nxb,a,exb,g)
   type(THalogenData), intent(in) :: halData
   integer, intent(in) :: n
   integer, intent(in) :: at(:)
   integer, intent(in) :: nxb
   integer, intent(in) :: xblist(:,:)
   real(wp), intent(in) :: xyz(:,:)
   real(wp), intent(inout) :: g(:,:)
   real(wp), intent(inout) :: exb
   real(wp), intent(in) :: a

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
   !$omp parallel do schedule(runtime) default(none) reduction(+:exb) &
   !$omp shared(nxb, xblist, at, halData, xyz, lj2, alp, a) &
   !$omp private(X, AA, B, ati, atj, cc, r0ax, dxa, dxb, dba, d2ax, &
   !$omp& d2bx, d2ab, rax, XY, aterm, t13, t14, TERM)
   do k=1,nxb
      X =xblist(1,k)
      AA=xblist(2,k)
      B=xblist(3,k)
      ati=at(X)
      atj=at(AA)
      cc=halData%bondStrength(ati)
      r0ax=halData%radScale*(halData%atomicRad(ati)+halData%atomicRad(atj))
      dxa=xyz(:,AA)-xyz(:,X)   ! acceptor - halogen
      dxb=xyz(:, B)-xyz(:,X)   ! neighbor - halogen
      dba=xyz(:,AA)-xyz(:, B)  ! acceptor - neighbor
      d2ax=sum(dxa*dxa)
      d2bx=sum(dxb*dxb)
      d2ab=sum(dba*dba)
      rax=sqrt(d2ax)
      ! angle part. term = cos angle B-X-A
      XY = SQRT(D2BX*D2AX)
      TERM = (D2BX+D2AX-D2AB) / XY
      aterm = (0.5_wp-0.25_wp*term)**alp
      t13 = r0ax/rax
      t14 = t13**a
      exb = exb +  aterm*cc*(t14-halData%dampingPar*t13**lj2) / (1.0_wp+t14)
   enddo

   ! analytic gradient
   !$omp parallel do schedule(runtime) default(none) reduction(+:g) &
   !$omp shared(nxb, xblist, at, halData, xyz, lj2, alp) &
   !$omp private(X, AA, B, ati, atj, cc, r0ax, dxa, dxb, dba, d2ax, &
   !$omp& d2bx, d2ab, rax, XY, aterm, numerator, denominator, termLJ, &
   !$omp& dtermlj, prefactor, dcosterm, t13, t14, rbx, TERM)
   do k=1,nxb
      X =xblist(1,k)
      AA=xblist(2,k)
      B=xblist(3,k)
      ati=at(X)
      atj=at(AA)
      cc=halData%bondStrength(ati)
      r0ax=halData%radScale*(halData%atomicRad(ati)+halData%atomicRad(atj))

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
      numerator = (t14*t14 - halData%dampingPar*t14)
      denominator = (1.0_wp + t14*t14)
      termLJ= numerator/denominator

      ! ----
      ! LJ derivative
      ! denominator part
      dtermlj=2.0_wp*lj2*numerator*t14*t14/(rax*denominator*denominator)
      ! numerator part
      dtermlj=dtermlj+lj2*t14*(halData%dampingPar - 2.0_wp*t14)/(rax*denominator)
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


end module xtb_xtb_halogen
