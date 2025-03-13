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

subroutine constralltors(n,at,xyz)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_scanparam
   use xtb_intmodes
   implicit none
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: xyz(3,n)

   real(wp) :: phi,valijkl,a1,a2
   integer  :: i,j,k,ii,jj,kk,ll,irm
   integer  :: nb(0:20,n),list(n,n)

   call neighbor(n,xyz,at,nb)
   list = 0
   do i=1,n
      k=nb(20,i)
      do j=1,k
         list(nb(j,i),i)=1
      enddo
   enddo

   irm=0
   k=0
   do i=1,n-1
      if(nb(20,i).le.1) cycle
      do j=i+1,n
         if(nb(20,j).le.1) cycle
         if(list(j,i).eq.1)then
            do ii=1,nb(20,i)
               kk=nb(ii,i)
               do jj=1,nb(20,j)
                  ll=nb(jj,j)
                  if(kk.eq.j) cycle
                  if(ll.eq.i) cycle
                  call BANGLE(XYZ,I,J,ll,A1)       ! avoid linear cas
                  call BANGLE(XYZ,J,I,KK,A2)
                  if(abs(pi-a1).lt.0.2.or.abs(pi-a2).lt.0.2)then
                     irm=irm+1
                     cycle
                  endif
                  k=k+1
                  atconstr(1,k)=kk
                  atconstr(2,k)=i
                  atconstr(3,k)=j
                  atconstr(4,k)=ll
                  phi=valijkl(n,xyz,kk,i,j,ll)
                  valconstr(k)=phi
                  !                    write(*,*) atconstr(1:4,k),phi*180./3.14159
               enddo
            enddo
         endif
      enddo
   enddo
   write(*,*) 'constraining ',k,' torsions'
   write(*,*) irm,' near linear torsions not included'
   nconstr=nconstr+k

end subroutine constralltors

subroutine constrallangles(n,at,xyz)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_scanparam
   implicit none
   integer n, at(n)
   real(wp) xyz(3,n)

   real(wp) phi,valijkl,a1,a2
   integer i,j,k,ii,jj,kk,ll,irm,nb(0:20,n),list(n,n)

   call neighbor(n,xyz,at,nb)
   list = 0
   do i=1,n
      k=nb(20,i)
      do j=1,k
         list(nb(j,i),i)=1
      enddo
   enddo

   irm=0
   k=0
   do i=1,n
      if(nb(20,i).le.1) cycle
      do ii=1,nb(20,i)
         kk=nb(ii,i)
         do jj=1,nb(20,i)
            ll=nb(jj,i)
            call bangl(xyz,kk,i,ll,phi)
            if(abs(pi-phi).lt.0.2)then
               irm=irm+1
               cycle
            endif
            k=k+1
            atconstr(1,k)=kk
            atconstr(2,k)=i
            atconstr(3,k)=ll
            valconstr(k)=phi
            !              write(*,*) atconstr(1:3,k),phi*180./3.14159
         enddo
      enddo
   enddo
   write(*,*) 'constraining ',k,' angles'
   write(*,*) irm,' near linear angles not included'
   nconstr=nconstr+k

end subroutine constrallangles

subroutine constrallbonds(nat,at,xyz)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autoaa
   use xtb_param_atomicrad, only : atomicRad
   use xtb_scanparam
   implicit none
   integer nat,at(nat)
   real(wp) xyz(3,nat)

   real(wp) rco,rij
   integer i,j

   do i = 1, nat
      do j =1, i
         if(i.eq.j) cycle
         rij=sqrt(sum((xyz(:,i)-xyz(:,j))**2))
         rco=(atomicRad(at(j))+atomicRad(at(i)))*autoaa
         if(0.52917726*rij.lt.1.2*rco)then
            nconstr = nconstr + 1
            atconstr(1,nconstr) = i
            atconstr(2,nconstr) = j
            valconstr(nconstr)  = rij
            !              write(*,*) at(i),at(j),rco,0.52917726*rij
         endif
      enddo
   enddo

   write(*,*) 'constraining ',nconstr,' bonds'

end subroutine constrallbonds

subroutine constrpot(nat,at,xyz,g,e)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_scanparam
   use xtb_splitparam
   implicit none
   integer nat,at(nat)
   real(wp) e
   real(wp) xyz(3,nat)
   real(wp) g(3,nat)

   integer i,j,k,l,kk
   real(wp)rij(3),dum,r0,ffc,vp(3),ra(3),rb(3)
   real(wp)va(3),vb(3),vc(3),vab(3),vcb(3),deda(3),dedc(3),dedb(3)
   real(wp)dda(3),ddb(3),ddc(3),ddd(3)
   real(wp)c0,rab2,rcb2,theta,dt,rmul1,rmul2,deddt,rp,cosa,ea,vlen,d
   real(wp)c1,phi0,phi,dphi1,x1cos,x1sin,dij,valijkl,ff,dum1,dum2
   real(wp)dx,dy,dz,termx,termy,termz,abcx,abcy,abcz,a,z0,expo,expom1

   if(nconstr.eq.0) return

   ffc=fcconstr

   do kk=1,nconstr
      ! CMA distance to valconstr (usually=0) for atoms in fragment1
      if(atconstr(1,kk).eq.-2)then
         rb=0
         do i=1,iatf1
            a=atmass(i)
            rb(1)=rb(1)+a*xyz(1,i)
            rb(2)=rb(2)+a*xyz(2,i)
            rb(3)=rb(3)+a*xyz(3,i)
         enddo
         z0=valconstr(kk)
         ra(1)=rb(1)/massf1 - z0
         ra(2)=rb(2)/massf1 - z0
         ra(3)=rb(3)/massf1 - z0
         d=ra(1)**2+ra(2)**2+ra(3)**2
         ffc=fcconstr/dble(iatf1)
         e=e+ffc*d
         do i=1,iatf1
            ff=2.0d0*ffc*atmass(i)/massf1
            g(1,i)=g(1,i)+ra(1)*ff
            g(2,i)=g(2,i)+ra(2)*ff
            g(3,i)=g(3,i)+ra(3)*ff
         enddo
      endif
      ! Z distance for atoms in fragment2
      if(atconstr(1,kk).eq.-1)then
         ffc=0.1*fcconstr/dble(iatf2)
         z0=valconstr(kk)
         do i=1,nat
            if(splitlist(i).eq.2)then
               dz=xyz(3,i)-z0
               e=e+ffc*dz**2
               g(3,i)=g(3,i)+2.0d0*ffc*dz
            endif
         enddo
      endif
      ! CMA distance (rcma in common block for average in MD)
      if(atconstr(1,kk).eq.0)then
         expo=6.0d0
         expom1=5.0d0
         ffc=0.1*fcconstr
         r0=valconstr (kk)
         call cmafrag(nat,at,xyz,ra,rb) ! rcma in common
         d=rcma-r0
         e=e+(1./expo)*ffc*d**expo
         do i=1,nat
            dx = (ra(2)-rb(2))**2+(ra(3)-rb(3))**2
            dy = (ra(1)-rb(1))**2+(ra(3)-rb(3))**2
            dz = (ra(2)-rb(2))**2+(ra(1)-rb(1))**2
            if(splitlist(i).eq.1)then      ! fragment 1
               a  = ffc*atmass(i)/massf1
               abcx=ra(1)-rb(1)
               abcy=ra(2)-rb(2)
               abcz=ra(3)-rb(3)
               termx=sqrt(abcx**2+dx)
               termy=sqrt(abcy**2+dy)
               termz=sqrt(abcz**2+dz)
               g(1,i)=g(1,i)+a*abcx*(termx-r0)**expom1/termx
               g(2,i)=g(2,i)+a*abcy*(termy-r0)**expom1/termy
               g(3,i)=g(3,i)+a*abcz*(termz-r0)**expom1/termz
            else                          ! fragment 2
               a  = ffc*atmass(i)/massf2
               abcx=rb(1)-ra(1)
               abcy=rb(2)-ra(2)
               abcz=rb(3)-ra(3)
               termx=sqrt(abcx**2+dx)
               termy=sqrt(abcy**2+dy)
               termz=sqrt(abcz**2+dz)
               g(1,i)=g(1,i)+a*abcx*(termx-r0)**expom1/termx
               g(2,i)=g(2,i)+a*abcy*(termy-r0)**expom1/termy
               g(3,i)=g(3,i)+a*abcz*(termz-r0)**expom1/termz
            endif
         enddo
      endif
      ! length (rcma in common block for average in MD)
      if(atconstr(3,kk).eq.0.and.atconstr(1,kk).ne.0&
         &.and.atconstr(2,kk).ne.0)then
         i= atconstr(1,kk)
         j= atconstr(2,kk)
         r0=valconstr (kk)
         rij=xyz(:,j)-xyz(:,i)
         rcma=sqrt(sum(rij*rij))
         d=rcma-r0
         !           e=e+ffc*d*d
         !           ff=ffc*2.0d0*d
         dum= d**springexpo
         dum2=d**(springexpo-1.)
         e=e+ffc*dum
         ff=ffc*springexpo*dum2
         if(kk.eq.1) pmf=ff
         dum=ff/rcma
         g(1,j)=g(1,j)+dum*rij(1)
         g(1,i)=g(1,i)-dum*rij(1)
         g(2,j)=g(2,j)+dum*rij(2)
         g(2,i)=g(2,i)-dum*rij(2)
         g(3,j)=g(3,j)+dum*rij(3)
         g(3,i)=g(3,i)-dum*rij(3)
      endif
      ! bend
      if (atconstr(4,kk).eq.0.and.atconstr(3,kk).gt.0) then
         i= atconstr(1,kk)
         j= atconstr(2,kk)
         k= atconstr(3,kk)
         c0 =valconstr(kk)
         va = xyz(1:3,i)
         vb = xyz(1:3,j)
         vc = xyz(1:3,k)
         vab = va-vb
         vcb = vc-vb
         rab2 = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
         rcb2 = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
         call crprod(vcb,vab,vp)
         rp = norm2(vp)+1.d-14
         call impsc(vab,vcb,cosa)
         cosa = dble(min(1.0d0,max(-1.0d0,cosa)))
         theta= dacos(cosa)
         dt  = theta - c0
         ea  = ffc* dt**2
         deddt = 2.d0 * ffc * dt
         e = e + ea
         call crprod(vab,vp,deda)
         rmul1 = -deddt / (rab2*rp)
         deda = deda*rmul1
         call crprod(vcb,vp,dedc)
         rmul2 =  deddt / (rcb2*rp)
         dedc = dedc*rmul2
         dedb = deda+dedc
         g(1:3,i) = g(1:3,i) + deda(1:3)
         g(1:3,j) = g(1:3,j) - dedb(1:3)
         g(1:3,k) = g(1:3,k) + dedc(1:3)
      endif
      ! dihed
      if(atconstr(4,kk).ne.0)then
         i= atconstr(1,kk)
         j= atconstr(2,kk)
         k= atconstr(3,kk)
         l= atconstr(4,kk)
         phi0 =valconstr(kk)
         phi=valijkl(nat,xyz,i,j,k,l)
         if(abs(phi-pi).lt.1.d-8.or.abs(phi).lt.1.d-8)phi=phi+1.d-8
         call dphidr(nat,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)
         dphi1=phi0-phi+pi
         x1cos=cos(dphi1)
         x1sin=sin(dphi1)
         dij  =ffc*x1sin
         g(1:3,i)=g(1:3,i)+dij*dda(1:3)
         g(1:3,j)=g(1:3,j)+dij*ddb(1:3)
         g(1:3,k)=g(1:3,k)+dij*ddc(1:3)
         g(1:3,l)=g(1:3,l)+dij*ddd(1:3)
         e=e+ffc*(1.+x1cos)
      endif

   enddo

end subroutine constrpot

subroutine constrhess(nat,at,xyz0,Hess)
   use xtb_mctc_accuracy, only : wp
   use xtb_scanparam
   use xtb_constrainpot
   use xtb_fixparam
   use xtb_metadynamic
   implicit none
   integer, intent(in)    :: nat
   integer, intent(in)    :: at(nat)
   real(wp),intent(inout) :: Hess((nat*3)*((nat*3)+1)/2)
   real(wp),intent(in)    :: xyz0(3,nat)

   integer ia,ic,ja,jc,i,j,k,n,n3,ii,jj
   real(wp) step,step2,e
   real(wp),allocatable :: gr(:,:),gl(:,:)
   real(wp),allocatable :: xyz(:,:)
   real(wp),allocatable :: h(:,:)
   allocate(gr(3,nat),gl(3,nat),xyz(3,nat),h(nat*3,nat*3))

   call qpothess2(potset%pos,nat,at,xyz0,hess)

   if(nconstr.eq.0 .and. potset%dist%n.eq.0 .and. potset%angle%n.eq.0&
      &.and. potset%dihedral%n.eq.0 ) return

   xyz=xyz0
   n =nat
   n3=3*nat
   step =0.0001
   step2=0.5d0/step
   e=0

   do ia = 1, n
      do ic = 1, 3
         ii = (ia-1)*3+ic
         xyz(ic,ia)=xyz(ic,ia)+step
         gr=0
         call constrpot(nat,at,xyz,gr,e)
         call constrain_dist    (potset%dist,    nat,at,xyz,gr,e)
         call constrain_angle   (potset%angle,   nat,at,xyz,gr,e)
         call constrain_dihedral(potset%dihedral,nat,at,xyz,gr,e)
         call metadynamic(metaset,nat,at,xyz,e,gr)
         xyz(ic,ia)=xyz(ic,ia)-2.*step
         gl=0
         call constrpot(nat,at,xyz,gl,e)
         call constrain_dist    (potset%dist,    nat,at,xyz,gl,e)
         call constrain_angle   (potset%angle,   nat,at,xyz,gl,e)
         call constrain_dihedral(potset%dihedral,nat,at,xyz,gl,e)
         call metadynamic(metaset,nat,at,xyz,e,gl)
         xyz(ic,ia)=xyz(ic,ia)+step
         do ja = 1, n
            do jc = 1, 3
               jj = (ja-1)*3 + jc
               h(jj,ii) =(gr(jc,ja) - gl(jc,ja)) * step2
            enddo
         enddo
      enddo
   enddo

   !     call prmat(6,h,3*n,3*n,'H')

   k=0
   do i=1,n3
      do j=1,i
         k=k+1
         hess(k)=hess(k)+0.5*(h(j,i)+h(i,j))
      enddo
   enddo

end subroutine

!  .....................................................................

Subroutine dthetadr(nat,xyz,j,k,i,theta,kijk,dei,dej,dek)
   use xtb_mctc_accuracy, only : wp

   !  .....................................................................
   !
   !   Input
   !   -----
   !                              __              __
   !    ra,rb          : distance ij and distnace ik
   !
   !    sinth          : sinus of angle j-i-k
   !
   !    costh          : cosinus of angle j-i-k
   !
   !
   !   Output
   !   ------
   !
   !    dei            : derivatives of the atom i
   !
   !    dej            : derivatives of the atom j
   !
   !    dek            : derivatives of the atom k
   !
   !  .....................................................................

   implicit none

   external vecnorm

   integer ic,nat,j,k,i

   real(wp) :: &
      &         ran,rbn,theta,xyz(3,nat),&
      &         ra(3),rb(3),sinth,costh,&
      &         dei(3),dej(3),dek(3),kijk,&
      &         fac,vecnorm,eps

   parameter (eps=1.d-14)

   ! ... dtheta/dri, dtheta/drj and dtheta/drk

   sinth=sin(theta)
   costh=cos(theta)

   do ic=1,3
      ra(ic)=xyz(ic,j)-xyz(ic,i)
      rb(ic)=xyz(ic,k)-xyz(ic,i)
   end do

   ran=vecnorm(ra,3,0)
   rbn=vecnorm(rb,3,0)
   if (abs(sinth).le.eps) then
      fac=sign(kijk,costh)/(ran*rbn)
   else
      fac=1.d0/(sinth*ran*rbn)
   end if

   do ic=1,3
      ! ... dedri
      dei(ic)=fac*( ra(ic)+rb(ic)-(rbn*rbn*ra(ic)+ran*ran*rb(ic))*&
         &                  costh/(ran*rbn))
      ! ... dedrj
      dej(ic)=fac*(ra(ic)*costh*rbn/ran-rb(ic))
      ! ... dedrk
      dek(ic)=fac*(rb(ic)*costh*ran/rbn-ra(ic))
   end do

End subroutine





Subroutine dphidr(nat,xyz,i,j,k,l,phi,&
      &                  dphidri,dphidrj,dphidrk,dphidrl)
   use xtb_mctc_accuracy, only : wp
   !     the torsion derivatives

   implicit none

   external vecnorm

   integer  ic,i,j,k,l,nat

   real(wp)&
      &         sinphi,cosphi,onenner,thab,thbc,&
      &         ra(3),rb(3),rc(3),rab(3),rac(3),rbc(3),rbb(3),&
      &         raa(3),rba(3),rapba(3),rapbb(3),rbpca(3),rbpcb(3),&
      &         rapb(3),rbpc(3),na(3),nb(3),nan,nbn,&
      &         dphidri(3),dphidrj(3),dphidrk(3),dphidrl(3),&
      &         xyz(3,nat),phi,vecnorm,nenner,eps,vz

   parameter (eps=1.d-14)

   cosphi=cos(phi)
   sinphi=sin(phi)
   do ic=1,3
      ra(ic)=xyz(ic,j)-xyz(ic,i)
      rb(ic)=xyz(ic,k)-xyz(ic,j)
      rc(ic)=xyz(ic,l)-xyz(ic,k)

      rapb(ic)=ra(ic)+rb(ic)
      rbpc(ic)=rb(ic)+rc(ic)
   end do

   call crossprod(ra,rb,na)
   call crossprod(rb,rc,nb)
   nan=vecnorm(na,3,0)
   nbn=vecnorm(nb,3,0)

   nenner=nan*nbn*sinphi
   if (abs(nenner).lt.eps) then
      dphidri=0
      dphidrj=0
      dphidrk=0
      dphidrl=0
      onenner=1.0d0/(nan*nbn)
   else
      onenner=1.d0/nenner
   endif

   call crossprod(na,rb,rab)
   call crossprod(nb,ra,rba)
   call crossprod(na,rc,rac)
   call crossprod(nb,rb,rbb)
   call crossprod(nb,rc,rbc)
   call crossprod(na,ra,raa)

   call crossprod(rapb,na,rapba)
   call crossprod(rapb,nb,rapbb)
   call crossprod(rbpc,na,rbpca)
   call crossprod(rbpc,nb,rbpcb)

   ! ... dphidri
   do ic=1,3
      dphidri(ic)=onenner*(cosphi*nbn/nan*rab(ic)-rbb(ic))

      ! ... dphidrj
      dphidrj(ic)=onenner*(cosphi*(nbn/nan*rapba(ic)&
         &                                +nan/nbn*rbc(ic))&
         &                        -(rac(ic)+rapbb(ic)))
      ! ... dphidrk
      dphidrk(ic)=onenner*(cosphi*(nbn/nan*raa(ic)&
         &                             +nan/nbn*rbpcb(ic))&
         &                        -(rba(ic)+rbpca(ic)))
      ! ... dphidrl
      dphidrl(ic)=onenner*(cosphi*nan/nbn*rbb(ic)-rab(ic))
   end do

End subroutine


Subroutine dphidrPBC(mode,nat,xyz,i,j,k,l,vTrR,vTrB,vTrC,phi,&
      &                  dphidri,dphidrj,dphidrk,dphidrl)
   use xtb_mctc_accuracy, only : wp
   !     the torsion derivatives

   implicit none

   external vecnorm

   integer  mode,ic,i,j,k,l,nat

   real(wp)&
      &         vTrR(3), vTrB(3), vTrC(3), & 
      &         sinphi,cosphi,onenner,thab,thbc,&
      &         ra(3),rb(3),rc(3),rab(3),rac(3),rbc(3),rbb(3),&
      &         raa(3),rba(3),rapba(3),rapbb(3),rbpca(3),rbpcb(3),&
      &         rapb(3),rbpc(3),na(3),nb(3),nan,nbn,&
      &         dphidri(3),dphidrj(3),dphidrk(3),dphidrl(3),&
      &         xyz(3,nat),phi,vecnorm,nenner,eps,vz

   parameter (eps=1.d-14)

   cosphi=cos(phi)
   sinphi=sin(phi)
   if(mode.eq.1)then
   do ic=1,3
      ra(ic)=xyz(ic,j)+vTrB(ic)-xyz(ic,i)-vTrR(ic)
      rb(ic)=xyz(ic,k)+vTrC(ic)-xyz(ic,j)-vTrB(ic)
      rc(ic)=xyz(ic,l)-xyz(ic,k)-vTrC(ic)

      rapb(ic)=ra(ic)+rb(ic)
      rbpc(ic)=rb(ic)+rc(ic)
   end do
   elseif(mode.eq.2) then
   do ic=1,3
     ra(ic) = -(xyz(ic,i)+vTrC(ic))+xyz(ic,j)
     rb(ic) = -xyz(ic,j)+(xyz(ic,k)+vTrR(ic))
     rc(ic) = -(xyz(ic,k)+vTrR(ic))+(xyz(ic,l)+vTrB(ic))

      rapb(ic)=ra(ic)+rb(ic)
      rbpc(ic)=rb(ic)+rc(ic)
   end do
   endif
   call crossprod(ra,rb,na)
   call crossprod(rb,rc,nb)
   nan=vecnorm(na,3,0)
   nbn=vecnorm(nb,3,0)

   nenner=nan*nbn*sinphi
   if (abs(nenner).lt.eps) then
      dphidri=0
      dphidrj=0
      dphidrk=0
      dphidrl=0
      if (abs(nan*nbn).gt.eps) then
         onenner=1.0d0/(nan*nbn)
      else
         onenner=0.0d0
      endif
   else
      onenner=1.d0/nenner
   endif

   call crossprod(na,rb,rab)
   call crossprod(nb,ra,rba)
   call crossprod(na,rc,rac)
   call crossprod(nb,rb,rbb)
   call crossprod(nb,rc,rbc)
   call crossprod(na,ra,raa)

   call crossprod(rapb,na,rapba)
   call crossprod(rapb,nb,rapbb)
   call crossprod(rbpc,na,rbpca)
   call crossprod(rbpc,nb,rbpcb)

   if (abs(onenner).gt.eps) then
      do ic=1,3
         ! ... dphidri
         dphidri(ic)=onenner*(cosphi*nbn/nan*rab(ic)-rbb(ic))
         ! ... dphidrj
         dphidrj(ic)=onenner*(cosphi*(nbn/nan*rapba(ic)&
            &                                +nan/nbn*rbc(ic))&
            &                        -(rac(ic)+rapbb(ic)))
         ! ... dphidrk
         dphidrk(ic)=onenner*(cosphi*(nbn/nan*raa(ic)&
            &                             +nan/nbn*rbpcb(ic))&
            &                        -(rba(ic)+rbpca(ic)))
         ! ... dphidrl
         dphidrl(ic)=onenner*(cosphi*nan/nbn*rbb(ic)-rab(ic))
      end do
   else
      dphidri=0.0d0
      dphidrj=0.0d0
      dphidrk=0.0d0
      dphidrl=0.0d0
   endif

End subroutine

!  .....................................................................

Subroutine domegadr&
      &           (nat,xyz,&
      &            i,j,k,l,omega,&
      &            domegadri,domegadrj,domegadrk,domegadrl)
   use xtb_mctc_accuracy, only : wp

   !     inversion derivatives
   !  .....................................................................

   implicit none

   external vecnorm

   integer  ic,i,j,k,l,nat

   real(wp)&
      &         omega,sinomega,&
      &         xyz(3,nat),onenner,vecnorm,rnn,rvn,&
      &         rn(3),rv(3),rd(3),re(3),rdme(3),rve(3),&
      &         rne(3),rdv(3),rdn(3),&
      &         rvdme(3),rndme(3),nenner,&
      &         domegadri(3),domegadrj(3),domegadrk(3),domegadrl(3),eps

   parameter (eps=1.d-14)

   sinomega=sin(omega)

   do ic=1,3
      rv(ic)=xyz(ic,l)-xyz(ic,i)
      rd(ic)=xyz(ic,k)-xyz(ic,j)
      re(ic)=xyz(ic,i)-xyz(ic,j)

      rdme(ic)=rd(ic)-re(ic)
   end do

   call crossprod(re,rd,rn)
   rvn=vecnorm(rv,3,0)
   rnn=vecnorm(rn,3,0)

   call crossprod(rv,re,rve)
   call crossprod(rn,re,rne)
   call crossprod(rd,rv,rdv)
   call crossprod(rd,rn,rdn)
   call crossprod(rv,rdme,rvdme)
   call crossprod(rn,rdme,rndme)

   nenner=rnn*rvn*cos(omega)
   if (abs(nenner).gt.eps) then
      onenner=1.d0/nenner
      do ic=1,3
         ! ... domega/dri
         domegadri(ic)=onenner*(rdv(ic)-rn(ic)-&
            &                       sinomega*(rvn/rnn*rdn(ic)-rnn/rvn*rv(ic)))

         ! ... domega/drj
         domegadrj(ic)=onenner*(rvdme(ic)-sinomega*rvn/rnn*rndme(ic))

         ! ... domega/drk
         domegadrk(ic)=onenner*(rve(ic)-sinomega*rvn/rnn*rne(ic))

         ! ... domega/drl
         domegadrl(ic)=onenner*(rn(ic)-sinomega*rnn/rvn*rv(ic))
      end do
   else
      do ic=1,3
         domegadri(ic)=0.d0
         domegadrj(ic)=0.d0
         domegadrk(ic)=0.d0
         domegadrl(ic)=0.d0
      end do
   end if

End subroutine

!  .....................................................................

Subroutine domegadrPBC&
      &           (nat,xyz,&
      &            i,j,k,l,vTr1,vTr2,vTr3,omega,&
      &            domegadri,domegadrj,domegadrk,domegadrl)
   use xtb_mctc_accuracy, only : wp

   !     inversion derivatives
   !  .....................................................................

   implicit none

   external vecnorm

   integer  ic,i,j,k,l,nat

   real(wp)&
      &         omega,sinomega,&
      &         vTr1(3), vTr2(3), vTr3(3), &
      &         xyz(3,nat),onenner,vecnorm,rnn,rvn,&
      &         rn(3),rv(3),rd(3),re(3),rdme(3),rve(3),&
      &         rne(3),rdv(3),rdn(3),&
      &         rvdme(3),rndme(3),nenner,&
      &         domegadri(3),domegadrj(3),domegadrk(3),domegadrl(3),eps

   parameter (eps=1.d-14)

   sinomega=sin(omega)

   do ic=1,3
      re(ic)= xyz(ic,i)-(xyz(ic,j)+vTr2(ic))            ! Vec central to 1st nb         
      rd(ic)=(xyz(ic,k)+vTr3(ic))-(xyz(ic,j)+vTr2(ic))  ! Vec 1st to 2nd nb
      rv(ic)=(xyz(ic,l)+vTr1(ic))-xyz(ic,i)             ! Vec central to 3rd nb 

      rdme(ic)=rd(ic)-re(ic)
   end do

   call crossprod(re,rd,rn)
   rvn=vecnorm(rv,3,0)
   rnn=vecnorm(rn,3,0)

   call crossprod(rv,re,rve)
   call crossprod(rn,re,rne)
   call crossprod(rd,rv,rdv)
   call crossprod(rd,rn,rdn)
   call crossprod(rv,rdme,rvdme)
   call crossprod(rn,rdme,rndme)

   nenner=rnn*rvn*cos(omega)
   if (abs(nenner).gt.eps) then
      onenner=1.d0/nenner
      do ic=1,3
         ! ... domega/dri
         domegadri(ic)=onenner*(rdv(ic)-rn(ic)-&
            &                       sinomega*(rvn/rnn*rdn(ic)-rnn/rvn*rv(ic)))

         ! ... domega/drj
         domegadrj(ic)=onenner*(rvdme(ic)-sinomega*rvn/rnn*rndme(ic))

         ! ... domega/drk
         domegadrk(ic)=onenner*(rve(ic)-sinomega*rvn/rnn*rne(ic))

         ! ... domega/drl
         domegadrl(ic)=onenner*(rn(ic)-sinomega*rnn/rvn*rv(ic))
      end do
   else
      do ic=1,3
         domegadri(ic)=0.d0
         domegadrj(ic)=0.d0
         domegadrk(ic)=0.d0
         domegadrl(ic)=0.d0
      end do
   end if

End subroutine

!  .....................................................................

real(wp) Function valijkl(nat,xyz,i,j,k,l)
   use xtb_mctc_accuracy, only : wp

   !  .....................................................................

   implicit none

   external vecnorm,valijk

   integer ic,i,j,k,l,nat

   real(wp) :: &
      &     xyz(3,nat),&
      &     eps,ra(3),rb(3),rc(3),na(3),nb(3),&
      &     rab,rbc,thab,thbc,valijk,&
      &     vecnorm,nan,nbn,rcn,snanb,deter,pi,test,&
      &     raabs,rbabs,rcabs

   parameter (eps=1.0d-14)
   data pi/3.1415926535897932384626433832795029d0/

   raabs=0.0d0
   rbabs=0.0d0
   rcabs=0.0d0
   ! ... get torsion coordinate
   do ic=1,3
      ra(ic)=xyz(ic,j)-xyz(ic,i)
      rb(ic)=xyz(ic,k)-xyz(ic,j)
      rc(ic)=xyz(ic,l)-xyz(ic,k)
      raabs=raabs+ra(ic)*ra(ic)
      rbabs=rbabs+rb(ic)*rb(ic)
      rcabs=rcabs+rc(ic)*rc(ic)
   end do
   raabs=sqrt(raabs)
   rbabs=sqrt(rbabs)
   rcabs=sqrt(rcabs)
   ! ... normalize vectors
   do ic=1,3
      ra(ic)=ra(ic)/raabs
      rb(ic)=rb(ic)/rbabs
      rc(ic)=rc(ic)/rcabs
   end do
   ! ... determinante of rb,ra,rc
   deter= ra(1)*(rb(2)*rc(3)-rb(3)*rc(2))&
      &      -ra(2)*(rb(1)*rc(3)-rb(3)*rc(1))&
      &      +ra(3)*(rb(1)*rc(2)-rb(2)*rc(1))

   call crossprod(-ra,rc,na)
   test=0.0d0
   do ic=1,3
      test=test+na(ic)*rb(ic)
   end do

   thab=valijk(nat,xyz,i,k,j)
   thbc=valijk(nat,xyz,j,l,k)
   call crossprod(ra,rb,na)
   call crossprod(rb,rc,nb)
   nan=vecnorm(na,3,1)
   nbn=vecnorm(nb,3,1)

   snanb=0.0d0
   do ic=1,3
      snanb=snanb+na(ic)*nb(ic)
   end do
   if (abs(abs(snanb)-1.d0).lt.eps) then
      snanb=sign(1.d0,snanb)
   end if

   valijkl=acos(snanb)
   if(test.lt.0.0d0) valijkl=-acos(snanb)

   ! the gradient dphir is only compatible with this subroutine
   ! if the statement below is commented out. If not, opt. and
   ! Hessian show large errors and imags. I don't understand
   ! this entirely but thats how it is.
   ! SG, Sat May 2411:41:42 CEST 2014

   !     if (deter.lt.0) then
   !        valijkl=2.d0*pi-valijkl
   !     end if

End function

!  .....................................................................

real(wp) Function valijk(nat,xyz,j,k,i)
   use xtb_mctc_accuracy, only : wp

   !  .....................................................................

   implicit none

   external vecnorm

   integer nat,j,k,i,ic

   real(wp) :: &
      &         ra(3),rb(3),rab,eps,&
      &         xyz(3,nat),vecnorm,ran,rbn

   parameter (eps=1.d-14)

   do ic=1,3
      ra(ic)=xyz(ic,j)-xyz(ic,i)
      rb(ic)=xyz(ic,k)-xyz(ic,i)
   end do

   ran=vecnorm(ra,3,1)
   rbn=vecnorm(rb,3,1)
   rab=0.d0
   do ic=1,3
      rab=rab+ra(ic)*rb(ic)
   end do

   if (abs(abs(rab)-1.d0).lt.eps) then
      rab=sign(1.d0,rab)
   end if
   valijk=acos(rab)

End function


!  .....................................................................

real(wp) Function valijkPBC(mode,nat,xyz,j,k,i,vTr1,vTr2,vTr3)
   use xtb_mctc_accuracy, only : wp

   !  .....................................................................

   implicit none

   external vecnorm

   integer mode,nat,j,k,i,ic

   real(wp) :: &
      &         ra(3),rb(3),rab,eps,&
      &         xyz(3,nat),vecnorm,ran,rbn,vTr1(3),vTr2(3),vTr3(3)

   parameter (eps=1.d-14)

   if (mode.eq.1) then ! here j=l,k=j,i=i are inserted, vTr1=vTrl, vTr2=vTrj
     do ic=1,3
       ra(ic)=(xyz(ic,j)+vTr1(ic))-xyz(ic,i)
       rb(ic)=(xyz(ic,k)+vTr2(ic))-xyz(ic,i)
     end do
   elseif(mode.eq.2) then ! here j=i,k=k,i=j are inserted vTr1=vTrk, vTr2=vTrj
     do ic=1,3
       ra(ic)= xyz(ic,j)-(xyz(ic,i)+vTr2(ic))
       rb(ic)=(xyz(ic,k)+vTr1(ic))-(xyz(ic,i)+vTr2(ic))
     end do
   elseif(mode.eq.3) then ! here j=R k=C i=B vTr1=vTrR vTr2=vTrC vTr3=vTrB
     do ic=1,3
       ra(ic)= (xyz(ic,j)+vTr1(ic))-(xyz(ic,i)+vTr3(ic)) ! R - B
       rb(ic)= (xyz(ic,k)+vTr2(ic))-(xyz(ic,i)+vTr3(ic)) ! C - B
     end do
   elseif(mode.eq.4) then ! here j=B k=H i=C vTr1=vTrB vTr2=vTrC
     do ic=1,3
       ra(ic)=(xyz(ic,j)+vTr1(ic))-(xyz(ic,i)+vTr2(ic))
       rb(ic)=(xyz(ic,k))         -(xyz(ic,i)+vTr2(ic))
     end do
   endif

   ran=vecnorm(ra,3,1)
   rbn=vecnorm(rb,3,1)
   rab=0.d0
   do ic=1,3
      rab=rab+ra(ic)*rb(ic)
   end do

   if (abs(abs(rab)-1.d0).lt.eps) then
      rab=sign(1.d0,rab)
   end if
   valijkPBC=acos(rab)

End function

!  .....................................................................

real(wp) Function omega (nat,xyz,i,j,k,l)
   use xtb_mctc_accuracy, only : wp

   !   Calculates the inversion angle
   !  .....................................................................

   implicit none

   external vecnorm

   integer ic,nat,i,j,k,l

   real(wp) :: &
      &        xyz(3,nat),&
      &        rd(3),re(3),rn(3),rv(3),rnv,&
      &        vecnorm,rkjn,rljn,rnn,rvn

   do ic=1,3
      re(ic)=xyz(ic,i)-xyz(ic,j)
      rd(ic)=xyz(ic,k)-xyz(ic,j)
      rv(ic)=xyz(ic,l)-xyz(ic,i)
   end do
   call crossprod(re,rd,rn)
   rnn=vecnorm(rn,3,1)
   rvn=vecnorm(rv,3,1)

   rnv=rn(1)*rv(1)+rn(2)*rv(2)+rn(3)*rv(3)
   omega=asin( rnv )

End function


real(wp) Function omegaPBC (nat,xyz,i,j,k,l,vTr1,vTr2,vTr3)
   use xtb_mctc_accuracy, only : wp

   !   Calculates the inversion angle
   !  .....................................................................

   implicit none

   external vecnorm

   integer ic,nat,i,j,k,l  
   
   real(wp) :: &
      &        xyz(3,nat),vTr1(3),vTr2(3),vTr3(3),&
      &        rd(3),re(3),rn(3),rv(3),rnv,&
      &        vecnorm,rkjn,rljn,rnn,rvn 
   ! out-of-plane case from ini; atoms and iTr's sorted by distance to atom i
   ! i=central, j=1st nb, k=2nd, l=3rd
   do ic=1,3
      re(ic)=xyz(ic,i)-(xyz(ic,j)+vTr2(ic))                ! Vec central to 1st nb         
      rd(ic)=(xyz(ic,k)+vTr3(ic))-(xyz(ic,j)+vTr2(ic))            ! Vec 1st to 2nd nb
      rv(ic)=(xyz(ic,l)+vTr1(ic))-xyz(ic,i) ! Vec central to 3rd nb 
   end do
   call crossprod(re,rd,rn)
   rnn=vecnorm(rn,3,1)
   rvn=vecnorm(rv,3,1)

   rnv=rn(1)*rv(1)+rn(2)*rv(2)+rn(3)*rv(3)
   omegaPBC=asin( rnv )

End function


!  .....................................................................

Subroutine crossprod(ra,rb,rab)
   use xtb_mctc_accuracy, only : wp

   implicit none

   real(wp) ra(3),rb(3),rab(3)

   rab(1)=ra(2)*rb(3)-ra(3)*rb(2)
   rab(2)=ra(3)*rb(1)-ra(1)*rb(3)
   rab(3)=ra(1)*rb(2)-ra(2)*rb(1)

End subroutine

!  .....................................................................

real(wp) Function  vecnorm (r,n,inorm)
   use xtb_mctc_accuracy, only : wp

   implicit none

   integer i,n,inorm

   real(wp) r(n),or,sp,rn

   sp=0.
   do i=1,n
      sp=sp+r(i)*r(i)
   end do
   rn=sqrt(sp)
   if (inorm.gt.0) then
      if (abs(rn).gt.1.d-14) then
         or=1.0d0/rn
         do i=1,n
            r(i)=or*r(i)
         end do
      end if
   end if
   vecnorm=rn

End function

!  .....................................................................

subroutine crprod(a,b,c)
   use xtb_mctc_accuracy, only : wp
   implicit double precision (a-h,o-z)
   dimension a(3),b(3),c(3)

   c(1)=a(2)*b(3)-a(3)*b(2)
   c(2)=a(3)*b(1)-a(1)*b(3)
   c(3)=a(1)*b(2)-a(2)*b(1)

   return
end subroutine crprod

!  .....................................................................

subroutine vsc1(a,scale,tol)
   use xtb_mctc_accuracy, only : wp
   implicit integer (i-n)
   implicit double precision (a-h,o-z)
   dimension a(3)

   rlen = norm2(a)

   if (rlen.gt.tol) then
      do i=1,3
         a(i) = a(i)*scale/rlen
      end do
   endif

   return
end subroutine vsc1

!  .....................................................................

pure subroutine impsc(a,b,c)
   use xtb_mctc_accuracy, only : wp
   implicit none
   real(wp),intent(in)  :: a(3),b(3)
   real(wp),intent(out) :: c
   integer  :: i
   real(wp) :: rimp,al,bl

   rimp = 0.0d0

   do i=1,3
      rimp = rimp + a(i)*b(i)
   end do

   al = norm2(a)
   bl = norm2(b)

   if (al.gt.0.0d0.and.bl.gt.0.0d0) then
      c = rimp/(al*bl)
   else
      c = 0.0d0
   endif

   return
end subroutine impsc

!  .....................................................................

pure subroutine bangl(xyz,i,j,k,angle)
   use xtb_mctc_accuracy, only : wp
   implicit none
   real(wp),intent(in)  :: xyz(3,*)
   real(wp),intent(out) :: angle
   integer, intent(in)  :: i,j,k

   real(wp) d2ij,d2jk,d2ik,xy,temp

   d2ij = sum((xyz(:,i)-xyz(:,j))**2)
   d2jk = sum((xyz(:,j)-xyz(:,k))**2)
   d2ik = sum((xyz(:,i)-xyz(:,k))**2)
   xy = sqrt(d2ij*d2jk)
   temp = 0.5d0 * (d2ij+d2jk-d2ik) / xy
   if (temp .gt. 1.0d0)  temp= 1.0d0
   if (temp .lt. -1.0d0) temp=-1.0d0
   angle = acos( temp )

end subroutine bangl

pure subroutine banglPBC(mode,xyz,i,j,k,iTr,iTr2,neigh,angle)
   use xtb_mctc_accuracy, only : wp
   use xtb_gfnff_neighbor
   implicit none
   real(wp),intent(in)  :: xyz(3,*)                  
   integer, intent(in)  :: mode,i,j,k,iTr,iTr2  ! j is in the middle
   type(TNeigh), intent(in) :: neigh
   real(wp),intent(out) :: angle   
                                  
   real(wp) d2ij,d2jk,d2ik,xy,temp,trV(3),trV2(3)
   trV =neigh%transVec(:,iTr)
   trV2=neigh%transVec(:,iTr2)
   if (mode.eq.1) then                                  
     d2ij = sum(((xyz(:,i)+trV)-xyz(:,j))**2)
     d2jk = sum((xyz(:,j)-(xyz(:,k)+trV2))**2)
     d2ik = sum(((xyz(:,i)+trV)-(xyz(:,k)+trV2))**2)  
   endif
   if (mode.eq.2) then
     d2ij = sum((xyz(:,i)-(xyz(:,j)+trV))**2)
     d2jk = sum(((xyz(:,j)+trV)-(xyz(:,k)+trV2))**2)
     d2ik = sum((xyz(:,i)-(xyz(:,k)+trV2))**2)  
   endif
   xy = sqrt(d2ij*d2jk)                
   temp = 0.5d0 * (d2ij+d2jk-d2ik) / xy  ! the angle is between side dij and djk
   if (temp .gt. 1.0d0)  temp= 1.0d0
   if (temp .lt. -1.0d0) temp=-1.0d0   
   angle = acos( temp )             
                                                     
end subroutine banglPBC
