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

module xtb_grad_core
   use xtb_mctc_accuracy, only : wp
   use xtb_xtb_data

   use xtb_mctc_la

!! ========================================================================
!  we use the complete xtb_scc_core, so we inherit all interfaces and all
!  parameters. Therefore, we don't need to declare them again, note that
!  by including the xtb_grad_core you also inherit the xtb_scc_core!
   use xtb_scc_core

   implicit none

   integer, private, parameter :: mmm(*)=(/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/)

contains

!! ========================================================================
!  convert H0/H1 from eV to Eh and calculate the energy weighted density
!  matrix
!! ========================================================================
subroutine prep_grad_conv(ndim,H0,H1,C,focc,emo,X)
   use xtb_mctc_convert, only : autoev,evtoau
   integer, intent(in)    :: ndim
   real(wp),intent(inout) :: H0(ndim*(ndim+1)/2)
   real(wp),intent(inout) :: H1(ndim*(ndim+1)/2)
   real(wp),intent(in)    :: C(ndim,ndim)
   real(wp),intent(in)    :: focc(ndim)
   real(wp),intent(in)    :: emo(ndim)
   real(wp),intent(out)   :: X(ndim,ndim)

   real(wp),allocatable :: temp(:)

   allocate( temp(ndim), source = 0.0_wp )
   H0 = H0*evtoau
   H1 = H1*evtoau

!  get energy weigthed density matrix
   temp = focc * emo*evtoau

   call dmat(ndim,temp,C,X)

   deallocate( temp )

end subroutine prep_grad_conv

!! ========================================================================
!  wave function terms
!! ========================================================================
subroutine poly_grad(hData,g,n,at,ndim,nmat2,matlist2,xyz,sqrab,P,S,aoat2,lao2,H0)
   use xtb_lin
   type(THamiltonianData), intent(in) :: hData
   real(wp),intent(inout) :: g(3,n)
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   integer, intent(in)    :: ndim
   integer, intent(in)    :: nmat2
   integer, intent(in)    :: matlist2(2,nmat2)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(in)    :: sqrab(n*(n+1)/2)
   real(wp),intent(in)    :: P(ndim,ndim)
   real(wp),intent(in)    :: S(ndim,ndim)
   integer, intent(in)    :: aoat2(ndim)
   integer, intent(in)    :: lao2(ndim)
   real(wp),intent(in)    :: H0(ndim*(ndim+1)/2)

   integer  :: i,j,kk,m
   integer  :: iat,jat,ishell,jshell,iZp,jZp
   real(wp) :: hji,rab2
   real(wp) :: h0s,h0sr
   real(wp) :: dum
   real(wp) :: drfdxyz(3)

!  using matlist2 => overlap between all atom pairs
   do m = 1, nmat2
      i = matlist2(1,m)
      j = matlist2(2,m)
      kk = j+i*(i-1)/2
      iat = aoat2(i)
      jat = aoat2(j)
      iZp = at(iat)
      jZp = at(jat)
      ishell=mmm(lao2(i))
      jshell=mmm(lao2(j))
      rab2 = sqrab(lin(jat,iat))
      hji = 2.0_wp*P(j,i)*S(j,i)
!     H0=S H (1+b x)
!     dH0/dx = H(1+bx)dS/dx + S H (b)
      call dshellPoly(hData%shellPoly(iShell,iZp),hData%shellPoly(jShell,jZp),&
         & hData%atomicRad(iZp),hData%atomicRad(jZp),rab2,xyz(:,iat),xyz(:,jat),&
         & dum,drfdxyz)
      h0s = H0(kk)/S(j,i) !H(1+bx)
      h0sr = h0s/dum      !H
      g(1,iat) = g(1,iat) + Hji*h0sr*drfdxyz(1)
      g(2,iat) = g(2,iat) + Hji*h0sr*drfdxyz(2)
      g(3,iat) = g(3,iat) + Hji*h0sr*drfdxyz(3)
      g(1,jat) = g(1,jat) - Hji*h0sr*drfdxyz(1)
      g(2,jat) = g(2,jat) - Hji*h0sr*drfdxyz(2)
      g(3,jat) = g(3,jat) - Hji*h0sr*drfdxyz(3)
   enddo

end subroutine poly_grad

!! ========================================================================
!  CN dependent part of GFN1 hamiltonian
!! ========================================================================
subroutine hcn_grad_gfn1(hData,g,n,at,ndim,nmat2,matlist2,xyz, &
   &                     kspd,kmagic,kenscal,kcnao,P,S,dcn, &
   &                     aoat2,lao2,valao2,hdiag2)
   use xtb_mctc_convert, only : autoev,evtoau
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   integer, intent(in)    :: ndim
   integer, intent(in)    :: nmat2
   integer, intent(in)    :: matlist2(2,nmat2)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(in)    :: kspd(6)
   real(wp),intent(in)    :: kmagic(4,4)
   real(wp),intent(in)    :: kenscal
   real(wp),intent(in)    :: kcnao(ndim)
   real(wp),intent(in)    :: P(ndim,ndim)
   real(wp),intent(in)    :: S(ndim,ndim)
   real(wp),intent(in)    :: dcn(3,n,n)
   integer, intent(in)    :: aoat2(ndim)
   integer, intent(in)    :: lao2(ndim)
   integer, intent(in)    :: valao2(ndim)
   real(wp),intent(in)    :: hdiag2(ndim)
   real(wp),intent(inout) :: g(3,n)

   real(wp),allocatable :: hcn(:)

   integer  :: i,j,m
   integer  :: iat,jat,iZp,jZp
   integer  :: ishell,jshell
   real(wp) :: hji
   real(wp) :: gtmp(3)
   real(wp) :: km
   real(wp) :: dum,dum1,dum2

   allocate( hcn(n), source = 0.0_wp )

!  CN dependent part of H0
   hcn=0.0_wp
!$omp parallel default(none) &
!$omp shared(nmat2,matlist2,aoat2,lao2,valao2,P,S,n,at,xyz,hdiag2,hData) &
!$omp shared(kspd,kmagic,kenscal,kcnao) &
!$omp private(i,m,j,jat,iat,dum1,dum2,km,dum,ishell,jshell,hji,iZp,jZp) &
!$omp reduction (+:hcn)
!$omp do
   do m=1,nmat2
      i=matlist2(1,m)
      j=matlist2(2,m)
      iat=aoat2(i)
      jat=aoat2(j)
      iZp = at(iat)
      jZp = at(jat)
      ishell=mmm(lao2(i))
      jshell=mmm(lao2(j))
      hji=P(j,i)*S(j,i)
      jat=aoat2(j)
      dum = shellPoly(hData%shellPoly(iShell, iZp), hData%shellPoly(jShell, jZp), &
         & hData%atomicRad(iZp), hData%atomicRad(jZp),xyz(:,iat),xyz(:,jat))
      call h0scal(n,at,i,j,ishell,jshell,iat,jat,valao2(i).ne.0,valao2(j).ne.0,  &
      &           kspd,kmagic,kenscal,km)
      dum1=hji*km*dum*hdiag2(i)*kcnao(i)*evtoau ! h independent part in H0
      dum2=hji*km*dum*hdiag2(j)*kcnao(j)*evtoau ! h independent part in H0
      hcn(jat)=hcn(jat)+dum2
      hcn(iat)=hcn(iat)+dum1
   enddo
!$omp end do
!$omp end parallel
!$omp parallel default(none) &
!$omp shared(ndim,aoat2,P,hdiag2,kcnao) &
!$omp private(i,iat,dum1) &
!$omp reduction (+:hcn)
!$omp do
   do i=1,ndim
      iat=aoat2(i)
      dum1=P(i,i)*hdiag2(i)*kcnao(i)*evtoau ! diagonal contribution
      hcn(iat)=hcn(iat)+dum1
   enddo
!$omp end do
!$omp end parallel

!  CN-level shift gradient
!$omp parallel default(none) &
!$omp shared(n,dcn,hcn) &
!$omp private(i,gtmp) shared ( g )
!$omp do
   do i=1,n
      call gemv('n',3,n,1.0d0,dcn(:,:,i),3,hcn,1,0.0d0,gtmp,1)
      g(1:3,i)=g(1:3,i)+gtmp(1:3)
   enddo
!$omp end do
!$omp end parallel

end subroutine hcn_grad_gfn1

!! ========================================================================
!  CN dependent part of GFN2 hamiltonian
!! ========================================================================
subroutine hcn_grad_gfn2(hData,g,n,at,ndim,nmat2,matlist2,xyz, &
   &                     kspd,kmagic,kenscal,kcnao,P,S,dcn, &
   &                     aoat2,lao2,valao2,hdiag2,aoexp)
   use xtb_mctc_convert, only : autoev,evtoau
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   integer, intent(in)    :: ndim
   integer, intent(in)    :: nmat2
   integer, intent(in)    :: matlist2(2,nmat2)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(in)    :: kspd(6)
   real(wp),intent(in)    :: kmagic(4,4)
   real(wp),intent(in)    :: kenscal
   real(wp),intent(in)    :: kcnao(ndim)
   real(wp),intent(in)    :: P(ndim,ndim)
   real(wp),intent(in)    :: S(ndim,ndim)
   real(wp),intent(in)    :: dcn(3,n,n)
   integer, intent(in)    :: aoat2(ndim)
   integer, intent(in)    :: lao2(ndim)
   integer, intent(in)    :: valao2(ndim)
   real(wp),intent(in)    :: hdiag2(ndim)
   real(wp),intent(in)    :: aoexp(ndim)
   real(wp),intent(inout) :: g(3,n)

   real(wp),allocatable :: hcn(:)

   integer  :: i,j,m
   integer  :: iat,jat,iZp,jZp
   integer  :: ishell,jshell
   real(wp) :: hji
   real(wp) :: gtmp(3)
   real(wp) :: km,fact
   real(wp) :: dum,dum1,dum2
   real(wp),parameter :: aot = -0.5_wp ! AO exponent dep. H0 scal

   allocate( hcn(n), source = 0.0_wp )

!  CN dependent part of H0
   hcn=0.0_wp
!$omp parallel default(none) &
!$omp shared(nmat2,matlist2,aoat2,lao2,valao2,P,S,n,at,xyz,hData) &
!$omp shared(kspd,kmagic,kenscal,aoexp,kcnao) &
!$omp private(i,m,j,jat,iat,dum1,dum2,km,dum,ishell,jshell,hji,fact,iZp,jZp) &
!$omp reduction (+:hcn)
!$omp do
   do m=1,nmat2
      i=matlist2(1,m)
      j=matlist2(2,m)
      iat=aoat2(i)
      jat=aoat2(j)
      iZp = at(iat)
      jZp = at(jat)
      ishell=mmm(lao2(i))
      jshell=mmm(lao2(j))
      hji=P(j,i)*S(j,i)
      dum = shellPoly(hData%shellPoly(iShell, iZp), hData%shellPoly(jShell, jZp), &
         & hData%atomicRad(iZp), hData%atomicRad(jZp),xyz(:,iat),xyz(:,jat))
      call h0scal(n,at,i,j,ishell,jshell,iat,jat,valao2(i).ne.0,valao2(j).ne.0,  &
      &               kspd,kmagic,kenscal,km)
      fact = 0.5_wp*(aoexp(i)+aoexp(j))/sqrt(aoexp(i)*aoexp(j))
      km = km*fact**aot
      dum1=hji*km*dum*kcnao(i)*evtoau ! h independent part in H0
      dum2=hji*km*dum*kcnao(j)*evtoau ! h independent part in H0
      hcn(jat)=hcn(jat)-dum2
      hcn(iat)=hcn(iat)-dum1
   enddo
!$omp end do
!$omp end parallel
!$omp parallel default(none) &
!$omp shared(ndim,aoat2,P,kcnao) &
!$omp private(i,iat,dum1) &
!$omp reduction (+:hcn)
!$omp do
   do i=1,ndim
      iat=aoat2(i)
      dum1=P(i,i)*kcnao(i)*evtoau ! diagonal contribution
      hcn(iat)=hcn(iat)-dum1
   enddo
!$omp end do
!$omp end parallel


!  CN-level shift gradient
!$omp parallel default(none) &
!$omp shared(n,dcn,hcn) &
!$omp private(i,gtmp) shared ( g )
!$omp do
   do i=1,n
      call gemv('n',3,n,1.0d0,dcn(:,:,i),3,hcn,1,0.0d0,gtmp,1)
      g(1:3,i)=g(1:3,i)+gtmp(1:3)
   enddo
!$omp end do
!$omp end parallel

end subroutine hcn_grad_gfn2

!! ========================================================================
!  derivative of the CM5 additional term for GBSA in GFN1
!! ========================================================================
subroutine cm5_grad_gfn1(g,n,q,fgb,fhb,dcm5a,lhb)
   use xtb_mctc_convert, only : autoev,evtoau
   real(wp),intent(inout) :: g(3,n)
   integer, intent(in)    :: n
   real(wp),intent(in)    :: q(n)
   real(wp),intent(inout) :: fgb(n,n)
   real(wp),intent(inout) :: fhb(n)
   real(wp),intent(in)    :: dcm5a(3,n,n)
   logical, intent(in)    :: lhb

   real(wp), allocatable  :: fgba(:)
   real(wp), allocatable  :: dcm5(:,:)
   integer :: iat,jat

   allocate( fgba(n), source = 0.0_wp )
   call dsymv('u',n,evtoau,fgb,n,q,1,0.0_wp,fgba,1)
   call dgemv('n',3*n,n,1.0_wp,dcm5a,3*n,fgba,1,1.0_wp,g,1)
   if (lhb) then
      fgba = evtoau * fhb * 2.0_wp * q
      call dgemv('n',3*n,n,1.0_wp,dcm5a,3*n,fgba,1,1.0_wp,g,1)
   end if

end subroutine cm5_grad_gfn1


!! ========================================================================
!  shellwise electrostatic gradient for GFN1
!! ========================================================================
subroutine shelles_grad_gfn1(g,jData,n,at,nshell,xyz,sqrab,ash,lsh,alphaj,qsh)
   use xtb_lin
   type(TCoulombData), intent(in) :: jData
   real(wp),intent(inout) :: g(3,n)
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   integer, intent(in) :: nshell
   real(wp),intent(in) :: xyz(3,n)
   real(wp),intent(in) :: sqrab(n*(n+1)/2)
   integer, intent(in) :: ash(nshell)
   integer, intent(in) :: lsh(nshell)
   real(wp),intent(in) :: alphaj
   real(wp),intent(in) :: qsh(nshell)

   integer  :: is,js,iat,jat,ati,atj
   real(wp) :: xa,ya,za,dx,dy,dz
   real(wp) :: gi,gj,r2,rr,yy,ff

   do is=1,nshell
      iat=ash(is)
      xa=xyz(1,iat)
      ya=xyz(2,iat)
      za=xyz(3,iat)
      ati=at(iat)
      gi=jData%chemicalHardness(ati)*(1.0d0+jData%shellHardness(lsh(is)+1,ati))
      do js=1,nshell
         jat=ash(js)
         if(jat.le.iat) cycle
         dx=xa-xyz(1,jat)
         dy=ya-xyz(2,jat)
         dz=za-xyz(3,jat)
         atj=at(jat)
         r2=sqrab(lin(jat,iat))
         gj=jData%chemicalHardness(atj)*(1.0d0+jData%shellHardness(lsh(js)+1,atj))
         rr=2.0d0/(1./gi+1./gj)
         rr=1.0d0/rr**alphaj
         ff=r2**(alphaj/2.0d0-1.0d0)* &
         &        (r2**(alphaj*0.5d0)+rr)**(-1.0d0/alphaj-1.0d0)
         yy=ff*qsh(is)*qsh(js)
         g(1,iat)=g(1,iat)-dx*yy
         g(2,iat)=g(2,iat)-dy*yy
         g(3,iat)=g(3,iat)-dz*yy
         g(1,jat)=g(1,jat)+dx*yy
         g(2,jat)=g(2,jat)+dy*yy
         g(3,jat)=g(3,jat)+dz*yy
      enddo
   enddo

end subroutine shelles_grad_gfn1

!! ========================================================================
!  shellwise electrostatic gradient for GFN2
!! ========================================================================
subroutine shelles_grad_gfn2(g,jData,n,at,nshell,xyz,sqrab,ash,lsh,qsh)
   use xtb_lin
   type(TCoulombData), intent(in) :: jData
   real(wp),intent(inout) :: g(3,n)
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   integer, intent(in) :: nshell
   real(wp),intent(in) :: xyz(3,n)
   real(wp),intent(in) :: sqrab(n*(n+1)/2)
   integer, intent(in) :: ash(nshell)
   integer, intent(in) :: lsh(nshell)
   real(wp),intent(in) :: qsh(nshell)

   integer  :: is,js,iat,jat,ati,atj
   real(wp) :: xa,ya,za,dx,dy,dz
   real(wp) :: gi,gj,r2,rr,yy

   do is=1,nshell
      iat=ash(is)
      xa=xyz(1,iat)
      ya=xyz(2,iat)
      za=xyz(3,iat)
      ati=at(iat)
      gi=jData%chemicalHardness(ati)*(1.0d0+jData%shellHardness(1+lsh(is),ati))
      do js=1,nshell
         jat=ash(js)
         if(jat.le.iat) cycle
         dx=xa-xyz(1,jat)
         dy=ya-xyz(2,jat)
         dz=za-xyz(3,jat)
         atj=at(jat)
         r2=sqrab(lin(jat,iat))
         gj=jData%chemicalHardness(atj)*(1.0d0+jData%shellHardness(1+lsh(js),atj))
         rr=0.5d0*(gi+gj)
         rr=1.0d0/rr**2
!        rr=1.0d0/(gi*gj) !NEWAV
         yy=qsh(is)*qsh(js)*(r2+rr)**(-1.5d0)
         g(1,iat)=g(1,iat)-dx*yy
         g(2,iat)=g(2,iat)-dy*yy
         g(3,iat)=g(3,iat)-dz*yy
         g(1,jat)=g(1,jat)+dx*yy
         g(2,jat)=g(2,jat)+dy*yy
         g(3,jat)=g(3,jat)+dz*yy
      enddo
   enddo

end subroutine shelles_grad_gfn2


!! ========================================================================
!  derivative of S(R) enhancement factor
!! ========================================================================
pure subroutine dshellPoly(iPoly,jPoly,iRad,jRad,rab2,xyz1,xyz2,rf,dxyz)
   real(wp), intent(in)  :: iPoly,jPoly
   real(wp), intent(in)  :: iRad,jRad
   real(wp), intent(in)  :: rab2
   real(wp), intent(out) :: dxyz(3),rf
   real(wp), intent(in)  :: xyz1(3),xyz2(3)
   real(wp) :: rab,k1,k2,rr,r,a,dum,rf1,rf2,dx,dy,dz
   real(wp) :: t14,t15,t17,t22,t20,t23,t10,t11,t35,t13

   a=0.5            ! R^a dependence 0.5 in GFN1

   dx=xyz1(1)-xyz2(1)
   dy=xyz1(2)-xyz2(2)
   dz=xyz1(3)-xyz2(3)

   rab=sqrt(rab2)

   ! this sloppy conv. factor has been used in development, keep it
   r=iRad+jRad

   rr=rab/r

   k1=iPoly*0.01_wp
   k2=jPoly*0.01_wp

   t14 = rr**a
   t15 = k1*t14
   t17 = 1.0_wp/rab2
   t22 = rr**a
   t23 = k2*t22
   rf=(1.0_wp+t15)*(1.0_wp+k2*t22)
   dxyz(:)=(t15*(1.0_wp+t23)+(1.0_wp+t15)*k2*t22)*a*t17*[dx,dy,dz]

end subroutine dshellPoly

end module xtb_grad_core
