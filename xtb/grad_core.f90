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

module grad_core
   use iso_fortran_env, only : wp => real64

   use mctc_la

!! ========================================================================
!  we use the complete scc_core, so we inherit all interfaces and all
!  parameters. Therefore, we don't need to declare them again, note that
!  by including the grad_core you also inherit the scc_core!
   use scc_core

   implicit none


contains

!! ========================================================================
!  convert H0/H1 from eV to Eh and calculate the energy weighted density
!  matrix
!! ========================================================================
subroutine prep_grad_conv(ndim,H0,H1,C,focc,emo,X)
   use mctc_econv, only : autoev,evtoau
   implicit none
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
subroutine poly_grad(g,n,at,ndim,nmat2,matlist2,xyz,sqrab,P,S,aoat2,lao2,H0)
   implicit none
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

   integer,external :: lin

   integer  :: i,j,kk,m
   integer  :: ati,atj,ishell,jshell
   real(wp) :: hji,rab2
   real(wp) :: h0s,h0sr
   real(wp) :: dum
   real(wp) :: drfdxyz(3)

!  using matlist2 => overlap between all atom pairs
   do m = 1, nmat2
      i = matlist2(1,m)
      j = matlist2(2,m)
      kk = j+i*(i-1)/2
      ati = aoat2(i)
      atj = aoat2(j)
      ishell=mmm(lao2(i))
      jshell=mmm(lao2(j))
      rab2 = sqrab(lin(atj,ati))
      hji = 2.0_wp*P(j,i)*S(j,i)
!     H0=S H (1+b x)
!     dH0/dx = H(1+bx)dS/dx + S H (b)
      call dhdr(n,at,xyz,ati,i,j,ati,atj,ishell,jshell,rab2,dum,drfdxyz)
      h0s = H0(kk)/S(j,i) !H(1+bx)
      h0sr = h0s/dum      !H
      g(1,ati) = g(1,ati) + Hji*h0sr*drfdxyz(1)
      g(2,ati) = g(2,ati) + Hji*h0sr*drfdxyz(2)
      g(3,ati) = g(3,ati) + Hji*h0sr*drfdxyz(3)
      g(1,atj) = g(1,atj) - Hji*h0sr*drfdxyz(1)
      g(2,atj) = g(2,atj) - Hji*h0sr*drfdxyz(2)
      g(3,atj) = g(3,atj) - Hji*h0sr*drfdxyz(3)
   enddo

end subroutine poly_grad

!! ========================================================================
!  CN dependent part of GFN1 hamiltonian
!! ========================================================================
subroutine hcn_grad_gfn1(g,n,at,ndim,nmat2,matlist2,xyz, &
   &                     kspd,kmagic,kenscal,kcnao,P,S,dcn, &
   &                     aoat2,lao2,valao2,hdiag2)
   use mctc_econv, only : autoev,evtoau
   use aoparam
   implicit none
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
   integer  :: ati,atj
   integer  :: ishell,jshell
   real(wp) :: hji
   real(wp) :: gtmp(3)
   real(wp) :: km
   real(wp) :: dum,dum1,dum2

   allocate( hcn(n), source = 0.0_wp )

!  CN dependent part of H0
   hcn=0.0_wp
!$omp parallel default(none) &
!$omp shared(nmat2,matlist2,aoat2,lao2,valao2,P,S,n,at,xyz,hdiag2) &
!$omp shared(kspd,kmagic,kenscal,kcnao) &
!$omp private(i,m,j,atj,ati,dum1,dum2,km,dum,ishell,jshell,hji) &
!$omp reduction (+:hcn)
!$omp do
   do m=1,nmat2
      i=matlist2(1,m)
      j=matlist2(2,m)
      ati=aoat2(i)
      atj=aoat2(j)
      ishell=mmm(lao2(i))
      jshell=mmm(lao2(j))
      hji=P(j,i)*S(j,i)
      atj=aoat2(j)
      dum=rfactor(ishell,jshell,at(ati),at(atj),  &
      &           xyz(:,ati),xyz(:,atj))
      call h0scal(n,at,i,j,ishell,jshell,ati,atj,valao2(i).ne.0,valao2(j).ne.0,  &
      &           kspd,kmagic,kenscal,km)
      dum1=hji*km*dum*hdiag2(i)*kcnao(i)*evtoau ! h independent part in H0
      dum2=hji*km*dum*hdiag2(j)*kcnao(j)*evtoau ! h independent part in H0
      hcn(atj)=hcn(atj)+dum2
      hcn(ati)=hcn(ati)+dum1
   enddo
!$omp end do
!$omp end parallel
!$omp parallel default(none) &
!$omp shared(ndim,aoat2,P,hdiag2,kcnao) &
!$omp private(i,ati,dum1) &
!$omp reduction (+:hcn)
!$omp do
   do i=1,ndim
      ati=aoat2(i)
      dum1=P(i,i)*hdiag2(i)*kcnao(i)*evtoau ! diagonal contribution
      hcn(ati)=hcn(ati)+dum1
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
subroutine hcn_grad_gfn2(g,n,at,ndim,nmat2,matlist2,xyz, &
   &                     kspd,kmagic,kenscal,kcnao,P,S,dcn, &
   &                     aoat2,lao2,valao2,hdiag2,aoexp)
   use mctc_econv, only : autoev,evtoau
   use aoparam
   implicit none
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
   integer  :: ati,atj
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
!$omp shared(nmat2,matlist2,aoat2,lao2,valao2,P,S,n,at,xyz) &
!$omp shared(kspd,kmagic,kenscal,aoexp,kcnao) &
!$omp private(i,m,j,atj,ati,dum1,dum2,km,dum,ishell,jshell,hji,fact) &
!$omp reduction (+:hcn)
!$omp do
   do m=1,nmat2
      i=matlist2(1,m)
      j=matlist2(2,m)
      ati=aoat2(i)
      atj=aoat2(j)
      ishell=mmm(lao2(i))
      jshell=mmm(lao2(j))
      hji=P(j,i)*S(j,i)
      dum=rfactor(ishell,jshell,at(ati),at(atj),  &
      &            xyz(:,ati),xyz(:,atj))
      call h0scal(n,at,i,j,ishell,jshell,ati,atj,valao2(i).ne.0,valao2(j).ne.0,  &
      &               kspd,kmagic,kenscal,km)
      fact = 0.5_wp*(aoexp(i)+aoexp(j))/sqrt(aoexp(i)*aoexp(j))
      km = km*fact**aot
      dum1=hji*km*dum*kcnao(i)*evtoau ! h independent part in H0
      dum2=hji*km*dum*kcnao(j)*evtoau ! h independent part in H0
      hcn(atj)=hcn(atj)-dum2
      hcn(ati)=hcn(ati)-dum1
   enddo
!$omp end do
!$omp end parallel
!$omp parallel default(none) &
!$omp shared(ndim,aoat2,P,kcnao) &
!$omp private(i,ati,dum1) &
!$omp reduction (+:hcn)
!$omp do
   do i=1,ndim
      ati=aoat2(i)
      dum1=P(i,i)*kcnao(i)*evtoau ! diagonal contribution
      hcn(ati)=hcn(ati)-dum1
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
   use mctc_econv, only : autoev,evtoau
   implicit none
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

   allocate( fgba(n),dcm5(3,n), source = 0.0_wp )

   fgb = fgb * evtoau
   fhb = fhb * evtoau

   do iat = 1, n
      do jat = 1, n
         fgba(iat)=fgba(iat)+q(jat)*fgb(jat,iat)
      enddo
   enddo
   do iat = 1, n
      do jat = 1, n
         g(1,iat)=g(1,iat)+fgba(jat)*dcm5a(1,jat,iat)
         g(2,iat)=g(2,iat)+fgba(jat)*dcm5a(2,jat,iat)
         g(3,iat)=g(3,iat)+fgba(jat)*dcm5a(3,jat,iat)
         dcm5(1,jat)=dcm5(1,jat)+dcm5a(1,jat,iat)
         dcm5(2,jat)=dcm5(2,jat)+dcm5a(2,jat,iat)
         dcm5(3,jat)=dcm5(3,jat)+dcm5a(3,jat,iat)
      enddo
   enddo
   do iat = 1, n
      g(1,iat)=g(1,iat)-fgba(iat)*dcm5(1,iat)
      g(2,iat)=g(2,iat)-fgba(iat)*dcm5(2,iat)
      g(3,iat)=g(3,iat)-fgba(iat)*dcm5(3,iat)
   enddo
   if(lhb) then
      do iat = 1, n
         fgba(iat)=fhb(iat)*2.0_wp*q(iat)
      enddo
      do iat = 1, n
         do jat = 1, n
            g(1,iat)=g(1,iat)+fgba(jat)*dcm5a(1,jat,iat)
            g(2,iat)=g(2,iat)+fgba(jat)*dcm5a(2,jat,iat)
            g(3,iat)=g(3,iat)+fgba(jat)*dcm5a(3,jat,iat)
         enddo
      enddo
      do iat = 1, n
         g(1,iat)=g(1,iat)-fgba(iat)*dcm5(1,iat)
         g(2,iat)=g(2,iat)-fgba(iat)*dcm5(2,iat)
         g(3,iat)=g(3,iat)-fgba(iat)*dcm5(3,iat)
      enddo
   endif

end subroutine cm5_grad_gfn1

!! ========================================================================
!  repulsion gradient of GFN1
!! ========================================================================
subroutine rep_grad_gfn1(g,ep,n,at,xyz,sqrab,kexp,rexp)
   use aoparam, only : rep
   implicit none
   real(wp),intent(inout) :: g(3,n)
   real(wp),intent(out)   :: ep
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(in)    :: sqrab(n*(n+1)/2)
   real(wp),intent(in)    :: kexp
   real(wp),intent(in)    :: rexp

   integer,external :: lin
   integer  :: iat,jat,ati,atj
   real(wp) :: t16,t19,t20,t22,t26,t27,t28,t39
   real(wp) :: dum
   real(wp) :: alpha,repab
   real(wp) :: xa,ya,za,dx,dy,dz
   real(wp) :: r2,rab

   ep = 0.0_wp
   do iat=1,n-1
      xa=xyz(1,iat)
      ya=xyz(2,iat)
      za=xyz(3,iat)
      ati=at(iat)
      do jat=iat+1,n
         r2=sqrab(lin(jat,iat))
         if(r2.gt.5000.0d0) cycle
         dx=xa-xyz(1,jat)
         dy=ya-xyz(2,jat)
         dz=za-xyz(3,jat)
         rab=sqrt(r2)
         atj=at(jat)
         alpha=sqrt(rep(1,ati)*rep(1,atj))
         repab=rep(2,ati)*rep(2,atj)
         t16 = rab**kexp
         t19 = 1/r2
         t26 = dexp(-alpha*t16)
         t27 = rab**rexp
         t28 = 1/t27
         ep  = ep + repab * t26 * t28 !energy
         t20 = 1/r2/rab
         t22 = 2.D0*dx
         t39 = -0.5D0*repab*alpha*t16*kexp*t19*t22*t26*t28 &
         &     -0.5D0*repab*t26*t28*rexp*t19*t22
         g(1,iat)=g(1,iat)+t39
         g(1,jat)=g(1,jat)-t39
         t22 = 2.D0*dy
         t39 = -0.5D0*repab*alpha*t16*kexp*t19*t22*t26*t28 &
         &     -0.5D0*repab*t26*t28*rexp*t19*t22
         g(2,iat)=g(2,iat)+t39
         g(2,jat)=g(2,jat)-t39
         t22 = 2.D0*dz
         t39 = -0.5D0*repab*alpha*t16*kexp*t19*t22*t26*t28 &
         &     -0.5D0*repab*t26*t28*rexp*t19*t22
         g(3,iat)=g(3,iat)+t39
         g(3,jat)=g(3,jat)-t39
      enddo
   enddo

end subroutine rep_grad_gfn1

!! ========================================================================
!  repulsion gradient of GFN2
!! ========================================================================
subroutine rep_grad_gfn2(g,ep,n,at,xyz,sqrab,rexp)
   use aoparam, only : rep
   implicit none
   real(wp),intent(inout) :: g(3,n)
   real(wp),intent(out)   :: ep
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(in)    :: sqrab(n*(n+1)/2)
   real(wp),intent(in)    :: rexp

   integer,external :: lin
   integer  :: iat,jat,ati,atj
   real(wp) :: t16,t19,t20,t22,t26,t27,t28,t39
   real(wp) :: kexpe
   real(wp) :: alpha,repab
   real(wp) :: xa,ya,za,xb,yb,zb,dx,dy,dz
   real(wp) :: r2,rab

   ep = 0.0_wp
   do iat=1,n-1
      xa=xyz(1,iat)
      ya=xyz(2,iat)
      za=xyz(3,iat)
      ati=at(iat)
      do jat=iat+1,n
         r2=sqrab(lin(jat,iat))
         if(r2.gt.5000.0d0) cycle
         xb=xyz(1,jat)
         yb=xyz(2,jat)
         zb=xyz(3,jat)
         dx=xa-xyz(1,jat)
         dy=ya-xyz(2,jat)
         dz=za-xyz(3,jat)
         rab=sqrt(r2)
         atj=at(jat)
         alpha=sqrt(rep(1,ati)*rep(1,atj))
         repab=rep(2,ati)*rep(2,atj)
         if(ati.le.2.and.atj.le.2) then
            kexpe=1.0_wp
            t16 = rab
         else
            kexpe=1.5_wp
            t16 = rab**kexpe
         endif
         t19 = 1/r2
         t26 = dexp(-alpha*t16)
         t27 = rab**rexp
         t28 = 1/t27
         ep  = ep + repab * t26 * t28 !energy
         t20 = 1/r2/rab
         t22 = 2.D0*dx
         t39 = -0.5D0*repab*alpha*t16*kexpe*t19*t22*t26*t28 &
         &            -0.5D0*repab*t26*t28*rexp*t19*t22
         g(1,iat)=g(1,iat)+t39
         g(1,jat)=g(1,jat)-t39
         t22 = 2.D0*dy
         t39 = -0.5D0*repab*alpha*t16*kexpe*t19*t22*t26*t28 &
         &            -0.5D0*repab*t26*t28*rexp*t19*t22
         g(2,iat)=g(2,iat)+t39
         g(2,jat)=g(2,jat)-t39
         t22 = 2.D0*dz
         t39 = -0.5D0*repab*alpha*t16*kexpe*t19*t22*t26*t28 &
         &            -0.5D0*repab*t26*t28*rexp*t19*t22
         g(3,iat)=g(3,iat)+t39
         g(3,jat)=g(3,jat)-t39
      enddo
   enddo

end subroutine rep_grad_gfn2

!! ========================================================================
!  shellwise electrostatic gradient for GFN1
!! ========================================================================
subroutine shelles_grad_gfn1(g,n,at,nshell,xyz,sqrab,ash,lsh,alphaj,qsh)
   use aoparam, only : lpar,gam
   implicit none
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

   integer,external :: lin
   integer  :: is,js,iat,jat,ati,atj
   real(wp) :: xa,ya,za,dx,dy,dz
   real(wp) :: gi,gj,r2,rr,yy,ff

   do is=1,nshell
      iat=ash(is)
      xa=xyz(1,iat)
      ya=xyz(2,iat)
      za=xyz(3,iat)
      ati=at(iat)
      gi=gam(ati)*(1.0d0+lpar(lsh(is),ati))
      do js=1,nshell
         jat=ash(js)
         if(jat.le.iat) cycle
         dx=xa-xyz(1,jat)
         dy=ya-xyz(2,jat)
         dz=za-xyz(3,jat)
         atj=at(jat)
         r2=sqrab(lin(jat,iat))
         gj=gam(atj)*(1.0d0+lpar(lsh(js),atj))
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
subroutine shelles_grad_gfn2(g,n,at,nshell,xyz,sqrab,ash,lsh,qsh)
   use aoparam, only : lpar,gam
   implicit none
   real(wp),intent(inout) :: g(3,n)
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   integer, intent(in) :: nshell
   real(wp),intent(in) :: xyz(3,n)
   real(wp),intent(in) :: sqrab(n*(n+1)/2)
   integer, intent(in) :: ash(nshell)
   integer, intent(in) :: lsh(nshell)
   real(wp),intent(in) :: qsh(nshell)

   integer,external :: lin
   integer  :: is,js,iat,jat,ati,atj
   real(wp) :: xa,ya,za,dx,dy,dz
   real(wp) :: gi,gj,r2,rr,yy

   do is=1,nshell
      iat=ash(is)
      xa=xyz(1,iat)
      ya=xyz(2,iat)
      za=xyz(3,iat)
      ati=at(iat)
      gi=gam(ati)*(1.0d0+lpar(lsh(is),ati))
      do js=1,nshell
         jat=ash(js)
         if(jat.le.iat) cycle
         dx=xa-xyz(1,jat)
         dy=ya-xyz(2,jat)
         dz=za-xyz(3,jat)
         atj=at(jat)
         r2=sqrab(lin(jat,iat))
         gj=gam(atj)*(1.0d0+lpar(lsh(js),atj))
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
!  dH0/dxyz part from S(R) enhancement
!! ========================================================================
pure subroutine dhdr(n,at,xyz,iat,i,j,ati,atj,ishell,jshell,rab2,rf,dHdxyz)
   implicit none
   integer,intent(in)  :: n
   integer,intent(in)  :: at(n)
   real(wp), intent(in)  :: xyz(3,n)
   integer,intent(in)  :: iat
   integer,intent(in)  :: i
   integer,intent(in)  :: j
   integer,intent(in)  :: ishell
   integer,intent(in)  :: jshell
   integer,intent(in)  :: ati
   integer,intent(in)  :: atj
   real(wp), intent(in)  :: rab2
   real(wp), intent(out) :: rf
   real(wp), intent(out) :: dHdxyz(3)

!  real(wp)  :: rfactor,r1,r2
!  real(wp)  :: xyzi(3),xyzj(3)

!  dHdxyz=0

   call drfactor(ishell,jshell,iat,at(ati),at(atj),rab2, &
   &              xyz(1,ati),xyz(1,atj),rf,dHdxyz)
   if(atj.eq.iat) dHdxyz=-dHdxyz

end subroutine dhdr

!! ========================================================================
!  derivative of S(R) enhancement factor
!! ========================================================================
pure subroutine drfactor(ish,jsh,iat,ati,atj,rab2,xyz1,xyz2,rf,dxyz)
   use mctc_econv
   use aoparam, only : rad,polyr
   implicit none
   integer,intent(in)  :: ati,atj,ish,jsh,iat
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
   r=(rad(ati)+rad(atj))*aatoau

   rr=rab/r

   k1=polyr(ish,ati)*0.01
   k2=polyr(jsh,atj)*0.01

   t14 = rr**a
   t15 = k1*t14
   t17 = 1/rab2
   t22 = rr**a
   t23 = k2*t22
   rf=(1.0d0+t15)*(1.0d0+k2*t22)
   t20 = 2.D0*dx
   dxyz(1)=.5*t15*a*t17*t20*(1.+t23)+.5*(1.+t15)*k2*t22*a*t17*t20
   t20 = 2.D0*dy
   dxyz(2)=.5*t15*a*t17*t20*(1.+t23)+.5*(1.+t15)*k2*t22*a*t17*t20
   t20 = 2.D0*dz
   dxyz(3)=.5*t15*a*t17*t20*(1.+t23)+.5*(1.+t15)*k2*t22*a*t17*t20

end subroutine drfactor

end module grad_core
