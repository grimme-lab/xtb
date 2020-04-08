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
subroutine prep_grad_conv(ndim,H0,C,focc,emo,X)
   use xtb_mctc_convert, only : autoev,evtoau
   integer, intent(in)    :: ndim
   real(wp),intent(inout) :: H0(ndim*(ndim+1)/2)
   real(wp),intent(in)    :: C(ndim,ndim)
   real(wp),intent(in)    :: focc(ndim)
   real(wp),intent(in)    :: emo(ndim)
   real(wp),intent(out)   :: X(ndim,ndim)

   real(wp),allocatable :: temp(:)

   allocate( temp(ndim), source = 0.0_wp )
   H0 = H0*evtoau

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
   integer  :: iat,jat,il,jl,iZp,jZp
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
      il=mmm(lao2(i))
      jl=mmm(lao2(j))
      rab2 = sqrab(lin(jat,iat))
      hji = 2.0_wp*P(j,i)*S(j,i)
!     H0=S H (1+b x)
!     dH0/dx = H(1+bx)dS/dx + S H (b)
      call dshellPoly(hData%shellPoly(il,iZp),hData%shellPoly(jl,jZp),&
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
!  CN dependent part of the xTB Hamiltonian
!! ========================================================================
subroutine hcn_grad(hData,g,n,at,ndim,nmat2,matlist2,xyz, &
   &                P,S,dcndr,selfEnergy,dSEdcn, &
   &                aoat2,lao2,valao2,aoexp,ao2sh)
   use xtb_mctc_convert, only : autoev,evtoau
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   integer, intent(in)    :: ndim
   integer, intent(in)    :: nmat2
   integer, intent(in)    :: matlist2(2,nmat2)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(in)    :: P(ndim,ndim)
   real(wp),intent(in)    :: S(ndim,ndim)
   real(wp),intent(in)    :: dcndr(3,n,n)
   real(wp),intent(in)    :: selfEnergy(:)
   real(wp),intent(in)    :: dSEdcn(:)
   integer, intent(in)    :: aoat2(ndim)
   integer, intent(in)    :: lao2(ndim)
   integer, intent(in)    :: valao2(ndim)
   integer, intent(in)    :: ao2sh(ndim)
   real(wp),intent(in)    :: aoexp(ndim)
   real(wp),intent(inout) :: g(3,n)

   real(wp),allocatable :: hcn(:)

   integer  :: i,j,m
   integer  :: iat,jat,iZp,jZp
   integer  :: il,jl,ish,jsh
   real(wp) :: psij,dHdSE
   real(wp) :: gtmp(3)
   real(wp) :: km,fact
   real(wp) :: shPoly,dum1,dum2
   real(wp),parameter :: aot = -0.5_wp ! AO exponent dep. H0 scal

   allocate( hcn(n), source = 0.0_wp )

!  CN dependent part of H0
   hcn=0.0_wp
!$omp parallel default(none) &
!$omp shared(nmat2,matlist2,aoat2,lao2,valao2,P,S,n,at,xyz,hData,ndim) &
!$omp shared(aoexp,ao2sh,dSEdcn) &
!$omp private(i,m,j,jat,iat,ish,jsh,dum1,dum2,km,shPoly,il,jl,psij,dHdSe,iZp,jZp) &
!$omp reduction(+:hcn)
!$omp do
   do m = 1, nmat2
      i = matlist2(1,m)
      j = matlist2(2,m)
      iat = aoat2(i)
      jat = aoat2(j)
      ish = ao2sh(i)
      jsh = ao2sh(j)
      iZp = at(iat)
      jZp = at(jat)
      il = mmm(lao2(i))
      jl = mmm(lao2(j))
      psij = P(j,i)*S(j,i)
      shPoly = shellPoly(hData%shellPoly(il, iZp), hData%shellPoly(jl, jZp), &
         & hData%atomicRad(iZp), hData%atomicRad(jZp),xyz(:,iat),xyz(:,jat))
      call h0scal(hData,n,at,i,j,il,jl,iat,jat,valao2(i).ne.0,valao2(j).ne.0,  &
      &               km)
      km = km*(2*sqrt(aoexp(i)*aoexp(j))/(aoexp(i)+aoexp(j)))**hData%wExp
      dHdSE = psij*km*shPoly*evtoau
      hcn(jat)=hcn(jat) + dHdSE*dSEdcn(jsh) ! h independent part in H0
      hcn(iat)=hcn(iat) + dHdSE*dSEdcn(ish) ! h independent part in H0
   end do
!$omp end do
!$omp do
   do i=1,ndim
      iat = aoat2(i)
      ish = ao2sh(i)
      dHdSE = P(i,i)*evtoau ! diagonal contribution
      hcn(iat)=hcn(iat) + dHdSE*dSEdcn(ish)
   enddo
!$omp end do
!$omp end parallel


   ! CN-level shift gradient
   call contract(dcndr, hcn, g, beta=1.0_wp)

end subroutine hcn_grad

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
   call dsymv('u',n,1.0_wp,fgb,n,q,1,0.0_wp,fgba,1)
   call dgemv('n',3*n,n,1.0_wp,dcm5a,3*n,fgba,1,1.0_wp,g,1)
   if (lhb) then
      fgba = fhb * 2.0_wp * q
      call dgemv('n',3*n,n,1.0_wp,dcm5a,3*n,fgba,1,1.0_wp,g,1)
   end if

end subroutine cm5_grad_gfn1


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
