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

!! ========================================================================
!  This module implements various model Hessians, the actual definition
!  of the model Hessian is given in the description of each header, the
!  implementation following tries to do its very best to transform it back
!  from redundant internal coordinates (where the model is defined) to
!  Cartesian coordinates which are then used in any imaginable fashion
!  later on in the actual optimization.                        - SAW190131
!! ========================================================================
module xtb_modelhessian
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   use xtb_chargemodel
   implicit none

   public :: mh_lindh
   public :: mh_lindh_d2
   public :: mh_swart
   private

!  van-der-Waals radii used in the D2 model
   real(wp),parameter :: vander(86) = aatoau * (/ &
     0.91_wp,0.92_wp, & ! H, He
     0.75_wp,1.28_wp,1.35_wp,1.32_wp,1.27_wp,1.22_wp,1.17_wp,1.13_wp, & ! Li-Ne
     1.04_wp,1.24_wp,1.49_wp,1.56_wp,1.55_wp,1.53_wp,1.49_wp,1.45_wp, & ! Na-Ar
     1.35_wp,1.34_wp, & ! K, Ca
     1.42_wp,1.42_wp,1.42_wp,1.42_wp,1.42_wp, & ! Sc-Zn
     1.42_wp,1.42_wp,1.42_wp,1.42_wp,1.42_wp, &
     1.50_wp,1.57_wp,1.60_wp,1.61_wp,1.59_wp,1.57_wp, & ! Ga-Kr
     1.48_wp,1.46_wp, & ! Rb, Sr
     1.49_wp,1.49_wp,1.49_wp,1.49_wp,1.49_wp, & ! Y-Cd
     1.49_wp,1.49_wp,1.49_wp,1.49_wp,1.49_wp, &
     1.52_wp,1.64_wp,1.71_wp,1.72_wp,1.72_wp,1.71_wp, & ! In-Xe
     2.00_wp,2.00_wp, &
     2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp, & ! La-Yb
     2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp, &
     2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp, & ! Lu-Hg
     2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp, &
     2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp,2.00_wp /) ! Tl-Rn
!  C6 coefficients used in the D2 model
   real(wp),parameter :: c6(86) = (/&
      0.14_wp, 0.08_wp, & ! H,He
      1.61_wp, 1.61_wp, 3.13_wp, 1.75_wp, 1.23_wp, 0.70_wp, 0.75_wp, 0.63_wp, &
      5.71_wp, 5.71_wp,10.79_wp, 9.23_wp, 7.84_wp, 5.57_wp, 5.07_wp, 4.61_wp, &
     10.80_wp,10.80_wp, & ! K,Ca
     10.80_wp,10.80_wp,10.80_wp,10.80_wp,10.80_wp, & ! Sc-Zn
     10.80_wp,10.80_wp,10.80_wp,10.80_wp,10.80_wp, &
     16.99_wp,17.10_wp,16.37_wp,12.64_wp,12.47_wp,12.01_wp, & ! Ga-Kr
     24.67_wp,24.67_wp, & ! Rb,Sr
     24.67_wp,24.67_wp,24.67_wp,24.67_wp,24.67_wp, & ! Y-Cd
     24.67_wp,24.67_wp,24.67_wp,24.67_wp,24.67_wp, &
     37.32_wp,38.71_wp,38.44_wp,31.74_wp,31.50_wp,29.99_wp, & ! In-Xe
     50.00_wp,50.00_wp, & ! Cs,Ba
     50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp, & ! La-Yb
     50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp, &
     50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp, & ! Lu-Hg
     50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp, &
     50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp,50.00_wp /) ! Tl-Rn

   real(wp),parameter :: min_fk = 1.0e-3_wp

contains

!! ========================================================================
!  Swart's Model Hessian augmented with D2
!! ------------------------------------------------------------------------
!  Implemented after:
!  M. Swart, F. M. Bickelhaupt, Int. J. Quantum Chem., 2006, 106, 2536–2544.
!  DOI:10.1002/qua.21049
!
!  gij = exp[-(Rij/Cij-1)]
!  kij   = rkr·gij
!  kijk  = rkf·gij·gjk
!  kijkl = rkt·gij·gjk·gkl
!
!  The proposed force constants by Swart are:
!  rkr = 0.35, rkf = 0.15, rkt = 0.005
!
!  This Hessian is additionally augmented with D2, please note that D2
!  is not implemented in atomic units and requires some magical conversion
!  factor somewhere hidden in the implementation below.
!! ------------------------------------------------------------------------
subroutine mh_swart(xyz,n,hess,at,modh)
   use xtb_mctc_constants
   use xtb_mctc_convert
   use xtb_mctc_param, only: rad => covalent_radius_2009

   use xtb_type_setvar
   use xtb_type_param

   implicit none

   integer, intent(in)  :: n
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(out) :: hess((3*n)*(3*n+1)/2)
   integer, intent(in)  :: at(n)
   type(modhess_setvar),intent(in) :: modh

   integer  :: n3
   real(wp),parameter :: rzero = 1.0e-10_wp
   logical, allocatable :: lcutoff(:,:)
   type(chrg_parameter) :: chrgeq
   real(wp) :: kd

   allocate( lcutoff(n,n), source=.false.)

   n3=3*n
   hess = 0.0d0

!  the dispersion force constant is used relative to the stretch force constant
   kd = modh%kd/modh%kr

   call mh_swart_stretch(n,at,xyz,hess,modh%kr,kd,modh%s6,rad,rad,lcutoff,modh%rcut)
   if (modh%kf.ne.0.0_wp) &
   call mh_swart_bend   (n,at,xyz,hess,modh%kf,kd,        rad,rad,lcutoff)
   if (modh%kt.ne.0.0_wp) &
   call mh_swart_torsion(n,at,xyz,hess,modh%kt,kd,        rad,rad,lcutoff)
   if (modh%ko.ne.0.0_wp) &
   call mh_swart_outofp (n,at,xyz,hess,modh%ko,kd,        rad,rad,lcutoff)
   if (modh%kq.ne.0.0_wp) then
      call new_charge_model_2019(chrgeq,n,at)
      call mh_eeq(n,at,xyz,0.0_wp,chrgeq,modh%kq,hess)
   endif

end subroutine mh_swart

pure subroutine mh_swart_stretch(n,at,xyz,hess,kr,kd,s6,rcov,rvdw,lcutoff,rcut)
   use xtb_mctc_constants
   use xtb_mctc_convert

   implicit none

   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: kr
   real(wp),intent(in)    :: kd
   real(wp),intent(in)    :: s6
   real(wp),intent(in)    :: rcov(:)
   real(wp),intent(in)    :: rvdw(:)
   logical, intent(out)   :: lcutoff(n,n)
   real(wp),intent(in)    :: rcut

   integer  :: i,j
   real(wp) :: xij,yij,zij,rij2,r0,d0
   real(wp) :: gmm
   real(wp) :: c6i,c6j,c6ij,rv
   real(wp) :: hxx,hxy,hxz,hyy,hyz,hzz
   real(wp) :: vdw(3,3)

!! ------------------------------------------------------------------------
!  Hessian for stretch
!! ------------------------------------------------------------------------
   stretch_iAt: do i = 1, n

      stretch_jAt: do j = 1, i-1

         ! save for later
         lcutoff(i,j) = rcutoff(xyz,i,j,rcut)
         lcutoff(j,i) = lcutoff(i,j)

         xij=xyz(1,i)-xyz(1,j)
         yij=xyz(2,i)-xyz(2,j)
         zij=xyz(3,i)-xyz(3,j)
         rij2 = xij**2 + yij**2 + zij**2
         r0 = rcov(at(i))+rcov(at(j))
         d0 = rvdw(at(i))+rvdw(at(j))

         !cccccc vdwx ccccccccccccccccccccccccccccccccc
         c6i=c6(at(i))
         c6j=c6(at(j))
         c6ij=sqrt(c6i*c6j)
         rv=(vander(at(i))+vander(at(j)))*aatoau

         call getvdwxx(xij, yij, zij, c6ij, s6, rv, vdw(1,1))
         call getvdwxy(xij, yij, zij, c6ij, s6, rv, vdw(1,2))
         call getvdwxy(xij, zij, yij, c6ij, s6, rv, vdw(1,3))
         call getvdwxx(yij, xij, zij, c6ij, s6, rv, vdw(2,2))
         call getvdwxy(yij, zij, xij, c6ij, s6, rv, vdw(2,3))
         call getvdwxx(zij, xij, yij, c6ij, s6, rv, vdw(3,3))
         !cccccc ende vdwx ccccccccccccccccccccccccccccccc

         gmm = kr*fk_swart(1.0_wp,  r0,rij2) &
            + kr*kd * fk_vdw(5.0_wp,d0,rij2)

         !gmm = max(gmm,min_fk)

         hxx=gmm*xij*xij/rij2-vdw(1,1)
         hxy=gmm*xij*yij/rij2-vdw(1,2)
         hxz=gmm*xij*zij/rij2-vdw(1,3)
         hyy=gmm*yij*yij/rij2-vdw(2,2)
         hyz=gmm*yij*zij/rij2-vdw(2,3)
         hzz=gmm*zij*zij/rij2-vdw(3,3)

         ! save diagonal elements for atom i
         hess(ind(1,i,1,i))=hess(ind(1,i,1,i))+hxx
         hess(ind(2,i,1,i))=hess(ind(2,i,1,i))+hxy
         hess(ind(2,i,2,i))=hess(ind(2,i,2,i))+hyy
         hess(ind(3,i,1,i))=hess(ind(3,i,1,i))+hxz
         hess(ind(3,i,2,i))=hess(ind(3,i,2,i))+hyz
         hess(ind(3,i,3,i))=hess(ind(3,i,3,i))+hzz
         ! save elements between atom i and atom j
         hess(ind(1,i,1,j))=hess(ind(1,i,1,j))-hxx
         hess(ind(1,i,2,j))=hess(ind(1,i,2,j))-hxy
         hess(ind(1,i,3,j))=hess(ind(1,i,3,j))-hxz
         hess(ind(2,i,1,j))=hess(ind(2,i,1,j))-hxy
         hess(ind(2,i,2,j))=hess(ind(2,i,2,j))-hyy
         hess(ind(2,i,3,j))=hess(ind(2,i,3,j))-hyz
         hess(ind(3,i,1,j))=hess(ind(3,i,1,j))-hxz
         hess(ind(3,i,2,j))=hess(ind(3,i,2,j))-hyz
         hess(ind(3,i,3,j))=hess(ind(3,i,3,j))-hzz
         ! save diagonal elements for atom j
         hess(ind(1,j,1,j))=hess(ind(1,j,1,j))+hxx
         hess(ind(2,j,1,j))=hess(ind(2,j,1,j))+hxy
         hess(ind(2,j,2,j))=hess(ind(2,j,2,j))+hyy
         hess(ind(3,j,1,j))=hess(ind(3,j,1,j))+hxz
         hess(ind(3,j,2,j))=hess(ind(3,j,2,j))+hyz
         hess(ind(3,j,3,j))=hess(ind(3,j,3,j))+hzz

      end do stretch_jAt
   end do stretch_iAt

end subroutine mh_swart_stretch

pure subroutine mh_swart_bend(n,at,xyz,hess,kf,kd,rcov,rvdw,lcutoff)
   use xtb_mctc_constants

   implicit none

   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: kf
   real(wp),intent(in)    :: kd
   real(wp),intent(in)    :: rcov(:)
   real(wp),intent(in)    :: rvdw(:)
   logical, intent(in)    :: lcutoff(n,n)

   integer  :: i,j,m,ic,jc,ii
   real(wp),parameter :: rzero = 1.0e-10_wp
   real(wp) :: xij,yij,zij,rij2,rrij,r1
   real(wp) :: xmi,ymi,zmi,rmi2,rmi,r0mi,d0mj,gmi
   real(wp) :: xmj,ymj,zmj,rmj2,rmj,r0mj,d0mi,gmj
   real(wp) :: test,gij,rl2,rl,rmidotrmj
   real(wp) :: sinphi,cosphi,costhetax,costhetay,costhetaz
   real(wp) :: alpha
   real(wp) :: si(3),sj(3),sm(3),x(2),y(2),z(2)

!! ------------------------------------------------------------------------
!  Hessian for bending
!! ------------------------------------------------------------------------
   bend_mAt: do m = 1, n
      bend_iAt: do i = 1, n
         if (i.eq.m) cycle bend_iAt
         if(lcutoff(i,m)) cycle bend_iAt

         xmi=(xyz(1,i)-xyz(1,m))
         ymi=(xyz(2,i)-xyz(2,m))
         zmi=(xyz(3,i)-xyz(3,m))
         rmi2 = xmi**2 + ymi**2 + zmi**2
         rmi=sqrt(rmi2)
         r0mi=rcov(at(m))+rcov(at(i))
         d0mi=rvdw(at(m))+rvdw(at(i))

         bend_jAt: do j = 1, i-1
            if (j.eq.m) cycle bend_jAt
            if(lcutoff(j,i)) cycle bend_jAt
            if(lcutoff(j,m)) cycle bend_jAt

            xmj=(xyz(1,j)-xyz(1,m))
            ymj=(xyz(2,j)-xyz(2,m))
            zmj=(xyz(3,j)-xyz(3,m))
            rmj2 = xmj**2 + ymj**2 + zmj**2
            rmj=sqrt(rmj2)
            r0mj=rcov(at(m))+rcov(at(j))
            d0mj=rvdw(at(m))+rvdw(at(j))

            ! test if zero angle
            test=xmi*xmj+ymi*ymj+zmi*zmj
            test=test/(rmi*rmj)
            if (abs(test-1.0_wp).lt.1.0e-12_wp) cycle bend_jAt

            xij=(xyz(1,j)-xyz(1,i))
            yij=(xyz(2,j)-xyz(2,i))
            zij=(xyz(3,j)-xyz(3,i))
            rij2 = xij**2 + yij**2 + zij**2
            rrij=sqrt(rij2)

            gmi = fk_swart(1.0_wp,r0mi,rmi2) &
                + 0.5_wp*kd * fk_vdw(5.0_wp,d0mi,rmi2)
            gmj = fk_swart(1.0_wp,r0mj,rmj2) &
                + 0.5_wp*kd * fk_vdw(5.0_wp,d0mj,rmj2)

            gij = kf*gmi*gmj

            rl2=(ymi*zmj-zmi*ymj)**2+(zmi*xmj-xmi*zmj)**2+(xmi*ymj-ymi*xmj)**2

            if(rl2.lt.1.e-14_wp) then
               rl=0.0_wp
            else
               rl=sqrt(rl2)
            end if

            !gij = max(gij,min_fk)

            if ((rmj.gt.rzero).and.(rmi.gt.rzero).and.(rrij.gt.rzero)) then
               sinphi=rl/(rmj*rmi)
               rmidotrmj=xmi*xmj+ymi*ymj+zmi*zmj
               cosphi=rmidotrmj/(rmj*rmi)
               ! none linear case
               if (sinphi.gt.rzero) then
                  si(1)=(xmi/rmi*cosphi-xmj/rmj)/(rmi*sinphi)
                  si(2)=(ymi/rmi*cosphi-ymj/rmj)/(rmi*sinphi)
                  si(3)=(zmi/rmi*cosphi-zmj/rmj)/(rmi*sinphi)
                  sj(1)=(cosphi*xmj/rmj-xmi/rmi)/(rmj*sinphi)
                  sj(2)=(cosphi*ymj/rmj-ymi/rmi)/(rmj*sinphi)
                  sj(3)=(cosphi*zmj/rmj-zmi/rmi)/(rmj*sinphi)
                  sm(1)=-si(1)-sj(1)
                  sm(2)=-si(2)-sj(2)
                  sm(3)=-si(3)-sj(3)
                  do ic=1,3
                     do jc=1,3
                        if (m.gt.i) then
                           hess(ind(ic,m,jc,i))=hess(ind(ic,m,jc,i)) &
                              +gij*sm(ic)*si(jc)
                        else
                           hess(ind(ic,i,jc,m))=hess(ind(ic,i,jc,m)) &
                              +gij*si(ic)*sm(jc)
                        end if
                        if (m.gt.j) then
                           hess(ind(ic,m,jc,j))=hess(ind(ic,m,jc,j)) &
                              +gij*sm(ic)*sj(jc)
                        else
                           hess(ind(ic,j,jc,m))=hess(ind(ic,j,jc,m)) &
                              +gij*sj(ic)*sm(jc)
                        end if
                        if (i.gt.j) then
                           hess(ind(ic,i,jc,j))=hess(ind(ic,i,jc,j)) &
                              +gij*si(ic)*sj(jc)
                        else
                           hess(ind(ic,j,jc,i))=hess(ind(ic,j,jc,i)) &
                              +gij*sj(ic)*si(jc)
                        end if
                     end do
                  end do
                  do ic=1,3
                     do jc=1,ic
                        hess(ind(ic,i,jc,i))=hess(ind(ic,i,jc,i))+gij*si(ic)*si(jc)
                        hess(ind(ic,m,jc,m))=hess(ind(ic,m,jc,m))+gij*sm(ic)*sm(jc)
                        hess(ind(ic,j,jc,j))=hess(ind(ic,j,jc,j))+gij*sj(ic)*sj(jc)
                     end do
                  end do
               else
                  ! linear case
                  if ((abs(ymi).gt.rzero).or.(abs(xmi).gt.rzero)) then
                     x(1)=-ymi
                     y(1)=xmi
                     z(1)=0.0_wp
                     x(2)=-xmi*zmi
                     y(2)=-ymi*zmi
                     z(2)=xmi*xmi+ymi*ymi
                  else
                     x(1)=1.0_wp
                     y(1)=0.0_wp
                     z(1)=0.0_wp
                     x(2)=0.0_wp
                     y(2)=1.0_wp
                     z(2)=0.0_wp
                  end if
                  do ii=1,2
                     r1=sqrt(x(ii)**2+y(ii)**2+z(ii)**2)
                     costhetax=x(ii)/r1
                     costhetay=y(ii)/r1
                     costhetaz=z(ii)/r1
                     si(1)=-costhetax/rmi
                     si(2)=-costhetay/rmi
                     si(3)=-costhetaz/rmi
                     sj(1)=-costhetax/rmj
                     sj(2)=-costhetay/rmj
                     sj(3)=-costhetaz/rmj
                     sm(1)=-(si(1)+sj(1))
                     sm(2)=-(si(2)+sj(2))
                     sm(3)=-(si(3)+sj(3))
                     !
                     do ic=1,3
                        do jc=1,3
                           if (m.gt.i) then
                              hess(ind(ic,m,jc,i))=hess(ind(ic,m,jc,i)) &
                                 +gij*sm(ic)*si(jc)
                           else
                              hess(ind(ic,i,jc,m))=hess(ind(ic,i,jc,m)) &
                                 +gij*si(ic)*sm(jc)
                           end if
                           if (m.gt.j) then
                              hess(ind(ic,m,jc,j))=hess(ind(ic,m,jc,j)) &
                                 +gij*sm(ic)*sj(jc)
                           else
                              hess(ind(ic,j,jc,m))=hess(ind(ic,j,jc,m)) &
                                 +gij*sj(ic)*sm(jc)
                           end if
                           if (i.gt.j) then
                              hess(ind(ic,i,jc,j))=hess(ind(ic,i,jc,j)) &
                                 +gij*si(ic)*sj(jc)
                           else
                              hess(ind(ic,j,jc,i))=hess(ind(ic,j,jc,i)) &
                                 +gij*sj(ic)*si(jc)
                           end if
                        end do
                     end do
                     do ic=1,3
                        do jc=1,ic
                           hess(ind(ic,i,jc,i))=hess(ind(ic,i,jc,i)) &
                              +gij*si(ic)*si(jc)
                           hess(ind(ic,m,jc,m))=hess(ind(ic,m,jc,m)) &
                              +gij*sm(ic)*sm(jc)
                           hess(ind(ic,j,jc,j))=hess(ind(ic,j,jc,j)) &
                              +gij*sj(ic)*sj(jc)
                        end do
                     end do
                  end do

               end if
            end if

         end do bend_jAt
      end do bend_iAt
   end do bend_mAt

end subroutine mh_swart_bend

pure subroutine mh_swart_torsion(n,at,xyz,hess,kt,kd,rcov,rvdw,lcutoff)
   use xtb_mctc_constants

   implicit none

   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: kt
   real(wp),intent(in)    :: kd
   real(wp),intent(in)    :: rcov(:)
   real(wp),intent(in)    :: rvdw(:)
   logical, intent(in)    :: lcutoff(n,n)

   integer  :: i,j,k,l,ic,jc,ij,kl
!  allow only angles in the range of 35-145
   real(wp),parameter :: a35 = (35.0d0/180.d0)* pi
   real(wp),parameter :: cosfi_max=cos(a35)
   real(wp) :: txyz(3,4),c(3,4)
   real(wp) :: rij(3),rij0,aij,rij2,d0ij,gij
   real(wp) :: rjk(3),rjk0,ajk,rjk2,d0jk,gjk
   real(wp) :: rkl(3),rkl0,akl,rkl2,d0kl,gkl
   real(wp) :: cosfi2,cosfi3,cosfi4
   real(wp) :: beta,tij,tau
   real(wp) :: si(3),sj(3),sk(3),sl(3)

!! ------------------------------------------------------------------------
!  Hessian for torsion
!! ------------------------------------------------------------------------
   torsion_jAt: do j = 1,n
      txyz(:,2)=xyz(:,j)
      torsion_kAt: do k = 1, n
         if (k.eq.j) cycle torsion_kAt
         if(lcutoff(k,j)) cycle torsion_kAt
         txyz(:,3) = xyz(:,k)
         torsion_iAt: do i = 1, n
            ij=n*(j-1)+i
            if (i.eq.j) cycle torsion_iAt
            if (i.eq.k) cycle torsion_iAt
            if(lcutoff(i,k)) cycle torsion_iAt
            if(lcutoff(i,j)) cycle torsion_iAt

            txyz(:,1)=xyz(:,i)
            torsion_lAt: do l = 1, n
               if (ij.le.kl) cycle torsion_lAt
               if (l.eq.i)   cycle torsion_lAt
               if (l.eq.j)   cycle torsion_lAt
               if (l.eq.k)   cycle torsion_lAt
!
               if(lcutoff(l,i)) cycle torsion_lAt
               if(lcutoff(l,k)) cycle torsion_lAt
               if(lcutoff(l,j)) cycle torsion_lAt

               txyz(:,4)=xyz(:,l)

               rij=xyz(:,i)-xyz(:,j)
               d0ij=rvdw(at(i))+rvdw(at(j))
               rij0=rcov(at(i))+rcov(at(j))

               rjk=xyz(:,j)-xyz(:,k)
               d0jk=rvdw(at(j))+rvdw(at(k))
               rjk0=rcov(at(j))+rcov(at(k))

               rkl=xyz(:,k)-xyz(:,l)
               d0kl=rvdw(at(k))+rvdw(at(l))
               rkl0=rcov(at(k))+rcov(at(l))

               rij2=sum(rij**2)
               rjk2=sum(rjk**2)
               rkl2=sum(rjk**2)

               cosfi2=dot_product(rij,rjk)/sqrt(rij2*rjk2)
               if (abs(cosfi2).gt.cosfi_max) cycle
               cosfi3=dot_product(rkl,rjk)/sqrt(rkl2*rjk2)
               if (abs(cosfi3).gt.cosfi_max) cycle

               gij = fk_swart(1.0_wp,rij0,rij2) &
                  + 0.5_wp*kd * fk_vdw(5.0_wp,d0ij,rij2)
               gjk = fk_swart(1.0_wp,rjk0,rjk2) &
                  + 0.5_wp*kd * fk_vdw(5.0_wp,d0jk,rjk2)
               gkl = fk_swart(1.0_wp,rkl0,rkl2) &
                  + 0.5_wp*kd * fk_vdw(5.0_wp,d0kl,rkl2)

               tij = kt * gij*gjk*gkl

               !tij = max(tij,10*min_fk)

               call trsn2(txyz,tau,c)
               si = c(:,1)
               sj = c(:,2)
               sk = c(:,3)
               sl = c(:,4)

               ! off diagonal block
               do ic=1,3
                  do jc=1,3
                     hess(ind(ic,i,jc,j))=hess(ind(ic,i,jc,j))+tij*si(ic) * sj(jc)
                     hess(ind(ic,i,jc,k))=hess(ind(ic,i,jc,k))+tij*si(ic) * sk(jc)
                     hess(ind(ic,i,jc,l))=hess(ind(ic,i,jc,l))+tij*si(ic) * sl(jc)
                     hess(ind(ic,j,jc,k))=hess(ind(ic,j,jc,k))+tij*sj(ic) * sk(jc)
                     hess(ind(ic,j,jc,l))=hess(ind(ic,j,jc,l))+tij*sj(ic) * sl(jc)
                     hess(ind(ic,k,jc,l))=hess(ind(ic,k,jc,l))+tij*sk(ic) * sl(jc)
                  end do
               end do

               ! diagonal block
               do ic=1,3
                  do jc=1,ic
                     hess(ind(ic,i,jc,i))=hess(ind(ic,i,jc,i))+tij*si(ic) * si(jc)
                     hess(ind(ic,j,jc,j))=hess(ind(ic,j,jc,j))+tij*sj(ic) * sj(jc)
                     hess(ind(ic,k,jc,k))=hess(ind(ic,k,jc,k))+tij*sk(ic) * sk(jc)
                     hess(ind(ic,l,jc,l))=hess(ind(ic,l,jc,l))+tij*sl(ic) * sl(jc)
                  end do
               end do

            end do torsion_lAt
         end do torsion_iAt
      end do torsion_kAt
   end do torsion_jAt

end subroutine mh_swart_torsion

pure subroutine mh_swart_outofp(n,at,xyz,hess,ko,kd,rcov,rvdw,lcutoff)
   use xtb_mctc_constants

   implicit none

   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: ko
   real(wp),intent(in)    :: kd
   real(wp),intent(in)    :: rcov(:)
   real(wp),intent(in)    :: rvdw(:)
   logical, intent(in)    :: lcutoff(n,n)

   integer  :: i,ir,j,jr,k,kr,l,lr,ic,jc
   real(wp) :: txyz(3,4),c(3,4)
   real(wp) :: rij(3),rij0,d0ij,rij2,gij
   real(wp) :: rik(3),rik0,d0ik,rik2,gik
   real(wp) :: ril(3),ril0,d0il,ril2,gil
   real(wp) :: cosfi2,cosfi3,cosfi4
   real(wp) :: beta,tij,tau
   real(wp) :: si(3),sj(3),sk(3),sl(3)

!! ------------------------------------------------------------------------
!  Hessian for out-of-plane
!! ------------------------------------------------------------------------
   outofplane_iAt: do i = 1, n
      txyz(:,4) = xyz(:,i)
      outofplane_jAt: do j = 1, n
         if (j.eq.i) cycle outofplane_jAt
         if(lcutoff(j,i)) cycle outofplane_jAt
         txyz(:,1) = xyz(:,j)
         outofplane_kAt: do k = 1, n
            if (i.eq.k) cycle outofplane_kAt
            if (j.eq.k) cycle outofplane_kat
            if(lcutoff(k,i)) cycle outofplane_kAt
            if(lcutoff(k,j)) cycle outofplane_kAt
            txyz(:,2) = xyz(:,k)
            outofplane_lAt: do l = 1, n
               txyz(:,3) = xyz(:,l)
               if (l.eq.i)   cycle outofplane_lAt
               if (l.eq.j)   cycle outofplane_lAt
               if (l.eq.k)   cycle outofplane_lAt
               if(lcutoff(l,i)) cycle outofplane_lAt
               if(lcutoff(l,k)) cycle outofplane_lAt
               if(lcutoff(l,j)) cycle outofplane_lAt

               rij=xyz(:,i)-xyz(:,j)
               rij0=rcov(at(i))+rcov(at(j))
               d0ij=rvdw(at(i))+rvdw(at(j))

               rik=xyz(:,i)-xyz(:,k)
               rik0=rcov(at(i))+rcov(at(k))
               d0ik=rvdw(at(i))+rvdw(at(k))

               ril=xyz(:,i)-xyz(:,l)
               ril0=rcov(at(i))+rcov(at(l))
               d0il=rvdw(at(i))+rvdw(at(l))

               rij2=sum(rij**2)
               rik2=sum(rik**2)
               ril2=sum(ril**2)

               cosfi2=dot_product(rij,rik)/sqrt(rij2*rik2)
               if (abs(abs(cosfi2)-1.0_wp).lt.1.0e-1_wp) cycle
               cosfi3=dot_product(rij,ril)/sqrt(rij2*ril2)
               if (abs(abs(cosfi3)-1.0_wp).lt.1.0e-1_wp) cycle
               cosfi4=dot_product(rik,ril)/sqrt(rik2*ril2)
               if (abs(abs(cosfi4)-1.0_wp).lt.1.0e-1_wp) cycle

               gij = fk_swart(1.0_wp,rij0,rij2) &
                  + 0.5_wp*kd * fk_vdw(5.0_wp,d0ij,rij2)
               gik = fk_swart(1.0_wp,rik0,rik2) &
                  + 0.5_wp*kd * fk_vdw(5.0_wp,d0ik,rik2)
               gil = fk_swart(1.0_wp,ril0,ril2) &
                  + 0.5_wp*kd * fk_vdw(5.0_wp,d0il,ril2)

               tij = ko * gij*gik*gil

               !tij = max(tij,10*min_fk)

               call outofp2(xyz,tau,c)
               If (abs(tau).gt.45.0d0*(pi/180.d0)) cycle

               si = c(:,4)
               sj = c(:,1)
               sk = c(:,2)
               sl = c(:,3)

               ! off diagonal block
               do ic=1,3
                  do jc=1,3
                     hess(ind(ic,i,jc,j))=hess(ind(ic,i,jc,j))+tij*si(ic) * sj(jc)
                     hess(ind(ic,i,jc,k))=hess(ind(ic,i,jc,k))+tij*si(ic) * sk(jc)
                     hess(ind(ic,i,jc,l))=hess(ind(ic,i,jc,l))+tij*si(ic) * sl(jc)
                     hess(ind(ic,j,jc,k))=hess(ind(ic,j,jc,k))+tij*sj(ic) * sk(jc)
                     hess(ind(ic,j,jc,l))=hess(ind(ic,j,jc,l))+tij*sj(ic) * sl(jc)
                     hess(ind(ic,k,jc,l))=hess(ind(ic,k,jc,l))+tij*sk(ic) * sl(jc)
                  end do
               end do

               ! diagonal block
               do ic=1,3
                  do jc=1,ic
                     hess(ind(ic,i,jc,i))=hess(ind(ic,i,jc,i))+tij*si(ic) * si(jc)
                     hess(ind(ic,j,jc,j))=hess(ind(ic,j,jc,j))+tij*sj(ic) * sj(jc)
                     hess(ind(ic,k,jc,k))=hess(ind(ic,k,jc,k))+tij*sk(ic) * sk(jc)
                     hess(ind(ic,l,jc,l))=hess(ind(ic,l,jc,l))+tij*sl(ic) * sl(jc)
                  end do
               end do

            enddo outofplane_lAt
         enddo outofplane_kAt
      enddo outofplane_jAt
   enddo outofplane_iAt

end subroutine mh_swart_outofp

!! ========================================================================
!  Lindh's Model Hessian augmented with D2
!! ------------------------------------------------------------------------
!  Implemented after:
!  Lindh, R., Bernhardsson, A., Karlström, G., & Malmqvist, P.-Å. (1995).
!  On the use of a Hessian model function in molecular geometry optimizations.
!  Chem. Phys. Lett., 241(4), 423–428. doi:10.1016/0009-2614(95)00646-l
!
!  gij = exp[αij(R²ref - R²ij)]
!  kij   = rkr·gij
!  kijk  = rkf·gij·gjk
!  kijkl = rkt·gij·gjk·gkl
!
!  Originally Lindh proposed (we tweaked those a little bit):
!  rkr = 0.45, rkf = 0.15, rkt = 0.005
!
!  the reference distances are divided by rows in the PSE:
!  rAv:        1        2        3       aAv:        1        2        3
!    1    1.3500   2.1000   2.5300         1    1.0000   0.3949   0.3949
!    2    2.1000   2.8700   3.4000         2    0.3949   0.2800   0.2800
!    3    2.5300   3.4000   3.4000         3    0.3949   0.2800   0.2800
!
!  This Hessian is additionally augmented with D2, please note that D2
!  is not implemented in atomic units and requires some magical conversion
!  factor somewhere hidden in the implementation below.
!! ------------------------------------------------------------------------
subroutine mh_lindh_d2(xyz,n,hess,at,modh)
   use xtb_mctc_constants
   use xtb_mctc_convert

   use xtb_type_setvar
   use xtb_type_param

   implicit none

   integer, intent(in)  :: n
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(out) :: hess((3*n)*(3*n+1)/2)
   integer, intent(in)  :: at(n)
   type(modhess_setvar),intent(in) :: modh

   real(wp),parameter :: rAv(3,3) = reshape( &
     (/1.3500_wp,2.1000_wp,2.5300_wp, &
       2.1000_wp,2.8700_wp,3.4000_wp, &
       2.5300_wp,3.4000_wp,3.4000_wp/), shape(rAv) )
   real(wp),parameter :: aAv(3,3) = reshape ( &
     (/1.0000_wp,0.3949_wp,0.3949_wp, &
       0.3949_wp,0.2800_wp,0.2800_wp, &
       0.3949_wp,0.2800_wp,0.2800_wp/), shape(aAv) )
   real(wp),parameter :: dAv(3,3) = reshape ( &
     (/0.0000_wp,0.0000_wp,0.0000_wp, &
       0.0000_wp,0.0000_wp,0.0000_wp, &
       0.0000_wp,0.0000_wp,0.0000_wp/), shape(aAv) )

   integer  :: n3
   real(wp) :: kd
   logical, allocatable :: lcutoff(:,:)
   type(chrg_parameter) :: chrgeq

   allocate( lcutoff(n,n), source=.false.)

   n3=3*n
   hess = 0.0d0

!  the dispersion force constant is used relative to the stretch force constant
   kd = modh%kd/modh%kr

   call mh_lindh_stretch(n,at,xyz,hess,modh%kr,kd,modh%s6,aav,rav,dav,lcutoff,modh%rcut)
   if (modh%kf.ne.0.0_wp) &
   call mh_lindh_bend   (n,at,xyz,hess,modh%kf,kd,        aav,rav,dav,lcutoff)
   if (modh%kt.ne.0.0_wp) &
   call mh_lindh_torsion(n,at,xyz,hess,modh%kt,kd,        aav,rav,dav,lcutoff)
   if (modh%ko.ne.0.0_wp) &
   call mh_lindh_outofp (n,at,xyz,hess,modh%ko,kd,        aav,rav,dav,lcutoff)
   if (modh%kq.ne.0.0_wp) then
      call new_charge_model_2019(chrgeq,n,at)
      call mh_eeq(n,at,xyz,0.0_wp,chrgeq,modh%kq,hess)
   endif

end subroutine mh_lindh_d2

!! ========================================================================
!  Lindh's Model Hessian updated around 2007
!! ------------------------------------------------------------------------
!  R. Lindh, personal communication.
!
!  gij = exp[αij(R²ref - R²ij)]
!  dij = exp[-4·(Rvdw - Rij)²]
!  kij   = rkr·gij + rkd·dij
!  kijk  = rkf·(gij+½·rkd/rkr·dij)·(gjk+½·rkd/rkr·djk)
!  kijkl = rkt·(gij+½·rkd/rkr·dij)·(gjk+½·rkd/rkr·djk)·(gkl+½·rkd/rkr·dkl)
!
!  parameters tweaked by R. Lindh in 2007:
!  rkr = 0.45, rkf = 0.10, rkt = 0.0025, rko = 0.16, rkd = 0.05
!
!  the reference distances are divided by rows in the PSE:
!  rAv:        1        2        3       aAv:        1        2        3
!    1    1.3500   2.1000   2.5300         1    1.0000   0.3949   0.3949
!    2    2.1000   2.8700   3.8000         2    0.3949   0.2800   0.1200
!    3    2.5300   3.8000   4.5000         3    0.3949   0.1200   0.0600
!
!  dAv:        1        2        3
!    1    0.0000   3.6000   3.6000
!    2    3.6000   5.3000   5.3000
!    3    3.6000   5.3000   5.3000
!
!! ------------------------------------------------------------------------
subroutine mh_lindh(xyz,n,hess,at,modh)
   use xtb_mctc_constants
   use xtb_mctc_convert

   use xtb_type_setvar
   use xtb_type_param

   implicit none

   integer, intent(in)  :: n
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(out) :: hess((3*n)*(3*n+1)/2)
   integer, intent(in)  :: at(n)
   type(modhess_setvar),intent(in) :: modh

   real(wp),parameter :: rAv(3,3) = reshape( &
     (/1.3500_wp,2.1000_wp,2.5300_wp, &
       2.1000_wp,2.8700_wp,3.8000_wp, &
       2.5300_wp,3.8000_wp,4.5000_wp/), shape(rAv) )
   real(wp),parameter :: aAv(3,3) = reshape ( &
     (/1.0000_wp,0.3949_wp,0.3949_wp, &
       0.3949_wp,0.2800_wp,0.1200_wp, &
       0.3949_wp,0.1200_wp,0.0600_wp/), shape(aAv) )
   real(wp),parameter :: dAv(3,3) = reshape ( &
     (/0.0000_wp,3.6000_wp,3.6000_wp, &
       3.6000_wp,5.3000_wp,5.3000_wp, &
       3.6000_wp,5.3000_wp,5.3000_wp/), shape(aAv) )

   integer  :: n3
   real(wp) :: kd
   logical, allocatable :: lcutoff(:,:)
   type(chrg_parameter) :: chrgeq

   allocate( lcutoff(n,n), source=.false.)

   n3=3*n
   hess = 0.0d0

!  the dispersion force constant is used relative to the stretch force constant
   kd = modh%kd/modh%kr

   call mh_lindh_stretch(n,at,xyz,hess,modh%kr,kd,modh%s6,aav,rav,dav,lcutoff,modh%rcut)
   if (modh%kf.ne.0.0_wp) &
   call mh_lindh_bend   (n,at,xyz,hess,modh%kf,kd,        aav,rav,dav,lcutoff)
   if (modh%kt.ne.0.0_wp) &
   call mh_lindh_torsion(n,at,xyz,hess,modh%kt,kd,        aav,rav,dav,lcutoff)
   if (modh%ko.ne.0.0_wp) &
   call mh_lindh_outofp (n,at,xyz,hess,modh%ko,0.0_wp,    aav,rav,dav,lcutoff)
   if (modh%kq.ne.0.0_wp) then
      call new_charge_model_2019(chrgeq,n,at)
      call mh_eeq(n,at,xyz,0.0_wp,chrgeq,modh%kq,hess)
   endif

end subroutine mh_lindh

pure subroutine mh_lindh_stretch(n,at,xyz,hess,kr,kd,s6,aav,rav,dav,lcutoff,rcut)
   use xtb_mctc_constants
   use xtb_mctc_convert

   implicit none

   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: kr
   real(wp),intent(in)    :: kd
   real(wp),intent(in)    :: s6
   real(wp),intent(in)    :: aav(3,3)
   real(wp),intent(in)    :: rav(3,3)
   real(wp),intent(in)    :: dav(3,3)
   logical, intent(out)   :: lcutoff(n,n)
   real(wp),intent(in)    :: rcut

   integer  :: i,ir,j,jr
   real(wp) :: xij,yij,zij,rij2,r0,d0
   real(wp) :: alpha,gmm
   real(wp) :: c6i,c6j,c6ij,rv
   real(wp) :: hxx,hxy,hxz,hyy,hyz,hzz
   real(wp) :: vdw(3,3)

!! ------------------------------------------------------------------------
!  Hessian for stretch
!! ------------------------------------------------------------------------
   stretch_iAt: do i = 1, n
      ir = itabrow(at(i))

      stretch_jAt: do j = 1, i-1
         jr=itabrow(at(j))

         ! save for later
         lcutoff(i,j) = rcutoff(xyz,i,j,rcut)
         lcutoff(j,i) = lcutoff(i,j)

         xij=xyz(1,i)-xyz(1,j)
         yij=xyz(2,i)-xyz(2,j)
         zij=xyz(3,i)-xyz(3,j)
         rij2 = xij**2 + yij**2 + zij**2
         r0 = rav(ir,jr)
         d0 = dav(ir,jr)
         alpha=aav(ir,jr)

         !cccccc vdwx ccccccccccccccccccccccccccccccccc
         c6i=c6(at(i))
         c6j=c6(at(j))
         c6ij=sqrt(c6i*c6j)
         rv=(vander(at(i))+vander(at(j)))*aatoau

         call getvdwxx(xij, yij, zij, c6ij, s6, rv, vdw(1,1))
         call getvdwxy(xij, yij, zij, c6ij, s6, rv, vdw(1,2))
         call getvdwxy(xij, zij, yij, c6ij, s6, rv, vdw(1,3))
         call getvdwxx(yij, xij, zij, c6ij, s6, rv, vdw(2,2))
         call getvdwxy(yij, zij, xij, c6ij, s6, rv, vdw(2,3))
         call getvdwxx(zij, xij, yij, c6ij, s6, rv, vdw(3,3))
         !cccccc ende vdwx ccccccccccccccccccccccccccccccc

         gmm = kr*fk_lindh(alpha,r0,rij2) &
            + kr*kd * fk_vdw(4.0_wp,d0,rij2)

         !gmm = max(gmm,min_fk)

         hxx=gmm*xij*xij/rij2-vdw(1,1)
         hxy=gmm*xij*yij/rij2-vdw(1,2)
         hxz=gmm*xij*zij/rij2-vdw(1,3)
         hyy=gmm*yij*yij/rij2-vdw(2,2)
         hyz=gmm*yij*zij/rij2-vdw(2,3)
         hzz=gmm*zij*zij/rij2-vdw(3,3)

         ! save diagonal elements for atom i
         hess(ind(1,i,1,i))=hess(ind(1,i,1,i))+hxx
         hess(ind(2,i,1,i))=hess(ind(2,i,1,i))+hxy
         hess(ind(2,i,2,i))=hess(ind(2,i,2,i))+hyy
         hess(ind(3,i,1,i))=hess(ind(3,i,1,i))+hxz
         hess(ind(3,i,2,i))=hess(ind(3,i,2,i))+hyz
         hess(ind(3,i,3,i))=hess(ind(3,i,3,i))+hzz
         ! save elements between atom i and atom j
         hess(ind(1,i,1,j))=hess(ind(1,i,1,j))-hxx
         hess(ind(1,i,2,j))=hess(ind(1,i,2,j))-hxy
         hess(ind(1,i,3,j))=hess(ind(1,i,3,j))-hxz
         hess(ind(2,i,1,j))=hess(ind(2,i,1,j))-hxy
         hess(ind(2,i,2,j))=hess(ind(2,i,2,j))-hyy
         hess(ind(2,i,3,j))=hess(ind(2,i,3,j))-hyz
         hess(ind(3,i,1,j))=hess(ind(3,i,1,j))-hxz
         hess(ind(3,i,2,j))=hess(ind(3,i,2,j))-hyz
         hess(ind(3,i,3,j))=hess(ind(3,i,3,j))-hzz
         ! save diagonal elements for atom j
         hess(ind(1,j,1,j))=hess(ind(1,j,1,j))+hxx
         hess(ind(2,j,1,j))=hess(ind(2,j,1,j))+hxy
         hess(ind(2,j,2,j))=hess(ind(2,j,2,j))+hyy
         hess(ind(3,j,1,j))=hess(ind(3,j,1,j))+hxz
         hess(ind(3,j,2,j))=hess(ind(3,j,2,j))+hyz
         hess(ind(3,j,3,j))=hess(ind(3,j,3,j))+hzz

      end do stretch_jAt
   end do stretch_iAt

end subroutine mh_lindh_stretch

pure subroutine mh_lindh_bend(n,at,xyz,hess,kf,kd,aav,rav,dav,lcutoff)
   use xtb_mctc_constants

   implicit none

   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: kf
   real(wp),intent(in)    :: kd
   real(wp),intent(in)    :: aav(3,3)
   real(wp),intent(in)    :: rav(3,3)
   real(wp),intent(in)    :: dav(3,3)
   logical, intent(in)    :: lcutoff(n,n)

   integer  :: i,ir,j,jr,m,mr,ic,jc,ii
   real(wp),parameter :: rzero = 1.0e-10_wp
   real(wp) :: xij,yij,zij,rij2,rrij,r1
   real(wp) :: xmi,ymi,zmi,rmi2,rmi,r0mi,ami,d0mj,gmi
   real(wp) :: xmj,ymj,zmj,rmj2,rmj,r0mj,amj,d0mi,gmj
   real(wp) :: test,gij,rl2,rl,rmidotrmj
   real(wp) :: sinphi,cosphi,costhetax,costhetay,costhetaz
   real(wp) :: alpha
   real(wp) :: si(3),sj(3),sm(3),x(2),y(2),z(2)

!! ------------------------------------------------------------------------
!  Hessian for bending
!! ------------------------------------------------------------------------
   bend_mAt: do m = 1, n
      mr=itabrow(at(m))
      bend_iAt: do i = 1, n
         if (i.eq.m) cycle bend_iAt
         ir=itabrow(at(i))
         if(lcutoff(i,m)) cycle bend_iAt

         xmi=(xyz(1,i)-xyz(1,m))
         ymi=(xyz(2,i)-xyz(2,m))
         zmi=(xyz(3,i)-xyz(3,m))
         rmi2 = xmi**2 + ymi**2 + zmi**2
         rmi=sqrt(rmi2)
         r0mi=rav(mr,ir)
         d0mi=dav(mr,ir)
         ami=aav(mr,ir)

         bend_jAt: do j = 1, i-1
            if (j.eq.m) cycle bend_jAt
            jr=itabrow(at(j))
            if(lcutoff(j,i)) cycle bend_jAt
            if(lcutoff(j,m)) cycle bend_jAt

            xmj=(xyz(1,j)-xyz(1,m))
            ymj=(xyz(2,j)-xyz(2,m))
            zmj=(xyz(3,j)-xyz(3,m))
            rmj2 = xmj**2 + ymj**2 + zmj**2
            rmj=sqrt(rmj2)
            r0mj=rav(mr,jr)
            d0mj=dav(mr,jr)
            amj=aav(mr,jr)

            ! test if zero angle
            test=xmi*xmj+ymi*ymj+zmi*zmj
            test=test/(rmi*rmj)
            if (abs(test-1.0_wp).lt.1.0e-12_wp) cycle bend_jAt

            xij=(xyz(1,j)-xyz(1,i))
            yij=(xyz(2,j)-xyz(2,i))
            zij=(xyz(3,j)-xyz(3,i))
            rij2 = xij**2 + yij**2 + zij**2
            rrij=sqrt(rij2)

            gmi = fk_lindh(ami,r0mi,rmi2) &
                + 0.5_wp*kd * fk_vdw(4.0_wp,d0mi,rmi2)
            gmj = fk_lindh(amj,r0mj,rmj2) &
                + 0.5_wp*kd * fk_vdw(4.0_wp,d0mj,rmj2)

            gij = kf*gmi*gmj

            rl2=(ymi*zmj-zmi*ymj)**2+(zmi*xmj-xmi*zmj)**2+(xmi*ymj-ymi*xmj)**2

            if(rl2.lt.1.e-14_wp) then
               rl=0.0_wp
            else
               rl=sqrt(rl2)
            end if

            !gij = max(gij,min_fk)

            if ((rmj.gt.rzero).and.(rmi.gt.rzero).and.(rrij.gt.rzero)) then
               sinphi=rl/(rmj*rmi)
               rmidotrmj=xmi*xmj+ymi*ymj+zmi*zmj
               cosphi=rmidotrmj/(rmj*rmi)
               ! none linear case
               if (sinphi.gt.rzero) then
                  si(1)=(xmi/rmi*cosphi-xmj/rmj)/(rmi*sinphi)
                  si(2)=(ymi/rmi*cosphi-ymj/rmj)/(rmi*sinphi)
                  si(3)=(zmi/rmi*cosphi-zmj/rmj)/(rmi*sinphi)
                  sj(1)=(cosphi*xmj/rmj-xmi/rmi)/(rmj*sinphi)
                  sj(2)=(cosphi*ymj/rmj-ymi/rmi)/(rmj*sinphi)
                  sj(3)=(cosphi*zmj/rmj-zmi/rmi)/(rmj*sinphi)
                  sm(1)=-si(1)-sj(1)
                  sm(2)=-si(2)-sj(2)
                  sm(3)=-si(3)-sj(3)
                  do ic=1,3
                     do jc=1,3
                        if (m.gt.i) then
                           hess(ind(ic,m,jc,i))=hess(ind(ic,m,jc,i)) &
                              +gij*sm(ic)*si(jc)
                        else
                           hess(ind(ic,i,jc,m))=hess(ind(ic,i,jc,m)) &
                              +gij*si(ic)*sm(jc)
                        end if
                        if (m.gt.j) then
                           hess(ind(ic,m,jc,j))=hess(ind(ic,m,jc,j)) &
                              +gij*sm(ic)*sj(jc)
                        else
                           hess(ind(ic,j,jc,m))=hess(ind(ic,j,jc,m)) &
                              +gij*sj(ic)*sm(jc)
                        end if
                        if (i.gt.j) then
                           hess(ind(ic,i,jc,j))=hess(ind(ic,i,jc,j)) &
                              +gij*si(ic)*sj(jc)
                        else
                           hess(ind(ic,j,jc,i))=hess(ind(ic,j,jc,i)) &
                              +gij*sj(ic)*si(jc)
                        end if
                     end do
                  end do
                  do ic=1,3
                     do jc=1,ic
                        hess(ind(ic,i,jc,i))=hess(ind(ic,i,jc,i))+gij*si(ic)*si(jc)
                        hess(ind(ic,m,jc,m))=hess(ind(ic,m,jc,m))+gij*sm(ic)*sm(jc)
                        hess(ind(ic,j,jc,j))=hess(ind(ic,j,jc,j))+gij*sj(ic)*sj(jc)
                     end do
                  end do
               else
                  ! linear case
                  if ((abs(ymi).gt.rzero).or.(abs(xmi).gt.rzero)) then
                     x(1)=-ymi
                     y(1)=xmi
                     z(1)=0.0_wp
                     x(2)=-xmi*zmi
                     y(2)=-ymi*zmi
                     z(2)=xmi*xmi+ymi*ymi
                  else
                     x(1)=1.0_wp
                     y(1)=0.0_wp
                     z(1)=0.0_wp
                     x(2)=0.0_wp
                     y(2)=1.0_wp
                     z(2)=0.0_wp
                  end if
                  do ii=1,2
                     r1=sqrt(x(ii)**2+y(ii)**2+z(ii)**2)
                     costhetax=x(ii)/r1
                     costhetay=y(ii)/r1
                     costhetaz=z(ii)/r1
                     si(1)=-costhetax/rmi
                     si(2)=-costhetay/rmi
                     si(3)=-costhetaz/rmi
                     sj(1)=-costhetax/rmj
                     sj(2)=-costhetay/rmj
                     sj(3)=-costhetaz/rmj
                     sm(1)=-(si(1)+sj(1))
                     sm(2)=-(si(2)+sj(2))
                     sm(3)=-(si(3)+sj(3))
                     !
                     do ic=1,3
                        do jc=1,3
                           if (m.gt.i) then
                              hess(ind(ic,m,jc,i))=hess(ind(ic,m,jc,i)) &
                                 +gij*sm(ic)*si(jc)
                           else
                              hess(ind(ic,i,jc,m))=hess(ind(ic,i,jc,m)) &
                                 +gij*si(ic)*sm(jc)
                           end if
                           if (m.gt.j) then
                              hess(ind(ic,m,jc,j))=hess(ind(ic,m,jc,j)) &
                                 +gij*sm(ic)*sj(jc)
                           else
                              hess(ind(ic,j,jc,m))=hess(ind(ic,j,jc,m)) &
                                 +gij*sj(ic)*sm(jc)
                           end if
                           if (i.gt.j) then
                              hess(ind(ic,i,jc,j))=hess(ind(ic,i,jc,j)) &
                                 +gij*si(ic)*sj(jc)
                           else
                              hess(ind(ic,j,jc,i))=hess(ind(ic,j,jc,i)) &
                                 +gij*sj(ic)*si(jc)
                           end if
                        end do
                     end do
                     do ic=1,3
                        do jc=1,ic
                           hess(ind(ic,i,jc,i))=hess(ind(ic,i,jc,i)) &
                              +gij*si(ic)*si(jc)
                           hess(ind(ic,m,jc,m))=hess(ind(ic,m,jc,m)) &
                              +gij*sm(ic)*sm(jc)
                           hess(ind(ic,j,jc,j))=hess(ind(ic,j,jc,j)) &
                              +gij*sj(ic)*sj(jc)
                        end do
                     end do
                  end do

               end if
            end if

         end do bend_jAt
      end do bend_iAt
   end do bend_mAt

end subroutine mh_lindh_bend

subroutine mh_lindh_torsion(n,at,xyz,hess,kt,kd,aav,rav,dav,lcutoff)
   use xtb_mctc_constants

   implicit none

   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: kt
   real(wp),intent(in)    :: kd
   real(wp),intent(in)    :: aav(3,3)
   real(wp),intent(in)    :: rav(3,3)
   real(wp),intent(in)    :: dav(3,3)
   logical, intent(in)    :: lcutoff(n,n)

   integer  :: i,ir,j,jr,k,kr,l,lr,ic,jc,ij,kl
!  allow only angles in the range of 35-145
   real(wp),parameter :: a35 = (35.0d0/180.d0)* pi
   real(wp),parameter :: cosfi_max=cos(a35)
   real(wp) :: txyz(3,4),c(3,4)
   real(wp) :: rij(3),rij0,aij,rij2,d0ij,gij
   real(wp) :: rjk(3),rjk0,ajk,rjk2,d0jk,gjk
   real(wp) :: rkl(3),rkl0,akl,rkl2,d0kl,gkl
   real(wp) :: cosfi2,cosfi3,cosfi4
   real(wp) :: beta,tij,tau,dum(3,4,3,4)
   real(wp) :: si(3),sj(3),sk(3),sl(3)

!! ------------------------------------------------------------------------
!  Hessian for torsion
!! ------------------------------------------------------------------------
   torsion_jAt: do j = 1,n
      jr=itabrow(at(j))
      txyz(:,2)=xyz(:,j)
      torsion_kAt: do k = 1, n
         if (k.eq.j) cycle torsion_kAt
         kr=itabrow(at(k))
         if(lcutoff(k,j)) cycle torsion_kAt
         txyz(:,3) = xyz(:,k)
         torsion_iAt: do i = 1, n
            ij=n*(j-1)+i
            if (i.eq.j) cycle torsion_iAt
            if (i.eq.k) cycle torsion_iAt
            ir=itabrow(at(i))
            if(lcutoff(i,k)) cycle torsion_iAt
            if(lcutoff(i,j)) cycle torsion_iAt

            txyz(:,1)=xyz(:,i)
            torsion_lAt: do l = 1, n
               kl=n*(k-1)+l
               if (ij.le.kl) cycle torsion_lAt
               if (l.eq.i)   cycle torsion_lAt
               if (l.eq.j)   cycle torsion_lAt
               if (l.eq.k)   cycle torsion_lAt
               lr=itabrow(at(l))
!
               if(lcutoff(l,i)) cycle torsion_lAt
               if(lcutoff(l,k)) cycle torsion_lAt
               if(lcutoff(l,j)) cycle torsion_lAt

               txyz(:,4)=xyz(:,l)

               rij=xyz(:,i)-xyz(:,j)
               d0ij=dav(ir,jr)
               rij0=rav(ir,jr)
               aij =aav(ir,jr)

               rjk=xyz(:,j)-xyz(:,k)
               d0jk=dav(jr,kr)
               rjk0=rav(jr,kr)
               ajk =aav(jr,kr)

               rkl=xyz(:,k)-xyz(:,l)
               d0kl=dav(kr,lr)
               rkl0=rav(kr,lr)
               akl =aav(kr,lr)

               rij2=sum(rij**2)
               rjk2=sum(rjk**2)
               rkl2=sum(rjk**2)

               cosfi2=dot_product(rij,rjk)/sqrt(rij2*rjk2)
               if (abs(cosfi2).gt.cosfi_max) cycle
               cosfi3=dot_product(rkl,rjk)/sqrt(rkl2*rjk2)
               if (abs(cosfi3).gt.cosfi_max) cycle

               gij = fk_lindh(aij,rij0,rij2) &
                  + 0.5_wp*kd * fk_vdw(4.0_wp,d0ij,rij2)
               gjk = fk_lindh(ajk,rjk0,rjk2) &
                  + 0.5_wp*kd * fk_vdw(4.0_wp,d0jk,rjk2)
               gkl = fk_lindh(akl,rkl0,rkl2) &
                  + 0.5_wp*kd * fk_vdw(4.0_wp,d0kl,rkl2)

               tij = kt * gij*gjk*gkl

               !tij = max(tij,10*min_fk)

               !call trsn2(txyz,tau,c)
               Call Trsn(txyz,4,Tau,C,.False.,.False.,'        ', &
      &                  Dum,.False.)
               si = c(:,1)
               sj = c(:,2)
               sk = c(:,3)
               sl = c(:,4)

               ! off diagonal block
               do ic=1,3
                  do jc=1,3
                     hess(ind(ic,i,jc,j))=hess(ind(ic,i,jc,j))+tij*si(ic) * sj(jc)
                     hess(ind(ic,i,jc,k))=hess(ind(ic,i,jc,k))+tij*si(ic) * sk(jc)
                     hess(ind(ic,i,jc,l))=hess(ind(ic,i,jc,l))+tij*si(ic) * sl(jc)
                     hess(ind(ic,j,jc,k))=hess(ind(ic,j,jc,k))+tij*sj(ic) * sk(jc)
                     hess(ind(ic,j,jc,l))=hess(ind(ic,j,jc,l))+tij*sj(ic) * sl(jc)
                     hess(ind(ic,k,jc,l))=hess(ind(ic,k,jc,l))+tij*sk(ic) * sl(jc)
                  end do
               end do

               ! diagonal block
               do ic=1,3
                  do jc=1,ic
                     hess(ind(ic,i,jc,i))=hess(ind(ic,i,jc,i))+tij*si(ic) * si(jc)
                     hess(ind(ic,j,jc,j))=hess(ind(ic,j,jc,j))+tij*sj(ic) * sj(jc)
                     hess(ind(ic,k,jc,k))=hess(ind(ic,k,jc,k))+tij*sk(ic) * sk(jc)
                     hess(ind(ic,l,jc,l))=hess(ind(ic,l,jc,l))+tij*sl(ic) * sl(jc)
                  end do
               end do

            end do torsion_lAt
         end do torsion_iAt
      end do torsion_kAt
   end do torsion_jAt

end subroutine mh_lindh_torsion

pure subroutine mh_lindh_outofp(n,at,xyz,hess,ko,kd,aav,rav,dav,lcutoff)
   use xtb_mctc_constants

   implicit none

   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: hess((3*n)*(3*n+1)/2)
   real(wp),intent(in)    :: ko
   real(wp),intent(in)    :: kd
   real(wp),intent(in)    :: aav(3,3)
   real(wp),intent(in)    :: rav(3,3)
   real(wp),intent(in)    :: dav(3,3)
   logical, intent(in)    :: lcutoff(n,n)

   integer  :: i,ir,j,jr,k,kr,l,lr,ic,jc
   real(wp) :: txyz(3,4),c(3,4)
   real(wp) :: rij(3),rij0,aij,rij2,gij,d0ij
   real(wp) :: rik(3),rik0,aik,rik2,gik,d0ik
   real(wp) :: ril(3),ril0,ail,ril2,gil,d0il
   real(wp) :: cosfi2,cosfi3,cosfi4
   real(wp) :: beta,tij,tau
   real(wp) :: si(3),sj(3),sk(3),sl(3)

!! ------------------------------------------------------------------------
!  Hessian for out-of-plane
!! ------------------------------------------------------------------------
   outofplane_iAt: do i = 1, n
      ir = itabrow(at(i))
      txyz(:,4) = xyz(:,i)
      outofplane_jAt: do j = 1, n
         if (j.eq.i) cycle outofplane_jAt
         if(lcutoff(j,i)) cycle outofplane_jAt
         jr = itabrow(at(j))
         txyz(:,1) = xyz(:,j)
         outofplane_kAt: do k = 1, n
            if (i.eq.k) cycle outofplane_kAt
            if (j.eq.k) cycle outofplane_kat
            if(lcutoff(k,i)) cycle outofplane_kAt
            if(lcutoff(k,j)) cycle outofplane_kAt
            kr = itabrow(at(k))
            txyz(:,2) = xyz(:,k)
            outofplane_lAt: do l = 1, n
               lr = itabrow(at(l))
               txyz(:,3) = xyz(:,l)
               if (l.eq.i) cycle outofplane_lAt
               if (l.eq.j) cycle outofplane_lAt
               if (l.eq.k) cycle outofplane_lAt
               if(lcutoff(l,i)) cycle outofplane_lAt
               if(lcutoff(l,k)) cycle outofplane_lAt
               if(lcutoff(l,j)) cycle outofplane_lAt

               rij=xyz(:,i)-xyz(:,j)
               d0ij=dav(ir,jr)
               rij0=rav(ir,jr)
               aij =aav(ir,jr)

               rik=xyz(:,i)-xyz(:,k)
               d0ik=dav(ir,kr)
               rik0=rav(ir,kr)
               aik =aav(ir,kr)

               ril=xyz(:,i)-xyz(:,l)
               d0il=dav(ir,lr)
               ril0=rav(ir,lr)
               ail =aav(ir,lr)

               rij2=sum(rij**2)
               rik2=sum(rik**2)
               ril2=sum(ril**2)

               cosfi2=dot_product(rij,rik)/sqrt(rij2*rik2)
               if (abs(abs(cosfi2)-1.0_wp).lt.1.0e-1_wp) cycle
               cosfi3=dot_product(rij,ril)/sqrt(rij2*ril2)
               if (abs(abs(cosfi3)-1.0_wp).lt.1.0e-1_wp) cycle
               cosfi4=dot_product(rik,ril)/sqrt(rik2*ril2)
               if (abs(abs(cosfi4)-1.0_wp).lt.1.0e-1_wp) cycle

               gij = fk_lindh(aij,rij0,rij2) &
                  + 0.5_wp*kd * fk_vdw(4.0_wp,d0ij,rij2)
               gik = fk_lindh(aik,rik0,rik2) &
                  + 0.5_wp*kd * fk_vdw(4.0_wp,d0ik,rik2)
               gil = fk_lindh(ail,ril0,ril2) &
                  + 0.5_wp*kd * fk_vdw(4.0_wp,d0il,ril2)

               tij = ko * gij*gik*gil

               !tij = max(tij,10*min_fk)

               call outofp2(xyz,tau,c)
               If (abs(tau).gt.45.0d0*(pi/180.d0)) cycle

               si = c(:,4)
               sj = c(:,1)
               sk = c(:,2)
               sl = c(:,3)

               ! off diagonal block
               do ic=1,3
                  do jc=1,3
                     hess(ind(ic,i,jc,j))=hess(ind(ic,i,jc,j))+tij*si(ic) * sj(jc)
                     hess(ind(ic,i,jc,k))=hess(ind(ic,i,jc,k))+tij*si(ic) * sk(jc)
                     hess(ind(ic,i,jc,l))=hess(ind(ic,i,jc,l))+tij*si(ic) * sl(jc)
                     hess(ind(ic,j,jc,k))=hess(ind(ic,j,jc,k))+tij*sj(ic) * sk(jc)
                     hess(ind(ic,j,jc,l))=hess(ind(ic,j,jc,l))+tij*sj(ic) * sl(jc)
                     hess(ind(ic,k,jc,l))=hess(ind(ic,k,jc,l))+tij*sk(ic) * sl(jc)
                  end do
               end do

               ! diagonal block
               do ic=1,3
                  do jc=1,ic
                     hess(ind(ic,i,jc,i))=hess(ind(ic,i,jc,i))+tij*si(ic) * si(jc)
                     hess(ind(ic,j,jc,j))=hess(ind(ic,j,jc,j))+tij*sj(ic) * sj(jc)
                     hess(ind(ic,k,jc,k))=hess(ind(ic,k,jc,k))+tij*sk(ic) * sk(jc)
                     hess(ind(ic,l,jc,l))=hess(ind(ic,l,jc,l))+tij*sl(ic) * sl(jc)
                  end do
               end do

            enddo outofplane_lAt
         enddo outofplane_kAt
      enddo outofplane_jAt
   enddo outofplane_iAt

end subroutine mh_lindh_outofp

pure elemental function ixyz(i,iatom)
   integer :: ixyz
   integer,intent(in) :: i,iatom
   ixyz = (iatom-1)*3 + i
end function ixyz
pure elemental function jnd(i,j)
   integer :: jnd
   integer,intent(in) :: i,j
   jnd = i*(i-1)/2 +j
end function jnd
pure elemental function ind(i,iatom,j,jatom)
   integer :: ind
   integer,intent(in) :: i,iatom,j,jatom
   ind=jnd(max(ixyz(i,iatom),ixyz(j,jatom)),min(ixyz(i,iatom),ixyz(j,jatom)))
end function ind

pure function rcutoff(xyz,katom,latom,rcut)
   implicit none
   logical  :: rcutoff
   real(wp),intent(in) :: xyz(3,*)
   real(wp),intent(in) :: rcut
   real(wp) :: rkl(3),rkl2
   integer, intent(in) :: katom,latom
   rcutoff=.false.
   rkl=xyz(:,kAtom)-xyz(:,lAtom)
   rkl2 = sum(rkl**2)
   if(rkl2.gt.rcut) rcutoff=.true.
end function rcutoff

pure elemental function itabrow(i)
   integer :: itabrow
   integer,intent(in) :: i

   itabrow=0
   if (i.gt. 0 .and. i.le. 2) then
      itabrow=1
   else if (i.gt. 2 .and. i.le.10) then
      itabrow=2
   else if (i.gt.10 .and. i.le.18) then
      itabrow=3
   else if (i.gt.18 .and. i.le.36) then
      itabrow=3
   else if (i.gt.36 .and. i.le.54) then
      itabrow=3
   else if (i.gt.54 .and. i.le.86) then
      itabrow=3
   else if (i.gt.86) then
      itabrow=3
   end if

   return
end function itabrow

pure subroutine outofp2(xyz,teta,bt)
   use xtb_mctc_constants
   implicit none
   real(wp),intent(out) :: teta
   real(wp),intent(out) :: bt(3,4)
   real(wp),intent(in)  :: xyz(3,4)
   real(wp) :: r1(3),r2(3),r3(3)
   real(wp) :: q41,q42,q43,e41(3),e42(3),e43(3)
   real(wp) :: cosfi1,fi1,dfi1,cosfi2,fi2,dfi2,cosfi3,fi3,dfi3
   real(wp) :: c14(3,3),br14(3,3)
   real(wp) :: r42(3),r43(3)
   integer  :: ix,iy,iz
!  4 -> 1 (bond)
   r1=xyz(:,1)-xyz(:,4)
   q41=norm2(r1)
   e41 = r1 / q41
!  4 -> 2 (bond in plane)
   r2=xyz(:,2)-xyz(:,4)
   q42=norm2(r2)
   e42 = r2 / q42
!  4 -> 3 (bond in plane)
   r3=xyz(:,3)-xyz(:,4)
   q43=norm2(r3)
   e43 = r3 / q43
!
!  get the angle between e43 and e42
!
   cosfi1 = dot_product(e43,e42)

   fi1=acos(cosfi1)
   dfi1 = 180.d0 * fi1 / pi
!
!  dirty exit! this happens when an earlier structure is ill defined.
!
   if (abs(fi1-pi).lt.1.0d-13) then
      teta=0.0_wp
      bt = 0.0_wp
      return
   end if
!
!  get the angle between e41 and e43
!
   cosfi2 = dot_product(e41,e43)

   fi2=acos(cosfi2)
   dfi2 = 180.d0 * fi2 / pi
!
!  get the angle between e41 and e42
!
   cosfi3 = dot_product(e41,e42)

   fi3=acos(cosfi3)
   dfi3 = 180.d0 * fi3 / pi
!
!  the first two centers are trivially
!
   c14(:,1) = xyz(:,1)
   c14(:,2) = xyz(:,4)
!
!  the 3rd is
!
   r42=xyz(:,2)-xyz(:,4)
   r43=xyz(:,3)-xyz(:,4)
   c14(1,3)=r42(2)*r43(3)-r42(3)*r43(2)
   c14(2,3)=r42(3)*r43(1)-r42(1)*r43(3)
   c14(3,3)=r42(1)*r43(2)-r42(2)*r43(1)
!
!  exit if 2-3-4 are collinear
!  (equivalent to the above check, but this is more concrete)
!
   if ((c14(1,3)**2+c14(2,3)**2+c14(3,3)**2).lt.1.0d-10) then
      teta=0.0d0
      bt = 0.0_wp
      return
   end if
   c14(1,3)=c14(1,3)+xyz(1,4)
   c14(2,3)=c14(2,3)+xyz(2,4)
   c14(3,3)=c14(3,3)+xyz(3,4)

   call bend2(c14,teta,br14)

   teta = teta - 0.5_wp*pi
!
!--compute the wdc matrix
!
   do ix = 1, 3
      iy = mod(ix+1, 4)+(ix+1)/4
      iz = mod(iy+1, 4)+(iy+1)/4

      bt(ix,1) = - br14(ix,1)
      bt(ix,2) =   r43(iz)*br14(iy,3) - r43(iy)*br14(iz,3)
      bt(ix,3) = - r42(iz)*br14(iy,3) + r42(iy)*br14(iz,3)

      bt(ix,4) = - (bt(ix,1)+bt(ix,2)+bt(ix,3))

   end do

   bt = -bt
end subroutine outofp2

pure subroutine trsn2(xyz,tau,bt)
   use xtb_mctc_constants
   implicit none
   real(wp),intent(out) :: bt(3,4)
   real(wp),intent(out) :: tau
   real(wp),intent(in)  :: xyz(3,4)
   real(wp) :: rij(3),rij1,brij(3,2)
   real(wp) :: rjk(3),rjk1,brjk(3,2)
   real(wp) :: rkl(3),rkl1,brkl(3,2)
   real(wp) :: bf2(3,3),fi2,sinfi2,cosfi2
   real(wp) :: bf3(3,3),fi3,sinfi3,cosfi3
   real(wp) :: costau,sintau
   integer  :: ix,iy,iz
   call strtch2(xyz(1,1),rij1,brij)
   call strtch2(xyz(1,2),rjk1,brjk)
   call strtch2(xyz(1,3),rkl1,brkl)
   call bend2(xyz(1,1),fi2,bf2)
   sinfi2=sin(fi2)
   cosfi2=cos(fi2)
   call bend2(xyz(1,2),fi3,bf3)
   sinfi3=sin(fi3)
   cosfi3=cos(fi3)
   costau = ( ( brij(2,1)*brjk(3,2) - brij(3,1)*brjk(2,2) ) * &
               ( brjk(2,1)*brkl(3,2) - brjk(3,1)*brkl(2,2) ) + &
               ( brij(3,1)*brjk(1,2) - brij(1,1)*brjk(3,2) ) * &
               ( brjk(3,1)*brkl(1,2) - brjk(1,1)*brkl(3,2) ) + &
               ( brij(1,1)*brjk(2,2) - brij(2,1)*brjk(1,2) ) * &
               ( brjk(1,1)*brkl(2,2) - brjk(2,1)*brkl(1,2) ) ) &
             / (sinfi2*sinfi3)
   sintau = ( brij(1,2) * (brjk(2,1)*brkl(3,2)-brjk(3,1)*brkl(2,2)) &
             + brij(2,2) * (brjk(3,1)*brkl(1,2)-brjk(1,1)*brkl(3,2)) &
             + brij(3,2) * (brjk(1,1)*brkl(2,2)-brjk(2,1)*brkl(1,2)) ) &
             / (sinfi2*sinfi3)
   tau = atan2(sintau,costau)
   if (abs(tau).eq.pi) tau=pi
   do ix = 1, 3
      iy=ix+1
      if (iy.gt.3) iy=iy-3
      iz=iy+1
      if (iz.gt.3) iz=iz-3
      bt(ix,1) = (brij(iy,2)*brjk(iz,2)-brij(iz,2)*brjk(iy,2)) &
         &           / (rij1*sinfi2**2)
      bt(ix,4) = (brkl(iy,1)*brjk(iz,1)-brkl(iz,1)*brjk(iy,1)) &
         &           / (rkl1*sinfi3**2)
      bt(ix,2) = -( (rjk1-rij1*cosfi2) * bt(ix,1) &
         &             +         rkl1*cosfi3  * bt(ix,4))/rjk1
      bt(ix,3) = - ( bt(ix,1)+bt(ix,2)+bt(ix,4))
   end do
end subroutine trsn2
pure subroutine strtch2(xyz,avst,b)
   implicit none
   real(wp),intent(out) :: b(3,2)
   real(wp),intent(in)  :: xyz(3,2)
   real(wp) :: r(3)
   real(wp) :: rr
   real(wp),intent(out) :: avst
   r=xyz(:,2)-xyz(:,1)
   rr=norm2(r)
   avst=rr
   b(:,1)=-r/rr
   b(:,2)=-b(:,1)
end subroutine strtch2
pure subroutine bend2(xyz,fir,bf)
   use xtb_mctc_constants
   implicit none
   real(wp),intent(out) :: bf(3,3)
   real(wp),intent(in)  :: xyz(3,3)
   real(wp) :: brij(3,2)
   real(wp) :: brjk(3,2)
   real(wp) :: co,crap
   real(wp),intent(out) :: fir
   real(wp) :: si
   real(wp) :: rij1,rjk1
   integer  :: i
   call strtch2(xyz(1,1),rij1,brij)
   call strtch2(xyz(1,2),rjk1,brjk)
   co=0.0_wp
   crap=0.0_wp
   do i = 1, 3
      co=co+brij(i,1)*brjk(i,2)
      crap=crap+(brjk(i,2)+brij(i,1))**2
   end do
   if (sqrt(crap).lt.1.0d-6) then
      fir=pi-asin(sqrt(crap))
      si=sqrt(crap)
   else
      fir=acos(co)
      si=sqrt(1.0_wp-co**2)
   end if
   if (abs(fir-pi).lt.1.0d-13) then
      fir=pi
      return
   end if
   do i = 1, 3
      bf(i,1)= (co*brij(i,1)-brjk(i,2))/(si*rij1)
      bf(i,3)= (co*brjk(i,2)-brij(i,1))/(si*rjk1)
      bf(i,2)=-(bf(i,1)+bf(i,3))
   end do
end subroutine bend2

pure elemental subroutine getvdwxy(rx,ry,rz, c66, s6,r0, vdw)
   !cc Ableitung nach rx und ry
   implicit none
   real(wp),intent(in)  :: rx,ry,rz,c66,s6,r0
   real(wp),intent(out) :: vdw
   real(wp) :: t1,t2,t3,t4,t5,t6,t7,t11,t12,t16,t17,t25,t26,t35
   real(wp) :: t40,t41,t43,t44,t56,avdw
   ! write(*,*) 's6:', s6
   avdw=20.0
   t1 = s6 * C66
   t2 = rx ** 2
   t3 = ry ** 2
   t4 = rz ** 2
   t5 = t2 + t3 + t4
   t6 = t5 ** 2
   t7 = t6 ** 2
   t11 = sqrt(t5)
   t12 = 0.1D1 / r0
   t16 = exp(-avdw * (t11 * t12 - 0.1D1))
   t17 = 0.1D1 + t16
   t25 = t17 ** 2
   t26 = 0.1D1 / t25
   t35 = 0.1D1 / t7
   t40 = avdw ** 2
   t41 = r0 ** 2
   t43 = t40 / t41
   t44 = t16 ** 2
   t56 = -0.48D2 * t1 / t7 / t5 / t17 * rx * ry + 0.13D2 * t1 / t11 / &
      t7 * t26 * rx * avdw * t12 * ry * t16 - 0.2D1 * t1 * t35 / t25 / &
      t17 * t43 * rx * t44 * ry + t1 * t35 * t26 * t43 * rx * ry * t16
   vdw=t56
end subroutine getvdwxy

pure elemental subroutine getvdwxx(rx, ry, rz, c66, s6, r0, vdw)
   !cc Ableitung nach rx und rx
   implicit none
   real(wp),intent(in)  :: rx,ry,rz,c66,s6,r0
   real(wp),intent(out) :: vdw
   real(wp) :: t1,t2,t3,t4,t5,t6,t7,t10,t11,t15,t16,t17,t24,t25,t29
   real(wp) :: t33,t41,t42,t44,t45,t62,avdw
   avdw=20.0
   ! write(*,*) 's6:', s6
   t1 = s6 * C66
   t2 = rx ** 2
   t3 = ry ** 2
   t4 = rz ** 2
   t5 = t2 + t3 + t4
   t6 = t5 ** 2
   t7 = t6 ** 2
   t10 = sqrt(t5)
   t11 = 0.1D1 / r0
   t15 = exp(-avdw * (t10 * t11 - 0.1D1))
   t16 = 0.1D1 + t15
   t17 = 0.1D1 / t16
   t24 = t16 ** 2
   t25 = 0.1D1 / t24
   t29 = t11 * t15
   t33 = 0.1D1 / t7
   t41 = avdw ** 2
   t42 = r0 ** 2
   t44 = t41 / t42
   t45 = t15 ** 2
   t62 = -0.48D2 * t1 / t7 / t5 * t17 * t2 + 0.13D2 * t1 / t10 / t7 * &
      t25 * t2 * avdw * t29 + 0.6D1 * t1 * t33 * t17 - 0.2D1 * t1 * t33 &
      / t24 / t16 * t44 * t2 * t45 - t1 / t10 / t6 / t5 * t25 * avdw * &
      t29 + t1 * t33 * t25 * t44 * t2 * t15
   vdw=t62
end subroutine getvdwxx

pure elemental function fk_lindh(alpha,r0,r2) result(gmm)
   implicit none
   real(wp),intent(in) :: alpha,r0,r2
   real(wp) :: gmm
   gmm = exp(alpha*(r0**2 - r2))
end function fk_lindh

pure elemental function fk_swart(alpha,r0,r2) result(gmm)
   implicit none
   real(wp),intent(in) :: alpha,r0,r2
   real(wp) :: gmm
   gmm = exp(-alpha*(sqrt(r2)/r0 - 1.0_wp))
end function fk_swart

pure elemental function fk_vdw(alpha,r0,r2) result(gmm)
   implicit none
   real(wp),intent(in) :: alpha,r0,r2
   real(wp) :: gmm
   gmm = exp(-alpha*(r0 - sqrt(r2))**2)
end function fk_vdw

subroutine mh_eeq(n,at,xyz,chrg,chrgeq,kq,hess)
   use xtb_type_param
   implicit none

!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: n                ! number of atoms
   integer, intent(in)    :: at(n)            ! ordinal numbers
   real(wp),intent(in)    :: xyz(3,n)         ! geometry
   real(wp),intent(in)    :: chrg             ! total charge
   real(wp),intent(in)    :: kq               ! scaling parameter
   type(chrg_parameter),intent(in) :: chrgeq  ! charge model
!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
   real(wp),intent(out)   :: hess((3*n)*(3*n+1)/2)
   real(wp),allocatable   :: hessian(:,:,:,:) ! molecular hessian of IES

!  π itself
   real(wp),parameter :: pi = 3.1415926535897932384626433832795029_wp
!  √π
   real(wp),parameter :: sqrtpi = sqrt(pi)
!  √(2/π)
   real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)
!
!! ------------------------------------------------------------------------
!  charge model
!! ------------------------------------------------------------------------
   integer  :: m ! dimension of the Lagrangian
   real(wp),allocatable :: Amat(:,:)
   real(wp),allocatable :: Xvec(:)
   real(wp),allocatable :: Ainv(:,:)
   real(wp),allocatable :: dAmat(:,:,:)
   real(wp),allocatable :: dqdr(:,:,:)

!! ------------------------------------------------------------------------
!  local variables
!! ------------------------------------------------------------------------
   integer  :: i,j,k,l
   real(wp) :: r,rij(3),r2
   real(wp) :: gamij,gamij2
   real(wp) :: arg,arg2,tmp,dtmp
   real(wp) :: lambda
   real(wp) :: es,expterm,erfterm
   real(wp) :: htmp,rxr(3,3)
   real(wp) :: rcovij,rr

!! ------------------------------------------------------------------------
!  scratch variables
!! ------------------------------------------------------------------------
   real(wp),allocatable :: alpha(:)
   real(wp),allocatable :: xtmp(:)
   real(wp),allocatable :: atmp(:,:)

!! ------------------------------------------------------------------------
!  Lapack work variables
!! ------------------------------------------------------------------------
   integer, allocatable :: ipiv(:)
   real(wp),allocatable :: temp(:)
   real(wp),allocatable :: work(:)
   integer  :: lwork
   integer  :: info
   real(wp) :: test(1)

!! ------------------------------------------------------------------------
!  initizialization
!! ------------------------------------------------------------------------
   m    = n+1
   allocate( ipiv(m), source = 0 )
   allocate( Amat(m,m), Xvec(m), alpha(n), dqdr(3,n,m), source = 0.0_wp )

!! ------------------------------------------------------------------------
!  set up the A matrix and X vector
!! ------------------------------------------------------------------------
!  αi -> alpha(i), ENi -> xi(i), κi -> kappa(i), Jii -> gam(i)
!  γij = 1/√(αi+αj)
!  Xi  = -ENi + κi·√CNi
!  Aii = Jii + 2/√π·γii
!  Aij = erf(γij·Rij)/Rij = 2/√π·F0(γ²ij·R²ij)
!! ------------------------------------------------------------------------
!  prepare some arrays
!$omp parallel default(none) &
!$omp shared(n,at,chrgeq) &
!$omp private(i) &
!$omp shared(Xvec,alpha)
!$omp do schedule(dynamic)
   do i = 1, n
      Xvec(i) = -chrgeq%en(i)
      alpha(i) = chrgeq%alpha(i)**2
   enddo
!$omp enddo
!$omp endparallel

!$omp parallel default(none) &
!$omp shared(n,at,xyz,chrgeq,alpha) &
!$omp private(i,j,r,gamij) &
!$omp shared(Amat)
!$omp do schedule(dynamic)
   ! prepare A matrix
   do i = 1, n
      ! EN of atom i
      do j = 1, i-1
         r = sqrt(sum((xyz(:,j) - xyz(:,i))**2))
         gamij = 1.0_wp/sqrt(alpha(i)+alpha(j))
         Amat(j,i) = erf(gamij*r)/r
         Amat(i,j) = Amat(j,i)
      enddo
      Amat(i,i) = chrgeq%gam(i) + sqrt2pi/sqrt(alpha(i))
   enddo
!$omp enddo
!$omp endparallel

!! ------------------------------------------------------------------------
!  solve the linear equations to obtain partial charges
!! ------------------------------------------------------------------------
   Amat(m,1:m) = 1.0_wp
   Amat(1:m,m) = 1.0_wp
   Amat(m,m  ) = 0.0_wp
   Xvec(m)     = chrg
   ! generate temporary copy
   allocate( Atmp(m,m), source = Amat )
   allocate( Xtmp(m),   source = Xvec )

   ! assume work space query, set best value to test after first dsysv call
   call dsysv('u', m, 1, Atmp, m, ipiv, Xtmp, m, test, -1, info)
   lwork = int(test(1))
   allocate( work(lwork), source = 0.0_wp )

   call dsysv('u',m,1,Atmp,m,ipiv,Xtmp,m,work,lwork,info)
   if(info > 0) call raise('E','(goedecker_solve) DSYSV failed')

   if(abs(sum(Xtmp(:n))-chrg) > 1.e-6_wp) &
      call raise('E','(goedecker_solve) charge constrain error')
   !print'(3f20.14)',Xtmp

!! ------------------------------------------------------------------------
!  calculate isotropic electrostatic (IES) energy
!! ------------------------------------------------------------------------
!  E = ∑i (ENi - κi·√CNi)·qi + ∑i (Jii + 2/√π·γii)·q²i
!      + ½ ∑i ∑j,j≠i qi·qj·2/√π·F0(γ²ij·R²ij)
!    = q·(½A·q - X)
!! ------------------------------------------------------------------------
!   work(:m) = Xvec
!   call dsymv('u',m,0.5_wp,Amat,m,Xtmp,1,-1.0_wp,work,1)
!   es = dot_product(Xtmp,work(:m))
!   energy = es + energy

!! ------------------------------------------------------------------------
!  calculate molecular gradient of the IES energy
!! ------------------------------------------------------------------------
!  dE/dRj -> g(:,j), ∂Xi/∂Rj -> -dcn(:,i,j), ½∂Aij/∂Rj -> dAmat(:,j,i)
!  dE/dR = (½∂A/∂R·q - ∂X/∂R)·q
!  ∂Aij/∂Rj = ∂Aij/∂Ri
!! ------------------------------------------------------------------------
   allocate( dAmat(3,n,m), source = 0.0_wp )
!$omp parallel default(none) &
!$omp shared(n,xyz,alpha,Amat,Xtmp) &
!$omp private(i,j,rij,r2,gamij,arg,dtmp) &
!$omp reduction(+:dAmat)
!$omp do schedule(dynamic)
   do i = 1, n
      do j = 1, i-1
         rij = xyz(:,i) - xyz(:,j)
         r2 = sum(rij**2)
         gamij = 1.0_wp/sqrt(alpha(i) + alpha(j))
         arg = gamij**2*r2
         dtmp = 2.0_wp*gamij*exp(-arg)/(sqrtpi*r2)-Amat(j,i)/r2
         dAmat(:,i,i) = +dtmp*rij*Xtmp(j) + dAmat(:,i,i)
         dAmat(:,j,j) = -dtmp*rij*Xtmp(i) + dAmat(:,j,j)
         dAmat(:,i,j) = +dtmp*rij*Xtmp(i)
         dAmat(:,j,i) = -dtmp*rij*Xtmp(j)
      enddo
   enddo
!$omp enddo
!$omp endparallel

!! ------------------------------------------------------------------------
!  invert the A matrix using a Bunch-Kaufman factorization
!  A⁻¹ = (L·D·L^T)⁻¹ = L^T·D⁻¹·L
!! ------------------------------------------------------------------------
   allocate( Ainv(m,m), source = Amat )

   ! assume work space query, set best value to test after first dsytrf call
   call dsytrf('L',m,Ainv,m,ipiv,test,-1,info)
   if (int(test(1)) > lwork) then
      deallocate(work)
      lwork=int(test(1))
      allocate( work(lwork), source = 0.0_wp )
   endif

   ! Bunch-Kaufman factorization A = L*D*L**T
   call dsytrf('L',m,Ainv,m,ipiv,work,lwork,info)
   if(info > 0)then
      call raise('E', '(goedecker_inversion) DSYTRF failed')
   endif

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Ainv matrix
   ! Ainv matrix is overwritten with lower triangular part of A⁻¹
   call dsytri('L',m,Ainv,m,ipiv,work,info)
   if (info > 0) then
      call raise('E', '(goedecker_inversion) DSYTRI failed')
   endif

   ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
   do i = 1, m
      do j = i+1, m
         Ainv(i,j)=Ainv(j,i)
      enddo
   enddo

!! ------------------------------------------------------------------------
!  calculate gradient of the partial charge w.r.t. the nuclear coordinates
!! ------------------------------------------------------------------------
   !call dsymm('r','l',3*n,m,-1.0_wp,Ainv,m,dAmat,3*n,1.0_wp,dqdr,3*n)
   call dgemm('n','n',3*n,m,m,-1.0_wp,dAmat,3*n,Ainv,m,1.0_wp,dqdr,3*n)
   !print'(/,"analytical gradient")'
   !print'(3f20.14)',dqdr(:,:,:n)

!! ------------------------------------------------------------------------
!  molecular Hessian calculation
!! ------------------------------------------------------------------------
   do i = 1, n
      do j = 1, i-1
         rij = xyz(:,j) - xyz(:,i)
         r2 = sum(rij**2)
         r = sqrt(r2)
         gamij = 1.0_wp/sqrt(alpha(i)+alpha(j))
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         erfterm = Xtmp(i)*Xtmp(j)*erf(arg)/r
         expterm = Xtmp(i)*Xtmp(j)*2*gamij*exp(-arg2)/sqrtpi
         ! ∂²(qAq)/(∂Ri∂Rj):
         ! ∂²(qAq)/(∂Xi∂Xi) = (1-3X²ij/R²ij-2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         !                  - (R²ij-3X²ij) erf[γij·Rij]/R⁵ij
         ! ∂²(qAq)/(∂Xi∂Xj) = (R²ij-3X²ij) erf[γij·Rij]/R⁵ij
         !                  - (1-3X²ij/R²ij-2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         ! ∂²(qAq)/(∂Xi∂Yi) = 3X²ij erf[γij·Rij]/R⁵ij
         !                  - (3X²ij/R²ij+2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         ! ∂²(qAq)/(∂Xi∂Yj) = (3X²ij/R²ij+2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         !                  - 3X²ij erf[γij·Rij]/R⁵ij
         rxr(1,1) = erfterm * ( 3*rij(1)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(1)**2/r2**2 + 2*gamij2*rij(1)**2/r2 - 1/r2 )
         rxr(2,2) = erfterm * ( 3*rij(2)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(2)**2/r2**2 + 2*gamij2*rij(2)**2/r2 - 1/r2 )
         rxr(3,3) = erfterm * ( 3*rij(3)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(3)**2/r2**2 + 2*gamij2*rij(3)**2/r2 - 1/r2 )
         rxr(2,1) = erfterm * 3*rij(2)*rij(1)/r2**2 &
                  - expterm * ( 3*rij(2)*rij(1)/r2**2 + 2*gamij2*rij(2)*rij(1)/r2 )
         rxr(3,1) = erfterm * 3*rij(3)*rij(1)/r2**2 &
                  - expterm * ( 3*rij(3)*rij(1)/r2**2 + 2*gamij2*rij(3)*rij(1)/r2 )
         rxr(3,2) = erfterm * 3*rij(3)*rij(2)/r2**2 &
                  - expterm * ( 3*rij(3)*rij(2)/r2**2 + 2*gamij2*rij(3)*rij(2)/r2 )

         do k = 1, m
            rxr(1,1) = rxr(1,1) + 0.5_wp*dqdr(1,i,k)*dAmat(1,j,k) &
                                + 0.5_wp*dqdr(1,j,k)*dAmat(1,i,k)
            rxr(2,1) = rxr(2,1) + 0.5_wp*dqdr(2,i,k)*dAmat(1,j,k) &
                                + 0.5_wp*dqdr(2,j,k)*dAmat(1,i,k)
            rxr(3,1) = rxr(3,1) + 0.5_wp*dqdr(3,i,k)*dAmat(1,j,k) &
                                + 0.5_wp*dqdr(3,j,k)*dAmat(1,i,k)
            rxr(2,2) = rxr(2,2) + 0.5_wp*dqdr(2,i,k)*dAmat(2,j,k) &
                                + 0.5_wp*dqdr(2,j,k)*dAmat(2,i,k)
            rxr(3,2) = rxr(3,2) + 0.5_wp*dqdr(3,i,k)*dAmat(2,j,k) &
                                + 0.5_wp*dqdr(3,j,k)*dAmat(2,i,k)
            rxr(3,3) = rxr(3,3) + 0.5_wp*dqdr(3,i,k)*dAmat(3,j,k) &
                                + 0.5_wp*dqdr(3,j,k)*dAmat(3,i,k)
         enddo
         ! symmetrize
         rxr(1,2) = rxr(2,1)
         rxr(1,3) = rxr(3,1)
         rxr(2,3) = rxr(3,2)

         ! save diagonal elements for atom i
         hess(ind(1,i,1,i))=hess(ind(1,i,1,i)) + kq*rxr(1,1)
         hess(ind(2,i,1,i))=hess(ind(2,i,1,i)) + kq*rxr(2,1)
         hess(ind(2,i,2,i))=hess(ind(2,i,2,i)) + kq*rxr(2,2)
         hess(ind(3,i,1,i))=hess(ind(3,i,1,i)) + kq*rxr(3,1)
         hess(ind(3,i,2,i))=hess(ind(3,i,2,i)) + kq*rxr(3,2)
         hess(ind(3,i,3,i))=hess(ind(3,i,3,i)) + kq*rxr(3,3)
         ! save elements between atom i and atom j
         hess(ind(1,i,1,j))=hess(ind(1,i,1,j)) - kq*rxr(1,1)
         hess(ind(1,i,2,j))=hess(ind(1,i,2,j)) - kq*rxr(2,1)
         hess(ind(1,i,3,j))=hess(ind(1,i,3,j)) - kq*rxr(3,1)
         hess(ind(2,i,1,j))=hess(ind(2,i,1,j)) - kq*rxr(2,1)
         hess(ind(2,i,2,j))=hess(ind(2,i,2,j)) - kq*rxr(2,2)
         hess(ind(2,i,3,j))=hess(ind(2,i,3,j)) - kq*rxr(3,2)
         hess(ind(3,i,1,j))=hess(ind(3,i,1,j)) - kq*rxr(3,1)
         hess(ind(3,i,2,j))=hess(ind(3,i,2,j)) - kq*rxr(3,2)
         hess(ind(3,i,3,j))=hess(ind(3,i,3,j)) - kq*rxr(3,3)
         ! save diagonal elements for atom j
         hess(ind(1,j,1,j))=hess(ind(1,j,1,j)) + kq*rxr(1,1)
         hess(ind(2,j,1,j))=hess(ind(2,j,1,j)) + kq*rxr(2,1)
         hess(ind(2,j,2,j))=hess(ind(2,j,2,j)) + kq*rxr(2,2)
         hess(ind(3,j,1,j))=hess(ind(3,j,1,j)) + kq*rxr(3,1)
         hess(ind(3,j,2,j))=hess(ind(3,j,2,j)) + kq*rxr(3,2)
         hess(ind(3,j,3,j))=hess(ind(3,j,3,j)) + kq*rxr(3,3)
      enddo
   enddo

   ! ∂²(qA)/(∂Ri∂q)·∂q/∂Rj
  ! hessian = hessian + reshape(matmul(reshape(dqdr,(/3*n,m/)),&
  !    transpose(reshape(dAmat,(/3*n,m/)))),(/3,n,3,n/))
   !call dgemm('n','t',3*n,m,3*n,+1.0_wp,dqdr,3*n,dAmat,3*n,1.0_wp,hessian,3*n)
   !call dgemm('n','t',3*n,m,3*n,+1.0_wp,dAmat,3*n,dqdr,3*n,1.0_wp,hessian,3*n)

end subroutine mh_eeq

end module xtb_modelhessian

      subroutine ddvopt(Cart,nAtoms,Hess,iANr,s6)
      Implicit Integer(i-n)
      Implicit Real*8 (a-h, o-z)
! include "common/real.inc" (molpro 2002.6)
      Real*8 Zero, One, Two, Three, Four, Five, Six, Seven, &
     &       Eight, RNine, Ten, Half, Pi, SqrtP2, TwoP34, &
     &       TwoP54, One2C2
      Parameter(Zero =0.0D0, One  =1.0D0, Two=2.0D0, Three=3.0D0, &
     &          Four =4.0D0, Five =5.0D0, Six=6.0D0, Seven=7.0D0, &
     &          Eight=8.0D0, rNine=9.0D0, Ten=1.0D1, Half=0.5D0, &
     &          Pi    =3.141592653589793D0, &
     &          SqrtP2=0.8862269254527579D0, &
     &          TwoP34=0.2519794355383808D0, &
     &          TwoP54=5.914967172795612D0, &
     &          One2C2=0.2662567690426443D-04)
! end: common/real.inc

      Real*8 Cart(3,nAtoms),rij(3),rjk(3),rkl(3), &
     &       Hess((3*nAtoms)*(3*nAtoms+1)/2),si(3),sj(3),sk(3), &
     &       sl(3),sm(3),x(2),y(2),z(2), &
     &       xyz(3,4), C(3,4), Dum(3,4,3,4)
      Integer   iANr(nAtoms)

      logical rcutoff

! include  "common/ddvdt.inc" (molpro 2002.6)
      Real*8 rAV(3,3), aAV(3,3), &
     &       B_Str(6), A_Bend(2), A_Trsn(2), A_StrH(2), &
     &       rkr, rkf, A_Str, RF_Const, &
     &       wthr

      Data rAv/1.3500d+00,2.1000d+00,2.5300d+00, &
     &         2.1000d+00,2.8700d+00,3.4000d+00, &
     &         2.5300d+00,3.4000d+00,3.4000d+00/
      Data aAv/1.0000d+00,0.3949d+00,0.3949d+00, &
     &         0.3949d+00,0.2800d+00,0.2800d+00, &
     &         0.3949d+00,0.2800d+00,0.2800d+00/
!org  Data rkr,rkf,rkt/0.4500D+00,0.1500D+00,0.5000D-02/
      Data rkr,rkf,rkt/0.4000D+00,0.1300D+00,0.7500D-02/
      Data A_Str/1.734d0/
      Data B_Str/-.244d0,0.352d0,1.085d0,0.660d0,1.522d0,2.068d0/
      Data A_Bend/0.160d0,0.250d0/
      Data A_Trsn/0.0023d0,0.07d0/
      Data A_StrH/0.3601d0,1.944d0/
      Data RF_Const/1.0D-2/
      Data wthr/0.2/
! end: "common/ddvdt.inc"

!cc VDWx-Parameters (Grimme) used for vdw-correction of model hessian
      real*8 alphavdw, damp, c6(100), c6k, c6l, c66, vander(100), &
     &   vdw(3,3), dr(3)
      integer kxyz, lxyz
      data vander &
! H, He
     &     /0.91d0,0.92d0, &
! Li-Ne
     &      0.75d0,1.28d0,1.35d0,1.32d0,1.27d0,1.22d0,1.17d0,1.13d0, &
! Na-Ar
     &      1.04d0,1.24d0,1.49d0,1.56d0,1.55d0,1.53d0,1.49d0,1.45d0, &
! K, Ca
     &      1.35d0,1.34d0, &
! Sc-Zn
     &      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0, &
     &      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0, &
! Ga-Kr
     &      1.50d0,1.57d0,1.60d0,1.61d0,1.59d0,1.57d0, &
! Rb, Sr
     &      1.48d0,1.46d0, &
! Y-Cd
     &      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0, &
     &      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0, &
! In, Sn, Sb, Te, I, Xe
     &      1.52d0,1.64d0,1.71d0,1.72d0,1.72d0,1.71d0, &
     &      46*2.d0/
!cc End: VDWx ccccccccccccccccc

!
!------- Statement functions
!
      ixyz(i,iAtom) = (iAtom-1)*3 + i
      Jnd(i,j) = i*(i-1)/2 +j
      Ind(i,iAtom,j,jAtom)=Jnd(Max(ixyz(i,iAtom),ixyz(j,jAtom)), &
     &                         Min(ixyz(i,iAtom),ixyz(j,jAtom)))
!end

!cc VDWx cccccccccccccccccccc
       c6 =   50.0
       c6( 1)=0.14
       c6( 2)=0.08
       c6( 3)=1.61
       c6( 4)=1.61
       c6( 5)=3.13
       c6( 6)=1.75
       c6( 7)=1.23
       c6( 8)=0.70
       c6( 9)=0.75
       c6(10)=0.63
       c6(11)=5.71
       c6(12)=5.71
       c6(13)=10.79
       c6(14)=9.23
       c6(15)=7.84
       c6(16)=5.57
       c6(17)=5.07
       c6(18)=4.61
       c6(19:30)=10.8
       c6(31)=16.99
       c6(32)=17.10
       c6(33)=16.37
       c6(34)=12.64
       c6(35)=12.47
       c6(36)=12.01
       c6(37:48)=24.67
       c6(49)=37.32
       c6(50)=38.71
       c6(51)=38.44
       c6(52)=31.74
       c6(53)=31.50
       c6(54)=29.99
!cc ENDE VDWx cccccccccccccccc

      bohr=0.52917726d0
      Fact=One
!hjw threshold reduced
      rZero=1.0d-10
      n3=3*nAtoms
      Hess = 0.0d0

!
!     Hessian for tension
!
 777  Do kAtom = 1, nAtoms
         kr=iTabRow(iANr(kAtom))
!        If (kr.eq.0) Go To 5

         Do lAtom = 1, kAtom-1
            lr=iTabRow(iANr(lAtom))
!           If (lr.eq.0) Go To 10
            xkl=Cart(1,kAtom)-Cart(1,lAtom)
            ykl=Cart(2,kAtom)-Cart(2,lAtom)
            zkl=Cart(3,kAtom)-Cart(3,lAtom)
            rkl2 = xkl**2 + ykl**2 + zkl**2
            r0=rAv(kr,lr)
            alpha=aAv(kr,lr)

!cccccc VDWx ccccccccccccccccccccccccccccccccc
            c6k=c6(iANr(katom))
            c6l=c6(iANr(latom))
            c66=sqrt(c6k*c6l)
            Rv=(vander(iANr(katom))+vander(iANr(latom)))/0.52917721

            call getvdwxx(xkl, ykl, zkl, c66, s6, Rv, vdw(1,1))
            call getvdwxy(xkl, ykl, zkl, c66, s6, Rv, vdw(1,2))
            call getvdwxy(xkl, zkl, ykl, c66, s6, Rv, vdw(1,3))
            call getvdwxx(ykl, xkl, zkl, c66, s6, Rv, vdw(2,2))
            call getvdwxy(ykl, zkl, xkl, c66, s6, Rv, vdw(2,3))
            call getvdwxx(zkl, xkl, ykl, c66, s6, Rv, vdw(3,3))
!cccccc Ende VDWx ccccccccccccccccccccccccccccccc

            gamma=rkr*Exp(alpha*r0**2)
! not better: *sqrt(abs(wb(kAtom,lAtom)))
            gmm=gamma*Exp(-alpha*rkl2)
            Hxx=gmm*xkl*xkl/rkl2-vdw(1,1)
            Hxy=gmm*xkl*ykl/rkl2-vdw(1,2)
            Hxz=gmm*xkl*zkl/rkl2-vdw(1,3)
            Hyy=gmm*ykl*ykl/rkl2-vdw(2,2)
            Hyz=gmm*ykl*zkl/rkl2-vdw(2,3)
            Hzz=gmm*zkl*zkl/rkl2-vdw(3,3)

!
            Hess(Ind(1,kAtom,1,kAtom))=Hess(Ind(1,kAtom,1,kAtom))+Hxx
            Hess(Ind(2,kAtom,1,kAtom))=Hess(Ind(2,kAtom,1,kAtom))+Hxy
            Hess(Ind(2,kAtom,2,kAtom))=Hess(Ind(2,kAtom,2,kAtom))+Hyy
            Hess(Ind(3,kAtom,1,kAtom))=Hess(Ind(3,kAtom,1,kAtom))+Hxz
            Hess(Ind(3,kAtom,2,kAtom))=Hess(Ind(3,kAtom,2,kAtom))+Hyz
            Hess(Ind(3,kAtom,3,kAtom))=Hess(Ind(3,kAtom,3,kAtom))+Hzz
!
            Hess(Ind(1,kAtom,1,lAtom))=Hess(Ind(1,kAtom,1,lAtom))-Hxx
            Hess(Ind(1,kAtom,2,lAtom))=Hess(Ind(1,kAtom,2,lAtom))-Hxy
            Hess(Ind(1,kAtom,3,lAtom))=Hess(Ind(1,kAtom,3,lAtom))-Hxz
            Hess(Ind(2,kAtom,1,lAtom))=Hess(Ind(2,kAtom,1,lAtom))-Hxy
            Hess(Ind(2,kAtom,2,lAtom))=Hess(Ind(2,kAtom,2,lAtom))-Hyy
            Hess(Ind(2,kAtom,3,lAtom))=Hess(Ind(2,kAtom,3,lAtom))-Hyz
            Hess(Ind(3,kAtom,1,lAtom))=Hess(Ind(3,kAtom,1,lAtom))-Hxz
            Hess(Ind(3,kAtom,2,lAtom))=Hess(Ind(3,kAtom,2,lAtom))-Hyz
            Hess(Ind(3,kAtom,3,lAtom))=Hess(Ind(3,kAtom,3,lAtom))-Hzz
!
            Hess(Ind(1,lAtom,1,lAtom))=Hess(Ind(1,lAtom,1,lAtom))+Hxx
            Hess(Ind(2,lAtom,1,lAtom))=Hess(Ind(2,lAtom,1,lAtom))+Hxy
            Hess(Ind(2,lAtom,2,lAtom))=Hess(Ind(2,lAtom,2,lAtom))+Hyy
            Hess(Ind(3,lAtom,1,lAtom))=Hess(Ind(3,lAtom,1,lAtom))+Hxz
            Hess(Ind(3,lAtom,2,lAtom))=Hess(Ind(3,lAtom,2,lAtom))+Hyz
            Hess(Ind(3,lAtom,3,lAtom))=Hess(Ind(3,lAtom,3,lAtom))+Hzz
!

 10         Continue
         End Do

 5       Continue
      End Do


!
!     Hessian for bending
!
      Do mAtom = 1, nAtoms
         mr=iTabRow(iANr(mAtom))
!        If (mr.eq.0) Go To 20
         Do iAtom = 1, nAtoms
           If (iAtom.eq.mAtom) Go To 30
           ir=iTabRow(iANr(iAtom))
!          If (ir.eq.0) Go To 30
           if(rcutoff(cart,iatom,matom)) cycle
!          if(wb(iatom,matom).lt.wthr) cycle
           Do jAtom = 1, iAtom-1
            If (jAtom.eq.mAtom) Go To 40
            jr=iTabRow(iANr(jAtom))
!           If (jr.eq.0) Go To 40
            if(rcutoff(cart,jatom,iatom)) cycle
            if(rcutoff(cart,jatom,matom)) cycle
!           if(wb(jatom,iatom).lt.wthr) cycle
!           if(wb(jatom,matom).lt.wthr) cycle

            xmi=(Cart(1,iAtom)-Cart(1,mAtom))
            ymi=(Cart(2,iAtom)-Cart(2,mAtom))
            zmi=(Cart(3,iAtom)-Cart(3,mAtom))
            rmi2 = xmi**2 + ymi**2 + zmi**2
            rmi=sqrt(rmi2)
            r0mi=rAv(mr,ir)
            ami=aAv(mr,ir)
!
            xmj=(Cart(1,jAtom)-Cart(1,mAtom))
            ymj=(Cart(2,jAtom)-Cart(2,mAtom))
            zmj=(Cart(3,jAtom)-Cart(3,mAtom))
            rmj2 = xmj**2 + ymj**2 + zmj**2
            rmj=sqrt(rmj2)
            r0mj=rAv(mr,jr)
            amj=aAv(mr,jr)
!
!---------- Test if zero angle
!
            Test=xmi*xmj+ymi*ymj+zmi*zmj
            Test=Test/(rmi*rmj)
            If (Test.eq.One) Go To 40
!
            xij=(Cart(1,jAtom)-Cart(1,iAtom))
            yij=(Cart(2,jAtom)-Cart(2,iAtom))
            zij=(Cart(3,jAtom)-Cart(3,iAtom))
            rij2 = xij**2 + yij**2 + zij**2
            rrij=sqrt(rij2)
!
            alpha=rkf*exp((ami*r0mi**2+amj*r0mj**2))
!
            r=sqrt(rmj2+rmi2)
            gij=alpha*exp(-(ami*rmi2+amj*rmj2))
!           Write (*,*) ' gij=',gij
            rL2=(ymi*zmj-zmi*ymj)**2+(zmi*xmj-xmi*zmj)**2+ &
     &         (xmi*ymj-ymi*xmj)**2
!hjw modified
            if(rL2.lt.1.d-14) then
              rL=0
            else
              rL=sqrt(rL2)
            end if
!
            if ((rmj.gt.rZero).and.(rmi.gt.rZero).and. &
     &                                (rrij.gt.rZero)) Then
              SinPhi=rL/(rmj*rmi)
              rmidotrmj=xmi*xmj+ymi*ymj+zmi*zmj
              CosPhi=rmidotrmj/(rmj*rmi)
!
!-------------None linear case
!
              If (SinPhi.gt.rZero) Then
!               Write (*,*) ' None linear case'
                si(1)=(xmi/rmi*cosphi-xmj/rmj)/(rmi*sinphi)
                si(2)=(ymi/rmi*cosphi-ymj/rmj)/(rmi*sinphi)
                si(3)=(zmi/rmi*cosphi-zmj/rmj)/(rmi*sinphi)
                sj(1)=(cosphi*xmj/rmj-xmi/rmi)/(rmj*sinphi)
                sj(2)=(cosphi*ymj/rmj-ymi/rmi)/(rmj*sinphi)
                sj(3)=(cosphi*zmj/rmj-zmi/rmi)/(rmj*sinphi)
                sm(1)=-si(1)-sj(1)
                sm(2)=-si(2)-sj(2)
                sm(3)=-si(3)-sj(3)
                Do icoor=1,3
                   Do jCoor=1,3
                    If (mAtom.gt.iAtom) Then
                       Hess(Ind(icoor,mAtom,jcoor,iAtom))= &
     &                        Hess(Ind(icoor,mAtom,jcoor,iAtom)) &
     &                        +gij*sm(icoor)*si(jcoor)
                    else
                      Hess(Ind(icoor,iAtom,jcoor,mAtom))= &
     &                        Hess(Ind(icoor,iAtom,jcoor,mAtom)) &
     &                        +gij*si(icoor)*sm(jcoor)
                    End If
                    If (mAtom.gt.jAtom) Then
                        Hess(Ind(icoor,mAtom,jcoor,jAtom))= &
     &                        Hess(Ind(icoor,mAtom,jcoor,jAtom)) &
     &                        +gij*sm(icoor)*sj(jcoor)
                    else
                      Hess(Ind(icoor,jAtom,jcoor,mAtom))= &
     &                        Hess(Ind(icoor,jAtom,jcoor,mAtom)) &
     &                        +gij*sj(icoor)*sm(jcoor)
                    End If
                    If (iAtom.gt.jAtom) Then
                        Hess(Ind(icoor,iAtom,jcoor,jAtom))= &
     &                        Hess(Ind(icoor,iAtom,jcoor,jAtom)) &
     &                        +gij*si(icoor)*sj(jcoor)
                     else
                        Hess(Ind(icoor,jAtom,jcoor,iAtom))= &
     &                        Hess(Ind(icoor,jAtom,jcoor,iAtom)) &
     &                        +gij*sj(icoor)*si(jcoor)
                     End If
                   End Do
                End Do
                Do icoor=1,3
                  Do jCoor=1,icoor
                    Hess(Ind(icoor,iAtom,jcoor,iAtom))= &
     &                        Hess(Ind(icoor,iAtom,jcoor,iAtom)) &
     &                        +gij*si(icoor)*si(jcoor)
                    Hess(Ind(icoor,mAtom,jcoor,mAtom))= &
     &                        Hess(Ind(icoor,mAtom,jcoor,mAtom)) &
     &                        +gij*sm(icoor)*sm(jcoor)
                    Hess(Ind(icoor,jAtom,jcoor,jAtom))= &
     &                        Hess(Ind(icoor,jAtom,jcoor,jAtom)) &
     &                        +gij*sj(icoor)*sj(jcoor)

!
                  End Do
                End Do
              Else
!
!----------------Linear case
!
!               Write (*,*) 'linear case'
                    if ((abs(ymi).gt.rZero).or. &
     &                 (abs(xmi).gt.rZero)) Then
                      x(1)=-ymi
                      y(1)=xmi
                      z(1)=Zero
                      x(2)=-xmi*zmi
                      y(2)=-ymi*zmi
                      z(2)=xmi*xmi+ymi*ymi
                    Else
                      x(1)=One
                      y(1)=Zero
                      z(1)=Zero
                      x(2)=Zero
                      y(2)=One
                      z(2)=Zero
                    End If
                    Do i=1,2
                     r1=sqrt(x(i)**2+y(i)**2+z(i)**2)
                     cosThetax=x(i)/r1
                     cosThetay=y(i)/r1
                     cosThetaz=z(i)/r1
                     si(1)=-cosThetax/rmi
                     si(2)=-cosThetay/rmi
                     si(3)=-cosThetaz/rmi
                     sj(1)=-cosThetax/rmj
                     sj(2)=-cosThetay/rmj
                     sj(3)=-cosThetaz/rmj
                     sm(1)=-(si(1)+sj(1))
                     sm(2)=-(si(2)+sj(2))
                     sm(3)=-(si(3)+sj(3))
!
                     Do icoor=1,3
                       Do jCoor=1,3
                        If (mAtom.gt.iAtom) Then
                          Hess(Ind(icoor,mAtom,jcoor,iAtom))= &
     &                        Hess(Ind(icoor,mAtom,jcoor,iAtom)) &
     &                         +gij*sm(icoor)*si(jcoor)
                        else
                           Hess(Ind(icoor,iAtom,jcoor,mAtom))= &
     &                        Hess(Ind(icoor,iAtom,jcoor,mAtom)) &
     &                         +gij*si(icoor)*sm(jcoor)
                        End If
                        If (mAtom.gt.jAtom) Then
                          Hess(Ind(icoor,mAtom,jcoor,jAtom))= &
     &                        Hess(Ind(icoor,mAtom,jcoor,jAtom)) &
     &                         +gij*sm(icoor)*sj(jcoor)
                        else
                          Hess(Ind(icoor,jAtom,jcoor,mAtom))= &
     &                        Hess(Ind(icoor,jAtom,jcoor,mAtom)) &
     &                         +gij*sj(icoor)*sm(jcoor)
                        End If
                        If (iAtom.gt.jAtom) Then
                           Hess(Ind(icoor,iAtom,jcoor,jAtom))= &
     &                        Hess(Ind(icoor,iAtom,jcoor,jAtom)) &
     &                         +gij*si(icoor)*sj(jcoor)
                        else
                           Hess(Ind(icoor,jAtom,jcoor,iAtom))= &
     &                        Hess(Ind(icoor,jAtom,jcoor,iAtom)) &
     &                         +gij*sj(icoor)*si(jcoor)
                        End If
                       End Do
                     End Do
                     Do icoor=1,3
                       Do jCoor=1,icoor
                         Hess(Ind(icoor,iAtom,jcoor,iAtom))= &
     &                        Hess(Ind(icoor,iAtom,jcoor,iAtom)) &
     &                         +gij*si(icoor)*si(jcoor)
                         Hess(Ind(icoor,mAtom,jcoor,mAtom))= &
     &                        Hess(Ind(icoor,mAtom,jcoor,mAtom)) &
     &                         +gij*sm(icoor)*sm(jcoor)
                         Hess(Ind(icoor,jAtom,jcoor,jAtom))= &
     &                         Hess(Ind(icoor,jAtom,jcoor,jAtom)) &
     &                         +gij*sj(icoor)*sj(jcoor)
                       End Do
                     End Do
                 End Do
              End If
            End If
!
 40         Continue
        End Do
 30     Continue
       End Do
 20    Continue
      End Do
!
!     Hessian for torsion
!
      Do jAtom = 1,nAtoms
        jr=iTabRow(iANr(jAtom))
!       If (jr.eq.0) Go To 444
!
        Call DCopy(3,Cart(1,jAtom),1,xyz(1,2),1)
!
        Do kAtom = 1, nAtoms
           If (kAtom.eq.jAtom) Go To 111
           kr=iTabRow(iANr(kAtom))
!          If (kr.eq.0) Go To 111

           if(rcutoff(cart,katom,jatom)) cycle
!          if(wb(katom,jatom).lt.wthr) cycle
!
           Call DCopy(3,Cart(1,kAtom),1,xyz(1,3),1)
!
           Do iAtom = 1, nAtoms
              ij_=nAtoms*(jAtom-1)+iAtom
              If (iAtom.eq.jAtom) Go To 333
              If (iAtom.eq.kAtom) Go To 333
              ir=iTabRow(iANr(iAtom))
!             If (ir.eq.0) Go To 333
!
              if(rcutoff(cart,iatom,katom)) cycle
              if(rcutoff(cart,iatom,jatom)) cycle
!             if(wb(iatom,katom).lt.wthr) cycle
!             if(wb(iatom,jatom).lt.wthr) cycle

              Call DCopy(3,Cart(1,iAtom),1,xyz(1,1),1)
!
              Do lAtom = 1, nAtoms
                 lk_=nAtoms*(kAtom-1)+lAtom
                 If (ij_.le.lk_) Go To 222
                 If (lAtom.eq.iAtom) Go To 222
                 If (lAtom.eq.jAtom) Go To 222
                 If (lAtom.eq.kAtom) Go To 222
                 lr=iTabRow(iANr(lAtom))
!                If (lr.eq.0) Go To 222
!
                 if(rcutoff(cart,latom,iatom)) cycle
                 if(rcutoff(cart,latom,katom)) cycle
                 if(rcutoff(cart,latom,jatom)) cycle
!                if(wb(latom,iatom).lt.wthr) cycle
!                if(wb(latom,katom).lt.wthr) cycle
!                if(wb(latom,jatom).lt.wthr) cycle

                 Call DCopy(3,Cart(1,lAtom),1,xyz(1,4),1)
!
                 rij(1)=Cart(1,iAtom)-Cart(1,jAtom)
                 rij(2)=Cart(2,iAtom)-Cart(2,jAtom)
                 rij(3)=Cart(3,iAtom)-Cart(3,jAtom)
                 rij0=rAv(ir,jr)**2
                 aij =aAv(ir,jr)
!
                 rjk(1)=Cart(1,jAtom)-Cart(1,kAtom)
                 rjk(2)=Cart(2,jAtom)-Cart(2,kAtom)
                 rjk(3)=Cart(3,jAtom)-Cart(3,kAtom)
                 rjk0=rAv(jr,kr)**2
                 ajk =aAv(jr,kr)
!
                 rkl(1)=Cart(1,kAtom)-Cart(1,lAtom)
                 rkl(2)=Cart(2,kAtom)-Cart(2,lAtom)
                 rkl(3)=Cart(3,kAtom)-Cart(3,lAtom)
                 rkl0=rAv(kr,lr)**2
                 akl =aAv(kr,lr)
!
                 rij2=rij(1)**2+rij(2)**2+rij(3)**2
                 rjk2=rjk(1)**2+rjk(2)**2+rjk(3)**2
                 rkl2=rkl(1)**2+rkl(2)**2+rkl(3)**2

!              Allow only angles in the range of 35-145
               A35 = (35.0D0/180.D0)* Pi
               CosFi_Max=Cos(A35)
               CosFi2=(rij(1)*rjk(1)+rij(2)*rjk(2)+rij(3)*rjk(3)) &
     &               /Sqrt(rij2*rjk2)
               If (Abs(CosFi2).gt.CosFi_Max) Go To 222
               CosFi3=(rkl(1)*rjk(1)+rkl(2)*rjk(2)+rkl(3)*rjk(3)) &
     &               /Sqrt(rkl2*rjk2)
               If (Abs(CosFi3).gt.CosFi_Max) Go To 222

               beta=rkt* &
     &                       exp( (aij*rij0+ajk*rjk0+akl*rkl0))
               tij=beta*exp(-(aij*rij2+ajk*rjk2+akl*rkl2))

               Call Trsn(xyz,4,Tau,C,.False.,.False.,'        ', &
     &                  Dum,.False.)
               Call DCopy(3,C(1,1),1,si,1)
               Call DCopy(3,C(1,2),1,sj,1)
               Call DCopy(3,C(1,3),1,sk,1)
               Call DCopy(3,C(1,4),1,sl,1)
!
!-------------Off diagonal block
!
              Do icoor=1,3
                Do jCoor=1,3
                 Hess(Ind(icoor,iAtom,jcoor,jAtom))= &
     &           Hess(Ind(icoor,iAtom,jcoor,jAtom)) &
     &            +tij*si(icoor) * sj(jcoor)
                 Hess(Ind(icoor,iAtom,jcoor,kAtom))= &
     &           Hess(Ind(icoor,iAtom,jcoor,kAtom)) &
     &            +tij*si(icoor) * sk(jcoor)
                 Hess(Ind(icoor,iAtom,jcoor,lAtom))= &
     &           Hess(Ind(icoor,iAtom,jcoor,lAtom)) &
     &            +tij*si(icoor) * sl(jcoor)
                 Hess(Ind(icoor,jAtom,jcoor,kAtom))= &
     &           Hess(Ind(icoor,jAtom,jcoor,kAtom)) &
     &            +tij*sj(icoor) * sk(jcoor)
                 Hess(Ind(icoor,jAtom,jcoor,lAtom))= &
     &           Hess(Ind(icoor,jAtom,jcoor,lAtom)) &
     &            +tij*sj(icoor) * sl(jcoor)
                 Hess(Ind(icoor,kAtom,jcoor,lAtom))= &
     &           Hess(Ind(icoor,kAtom,jcoor,lAtom)) &
     &            +tij*sk(icoor) * sl(jcoor)

                End Do
              End Do
!
!-------------Diagonal block
!
              Do icoor=1,3
                Do jCoor=1,icoor
                 Hess(Ind(icoor,iAtom,jcoor,iAtom))= &
     &           Hess(Ind(icoor,iAtom,jcoor,iAtom)) &
     &            +tij*si(icoor) * si(jcoor)
                 Hess(Ind(icoor,jAtom,jcoor,jAtom))= &
     &           Hess(Ind(icoor,jAtom,jcoor,jAtom)) &
     &            +tij*sj(icoor) * sj(jcoor)
                 Hess(Ind(icoor,kAtom,jcoor,kAtom))= &
     &           Hess(Ind(icoor,kAtom,jcoor,kAtom)) &
     &            +tij*sk(icoor) * sk(jcoor)
                 Hess(Ind(icoor,lAtom,jcoor,lAtom))= &
     &           Hess(Ind(icoor,lAtom,jcoor,lAtom)) &
     &            +tij*sl(icoor) * sl(jcoor)

!
                 End Do
               End Do
222            Continue
             End Do        ! lAtom
333          Continue
           End Do          ! iAtom
111        Continue
        End Do             ! kAtom
444     Continue
      End Do               ! jAtom
      Return
      End

      subroutine gff_ddvopt(Cart,nAtoms,Hess,at,s6,param,topo,neigh)
      use xtb_gfnff_data, only : TGFFData
      use xtb_gfnff_topology, only : TGFFTopology
      use xtb_gfnff_neighbor, only : TNeigh
      use xtb_type_timer
      Implicit Integer (i-n)
      Implicit Real*8 (a-h, o-z)
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(in) :: topo
      type(TNeigh), intent(in) :: neigh

      Integer at(nAtoms)
      Real*8 Cart(3,nAtoms)
      Real*8 Hess((3*nAtoms)*(3*nAtoms+1)/2)
      Real*8 s6
!     charges qa taken from gffcom

      Real*8 rij(3),rjk(3),rkl(3),&
     &       si(3),sj(3),sk(3),&
     &       sl(3),sm(3),x(2),y(2),z(2),&
     &       xyz(3,4), C(3,4), Dum(3,4,3,4) 
      Real*8 rAV(3,3), aAV(3,3),rkr, rkf
      Integer iANr(nAtoms)
      logical profile

      Data rAv/1.3500d+00,2.1000d+00,2.5300d+00,&
     &         2.1000d+00,2.8700d+00,3.4000d+00,&
     &         2.5300d+00,3.4000d+00,3.4000d+00/ 
      Data aAv/1.0000d+00,0.3949d+00,0.3949d+00,&
     &         0.3949d+00,0.2800d+00,0.2800d+00,&
     &         0.3949d+00,0.2800d+00,0.2800d+00/ 
!     Data rkr,rkf,rkt/0.4500D+00,0.1500D+00,0.5000D-02/ ! org
      Data rkr,rkf,rkt/0.4500D+00,0.3000D+00,0.75000  /  ! adjusted to account for redundant angles in old version

!cc VDWx-Parameters (Grimme) used for vdw-correction of model hessian 
      real*8 c6(55),c66
! NOT USED ANYMORE (BJ instead)
!     data vander &
! H, He
!    &     /0.91d0,0.92d0,&
! Li-Ne
!    &      0.75d0,1.28d0,1.35d0,1.32d0,1.27d0,1.22d0,1.17d0,1.13d0,&
! Na-Ar
!    &      1.04d0,1.24d0,1.49d0,1.56d0,1.55d0,1.53d0,1.49d0,1.45d0,&
! K, Ca
!    &      1.35d0,1.34d0,&
! Sc-Zn
!    &      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,&
!    &      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,&
! Ga-Kr
!    &      1.50d0,1.57d0,1.60d0,1.61d0,1.59d0,1.57d0,&
! Rb, Sr
!    &      1.48d0,1.46d0,&
! Y-Cd
!    &      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,&
!    &      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,&
! In, Sn, Sb, Te, I, Xe
!    &      1.52d0,1.64d0,1.71d0,1.72d0,1.72d0,1.71d0/ 

      Real*8 Zero, One, Two, Three, Four, Five, Six, Seven,&
     &       Eight, RNine, Ten, Half, Pi, SqrtP2, TwoP34,&
     &       TwoP54, One2C2
      Parameter(Zero =0.0D0, One  =1.0D0, Two=2.0D0, Three=3.0D0,&
     &          Four =4.0D0, Five =5.0D0, Six=6.0D0, Seven=7.0D0,&
     &          Eight=8.0D0, rNine=9.0D0, Ten=1.0D1, Half=0.5D0,&
     &          Pi    =3.141592653589793D0,&
     &          SqrtP2=0.8862269254527579D0,&
     &          TwoP34=0.2519794355383808D0,&
     &          TwoP54=5.914967172795612D0,&
     &          One2C2=0.2662567690426443D-04)

      type(tb_timer) :: timer
!     inline fct
      lina(i,j)=min(i,j)+max(i,j)*(max(i,j)-1)/2        
      ixyz(i,iAtom) = (iAtom-1)*3 + i
      Jnd(i,j) = i*(i-1)/2 +j
      Ind(i,iAtom,j,jAtom)=Jnd(Max(ixyz(i,iAtom),ixyz(j,jAtom)),&
     &                         Min(ixyz(i,iAtom),ixyz(j,jAtom)))

!     map heavy atoms to Z<=54
      do i=1,nAtoms
         iANr(i)=at(i)
         if(at(i).gt.54) iANr(i)=at(i) - 18
         if(at(i).gt.72) iANr(i)=at(i) - 32
         if(at(i).gt.56.and.at(i).lt.72) iANr(i)=39 ! set LNs to Y
         if(at(i).gt.86) iANr(i)=55 ! use same c6 for all atom types > 86
      enddo

      profile=.false.

      if (profile) call timer%new(3,.false.)
      if (profile) call timer%measure(1,'bond/vdw/qq')
!     D2 C6
      c6( 1)=0.14; c6( 2)=0.08; c6( 3)=1.61; c6( 4)=1.61; c6( 5)=3.13
      c6( 6)=1.75; c6( 7)=1.23; c6( 8)=0.70; c6( 9)=0.75; c6(10)=0.63
      c6(11)=5.71; c6(12)=5.71; c6(13)=10.79; c6(14)=9.23; c6(15)=7.84
      c6(16)=5.57; c6(17)=5.07; c6(18)=4.61; c6(19:30)=10.8; c6(31)=16.99
      c6(32)=17.10; c6(33)=16.37; c6(34)=12.64; c6(35)=12.47; c6(36)=12.01
      c6(37:48)=24.67; c6(49)=37.32; c6(50)=38.71; c6(51)=38.44
      c6(52)=31.74; c6(53)=31.50; c6(54)=29.99
      c6(55)=40.0

!hjw threshold reduced
      rZero=1.0d-10
      n3=3*nAtoms
     Hess = 0.0d0


!
!     Hessian for tension
!
      do ibond=1,neigh%nbond
         kAtom=neigh%blist(1,ibond)
         lAtom=neigh%blist(2,ibond)
         kr=iTabRow(iANr(kAtom))
         lr=iTabRow(iANr(lAtom))
            xkl=Cart(1,kAtom)-Cart(1,lAtom)
            ykl=Cart(2,kAtom)-Cart(2,lAtom)
            zkl=Cart(3,kAtom)-Cart(3,lAtom)
            rkl2 = xkl**2 + ykl**2 + zkl**2
            r0=rAv(kr,lr)
            alpha=aAv(kr,lr)
            gmm=Exp(-alpha*rkl2)*rkr*Exp(alpha*r0**2)
            Hxx=gmm*xkl*xkl/rkl2
            Hxy=gmm*xkl*ykl/rkl2
            Hxz=gmm*xkl*zkl/rkl2
            Hyy=gmm*ykl*ykl/rkl2
            Hyz=gmm*ykl*zkl/rkl2
            Hzz=gmm*zkl*zkl/rkl2
            Hess(Ind(1,kAtom,1,kAtom))=Hess(Ind(1,kAtom,1,kAtom))+Hxx
            Hess(Ind(2,kAtom,1,kAtom))=Hess(Ind(2,kAtom,1,kAtom))+Hxy
            Hess(Ind(2,kAtom,2,kAtom))=Hess(Ind(2,kAtom,2,kAtom))+Hyy
            Hess(Ind(3,kAtom,1,kAtom))=Hess(Ind(3,kAtom,1,kAtom))+Hxz
            Hess(Ind(3,kAtom,2,kAtom))=Hess(Ind(3,kAtom,2,kAtom))+Hyz
            Hess(Ind(3,kAtom,3,kAtom))=Hess(Ind(3,kAtom,3,kAtom))+Hzz
            Hess(Ind(1,kAtom,1,lAtom))=Hess(Ind(1,kAtom,1,lAtom))-Hxx
            Hess(Ind(1,kAtom,2,lAtom))=Hess(Ind(1,kAtom,2,lAtom))-Hxy
            Hess(Ind(1,kAtom,3,lAtom))=Hess(Ind(1,kAtom,3,lAtom))-Hxz
            Hess(Ind(2,kAtom,1,lAtom))=Hess(Ind(2,kAtom,1,lAtom))-Hxy
            Hess(Ind(2,kAtom,2,lAtom))=Hess(Ind(2,kAtom,2,lAtom))-Hyy
            Hess(Ind(2,kAtom,3,lAtom))=Hess(Ind(2,kAtom,3,lAtom))-Hyz
            Hess(Ind(3,kAtom,1,lAtom))=Hess(Ind(3,kAtom,1,lAtom))-Hxz
            Hess(Ind(3,kAtom,2,lAtom))=Hess(Ind(3,kAtom,2,lAtom))-Hyz
            Hess(Ind(3,kAtom,3,lAtom))=Hess(Ind(3,kAtom,3,lAtom))-Hzz
            Hess(Ind(1,lAtom,1,lAtom))=Hess(Ind(1,lAtom,1,lAtom))+Hxx
            Hess(Ind(2,lAtom,1,lAtom))=Hess(Ind(2,lAtom,1,lAtom))+Hxy
            Hess(Ind(2,lAtom,2,lAtom))=Hess(Ind(2,lAtom,2,lAtom))+Hyy
            Hess(Ind(3,lAtom,1,lAtom))=Hess(Ind(3,lAtom,1,lAtom))+Hxz
            Hess(Ind(3,lAtom,2,lAtom))=Hess(Ind(3,lAtom,2,lAtom))+Hyz
            Hess(Ind(3,lAtom,3,lAtom))=Hess(Ind(3,lAtom,3,lAtom))+Hzz
      End do
 
!
!     Hessian for Coulomb+dispersion
!

      do kAtom=1,nAtoms
         do lAtom=1,kAtom-1
            xkl=Cart(1,kAtom)-Cart(1,lAtom)
            ykl=Cart(2,kAtom)-Cart(2,lAtom)
            zkl=Cart(3,kAtom)-Cart(3,lAtom)
            rkl2 = xkl**2 + ykl**2 + zkl**2
            if(rkl2.gt.1600.d0) cycle
            c66=-s6*sqrt(c6(iANr(katom))*c6(iANr(latom))) ! D2 value
            rr=sqrt(rkl2)
            rr3=rr*rkl2
            r02=param%d3r0(lina(at(katom),at(latom)))
            rrpa=rr+sqrt(r02) ! qq damping with BJ radius
            cqq=2.0d0*topo%qa(kAtom)*topo%qa(lAtom) ! a bit upscaled
            call getqqxx(xkl,    cqq,c66,rr,rkl2,rr3,rrpa,r02,hxx)
            call getqqxy(xkl,ykl,cqq,c66,rr,rkl2,rr3,rrpa,r02,hxy)
            call getqqxy(xkl,zkl,cqq,c66,rr,rkl2,rr3,rrpa,r02,hxz)
            call getqqxx(ykl,    cqq,c66,rr,rkl2,rr3,rrpa,r02,hyy)
            call getqqxy(ykl,zkl,cqq,c66,rr,rkl2,rr3,rrpa,r02,hyz)
            call getqqxx(zkl,    cqq,c66,rr,rkl2,rr3,rrpa,r02,hzz)
            Hess(Ind(1,kAtom,1,kAtom))=Hess(Ind(1,kAtom,1,kAtom))+Hxx
            Hess(Ind(2,kAtom,1,kAtom))=Hess(Ind(2,kAtom,1,kAtom))+Hxy
            Hess(Ind(2,kAtom,2,kAtom))=Hess(Ind(2,kAtom,2,kAtom))+Hyy
            Hess(Ind(3,kAtom,1,kAtom))=Hess(Ind(3,kAtom,1,kAtom))+Hxz
            Hess(Ind(3,kAtom,2,kAtom))=Hess(Ind(3,kAtom,2,kAtom))+Hyz
            Hess(Ind(3,kAtom,3,kAtom))=Hess(Ind(3,kAtom,3,kAtom))+Hzz
            Hess(Ind(1,kAtom,1,lAtom))=Hess(Ind(1,kAtom,1,lAtom))-Hxx
            Hess(Ind(1,kAtom,2,lAtom))=Hess(Ind(1,kAtom,2,lAtom))-Hxy
            Hess(Ind(1,kAtom,3,lAtom))=Hess(Ind(1,kAtom,3,lAtom))-Hxz
            Hess(Ind(2,kAtom,1,lAtom))=Hess(Ind(2,kAtom,1,lAtom))-Hxy
            Hess(Ind(2,kAtom,2,lAtom))=Hess(Ind(2,kAtom,2,lAtom))-Hyy
            Hess(Ind(2,kAtom,3,lAtom))=Hess(Ind(2,kAtom,3,lAtom))-Hyz
            Hess(Ind(3,kAtom,1,lAtom))=Hess(Ind(3,kAtom,1,lAtom))-Hxz
            Hess(Ind(3,kAtom,2,lAtom))=Hess(Ind(3,kAtom,2,lAtom))-Hyz
            Hess(Ind(3,kAtom,3,lAtom))=Hess(Ind(3,kAtom,3,lAtom))-Hzz
            Hess(Ind(1,lAtom,1,lAtom))=Hess(Ind(1,lAtom,1,lAtom))+Hxx
            Hess(Ind(2,lAtom,1,lAtom))=Hess(Ind(2,lAtom,1,lAtom))+Hxy
            Hess(Ind(2,lAtom,2,lAtom))=Hess(Ind(2,lAtom,2,lAtom))+Hyy
            Hess(Ind(3,lAtom,1,lAtom))=Hess(Ind(3,lAtom,1,lAtom))+Hxz
            Hess(Ind(3,lAtom,2,lAtom))=Hess(Ind(3,lAtom,2,lAtom))+Hyz
            Hess(Ind(3,lAtom,3,lAtom))=Hess(Ind(3,lAtom,3,lAtom))+Hzz
         End Do
      End Do

      if (profile) call timer%measure(1)

!
!     Hessian for bending
!
      if (profile) call timer%measure(2,'bend')
      do iangl=1,topo%nangl
         mAtom=topo%alist(1,iangl)
         iAtom=topo%alist(2,iangl)
         jAtom=topo%alist(3,iangl)
         mr=iTabRow(iANr(mAtom))
         ir=iTabRow(iANr(iAtom))
         jr=iTabRow(iANr(jAtom))
 
            xmi=(Cart(1,iAtom)-Cart(1,mAtom))
            ymi=(Cart(2,iAtom)-Cart(2,mAtom))
            zmi=(Cart(3,iAtom)-Cart(3,mAtom))
            rmi2 = xmi**2 + ymi**2 + zmi**2
            rmi=sqrt(rmi2)
            r0mi=rAv(mr,ir)
            ami=aAv(mr,ir)
            xmj=(Cart(1,jAtom)-Cart(1,mAtom))
            ymj=(Cart(2,jAtom)-Cart(2,mAtom))
            zmj=(Cart(3,jAtom)-Cart(3,mAtom))
            rmj2 = xmj**2 + ymj**2 + zmj**2
            rmj=sqrt(rmj2)
            r0mj=rAv(mr,jr)
            amj=aAv(mr,jr)
            xij=(Cart(1,jAtom)-Cart(1,iAtom))
            yij=(Cart(2,jAtom)-Cart(2,iAtom))
            zij=(Cart(3,jAtom)-Cart(3,iAtom))
            rij2 = xij**2 + yij**2 + zij**2
            rrij=sqrt(rij2)

            alpha=rkf*exp((ami*r0mi**2+amj*r0mj**2))
            r=sqrt(rmj2+rmi2)
            gij=alpha*exp(-(ami*rmi2+amj*rmj2))
            rL2=(ymi*zmj-zmi*ymj)**2+(zmi*xmj-xmi*zmj)**2+&
     &         (xmi*ymj-ymi*xmj)**2
!hjw modified
            if(rL2.lt.1.d-14) then
              rL=0
            else
              rL=sqrt(rL2)
            end if
            if ((rmj.gt.rZero).and.(rmi.gt.rZero).and.&
     &                                (rrij.gt.rZero)) Then
              SinPhi=rL/(rmj*rmi)
              rmidotrmj=xmi*xmj+ymi*ymj+zmi*zmj
              CosPhi=rmidotrmj/(rmj*rmi)
!-------------None linear case
              If (SinPhi.gt.rZero) Then
                si(1)=(xmi/rmi*cosphi-xmj/rmj)/(rmi*sinphi)
                si(2)=(ymi/rmi*cosphi-ymj/rmj)/(rmi*sinphi)
                si(3)=(zmi/rmi*cosphi-zmj/rmj)/(rmi*sinphi)
                sj(1)=(cosphi*xmj/rmj-xmi/rmi)/(rmj*sinphi)
                sj(2)=(cosphi*ymj/rmj-ymi/rmi)/(rmj*sinphi)
                sj(3)=(cosphi*zmj/rmj-zmi/rmi)/(rmj*sinphi)
                sm(1)=-si(1)-sj(1)
                sm(2)=-si(2)-sj(2)
                sm(3)=-si(3)-sj(3)
                Do icoor=1,3
                   Do jCoor=1,3
                    If (mAtom.gt.iAtom) Then
                       Hess(Ind(icoor,mAtom,jcoor,iAtom))=&
     &                        Hess(Ind(icoor,mAtom,jcoor,iAtom))&
     &                        +gij*sm(icoor)*si(jcoor)
                    else
                      Hess(Ind(icoor,iAtom,jcoor,mAtom))=&
     &                        Hess(Ind(icoor,iAtom,jcoor,mAtom))&
     &                        +gij*si(icoor)*sm(jcoor)
                    End If
                    If (mAtom.gt.jAtom) Then
                        Hess(Ind(icoor,mAtom,jcoor,jAtom))=&
     &                        Hess(Ind(icoor,mAtom,jcoor,jAtom))&
     &                        +gij*sm(icoor)*sj(jcoor)
                    else
                      Hess(Ind(icoor,jAtom,jcoor,mAtom))=&
     &                        Hess(Ind(icoor,jAtom,jcoor,mAtom))&
     &                        +gij*sj(icoor)*sm(jcoor)
                    End If
                    If (iAtom.gt.jAtom) Then
                        Hess(Ind(icoor,iAtom,jcoor,jAtom))=&
     &                        Hess(Ind(icoor,iAtom,jcoor,jAtom))&
     &                        +gij*si(icoor)*sj(jcoor)
                     else
                        Hess(Ind(icoor,jAtom,jcoor,iAtom))=&
     &                        Hess(Ind(icoor,jAtom,jcoor,iAtom))&
     &                        +gij*sj(icoor)*si(jcoor)
                     End If
                   End Do
                End Do
                Do icoor=1,3
                  Do jCoor=1,icoor
                    Hess(Ind(icoor,iAtom,jcoor,iAtom))=&
     &                        Hess(Ind(icoor,iAtom,jcoor,iAtom))&
     &                        +gij*si(icoor)*si(jcoor)
                    Hess(Ind(icoor,mAtom,jcoor,mAtom))=&
     &                        Hess(Ind(icoor,mAtom,jcoor,mAtom))&
     &                        +gij*sm(icoor)*sm(jcoor)
                    Hess(Ind(icoor,jAtom,jcoor,jAtom))=&
     &                        Hess(Ind(icoor,jAtom,jcoor,jAtom))&
     &                        +gij*sj(icoor)*sj(jcoor)

                  End Do
                End Do
              Else
!----------------Linear case
                    if ((abs(ymi).gt.rZero).or.(abs(xmi).gt.rZero)) Then
                      x(1)=-ymi
                      y(1)=xmi
                      z(1)=Zero
                      x(2)=-xmi*zmi
                      y(2)=-ymi*zmi
                      z(2)=xmi*xmi+ymi*ymi
                    Else
                      x(1)=One
                      y(1)=Zero
                      z(1)=Zero
                      x(2)=Zero
                      y(2)=One
                      z(2)=Zero
                    End If
                    Do i=1,2
                     r1=sqrt(x(i)**2+y(i)**2+z(i)**2)
                     cosThetax=x(i)/r1
                     cosThetay=y(i)/r1
                     cosThetaz=z(i)/r1
                     si(1)=-cosThetax/rmi
                     si(2)=-cosThetay/rmi
                     si(3)=-cosThetaz/rmi
                     sj(1)=-cosThetax/rmj
                     sj(2)=-cosThetay/rmj
                     sj(3)=-cosThetaz/rmj
                     sm(1)=-(si(1)+sj(1))
                     sm(2)=-(si(2)+sj(2))
                     sm(3)=-(si(3)+sj(3))
                     Do icoor=1,3
                       Do jCoor=1,3
                        If (mAtom.gt.iAtom) Then
                          Hess(Ind(icoor,mAtom,jcoor,iAtom))=&
     &                        Hess(Ind(icoor,mAtom,jcoor,iAtom))&
     &                         +gij*sm(icoor)*si(jcoor)
                        else
                           Hess(Ind(icoor,iAtom,jcoor,mAtom))=&
     &                        Hess(Ind(icoor,iAtom,jcoor,mAtom))&
     &                         +gij*si(icoor)*sm(jcoor)
                        End If
                        If (mAtom.gt.jAtom) Then
                          Hess(Ind(icoor,mAtom,jcoor,jAtom))=&
     &                        Hess(Ind(icoor,mAtom,jcoor,jAtom))&
     &                         +gij*sm(icoor)*sj(jcoor)
                        else
                          Hess(Ind(icoor,jAtom,jcoor,mAtom))=&
     &                        Hess(Ind(icoor,jAtom,jcoor,mAtom))&
     &                         +gij*sj(icoor)*sm(jcoor)
                        End If
                        If (iAtom.gt.jAtom) Then
                           Hess(Ind(icoor,iAtom,jcoor,jAtom))=&
     &                        Hess(Ind(icoor,iAtom,jcoor,jAtom))&
     &                         +gij*si(icoor)*sj(jcoor)
                        else
                           Hess(Ind(icoor,jAtom,jcoor,iAtom))=&
     &                        Hess(Ind(icoor,jAtom,jcoor,iAtom))&
     &                         +gij*sj(icoor)*si(jcoor)
                        End If
                       End Do
                     End Do
                     Do icoor=1,3
                       Do jCoor=1,icoor
                         Hess(Ind(icoor,iAtom,jcoor,iAtom))=&
     &                        Hess(Ind(icoor,iAtom,jcoor,iAtom))&
     &                         +gij*si(icoor)*si(jcoor)
                         Hess(Ind(icoor,mAtom,jcoor,mAtom))=&
     &                        Hess(Ind(icoor,mAtom,jcoor,mAtom))&
     &                         +gij*sm(icoor)*sm(jcoor)
                         Hess(Ind(icoor,jAtom,jcoor,jAtom))=&
     &                         Hess(Ind(icoor,jAtom,jcoor,jAtom))&
     &                         +gij*sj(icoor)*sj(jcoor)
                       End Do
                     End Do
                 End Do
              End If
            End If
      End Do ! bend list

      if (profile) call timer%measure(2)

!
!     Hessian for torsion
!
      if (profile) call timer%measure(3,'torsion')
      do itors=1,topo%ntors
         iAtom=topo%tlist(3,itors)
         jAtom=topo%tlist(1,itors)
         kAtom=topo%tlist(2,itors)
         lAtom=topo%tlist(4,itors)

                 ir=iTabRow(iANr(iAtom))
                 jr=iTabRow(iANr(jAtom))
                 kr=iTabRow(iANr(kAtom))
                 lr=iTabRow(iANr(lAtom))
                 Call DCopy(3,Cart(1,iAtom),1,xyz(1,1),1)
                 Call DCopy(3,Cart(1,jAtom),1,xyz(1,2),1)
                 Call DCopy(3,Cart(1,kAtom),1,xyz(1,3),1)
                 Call DCopy(3,Cart(1,lAtom),1,xyz(1,4),1)
                 rij(1)=Cart(1,iAtom)-Cart(1,jAtom)
                 rij(2)=Cart(2,iAtom)-Cart(2,jAtom)
                 rij(3)=Cart(3,iAtom)-Cart(3,jAtom)
                 rij0=rAv(ir,jr)**2
                 aij =aAv(ir,jr)
                 rjk(1)=Cart(1,jAtom)-Cart(1,kAtom)
                 rjk(2)=Cart(2,jAtom)-Cart(2,kAtom)
                 rjk(3)=Cart(3,jAtom)-Cart(3,kAtom)
                 rjk0=rAv(jr,kr)**2
                 ajk =aAv(jr,kr)
                 rkl(1)=Cart(1,kAtom)-Cart(1,lAtom)
                 rkl(2)=Cart(2,kAtom)-Cart(2,lAtom)
                 rkl(3)=Cart(3,kAtom)-Cart(3,lAtom)
                 rkl0=rAv(kr,lr)**2
                 akl =aAv(kr,lr)
                 rij2=rij(1)**2+rij(2)**2+rij(3)**2
                 rjk2=rjk(1)**2+rjk(2)**2+rjk(3)**2
                 rkl2=rkl(1)**2+rkl(2)**2+rkl(3)**2

               beta=rkt*exp( (aij*rij0+ajk*rjk0+akl*rkl0))
               tij=beta*exp(-(aij*rij2+ajk*rjk2+akl*rkl2))

               Call Trsn(xyz,4,Tau,C,.False.,.False.,'        ',Dum,.False.)
               Call DCopy(3,C(1,1),1,si,1)
               Call DCopy(3,C(1,2),1,sj,1)
               Call DCopy(3,C(1,3),1,sk,1)
               Call DCopy(3,C(1,4),1,sl,1)
!-------------Off diagonal block
              Do icoor=1,3
                Do jCoor=1,3
                 Hess(Ind(icoor,iAtom,jcoor,jAtom))=&
     &           Hess(Ind(icoor,iAtom,jcoor,jAtom))&
     &            +tij*si(icoor) * sj(jcoor)
                 Hess(Ind(icoor,iAtom,jcoor,kAtom))=&
     &           Hess(Ind(icoor,iAtom,jcoor,kAtom))&
     &            +tij*si(icoor) * sk(jcoor)
                 Hess(Ind(icoor,iAtom,jcoor,lAtom))=&
     &           Hess(Ind(icoor,iAtom,jcoor,lAtom))&
     &            +tij*si(icoor) * sl(jcoor)
                 Hess(Ind(icoor,jAtom,jcoor,kAtom))=&
     &           Hess(Ind(icoor,jAtom,jcoor,kAtom))&
     &            +tij*sj(icoor) * sk(jcoor)
                 Hess(Ind(icoor,jAtom,jcoor,lAtom))=&
     &           Hess(Ind(icoor,jAtom,jcoor,lAtom))&
     &            +tij*sj(icoor) * sl(jcoor)
                 Hess(Ind(icoor,kAtom,jcoor,lAtom))=&
     &           Hess(Ind(icoor,kAtom,jcoor,lAtom))&
     &            +tij*sk(icoor) * sl(jcoor)

                End Do
              End Do
!-------------Diagonal block
              Do icoor=1,3
                Do jCoor=1,icoor
                 Hess(Ind(icoor,iAtom,jcoor,iAtom))=&
     &           Hess(Ind(icoor,iAtom,jcoor,iAtom))&
     &            +tij*si(icoor) * si(jcoor)
                 Hess(Ind(icoor,jAtom,jcoor,jAtom))=&
     &           Hess(Ind(icoor,jAtom,jcoor,jAtom))&
     &            +tij*sj(icoor) * sj(jcoor)
                 Hess(Ind(icoor,kAtom,jcoor,kAtom))=&
     &           Hess(Ind(icoor,kAtom,jcoor,kAtom))&
     &            +tij*sk(icoor) * sk(jcoor)
                 Hess(Ind(icoor,lAtom,jcoor,lAtom))=&
     &           Hess(Ind(icoor,lAtom,jcoor,lAtom))&
     &            +tij*sl(icoor) * sl(jcoor)

                 End Do
               End Do
 
      enddo ! tors list

      if (profile) call timer%measure(3)
      if (profile) call timer%write(6,'modhes')
      Return
      End

! Coulomb + dispersion Hessian (SG, 5/19)
! d/(dxdx)  q*q/(damp+r),  rpa=damp+r
! d/(dxdx)   c6/(r02^3+r^6)
      subroutine getqqxx(dx,qq,c6,r,r2,r3,rpa,r02,d2)      
      implicit none
      real*8 dx,qq,c6,r,r2,r3,rpa,r02,d2
      real*8 rpa2,dx2,ar6,r6,r8

!     Coulomb part
      rpa2=rpa**2
      dx2 = dx**2
      d2  = qq*(2.d0*dx2/(r2*rpa*rpa2) + dx2/(r3*rpa2) - 1.d0/(r*rpa2))
!     disp part
      r6  = r3*r3
      r8  = r6*r2
      ar6 = r02**3+r6
      d2  = d2 + c6*(dx2*72.d0*r8/ar6**3 - dx2*24.d0*r2/ar6**2 &
     &         - 6.d0*r2*r2/ar6**2)
      end

! d/(dxdy)  q*q/(damp+r),  rpa=damp+r
! d/(dxdy)   c6/(r02^3+r^6)
      subroutine getqqxy(dx,dy,qq,c6,r,r2,r3,rpa,r02,d2)      
      implicit none
      real*8 dx,dy,qq,c6,r,r2,r3,rpa,r02,d2
      real*8 rpa2,r6,r8,ar6

!     Coulomb part
      rpa2=rpa**2
      d2  =qq*(2.d0*dx*dy/(r2*rpa*rpa2) + dx*dy/(r3*rpa2))
!     disp part
      r6  = r3*r3
      r8  = r6*r2
      ar6 = r02**3+r6
      d2  = d2 + c6*(dx*dy*72.d0*r8/ar6**3 - dx*dy*24.d0*r2/ar6**2)

      end
