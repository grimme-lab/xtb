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

!! ========================================================================
!  GENERAL FUNCTIONS FOR CORE FUNCTIONALITIES OF THE SCC
!! ------------------------------------------------------------------------
!  GFN1:
!  -> build_h0_gfn1
!  GFN2:
!  -> build_h0_gfn2
!! ========================================================================
module scc_core
   use iso_fortran_env, only : wp => real64
   use mctc_la, only : sygvd,gemm,symm
   implicit none

   integer, parameter :: mmm(*)=(/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/)
 
contains

!! ========================================================================
!  build GFN1 core Hamiltonian
!! ========================================================================
subroutine build_h0_gfn1(H0,n,at,ndim,nmat,matlist,kspd,kmagic,kenscal, &
   &                     xyz,cn,kcnao,S,aoat2,lao2,valao2,hdiag2)
   implicit none
   real(wp),intent(out) :: H0(ndim*(ndim+1)/2)
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: ndim
   integer, intent(in)  :: nmat
   integer, intent(in) :: matlist(2,nmat)
   real(wp),intent(in)  :: kspd(6)
   real(wp),intent(in)  :: kmagic(4,4)
   real(wp),intent(in)  :: kenscal
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: cn(n)
   real(wp),intent(in)  :: kcnao(ndim)
   real(wp),intent(in)  :: S(ndim,ndim)
   integer, intent(in)  :: aoat2(ndim)
   integer, intent(in)  :: lao2(ndim)
   integer, intent(in)  :: valao2(ndim)
   real(wp),intent(in)  :: hdiag2(ndim)

   integer  :: i,j,k,m
   integer  :: iat,jat,ishell,jshell
   real(wp) :: hdii,hdjj,hav
   real(wp) :: km

   H0=0.0_wp
   do m=1,nmat
      i=matlist(1,m)
      j=matlist(2,m)
      k=j+i*(i-1)/2
      iat=aoat2(i)
      jat=aoat2(j)
      ishell=mmm(lao2(i))
      jshell=mmm(lao2(j))
      hdii=hdiag2(i)
      hdii=hdii*(1.0d0+kcnao(i)*cn(iat))  ! CN dependent shift
      hdjj=hdiag2(j)
      hdjj=hdjj*(1.0d0+kcnao(j)*cn(jat))  ! CN dependent shift
      call h0scal(n,at,i,j,ishell,jshell,iat,jat,valao2(i).ne.0,valao2(j).ne.0, &
      &           kspd,kmagic,kenscal,km)
      hav=0.5d0*(hdii+hdjj)* &
      &      rfactor(ishell,jshell,at(iat),at(jat),xyz(:,iat),xyz(:,jat))
      H0(k)=S(j,i)*km*hav
   enddo
!  diagonal
   k=0
   do i=1,ndim
      k=k+i
      iat=aoat2(i)
      ishell=mmm(lao2(i))
      H0(k)=hdiag2(i)*(1.0d0+kcnao(i)*cn(iat))  ! CN dependent shift
   enddo

end subroutine build_h0_gfn1

!! ========================================================================
!  build GFN2 core Hamiltonian
!! ========================================================================
subroutine build_h0_gfn2(H0,n,at,ndim,nmat,matlist,kspd,kmagic,kenscal, &
   &                     xyz,cn,kcnao,S,aoat2,lao2,valao2,hdiag2,aoexp)
   implicit none
   real(wp),intent(out) :: H0(ndim*(ndim+1)/2)
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: ndim
   integer, intent(in)  :: nmat
   integer, intent(in)  :: matlist(2,nmat)
   real(wp),intent(in)  :: kspd(6)
   real(wp),intent(in)  :: kmagic(4,4)
   real(wp),intent(in)  :: kenscal
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: cn(n)
   real(wp),intent(in)  :: kcnao(ndim)
   real(wp),intent(in)  :: S(ndim,ndim)
   integer, intent(in)  :: aoat2(ndim)
   integer, intent(in)  :: lao2(ndim)
   integer, intent(in)  :: valao2(ndim)
   real(wp),intent(in)  :: hdiag2(ndim)
   real(wp),intent(in)  :: aoexp(ndim)

   integer  :: i,j,k,m
   integer  :: iat,jat,ishell,jshell
   real(wp) :: hdii,hdjj,hav
   real(wp) :: km
   real(wp),parameter :: aot = -0.5d0 ! AO exponent dep. H0 scal

   H0=0.0_wp

   do m=1,nmat
      i=matlist(1,m)
      j=matlist(2,m)
      k=j+i*(i-1)/2
      iat=aoat2(i)
      jat=aoat2(j)
      ishell=mmm(lao2(i))
      jshell=mmm(lao2(j))
      hdii=hdiag2(i)
      hdii=hdii-kcnao(i)*cn(iat)  ! CN dependent shift
      hdjj=hdiag2(j)
      hdjj=hdjj-kcnao(j)*cn(jat)  ! CN dependent shift
      call h0scal(n,at,i,j,ishell,jshell,iat,jat,valao2(i).ne.0,valao2(j).ne.0, &
      &           kspd,kmagic,kenscal,km)
      km=km*(0.5*((aoexp(i)+aoexp(j))/(aoexp(i)*aoexp(j))**0.5))**aot
      hav=0.5d0*(hdii+hdjj)* &
      &      rfactor(ishell,jshell,at(iat),at(jat),xyz(:,iat),xyz(:,jat))
      H0(k)=S(j,i)*km*hav
   enddo
!  diagonal
   k=0
   do i=1,ndim
      k=k+i
      iat=aoat2(i)
      ishell=mmm(lao2(i))
      H0(k)=hdiag2(i)-kcnao(i)*cn(iat)  ! CN dependent shift
   enddo

end subroutine build_h0_gfn2

!! ========================================================================
!  build GFN1 Fockian
!! ========================================================================
subroutine build_h1_gfn1(n,at,ndim,nshell,nmat,matlist,H,H1,H0,S,ves,q, &
                         cm5,fgb,fhb,aoat2,ao2sh)
   use mctc_econv, only : autoev,evtoau
   use aoparam,  only : gam3
   use gbobc, only : lgbsa
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: ndim
   integer, intent(in)  :: nshell
   integer, intent(in)  :: nmat
   integer, intent(in)  :: matlist(2,nmat)
   real(wp),intent(in)  :: H0(ndim*(1+ndim)/2)
   real(wp),intent(in)  :: S(ndim,ndim)
   real(wp),intent(in)  :: ves(nshell)
   real(wp),intent(in)  :: q(n)
   real(wp),intent(in)  :: cm5(n)
   real(wp),intent(in)  :: fgb(n,n)
   real(wp),intent(in)  :: fhb(n)
   integer, intent(in)  :: aoat2(ndim)
   integer, intent(in)  :: ao2sh(ndim)
   real(wp),intent(out) :: H(ndim,ndim)
   real(wp),intent(out) :: H1(ndim*(1+ndim)/2)

   integer  :: m,i,j,k
   integer  :: ishell,jshell
   integer  :: ii,jj,kk
   real(wp) :: dum
   real(wp) :: eh1,t8,t9,tgb

   H = 0.0_wp
   H1 = 0.0_wp

   do m = 1, nmat
      i = matlist(1,m)
      j = matlist(2,m)
      k = j+i*(i-1)/2
      ishell = ao2sh(i)
      jshell = ao2sh(j)
      dum = S(j,i)
!     SCC terms
!     2nd order ES term (optional: including point charge potential)
      eh1 = ves(ishell) + ves(jshell)
!     3rd order and set-up of H
      ii = aoat2(i)
      jj = aoat2(j)
      dum = S(j,i)
!     third-order diagonal term, unscreened
      t8 = q(ii)**2 * gam3(at(ii))
      t9 = q(jj)**2 * gam3(at(jj))
      eh1 = eh1 + autoev*(t8+t9)
      H1(k) = -dum*eh1*0.5_wp
      H(j,i) = H0(k)+H1(k)
      H(i,j) = H(j,i)
   enddo
!  add the gbsa SCC term
   if (lgbsa) then
!     hbpow=2.d0*c3-1.d0
      do m=1,nmat
         i=matlist(1,m)
         j=matlist(2,m)
         k=j+i*(i-1)/2
         ii=aoat2(i)
         jj=aoat2(j)
         dum=S(j,i)
!        GBSA SCC terms
         eh1=0.0_wp
         do kk=1,n
            eh1=eh1+cm5(kk)*(fgb(kk,ii)+fgb(kk,jj))
         enddo
         t8=fhb(ii)*cm5(ii)+fhb(jj)*cm5(jj)
         tgb=-dum*(0.5_wp*eh1+t8)
         H1(k)=H1(k)+tgb
         H(j,i)=H(j,i)+tgb
         H(i,j)=H(j,i)
      enddo
   endif

end subroutine build_h1_gfn1

!! ========================================================================
!  build GFN1 Fockian
!! ========================================================================
subroutine build_h1_gfn2(n,at,ndim,nshell,nmat,ndp,nqp,matlist,mdlst,mqlst,&
                         H,H1,H0,S,dpint,qpint,ves,vs,vd,vq,q,qsh,gam3sh, &
                         hdisp,fgb,fhb,aoat2,ao2sh)
   use mctc_econv, only : autoev,evtoau
   use gbobc, only : lgbsa
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: ndim
   integer, intent(in)  :: nshell
   integer, intent(in)  :: nmat
   integer, intent(in)  :: ndp
   integer, intent(in)  :: nqp
   integer, intent(in)  :: matlist(2,nmat)
   integer, intent(in)  :: mdlst(2,ndp)
   integer, intent(in)  :: mqlst(2,nqp)
   real(wp),intent(in)  :: H0(ndim*(1+ndim)/2)
   real(wp),intent(in)  :: S(ndim,ndim)
   real(wp),intent(in)  :: dpint(3,ndim*(1+ndim)/2)
   real(wp),intent(in)  :: qpint(6,ndim*(1+ndim)/2)
   real(wp),intent(in)  :: ves(nshell)
   real(wp),intent(in)  :: vs(n)
   real(wp),intent(in)  :: vd(3,n)
   real(wp),intent(in)  :: vq(6,n)
   real(wp),intent(in)  :: q(n)
   real(wp),intent(in)  :: qsh(nshell)
   real(wp),intent(in)  :: gam3sh(nshell)
   real(wp),intent(in)  :: hdisp(n)
   real(wp),intent(in)  :: fgb(n,n)
   real(wp),intent(in)  :: fhb(n)
   integer, intent(in)  :: aoat2(ndim)
   integer, intent(in)  :: ao2sh(ndim)
   real(wp),intent(out) :: H(ndim,ndim)
   real(wp),intent(out) :: H1(ndim*(1+ndim)/2)

   integer, external :: lin
   integer  :: m,i,j,k,l
   integer  :: ii,jj,kk
   integer  :: ishell,jshell
   real(wp) :: dum,eh1,t8,t9,tgb

   H1=0.0_wp
   H =0.0_wp
! --- set up of Fock matrix
!  overlap dependent terms
!  on purpose, vs is NOT added to H1 (gradient is calculated separately)
   do m=1,nmat
      i=matlist(1,m)
      j=matlist(2,m)
      k=j+i*(i-1)/2
      ii=aoat2(i)
      jj=aoat2(j)
      ishell=ao2sh(i)
      jshell=ao2sh(j)
      dum=S(j,i)
!     SCC terms
!     2nd order ES term (optional: including point charge potential)
      eh1=ves(ishell)+ves(jshell)
!     SIE term
!     t6=gsie(at(ii))*pi*sin(2.0d0*pi*q(ii))
!     t7=gsie(at(jj))*pi*sin(2.0d0*pi*q(jj))
!     eh1=eh1+autoev*(t6+t7)*gscal
!     third order term
      t8=qsh(ishell)**2*gam3sh(ishell)
      t9=qsh(jshell)**2*gam3sh(jshell)
      eh1=eh1+autoev*(t8+t9)
      H1(k)=-dum*eh1*0.5d0
! SAW start - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1801
!     Dispersion contribution to Hamiltionian
      H1(k)=H1(k) - 0.5d0*dum*autoev*(hdisp(ii)+hdisp(jj))
! SAW end - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1801
      H(j,i)=H0(k)+H1(k)
!     CAMM potential
      eh1=0.50d0*dum*(vs(ii)+vs(jj))*autoev
      H(j,i)=H(j,i)+eh1
      H(i,j)=H(j,i)
   enddo
!  quadrupole-dependent terms
   do m=1,nqp
      i=mqlst(1,m)
      j=mqlst(2,m)
      ii=aoat2(i)
      jj=aoat2(j)
      k=lin(j,i)
      eh1=0.0d0
      ! note: these come in the following order
      ! xx, yy, zz, xy, xz, yz
      do l=1,6
         eh1=eh1+qpint(l,k)*(vq(l,ii)+vq(l,jj))
      enddo
      eh1=0.50d0*eh1*autoev
!     purposely, do NOT add dip/qpole-int terms onto H1
!     (due to gradient computation later on)
      H(i,j)=H(i,j)+eh1
      H(j,i)=H(i,j)
   enddo
!  dipolar terms
   do m=1,ndp
      i=mdlst(1,m)
      j=mdlst(2,m)
      k=lin(j,i)
      ii=aoat2(i)
      jj=aoat2(j)
      eh1=0.0d0
      do l=1,3
         eh1=eh1+dpint(l,k)*(vd(l,ii)+vd(l,jj))
      enddo
      eh1=0.50d0*eh1*autoev
!     purposely, do NOT add dip/qpole-int terms onto H1
!     (due to gradient computation later on)
      H(i,j)=H(i,j)+eh1
      H(j,i)=H(i,j)
   enddo
!                                          call timing(t2,w2)
!                           call prtime(6,t2-t1,w2-w1,'Fmat')
!  add the gbsa SCC term
   if(lgbsa) then
!     hbpow=2.d0*c3-1.d0
      do m=1,nmat
         i=matlist(1,m)
         j=matlist(2,m)
         k=j+i*(i-1)/2
         ii=aoat2(i)
         jj=aoat2(j)
         dum=S(j,i)
!        GBSA SCC terms
         eh1=0.0d0
         do kk=1,n
            eh1=eh1+q(kk)*(fgb(kk,ii)+fgb(kk,jj))
         enddo
         t8=fhb(ii)*q(ii)+fhb(jj)*q(jj)
         tgb=-dum*(0.5d0*eh1+t8)
         H1(k)=H1(k)+tgb
         H(j,i)=H(j,i)+tgb
         H(i,j)=H(j,i)
      enddo
   endif

end subroutine build_h1_gfn2

!! ========================================================================
!  self consistent charge iterator for GFN1 Hamiltonian
!! ========================================================================
subroutine scc_gfn1(iunit,n,nel,nopen,ndim,nmat,nshell, &
   &                at,matlist,aoat2,ao2sh, &
   &                q,qq,qlmom,qsh,zsh, &
   &                gbsa,fgb,fhb,cm5,cm5a,gborn, &
   &                broy,broydamp,damp0, &
   &                pcem,ves,vpc, &
   &                et,focc,focca,foccb,efa,efb, &
   &                eel,ees,epcem,egap,emo,ihomo,ihomoa,ihomob, &
   &                H0,H1,H,S,X,P,jab, &
   &                maxiter,startpdiag,scfconv,qconv, &
   &                minpr,pr, &
   &                fail,jter)
   use mctc_econv, only : autoev,evtoau

   use aoparam,  only : gam3

   use gbobc, only : lgbsa,lhb,tb_solvent
   use embedding, only : electro_pcem

   implicit none

   integer, intent(in)  :: iunit

   integer, intent(in)  :: n
   integer, intent(in)  :: nel
   integer, intent(in)  :: nopen
   integer, intent(in)  :: ndim
   integer, intent(in)  :: nmat
   integer, intent(in)  :: nshell
!! ------------------------------------------------------------------------
!  general options for the iterator
   integer, intent(in)  :: maxiter
   integer, intent(in)  :: startpdiag
   real(wp),intent(in)  :: scfconv
   real(wp),intent(in)  :: qconv
   logical, intent(in)  :: minpr
   logical, intent(in)  :: pr
   logical, intent(out) :: fail
!! ------------------------------------------------------------------------
   integer, intent(in)  :: at(n)
!! ------------------------------------------------------------------------
   integer, intent(in)  :: matlist(2,nmat)
   integer, intent(in)  :: aoat2(ndim)
   integer, intent(in)  :: ao2sh(ndim)
!! ------------------------------------------------------------------------
!  a bunch of charges
   real(wp),intent(inout) :: q(n)
   real(wp),intent(inout) :: qq(n)
   real(wp),intent(inout) :: qlmom(3,n)
   real(wp),intent(inout) :: qsh(nshell)
   real(wp),intent(in)    :: zsh(nshell)
!! ------------------------------------------------------------------------
!  continuum solvation model GBSA
   type(tb_solvent),intent(inout) :: gbsa
   real(wp),intent(inout) :: fgb(n,n)
   real(wp),intent(inout) :: fhb(n)
   real(wp),intent(in)    :: cm5a(n)
   real(wp),intent(inout) :: cm5(n)
   real(wp),intent(inout) :: gborn
!! ------------------------------------------------------------------------
!  point charge embedding potentials
   logical, intent(in)    :: pcem
   real(wp),intent(inout) :: ves(nshell)
   real(wp),intent(inout) :: vpc(nshell)
!! ------------------------------------------------------------------------
!  Fermi-smearing
   real(wp),intent(in)    :: et
   real(wp),intent(inout) :: focc(ndim)
   real(wp),intent(inout) :: foccb(ndim),focca(ndim)
   real(wp),intent(inout) :: efa,efb
!! ------------------------------------------------------------------------
!  Convergence accelerators, a simple damping as well as a Broyden mixing
!  are available. The Broyden mixing is used by default seems reliable.
   real(wp),intent(in)    :: damp0
   real(wp)               :: damp
!  Broyden
   logical, intent(in)    :: broy
   real(wp),intent(inout) :: broydamp
   real(wp)               :: omegap
   real(wp),allocatable   :: df(:,:)
   real(wp),allocatable   :: u(:,:)
   real(wp),allocatable   :: a(:,:)
   real(wp),allocatable   :: q_in(:)
   real(wp),allocatable   :: dq(:)
   real(wp),allocatable   :: qlast_in(:)
   real(wp),allocatable   :: dqlast(:)
   real(wp),allocatable   :: omega(:)
!! ------------------------------------------------------------------------
!  results of the SCC iterator
   real(wp),intent(out)   :: eel
   real(wp),intent(out)   :: epcem
   real(wp),intent(out)   :: ees
   real(wp),intent(out)   :: egap
   real(wp),intent(out)   :: emo(ndim)
   integer, intent(inout) :: ihomoa
!! ------------------------------------------------------------------------
   real(wp),intent(in)    :: H0(ndim*(ndim+1)/2)
   real(wp),intent(out)   :: H1(ndim*(ndim+1)/2)
   real(wp),intent(out)   :: H(ndim,ndim)
   real(wp),intent(inout) :: P(ndim,ndim)
   real(wp),intent(inout) :: X(ndim,ndim)
   real(wp),intent(in)    :: S(ndim,ndim)
   real(wp),intent(inout) :: jab(nshell,nshell)

   integer, intent(inout) :: jter
!! ------------------------------------------------------------------------
!  local variables
   integer,external :: lin
   integer  :: i,ii,j,jj,k,kk,l,m
   integer  :: ishell,jshell
   integer  :: ihomo,ihomob
   real(wp) :: t8,t9
   real(wp) :: eh1,dum,tgb
   real(wp) :: eold
   real(wp) :: ga,gb
   real(wp) :: rmsq
   real(wp) :: nfoda,nfodb
   logical  :: fulldiag
   logical  :: lastdiag
   integer  :: iter
   integer  :: thisiter
   logical  :: converged
   logical  :: econverged
   logical  :: qconverged

   converged = .false.
   lastdiag = .false.
   ! number of iterations for this iterator
   thisiter = maxiter - jter

   damp = damp0
!  broyden data storage and init
   allocate( df(thisiter,nshell),u(thisiter,nshell), &
   &         a(thisiter,thisiter),dq(nshell),dqlast(nshell), &
   &         qlast_in(nshell),omega(thisiter),q_in(nshell), &
   &         source = 0.0_wp )

!! ------------------------------------------------------------------------
!  iteration entry point
   scc_iterator: do iter = 1, thisiter
!! ------------------------------------------------------------------------
!  build the Fockian from current ES potential and partial charges
!  includes GBSA contribution to Fockian
   call build_H1_gfn1(n,at,ndim,nshell,nmat,matlist,H,H1,H0,S,ves,q, &
                      cm5,fgb,fhb,aoat2,ao2sh)

!! ------------------------------------------------------------------------
!  solve HC=SCemo(X,P are scratch/store)
!  solution is in H(=C)/emo
!! ------------------------------------------------------------------------
   fulldiag=.false.
   if (iter.lt.startpdiag) fulldiag=.true.
   if (lastdiag )          fulldiag=.true.
   call solve(fulldiag,ndim,ihomo,scfconv,H,S,X,P,emo,fail)

   if (fail) then
      eel = 1.e+99_wp
      return
   endif

   if ((ihomo+1.le.ndim).and.(ihomo.ge.1)) egap = emo(ihomo+1)-emo(ihomo)
!  automatic reset to small value
   if ((egap.lt.0.1_wp).and.(iter.eq.0)) broydamp = 0.03_wp

!! ------------------------------------------------------------------------
!  Fermi smearing
   if (et.gt.0.1_wp) then
!     convert restricted occ first to alpha/beta
      if(nel.gt.0) then
         call occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
      else
         focca=0.0_wp
         foccb=0.0_wp
         ihomoa=0
         ihomob=0
      endif
      if(ihomoa+1.le.ndim) then
         call fermismear(.false.,ndim,ihomoa,et,emo,focca,nfoda,efa,ga)
      endif
      if(ihomob+1.le.ndim) then
         call fermismear(.false.,ndim,ihomob,et,emo,foccb,nfodb,efb,gb)
      endif
      focc = focca + foccb
   else
      ga = 0.0_wp
      gb = 0.0_wp
   endif
!! ------------------------------------------------------------------------

!  save q
   q_in(1:nshell)=qsh(1:nshell)
   k=nshell

!  density matrix
   call dmat(ndim,focc,H,P)

!  new q
   call mpopsh (n,ndim,nshell,ao2sh,S,P,qsh)
   qsh = zsh - qsh

!  qat from qsh
   call qsh2qat(n,at,nshell,qsh,q)

   eold=eel
   call electro(n,at,ndim,nshell,jab,H0,P,q,qsh,ees,eel)

!  point charge contribution
   if (pcem) call electro_pcem(nshell,qsh,Vpc,epcem,eel)

!  new cm5 charges and gborn energy
   if(lgbsa) then
      cm5=q+cm5a
      call electro_gbsa(n,at,fgb,fhb,cm5,gborn,eel)
   endif

!  ad el. entropies*T
   eel=eel+ga+gb

!! ------------------------------------------------------------------------
!  check for energy convergence
   econverged = abs(eel - eold) < scfconv
!! ------------------------------------------------------------------------

   dq(1:nshell)=qsh(1:nshell)-q_in(1:nshell)
   rmsq=sum(dq(1:nshell)**2)/dble(n)
   rmsq=sqrt(rmsq)

!! ------------------------------------------------------------------------
!  end of SCC convergence part
   qconverged = rmsq < qconv
!! ------------------------------------------------------------------------

!  SCC convergence acceleration
   if (.not.broy) then

!     simple damp
      if(iter.gt.0) then
         omegap=egap
         ! monopoles only
         do i=1,nshell
            qsh(i)=damp*qsh(i)+(1.0_wp-damp)*q_in(i)
         enddo
         if(eel-eold.lt.0) then
            damp=damp*1.15_wp
         else
            damp=damp0
         endif
         damp=min(damp,1.0_wp)
         if (egap.lt.1.0_wp) damp=min(damp,0.5_wp)
      endif

   else

!     Broyden mixing
      omegap=0.0_wp
      call broyden(nshell,q_in,qlast_in,dq,dqlast, &
      &            iter,thisiter,broydamp,omega,df,u,a)
      qsh(1:nshell)=q_in(1:nshell)
      if(iter.gt.1) omegap=omega(iter-1)
   endif ! Broyden?

   call qsh2qat(n,at,nshell,qsh,q) !new qat

   if(minpr) write(iunit,'(i4,F15.7,E14.6,E11.3,f8.2,2x,f8.1,l3)') &
   &         iter+jter,eel,eel-eold,rmsq,egap,omegap,fulldiag
   qq=q

   if(lgbsa) cm5=q+cm5a
!  set up ES potential
   if(pcem) then
      ves(1:nshell)=Vpc(1:nshell)
   else
      ves=0.0_wp
   endif
   call setespot(nshell,qsh,jab,ves)

!! ------------------------------------------------------------------------
   if (econverged.and.qconverged) then
      converged = .true.
      if (lastdiag) exit scc_iterator
      lastdiag = .true.
   endif
!! ------------------------------------------------------------------------

   enddo scc_iterator

   jter = jter + min(iter,maxiter-jter)
   fail = .not.converged

end subroutine scc_gfn1

!! ========================================================================
!  self consistent charge iterator for GFN2 Hamiltonian
!! ========================================================================
subroutine scc_gfn2(iunit,n,nel,nopen,ndim,ndp,nqp,nmat,nshell, &
   &                at,matlist,mdlst,mqlst,aoat2,ao2sh, &
   &                q,dipm,qp,qq,qlmom,qsh,zsh, &
   &                xyz,vs,vd,vq,gab3,gab5,gscal, &
   &                gbsa,fgb,fhb,cm5,cm5a,gborn, &
   &                newdisp,dispdim,g_a,g_c,gw,wdispmat,hdisp, &
   &                broy,broydamp,damp0, &
   &                pcem,ves,vpc, &
   &                et,focc,focca,foccb,efa,efb, &
   &                eel,ees,eaes,epol,ed,epcem,egap,emo,ihomo,ihomoa,ihomob, &
   &                H0,H1,H,S,dpint,qpint,X,P,jab,gam3sh, &
   &                maxiter,startpdiag,scfconv,qconv, &
   &                minpr,pr, &
   &                fail,jter)
   use mctc_econv, only : autoev,evtoau

   use aoparam,  only : gam3

   use gbobc,  only : lgbsa,lhb,tb_solvent
   use dftd4,  only : disppot,edisp_scc
   use aespot, only : gfn2broyden_diff,gfn2broyden_out,gfn2broyden_save, &
   &                  mmompop,aniso_electro,setvsdq
   use embedding, only : electro_pcem

   implicit none

   integer, intent(in)  :: iunit

   integer, intent(in)  :: n
   integer, intent(in)  :: nel
   integer, intent(in)  :: nopen
   integer, intent(in)  :: ndim
   integer, intent(in)  :: ndp
   integer, intent(in)  :: nqp
   integer, intent(in)  :: nmat
   integer, intent(in)  :: nshell
!! ------------------------------------------------------------------------
!  general options for the iterator
   integer, intent(in)  :: maxiter
   integer, intent(in)  :: startpdiag
   real(wp),intent(in)  :: scfconv
   real(wp),intent(in)  :: qconv
   logical, intent(in)  :: minpr
   logical, intent(in)  :: pr
   logical, intent(out) :: fail
!! ------------------------------------------------------------------------
   integer, intent(in)  :: at(n)
!! ------------------------------------------------------------------------
   integer, intent(in)  :: matlist(2,nmat)
   integer, intent(in)  :: mdlst(2,ndp)
   integer, intent(in)  :: mqlst(2,nqp)
   integer, intent(in)  :: aoat2(ndim)
   integer, intent(in)  :: ao2sh(ndim)
!! ------------------------------------------------------------------------
!  a bunch of charges and CAMMs
   real(wp),intent(inout) :: q(n)
   real(wp),intent(inout) :: dipm(3,n)
   real(wp),intent(inout) :: qp(6,n)
   real(wp),intent(inout) :: qq(n)
   real(wp),intent(inout) :: qlmom(3,n)
   real(wp),intent(inout) :: qsh(nshell)
   real(wp),intent(in)    :: zsh(nshell)
!! ------------------------------------------------------------------------
!  anisotropic electrostatic
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: vs(n)
   real(wp),intent(inout) :: vd(3,n)
   real(wp),intent(inout) :: vq(6,n)
   real(wp),intent(inout) :: gab3(n*(n+1)/2)
   real(wp),intent(inout) :: gab5(n*(n+1)/2)
   real(wp),intent(in)    :: gscal
!! ------------------------------------------------------------------------
!  continuum solvation model GBSA
   type(tb_solvent),intent(inout) :: gbsa
   real(wp),intent(inout) :: fgb(n,n)
   real(wp),intent(inout) :: fhb(n)
   real(wp),intent(in)    :: cm5a(n)
   real(wp),intent(inout) :: cm5(n)
   real(wp),intent(inout) :: gborn
!! ------------------------------------------------------------------------
!  selfconsistent DFT-D4 dispersion correction
   logical, intent(in)    :: newdisp
   integer, intent(in)    :: dispdim
   real(wp),intent(in)    :: g_a,g_c
   real(wp),intent(in)    :: gw(dispdim)
   real(wp),intent(in)    :: wdispmat(dispdim,dispdim)
   real(wp),intent(inout) :: hdisp(n)
!! ------------------------------------------------------------------------
!  point charge embedding potentials
   logical, intent(in)    :: pcem
   real(wp),intent(inout) :: ves(nshell)
   real(wp),intent(inout) :: vpc(nshell)
!! ------------------------------------------------------------------------
!  Fermi-smearing
   real(wp),intent(in)    :: et
   real(wp),intent(inout) :: focc(ndim)
   real(wp),intent(inout) :: foccb(ndim),focca(ndim)
   real(wp),intent(inout) :: efa,efb
!! ------------------------------------------------------------------------
!  Convergence accelerators, a simple damping as well as a Broyden mixing
!  are available. The Broyden mixing is used by default seems reliable.
   real(wp),intent(in)    :: damp0
   real(wp)               :: damp
!  Broyden
   integer                :: nbr
   logical, intent(in)    :: broy
   real(wp),intent(inout) :: broydamp
   real(wp)               :: omegap
   real(wp),allocatable   :: df(:,:)
   real(wp),allocatable   :: u(:,:)
   real(wp),allocatable   :: a(:,:)
   real(wp),allocatable   :: q_in(:)
   real(wp),allocatable   :: dq(:)
   real(wp),allocatable   :: qlast_in(:)
   real(wp),allocatable   :: dqlast(:)
   real(wp),allocatable   :: omega(:)
!! ------------------------------------------------------------------------
!  results of the SCC iterator
   real(wp),intent(out)   :: eel
   real(wp),intent(out)   :: epcem
   real(wp),intent(out)   :: ees
   real(wp),intent(out)   :: eaes
   real(wp),intent(out)   :: epol
   real(wp),intent(out)   :: ed
   real(wp),intent(out)   :: egap
   real(wp),intent(out)   :: emo(ndim)
   integer, intent(inout) :: ihomo
   integer, intent(inout) :: ihomoa
   integer, intent(inout) :: ihomob
!! ------------------------------------------------------------------------
   real(wp),intent(in)    :: H0(ndim*(ndim+1)/2)
   real(wp),intent(out)   :: H1(ndim*(ndim+1)/2)
   real(wp),intent(out)   :: H(ndim,ndim)
   real(wp),intent(inout) :: P(ndim,ndim)
   real(wp),intent(inout) :: X(ndim,ndim)
   real(wp),intent(in)    :: S(ndim,ndim)
   real(wp),intent(in)    :: dpint(3,ndim*(ndim+1)/2)
   real(wp),intent(in)    :: qpint(6,ndim*(ndim+1)/2)
   real(wp),intent(inout) :: jab(nshell,nshell)
   real(wp),intent(in)    :: gam3sh(nshell)

   integer, intent(inout) :: jter
!! ------------------------------------------------------------------------
!  local variables
   integer,external :: lin
   integer  :: i,ii,j,jj,k,kk,l,m
   integer  :: ishell,jshell
   real(wp) :: t8,t9
   real(wp) :: eh1,dum,tgb
   real(wp) :: eold
   real(wp) :: ga,gb
   real(wp) :: rmsq
   real(wp) :: nfoda,nfodb
   logical  :: fulldiag
   logical  :: lastdiag
   integer  :: iter
   integer  :: thisiter
   logical  :: converged
   logical  :: econverged
   logical  :: qconverged

   converged = .false.
   lastdiag = .false.
   ! number of iterations for this iterator
   thisiter = maxiter - jter

   damp = damp0
   nbr = nshell + 9*n
!  broyden data storage and init
   allocate( df(thisiter,nbr),u(thisiter,nbr),a(thisiter,thisiter), &
   &         dq(nbr),dqlast(nbr),qlast_in(nbr),omega(thisiter), &
   &         q_in(nbr), source = 0.0_wp )

!! ------------------------------------------------------------------------
!  Iteration entry point
   scc_iterator: do iter = 1, thisiter
!! ------------------------------------------------------------------------
   call build_h1_gfn2(n,at,ndim,nshell,nmat,ndp,nqp,matlist,mdlst,mqlst,&
                      H,H1,H0,S,dpint,qpint,ves,vs,vd,vq,q,qsh,gam3sh, &
                      hdisp,fgb,fhb,aoat2,ao2sh)

!! ------------------------------------------------------------------------
!  solve HC=SCemo(X,P are scratch/store)
!  solution is in H(=C)/emo
!! ------------------------------------------------------------------------
   fulldiag=.false.
   if(iter.lt.startpdiag) fulldiag=.true.
   if(lastdiag )          fulldiag=.true.
!                                            call timing(t1,w1)
   call solve(fulldiag,ndim,ihomo,scfconv,H,S,X,P,emo,fail)
!                                            call timing(t2,w2)
!                            call prtime(6,t2-t1,w2-w1,'diag')

   if(fail)then
      eel=1.d+99
      return
   endif

   if(ihomo+1.le.ndim.and.ihomo.ge.1)egap=emo(ihomo+1)-emo(ihomo)
!  automatic reset to small value
   if(egap.lt.0.1.and.iter.eq.0) broydamp=0.03

!  Fermi smearing
   if(et.gt.0.1)then
!     convert restricted occ first to alpha/beta
      if(nel.gt.0) then
         call occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
      else
         focca=0.0d0
         foccb=0.0d0
         ihomoa=0
         ihomob=0
      endif
      if (ihomoa+1.le.ndim) then
         call fermismear(.false.,ndim,ihomoa,et,emo,focca,nfoda,efa,ga)
      endif
      if (ihomob+1.le.ndim) then
         call fermismear(.false.,ndim,ihomob,et,emo,foccb,nfodb,efb,gb)
      endif
      focc = focca + foccb
   else
      ga = 0.0_wp
      gb = 0.0_wp
   endif

!  save q
   q_in(1:nshell)=qsh(1:nshell)
   k=nshell
   call gfn2broyden_save(n,k,nbr,dipm,qp,q_in)

!  density matrix
   call dmat(ndim,focc,H,P)

!  new q
   call mpopsh (n,ndim,nshell,ao2sh,S,P,qsh)
   qsh = zsh - qsh

!  qat from qsh
   call qsh2qat(n,at,nshell,qsh,q)

   eold=eel
   call electro2(n,at,ndim,nshell,jab,H0,P,q, &
   &                gam3sh,qsh,gscal,ees,eel)
!  multipole electrostatic
   call mmompop(n,ndim,aoat2,xyz,p,s,dpint,qpint,dipm,qp)
!  call scalecamm(n,at,dipm,qp)
!  DEBUG option: check, whether energy from Fock matrix coincides
!                w/ energy routine
!  include 'cammcheck.inc'
!  evaluate energy
   call aniso_electro(n,at,xyz,q,dipm,qp,gab3,gab5,eaes,epol)
   eel=eel+eaes+epol
! SAW start - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1804
   if (newdisp) then
      ed = edisp_scc(n,dispdim,at,q,g_a,g_c,wdispmat,gw)
      eel = eel + ed
   endif
! SAW end - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1804

!  point charge contribution
   if(pcem) call electro_pcem(nshell,qsh,Vpc,epcem,eel)

!!  new cm5 charges and gborn energy
   if(lgbsa) then
      cm5=q+cm5a
      call electro_gbsa(n,at,fgb,fhb,cm5,gborn,eel)
   endif

!  ad el. entropies*T
   eel=eel+ga+gb

!! ------------------------------------------------------------------------
!  check for energy convergence
   econverged = abs(eel - eold) < scfconv
!! ------------------------------------------------------------------------

   dq(1:nshell)=qsh(1:nshell)-q_in(1:nshell)
   k=nshell
   call gfn2broyden_diff(n,k,nbr,dipm,qp,q_in,dq) ! CAMM case
   rmsq=sum(dq(1:nbr)**2)/dble(n)
   rmsq=sqrt(rmsq)

!! ------------------------------------------------------------------------
!  end of SCC convergence part
   qconverged = rmsq < qconv
!! ------------------------------------------------------------------------

!  SCC convergence acceleration
   if(.not.broy)then

!      simple damp
      if(iter.gt.0) then
         omegap=egap
         ! monopoles only
         do i=1,nshell
            qsh(i)=damp*qsh(i)+(1.0d0-damp)*q_in(i)
         enddo
         ! CAMM
         k=nshell
         do i=1,n
            do j=1,3
               k=k+1
               dipm(j,i)=damp*dipm(j,i)+(1.0d0-damp)*q_in(k)
            enddo
            do j=1,6
               k=k+1
               qp(j,i)=damp*qp(j,i)+(1.0d0-damp)*q_in(k)
            enddo
         enddo
         if(eel-eold.lt.0) then
            damp=damp*1.15
         else
            damp=damp0
         endif
         damp=min(damp,1.0)
         if(egap.lt.1.0)damp=min(damp,0.5)
      endif

   else

!     Broyden mixing
      omegap=0.0d0
      call broyden(nbr,q_in,qlast_in,dq,dqlast, &
     &             iter,thisiter,broydamp,omega,df,u,a)
      qsh(1:nshell)=q_in(1:nshell)
      k=nshell
      call gfn2broyden_out(n,k,nbr,q_in,dipm,qp) ! CAMM case
      if(iter.gt.1) omegap=omega(iter-1)
   endif ! Broyden?

   call qsh2qat(n,at,nshell,qsh,q) !new qat

! SAW start - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1801
   if(newdisp) call disppot(n,dispdim,at,q,g_a,g_c,wdispmat,gw,hdisp)
! SAW end - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 1801

   if(minpr)write(iunit,'(i4,F15.7,E14.6,E11.3,f8.2,2x,f8.1,l3)') &
   &  iter+jter,eel,eel-eold,rmsq,egap,omegap,fulldiag
   qq=q

   if(lgbsa) cm5=q+cm5a
!  set up ES potential
   if(pcem) then
      ves(1:nshell)=Vpc(1:nshell)
   else
      ves=0.0d0
   endif
   call setespot(nshell,qsh,jab,ves)
!  compute potential intermediates
   call setvsdq(n,at,xyz,q,dipm,qp,gab3,gab5,vs,vd,vq)

!  end of SCC convergence part

!! ------------------------------------------------------------------------
   if (econverged.and.qconverged) then
      converged = .true.
      if (lastdiag) exit scc_iterator
      lastdiag = .true.
   endif
!! ------------------------------------------------------------------------

   enddo scc_iterator

   jter = jter + min(iter,thisiter)
   fail = .not.converged

end subroutine scc_gfn2

!! ========================================================================
!  H0 off-diag scaling
!! ========================================================================
subroutine h0scal(n,at,i,j,ishell,jshell,iat,jat,valaoi,valaoj,kspd,kmagic, &
   &              kenscal,km)
   use aoparam,  only : kpair,en
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: i
   integer, intent(in)  :: j
   integer, intent(in)  :: ishell
   integer, intent(in)  :: jshell
   integer, intent(in)  :: iat
   integer, intent(in)  :: jat
   logical, intent(in)  :: valaoi
   logical, intent(in)  :: valaoj
   real(wp),intent(in)  :: kspd(6)
   real(wp),intent(in)  :: kmagic(4,4)
   real(wp),intent(in)  :: kenscal
   real(wp),intent(out) :: km
   integer  :: ii,jj
   real(wp) :: den

   km = 0.0_wp

!  valence
   if(valaoi.and.valaoj) then
      ii=at(iat)
      jj=at(jat)
      den=(en(ii)-en(jj))**2
      km=kmagic(jshell,ishell)*(1.0d0-kenscal*0.01*den)*kpair(ii,jj)
      return
   endif

!  "DZ" functions (on H for GFN or 3S for EA calc on all atoms)
   if((.not.valaoi).and.(.not.valaoj)) then
      km=kspd(6)
      return
   endif
   if(.not.valaoi.and.valaoj) then
      km=0.5*(kspd(jshell)+kspd(6))
      return
   endif
   if(.not.valaoj.and.valaoi) then
      km=0.5*(kspd(ishell)+kspd(6))
   endif


end subroutine h0scal

!! ========================================================================
!  total energy for GFN1
!! ========================================================================
pure subroutine electro(n,at,nbf,nshell,gab,H0,P,dq,dqsh,es,scc)
   use mctc_econv, only : evtoau
   use aoparam, only : gam3
   implicit none
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   integer, intent(in) :: nbf
   integer, intent(in) :: nshell
   real(wp),intent(in)  :: H0(nbf*(nbf+1)/2)
   real(wp),intent(in)  :: P (nbf,nbf)
   real(wp),intent(in)  :: gab(nshell,nshell)
   real(wp),intent(in)  :: dq(n)
   real(wp),intent(in)  :: dqsh(nshell)
   real(wp),intent(out) :: es
   real(wp),intent(out) :: scc
   real(wp) :: ehb ! not used

   integer  :: i,j,k
   real(wp) :: h,t

!  second order non-diagonal
   es =0.0d0
   do i=1,nshell-1
      do j=i+1,nshell
         es =es + dqsh(i)*dqsh(j)*gab(j,i)
      enddo
   enddo

   es=es*2.0_wp

!  second-order diagonal term
   do i=1,nshell
      es =es + dqsh(i)*dqsh(i)*gab(i,i)
   enddo

   t=0.0_wp
   do i=1,n
!     third-order diagonal term
      t = t + gam3(at(i))*dq(i)**3
   enddo

!  ES energy in Eh (gam3 in Eh)
   es=0.50_wp*es*evtoau+t/3.0_wp

!  H0 part
   k=0
   h=0.0_wp
   do i=1,nbf
      do j=1,i-1
         k=k+1
         h=h+P(j,i)*H0(k)
      enddo
      k=k+1
      h=h+P(i,i)*H0(k)*0.5_wp
   enddo

!  Etotal in Eh
   scc = es + 2.0_wp*h*evtoau

end subroutine electro

!! ========================================================================
!  total energy for GFN2
!! ========================================================================
pure subroutine electro2(n,at,nbf,nshell,gab,H0,P,q,  &
   &                     gam3sh,dqsh,gscal,es,scc)
   use mctc_constants, only : pi
   use mctc_econv, only : evtoau
   implicit none
   integer,intent(in)  :: n
   integer,intent(in)  :: at(n)
   integer,intent(in)  :: nbf
   integer,intent(in)  :: nshell
   real(wp), intent(in)  :: H0(nbf*(nbf+1)/2)
   real(wp), intent(in)  :: P (nbf,nbf)
   real(wp), intent(in)  :: q (n) ! not used
   real(wp), intent(in)  :: gscal ! not used
   real(wp), intent(in)  :: gab(nshell,nshell)
   real(wp), intent(in)  :: gam3sh(nshell)
   real(wp), intent(in)  :: dqsh(nshell)
   real(wp), intent(out) :: scc
   real(wp), intent(out) :: es
   real(wp) :: ehb ! not used

   integer :: i,j,k
   real(wp)  :: h,t,esie

!  second order non-diagonal
   es =0.0d0
   do i=1,nshell-1
      do j=i+1,nshell
         es =es + dqsh(i)*dqsh(j)*gab(j,i)
      enddo
   enddo

   es=es*2.0d0

!  second-order diagonal term
   do i=1,nshell
      es =es + dqsh(i)*dqsh(i)*gab(i,i)
   enddo

   t=0.0d0
   do i=1,nshell
!     third-order diagonal term
      t = t + gam3sh(i)*dqsh(i)**3
   enddo

!  SIE diagonal term
!  esie=0
!  do i=1,n
!     esie=esie+gsie(at(i))*(sin(pi*q(i)))**2
!  enddo

!  ES energy in Eh (gam3 in Eh)
   es=0.50d0*es*evtoau+t/3.0d0

!  H0 part
   k=0
   h=0.0d0
   do i=1,nbf
      do j=1,i-1
         k=k+1
         h=h+P(j,i)*H0(k)
      enddo
      k=k+1
      h=h+P(i,i)*H0(k)*0.5d0
   enddo

!  Etotal in Eh
   scc = es + 2.0d0*h*evtoau

end subroutine electro2


!! ========================================================================
!  GBSA related subroutine
!! ========================================================================
pure subroutine electro_gbsa(n,at,gab,fhb,dqsh,es,scc)
   use mctc_econv, only : evtoau
   use gbobc, only: lhb
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: gab(n,n)
   real(wp),intent(in)  :: fhb(n)
   real(wp),intent(in)  :: dqsh(n)
   real(wp),intent(out) :: es
   real(wp),intent(inout) :: scc

   integer :: i,j,k
   real(wp)  :: h,t
   real(wp)  :: ehb

!  second order non-diagonal
   es =0
   do i=1,n-1
      do j=i+1,n
         es =es + dqsh(i)*dqsh(j)*gab(j,i)
      enddo
   enddo

   es=es*2.0d0

!  second-order diagonal term + HB contribution
   do i=1,n
      es =es + dqsh(i)*dqsh(i)*gab(i,i)
   enddo

!  HB energy
   ehb=0.d0
   if(lhb)then
      do i = 1, n
!        ehb = ehb + fhb(i)*(dqsh(i)**2)**c3
         ehb = ehb + fhb(i)*(dqsh(i)**2)
      enddo
   endif

!  ES energy in Eh
   es=0.5*es*evtoau

!  HB energy in Eh
   ehb=ehb*evtoau

!  Etotal in Eh
   scc = scc + es + ehb

end subroutine electro_gbsa

!! ========================================================================
!  S(R) enhancement factor
!! ========================================================================
   pure function rfactor(ish,jsh,ati,atj,xyz1,xyz2)
   use aoparam, only : rad,polyr
   use mctc_econv, only : aatoau
   implicit none
   integer,intent(in) :: ati,atj,ish,jsh
   real(wp), intent(in) :: xyz1(3),xyz2(3)
   real(wp) :: rfactor
   real(wp) :: rab,k1,rr,r,rf1,rf2,dx,dy,dz,a

   a=0.5           ! R^a dependence 0.5 in GFN1

   dx=xyz1(1)-xyz2(1)
   dy=xyz1(2)-xyz2(2)
   dz=xyz1(3)-xyz2(3)

   rab=sqrt(dx**2+dy**2+dz**2)

   ! this sloppy conv. factor has been used in development, keep it
   rr=(rad(ati)+rad(atj))*aatoau

   r=rab/rr

   k1=polyr(ish,ati)
   rf1=1.0d0+0.01*k1*r**a
   k1=polyr(jsh,atj)
   rf2=1.0d0+0.01*k1*r**a

   rfactor= rf1*rf2

end function rfactor

!! ========================================================================
!  set up Coulomb potential due to 2nd order fluctuation
!! ========================================================================
pure subroutine setespot(nshell,qsh,jab,ves)
   implicit none
   integer, intent(in) :: nshell
   real(wp),intent(in) ::  qsh(nshell),jab(nshell,nshell)
!  ves possibly already contains with PC-potential
   real(wp),intent(inout) :: ves(nshell)
   real(wp) :: qshi,vesi
   integer  :: i,j,k
   do i=1,nshell
      qshi=qsh(i)
      vesi=0.0_wp
      do j=1,i-1
         ves(j)=ves(j)+qshi*jab(j,i)
         vesi=vesi+qsh(j)*jab(j,i)
      enddo
      vesi=vesi+qshi*jab(i,i)
      ves(i)=ves(i)+vesi
   enddo
end subroutine setespot

pure subroutine jpot_gfn1(nat,nshell,ash,lsh,at,sqrab,alphaj,jab)
   use mctc_econv
   use aoparam
   use lin_mod
   implicit none
   integer, intent(in) :: nat
   integer, intent(in) :: nshell
   integer, intent(in) :: ash(nshell)
   integer, intent(in) :: lsh(nshell)
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: sqrab(nat*(nat+1)/2)
   real(wp),intent(in) :: alphaj
   real(wp),intent(inout) :: jab(nshell,nshell)

   integer  :: is,iat,ati,js,jat,atj,k
   real(wp) :: gi,gj,xj,rab

   do is=1,nshell
      iat=ash(is)
      ati=at(iat)
      gi=gam(ati)*(1.0_wp+lpar(lsh(is),ati))
      do js=1,is
         jat=ash(js)
         atj=at(jat)
         k=lin(jat,iat)
         gj=gam(atj)*(1.0_wp+lpar(lsh(js),atj))
         xj=2.0_wp/(1./gi+1./gj)
         if(is.eq.js)then
            jab(is,js)=xj*autoev
         else
            rab=sqrt(sqrab(k))
            jab(js,is)=autoev/(rab**alphaj &
               &  + 1._wp/xj**alphaj)**(1._wp/alphaj)
            jab(is,js)=jab(js,is)
         endif
      enddo
   enddo

end subroutine jpot_gfn1

pure subroutine jpot_gfn2(nat,nshell,ash,lsh,at,sqrab,jab)
   use mctc_econv
   use aoparam
   use lin_mod
   implicit none
   integer, intent(in) :: nat
   integer, intent(in) :: nshell
   integer, intent(in) :: ash(nshell)
   integer, intent(in) :: lsh(nshell)
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: sqrab(nat*(nat+1)/2)
   real(wp),intent(inout) :: jab(nshell,nshell)

   integer  :: is,iat,ati,js,jat,atj,k
   real(wp) :: gi,gj,xj,rab

   do is=1,nshell
      iat=ash(is)
      ati=at(iat)
      gi=gam(ati)*(1.0_wp+lpar(lsh(is),ati))
      do js=1,is-1
         jat=ash(js)
         atj=at(jat)
         k=lin(jat,iat)
         gj=gam(atj)*(1.0_wp+lpar(lsh(js),atj))
         xj=0.5_wp*(gi+gj)
         jab(js,is)=autoev/sqrt(sqrab(k)+1._wp/xj**2)
         ! jab(js,is)=autoev/sqrt(sqrab(k)+1._wp/(gi*gj))  ! NEWAV
         jab(is,js)=jab(js,is)
      enddo
      jab(is,is)=autoev*gi
   enddo

end subroutine jpot_gfn2

!! ========================================================================
!  eigenvalue solver single-precision
!! ========================================================================
subroutine solve4(full,ndim,ihomo,acc,H,S,X,P,e,fail)
   use iso_fortran_env, sp => real32
   implicit none
   integer, intent(in)   :: ndim
   logical, intent(in)   :: full
   integer, intent(in)   :: ihomo
   real(wp),intent(inout):: H(ndim,ndim)
   real(wp),intent(in)   :: S(ndim,ndim)
   real(wp),intent(out)  :: X(ndim,ndim)
   real(wp),intent(out)  :: P(ndim,ndim)
   real(wp),intent(out)  :: e(ndim)
   real(wp),intent(in)   :: acc
   logical, intent(out)  :: fail

   integer i,j,info,lwork,liwork,nfound,iu,nbf
   integer, allocatable :: iwork(:),ifail(:)
   real(wp),allocatable :: aux  (:)
   real(wp) w0,w1,t0,t1

   real(sp),allocatable :: H4(:,:)
   real(sp),allocatable :: S4(:,:)
   real(sp),allocatable :: X4(:,:)
   real(sp),allocatable :: P4(:,:)
   real(sp),allocatable :: e4(:)
   real(sp),allocatable :: aux4(:)


   allocate(H4(ndim,ndim),S4(ndim,ndim))
   allocate(X4(ndim,ndim),P4(ndim,ndim),e4(ndim))
 
   H4 = H
   S4 = S

   fail =.false.
!  standard first full diag call
   if(full) then
!                                                     call timing(t0,w0)
!     if(ndim.gt.0)then
!     USE DIAG IN NON-ORTHORGONAL BASIS
      allocate (aux4(1),iwork(1),ifail(ndim))
      P4 = s4
      call sygvd(1,'v','u',ndim,h4,ndim,p4,ndim,e4,aux4, &!workspace query
     &           -1,iwork,liwork,info)
      lwork=int(aux4(1))
      liwork=iwork(1)
      deallocate(aux4,iwork)
      allocate (aux4(lwork),iwork(liwork))              !do it
      call sygvd(1,'v','u',ndim,h4,ndim,p4,ndim,e4,aux4, &
     &           lwork,iwork,liwork,info)
      if(info.ne.0) then
         fail=.true.
         return
      endif
      X4 = H4 ! save
      deallocate(aux4,iwork,ifail)

!     else
!        USE DIAG IN ORTHOGONAL BASIS WITH X=S^-1/2 TRAFO
!        nbf = ndim
!        lwork  = 1 + 6*nbf + 2*nbf**2
!        allocate (aux(lwork))
!        call gemm('N','N',nbf,nbf,nbf,1.0d0,H,nbf,X,nbf,0.0d0,P,nbf)
!        call gemm('T','N',nbf,nbf,nbf,1.0d0,X,nbf,P,nbf,0.0d0,H,nbf)
!        call SYEV('V','U',nbf,H,nbf,e,aux,lwork,info)
!        if(info.ne.0) error stop 'diag error'
!        call gemm('N','N',nbf,nbf,nbf,1.0d0,X,nbf,H,nbf,0.0d0,P,nbf)
!        H = P
!        deallocate(aux)
!     endif
!                                                     call timing(t1,w1)
!                                    call prtime(6,t1-t0,w1-w0,'dsygvd')

   else
!                                                     call timing(t0,w0)
!     go to MO basis using trafo(X) from first iteration (=full diag)
!      call gemm('N','N',ndim,ndim,ndim,1.d0,H4,ndim,X4,ndim,0.d0,P4,ndim)
!      call gemm('T','N',ndim,ndim,ndim,1.d0,X4,ndim,P4,ndim,0.d0,H4,ndim)
!                                                     call timing(t1,w1)
!                       call prtime(6,1.5*(t1-t0),1.5*(w1-w0),'3xdgemm')
!                                                     call timing(t0,w0)
!      call pseudodiag(ndim,ihomo,H4,e4)
!                                                     call timing(t1,w1)
!                                call prtime(6,t1-t0,w1-w0,'pseudodiag')

!     C = X C', P=scratch
!      call gemm('N','N',ndim,ndim,ndim,1.d0,X4,ndim,H4,ndim,0.d0,P4,ndim)
!     save and output MO matrix in AO basis
!      H4 = P4
   endif

   H = H4
   P = P4
   X = X4
   e = e4

   deallocate(e4,P4,X4,S4,H4)

end subroutine solve4

!! ========================================================================
!  eigenvalue solver
!! ========================================================================
subroutine solve(full,ndim,ihomo,acc,H,S,X,P,e,fail)
   implicit none
   integer, intent(in)   :: ndim
   logical, intent(in)   :: full
   integer, intent(in)   :: ihomo
   real(wp),intent(inout):: H(ndim,ndim)
   real(wp),intent(in)   :: S(ndim,ndim)
   real(wp),intent(out)  :: X(ndim,ndim)
   real(wp),intent(out)  :: P(ndim,ndim)
   real(wp),intent(out)  :: e(ndim)
   real(wp),intent(in)   :: acc
   logical, intent(out)  :: fail

   integer i,j,info,lwork,liwork,nfound,iu,nbf
   integer, allocatable :: iwork(:),ifail(:)
   real(wp),allocatable :: aux  (:)
   real(wp) w0,w1,t0,t1

   fail =.false.

!  standard first full diag call
   if(full) then
!                                                     call timing(t0,w0)
!     if(ndim.gt.0)then
!     USE DIAG IN NON-ORTHORGONAL BASIS
      allocate (aux(1),iwork(1),ifail(ndim))
      P = s
      call sygvd(1,'v','u',ndim,h,ndim,p,ndim,e,aux, &!workspace query
     &           -1,iwork,liwork,info)
      lwork=int(aux(1))
      liwork=iwork(1)
      deallocate(aux,iwork)
      allocate (aux(lwork),iwork(liwork))              !do it
      call sygvd(1,'v','u',ndim,h,ndim,p,ndim,e,aux, &
     &           lwork,iwork,liwork,info)
      !write(*,*)'SYGVD INFO', info
      if(info.ne.0) then
         fail=.true.
         return
      endif
      X = H ! save
      deallocate(aux,iwork,ifail)

!     else
!        USE DIAG IN ORTHOGONAL BASIS WITH X=S^-1/2 TRAFO
!        nbf = ndim
!        lwork  = 1 + 6*nbf + 2*nbf**2
!        allocate (aux(lwork))
!        call gemm('N','N',nbf,nbf,nbf,1.0d0,H,nbf,X,nbf,0.0d0,P,nbf)
!        call gemm('T','N',nbf,nbf,nbf,1.0d0,X,nbf,P,nbf,0.0d0,H,nbf)
!        call SYEV('V','U',nbf,H,nbf,e,aux,lwork,info)
!        if(info.ne.0) error stop 'diag error'
!        call gemm('N','N',nbf,nbf,nbf,1.0d0,X,nbf,H,nbf,0.0d0,P,nbf)
!        H = P
!        deallocate(aux)
!     endif
!                                                     call timing(t1,w1)
!                                    call prtime(6,t1-t0,w1-w0,'dsygvd')

   else
!                                                     call timing(t0,w0)
!     go to MO basis using trafo(X) from first iteration (=full diag)
      call gemm('N','N',ndim,ndim,ndim,1.d0,H,ndim,X,ndim,0.d0,P,ndim)
      call gemm('T','N',ndim,ndim,ndim,1.d0,X,ndim,P,ndim,0.d0,H,ndim)
!                                                     call timing(t1,w1)
!                       call prtime(6,1.5*(t1-t0),1.5*(w1-w0),'3xdgemm')
!                                                     call timing(t0,w0)
      call pseudodiag(ndim,ihomo,H,e)
!                                                     call timing(t1,w1)
!                                call prtime(6,t1-t0,w1-w0,'pseudodiag')

!     C = X C', P=scratch
      call gemm('N','N',ndim,ndim,ndim,1.d0,X,ndim,H,ndim,0.d0,P,ndim)
!     save and output MO matrix in AO basis
      H = P
   endif

end subroutine solve

subroutine fermismear(prt,norbs,nel,t,eig,occ,fod,e_fermi,s)
   use mctc_econv, only : autoev
   use mctc_constants, only : kB
   implicit none
   integer, intent(in)  :: norbs
   integer, intent(in)  :: nel
   real(wp),intent(in)  :: eig(norbs)
   real(wp),intent(out) :: occ(norbs)
   real(wp),intent(in)  :: t
   real(wp),intent(out) :: fod
   real(wp),intent(out) :: e_fermi
   logical, intent(in)  :: prt

   real(wp) :: boltz,bkt,occt,total_number,thr
   real(wp) :: total_dfermi,dfermifunct,fermifunct,s,change_fermi

   parameter (boltz = kB*autoev)
   parameter (thr   = 1.d-9)
   integer :: ncycle,i,j,m,k,i1,i2

   bkt = boltz*t

   e_fermi = 0.5*(eig(nel)+eig(nel+1))
   occt=nel

   do ncycle = 1, 200  ! this loop would be possible instead of gotos
      total_number = 0.0
      total_dfermi = 0.0
      do i = 1, norbs
         fermifunct = 0.0
         if((eig(i)-e_fermi)/bkt.lt.50) then
            fermifunct = 1.0/(exp((eig(i)-e_fermi)/bkt)+1.0)
            dfermifunct = exp((eig(i)-e_fermi)/bkt) / &
            &       (bkt*(exp((eig(i)-e_fermi)/bkt)+1.0)**2)
         else
            dfermifunct = 0.0
         end if
         occ(i) = fermifunct
         total_number = total_number + fermifunct
         total_dfermi = total_dfermi + dfermifunct
      end do
      change_fermi = (occt-total_number)/total_dfermi
      e_fermi = e_fermi+change_fermi
      if (abs(occt-total_number).le.thr) exit
   enddo

   fod=0
   s  =0
   do i=1,norbs
      if(occ(i).gt.thr.and.1.0d00-occ(i).gt.thr) &
      &   s=s+occ(i)*log(occ(i))+(1.0d0-occ(i))*log(1.0d00-occ(i))
      if (eig(i).lt.e_fermi) then
         fod=fod+1.0d0-occ(i)
      else
         fod=fod+      occ(i)
      endif
   enddo
   s=s*kB*t

   if (prt) then
      write(*,'('' t,e(fermi),nfod : '',2f10.3,f10.6)') t,e_fermi,fod
   endif

end subroutine fermismear

subroutine occ(ndim,nel,nopen,ihomo,focc)
   implicit none
   integer  :: nel
   integer  :: nopen
   integer  :: ndim
   integer  :: ihomo
   real(wp) :: focc(ndim)
   integer  :: i,na,nb

   focc=0
!  even nel      
   if(mod(nel,2).eq.0)then
      ihomo=nel/2
      do i=1,ihomo 
         focc(i)=2.0d0
      enddo
      if(2*ihomo.ne.nel) then
         ihomo=ihomo+1
         focc(ihomo)=1.0d0
         if(nopen.eq.0)nopen=1
      endif
      if(nopen.gt.1)then
         do i=1,nopen/2
            focc(ihomo-i+1)=focc(ihomo-i+1)-1.0
            focc(ihomo+i)=focc(ihomo+i)+1.0
         enddo
      endif
!  odd nel      
   else
      na=nel/2+(nopen-1)/2+1
      nb=nel/2-(nopen-1)/2
      do i=1,na             
         focc(i)=focc(i)+1.
      enddo
      do i=1,nb             
         focc(i)=focc(i)+1.
      enddo
   endif

   do i=1,ndim
      if(focc(i).gt.0.99) ihomo=i
   enddo

end subroutine occ

subroutine occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
   implicit none
   integer  :: nel
   integer  :: nopen
   integer  :: ndim
   integer  :: ihomoa
   integer  :: ihomob
   real(wp) :: focca(ndim)
   real(wp) :: foccb(ndim)
   integer  :: focc(ndim)
   integer  :: i,na,nb,ihomo

   focc=0
   focca=0
   foccb=0
!  even nel      
   if(mod(nel,2).eq.0)then
      ihomo=nel/2
      do i=1,ihomo 
         focc(i)=2
      enddo
      if(2*ihomo.ne.nel) then
         ihomo=ihomo+1
         focc(ihomo)=1
         if(nopen.eq.0)nopen=1
      endif
      if(nopen.gt.1)then
         do i=1,nopen/2
            focc(ihomo-i+1)=focc(ihomo-i+1)-1
            focc(ihomo+i)=focc(ihomo+i)+1
         enddo
      endif
!  odd nel      
   else
      na=nel/2+(nopen-1)/2+1
      nb=nel/2-(nopen-1)/2
      do i=1,na             
         focc(i)=focc(i)+1
      enddo
      do i=1,nb             
         focc(i)=focc(i)+1
      enddo
   endif

   do i=1,ndim
      if(focc(i).eq.2)then
         focca(i)=1.0d0
         foccb(i)=1.0d0
      endif
      if(focc(i).eq.1)focca(i)=1.0d0
   enddo

   ihomoa=0
   ihomob=0
   do i=1,ndim
      if(focca(i).gt.0.99) ihomoa=i
      if(foccb(i).gt.0.99) ihomob=i
   enddo

   if(ihomoa.lt.1) call raise('E','internal error in occu',1)
end subroutine occu


subroutine epart
   use aoparam
   use mctc_econv, only : evtoau
   implicit none
   real(wp) :: h,e,p

   p=0.5*gam(1)+gam3(1)/3.0_wp

   h=ao_lev(1,1)*evtoau

   e=h-p

   write(output_unit,'(a)')
   write(output_unit,'(''H atom         energy :'',F10.6)') h
   write(output_unit,'(''proton   (self)energy :'',F10.6)') p
   write(output_unit,'(''electron (self)energy :'',F10.6)') e
   write(output_unit,'(''these values must be considered in IP/EA/PA calc.!'')')

end subroutine epart

subroutine eself(n,at,z)
   use aoparam
   use mctc_econv, only : evtoau
   implicit none
   integer, intent(in) :: n,at(n)
   real(wp),intent(in) :: z(n)

   real(wp) :: e
   integer  :: i,ii

   e=0.0_wp
   do i=1,n
      ii=at(i)
      e=e+0.5_wp*z(i)**2*gam(ii)+z(i)**3*gam3(ii)/3.0_wp
   enddo

   write(output_unit,'(''NO ELECTRONS!'')')
   write(output_unit,'(''molecular/atomic self energy :'',F10.6)') e

end subroutine eself

!ccccccccccccccccccccccccccccccccccccccccccccc
! density matrix
! C: MO coefficient
! X: scratch
! P  dmat
!ccccccccccccccccccccccccccccccccccccccccccccc

subroutine dmat(ndim,focc,C,P)
   use iso_fortran_env, only : wp => real64
   use mctc_la, only : gemm
   implicit none
   integer, intent(in)  :: ndim
   real(wp),intent(in)  :: focc(*)
   real(wp),intent(in)  :: C(ndim,ndim)
   real(wp),intent(out) :: P(ndim,ndim)
   integer :: i,m
   real(wp),allocatable :: Ptmp(:,:)

   allocate( Ptmp(ndim,ndim), source = 0.0_wp )

   do m=1,ndim
      do i=1,ndim
         Ptmp(i,m)=C(i,m)*focc(m)
      enddo
   enddo
   call gemm('n','t',ndim,ndim,ndim,1.0_wp,C,ndim,Ptmp,ndim,0.0_wp,P,ndim)

   deallocate(Ptmp)

end subroutine dmat

!ccccccccccccccccccccccccccccccccccccccccccccc
!            Wiberg BOs
!ccccccccccccccccccccccccccccccccccccccccccccc

subroutine wiberg(n,ndim,at,xyz,P,S,wb,ex,pr,fila2)
   use iso_fortran_env, only : wp => real64
   use mctc_la, only : gemm
   implicit none
   integer, intent(in)  :: n,ndim,at(n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: P(ndim,ndim)
   real(wp),intent(in)  :: S(ndim,ndim)
   real(wp),intent(out) :: wb (n,n)
   logical, intent(in)  :: ex,pr
   integer, intent(in)  :: fila2(:,:)

   real(wp),allocatable ::Ptmp(:,:)
   real(wp) xsum,rab
!  real(wp) wbr(n,n)
!  real(wp) ddum(3,n,n),cdum(n)
   real(wp),allocatable :: wbr(:,:)
   real(wp),allocatable :: ddum(:,:,:),cdum(:)
!  integer i,j,k,m,ibmax,imem(n)
   integer i,j,k,m,ibmax
   integer, allocatable :: imem(:)
   character(2) asym
   integer :: ich

   allocate(Ptmp(ndim,ndim))
   call GEMM('N','N',ndim,ndim,ndim,1.0d0,P,ndim,S,ndim,0.0d0,Ptmp,ndim)
   call open_file(ich,'wbo','w')
   wb=0
   do i=1,n
      do j=1,i-1
         xsum=0.0
         rab=(xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2
         if(rab.lt.100.0)then
            do k=fila2(1,i),fila2(2,i)   ! AOs on atom i
               do m=fila2(1,j),fila2(2,j) ! AOs on atom j
                  xsum=xsum+Ptmp(k,m)*Ptmp(m,k)
               enddo
            enddo
         endif
         wb(i,j)=xsum
         wb(j,i)=xsum
         if(abs(xsum).gt.0.3) write(ich,*)i,j,xsum
      enddo
   enddo
   call close_file(ich)
   deallocate(Ptmp)

   if (pr) then
      allocate( wbr(n,n), source = 0.0_wp )
      allocate( imem(n),  source = 0 )
      wbr = wb

      write(*,*)
      write(*,*)'Wiberg/Mayer (AO) data. WBOs > 0.1 written file <wbo>'
      write(*,*)'largest (>0.05) Wiberg bond orders for each atom'
      write(*,*)'          total WBO             WBO to atom ...'
      do i=1,n
         do j=1,n
            imem(j)=j
         enddo
         call wibsort(n,i,imem,wbr)
         ibmax=0
         xsum =0.0_wp
         do j=1,n
            if(wbr(i,j).gt.0.05_wp)ibmax=j
            xsum=xsum+wbr(i,j)
         enddo
         write(*,'(i6,a4,1x,f6.3,9(4x,a2,i4,f6.3))') &
            & i,asym(at(i)),xsum, &
            & (asym(at(imem(j))),imem(j),wbr(i,j),j=1,ibmax) 
         !        if(fit) write(112,'(2F22.12)') xsum,cdum(i)
      enddo

      deallocate(wbr,imem)
   endif

end subroutine wiberg

subroutine wiberg_nosort(n,ndim,at,xyz,P,S,wb,ex,fila2)
   use iso_fortran_env, only : wp => real64
   use mctc_la, only : gemm
   implicit none
   integer n,ndim,at(n)
   real(wp) xyz(3,n)
   real(wp) P(ndim,ndim)
   real(wp) S(ndim,ndim)
   real(wp) wb (n,n)
   logical ex
   integer, intent(in) :: fila2(:,:)

   real(wp),allocatable ::Ptmp(:,:)
   real(wp) xsum,rab
   real(wp) wbr(n,n)
   integer i,j,k,m,ibmax,imem(n)
   character(2) asym
   integer :: ich

   allocate(Ptmp(ndim,ndim))
   call GEMM('N','N',ndim,ndim,ndim,1.0d0,P,ndim,S,ndim,0.0d0,Ptmp,ndim)
   call open_file(ich,'wbo','w')
   wb=0
   do i=1,n
      do j=1,i-1
         xsum=0.0
         rab=(xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2
         if(rab.lt.100.0)then
            do k=fila2(1,i),fila2(2,i)
               do m=fila2(1,j),fila2(2,j)
                  xsum=xsum+Ptmp(k,m)*Ptmp(m,k)
               enddo
            enddo
         endif
         wb(i,j)=xsum
         wb(j,i)=xsum
         if(abs(xsum).gt.0.3) write(ich,*)i,j,xsum
      enddo
   enddo
   call close_file(ich)
   deallocate(Ptmp)

end subroutine wiberg_nosort

!cccccccccccccccccccccccccccccccccccccccc

SUBROUTINE wibsort(ncent,imo,imem,qmo)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer  :: ncent
   integer  :: imo
   real(wp) :: qmo(ncent,ncent)
   integer  :: imem(ncent)
   integer  :: ii,i,j,k,ihilf
   real(wp) :: pp

   do ii = 2,ncent
      i = ii - 1
      k = i
      pp= qmo(imo,i)
      do j = ii, ncent
         if (qmo(imo,j) .lt. pp) cycle
         k = j
         pp=qmo(imo,j)
      enddo
      if (k .eq. i) cycle
      qmo(imo,k) = qmo(imo,i)
      qmo(imo,i) = pp

      ihilf=imem(i)
      imem(i)=imem(k)
      imem(k)=ihilf
   enddo

end SUBROUTINE wibsort

!cccccccccccccccccccccccccccccccccccccccc
!c Mulliken pop + AO pop
!cccccccccccccccccccccccccccccccccccccccc

subroutine mpopall(n,nao,aoat,S,P,qao,q)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer nao,n,aoat(nao)
   real(wp)  S (nao,nao)
   real(wp)  P (nao,nao)
   real(wp)  qao(nao),q(n),ps

   integer i,j,ii,jj,ij,is,js

   q  = 0
   qao= 0
   do i=1,nao
      ii=aoat(i)
      do j=1,i-1
         jj=aoat(j)
         ps=p(j,i)*s(j,i)
         q(ii)=q(ii)+ps
         q(jj)=q(jj)+ps
         qao(i)=qao(i)+ps
         qao(j)=qao(j)+ps
      enddo
      ps=p(i,i)*s(i,i)
      q(ii)=q(ii)+ps
      qao(i)=qao(i)+ps
   enddo

end subroutine mpopall

!cccccccccccccccccccccccccccccccccccccccc
!c Mulliken pop
!cccccccccccccccccccccccccccccccccccccccc

subroutine mpop0(n,nao,aoat,S,P,q)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer nao,n,aoat(nao)
   real(wp)  S (nao,nao)
   real(wp)  P (nao,nao)
   real(wp)  q(n),ps

   integer i,j,ii,jj,ij,is,js

   q = 0
   do i=1,nao
      ii=aoat(i)
      do j=1,i-1
         jj=aoat(j)
         ps=p(j,i)*s(j,i)
         q(ii)=q(ii)+ps
         q(jj)=q(jj)+ps
      enddo
      ps=p(i,i)*s(i,i)
      q(ii)=q(ii)+ps
   enddo

end subroutine mpop0

!cccccccccccccccccccccccccccccccccccccccc
!c Mulliken AO pop
!cccccccccccccccccccccccccccccccccccccccc

subroutine mpopao(n,nao,S,P,qao)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer nao,n
   real(wp)  S (nao,nao)
   real(wp)  P (nao,nao)
   real(wp)  qao(nao),ps

   integer i,j

   qao = 0
   do i=1,nao
      do j=1,i-1
         ps=p(j,i)*s(j,i)
         qao(i)=qao(i)+ps
         qao(j)=qao(j)+ps
      enddo
      ps=p(i,i)*s(i,i)
      qao(i)=qao(i)+ps
   enddo

end subroutine mpopao

!cccccccccccccccccccccccccccccccccccccccc
!c Mulliken pop
!cccccccccccccccccccccccccccccccccccccccc

subroutine mpop(n,nao,aoat,lao,S,P,q,ql)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer nao,n,aoat(nao),lao(nao)
   real(wp)  S (nao,nao)
   real(wp)  P (nao,nao)
   real(wp)  q(n),ps
   real(wp)  ql(3,n)

   integer i,j,ii,jj,ij,is,js,mmm(20)
   data    mmm/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/

   ql= 0
   q = 0
   do i=1,nao
      ii=aoat(i)
      is=mmm(lao(i))
      do j=1,i-1
         jj=aoat(j)
         js=mmm(lao(j))
         ps=p(j,i)*s(j,i)
         q(ii)=q(ii)+ps
         q(jj)=q(jj)+ps
         ql(is,ii)=ql(is,ii)+ps
         ql(js,jj)=ql(js,jj)+ps
      enddo
      ps=p(i,i)*s(i,i)
      q(ii)=q(ii)+ps
      ql(is,ii)=ql(is,ii)+ps
   enddo

end subroutine mpop

!cccccccccccccccccccccccccccccccccccccccc
!c Mulliken pop shell wise
!cccccccccccccccccccccccccccccccccccccccc

subroutine mpopsh(n,nao,nshell,ao2sh,S,P,qsh)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer nao,n,nshell,ao2sh(nao)
   real(wp)  S (nao,nao)
   real(wp)  P (nao,nao)
   real(wp)  qsh(nshell),ps

   integer i,j,ii,jj,ij

   qsh=0
   do i=1,nao
      ii =ao2sh(i)
      do j=1,i-1
         jj =ao2sh(j)
         ps=p(j,i)*s(j,i)
         qsh(ii)=qsh(ii)+ps
         qsh(jj)=qsh(jj)+ps
      enddo
      ps=p(i,i)*s(i,i)
      qsh(ii)=qsh(ii)+ps
   enddo

end subroutine mpopsh

subroutine qsh2qat(n,at,nshell,qsh,q)
   use iso_fortran_env, only : wp => real64
   use aoparam
   implicit none
   integer,intent(in) :: n,nshell,at(n)
   real(wp),intent(in) :: qsh(nshell)
   real(wp),intent(out) :: q(n)

   integer i,mi,k

   k=0
   do i=1,n
      q(i)=0
      do mi=1,ao_n(at(i))
         k=k+1
         q(i)=q(i)+qsh(k)
      enddo
   enddo

end subroutine qsh2qat


!cccccccccccccccccccccccccccccccccccccccc
!c Loewdin pop
!cccccccccccccccccccccccccccccccccccccccc

subroutine lpop(n,nao,aoat,lao,occ,C,f,q,ql)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer nao,n,aoat(nao),lao(nao)
   real(wp)  C (nao,nao)
   real(wp)  occ(nao)
   real(wp)  q(n)
   real(wp)  ql(3,n)
   real(wp)  f

   integer i,j,ii,jj,js,mmm(20)
   data    mmm/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/
   real(wp)  cc

   do i=1,nao
      if(occ(i).lt.1.d-8) cycle
      do j=1,nao
         cc=f*C(j,i)*C(j,i)*occ(i)
         jj=aoat(j)
         js=mmm(lao(j))
         q(jj)=q(jj)+cc
         ql(js,jj)=ql(js,jj)+cc
      enddo
   enddo

end subroutine lpop

!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c atomic valence shell pops and total atomic energy
!cccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine iniqshell(n,at,z,nshell,q,qsh,gfn_method)
   use iso_fortran_env, only : wp => real64
   use aoparam
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: nshell
   integer, intent(in)  :: gfn_method
   real(wp),intent(in)  :: z(n)
   real(wp),intent(in)  :: q(n)
   real(wp),intent(out) :: qsh(nshell)
   real(wp) :: zshell
   real(wp) :: ntot,fracz
   real(wp) :: iox(86,0:2,2) ! GFN1 on 1
   integer  :: i,j,k,m,l,ll(0:3),iat,lll,iver
   data ll /1,3,5,7/

   !     H           Initial     Orbital Occupancies                     He
   !     Li Be                                            B  C  N  O  F  Ne
   !     Na Mg                                            Al Si P  S  Cl Ar
   !     K  Ca Sc            Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
   !     Rb Sr Y             Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
   !     Cs Ba La Ce-Lu      Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
   !                                      spd shell

   ! GFN1
   data iox / &
      &1,                                                          2, & !He
      &1,2,                                         2, 2, 2, 2, 2, 2, & !Ne
      &1,2,                                         2, 2, 2, 2, 2, 2, & !Ar
      &1,2,2,             2, 2, 2, 2, 2, 2, 2, 2, 2,2, 2, 2, 2, 2, 2, & !Kr
      &1,2,2,             2, 2, 2, 2, 2, 2, 2, 2, 2,2, 2, 2, 2, 2, 2, & !Xe
      &1,2,2, 14*2       ,2, 2, 2, 2, 2, 2, 2, 2, 2,2, 2, 2, 2, 2, 2, & !Rn
      !
      &0,                                                          0, & !He
      &0,0,                                         1, 2, 3, 4, 5, 6, & !Ne
      &0,0,                                         1, 2, 3, 4, 5, 6, & !Ar
      &0,0,0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, & !Kr
      &0,0,0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, & !Xe
      &0,0,0,  14*0,   0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, & !Rn
      !
      & 0,                                                         0, & !He
      & 0,0,                                        0, 0, 0, 0, 0, 0, & !Ne
      & 0,0,                                        0, 0, 0, 0, 0, 0, & !Ar
      & 0,0,1,          2, 3, 4, 5, 6, 7, 8,  9,10, 0, 0, 0, 0, 0, 0, & !Kr
      & 0,0,1,          2, 3, 4, 5, 6, 7, 8,  9,10, 0, 0, 0, 0, 0, 0, & !Xe
      & 0,0,1, 14*1,    2, 3, 4, 5, 6, 7, 8,  9,10, 0, 0, 0, 0, 0, 0, & !Rn
      ! GFN2
      &1,                                                            2, & !He
      &1,2,                                        2  ,1  ,1.5,  2,2,2, & !Ne
      &1,2,                                        2  ,1.5,1.5,  2,2,2, & !Ar
      &1,2,1,             1, 1, 1, 1, 1, 1, 1,1,2, 2,  1.5,1.5,  2,2,2, & !Kr
      &1,2,1,             1, 1, 1, 1, 1, 1, 1,1,2, 2,  2,  2 ,   2,2,2, & !Xe
      &1,2,1, 14*1       ,1, 1, 1, 1, 1, 1, 1,1,2, 2,  2,  2 ,   2,2,2, & !Rn
      !
      &0,                                                            0, & !He
      &0,0,                                        1  , 3  ,3.5, 4,5,6, & !Ne
      &0,0,                                        1  , 2.5,3.5, 4,5,6, & !Ar
      &0,0,1,            1, 1, 1, 1, 1, 1, 1, 0,0, 1,   2.5,3.5, 4,5,6, & !Kr
      &0,0,1,            1, 1, 1, 1, 1, 1, 1, 0,0, 1,   2,  3,   4,5,6, & !Xe
      &0,0,1,  14*1,     1, 1, 1, 1, 1, 1, 1, 0,0, 1,   2,  3,   4,5,6, & !Rn
      !
      & 0,                                                          0, & !He
      & 0,0,                                        0, 0 , 0, 0, 0, 0, & !Ne
      & 0,0,                                        0, 0 , 0, 0, 0, 0, & !Ar
      & 0,0,1,          2, 3, 4, 5, 6, 7, 8, 10,10, 0, 0 , 0, 0, 0, 0, & !Kr
      & 0,0,1,          2, 3, 4, 5, 6, 7, 8, 10,10, 0, 0 , 0, 0, 0, 0, & !Xe
      & 0,0,1, 14*1,    2, 3, 4, 5, 6, 7, 8, 10,10, 0, 0 , 0, 0, 0, 0  & !Rn
      & /

   qsh = 0.0_wp

   iver=1
   if(gfn_method.gt.1) iver = 2

   k=0
   do i=1,n
      iat=at(i)
      ntot=-1.d-6
      do m=1,ao_n(iat)
         l=ao_l(m,iat)
         k=k+1
         zshell=iox(iat,l,iver)
         ntot=ntot+zshell
         if(ntot.gt.z(i)) zshell=0
         fracz=zshell/z(i)
         qsh(k)=fracz*q(i)
      enddo
   enddo
   if(k.ne.nshell) error stop 'internal setzshell error 1'

end subroutine iniqshell


subroutine setzshell(n,at,nshell,z,zsh,e,gfn_method)
   use iso_fortran_env, only : wp => real64
   use aoparam
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: nshell
   integer, intent(in)  :: gfn_method
   real(wp),intent(in)  :: z(n)
   real(wp),intent(out) :: zsh(nshell)
!   integer, intent(out) :: ash(nshell)
!   integer, intent(out) :: lsh(nshell)
   real(wp),intent(out) :: e

   real(wp)  ntot,fracz
   real(wp)  iox(86,0:2,2) ! GFN1 on 1
   integer i,j,k,m,l,ll(0:3),iat,lll,iver
   data ll /1,3,5,7/

   !     H           Initial     Orbital Occupancies                     He
   !     Li Be                                            B  C  N  O  F  Ne
   !     Na Mg                                            Al Si P  S  Cl Ar
   !     K  Ca Sc            Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
   !     Rb Sr Y             Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
   !     Cs Ba La Ce-Lu      Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
   !                                      spd shell

   ! GFN1
   data iox / &
      &1,                                                          2, & !He
      &1,2,                                         2, 2, 2, 2, 2, 2, & !Ne
      &1,2,                                         2, 2, 2, 2, 2, 2, & !Ar
      &1,2,2,             2, 2, 2, 2, 2, 2, 2, 2, 2,2, 2, 2, 2, 2, 2, & !Kr
      &1,2,2,             2, 2, 2, 2, 2, 2, 2, 2, 2,2, 2, 2, 2, 2, 2, & !Xe
      &1,2,2, 14*2       ,2, 2, 2, 2, 2, 2, 2, 2, 2,2, 2, 2, 2, 2, 2, & !Rn
      !
      &0,                                                          0, & !He
      &0,0,                                         1, 2, 3, 4, 5, 6, & !Ne
      &0,0,                                         1, 2, 3, 4, 5, 6, & !Ar
      &0,0,0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, & !Kr
      &0,0,0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, & !Xe
      &0,0,0,  14*0,   0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, & !Rn
      !
      & 0,                                                         0, & !He
      & 0,0,                                        0, 0, 0, 0, 0, 0, & !Ne
      & 0,0,                                        0, 0, 0, 0, 0, 0, & !Ar
      & 0,0,1,          2, 3, 4, 5, 6, 7, 8,  9,10, 0, 0, 0, 0, 0, 0, & !Kr
      & 0,0,1,          2, 3, 4, 5, 6, 7, 8,  9,10, 0, 0, 0, 0, 0, 0, & !Xe
      & 0,0,1, 14*1,    2, 3, 4, 5, 6, 7, 8,  9,10, 0, 0, 0, 0, 0, 0, & !Rn
      ! GFN2
      &1,                                                            2, & !He
      &1,2,                                        2  ,1  ,1.5,  2,2,2, & !Ne
      &1,2,                                        2  ,1.5,1.5,  2,2,2, & !Ar
      &1,2,1,             1, 1, 1, 1, 1, 1, 1,1,2, 2,  1.5,1.5,  2,2,2, & !Kr
      &1,2,1,             1, 1, 1, 1, 1, 1, 1,1,2, 2,  2,  2 ,   2,2,2, & !Xe
      &1,2,1, 14*1       ,1, 1, 1, 1, 1, 1, 1,1,2, 2,  2,  2 ,   2,2,2, & !Rn
      !
      &0,                                                            0, & !He
      &0,0,                                        1  , 3  ,3.5, 4,5,6, & !Ne
      &0,0,                                        1  , 2.5,3.5, 4,5,6, & !Ar
      &0,0,1,            1, 1, 1, 1, 1, 1, 1, 0,0, 1,   2.5,3.5, 4,5,6, & !Kr
      &0,0,1,            1, 1, 1, 1, 1, 1, 1, 0,0, 1,   2,  3,   4,5,6, & !Xe
      &0,0,1,  14*1,     1, 1, 1, 1, 1, 1, 1, 0,0, 1,   2,  3,   4,5,6, & !Rn
      !
      & 0,                                                          0, & !He
      & 0,0,                                        0, 0 , 0, 0, 0, 0, & !Ne
      & 0,0,                                        0, 0 , 0, 0, 0, 0, & !Ar
      & 0,0,1,          2, 3, 4, 5, 6, 7, 8, 10,10, 0, 0 , 0, 0, 0, 0, & !Kr
      & 0,0,1,          2, 3, 4, 5, 6, 7, 8, 10,10, 0, 0 , 0, 0, 0, 0, & !Xe
      & 0,0,1, 14*1,    2, 3, 4, 5, 6, 7, 8, 10,10, 0, 0 , 0, 0, 0, 0  & !Rn
      & /

   iver=1
   if(gfn_method.gt.1) iver = 2

   k=0
   e=0.0_wp
   do i=1,n
      iat=at(i)
      ntot=-1.d-6
      do m=1,ao_n(iat)
         l=ao_l(m,iat)
         k=k+1
         zsh(k)=iox(iat,l,iver)
!         lsh(k)=l
!         ash(k)=i
         ntot=ntot+zsh(k)
         if(ntot.gt.z(i)) zsh(k)=0
         e=e+ao_lev(m,iat)*zsh(k)
      enddo
   enddo
   if(abs(sum(z)-sum(zsh)).gt.1.d-4) then
      write(*,*) i,sum(z),sum(zsh)
      error stop 'internal setzshell error 2'
   endif

end subroutine setzshell

!cccccccccccccccccccccccccccccccccccccccccccccccccccc

end module scc_core
