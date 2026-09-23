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

!> general functions for core functionalities of the SCC
module xtb_scc_core
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_la, only : contract
   use xtb_mctc_lapack, only : lapack_sygvd
   use xtb_mctc_blas, only : blas_gemm, mctc_symv, mctc_gemm
   use xtb_mctc_lapack_eigensolve, only : TEigenSolver
   use xtb_type_environment, only : TEnvironment
   use xtb_type_solvation, only : TSolvation
   use xtb_solv_gbsa, only : TBorn
   use xtb_xtb_data
   use xtb_xtb_coulomb
   use xtb_xtb_dispersion
   use xtb_xtb_multipole
   use xtb_broyden
   implicit none
   private

   public :: build_h0, scc, electro, solve, solve4
   public :: fermismear, occ, occu, dmat, get_unrestricted_wiberg
   public :: get_wiberg, mpopall, mpop0, mpopao, mpop, mpopsh, qsh2qat, lpop
   public :: iniqshell, setzshell
   public :: shellPoly, h0scal


   integer, private, parameter :: mmm(20)=(/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/)


contains

!! ========================================================================
!  build GFN2 core Hamiltonian
!! ========================================================================
subroutine build_h0(hData,H0,n,at,ndim,nmat,matlist, &
   &                xyz,selfEnergy,S,aoat2,lao2,valao2,aoexp,ao2sh)
   type(THamiltonianData), intent(in) :: hData
   real(wp),intent(out) :: H0(ndim*(ndim+1)/2)
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: ndim
   integer, intent(in)  :: nmat
   integer, intent(in)  :: matlist(2,nmat)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: selfEnergy(:)
   real(wp),intent(in)  :: S(ndim,ndim)
   integer, intent(in)  :: aoat2(ndim)
   integer, intent(in)  :: lao2(ndim)
   integer, intent(in)  :: valao2(ndim)
   integer, intent(in)  :: ao2sh(ndim)
   real(wp),intent(in)  :: aoexp(ndim)

   integer  :: i,j,k,m
   integer  :: iat,jat,ish,jsh,il,jl,iZp,jZp
   real(wp) :: hdii,hdjj,hav
   real(wp) :: km

   H0=0.0_wp

   do m = 1, nmat
      i = matlist(1,m)
      j = matlist(2,m)
      k = j+i*(i-1)/2
      iat = aoat2(i)
      jat = aoat2(j)
      ish = ao2sh(i)
      jsh = ao2sh(j)
      iZp = at(iat)
      jZp = at(jat)
      il = mmm(lao2(i))
      jl = mmm(lao2(j))
      hdii = selfEnergy(ish)
      hdjj = selfEnergy(jsh)
      call h0scal(hData,il,jl,izp,jzp,valao2(i).ne.0,valao2(j).ne.0, &
      &           km)
      km = km*(2*sqrt(aoexp(i)*aoexp(j))/(aoexp(i)+aoexp(j)))**hData%wExp
      hav = 0.5_wp*(hdii+hdjj)* &
      &      shellPoly(hData%shellPoly(il, iZp), hData%shellPoly(jl, jZp), &
      &                hData%atomicRad(iZp), hData%atomicRad(jZp),xyz(:,iat),xyz(:,jat))
      H0(k) = S(j,i)*km*hav
   enddo
!  diagonal
   k=0
   do i=1,ndim
      k=k+i
      iat = aoat2(i)
      ish = ao2sh(i)
      il = mmm(lao2(i))
      H0(k) = selfEnergy(ish)
   enddo

end subroutine build_h0

!> build isotropic H1/Fockian
subroutine buildIsotropicH1(n, at, ndim, nshell, nmat, matlist, H, &
      & H0, S, shellShift, aoat2, ao2sh)
   use xtb_mctc_convert, only : autoev,evtoau
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: ndim
   integer, intent(in)  :: nshell
   integer, intent(in)  :: nmat
   integer, intent(in)  :: matlist(2,nmat)
   real(wp),intent(in)  :: H0(ndim*(1+ndim)/2)
   real(wp),intent(in)  :: S(ndim,ndim)
   real(wp),intent(in)  :: shellShift(nshell)
   integer, intent(in)  :: aoat2(ndim)
   integer, intent(in)  :: ao2sh(ndim)
   real(wp),intent(out) :: H(ndim,ndim)

   integer  :: m,i,j,k
   integer  :: ishell,jshell
   integer  :: ii,jj,kk
   real(wp) :: dum
   real(wp) :: eh1,t8,t9,tgb,h1

   H = 0.0_wp

   !$omp parallel do default(none) &
   !$omp private(m, i, j, k, ishell, jshell, eh1, H1) &
   !$omp shared(H, H0, S, matlist, nmat, ao2sh, shellShift) &
   !$omp schedule(static)
   do m = 1, nmat
      i = matlist(1,m)
      j = matlist(2,m)
      k = j+i*(i-1)/2
      ishell = ao2sh(i)
      jshell = ao2sh(j)
      ! SCC terms
      eh1 = autoev*(shellShift(ishell) + shellShift(jshell))
      H1 = -S(j,i)*eh1*0.5_wp
      H(j,i) = H0(k) + H1
      H(i,j) = H(j,i)
   enddo

end subroutine buildIsotropicH1

!> build isotropic & anisotropic H1/Fockian
subroutine buildIsoAnisotropicH1(n,at,ndim,nshell,nmat,ndp,nqp,matlist,mdlst,mqlst,&
                         H,H0,S,shellShift,dpint,qpint,vs,vd,vq,aoat2,ao2sh)
   use xtb_mctc_convert, only : autoev,evtoau
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
   real(wp),intent(in)  :: shellShift(nshell)
   real(wp),intent(in)  :: dpint(3,ndim,ndim)
   real(wp),intent(in)  :: qpint(6,ndim,ndim)
   real(wp),intent(in)  :: vs(n)
   real(wp),intent(in)  :: vd(3,n)
   real(wp),intent(in)  :: vq(6,n)
   integer, intent(in)  :: aoat2(ndim)
   integer, intent(in)  :: ao2sh(ndim)
   real(wp),intent(inout) :: H(ndim,ndim)

   integer, external :: lin
   integer  :: m,i,j,k,l
   integer  :: ii,jj,kk
   integer  :: ishell,jshell
   real(wp) :: dum,eh1,t8,t9,tgb

   !$omp parallel default(none) &
   !$omp private(m, i, j, k, l, ii, jj, dum, eh1) &
   !$omp shared(matlist, mdlst, mqlst, ao2sh, aoat2, nmat, ndp, nqp) &
   !$omp shared(S, H, H0, shellShift, vs, vd, vq, dpint, qpint)

   !> overlap dependent terms
   !$omp do schedule(static)
   do m=1,nmat
      i=matlist(1,m)
      j=matlist(2,m)
      k=j+i*(i-1)/2
      dum = S(j,i)*autoev*0.5_wp

      ii = ao2sh(i)
      jj = ao2sh(j)
      ! SCC terms (isotropic; must be first!)
      eh1 = -dum*(shellShift(ii) + shellShift(jj))
      H(j,i) = H0(k) + eh1

      ii=aoat2(i)
      jj=aoat2(j)
      ! CAMM potential
      eh1=dum*(vs(ii)+vs(jj))
      H(j,i)=H(j,i)+eh1

      H(i,j)=H(j,i)
   enddo

   !> dipolar terms
   !$omp do schedule(static)
   do m=1,ndp
      i=mdlst(1,m)
      j=mdlst(2,m)
      k=lin(j,i)
      ii=aoat2(i)
      jj=aoat2(j)
      eh1=0.0_wp
      do l=1,3
         eh1=eh1+dpint(l,i,j)*(vd(l,ii)+vd(l,jj))
      enddo
      eh1=0.50_wp*eh1*autoev
      H(i,j)=H(i,j)+eh1
      H(j,i)=H(i,j)
   enddo

   !> quadrupole-dependent terms
   !$omp do schedule(static)
   do m=1,nqp
      i=mqlst(1,m)
      j=mqlst(2,m)
      ii=aoat2(i)
      jj=aoat2(j)
      k=lin(j,i)
      eh1=0.0_wp
      ! note: these come in the following order
      ! xx, yy, zz, xy, xz, yz
      do l=1,6
         eh1=eh1+qpint(l,i,j)*(vq(l,ii)+vq(l,jj))
      enddo
      eh1=0.50_wp*eh1*autoev
      H(i,j)=H(i,j)+eh1
      H(j,i)=H(i,j)
   enddo

   !$omp end parallel

end subroutine buildIsoAnisotropicH1


!> self consistent charge iterator
subroutine scc(env,xtbData,solver,n,nel,nopen,ndim,ndp,nqp,nmat,nshell, &
      &        at,matlist,mdlst,mqlst,aoat2,ao2sh,ash, &
      &        q,dipm,qp,qq,qlmom,qsh,zsh, &
      &        xyz,aes, &
      &        cm5,cm5a,gborn,solvation, &
      &        scD4, &
      &        broy,broydamp,damp0, &
      &        pcem,shellShift,externShift, &
      &        et,focc,focca,foccb,efa,efb, &
      &        eel,ees,eaes,epol,ed,epcem,egap,emo,ihomo,ihomoa,ihomob, &
      &        H0,H,S,dpint,qpint,X,P,ies, &
      &        maxiter,startpdiag,scfconv,qconv, &
      &        minpr,pr, &
      &        fail,jter)
   use xtb_mctc_convert, only : autoev,evtoau
   use xtb_mctc_lapack_trf, only : mctc_potrf

   use xtb_disp_dftd4,  only: disppot,edisp_scc
   use xtb_aespot, only : gfn2broyden_diff,gfn2broyden_out,gfn2broyden_save, &
   &                  mmompop,aniso_electro,setvsdq
   use xtb_embedding, only : electro_pcem

   character(len=*), parameter :: source = 'scc_core'

   type(TEnvironment), intent(inout) :: env

   type(TxTBData), intent(in) :: xtbData
   type(TEigenSolver), intent(inout) :: solver

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
   integer, intent(in)  :: ash(:)
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
   type(TxTBMultipole), intent(in), optional :: aes
   real(wp),intent(in)   :: xyz(3,n)
   real(wp), allocatable :: vs(:)
   real(wp), allocatable :: vd(:, :)
   real(wp), allocatable :: vq(:, :)
!! ------------------------------------------------------------------------
!  continuum solvation model GBSA
   class(TSolvation), allocatable, intent(inout) :: solvation
   real(wp),intent(in)    :: cm5a(n)
   real(wp),intent(inout) :: cm5(n)
   real(wp),intent(inout) :: gborn
!! ------------------------------------------------------------------------
!  selfconsistent DFT-D4 dispersion correction
   type(TxTBDispersion), intent(inout), optional :: scD4
!! ------------------------------------------------------------------------
!  point charge embedding potentials
   logical, intent(in)    :: pcem
   real(wp),intent(inout) :: shellShift(nshell)
   real(wp),intent(inout) :: externShift(nshell)
   real(wp), allocatable :: atomicShift(:)
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
!  Factorized overlap to avoid multiple factorizations
   real(wp), allocatable :: S_factorized(:,:)
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
   real(wp),intent(out)   :: H(ndim,ndim)
   real(wp),intent(inout) :: P(ndim,ndim)
   real(wp),intent(inout) :: X(ndim,ndim)
   real(wp),intent(in)    :: S(ndim,ndim)
   real(wp),intent(in)    :: dpint(3,ndim,ndim)
   real(wp),intent(in)    :: qpint(6,ndim,ndim)
   type(TxTBCoulomb), intent(inout) :: ies

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
   ! Response correction
   integer, parameter :: maxnewton = 3
   integer :: nnewton
   logical :: use_newton, newton_ok
   real(wp), allocatable :: newton_step(:), newton_base(:)

   allocate(S_factorized(ndim, ndim), source = 0.0_wp )
   S_factorized = S
   call mctc_potrf(env, S_factorized)

   converged = .false.
   lastdiag = .false.
   ! number of iterations for this iterator
   thisiter = maxiter - jter

   ! Setup response solver to accelerate SCC when qconv is tiny.
   nnewton = 0
   use_newton = qconv <= 1.0e-5_wp
   if (allocated(solvation)) then
      select type(solvation)
      type is (TBorn)
         use_newton = use_newton .or. solvation%alpbet > 0.0_wp
      class default
         use_newton = .false.
      end select
   end if

   damp = damp0
   if (present(aes)) then
      nbr = nshell + 9*n
      allocate(vs(n), vd(3, n), vq(6, n))
   else
      nbr = nshell
   end if
!  broyden data storage and init
   allocate(newton_step(nbr), newton_base(nbr))
   allocate( df(thisiter,nbr),u(thisiter,nbr),a(thisiter,thisiter), &
   &         dq(nbr),dqlast(nbr),qlast_in(nbr),omega(thisiter), &
   &         q_in(nbr),atomicShift(n), source = 0.0_wp )

!! ------------------------------------------------------------------------
!  Iteration entry point
   scc_iterator: do iter = 1, thisiter

   ! set up ES potential
   atomicShift(:) = 0.0_wp
   shellShift(:) = externShift
   call ies%addShift(q, qsh, atomicShift, shellShift)
   ! compute potential intermediates
   if (present(aes)) then
      call setvsdq(aes,n,at,xyz,q,dipm,qp,aes%gab3,aes%gab5,vs,vd,vq)
   end if
   ! Solvation contributions
   if (allocated(solvation)) then
      cm5(:) = q + cm5a
      call solvation%addShift(env, cm5, qsh, atomicShift, shellShift)
   end if
   ! self consistent dispersion contributions
   if (present(scD4)) then
      call scD4%addShift(at, q, atomicShift)
   end if
   ! expand all atomic potentials to shell resolved potentials
   call addToShellShift(ash, atomicShift, shellShift)

   ! build the charge dependent Hamiltonian
   if (present(aes)) then
      call buildIsoAnisotropicH1(n,at,ndim,nshell,nmat,ndp,nqp,matlist,mdlst,mqlst,&
         & H,H0,S,shellShift,dpint,qpint,vs,vd,vq,aoat2,ao2sh)
   else
      call buildIsotropicH1(n,at,ndim,nshell,nmat,matlist,H,H0,S, &
         & shellShift,aoat2,ao2sh)
   end if

   ! ------------------------------------------------------------------------
   ! solve HC=SCemo(X,P are scratch/store)
   ! solution is in H(=C)/emo
   ! ------------------------------------------------------------------------
   fulldiag=.false.
   if(iter.lt.startpdiag) fulldiag=.true.
   if(lastdiag )          fulldiag=.true.

   !call solve(fulldiag,ndim,ihomo,scfconv,H,S,X,P,emo,fail)

   call solver%fact_solve(env, H, S_factorized, emo)
   call env%check(fail)
   if(fail)then
      call env%error("Diagonalization of Hamiltonian failed", source)
      return
   endif

   if(ihomo+1.le.ndim.and.ihomo.ge.1)egap=emo(ihomo+1)-emo(ihomo)
   ! automatic reset to small value
   if(egap.lt.0.1_wp.and.iter.eq.0) broydamp=0.03_wp

   ! Fermi smearing
   if(et.gt.0.1_wp)then
      ! convert restricted occ first to alpha/beta
      if(nel.gt.0) then
         call occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
      else
         focca=0.0_wp
         foccb=0.0_wp
         ihomoa=0
         ihomob=0
      endif
      call fermismear(.false.,ndim,ihomoa,et,emo,focca,nfoda,efa,ga)
      call fermismear(.false.,ndim,ihomob,et,emo,foccb,nfodb,efb,gb)
      focc = focca + foccb
   else
      ga = 0.0_wp
      gb = 0.0_wp
   endif

   ! save q
   q_in(1:nshell)=qsh(1:nshell)
   if (present(aes)) then
      k=nshell
      call gfn2broyden_save(n,k,nbr,dipm,qp,q_in)
   end if

   ! density matrix
   call dmat(ndim,focc,H,P)

   ! new q
   call mpopsh(n,ndim,nshell,ao2sh,S,P,qsh)
   qsh = zsh - qsh

   ! qat from qsh
   call qsh2qat(ash,qsh,q)

   eold=eel
   call electro(n,at,ndim,nshell,ies,H0,P,q,qsh,ees,eel)
   ! multipole electrostatic
   if (present(aes)) then
      call mmompop(n,ndim,aoat2,xyz,p,s,dpint,qpint,dipm,qp)
      ! evaluate energy
      call aniso_electro(aes,n,at,xyz,q,dipm,qp,aes%gab3,aes%gab5,eaes,epol)
      eel=eel+eaes+epol
   end if

   ! Self consistent dispersion
   if (present(scD4)) then
      call scD4%getEnergy(at, q, ed)
      eel = eel + ed
   endif

   ! point charge contribution
   if(pcem) then
      call electro_pcem(nshell,qsh,externShift,epcem,eel)
   end if

   ! new cm5 charges and gborn energy
   if (allocated(solvation)) then
      cm5=q+cm5a
      call solvation%getEnergy(env, cm5, qsh, gborn)
      eel = eel + gborn
   end if

   ! add el. entropies*T
   eel=eel+ga+gb

   ! ------------------------------------------------------------------------
   ! check for energy convergence
   econverged = abs(eel - eold) < scfconv
   ! ------------------------------------------------------------------------

   dq(1:nshell)=qsh(1:nshell)-q_in(1:nshell)
   if (present(aes)) then
      k=nshell
      call gfn2broyden_diff(n,k,nbr,dipm,qp,q_in,dq) ! CAMM case
   end if
   rmsq=sum(dq(1:nbr)**2)/dble(n)
   rmsq=sqrt(rmsq)

   ! ------------------------------------------------------------------------
   ! end of SCC convergence part
   qconverged = rmsq < qconv
   ! ------------------------------------------------------------------------

   ! SCC convergence acceleration
   if(.not.broy)then

      ! simple damp
      if(iter.gt.0) then
         omegap=egap
         ! monopoles only
         do i=1,nshell
            qsh(i)=damp*qsh(i)+(1.0_wp-damp)*q_in(i)
         enddo
         if (present(aes)) then
            ! CAMM
            k=nshell
            do i=1,n
               do j=1,3
                  k=k+1
                  dipm(j,i)=damp*dipm(j,i)+(1.0_wp-damp)*q_in(k)
               enddo
               do j=1,6
                  k=k+1
                  qp(j,i)=damp*qp(j,i)+(1.0_wp-damp)*q_in(k)
               enddo
            enddo
         end if
         if(eel-eold.lt.0.0_wp) then
            damp=damp*1.15_wp
         else
            damp=damp0
         endif
         damp=min(damp,1.0_wp)
         if(egap.lt.1.0_wp)damp=min(damp,0.5_wp)
      endif

   else

      ! Broyden mixing
      omegap=0.0_wp

      newton_ok = .false.
      newton_base(:nbr) = q_in(:nbr)
      if (use_newton .and. nnewton < maxnewton .and. rmsq < 1.0e-3_wp .and. .not.qconverged) then
         call scc_coupled_newton(env, minpr, &
                                 n, nbr, nshell, nmat, ndp, nqp, &
                                 q_in, dq, &
                                 et, &
                                 focc, focca, foccb, &
                                 emo, &
                                 newton_step(:nbr), newton_ok, ies, aes, solvation, &
                                 xyz, H, S, &
                                 dpint, qpint, &
                                 matlist, mdlst, mqlst, &
                                 at, ash, aoat2, ao2sh)
      end if

      call broyden(nbr,q_in,qlast_in,dq,dqlast,iter,thisiter,broydamp,omega,df,u,a)

      if (newton_ok) then
         q_in = newton_base(:nbr) + newton_step(:nbr)
         nnewton = nnewton + 1
      end if

      qsh(1:nshell)=q_in(1:nshell)
      if (present(aes)) then
         k=nshell
         call gfn2broyden_out(n,k,nbr,q_in,dipm,qp) ! CAMM case
      end if
      if(iter.gt.1) omegap=omega(iter-1)
   endif ! Broyden?

   call qsh2qat(ash, qsh, q) !new qat

   if(allocated(solvation)) cm5 = q+cm5a

   if(minpr)write(env%unit,'(i4,F15.7,E14.6,E11.3,f8.2,2x,f8.1,l3)') &
   &  iter+jter,eel,eel-eold,rmsq,egap,omegap,fulldiag
   qq=q

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

end subroutine scc

! Solve (I - d q_out / d q_in) delta_q = scc_residual.
subroutine scc_coupled_newton(env, minpr, &
                              nat, nbr, nshell, nmat, ndp, nqp, &
                              q0, residual, &
                              et, &
                              focc, focca, foccb, &
                              emo, &
                              step, ok, &
                              ies, aes, solvation, &
                              xyz, H, S, &
                              dpint, qpint, &
                              matlist, mdlst, mqlst, &
                              at, ash, aoat2, ao2sh)
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use xtb_mctc_blas, only : mctc_dot, mctc_gemv
   use xtb_mctc_constants, only : kB
   use xtb_mctc_convert, only : autoev

   type(TEnvironment), intent(inout) :: env
   logical,  intent(in)  :: minpr
   integer,  intent(in)  :: nat, nbr, nshell, nmat, ndp, nqp
   real(wp), intent(in)  :: q0(:), residual(:)
   real(wp), intent(in)  :: et
   real(wp), intent(in)  :: focc(:), focca(:), foccb(:)
   real(wp), intent(in)  :: emo(:)
   real(wp), intent(out) :: step(:)
   logical,  intent(out) :: ok
   ! bypass variables
   type(TxTBCoulomb),   intent(inout)              :: ies
   type(TxTBMultipole), intent(in),    optional    :: aes
   class(TSolvation), target, intent(inout), allocatable :: solvation
   real(wp), intent(in)  :: xyz(:,:), H(:,:), S(:,:)
   real(wp), intent(in)  :: dpint(:,:,:), qpint(:,:,:)
   integer,  intent(in)  :: matlist(:,:), mdlst(:,:), mqlst(:,:)
   integer,  intent(in)  :: at(:), ash(:), aoat2(:), ao2sh(:)
   ! parameter
   integer, parameter :: maxk = 20
   ! local variables
   real(wp) :: beta, normw0, normw, gap, tolerance, tmp
   integer :: col, row, mo, mj, nk, nocc, ndim
   ! dynamic memory
   real(wp), allocatable :: v(:,:), hh(:,:), cs(:), sn(:), rhs(:), y(:), w(:)
   real(wp), allocatable :: fdiv(:,:), fprime(:,:)
   real(wp), allocatable :: hsp(:,:), val(:)

   ok = .false.
   step = 0.0_wp

   beta = norm2(residual)
   if (beta == 0.0_wp) return

   nocc = count(focc > 1.0e-14_wp)
   ndim = size(focc, 1)

   allocate(v(nbr,maxk+1), &
            hh(maxk+1,maxk), &
            cs(maxk), sn(maxk), &
            rhs(maxk+1), &
            y(maxk), &
            w(nbr), &
            fdiv(ndim,nocc), &
            fprime(ndim,2), &
            hsp(ndim,ndim), &
            val(maxk+1), &
            source = 0.0_8)

   if (et > 0.1_wp) then
      fprime(:,1) = -focca * (1.0_wp - focca) / (kB*autoev*et)
      fprime(:,2) = -foccb * (1.0_wp - foccb) / (kB*autoev*et)
   end if
   do mj = 1, nocc
      do mo = 1, ndim
         gap = emo(mo)-emo(mj)
         if (abs(gap) > 1.0e-8_wp) then
            fdiv(mo,mj) = (focc(mo) - focc(mj)) / gap
         else
            fdiv(mo,mj) = 0.5_wp * sum(fprime(mo,:) + fprime(mj,:))
         end if
      end do
   end do

   ! tighten inexact Newton solve near fixed point
   tolerance = max(1.0e-4_wp, min(1.0e-2_wp, sqrt(beta))) * beta

   rhs(1) = beta
   v(:,1) = residual / beta

   do col = 1, maxk
      call scc_jvp(env, ies, aes, solvation, &
                   nat, nbr, nshell, nmat, ndp, nqp, &
                   q0, v(:,col), &
                   fdiv, fprime, &
                   et, focc, focca, foccb, emo, &
                   xyz, H, S, &
                   dpint, qpint, &
                   matlist, mdlst, mqlst, &
                   at, ash, aoat2, ao2sh, &
                   w) ! long call

      ! DGKS correction
      normw0 = norm2(w)

      call mctc_gemv(v, w, val, alpha=1.0_wp, beta=0.0_wp, trans='T')  ! val = V^T w
      hh(1:col,col) = hh(1:col,col) + val(1:col)
      call mctc_gemv(v, val, w, alpha=-1.0_wp, beta=1.0_wp, trans='N')  ! w = w - V val
      normw = norm2(w)

      if (normw < 0.7071_wp * normw0) then
         call mctc_gemv(v, w, val, alpha=1.0_wp, beta=0.0_wp, trans='T')  ! val = V^T w
         hh(1:col,col) = hh(1:col,col) + val(1:col)
         call mctc_gemv(v, val, w, alpha=-1.0_wp, beta=1.0_wp, trans='N')  ! w = w - V val
         normw = norm2(w)
      end if

      ! rotation
      hh(col+1,col) = normw
      if (normw > 1.0e-14_wp) then
        v(:,col+1) = w / normw
      end if
      do row = 1, col-1
         tmp = cs(row) * hh(row,col) + sn(row) * hh(row+1,col)
         hh(row+1,col) = -sn(row) * hh(row,col) + cs(row) * hh(row+1,col)
         hh(row,col) = tmp
      end do
      tmp = hypot(hh(col,col), hh(col+1,col))
      if (tmp < 1.0e-14_wp) return

      cs(col) = hh(col,col) / tmp
      sn(col) = hh(col+1,col) / tmp
      hh(col,col) = tmp
      hh(col+1,col) = 0.0_wp
      rhs(col+1) = -sn(col) * rhs(col)
      rhs(col) = cs(col) * rhs(col)
      if (abs(rhs(col+1)) < tolerance .or. normw < 1.0e-14_wp) exit
   end do

   nk = min(col, maxk)
   if (abs(rhs(nk+1)) > tolerance) return
   y(:nk) = rhs(:nk)
   do row = nk, 1, -1
      y(row) = (y(row) - dot_product(hh(row,row+1:nk), y(row+1:nk))) / hh(row,row)
   end do
   call mctc_gemv(v(:,:nk), y(:nk), step)

   ok = all(ieee_is_finite(step)) .and. norm2(step) < 100.0_wp*beta
   if (minpr) write(env%unit, '(a,i4,a,es10.2)') ' Newton response: ', nk, &
      & ' products, relative residual ', abs(rhs(nk+1)) / beta

end subroutine scc_coupled_newton

! Apply the SCC residual Jacobian
subroutine scc_jvp(env, ies, aes, solvation, &
                   nat, nbr, nshell, nmat, ndp, nqp, &
                   q0, direction, &
                   fdiv, fprime, &
                   et, focc, focca, foccb, emo, &
                   xyz, H, S, &
                   dpint, qpint, &
                   matlist, mdlst, mqlst, &
                   at, ash, aoat2, ao2sh, &
                   jv)
   use xtb_aespot, only : gfn2broyden_out, gfn2broyden_save, mmompop, setvsdq
   use xtb_mctc_blas, only : mctc_dot, mctc_gemv
   type(TEnvironment),  intent(inout)              :: env
   type(TxTBCoulomb),   intent(inout)              :: ies
   type(TxTBMultipole), intent(in),    optional    :: aes
   class(TSolvation),   intent(inout), allocatable :: solvation
   integer,  intent(in)  :: nat, nbr, nshell, nmat, ndp, nqp
   real(wp), intent(in)  :: q0(:), direction(:)
   real(wp), intent(in)  :: fdiv(:,:), fprime(:,:)
   real(wp), intent(in)  :: et, focc(:), focca(:), foccb(:), emo(:)
   real(wp), intent(in)  :: xyz(:,:), H(:,:), S(:,:)
   real(wp), intent(in)  :: dpint(:,:,:), qpint(:,:,:)
   integer, intent(in)   :: matlist(:,:), mdlst(:,:), mqlst(:,:)
   integer, intent(in)   :: at(:), ash(:), aoat2(:), ao2sh(:)
   real(wp), intent(out) :: jv(:)
   ! local variables
   integer :: nocc, nvirt, ndim
   real(wp) :: mu, weight, gap
   integer :: mo, mj, spin, offset
   logical :: ov_response
   ! dynamic memory
   real(wp), allocatable :: dh(:,:), dm(:,:), scratch(:,:), work(:)
   real(wp), allocatable :: qat0(:), dqat(:), ashift(:), sshift(:)
   real(wp), allocatable :: ddip(:,:), dquad(:,:), ds(:), dd(:,:), dqq(:,:)
   real(wp), allocatable :: hzero(:), de(:)

   nocc = count(focc > 1.0e-14_wp)
   ndim = size(focc, 1)
   nvirt = ndim - nocc

   allocate(ddip(3,nat), dquad(6,nat), &
            ds(nat), dd(3,nat), dqq(6,nat), &
            hzero(ndim*(ndim+1)/2), de(ndim), &
            dh(ndim,ndim), &
            dm(ndim,nocc), &
            scratch(ndim,nocc), &
            work(ndim), &
            qat0(nat), dqat(nat), &
            ashift(nat), sshift(nshell), &
            source = 0.0_wp)

   call qsh2qat(ash, q0(:nshell), qat0)
   call qsh2qat(ash, direction(:nshell), dqat)
   call mctc_symv(ies%jmat, direction(:nshell), sshift)
   if (allocated(ies%thirdOrder%atomicGam)) then
      ashift = 2.0_wp * ies%thirdOrder%atomicGam * qat0 * dqat
   end if
   if (allocated(ies%thirdOrder%shellGam)) then
      sshift = sshift + 2.0_wp * ies%thirdOrder%shellGam * q0(:nshell) * direction(:nshell)
   end if
   if (allocated(solvation)) then
      call solvation%addShift(env, dqat, direction(:nshell), ashift, sshift)
   end if

   ! TODO: include self-consistent D4 contribution

   call addToShellShift(ash, ashift, sshift)

   if (present(aes)) then
      offset = nshell
      call gfn2broyden_out(nat, offset, nbr, direction, ddip, dquad)
      call setvsdq(aes, nat, at, xyz, dqat, ddip, dquad, aes%gab3, aes%gab5, ds, dd, dqq)
      call buildIsoAnisotropicH1(nat, at, ndim, nshell, nmat, ndp, nqp, matlist, mdlst, mqlst, &
                                 dh, hzero, S, sshift, dpint, qpint, ds, dd, dqq, aoat2, ao2sh)
   else
      call buildIsotropicH1(nat, at, ndim, nshell, nmat, matlist, dh, hzero, S, sshift, aoat2, ao2sh)
   end if

   call mctc_gemm(dh, H(:,:nocc), scratch)

   ov_response = et <= 0.1_wp
   if (ov_response) then
      ! response at zero temperature:
      !   dP = C_v [f_v-f_o]/(e_v-e_o) (C_v^T dH C_o) C_o^T + response_H.c.
      call mctc_gemm(H(:,nocc+1:), scratch, dm(:nvirt,:), transa='t')

      do mo = 1, nvirt
         do mj = 1, nocc
            gap = emo(nocc+mo) - emo(mj)
            if (abs(gap) > 1.0e-8_wp) then
               dm(mo,mj) = dm(mo,mj) * (focc(nocc+mo) - focc(mj)) / gap
            else
               dm(mo,mj) = 0.0_wp
            end if
         end do
      end do

      call mctc_gemm(H(:,nocc+1:), dm(:nvirt,:), scratch)
      call mctc_gemm(scratch, H(:,:nocc), dh, transb='t')
   else
      ! Finite-temperature response
      call mctc_gemm(H, scratch, dm, transa='t')

      de = 0.0_wp
      do mo = 1, ndim
         call mctc_gemv(dh, H(:,mo), work)
         de(mo) = mctc_dot(H(:,mo), work)
      end do
      do mo = 1, ndim
         do mj = 1, nocc
            gap = emo(mo)-emo(mj)
            if (mo /= mj .and. abs(gap) > 1.0e-8_wp) then
               dm(mo,mj) = dm(mo,mj) * (focca(mo)-focca(mj) + &
                                        foccb(mo)-foccb(mj)) / gap
            else
               dm(mo,mj) = 0.0_wp
            end if
         end do
      end do
      dm(:nocc,:) = 0.5_wp * dm(:nocc,:)

      do spin = 1, 2
         weight = sum(fprime(:,spin))
         if (abs(weight) < 1.0e-14_wp) cycle
         mu = mctc_dot(fprime(:,spin), de) / weight
         do mo = 1, nocc
            dm(mo,mo) = dm(mo,mo) + 0.5_wp * fprime(mo,spin) * (de(mo) - mu)
         end do
      end do

      call mctc_gemm(H, dm, scratch)
      call mctc_gemm(scratch, H(:,:nocc), dh, transb='t')
   end if

   do mj = 1, ndim
      dh(mj,mj) = 2.0_wp * dh(mj,mj)
      do mo = 1, mj-1
         dh(mo,mj) = dh(mo,mj) + dh(mj,mo)
         dh(mj,mo) = dh(mo,mj)
      enddo
   enddo
   call mpopsh(nat, ndim, nshell, ao2sh, S, dh, jv)
   jv = -jv

   if (present(aes)) then
      call mmompop(nat, ndim, aoat2, xyz, dh, S, dpint, qpint, ddip, dquad)
      offset = nshell
      call gfn2broyden_save(nat, offset, nbr, ddip, dquad, jv)
   end if

   jv = direction - jv
end subroutine scc_jvp

!> H0 off-diag scaling
subroutine h0scal(hData,il,jl,izp,jzp,valaoi,valaoj,km)
  !$acc routine seq
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in)  :: il
   integer, intent(in)  :: jl
   integer, intent(in)  :: izp
   integer, intent(in)  :: jzp
   logical, intent(in)  :: valaoi
   logical, intent(in)  :: valaoj
   real(wp),intent(out) :: km
   real(wp) :: den, enpoly

   km = 0.0_wp

!  valence
   if(valaoi.and.valaoj) then
      den=(hData%electronegativity(izp)-hData%electronegativity(jzp))**2
      enpoly = (1.0_wp+hData%enScale(jl-1,il-1)*den*(1.0_wp+hData%enScale4*den))
      km=hData%kScale(jl-1,il-1)*enpoly*hData%pairParam(izp,jzp)
      return
   endif

!  "DZ" functions (on H for GFN or 3S for EA calc on all atoms)
   if((.not.valaoi).and.(.not.valaoj)) then
      km=hData%kDiff
      return
   endif
   if(.not.valaoi.and.valaoj) then
      km=0.5_wp*(hData%kScale(jl-1,jl-1)+hData%kDiff)
      return
   endif
   if(.not.valaoj.and.valaoi) then
      km=0.5_wp*(hData%kScale(il-1,il-1)+hData%kDiff)
   endif


end subroutine h0scal


!! ========================================================================
!  total energy for GFN1
!! ========================================================================
pure subroutine electro(n,at,nbf,nshell,ies,H0,P,dq,dqsh,es,scc)
   use xtb_mctc_convert, only : evtoau
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   integer, intent(in) :: nbf
   integer, intent(in) :: nshell
   real(wp),intent(in)  :: H0(nbf*(nbf+1)/2)
   real(wp),intent(in)  :: P (nbf,nbf)
   type(TxTBCoulomb), intent(inout) :: ies
   real(wp),intent(in)  :: dq(n)
   real(wp),intent(in)  :: dqsh(nshell)
   real(wp),intent(out) :: es
   real(wp),intent(out) :: scc
   real(wp) :: ehb ! not used

   integer  :: i,j,k
   real(wp) :: h,t

   ! ES energy in Eh
   call ies%getEnergy(dq, dqsh, es)

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
!  S(R) enhancement factor
!! ========================================================================
pure function shellPoly(iPoly,jPoly,iRad,jRad,xyz1,xyz2)
   use xtb_mctc_convert, only : aatoau
   !$acc routine seq
   real(wp), intent(in) :: iPoly,jPoly
   real(wp), intent(in) :: iRad,jRad
   real(wp), intent(in) :: xyz1(3),xyz2(3)
   real(wp) :: shellPoly
   real(wp) :: rab,k1,rr,r,rf1,rf2,dx,dy,dz
   real(wp), parameter :: a = 0.5_wp ! R^a dependence 0.5 in GFN1

   dx=xyz1(1)-xyz2(1)
   dy=xyz1(2)-xyz2(2)
   dz=xyz1(3)-xyz2(3)

   rab=sqrt(dx**2+dy**2+dz**2)

   ! this sloppy conv. factor has been used in development, keep it
   rr=jRad+iRad

   r=rab/rr

   rf1=1.0_wp+0.01_wp*iPoly*sqrt(r)
   rf2=1.0_wp+0.01_wp*jPoly*sqrt(r)

   shellPoly= rf1*rf2

end function shellPoly

!! ========================================================================
!  set up Coulomb potential due to 2nd order fluctuation
!! ========================================================================
pure subroutine setespot(nshell, qsh, jmat, shellShift)

   !> Dimension of the Coulomb matrix
   integer, intent(in) :: nshell

   !> Shell resolved partial charges
   real(wp),intent(in) :: qsh(:)

   !> Coulomb matrix
   real(wp),intent(in) :: jmat(:, :)

   !> shell-resolved potential shift
   real(wp),intent(inout) :: shellShift(:)

   !call blas_symv('l', nshell, 1.0_wp, jmat, nshell, qsh, 1, 1.0_wp, shellShift, 1)
   call mctc_symv(jmat, qsh, shellShift, beta=1.0_wp)

end subroutine setespot


pure subroutine addToShellShift(shellAtom, atomicShift, shellShift)

   !> Atom each shell is centered on
   integer, intent(in) :: shellAtom(:)

   !> Atomic potential shift
   real(wp), intent(in) :: atomicShift(:)

   !> Shell-resolved potential shift
   real(wp), intent(inout) :: shellShift(:)

   integer :: ii

   do ii = 1, size(shellShift, dim=1)
      shellShift(ii) = shellShift(ii) + atomicShift(shellAtom(ii))
   end do

end subroutine addToShellShift


!! ========================================================================
!  eigenvalue solver single-precision
!! ========================================================================
subroutine solve4(full,ndim,ihomo,acc,H,S,X,P,e,fail)
   use xtb_mctc_accuracy, only : sp
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
      call lapack_sygvd(1,'v','u',ndim,h4,ndim,p4,ndim,e4,aux4, &!workspace query
     &           -1,iwork,liwork,info)
      lwork=int(aux4(1))
      liwork=iwork(1)
      deallocate(aux4,iwork)
      allocate (aux4(lwork),iwork(liwork))              !do it
      call lapack_sygvd(1,'v','u',ndim,h4,ndim,p4,ndim,e4,aux4, &
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
!        call blas_gemm('N','N',nbf,nbf,nbf,1.0d0,H,nbf,X,nbf,0.0d0,P,nbf)
!        call blas_gemm('T','N',nbf,nbf,nbf,1.0d0,X,nbf,P,nbf,0.0d0,H,nbf)
!        call SYEV('V','U',nbf,H,nbf,e,aux,lwork,info)
!        if(info.ne.0) error stop 'diag error'
!        call blas_gemm('N','N',nbf,nbf,nbf,1.0d0,X,nbf,H,nbf,0.0d0,P,nbf)
!        H = P
!        deallocate(aux)
!     endif
!                                                     call timing(t1,w1)
!                                    call prtime(6,t1-t0,w1-w0,'dsygvd')

   else
!                                                     call timing(t0,w0)
!     go to MO basis using trafo(X) from first iteration (=full diag)
!      call blas_gemm('N','N',ndim,ndim,ndim,1.d0,H4,ndim,X4,ndim,0.d0,P4,ndim)
!      call blas_gemm('T','N',ndim,ndim,ndim,1.d0,X4,ndim,P4,ndim,0.d0,H4,ndim)
!                                                     call timing(t1,w1)
!                       call prtime(6,1.5*(t1-t0),1.5*(w1-w0),'3xdgemm')
!                                                     call timing(t0,w0)
!      call pseudodiag(ndim,ihomo,H4,e4)
!                                                     call timing(t1,w1)
!                                call prtime(6,t1-t0,w1-w0,'pseudodiag')

!     C = X C', P=scratch
!      call blas_gemm('N','N',ndim,ndim,ndim,1.d0,X4,ndim,H4,ndim,0.d0,P4,ndim)
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
      call lapack_sygvd(1,'v','u',ndim,h,ndim,p,ndim,e,aux, &!workspace query
     &           -1,iwork,liwork,info)
      lwork=int(aux(1))
      liwork=iwork(1)
      deallocate(aux,iwork)
      allocate (aux(lwork),iwork(liwork))              !do it
      call lapack_sygvd(1,'v','u',ndim,h,ndim,p,ndim,e,aux, &
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
!        call blas_gemm('N','N',nbf,nbf,nbf,1.0d0,H,nbf,X,nbf,0.0d0,P,nbf)
!        call blas_gemm('T','N',nbf,nbf,nbf,1.0d0,X,nbf,P,nbf,0.0d0,H,nbf)
!        call SYEV('V','U',nbf,H,nbf,e,aux,lwork,info)
!        if(info.ne.0) error stop 'diag error'
!        call blas_gemm('N','N',nbf,nbf,nbf,1.0d0,X,nbf,H,nbf,0.0d0,P,nbf)
!        H = P
!        deallocate(aux)
!     endif
!                                                     call timing(t1,w1)
!                                    call prtime(6,t1-t0,w1-w0,'dsygvd')

   else
!                                                     call timing(t0,w0)
!     go to MO basis using trafo(X) from first iteration (=full diag)
      call blas_gemm('N','N',ndim,ndim,ndim,1.0_wp,H,ndim,X,ndim,0.0_wp,P,ndim)
      call blas_gemm('T','N',ndim,ndim,ndim,1.0_wp,X,ndim,P,ndim,0.0_wp,H,ndim)
!                                                     call timing(t1,w1)
!                       call prtime(6,1.5*(t1-t0),1.5*(w1-w0),'3xdgemm')
!                                                     call timing(t0,w0)
      call pseudodiag(ndim,ihomo,H,e)
!                                                     call timing(t1,w1)
!                                call prtime(6,t1-t0,w1-w0,'pseudodiag')

!     C = X C', P=scratch
      call blas_gemm('N','N',ndim,ndim,ndim,1.0_wp,X,ndim,H,ndim,0.0_wp,P,ndim)
!     save and output MO matrix in AO basis
      H = P
   endif

end subroutine solve


subroutine fermismear(prt,norbs,nel,t,eig,occ,fod,e_fermi,s)
   use xtb_mctc_convert, only : autoev
   use xtb_mctc_constants, only : kB
   integer, intent(in)  :: norbs
   integer, intent(in)  :: nel
   real(wp),intent(in)  :: eig(norbs)
   real(wp),intent(out) :: occ(norbs)
   real(wp),intent(in)  :: t
   real(wp),intent(out) :: fod
   real(wp),intent(out) :: e_fermi
   logical, intent(in)  :: prt

   real(wp) :: bkt,occt,total_number
   real(wp) :: total_dfermi,dfermifunct,fermifunct,s,change_fermi
   real(wp), parameter :: boltz = kB*autoev
   real(wp), parameter :: thr   = 1e-9_wp
   real(wp), parameter :: sqrttiny = sqrt(tiny(1.0_wp))
   integer :: ncycle,i,j,m,k,i1,i2

   bkt = boltz*t


   ! First we need a good guess for the Fermi level
   if(nel+1 .gt. norbs) then 
      ! some atoms (e.g., He) do not have a LUMO because of the valence basis and
      ! the LUMO index becomes larger than No. MOs
      e_fermi = eig(nel)
   else if (nel .eq. 0) then
      ! without electrons the Fermi energy is the energy of the LUMO
      ! i.e. the lowest orbital
      e_fermi = eig(nel+1)
   else
      ! In all other cases the Fermi energy starts as the midpoint between HOMO and LUMO
      e_fermi = 0.5_wp*(eig(nel)+eig(nel+1))

      ! With this we can refine it to meet the definition at the current temperture
      occt=nel
      do ncycle = 1, 200  ! this loop would be possible instead of gotos
         total_number = 0.0_wp
         total_dfermi = 0.0_wp
         do i = 1, norbs
            fermifunct = 0.0_wp
            if((eig(i)-e_fermi)/bkt.lt.50) then
               fermifunct = 1.0_wp/(exp((eig(i)-e_fermi)/bkt)+1.0_wp)
               dfermifunct = exp((eig(i)-e_fermi)/bkt) / &
               &       (bkt*(exp((eig(i)-e_fermi)/bkt)+1.0_wp)**2)
            else
               dfermifunct = 0.0_wp
            end if
            occ(i) = fermifunct
            total_number = total_number + fermifunct
            total_dfermi = total_dfermi + dfermifunct
         end do
         if (total_dfermi > sqrttiny) then
            change_fermi = (occt-total_number)/total_dfermi
         else
            change_fermi = 0.0_wp
         end if
         change_fermi = (occt-total_number)/total_dfermi
         e_fermi = e_fermi+change_fermi
         if (abs(occt-total_number).le.thr) exit
      enddo
   end if 

   fod=0
   s  =0
   do i=1,norbs
      if(occ(i).gt.thr.and.1.0_wp-occ(i).gt.thr) &
      &   s=s+occ(i)*log(occ(i))+(1.0_wp-occ(i))*log(1.0_wp-occ(i))
      if (eig(i).lt.e_fermi) then
         fod=fod+1.0_wp-occ(i)
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
         focc(i)=2.0_wp
      enddo
      if(2*ihomo.ne.nel) then
         ihomo=ihomo+1
         focc(ihomo)=1.0_wp
         if(nopen.eq.0)nopen=1
      endif
      if(nopen.gt.1)then
         do i=1,nopen/2
            focc(ihomo-i+1)=focc(ihomo-i+1)-1.0_wp
            focc(ihomo+i)=focc(ihomo+i)+1.0_wp
         enddo
      endif
!  odd nel
   else
      na=nel/2+(nopen-1)/2+1
      nb=nel/2-(nopen-1)/2
      do i=1,na
         focc(i)=focc(i)+1.0_wp
      enddo
      do i=1,nb
         focc(i)=focc(i)+1.0_wp
      enddo
   endif

   do i=1,ndim
      if(focc(i).gt.0.99_wp) ihomo=i
   enddo

end subroutine occ


subroutine occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
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
         focc(i)=2.0_wp
      enddo
      if(2*ihomo.ne.nel) then
         ihomo=ihomo+1
         focc(ihomo)=1
         if(nopen.eq.0)nopen=1
      endif
      if(nopen.gt.1)then
         do i=1,nopen/2
            focc(ihomo-i+1)=focc(ihomo-i+1)-1.0_wp
            focc(ihomo+i)=focc(ihomo+i)+1.0_wp
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
         focca(i)=1.0_wp
         foccb(i)=1.0_wp
      endif
      if(focc(i).eq.1)focca(i)=1.0_wp
   enddo

   ihomoa=0
   ihomob=0
   do i=1,ndim
      if(focca(i).gt.0.99_wp) ihomoa=i
      if(foccb(i).gt.0.99_wp) ihomob=i
   enddo

end subroutine occu


!> density matrix
! C: MO coefficient
! X: scratch
! P  dmat
subroutine dmat(ndim,focc,C,P)
   integer, intent(in)  :: ndim
   real(wp),intent(in)  :: focc(:)
   real(wp),intent(in)  :: C(:,:)
   real(wp),intent(out) :: P(:,:)
   integer :: i,m
   real(wp),allocatable :: Ptmp(:,:)

   allocate(Ptmp(ndim,ndim))
   ! acc enter data create(Ptmp(:,:)) copyin(C(:, :), focc(:), P(:, :))
   ! acc kernels default(present)
   Ptmp = 0.0_wp
   ! acc end kernels

   ! acc parallel
   ! acc loop gang collapse(2)
   do m=1,ndim
      do i=1,ndim
         Ptmp(i,m)=C(i,m)*focc(m)
      enddo
   enddo
   ! acc end parallel
   ! acc update host(Ptmp)
   call mctc_gemm(C, Ptmp, P, transb='t')
   ! acc exit data copyout(P(:,:)) delete(C(:,:), focc(:), Ptmp(:, :))

   deallocate(Ptmp)

end subroutine dmat

! Reference: I. Mayer, "Simple Theorems, Proofs, and Derivations in Quantum Chemistry", formula (7.35)
subroutine get_wiberg(n,ndim,at,xyz,P,S,wb,fila2)
   integer, intent(in)  :: n,ndim,at(n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: P(ndim,ndim)
   real(wp),intent(in)  :: S(ndim,ndim)
   real(wp),intent(out) :: wb (n,n)
   integer, intent(in)  :: fila2(:,:)

   real(wp),allocatable :: Ptmp(:,:)
   real(wp) xsum,rab
   integer i,j,k,m

   allocate(Ptmp(ndim,ndim))
   call blas_gemm('N','N',ndim,ndim,ndim,1.0d0,P,ndim,S,ndim,0.0d0,Ptmp,ndim)
   wb = 0
   !$omp parallel do default(none) &
   !$omp private(i,j,k,m,xsum,rab) &
   !$omp shared(n,xyz,fila2,Ptmp,wb) &
   !$omp schedule(dynamic,32) collapse(2)
   do i = 1, n
      do j = 1, n
         if (j >= i) cycle
         xsum = 0.0_wp
         rab = sum((xyz(:,i) - xyz(:,j))**2)
         if(rab < 100.0_wp)then
            do k = fila2(1,i), fila2(2,i) ! AOs on atom i
               do m = fila2(1,j), fila2(2,j) ! AOs on atom j
                  xsum = xsum + Ptmp(k,m)*Ptmp(m,k)
               enddo
            enddo
         endif
         wb(i,j) = xsum
         wb(j,i) = xsum
      enddo
   enddo
   deallocate(Ptmp)

end subroutine get_wiberg

! Reference: I. Mayer, "Simple Theorems, Proofs, and Derivations in Quantum Chemistry", formula (7.36)
subroutine get_unrestricted_wiberg(n,ndim,at,xyz,Pa,Pb,S,wb,fila2)
   integer, intent(in)  :: n,ndim,at(n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: Pa(ndim,ndim)
   real(wp),intent(in)  :: Pb(ndim,ndim)
   real(wp),intent(in)  :: S(ndim,ndim)
   real(wp),intent(out) :: wb (n,n)
   integer, intent(in)  :: fila2(:,:)

   real(wp),allocatable :: Ptmp_a(:,:)
   real(wp),allocatable :: Ptmp_b(:,:)
   real(wp) xsum,rab
   integer i,j,k,m

   allocate(Ptmp_a(ndim,ndim))
   allocate(Ptmp_b(ndim,ndim))

   ! P^(alpha) * S !
   call blas_gemm('N','N',ndim,ndim,ndim,1.0_wp,Pa,ndim,S,ndim,0.0_wp,Ptmp_a,ndim)
   
   ! P^(beta) * S !
   call blas_gemm('N','N',ndim,ndim,ndim,1.0_wp,Pb,ndim,S,ndim,0.0_wp,Ptmp_b,ndim)
   
   wb = 0
   !$omp parallel do default(none) &
   !$omp private(i,j,k,m,xsum,rab) &
   !$omp shared(n,xyz,fila2,Ptmp_a,Ptmp_b,wb) &
   !$omp schedule(dynamic,32) collapse(2)
   do i = 1, n
      do j = 1, n
         if (j >= i) cycle
         xsum = 0.0_wp
         rab = sum((xyz(:,i) - xyz(:,j))**2)
         if(rab < 100.0_wp)then
            do k = fila2(1,i), fila2(2,i) ! AOs on atom i
               do m = fila2(1,j), fila2(2,j) ! AOs on atom j
                  xsum = xsum + Ptmp_a(k,m)*Ptmp_a(m,k) + Ptmp_b(k,m)*Ptmp_b(m,k) 
               enddo
            enddo
         endif
         wb(i,j) = 2.0_wp*xsum
         wb(j,i) = 2.0_wp*xsum
      enddo
   enddo
   deallocate(Ptmp_a)
   deallocate(Ptmp_b)

end subroutine get_unrestricted_wiberg

!> Mulliken pop + AO pop
subroutine mpopall(n,nao,aoat,S,P,qao,q)
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


!> Mulliken pop
subroutine mpop0(n,nao,aoat,S,P,q)
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


!> Mulliken AO pop
subroutine mpopao(n,nao,S,P,qao)
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


!> Mulliken pop
subroutine mpop(n,nao,aoat,lao,S,P,q,ql)
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


!> Mulliken pop shell wise
subroutine mpopsh(n,nao,nshell,ao2sh,S,P,qsh)
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


subroutine qsh2qat(ash,qsh,qat)
   integer, intent(in) :: ash(:)
   real(wp), intent(in) :: qsh(:)
   real(wp), intent(out) :: qat(:)

   integer :: iSh

   qat(:) = 0.0_wp
   do iSh = 1, size(qsh)
      qat(ash(iSh)) = qat(ash(iSh)) + qsh(iSh)
   enddo

end subroutine qsh2qat


!> Loewdin pop
subroutine lpop(n,nao,aoat,lao,occ,C,f,q,ql)
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
      if(occ(i).lt.1.0e-8_wp) cycle
      do j=1,nao
         cc=f*C(j,i)*C(j,i)*occ(i)
         jj=aoat(j)
         js=mmm(lao(j))
         q(jj)=q(jj)+cc
         ql(js,jj)=ql(js,jj)+cc
      enddo
   enddo

end subroutine lpop


!> atomic valence shell pops and total atomic energy
subroutine iniqshell(xtbData,n,at,z,nshell,q,qsh,gfn_method)
   type(TxTBData), intent(in) :: xtbData
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   integer, intent(in)  :: nshell
   integer, intent(in)  :: gfn_method
   real(wp),intent(in)  :: z(n)
   real(wp),intent(in)  :: q(n)
   real(wp),intent(out) :: qsh(nshell)
   real(wp) :: zshell
   real(wp) :: ntot,fracz
   integer  :: i,j,k,m,l,ll(0:3),iat,lll,iver
   data ll /1,3,5,7/

   qsh = 0.0_wp

   k=0
   do i=1,n
      iat=at(i)
      ntot=-1.0e-6_wp
      do m=1,xtbData%nShell(iat)
         l=xtbData%hamiltonian%angShell(m,iat)
         k=k+1
         zshell=xtbData%hamiltonian%referenceOcc(m,iat)
         ntot=ntot+zshell
         if(ntot.gt.z(i)) zshell=0.0_wp
         fracz=zshell/z(i)
         qsh(k)=fracz*q(i)
      enddo
   enddo

end subroutine iniqshell


subroutine setzshell(xtbData,n,at,nshell,z,zsh,e,gfn_method)
   type(TxTBData), intent(in) :: xtbData
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
   integer i,j,k,m,l,ll(0:3),iat,lll,iver
   data ll /1,3,5,7/

   k=0
   e=0.0_wp
   do i=1,n
      iat=at(i)
      ntot=-1.0e-6_wp
      do m=1,xtbData%nShell(iat)
         l=xtbData%hamiltonian%angShell(m,iat)
         k=k+1
         zsh(k)=xtbData%hamiltonian%referenceOcc(m,iat)
!         lsh(k)=l
!         ash(k)=i
         ntot=ntot+zsh(k)
         if(ntot.gt.z(i)) zsh(k)=0.0_wp
         e=e+xtbData%hamiltonian%selfEnergy(m,iat)*zsh(k)
      enddo
   enddo

end subroutine setzshell


end module xtb_scc_core
