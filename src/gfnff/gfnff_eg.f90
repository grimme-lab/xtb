! This file is part of xtb.
!
! Copyright (C) 2019-2020 Stefan Grimme
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

module xtb_gfnff_eg
   use xtb_gfnff_ini2
   use xtb_gfnff_data, only : TGFFData
   use xtb_gfnff_neighbourlist, only : TGFFNeighbourList, new
   use xtb_gfnff_topology, only : TGFFTopology
   use xtb_solv_gbsa, only : TBorn
   use xtb_type_environment, only : TEnvironment
   use xtb_mctc_lapack, only : mctc_sytrf, mctc_sytrs
   use xtb_mctc_blas, only : mctc_gemv
   use xtb_param_sqrtzr4r2, only : sqrtZr4r2
   implicit none
   private
   public :: gfnff_eg, gfnff_dlogcoord

contains

!---------------------------------------------------
! GFN-FF
! energy and analytical gradient for given xyz and
! charge ichrg
! requires D3 ini (rcov,r2r4,copyc6) as well as
! gfnff_ini call
!
! the total energy is
! ees + edisp + erep + ebond + eangl + etors + ehb + exb + ebatm + eext
!
! uses EEQ charge and D3 routines
! basic trigonometry for bending and torsion angles
! taken slightly modified from QMDFF code
! repulsion and rabguess from xtb GFN0 part
!
! requires setup of
!     integer,allocatable :: blist(:,:)
!     integer,allocatable :: alist(:,:)
!     integer,allocatable :: tlist(:,:)
!     integer,allocatable ::b3list(:,:)
!     real(wp),allocatable:: vbond(:,:)
!     real(wp),allocatable:: vangl(:,:)
!     real(wp),allocatable:: vtors(:,:)
!     chi,gam,alp,cnf
!     repa,repz,alphanb
!
!---------------------------------------------------

   subroutine gfnff_eg(env,pr,n,ichrg,at,xyz,makeq,g,etot,res_gff, &
         & param,topo,nlist,solvation,update,version,accuracy)
      use xtb_mctc_accuracy, only : wp
      use xtb_gfnff_param, only : efield, gffVersion, gfnff_thresholds
      use xtb_type_data
      use xtb_type_timer
      use xtb_gfnff_gdisp0
      use xtb_mctc_constants
      implicit none
      character(len=*), parameter :: source = 'gfnff_eg'
      type(TEnvironment), intent(inout) :: env
      type(scc_results),intent(out) :: res_gff
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(in) :: topo
      type(TGFFNeighbourList), intent(inout) :: nlist
      type(TBorn), allocatable, intent(inout) :: solvation
      logical, intent(in) :: update
      integer, intent(in) :: version
      real(wp), intent(in) :: accuracy
      integer n
      integer ichrg
      integer at(n)
      real*8 xyz(3,n)
      real*8 g  (3,n)
      real*8 etot
      logical pr
      logical makeq

      real*8 edisp,ees,ebond,eangl,etors,erep,ehb,exb,ebatm,eext
      real*8 :: gsolv, gborn, ghb, gsasa, gshift

      integer i,j,k,l,m,ij,nd3
      integer ati,atj,iat,jat
      integer hbA,hbB,nbk,nbnbk
      integer lin
      logical ex, require_update
      integer nhb1, nhb2, nxb
      real*8  r2,rab,qq0,erff,dd,dum1,r3(3),t8,dum,t22,t39
      real*8  dx,dy,dz,yy,t4,t5,t6,alpha,t20
      real*8  repab,t16,t19,t26,t27,xa,ya,za,cosa,de,t28
      real*8  gammij,eesinf,etmp,phi,valijklff
      real*8  omega,rn,dr,g3tmp(3,3),g4tmp(3,4)
      real*8 rij,drij(3,n)

      real*8, allocatable :: grab0(:,:,:), rab0(:), eeqtmp(:,:)
      real*8, allocatable :: cn(:), dcn(:,:,:), qtmp(:)
      real*8, allocatable :: hb_cn(:), hb_dcn(:,:,:)
      real*8, allocatable :: sqrab(:), srab(:)
      real*8, allocatable :: g5tmp(:,:)
      integer,allocatable :: d3list(:,:)
      type(tb_timer) :: timer
      real(wp) :: dispthr, cnthr, repthr, hbthr1, hbthr2

      call gfnff_thresholds(accuracy, dispthr, cnthr, repthr, hbthr1, hbthr2)

      g  =  0
      exb = 0
      ehb = 0
      erep= 0
      ees = 0
      edisp=0
      ebond=0
      eangl=0
      etors=0
      ebatm=0
      eext =0

      gsolv = 0.0d0
      gsasa = 0.0d0
      gborn = 0.0d0
      ghb = 0.0d0
      gshift = 0.0d0

      allocate(sqrab(n*(n+1)/2),srab(n*(n+1)/2),qtmp(n),g5tmp(3,n), &
     &         eeqtmp(2,n*(n+1)/2),d3list(2,n*(n+1)/2),dcn(3,n,n),cn(n), &
     &         hb_dcn(3,n,n),hb_cn(n))

      if (pr) call timer%new(10 + count([allocated(solvation)]),.false.)

      if (pr) call timer%measure(1,'distance/D3 list')
      nd3=0
      do i=1,n
         ij = i*(i-1)/2
         do j=1,i-1
            k = ij+j
            sqrab(k)=(xyz(1,i)-xyz(1,j))**2+&
     &               (xyz(2,i)-xyz(2,j))**2+&
     &               (xyz(3,i)-xyz(3,j))**2
            if(sqrab(k).lt.dispthr)then
               nd3=nd3+1
               d3list(1,nd3)=i
               d3list(2,nd3)=j
            endif
            srab (k)=sqrt(sqrab(k))
         enddo
! The loop above only runs over the off diagonal elements
! This initializes the unitialized diagonal to zero but does not
! add it to the dispersion list.
         sqrab(ij + i) = 0.0d0
         srab(ij + i) = 0.0d0
      enddo
      if (pr) call timer%measure(1)

!!!!!!!!!!!!
! Setup HB
!!!!!!!!!!!!

      if (pr) call timer%measure(10,'HB/XB (incl list setup)')
      if (allocated(nlist%q)) then
         nlist%initialized = size(nlist%q) == n
      end if
      call gfnff_hbset0(n,at,xyz,sqrab,topo,nhb1,nhb2,nxb,hbthr1,hbthr2)
      nlist%initialized = nlist%initialized .and. nhb1 <= nlist%nhb1 &
         & .and. nhb2 <= nlist%nhb2 .and. nxb <= nlist%nxb
      require_update = .not.nlist%initialized
      if (.not.nlist%initialized) then
         if (pr) then
            write(env%unit,'(10x,"Number of HB bonds (bound hydrogen)",5x,i0,x,i0,x,i0)') &
               & nhb1
            write(env%unit,'(10x,"Number of HB bonds (unbound hydrogen)",3x,i0,x,i0,x,i0)') &
               & nhb2
            write(env%unit,'(10x,"Number of XB bonds",22x,i0,x,i0,x,i0)') &
               & nxb
         end if
         call new(nlist, n, 5*nhb1, 5*nhb2, 3*nxb)
         nlist%hbrefgeo(:, :) = xyz
      end if
      if (update .or. require_update) then
         call gfnff_hbset(n,at,xyz,sqrab,topo,nlist,hbthr1,hbthr2)
      end if
      if (pr) call timer%measure(10)

!!!!!!!!!!!!!
! Setup
! GBSA
!!!!!!!!!!!!!

   if (allocated(solvation)) then
      call timer%measure(11, "GBSA")
      call solvation%update(env, at, xyz)
      call timer%measure(11)
   endif


!!!!!!!!!!!!!
! REP part
! non-bonded
!!!!!!!!!!!!!

      if (pr) call timer%measure(2,'non bonded repulsion')
      !$omp parallel do default(none) reduction(+:erep, g) &
      !$omp shared(n, at, xyz, srab, sqrab, repthr, topo, param) &
      !$omp private(iat, jat, m, ij, ati, atj, rab, r2, r3, t8, t16, t19, t26, t27)
      do iat=1,n
         m=iat*(iat-1)/2
         do jat=1,iat-1
         ij=m+jat
         r2=sqrab(ij)
         if(r2.gt.repthr)   cycle ! cut-off
         if(topo%bpair(ij).eq.1) cycle ! list avoided because of memory
         ati=at(iat)
         atj=at(jat)
         rab=srab(ij)
         t16=r2**0.75
         t19=t16*t16
         t8 =t16*topo%alphanb(ij)
         t26=exp(-t8)*param%repz(ati)*param%repz(atj)*param%repscaln
         erep=erep+t26/rab !energy
         t27=t26*(1.5d0*t8+1.0d0)/t19
         r3 =(xyz(:,iat)-xyz(:,jat))*t27
         g(:,iat)=g(:,iat)-r3
         g(:,jat)=g(:,jat)+r3
         enddo
      enddo
      !$omp end parallel do
      if (pr) call timer%measure(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! just a extremely crude mode for 2D-3D conversion
! i.e. an harmonic potential with estimated Re
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (version == gffVersion%harmonic2020) then
      ebond=0
      !$omp parallel do default(none) reduction(+:ebond, g) &
      !$omp shared(topo, param, xyz, at) private(i, iat, jat, rab, r2, r3, rn, dum)
      do i=1,topo%nbond
         iat=topo%blist(1,i)
         jat=topo%blist(2,i)
         r3 =xyz(:,iat)-xyz(:,jat)
         rab=sqrt(sum(r3*r3))
         rn=0.7*(param%rcov(at(iat))+param%rcov(at(jat)))
         r2=rn-rab
         ebond=ebond+0.1d0*r2**2  ! fixfc = 0.1
         dum=0.1d0*2.0d0*r2/rab
         g(:,jat)=g(:,jat)+dum*r3
         g(:,iat)=g(:,iat)-dum*r3
      enddo
      !$omp end parallel do
      etot = ebond + erep
      return
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! erf CN and gradient for disp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (pr) call timer%measure(3,'dCN')
      call gfnff_dlogcoord(n,at,xyz,srab,cn,dcn,cnthr,param) ! new erf used in GFN0
      if (sum(topo%nr_hb).gt.0) call dncoord_erf(n,at,xyz,param%rcov,hb_cn,hb_dcn,900.0d0,topo) ! HB erf CN
      if (pr) call timer%measure(3)

!!!!!!
! EEQ
!!!!!!

      if (pr) call timer%measure(4,'EEQ energy and q')
      call goed_gfnff(env,accuracy.gt.1,n,at,sqrab,srab,&         ! modified version
     &                dfloat(ichrg),eeqtmp,cn,nlist%q,ees,solvation,param,topo)  ! without dq/dr
      if (pr) call timer%measure(4)

!!!!!!!!
!D3(BJ)
!!!!!!!!

      if (pr) call timer%measure(5,'D3')
      if(nd3.gt.0) then
         call d3_gradient(topo%dispm, n, at, xyz, nd3, d3list, topo%zetac6, &
            & param%d3r0, sqrtZr4r2, 4.0d0, param%dispscale, cn, dcn, edisp, g)
      endif
      deallocate(d3list)
      if (pr) call timer%measure(5)

!!!!!!!!
! ES part
!!!!!!!!
      if (pr) call timer%measure(6,'EEQ gradient')
      !$omp parallel do default(none) reduction (+:g) &
      !$omp shared(topo,nlist,n,sqrab,srab,eeqtmp,xyz,at) &
      !$omp private(i,j,k,ij,r3,r2,rab,gammij,erff,dd)
      do i=1,n
         k = i*(i-1)/2
         do j=1,i-1
            ij = k+j
            r2 =sqrab(ij)
            rab= srab(ij)
            gammij=eeqtmp(1,ij)
            erff  =eeqtmp(2,ij)
            dd=(2.0d0*gammij*exp(-gammij**2*r2) &
               & /(sqrtpi*r2)-erff/(rab*r2))*nlist%q(i)*nlist%q(j)
            r3=(xyz(:,i)-xyz(:,j))*dd
            g(:,i)=g(:,i)+r3
            g(:,j)=g(:,j)-r3
         enddo
      enddo
      !$omp end parallel do
      if(.not.pr) deallocate(eeqtmp)

      if (allocated(solvation)) then
         call timer%measure(11, "GBSA")
         call solvation%addGradient(env, at, xyz, nlist%q, nlist%q, g)
         call solvation%getEnergyParts(env, nlist%q, nlist%q, gborn, ghb, gsasa, &
            & gshift)
         gsolv = gsasa + gborn + ghb + gshift
         call timer%measure(11)
      else
         gborn = 0.0d0
         ghb = 0.0d0
      endif

      do i=1,n
         qtmp(i)=nlist%q(i)*param%cnf(at(i))/(2.0d0*sqrt(cn(i))+1.d-16)
      enddo

      call mctc_gemv(dcn, qtmp, g, alpha=-1.0_wp, beta=1.0_wp)
      if (pr) call timer%measure(6)

!!!!!!!!!!!!!!!!!!
! SRB bonded part
!!!!!!!!!!!!!!!!!!

      if (pr) call timer%measure(7,'bonds')
      if(topo%nbond.gt.0)then
      allocate(grab0(3,n,topo%nbond),rab0(topo%nbond))
      rab0(:)=topo%vbond(1,:) ! shifts
      call gfnffdrab(n,at,xyz,cn,dcn,topo%nbond,topo%blist,rab0,grab0)
      deallocate(dcn)

      !$omp parallel do default(none) reduction(+:g, ebond) &
      !$omp shared(grab0, topo, param, rab0, srab, xyz, at, hb_cn, hb_dcn, n) &
      !$omp private(i, k, iat, jat, ij, rab, rij, drij, t8, dr, dum, yy, &
      !$omp& dx, dy, dz, t4, t5, t6, ati, atj)
      do i=1,topo%nbond
         iat=topo%blist(1,i)
         jat=topo%blist(2,i)
         ati=at(iat)
         atj=at(jat)
         ij=iat*(iat-1)/2+jat
         rab=srab(ij)
         rij=rab0(i)
         drij=grab0(:,:,i)
         if (topo%nr_hb(i).ge.1) then
            call egbond_hb(i,iat,jat,rab,rij,drij,hb_cn,hb_dcn,n,at,xyz,ebond,g,param,topo)
         else
            call egbond(i,iat,jat,rab,rij,drij,n,at,xyz,ebond,g,topo)
         end if
      enddo
      !$omp end parallel do

      deallocate(hb_dcn)

!!!!!!!!!!!!!!!!!!
! bonded REP
!!!!!!!!!!!!!!!!!!

      !$omp parallel do default(none) reduction(+:erep, g) &
      !$omp shared(topo, param, at, sqrab, srab, xyz) &
      !$omp private(i, iat, jat, ij, xa, ya, za, dx, dy, dz, r2, rab, ati, atj, &
      !$omp& alpha, repab, t16, t19, t26, t27)
      do i=1,topo%nbond
         iat=topo%blist(1,i)
         jat=topo%blist(2,i)
         ij=iat*(iat-1)/2+jat
         xa=xyz(1,iat)
         ya=xyz(2,iat)
         za=xyz(3,iat)
         dx=xa-xyz(1,jat)
         dy=ya-xyz(2,jat)
         dz=za-xyz(3,jat)
         r2=sqrab(ij)
         rab=srab(ij)
         ati=at(iat)
         atj=at(jat)
         alpha=sqrt(param%repa(ati)*param%repa(atj))
         repab=param%repz(ati)*param%repz(atj)*param%repscalb
         t16=r2**0.75d0
         t19=t16*t16
         t26=exp(-alpha*t16)*repab
         erep=erep+t26/rab !energy
         t27=t26*(1.5d0*alpha*t16+1.0d0)/t19
         g(1,iat)=g(1,iat)-dx*t27
         g(2,iat)=g(2,iat)-dy*t27
         g(3,iat)=g(3,iat)-dz*t27
         g(1,jat)=g(1,jat)+dx*t27
         g(2,jat)=g(2,jat)+dy*t27
         g(3,jat)=g(3,jat)+dz*t27
      enddo
      !$omp end parallel do
      endif
      if (pr) call timer%measure(7)

!!!!!!!!!!!!!!!!!!
! bend
!!!!!!!!!!!!!!!!!!

      if (pr) call timer%measure(8,'bend and torsion')
      if(topo%nangl.gt.0)then
         !$omp parallel do default(none) reduction (+:eangl, g) &
         !$omp shared(n, at, xyz, topo, param) &
         !$omp private(m, j, i, k, etmp, g3tmp)
         do m=1,topo%nangl
            j = topo%alist(1,m)
            i = topo%alist(2,m)
            k = topo%alist(3,m)
            call egbend(m,j,i,k,n,at,xyz,etmp,g3tmp,param,topo)
            g(1:3,j)=g(1:3,j)+g3tmp(1:3,1)
            g(1:3,i)=g(1:3,i)+g3tmp(1:3,2)
            g(1:3,k)=g(1:3,k)+g3tmp(1:3,3)
            eangl=eangl+etmp
         enddo
         !$omp end parallel do
      endif

!!!!!!!!!!!!!!!!!!
! torsion
!!!!!!!!!!!!!!!!!!

      if(topo%ntors.gt.0)then
         !$omp parallel do default(none) reduction(+:etors, g) &
         !$omp shared(param, topo, n, at, xyz) &
         !$omp private(m, i, j, k, l, etmp, g4tmp)
         do m=1,topo%ntors
            i=topo%tlist(1,m)
            j=topo%tlist(2,m)
            k=topo%tlist(3,m)
            l=topo%tlist(4,m)
            call egtors(m,i,j,k,l,n,at,xyz,etmp,g4tmp,param,topo)
            g(1:3,i)=g(1:3,i)+g4tmp(1:3,1)
            g(1:3,j)=g(1:3,j)+g4tmp(1:3,2)
            g(1:3,k)=g(1:3,k)+g4tmp(1:3,3)
            g(1:3,l)=g(1:3,l)+g4tmp(1:3,4)
            etors=etors+etmp
         enddo
         !$omp end parallel do
      endif
      if (pr) call timer%measure(8)

!!!!!!!!!!!!!!!!!!
! BONDED ATM
!!!!!!!!!!!!!!!!!!

      if (pr) call timer%measure(9,'bonded ATM')
      if(topo%nbatm.gt.0) then
         !$omp parallel do default(none) reduction(+:ebatm, g) &
         !$omp shared(n, at, xyz, srab, sqrab, topo, param) &
         !$omp private(i, j, k, l, etmp, g3tmp)
         do i=1,topo%nbatm
            j=topo%b3list(1,i)
            k=topo%b3list(2,i)
            l=topo%b3list(3,i)
            call batmgfnff_eg(n,j,k,l,at,xyz,topo%qa,sqrab,srab,etmp,g3tmp,param)
            g(1:3,j)=g(1:3,j)+g3tmp(1:3,1)
            g(1:3,k)=g(1:3,k)+g3tmp(1:3,2)
            g(1:3,l)=g(1:3,l)+g3tmp(1:3,3)
            ebatm=ebatm+etmp
         enddo
         !$omp end parallel do
      endif
      if (pr) call timer%measure(9)

!!!!!!!!!!!!!!!!!!
! EHB
!!!!!!!!!!!!!!!!!!

      if (pr) call timer%measure(10,'HB/XB (incl list setup)')

      if(nlist%nhb1.gt.0) then
         !$omp parallel do default(none) reduction(+:ehb, g) &
         !$omp shared(topo, nlist, param, n, at, xyz, sqrab, srab) &
         !$omp private(i, j, k, l, etmp, g3tmp)
         do i=1,nlist%nhb1
            j=nlist%hblist1(1,i)
            k=nlist%hblist1(2,i)
            l=nlist%hblist1(3,i)
            call abhgfnff_eg1(n,j,k,l,at,xyz,topo%qa,sqrab,srab,etmp,g3tmp,param,topo)
            g(1:3,j)=g(1:3,j)+g3tmp(1:3,1)
            g(1:3,k)=g(1:3,k)+g3tmp(1:3,2)
            g(1:3,l)=g(1:3,l)+g3tmp(1:3,3)
            ehb=ehb+etmp
            nlist%hbe1(i)=etmp ! HB energies  
         enddo
         !$omp end parallel do
      endif


      if(nlist%nhb2.gt.0) then
         !$omp parallel do default(none) reduction(+:ehb, g) &
         !$omp shared(topo, nlist, param, n, at, xyz, sqrab, srab) &
         !$omp private(i, j, k, l, nbk, nbnbk, etmp, g5tmp)
         do i=1,nlist%nhb2
            j=nlist%hblist2(1,i)
            k=nlist%hblist2(2,i)
            l=nlist%hblist2(3,i)
            ! prepare variables for check in carbonyl/nitro case
            ! Number neighbors of C/N should be > 1 for carbonyl/nitro
            nbk=topo%nb(1,k) ! index of first neighbor of k
            ! get number neighbors of neighbor of k
            nbnbk=0
            if(nbk.ne.0) nbnbk=topo%nb(20,nbk) 
            !Carbonyl case R-C=O...H_A
            if(at(k).eq.8.and.topo%nb(20,k).eq.1.and.at(topo%nb(1,k)).eq.6 &
                    & .and.nbnbk.gt.1) then
               call abhgfnff_eg3(n,j,k,l,at,xyz,topo%qa,sqrab,srab, &
                  & etmp,g5tmp,param,topo)
               !Nitro case R-N=O...H_A
            else if(at(k).eq.8.and.topo%nb(20,k).eq.1.and.at(topo%nb(1,k)).eq.7 &
                    &  .and.nbnbk.gt.1) then
               call abhgfnff_eg3(n,j,k,l,at,xyz,topo%qa,sqrab,srab, &
                  & etmp,g5tmp,param,topo)
               !N hetero aromat
            else if(at(k).eq.7.and.topo%nb(20,k).eq.2) then
               call abhgfnff_eg2_rnr(n,j,k,l,at,xyz,topo%qa,sqrab,srab, &
                  & etmp,g5tmp,param,topo)
            else
               !Default
               call abhgfnff_eg2new(n,j,k,l,at,xyz,topo%qa,sqrab,srab, &
                  & etmp,g5tmp,param,topo)
            end if
            g=g+g5tmp
            ehb=ehb+etmp
            nlist%hbe2(i)=etmp
         enddo
         !$omp end parallel do
      endif

!!!!!!!!!!!!!!!!!!
! EXB
!!!!!!!!!!!!!!!!!!

      if(nlist%nxb.gt.0) then
         !$omp parallel do default(none) reduction(+:exb, g) &
         !$omp shared(topo, nlist, param, n, at, xyz) &
         !$omp private(i, j, k, l, etmp, g3tmp)
         do i=1,nlist%nxb
            j=nlist%hblist3(1,i)
            k=nlist%hblist3(2,i)
            l=nlist%hblist3(3,i)
            call rbxgfnff_eg(n,j,k,l,at,xyz,topo%qa,etmp,g3tmp,param)
            g(1:3,j)=g(1:3,j)+g3tmp(1:3,1)
            g(1:3,k)=g(1:3,k)+g3tmp(1:3,2)
            g(1:3,l)=g(1:3,l)+g3tmp(1:3,3)
            exb=exb+etmp
            nlist%hbe3(i)=etmp
         enddo
         !$omp end parallel do
      endif
      if (pr) call timer%measure(10)

!!!!!!!!!!!!!!!!!!
! external stuff
!!!!!!!!!!!!!!!!!!

      if(sum(abs(efield)).gt.1d-6)then
         do i=1,n
            r3(:) =-nlist%q(i)*efield(:)
            g(:,i)= g(:,i) + r3(:)
            eext = eext + r3(1)*(xyz(1,i)-topo%xyze0(1,i))+&
     &                    r3(2)*(xyz(2,i)-topo%xyze0(2,i))+&
     &                    r3(3)*(xyz(3,i)-topo%xyze0(3,i))
         enddo
      endif

!!!!!!!!!!!!!!!!!!
! total energy
!!!!!!!!!!!!!!!!!!
      etot = ees + edisp + erep + ebond &
     &           + eangl + etors + ehb + exb + ebatm + eext &
     &           + gsolv

!!!!!!!!!!!!!!!!!!
! printout
!!!!!!!!!!!!!!!!!!
      if (pr) then
        call timer%write(6,'E+G')
        if(abs(sum(nlist%q)-ichrg).gt.1.d-1) then ! check EEQ only once
          write(env%unit,*) nlist%q
          write(env%unit,*) sum(nlist%q),ichrg
          call env%error('EEQ charge constrain error', source)
          return
        endif
        r3 = 0
        do i=1,n
           r3(:) = r3(:)+nlist%q(i)*xyz(:,i)
        enddo

!       just for fit De calc
        sqrab = 1.d+12
        srab = 1.d+6
        cn = 0
!       asymtotically for R=inf, Etot is the SIE contaminted EES
!       which is computed here to get the atomization energy De,n,at(n)
        call goed_gfnff(env,.true.,n,at,sqrab,srab,dfloat(ichrg),eeqtmp,cn,qtmp,eesinf,solvation,param,topo)
        de=-(etot - eesinf)
      endif
!     write resusts to res type
      res_gff%e_total = etot
      res_gff%gnorm   = sqrt(sum(g**2))
      res_gff%e_bond  = ebond
      res_gff%e_angl  = eangl
      res_gff%e_tors  = etors
      res_gff%e_es    = ees
      res_gff%e_rep   = erep
      res_gff%e_disp  = edisp
      res_gff%e_hb    = ehb
      res_gff%e_xb    = exb
      res_gff%e_batm  = ebatm
      res_gff%e_ext   = eext
      res_gff%g_hb    = ghb
      res_gff%g_born  = gborn
      res_gff%g_solv  = gsolv
      res_gff%g_shift = gshift
      res_gff%g_sasa  = gsasa
      call mctc_gemv(xyz, nlist%q, res_gff%dipole)

   end subroutine gfnff_eg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine egbond(i,iat,jat,rab,rij,drij,n,at,xyz,e,g,topo)
      implicit none
      !Dummy
      type(TGFFTopology), intent(in) :: topo
      integer,intent(in)   :: i
      integer,intent(in)   :: n
      integer,intent(in)   :: iat
      integer,intent(in)   :: jat
      integer,intent(in)   :: at(n)
      real*8,intent(in)    :: rab
      real*8,intent(in)    :: rij
      real*8,intent(in)    :: drij(3,n)
      real*8,intent(in)    :: xyz(3,n)
      real*8,intent(inout) :: e
      real*8,intent(inout) :: g(3,n)
      !Stack
      integer j,k
      real*8 dr,dum
      real*8 dx,dy,dz
      real*8 yy
      real*8 t4,t5,t6,t8

         t8 =topo%vbond(2,i)
         dr =rab-rij
         dum=topo%vbond(3,i)*exp(-t8*dr**2)
         e=e+dum                      ! bond energy
         yy=2.0d0*t8*dr*dum
         dx=xyz(1,iat)-xyz(1,jat)
         dy=xyz(2,iat)-xyz(2,jat)
         dz=xyz(3,iat)-xyz(3,jat)
         t4=-yy*(dx/rab-drij(1,iat))
         t5=-yy*(dy/rab-drij(2,iat))
         t6=-yy*(dz/rab-drij(3,iat))
         g(1,iat)=g(1,iat)+t4-drij(1,iat)*yy ! to avoid if in loop below
         g(2,iat)=g(2,iat)+t5-drij(2,iat)*yy
         g(3,iat)=g(3,iat)+t6-drij(3,iat)*yy
         t4=-yy*(-dx/rab-drij(1,jat))
         t5=-yy*(-dy/rab-drij(2,jat))
         t6=-yy*(-dz/rab-drij(3,jat))
         g(1,jat)=g(1,jat)+t4-drij(1,jat)*yy ! to avoid if in loop below
         g(2,jat)=g(2,jat)+t5-drij(2,jat)*yy
         g(3,jat)=g(3,jat)+t6-drij(3,jat)*yy
         do k=1,n !3B gradient
            g(:,k)=g(:,k)+drij(:,k)*yy
         enddo

      end subroutine egbond

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine egbond_hb(i,iat,jat,rab,rij,drij,hb_cn,hb_dcn,n,at,xyz,e,g,param,topo)
      implicit none
      !Dummy
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(in) :: topo
      integer,intent(in)   :: i
      integer,intent(in)   :: n
      integer,intent(in)   :: iat
      integer,intent(in)   :: jat
      integer,intent(in)   :: at(n)
      real*8,intent(in)    :: rab
      real*8,intent(in)    :: rij
      real*8,intent(in)    :: drij(3,n)
      real*8,intent(in)    :: xyz(3,n)
      real*8,intent(in)    :: hb_cn(n)
      real*8,intent(in)    :: hb_dcn(3,n,n)
      real*8,intent(inout) :: e
      real*8,intent(inout) :: g(3,n)
      !Stack
      integer j,k
      integer jA,jH
      integer hbH,hbB,hbA
      real*8 dr,dum
      real*8 dx,dy,dz
      real*8 yy,zz
      real*8 t1,t4,t5,t6,t8

         if (at(iat).eq.1) then
           hbH=iat
           hbA=jat
         else if (at(jat).eq.1) then
           hbH=jat
           hbA=iat
         else
           write(*,'(10x,"No H-atom found in this bond ",i0,1x,i0)') iat,jat
           return
         end if

         t1=1.0-param%vbond_scale
         t8 =(-t1*hb_cn(hbH)+1.0)*topo%vbond(2,i)
         dr =rab-rij
         dum=topo%vbond(3,i)*exp(-t8*dr**2)
         e=e+dum                      ! bond energy
         yy=2.0d0*t8*dr*dum
         dx=xyz(1,iat)-xyz(1,jat)
         dy=xyz(2,iat)-xyz(2,jat)
         dz=xyz(3,iat)-xyz(3,jat)
         t4=-yy*(dx/rab-drij(1,iat))
         t5=-yy*(dy/rab-drij(2,iat))
         t6=-yy*(dz/rab-drij(3,iat))
         g(1,iat)=g(1,iat)+t4-drij(1,iat)*yy ! to avoid if in loop below
         g(2,iat)=g(2,iat)+t5-drij(2,iat)*yy
         g(3,iat)=g(3,iat)+t6-drij(3,iat)*yy
         t4=-yy*(-dx/rab-drij(1,jat))
         t5=-yy*(-dy/rab-drij(2,jat))
         t6=-yy*(-dz/rab-drij(3,jat))
         g(1,jat)=g(1,jat)+t4-drij(1,jat)*yy ! to avoid if in loop below
         g(2,jat)=g(2,jat)+t5-drij(2,jat)*yy
         g(3,jat)=g(3,jat)+t6-drij(3,jat)*yy
         do k=1,n !3B gradient
            g(:,k)=g(:,k)+drij(:,k)*yy
         end do
         zz=dum*topo%vbond(2,i)*dr**2*t1
         do j=1,topo%bond_hb_nr !CN gradient
            jH = topo%bond_hb_AH(2,j)
            jA = topo%bond_hb_AH(1,j)
            if (jH.eq.hbH.and.jA.eq.hbA) then
               g(:,hbH)=g(:,hbH)+hb_dcn(:,hbH,hbH)*zz
               do k=1,topo%bond_hb_Bn(j)
                  hbB = topo%bond_hb_B(k,j)
                  g(:,hbB)=g(:,hbB)-hb_dcn(:,hbB,hbH)*zz
               end do
            end if
         end do

      end subroutine egbond_hb

      subroutine dncoord_erf(nat,at,xyz,rcov,cn,dcn,thr,topo)
      use xtb_mctc_accuracy, only : wp
      implicit none
      !Dummy
      type(TGFFTopology), intent(in) :: topo
      integer,intent(in)   :: nat
      integer,intent(in)   :: at(nat)
      real(wp),intent(in)  :: xyz(3,nat)
      real(wp),intent(in)  :: rcov(:)
      real(wp),intent(out) :: cn(nat)
      real(wp),intent(out) :: dcn(3,nat,nat)
      real(wp),intent(in),optional :: thr
      real(wp) :: cn_thr
      !Stack
      integer  :: i, j
      integer  :: lin,linAH
      integer  :: iat, jat
      integer  :: iA,jA,jH
      integer  :: ati, atj
      real(wp) :: r, r2, rij(3)
      real(wp) :: rcovij
      real(wp) :: dtmp, tmp
      real(wp),parameter :: hlfosqrtpi = 1.0_wp/1.77245385091_wp
      real(wp),parameter :: kn=27.5_wp
      real(wp),parameter :: rcov_scal=1.78

         cn  = 0._wp
         dcn = 0._wp

         do i = 1,topo%bond_hb_nr
            iat = topo%bond_hb_AH(2,i)
            ati = at(iat)
            iA  = topo%bond_hb_AH(1,i)
            do j = 1, topo%bond_hb_Bn(i)
               jat = topo%bond_hb_B(j,i)
               atj = at(jat)
               rij = xyz(:,jat) - xyz(:,iat)
               r2  = sum( rij**2 )
               if (r2.gt.thr) cycle
               r = sqrt(r2)
               rcovij=rcov_scal*(rcov(ati)+rcov(atj))
               tmp = 0.5_wp * (1.0_wp + erf(-kn*(r-rcovij)/rcovij))
               dtmp =-hlfosqrtpi*kn*exp(-kn**2*(r-rcovij)**2/rcovij**2)/rcovij
               cn(iat) = cn(iat) + tmp
               cn(jat) = cn(jat) + tmp
               dcn(:,jat,jat)= dtmp*rij/r + dcn(:,jat,jat)
               dcn(:,iat,jat)= dtmp*rij/r
               dcn(:,jat,iat)=-dtmp*rij/r
               dcn(:,iat,iat)=-dtmp*rij/r + dcn(:,iat,iat)
            end do
         end do

      end subroutine dncoord_erf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine egbend(m,j,i,k,n,at,xyz,e,g,param,topo)
      use xtb_mctc_constants
      implicit none
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(in) :: topo
      integer m,n,at(n)
      integer i,j,k
      real*8 xyz(3,n),g(3,3),e

      real*8  c0,kijk,va(3),vb(3),vc(3),cosa
      real*8  dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
      real*8  term1(3),term2(3),rab2,vab(3),vcb(3),rp
      real*8  rcb2,damp,dampij,damp2ij,dampjk,damp2jk
      real*8  theta,deda(3),vlen,vp(3),et,dij,c1
      real*8  term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
      real*8  ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
      real*8  omega,rij,rijk,phi0,rkl,rjk,dampkl,damp2kl
      real*8  dampjl,damp2jl,valijklff,rn

         c0  =topo%vangl(1,m)
         kijk=topo%vangl(2,m)
         va(1:3) = xyz(1:3,i)
         vb(1:3) = xyz(1:3,j)
         vc(1:3) = xyz(1:3,k)
         call vsub(va,vb,vab,3)
         call vsub(vc,vb,vcb,3)
         rab2 = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
         rcb2 = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
         call crprod(vcb,vab,vp)
         rp = vlen(vp)+1.d-14
         call impsc(vab,vcb,cosa)
         cosa = dble(min(1.0d0,max(-1.0d0,cosa)))
         theta= dacos(cosa)

         call gfnffdampa(at(i),at(j),rab2,dampij,damp2ij,param)
         call gfnffdampa(at(k),at(j),rcb2,dampjk,damp2jk,param)
         damp=dampij*dampjk

         if(pi-c0.lt.1.d-6)then ! linear
         dt  = theta - c0
         ea  = kijk * dt**2
         deddt = 2.d0 * kijk * dt
         else
         ea=kijk*(cosa-cos(c0))**2
         deddt=2.0d0*kijk*sin(theta)*(cos(c0)-cosa)
         endif

         e = ea * damp
         call crprod(vab,vp,deda)
         rmul1 = -deddt / (rab2*rp)
         deda = deda*rmul1
         call crprod(vcb,vp,dedc)
         rmul2 =  deddt / (rcb2*rp)
         dedc = dedc*rmul2
         dedb = deda+dedc
         term1(1:3)=ea*damp2ij*dampjk*vab(1:3)
         term2(1:3)=ea*damp2jk*dampij*vcb(1:3)
         g(1:3,1) = -dedb(1:3)*damp-term1(1:3)-term2(1:3)
         g(1:3,2) =  deda(1:3)*damp+term1(1:3)
         g(1:3,3) =  dedc(1:3)*damp+term2(1:3)

      end subroutine egbend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine egbend_nci_mul(j,i,k,c0,fc,n,at,xyz,e,g)
      use xtb_mctc_constants
      implicit none
      !Dummy
      integer n,at(n)
      integer i,j,k
      real*8  c0,fc
      real*8  xyz(3,n),g(3,3),e
      !Stack
      real*8  kijk,va(3),vb(3),vc(3),cosa
      real*8  dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
      real*8  term1(3),term2(3),rab2,vab(3),vcb(3),rp
      real*8  rcb2,damp,dampij,damp2ij,dampjk,damp2jk
      real*8  theta,deda(3),vlen,vp(3),et,dij,c1
      real*8  term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
      real*8  ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
      real*8  omega,rij,rijk,phi0,rkl,rjk,dampkl,damp2kl
      real*8  dampjl,damp2jl,valijklff,rn

         kijk=fc/(cos(0.0d0)-cos(c0))**2
         va(1:3) = xyz(1:3,i)
         vb(1:3) = xyz(1:3,j)
         vc(1:3) = xyz(1:3,k)
         call vsub(va,vb,vab,3)
         call vsub(vc,vb,vcb,3)
         rab2 = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
         rcb2 = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
         call crprod(vcb,vab,vp)
         rp = vlen(vp)+1.d-14
         call impsc(vab,vcb,cosa)
         cosa = dble(min(1.0d0,max(-1.0d0,cosa)))
         theta= dacos(cosa)

         if(pi-c0.lt.1.d-6)then     ! linear
         dt  = theta - c0
         ea  = kijk * dt**2
         deddt = 2.d0 * kijk * dt
         else
         ea=kijk*(cosa-cos(c0))**2  ! not linear
         deddt=2.0d0*kijk*sin(theta)*(cos(c0)-cosa)
         endif

         e = (1.0d0-ea)
         call crprod(vab,vp,deda)
         rmul1 = -deddt / (rab2*rp)
         deda = deda*rmul1
         call crprod(vcb,vp,dedc)
         rmul2 =  deddt / (rcb2*rp)
         dedc = dedc*rmul2
         dedb = deda+dedc
         g(1:3,1) =  dedb(1:3)
         g(1:3,2) = -deda(1:3)
         g(1:3,3) = -dedc(1:3)

      end subroutine egbend_nci_mul

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine egbend_nci(j,i,k,c0,kijk,n,at,xyz,e,g,param)
      use xtb_mctc_constants
      implicit none
      !Dummy
      type(TGFFData), intent(in) :: param
      integer n,at(n)
      integer i,j,k
      real*8 c0,kijk
      real*8 xyz(3,n),g(3,3),e
      !Stack
      real*8  va(3),vb(3),vc(3),cosa
      real*8  dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
      real*8  term1(3),term2(3),rab2,vab(3),vcb(3),rp
      real*8  rcb2,damp,dampij,damp2ij,dampjk,damp2jk
      real*8  theta,deda(3),vlen,vp(3),et,dij,c1
      real*8  term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
      real*8  ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
      real*8  omega,rij,rijk,phi0,rkl,rjk,dampkl,damp2kl
      real*8  dampjl,damp2jl,valijklff,rn

         va(1:3) = xyz(1:3,i)
         vb(1:3) = xyz(1:3,j)
         vc(1:3) = xyz(1:3,k)
         call vsub(va,vb,vab,3)
         call vsub(vc,vb,vcb,3)
         rab2 = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
         rcb2 = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
         call crprod(vcb,vab,vp)
         rp = vlen(vp)+1.d-14
         call impsc(vab,vcb,cosa)
         cosa = dble(min(1.0d0,max(-1.0d0,cosa)))
         theta= dacos(cosa)

         call gfnffdampa_nci(at(i),at(j),rab2,dampij,damp2ij,param)
         call gfnffdampa_nci(at(k),at(j),rcb2,dampjk,damp2jk,param)
         damp=dampij*dampjk

         if(pi-c0.lt.1.d-6)then ! linear
         dt  = theta - c0
         ea  = kijk * dt**2
         deddt = 2.d0 * kijk * dt
         else
         ea=kijk*(cosa-cos(c0))**2
         deddt=2.0d0*kijk*sin(theta)*(cos(c0)-cosa)
         endif

         e = ea * damp
         call crprod(vab,vp,deda)
         rmul1 = -deddt / (rab2*rp)
         deda = deda*rmul1
         call crprod(vcb,vp,dedc)
         rmul2 =  deddt / (rcb2*rp)
         dedc = dedc*rmul2
         dedb = deda+dedc
         term1(1:3)=ea*damp2ij*dampjk*vab(1:3)
         term2(1:3)=ea*damp2jk*dampij*vcb(1:3)
         g(1:3,1) = -dedb(1:3)*damp-term1(1:3)-term2(1:3)
         g(1:3,2) =  deda(1:3)*damp+term1(1:3)
         g(1:3,3) =  dedc(1:3)*damp+term2(1:3)

      end subroutine egbend_nci

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine egtors(m,i,j,k,l,n,at,xyz,e,g,param,topo)
      use xtb_mctc_constants
      implicit none
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(in) :: topo
      integer m,n,at(n)
      integer i,j,k,l
      real*8 xyz(3,n),g(3,4),e

      real*8  c0,kijk,va(3),vb(3),vc(3),cosa
      real*8  dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
      real*8  term1(3),term2(3),rab2,vab(3),vcb(3),rp
      real*8  rcb2,damp,dampij,damp2ij,dampjk,damp2jk
      real*8  theta,deda(3),vlen,vp(3),et,dij,c1
      real*8  term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
      real*8  ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
      real*8  omega,rij,rijk,phi0,rkl,rjk,dampkl,damp2kl
      real*8  dampjl,damp2jl,valijklff,rn

         rn=dble(topo%tlist(5,m))
         phi0 =topo%vtors(1,m)
         if(topo%tlist(5,m).gt.0)then
         vab(1:3) = xyz(1:3,i)-xyz(1:3,j)
         vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
         vdc(1:3) = xyz(1:3,k)-xyz(1:3,l)
         rij = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
         rjk = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
         rkl = vdc(1)*vdc(1) + vdc(2)*vdc(2) + vdc(3)*vdc(3)
         call gfnffdampt(at(i),at(j),rij,dampij,damp2ij,param)
         call gfnffdampt(at(k),at(j),rjk,dampjk,damp2jk,param)
         call gfnffdampt(at(k),at(l),rkl,dampkl,damp2kl,param)
         damp= dampjk*dampij*dampkl
         phi=valijklff(n,xyz,i,j,k,l)
         call dphidr  (n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)
         dphi1=phi-phi0
         c1=rn*dphi1+pi
         x1cos=cos(c1)
         x1sin=sin(c1)
         et =(1.+x1cos)*topo%vtors(2,m)
         dij=-rn*x1sin*topo%vtors(2,m)*damp
         term1(1:3)=et*damp2ij*dampjk*dampkl*vab(1:3)
         term2(1:3)=et*damp2jk*dampij*dampkl*vcb(1:3)
         term3(1:3)=et*damp2kl*dampij*dampjk*vdc(1:3)
         g(1:3,1)=dij*dda(1:3)+term1
         g(1:3,2)=dij*ddb(1:3)-term1+term2
         g(1:3,3)=dij*ddc(1:3)+term3-term2
         g(1:3,4)=dij*ddd(1:3)-term3
         e=et*damp
         else
         vab(1:3) = xyz(1:3,j)-xyz(1:3,i)
         vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
         vdc(1:3) = xyz(1:3,j)-xyz(1:3,l)
         rij = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
         rjk = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
         rjl = vdc(1)*vdc(1) + vdc(2)*vdc(2) + vdc(3)*vdc(3)
         call gfnffdampt(at(i),at(j),rij,dampij,damp2ij,param)
         call gfnffdampt(at(k),at(j),rjk,dampjk,damp2jk,param)
         call gfnffdampt(at(j),at(l),rjl,dampjl,damp2jl,param)
         damp= dampjk*dampij*dampjl
         phi=omega(n,xyz,i,j,k,l)
         call domegadr(n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)
         if(topo%tlist(5,m).eq.0)then  ! phi0=0 case
         dphi1=phi-phi0
         c1=dphi1+pi
         x1cos=cos(c1)
         x1sin=sin(c1)
         et   =(1.+x1cos)*topo%vtors(2,m)
         dij  =-x1sin*topo%vtors(2,m)*damp
         else                     ! double min at phi0,-phi0
         et =   topo%vtors(2,m)*(cos(phi) -cos(phi0))**2
         dij=2.*topo%vtors(2,m)* sin(phi)*(cos(phi0)-cos(phi))*damp
         endif
         term1(1:3)=et*damp2ij*dampjk*dampjl*vab(1:3)
         term2(1:3)=et*damp2jk*dampij*dampjl*vcb(1:3)
         term3(1:3)=et*damp2jl*dampij*dampjk*vdc(1:3)
         g(1:3,1)=dij*dda(1:3)-term1
         g(1:3,2)=dij*ddb(1:3)+term1+term2+term3
         g(1:3,3)=dij*ddc(1:3)-term2
         g(1:3,4)=dij*ddd(1:3)-term3
         e=et*damp
         endif

      end subroutine egtors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !torsion without distance damping!!! damping is inherint in the HB term
      subroutine egtors_nci_mul(i,j,k,l,rn,phi0,tshift,n,at,xyz,e,g)
      use xtb_mctc_constants
      implicit none
      !Dummy
      integer n,at(n)
      integer i,j,k,l
      integer rn
      real*8 phi0,tshift
      real*8 xyz(3,n),g(3,4),e
      !Stack
      real*8  c0,fc,kijk,va(3),vb(3),vc(3),cosa
      real*8  dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
      real*8  term1(3),term2(3),rab2,vab(3),vcb(3),rp
      real*8  rcb2,damp,dampij,damp2ij,dampjk,damp2jk
      real*8  theta,deda(3),vlen,vp(3),et,dij,c1
      real*8  term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
      real*8  ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
      real*8  omega,rij,rijk,rkl,rjk,dampkl,damp2kl
      real*8  dampjl,damp2jl,valijklff

         fc=(1.0d0-tshift)/2.0d0
         vab(1:3) = xyz(1:3,i)-xyz(1:3,j)
         vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
         vdc(1:3) = xyz(1:3,k)-xyz(1:3,l)
         rij = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
         rjk = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
         rkl = vdc(1)*vdc(1) + vdc(2)*vdc(2) + vdc(3)*vdc(3)
         phi=valijklff(n,xyz,i,j,k,l)
         call dphidr  (n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)
         dphi1=phi-phi0
         c1=rn*dphi1+pi
         x1cos=cos(c1)
         x1sin=sin(c1)
         et =(1.+x1cos)*fc+tshift
         dij=-rn*x1sin*fc
         g(1:3,1)=dij*dda(1:3)
         g(1:3,2)=dij*ddb(1:3)
         g(1:3,3)=dij*ddc(1:3)
         g(1:3,4)=dij*ddd(1:3)
         e=et !*damp
      end subroutine egtors_nci_mul

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine egtors_nci(i,j,k,l,rn,phi0,fc,n,at,xyz,e,g,param)
      use xtb_mctc_constants
      implicit none
      !Dummy
      type(TGFFData), intent(in) :: param
      integer n,at(n)
      integer i,j,k,l
      integer rn
      real*8 phi0,fc
      real*8 xyz(3,n),g(3,4),e
      !Stack
      real*8  c0,kijk,va(3),vb(3),vc(3),cosa
      real*8  dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
      real*8  term1(3),term2(3),rab2,vab(3),vcb(3),rp
      real*8  rcb2,damp,dampij,damp2ij,dampjk,damp2jk
      real*8  theta,deda(3),vlen,vp(3),et,dij,c1
      real*8  term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
      real*8  ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
      real*8  omega,rij,rijk,rkl,rjk,dampkl,damp2kl
      real*8  dampjl,damp2jl,valijklff

         vab(1:3) = xyz(1:3,i)-xyz(1:3,j)
         vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
         vdc(1:3) = xyz(1:3,k)-xyz(1:3,l)
         rij = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
         rjk = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
         rkl = vdc(1)*vdc(1) + vdc(2)*vdc(2) + vdc(3)*vdc(3)
         call gfnffdampt_nci(at(i),at(j),rij,dampij,damp2ij,param)
         call gfnffdampt_nci(at(k),at(j),rjk,dampjk,damp2jk,param)
         call gfnffdampt_nci(at(k),at(l),rkl,dampkl,damp2kl,param)
         damp= dampjk*dampij*dampkl
         phi=valijklff(n,xyz,i,j,k,l)
         call dphidr  (n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)
         dphi1=phi-phi0
         c1=rn*dphi1+pi
         x1cos=cos(c1)
         x1sin=sin(c1)
         et =(1.+x1cos)*fc
         dij=-rn*x1sin*fc*damp
         term1(1:3)=et*damp2ij*dampjk*dampkl*vab(1:3)
         term2(1:3)=et*damp2jk*dampij*dampkl*vcb(1:3)
         term3(1:3)=et*damp2kl*dampij*dampjk*vdc(1:3)
         g(1:3,1)=dij*dda(1:3)+term1
         g(1:3,2)=dij*ddb(1:3)-term1+term2
         g(1:3,3)=dij*ddc(1:3)+term3-term2
         g(1:3,4)=dij*ddd(1:3)-term3
         e=et*damp
      end subroutine egtors_nci

!cccccccccccccccccccccccccccccccccccccccccccccc
! damping of bend and torsion for long
! bond distances to allow proper dissociation
!cccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gfnffdampa(ati,atj,r2,damp,ddamp,param)
      implicit none
      type(TGFFData), intent(in) :: param
      integer ati,atj
      real*8 r2,damp,ddamp,rr,rcut
      rcut =param%atcuta*(param%rcov(ati)+param%rcov(atj))**2
      rr   =(r2/rcut)**2
      damp = 1.0d0/(1.0d0+rr)
      ddamp=-2.d0*2*rr/(r2*(1.0d0+rr)**2)
      end subroutine gfnffdampa

      subroutine gfnffdampt(ati,atj,r2,damp,ddamp,param)
      implicit none
      type(TGFFData), intent(in) :: param
      integer ati,atj
      real*8 r2,damp,ddamp,rr,rcut
      rcut =param%atcutt*(param%rcov(ati)+param%rcov(atj))**2
      rr   =(r2/rcut)**2
      damp = 1.0d0/(1.0d0+rr)
      ddamp=-2.d0*2*rr/(r2*(1.0d0+rr)**2)
      end subroutine gfnffdampt

      subroutine gfnffdampa_nci(ati,atj,r2,damp,ddamp,param)
      implicit none
      type(TGFFData), intent(in) :: param
      integer ati,atj
      real*8 r2,damp,ddamp,rr,rcut
      rcut =param%atcuta_nci*(param%rcov(ati)+param%rcov(atj))**2
      rr   =(r2/rcut)**2
      damp = 1.0d0/(1.0d0+rr)
      ddamp=-2.d0*2*rr/(r2*(1.0d0+rr)**2)
      end subroutine gfnffdampa_nci

      subroutine gfnffdampt_nci(ati,atj,r2,damp,ddamp,param)
      implicit none
      type(TGFFData), intent(in) :: param
      integer ati,atj
      real*8 r2,damp,ddamp,rr,rcut
      rcut =param%atcutt_nci*(param%rcov(ati)+param%rcov(atj))**2
      rr   =(r2/rcut)**2
      damp = 1.0d0/(1.0d0+rr)
      ddamp=-2.d0*2*rr/(r2*(1.0d0+rr)**2)
      end subroutine gfnffdampt_nci

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ref.: S. Alireza Ghasemi, Albert Hofstetter, Santanu Saha, and Stefan Goedecker
!       PHYSICAL REVIEW B 92, 045131 (2015)
!       Interatomic potentials for ionic systems with density functional accuracy
!       based on charge densities obtained by a neural network
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine goed_gfnff(env,single,n,at,sqrab,r,chrg,eeqtmp,cn,q,es,gbsa,param,topo)
      use xtb_mctc_accuracy, only : wp, sp
      use xtb_mctc_la
      implicit none
      character(len=*), parameter :: source = 'gfnff_eg_goed'
      type(TEnvironment), intent(inout) :: env
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(in) :: topo
      logical, intent(in)  :: single     ! real*4 flag for solver
      integer, intent(in)  :: n          ! number of atoms
      integer, intent(in)  :: at(n)      ! ordinal numbers
      real(wp),intent(in)  :: sqrab(n*(n+1)/2)   ! squared dist
      real(wp),intent(in)  :: r(n*(n+1)/2)       ! dist
      real(wp),intent(in)  :: chrg       ! total charge on system
      real(wp),intent(in)  :: cn(n)      ! CN
      real(wp),intent(out) :: q(n)       ! output charges
      real(wp),intent(out) :: es         ! ES energy
      real(wp),intent(out) :: eeqtmp(2,n*(n+1)/2)    ! intermediates
      type(TBorn), allocatable, intent(in) :: gbsa

!  local variables
      integer  :: m,i,j,k,ii,ij
      integer,allocatable :: ipiv(:)
      real(wp) :: gammij,tsqrt2pi,r2,tmp
      real(wp),allocatable :: A (:,:),x (:)
      real(sp),allocatable :: A4(:,:),x4(:)
!  parameter
      parameter (tsqrt2pi = 0.797884560802866_wp)
      logical :: exitRun

      m=n+topo%nfrag ! # atoms + chrg constrain + frag constrain

      allocate(A(m,m),x(m))
!  setup RHS
      do i=1,n
         x(i) = topo%chieeq(i) + param%cnf(at(i))*sqrt(cn(i))
      enddo

      A = 0
!  setup A matrix
!$omp parallel default(none) &
!$omp shared(topo,n,sqrab,r,eeqtmp,A,at) &
!$omp private(i,j,k,ij,gammij,tmp)
!$omp do schedule(dynamic)
      do i=1,n
      A(i,i)=tsqrt2pi/sqrt(topo%alpeeq(i))+topo%gameeq(i) ! J of i
      k = i*(i-1)/2
      do j=1,i-1
         ij = k+j
         gammij=1./sqrt(topo%alpeeq(i)+topo%alpeeq(j)) ! squared above
         tmp = erf(gammij*r(ij))
         eeqtmp(1,ij)=gammij
         eeqtmp(2,ij)=tmp
         A(j,i) = tmp/r(ij)
         A(i,j) = A(j,i)
      enddo
      enddo
!$omp enddo
!$omp end parallel

!  fragment charge constrain
      do i=1,topo%nfrag
        x(n+i)=topo%qfrag(i)
        do j=1,n
         if(topo%fraglist(j).eq.i) then
            A(n+i,j)=1
            A(j,n+i)=1
         endif
        enddo
      enddo

      if (allocated(gbsa)) then
         A(:n, :n) = A(:n, :n) + gbsa%bornMat(:, :)
      end if

!     call prmat(6,A,m,m,'A eg')
      allocate(ipiv(m))

      if(single) then
         allocate(A4(m,m),x4(m))
         A4=A
         x4=x
         deallocate(A,x)
         call mctc_sytrf(env, a4, ipiv)
         call mctc_sytrs(env, a4, x4, ipiv)
         q(1:n)=x4(1:n)
         deallocate(A4,x4)
      else
         call mctc_sytrf(env, a, ipiv)
         call mctc_sytrs(env, a, x, ipiv)
         q(1:n)=x(1:n)
         deallocate(A,x)
      endif

      call env%check(exitRun)
      if(exitRun) then
         call env%error('Solving linear equations failed', source)
         return
      end if

      if(n.eq.1) q(1)=chrg

!  energy
      es = 0.0_wp
      do i=1,n
      ii = i*(i-1)/2
      do j=1,i-1
         ij = ii+j
         tmp   =eeqtmp(2,ij)
         es = es + q(i)*q(j)*tmp/r(ij)
      enddo
      es = es - q(i)*(topo%chieeq(i) + param%cnf(at(i))*sqrt(cn(i))) &
     &        + q(i)*q(i)*0.5d0*(topo%gameeq(i)+tsqrt2pi/sqrt(topo%alpeeq(i)))
      enddo

      !work = x
      !call dsymv('u', n, 0.5d0, A, m, q, 1, -1.0_wp, work, 1)
      !es = ddot(n, q, 1, work, 1)

!     deallocate(cn)

      end subroutine goed_gfnff

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! HB energy and analytical gradient
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!subroutine for case 1: A...H...B
subroutine abhgfnff_eg1(n,A,B,H,at,xyz,q,sqrab,srab,energy,gdr,param,topo)
      implicit none
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(in) :: topo
      integer A,B,H,n,at(n)
      real*8 xyz(3,n),energy,gdr(3,3)
      real*8 q(n)
      real*8 sqrab(n*(n+1)/2)   ! squared dist
      real*8 srab(n*(n+1)/2)    ! dist

      real*8 outl,dampl,damps,rdamp,damp,dd24a,dd24b
      real*8 ratio1,ratio2,ratio3
      real*8 xm,ym,zm
      real*8 rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
      real*8 drah(3),drbh(3),drab(3),drm(3)
      real*8 dg(3),dga(3),dgb(3),dgh(3)
      real*8 ga(3),gb(3),gh(3)
      real*8 gi,denom,ratio,tmp,qhoutl,radab,rahprbh
      real*8 ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo
      real*8 bas,aci
      real*8 eabh
      real*8 aterm,rterm,dterm,sterm
      real*8 qa,qb,qh
      real*8 ca(2),cb(2)
      real*8 gqa,gqb,gqh
      real*8 caa,cbb
      real*8 shortcut

      integer i,j,ij,lina
      lina(i,j)=min(i,j)+max(i,j)*(max(i,j)-1)/2

      gdr  = 0
      energy=0

      call hbonds(A,B,ca,cb,param,topo)

!     A-B distance
      ij=lina(A,B)
      rab2=sqrab(ij)
      rab =srab (ij)

!     A-H distance
      ij=lina(A,H)
      rah2= sqrab(ij)
      rah = srab (ij)

!     B-H distance
      ij=lina(B,H)
      rbh2= sqrab(ij)
      rbh = srab (ij)

      rahprbh=rah+rbh+1.d-12
      radab=param%rad(at(A))+param%rad(at(B))

!     out-of-line damp
      expo=(param%hbacut/radab)*(rahprbh/rab-1.d0)
      if(expo.gt.15.0d0) return ! avoid overflow
      ratio2=exp(expo)
      outl=2.d0/(1.d0+ratio2)

!     long damping
      ratio1=(rab2/param%hblongcut)**param%hbalp
      dampl=1.d0/(1.d0+ratio1)

!     short damping
      shortcut=param%hbscut*radab
      ratio3=(shortcut/rab2)**param%hbalp
      damps=1.d0/(1.d0+ratio3)

      damp = damps*dampl
      rdamp = damp/rab2/rab

!     hydrogen charge scaled term
      ex1h=exp(param%hbst*q(H))
      ex2h=ex1h+param%hbsf
      qh=ex1h/ex2h

!     hydrogen charge scaled term
      ex1a=exp(-param%hbst*q(A))
      ex2a=ex1a+param%hbsf
      qa=ex1a/ex2a

!     hydrogen charge scaled term
      ex1b=exp(-param%hbst*q(B))
      ex2b=ex1b+param%hbsf
      qb=ex1b/ex2b

!     donor-acceptor term
      rah4 = rah2*rah2
      rbh4 = rbh2*rbh2
      denom = 1.d0/(rah4+rbh4)

      caa=qa*ca(1)
      cbb=qb*cb(1)
      qhoutl=qh*outl

      bas = (caa*rah4 + cbb*rbh4)*denom
      aci = (cb(2)*rah4+ca(2)*rbh4)*denom

!     energy
      rterm  = -aci*rdamp*qhoutl
      energy =  bas*rterm

!     gradient
      drah(1:3)=xyz(1:3,A)-xyz(1:3,H)
      drbh(1:3)=xyz(1:3,B)-xyz(1:3,H)
      drab(1:3)=xyz(1:3,A)-xyz(1:3,B)

      aterm= -aci*bas*rdamp*qh
      sterm= -rdamp*bas*qhoutl
      dterm= -aci*bas*qhoutl

      tmp=denom*denom*4.0d0
      dd24a=rah2*rbh4*tmp
      dd24b=rbh2*rah4*tmp

!     donor-acceptor part: bas
      gi = (caa-cbb)*dd24a*rterm
      ga(1:3) = gi*drah(1:3)
      gi = (cbb-caa)*dd24b*rterm
      gb(1:3) = gi*drbh(1:3)
      gh(1:3) = -ga(1:3)-gb(1:3)

!     donor-acceptor part: aci
      gi = (cb(2)-ca(2))*dd24a
      dga(1:3) = gi*drah(1:3)*sterm
      ga(1:3) = ga(1:3) + dga(1:3)

      gi = (ca(2)-cb(2))*dd24b
      dgb(1:3) = gi*drbh(1:3)*sterm
      gb(1:3) = gb(1:3) + dgb(1:3)

      dgh(1:3) = -dga(1:3)-dgb(1:3)
      gh(1:3) = gh(1:3) + dgh(1:3)

!     damping part rab
      gi = rdamp*(-(2.d0*param%hbalp*ratio1/(1+ratio1))+(2.d0*param%hbalp*ratio3/(1+ratio3))-3.d0)/rab2
      dg(1:3) = gi*drab(1:3)*dterm
      ga(1:3) = ga(1:3) + dg(1:3)
      gb(1:3) = gb(1:3) - dg(1:3)

!     out of line term: rab
      gi = aterm*2.d0*ratio2*expo*rahprbh/(1+ratio2)**2/(rahprbh-rab)/rab2
      dg(1:3) = gi*drab(1:3)
      ga(1:3) = ga(1:3) + dg(1:3)
      gb(1:3) = gb(1:3) - dg(1:3)

!     out of line term: rah,rbh
      tmp= -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
      dga(1:3) = drah(1:3)*tmp/rah
      ga(1:3) = ga(1:3) + dga(1:3)
      dgb(1:3) = drbh(1:3)*tmp/rbh
      gb(1:3) = gb(1:3) + dgb(1:3)
      dgh(1:3) = -dga(1:3)-dgb(1:3)
      gh(1:3) = gh(1:3) + dgh(1:3)

!     move gradients into place
      gdr(1:3,1) = ga(1:3)
      gdr(1:3,2) = gb(1:3)
      gdr(1:3,3) = gh(1:3)

end subroutine abhgfnff_eg1

!subroutine for case 2: A-H...B including orientation of neighbors at B
subroutine abhgfnff_eg2new(n,A,B,H,at,xyz,q,sqrab,srab,energy,gdr,param,topo)
      implicit none
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(in) :: topo
      integer A,B,H,n,at(n)
      real*8 xyz(3,n),energy,gdr(3,n)
      real*8 q(n)
      real*8 sqrab(n*(n+1)/2)   ! squared dist
      real*8 srab(n*(n+1)/2)    ! dist

      real*8 outl,dampl,damps,rdamp,damp
      real*8 ddamp,rabdamp,rbhdamp
      real*8 ratio1,ratio2,ratio2_nb(topo%nb(20,B)),ratio3
      real*8 xm,ym,zm
      real*8 rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
      real*8 ranb(topo%nb(20,B)),ranb2(topo%nb(20,B)),rbnb(topo%nb(20,B)),rbnb2(topo%nb(20,B))
      real*8 drah(3),drbh(3),drab(3),drm(3)
      real*8 dranb(3,topo%nb(20,B)),drbnb(3,topo%nb(20,B))
      real*8 dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
      real*8 ga(3),gb(3),gh(3),gnb(3,topo%nb(20,B))
      real*8 denom,ratio,qhoutl,radab
      real*8 gi,gi_nb(topo%nb(20,B))
      real*8 tmp1,tmp2(topo%nb(20,B))
      real*8 rahprbh,ranbprbnb(topo%nb(20,B))
      real*8 ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_nb(topo%nb(20,B))
      real*8 eabh
      real*8 aterm,dterm,nbterm
      real*8 qa,qb,qh
      real*8 ca(2),cb(2)
      real*8 gqa,gqb,gqh
      real*8 shortcut
      real*8 const
      real*8 outl_nb(topo%nb(20,B)),outl_nb_tot
      real*8 hbnbcut_save
      logical mask_nb(topo%nb(20,B))

!     proportion between Rbh und Rab distance dependencies
      real*8 :: p_bh
      real*8 :: p_ab

      integer i,j,ij,lina,nbb
      lina(i,j)=min(i,j)+max(i,j)*(max(i,j)-1)/2

      p_bh=1.d0+param%hbabmix
      p_ab=    -param%hbabmix

      gdr    = 0
      energy = 0

      call hbonds(A,B,ca,cb,param,topo)

      nbb=topo%nb(20,B)
!     Neighbours of B
      do i=1,nbb
!        compute distances
         dranb(1:3,i) = xyz(1:3,A) - xyz(1:3,topo%nb(i,B))
         drbnb(1:3,i) = xyz(1:3,B) - xyz(1:3,topo%nb(i,B))
!        A-nb(B) distance
         ranb2(i) = sum(dranb(1:3,i)**2)
         ranb(i)  = sqrt(ranb2(i))
!        B-nb(B) distance
         rbnb2(i) = sum(drbnb(1:3,i)**2)
         rbnb(i)  = sqrt(rbnb2(i))
      end do

!     A-B distance
      ij=lina(A,B)
      rab2=sqrab(ij)
      rab =srab (ij)

!     A-H distance
      ij=lina(A,H)
      rah2= sqrab(ij)
      rah = srab (ij)

!     B-H distance
      ij=lina(B,H)
      rbh2= sqrab(ij)
      rbh = srab (ij)

      rahprbh=rah+rbh+1.d-12
      radab=param%rad(at(A))+param%rad(at(B))

!     out-of-line damp: A-H...B
      expo=(param%hbacut/radab)*(rahprbh/rab-1.d0)
      if(expo.gt.15.0d0) return ! avoid overflow
      ratio2=exp(expo)
      outl=2.d0/(1.d0+ratio2)

!     out-of-line damp: A...nb(B)-B
      if(at(B).eq.7.and.topo%nb(20,B).eq.1) then
        hbnbcut_save = 2.0
      else
        hbnbcut_save = param%hbnbcut
      end if
      do i=1,nbb
         ranbprbnb(i)=ranb(i)+rbnb(i)+1.d-12
         expo_nb(i)=(hbnbcut_save/radab)*(ranbprbnb(i)/rab-1.d0)
         ratio2_nb(i)=exp(-expo_nb(i))**(1.0)
         outl_nb(i)=( 2.d0/(1.d0+ratio2_nb(i)) ) - 1.0d0
      end do
      outl_nb_tot = product(outl_nb)

!     long damping
      ratio1=(rab2/param%hblongcut)**param%hbalp
      dampl=1.d0/(1.d0+ratio1)

!     short damping
      shortcut=param%hbscut*radab
      ratio3=(shortcut/rab2)**param%hbalp
      damps=1.d0/(1.d0+ratio3)

      damp  = damps*dampl
      ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
      rbhdamp = damp * ( (p_bh/rbh2/rbh) )
      rabdamp = damp * ( (p_ab/rab2/rab) )
      rdamp   = rbhdamp + rabdamp

!     hydrogen charge scaled term
      ex1h=exp(param%hbst*q(H))
      ex2h=ex1h+param%hbsf
      qh=ex1h/ex2h

!     hydrogen charge scaled term
      ex1a=exp(-param%hbst*q(A))
      ex2a=ex1a+param%hbsf
      qa=ex1a/ex2a

!     hydrogen charge scaled term
      ex1b=exp(-param%hbst*q(B))
      ex2b=ex1b+param%hbsf
      qb=ex1b/ex2b

      qhoutl=qh*outl*outl_nb_tot

!     constant values, no gradient
      const = ca(2)*qa*cb(1)*qb*param%xhaci_globabh

!     energy
      energy = -rdamp*qhoutl*const

!     gradient
      drah(1:3)=xyz(1:3,A)-xyz(1:3,H)
      drbh(1:3)=xyz(1:3,B)-xyz(1:3,H)
      drab(1:3)=xyz(1:3,A)-xyz(1:3,B)

      aterm  = -rdamp*qh*outl_nb_tot*const
      nbterm = -rdamp*qh*outl*const
      dterm  = -qhoutl*const

!------------------------------------------------------------------------------
!     damping part: rab
      gi = ( (rabdamp+rbhdamp)*ddamp-3.d0*rabdamp ) / rab2
      gi = gi*dterm
      dg(1:3) = gi*drab(1:3)
      ga(1:3) =  dg(1:3)
      gb(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!     damping part: rbh
      gi = -3.d0*rbhdamp/rbh2
      gi = gi*dterm
      dg(1:3) = gi*drbh(1:3)
      gb(1:3) = gb(1:3)+dg(1:3)
      gh(1:3) =       - dg(1:3)

!------------------------------------------------------------------------------
!     angular A-H...B term
!------------------------------------------------------------------------------
!     out of line term: rab
      tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
      gi   = -tmp1 *rahprbh/rab2
      dg(1:3) = gi*drab(1:3)
      ga(1:3) = ga(1:3) + dg(1:3)
      gb(1:3) = gb(1:3) - dg(1:3)

!     out of line term: rah,rbh
      gi = tmp1/rah
      dga(1:3) = gi*drah(1:3)
      ga(1:3)  = ga(1:3) + dga(1:3)
      gi = tmp1/rbh
      dgb(1:3) = gi*drbh(1:3)
      gb(1:3)  = gb(1:3) + dgb(1:3)
      dgh(1:3) = -dga(1:3)-dgb(1:3)
      gh(1:3)  =  gh(1:3) + dgh(1:3)

!------------------------------------------------------------------------------
!     angular A...nb(B)-B term
!------------------------------------------------------------------------------
!     out of line term: rab
      mask_nb=.true.
      do i=1,nbb
         mask_nb(i)=.false.
         tmp2(i)  = 2.d0*nbterm*product(outl_nb,mask_nb)*ratio2_nb(i)*expo_nb(i)/&
                  & (1+ratio2_nb(i))**2/(ranbprbnb(i)-rab)
         gi_nb(i) = -tmp2(i) *ranbprbnb(i)/rab2
         dg(1:3)  = gi_nb(i)*drab(1:3)
         ga(1:3)  = ga(1:3) + dg(1:3)
         gb(1:3)  = gb(1:3) - dg(1:3)
         mask_nb=.true.
      end do

!     out of line term: ranb,rbnb
      do i=1,nbb
         gi_nb(i)   = tmp2(i)/ranb(i)
         dga(1:3)   = gi_nb(i)*dranb(1:3,i)
         ga(1:3)    = ga(1:3) + dga(1:3)
         gi_nb(i)   = tmp2(i)/rbnb(i)
         dgb(1:3)   = gi_nb(i)*drbnb(1:3,i)
         gb(1:3)    = gb(1:3) + dgb(1:3)
         dgnb(1:3)  = -dga(1:3)-dgb(1:3)
         gnb(1:3,i) = dgnb(1:3)
      end do

!------------------------------------------------------------------------------
      if(nbb.lt.1) then
         gdr(1:3,A) = gdr(1:3,A) + ga(1:3)
         gdr(1:3,B) = gdr(1:3,B) + gb(1:3)
         gdr(1:3,H) = gdr(1:3,H) + gh(1:3)
         return
      endif

!------------------------------------------------------------------------------
!     move gradients into place
      gdr(1:3,A) = gdr(1:3,A) + ga(1:3)
      gdr(1:3,B) = gdr(1:3,B) + gb(1:3)
      gdr(1:3,H) = gdr(1:3,H) + gh(1:3)
      do i=1,nbb
         gdr(1:3,topo%nb(i,B)) = gdr(1:3,topo%nb(i,B)) + gnb(1:3,i)
      end do

end subroutine abhgfnff_eg2new

!subroutine for case 2: A-H...B including LP position
subroutine abhgfnff_eg2_rnr(n,A,B,H,at,xyz,q,sqrab,srab,energy,gdr,param,topo)
      implicit none
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(in) :: topo
      integer A,B,H,n,at(n)
      real*8 xyz(3,n),energy,gdr(3,n)
      real*8 q(n)
      real*8 sqrab(n*(n+1)/2)   ! squared dist
      real*8 srab(n*(n+1)/2)    ! dist

      real*8 outl,dampl,damps,rdamp,damp
      real*8 ddamp,rabdamp,rbhdamp
      real*8 ratio1,ratio2,ratio2_lp,ratio2_nb(topo%nb(20,B)),ratio3
      real*8 xm,ym,zm
      real*8 rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
      real*8 ranb(topo%nb(20,B)),ranb2(topo%nb(20,B)),rbnb(topo%nb(20,B)),rbnb2(topo%nb(20,B))
      real*8 drah(3),drbh(3),drab(3),drm(3),dralp(3),drblp(3)
      real*8 dranb(3,topo%nb(20,B)),drbnb(3,topo%nb(20,B))
      real*8 dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
      real*8 ga(3),gb(3),gh(3),gnb(3,topo%nb(20,B)),gnb_lp(3),glp(3)
      real*8 denom,ratio,qhoutl,radab
      real*8 gi,gi_nb(topo%nb(20,B))
      real*8 tmp1,tmp2(topo%nb(20,B)),tmp3
      real*8 rahprbh,ranbprbnb(topo%nb(20,B))
      real*8 ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_lp,expo_nb(topo%nb(20,B))
      real*8 eabh
      real*8 aterm,dterm,nbterm,lpterm
      real*8 qa,qb,qh
      real*8 ca(2),cb(2)
      real*8 gqa,gqb,gqh
      real*8 shortcut
      real*8 const
      real*8 outl_nb(topo%nb(20,B)),outl_nb_tot,outl_lp
      real*8 vector(3),vnorm
      real*8 gii(3,3)
      real*8 unit_vec(3)
      real*8 drnb(3,topo%nb(20,B))
      real*8 lp(3)   !lonepair position
      real*8 lp_dist !distance parameter between B and lonepair
      real*8 ralp,ralp2,rblp,rblp2,ralpprblp
      logical mask_nb(topo%nb(20,B))

!     proportion between Rbh und Rab distance dependencies
      real*8 :: p_bh
      real*8 :: p_ab
!     lone-pair out-of-line damping
      real*8 hblpcut

      integer i,j,ij,lina,nbb
      lina(i,j)=min(i,j)+max(i,j)*(max(i,j)-1)/2

      p_bh=1.d0+param%hbabmix
      p_ab=    -param%hbabmix

      gdr    = 0
      energy = 0
      vector = 0
      lp_dist = 0.50-0.018*param%repz(at(B))
      hblpcut=56

      call hbonds(A,B,ca,cb,param,topo)

      nbb=topo%nb(20,B)
!     Neighbours of B
      do i=1,nbb
!        compute distances
         dranb(1:3,i) = xyz(1:3,A) - xyz(1:3,topo%nb(i,B))
         drbnb(1:3,i) = xyz(1:3,B) - xyz(1:3,topo%nb(i,B))
!        A-nb(B) distance
         ranb2(i) = sum(dranb(1:3,i)**2)
         ranb(i)  = sqrt(ranb2(i))
!        B-nb(B) distance
         rbnb2(i) = sum(drbnb(1:3,i)**2)
         rbnb(i)  = sqrt(rbnb2(i))
      end do

!     Neighbours of B
      do i=1,nbb
         drnb(1:3,i)=xyz(1:3,topo%nb(i,B))-xyz(1:3,B)
         vector = vector + drnb(1:3,i)
      end do

      vnorm = norm2(vector)
!     lonepair coordinates
      if(vnorm.gt.1.d-10) then
      lp = xyz(1:3,B) - lp_dist * ( vector / vnorm )
      else
      lp = xyz(1:3,B)
      nbb = 0
      endif

!     A-B distance
      ij=lina(A,B)
      rab2=sqrab(ij)
      rab =srab (ij)

!     A-H distance
      ij=lina(A,H)
      rah2= sqrab(ij)
      rah = srab (ij)

!     B-H distance
      ij=lina(B,H)
      rbh2= sqrab(ij)
      rbh = srab (ij)

      rahprbh=rah+rbh+1.d-12
      radab=param%rad(at(A))+param%rad(at(B))

!     out-of-line damp: A-H...B
      expo=(param%hbacut/radab)*(rahprbh/rab-1.d0)
      if(expo.gt.15.0d0) return ! avoid overflow
      ratio2=exp(expo)
      outl=2.d0/(1.d0+ratio2)

!     out-of-line damp: A...LP-B
      rblp2 = sum( ( xyz(1:3,B) - lp(1:3) )**2 )
      rblp  = sqrt(rblp2)
      ralp2 = sum( ( xyz(1:3,A) - lp(1:3) )**2 )
      ralp  = sqrt(ralp2)
      ralpprblp=ralp+rblp+1.d-12
      expo_lp=(hblpcut/radab)*(ralpprblp/rab-1.d0)
      ratio2_lp=exp(expo_lp)
      outl_lp=2.d0/(1.d0+ratio2_lp)

!     out-of-line damp: A...nb(B)-B
      do i=1,nbb
         ranbprbnb(i)=ranb(i)+rbnb(i)+1.d-12
         expo_nb(i)=(param%hbnbcut/radab)*(ranbprbnb(i)/rab-1.d0)
         ratio2_nb(i)=exp(-expo_nb(i))**(1.0)
         outl_nb(i)=( 2.d0/(1.d0+ratio2_nb(i)) ) - 1.0d0
      end do
      outl_nb_tot = product(outl_nb)

!     long damping
      ratio1=(rab2/param%hblongcut)**param%hbalp
      dampl=1.d0/(1.d0+ratio1)

!     short damping
      shortcut=param%hbscut*radab
      ratio3=(shortcut/rab2)**param%hbalp
      damps=1.d0/(1.d0+ratio3)

      damp  = damps*dampl
      ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
      rbhdamp = damp * ( (p_bh/rbh2/rbh) )
      rabdamp = damp * ( (p_ab/rab2/rab) )
      rdamp   = rbhdamp + rabdamp

!     hydrogen charge scaled term
      ex1h=exp(param%hbst*q(H))
      ex2h=ex1h+param%hbsf
      qh=ex1h/ex2h

!     hydrogen charge scaled term
      ex1a=exp(-param%hbst*q(A))
      ex2a=ex1a+param%hbsf
      qa=ex1a/ex2a

!     hydrogen charge scaled term
      ex1b=exp(-param%hbst*q(B))
      ex2b=ex1b+param%hbsf
      qb=ex1b/ex2b

      qhoutl=qh*outl*outl_nb_tot*outl_lp

!     constant values, no gradient
      const = ca(2)*qa*cb(1)*qb*param%xhaci_globabh

!     energy
      energy = -rdamp*qhoutl*const

!     gradient
      drah(1:3)=xyz(1:3,A)-xyz(1:3,H)
      drbh(1:3)=xyz(1:3,B)-xyz(1:3,H)
      drab(1:3)=xyz(1:3,A)-xyz(1:3,B)
      dralp(1:3)=xyz(1:3,A)-lp(1:3)
      drblp(1:3)=xyz(1:3,B)-lp(1:3)

      aterm  = -rdamp*qh*outl_nb_tot*outl_lp*const
      nbterm = -rdamp*qh*outl*outl_lp*const
      lpterm = -rdamp*qh*outl*outl_nb_tot*const
      dterm  = -qhoutl*const

!------------------------------------------------------------------------------
!     damping part: rab
      gi = ( (rabdamp+rbhdamp)*ddamp-3.d0*rabdamp ) / rab2
      gi = gi*dterm
      dg(1:3) = gi*drab(1:3)
      ga(1:3) =  dg(1:3)
      gb(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!     damping part: rbh
      gi = -3.d0*rbhdamp/rbh2
      gi = gi*dterm
      dg(1:3) = gi*drbh(1:3)
      gb(1:3) = gb(1:3)+dg(1:3)
      gh(1:3) =       - dg(1:3)

!------------------------------------------------------------------------------
!     angular A-H...B term
!------------------------------------------------------------------------------
!     out of line term: rab
      tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
      gi   = -tmp1 *rahprbh/rab2
      dg(1:3) = gi*drab(1:3)
      ga(1:3) = ga(1:3) + dg(1:3)
      gb(1:3) = gb(1:3) - dg(1:3)

!     out of line term: rah,rbh
      gi = tmp1/rah
      dga(1:3) = gi*drah(1:3)
      ga(1:3)  = ga(1:3) + dga(1:3)
      gi = tmp1/rbh
      dgb(1:3) = gi*drbh(1:3)
      gb(1:3)  = gb(1:3) + dgb(1:3)
      dgh(1:3) = -dga(1:3)-dgb(1:3)
      gh(1:3)  =  gh(1:3) + dgh(1:3)

!!------------------------------------------------------------------------------
!!     angular A...LP-B term
!!------------------------------------------------------------------------------
!     out of line term: rab
      tmp3 = -2.d0*lpterm*ratio2_lp*expo_lp/(1+ratio2_lp)**2/(ralpprblp-rab)
      gi   = -tmp3 *ralpprblp/rab2
      dg(1:3) = gi*drab(1:3)
      ga(1:3) = ga(1:3) + dg(1:3)
      gb(1:3) = gb(1:3) - dg(1:3)

!     out of line term: ralp,rblp
      gi = tmp3/ralp
      dga(1:3) = gi*dralp(1:3)
      ga(1:3)  = ga(1:3) + dga(1:3)
      gi = tmp3/(rblp+1.0d-12)
      dgb(1:3) = gi*drblp(1:3)
      gb(1:3)  = gb(1:3) - dga(1:3)
      glp(1:3) = -dga(1:3)!-dgb(1:3)

!     neighbor part: LP
      unit_vec=0
      do i=1,3
         unit_vec(i)=-1
         gii(1:3,i) = -lp_dist * dble(nbb) * ( unit_vec/vnorm + (vector*vector(i)/sum(vector**2)**(1.5d0)) )
         unit_vec=0
      end do
      gnb_lp=matmul(gii,glp)

!------------------------------------------------------------------------------
!     angular A...nb(B)-B term
!------------------------------------------------------------------------------
!     out of line term: rab
      mask_nb=.true.
      do i=1,nbb
         mask_nb(i)=.false.
         tmp2(i)  = 2.d0*nbterm*product(outl_nb,mask_nb)*ratio2_nb(i)*expo_nb(i)/&
                  & (1+ratio2_nb(i))**2/(ranbprbnb(i)-rab)
         gi_nb(i) = -tmp2(i) *ranbprbnb(i)/rab2
         dg(1:3)  = gi_nb(i)*drab(1:3)
         ga(1:3)  = ga(1:3) + dg(1:3)
         gb(1:3)  = gb(1:3) - dg(1:3)
         mask_nb=.true.
      end do

!     out of line term: ranb,rbnb
      do i=1,nbb
         gi_nb(i)   = tmp2(i)/ranb(i)
         dga(1:3)   = gi_nb(i)*dranb(1:3,i)
         ga(1:3)    = ga(1:3) + dga(1:3)
         gi_nb(i)   = tmp2(i)/rbnb(i)
         dgb(1:3)   = gi_nb(i)*drbnb(1:3,i)
         gb(1:3)    = gb(1:3) + dgb(1:3)
         dgnb(1:3)  = -dga(1:3)-dgb(1:3)
         gnb(1:3,i) = dgnb(1:3)
      end do

!------------------------------------------------------------------------------
      if(nbb.lt.1) then
         gdr(1:3,A) = gdr(1:3,A) + ga(1:3)
         gdr(1:3,B) = gdr(1:3,B) + gb(1:3)
         gdr(1:3,H) = gdr(1:3,H) + gh(1:3)
         return
      endif

!------------------------------------------------------------------------------
!     move gradients into place
      gdr(1:3,A) = gdr(1:3,A) + ga(1:3)
      gdr(1:3,B) = gdr(1:3,B) + gb(1:3) + gnb_lp(1:3)
      gdr(1:3,H) = gdr(1:3,H) + gh(1:3)
      do i=1,nbb
         gdr(1:3,topo%nb(i,B)) = gdr(1:3,topo%nb(i,B)) + gnb(1:3,i) - gnb_lp(1:3)/dble(nbb)
      end do

end subroutine abhgfnff_eg2_rnr

!subroutine for case 3: A-H...B, B is 0=C including two in plane LPs at B
!this is the multiplicative version of incorporationg etors and ebend
!equal to abhgfnff_eg2_new multiplied by etors and eangl
subroutine abhgfnff_eg3(n,A,B,H,at,xyz,q,sqrab,srab,energy,gdr,param,topo)
      use xtb_mctc_constants
      implicit none
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(in) :: topo
      integer A,B,H,n,at(n)
      real*8 xyz(3,n),energy,gdr(3,n)
      real*8 q(n)
      real*8 sqrab(n*(n+1)/2)   ! squared dist
      real*8 srab(n*(n+1)/2)    ! dist

      real*8 outl,dampl,damps,rdamp,damp
      real*8 ddamp,rabdamp,rbhdamp
      real*8 ratio1,ratio2,ratio2_nb(topo%nb(20,B)),ratio3
      real*8 xm,ym,zm
      real*8 rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
      real*8 ranb(topo%nb(20,B)),ranb2(topo%nb(20,B)),rbnb(topo%nb(20,B)),rbnb2(topo%nb(20,B))
      real*8 drah(3),drbh(3),drab(3),drm(3)
      real*8 dranb(3,topo%nb(20,B)),drbnb(3,topo%nb(20,B))
      real*8 dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
      real*8 ga(3),gb(3),gh(3),gnb(3,topo%nb(20,B))
      real*8 phi,phi0,r0,t0,fc,tshift,bshift
      real*8 eangl,etors,gangl(3,n),gtors(3,n)
      real*8 etmp(20),g3tmp(3,3),g4tmp(3,4,20)
      real*8 ratio,qhoutl,radab
      real*8 gi,gi_nb(topo%nb(20,B))
      real*8 tmp1,tmp2(topo%nb(20,B))
      real*8 rahprbh,ranbprbnb(topo%nb(20,B))
      real*8 ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_nb(topo%nb(20,B))
      real*8 eabh
      real*8 aterm,dterm,nbterm,bterm,tterm
      real*8 qa,qb,qh
      real*8 ca(2),cb(2)
      real*8 gqa,gqb,gqh
      real*8 shortcut
      real*8 tlist(5,topo%nb(20,topo%nb(1,B)))
      real*8 vtors(2,topo%nb(20,topo%nb(1,B)))
      real*8 valijklff
      real*8 const
      real*8 outl_nb(topo%nb(20,B)),outl_nb_tot
      logical mask_nb(topo%nb(20,B)),t_mask(20)

!     proportion between Rbh und Rab distance dependencies
      real*8 :: p_bh
      real*8 :: p_ab

      integer C,D
      integer i,j,ii,jj,kk,ll,ij,lina
      integer nbb,nbc
      integer ntors,rn

      lina(i,j)=min(i,j)+max(i,j)*(max(i,j)-1)/2

      p_bh=1.d0+param%hbabmix
      p_ab=    -param%hbabmix

      gdr    = 0
      energy = 0
      etors  = 0
      gtors  = 0
      eangl  = 0
      gangl  = 0

      call hbonds(A,B,ca,cb,param,topo)

      !Determine all neighbors for torsion term
      !  A
      !   \         tors:
      !    H        ll
      !     :        \
      !      O        jj
      !      ||       |
      !      C        kk
      !     / \       \
      !    R1  R2      ii
      !------------------------------------------
      nbb=topo%nb(20,B)
      C=topo%nb(nbb,B)
      nbc=topo%nb(20,C)
      ntors=nbc-nbb

      nbb=topo%nb(20,B)
!     Neighbours of B
      do i=1,nbb
!        compute distances
         dranb(1:3,i) = xyz(1:3,A) - xyz(1:3,topo%nb(i,B))
         drbnb(1:3,i) = xyz(1:3,B) - xyz(1:3,topo%nb(i,B))
!        A-nb(B) distance
         ranb2(i) = sum(dranb(1:3,i)**2)
         ranb(i)  = sqrt(ranb2(i))
!        B-nb(B) distance
         rbnb2(i) = sum(drbnb(1:3,i)**2)
         rbnb(i)  = sqrt(rbnb2(i))
      end do

!     A-B distance
      ij=lina(A,B)
      rab2=sqrab(ij)
      rab =srab (ij)

!     A-H distance
      ij=lina(A,H)
      rah2= sqrab(ij)
      rah = srab (ij)

!     B-H distance
      ij=lina(B,H)
      rbh2= sqrab(ij)
      rbh = srab (ij)

      rahprbh=rah+rbh+1.d-12
      radab=param%rad(at(A))+param%rad(at(B))

!     out-of-line damp: A-H...B
      expo=(param%hbacut/radab)*(rahprbh/rab-1.d0)
      if(expo.gt.15.0d0) return ! avoid overflow
      ratio2=exp(expo)
      outl=2.d0/(1.d0+ratio2)

!     out-of-line damp: A...nb(B)-B
      do i=1,nbb
         ranbprbnb(i)=ranb(i)+rbnb(i)+1.d-12
         expo_nb(i)=(param%hbnbcut/radab)*(ranbprbnb(i)/rab-1.d0)
         ratio2_nb(i)=exp(-expo_nb(i))
         outl_nb(i)=( 2.d0/(1.d0+ratio2_nb(i)) ) - 1.0d0
      end do
      outl_nb_tot = product(outl_nb)

!     long damping
      ratio1=(rab2/param%hblongcut)**param%hbalp
      dampl=1.d0/(1.d0+ratio1)

!     short damping
      shortcut=param%hbscut*radab
      ratio3=(shortcut/rab2)**6
      damps=1.d0/(1.d0+ratio3)

      damp  = damps*dampl
      ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
      rbhdamp = damp * ( (p_bh/rbh2/rbh) )
      rabdamp = damp * ( (p_ab/rab2/rab) )
      rdamp   = rbhdamp + rabdamp

      !Set up torsion paramter
      j = 0
      do i = 1,nbc
         if( topo%nb(i,C) == B ) cycle
         j = j + 1
         tlist(1,j)=topo%nb(i,C)
         tlist(2,j)=B
         tlist(3,j)=C
         tlist(4,j)=H
         tlist(5,j)=2
         t0=180
         phi0=t0*pi/180.
         vtors(1,j)=phi0-(pi/2.0)
         vtors(2,j)=param%tors_hb
      end do

      !Calculate etors
      do i = 1,ntors
         ii =tlist(1,i)
         jj =tlist(2,i)
         kk =tlist(3,i)
         ll =tlist(4,i)
         rn =tlist(5,i)
         phi0=vtors(1,i)
         tshift=vtors(2,i)
         phi=valijklff(n,xyz,ii,jj,kk,ll)
         call egtors_nci_mul(ii,jj,kk,ll,rn,phi0,tshift,n,at,xyz,etmp(i),g4tmp(:,:,i))
      end do
      etors=product(etmp(1:ntors))
      !Calculate gtors
      t_mask=.true.
      do i = 1,ntors
         t_mask(i)=.false.
         ii =tlist(1,i)
         jj =tlist(2,i)
         kk =tlist(3,i)
         ll =tlist(4,i)
         gtors(1:3,ii) = gtors(1:3,ii)+g4tmp(1:3,1,i)*product(etmp(1:ntors),t_mask(1:ntors))
         gtors(1:3,jj) = gtors(1:3,jj)+g4tmp(1:3,2,i)*product(etmp(1:ntors),t_mask(1:ntors))
         gtors(1:3,kk) = gtors(1:3,kk)+g4tmp(1:3,3,i)*product(etmp(1:ntors),t_mask(1:ntors))
         gtors(1:3,ll) = gtors(1:3,ll)+g4tmp(1:3,4,i)*product(etmp(1:ntors),t_mask(1:ntors))
         t_mask=.true.
      end do

      !Calculate eangl + gangl
      r0=120
      phi0=r0*pi/180.
      bshift=param%bend_hb
      fc=1.0d0-bshift
      call bangl(xyz,kk,jj,ll,phi)
      call egbend_nci_mul(jj,kk,ll,phi0,fc,n,at,xyz,eangl,g3tmp)
      gangl(1:3,jj)=gangl(1:3,jj)+g3tmp(1:3,1)
      gangl(1:3,kk)=gangl(1:3,kk)+g3tmp(1:3,2)
      gangl(1:3,ll)=gangl(1:3,ll)+g3tmp(1:3,3)

!     hydrogen charge scaled term
      ex1h=exp(param%hbst*q(H))
      ex2h=ex1h+param%hbsf
      qh=ex1h/ex2h

!     hydrogen charge scaled term
      ex1a=exp(-param%hbst*q(A))
      ex2a=ex1a+param%hbsf
      qa=ex1a/ex2a

!     hydrogen charge scaled term
      ex1b=exp(-param%hbst*q(B))
      ex2b=ex1b+param%hbsf
      qb=ex1b/ex2b

      qhoutl=qh*outl*outl_nb_tot

!     constant values, no gradient
      const = ca(2)*qa*cb(1)*qb*param%xhaci_coh

!     energy
      energy = -rdamp*qhoutl*eangl*etors*const

!     gradient
      drah(1:3)=xyz(1:3,A)-xyz(1:3,H)
      drbh(1:3)=xyz(1:3,B)-xyz(1:3,H)
      drab(1:3)=xyz(1:3,A)-xyz(1:3,B)

      aterm  = -rdamp*qh*outl_nb_tot*eangl*etors*const
      nbterm = -rdamp*qh*outl*eangl*etors*const
      dterm  = -qhoutl*eangl*etors*const
      tterm  = -rdamp*qhoutl*eangl*const
      bterm  = -rdamp*qhoutl*etors*const

!------------------------------------------------------------------------------
!     damping part: rab
      gi = ( (rabdamp+rbhdamp)*ddamp-3.d0*rabdamp ) / rab2
      gi = gi*dterm
      dg(1:3) = gi*drab(1:3)
      ga(1:3) =  dg(1:3)
      gb(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!     damping part: rbh
      gi = -3.d0*rbhdamp/rbh2
      gi = gi*dterm
      dg(1:3) = gi*drbh(1:3)
      gb(1:3) = gb(1:3)+dg(1:3)
      gh(1:3) =       - dg(1:3)

!------------------------------------------------------------------------------
!     angular A-H...B term
!------------------------------------------------------------------------------
!     out of line term: rab
      tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
      gi   = -tmp1 *rahprbh/rab2
      dg(1:3) = gi*drab(1:3)
      ga(1:3) = ga(1:3) + dg(1:3)
      gb(1:3) = gb(1:3) - dg(1:3)

!     out of line term: rah,rbh
      gi = tmp1/rah
      dga(1:3) = gi*drah(1:3)
      ga(1:3)  = ga(1:3) + dga(1:3)
      gi = tmp1/rbh
      dgb(1:3) = gi*drbh(1:3)
      gb(1:3)  = gb(1:3) + dgb(1:3)
      dgh(1:3) = -dga(1:3)-dgb(1:3)
      gh(1:3)  =  gh(1:3) + dgh(1:3)

!------------------------------------------------------------------------------
!     angular A...nb(B)-B term
!------------------------------------------------------------------------------
!     out of line term: rab
      mask_nb=.true.
      do i=1,nbb
         mask_nb(i)=.false.
         tmp2(i)  = 2.d0*nbterm*product(outl_nb,mask_nb)*ratio2_nb(i)*expo_nb(i)/&
                  & (1+ratio2_nb(i))**2/(ranbprbnb(i)-rab)
         gi_nb(i) = -tmp2(i) *ranbprbnb(i)/rab2
         dg(1:3)  = gi_nb(i)*drab(1:3)
         ga(1:3)  = ga(1:3) + dg(1:3)
         gb(1:3)  = gb(1:3) - dg(1:3)
         mask_nb=.true.
      end do

!     out of line term: ranb,rbnb
      do i=1,nbb
         gi_nb(i)   = tmp2(i)/ranb(i)
         dga(1:3)   = gi_nb(i)*dranb(1:3,i)
         ga(1:3)    = ga(1:3) + dga(1:3)
         gi_nb(i)   = tmp2(i)/rbnb(i)
         dgb(1:3)   = gi_nb(i)*drbnb(1:3,i)
         gb(1:3)    = gb(1:3) + dgb(1:3)
         dgnb(1:3)  = -dga(1:3)-dgb(1:3)
         gnb(1:3,i) = dgnb(1:3)
      end do

!------------------------------------------------------------------------------
      if(nbb.lt.1) then
         gdr(1:3,A) = gdr(1:3,A) + ga(1:3)
         gdr(1:3,B) = gdr(1:3,B) + gb(1:3)
         gdr(1:3,H) = gdr(1:3,H) + gh(1:3)
         return
      endif

!------------------------------------------------------------------------------
!     torsion term H...B=C<R1,R2
!------------------------------------------------------------------------------
      do i = 1,ntors
         ii =tlist(1,i)
         gdr(1:3,ii) = gdr(1:3,ii)+gtors(1:3,ii)*tterm
      end do
      gdr(1:3,jj) = gdr(1:3,jj)+gtors(1:3,jj)*tterm
      gdr(1:3,kk) = gdr(1:3,kk)+gtors(1:3,kk)*tterm
      gdr(1:3,ll) = gdr(1:3,ll)+gtors(1:3,ll)*tterm

!------------------------------------------------------------------------------
!     angle term H...B=C
!------------------------------------------------------------------------------
      gdr(1:3,jj) = gdr(1:3,jj)+gangl(1:3,jj)*bterm
      gdr(1:3,kk) = gdr(1:3,kk)+gangl(1:3,kk)*bterm
      gdr(1:3,ll) = gdr(1:3,ll)+gangl(1:3,ll)*bterm

!------------------------------------------------------------------------------
!     move gradients into place
      gdr(1:3,A) = gdr(1:3,A) + ga(1:3)
      gdr(1:3,B) = gdr(1:3,B) + gb(1:3)
      gdr(1:3,H) = gdr(1:3,H) + gh(1:3)
      do i=1,nbb
         gdr(1:3,topo%nb(i,B)) = gdr(1:3,topo%nb(i,B)) + gnb(1:3,i)
      end do

end subroutine abhgfnff_eg3

!subroutine for case 3: A-H...B, B is 0=C including two in plane LPs at B
!this is the multiplicative version of incorporationg etors and ebend without neighbor LP
subroutine abhgfnff_eg3_mul(n,A,B,H,at,xyz,q,sqrab,srab,energy,gdr,param,topo)
      use xtb_mctc_constants
      implicit none
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(in) :: topo
      integer A,B,H,C,D,n,at(n)
      real*8 xyz(3,n),energy,gdr(3,n)
      real*8 q(n)
      real*8 sqrab(n*(n+1)/2)   ! squared dist
      real*8 srab(n*(n+1)/2)    ! dist

      real*8 outl,dampl,damps,rdamp,damp
      real*8 ddamp,rabdamp,rbhdamp
      real*8 ratio1,ratio2,ratio2_nb(topo%nb(20,B)),ratio3
      real*8 xm,ym,zm
      real*8 rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
      real*8 ranb(topo%nb(20,B)),ranb2(topo%nb(20,B)),rbnb(topo%nb(20,B)),rbnb2(topo%nb(20,B))
      real*8 drah(3),drbh(3),drab(3),drm(3)
      real*8 dranb(3,topo%nb(20,B)),drbnb(3,topo%nb(20,B))
      real*8 dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
      real*8 ga(3),gb(3),gh(3),gnb(3,topo%nb(20,B))
      real*8 phi,phi0,r0,fc,tshift,bshift
      real*8 eangl,etors,gangl(3,n),gtors(3,n)
      real*8 etmp,g3tmp(3,3),g4tmp(3,4)
      real*8 denom,ratio,qhoutl,radab
      real*8 gi,gi_nb(topo%nb(20,B))
      real*8 tmp1,tmp2(topo%nb(20,B))
      real*8 rahprbh,ranbprbnb(topo%nb(20,B))
      real*8 ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_nb(topo%nb(20,B))
      real*8 eabh
      real*8 aterm,dterm,bterm,tterm
      real*8 qa,qb,qh
      real*8 ca(2),cb(2)
      real*8 gqa,gqb,gqh
      real*8 shortcut
      real*8 const
      real*8 tlist(5,topo%nb(20,topo%nb(1,B)))
      real*8 vtors(2,topo%nb(20,topo%nb(1,B)))
      real*8 valijklff
      logical mask_nb(topo%nb(20,B))

!     proportion between Rbh und Rab distance dependencies
      real*8 :: p_bh
      real*8 :: p_ab

      integer i,j,ii,jj,kk,ll,ij,lina
      integer nbb,nbc
      integer ntors,rn

      lina(i,j)=min(i,j)+max(i,j)*(max(i,j)-1)/2

      p_bh=1.d0+param%hbabmix
      p_ab=    -param%hbabmix

      gdr    = 0
      energy = 0
      etors  = 0
      gtors  = 0
      eangl  = 0
      gangl  = 0

      call hbonds(A,B,ca,cb,param,topo)

      !Determine all neighbors for torsion term
      !  A
      !   \         tors:
      !    H        ll
      !     :        \
      !      O        jj
      !      ||       |
      !      C        kk
      !     / \       \
      !    R1  R2      ii
      !------------------------------------------
      nbb=topo%nb(20,B)
      C=topo%nb(nbb,B)
      nbc=topo%nb(20,C)
      ntors=nbc-nbb

      !A-B distance
      ij=lina(A,B)
      rab2=sqrab(ij)
      rab =srab (ij)

      !A-H distance
      ij=lina(A,H)
      rah2= sqrab(ij)
      rah = srab (ij)

      !B-H distance
      ij=lina(B,H)
      rbh2= sqrab(ij)
      rbh = srab (ij)

      rahprbh=rah+rbh+1.d-12
      radab=param%rad(at(A))+param%rad(at(B))

!     out-of-line damp: A-H...B
      expo=(param%hbacut/radab)*(rahprbh/rab-1.d0)
      if(expo.gt.15.0d0) return ! avoid overflow
      ratio2=exp(expo)
      outl=2.d0/(1.d0+ratio2)

!     long damping
      ratio1=(rab2/param%hblongcut)**param%hbalp
      dampl=1.d0/(1.d0+ratio1)

!     short damping
      shortcut=param%hbscut*radab
      ratio3=(shortcut/rab2)**param%hbalp
      damps=1.d0/(1.d0+ratio3)

      damp  = damps*dampl
      ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
      rbhdamp = damp * ( (p_bh/rbh2/rbh) )
      rabdamp = damp * ( (p_ab/rab2/rab) )
      rdamp   = rbhdamp + rabdamp

      !Set up torsion paramter
      j = 0
      do i = 1,nbc
         if( topo%nb(i,C) == B ) cycle
         j = j + 1
         tlist(1,j)=topo%nb(i,C)
         tlist(2,j)=B
         tlist(3,j)=C
         tlist(4,j)=H
         tlist(5,j)=2
         vtors(1,j)=pi/2.0
         vtors(2,j)=0.70
      end do

      !Calculate etors
      do i = 1,ntors
         ii =tlist(1,i)
         jj =tlist(2,i)
         kk =tlist(3,i)
         ll =tlist(4,i)
         rn =tlist(5,i)
         phi0=vtors(1,i)
         tshift=vtors(2,i)
         phi=valijklff(n,xyz,ii,jj,kk,ll)
         call egtors_nci_mul(ii,jj,kk,ll,rn,phi0,tshift,n,at,xyz,etmp,g4tmp)
         gtors(1:3,ii) = gtors(1:3,ii)+g4tmp(1:3,1)
         gtors(1:3,jj) = gtors(1:3,jj)+g4tmp(1:3,2)
         gtors(1:3,kk) = gtors(1:3,kk)+g4tmp(1:3,3)
         gtors(1:3,ll) = gtors(1:3,ll)+g4tmp(1:3,4)
         etors = etors + etmp
      end do
      etors=etors/ntors

      r0=120
      phi0=r0*pi/180.
      bshift=0.1
      fc=1.0d0-bshift
      call bangl(xyz,kk,jj,ll,phi)
      call egbend_nci_mul(jj,kk,ll,phi0,fc,n,at,xyz,etmp,g3tmp)
      gangl(1:3,jj)=gangl(1:3,jj)+g3tmp(1:3,1)
      gangl(1:3,kk)=gangl(1:3,kk)+g3tmp(1:3,2)
      gangl(1:3,ll)=gangl(1:3,ll)+g3tmp(1:3,3)
      eangl=eangl+etmp

!     hydrogen charge scaled term
      ex1h=exp(param%hbst*q(H))
      ex2h=ex1h+param%hbsf
      qh=ex1h/ex2h

!     hydrogen charge scaled term
      ex1a=exp(-param%hbst*q(A))
      ex2a=ex1a+param%hbsf
      qa=ex1a/ex2a

!     hydrogen charge scaled term
      ex1b=exp(-param%hbst*q(B))
      ex2b=ex1b+param%hbsf
      qb=ex1b/ex2b

!     max distance to neighbors excluded, would lead to linear C=O-H
      qhoutl=qh*outl

!     constant values, no gradient
      const = ca(2)*qa*cb(1)*qb*param%xhaci_globabh

!     energy
      energy = -rdamp*qhoutl*const*eangl*etors

!     gradient
      drah(1:3)=xyz(1:3,A)-xyz(1:3,H)
      drbh(1:3)=xyz(1:3,B)-xyz(1:3,H)
      drab(1:3)=xyz(1:3,A)-xyz(1:3,B)

      aterm  = -rdamp*qh*etors*eangl*const
      dterm  = -qhoutl*etors*eangl*const
      tterm  = -rdamp*qhoutl*eangl*const/ntors
      bterm  = -rdamp*qhoutl*etors*const

!------------------------------------------------------------------------------
!     damping part: rab
      gi = ( (rabdamp+rbhdamp)*ddamp-3.d0*rabdamp ) / rab2
      gi = gi*dterm
      dg(1:3) = gi*drab(1:3)
      ga(1:3) =  dg(1:3)
      gb(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!     damping part: rbh
      gi = -3.d0*rbhdamp/rbh2
      gi = gi*dterm
      dg(1:3) = gi*drbh(1:3)
      gb(1:3) = gb(1:3)+dg(1:3)
      gh(1:3) =       - dg(1:3)

!------------------------------------------------------------------------------
!     angular A-H...B term
!------------------------------------------------------------------------------
!     out of line term: rab
      tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
      gi   = -tmp1 *rahprbh/rab2
      dg(1:3) = gi*drab(1:3)
      ga(1:3) = ga(1:3) + dg(1:3)
      gb(1:3) = gb(1:3) - dg(1:3)

!     out of line term: rah,rbh
      gi = tmp1/rah
      dga(1:3) = gi*drah(1:3)
      ga(1:3)  = ga(1:3) + dga(1:3)
      gi = tmp1/rbh
      dgb(1:3) = gi*drbh(1:3)
      gb(1:3)  = gb(1:3) + dgb(1:3)
      dgh(1:3) = -dga(1:3)-dgb(1:3)
      gh(1:3)  =  gh(1:3) + dgh(1:3)


!------------------------------------------------------------------------------
!     torsion term H...B=C<R1,R2
!------------------------------------------------------------------------------
      do i = 1,ntors
         ii =tlist(1,i)
         gdr(1:3,ii) = gdr(1:3,ii)+gtors(1:3,ii)*tterm
      end do
      gdr(1:3,jj) = gdr(1:3,jj)+gtors(1:3,jj)*tterm
      gdr(1:3,kk) = gdr(1:3,kk)+gtors(1:3,kk)*tterm
      gdr(1:3,ll) = gdr(1:3,ll)+gtors(1:3,ll)*tterm

!------------------------------------------------------------------------------
!     angle term H...B=C
!------------------------------------------------------------------------------
      gdr(1:3,jj) = gdr(1:3,jj)+gangl(1:3,jj)*bterm
      gdr(1:3,kk) = gdr(1:3,kk)+gangl(1:3,kk)*bterm
      gdr(1:3,ll) = gdr(1:3,ll)+gangl(1:3,ll)*bterm

!------------------------------------------------------------------------------
!     move gradients into place
      gdr(1:3,A) = gdr(1:3,A) + ga(1:3)
      gdr(1:3,B) = gdr(1:3,B) + gb(1:3)
      gdr(1:3,H) = gdr(1:3,H) + gh(1:3)

end subroutine abhgfnff_eg3_mul

!subroutine for case 3: A-H...B, B is 0=C including two in plane LPs at B
!this is the additive version of incorporationg etors and ebend
subroutine abhgfnff_eg3_add(n,A,B,H,at,xyz,q,sqrab,srab,energy,gdr,param,topo)
      use xtb_mctc_constants
      implicit none
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(in) :: topo
      integer A,B,H,C,D,n,at(n)
      real*8 xyz(3,n),energy,gdr(3,n)
      real*8 q(n)
      real*8 sqrab(n*(n+1)/2)   ! squared dist
      real*8 srab(n*(n+1)/2)    ! dist

      real*8 outl,dampl,damps,rdamp,damp
      real*8 ddamp,rabdamp,rbhdamp
      real*8 ratio1,ratio2,ratio2_nb(topo%nb(20,B)),ratio3
      real*8 xm,ym,zm
      real*8 rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
      real*8 ranb(topo%nb(20,B)),ranb2(topo%nb(20,B)),rbnb(topo%nb(20,B)),rbnb2(topo%nb(20,B))
      real*8 drah(3),drbh(3),drab(3),drm(3)
      real*8 dranb(3,topo%nb(20,B)),drbnb(3,topo%nb(20,B))
      real*8 dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
      real*8 ga(3),gb(3),gh(3),gnb(3,topo%nb(20,B))
      real*8 phi,phi0,r0,fc
      real*8 eangl,etors
      real*8 etmp,g3tmp(3,3),g4tmp(3,4)
      real*8 denom,ratio,qhoutl,radab
      real*8 gi,gi_nb(topo%nb(20,B))
      real*8 tmp1,tmp2(topo%nb(20,B))
      real*8 rahprbh,ranbprbnb(topo%nb(20,B))
      real*8 ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_nb(topo%nb(20,B))
      real*8 eabh
      real*8 aterm,dterm,nbterm
      real*8 qa,qb,qh
      real*8 ca(2),cb(2)
      real*8 gqa,gqb,gqh
      real*8 shortcut
      real*8 const
      real*8 outl_nb(topo%nb(20,B)),outl_nb_tot
      real*8 tlist(5,topo%nb(20,topo%nb(1,B)))
      real*8 vtors(2,topo%nb(20,topo%nb(1,B)))
      real*8 valijklff
      logical mask_nb(topo%nb(20,B))

!     proportion between Rbh und Rab distance dependencies
      real*8 :: p_bh
      real*8 :: p_ab

      integer i,j,ii,jj,kk,ll,ij,lina
      integer nbb,nbc
      integer ntors,rn

      lina(i,j)=min(i,j)+max(i,j)*(max(i,j)-1)/2

      p_bh=1.d0+param%hbabmix
      p_ab=    -param%hbabmix

      gdr    = 0
      energy = 0
      etors  = 0
      eangl  = 0

      call hbonds(A,B,ca,cb,param,topo)

      !Determine all neighbors for torsion term
      !  A
      !   \         tors:
      !    H        ll
      !     :        \
      !      O        jj
      !      ||       |
      !      C        kk
      !     / \       \
      !    R1  R2      ii
      !------------------------------------------
      nbb=topo%nb(20,B)
      C=topo%nb(nbb,B)
      nbc=topo%nb(20,C)
      ntors=nbc-nbb

      !A-B distance
      ij=lina(A,B)
      rab2=sqrab(ij)
      rab =srab (ij)

      !A-H distance
      ij=lina(A,H)
      rah2= sqrab(ij)
      rah = srab (ij)

      !B-H distance
      ij=lina(B,H)
      rbh2= sqrab(ij)
      rbh = srab (ij)

      rahprbh=rah+rbh+1.d-12
      radab=param%rad(at(A))+param%rad(at(B))

!     out-of-line damp: A-H...B
      expo=(param%hbacut/radab)*(rahprbh/rab-1.d0)
      if(expo.gt.15.0d0) return ! avoid overflow
      ratio2=exp(expo)
      outl=2.d0/(1.d0+ratio2)

!     long damping
      ratio1=(rab2/param%hblongcut)**param%hbalp
      dampl=1.d0/(1.d0+ratio1)

!     short damping
      shortcut=param%hbscut*radab
      ratio3=(shortcut/rab2)**param%hbalp
      damps=1.d0/(1.d0+ratio3)

      damp  = damps*dampl
      ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
      rbhdamp = damp * ( (p_bh/rbh2/rbh) )
      rabdamp = damp * ( (p_ab/rab2/rab) )
      rdamp   = rbhdamp + rabdamp

      !Set up torsion paramter
      j = 0
      do i = 1,nbc
         if( topo%nb(i,C) == B ) cycle
         j = j + 1
         tlist(1,j)=topo%nb(i,C)
         tlist(2,j)=B
         tlist(3,j)=C
         tlist(4,j)=H
         tlist(5,j)=2
         vtors(1,j)=pi
         vtors(2,j)=0.30
      end do

      !Calculate etors
      !write(*,*)'torsion atoms       nrot      phi0   phi      FC'
      do i = 1,ntors
         ii =tlist(1,i)
         jj =tlist(2,i)
         kk =tlist(3,i)
         ll =tlist(4,i)
         rn=tlist(5,i)
         phi0=vtors(1,i)
         fc=vtors(2,i)
         phi=valijklff(n,xyz,ii,jj,kk,ll)
         call egtors_nci(ii,jj,kk,ll,rn,phi0,fc,n,at,xyz,etmp,g4tmp,param)
         gdr(1:3,ii) = gdr(1:3,ii)+g4tmp(1:3,1)
         gdr(1:3,jj) = gdr(1:3,jj)+g4tmp(1:3,2)
         gdr(1:3,kk) = gdr(1:3,kk)+g4tmp(1:3,3)
         gdr(1:3,ll) = gdr(1:3,ll)+g4tmp(1:3,4)
         etors = etors + etmp
      end do

      !Calculate eangl + gangl
      write(*,*)'angle atoms          phi0   phi      FC'
      r0=120
      phi0=r0*pi/180.
      fc=0.20
      call bangl(xyz,kk,jj,ll,phi)
      write(*,'(3i5,2x,3f8.3)') &
      &   jj,kk,ll,phi0*180./pi,phi*180./pi,fc
      call egbend_nci(jj,kk,ll,phi0,fc,n,at,xyz,etmp,g3tmp,param)
      gdr(1:3,jj)=gdr(1:3,jj)+g3tmp(1:3,1)
      gdr(1:3,kk)=gdr(1:3,kk)+g3tmp(1:3,2)
      gdr(1:3,ll)=gdr(1:3,ll)+g3tmp(1:3,3)
      eangl=eangl+etmp

!     hydrogen charge scaled term
      ex1h=exp(param%hbst*q(H))
      ex2h=ex1h+param%hbsf
      qh=ex1h/ex2h

!     hydrogen charge scaled term
      ex1a=exp(-param%hbst*q(A))
      ex2a=ex1a+param%hbsf
      qa=ex1a/ex2a

!     hydrogen charge scaled term
      ex1b=exp(-param%hbst*q(B))
      ex2b=ex1b+param%hbsf
      qb=ex1b/ex2b

!     max distance to neighbors excluded, would lead to linear C=O-H
      qhoutl=qh*outl

!     constant values, no gradient
      const = ca(2)*qa*cb(1)*qb*param%xhaci_globabh

!     energy
      energy = -rdamp*qhoutl*const+etors+eangl

!     gradient
      drah(1:3)=xyz(1:3,A)-xyz(1:3,H)
      drbh(1:3)=xyz(1:3,B)-xyz(1:3,H)
      drab(1:3)=xyz(1:3,A)-xyz(1:3,B)

      aterm  = -rdamp*qh*const
      dterm  = -qhoutl*const

!------------------------------------------------------------------------------
!     damping part: rab
      gi = ( (rabdamp+rbhdamp)*ddamp-3.d0*rabdamp ) / rab2
      gi = gi*dterm
      dg(1:3) = gi*drab(1:3)
      ga(1:3) =  dg(1:3)
      gb(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!     damping part: rbh
      gi = -3.d0*rbhdamp/rbh2
      gi = gi*dterm
      dg(1:3) = gi*drbh(1:3)
      gb(1:3) = gb(1:3)+dg(1:3)
      gh(1:3) =       - dg(1:3)

!------------------------------------------------------------------------------
!     angular A-H...B term
!------------------------------------------------------------------------------
!     out of line term: rab
      tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
      gi   = -tmp1 *rahprbh/rab2
      dg(1:3) = gi*drab(1:3)
      ga(1:3) = ga(1:3) + dg(1:3)
      gb(1:3) = gb(1:3) - dg(1:3)

!     out of line term: rah,rbh
      gi = tmp1/rah
      dga(1:3) = gi*drah(1:3)
      ga(1:3)  = ga(1:3) + dga(1:3)
      gi = tmp1/rbh
      dgb(1:3) = gi*drbh(1:3)
      gb(1:3)  = gb(1:3) + dgb(1:3)
      dgh(1:3) = -dga(1:3)-dgb(1:3)
      gh(1:3)  =  gh(1:3) + dgh(1:3)

!------------------------------------------------------------------------------
!     move gradients into place
      gdr(1:3,A) = gdr(1:3,A) + ga(1:3)
      gdr(1:3,B) = gdr(1:3,B) + gb(1:3)
      gdr(1:3,H) = gdr(1:3,H) + gh(1:3)

end subroutine abhgfnff_eg3_add

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! XB energy and analytical gradient
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine rbxgfnff_eg(n,A,B,X,at,xyz,q,energy,gdr,param)
      implicit none
      type(TGFFData), intent(in) :: param
      integer               :: A,B,X,n,at(n)
      real*8                :: xyz(3,n)
      real*8,intent(inout)  :: energy,gdr(3,3)
      real*8                :: q(n)

      real*8 outl,dampl,damps,rdamp,damp
      real*8 ratio1,ratio2,ratio3
      real*8 rab,rax,rbx,rab2,rax2,rbx2,rax4,rbx4
      real*8 drax(3),drbx(3),drab(3),drm(3)
      real*8 dg(3),dga(3),dgb(3),dgx(3)
      real*8 gi,ga(3),gb(3),gx(3)
      real*8 ex1_a,ex2_a,ex1_b,ex2_b,ex1_x,ex2_x,expo
      real*8 aterm,dterm
      real*8 qa,qb,qx
      real*8 cx,cb
      real*8 gqa,gqb,gqx
      real*8 shortcut,const

      integer i,j

      gdr =  0
      energy=0

      cb = 1.!param%xhbas(at(B))
      cx = param%xbaci(at(X))

!     compute distances
      drax(1:3) = xyz(1:3,A)-xyz(1:3,X)
      drbx(1:3) = xyz(1:3,B)-xyz(1:3,X)
      drab(1:3) = xyz(1:3,A)-xyz(1:3,B)

!     A-B distance
      rab2 = sum(drab**2)
      rab  = sqrt(rab2)

!     A-X distance
      rax2 = sum(drax**2)
      rax  = sqrt(rax2)+1.d-12

!     B-X distance
      rbx2 = sum(drbx**2)
      rbx  = sqrt(rbx2)+1.d-12

!     out-of-line damp
      expo   = param%xbacut*((rax+rbx)/rab-1.d0)
      if(expo.gt.15.0d0) return ! avoid overflow
      ratio2 = exp(expo)
      outl   = 2.d0/(1.d0+ratio2)

!     long damping
      ratio1 = (rbx2/param%hblongcut_xb)**param%hbalp
      dampl  = 1.d0/(1.d0+ratio1)

!     short damping
      shortcut = param%xbscut*(param%rad(at(A))+param%rad(at(B)))
      ratio3   = (shortcut/rbx2)**param%hbalp
      damps    = 1.d0/(1.d0+ratio3)

      damp  = damps*dampl
      rdamp = damp/rbx2/rbx ! **2

!     halogen charge scaled term
      ex1_x = exp(param%xbst*q(X))
      ex2_x = ex1_x+param%xbsf
      qx    = ex1_x/ex2_x

!     donor charge scaled term
      ex1_b = exp(-param%xbst*q(B))
      ex2_b = ex1_b+param%xbsf
      qb    = ex1_b/ex2_b

!     constant values, no gradient
      const = cb*qb*cx*qx

!     r^3 only sligxtly better than r^4
      aterm = -rdamp*const
      dterm = -outl*const
      energy= -rdamp*outl*const

!     damping part rab
      gi = rdamp*(-(2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3&
     &     /(1.d0+ratio3))-3.d0)/rbx2   ! 4,5,6 instead of 3.
      gi = gi*dterm
      dg(1:3) = gi*drbx(1:3)
      gb(1:3) =  dg(1:3)
      gx(1:3) = -dg(1:3)

!     out of line term: rab
      gi =2.d0*ratio2*expo*(rax+rbx)/(1.d0+ratio2)**2/(rax+rbx-rab)/rab2
      gi = gi*aterm
      dg(1:3) = gi*drab(1:3)
      ga(1:3) =         + dg(1:3)
      gb(1:3) = gb(1:3) - dg(1:3)

!     out of line term: rax,rbx
      gi = -2.d0*ratio2*expo/(1.d0+ratio2)**2/(rax+rbx-rab)/rax
      gi = gi*aterm
      dga(1:3) = gi*drax(1:3)
      ga(1:3) = ga(1:3) + dga(1:3)
      gi = -2.d0*ratio2*expo/(1.d0+ratio2)**2/(rax+rbx-rab)/rbx
      gi = gi*aterm
      dgb(1:3) = gi*drbx(1:3)
      gb(1:3) = gb(1:3) + dgb(1:3)
      dgx(1:3) = -dga(1:3) - dgb(1:3)
      gx(1:3) = gx(1:3) + dgx(1:3)

!     move gradients into place

      gdr(1:3,1) = ga(1:3)
      gdr(1:3,2) = gb(1:3)
      gdr(1:3,3) = gx(1:3)

      return
      end subroutine rbxgfnff_eg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! taken from D3 ATM code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine batmgfnff_eg(n,iat,jat,kat,at,xyz,q,sqrab,srab,energy,g,param)
      implicit none
      type(TGFFData), intent(in) :: param
      integer, intent(in) :: iat,jat,kat,n,at(n)
      real*8, intent(in) :: xyz(3,n),q(n)
      real*8, intent(out) :: energy,g(3,3)
      real*8, intent(in) :: sqrab(n*(n+1)/2)   ! squared dist
      real*8, intent(in) :: srab (n*(n+1)/2)   ! dist

      real*8 r2ij,r2jk,r2ik,c9,mijk,imjk,ijmk,rijk3,ang,angr9,rav3
      real*8 rij(3),rik(3),rjk(3),drij,drik,drjk,dang,ff,fi,fj,fk,fqq
      parameter (fqq=3.0d0)
      integer linij,linik,linjk,lina,i,j
      lina(i,j)=min(i,j)+max(i,j)*(max(i,j)-1)/2

      fi=(1.d0-fqq*q(iat))
      fi=min(max(fi,-4.0d0),4.0d0)
      fj=(1.d0-fqq*q(jat))
      fj=min(max(fj,-4.0d0),4.0d0)
      fk=(1.d0-fqq*q(kat))
      fk=min(max(fk,-4.0d0),4.0d0)
      ff=fi*fj*fk ! charge term
      c9=ff*param%zb3atm(at(iat))*param%zb3atm(at(jat))*param%zb3atm(at(kat)) ! strength of interaction
      linij=lina(iat,jat)
      linik=lina(iat,kat)
      linjk=lina(jat,kat)
      r2ij=sqrab(linij)
      r2jk=sqrab(linjk)
      r2ik=sqrab(linik)
      mijk=-r2ij+r2jk+r2ik
      imjk= r2ij-r2jk+r2ik
      ijmk= r2ij+r2jk-r2ik
      rijk3=r2ij*r2jk*r2ik
      rav3=rijk3**1.5 ! R^9
      ang=0.375d0*ijmk*imjk*mijk/rijk3
      angr9=(ang +1.0d0)/rav3
      energy=c9*angr9 ! energy

!     derivatives of each part w.r.t. r_ij,jk,ik
              dang=-0.375d0*(r2ij**3+r2ij**2*(r2jk+r2ik) &
     &             +r2ij*(3.0d0*r2jk**2+2.0*r2jk*r2ik+3.0*r2ik**2) &
     &             -5.0*(r2jk-r2ik)**2*(r2jk+r2ik)) &
     &             /(srab(linij)*rijk3*rav3)
              drij=-dang*c9
              dang=-0.375d0*(r2jk**3+r2jk**2*(r2ik+r2ij) &
     &             +r2jk*(3.0d0*r2ik**2+2.0*r2ik*r2ij+3.0*r2ij**2) &
     &             -5.0*(r2ik-r2ij)**2*(r2ik+r2ij)) &
     &             /(srab(linjk)*rijk3*rav3)
              drjk=-dang*c9
              dang=-0.375d0*(r2ik**3+r2ik**2*(r2jk+r2ij) &
     &             +r2ik*(3.0d0*r2jk**2+2.0*r2jk*r2ij+3.0*r2ij**2) &
     &             -5.0*(r2jk-r2ij)**2*(r2jk+r2ij)) &
     &             /(srab(linik)*rijk3*rav3)
              drik=-dang*c9

      rij=xyz(:,jat)-xyz(:,iat)
      rik=xyz(:,kat)-xyz(:,iat)
      rjk=xyz(:,kat)-xyz(:,jat)
      g(:,1  )=         drij*rij/srab(linij)
      g(:,1  )=g(:,1  )+drik*rik/srab(linik)
      g(:,2  )=         drjk*rjk/srab(linjk)
      g(:,2  )=g(:,2  )-drij*rij/srab(linij)
      g(:,3  )=        -drik*rik/srab(linik)
      g(:,3  )=g(:,3  )-drjk*rjk/srab(linjk)

      end subroutine batmgfnff_eg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CN routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> logCN derivative saved in dlogCN array
subroutine gfnff_dlogcoord(n,at,xyz,rab,logCN,dlogCN,thr2,param)
use xtb_mctc_accuracy, only : wp
      implicit none
      type(TGFFData), intent(in) :: param
      integer, intent(in)  :: n
      integer, intent(in)  :: at(n)
      real(wp),intent(in)  :: xyz(3,n)
      real(wp),intent(in)  :: rab(n*(n+1)/2)
      real(wp),intent(out) :: logCN(n)
      real(wp),intent(out) :: dlogCN(3,n,n)
      real(wp),intent(in)  :: thr2
      real(wp)             :: cn(n)
      !> counter
      integer  :: i,j,ij,ii
      !> distances
      real(wp) :: r
      real(wp) :: r0
      real(wp) :: dr
      real(wp), dimension(3) :: rij
      !> treshold
      real(wp) :: thr
      !> cn parameter
      real(wp) :: erfCN
      real(wp) :: dlogdcni
      real(wp) :: dlogdcnj
      real(wp) :: derivative
      !> local parameter
      real(wp),parameter :: kn = -7.5_wp
      !> initialize vaiiables
      cn     = 0.0_wp
      logCN  = 0.0_wp
      dlogCN = 0.0_wp
      dlogdcni = 0.0_wp
      dlogdcnj = 0.0_wp
      derivative = 0.0_wp
      r   = 0.0_wp
      r0  = 0.0_wp
      rij = 0.0_wp
      thr = sqrt(thr2)
      !> create error function CN
      do i = 2, n
         ii=i*(i-1)/2
         do j = 1, i-1
            ij = ii+j
            r = rab(ij)
            if (r.gt.thr) cycle
            r0=(param%rcov(at(i))+param%rcov(at(j)))
            dr = (r-r0)/r0
            !> hier kommt die CN funktion hin
!           erfCN = create_expCN(16.0d0,r,r0)
            erfCN = 0.5_wp * (1.0_wp + erf(kn*dr))
            cn(i) = cn(i) + erfCN
            cn(j) = cn(j) + erfCN
         enddo
      enddo
      !> create cutted logarithm CN + derivatives
      do i = 1, n
         ii=i*(i-1)/2
         logCN(i) = create_logCN(cn(i),param)
         !> get dlogCN/dCNi
         dlogdcni = create_dlogCN(cn(i),param)
         do j = 1, i-1
            ij = ii+j
            !> get dlogCN/dCNj
            dlogdcnj = create_dlogCN(cn(j),param)
            r = rab(ij)
            if (r.gt.thr) cycle
            r0 = (param%rcov(at(i)) + param%rcov(at(j)))
            !> get derfCN/dRij
            derivative = create_derfCN(kn,r,r0)
!           derivative = create_dexpCN(16.0d0,r,r0)
            rij  = derivative*(xyz(:,j) - xyz(:,i))/r
            !> project rij gradient onto cartesians
            dlogCN(:,j,j)= dlogdcnj*rij + dlogCN(:,j,j)
            dlogCN(:,i,j)=-dlogdcnj*rij
            dlogCN(:,j,i)= dlogdcni*rij
            dlogCN(:,i,i)=-dlogdcni*rij + dlogCN(:,i,i)
         enddo
      enddo

    contains

pure elemental function create_logCN(cn,param) result(count)
use xtb_mctc_accuracy, only : wp
      type(TGFFData), intent(in) :: param
   real(wp), intent(in) :: cn
   real(wp) :: count
   count = log(1 + exp(param%cnmax)) - log(1 + exp(param%cnmax - cn) )
end function create_logCN

pure elemental function create_dlogCN(cn,param) result(count)
use xtb_mctc_accuracy, only : wp
      type(TGFFData), intent(in) :: param
   real(wp), intent(in) :: cn
   real(wp) :: count
   count = exp(param%cnmax)/(exp(param%cnmax) + exp(cn))
end function create_dlogCN

pure elemental function create_erfCN(k,r,r0) result(count)
use xtb_mctc_accuracy, only : wp
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   real(wp) :: dr
   dr = (r -r0)/r0
   count = 0.5_wp * (1.0_wp + erf(kn*dr))
end function create_erfCN

pure elemental function create_derfCN(k,r,r0) result(count)
use xtb_mctc_accuracy, only : wp
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp), parameter :: sqrtpi = 1.77245385091_wp
   real(wp) :: count
   real(wp) :: dr
   dr = (r -r0)/r0
   count = k/sqrtpi*exp(-k**2*dr*dr)/r0
end function create_derfCN

pure elemental function create_expCN(k,r,r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   count =1.0_wp/(1.0_wp+exp(-k*(r0/r-1.0_wp)))
end function create_expCN

pure elemental function create_dexpCN(k,r,r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   real(wp) :: expterm
   expterm=exp(-k*(r0/r-1._wp))
   count = (-k*r0*expterm)/(r**2*((expterm+1._wp)**2))
end function create_dexpCN

end subroutine gfnff_dlogcoord

end module xtb_gfnff_eg
