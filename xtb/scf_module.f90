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

module scf_module
! ========================================================================
   use iso_fortran_env, wp => real64
   implicit none

   logical,private,parameter :: profile = .true.

contains

subroutine scf(iunit,mol,wfn,basis,param,pcem, &
&              egap,et,maxiter,prlevel,restart,grd,acc,eel,g, &
&              res)

   use mctc_econv, only : autoev,evtoau

! ========================================================================
!  type definitions
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data
   use tbdef_timer
   use tbdef_pcem

! ========================================================================
!  global storage
   use aoparam
   use setparam

! ========================================================================
   use scc_core

! ========================================================================
   use aespot,    only : ovlp2,sdqint,setdqlist,get_radcn,setvsdq, &
   &                     mmomgabzero,mmompop,molmom
   use gbobc,     only : lgbsa,lhb,tb_solvent,gshift, &
   &                     new_gbsa,deallocate_gbsa, &
   &                     update_nnlist_gbsa,compute_fgb,compute_brad_sasa
   use dftd4,     only : build_wdispmat,d4dim,d4,disppot,p_refq_gfn2xtb, &
   &                     mdisp,prmolc6,edisp_scc,edisp,abcappr
   use ncoord,    only : dncoord_gfn,ncoord_d4,dncoord_d3
   use embedding, only : read_pcem,jpot_pcem_gfn1,jpot_pcem_gfn2

   use readin

   implicit none

! ========================================================================
   integer, intent(in)    :: iunit
   type(tb_molecule), intent(in) :: mol
   type(tb_wavefunction),intent(inout) :: wfn
   type(tb_basisset),intent(in) :: basis
   type(scc_parameter),intent(in) :: param
   type(tb_pcem),intent(inout) :: pcem
   real(wp),intent(inout) :: egap
   real(wp),intent(in)    :: et
   integer, intent(in)    :: maxiter
   integer, intent(in)    :: prlevel
   logical, intent(in)    :: restart
   logical, intent(in)    :: grd
   real(wp),intent(in)    :: acc
   real(wp),intent(inout) :: eel
   real(wp),intent(inout) :: g(3,mol%n)
   type(scc_results),intent(out) :: res

! ========================================================================
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: X(:,:)
   real(wp),allocatable :: S(:,:)
   real(wp),allocatable :: S12(:,:)
   real(wp),allocatable :: H0(:)
   real(wp),allocatable :: H1(:)
   real(wp),allocatable :: jab(:,:)
   real(wp),allocatable :: ves(:) ! shell ES potential
   real(wp),allocatable :: tmp(:)
   real(wp),allocatable :: zsh(:)
   real(wp),allocatable :: qq(:)
   real(wp),allocatable :: qlmom(:,:)
   real(wp),allocatable :: cm5(:)
   real(wp),allocatable :: kcnao(:)
   real(wp),allocatable :: Xcao(:,:)
!  AES stuff
   real(wp),allocatable  :: dpint(:,:),qpint(:,:)
   real(wp),allocatable  :: gab3(:),gab5(:)
   real(wp),allocatable  :: vs(:),vq(:,:),vd(:,:)
   real(wp),allocatable  :: gam3sh(:)
   real(wp),allocatable  :: radcn(:) ! CBNEW

! ========================================================================
   real(wp) :: embd
   integer  :: dispdim
   real(wp),allocatable :: c6abns(:,:)
   real(wp),allocatable :: wdispmat(:,:)
   real(wp),allocatable :: hdisp(:)
   real(wp),allocatable :: covcn(:)
   real(wp),allocatable :: gw(:)
   integer  :: mbd
   parameter(mbd=3)
   real(wp) :: molc6,molc8,molpol
   real(wp),allocatable :: aw(:,:)
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: rvol(:)

! ========================================================================
   integer,allocatable :: mdlst(:,:),mqlst(:,:)
   integer :: ndp,nqp

   integer,allocatable :: matlist (:,:)
   integer,allocatable :: matlist2(:,:)
   integer,allocatable :: xblist(:,:)
   real(wp),allocatable :: sqrab(:)
   real(wp),allocatable :: dcn(:,:,:)

   real(wp) :: dipol(3),dip,gsolv,eat,hlgap,efix
   real(wp) :: temp,xsum,eh1,rab,eold,dum,xx,r2
   real(wp) :: t0,t1,sum,rr,scal,hav,alpha,ep,dx,dy,dz,dum1,r0i,r0j
   real(wp) :: efa,efb,nfoda,nfodb,hdii,hdjj,qconv,ff
   real(wp) :: x1,x2,ed,intcut,neglect,ga,gb,ehb,h0s,hmat,rab2
   real(wp) :: h0sr,scfconv,rmsq,dum2,drfdxyz(3),yy,tex,rav,tab,ljexp
   real(wp) :: ees,xa,xb,ya,yb,za,zb,repab,jmpol,kdampxb,vvxb,exb
   real(wp) :: w2,xj,gi,gj,rexp,eatoms,kexp,neglect2,gnorm
   real(wp) :: w0,w1,damp,damp0,gtmp(3),exc
   real(wp) :: d3atm,esave
   real(wp) :: eaes,t6,t7,pi,epol,aot,kexpe
   parameter (pi =  3.14159265358979_wp)

!  some parameter defaults which are not fitted
   data kexp /1.50_wp/ ! rep exp for exp(-R**kexp)
   data rexp /1.00_wp/ ! rep exp in 1/R**rexp
   data d3atm/0.00_wp/ ! ATM scal, zero in GFN1
   data ljexp/12.0_wp/ ! XB parameter for damped LJ in GFN1
   data aot  /-0.5_wp/ ! AO exponent dep. H0 scal

   integer :: ich ! file handle
   integer :: npr,ii,jj,kk,i,j,k,m,iat,jat,mi,jter,atj,kkk,mj,mm
   integer :: ishell,jshell,np,ia,ndimv,l,nmat,nmat2
   integer :: ll,i1,i2,nn,ati,nxb,lin,startpdiag,lladr(4)
   integer :: is,js,lladr2(0:3),tmmetal
   data    lladr  /1,3,6,10/
   data    lladr2 /1,3,5,7/

   character(len=2),external :: asym
   character(len=128) :: atmp,ftmp
   logical :: ex,minpr,pr,fulldiag,xbond,lastdiag,iniqsh,fail,early3d

!  GBSA stuff
   type(tb_solvent) :: gbsa
   real(wp),allocatable :: fgb(:,:)
   real(wp),allocatable :: fhb(:)
   real(wp) :: gborn,tgb
!  for the CM5 charges
   real(wp),allocatable :: cm5a(:)
   real(wp),allocatable :: dcm5a(:,:,:)
   real(wp),allocatable :: fgba(:),dcm5(:,:)
   real(wp) :: hbpow

!  point charge embedding stuff
   logical  :: lpcem
   real(wp),allocatable :: Vpc(:)
   real(wp) :: epcem

   type(tb_timer) :: timer
   character(len=*),parameter :: scifmt = &
      '(10x,":",2x,a,e22.7,1x,a,1x,":")'
   character(len=*),parameter :: dblfmt = &
      '(10x,":",2x,a,f18.7,5x,a,1x,":")'
   character(len=*),parameter :: intfmt = &
      '(10x,":",2x,a,i18,      10x,":")'
   character(len=*),parameter :: chrfmt = &
      '(10x,":",2x,a,a18,      10x,":")'



!  broyden stuff
   logical  :: broy

! ------------------------------------------------------------------------
!  initialization
! ------------------------------------------------------------------------
   if (profile) call timer%new(7,.false.)
   if (profile) call timer%measure(1,"SCC setup")
   rmsq  =1.e+42_wp
   lastdiag=.false.
   lpcem = pcem%n > 0
   jter=0

   dipol= 0.0_wp
   eaes = 0.0_wp
   gsolv= 0.0_wp
   epol = 0.0_wp
   ees  = 0.0_wp
   epcem= 0.0_wp
   ga   = 0.0_wp
   gb   = 0.0_wp
   eel  = 0.0_wp
   ep   = 0.0_wp
   ed   = 0.0_wp
   embd = 0.0_wp
   eat  = 0.0_wp
   egap = 0.0_wp
   molpol = 0.0_wp

   pr   = prlevel.gt.1
   minpr= prlevel.gt.0

!  numerical stuff
!  primitive cut-off
   intcut=25.0_wp-10.0*log10(acc)
   intcut=max(20.0_wp,intcut)
!  integral neglect threshold (the distance criterion used in ovlp limits S to about >1.d-9
   neglect =10.0d-9*acc
   neglect2=neglect*10.0_wp
!  exit if E < scfconv, Hessian sensitive to this
   scfconv=1.d-6*acc
!  old/new q mixing start
   damp0=0.20
!  conv threshold in dRMS q
   qconv=2.d-5*acc
   if(gfn_method.gt.1) qconv=1.d-4*acc

   broy = mol%n > 1
!  when to start pseudodiag (=1 after one full diag)
   startpdiag=1000 !large number=never

!  do the first SCC by full diag
   if(egap.eq.0) startpdiag=1000

!ccccccccccccccccccc
! note: H is in eV!
!ccccccccccccccccccc

!  # atom arrays
   allocate(qq(mol%n),qlmom(3,mol%n),cm5(mol%n),sqrab(mol%n*(mol%n+1)/2),dcn(3,mol%n,mol%n),cn(mol%n))

!  initialize the GBSA module (GBSA works with CM5 charges)
   if(lgbsa) then
      call new_gbsa(gbsa,mol%n,mol%at)
      allocate(fgb(mol%n,mol%n),fhb(mol%n),cm5a(mol%n),dcm5a(3,mol%n,mol%n))
      gborn=0._wp
      gbsa%gsasa=0._wp
      gbsa%ghb=0._wp
      qq=0._wp
!     initialize the neighbor list
      call update_nnlist_gbsa(gbsa,mol%xyz,.false.)
      ! compute Born radii
      call compute_brad_sasa(gbsa,mol%xyz)
!     initialize the fgb matrix (dielectric screening of the Coulomb potential)
      call compute_fgb(gbsa,fgb,fhb)
!     initialize the CM5 charges computation
      if(gfn_method.gt.1) then !GFN2 does not use CM5 charges
        cm5=wfn%q
        cm5a=0.d0
        dcm5a=0.d0
      else
        call prepcm5(mol%n,mol%at,mol%xyz,qq,cm5,cm5a,dcm5a)
      endif
   endif

   allocate(H0(basis%nao*(basis%nao+1)/2), &
   &        S(basis%nao,basis%nao),tmp(basis%nao), &
   &        X(basis%nao,basis%nao),H1(basis%nao*(basis%nao+1)/2), &
   &        kcnao(basis%nao),ves(basis%nshell), &
   &        zsh(basis%nshell),&
   &        jab(basis%nshell,basis%nshell), &
   &        matlist (2,basis%nao*(basis%nao+1)/2), &
   &        matlist2(2,basis%nao*(basis%nao+1)/2-basis%nao))

   call setzshell(mol%n,mol%at,basis%nshell,mol%z,zsh,eatoms,gfn_method)

   ! fill levels
   if(wfn%nel.ne.0) then
      call occu(basis%nao,wfn%nel,wfn%nopen,wfn%ihomoa,wfn%ihomob,wfn%focca,wfn%foccb)
      wfn%focc = wfn%focca + wfn%foccb
      wfn%ihomo=wfn%ihomoa
   else
      wfn%focc=0.0_wp
      wfn%ihomo=0
      wfn%ihomoa=0
      wfn%nopen=0
   endif

!  distances and XB set
   nxb = 0
   k = 0
   do i=1,mol%n
      do j=1,i
!        save distances
         k=k+1
         sqrab(k)=(mol%xyz(1,i)-mol%xyz(1,j))**2 &
         &       +(mol%xyz(2,i)-mol%xyz(2,j))**2 &
         &       +(mol%xyz(3,i)-mol%xyz(3,j))**2
!        XB
         if(xbond(mol%at(i),mol%at(j)).and.sqrab(k).lt.400.) nxb=nxb+1
      enddo
   enddo

   ! XB part
   if(gfn_method.lt.2) then ! not used in GFN2
      allocate(xblist(3,nxb+1))
      nxb = 0
      k = 0
      do i=1,mol%n
         ati=mol%at(i)
         do j=1,i
            atj=mol%at(j)
            k=k+1
            if(xbond(ati,atj).and.sqrab(k).lt.400) then
               nxb=nxb+1
               dum1=1.d+42
               if(ati.eq.17.or.ati.eq.35.or.ati.eq.53.or.ati.eq.85)then
                  xblist(1,nxb)=i
                  xblist(2,nxb)=j
                  do m=1,mol%n
                     if(m.ne.i.and.sqrab(lin(m,i)).lt.dum1)then
                        dum1=sqrab(lin(m,i))
                        xblist(3,nxb)=m
                     endif
                  enddo
               endif
               if(atj.eq.17.or.atj.eq.35.or.atj.eq.53.or.atj.eq.85)then
                  xblist(1,nxb)=j
                  xblist(2,nxb)=i
                  do m=1,mol%n
                     if(m.ne.j.and.sqrab(lin(m,j)).lt.dum1)then
                        dum1=sqrab(lin(m,j))
                        xblist(3,nxb)=m
                     endif
                  enddo
               endif
            endif
         enddo
      enddo
   endif
   ! ldep J potentials (in eV) for SCC
   if(gfn_method.eq.1)then
      call jpot_gfn1(mol%n,basis%nshell,basis%ash,basis%lsh,mol%at,sqrab,param%alphaj,jab)
   else !GFN2
      call jpot_gfn2(mol%n,basis%nshell,basis%ash,basis%lsh,mol%at,sqrab,jab)
   endif

!  J potentials including the point charge stuff
   if(lpcem)then
      allocate( Vpc(basis%nshell), source = 0.0_wp )
      if (gfn_method.eq.1)then
         call jpot_pcem_gfn1(mol%n,pcem,basis%nshell,mol%at,mol%xyz,basis%ash,basis%lsh,param%alphaj,Vpc)
      else ! GFN2
         call jpot_pcem_gfn2(mol%n,pcem,basis%nshell,mol%at,mol%xyz,basis%ash,basis%lsh,Vpc)
      endif
   endif

   ! set 3rd order shell gammas
   if(gfn_method.gt.1) then
      allocate(gam3sh(basis%nshell),source = 0.0_wp)
      do is=1,basis%nshell
         iat=basis%ash(is)
         ati=mol%at(iat)
         dum=param%gam3l(basis%lsh(is))  ! sp or d-pol
         if ((basis%lsh(is).eq.2).and.(tmmetal(ati).ge.1)) dum=param%gam3l(3) ! d-val
         gam3sh(is)=gam3(ati)*dum
      enddo
   endif

   if (prlevel > 1) then
      write(iunit,'(/,10x,51("."))')
      write(iunit,'(10x,":",22x,a,22x,":")') "SETUP"
      write(iunit,'(10x,":",49("."),":")')
      write(iunit,intfmt) "# basis functions  ",basis%nbf
      write(iunit,intfmt) "# atomic orbitals  ",basis%nao
      write(iunit,intfmt) "# shells           ",basis%nshell
      write(iunit,intfmt) "# electrons        ",wfn%nel
      if (gfn_method.eq.1) &
      write(iunit,intfmt) "# halogen bonds    ",nxb
      write(iunit,intfmt) "max. iterations    ",maxiter
      if (gfn_method.eq.2) &
      write(iunit,chrfmt) "Hamiltonian        ","GFN2-xTB"
      if (gfn_method.eq.1) &
      write(iunit,chrfmt) "Hamiltonian        ","GFN1-xTB"
      write(iunit,chrfmt) "restarted?         ",bool2string(restart)
      write(iunit,chrfmt) "GBSA solvation     ",bool2string(lgbsa)
      write(iunit,chrfmt) "PC potential       ",bool2string(lpcem)
      if (lpcem) then
         write(iunit,intfmt) "-> # point charges ",pcem%n
         write(iunit,dblfmt) "-> sum of PC       ",sum(pcem%q),"e   "
      endif
      write(iunit,dblfmt) "electronic temp.   ",et,      "K   "
      write(iunit,dblfmt) "accuracy           ",acc,     "    "
      write(iunit,scifmt) "-> integral cutoff ",intcut,  "    "
      write(iunit,scifmt) "-> integral neglect",neglect, "    "
      write(iunit,scifmt) "-> SCF convergence ",scfconv, "Eh  "
      write(iunit,scifmt) "-> wf. convergence ",qconv,   "e   "
      write(iunit,dblfmt) "Broyden damping    ",broydamp,"    "
      write(iunit,'(10x,51("."))')
   endif

   qq    =wfn%q
   damp  =damp0

   if (profile) call timer%measure(1)
   if (profile) call timer%measure(2,"Dispersion")

   if ((gfn_method.gt.1).and.newdisp) then
      call d4dim(mol%n,mol%at,dispdim)
      allocate( gw(dispdim), &
      &         c6abns(dispdim,dispdim), &
      &         wdispmat(dispdim,dispdim), &
      &         covcn(mol%n), &
      &         hdisp(mol%n), &
      &         source=0.0_wp )
      call ncoord_d4(mol%n,mol%at,mol%xyz,covcn,thr=1600.0_wp)
      call d4(mol%n,dispdim,mol%at,param%wf,param%g_a,param%g_c,covcn,gw,c6abns)
      call build_wdispmat(mol%n,dispdim,mol%at,mol%xyz,param%disp,c6abns,gw, &
      &                    wdispmat)
   else
      allocate( hdisp(mol%n), source=0.0_wp )
      ! D3 part first because we need CN
      call dncoord_d3(mol%n,mol%at,mol%xyz,cn,dcn)
   endif

   if (profile) call timer%measure(2)
   if (profile) call timer%measure(3,"integral evaluation")

   allocate(dpint(3,basis%nao*(basis%nao+1)/2), &
      &     qpint(6,basis%nao*(basis%nao+1)/2), &
      &     source = 0.0_wp)
   ! compute integrals and prescreen to set up list arrays
   call sdqint(mol%n,mol%at,basis%nbf,basis%nao,mol%xyz,neglect,ndp,nqp,intcut, &
      &        basis%caoshell,basis%saoshell,basis%nprim,basis%primcount, &
      &        basis%alp,basis%cont,S,dpint,qpint)

   ! prepare aes stuff
   if(gfn_method.gt.1) then
!     CN/dCN replaced by special smoother and faster decaying function
      call dncoord_gfn(mol%n,mol%at,mol%xyz,cn,dcn)

!     allocate arrays for lists and fill (to exploit sparsity)
      allocate(mdlst(2,ndp),mqlst(2,nqp))
      call setdqlist(basis%nao,ndp,nqp,neglect,dpint,qpint,mdlst,mqlst)
!     set up 1/R^n * damping function terms
      ii=mol%n*(mol%n+1)/2
      allocate(gab3(ii),gab5(ii),radcn(mol%n))
      call get_radcn(mol%n,mol%at,cn,param%cn_shift,param%cn_expo,param%cn_rmax,radcn)
      call mmomgabzero(mol%n,mol%at,mol%xyz,param%xbrad,param%xbdamp,radcn,gab3,gab5) ! zero damping, xbrad=kdmp3,xbdamp=kdmp5
!     allocate CAMM arrays
   endif

   if (profile) call timer%measure(3)
   if (profile) call timer%measure(1)

   if(gfn_method.gt.1) then
!     if no CAMMs were read, get them from P (e.g., in geometry opts., MD runs)
      if(.not.restart) &
      &  call mmompop(mol%n,basis%nao,basis%aoat2,mol%xyz,wfn%p,s,dpint,qpint, &
      &               wfn%dipm,wfn%qp)

!     scale CAMMs before setting up the potential
!     call scalecamm(mol%n,mol%at,dipm,qp)
!     compute intermediates for potential
      allocate(vs(mol%n),vq(6,mol%n),vd(3,mol%n))
      vs=0.0_wp
      vd=0.0_wp
      vq=0.0_wp
      call setvsdq(mol%n,mol%at,mol%xyz,wfn%q,wfn%dipm,wfn%qp,gab3,gab5,vs,vd,vq)
   endif

   if (gfn_method.gt.1) &
   & call disppot(mol%n,dispdim,mol%at,wfn%q,param%g_a,param%g_c,wdispmat,gw,hdisp)

   if(lgbsa) cm5=wfn%q+cm5a

   ! set up first ES potential
   if(lpcem) then
      ves(1:basis%nshell)=Vpc(1:basis%nshell)
   else
      ves=0.0_wp
   endif
   call setespot(basis%nshell,wfn%qsh,jab,ves)

   ! prepare matrix indices
   nmat =0
   nmat2=0
   do ii=1,basis%nao
      iat=basis%aoat2(ii)
      do jj=1,ii-1
         jat=basis%aoat2(jj)
         if(abs(S(jj,ii)).lt.neglect) then
            S(jj,ii)=0.0_wp
            S(ii,jj)=0.0_wp
            cycle
         endif
         nmat=nmat+1
         matlist(1,nmat)=ii
         matlist(2,nmat)=jj
         if(iat.ne.jat)then
            nmat2=nmat2+1
            matlist2(1,nmat2)=ii
            matlist2(2,nmat2)=jj
         endif
      enddo
      ! CB: moved this here so j/i indices from matlist come in a reasonable order
      ! also setup CN dep. stuff
      nmat=nmat+1
      matlist(1,nmat)=ii
      matlist(2,nmat)=ii
      ishell=mmm(basis%lao2(ii))
      if(gfn_method.lt.2) then
         kcnao(ii)=param%kcnsh(ishell)
         if(metal(mol%at(iat)).eq.1) kcnao(ii)=0.0_wp  ! CN dep. bad for metals
         if(early3d(mol%at(iat))) then
            kcnao(ii)=param%kcnsh(ishell)
            if(ishell.eq.3) kcnao(ii)=param%kcnsh(4) ! fix problems with too low-coord CP rings
         endif
      else
         kcnao(ii)=kcnat(ishell-1,mol%at(iat))  ! clean GFN2 version
      endif
   enddo

   if (profile) call timer%measure(1)
   if (profile) call timer%measure(6,"classical contributions")

   ! this is the classical part of the energy/gradient
   ! dispersion/XB/repulsion for GFN1-xTB
   ! only repulsion for GFN2-xTB
   call cls_grad(mol%n,mol%at,mol%xyz,sqrab,param,rexp,kexp,nxb,ljexp,xblist, &
      &          ed,exb,ep,g,prlevel)

   if (profile) call timer%measure(6)
   if (profile) call timer%measure(4,"zeroth order Hamiltonian")

   if(pr)then
      write(iunit,'(a)')
      write(iunit,*) 'iter      E             dE          RMSdq', &
      &'      gap      omega  full diag'
   endif

! ========================================================================
!  do H0 once
! ========================================================================
   H0=0
   if(gfn_method.eq.1)then
      call build_h0_gfn1(H0,mol%n,mol%at,basis%nao,nmat,matlist, &
      &                  param%kspd,param%kmagic,param%kenscal, &
      &                  mol%xyz,cn,kcnao,S,basis%aoat2,basis%lao2,basis%valao2,basis%hdiag2)
   else
      call build_h0_gfn2(H0,mol%n,mol%at,basis%nao,nmat,matlist, &
      &                  param%kspd,param%kmagic,param%kenscal, &
      &                  mol%xyz,cn,kcnao,S,basis%aoat2,basis%lao2,basis%valao2,basis%hdiag2,basis%aoexp)
   endif
! ========================================================================

   ! first order energy for given geom. and density, i.e. skip SCC and grad
   if(maxiter.eq.0) then
      call qsh2qat(mol%n,mol%at,basis%nshell,wfn%qsh,wfn%q)
      if(gfn_method.gt.1) then
         call electro2(mol%n,mol%at,basis%nao,basis%nshell,jab,H0,wfn%P, &
         &             wfn%q,gam3sh,wfn%qsh,param%gscal,ees,eel)
      else
         call electro(mol%n,mol%at,basis%nao,basis%nshell,jab,H0,wfn%P,wfn%q,wfn%qsh,ees,eel)
      endif
      if(lgbsa) then
         cm5=wfn%q+cm5a
         call electro_gbsa(mol%n,mol%at,fgb,fhb,cm5,gborn,eel)
      endif
      goto 9999
   endif
   if (profile) call timer%measure(4)

! ========================================================================
!  SCC iterations
! ========================================================================
   if (profile) call timer%measure(5,"iterations")
   if (gfn_method.eq.1) then
      call scc_gfn1(iunit,mol%n,wfn%nel,wfn%nopen,basis%nao,nmat,basis%nshell, &
      &             mol%at,matlist,basis%aoat2,basis%ao2sh, &
      &             wfn%q,qq,qlmom,wfn%qsh,zsh, &
      &             gbsa,fgb,fhb,cm5,cm5a,gborn, &
      &             broy,broydamp,damp0, &
      &             lpcem,ves,vpc, &
      &             et,wfn%focc,wfn%focca,wfn%foccb,wfn%efa,wfn%efb, &
      &             eel,ees,epcem,egap,wfn%emo,wfn%ihomo,wfn%ihomoa,wfn%ihomob, &
      &             H0,H1,wfn%C,S,X,wfn%P,jab, &
      &             maxiter,startpdiag,scfconv,qconv, &
      &             minpr,pr, &
      &             fail,jter)
   else
      call scc_gfn2(iunit,mol%n,wfn%nel,wfn%nopen,basis%nao,ndp,nqp,nmat,basis%nshell, &
      &             mol%at,matlist,mdlst,mqlst,basis%aoat2,basis%ao2sh, &
      &             wfn%q,wfn%dipm,wfn%qp,qq,qlmom,wfn%qsh,zsh, &
      &             mol%xyz,vs,vd,vq,gab3,gab5,param%gscal, &
      &             gbsa,fgb,fhb,cm5,cm5a,gborn, &
      &             newdisp,dispdim,param%g_a,param%g_c,gw,wdispmat,hdisp, &
      &             broy,broydamp,damp0, &
      &             lpcem,ves,vpc, &
      &             et,wfn%focc,wfn%focca,wfn%foccb,wfn%efa,wfn%efb, &
      &             eel,ees,eaes,epol,ed,epcem,egap, &
      &             wfn%emo,wfn%ihomo,wfn%ihomoa,wfn%ihomob, &
      &             H0,H1,wfn%C,S,dpint,qpint,X,wfn%P,jab,gam3sh, &
      &             maxiter,startpdiag,scfconv,qconv, &
      &             minpr,pr, &
      &             fail,jter)
   endif
! ========================================================================
   ! free some memory (this stuff is not needed for gradients)
   deallocate(ves)
   if (allocated(vpc)) deallocate(vpc)
   if(allocated(wdispmat)) deallocate( wdispmat )

   9999  continue

! ------------------------------------------------------------------------
!  check for convergence, only do this if printlevel is maximal (WHY?)
   res % converged = .not. fail
   if (pr) then
      if (fail) then
         call touch_file('.sccnotconverged')
         call raise('S',"SCC is not converged properly!",1)
         write(iunit,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***")') &
            "convergence criteria cannot be satisfied within",jter,"iterations"
      else
         write(iunit,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***")') &
            "convergence criteria satisfied after",jter,"iterations"
      endif
   endif

   if (profile) call timer%measure(5)
   if (.not.pr.and.profile.and.minpr) &
      call timer%write_timing(iunit,5,'SCC iter.')
   if (profile) call timer%measure(6,"molecular gradient")

!  print'("Entering gradient calculation")'
! ========================================================================
!  GRADIENT (now 100% analytical (that's not true!))
! ========================================================================

   call scf_grad(mol%n,mol%at,nmat2,matlist2, &
        &        H0,H1,S, &
        &        mol%xyz,sqrab,wfn,basis, &
        &        param,kcnao, &
        &        dispdim,c6abns,mbd, &
        &        intcut, &
        &        gab3,gab5,radcn, &
        &        lpcem,pcem, &
        &        gbsa,gborn,fgb,fhb,cm5a,dcm5a,gsolv, &
        &        eel,ed,embd, &
        &        g,prlevel)
!  print'("Finished gradient calculation")'

!  calculate the norm for printout
   gnorm = sqrt(sum( g**2 ))
! ========================================================================
!  clear some space
   if (gfn_method.gt.1) then
      deallocate(radcn,gab3,gab5,mdlst,mqlst,vs,vd,vq,gam3sh) !CBNEW
   endif

   if (profile) call timer%measure(6)
   if (.not.pr.and.profile.and.minpr) &
      call timer%write_timing(iunit,6,'gradient')
   if (profile) call timer%measure(7,"printout")

! ========================================================================
!  PROPERTIES & PRINTOUT
! ========================================================================
   printing: if (pr) then
! ------------------------------------------------------------------------
!     print orbital energies and occupation numbers
      if (pr_eig) then
         !call preig(iunit,wfn%focc,1.0_wp,wfn%emo, &
                    !max(wfn%ihomoa-12,1),min(wfn%ihomoa+11,basis%nao))
         call print_orbital_eigenvalues(iunit,wfn,5)
      endif

! ------------------------------------------------------------------------
!     get Wiberg bond orders
      call wiberg(mol%n,basis%nao,mol%at,mol%xyz,wfn%P,S,wfn%wbo,.false.,.false., &
         &        basis%fila2)

! ------------------------------------------------------------------------
!     HOMO-LUMO excitation properties if  UHF=2        
      if (wfn%nopen.eq.2) then
         call hlex(mol%n,mol%at,basis%nbf,basis%nao,wfn%ihomoa,mol%xyz,wfn%focc,S,wfn%C,wfn%emo,basis)
      endif

! ------------------------------------------------------------------------
!     LMO /xTB-IFF
      if (pr_lmo) then
         tmp=wfn%emo*evtoau
         call local(mol%n,mol%at,basis%nbf,basis%nao,wfn%ihomoa,mol%xyz,mol%z,wfn%focc,S,wfn%P,wfn%C,tmp,wfn%q,eel,lgbsa,basis)
      endif

! ------------------------------------------------------------------------
!  exchange energy correction ala sTDA
      if (wfn%nopen.ge.2) then
         call exch(mol%n,mol%at,basis%nao,wfn%nopen,wfn%ihomoa,mol%xyz,wfn%focc,S,wfn%C,exc,basis%aoat)
         write(iunit,'(''open-shell EX :'',F16.7)') -exc
         write(iunit,'(''corrected Etot:'',F16.7, &
         &   '' (not used further except for this printout!)'')') eel - exc
      endif

   endif printing

! ------------------------------------------------------------------------
!  dipole calculation (always done because its free)
   if (gfn_method.lt.2) then
      call mmompop(mol%n,basis%nao,basis%aoat2,mol%xyz,wfn%p,s,dpint,qpint, &
         &         wfn%dipm,wfn%qp)
   endif

   call calc_dipole(mol%n,mol%at,mol%xyz,mol%z,basis%nao,wfn%P,dpint,dip,dipol)

   if (profile) call timer%measure(7)

!  END OF PROPERTY & PRINTOUT BLOCK
! ========================================================================

! ========================================================================
!  SAVE FOR FINAL PRINTOUT
! ========================================================================
   dum=eel
!  Etot
   eel = eel + ep + exb + embd
   eat = eatoms*evtoau - eel
   if(gfn_method.lt.2) then
      eel = eel + ed
   endif
   res%e_elec  = dum
   res%e_atom  = eat
   res%e_rep   = ep
   res%e_es    = ees
   res%e_aes   = eaes
   res%e_axc   = epol
   res%e_disp  = ed+embd
   res%e_total = eel
   res%hl_gap  = egap
   res%dipole  = dipol
   if (lgbsa) then
      res%g_solv  = gsolv
      res%g_born  = gborn
      res%g_sasa  = gbsa%gsasa
      res%g_hb    = gbsa%ghb
      res%g_shift = gshift
   endif
   res%gnorm = norm2(g)

   if (profile.and.pr) call timer%write(iunit,'SCC')

! ========================================================================
   if (profile) call timer%deallocate

   deallocate(S,H0,tmp,X,H1,kcnao,zsh,jab,matlist,matlist2,dpint)
   call deallocate_gbsa(gbsa)
end subroutine scf

subroutine scf_grad(n,at,nmat2,matlist2, &
      &             H0,H1,S, &
      &             xyz,sqrab,wfn,basis, &
      &             param,kcnao, &
      &             dispdim,c6abns,mbd, &
      &             intcut, &
      &             gab3,gab5,radcn, &
      &             lpcem,pcem, &
      &             gbsa,gborn,fgb,fhb,cm5a,dcm5a,gsolv, &
      &             eel,ed,embd, &
      &             g,printlvl)

! ========================================================================
!  type definitions
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_pcem

! ========================================================================
!  global storage
   use aoparam
   use setparam

! ========================================================================
!  interfaces
   use scc_core
   use grad_core

   use aespot,    only : ddqint,dradcn,aniso_grad,setdvsdq,dsint
   use gbobc,     only : lgbsa,lhb,tb_solvent,gshift,compute_gb_egrad
   use dftd4,     only : dispgrad,d4_numgrad
   use ncoord,    only : dncoord_gfn,dncoord_d3
   use embedding, only : pcem_grad_gfn1,pcem_grad_gfn2

   implicit none

   type(tb_wavefunction),intent(in) :: wfn
   type(tb_basisset),    intent(in) :: basis
   type(scc_parameter),  intent(in) :: param
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   integer, intent(in)    :: nmat2
   integer,intent(in) :: matlist2(2,nmat2)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(in)    :: sqrab(n*(n+1)/2)
   real(wp),intent(inout) :: H0(basis%nao*(basis%nao+1)/2)
   real(wp),intent(inout) :: H1(basis%nao*(basis%nao+1)/2)
   real(wp),intent(inout) :: g(3,n)
   real(wp),intent(in)    :: S(basis%nao,basis%nao)
   real(wp),intent(in)    :: kcnao(basis%nao)
   real(wp),intent(in)    :: gab3(n*(n+1)/2)
   real(wp),intent(in)    :: gab5(n*(n+1)/2)
   real(wp),intent(in)    :: intcut
   real(wp),intent(in)    :: radcn(n)
   integer, intent(in)    :: dispdim
   real(wp),intent(in)    :: c6abns(dispdim,dispdim)
   integer, intent(in)    :: mbd
   real(wp),intent(inout) :: ed
   real(wp),intent(out)   :: embd
   real(wp),intent(inout) :: eel
   type(tb_solvent),intent(inout) :: gbsa
   real(wp),intent(inout) :: gborn
   real(wp)               :: ghb
   real(wp),intent(inout) :: gsolv
   real(wp),intent(in)    :: cm5a(n)
   real(wp),intent(in)    :: dcm5a(3,n,n)
   real(wp),intent(inout) :: fgb(n,n)
   real(wp),intent(inout) :: fhb(n)
   type(tb_pcem),intent(inout) :: pcem
   logical, intent(in)    :: lpcem
   integer, intent(in)    :: printlvl

   integer :: m,i,j,kk
   real(wp),allocatable :: qq(:)
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: dcn(:,:,:)
   real(wp),allocatable :: X(:,:)
   real(wp),allocatable :: H(:,:)
   real(wp),allocatable :: vs(:),vd(:,:),vq(:,:)
   logical :: pr,minpr

!  print'("Entered gradient calculation")'

   minpr = printlvl.gt.0
   pr = printlvl.gt.1

!  print'("Allocating local memory")'
   allocate( cn(n), source = 0.0_wp )
   allocate( dcn(3,n,n), source = 0.0_wp )
   allocate( H(basis%nao,basis%nao), source = 0.0_wp )
   allocate( X(basis%nao,basis%nao), source = 0.0_wp )
!  print'("Allocated local memory")'

!  get energy weighted density matrix and convert from eV to Eh
!  print'("Getting energy weighted density matrix")'
   call prep_grad_conv(basis%nao,H0,H1,wfn%C,wfn%focc,wfn%emo,X)

!  wave function terms
!  print'("Calculating polynomial derivatives")'
   call poly_grad(g,n,at,basis%nao,nmat2,matlist2,xyz,sqrab,wfn%P,S,basis%aoat2,basis%lao2,H0)

!  CN dependent part
!  print'("Calculating CN dependent derivatives")'
   if (gfn_method.gt.1) then
      call dncoord_gfn(n,at,xyz,cn,dcn)
      call hcn_grad_gfn2(g,n,at,basis%nao,nmat2,matlist2,xyz, &
           &             param%kspd,param%kmagic,param%kenscal,kcnao,wfn%P,S,dcn, &
           &             basis%aoat2,basis%lao2,basis%valao2,basis%hdiag2,basis%aoexp)
   else
      call dncoord_d3(n,at,xyz,cn,dcn)
      call hcn_grad_gfn1(g,n,at,basis%nao,nmat2,matlist2,xyz, &
           &             param%kspd,param%kmagic,param%kenscal,kcnao,wfn%P,S,dcn, &
           &             basis%aoat2,basis%lao2,basis%valao2,basis%hdiag2)
   endif

!  preccalc
!  print'("Resetting the Hamiltonian")'
   do m=1,nmat2
      i=matlist2(1,m)
      j=matlist2(2,m)
      kk=j+i*(i-1)/2
      H(j,i)=(H1(kk)+H0(kk))*wfn%P(j,i)/S(j,i)-X(j,i)
      H(i,j)=H(j,i)
   enddo

!  multipole gradient stuff
!  print'("Calculating multipole gradient")'
   if (gfn_method.gt.1) then
      allocate( vs(n),vd(3,n),vq(6,n), source = 0.0_wp )
!     VS, VD, VQ-dependent potentials are changed w.r.t. SCF,
!     since moment integrals are now computed with origin at
!     respective atoms
      call setdvsdq(n,at,xyz,wfn%q,wfn%dipm,wfn%qp,gab3,gab5,vs,vd,vq)
      call ddqint(intcut,n,basis%nao,basis%nbf,at,xyz, &
         &        basis%caoshell,basis%saoshell,basis%nprim,basis%primcount, &
         &        basis%alp,basis%cont,wfn%p,vs,vd,vq,H,g)

! WARNING: dcn is overwritten on output and now dR0A/dXC, 
!          and index i & j are flipped
      call dradcn(n,at,cn,param%cn_shift,param%cn_expo,param%cn_rmax,dcn)
      call aniso_grad(n,at,xyz,wfn%q,wfn%dipm,wfn%qp,param%xbrad,param%xbdamp, &
           &          radcn,dcn,gab3,gab5,g)

   else
!     wave function terms 2/overlap dependent parts of H
      call dsint(intcut,n,basis%nao,basis%nbf,at,xyz,sqrab, &
         &        basis%caoshell,basis%saoshell,basis%nprim,basis%primcount, &
         &        basis%alp,basis%cont,H,g)
   endif

!  dispersion (DFT-D type correction)
!  print'("Calculating dispersion gradient")'
   if ((gfn_method.gt.1).and.newdisp) then
      call dispgrad(n,dispdim,at,wfn%q,xyz, &
           &        param%disp,param%wf,param%g_a,param%g_c, &
           &        c6abns,mbd,g,embd)
      embd = embd-ed
   endif

! GBSA
! start GBSA gradient
!  print'("Calculating GBSA gradient")'
   if (lgbsa) then
      if (gfn_method.gt.1) then
         call compute_gb_egrad(gbsa,wfn%q,gborn,ghb,g,minpr)
      else
         allocate( qq(n), source = wfn%q )
         qq=qq+cm5a
         call compute_gb_egrad(gbsa,qq,gborn,ghb,g,minpr)
         call cm5_grad_gfn1(g,n,qq,fgb,fhb,dcm5a,lhb)
      endif
!     solvation energy
      gbsa%gborn = gborn
      gbsa%ghb = ghb
      gsolv=gborn+gbsa%gsasa+gbsa%ghb+gshift
      eel=eel+gbsa%gsasa+gshift
   endif

!  print'("Calculating shell es and repulsion gradient")'
   if(gfn_method.eq.1)then
      call shelles_grad_gfn1(g,n,at,basis%nshell,xyz,sqrab, &
         &                 basis%ash,basis%lsh,param%alphaj,wfn%qsh)
   else ! GFN2
      call shelles_grad_gfn2(g,n,at,basis%nshell,xyz,sqrab,basis%ash,basis%lsh,wfn%qsh)
   endif

! --- ES point charge embedding
!  print'("Calculating embedding gradient")'
   if (lpcem) then
      if (gfn_method.eq.1) then
         call pcem_grad_gfn1(g,pcem%grd,n,pcem,at,basis%nshell,xyz, &
            &                basis%ash,basis%lsh,param%alphaj,wfn%qsh)
      else
         call pcem_grad_gfn2(g,pcem%grd,n,pcem,at,basis%nshell,xyz, &
            &                basis%ash,basis%lsh,wfn%qsh)
      endif
   endif

!  print'("Finished in gradient subroutine")'

end subroutine scf_grad

subroutine cls_grad(n,at,xyz,sqrab, &
      &             param,rexp,kexp, &
      &             nxb,ljexp,xblist, &
      &             ed,exb,ep, &
      &             g,printlvl)

! ========================================================================
!  type definitions
   use tbdef_param

! ========================================================================
!  global storage
   use aoparam
   use setparam

! ========================================================================
!  interfaces
   use scc_core
   use grad_core

   implicit none

   type(scc_parameter),  intent(in) :: param
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(in)    :: sqrab(n*(n+1)/2)
   real(wp),intent(inout) :: g(3,n)
   real(wp),intent(inout) :: ed
   integer, intent(in)    :: nxb
   integer, intent(in)    :: xblist(3,nxb+1)
   real(wp),intent(in)    :: ljexp
   real(wp),intent(inout) :: exb
   real(wp),intent(inout) :: ep
   real(wp),intent(inout) :: rexp
   real(wp),intent(inout) :: kexp
   integer, intent(in)    :: printlvl

   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: dcn(:,:,:)
   logical :: pr,minpr

!  print'("Entered gradient calculation")'

   minpr = printlvl.gt.0
   pr = printlvl.gt.1

!  print'("Allocating local memory")'
   allocate( cn(n), source = 0.0_wp )
   allocate( dcn(3,n,n), source = 0.0_wp )

!  dispersion (DFT-D type correction)
!  print'("Calculating dispersion gradient")'
   if (gfn_method.eq.1) then
      call gdisp(n,at,xyz,param%disp%a1,param%disp%a2,param%disp%s8,param%disp%s9, &
      &   ed,g,cn,dcn)
   endif

! XB or gCP ----------------------------------------------------
!  print'("Calculating xbond gradient")'
   exb=0.0_wp
   if(gfn_method.lt.2) then
      call xbpot(n,at,xyz,sqrab,xblist,nxb,param%xbdamp,param%xbrad,ljexp,exb,g)
   endif

!  print'("Calculating shell es and repulsion gradient")'
   if(gfn_method.eq.1)then
      call rep_grad_gfn1(g,ep,n,at,xyz,sqrab,kexp,rexp)
   else ! GFN2
      call rep_grad_gfn2(g,ep,n,at,xyz,sqrab,rexp)
   endif

end subroutine cls_grad

end module scf_module
