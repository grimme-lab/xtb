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

module xtb_scf
! ========================================================================
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autoev,evtoau
   use xtb_mctc_la, only : contract
   use xtb_coulomb_klopmanohno
   use xtb_disp_coordinationnumber, only : getCoordinationNumber, cnType
   use xtb_disp_dftd3, only : d3_gradient
   use xtb_type_basisset
   use xtb_type_coulomb, only : TCoulomb
   use xtb_type_data
   use xtb_type_environment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_timer
   use xtb_type_wavefunction
   use xtb_xtb_data
   use xtb_xtb_halogen
   use xtb_xtb_repulsion
   use xtb_xtb_coulomb
   use xtb_xtb_dispersion
   use xtb_xtb_hamiltonian, only : getSelfEnergy, build_SDQH0, build_dSDQH0, &
      & count_dpint, count_qpint
   use xtb_xtb_multipole
   use xtb_paramset, only : tmmetal
   use xtb_scc_core
   use xtb_grad_core
   use xtb_hlex
   use xtb_local
   use xtb_dipole
   implicit none
   private

   public :: scf

   logical, parameter :: profile = .true.

   integer, parameter :: mmm(*)=(/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/)

   integer, parameter :: metal(1:86) = [&
      & 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
      & 1, 0, 0, 0, 0, 0]

contains

subroutine scf(env, mol, wfn, basis, pcem, xtbData, &
      & egap, et, maxiter, prlevel, restart, grd, acc, energy, gradient, res)


! ========================================================================
!  global storage
   use xtb_setparam

! ========================================================================
! ========================================================================
   use xtb_aespot,    only : setdqlist,get_radcn,setvsdq, &
   &                     mmomgabzero,mmompop,molmom
   use xtb_solv_gbobc,     only : lgbsa,lhb,TSolvent,gshift, &
   &                     new_gbsa,deallocate_gbsa, &
   &                     update_nnlist_gbsa,compute_fgb,compute_brad_sasa
   use xtb_disp_dftd4, only: build_wdispmat,d4dim,d4,disppot,p_refq_gfn2xtb, &
   &                     mdisp,prmolc6,edisp_scc,edisp,abcappr,d4_gradient
   use xtb_disp_ncoord,    only : dncoord_gfn,ncoord_d4,dncoord_d3
   use xtb_embedding, only : read_pcem,jpot_pcem_gfn1,jpot_pcem_gfn2
   use xtb_aespot,    only : dradcn,aniso_grad,setdvsdq,dsint
   use xtb_solv_gbobc,     only : lgbsa,lhb,TSolvent,gshift,compute_gb_egrad
   use xtb_disp_dftd4, only: dispgrad
   use xtb_disp_ncoord,    only : dncoord_gfn,dncoord_d3
   use xtb_embedding, only : pcem_grad_gfn1,pcem_grad_gfn2

   use xtb_readin

   implicit none

   character(len=*), parameter :: source = 'scf'

! ========================================================================
   type(TEnvironment), intent(inout)    :: env
   type(TMolecule), intent(in) :: mol
   type(TWavefunction),intent(inout) :: wfn
   type(TBasisset),intent(in) :: basis
   type(tb_pcem),intent(inout) :: pcem
   real(wp),intent(inout) :: egap
   real(wp),intent(in)    :: et
   integer, intent(in)    :: maxiter
   integer, intent(in)    :: prlevel
   logical, intent(in)    :: restart
   logical, intent(in)    :: grd
   real(wp),intent(in)    :: acc
   real(wp),intent(inout) :: energy
   real(wp),intent(inout) :: gradient(3,mol%n)
   real(wp) :: sigma(3,3)
   type(scc_results),intent(out) :: res
   type(TxTBData), intent(in) :: xtbData

! ========================================================================
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: X(:,:)
   real(wp),allocatable :: S(:,:)
   real(wp),allocatable :: S12(:,:)
   real(wp),allocatable :: H0(:)
   real(wp),allocatable :: ves(:) ! shell ES potential
   real(wp),allocatable :: tmp(:)
   real(wp),allocatable :: zsh(:)
   real(wp),allocatable :: qq(:)
   real(wp),allocatable :: qlmom(:,:)
   real(wp),allocatable :: cm5(:)
   real(wp),allocatable :: Xcao(:,:)
   real(wp),allocatable :: selfEnergy(:, :)
   real(wp),allocatable :: dSEdcn(:, :)
   integer :: nid
   integer, allocatable :: idnum(:)
   type(TxTBCoulomb) :: ies
   type(TKlopmanOhno) :: coulomb
   real(wp), allocatable :: djdr(:, :, :)
   real(wp), allocatable :: djdtr(:, :)
   real(wp), allocatable :: djdL(:, :, :)
!  AES stuff
   type(TxTBMultipole), allocatable :: aes
   real(wp),allocatable  :: dpint(:,:,:),qpint(:,:,:)
   real(wp),allocatable  :: radcn(:) ! CBNEW

! ========================================================================
   type(TxTBDispersion), allocatable :: scD4
   real(wp) :: embd
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
   real(wp),allocatable :: dcndr(:,:,:)
   real(wp),allocatable :: dcndL(:,:,:)
   real(wp),allocatable :: trans(:,:)

   real(wp) :: dipol(3),dip,gsolv,eat,hlgap,efix
   real(wp) :: temp,xsum,eh1,rab,eold,dum,xx,r2
   real(wp) :: t0,t1,sum,rr,hav,alpha,ep,dx,dy,dz,dum1,r0i,r0j
   real(wp) :: efa,efb,nfoda,nfodb,hdii,hdjj,qconv,ff
   real(wp) :: x1,x2,ed,intcut,neglect,ga,gb,ehb,h0s,hmat,rab2
   real(wp) :: h0sr,scfconv,rmsq,dum2,drfdxyz(3),yy,tex,rav,tab,ljexp
   real(wp) :: ees,xa,xb,ya,yb,za,zb,repab,jmpol,kdampxb,vvxb,exb
   real(wp) :: w2,xj,gi,gj,eatoms,neglect2
   real(wp) :: w0,w1,damp,damp0,gtmp(3),exc
   real(wp) :: esave,eel
   real(wp) :: eaes,t6,t7,pi,epol,kexpe
   parameter (pi =  3.14159265358979_wp)

!  some parameter defaults which are not fitted
   data ljexp/12.0_wp/ ! XB parameter for damped LJ in GFN1

   integer :: ich ! file handle
   integer :: npr,ii,jj,kk,i,j,k,m,iat,jat,mi,jter,atj,kkk,mj,mm
   integer :: ishell,jshell,np,ia,ndimv,l,nmat,nmat2
   integer :: ll,i1,i2,nn,ati,nxb,lin,startpdiag
   integer :: is,js,gav

   character(len=128) :: atmp,ftmp
   logical :: ex,minpr,pr,fulldiag,lastdiag,iniqsh,fail

!  GBSA stuff
   type(TSolvent) :: gbsa
   real(wp),allocatable :: fgb(:,:)
   real(wp),allocatable :: fhb(:)
   real(wp) :: gborn,ghb
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

   logical :: exitRun

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
   if(allocated(xtbData%multipole)) qconv=1.d-4*acc

   broy = mol%n > 1
!  when to start pseudodiag (=1 after one full diag)
   startpdiag=1000 !large number=never

!  do the first SCC by full diag
   if(egap.eq.0) startpdiag=1000

!ccccccccccccccccccc
! note: H is in eV!
!ccccccccccccccccccc
   !> update atomic Mulliken charges
   call qsh2qat(basis%ash,wfn%qsh,wfn%q)

!  # atom arrays
   allocate(qq(mol%n),qlmom(3,mol%n),cm5(mol%n),sqrab(mol%n*(mol%n+1)/2))
   allocate(dcndr(3,mol%n,mol%n),cn(mol%n),dcndL(3,3,mol%n))
   allocate(trans(3,1))
   trans(:, :) = 0.0_wp

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
      if(gfn_method > 1) then !GFN2 does not use CM5 charges
        cm5=wfn%q
        cm5a=0.d0
        dcm5a=0.d0
      else
        call calc_cm5(mol%n,mol%at,mol%xyz,cm5a,dcm5a)
        cm5 = wfn%q + cm5a
      endif
   endif

   allocate(H0(basis%nao*(basis%nao+1)/2), &
   &        S(basis%nao,basis%nao),tmp(basis%nao), &
   &        X(basis%nao,basis%nao), &
   &        ves(basis%nshell), &
   &        zsh(basis%nshell),&
   &        matlist (2,basis%nao*(basis%nao+1)/2), &
   &        matlist2(2,basis%nao*(basis%nao+1)/2-basis%nao))
   allocate(selfEnergy(maxval(xtbData%nshell), mol%n))
   allocate(dSEdcn(maxval(xtbData%nshell), mol%n))

   call setzshell(xtbData,mol%n,mol%at,basis%nshell,mol%z,zsh,eatoms,gfn_method)

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
   if(allocated(xtbData%halogen)) then ! not used in GFN2
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

   ! setup isotropic electrostatics
   call init(ies, xtbData%coulomb, xtbData%nshell, mol%at)

   nid = maxval(mol%id)
   allocate(idnum(nid))
   do ii = 1, nId
      jat = 0
      do iat = 1, mol%n
         if (mol%id(iat) == ii) then
            jat = iat
            exit
         end if
      end do
      idnum(ii) = mol%at(jat)
   end do

   if(gfn_method == 1)then
      gav = gamAverage%harmonic
   else !GFN2
      gav = gamAverage%arithmetic
   endif
   call init(coulomb, env, mol, gav, xtbData%coulomb%shellHardness, &
      & xtbData%coulomb%gExp, num=idnum, nshell=xtbData%nShell)

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Setup of Coulomb evaluator failed", source)
      return
   end if

   call coulomb%getCoulombMatrix(mol, ies%jmat)

!  J potentials including the point charge stuff
   allocate(Vpc(basis%nshell))
   vpc(:) = 0.0_wp
   if(lpcem)then
      if (gfn_method == 1)then
         call jpot_pcem_gfn1(xtbData%coulomb,mol%n,pcem,xtbData%nshell,mol%at, &
            & mol%xyz,xtbData%coulomb%gExp,Vpc)
      else ! GFN2
         call jpot_pcem_gfn2(xtbData%coulomb,mol%n,pcem,xtbData%nshell,mol%at, &
            & mol%xyz,Vpc)
      endif
   endif

   if (prlevel > 1) then
      write(env%unit,'(/,10x,51("."))')
      write(env%unit,'(10x,":",22x,a,22x,":")') "SETUP"
      write(env%unit,'(10x,":",49("."),":")')
      write(env%unit,intfmt) "# basis functions  ",basis%nbf
      write(env%unit,intfmt) "# atomic orbitals  ",basis%nao
      write(env%unit,intfmt) "# shells           ",basis%nshell
      write(env%unit,intfmt) "# electrons        ",wfn%nel
      if (allocated(xtbData%halogen)) &
      write(env%unit,intfmt) "# halogen bonds    ",nxb
      write(env%unit,intfmt) "max. iterations    ",maxiter
      if (gfn_method == 2) &
      write(env%unit,chrfmt) "Hamiltonian        ","GFN2-xTB"
      if (gfn_method == 1) &
      write(env%unit,chrfmt) "Hamiltonian        ","GFN1-xTB"
      write(env%unit,chrfmt) "restarted?         ",bool2string(restart)
      write(env%unit,chrfmt) "GBSA solvation     ",bool2string(lgbsa)
      write(env%unit,chrfmt) "PC potential       ",bool2string(lpcem)
      if (lpcem) then
         write(env%unit,intfmt) "-> # point charges ",pcem%n
         write(env%unit,dblfmt) "-> sum of PC       ",sum(pcem%q),"e   "
      endif
      write(env%unit,dblfmt) "electronic temp.   ",et,      "K   "
      write(env%unit,dblfmt) "accuracy           ",acc,     "    "
      write(env%unit,scifmt) "-> integral cutoff ",intcut,  "    "
      write(env%unit,scifmt) "-> integral neglect",neglect, "    "
      write(env%unit,scifmt) "-> SCF convergence ",scfconv, "Eh  "
      write(env%unit,scifmt) "-> wf. convergence ",qconv,   "e   "
      write(env%unit,dblfmt) "Broyden damping    ",broydamp,"    "
      write(env%unit,'(10x,51("."))')
   endif

   qq    =wfn%q
   damp  =damp0

   if (profile) call timer%measure(1)
   if (profile) call timer%measure(2,"Dispersion")

   ! ========================================================================
   ! Coordination number
   if (gfn_method == 1) then
      ! D3 part first because we need CN
      call getCoordinationNumber(mol, trans, 40.0_wp, cnType%exp, cn, dcndr, dcndL)
   else
      ! CN/dCN replaced by special smoother and faster decaying function
      call getCoordinationNumber(mol, trans, 40.0_wp, cnType%gfn, cn, dcndr, dcndL)
   endif

   ! ------------------------------------------------------------------------
   ! dispersion (DFT-D type correction)
   if (gfn_method == 1) then
      call d3_gradient &
         & (mol, trans, xtbData%dispersion%dpar, 4.0_wp, 60.0_wp, &
         &  cn, dcndr, dcndL, ed, gradient, sigma)
   else
      allocate(scD4)
      call init(scD4, xtbData%dispersion, mol)
   endif

   if (profile) call timer%measure(2)
   if (profile) call timer%measure(6,"classical contributions")

   ! ========================================================================
   ! Halogen bond correction
   exb=0.0_wp
   if (allocated(xtbData%halogen)) then
      call xbpot(xtbData%halogen,mol%n,mol%at,mol%xyz,xblist,nxb,&
         & ljexp,exb,gradient)
   end if

   ! ------------------------------------------------------------------------
   ! Repulsion energy
   ep = 0.0_wp
   call repulsionEnGrad(mol, xtbData%repulsion, trans, 40.0_wp, &
      & ep, gradient, sigma)

   if (profile) call timer%measure(6)
   if (profile) call timer%measure(3,"integral evaluation")

   ! ========================================================================
   ! Overlap integrals
   allocate(dpint(3,basis%nao,basis%nao), &
      &     qpint(6,basis%nao,basis%nao), &
      &     source = 0.0_wp)

   call getSelfEnergy(xtbData%hamiltonian, xtbData%nShell, mol%at, cn=cn, &
      & selfEnergy=selfEnergy, dSEdcn=dSEdcn)
   ! compute integrals and prescreen to set up list arrays
   call build_SDQH0(xtbData%nShell, xtbData%hamiltonian, mol%n, mol%at, basis%nbf, &
      & basis%nao, mol%xyz, selfEnergy, intcut, basis%caoshell, basis%saoshell, &
      & basis%nprim, basis%primcount, basis%alp, basis%cont, S, dpint, qpint, H0)
   call count_dpint(ndp, dpint, neglect)
   call count_qpint(nqp, qpint, neglect)

   ! prepare aes stuff
   if (allocated(xtbData%multipole)) then
      allocate(aes)
      call init(aes, xtbData%multipole)

      ! allocate arrays for lists and fill (to exploit sparsity)
      allocate(mdlst(2,ndp),mqlst(2,nqp))
      call setdqlist(basis%nao,ndp,nqp,neglect,dpint,qpint,mdlst,mqlst)
      ! set up 1/R^n * damping function terms
      ii=mol%n*(mol%n+1)/2
      allocate(aes%gab3(ii),aes%gab5(ii),radcn(mol%n))
      call get_radcn(xtbData%multipole,mol%n,mol%at,cn,aes%cnShift, &
         & aes%cnExp,aes%cnRMax,radcn)
      call mmomgabzero(mol%n,mol%at,mol%xyz,aes%dipDamp, &
         & aes%quadDamp,radcn,aes%gab3,aes%gab5) ! zero damping
   end if

   ! ------------------------------------------------------------------------
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
   enddo

   if (profile) call timer%measure(3)
   if (profile) call timer%measure(5,"iterations")

   ! ========================================================================
   ! SCC iterations

   if(pr)then
      write(env%unit,'(a)')
      write(env%unit,*) 'iter      E             dE          RMSdq', &
      &'      gap      omega  full diag'
   endif

   call scc(env,xtbData,mol%n,wfn%nel,wfn%nopen,basis%nao,ndp,nqp,nmat,basis%nshell, &
      &     mol%at,matlist,mdlst,mqlst,basis%aoat2,basis%ao2sh,basis%ash, &
      &     wfn%q,wfn%dipm,wfn%qp,qq,qlmom,wfn%qsh,zsh, &
      &     mol%xyz,aes, &
      &     gbsa,fgb,fhb,cm5,cm5a,gborn, &
      &     scD4, &
      &     broy,broydamp,damp0, &
      &     lpcem,ves,vpc, &
      &     et,wfn%focc,wfn%focca,wfn%foccb,wfn%efa,wfn%efb, &
      &     eel,ees,eaes,epol,ed,epcem,egap, &
      &     wfn%emo,wfn%ihomo,wfn%ihomoa,wfn%ihomob, &
      &     H0,wfn%C,S,dpint,qpint,X,wfn%P,ies, &
      &     maxiter,startpdiag,scfconv,qconv, &
      &     minpr,pr, &
      &     fail,jter)

   ! check if something terrible happend in the SCC
   call env%check(exitRun)
   if (exitRun) then
      call env%error("Self consistent charge iterator terminated", source)
      return
   end if

   ! check for convergence, only do this if printlevel is maximal (WHY?)
   res % converged = .not. fail
   if (fail) then
      call env%warning("Self consistent charge iterator did not converge", source)
   end if
   if (pr) then
      if (fail) then
         call touch_file('.sccnotconverged')
         write(env%unit,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***")') &
            "convergence criteria cannot be satisfied within",jter,"iterations"
      else
         write(env%unit,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***")') &
            "convergence criteria satisfied after",jter,"iterations"
      endif
   endif

   if (profile) call timer%measure(5)
   if (.not.pr.and.profile.and.minpr) &
      call timer%write_timing(env%unit,5,'SCC iter.')
   if (profile) call timer%measure(6,"molecular gradient")

   ! ========================================================================
   ! GRADIENT (finally 100% analytical)
   call scf_grad(mol, nmat2, matlist2, H0, ves, S, xtbData, trans, selfEnergy, &
      & dSEdcn, sqrab, wfn, basis, intcut, aes, radcn, gradient, prlevel)

   ! ------------------------------------------------------------------------
   ! dispersion (DFT-D type correction)
   if (allocated(scD4)) then
      call getCoordinationNumber(mol, trans, 40.0_wp, cnType%cov, &
         & cn, dcndr, dcndL)
      call d4_gradient(mol, trans, xtbData%dispersion%dpar, scD4%g_a, scD4%g_c, &
         &  scD4%wf, 60.0_wp, 40.0_wp, cn, dcndr, dcndL, wfn%q, &
         &  energy=dum, gradient=gradient, sigma=sigma, e3=embd)
   endif

   ! ------------------------------------------------------------------------
   ! Solvation contributions from GBSA
   if (lgbsa) then
      if (gfn_method > 1) then
         call compute_gb_egrad(gbsa, wfn%q, gborn, ghb, gradient, minpr)
      else
         cm5(:)=wfn%q+cm5a
         call compute_gb_egrad(gbsa, cm5, gborn, ghb, gradient, minpr)
         call cm5_grad_gfn1(gradient, mol%n, cm5, fgb, fhb, dcm5a, lhb)
      endif
!     solvation energy
      gbsa%gborn = gborn
      gbsa%ghb = ghb
      gsolv = gborn+gbsa%gsasa+gbsa%ghb+gshift
      eel = eel+gbsa%gsasa+gshift
   endif

   ! ------------------------------------------------------------------------
   ! Derivative of electrostatic energy
   allocate(djdr(3, mol%n, basis%nshell))
   allocate(djdtr(3, basis%nshell))
   allocate(djdL(3, 3, basis%nshell))
   call coulomb%getCoulombDerivs(mol, wfn%qsh, djdr, djdtr, djdL)
   call contract(djdr, wfn%qsh, gradient, beta=1.0_wp)
   !call contract(djdL, wfn%qsh, sigma, beta=1.0_wp)

   ! ------------------------------------------------------------------------
   ! ES point charge embedding
   if (lpcem) then
      if (gfn_method == 1) then
         call pcem_grad_gfn1(xtbData%coulomb,gradient,pcem%grd,mol%n,pcem,mol%at, &
            & xtbData%nshell,mol%xyz,xtbData%coulomb%gExp,wfn%qsh)
      else
         call pcem_grad_gfn2(xtbData%coulomb,gradient,pcem%grd,mol%n,pcem,mol%at, &
            & xtbData%nshell,mol%xyz,wfn%qsh)
      endif
   endif

   if (profile) call timer%measure(6)
   if (.not.pr.and.profile.and.minpr) &
      call timer%write_timing(env%unit,6,'gradient')
   if (profile) call timer%measure(7,"printout")

   ! ========================================================================
   ! PROPERTIES & PRINTOUT

   printing: if (pr) then
      ! ---------------------------------------------------------------------
      ! print orbital energies and occupation numbers
      if (pr_eig) then
         call print_orbital_eigenvalues(env%unit, wfn, 5)
      endif

      ! ---------------------------------------------------------------------
      ! HOMO-LUMO excitation properties if UHF=2
      if (wfn%nopen.eq.2) then
         call hlex(mol%n, mol%at, basis%nbf, basis%nao, wfn%ihomoa, mol%xyz, &
            & wfn%focc, S, wfn%C, wfn%emo, basis)
      endif

      ! ---------------------------------------------------------------------
      ! LMO /xTB-IFF
      if (pr_lmo) then
         tmp=wfn%emo*evtoau
         call local(mol%n, mol%at, basis%nbf, basis%nao, wfn%ihomoa, mol%xyz, &
            & mol%z, wfn%focc, S, wfn%P, wfn%C, tmp, wfn%q, eel, lgbsa, basis)
      endif

   endif printing

   ! ------------------------------------------------------------------------
   ! get Wiberg bond orders
   call get_wiberg(mol%n, basis%nao, mol%at, mol%xyz, wfn%P, S, wfn%wbo, &
      & basis%fila2)

   ! ------------------------------------------------------------------------
   ! dipole calculation (always done because its free)
   if (.not.allocated(xtbData%multipole)) then
      call mmompop(mol%n,basis%nao,basis%aoat2,mol%xyz,wfn%p,s,dpint,qpint, &
         &         wfn%dipm,wfn%qp)
   endif
   call calc_dipole(mol%n,mol%at,mol%xyz,mol%z,basis%nao,wfn%P,dpint,dip,dipol)

   if (profile) call timer%measure(7)


   ! ========================================================================
   ! SAVE FOR FINAL PRINTOUT
   ! total energy
   energy = eel + ep + exb + embd
   eat = eatoms*evtoau - energy
   if (.not.allocated(scD4)) then
      energy = energy + ed
   endif
   res%e_elec  = eel
   res%e_atom  = eat
   res%e_rep   = ep
   res%e_es    = ees
   res%e_aes   = eaes
   res%e_axc   = epol
   res%e_disp  = ed+embd
   res%e_total = energy
   res%hl_gap  = egap
   res%dipole  = dipol
   if (allocated(xtbData%halogen)) res%e_xb = exb
   if (lgbsa) then
      res%g_solv  = gsolv
      res%g_born  = gborn
      res%g_sasa  = gbsa%gsasa
      res%g_hb    = gbsa%ghb
      res%g_shift = gshift
   endif
   res%gnorm = norm2(gradient)

   if (profile.and.pr) call timer%write(env%unit,'SCC')

! ========================================================================
   if (profile) call timer%deallocate

   call deallocate_gbsa(gbsa)
end subroutine scf

subroutine scf_grad(mol, nmat2, matlist2, H0, shellShift, S, xtbData, trans, &
      & selfEnergy, dSEdcn, sqrab, wfn, basis, intcut, aes, radcn, gradient, &
      & printlvl)

! ========================================================================
!  type definitions
   use xtb_type_wavefunction
   use xtb_type_basisset
   use xtb_type_param
   use xtb_type_pcem

! ========================================================================
!  global storage
   use xtb_setparam

! ========================================================================
!  interfaces
   use xtb_scc_core
   use xtb_grad_core

   use xtb_aespot,    only : dradcn,aniso_grad,setdvsdq,dsint
   use xtb_solv_gbobc,     only : lgbsa,lhb,TSolvent,gshift,compute_gb_egrad
   use xtb_disp_dftd4, only: d4_gradient
   use xtb_disp_ncoord,    only : dncoord_gfn,dncoord_d3
   use xtb_embedding, only : pcem_grad_gfn1,pcem_grad_gfn2

   implicit none

   type(TMolecule), intent(in) :: mol
   type(TWavefunction),intent(in) :: wfn
   type(TBasisset),    intent(in) :: basis
   type(TxTBData), intent(in) :: xtbData
   integer, intent(in)    :: nmat2
   integer,intent(in) :: matlist2(2,nmat2)
   real(wp),intent(in)    :: selfEnergy(:, :)
   real(wp),intent(in)    :: dSEdcn(:, :)
   real(wp),intent(in) :: trans(:, :)
   real(wp),intent(in)    :: sqrab(:)
   real(wp),intent(inout) :: H0(basis%nao*(basis%nao+1)/2)
   real(wp),intent(in) :: shellShift(:)
   real(wp),intent(inout) :: gradient(:, :)
   real(wp),intent(in)    :: S(basis%nao,basis%nao)
   type(TxTBMultipole), intent(in), optional :: aes
   real(wp),intent(in)    :: intcut
   real(wp),intent(in)    :: radcn(:)
   integer, intent(in)    :: printlvl

   integer :: m,i,j,kk,ish,jsh,iat,izp
   real(wp) :: h1
   real(wp),allocatable :: qq(:)
   real(wp),allocatable :: cn(:), dcndr(:,:,:), dcndL(:,:,:)
   real(wp),allocatable :: Pew(:,:)
   real(wp),allocatable :: ves(:, :)
   real(wp),allocatable :: dhdcn(:)
   real(wp),allocatable :: vs(:),vd(:,:),vq(:,:)
   logical :: pr,minpr

   minpr = printlvl.gt.0
   pr = printlvl.gt.1

   allocate(cn(mol%n), source=0.0_wp)
   allocate(dcndr(3,mol%n,mol%n), source=0.0_wp)
   allocate(dcndL(3,3,mol%n), source=0.0_wp)
   allocate(Pew(basis%nao,basis%nao), source=0.0_wp)
   allocate(ves(maxval(xtbData%nShell), mol%n))
   allocate(dhdcn(mol%n))

   ! get energy weighted density matrix and convert from eV to Eh
   call prep_grad_conv(basis%nao, wfn%C, wfn%focc, wfn%emo, Pew)

   ! CN dependent part
   if (gfn_method == 1) then
      call getCoordinationNumber(mol, trans, 40.0_wp, cnType%exp, cn, dcndr, dcndL)
   else
      call getCoordinationNumber(mol, trans, 40.0_wp, cnType%gfn, cn, dcndr, dcndL)
   end if

   ves(:, :) = 0.0_wp
   i = 0
   do iat = 1, mol%n
      izp = mol%at(iat)
      do ish = 1, xtbData%nShell(izp)
         ves(ish, iat) = shellShift(i+ish)
      end do
      i = i+xtbData%nShell(izp)
   end do

   allocate( vs(mol%n),vd(3,mol%n),vq(6,mol%n), source = 0.0_wp )
   if (allocated(xtbData%multipole)) then
      ! VS, VD, VQ-dependent potentials are changed w.r.t. SCF,
      ! since moment integrals are now computed with origin at
      ! respective atoms
      call setdvsdq(xtbData%multipole, mol%n, mol%at, mol%xyz, wfn%q, wfn%dipm, &
         & wfn%qp, aes%gab3, aes%gab5, vs, vd, vq)
   end if

   dhdcn(:) = 0.0_wp
   call build_dSDQH0(xtbData%nShell, xtbData%hamiltonian, selfEnergy, dSEdcn, &
      & intcut, mol%n, basis%nao, basis%nbf, mol%at, mol%xyz, basis%caoshell, &
      & basis%saoshell, basis%nprim, basis%primcount, basis%alp, basis%cont, &
      & wfn%p, Pew, ves, vs, vd, vq, dhdcn, gradient)
   ! setup CN gradient
   call contract(dcndr, dhdcn, gradient, beta=1.0_wp)

   ! multipole gradient stuff
   if (allocated(xtbData%multipole)) then
      ! VS, VD, VQ-dependent potentials are changed w.r.t. SCF,
      ! since moment integrals are now computed with origin at
      ! respective atoms
      call setdvsdq(xtbData%multipole, mol%n, mol%at, mol%xyz, wfn%q, wfn%dipm, &
         & wfn%qp, aes%gab3, aes%gab5, vs, vd, vq)

      ! WARNING: dcndr is overwritten on output and now dR0A/dXC
      call dradcn(xtbData%multipole, mol%n, mol%at, cn, aes%cnShift, &
         & aes%cnExp, aes%cnRMax, dcndr)
      call aniso_grad(mol%n, mol%at, mol%xyz, wfn%q, wfn%dipm, wfn%qp, &
         & aes%dipDamp, aes%quadDamp, radcn, dcndr, aes%gab3, aes%gab5, gradient)
   end if

end subroutine scf_grad


pure elemental function early3d(i) result(bool)
   implicit none
   integer,intent(in) :: i
   logical :: bool
   bool = .false.
   if ((i.ge.21).and.(i.le.24)) bool = .true.
end function early3d


!> X..Y bond? (X=halogen but not F, Y=N,O,P,S) !CB: Cl is formally included, but
!  has a parameter of zero in GFN1-xTB
pure elemental function xbond(ati,atj) result(bool)
   integer,intent(in) :: ati,atj
   logical :: bool
   logical :: lx1,lx2,ly1,ly2

   bool=.false.
   lx1=.false.
   lx2=.false.
   ly1=.false.
   ly2=.false.

   if(ati.eq.17.or.ati.eq.35.or.ati.eq.53.or.ati.eq.85) lx1=.true.
   if(atj.eq.17.or.atj.eq.35.or.atj.eq.53.or.atj.eq.85) lx2=.true.
   if(ati.eq. 7.or.ati.eq. 8.or.ati.eq.15.or.ati.eq.16) ly1=.true.
   if(atj.eq. 7.or.atj.eq. 8.or.atj.eq.15.or.atj.eq.16) ly2=.true.

   if(lx1.and.ly2) bool=.true.
   if(lx2.and.ly1) bool=.true.

end function xbond


end module xtb_scf
