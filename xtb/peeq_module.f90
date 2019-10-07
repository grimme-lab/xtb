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

module peeq_module
! ------------------------------------------------------------------------
!  PEEQ method developed and implemented by E. Caldeweyher
!   01/19
!  with help of (alphabetically sorted):
!  S. Ehlert, S. Grimme, P. Pracht
! ------------------------------------------------------------------------
   use iso_fortran_env, wp => real64
   use mctc_econv
   use mctc_la
   implicit none

   !> print the different gradient contributions (for debugging)
   logical,parameter,private :: gpr =.false.

   !> profiling
   logical,parameter,private :: profile = .true.

!--- Parameter for the hydrogen bond correction
!-- cut-offs tested on 1000 atom system riva+H2O yielding sub mEh accuracy for EHB
   !> save R^2 ABH triangle cut-off
   real(wp),parameter,private :: hbthr  = 400.0_wp
   !> hbthr*2
   real(wp),parameter,private :: hbthr3 = 800.0_wp
   !> neighbor list R^2 cut-off for LP grad part
   real(wp),parameter,private :: hbthr2 = 150.0_wp

   !> about 8 Bohr smooth cut-off in rAB
   real(wp),parameter,private :: longcut=56.0_wp
   !> smoothness of damping long range
   real(wp),parameter,private :: hbalp=6.0_wp
   !> smoothness of damping short range
   real(wp),parameter,private :: hbalp1=4.0_wp
   !> smoothness of damping out-of-lin
   real(wp),parameter,private :: hbalp2=1.5_wp
   !> shift of A,B charge scaling fct
   real(wp),parameter,private :: sf_ab=5.0_wp
   !> steepness of A,B charge scaling fct
   real(wp),parameter,private :: st_ab=15.0_wp
   !> shift of H charge scaling fct
   real(wp),parameter,private :: sf_h=15.0_wp
   !> steepness of H charge scaling fct
   real(wp),parameter,private :: st_h=15.0_wp

   !> smoothing factor for periodic CN cutoff
   real(wp),parameter,private :: beta = 7.5_wp

   !> mapping from l-quantum numbers to shell types
   integer,parameter,private :: shell(*) = [1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4]

contains

subroutine peeq &
      (iunit,mol,wfn,basis,param,egap,et,prlevel,grd,ccm,acc,etot,g,sigma,res)

! ------------------------------------------------------------------------
!  Class definitions
! ------------------------------------------------------------------------
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data
   use tbdef_timer

! ------------------------------------------------------------------------
!  Global storage
! ------------------------------------------------------------------------
   use aoparam

! ------------------------------------------------------------------------
!  Get interfaces
! ------------------------------------------------------------------------
   use readin
   use aespot
   use scc_core
   use grad_core
   use eeq_model
   use ncoord
   use lidep
   use gbobc
   use pbc

   implicit none

! ------------------------------------------------------------------------
!  INPUT
! ------------------------------------------------------------------------
   integer, intent(in)            :: iunit
   type(tb_molecule),  intent(in) :: mol     !< molecular structure infomation
   type(tb_basisset),  intent(in) :: basis   !< basis set
   type(scc_parameter),intent(in) :: param   !< method parameters
   real(wp),intent(in)            :: et      !< electronic temperature
   integer, intent(in)            :: prlevel !< amount of printout
   logical, intent(in)            :: grd     !< toggles gradient calculation
   real(wp),intent(in)            :: acc     !< numerical accuracy
   logical, intent(in)            :: ccm     !< use cyclic cluster model

! ------------------------------------------------------------------------
!  INPUT/OUTPUT
! ------------------------------------------------------------------------
   real(wp),intent(inout)                    :: etot !< total energy
   real(wp),intent(inout)                    :: egap !< HOMO-LUMO gap
   real(wp),intent(inout),dimension(3,mol%n) :: g    !< molecular gradient
   type(tb_wavefunction),intent(inout)       :: wfn  !< TB-wavefunction
! ------------------------------------------------------------------------
!  OUTPUT
! ------------------------------------------------------------------------
   type(scc_results),intent(out) :: res !< bundles all calculation informations

! ------------------------------------------------------------------------
   type(chrg_parameter) :: chrgeq
   real(wp),allocatable,dimension(:)     :: cn
   real(wp),allocatable,dimension(:)     :: sqrab
   real(wp),allocatable,dimension(:,:,:) :: dcndr
   real(wp),allocatable,dimension(:,:,:) :: dcndL
   real(wp),allocatable,dimension(:,:)   :: X
   real(wp),allocatable,dimension(:,:)   :: S
   real(wp),allocatable,dimension(:,:)   :: H
   real(wp),allocatable,dimension(:)     :: H0
   real(wp),allocatable,dimension(:)     :: H1
   real(wp),allocatable,dimension(:)     :: zsh
   real(wp),allocatable,dimension(:)     :: kcnao
   real(wp),allocatable,dimension(:)     :: kqao
   integer :: rep_cn(3)

! ------------------------------------------------------------------------
   ! local variable (without comment)
   integer :: ii,jj,i,j,k,m,iat,jat,atj,kk
   integer :: il,jl

   integer,parameter :: lladr(4)    = [1,3,6,10]
   integer,parameter :: lladr2(0:3) = [1,3,5,7 ]

   real(wp) :: dipol(3)
   real(wp) :: eatoms,eat
   real(wp) :: dum

   character(len=2), external :: asym
   logical :: ex
   logical :: minpr
   logical :: pr
   logical :: fail

   ! timer
   type(tb_timer) :: timer
   ! Fermi smearing variables
   real(wp) :: efa,efb,nfoda,nfodb,ga,gb
   ! numerical thresholds
   real(wp) :: intcut,neglect,neglect2,scfconv
   integer  :: intrep(3)
   ! energies
   real(wp) :: eel,ed,ees,exb,ep

! ---------------------------------------
!  PEEQ information
! ---------------------------------------
   real(wp) , allocatable, dimension(:)     :: rab0
   real(wp) , allocatable, dimension(:,:,:) :: grab0
   real(wp) , allocatable, dimension(:)     :: lcn
   real(wp) , allocatable, dimension(:,:,:) :: dlcn
   integer(int16), allocatable, dimension(:,:)   :: hblist
   integer(int16), allocatable, dimension(:,:)   :: nlist
   real(wp),               dimension(2)     :: hbpi
   real(wp),               dimension(2)     :: hbpj
   real(wp),               dimension(4,4)   :: ken
   real(wp) :: esrb
   real(wp) :: ehb
   real(wp) :: den2
   real(wp) :: den4
   real(wp) :: ken2
   real(wp) :: ken4
   integer  :: nhb
   integer  :: nh
   integer  :: ich
   logical :: debug

! ---------------------------------------
!  PEEQ q dependent H0
! ---------------------------------------
   real(wp),allocatable, dimension(:)     :: qeeq
   real(wp),allocatable, dimension(:,:,:) :: dqdr
   real(wp),allocatable, dimension(:,:,:) :: dqdL
   real(wp),allocatable, dimension(:)     :: dHdq
   real(wp),allocatable, dimension(:)     :: dHdcn
   real(wp),allocatable, dimension(:,:)   :: pew
   real(wp),allocatable, dimension(:)     :: tmp

! ---------------------------------------
!  PEEQ WSC information
! ---------------------------------------
   real(wp),dimension(3,3),intent(inout)     :: sigma
   real(wp),dimension(3,3,5)                 :: sigma_tmp

! For eigenvalues of S via dsyev routine
   real(wp),allocatable,dimension(:,:)       :: Stmp
   real(wp),allocatable,dimension(:)         :: ev_s
   integer                                   :: info
   integer                                   :: lwork
   real(wp),allocatable,dimension(:)         :: aux
   real(wp) :: test(1)

   integer :: xdim ! new dimension of the cut Overlap (lin.depend.)

   character(len=*),parameter :: scifmt = &
      '(10x,":",2x,a,e22.7,1x,a,1x,":")'
   character(len=*),parameter :: dblfmt = &
      '(10x,":",2x,a,f18.7,5x,a,1x,":")'
   character(len=*),parameter :: intfmt = &
      '(10x,":",2x,a,i18,      10x,":")'
   character(len=*),parameter :: chrfmt = &
      '(10x,":",2x,a,a18,      10x,":")'

! ---------------------------------------
!  EEQ/GBSA information
! ---------------------------------------
   type(tb_solvent) :: gbsa
   real(wp) :: gsolv

   associate( nao => basis%nao, &
         &    nbf => basis%nbf, &
         &    nshell => basis%nshell, &
         &    naop => basis%nao*(basis%nao+1)/2, &
         &    nbfp => basis%nbf*(basis%nbf+1)/2)

   if (profile) then
      if (lgbsa) then
         call timer%new(9,.false.)
      else
         call timer%new(8,.false.)
      endif
   endif
   if (profile) call timer%measure(1,"EHT setup")

   etot  = 0.0_wp
   dipol = 0.0_wp
   exb   = 0.0_wp
   ga    = 0.0_wp
   gb    = 0.0_wp
   eel   = 0.0_wp
   ep    = 0.0_wp
   ehb   = 0.0_wp
   esrb  = 0.0_wp
   ed    = 0.0_wp
   ees   = 0.0_wp
   sigma_tmp = 0.0_wp
   g = 0.0_wp

   !debug = prlevel.gt.2
   debug =.false.
   pr    = prlevel.gt.1
   minpr = prlevel.gt.0

! ---------------------------------------
!  Numerical stuff and cutoffs
! ---------------------------------------
   intcut=25.0_wp-10.0_wp*log10(acc)
   intcut=max(20.0_wp,intcut)
   neglect =10.0e-9_wp*acc
   neglect2=neglect*10.0_wp
   scfconv=1.e-6_wp*acc

   if (mol%npbc > 0) call get_realspace_cutoff(mol%lattice,800.0_wp,intrep)

! ---------------------------------------
!  IMPORTANT FACT: H is given in eV
! ---------------------------------------

! ---------------------------------------
!  Get memory
! ---------------------------------------
   allocate(sqrab(mol%n*(mol%n+1)/2));   sqrab = 0.0_wp
   allocate(cn(mol%n));                     cn = 0.0_wp
   allocate(dcndr(3,mol%n,mol%n));       dcndr = 0.0_wp
   allocate(dcndL(3,3,mol%n));           dcndL = 0.0_wp
   allocate(H(nao,nao));                     H = 0.0_wp
   allocate(X(nao,nao));                     X = 0.0_wp
   allocate(H0(naop));                      H0 = 0.0_wp
   allocate(S(nao,nao));                     S = 0.0_wp
   allocate(H1(naop));                      H1 = 0.0_wp
   allocate(kcnao(nao));                 kcnao = 0.0_wp
   allocate(kqao(nao));                   kqao = 0.0_wp
   allocate(zsh(nshell));                  zsh = 0.0_wp
   allocate(qeeq(mol%n));           qeeq = 0.0_wp
   allocate(dqdr(3,mol%n,mol%n+1)); dqdr = 0.0_wp
   allocate(dqdL(3,3,mol%n+1));     dqdL = 0.0_wp

! ---------------------------------------
!  Fill levels
! ---------------------------------------
   call setzshell(mol%n,mol%at,nshell,mol%z,zsh,eatoms,0)
   if(wfn%nel.ne.0) then
      call occu(nao,wfn%nel,wfn%nopen,wfn%ihomoa,wfn%ihomob,wfn%focca,wfn%foccb)
      wfn%focc = wfn%focca + wfn%foccb
      wfn%ihomo=wfn%ihomoa
   else
      wfn%focc=0.0_wp
      wfn%ihomo=0
      wfn%ihomoa=0
      wfn%nopen=0
   endif

! ---------------------------------------
!  Calculate distances
! ---------------------------------------
   k = 0
   do i=1,mol%n
      do j=1,i
         k=k+1
         sqrab(k)=(mol%xyz(1,i)-mol%xyz(1,j))**2 &
         &       +(mol%xyz(2,i)-mol%xyz(2,j))**2 &
         &       +(mol%xyz(3,i)-mol%xyz(3,j))**2
      enddo
   enddo

   if (prlevel > 1) then
      write(iunit,'(/,10x,51("."))')
      write(iunit,'(10x,":",22x,a,22x,":")') "SETUP"
      write(iunit,'(10x,":",49("."),":")')
      write(iunit,intfmt) "# basis functions  ",basis%nbf
      write(iunit,intfmt) "# atomic orbitals  ",basis%nao
      write(iunit,intfmt) "# shells           ",basis%nshell
      write(iunit,intfmt) "# electrons        ",wfn%nel
      if (mol%npbc > 0) &
      write(iunit,chrfmt) "PBC by CCM         ",bool2string(ccm)
      write(iunit,dblfmt) "electronic temp.   ",et,      "K   "
      write(iunit,dblfmt) "accuracy           ",acc,     "    "
      write(iunit,scifmt) "-> integral cutoff ",intcut,  "    "
      write(iunit,scifmt) "-> integral neglect",neglect, "    "
      if (mol%npbc > 0 .and. .not.ccm) then
      write(iunit,intfmt) "# images for ints. ",product(2*intrep+1)
      endif
      write(iunit,'(10x,51("."))')
   endif

! ---------------------------------------
!  set up parameters (TODO - move to gfn_paramset.f90)
! ---------------------------------------
   ken(1,1)=0.01_wp*param%kspd(4)
   ken(2,2)=0.01_wp*param%kspd(6)
   ken(3,3)=0.01_wp*param%kenscal
   do i=1,3
     do j=1,3
     ken(i,j)=(ken(i,i)+ken(j,j))*0.5_wp
     enddo
   enddo

   if (profile) call timer%measure(1)
   if (profile) call timer%measure(2,"Coordination number")
! ---------------------------------------
!  Get CN(1:n) + dcndr(3,1:n,1:n) under pbc
! ---------------------------------------
   if (mol%npbc > 0) then
      call get_erf_cn(mol,cn,dcndr,dcndL,thr=900.0_wp)
      call dncoord_logcn(mol%n,cn,dcndr,dcndL,cn_max=8.0_wp)
   else
      call get_erf_cn(mol,cn,dcndr,thr=900.0_wp)
      call dncoord_logcn(mol%n,cn,dcndr,cn_max=8.0_wp)
   endif

   if (profile) call timer%measure(2)
! ---------------------------------------
!  Get EEQ charges q(1:n) + dqdr(3,1:n,1:n) under pbc
! ---------------------------------------

   if (lgbsa) then
      if (profile) call timer%measure(9,"GBSA setup")
      call new_gbsa(gbsa,mol%n,mol%at)
      call update_nnlist_gbsa(gbsa,mol%xyz,.false.)
      ! compute Born radii
      call compute_brad_sasa(gbsa,mol%xyz)
      ! add SASA term to energy and gradient
      ees = gbsa%gsasa
      gsolv = gbsa%gsasa
      g = g + gbsa%dsdr
      if (profile) call timer%measure(9)
   endif


   if (profile) call timer%measure(3,"EEQ model density")
   ! names DO NOT corresponds to content of variables, obviously...
   call gfn0_charge_model(chrgeq,mol%n,mol%at,dpolc,gam,cxb,alp0)
   ! initialize electrostatic energy
   if (lgbsa) then
      call eeq_chrgeq(mol,chrgeq,gbsa,cn,dcndr,qeeq,dqdr, &
         &            ees,gsolv,g,.false.,.true.,.true.)
   else
      call eeq_chrgeq(mol,chrgeq,cn,dcndr,dcndL,qeeq,dqdr,dqdL, &
         &            ees,g,sigma_tmp(:,:,1),&
         &            .false.,.true.,.true.)
   endif

   wfn%q = qeeq

   if (profile) call timer%measure(3)
   if (profile) call timer%measure(4,"D4 Dispersion")
! ----------------------------------------
!  D4 dispersion energy + gradient (2B) under pbc
! ----------------------------------------
   call ddisp_peeq(mol,param,cn,dcndr,dcndL,grd,ed,g,sigma_tmp(:,:,2))

   if (profile) call timer%measure(4)
   if (profile) call timer%measure(5,"Integral evaluation")
! ---------------------------------------
!  Build AO overlap S and H0 integrals under pbc
! ---------------------------------------
   if (mol%npbc > 0) then
      if (ccm) then
         call ccm_build_SH0(mol%n,mol%at,basis,nbf,nao,mol%xyz,mol%lattice, &
            &               wfn%q,cn,intcut, &
            &               param%kmagic,ken,param%alphaj,param%kcnsh, &
            &               param%xbdamp,s,h0,mol%wsc)
      else
         call pbc_build_SH0(mol%n,mol%at,basis,nbf,nao,mol%xyz,mol%lattice,intrep,&
            &               wfn%q,cn,intcut, &
            &               param%kmagic,ken,param%alphaj,param%kcnsh, &
            &               param%xbdamp,s,h0)
      endif
   else
      call mol_build_SH0(mol%n,mol%at,basis,nbf,nao,mol%xyz,wfn%q,cn,intcut, &
         &               param%kmagic,ken,param%alphaj,param%kcnsh, &
         &               param%xbdamp,s,h0)
   endif

   if (profile) call timer%measure(5)
   if (profile) call timer%measure(6,"Cholesky factorization")
  
! ---------------------------------------
!  Check for near linear dependencies via Cholesky decomposition
! ---------------------------------------
   call cholesky(iunit,pr,nao,S,orthog) ! S is not modified
   if(orthog)then
      !if (profile) call timer%measure(7,"Canonical orthogonalization")
      call renorm(nao,S) ! S is renormalized
      call canorthog(iunit,nao,S,X,xdim,pr,fail)
      !if (profile) call timer%measure(7)
   endif

   if (profile) call timer%measure(6)
   if (profile) call timer%measure(7,"Zeroth order Hamiltonian")
! ---------------------------------------
!  Setup H0 under pbc
! ---------------------------------------
!  H0 already initialized
   do i = 1, nao
      do j = 1, i
         k = j+i*(i-1)/2
         H(j,i) = H0(k)
         H(i,j) = H(j,i)
      enddo
   enddo

   if(.not.orthog)then
      call solve(.true.,nao,wfn%ihomo,scfconv,H,S,X,wfn%P,wfn%emo,fail)
   else
      call orthgsolve2(.true.,nao,xdim,wfn%ihomo,scfconv,H,S,X,wfn%P,wfn%emo,fail)
   endif
   ! for geometry optimization
   res % converged = .not. fail

   ! save eigenvectors
   wfn%C = H

   if (pr.and.fail) call raise('E',"Diagonalization failed!",1)

! ---------------------------------------
!  Fermi smearing
! ---------------------------------------
   if(et.gt.0.1_wp)then
      if (wfn%ihomoa+1.le.nao) then
         call fermismear(.false.,nao,wfn%ihomoa,et,wfn%emo,wfn%focca,nfoda,efa,ga)
      endif
      if (wfn%ihomob+1.le.nao) then
         call fermismear(.false.,nao,wfn%ihomob,et,wfn%emo,wfn%foccb,nfodb,efb,gb)
      endif
      wfn%focc = wfn%focca + wfn%foccb
   endif
   ! create density matrix = save in wfn%P
   call dmat(nao,wfn%focc,wfn%C,wfn%P)

   call mpopsh(mol%n,nao,nshell,basis%ao2sh,S,wfn%P,wfn%qsh)

! ---------------------------------------
!  emo: Eigenvalues from H0
! ---------------------------------------
   eel = sum(wfn%focc*wfn%emo)*evtoau + ga + gb

   if (profile) call timer%measure(7)
   if (.not.pr.and.profile.and.minpr) &
      call timer%write_timing(iunit,7,"Diagonalization")
   if (profile) call timer%measure(8,"Gradient calculation")
! ======================================================================
!  GRADIENT (100% analytical)
! ======================================================================
   ! repulsion energy + gradient
   !g = 0.0_wp; sigma = 0.0_wp
   call drep_grad(mol,param,ep,g,sigma_tmp(:,:,3))
   ! short ranged bond energy + gradient
   call dsrb_grad(mol,param,cn,dcndr,dcndL,esrb,g,sigma_tmp(:,:,4)) ! WRONG
   !etot = ep + esrb; return
   ! h0 gradient
   allocate( dHdcn(mol%n), dHdq(mol%n), pew(nao,nao), tmp(nao), &
      &      source = 0.0_wp )
   tmp = wfn%focc*wfn%emo*evtoau
   ! setup energy weighted density matrix = pew
   call dmat(nao,tmp,wfn%C,pew)
   if (mol%npbc > 0) then
      if (ccm) then
         call ccm_build_dSH0(mol%n,basis,intcut,nao,nbf,mol%at,mol%xyz, &
            &                mol%lattice,wfn%q,cn, &
            &                wfn%P,Pew,g,sigma,dhdcn,dhdq,param%kmagic,ken, &
            &                param%alphaj,param%kcnsh,param%xbdamp,mol%wsc)
      else
         call pbc_build_dSH0(mol%n,basis,intcut,nao,nbf,mol%at,mol%xyz, &
            &                mol%lattice,intrep,wfn%q,cn, &
            &                wfn%P,Pew,g,sigma,dhdcn,dhdq,param%kmagic,ken, &
            &                param%alphaj,param%kcnsh,param%xbdamp)
      endif
      ! setup CN sigma
      call dgemv('n',9,mol%n,1.0_wp,dcndL,9,dhdcn,1,1.0_wp,sigma,1)
      ! setup  q sigma
      call dgemv('n',9,mol%n,1.0_wp,dqdL, 9, dhdq,1,1.0_wp,sigma,1)
   else
      call mol_build_dSH0(mol%n,basis,intcut,nao,nbf,mol%at,mol%xyz,wfn%q,cn, &
         &                wfn%P,Pew,g,sigma,dhdcn,dhdq,param%kmagic,ken, &
         &                param%alphaj,param%kcnsh,param%xbdamp)
   endif
   ! setup CN gradient
   call dgemv('n',3*mol%n,mol%n,-1.0_wp,dcndr,3*mol%n,dhdcn,1,1.0_wp,g,1)
   ! setup  q gradient
   call dgemv('n',3*mol%n,mol%n, 1.0_wp,dqdr, 3*mol%n, dhdq,1,1.0_wp,g,1)

   if (profile) call timer%measure(8)

!  calculate the norm for printout
   res%gnorm = sqrt(sum( g**2 ))

! ---------------------------------------
!  Properties + printout
! ---------------------------------------
   printing: if (pr) then
      ! print orbital energies and occupation numbers
      !call preig(6,wfn%focc,1.0_wp,wfn%emo, &
      !     max(wfn%ihomoa-12,1),min(wfn%ihomoa+11,nao))
      if(.not.orthog)then
        call print_orbital_eigenvalues(iunit,wfn,5)
        if ((wfn%ihomo+1.le.nao).and.(wfn%ihomo.ge.1)) &
           egap = wfn%emo(wfn%ihomo+1)-wfn%emo(wfn%ihomo)
      else
        wfn%nao = xdim
        call print_orbital_eigenvalues(iunit,wfn,5)
        if ((wfn%ihomo+1.le.xdim).and.(wfn%ihomo.ge.1)) &
           egap = wfn%emo(wfn%ihomo+1)-wfn%emo(wfn%ihomo)
        wfn%nao = nao
      endif
     

   endif printing

   sigma = sigma + sum(sigma_tmp, dim=3)
   if (debug) then
      write(iunit,'("ij",2x,6a12)') "total", "EEQ", "D4", "rep", "SRB", "H0"
      write(iunit,'("xx:",1x,6f12.7)') &
         sigma(1,1), sigma_tmp(1,1,1), sigma_tmp(1,1,2), sigma_tmp(1,1,3), sigma_tmp(1,1,4), sigma_tmp(1,1,5)
      write(iunit,'("xy:",1x,6f12.7)') &                                                                  
         sigma(1,2), sigma_tmp(1,2,1), sigma_tmp(1,2,2), sigma_tmp(1,2,3), sigma_tmp(1,2,4), sigma_tmp(1,2,5)
      write(iunit,'("xz:",1x,6f12.7)') &                                                                  
         sigma(1,3), sigma_tmp(1,3,1), sigma_tmp(1,3,2), sigma_tmp(1,3,3), sigma_tmp(1,3,4), sigma_tmp(1,3,5)
      write(iunit,'("yx:",1x,6f12.7)') &                                                                  
         sigma(2,1), sigma_tmp(2,1,1), sigma_tmp(2,1,2), sigma_tmp(2,1,3), sigma_tmp(2,1,4), sigma_tmp(2,1,5)
      write(iunit,'("yy:",1x,6f12.7)') &                                                                  
         sigma(2,2), sigma_tmp(2,2,1), sigma_tmp(2,2,2), sigma_tmp(2,2,3), sigma_tmp(2,2,4), sigma_tmp(2,2,5)
      write(iunit,'("yz:",1x,6f12.7)') &                                                                  
         sigma(2,3), sigma_tmp(2,3,1), sigma_tmp(2,3,2), sigma_tmp(2,3,3), sigma_tmp(2,3,4), sigma_tmp(2,3,5)
      write(iunit,'("zx:",1x,6f12.7)') &                                                                  
         sigma(3,1), sigma_tmp(3,1,1), sigma_tmp(3,1,2), sigma_tmp(3,1,3), sigma_tmp(3,1,4), sigma_tmp(3,1,5)
      write(iunit,'("zy:",1x,6f12.7)') &                                                                  
         sigma(3,2), sigma_tmp(3,2,1), sigma_tmp(3,2,2), sigma_tmp(3,2,3), sigma_tmp(3,2,4), sigma_tmp(3,2,5)
      write(iunit,'("zz:",1x,6f12.7)') &                                                                  
         sigma(3,3), sigma_tmp(3,3,1), sigma_tmp(3,3,2), sigma_tmp(3,3,3), sigma_tmp(3,3,4), sigma_tmp(3,3,5)
   endif

! ---------------------------------------
!  Save all the energy contributions
! ---------------------------------------
   etot = eel + ees + ep + exb + esrb
   eat = eatoms*evtoau - etot
   etot = etot + ed

   res%e_elec  = eel
   res%e_atom  = eat
   res%e_rep   = ep
   res%e_es    = ees
   res%e_aes   = 0.0_wp
   res%e_axc   = 0.0_wp
   res%e_xb    = esrb
   res%e_disp  = ed
   res%g_hb    = ehb
   res%e_total = etot
   res%hl_gap  = egap
   if (lgbsa) then
      res%g_solv  = gsolv
      !res%g_born  = gborn    ! not returned
      res%g_sasa  = gbsa%gsasa
      !res%g_hb    = gbsa%ghb ! not returned
      res%g_shift = gshift
   endif
   ! do NOT calculate the dipole moment from the density, because it's really bad
   res%dipole  = matmul(mol%xyz,wfn%q)
   res%g_solv  = 0.0_wp

   if (profile.and.pr) call timer%write(iunit,'EHT')

end associate

! ---------------------------------------
! Free memory
! ---------------------------------------
   if (profile) call timer%deallocate
   if(allocated(S))         deallocate(S)
   if(allocated(H))         deallocate(H)
   if(allocated(H0))        deallocate(H0)
   if(allocated(X))         deallocate(X)
   if(allocated(H1))        deallocate(H1)
   if(allocated(kcnao))     deallocate(kcnao)
   if(allocated(zsh))       deallocate(zsh)
   if(allocated(dqdr))      deallocate(dqdr)

end subroutine peeq

! -----------------------------------------------------------------------
!  Calculate D4 dispersion gradient
! -----------------------------------------------------------------------
subroutine ddisp_peeq(mol,param,cn,dcndr,dcndL,grd,ed,gd,sigma)
  use iso_fortran_env, only : wp => real64
! -----------------------------------------------------------------------
!  Type definitions
! -----------------------------------------------------------------------
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data
! -----------------------------------------------------------------------
!  DFT-D4 definitions and PBC definitions
! -----------------------------------------------------------------------
   use dftd4
   use eeq_model
   use pbc,    only : get_realspace_cutoff
   use ncoord
   implicit none

! -----------------------------------------------------------------------
!  Intent IN
! -----------------------------------------------------------------------
   type(tb_molecule),           intent(in)     :: mol
   type(scc_parameter),         intent(in)     :: param
   type(dftd_parameter)                        :: par
   type(chrg_parameter)                        :: chrgeq
   ! EEQ partial charges and derivatives
   real(wp), dimension(mol%n),        intent(in) :: cn
   real(wp), dimension(3,mol%n,mol%n),intent(in) :: dcndr
   real(wp), dimension(3,3,mol%n),    intent(in) :: dcndL

! -----------------------------------------------------------------------
!  Intent INOUT
! -----------------------------------------------------------------------
! dispersion energy and derivative
   real(wp),                     intent(inout) :: ed
   real(wp), dimension(3,mol%n), intent(inout) :: gd
   real(wp), dimension(3,3),     intent(inout) :: sigma

!--- allocatables
   real(wp),allocatable, dimension(:,:)   :: sdum
   real(wp),allocatable, dimension(:,:)   :: gdum
   real(wp),allocatable, dimension(:)     :: q
   real(wp),allocatable, dimension(:,:,:) :: dqdr
   real(wp),allocatable, dimension(:,:,:) :: dqdL
   real(wp),allocatable, dimension(:)     :: covcn
   real(wp),allocatable, dimension(:,:,:) :: dcovcndr
   real(wp),allocatable, dimension(:,:,:) :: dcovcndL
   real(wp),allocatable, dimension(:,:)   :: c6abns
   real(wp),allocatable, dimension(:)     :: gw
   real(wp),allocatable, dimension(:,:)   :: gdummy

! -----------------------------------------------------------------------
!  Variables
! -----------------------------------------------------------------------
   integer               :: i,j,k,l,m
   integer               :: ndim
   integer               :: mbd
   logical, intent(in)   :: grd
   integer, dimension(3) :: rep_vdw,rep_cn
 ! real space cutoffs
   real(wp), parameter   :: cn_thr = 1600.0_wp
   real(wp), parameter   :: crit_vdw         = 4000.0_wp

 ! damping variable
   real(wp) :: edum
   real(wp) :: t6
   real(wp) :: t8
   real(wp) :: expterm

! -----------------------------------------------------------------------
!  Initialization
! -----------------------------------------------------------------------
   ed  = 0.0_wp
   mbd = 0

! -----------------------------------------------------------------------
!  Get ndim
! -----------------------------------------------------------------------
   call d4dim(mol%n,mol%at,ndim)
   if (mol%npbc > 0) &
   call get_realspace_cutoff(mol%lattice,crit_vdw,rep_vdw)

! -----------------------------------------------------------------------
!  Get memory
! -----------------------------------------------------------------------
  allocate(c6abns(ndim,ndim));   c6abns = 0.0_wp
  allocate(gw(ndim));                gw = 0.0_wp
  allocate(q(mol%n));                 q = 0.0_wp
  allocate(dqdr(3,mol%n,mol%n+1)); dqdr = 0.0_wp
  allocate(dqdL(3,3,mol%n+1));     dqdL = 0.0_wp
  allocate(covcn(mol%n));               covcn = 0.0_wp
  allocate(dcovcndr(3,mol%n,mol%n)); dcovcndr = 0.0_wp
  allocate(dcovcndL(3,3,mol%n));     dcovcndL = 0.0_wp

  ! get D4(EEQ) charges
  call new_charge_model_2019(chrgeq,mol%n,mol%at)
  ! neither sdum nor gdum need to be dummy allocated, since lgrad = .false. is set
  ! eeq_chrgeq must not reference sdum/gdum and we can avoid the dummy allocate
  call eeq_chrgeq(mol,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,edum,gdum,sdum, &
     &            .false.,.false.,.true.)

  if (mol%npbc > 0) then
     call get_d4_cn(mol,covcn,dcovcndr,dcovcndL,thr=cn_thr)
  !else
     !call dncoord_d4(mol%n,mol%at,mol%xyz,covcn,dcovcndr,cn_thr)
  endif

  ! setup c6abns with diagonal terms: i interaction with its images
  call pbc_d4(mol%n,ndim,mol%at,param%wf,param%g_a,param%g_c,covcn,gw,c6abns)

! -----------------------------------------------------------------------
!  Set dispersion parameters and calculate Edisp or Gradient
! -----------------------------------------------------------------------
   if (mol%npbc > 0) then
      call dispgrad_3d(mol,ndim,q,covcn,dcovcndr,dcovcndL,rep_vdw,rep_vdw,crit_vdw,crit_vdw, &
         &             param%disp,param%wf,param%g_a,param%g_c,c6abns,mbd, &
         &             gd,sigma,ed,dqdr,dqdL)
   else
      call dispgrad(mol%n,ndim,mol%at,q,mol%xyz, &
         &          param%disp,param%wf,param%g_a,param%g_c, &
         &          c6abns,mbd, &
         &          gd,ed,dqdr)
   endif

! -----------------------------------------------------------------------
!  Free willy (no this is not ORCA)
! -----------------------------------------------------------------------
  if(allocated(c6abns))  deallocate(c6abns)
  if(allocated(gw))      deallocate(gw)
  if(allocated(covcn))    deallocate(covcn)
  if(allocated(dcovcndr)) deallocate(dcovcndr)
  if(allocated(dcovcndL)) deallocate(dcovcndL)
  if(allocated(q))       deallocate(q)
  if(allocated(dqdr))    deallocate(dqdr)
  if(allocated(dqdL))    deallocate(dqdL)
end subroutine ddisp_peeq

! repulsion
pure subroutine drep_grad(mol,param,erep,g,sigma)
   use aoparam,   only : rep,en
   use tbdef_molecule
   use tbdef_param
   use pbc_tools
   use pbc
   implicit none
   ! intent in
   type(scc_parameter), intent(in) :: param
   type(tb_molecule),   intent(in) :: mol
   ! intent inout
   real(wp), intent(inout), dimension(3,mol%n) :: g
   real(wp), intent(inout), dimension(3,3)     :: sigma
   ! intent out
   real(wp), intent(out)                   :: erep
   ! local variables
   integer :: i,j,k,lin
   integer :: iat,jat,ati,atj
   real(wp), dimension(3) :: ri,dr
   real(wp), dimension(3) :: rij
   real(wp)               :: r2,r,r2top34,r2top34top2
   real(wp)               :: den2,den4
   real(wp)               :: alpha
   real(wp)               :: repab
   real(wp)               :: expterm
   real(wp)               :: dtmp
   real(wp), parameter    :: rthr = 1600.0_wp
   real(wp)               :: w,t(3)
   integer                :: latrep(3),tx,ty,tz
   call get_realspace_cutoff(mol%lattice,rthr,latrep)
   w = 1.0_wp
   ! initialize
   erep = 0.0_wp
   if (mol%npbc > 0) then
      do i = 1, mol%n
         ati=mol%at(i)
         do j = 1, i
            do concurrent(tx = -latrep(1):latrep(1), &
                  &       ty = -latrep(2):latrep(2), &
                  &       tz = -latrep(3):latrep(3), &
                  &       i.ne.j .or. (tx.ne.0 .or. ty.ne.0 .or. tz.ne.0))
               t = [tx,ty,tz]
               rij = mol%xyz(:,i) - (mol%xyz(:,j) + matmul(mol%lattice,t))
               r2 = sum(rij**2)
               if(r2.gt.rthr) cycle
               atj=mol%at(j)
               r = sqrt(r2)
               den2 = (en(ati) - en(atj))**2
               den4 = den2**2
               alpha=sqrt(rep(1,ati)*rep(1,atj))&
                  *(1.0_wp+(0.01_wp*den2+0.01_wp*den4)*param%xbrad)
               repab = rep(2,ati)*rep(2,atj)
               r2top34 = r2**0.75_wp
               r2top34top2 = r2top34**2
               expterm = exp(-alpha*r2top34)*repab
               ! save repulsion energy
               erep = erep + expterm/r * w
               ! save repulsion gradient
               dtmp = expterm*(1.5_wp*alpha*r2top34 + 1)/r2top34top2 * w
               g(:,i) = g(:,i) - dtmp*rij
               g(:,j) = g(:,j) + dtmp*rij
               sigma = sigma - dtmp*outer_prod_3x3(rij,rij)
            enddo ! k WSC partner
         enddo ! j atom
      enddo ! i atom
   else
      do i = 1, mol%n
         ati=mol%at(i)
         do j = 1, i-1
            rij = mol%xyz(:,i) - mol%xyz(:,j)
            r2 = sum(rij**2)
            if(r2.gt.rthr) cycle
            atj=mol%at(j)
            r = sqrt(r2)
            den2 = (en(ati) - en(atj))**2
            den4 = den2**2
            alpha=sqrt(rep(1,ati)*rep(1,atj))&
               *(1.0_wp+(0.01_wp*den2+0.01_wp*den4)*param%xbrad)
            repab = rep(2,ati)*rep(2,atj)
            r2top34 = r2**0.75_wp
            r2top34top2 = r2top34**2
            expterm = exp(-alpha*r2top34)*repab
            ! save repulsion energy
            erep = erep + expterm/r
            ! save repulsion gradient
            dtmp = expterm*(1.5_wp*alpha*r2top34 + 1)&
               /r2top34top2
            g(:,i) = g(:,i) - dtmp*rij
            g(:,j) = g(:,j) + dtmp*rij
         enddo ! j atom
      enddo ! i atom
   endif

end subroutine drep_grad

! short-ranged bond correction
pure subroutine dsrb_grad(mol,param,cn,dcndr,dcndL,esrb,g,sigma)

   use tbdef_param
   use tbdef_molecule

   use aoparam,   only : rep

   use approxrab
   use ncoord
   use pbc_tools
   use pbc

   implicit none

   ! intent in
   type(tb_molecule),                intent(in) :: mol
   type(scc_parameter),              intent(in) :: param
   ! intent inout
   real(wp),     dimension(3,mol%n), intent(inout) :: g
   real(wp),     dimension(3,3),     intent(inout) :: sigma
   ! intent out
   real(wp), intent(out) :: esrb
   ! local variables
   real(wp), dimension(:),     intent(in) :: cn
   real(wp), dimension(:,:,:), intent(in) :: dcndr
   real(wp), dimension(:,:,:), intent(in) :: dcndL
   integer :: i,j,k,wsAt
   integer :: iat,jat,ati,atj
   integer :: nsrb,lin
   integer :: dx,dy,dz
   real(wp), dimension(4) :: kcn
   real(wp) :: xbrad
   real(wp) :: gscal
   real(wp) :: den
   real(wp) :: expterm
   ! distances
   real(wp), dimension(3) :: rij
   real(wp)               :: r2,r,dr,rab
   real(wp)               :: dtmp,pre
   integer                :: tx,ty,tz,latrep(3)
   real(wp)               :: w,t(3)
   ! allocatables
   real(wp), allocatable, dimension(:,:,:) :: drab0dr
   real(wp), allocatable, dimension(:,:,:) :: drab0dL
   real(wp), allocatable, dimension(:)     :: rab0
   integer,  allocatable, dimension(:,:)   :: srblist
   ! initialize
   kcn=param%kcnsh
   xbrad=param%xbrad
   gscal=param%gscal
   esrb = 0.0_wp

   if (mol%npbc > 0) call get_realspace_cutoff(mol%lattice,200.0_wp,latrep)
   w = 1.0_wp

   call build_srblist(mol,nsrb,srblist)
   if (nsrb.eq.0) return

   ! get memory
   allocate( drab0dr(3,mol%n,nsrb), source = 0.0_wp )
   allocate( drab0dL(3,3,nsrb),     source = 0.0_wp )
   allocate( rab0(nsrb),            source = 0.0_wp )
   ! get approximated distances rab and gradients
   periodic: if (mol%npbc > 0) then
      call pbc_approx_rab(mol%n,mol%at,mol%xyz,cn,dcndr,dcndL,nsrb,srblist,kcn(2), &
         &            rab0,drab0dr,drab0dL)
      do i = 1, nsrb
         iat = srblist(1,i)
         jat = srblist(2,i)
         ati = mol%at(iat)
         atj = mol%at(jat)
         den = en(ati) - en(atj)
         pre = kcn(4)*(1.0_wp + gscal*den**2)
         do concurrent(tx = -latrep(1):latrep(1), &
               &       ty = -latrep(2):latrep(2), &
               &       tz = -latrep(3):latrep(3))
            t = [tx,ty,tz]
            rij = mol%xyz(:,iat) - (mol%xyz(:,jat) + matmul(mol%lattice,t))
            rab = norm2(rij)
            dr = rab - rab0(i)
            expterm = kcn(3)*exp(-pre*dr**2)
            ! save SRB energy
            esrb = esrb + expterm * w
            dtmp = 2.0_wp*pre*dr*expterm * w
            g(:,iat) = g(:,iat) - dtmp*rij/rab
            g(:,jat) = g(:,jat) + dtmp*rij/rab
            ! three body gradient
            g(:,:) = g(:,:) - dtmp*drab0dr(:,:,i)
            sigma = sigma + dtmp*(drab0dL(:,:,i) - outer_prod_3x3(rij,rij)/rab)
         enddo ! rep
      enddo ! i

   else
      call approx_rab(mol%n,mol%at,mol%xyz,cn,dcndr,nsrb,srblist,kcn(2),rab0,drab0dr)
      do i = 1, nsrb
         iat = srblist(1,i)
         jat = srblist(2,i)
         ati = mol%at(iat)
         atj = mol%at(jat)
         den = en(ati) - en(atj)
         pre = kcn(4)*(1.0_wp + gscal*den**2)
         rij = mol%xyz(:,iat) - mol%xyz(:,jat)
         rab = norm2(rij)
         dr = rab - rab0(i)
         expterm = kcn(3)*exp(-pre*dr**2)
         ! save SRB energy
         esrb = esrb + expterm
         dtmp = 2.0_wp*pre*dr*expterm
         g(:,iat) = g(:,iat) - dtmp*rij/rab
         g(:,jat) = g(:,jat) + dtmp*rij/rab
         ! three body gradient
         g(:,:) = g(:,:) - dtmp*drab0dr(:,:,i)
      enddo ! i

   endif periodic

   ! free memory
   if(allocated(drab0dr))   deallocate(drab0dr)
   if(allocated(drab0dl))   deallocate(drab0dl)
   if(allocated(rab0))    deallocate(rab0)
   if(allocated(srblist)) deallocate(srblist)

end subroutine dsrb_grad

pure subroutine build_srblist(mol,nsrb,srblist)
   use tbdef_molecule
   implicit none
   type(tb_molecule), intent(in) :: mol
   integer, intent(out) :: nsrb
   integer, allocatable, intent(out) :: srblist(:,:)
   integer  :: i,j,k
   real(wp) :: r2
   ! cutoff
   real(wp), parameter    :: srb_cut = 200.0_wp
   nsrb = 0
   do i = 1, mol%n-1
      if (srbatom(mol%at(i))) cycle    ! i i case
      do j = i+1, mol%n
         if (srbatom(mol%at(j))) cycle ! j j case
         if (mol%at(i).eq.mol%at(j)) cycle ! hetero only
         r2  = mol%dist(j,i)**2 ! get minimum image distance
         ! increment array size according to SRB cutoff criterium
         ! default: srb_cut = 200.0_wp
         if(r2.lt.srb_cut) nsrb=nsrb+1
      enddo ! j
   enddo ! i
   ! nsrb = 0: no pairs for SRB correction, nothing to do here
   if(nsrb.eq.0) return
   ! get memory
   allocate( srblist(2,nsrb), source = 0 )
   ! setup srblist
   k=0
   do i = 1, mol%n-1
      if (srbatom(mol%at(i))) cycle    ! i i case
      do j = i+1, mol%n
         if (srbatom(mol%at(j))) cycle ! j j case
         if (mol%at(i).eq.mol%at(j)) cycle ! hetero only
         r2  = mol%dist(j,i)**2 ! get minimum image distance
         if(r2.gt.srb_cut) cycle
         k = k + 1
         srblist(1,k) = i
         srblist(2,k) = j
      enddo ! j
   enddo ! i
   nsrb = k

contains

pure elemental function srbatom(iat) result(lsrb)
   integer, intent(in) :: iat
   logical :: lsrb
   lsrb = iat < 5 .or. iat > 9
end function srbatom

end subroutine build_srblist

! ------------------------------------------------------------------------
!  Calculate the periodic AO overlap matrix
! ------------------------------------------------------------------------
subroutine mol_build_SH0(nat,at,basis,nbf,nao,xyz,q,cn,intcut, &
      &                  kmagic,ken,alphaj,kcn,xbdamp,sint,h0)

   use tbdef_basisset

   use aoparam

   use lin_mod, only : lin
   use intgrad
   use scc_core

   implicit none

   type(tb_basisset), intent(in) :: basis

   integer, intent(in)  :: nat
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nao
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: q(nat)
   real(wp),intent(in)  :: cn(nat)
   real(wp),intent(in)  :: intcut

   real(wp),intent(in)  :: kmagic(:,:)
   real(wp),intent(in)  :: ken(:,:)
   real(wp),intent(in)  :: kcn(:)
   real(wp),intent(in)  :: alphaj
   real(wp),intent(in)  :: xbdamp

   real(wp),intent(out) :: sint(nao,nao)  ! overlap matrix S
   real(wp),intent(out) :: h0(nao*(1+nao)/2)

   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij
   real(wp) tmp1,tmp2,tmp3,tmp4,step,step2
   real(wp) dx,dy,dz,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,f,ci,cc,cj,alpi,rab2,ab,est

   real(wp) :: den,den2,den4,enpoly
   real(wp) :: zi,zj,zetaij
   real(wp) :: hii,hjj,hav
   real(wp) :: shpoly,km
   logical  :: valaoi,valaoj

   real(wp)  ri(3),rj(3),f1,f2
   real(wp), parameter :: point(3) = 0.0_wp
   real(wp) dtmp(3),qtmp(6),ss(6,6),dd(3,6,6),qq(6,6,6),tmp(6,6)
   real(wp) sstmp(6,6) !TESTST
   integer ip,jp,iat,jat,ati,atj,ish,jsh,icao,jcao,iao,jao,jshmax,il,jl
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer itt(0:3)

   integer kat

   parameter(itt  =(/0,1,4,10/))
   real(wp) :: saw(10)

   sint = 0.0_wp
   h0   = 0.0_wp

   !$omp parallel default(none) &
   !$omp private(iat,jat,ij,ati,cc,ci,rab2,atj,ish,ishtyp,valaoi,valaoj, &
   !$omp&        ri,rj,icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp, &
   !$omp&        ss,saw,est,alpi,alpj,ab,iprim,jprim,ip,jp,km,shpoly, &
   !$omp&        mli,mlj,tmp,zi,zj,zetaij,enpoly,iao,jao, &
   !$omp&        ii,jj,k,den,den2,den4,i,j,il,jl,hii,hjj,hav) &
   !$omp reduction (+:sint,h0) &
   !$omp shared(basis,at,ao_n,ao_l,xyz,intcut,nat,kqat,kcnat,cn,q,en, &
   !$omp        ken,gam3,xbdamp,kcn,kmagic,kpair,alphaj)
   !$omp do schedule(runtime)
   do iat = 1, nat
      ri  = xyz(:,iat)
      ati = at(iat)
      do jat = 1, iat-1
         atj = at(jat)

         den = en(ati) - en(atj)
         den2 =  den**2
         den4 = den2**2

         ishells: do ish = 1, ao_n(ati)
            ishtyp = ao_l(ish,ati)
            icao = basis%caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = ao_n(atj)
            if (iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jshtyp = ao_l(jsh,atj)
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               ! get indices
               i = 1+basis%saoshell(ish,iat)
               j = 1+basis%saoshell(jsh,jat)
               il = shell(basis%lao2(i)) 
               jl = shell(basis%lao2(j))
               ! diagonals are the same for all H0 elements
               hii = basis%hdiag2(i) - kcnat(il-1,ati)*cn(iat) &
                  &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
               hjj = basis%hdiag2(j) - kcnat(jl-1,atj)*cn(jat) &
                  &  - kqat(jl,atj)*q(jat) - gam3(atj)*q(jat)**2

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + ken(jl,il)*den2 &
                  &      + xbdamp*ken(jl,il)*den4)

               ! we scale the two shells depending on their exponent
               zi = basis%aoexp(i)
               zj = basis%aoexp(j)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = kmagic(jl,il) * kpair(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = basis%valao2(i).eq.0
               valaoj = basis%valao2(j).eq.0
               ! and scale appropiately
               if (valaoi) then
                  if (valaoj) then
                     km = km * kcn(1)
                  else
                     km = km * alphaj
                  endif
               else
                  if (valaoj) km = km * alphaj
               endif

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * km * (hii + hjj)

               ss = 0.0_wp
               rj = xyz(:,jat)

               ! distance dependent polynomial
               shpoly = rfactor(il,jl,ati,atj,ri,rj)

               ! get overlap integral
               call get_overlap(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj,point, &
                  &             intcut,basis%nprim,basis%primcount, &
                  &             basis%alp,basis%cont,ss)

               !transform from CAO to SAO
               call dtrf2(ss,ishtyp,jshtyp)

               do ii = 1, llao2(ishtyp)
                  iao = ii+basis%saoshell(ish,iat)
                  do jj = 1, llao2(jshtyp)
                     jao = jj+basis%saoshell(jsh,jat)
                     if (jao > iao) cycle
                     ij = lin(iao,jao)
                     ! add all WSC images
                     sint(iao,jao) = sint(iao,jao) + ss(jj,ii)
                     sint(jao,iao) = sint(jao,iao) + ss(jj,ii)
                     ! Hamiltonian
                     H0(ij) = H0(ij) + hav * shpoly * ss(jj,ii)
                  enddo
               enddo

            enddo jshells
         enddo ishells
      enddo
   enddo
   !$OMP END DO
   !$OMP END PARALLEL

   ! diagonal elements
   do i = 1, nao
      sint(i,i)=1.0_wp+sint(i,i)

      ii  = i*(1+i)/2 ! H0 is packed, note i*(i-1)/2+i = i*(1+i)/2
      iat = basis%aoat2(i)
      ati = at(iat)
      il  = shell(basis%lao2(i))

      ! calculate environment dependent shift
      hii = basis%hdiag2(i) - kcnat(il-1,ati)*cn(iat) &
         &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
      H0(ii) = hii

   end do

end subroutine mol_build_SH0

! ------------------------------------------------------------------------
!  Calculate the periodic AO overlap matrix
! ------------------------------------------------------------------------
subroutine ccm_build_SH0(nat,at,basis,nbf,nao,xyz,lattice,q,cn,intcut, &
      &                  kmagic,ken,alphaj,kcn,xbdamp,sint,h0,wsc)

   use tbdef_basisset
   use tbdef_wsc

   use aoparam

   use lin_mod, only : lin
   use intgrad
   use scc_core

   implicit none

   type(tb_wsc), intent(in) :: wsc
   type(tb_basisset), intent(in) :: basis

   integer, intent(in)  :: nat
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nao
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: lattice(3,3)
   real(wp),intent(in)  :: q(nat)
   real(wp),intent(in)  :: cn(nat)
   real(wp),intent(in)  :: intcut

   real(wp),intent(in)  :: kmagic(:,:)
   real(wp),intent(in)  :: ken(:,:)
   real(wp),intent(in)  :: kcn(:)
   real(wp),intent(in)  :: alphaj
   real(wp),intent(in)  :: xbdamp

   real(wp),intent(out) :: sint(nao,nao)  ! overlap matrix S
   real(wp),intent(out) :: h0(nao*(1+nao)/2)

   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij
   real(wp) tmp1,tmp2,tmp3,tmp4,step,step2
   real(wp) dx,dy,dz,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,f,ci,cc,cj,alpi,rab2,ab,est

   real(wp) :: den,den2,den4,enpoly
   real(wp) :: zi,zj,zetaij
   real(wp) :: hii,hjj,hav
   real(wp) :: shpoly,km
   logical  :: valaoi,valaoj

   real(wp)  ri(3),rj(3),t(3),f1,f2
   real(wp),parameter ::point(3) = 0.0_wp
   real(wp) dtmp(3),qtmp(6),ss(6,6),dd(3,6,6),qq(6,6,6),tmp(6,6)
   real(wp) sstmp(6,6) !TESTST
   integer ip,jp,iat,jat,ati,atj,ish,jsh,icao,jcao,iao,jao,jshmax,il,jl
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer itt(0:3)

   integer kat

   parameter(itt  =(/0,1,4,10/))
   real(wp) :: saw(10)

   sint = 0.0_wp
   h0   = 0.0_wp

   !$omp parallel default(none) &
   !$omp private(iat,jat,ij,ati,cc,ci,rab2,atj,ish,ishtyp,valaoi,valaoj, &
   !$omp&        ri,rj,icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp, &
   !$omp&        ss,saw,est,alpi,alpj,ab,iprim,jprim,ip,jp,km,shpoly, &
   !$omp&        mli,mlj,tmp,zi,zj,zetaij,enpoly,iao,jao, &
   !$omp&        ii,jj,k,den,den2,den4,i,j,il,jl,hii,hjj,hav,t) &
   !$omp reduction (+:sint,h0) &
   !$omp shared(wsc,basis,at,ao_n,ao_l,xyz,lattice,intcut,nat,kqat,kcnat,cn,q,en, &
   !$omp        ken,gam3,xbdamp,kcn,kmagic,kpair,alphaj)
   !$omp do schedule(runtime)
   do iat = 1, nat
      ri  = xyz(:,iat)
      ati = at(iat)
      do jat = 1, iat-1
         atj = at(jat)

         den = en(ati) - en(atj)
         den2 =  den**2
         den4 = den2**2

         ishells: do ish = 1, ao_n(ati)
            ishtyp = ao_l(ish,ati)
            icao = basis%caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = ao_n(atj)
            if (iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jshtyp = ao_l(jsh,atj)
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               ! get indices
               i = 1+basis%saoshell(ish,iat)
               j = 1+basis%saoshell(jsh,jat)
               il = shell(basis%lao2(i)) 
               jl = shell(basis%lao2(j))
               ! diagonals are the same for all H0 elements
               hii = basis%hdiag2(i) - kcnat(il-1,ati)*cn(iat) &
                  &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
               hjj = basis%hdiag2(j) - kcnat(jl-1,atj)*cn(jat) &
                  &  - kqat(jl,atj)*q(jat) - gam3(atj)*q(jat)**2

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + ken(jl,il)*den2 &
                  &      + xbdamp*ken(jl,il)*den4)

               ! we scale the two shells depending on their exponent
               zi = basis%aoexp(i)
               zj = basis%aoexp(j)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = kmagic(jl,il) * kpair(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = basis%valao2(i).eq.0
               valaoj = basis%valao2(j).eq.0
               ! and scale appropiately
               if (valaoi) then
                  if (valaoj) then
                     km = km * kcn(1)
                  else
                     km = km * alphaj
                  endif
               else
                  if (valaoj) km = km * alphaj
               endif

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * km * (hii + hjj)

               wscAt: do kat = 1, wsc%itbl(jat,iat)
                  ss = 0.0_wp
                  t = wsc%lattr(:,kat,jat,iat)
                  rj = xyz(:,jat) + matmul(lattice,t)

                  ! distance dependent polynomial
                  shpoly = rfactor(il,jl,ati,atj,ri,rj)

                  ! get overlap integral
                  call get_overlap(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj,point, &
                     &             intcut,basis%nprim,basis%primcount, &
                     &             basis%alp,basis%cont,ss)

                  !transform from CAO to SAO
                  call dtrf2(ss,ishtyp,jshtyp)

                  do ii = 1, llao2(ishtyp)
                     iao = ii+basis%saoshell(ish,iat)
                     do jj = 1, llao2(jshtyp)
                        jao = jj+basis%saoshell(jsh,jat)
                        if (jao > iao) cycle
                        ij = lin(iao,jao)
                        ! add all WSC images
                        sint(iao,jao) = sint(iao,jao) + ss(jj,ii)*wsc%w(jat,iat)
                        sint(jao,iao) = sint(jao,iao) + ss(jj,ii)*wsc%w(jat,iat)
                        ! Hamiltonian
                        H0(ij) = H0(ij) + hav * shpoly * ss(jj,ii) * wsc%w(jat,iat)
                     enddo
                  enddo

               enddo wscAt
            enddo jshells
         enddo ishells
      enddo
   enddo
   !$OMP END DO
   !$OMP END PARALLEL

   ! diagonal elements
   do i = 1, nao
      sint(i,i)=1.0_wp+sint(i,i)

      ii  = i*(1+i)/2 ! H0 is packed, note i*(i-1)/2+i = i*(1+i)/2
      iat = basis%aoat2(i)
      ati = at(iat)
      il  = shell(basis%lao2(i))

      ! calculate environment dependent shift
      hii = basis%hdiag2(i) - kcnat(il-1,ati)*cn(iat) &
         &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
      H0(ii) = hii

   end do

end subroutine ccm_build_SH0

! ------------------------------------------------------------------------
!  Calculate the periodic AO overlap matrix
! ------------------------------------------------------------------------
subroutine pbc_build_SH0(nat,at,basis,nbf,nao,xyz,lat,latrep,q,cn,intcut, &
      &                  kmagic,ken,alphaj,kcn,xbdamp,sint,h0)

   use tbdef_basisset

   use aoparam

   use lin_mod, only : lin
   use intgrad
   use scc_core

   implicit none

   type(tb_basisset), intent(in) :: basis

   integer, intent(in)  :: nat
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nao
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: lat(3,3)
   integer, intent(in)  :: latrep(3)
   real(wp),intent(in)  :: q(nat)
   real(wp),intent(in)  :: cn(nat)
   real(wp),intent(in)  :: intcut

   real(wp),intent(in)  :: kmagic(:,:)
   real(wp),intent(in)  :: ken(:,:)
   real(wp),intent(in)  :: kcn(:)
   real(wp),intent(in)  :: alphaj
   real(wp),intent(in)  :: xbdamp

   real(wp),intent(out) :: sint(nao,nao)  ! overlap matrix S
   real(wp),intent(out) :: h0(nao*(1+nao)/2)

   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij
   real(wp) tmp1,tmp2,tmp3,tmp4,step,step2
   real(wp) dx,dy,dz,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,f,ci,cc,cj,alpi,rab2,ab,est

   real(wp) :: den,den2,den4,enpoly
   real(wp) :: zi,zj,zetaij
   real(wp) :: hii,hjj,hav
   real(wp) :: shpoly,km
   logical  :: valaoi,valaoj

   real(wp)  ri(3),rj(3),f1,f2
   real(wp), parameter :: point(3) = 0.0_wp
   real(wp) dtmp(3),qtmp(6),ss(6,6),dd(3,6,6),qq(6,6,6),tmp(6,6)
   real(wp) sstmp(6,6) !TESTST
   integer ip,jp,iat,jat,ati,atj,ish,jsh,icao,jcao,iao,jao,jshmax,il,jl
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer itt(0:3)

   integer kat

   parameter(itt  =(/0,1,4,10/))
   real(wp) :: saw(10)
   integer  :: tx,ty,tz
   real(wp) :: t(3),w

   w = 1.0_wp !/ real(product(2*latrep+1),wp)

   sint = 0.0_wp
   h0   = 0.0_wp

   !$omp parallel default(none) &
   !$omp private(iat,jat,ij,ati,cc,ci,rab2,atj,ish,ishtyp,valaoi,valaoj, &
   !$omp&        ri,rj,icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp, &
   !$omp&        ss,saw,est,alpi,alpj,ab,iprim,jprim,ip,jp,km,shpoly, &
   !$omp&        mli,mlj,tmp,zi,zj,zetaij,enpoly,iao,jao, &
   !$omp&        ii,jj,k,den,den2,den4,i,j,il,jl,hii,hjj,hav,tx,ty,tz,t) &
   !$omp reduction (+:sint,h0) &
   !$omp shared(basis,at,ao_n,ao_l,xyz,intcut,nat,kqat,kcnat,cn,q,en, &
   !$omp        ken,gam3,xbdamp,kcn,kmagic,kpair,alphaj,lat,latrep,w)
   !$omp do schedule(runtime)
   do iat = 1, nat
      ri  = xyz(:,iat)
      ati = at(iat)
      do jat = 1, iat-1
         atj = at(jat)

         den = en(ati) - en(atj)
         den2 =  den**2
         den4 = den2**2

         ishells: do ish = 1, ao_n(ati)
            ishtyp = ao_l(ish,ati)
            icao = basis%caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = ao_n(atj)
            if (iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jshtyp = ao_l(jsh,atj)
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               ! get indices
               i = 1+basis%saoshell(ish,iat)
               j = 1+basis%saoshell(jsh,jat)
               il = shell(basis%lao2(i)) 
               jl = shell(basis%lao2(j))
               ! diagonals are the same for all H0 elements
               hii = basis%hdiag2(i) - kcnat(il-1,ati)*cn(iat) &
                  &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
               hjj = basis%hdiag2(j) - kcnat(jl-1,atj)*cn(jat) &
                  &  - kqat(jl,atj)*q(jat) - gam3(atj)*q(jat)**2

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + ken(jl,il)*den2 &
                  &      + xbdamp*ken(jl,il)*den4)

               ! we scale the two shells depending on their exponent
               zi = basis%aoexp(i)
               zj = basis%aoexp(j)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = kmagic(jl,il) * kpair(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = basis%valao2(i).eq.0
               valaoj = basis%valao2(j).eq.0
               ! and scale appropiately
               if (valaoi) then
                  if (valaoj) then
                     km = km * kcn(1)
                  else
                     km = km * alphaj
                  endif
               else
                  if (valaoj) km = km * alphaj
               endif

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * km * (hii + hjj)

               do concurrent(tx = -latrep(1):latrep(1), &
                     &       ty = -latrep(2):latrep(2), &
                     &       tz = -latrep(3):latrep(3))
                  t = [tx,ty,tz]
                  rj = xyz(:,jat) + matmul(lat,t)
                  ss = 0.0_wp

                  ! distance dependent polynomial
                  shpoly = rfactor(il,jl,ati,atj,ri,rj)

                  ! get overlap integral
                  call get_overlap(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj,point, &
                     &             intcut,basis%nprim,basis%primcount, &
                     &             basis%alp,basis%cont,ss)

                  !transform from CAO to SAO
                  call dtrf2(ss,ishtyp,jshtyp)

                  do ii = 1, llao2(ishtyp)
                     iao = ii+basis%saoshell(ish,iat)
                     do jj = 1, llao2(jshtyp)
                        jao = jj+basis%saoshell(jsh,jat)
                        if (jao > iao) cycle
                        ij = lin(iao,jao)
                        ! add all WSC images
                        sint(iao,jao) = sint(iao,jao) + ss(jj,ii)*w
                        sint(jao,iao) = sint(jao,iao) + ss(jj,ii)*w
                        ! Hamiltonian
                        H0(ij) = H0(ij) + hav * shpoly * ss(jj,ii)*w
                     enddo
                  enddo

               enddo
            enddo jshells
         enddo ishells
      enddo
   enddo
   !$OMP END DO
   !$OMP END PARALLEL

   ! diagonal elements
   do i = 1, nao
      sint(i,i)=1.0_wp+sint(i,i)

      ii  = i*(1+i)/2 ! H0 is packed, note i*(i-1)/2+i = i*(1+i)/2
      iat = basis%aoat2(i)
      ati = at(iat)
      il  = shell(basis%lao2(i))

      ! calculate environment dependent shift
      hii = basis%hdiag2(i) - kcnat(il-1,ati)*cn(iat) &
         &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
      H0(ii) = hii

   end do

end subroutine pbc_build_SH0

pure subroutine get_overlap(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj,point,intcut, &
      &                nprim,primcount,alp,cont,sint)
   use intgrad
   implicit none
   integer, intent(in)  :: icao
   integer, intent(in)  :: jcao
   integer, intent(in)  :: naoi
   integer, intent(in)  :: naoj
   integer, intent(in)  :: iptyp
   integer, intent(in)  :: jptyp
   real(wp),intent(in)  :: ri(3)
   real(wp),intent(in)  :: rj(3)
   real(wp),intent(in)  :: point(3)
   real(wp),intent(in)  :: intcut
   real(wp),intent(out) :: sint(:,:)

   integer, intent(in)  :: nprim(:)
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)

   integer  :: ip,iprim,mli,jp,jprim,mlj
   real(wp) :: rij(3),rij2,alpi,alpj,ci,cj,cc
   real(wp) :: ab,est,saw(10)

   real(wp),parameter :: max_r2 = 2000.0_wp

   sint = 0.0_wp

   rij = ri - rj
   rij2 = rij(1)**2 + rij(2)**2 + rij(3)**2

   if(rij2.gt.max_r2) return

   do ip = 1, nprim(icao+1)
      iprim = ip + primcount(icao+1)
      ! exponent the same for each l component
      alpi = alp(iprim)
      do jp = 1, nprim(jcao+1)
         jprim=jp+primcount(jcao+1)
         ! exponent the same for each l component
         alpj=alp(jprim)
         ab=1.0_wp/(alpi+alpj)
         est=rij2*alpi*alpj*ab
         if(est.gt.intcut) cycle
         ! now compute integrals
         do mli=1,naoi
            iprim = ip + primcount(icao+mli)
            ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
            ci = cont(iprim)
            do mlj=1,naoj
               jprim = jp + primcount(jcao+mlj)
               saw = 0.0_wp
               ! prim-prim  integrals
               call build_sdq_ints(ri,rj,point,alpi,alpj,iptyp+mli,jptyp+mlj,saw)
               cc = cont(jprim)*ci
               sint(mlj,mli) = sint(mlj,mli)+saw(1)*cc! pbc_w(jat,iat)
            enddo ! mlj
         enddo ! mli
      enddo ! jp
   enddo ! ip

end subroutine get_overlap

pure subroutine get_grad_overlap(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj,point,intcut, &
      &                     nprim,primcount,alp,cont,sdq,sdqg)
   use intgrad
   implicit none
   integer, intent(in)  :: icao
   integer, intent(in)  :: jcao
   integer, intent(in)  :: naoi
   integer, intent(in)  :: naoj
   integer, intent(in)  :: iptyp
   integer, intent(in)  :: jptyp
   real(wp),intent(in)  :: ri(3)
   real(wp),intent(in)  :: rj(3)
   real(wp),intent(in)  :: point(3)
   real(wp),intent(in)  :: intcut
   real(wp),intent(out) :: sdq(:,:)
   real(wp),intent(out) :: sdqg(:,:,:)

   integer, intent(in)  :: nprim(:)
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)

   integer  :: ip,iprim,mli,jp,jprim,mlj
   real(wp) :: rij(3),rij2,alpi,alpj,ci,cj,cc
   real(wp) :: ab,est,saw,sawg(3)

   real(wp),parameter :: max_r2 = 2000.0_wp

   sdqg = 0.0_wp
   sdq  = 0.0_wp

   rij = ri - rj
   rij2 = rij(1)**2 + rij(2)**2 + rij(3)**2

   if(rij2.gt.max_r2) return

   do ip=1,nprim(icao+1)
      iprim=ip+primcount(icao+1)
      ! exponent the same for each l component
      alpi=alp(iprim)
      do jp=1,nprim(jcao+1)
         jprim=jp+primcount(jcao+1)
         ! exponent the same for each l component
         alpj=alp(jprim)
         est=alpi*alpj*rij2/(alpi+alpj)
         if(est.gt.intcut) cycle
         !--------------- compute gradient ----------
         ! now compute integrals  for different components of i(e.g., px,py,pz)
         do mli=1,naoi
            iprim=ip+primcount(icao+mli)
            ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
            ci=cont(iprim)
            do mlj=1,naoj
               jprim=jp+primcount(jcao+mlj)
               cc=cont(jprim)*ci
               saw=0;sawg=0
               call build_ds_ints(ri,rj,alpi,alpj,iptyp+mli,jptyp+mlj,saw,sawg)
               sdq(mlj,mli) = sdq(mlj,mli)+saw*cc
               sdqg(:,mlj,mli) = sdqg(:,mlj,mli) &
                  & + sawg(:)*cc
            enddo ! mlj : Cartesian component of j prims
         enddo  ! mli : Cartesian component of i prims
      enddo ! jp : loop over j prims
   enddo  ! ip : loop over i prims
end subroutine get_grad_overlap

! ------------------------------------------------------------------------
!  Calculate the gradient resulting from a periodic AO overlap matrix
! ------------------------------------------------------------------------
subroutine mol_build_dSH0(nat,basis,thr,nao,nbf,at,xyz,q,cn,P,Pew,g,sigma, &
      &                   dHdcn,dHdq,kmagic,ken,alphaj,kcn,xbdamp)
   use mctc_constants, only : pi
   use mctc_econv

   use tbdef_wsc
   use tbdef_basisset

   use aoparam

   use intgrad
   use grad_core

   implicit none

   ! intent in
   integer, intent(in)      :: nat
   type(tb_basisset), intent(in) :: basis
   real(wp),intent(in)      :: thr
   integer, intent(in)      :: nao
   integer, intent(in)      :: nbf
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   real(wp),intent(in) :: q(nat)
   real(wp),intent(in) :: cn(nat)
   real(wp),intent(in) :: P(nao,nao)
   real(wp),intent(in) :: Pew(nao,nao)

   real(wp),intent(in)  :: kmagic(:,:)
   real(wp),intent(in)  :: ken(:,:)
   real(wp),intent(in)  :: kcn(:)
   real(wp),intent(in)  :: alphaj
   real(wp),intent(in)  :: xbdamp
   ! intent inout
   real(wp),intent(inout) :: dHdq(nat)
   real(wp),intent(inout) :: dHdcn(nat)
   real(wp),intent(inout) :: g(3,nat)
   real(wp),intent(inout) :: sigma(3,3)

   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

   integer itt(0:3)
   parameter (itt=(/0,1,4,10/))
   real(wp) tmp1,tmp2,tmp3,step,step2,step3,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4
   real(wp) thr2,f,ci,cc,cj,alpi,rab2,ab,est
   real(wp) f1,f2,tmp(6,6),rij(3),ri(3),rj(3),rij2
   real(wp) stmp,ral(3,3),rar(3,3),rbl(3,3),pre
   real(wp) dtmp,qtmp,rbr(3,3),r2l(3),r2r(3),qqa(6,6,6,3)
   real(wp)  ss(6,6,3),ddc(3,6,6,3),qqc(6,6,6,3),dda(3,6,6,3)
   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij,lin,jshmax
   integer ip,jp,iat,jat,kat,ati,atj,ish,jsh,icao,jcao,iao,jao,ixyz,jxyz
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   real(wp) :: dumdum(3),dum,sdq(6,6),sdqg(3,6,6)
   real(wp),parameter :: point(3) = 0.0_wp

   integer  :: il,jl
   real(wp) :: Hij,Pii,Pij,HPij,H0sr
   real(wp) :: den,den2,den4,enpoly
   real(wp) :: zi,zj,zetaij
   real(wp) :: hii,hjj,hav
   real(wp) :: shpoly,dshpoly(3),km
   logical  :: valaoi,valaoj
   real(wp) :: g_xyz(3)

   !$omp parallel default(none) &
   !$omp private(iat,jat,ixyz,ati,cc,ci,rab2,atj,ish,ishtyp,valaoi,valaoj,g_xyz, &
   !$omp&        ri,rj,icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp,rij2, &
   !$omp&        sdq,sdqg,est,alpi,alpj,ab,iprim,jprim,ip,jp,km,shpoly,dshpoly, &
   !$omp&        mli,mlj,dum,dumdum,tmp,stmp,rij,zi,zj,zetaij,enpoly,iao,jao, &
   !$omp&        ii,jj,k,den,den2,den4,i,j,il,jl,hii,hjj,hav,Pij,Hij,HPij,H0sr) &
   !$omp reduction (+:g,sigma,dhdcn,dhdq) &
   !$omp shared(basis,at,ao_n,ao_l,xyz,thr,nat,kqat,kcnat,cn,q,en, &
   !$omp        ken,gam3,xbdamp,kcn,kmagic,kpair,alphaj,P,Pew)
   !$omp do schedule(runtime)
   do iat = 1, nat
      ri  = xyz(:,iat)
      ati = at(iat)
      do jat = 1, iat-1
         atj = at(jat)

         den = en(ati) - en(atj)
         den2 =  den**2
         den4 = den2**2

         ishells: do ish = 1, ao_n(ati)
            ishtyp = ao_l(ish,ati)
            icao = basis%caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = ao_n(atj)
            if(iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jshtyp = ao_l(jsh,atj)
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               ! get indices
               i = 1+basis%saoshell(ish,iat)
               j = 1+basis%saoshell(jsh,jat)
               il = shell(basis%lao2(i)) 
               jl = shell(basis%lao2(j))
               ! diagonals are the same for all H0 elements
               hii = basis%hdiag2(i) - kcnat(il-1,ati)*cn(iat) &
                  &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
               hjj = basis%hdiag2(j) - kcnat(jl-1,atj)*cn(jat) &
                  &  - kqat(jl,atj)*q(jat) - gam3(atj)*q(jat)**2

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + ken(jl,il)*den2 &
                  &      + xbdamp*ken(jl,il)*den4)

               ! we scale the two shells depending on their exponent
               zi = basis%aoexp(i)
               zj = basis%aoexp(j)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = kmagic(jl,il) * kpair(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = basis%valao2(i).eq.0
               valaoj = basis%valao2(j).eq.0
               ! and scale appropiately
               if (valaoi) then
                  if (valaoj) then
                     km = km * kcn(1)
                  else
                     km = km * alphaj
                  endif
               else
                  if (valaoj) km = km * alphaj
               endif

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * km * (hii + hjj)

               rj = xyz(:,jat)
               rij = ri - rj
               rij2 = sum(rij**2)

               ! distance dependent polynomial
               !shpoly = rfactor(il,jl,ati,atj,ri,rj)
               call drfactor(il,jl,0,ati,atj,rij2,ri,rj,shpoly,dshpoly)

               call get_grad_overlap(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj, &
                  &                  point,thr,basis%nprim,basis%primcount, &
                  &                  basis%alp,basis%cont,sdq,sdqg)

               call dtrf2(sdq,ishtyp,jshtyp)
               do ixyz = 1, 3
                  ! transform from CAO to SAO, only transform the gradient
                  tmp = sdqg(ixyz,:,:)
                  call dtrf2(tmp,ishtyp,jshtyp)
                  sdqg(ixyz,:,:) = tmp
               enddo

               do ii = 1, llao2(ishtyp)
                  iao = ii+basis%saoshell(ish,iat)
                  do jj = 1, llao2(jshtyp)
                     jao = jj+basis%saoshell(jsh,jat)

                     Pij = P(jao,iao) * evtoau

                     ! Hamiltonian element without overlap
                     Hij  = hav * shpoly
                     HPij = Hij * Pij

                     ! overlap contribution
                     stmp = 2*(HPij - Pew(jao,iao))

                     ! distance dependent polynomial contribution
                     H0sr = 2*HPij*sdq(jj,ii)
                     g_xyz = stmp*sdqg(:,jj,ii) + H0sr*dshpoly/shpoly

                     ! distribute to atom i
                     g(:,iat) = g(:,iat) + g_xyz

                     ! distribute to atom j
                     g(:,jat) = g(:,jat) - g_xyz

                     ! reduce to strain
                     sigma(:,1) = sigma(:,1) + g_xyz*rij(1)
                     sigma(:,2) = sigma(:,2) + g_xyz*rij(2)
                     sigma(:,3) = sigma(:,3) + g_xyz*rij(3)

                     ! Hamiltonian without Hav
                     HPij = km * shpoly * Pij * sdq(jj,ii)
                     ! save dE/dCN for CNi
                     dhdcn(iat) = dhdcn(iat) - HPij*kcnat(il-1,ati)
                     ! save dE/dCN for CNj
                     dhdcn(jat) = dhdcn(jat) - HPij*kcnat(jl-1,atj)

                     ! save dE/dq for qi
                     dhdq(iat) = dhdq(iat) - HPij*kqat(il,ati) &
                        &                  - HPij*gam3(ati)*2*q(iat)
                     ! save dE/dq for qj
                     dhdq(jat) = dhdq(jat) - HPij*kqat(jl,atj) &
                        &                  - HPij*gam3(atj)*2*q(jat)

                  enddo
               enddo

            enddo jshells
         enddo  ishells

      enddo ! jat
   enddo  ! iat
   !$omp end do
   !$omp end parallel

   ! diagonal contributions
   do i = 1, nao
      ii  = i*(1+i)/2 ! H0 is packed, note i*(i-1)/2+i = i*(1+i)/2
      iat = basis%aoat2(i)
      ati = at(iat)
      il  = shell(basis%lao2(i))

      Pii = P(i,i)*evtoau

      ! save dE/dCN for CNi
      dhdcn(iat) = dhdcn(iat) - Pii*kcnat(il-1,ati)
      ! save dE/dq for qi
      dhdq(iat) = dhdq(iat) - Pii*kqat(il,ati) - Pii*gam3(ati)*2*q(iat)
   enddo

end subroutine mol_build_dSH0

! ------------------------------------------------------------------------
!  Calculate the gradient resulting from a periodic AO overlap matrix
! ------------------------------------------------------------------------
subroutine ccm_build_dSH0(nat,basis,thr,nao,nbf,at,xyz,lattice,q,cn,P,Pew,g,sigma,&
      &                   dHdcn,dHdq,kmagic,ken,alphaj,kcn,xbdamp,wsc)
   use mctc_constants, only : pi
   use mctc_econv

   use tbdef_wsc
   use tbdef_basisset

   use aoparam

   use intgrad
   use grad_core

   implicit none

   ! intent in
   integer, intent(in)      :: nat
   type(tb_wsc), intent(in) :: wsc
   type(tb_basisset), intent(in) :: basis
   real(wp),intent(in)      :: thr
   integer, intent(in)      :: nao
   integer, intent(in)      :: nbf
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   real(wp),intent(in) :: lattice(3,3)
   real(wp),intent(in) :: q(nat)
   real(wp),intent(in) :: cn(nat)
   real(wp),intent(in) :: P(nao,nao)
   real(wp),intent(in) :: Pew(nao,nao)

   real(wp),intent(in)  :: kmagic(:,:)
   real(wp),intent(in)  :: ken(:,:)
   real(wp),intent(in)  :: kcn(:)
   real(wp),intent(in)  :: alphaj
   real(wp),intent(in)  :: xbdamp
   ! intent inout
   real(wp),intent(inout) :: dHdq(nat)
   real(wp),intent(inout) :: dHdcn(nat)
   real(wp),intent(inout) :: g(3,nat)
   real(wp),intent(inout) :: sigma(3,3)

   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

   integer itt(0:3)
   parameter (itt=(/0,1,4,10/))
   real(wp) tmp1,tmp2,tmp3,step,step2,step3,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4
   real(wp) thr2,f,ci,cc,cj,alpi,rab2,ab,est
   real(wp) f1,f2,tmp(6,6),rij(3),ri(3),rj(3),rij2,t(3)
   real(wp) stmp,ral(3,3),rar(3,3),rbl(3,3),pre
   real(wp) dtmp,qtmp,rbr(3,3),r2l(3),r2r(3),qqa(6,6,6,3)
   real(wp)  ss(6,6,3),ddc(3,6,6,3),qqc(6,6,6,3),dda(3,6,6,3)
   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij,lin,jshmax
   integer ip,jp,iat,jat,kat,ati,atj,ish,jsh,icao,jcao,iao,jao,ixyz,jxyz
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   real(wp) :: dumdum(3),dum,sdq(6,6),sdqg(3,6,6)
   real(wp),parameter :: point(3) = 0.0_wp

   integer  :: il,jl
   real(wp) :: Hij,Pii,Pij,HPij,H0sr
   real(wp) :: den,den2,den4,enpoly
   real(wp) :: zi,zj,zetaij
   real(wp) :: hii,hjj,hav
   real(wp) :: shpoly,dshpoly(3),km
   logical  :: valaoi,valaoj
   real(wp) :: g_xyz(3)

   !$omp parallel default(none) &
   !$omp private(iat,jat,ixyz,ati,cc,ci,rab2,atj,ish,ishtyp,valaoi,valaoj,g_xyz, &
   !$omp&        ri,rj,icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp,rij2, &
   !$omp&        sdq,sdqg,est,alpi,alpj,ab,iprim,jprim,ip,jp,km,shpoly,dshpoly, &
   !$omp&        mli,mlj,dum,dumdum,tmp,stmp,rij,zi,zj,zetaij,enpoly,iao,jao,t, &
   !$omp&        ii,jj,k,den,den2,den4,i,j,il,jl,hii,hjj,hav,Pij,Hij,HPij,H0sr) &
   !$omp reduction (+:g,sigma,dhdcn,dhdq) &
   !$omp shared(wsc,basis,at,ao_n,ao_l,xyz,lattice,thr,nat,kqat,kcnat,cn,q,en, &
   !$omp        ken,gam3,xbdamp,kcn,kmagic,kpair,alphaj,P,Pew)
   !$omp do schedule(runtime)
   do iat = 1, nat
      ri  = xyz(:,iat)
      ati = at(iat)
      do jat = 1, iat-1
         atj = at(jat)

         den = en(ati) - en(atj)
         den2 =  den**2
         den4 = den2**2

         ishells: do ish = 1, ao_n(ati)
            ishtyp = ao_l(ish,ati)
            icao = basis%caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = ao_n(atj)
            if(iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jshtyp = ao_l(jsh,atj)
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               ! get indices
               i = 1+basis%saoshell(ish,iat)
               j = 1+basis%saoshell(jsh,jat)
               il = shell(basis%lao2(i)) 
               jl = shell(basis%lao2(j))
               ! diagonals are the same for all H0 elements
               hii = basis%hdiag2(i) - kcnat(il-1,ati)*cn(iat) &
                  &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
               hjj = basis%hdiag2(j) - kcnat(jl-1,atj)*cn(jat) &
                  &  - kqat(jl,atj)*q(jat) - gam3(atj)*q(jat)**2

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + ken(jl,il)*den2 &
                  &      + xbdamp*ken(jl,il)*den4)

               ! we scale the two shells depending on their exponent
               zi = basis%aoexp(i)
               zj = basis%aoexp(j)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = kmagic(jl,il) * kpair(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = basis%valao2(i).eq.0
               valaoj = basis%valao2(j).eq.0
               ! and scale appropiately
               if (valaoi) then
                  if (valaoj) then
                     km = km * kcn(1)
                  else
                     km = km * alphaj
                  endif
               else
                  if (valaoj) km = km * alphaj
               endif

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * km * (hii + hjj)

               wscAt: do kat = 1,wsc%itbl(jat,iat)
                  t = wsc%lattr(:,kat,jat,iat)
                  rj = xyz(:,jat) + matmul(lattice,t)
                  rij = ri - rj
                  rij2 = sum(rij**2)

                  ! distance dependent polynomial
                  !shpoly = rfactor(il,jl,ati,atj,ri,rj)
                  call drfactor(il,jl,0,ati,atj,rij2,ri,rj,shpoly,dshpoly)

                  call get_grad_overlap(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj, &
                     &                  point,thr,basis%nprim,basis%primcount, &
                     &                  basis%alp,basis%cont,sdq,sdqg)

                  call dtrf2(sdq,ishtyp,jshtyp)
                  do ixyz = 1, 3
                     ! transform from CAO to SAO, only transform the gradient
                     tmp = sdqg(ixyz,:,:)
                     call dtrf2(tmp,ishtyp,jshtyp)
                     sdqg(ixyz,:,:) = tmp
                  enddo

                  do ii = 1, llao2(ishtyp)
                     iao = ii+basis%saoshell(ish,iat)
                     do jj = 1, llao2(jshtyp)
                        jao = jj+basis%saoshell(jsh,jat)

                        Pij = P(jao,iao) * evtoau

                        ! Hamiltonian element without overlap
                        Hij  = hav * shpoly
                        HPij = Hij * Pij

                        ! overlap contribution
                        stmp = 2*wsc%w(jat,iat) * (HPij - Pew(jao,iao))

                        ! distance dependent polynomial contribution
                        H0sr = 2*HPij*sdq(jj,ii)*wsc%w(jat,iat)
                        g_xyz = stmp*sdqg(:,jj,ii) + H0sr*dshpoly/shpoly

                        ! distribute to atom i
                        g(:,iat) = g(:,iat) + g_xyz

                        ! distribute to atom j
                        g(:,jat) = g(:,jat) - g_xyz

                        ! reduce to strain
                        sigma(:,1) = sigma(:,1) + g_xyz*rij(1)
                        sigma(:,2) = sigma(:,2) + g_xyz*rij(2)
                        sigma(:,3) = sigma(:,3) + g_xyz*rij(3)

                        ! Hamiltonian without Hav
                        HPij = km * shpoly * Pij * sdq(jj,ii) * wsc%w(jat,iat)
                        ! save dE/dCN for CNi
                        dhdcn(iat) = dhdcn(iat) - HPij*kcnat(il-1,ati)
                        ! save dE/dCN for CNj
                        dhdcn(jat) = dhdcn(jat) - HPij*kcnat(jl-1,atj)

                        ! save dE/dq for qi
                        dhdq(iat) = dhdq(iat) - HPij*kqat(il,ati) &
                           &                  - HPij*gam3(ati)*2*q(iat)
                        ! save dE/dq for qj
                        dhdq(jat) = dhdq(jat) - HPij*kqat(jl,atj) &
                           &                  - HPij*gam3(atj)*2*q(jat)

                     enddo
                  enddo

               enddo wscAt

            enddo jshells
         enddo  ishells

      enddo ! jat
   enddo  ! iat
   !$omp end do
   !$omp end parallel

   ! diagonal contributions
   do i = 1, nao
      ii  = i*(1+i)/2 ! H0 is packed, note i*(i-1)/2+i = i*(1+i)/2
      iat = basis%aoat2(i)
      ati = at(iat)
      il  = shell(basis%lao2(i))

      Pii = P(i,i)*evtoau

      ! save dE/dCN for CNi
      dhdcn(iat) = dhdcn(iat) - Pii*kcnat(il-1,ati)
      ! save dE/dq for qi
      dhdq(iat) = dhdq(iat) - Pii*kqat(il,ati) - Pii*gam3(ati)*2*q(iat)
   enddo

end subroutine ccm_build_dSH0

! ------------------------------------------------------------------------
!  Calculate the gradient resulting from a periodic AO overlap matrix
! ------------------------------------------------------------------------
subroutine pbc_build_dSH0(nat,basis,thr,nao,nbf,at,xyz,lat,latrep,q,cn,P,Pew,g,sigma, &
      &                   dHdcn,dHdq,kmagic,ken,alphaj,kcn,xbdamp)
   use mctc_constants, only : pi
   use mctc_econv

   use tbdef_wsc
   use tbdef_basisset

   use aoparam

   use intgrad
   use grad_core

   implicit none

   ! intent in
   integer, intent(in)      :: nat
   type(tb_basisset), intent(in) :: basis
   real(wp),intent(in)      :: thr
   integer, intent(in)      :: nao
   integer, intent(in)      :: nbf
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   real(wp),intent(in) :: lat(3,nat)
   integer, intent(in) :: latrep(3)
   real(wp),intent(in) :: q(nat)
   real(wp),intent(in) :: cn(nat)
   real(wp),intent(in) :: P(nao,nao)
   real(wp),intent(in) :: Pew(nao,nao)

   real(wp),intent(in)  :: kmagic(:,:)
   real(wp),intent(in)  :: ken(:,:)
   real(wp),intent(in)  :: kcn(:)
   real(wp),intent(in)  :: alphaj
   real(wp),intent(in)  :: xbdamp
   ! intent inout
   real(wp),intent(inout) :: dHdq(nat)
   real(wp),intent(inout) :: dHdcn(nat)
   real(wp),intent(inout) :: g(3,nat)
   real(wp),intent(inout) :: sigma(3,3)

   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

   integer itt(0:3)
   parameter (itt=(/0,1,4,10/))
   real(wp) tmp1,tmp2,tmp3,step,step2,step3,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4
   real(wp) thr2,f,ci,cc,cj,alpi,rab2,ab,est
   real(wp) f1,f2,tmp(6,6),rij(3),ri(3),rj(3),rij2
   real(wp) stmp,ral(3,3),rar(3,3),rbl(3,3),pre
   real(wp) dtmp,qtmp,rbr(3,3),r2l(3),r2r(3),qqa(6,6,6,3)
   real(wp)  ss(6,6,3),ddc(3,6,6,3),qqc(6,6,6,3),dda(3,6,6,3)
   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij,lin,jshmax
   integer ip,jp,iat,jat,kat,ati,atj,ish,jsh,icao,jcao,iao,jao,ixyz,jxyz
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   real(wp) :: dumdum(3),dum,sdq(6,6),sdqg(3,6,6)
   real(wp),parameter :: point(3) = 0.0_wp

   integer  :: il,jl
   real(wp) :: Hij,Pii,Pij,HPij,H0sr
   real(wp) :: den,den2,den4,enpoly
   real(wp) :: zi,zj,zetaij
   real(wp) :: hii,hjj,hav
   real(wp) :: shpoly,dshpoly(3),km
   logical  :: valaoi,valaoj
   real(wp) :: g_xyz(3)
   integer  :: tx,ty,tz
   real(wp) :: t(3),w

   w = 1.0_wp !/ real(product(2*latrep+1),wp)

   !$omp parallel default(none) &
   !$omp private(iat,jat,ixyz,ati,cc,ci,rab2,atj,ish,ishtyp,valaoi,valaoj,g_xyz, &
   !$omp&        ri,rj,icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp,rij2, &
   !$omp&        sdq,sdqg,est,alpi,alpj,ab,iprim,jprim,ip,jp,km,shpoly,dshpoly, &
   !$omp&        mli,mlj,dum,dumdum,tmp,stmp,rij,zi,zj,zetaij,enpoly,iao,jao, &
   !$omp&        ii,jj,k,den,den2,den4,i,j,il,jl,hii,hjj,hav,Pij,Hij,HPij,H0sr, &
   !$omp&        tx,ty,tz,t) &
   !$omp reduction (+:g,sigma,dhdcn,dhdq) &
   !$omp shared(basis,at,ao_n,ao_l,xyz,thr,nat,kqat,kcnat,cn,q,en, &
   !$omp        ken,gam3,xbdamp,kcn,kmagic,kpair,alphaj,P,Pew,w,lat,latrep)
   !$omp do schedule(runtime)
   do iat = 1, nat
      ri  = xyz(:,iat)
      ati = at(iat)
      do jat = 1, iat-1
         atj = at(jat)

         den = en(ati) - en(atj)
         den2 =  den**2
         den4 = den2**2

         ishells: do ish = 1, ao_n(ati)
            ishtyp = ao_l(ish,ati)
            icao = basis%caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = ao_n(atj)
            if(iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jshtyp = ao_l(jsh,atj)
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               ! get indices
               i = 1+basis%saoshell(ish,iat)
               j = 1+basis%saoshell(jsh,jat)
               il = shell(basis%lao2(i)) 
               jl = shell(basis%lao2(j))
               ! diagonals are the same for all H0 elements
               hii = basis%hdiag2(i) - kcnat(il-1,ati)*cn(iat) &
                  &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
               hjj = basis%hdiag2(j) - kcnat(jl-1,atj)*cn(jat) &
                  &  - kqat(jl,atj)*q(jat) - gam3(atj)*q(jat)**2

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + ken(jl,il)*den2 &
                  &      + xbdamp*ken(jl,il)*den4)

               ! we scale the two shells depending on their exponent
               zi = basis%aoexp(i)
               zj = basis%aoexp(j)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = kmagic(jl,il) * kpair(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = basis%valao2(i).eq.0
               valaoj = basis%valao2(j).eq.0
               ! and scale appropiately
               if (valaoi) then
                  if (valaoj) then
                     km = km * kcn(1)
                  else
                     km = km * alphaj
                  endif
               else
                  if (valaoj) km = km * alphaj
               endif

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * km * (hii + hjj)

               do concurrent(tx = -latrep(1):latrep(1), &
                     &       ty = -latrep(2):latrep(2), &
                     &       tz = -latrep(3):latrep(3))
                  t = [tx,ty,tz]
                  rj = xyz(:,jat) + matmul(lat,t)
                  rij = ri - rj
                  rij2 = sum(rij**2)

                  ! distance dependent polynomial
                  !shpoly = rfactor(il,jl,ati,atj,ri,rj)
                  call drfactor(il,jl,0,ati,atj,rij2,ri,rj,shpoly,dshpoly)

                  call get_grad_overlap(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj, &
                     &                  point,thr,basis%nprim,basis%primcount, &
                     &                  basis%alp,basis%cont,sdq,sdqg)

                  call dtrf2(sdq,ishtyp,jshtyp)
                  do ixyz = 1, 3
                     ! transform from CAO to SAO, only transform the gradient
                     tmp = sdqg(ixyz,:,:)
                     call dtrf2(tmp,ishtyp,jshtyp)
                     sdqg(ixyz,:,:) = tmp
                  enddo

                  do ii = 1, llao2(ishtyp)
                     iao = ii+basis%saoshell(ish,iat)
                     do jj = 1, llao2(jshtyp)
                        jao = jj+basis%saoshell(jsh,jat)

                        Pij = P(jao,iao) * evtoau

                        ! Hamiltonian element without overlap
                        Hij  = hav * shpoly
                        HPij = Hij * Pij

                        ! overlap contribution
                        stmp = 2*w * (HPij - Pew(jao,iao))

                        ! distance dependent polynomial contribution
                        H0sr = 2*HPij*sdq(jj,ii)*w
                        g_xyz = stmp*sdqg(:,jj,ii) + H0sr*dshpoly/shpoly

                        ! distribute to atom i
                        g(:,iat) = g(:,iat) + g_xyz

                        ! distribute to atom j
                        g(:,jat) = g(:,jat) - g_xyz

                        ! reduce to strain
                        sigma(:,1) = sigma(:,1) + g_xyz*rij(1)
                        sigma(:,2) = sigma(:,2) + g_xyz*rij(2)
                        sigma(:,3) = sigma(:,3) + g_xyz*rij(3)

                        ! Hamiltonian without Hav
                        HPij = km * shpoly * Pij * sdq(jj,ii) * w
                        ! save dE/dCN for CNi
                        dhdcn(iat) = dhdcn(iat) - HPij*kcnat(il-1,ati)
                        ! save dE/dCN for CNj
                        dhdcn(jat) = dhdcn(jat) - HPij*kcnat(jl-1,atj)

                        ! save dE/dq for qi
                        dhdq(iat) = dhdq(iat) - HPij*kqat(il,ati) &
                           &                  - HPij*gam3(ati)*2*q(iat)
                        ! save dE/dq for qj
                        dhdq(jat) = dhdq(jat) - HPij*kqat(jl,atj) &
                           &                  - HPij*gam3(atj)*2*q(jat)

                     enddo
                  enddo

               enddo

            enddo jshells
         enddo  ishells

      enddo ! jat
   enddo  ! iat
   !$omp end do
   !$omp end parallel

   ! diagonal contributions
   do i = 1, nao
      ii  = i*(1+i)/2 ! H0 is packed, note i*(i-1)/2+i = i*(1+i)/2
      iat = basis%aoat2(i)
      ati = at(iat)
      il  = shell(basis%lao2(i))

      Pii = P(i,i)*evtoau

      ! save dE/dCN for CNi
      dhdcn(iat) = dhdcn(iat) - Pii*kcnat(il-1,ati)
      ! save dE/dq for qi
      dhdq(iat) = dhdq(iat) - Pii*kqat(il,ati) - Pii*gam3(ati)*2*q(iat)
   enddo

end subroutine pbc_build_dSH0

end module peeq_module

