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

module xtb_peeq
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   use xtb_mctc_la
   use xtb_xtb_data
   use xtb_chargemodel
   use xtb_type_latticepoint, only : TLatticePoint, init
   use xtb_disp_coordinationnumber, only : getCoordinationNumber, cnType, &
      & cutCoordinationNumber
   use xtb_disp_dftd4, only : d4_gradient
   use xtb_disp_encharges, only : getENCharges
   use xtb_coulomb_gaussian
   use xtb_intgrad, only : get_grad_overlap, get_overlap
   use xtb_xtb_eeq
   use xtb_xtb_hamiltonian, only : getSelfEnergy
   use xtb_type_wsc, only : tb_wsc
   implicit none
   private

   public :: peeq

   !> print the different gradient contributions (for debugging)
   logical,parameter,private :: gpr =.false.

   !> profiling
   logical,parameter,private :: profile = .true.

contains

subroutine peeq &
      (env,mol,wfn,basis,xtbData,gbsa,egap,et,prlevel,grd,ccm,acc,etot,gradient,sigma,res)

! ------------------------------------------------------------------------
!  Class definitions
! ------------------------------------------------------------------------
   use xtb_solv_gbsa, only : TBorn
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_wavefunction
   use xtb_type_basisset
   use xtb_type_param
   use xtb_type_data
   use xtb_type_timer

! ------------------------------------------------------------------------
!  Get interfaces
! ------------------------------------------------------------------------
   use xtb_readin
   use xtb_aespot
   use xtb_scc_core
   use xtb_grad_core
   use xtb_eeq
   use xtb_disp_ncoord
   use xtb_lineardep
   use xtb_pbc

   implicit none

   character(len=*), parameter :: source = 'peeq'

! ------------------------------------------------------------------------
!  INPUT
! ------------------------------------------------------------------------
   type(TEnvironment), intent(inout)    :: env
   type(TMolecule),  intent(in) :: mol     !< molecular structure infomation
   type(TBasisset),  intent(in) :: basis   !< basis set
   type(TxTBData), intent(in) :: xtbData
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
   real(wp),intent(inout),dimension(3,mol%n) :: gradient    !< molecular gradient
   type(TWavefunction),intent(inout)       :: wfn  !< TB-wavefunction
   type(TBorn),allocatable,intent(inout) :: gbsa
! ------------------------------------------------------------------------
!  OUTPUT
! ------------------------------------------------------------------------
   type(scc_results),intent(out) :: res !< bundles all calculation informations

! ------------------------------------------------------------------------
   type(chrg_parameter) :: chrgeq
   real(wp), allocatable :: cn(:), dcndr(:,:,:), dcndL(:,:,:)
   real(wp), allocatable :: ccn(:), dccndr(:,:,:), dccndL(:,:,:)
   real(wp), allocatable :: X    (:,:)
   real(wp), allocatable :: S    (:,:)
   real(wp), allocatable :: H    (:,:)
   real(wp), allocatable :: H0   (:)
   real(wp), allocatable :: H1   (:)
   real(wp), allocatable :: zsh  (:)
   real(wp), allocatable :: kcnao(:)
   real(wp), allocatable :: kqao (:)
   real(wp), allocatable :: selfEnergy(:, :)
   real(wp), allocatable :: dSEdcn(:, :)
   real(wp), allocatable :: dSEdq(:, :)
   real(wp), allocatable :: Pa   (:,:)
   real(wp), allocatable :: Pb   (:,:)
   integer :: rep_cn(3)
   integer :: nid
   integer, allocatable :: idnum(:)
   real(wp),allocatable :: chargeWidth(:, :)
   type(TGaussianSmeared) :: coulomb
   type(TENEquilibration) :: eeq

! ------------------------------------------------------------------------
   ! local variable (without comment)
   integer :: ii,jj,i,j,k,m,iat,jat,ati,atj,kk
   integer :: il,jl

   integer,parameter :: lladr(4)    = [1,3,6,10]
   integer,parameter :: lladr2(0:3) = [1,3,5,7 ]

   real(wp) :: dipol(3)
   real(wp) :: eatoms,eat
   real(wp) :: dum

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
   ! energies
   real(wp) :: eel,ed,ees,exb,ep

! ---------------------------------------
!  PEEQ information
! ---------------------------------------
   real(wp) , allocatable, dimension(:)     :: rab0
   real(wp) , allocatable, dimension(:,:,:) :: grab0
   real(wp) , allocatable, dimension(:)     :: lcn
   real(wp) , allocatable, dimension(:,:,:) :: dlcn
   real(wp),               dimension(2)     :: hbpi
   real(wp),               dimension(2)     :: hbpj
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
   logical :: exitRun

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
   type(TLatticePoint) :: latp
   real(wp), allocatable :: trans(:, :)

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

   interface
      subroutine generate_wsc(mol,wsc)
         import :: TMolecule, tb_wsc
         type(TMolecule), intent(in) :: mol
         type(tb_wsc),    intent(inout) :: wsc
      end subroutine generate_wsc
   end interface

   type(tb_wsc) :: wsc

! ---------------------------------------
!  EEQ/GBSA information
! ---------------------------------------
   real(wp) :: gsolv

   associate( nao => basis%nao, &
         &    nbf => basis%nbf, &
         &    nshell => basis%nshell, &
         &    naop => basis%nao*(basis%nao+1)/2, &
         &    nbfp => basis%nbf*(basis%nbf+1)/2)

   if (ccm) &
      call generate_wsc(mol, wsc)

   if (profile) then
      if (allocated(gbsa)) then
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
   gradient = 0.0_wp

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

   call init(latp, env, mol, 60.0_wp)

! ---------------------------------------
!  IMPORTANT FACT: H is given in eV
! ---------------------------------------

! ---------------------------------------
!  Get memory
! ---------------------------------------
   allocate(cn(mol%n));                     cn = 0.0_wp
   allocate(dcndr(3,mol%n,mol%n));       dcndr = 0.0_wp
   allocate(dcndL(3,3,mol%n));           dcndL = 0.0_wp
   allocate(ccn(mol%n));                   ccn = 0.0_wp
   allocate(dccndr(3,mol%n,mol%n));     dccndr = 0.0_wp
   allocate(dccndL(3,3,mol%n));         dccndL = 0.0_wp
   allocate(H(nao,nao));                     H = 0.0_wp
   allocate(X(nao,nao));                     X = 0.0_wp
   allocate(H0(naop));                      H0 = 0.0_wp
   allocate(S(nao,nao));                     S = 0.0_wp
   allocate(H1(naop));                      H1 = 0.0_wp
   allocate(kcnao(nao));                 kcnao = 0.0_wp
   allocate(kqao(nao));                   kqao = 0.0_wp
   allocate(zsh(nshell));                  zsh = 0.0_wp
   allocate(selfEnergy(maxval(xtbData%nshell), mol%n))
   allocate(dSEdcn(maxval(xtbData%nshell), mol%n))
   allocate(dSEdq(maxval(xtbData%nshell), mol%n))
   allocate(qeeq(mol%n));           qeeq = 0.0_wp
   allocate(dqdr(3,mol%n,mol%n)); dqdr = 0.0_wp
   allocate(dqdL(3,3,mol%n));     dqdL = 0.0_wp

   wfn%nel = nint(sum(mol%z) - mol%chrg)
   wfn%nopen = mol%uhf
   if(wfn%nopen == 0 .and. mod(wfn%nel,2) /= 0) wfn%nopen=1

! ---------------------------------------
!  Fill levels
! ---------------------------------------
   call setzshell(xtbData,mol%n,mol%at,nshell,mol%z,zsh,eatoms,0)
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

   if (prlevel > 1) then
      write(env%unit,'(/,10x,51("."))')
      write(env%unit,'(10x,":",22x,a,22x,":")') "SETUP"
      write(env%unit,'(10x,":",49("."),":")')
      write(env%unit,intfmt) "# basis functions  ",basis%nbf
      write(env%unit,intfmt) "# atomic orbitals  ",basis%nao
      write(env%unit,intfmt) "# shells           ",basis%nshell
      write(env%unit,intfmt) "# electrons        ",wfn%nel
      if (mol%npbc > 0) &
      write(env%unit,chrfmt) "PBC by CCM         ",bool2string(ccm)
      write(env%unit,dblfmt) "electronic temp.   ",et,      "K   "
      write(env%unit,dblfmt) "accuracy           ",acc,     "    "
      write(env%unit,scifmt) "-> integral cutoff ",intcut,  "    "
      write(env%unit,scifmt) "-> integral neglect",neglect, "    "
      write(env%unit,'(10x,51("."))')
   endif

   if (profile) call timer%measure(1)
   if (profile) call timer%measure(2,"Coordination number")
! ---------------------------------------
!  Get CN(1:n) + dcndr(3,1:n,1:n) under pbc
! ---------------------------------------
   call latp%getLatticePoints(trans, 40.0_wp)
   call getCoordinationNumber(mol, trans, 40.0_wp, cnType%erf, &
      & cn, dcndr, dcndL)
   call cutCoordinationNumber(mol%n, cn, dcndr, dcndL, maxCN=8.0_wp)

   if (profile) call timer%measure(2)
   if (profile) call timer%measure(4,"D4 Dispersion")
! ----------------------------------------
!  D4 dispersion energy + gradient (2B) under pbc
! ----------------------------------------
   call getENCharges(env, mol, cn, dcndr, dcndL, qeeq, dqdr, dqdL)
   call env%check(exitRun)
   if (exitRun) then
      call env%error("Could not get EN charges for D4 dispersion", source)
      return
   end if

   call getCoordinationNumber(mol, trans, 40.0_wp, cnType%cov, &
      & ccn, dccndr, dccndL)
   call latp%getLatticePoints(trans, 60.0_wp)
   call d4_gradient(mol, xtbData%dispersion%dispm, trans, xtbData%dispersion%dpar, &
      & xtbData%dispersion%g_a, xtbData%dispersion%g_c, xtbData%dispersion%wf, &
      & 60.0_wp, ccn, dccndr, dccndL, qeeq, dqdr, dqdL, ed, gradient, sigma)

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Evaluation of dispersion energy failed", source)
      return
   end if

   ! better save than sorry, delete the D4-EEQ charges
   qeeq(:) = 0.0_wp
   dqdr(:, :, :) = 0.0_wp
   dqdL(:, :, :) = 0.0_wp

   if (profile) call timer%measure(4)
! ---------------------------------------
!  Get EEQ charges q(1:n) + dqdr(3,1:n,1:n) under pbc
! ---------------------------------------

   if (allocated(gbsa)) then
      if (profile) call timer%measure(9,"GBSA setup")
      ! compute Born radii
      call gbsa%update(env, mol%at, mol%xyz)
      ! add SASA term to energy and gradient
      ees = gbsa%gsasa + gbsa%gshift
      gsolv = gbsa%gsasa
      gradient = gradient + gbsa%dsdr
      if (profile) call timer%measure(9)
   endif


   if (profile) call timer%measure(3,"EEQ model density")
   ! names DO NOT corresponds to content of variables, obviously...
   call gfn0_charge_model(chrgeq,mol%n,mol%at,xtbData%coulomb)
   ! initialize electrostatic energy
   if (allocated(gbsa)) then
      call eeq_chrgeq(mol,env,chrgeq,gbsa,cn,dcndr,qeeq,dqdr, &
         &            ees,gsolv,gradient,.false.,.true.,.true.)
   else
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
      allocate(chargeWidth(1, nid))
      do ii = 1, nid
         ati = idnum(ii)
         chargeWidth(1, ii) = xtbData%coulomb%chargeWidth(ati)
      end do
      call init(coulomb, env, mol, chargeWidth)
      call init(eeq, env, xtbData%coulomb%electronegativity, &
         & xtbData%coulomb%kcn, xtbData%coulomb%chemicalHardness, num=idnum)
      call eeq%chargeEquilibration(env, mol, coulomb, cn, dcndr, dcndL, &
         & ees, gradient, sigma, qat=qeeq, dqdr=dqdr, dqdL=dqdL)
   endif

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Electronegativity equilibration failed", source)
      return
   end if

   wfn%q = qeeq

   if (profile) call timer%measure(3)
   if (profile) call timer%measure(5,"Integral evaluation")
! ---------------------------------------
!  Build AO overlap S and H0 integrals under pbc
! ---------------------------------------
   call getSelfEnergy(xtbData%hamiltonian, xtbData%nShell, mol%at, cn, wfn%q, &
      & selfEnergy, dSEdcn, dSEdq)
   if (ccm) then
      call ccm_build_SH0(xtbData%nShell, xtbData%hamiltonian, selfEnergy, &
         & mol%n, mol%at, basis, nbf, nao, mol%xyz, mol%lattice, intcut, &
         & s, h0, wsc)
   else
      call latp%getLatticePoints(trans, sqrt(800.0_wp))
      call pbc_build_SH0(xtbData%nShell, xtbData%hamiltonian, selfEnergy, &
         & mol%n, mol%at, basis, nbf, nao, mol%xyz, trans, intcut, s, h0)
   endif

   if (profile) call timer%measure(5)
   if (profile) call timer%measure(6,"Cholesky factorization")

! ---------------------------------------
!  Check for near linear dependencies via Cholesky decomposition
! ---------------------------------------
   call cholesky(env%unit,pr,nao,S,orthog) ! S is not modified
   if(orthog)then
      !if (profile) call timer%measure(7,"Canonical orthogonalization")
      call renorm(nao,S) ! S is renormalized
      call canorthog(env%unit,nao,S,X,xdim,pr,fail)
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

   if (fail) then
      call env%error("Diagonalization of Hamiltonian failed", source)
      return
   end if

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
      call timer%write_timing(env%unit,7,"Diagonalization")
   if (profile) call timer%measure(8,"Gradient calculation")
! ======================================================================
!  GRADIENT (100% analytical)
! ======================================================================
   ! repulsion energy + gradient
   !gradient = 0.0_wp; sigma = 0.0_wp
   call latp%getLatticePoints(trans, 40.0_wp)
   call drep_grad(xtbData%repulsion,mol,trans,ep,gradient,sigma)
   ! short ranged bond energy + gradient
   if (allocated(xtbData%srb)) then
      call latp%getLatticePoints(trans, sqrt(200.0_wp))
      call dsrb_grad(mol,xtbData%srb,cn,dcndr,dcndL,trans,esrb,gradient,sigma)
   end if
   !etot = ep + esrb; return
   ! h0 gradient
   allocate( dHdcn(mol%n), dHdq(mol%n), pew(nao,nao), tmp(nao), &
      &      source = 0.0_wp )
   tmp = wfn%focc*wfn%emo*evtoau
   ! setup energy weighted density matrix = pew
   call dmat(nao,tmp,wfn%C,pew)
   if (ccm) then
      call ccm_build_dSH0(xtbData%nShell, xtbData%hamiltonian, selfEnergy, &
         & dSEdcn, dSEdq, mol%n, basis, intcut, nao, nbf, mol%at, mol%xyz, &
         & mol%lattice, wfn%P, Pew, gradient, sigma, dhdcn, dhdq, wsc)
   else
      call latp%getLatticePoints(trans, sqrt(800.0_wp))
      call pbc_build_dSH0(xtbData%nShell, xtbData%hamiltonian, selfEnergy, &
         & dSEdcn, dSEdq, mol%n, basis, intcut, nao, nbf, mol%at, mol%xyz, &
         & trans, wfn%P, Pew, gradient, sigma, dhdcn, dhdq)
   endif
   if (mol%npbc > 0) then
      ! setup CN sigma
      call dgemv('n',9,mol%n,1.0_wp,dcndL,9,dhdcn,1,1.0_wp,sigma,1)
      ! setup  q sigma
      call dgemv('n',9,mol%n,1.0_wp,dqdL, 9, dhdq,1,1.0_wp,sigma,1)
   endif
   ! setup CN gradient
   call dgemv('n',3*mol%n,mol%n, 1.0_wp,dcndr,3*mol%n,dhdcn,1,1.0_wp,gradient,1)
   ! setup  q gradient
   call dgemv('n',3*mol%n,mol%n, 1.0_wp,dqdr, 3*mol%n, dhdq,1,1.0_wp,gradient,1)

   if (profile) call timer%measure(8)

!  calculate the norm for printout
   res%gnorm = sqrt(sum( gradient**2 ))

! ---------------------------------------
!  Properties + printout
! ---------------------------------------
   printing: if (pr) then
      ! print orbital energies and occupation numbers
      !call preig(6,wfn%focc,1.0_wp,wfn%emo, &
      !     max(wfn%ihomoa-12,1),min(wfn%ihomoa+11,nao))
      if(.not.orthog)then
        call print_orbital_eigenvalues(env%unit,wfn,5)
        if ((wfn%ihomo+1.le.nao).and.(wfn%ihomo.ge.1)) &
           egap = wfn%emo(wfn%ihomo+1)-wfn%emo(wfn%ihomo)
      else
        wfn%nao = xdim
        call print_orbital_eigenvalues(env%unit,wfn,5)
        if ((wfn%ihomo+1.le.xdim).and.(wfn%ihomo.ge.1)) &
           egap = wfn%emo(wfn%ihomo+1)-wfn%emo(wfn%ihomo)
        wfn%nao = nao
      endif


   endif printing
   
   !--------------------------!
   ! Wiberg-Mayer bond orders !
   !--------------------------!

   ! closed-shell !
   if (wfn%nopen == 0) then
      
      call get_wiberg(mol%n,basis%nao,mol%at,mol%xyz,wfn%P,S,wfn%wbo,basis%fila2)
   
   ! (restricted) open-shell !
   else if (wfn%nopen > 0) then   
         
      allocate(Pa(basis%nao,basis%nao))
      allocate(Pb(basis%nao,basis%nao))
      
      ! obtain alpha and beta spin densities !
      call dmat(basis%nao, wfn%focca, wfn%C, Pa) 
      call dmat(basis%nao, wfn%foccb, wfn%C, Pb) 
      
      call get_unrestricted_wiberg(mol%n, basis%nao, mol%at, mol%xyz, Pa, Pb ,S, wfn%wbo, &
         & basis%fila2)

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
   if (allocated(gbsa)) then
      res%g_solv  = gsolv
      !res%g_born  = gborn    ! not returned
      res%g_sasa  = gbsa%gsasa
      !res%g_hb    = gbsa%ghb ! not returned
      res%g_shift = gbsa%gshift
   endif
   ! do NOT calculate the dipole moment from the density, because it's really bad
   res%dipole  = matmul(mol%xyz,wfn%q)
   res%g_solv  = 0.0_wp

   if (profile.and.pr) call timer%write(env%unit,'EHT')

end associate

   if (profile) call timer%deallocate

end subroutine peeq

! repulsion
pure subroutine drep_grad(repData,mol,trans,erep,gradient,sigma)
   use xtb_type_molecule
   use xtb_type_param
   use xtb_pbc_tools
   use xtb_pbc
   implicit none
   type(TRepulsionData), intent(in) :: repData
   type(TMolecule), intent(in) :: mol
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(out) :: erep
   integer :: i,j,k,lin
   integer :: iat,jat,ati,atj
   real(wp) :: ri(3),dr(3)
   real(wp) :: rij(3)
   real(wp) :: r2,r,r2top34,r2top34top2
   real(wp) :: den2,den4
   real(wp) :: alpha
   real(wp) :: repab
   real(wp) :: expterm
   real(wp) :: dtmp
   real(wp), parameter :: rthr = 1600.0_wp
   real(wp) :: w,t(3)
   integer  :: latrep(3),tx,ty,tz,itr
   call get_realspace_cutoff(mol%lattice,rthr,latrep)
   w = 1.0_wp
   ! initialize
   erep = 0.0_wp
   do i = 1, mol%n
      ati=mol%at(i)
      do j = 1, i
         do itr = 1, size(trans, dim=2)
            t = trans(:, itr)
            rij = mol%xyz(:,i) - mol%xyz(:,j) + t
            r2 = sum(rij**2)
            if(r2.gt.rthr .or. r2.lt.1.0e-6_wp) cycle
            atj=mol%at(j)
            r = sqrt(r2)
            den2 = (repData%electronegativity(ati) - repData%electronegativity(atj))**2
            den4 = den2**2
            alpha=sqrt(repData%alpha(ati)*repData%alpha(atj))&
               *(1.0_wp+(0.01_wp*den2+0.01_wp*den4)*repData%enScale)
            repab = repData%zeff(ati)*repData%zeff(atj)
            r2top34 = r2**0.75_wp
            r2top34top2 = r2top34**2
            expterm = exp(-alpha*r2top34)*repab
            ! save repulsion energy
            erep = erep + expterm/r * w
            ! save repulsion gradient
            dtmp = expterm*(1.5_wp*alpha*r2top34 + 1)/r2top34top2 * w
            gradient(:,i) = gradient(:,i) - dtmp*rij
            gradient(:,j) = gradient(:,j) + dtmp*rij
            sigma = sigma - dtmp*outer_prod_3x3(rij,rij)
         enddo ! k WSC partner
      enddo ! j atom
   enddo ! i atom

end subroutine drep_grad

! short-ranged bond correction
subroutine dsrb_grad(mol,srb,cn,dcndr,dcndL,trans,esrb,gradient,sigma)

   use xtb_type_param
   use xtb_type_molecule

   use xtb_approxrab
   use xtb_disp_ncoord
   use xtb_pbc_tools
   use xtb_pbc

   implicit none

   type(TMolecule), intent(in) :: mol
   type(TShortRangeData), intent(in) :: srb
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(out) :: esrb
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:,:,:)
   real(wp), intent(in) :: dcndL(:,:,:)

   integer :: i,j,k,wsAt
   integer :: iat,jat,ati,atj
   integer :: nsrb,itr
   real(wp) :: den
   real(wp) :: expterm
   ! distances
   real(wp) :: rij(3)
   real(wp) :: r2,r,dr,rab
   real(wp) :: dtmp,pre
   real(wp) :: w
   ! allocatables
   real(wp), allocatable :: dEdr0(:)
   real(wp), allocatable :: drab0dr(:,:,:)
   real(wp), allocatable :: drab0dL(:,:,:)
   real(wp), allocatable :: rab0(:)
   integer,  allocatable :: srblist(:,:)

   ! initialize
   esrb = 0.0_wp

   w = 1.0_wp

   call build_srblist(mol,nsrb,srblist)
   if (nsrb.eq.0) return

   ! get memory
   allocate( dEdr0(nsrb), source = 0.0_wp )
   allocate( drab0dr(3,mol%n,nsrb), source = 0.0_wp )
   allocate( drab0dL(3,3,nsrb),     source = 0.0_wp )
   allocate( rab0(nsrb),            source = 0.0_wp )
   ! get approximated distances rab and gradients
   call pbc_approx_rab(mol%n,mol%at,mol%xyz,cn,dcndr,dcndL,nsrb,srblist,srb%shift, &
      &                rab0,drab0dr,drab0dL)
   do i = 1, nsrb
      iat = srblist(1,i)
      jat = srblist(2,i)
      ati = mol%at(iat)
      atj = mol%at(jat)
      den = en(ati) - en(atj)
      pre = srb%steepness*(1.0_wp + srb%enScale*den**2)
      do itr = 1, size(trans, dim=2)
         rij = mol%xyz(:,iat) - (mol%xyz(:,jat) + trans(:, itr))
         rab = norm2(rij)
         dr = rab - rab0(i)
         expterm = srb%prefactor*exp(-pre*dr**2)
         ! save SRB energy
         esrb = esrb + expterm * w
         dtmp = 2.0_wp*pre*dr*expterm * w
         gradient(:,iat) = gradient(:,iat) - dtmp*rij/rab
         gradient(:,jat) = gradient(:,jat) + dtmp*rij/rab
         ! three body gradient
         dEdr0(i) = dEdr0(i) + dtmp
         sigma = sigma - dtmp*spread(rij, 1, 3)*spread(rij, 2, 3)/rab
      enddo ! rep
   enddo ! i

   call contract(drab0dr, dEdr0, gradient, beta=1.0_wp)
   call contract(drab0dL, dEdr0, sigma, beta=1.0_wp)

end subroutine dsrb_grad

pure subroutine build_srblist(mol,nsrb,srblist)
   use xtb_type_molecule
   implicit none
   type(TMolecule), intent(in) :: mol
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
subroutine ccm_build_SH0(nShell, hData, selfEnergy, nat, at, basis, nbf, nao, &
      & xyz, lattice, intcut, sint, h0, wsc)

   use xtb_type_basisset
   use xtb_type_wsc

   use xtb_lin, only : lin
   use xtb_intgrad
   use xtb_scc_core

   implicit none

   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
   type(tb_wsc), intent(in) :: wsc
   type(TBasisset), intent(in) :: basis
   real(wp), intent(in) :: selfEnergy(:, :)

   integer, intent(in)  :: nat
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nao
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: lattice(3,3)
   real(wp),intent(in)  :: intcut

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

   !$omp parallel do default(none) schedule(dynamic) &
   !$omp private(iat,jat,ij,ati,cc,ci,rab2,atj,ish,ishtyp,valaoi,valaoj, &
   !$omp&        ri,rj,icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp, &
   !$omp&        ss,saw,est,alpi,alpj,ab,iprim,jprim,ip,jp,km,shpoly, &
   !$omp&        mli,mlj,tmp,zi,zj,zetaij,enpoly,iao,jao, &
   !$omp&        ii,jj,k,den,den2,den4,i,j,il,jl,hii,hjj,hav,t) &
   !$omp shared(sint,h0) &
   !$omp shared(wsc,basis,at,nShell,hData,xyz,lattice,intcut,nat,selfEnergy)
   do iat = 1, nat
      ri  = xyz(:,iat)
      ati = at(iat)
      do jat = 1, iat-1
         atj = at(jat)

         den = hData%electronegativity(ati) - hData%electronegativity(atj)
         den2 =  den**2
         den4 = den2**2

         ishells: do ish = 1, nShell(ati)
            ishtyp = hData%angShell(ish,ati)
            icao = basis%caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = nShell(atj)
            if (iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jshtyp = hData%angShell(jsh,atj)
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               il = ishtyp+1
               jl = jshtyp+1
               ! diagonals are the same for all H0 elements
               hii = selfEnergy(ish, iat)
               hjj = selfEnergy(jsh, jat)

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + hData%enScale(jl-1,il-1)*den2 &
                  & + hData%enScale4*hData%enScale(jl-1,il-1)*den4)

               ! we scale the two shells depending on their exponent
               zi = hData%slaterExponent(ish, ati)
               zj = hData%slaterExponent(jsh, atj)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = hData%kScale(jl-1,il-1) * hData%pairParam(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = hData%valenceShell(ish, ati).eq.0
               valaoj = hData%valenceShell(jsh, atj).eq.0
               ! and scale appropiately
               if (valaoi) then
                  if (valaoj) then
                     km = 0.0_wp
                  else
                     km = km * hData%kDiff
                  endif
               else
                  if (valaoj) km = km * hData%kDiff
               endif

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * km * (hii + hjj)

               wscAt: do kat = 1, wsc%itbl(jat,iat)
                  ss = 0.0_wp
                  t = wsc%lattr(:,kat,jat,iat)
                  rj = xyz(:,jat) + matmul(lattice,t)

                  ! distance dependent polynomial
                  shpoly = shellPoly(hData%shellPoly(il,ati),hData%shellPoly(jl,atj),&
                     &               hData%atomicRad(ati),hData%atomicRad(atj),ri,rj)

                  ! get overlap integral
                  call get_overlap(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj,point, &
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
                        !sint(iao,jao) = sint(iao,jao) + ss(jj,ii)*wsc%w(jat,iat)
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
   !$omp parallel do default(none) shared(nao, sint) private(iao, jao)
   do iao = 1, nao
      do jao = 1, iao - 1
         sint(iao, jao) = sint(jao, iao)
      end do
   end do

   ! diagonal elements
   do iat = 1, nat
      ati = at(iat)
      do ish = 1, nShell(ati)
         iao = 1+basis%saoshell(ish,iat)
         ishtyp = hData%angShell(ish,ati)
         il = ishtyp + 1
         do iao = 1, llao2(ishtyp)
            i = iao+basis%saoshell(ish,iat)
            sint(i,i)=1.0_wp+sint(i,i)

            ii  = i*(1+i)/2 ! H0 is packed, note i*(i-1)/2+i = i*(1+i)/2

            ! calculate environment dependent shift
            hii = selfEnergy(ish, iat)
            H0(ii) = hii
         end do
      end do
   end do

end subroutine ccm_build_SH0

! ------------------------------------------------------------------------
!  Calculate the periodic AO overlap matrix
! ------------------------------------------------------------------------
subroutine pbc_build_SH0(nShell, hData, selfEnergy, nat, at, basis, nbf, nao, &
      & xyz, trans, intcut, sint, h0)

   use xtb_type_basisset

   use xtb_lin, only : lin
   use xtb_intgrad
   use xtb_scc_core

   implicit none

   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
   type(TBasisset), intent(in) :: basis
   real(wp), intent(in) :: selfEnergy(:, :)

   integer, intent(in)  :: nat
   integer, intent(in)  :: nbf
   integer, intent(in)  :: nao
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: trans(:,:)
   real(wp),intent(in)  :: intcut

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
   integer  :: itr
   real(wp) :: t(3),w

   w = 1.0_wp !/ real(product(2*latrep+1),wp)

   sint = 0.0_wp
   h0   = 0.0_wp

   !$omp parallel do default(none) schedule(dynamic) &
   !$omp private(iat,jat,ij,ati,cc,ci,rab2,atj,ish,ishtyp,valaoi,valaoj, &
   !$omp&        ri,rj,icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp, &
   !$omp&        ss,saw,est,alpi,alpj,ab,iprim,jprim,ip,jp,km,shpoly, &
   !$omp&        mli,mlj,tmp,zi,zj,zetaij,enpoly,iao,jao, &
   !$omp&        ii,jj,k,den,den2,den4,i,j,il,jl,hii,hjj,hav,itr,t) &
   !$omp shared(sint,h0) &
   !$omp shared(basis,at,nShell,hData,xyz,intcut,nat,trans,w,selfEnergy)
   do iat = 1, nat
      ri  = xyz(:,iat)
      ati = at(iat)
      do jat = 1, iat-1
         atj = at(jat)

         den = hData%electronegativity(ati) - hData%electronegativity(atj)
         den2 =  den**2
         den4 = den2**2

         ishells: do ish = 1, nShell(ati)
            ishtyp = hData%angShell(ish,ati)
            icao = basis%caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = nShell(atj)
            if (iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jshtyp = hData%angShell(jsh,atj)
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               il = ishtyp+1
               jl = jshtyp+1
               ! diagonals are the same for all H0 elements
               hii = selfEnergy(ish, iat)
               hjj = selfEnergy(jsh, jat)

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + hData%enScale(jl-1,il-1)*den2 &
                  & + hData%enScale4*hData%enScale(jl-1,il-1)*den4)

               ! we scale the two shells depending on their exponent
               zi = hData%slaterExponent(ish, ati)
               zj = hData%slaterExponent(jsh, atj)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = hData%kScale(jl-1,il-1) * hData%pairParam(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = hData%valenceShell(ish, ati).eq.0
               valaoj = hData%valenceShell(jsh, atj).eq.0
               ! and scale appropiately
               if (valaoi) then
                  if (valaoj) then
                     km = 0.0_wp
                  else
                     km = km * hData%kDiff
                  endif
               else
                  if (valaoj) km = km * hData%kDiff
               endif

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * km * (hii + hjj)

               do itr = 1, size(trans, dim=2)
                  t = trans(:, itr)
                  rj = xyz(:,jat) + t
                  ss = 0.0_wp

                  ! distance dependent polynomial
                  shpoly = shellPoly(hData%shellPoly(il,ati),hData%shellPoly(jl,atj),&
                     &               hData%atomicRad(ati),hData%atomicRad(atj),ri,rj)

                  ! get overlap integral
                  call get_overlap(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj,point, &
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
   !$omp parallel do default(none) shared(nao, sint) private(iao, jao)
   do iao = 1, nao
      do jao = 1, iao - 1
         sint(iao, jao) = sint(jao, iao)
      end do
   end do

   ! diagonal elements
   do iat = 1, nat
      ati = at(iat)
      do ish = 1, nShell(ati)
         iao = 1+basis%saoshell(ish,iat)
         ishtyp = hData%angShell(ish,ati)
         il = ishtyp + 1
         do iao = 1, llao2(ishtyp)
            i = iao+basis%saoshell(ish,iat)
            sint(i,i)=1.0_wp+sint(i,i)

            ii  = i*(1+i)/2 ! H0 is packed, note i*(i-1)/2+i = i*(1+i)/2

            ! calculate environment dependent shift
            hii = selfEnergy(ish, iat)
            H0(ii) = hii
         end do
      end do
   end do

end subroutine pbc_build_SH0

! ------------------------------------------------------------------------
!  Calculate the gradient resulting from a periodic AO overlap matrix
! ------------------------------------------------------------------------
subroutine ccm_build_dSH0(nShell, hData, selfEnergy, dSEdcn, dSEdq, nat, basis, &
      & thr, nao, nbf, at, xyz, lattice, P, Pew, g, sigma, dHdcn, dHdq, wsc)
   use xtb_mctc_constants, only : pi
   use xtb_mctc_convert

   use xtb_type_wsc
   use xtb_type_basisset

   use xtb_intgrad
   use xtb_grad_core

   implicit none

   ! intent in
   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
   real(wp), intent(in) :: selfEnergy(:, :)
   real(wp), intent(in) :: dSEdcn(:, :)
   real(wp), intent(in) :: dSEdq(:, :)
   integer, intent(in)      :: nat
   type(tb_wsc), intent(in) :: wsc
   type(TBasisset), intent(in) :: basis
   real(wp),intent(in)      :: thr
   integer, intent(in)      :: nao
   integer, intent(in)      :: nbf
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   real(wp),intent(in) :: lattice(3,3)
   real(wp),intent(in) :: P(nao,nao)
   real(wp),intent(in) :: Pew(nao,nao)
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
   !$omp shared(wsc,basis,at,nShell,hData,xyz,lattice,thr,nat,selfEnergy, &
   !$omp        dSEdcn,dSEdq,P,Pew)
   !$omp do schedule(runtime)
   do iat = 1, nat
      ri  = xyz(:,iat)
      ati = at(iat)
      do jat = 1, iat-1
         atj = at(jat)

         den = hData%electronegativity(ati) - hData%electronegativity(atj)
         den2 =  den**2
         den4 = den2**2

         ishells: do ish = 1, nShell(ati)
            ishtyp = hData%angShell(ish,ati)
            icao = basis%caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = nShell(atj)
            if(iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jshtyp = hData%angShell(jsh,atj)
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               il = ishtyp+1
               jl = jshtyp+1
               ! diagonals are the same for all H0 elements
               hii = selfEnergy(ish, iat)
               hjj = selfEnergy(jsh, jat)

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + hData%enScale(jl-1,il-1)*den2 &
                  & + hData%enScale4*hData%enScale(jl-1,il-1)*den4)

               ! we scale the two shells depending on their exponent
               zi = hData%slaterExponent(ish, ati)
               zj = hData%slaterExponent(jsh, atj)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = hData%kScale(jl-1,il-1) * hData%pairParam(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = hData%valenceShell(ish, ati).eq.0
               valaoj = hData%valenceShell(jsh, atj).eq.0
               ! and scale appropiately
               if (valaoi) then
                  if (valaoj) then
                     km = 0.0_wp
                  else
                     km = km * hData%kDiff
                  endif
               else
                  if (valaoj) km = km * hData%kDiff
               endif

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * km * (hii + hjj)

               wscAt: do kat = 1,wsc%itbl(jat,iat)
                  t = wsc%lattr(:,kat,jat,iat)
                  rj = xyz(:,jat) + matmul(lattice,t)
                  rij = ri - rj
                  rij2 = sum(rij**2)

                  ! distance dependent polynomial
                  call dshellPoly(hData%shellPoly(il,ati),hData%shellPoly(jl,atj),&
                     & hData%atomicRad(ati),hData%atomicRad(atj),rij2,ri,rj,&
                     & shpoly,dshpoly)

                  call get_grad_overlap(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj, &
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
                        dhdcn(iat) = dhdcn(iat) + HPij*dSEdcn(ish, iat)
                        ! save dE/dCN for CNj
                        dhdcn(jat) = dhdcn(jat) + HPij*dSEdcn(jsh, jat)

                        ! save dE/dq for qi
                        dhdq(iat) = dhdq(iat) + HPij*dSEdq(ish, iat)
                        ! save dE/dq for qj
                        dhdq(jat) = dhdq(jat) + HPij*dSEdq(jsh, jat)

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
   do iat = 1, nat
      ati = at(iat)
      do ish = 1, nShell(ati)
         iao = 1+basis%saoshell(ish,iat)
         ishtyp = hData%angShell(ish,ati)
         il = ishtyp + 1
         do iao = 1, llao2(ishtyp)
            i = iao+basis%saoshell(ish,iat)
            ii = i*(1+i)/2 ! H0 is packed, note i*(i-1)/2+i = i*(1+i)/2

            Pii = P(i,i)*evtoau

            ! save dE/dCN for CNi
            dhdcn(iat) = dhdcn(iat) + Pii*dSEdcn(ish, iat)
            ! save dE/dq for qi
            dhdq(iat) = dhdq(iat) + Pii*dSEdq(ish, iat)
         end do
      end do
   end do

end subroutine ccm_build_dSH0

! ------------------------------------------------------------------------
!  Calculate the gradient resulting from a periodic AO overlap matrix
! ------------------------------------------------------------------------
subroutine pbc_build_dSH0(nShell, hData, selfEnergy, dSEdcn, dSEdq, nat, basis, &
      & thr, nao, nbf, at, xyz, trans, P, Pew, g, sigma, dHdcn, dHdq)
   use xtb_mctc_constants, only : pi
   use xtb_mctc_convert

   use xtb_type_wsc
   use xtb_type_basisset

   use xtb_intgrad
   use xtb_grad_core

   implicit none

   ! intent in
   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
   real(wp), intent(in) :: selfEnergy(:, :)
   real(wp), intent(in) :: dSEdcn(:, :)
   real(wp), intent(in) :: dSEdq(:, :)
   integer, intent(in)      :: nat
   type(TBasisset), intent(in) :: basis
   real(wp),intent(in)      :: thr
   integer, intent(in)      :: nao
   integer, intent(in)      :: nbf
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   real(wp),intent(in) :: trans(:,:)
   real(wp),intent(in) :: P(nao,nao)
   real(wp),intent(in) :: Pew(nao,nao)
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

   integer  :: il,jl,itr
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
   !$omp&        itr,t) &
   !$omp reduction (+:g,sigma,dhdcn,dhdq) &
   !$omp shared(basis,at,nShell,hData,xyz,thr,nat,trans,P,Pew,w,selfEnergy, &
   !$omp& dSEdcn,dSEdq)
   !$omp do schedule(runtime)
   do iat = 1, nat
      ri  = xyz(:,iat)
      ati = at(iat)
      do jat = 1, iat-1
         atj = at(jat)

         den = hData%electronegativity(ati) - hData%electronegativity(atj)
         den2 =  den**2
         den4 = den2**2

         ishells: do ish = 1, nShell(ati)
            ishtyp = hData%angShell(ish,ati)
            icao = basis%caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = nShell(atj)
            if(iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jshtyp = hData%angShell(jsh,atj)
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               il = ishtyp+1
               jl = jshtyp+1
               ! diagonals are the same for all H0 elements
               hii = selfEnergy(ish, iat)
               hjj = selfEnergy(jsh, jat)

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + hData%enScale(jl-1,il-1)*den2 &
                  & + hData%enScale4*hData%enScale(jl-1,il-1)*den4)

               ! we scale the two shells depending on their exponent
               zi = hData%slaterExponent(ish, ati)
               zj = hData%slaterExponent(jsh, atj)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = hData%kScale(jl-1,il-1) * hData%pairParam(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = hData%valenceShell(ish, ati).eq.0
               valaoj = hData%valenceShell(jsh, atj).eq.0
               ! and scale appropiately
               if (valaoi) then
                  if (valaoj) then
                     km = 0.0_wp
                  else
                     km = km * hData%kDiff
                  endif
               else
                  if (valaoj) km = km * hData%kDiff
               endif

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * km * (hii + hjj)

               do itr = 1, size(trans, dim=2)
                  t = trans(:, itr)
                  rj = xyz(:,jat) + t
                  rij = ri - rj
                  rij2 = sum(rij**2)

                  ! distance dependent polynomial
                  call dshellPoly(hData%shellPoly(il,ati),hData%shellPoly(jl,atj),&
                     & hData%atomicRad(ati),hData%atomicRad(atj),rij2,ri,rj,&
                     & shpoly,dshpoly)

                  call get_grad_overlap(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj, &
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
                        dhdcn(iat) = dhdcn(iat) + HPij*dSEdcn(ish, iat)
                        ! save dE/dCN for CNj
                        dhdcn(jat) = dhdcn(jat) + HPij*dSEdcn(jsh, jat)

                        ! save dE/dq for qi
                        dhdq(iat) = dhdq(iat) + HPij*dSEdq(ish, iat)
                        ! save dE/dq for qj
                        dhdq(jat) = dhdq(jat) + HPij*dSEdq(jsh, jat)

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
   do iat = 1, nat
      ati = at(iat)
      do ish = 1, nShell(ati)
         iao = 1+basis%saoshell(ish,iat)
         ishtyp = hData%angShell(ish,ati)
         il = ishtyp + 1
         do iao = 1, llao2(ishtyp)
            i = iao+basis%saoshell(ish,iat)
            ii = i*(1+i)/2 ! H0 is packed, note i*(i-1)/2+i = i*(1+i)/2

            Pii = P(i,i)*evtoau

            ! save dE/dCN for CNi
            dhdcn(iat) = dhdcn(iat) + Pii*dSEdcn(ish, iat)
            ! save dE/dq for qi
            dhdq(iat) = dhdq(iat) + Pii*dSEdq(ish, iat)
         end do
      end do
   end do

end subroutine pbc_build_dSH0


end module xtb_peeq
