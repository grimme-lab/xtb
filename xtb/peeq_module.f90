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

module peeq_module
   use iso_fortran_env, wp => real64
   use mctc_econv
   use mctc_la
   implicit none

   !> profiling
   logical,parameter,private :: profile = .true.

   !> Numerical thresholds for GFN0-xTB
   type :: eeq_thresholds
      !> primitive cut-off
      real(wp) :: intcut
      !> integral neglect threshold
      real(wp) :: neglect
      !> coordination number cutoff
      real(wp) :: cutoff_cn
      !> integral distance cutoff
      real(wp) :: cutoff_ham
      !> repulsion polynomial cutoff
      real(wp) :: cutoff_rep
      !> short range basis cutoff
      real(wp) :: cutoff_srb
      !> long-ranged dispersion cutoff
      real(wp) :: cutoff_disp
   end type eeq_thresholds

   interface eeq_thresholds
      module procedure :: new_thresholds
   end interface

   interface maxval
      module procedure :: eeqthr_maxval
   end interface

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
   use tbdef_neighbourlist

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
   type(tb_neighbourlist) :: neighlist
   type(tb_neighlist_generator) :: neighgen
   integer, allocatable :: neighs(:)

   type(chrg_parameter) :: chrgeq
   real(wp),allocatable,dimension(:)     :: cn
   real(wp),allocatable,dimension(:,:,:) :: dcndr
   real(wp),allocatable,dimension(:,:,:) :: dcndL
   real(wp),allocatable,dimension(:,:)   :: X
   real(wp),allocatable,dimension(:,:)   :: S
   real(wp),allocatable,dimension(:,:)   :: H
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
   type(eeq_thresholds) :: thr
   ! energies
   real(wp) :: eel,ed,ees,exb,ep

   ! PEEQ information
   real(wp) :: ken(4, 4)
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

   if (profile) then
      if (lgbsa) then
         call timer%new(10,.false.)
      else
         call timer%new(9,.false.)
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

   ! Numerical stuff and cutoffs
   thr = eeq_thresholds(acc)

! ---------------------------------------
!  IMPORTANT FACT: H is given in eV
! ---------------------------------------

! ---------------------------------------
!  Get memory
! ---------------------------------------
   allocate(cn(mol%n), source=0.0_wp)
   allocate(dcndr(3,mol%n,mol%n), source=0.0_wp)
   allocate(dcndL(3,3,mol%n), source=0.0_wp)
   allocate(H(basis%nao,basis%nao), source=0.0_wp)
   allocate(X(basis%nao,basis%nao), source=0.0_wp)
   allocate(S(basis%nao,basis%nao), source=0.0_wp)
   allocate(kcnao(basis%nao), source=0.0_wp)
   allocate(kqao(basis%nao), source=0.0_wp)
   allocate(zsh(basis%nshell), source=0.0_wp)
   allocate(qeeq(mol%n), source=0.0_wp)
   allocate(dqdr(3,mol%n,mol%n+1), source=0.0_wp)
   allocate(dqdL(3,3,mol%n+1), source=0.0_wp)
   allocate(neighs(len(mol)), source=0)

! ---------------------------------------
!  Fill levels
! ---------------------------------------
   call setzshell(mol%n,mol%at,basis%nshell,mol%z,zsh,eatoms,0)
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

! ---------------------------------------
!  Calculate distances
! ---------------------------------------
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
      write(iunit,scifmt) "-> integral cutoff ",thr%intcut,  "    "
      write(iunit,scifmt) "-> integral neglect",thr%neglect, "    "
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
   if (profile) call timer%measure(9,"Neighbourlist")

   call neighgen%new(maxval(thr), len(mol), mol%xyz, mol%lattice, mol%npbc > 0)
   call neighlist%new(neighgen)

   if (profile) call timer%measure(9)
   if (profile) call timer%measure(2,"Coordination number")
! ---------------------------------------
!  Get CN(1:n) + dcndr(3,1:n,1:n) under pbc
! ---------------------------------------
   call neighlist%get_neighs(neighs, cutoff=thr%cutoff_cn)
   call get_coordination_number(mol, neighs, neighlist, tb_cn_type%erf, &
      &                         cn, dcndr, dcndL)
   call dncoord_logcn(len(mol), cn, dcndr, dcndL, cn_max=8.0_wp)

   if (profile) call timer%measure(2)
! ---------------------------------------
!  Get EEQ charges q(1:n) + dqdr(3,1:n,1:n) under pbc
! ---------------------------------------

   if (lgbsa) then
      if (profile) call timer%measure(10,"GBSA setup")
      call new_gbsa(gbsa,mol%n,mol%at)
      call update_nnlist_gbsa(gbsa,mol%xyz,.false.)
      ! compute Born radii
      call compute_brad_sasa(gbsa,mol%xyz)
      ! add SASA term to energy and gradient
      ees = gbsa%gsasa
      gsolv = gbsa%gsasa
      g = g + gbsa%dsdr
      if (profile) call timer%measure(10)
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
   call ddisp_peeq(mol,neighlist,thr,param,cn,dcndr,dcndL,grd,ed,g,sigma_tmp(:,:,2))

   if (profile) call timer%measure(4)
   if (profile) call timer%measure(5,"Integral evaluation")
! ---------------------------------------
!  Build AO overlap S and H0 integrals under pbc
! ---------------------------------------
   call neighlist%get_neighs(neighs, cutoff=thr%cutoff_ham)
   if (mol%npbc > 0 .and. ccm) then
      call ccm_build_SH0(mol%n,mol%at,basis,mol%xyz,mol%lattice, &
         &               wfn%q,cn,thr%intcut, &
         &               param%kmagic,ken,param%alphaj,param%kcnsh, &
         &               param%xbdamp,s,h,mol%wsc)
   else
      call build_SH0(mol,neighs,neighlist,basis,wfn%q,cn,thr%intcut, &
         &           param%kmagic,ken,param%alphaj,param%kcnsh, &
         &           param%xbdamp,s,h)
   endif

   if (profile) call timer%measure(5)
   if (profile) call timer%measure(6,"Cholesky factorization")

! ---------------------------------------
!  Check for near linear dependencies via Cholesky decomposition
! ---------------------------------------
   call cholesky(iunit,pr,basis%nao,S,orthog) ! S is not modified
   if(orthog)then
      !if (profile) call timer%measure(7,"Canonical orthogonalization")
      call renorm(basis%nao,S) ! S is renormalized
      call canorthog(iunit,basis%nao,S,X,xdim,pr,fail)
      !if (profile) call timer%measure(7)
   endif

   if (profile) call timer%measure(6)
   if (profile) call timer%measure(7,"Zeroth order Hamiltonian")

   if(.not.orthog)then
      call solve(.true.,basis%nao,wfn%ihomo,H,S,X,wfn%P,wfn%emo,fail)
   else
      call orthgsolve2(.true.,basis%nao,xdim,wfn%ihomo,H,S,X,wfn%P,wfn%emo,fail)
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
      if (wfn%ihomoa+1.le.basis%nao) then
         call fermismear(.false.,basis%nao,wfn%ihomoa,et,wfn%emo,wfn%focca,nfoda,efa,ga)
      endif
      if (wfn%ihomob+1.le.basis%nao) then
         call fermismear(.false.,basis%nao,wfn%ihomob,et,wfn%emo,wfn%foccb,nfodb,efb,gb)
      endif
      wfn%focc = wfn%focca + wfn%foccb
   endif
   ! create density matrix = save in wfn%P
   call dmat(basis%nao,wfn%focc,wfn%C,wfn%P)

   call mpopsh(mol%n,basis%nao,basis%nshell,basis%ao2sh,S,wfn%P,wfn%qsh)

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
   call neighlist%get_neighs(neighs, cutoff=thr%cutoff_rep)
   call drep_grad(mol, neighs, neighlist, param, ep, g, sigma_tmp(:,:,3))
   ! short ranged bond energy + gradient
   call neighlist%get_neighs(neighs, cutoff=thr%cutoff_srb)
   call dsrb_grad(mol, neighs, neighlist, param, cn, dcndr, dcndL, &
      &           esrb, g, sigma_tmp(:,:,4))
   !etot = ep + esrb; return
   ! h0 gradient
   allocate( dHdcn(mol%n), dHdq(mol%n), pew(basis%nao,basis%nao), tmp(basis%nao), &
      &      source = 0.0_wp )
   tmp = wfn%focc*wfn%emo*evtoau
   ! setup energy weighted density matrix = pew
   call dmat(basis%nao,tmp,wfn%C,pew)
   call neighlist%get_neighs(neighs, cutoff=thr%cutoff_ham)
   if (mol%npbc > 0 .and. ccm) then
      call ccm_build_dSH0(mol%n,basis,thr%intcut,mol%at,mol%xyz, &
         &                mol%lattice,wfn%q,cn, &
         &                wfn%P,Pew,g,sigma,dhdcn,dhdq,param%kmagic,ken, &
         &                param%alphaj,param%kcnsh,param%xbdamp,mol%wsc)
   else
      call build_dSH0(mol,basis,neighs,neighlist,thr%intcut,wfn%q,cn, &
         &            wfn%P,Pew,g,sigma,dhdcn,dhdq,param%kmagic,ken, &
         &            param%alphaj,param%kcnsh,param%xbdamp)
   endif
   ! setup CN gradient
   call dgemv('n',3*mol%n,mol%n,-1.0_wp,dcndr,3*mol%n,dhdcn,1,1.0_wp,g,1)
   ! setup  q gradient
   call dgemv('n',3*mol%n,mol%n, 1.0_wp,dqdr, 3*mol%n, dhdq,1,1.0_wp,g,1)
   if (mol%npbc > 0) then
      ! setup CN sigma
      call dgemv('n',9,mol%n,1.0_wp,dcndL,9,dhdcn,1,1.0_wp,sigma,1)
      ! setup  q sigma
      call dgemv('n',9,mol%n,1.0_wp,dqdL, 9, dhdq,1,1.0_wp,sigma,1)
   endif

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
        if ((wfn%ihomo+1.le.basis%nao).and.(wfn%ihomo.ge.1)) &
           egap = wfn%emo(wfn%ihomo+1)-wfn%emo(wfn%ihomo)
      else
        wfn%nao = xdim
        call print_orbital_eigenvalues(iunit,wfn,5)
        if ((wfn%ihomo+1.le.xdim).and.(wfn%ihomo.ge.1)) &
           egap = wfn%emo(wfn%ihomo+1)-wfn%emo(wfn%ihomo)
        wfn%nao = basis%nao
      endif

   endif printing

! ------------------------------------------------------------------------
!  get Wiberg bond orders
   call get_wiberg(mol%n,basis%nao,mol%at,mol%xyz,wfn%P,S,wfn%wbo,basis%fila2)

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
   ! TODO: maybe use density matrix and overlap integrals instead
   res%dipole  = matmul(mol%xyz,wfn%q)
   res%g_solv  = 0.0_wp

   if (profile.and.pr) call timer%write(iunit,'EHT')

! ---------------------------------------
! Free memory
! ---------------------------------------
   if (profile) call timer%deallocate
   if(allocated(S))         deallocate(S)
   if(allocated(H))         deallocate(H)
   if(allocated(X))         deallocate(X)
   if(allocated(kcnao))     deallocate(kcnao)
   if(allocated(zsh))       deallocate(zsh)
   if(allocated(dqdr))      deallocate(dqdr)

end subroutine peeq

! -----------------------------------------------------------------------
!  Calculate D4 dispersion gradient
! -----------------------------------------------------------------------
subroutine ddisp_peeq(mol,neighlist,thr,param,cn,dcndr,dcndL,grd,ed,gd,sigma)
   use iso_fortran_env, only : wp => real64
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data
   use tbdef_neighbourlist
   use tbmod_dftd4
   use eeq_model
   use ncoord

   type(tb_molecule), intent(in) :: mol
   type(scc_parameter), intent(in) :: param
   type(tb_neighbourlist), intent(in) :: neighlist
   type(eeq_thresholds), intent(in) :: thr
   type(chrg_parameter) :: chrgeq
   ! EEQ partial charges and derivatives
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)

   ! dispersion energy and derivative
   real(wp), intent(inout) :: ed
   real(wp), intent(inout) :: gd(:, :)
   real(wp), intent(inout) :: sigma(:, :)

   integer, allocatable :: neighs(:)
   real(wp), allocatable :: sdum(:, :)
   real(wp), allocatable :: gdum(:, :)
   real(wp), allocatable :: q(:)
   real(wp), allocatable :: dqdr(:, :, :)
   real(wp), allocatable :: dqdL(:, :, :)
   real(wp), allocatable :: covcn(:)
   real(wp), allocatable :: dcovcndr(:, :, :)
   real(wp), allocatable :: dcovcndL(:, :, :)
   real(wp), allocatable :: c6abns(:, :)
   real(wp), allocatable :: gw(:)
   real(wp), allocatable :: gdummy(:, :)

! -----------------------------------------------------------------------
!  Variables
! -----------------------------------------------------------------------
   integer :: i,j,k,l,m
   integer :: ndim
   integer :: mbd
   logical, intent(in) :: grd
   integer, dimension(3) :: rep_vdw(3), rep_cn(3)
 ! real space cutoffs
   real(wp), parameter :: cn_thr = 1600.0_wp
   real(wp), parameter :: crit_vdw = 4000.0_wp

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

! -----------------------------------------------------------------------
!  Get memory
! -----------------------------------------------------------------------
  allocate(neighs(len(mol)), source=0)
  allocate(c6abns(ndim, ndim), source=0.0_wp)
  allocate(gw(ndim), source=0.0_wp)
  allocate(q(mol%n), source=0.0_wp)
  allocate(dqdr(3, mol%n, mol%n+1), source=0.0_wp)
  allocate(dqdL(3, 3, mol%n+1), source=0.0_wp)
  allocate(covcn(mol%n), source=0.0_wp)
  allocate(dcovcndr(3, mol%n, mol%n), source=0.0_wp)
  allocate(dcovcndL(3, 3, mol%n), source=0.0_wp)

  ! get D4(EEQ) charges
  call new_charge_model_2019(chrgeq,mol%n,mol%at)
  ! neither sdum nor gdum need to be dummy allocated, since lgrad = .false. is set
  ! eeq_chrgeq must not reference sdum/gdum and we can avoid the dummy allocate
  call eeq_chrgeq(mol,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,edum,gdum,sdum, &
     &            .false.,.false.,.true.)

  call neighlist%get_neighs(neighs, cutoff=thr%cutoff_cn)
  call get_coordination_number(mol, neighs, neighlist, tb_cn_type%cov, &
     &                         covcn, dcovcndr, dcovcndL)
  dcovcndr = -dcovcndr ! FIXME

  call neighlist%get_neighs(neighs, cutoff=thr%cutoff_disp)
  call d4_gradient(mol, neighs, neighs, neighlist, param%disp, param%g_a, &
     &             param%g_c, param%wf, covcn, dcovcndr, dcovcndL, q, dqdr, dqdL, &
     &             ed, gd, sigma)

end subroutine ddisp_peeq

! repulsion
subroutine drep_grad(mol, neighs, neighlist, param, energy, gradient, sigma)
   use aoparam, only : rep, en
   use tbdef_molecule
   use tbdef_neighbourlist
   use tbdef_param
   implicit none
   ! intent in
   type(scc_parameter), intent(in) :: param
   type(tb_molecule), intent(in) :: mol
   integer, intent(in) :: neighs(:)
   type(tb_neighbourlist), intent(in) :: neighlist
   ! intent inout
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   ! intent out
   real(wp), intent(out) :: energy
   ! local variables
   integer :: iat, jat, ati, atj, ij, img
   real(wp) :: r2, r1, rij(3), r2top34,r2top34top2
   real(wp) :: den2, den4, alpha, repab, expterm
   real(wp) :: dE, dG(3), dS(3, 3)
   real(wp), allocatable :: energies(:)
   ! initialize
   energy = 0.0_wp
   allocate(energies(len(mol)), source=0.0_wp)
   !$omp parallel do default(none) &
   !$omp reduction(+:energies, gradient, sigma) &
   !$omp shared(mol, neighlist, neighs, param, en, rep) &
   !$omp private(r2, rij, r1, den2, den4, alpha, repab, r2top34, r2top34top2, &
   !$omp&        expterm, dE, dG, dS, ij, img, jat, ati, atj)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      do ij = 1, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dists2(ij, iat)
         rij = mol%xyz(:, iat) - neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         r1 = sqrt(r2)
         den2 = (en(ati) - en(atj))**2
         den4 = den2**2
         alpha = sqrt(rep(1, ati)*rep(1, atj)) &
            & * (1.0_wp+(0.01_wp*den2+0.01_wp*den4)*param%xbrad)
         repab = rep(2, ati)*rep(2, atj)
         r2top34 = r2**0.75_wp
         r2top34top2 = r2top34**2
         expterm = exp(-alpha*r2top34)*repab
         ! save repulsion energy
         dE = expterm/r1
         dG = expterm*(1.5_wp*alpha*r2top34 + 1)/r2top34top2 * rij
         dS = spread(dG, 1, 3) * spread(rij, 2, 3) * 0.5_wp
         energies(iat) = energies(iat) + 0.5_wp * dE
         sigma = sigma - dS
         if (iat /= jat) then
            energies(jat) = energies(jat) + 0.5_wp * dE
            sigma = sigma - dS
         endif
         ! save repulsion gradient
         gradient(:, iat) = gradient(:, iat) - dG
         gradient(:, jat) = gradient(:, jat) + dG
      enddo
   enddo
   !$omp end parallel do
   energy = sum(energies)

end subroutine drep_grad

! short-ranged bond correction
pure subroutine dsrb_grad(mol, neighs, neighlist, param, cn, dcndr, dcndL, &
      &                   esrb, g, sigma)

   use tbdef_param
   use tbdef_molecule
   use tbdef_neighbourlist

   use aoparam,   only : rep

   use approxrab
   use ncoord
   use pbc_tools
   use pbc

   implicit none

   ! intent in
   type(tb_molecule), intent(in) :: mol
   type(scc_parameter), intent(in) :: param
   type(tb_neighbourlist), intent(in) :: neighlist
   integer, intent(in) :: neighs(:)
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)
   ! intent inout
   real(wp), intent(inout) :: g(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   ! intent out
   real(wp), intent(out) :: esrb
   ! local variables
   integer :: i,j,k,wsAt
   integer :: iat,jat,ati,atj
   integer :: nsrb,lin
   integer :: dx,dy,dz
   real(wp), dimension(4) :: kcn
   real(wp) :: gscal
   real(wp) :: den
   real(wp) :: expterm
   ! distances
   real(wp) :: rij(3)
   real(wp) :: r2,r,dr,rab
   real(wp) :: dtmp,pre
   real(wp) :: w,t(3)
   integer :: tx,ty,tz,latrep(3)
   ! allocatables
   real(wp), allocatable :: drab0dr(:, :, :)
   real(wp), allocatable :: drab0dL(:, :, :)
   real(wp), allocatable :: rab0(:)
   integer, allocatable :: srblist(:, :)
   ! initialize
   kcn=param%kcnsh
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
subroutine build_SH0(mol, neighs, neighlist, basis,q,cn,intcut, &
      &              kmagic, ken, alphaj, kcn, xbdamp, sint, h0)

   use tbdef_molecule
   use tbdef_basisset
   use tbdef_neighbourlist

   use aoparam

   use lin_mod, only : lin
   use intgrad
   use scc_core

   implicit none

   type(tb_molecule), intent(in) :: mol
   type(tb_basisset), intent(in) :: basis
   type(tb_neighbourlist), intent(in) :: neighlist
   integer, intent(in) :: neighs(:)

   real(wp),intent(in)  :: q(:)
   real(wp),intent(in)  :: cn(:)
   real(wp),intent(in)  :: intcut

   real(wp),intent(in)  :: kmagic(:, :)
   real(wp),intent(in)  :: ken(:, :)
   real(wp),intent(in)  :: kcn(:)
   real(wp),intent(in)  :: alphaj
   real(wp),intent(in)  :: xbdamp

   real(wp),intent(out) :: sint(:, :)  ! overlap matrix S
   real(wp),intent(out) :: h0(:, :)

   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

   integer k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij,ijao,img,io,jo
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
   integer iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer itt(0:3)

   integer kat

   parameter(itt  =(/0,1,4,10/))
   real(wp) :: saw(10)

   sint = 0.0_wp
   h0   = 0.0_wp

   ! diagonal elements
   do iao = 1, basis%nao
      sint(iao,iao)=1.0_wp+sint(iao,iao)

      ii  = iao*(1+iao)/2 ! H0 is packed, note iao*(iao-1)/2+iao = iao*(1+iao)/2
      iat = basis%aoat2(iao)
      ati = mol%at(iat)
      ish = basis%ao2sh(iao)
      il  = basis%lsh(ish)+1

      ! calculate environment dependent shift
      hii = basis%level(ish) - kcnat(il-1,ati)*cn(iat) &
         &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
      H0(iao,iao) = hii

   end do

   !$omp parallel do default(none) reduction (+:sint,h0) &
   !$omp shared(basis,mol,neighlist,neighs,ao_n,ao_l,intcut,kqat,kcnat,cn,q,en, &
   !$omp&       ken,gam3,xbdamp,kcn,kmagic,kpair,alphaj) &
   !$omp private(iat,jat,ij,ijao,img,ati,cc,ci,rab2,atj,ish,valaoi,valaoj, &
   !$omp&        ri,rj,icao,naoi,iptyp,jsh,jshmax,jcao,naoj,jptyp, &
   !$omp&        ss,saw,est,alpi,alpj,ab,iprim,jprim,ip,jp,km,shpoly,io,jo, &
   !$omp&        mli,mlj,tmp,zi,zj,zetaij,enpoly,iao,jao, &
   !$omp&        ii,jj,k,den,den2,den4,il,jl,hii,hjj,hav)
   do iat = 1, len(mol)
      ri = mol%xyz(:, iat)
      ati = mol%at(iat)
      io = basis%shells(1, iat)-1
      do ij = 1, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         !r2 = neighlist%dists2(ij, iat)
         rj = neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         !r1 = sqrt(r2)

         if (iat == jat) cycle ! FIXME

         den = en(ati) - en(atj)
         den2 =  den**2
         den4 = den2**2

         jo = basis%shells(1, jat)-1
         ishells: do ish = 1, ao_n(ati)
            il = basis%lsh(ish+io)+1
            icao = basis%caoshell(ish,iat)
            naoi = llao(il-1)
            iptyp = itt(il-1)
            jshmax = ao_n(atj)
            if(iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jl = basis%lsh(jsh+jo)+1
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jl-1)
               jptyp = itt(jl-1)

               ! diagonals are the same for all H0 elements
               hii = basis%level(ish+io) - kcnat(il-1,ati)*cn(iat) &
                  &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
               hjj = basis%level(jsh+jo) - kcnat(jl-1,atj)*cn(jat) &
                  &  - kqat(jl,atj)*q(jat) - gam3(atj)*q(jat)**2

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + ken(jl,il)*den2 &
                  &      + xbdamp*ken(jl,il)*den4)

               ! we scale the two shells depending on their exponent
               zi = basis%zeta(ish+io)
               zj = basis%zeta(jsh+jo)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = kmagic(jl,il) * kpair(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = basis%valsh(ish+io).eq.0
               valaoj = basis%valsh(jsh+jo).eq.0
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

               ! distance dependent polynomial
               shpoly = rfactor(il,jl,ati,atj,ri,rj)

               ! get overlap integral
               call get_overlap(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj,point, &
                  &             intcut,basis%nprim,basis%primcount, &
                  &             basis%alp,basis%cont,ss)

               !transform from CAO to SAO
               call dtrf2(ss,il-1,jl-1)

               do ii = 1, llao2(il-1)
                  iao = ii+basis%saoshell(ish,iat)
                  do jj = 1, llao2(jl-1)
                     jao = jj+basis%saoshell(jsh,jat)
                     ijao = lin(iao,jao)
                     ! add all WSC images
                     sint(iao,jao) = sint(iao,jao) + ss(jj,ii)
                     sint(jao,iao) = sint(jao,iao) + ss(jj,ii)
                     ! Hamiltonian
                     H0(iao,jao) = H0(iao,jao) + hav * shpoly * ss(jj,ii)
                     H0(jao,iao) = H0(jao,iao) + hav * shpoly * ss(jj,ii)
                  enddo
               enddo

            enddo jshells
         enddo ishells
      enddo
   enddo
   !$omp end parallel do

end subroutine build_SH0

! ------------------------------------------------------------------------
!  Calculate the periodic AO overlap matrix
! ------------------------------------------------------------------------
subroutine ccm_build_SH0(nat,at,basis,xyz,lattice,q,cn,intcut, &
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

   real(wp),intent(out) :: sint(:,:)  ! overlap matrix S
   real(wp),intent(out) :: h0(:,:)

   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij,io,jo
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
   integer iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer itt(0:3)

   integer kat

   parameter(itt  =(/0,1,4,10/))
   real(wp) :: saw(10)

   sint = 0.0_wp
   h0   = 0.0_wp

   ! diagonal elements
   do iao = 1, basis%nao
      sint(iao,iao)=1.0_wp+sint(iao,iao)

      ii  = iao*(1+iao)/2 ! H0 is packed, note iao*(iao-1)/2+iao = iao*(1+iao)/2
      iat = basis%aoat2(iao)
      ati = at(iat)
      ish = basis%ao2sh(iao)
      il  = basis%lsh(ish)+1

      ! calculate environment dependent shift
      hii = basis%level(ish) - kcnat(il-1,ati)*cn(iat) &
         &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
      H0(iao,iao) = hii

   end do

   !$omp parallel default(none) &
   !$omp private(iat,jat,ij,ati,cc,ci,rab2,atj,ish,valaoi,valaoj, &
   !$omp&        ri,rj,icao,naoi,iptyp,jsh,jshmax,jcao,naoj,jptyp, &
   !$omp&        ss,saw,est,alpi,alpj,ab,iprim,jprim,ip,jp,km,shpoly, &
   !$omp&        mli,mlj,tmp,zi,zj,zetaij,enpoly,iao,jao,io,jo, &
   !$omp&        ii,jj,k,den,den2,den4,il,jl,hii,hjj,hav,t) &
   !$omp reduction (+:sint,h0) &
   !$omp shared(wsc,basis,at,ao_n,ao_l,xyz,lattice,intcut,nat,kqat,kcnat,cn,q,en, &
   !$omp        ken,gam3,xbdamp,kcn,kmagic,kpair,alphaj)
   !$omp do
   do iat = 1, nat
      ri  = xyz(:,iat)
      ati = at(iat)
      io = basis%shells(1, iat)-1
      do jat = 1, iat-1
         atj = at(jat)

         den = en(ati) - en(atj)
         den2 =  den**2
         den4 = den2**2

         jo = basis%shells(1, jat)-1
         ishells: do ish = 1, ao_n(ati)
            il = basis%lsh(ish+io)+1
            icao = basis%caoshell(ish,iat)
            naoi = llao(il-1)
            iptyp = itt(il-1)
            jshmax = ao_n(atj)
            if (iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jl = basis%lsh(jsh+jo)+1
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jl-1)
               jptyp = itt(jl-1)

               ! diagonals are the same for all H0 elements
               hii = basis%level(ish+io) - kcnat(il-1,ati)*cn(iat) &
                  &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
               hjj = basis%level(jsh+jo) - kcnat(jl-1,atj)*cn(jat) &
                  &  - kqat(jl,atj)*q(jat) - gam3(atj)*q(jat)**2

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + ken(jl,il)*den2 &
                  &      + xbdamp*ken(jl,il)*den4)

               ! we scale the two shells depending on their exponent
               zi = basis%zeta(ish+io)
               zj = basis%zeta(jsh+jo)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = kmagic(jl,il) * kpair(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = basis%valsh(ish+io).eq.0
               valaoj = basis%valsh(jsh+jo).eq.0
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
                  call dtrf2(ss,il-1,jl-1)

                  do ii = 1, llao2(il-1)
                     iao = ii+basis%saoshell(ish,iat)
                     do jj = 1, llao2(jl-1)
                        jao = jj+basis%saoshell(jsh,jat)
                        if (jao > iao) cycle
                        ! add all WSC images
                        sint(iao,jao) = sint(iao,jao) + ss(jj,ii)*wsc%w(jat,iat)
                        sint(jao,iao) = sint(jao,iao) + ss(jj,ii)*wsc%w(jat,iat)
                        ! Hamiltonian
                        H0(iao,jao) = H0(iao,jao) + hav * shpoly * ss(jj,ii)*wsc%w(jat,iat)
                        H0(jao,iao) = H0(jao,iao) + hav * shpoly * ss(jj,ii)*wsc%w(jat,iat)
                     enddo
                  enddo

               enddo wscAt
            enddo jshells
         enddo ishells
      enddo
   enddo
   !$OMP END DO
   !$OMP END PARALLEL

end subroutine ccm_build_SH0

! ------------------------------------------------------------------------
!  Calculate the gradient resulting from a periodic AO overlap matrix
! ------------------------------------------------------------------------
subroutine build_dSH0(mol,basis,neighs,neighlist,thr,q,cn,P,Pew,g,sigma, &
      &               dHdcn,dHdq,kmagic,ken,alphaj,kcn,xbdamp)
   use mctc_constants, only : pi
   use mctc_econv

   use tbdef_molecule
   use tbdef_wsc
   use tbdef_basisset
   use tbdef_neighbourlist

   use aoparam

   use intgrad
   use grad_core

   implicit none

   ! intent in
   type(tb_molecule), intent(in) :: mol
   type(tb_basisset), intent(in) :: basis
   type(tb_neighbourlist), intent(in) :: neighlist
   integer, intent(in) :: neighs(:)
   real(wp),intent(in) :: thr
   real(wp),intent(in) :: q(:)
   real(wp),intent(in) :: cn(:)
   real(wp),intent(in) :: P(:, :)
   real(wp),intent(in) :: Pew(:, :)

   real(wp),intent(in)  :: kmagic(:, :)
   real(wp),intent(in)  :: ken(:, :)
   real(wp),intent(in)  :: kcn(:)
   real(wp),intent(in)  :: alphaj
   real(wp),intent(in)  :: xbdamp
   ! intent inout
   real(wp),intent(inout) :: dHdq(:)
   real(wp),intent(inout) :: dHdcn(:)
   real(wp),intent(inout) :: g(:, :)
   real(wp),intent(inout) :: sigma(:, :)

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
   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,lin,jshmax,io,jo
   integer ip,jp,iat,jat,kat,ati,atj,ish,jsh,icao,jcao,iao,jao,ixyz,jxyz
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   real(wp) :: dumdum(3),dum,sdq(6,6),sdqg(3,6,6)
   real(wp),parameter :: point(3) = 0.0_wp

   integer  :: il,jl, ij, img
   real(wp) :: Hij,Pii,Pij,HPij,H0sr
   real(wp) :: den,den2,den4,enpoly
   real(wp) :: zi,zj,zetaij
   real(wp) :: hii,hjj,hav
   real(wp) :: shpoly,dshpoly(3),km
   logical  :: valaoi,valaoj
   real(wp) :: g_xyz(3)

   !$omp parallel do default(none) &
   !$omp reduction(+:g,sigma,dhdcn,dhdq) &
   !$omp private(iat,jat,ixyz,ati,cc,ci,rab2,atj,ish,ishtyp,valaoi,valaoj,g_xyz, &
   !$omp&        ri,rj,icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp,rij2, &
   !$omp&        sdq,sdqg,est,alpi,alpj,ab,iprim,jprim,ip,jp,km,shpoly,dshpoly, &
   !$omp&        mli,mlj,dum,dumdum,tmp,stmp,rij,zi,zj,zetaij,enpoly,iao,jao,img, &
   !$omp&        ii,jj,k,den,den2,den4,io,jo,il,jl,hii,hjj,hav,Pij,Hij,HPij,H0sr) &
   !$omp shared(mol,neighs,neighlist,basis,ao_n,ao_l,thr,kqat,kcnat,cn,q,en, &
   !$omp        ken,gam3,xbdamp,kcn,kmagic,kpair,alphaj,P,Pew)
   do iat = 1, len(mol)
      ri  = mol%xyz(:,iat)
      ati = mol%at(iat)
      io = basis%shells(1, iat)-1
      do ij = 1, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         rij2 = neighlist%dists2(ij, iat)
         rj = neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         rij = ri - rj

         if (iat == jat) cycle ! FIXME

         den = en(ati) - en(atj)
         den2 =  den**2
         den4 = den2**2

         jo = basis%shells(1, jat)-1
         ishells: do ish = 1, ao_n(ati)
            il = basis%lsh(ish+io)+1
            icao = basis%caoshell(ish,iat)
            naoi = llao(il-1)
            iptyp = itt(il-1)
            jshmax = ao_n(atj)
            if(iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jl = basis%lsh(jsh+jo)+1
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jl-1)
               jptyp = itt(jl-1)

               ! diagonals are the same for all H0 elements
               hii = basis%level(ish+io) - kcnat(il-1,ati)*cn(iat) &
                  &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
               hjj = basis%level(jsh+jo) - kcnat(jl-1,atj)*cn(jat) &
                  &  - kqat(jl,atj)*q(jat) - gam3(atj)*q(jat)**2

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + ken(jl,il)*den2 &
                  &      + xbdamp*ken(jl,il)*den4)

               ! we scale the two shells depending on their exponent
               zi = basis%zeta(ish+io)
               zj = basis%zeta(jsh+jo)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = kmagic(jl,il) * kpair(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = basis%valsh(ish+io).eq.0
               valaoj = basis%valsh(jsh+jo).eq.0
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

               ! distance dependent polynomial
               !shpoly = rfactor(il,jl,ati,atj,ri,rj)
               call drfactor(il,jl,0,ati,atj,rij2,ri,rj,shpoly,dshpoly)

               call get_grad_overlap(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj, &
                  &                  point,thr,basis%nprim,basis%primcount, &
                  &                  basis%alp,basis%cont,sdq,sdqg)

               call dtrf2(sdq,il-1,jl-1)
               do ixyz = 1, 3
                  ! transform from CAO to SAO, only transform the gradient
                  tmp = sdqg(ixyz,:,:)
                  call dtrf2(tmp,il-1,jl-1)
                  sdqg(ixyz,:,:) = tmp
               enddo

               do ii = 1, llao2(il-1)
                  iao = ii+basis%saoshell(ish,iat)
                  do jj = 1, llao2(jl-1)
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
                     sigma = sigma + spread(g_xyz, 1, 3) * spread(rij, 2, 3)

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
   !$omp end parallel do

   ! diagonal contributions
   do iao = 1, basis%nao
      ii  = iao*(1+iao)/2 ! H0 is packed, note iao*(iao-1)/2+iao = iao*(1+iao)/2
      iat = basis%aoat2(iao)
      ati = mol%at(iat)
      ish = basis%ao2sh(iao)
      il  = basis%lsh(ish)+1

      Pii = P(iao,iao)*evtoau

      ! save dE/dCN for CNi
      dhdcn(iat) = dhdcn(iat) - Pii*kcnat(il-1,ati)
      ! save dE/dq for qi
      dhdq(iat) = dhdq(iat) - Pii*kqat(il,ati) - Pii*gam3(ati)*2*q(iat)
   enddo

end subroutine build_dSH0

! ------------------------------------------------------------------------
!  Calculate the gradient resulting from a periodic AO overlap matrix
! ------------------------------------------------------------------------
subroutine ccm_build_dSH0(nat,basis,thr,at,xyz,lattice,q,cn,P,Pew,g,sigma,&
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
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   real(wp),intent(in) :: lattice(3,3)
   real(wp),intent(in) :: q(nat)
   real(wp),intent(in) :: cn(nat)
   real(wp),intent(in) :: P(:,:)
   real(wp),intent(in) :: Pew(:,:)

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
   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij,lin,jshmax,io,jo
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
   !$omp&        ii,jj,k,den,den2,den4,io,jo,il,jl,hii,hjj,hav,Pij,Hij,HPij,H0sr) &
   !$omp reduction (+:g,sigma,dhdcn,dhdq) &
   !$omp shared(wsc,basis,at,ao_n,ao_l,xyz,lattice,thr,nat,kqat,kcnat,cn,q,en, &
   !$omp        ken,gam3,xbdamp,kcn,kmagic,kpair,alphaj,P,Pew)
   !$omp do
   do iat = 1, nat
      ri  = xyz(:,iat)
      ati = at(iat)
      io = basis%shells(1, iat)-1
      do jat = 1, iat-1
         atj = at(jat)

         den = en(ati) - en(atj)
         den2 =  den**2
         den4 = den2**2

         jo = basis%shells(1, jat)-1
         ishells: do ish = 1, ao_n(ati)
            il = basis%lsh(ish+io)+1
            icao = basis%caoshell(ish,iat)
            naoi = llao(il-1)
            iptyp = itt(il-1)
            jshmax = ao_n(atj)
            if(iat == jat) jshmax = ish

            jshells: do jsh = 1, jshmax
               jl = basis%lsh(jsh+jo)+1
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jl-1)
               jptyp = itt(jl-1)

               ! diagonals are the same for all H0 elements
               hii = basis%level(ish+io) - kcnat(il-1,ati)*cn(iat) &
                  &  - kqat(il,ati)*q(iat) - gam3(ati)*q(iat)**2
               hjj = basis%level(jsh+jo) - kcnat(jl-1,atj)*cn(jat) &
                  &  - kqat(jl,atj)*q(jat) - gam3(atj)*q(jat)**2

               ! evaluate the EN polynom for this shells
               enpoly = (1.0_wp + ken(jl,il)*den2 &
                  &      + xbdamp*ken(jl,il)*den4)

               ! we scale the two shells depending on their exponent
               zi = basis%zeta(ish+io)
               zj = basis%zeta(jsh+jo)
               zetaij = 2 * sqrt(zi*zj)/(zi+zj)

               ! now do the real magic (called EHT enhancement factor)
               km = kmagic(jl,il) * kpair(ati,atj) * zetaij * enpoly

               ! check for valence orbitals
               valaoi = basis%valsh(ish+io).eq.0
               valaoj = basis%valsh(jsh+jo).eq.0
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

                  call dtrf2(sdq,il-1,jl-1)
                  do ixyz = 1, 3
                     ! transform from CAO to SAO, only transform the gradient
                     tmp = sdqg(ixyz,:,:)
                     call dtrf2(tmp,il-1,jl-1)
                     sdqg(ixyz,:,:) = tmp
                  enddo

                  do ii = 1, llao2(il-1)
                     iao = ii+basis%saoshell(ish,iat)
                     do jj = 1, llao2(jl-1)
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
   do iao = 1, basis%nao
      ii  = iao*(1+iao)/2 ! H0 is packed, note iao*(iao-1)/2+iao = iao*(1+iao)/2
      iat = basis%aoat2(iao)
      ati = at(iat)
      ish = basis%ao2sh(iao)
      il  = basis%lsh(ish)+1

      Pii = P(iao,iao)*evtoau

      ! save dE/dCN for CNi
      dhdcn(iat) = dhdcn(iat) - Pii*kcnat(il-1,ati)
      ! save dE/dq for qi
      dhdq(iat) = dhdq(iat) - Pii*kqat(il,ati) - Pii*gam3(ati)*2*q(iat)
   enddo

end subroutine ccm_build_dSH0

real(wp) elemental function eeqthr_maxval(self) result(cutoff)
   type(eeq_thresholds), intent(in) :: self
   cutoff = max(self%cutoff_disp, &
      &         self%cutoff_cn, &
      &         self%cutoff_ham, &
      &         self%cutoff_srb, &
      &         self%cutoff_rep)
end function eeqthr_maxval

!> Generate a set of numerical thresholds from a given accuracy value.
type(eeq_thresholds) elemental function new_thresholds(accuracy) result(thr)
   real(wp), intent(in) :: accuracy
   !> primitive cut-off
   thr%intcut = max(20.0_wp, 25.0_wp - 10.0_wp*log10(accuracy))
   !> integral neglect threshold
   thr%neglect = 10.0e-9_wp * accuracy
   !> coordination number cutoff
   thr%cutoff_cn = 30.0_wp - log10(accuracy)*10.0_wp
   !> integral distance cutoff
   thr%cutoff_ham = sqrt(800.0_wp) - log10(accuracy)*10.0_wp
   !> repulsion polynomial cutoff
   thr%cutoff_rep = 40.0_wp - log10(accuracy)*20.0_wp
   !> short range basis cutoff
   thr%cutoff_srb = 14.0_wp - log10(accuracy)*10.0_wp
   !> long-ranged dispersion cutoff
   thr%cutoff_disp = sqrt(4000.0_wp) - log10(accuracy)*30.0_wp
end function new_thresholds

end module peeq_module

