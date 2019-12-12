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

!> Implementation of the self-consistent charge extended tight binding Hamiltonians
module scf_module
! ========================================================================
   use iso_fortran_env, wp => real64
   implicit none
   public :: scf, scf_options
   private

   !> Profile the steps in the evaluation by saving timings for each step.
   logical, parameter :: profile = .true.

   !> options for the SCF procedure
   type :: scf_options
      integer :: prlevel
      integer :: maxiter = 250
      logical :: lqpc = .true.
      logical :: ccm = .false.
      real(wp) :: cf
      real(wp) :: etemp = 300.0_wp
      real(wp) :: accuracy = 1.0_wp
   end type scf_options

   !> Numerical thresholds
   type :: scf_thresholds
      !> primitive cut-off
      real(wp) :: intcut
      !> integral neglect threshold
      real(wp) :: neglect
      !> convergence threshold of the energy
      real(wp) :: scfconv
      !> convergence threshold for the charges (and moments)
      real(wp) :: qconv
      real(wp) :: cutoff_disp
      real(wp) :: cutoff_cn
      real(wp) :: cutoff_ham
      real(wp) :: cutoff_rep
      real(wp) :: cutoff_atm
   end type scf_thresholds

   interface scf_thresholds
      module procedure :: new_thresholds
   end interface

   interface maxval
      module procedure :: scfthr_maxval
   end interface

   interface
      module subroutine build_SDH0(mol, neighs, neighlist, basis, param, intcut, &
            &                      cn, kcnsh, sint, dpint, H0)
         use tbdef_molecule
         use tbdef_neighbourlist
         use tbdef_basisset
         use tbdef_param
         type(tb_molecule), intent(in) :: mol
         type(tb_neighbourlist), intent(in) :: neighlist
         type(tb_basisset), intent(in) :: basis
         type(scc_parameter), intent(in) :: param
         integer, intent(in) :: neighs(:)
         real(wp), intent(in) :: intcut
         real(wp), intent(in) :: cn(:)
         real(wp), intent(in) :: kcnsh(:)
         real(wp), intent(out) :: H0(:)
         real(wp), intent(out) :: sint(:,:)
         real(wp), intent(out) :: dpint(:,:,:)
      end subroutine build_SDH0
      module subroutine build_SDQH0(mol, neighs, neighlist, basis, param, intcut, &
            &                       cn, sint, dpint, qpint, H0)
         use tbdef_molecule
         use tbdef_neighbourlist
         use tbdef_basisset
         use tbdef_param
         type(tb_molecule), intent(in) :: mol
         type(tb_neighbourlist), intent(in) :: neighlist
         type(tb_basisset), intent(in) :: basis
         type(scc_parameter), intent(in) :: param
         integer, intent(in) :: neighs(:)
         real(wp), intent(in) :: intcut
         real(wp), intent(in) :: cn(:)
         real(wp), intent(out) :: H0(:)
         real(wp), intent(out) :: sint(:,:)
         real(wp), intent(out) :: dpint(:,:,:)
         real(wp), intent(out) :: qpint(:,:,:)
      end subroutine build_SDQH0
      module subroutine build_dSDH0(mol, neighs, neighlist, basis, param, intcut, &
            &                       cn, kcnsh, P, Pew, ves, dhdcn, gradient, sigma)
         use tbdef_molecule
         use tbdef_neighbourlist
         use tbdef_basisset
         use tbdef_param
         type(tb_molecule), intent(in) :: mol
         type(tb_neighbourlist), intent(in) :: neighlist
         type(tb_basisset), intent(in) :: basis
         type(scc_parameter), intent(in) :: param
         integer, intent(in) :: neighs(:)
         real(wp), intent(in) :: intcut
         real(wp), intent(in) :: cn(:)
         real(wp), intent(in) :: kcnsh(:)
         real(wp), intent(in) :: P(:, :)
         real(wp), intent(in) :: Pew(:, :)
         real(wp), intent(in) :: ves(:)
         real(wp), intent(out) :: dhdcn(:)
         real(wp), intent(inout) :: gradient(:,:)
         real(wp), intent(inout) :: sigma(:,:)
      end subroutine build_dSDH0
      module subroutine build_dSDQH0(mol, neighs, neighlist, basis, param, &
            &                        intcut, cn, P, Pew, ves, vs, vd, vq, dhdcn, &
            &                        gradient, sigma)
         use tbdef_molecule
         use tbdef_neighbourlist
         use tbdef_basisset
         use tbdef_param
         type(tb_molecule), intent(in) :: mol
         type(tb_neighbourlist), intent(in) :: neighlist
         type(tb_basisset), intent(in) :: basis
         type(scc_parameter), intent(in) :: param
         integer, intent(in) :: neighs(:)
         real(wp), intent(in) :: intcut
         real(wp), intent(in) :: cn(:)
         real(wp), intent(in) :: ves(:)
         real(wp), intent(in) :: vs(:)
         real(wp), intent(in) :: vd(:,:)
         real(wp), intent(in) :: vq(:,:)
         real(wp), intent(in) :: P(:,:)
         real(wp), intent(in) :: Pew(:,:)
         real(wp), intent(out) :: dhdcn(:)
         real(wp), intent(inout) :: gradient(:,:)
         real(wp), intent(inout) :: sigma(:,:)
      end subroutine build_dSDQH0
   end interface

contains

!> Wrapper for the implemented self-consistent charge Hamiltonians.
subroutine scf(iunit, mol, wfn, basis, param, pcem, opt, &
      &        egap, etotal, g, sigma, res)
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data
   use tbdef_timer
   use tbdef_pcem
   use tbdef_neighbourlist
   use setparam
   !> Unit identifier for all IO, only used for prlevel > 0.
   integer, intent(in) :: iunit
   !> Molecular structure information.
   type(tb_molecule), intent(in) :: mol
   !> Tight binding wavefunction, contains previous wavefunction for restart
   !  and is updated after this run.
   type(tb_wavefunction), intent(inout) :: wfn
   !> Atomic orbital basis set.
   type(tb_basisset), intent(in) :: basis
   !> Global method parameters.
   type(scc_parameter), intent(in) :: param
   !> External potentials, like point charge fields.
   type(tb_pcem), intent(inout) :: pcem
   !> HOMO-LUMO gap in eV.
   real(wp), intent(inout) :: egap
   !> Total energy in Eh.
   real(wp), intent(inout) :: etotal
   !> Gradient in Eh/a0.
   real(wp), intent(inout) :: g(:, :)
   !> Stress tensor in Eh (volume contribution not included).
   real(wp), intent(inout) :: sigma(:, :)
   !> Detailed results on the calculation.
   type(scc_results), intent(out) :: res
   !> Options for the SCF
   type(scf_options), intent(in) :: opt
   select case(gfn_method)
   case(1)
      call scf1(iunit, mol, wfn, basis, param, pcem, opt, &
         &      egap, etotal, g, sigma, res)
   case(2)
      call scf2(iunit, mol, wfn, basis, param, pcem, opt, &
         &      egap, etotal, g, sigma, res)
   end select
end subroutine scf

!> Implementation of the GFN1-xTB Hamiltonian and related parametrisations.
subroutine scf1(iunit, mol, wfn, basis, param, pcem, opt, &
      &         egap, etotal, g, sigma, res)

   use mctc_econv, only : autoev,evtoau

   ! type definitions
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data
   use tbdef_timer
   use tbdef_pcem
   use tbdef_neighbourlist

   ! global storage
   use aoparam
   use setparam

   ! interfaces
   use scc_core
   use grad_core
   use tbmod_dftd3
   use aespot
   use gbobc
   use ncoord, only: get_coordination_number, tb_cn_type
   use embedding
   use readin
   use lidep

   !> Unit identifier for all IO, only used for prlevel > 0.
   integer, intent(in) :: iunit
   !> Molecular structure information.
   type(tb_molecule), intent(in) :: mol
   !> Tight binding wavefunction, contains previous wavefunction for restart
   !  and is updated after this run.
   type(tb_wavefunction), intent(inout) :: wfn
   !> Atomic orbital basis set.
   type(tb_basisset), intent(in) :: basis
   !> Global method parameters.
   type(scc_parameter), intent(in) :: param
   !> External potentials, like point charge fields.
   type(tb_pcem), intent(inout) :: pcem
   !> Options for the SCF.
   type(scf_options), intent(in) :: opt
   !> HOMO-LUMO gap in eV.
   real(wp), intent(inout) :: egap
   !> Total energy in Eh.
   real(wp), intent(inout) :: etotal
   !> Gradient in Eh/a0.
   real(wp), intent(inout) :: g(:, :)
   !> Detailed results on the calculation.
   type(scc_results), intent(out) :: res
   !> Stress tensor in Eh (volume contribution not included).
   real(wp), intent(inout) :: sigma(:, :)
   !> Static neighbourlist including all images up to a certain cutoff radius.
   type(tb_neighbourlist) :: neighlist
   !> Generator for neighbourlist.
   type(tb_neighlist_generator) :: neighgen
   !> Number of neighbours for each atom up to a certain cutoff radius.
   integer, allocatable :: neighs(:)

   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: X(:,:)
   real(wp),allocatable :: S(:,:)
   real(wp),allocatable :: H(:,:)
   real(wp),allocatable :: H0(:)
   real(wp),allocatable :: H1(:)
   real(wp),allocatable :: jab(:,:)
   real(wp),allocatable :: ves(:) ! shell ES potential
   real(wp),allocatable :: tmp(:)
   real(wp),allocatable :: zsh(:)
   real(wp),allocatable :: qq(:)
   real(wp),allocatable :: qlmom(:,:)
   real(wp),allocatable :: cm5(:)
   real(wp),allocatable :: kcnsh(:)
   real(wp), allocatable :: gam2sh(:)
!  AES stuff
   integer :: ndp,nqp
   real(wp),allocatable  :: dpint(:,:,:),qpint(:,:,:)

   integer,allocatable :: matlist (:,:)
   integer,allocatable :: matlist2(:,:)
   integer,allocatable :: xblist(:,:)
   real(wp),allocatable :: sqrab(:)
   real(wp),allocatable :: dHdcn(:)
   real(wp),allocatable :: dcndr(:,:,:)
   real(wp),allocatable :: dcndL(:,:,:)

   real(wp) :: eh1
   real(wp) :: dipol(3),dip,gsolv,eat
   real(wp) :: eelec,ep,dum1
   real(wp) :: ed
   real(wp) :: rmsq
   real(wp) :: ees,repab,jmpol,kdampxb,vvxb,exb
   real(wp) :: eatoms,gnorm
   real(wp) :: damp,damp0,exc
   real(wp) :: eaes,epol,kexpe
   type(scf_thresholds) :: thr

!  some parameter defaults which are not fitted
   !> rep exp for exp(-R**kexp)
   real(wp), parameter :: kexp = 1.50_wp
   !> rep exp in 1/R**rexp
   real(wp), parameter :: rexp = 1.00_wp
   !> XB parameter for damped LJ in GFN1
   real(wp), parameter ::  ljexp = 12.0_wp

   integer, external :: lin, tmmetal
   integer :: ii,jj,kk,i,j,k,m,iat,jat,mi,jter,atj,kkk,mj,mm,ish,jsh
   integer :: ishell,jshell,nmat,nmat2
   integer :: ati,nxb,startpdiag
   integer, parameter :: lladr(4) = (/1,3,6,10/)
   integer, parameter :: lladr2(0:3) = (/1,3,5,7/)

   logical, external :: xbond, early3d
   logical :: minpr,pr,lastdiag,iniqsh,fail

!  GBSA stuff
   type(tb_solvent) :: gbsa
   real(wp),allocatable :: fgb(:,:)
   real(wp),allocatable :: fhb(:)
   real(wp) :: gborn,ghb
!  for the CM5 charges
   real(wp),allocatable :: cm5a(:)
   real(wp),allocatable :: dcm5a(:,:,:)
   real(wp),allocatable :: vborn(:)

!  point charge embedding stuff
   logical  :: lpcem
   real(wp),allocatable :: vpc(:)
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
   logical :: broy

! ------------------------------------------------------------------------
!  initialization
! ------------------------------------------------------------------------
   if (profile) call timer%new(8,.false.)
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
   etotal  = 0.0_wp
   ep   = 0.0_wp
   ed   = 0.0_wp
   eat  = 0.0_wp
   egap = 0.0_wp

   pr   = opt%prlevel.gt.1
   minpr= opt%prlevel.gt.0

   thr = scf_thresholds(opt%accuracy, 1)

   !> old/new q mixing start
   damp0=0.20
   broy = mol%n > 1
   !> when to start pseudodiag (=1 after one full diag)
   startpdiag=1000 !large number=never

!ccccccccccccccccccc
! note: H is in eV!
!ccccccccccccccccccc

   ! # atom arrays
   allocate(qq(mol%n),qlmom(3,mol%n),cm5(mol%n),sqrab(mol%n*(mol%n+1)/2), &
      &     dcndr(3,mol%n,mol%n),cn(mol%n),dcndL(3,3,mol%n),vborn(mol%n), &
      &     dHdcn(mol%n))

   ! initialize the GBSA module (GBSA works with CM5 charges)
   if(lgbsa) then
      call new_gbsa(gbsa,mol%n,mol%at)
      allocate(fgb(mol%n,mol%n),fhb(mol%n),cm5a(mol%n),dcm5a(3,mol%n,mol%n))
      gborn=0._wp
      gbsa%gsasa=0._wp
      gbsa%ghb=0._wp
      qq=0._wp
      ! initialize the neighbor list
      call update_nnlist_gbsa(gbsa,mol%xyz,.false.)
      ! compute Born radii
      call compute_brad_sasa(gbsa,mol%xyz)
      ! initialize the fgb matrix (dielectric screening of the Coulomb potential)
      call compute_fgb(gbsa,fgb,fhb)
      ! initialize the CM5 charges computation
      call calc_cm5(mol%n,mol%at,mol%xyz,qq,cm5,cm5a,dcm5a)
   endif

   allocate(H0(basis%nao*(basis%nao+1)/2), &
      &     S(basis%nao,basis%nao),tmp(basis%nao), &
      &     X(basis%nao,basis%nao),H1(basis%nao*(basis%nao+1)/2), &
      &     kcnsh(basis%nshell),ves(basis%nshell), &
      &     zsh(basis%nshell),&
      &     jab(basis%nshell,basis%nshell), &
      &     matlist (2,basis%nao*(basis%nao+1)/2), &
      &     matlist2(2,basis%nao*(basis%nao+1)/2-basis%nao))

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
   do iat=1,mol%n
      do jat=1,iat
!        save distances
         k=k+1
         sqrab(k)=(mol%xyz(1,iat)-mol%xyz(1,jat))**2 &
         &       +(mol%xyz(2,iat)-mol%xyz(2,jat))**2 &
         &       +(mol%xyz(3,iat)-mol%xyz(3,jat))**2
!        XB
         if(xbond(mol%at(iat),mol%at(jat)).and.sqrab(k).lt.400.) nxb=nxb+1
      enddo
   enddo

   ! XB part
   allocate(xblist(3,nxb+1))
   nxb = 0
   k = 0
   do iat=1,mol%n
      ati=mol%at(iat)
      do jat=1,iat
         atj=mol%at(jat)
         k=k+1
         if(xbond(ati,atj).and.sqrab(k).lt.400) then
            nxb=nxb+1
            dum1=1.d+42
            if(ati.eq.17.or.ati.eq.35.or.ati.eq.53.or.ati.eq.85)then
               xblist(1,nxb)=iat
               xblist(2,nxb)=jat
               do m=1,mol%n
                  if(m.ne.iat.and.sqrab(lin(m,iat)).lt.dum1)then
                     dum1=sqrab(lin(m,iat))
                     xblist(3,nxb)=m
                  endif
               enddo
            endif
            if(atj.eq.17.or.atj.eq.35.or.atj.eq.53.or.atj.eq.85)then
               xblist(1,nxb)=jat
               xblist(2,nxb)=iat
               do m=1,mol%n
                  if(m.ne.jat.and.sqrab(lin(m,jat)).lt.dum1)then
                     dum1=sqrab(lin(m,jat))
                     xblist(3,nxb)=m
                  endif
               enddo
            endif
         endif
      enddo
   enddo

   ! ldep J potentials (in Eh) for SCC
   allocate(gam2sh(basis%nshell), source=0.0_wp)
   do ishell = 1, basis%nshell
      iat=basis%ash(ishell)
      ati=mol%at(iat)
      gam2sh(ishell) = gam(ati)*(1.0_wp+lpar(basis%lsh(ishell),ati))
   enddo
   call get_gfn_coulomb_matrix(mol, basis%nshell, basis%ash, gam2sh, gfn_method, &
      &                        opt%cf, opt%lqpc, jab)

   ! J potentials including the point charge field
   if(lpcem)then
      allocate(vpc(basis%nshell), source=0.0_wp)
      call jpot_pcem_gfn1(mol%n,pcem,basis%nshell,mol%at,mol%xyz,basis%ash, &
         &                basis%lsh,param%alphaj,vpc)
   endif

   if (opt%prlevel > 1) then
      write(iunit,'(/,10x,51("."))')
      write(iunit,'(10x,":",22x,a,22x,":")') "SETUP"
      write(iunit,'(10x,":",49("."),":")')
      write(iunit,intfmt) "# basis functions  ",basis%nbf
      write(iunit,intfmt) "# atomic orbitals  ",basis%nao
      write(iunit,intfmt) "# shells           ",basis%nshell
      write(iunit,intfmt) "# electrons        ",wfn%nel
      write(iunit,intfmt) "# halogen bonds    ",nxb
      write(iunit,intfmt) "max. iterations    ",opt%maxiter
      write(iunit,chrfmt) "Hamiltonian        ","GFN1-xTB"
      write(iunit,chrfmt) "GBSA solvation     ",bool2string(lgbsa)
      write(iunit,chrfmt) "PC potential       ",bool2string(lpcem)
      if (lpcem) then
         write(iunit,intfmt) "-> # point charges ",pcem%n
         write(iunit,dblfmt) "-> sum of PC       ",sum(pcem%q),"e   "
      endif
      write(iunit,dblfmt) "electronic temp.   ",opt%etemp,"K   "
      write(iunit,dblfmt) "accuracy           ",opt%accuracy,"    "
      write(iunit,scifmt) "-> integral cutoff ",thr%intcut,  "    "
      write(iunit,scifmt) "-> integral neglect",thr%neglect, "    "
      write(iunit,scifmt) "-> SCF convergence ",thr%scfconv, "Eh  "
      write(iunit,scifmt) "-> wf. convergence ",thr%qconv,   "e   "
      write(iunit,dblfmt) "Broyden damping    ",broydamp,"    "
      write(iunit,'(10x,51("."))')
   endif

   qq = wfn%q
   damp = damp0

   do ii = 1, basis%nshell
      iat=basis%ash(ii)
      ishell=basis%lsh(ii)+1
      kcnsh(ii)=basis%level(ii)*param%kcnsh(ishell)
      if(metal(mol%at(iat)).eq.1) kcnsh(ii)=0.0_wp  ! CN dep. bad for metals
      if(early3d(mol%at(iat))) then
         kcnsh(ii)=basis%level(ii)*param%kcnsh(ishell)
         ! fix problems with too low-coord CP rings
         if(ishell.eq.3) kcnsh(ii)=basis%level(ii)*param%kcnsh(4)
      endif
   enddo

   if (profile) call timer%measure(1)
   if (profile) call timer%measure(8,"Neighbourlist")

   call neighgen%new(maxval(thr), len(mol), mol%xyz, mol%lattice, mol%npbc > 0)
   call neighlist%new(neighgen)
   allocate(neighs(len(mol)), source=0)

   if (profile) call timer%measure(8)
   if (profile) call timer%measure(2,"Dispersion")

   call neighlist%get_neighs(neighs, cutoff=thr%cutoff_cn)
   call get_coordination_number(mol, neighs, neighlist, tb_cn_type%exp, &
      &                         cn, dcndr, dcndL)
   call neighlist%get_neighs(neighs, cutoff=thr%cutoff_disp)
   call d3_gradient(mol, neighs, neighlist, param%disp, 4.0_wp, &
      &             cn, dcndr, dcndL, ed, g, sigma)

   if (profile) call timer%measure(2)
   if (profile) call timer%measure(6,"classical contributions")

   exb = 0.0_wp
   call xbpot(mol%n,mol%at,mol%xyz,sqrab,xblist,nxb,param%xbdamp,param%xbrad, &
      &       ljexp,exb,g)

   call neighlist%get_neighs(neighs, cutoff=thr%cutoff_rep)
   call rep_grad_gfn1(mol, neighs, neighlist, kexp, rexp, ep, g, sigma)

   if (profile) call timer%measure(6)
   if (profile) call timer%measure(3,"integral evaluation")

   allocate(dpint(3,basis%nao,basis%nao), source = 0.0_wp)
   call neighlist%get_neighs(neighs, cutoff=thr%cutoff_ham)
   ! compute integrals and prescreen to set up list arrays
   call build_SDH0(mol, neighs, neighlist, basis, param, thr%intcut, cn, kcnsh, &
      &            S, dpint, H0)

   call cholesky(iunit,pr,basis%nao,S,orthog) ! S is not modified
   if (orthog) then
      call prmat(6, S, basis%nao, basis%nao, 'overlap')
      call raise('E', 'System is linear dependent', 1)
   endif

   if (profile) call timer%measure(3)
   if (profile) call timer%measure(1)

   ! set first electrostatic potential (start with third order on-side)
   vborn = wfn%q**2 * gam3(mol%at) * autoev ! abuse of vborn
   ves = vborn(basis%ash)
   ! add Born shifts to potential
   if (lgbsa) then
      cm5 = wfn%q + cm5a
      call get_born_shift(mol%n, cm5, fgb, fhb, vborn)
      ves = ves + vborn(basis%ash)
   else
      vborn = 0.0_wp
   endif
   ! add external potential
   if (lpcem) then
      ves = ves + vpc
   endif
   ! get second order electrostatic potential
   call get_charge_shift(basis%nshell,wfn%qsh,jab,ves)

   ! prepare matrix indices
   nmat =0
   nmat2=0
   do ii=1,basis%nao
      iat=basis%aoat2(ii)
      do jj=1,ii-1
         jat=basis%aoat2(jj)
         if(abs(S(jj,ii)).lt.thr%neglect) then
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
   enddo

   if (profile) call timer%measure(1)
   if (profile) call timer%measure(4,"zeroth order Hamiltonian")

   if(pr)then
      write(iunit,'(a)')
      write(iunit,*) 'iter      E             dE          RMSdq', &
      &'      gap      omega  full diag'
   endif

   ! first order energy for given geom. and density, i.e. skip SCC and grad
   if(opt%maxiter.eq.0) then
      call qsh2qat(mol%n,mol%at,basis%nshell,wfn%qsh,wfn%q)
      call electro(mol%n,mol%at,basis%nao,basis%nshell,jab,H0,wfn%P, &
         &         wfn%q,wfn%qsh,ees,etotal)
      if(lgbsa) then
         cm5=wfn%q+cm5a
         call electro_gbsa(mol%n,fgb,fhb,cm5,gborn,etotal)
      endif
      goto 9999
   endif
   if (profile) call timer%measure(4)

! ========================================================================
!  SCC iterations
! ========================================================================
   if (profile) call timer%measure(5,"iterations")
   call scc_gfn1(iunit,mol%n,wfn%nel,wfn%nopen,basis%nao,nmat,basis%nshell, &
      &          mol%at,matlist,basis%aoat2,basis%ao2sh,basis%ash, &
      &          wfn%q,qq,qlmom,wfn%qsh,zsh, &
      &          gbsa,fgb,fhb,cm5,cm5a,gborn,vborn, &
      &          broy,broydamp,damp0, &
      &          lpcem,ves,vpc, &
      &          opt%etemp,wfn%focc,wfn%focca,wfn%foccb,wfn%efa,wfn%efb, &
      &          etotal,ees,epcem,egap,wfn%emo,wfn%ihomo,wfn%ihomoa,wfn%ihomob, &
      &          H0,H1,wfn%C,S,X,wfn%P,jab, &
      &          opt%maxiter,startpdiag,thr%scfconv,thr%qconv, &
      &          minpr,pr, &
      &          fail,jter)
! ========================================================================
   ! free some memory (this stuff is not needed for gradients)
   if (allocated(vpc)) deallocate(vpc)

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

   ! ========================================================================
   !  GRADIENT (now 100% analytical (that's not true!))
   ! ========================================================================

   ! get energy weighted density matrix
   tmp = wfn%focc * wfn%emo*evtoau
   call dmat(basis%nao,tmp,wfn%C,X)

   !call neighlist%get_neighs(neighs, cutoff=thr%cutoff_ham)
   ! wave function terms/overlap dependent parts of H
   call build_dSDH0(mol, neighs, neighlist, basis, param, thr%intcut, cn, kcnsh, &
      &             wfn%P, X, ves, dhdcn, g, sigma)
   ! CN-level shift gradient
   call dgemv('n', 3*mol%n, mol%n,-1.0_wp, dcndr, 3*mol%n, dhdcn, 1, 1.0_wp, g, 1)
   call dgemv('n', 9, mol%n, 1.0_wp, dcndL, 9, dhdcn, 1, 1.0_wp, sigma, 1)

   ! GBSA
   if (lgbsa) then
      qq=qq+cm5a
      call compute_gb_egrad(gbsa,qq,gborn,ghb,g,minpr)
      call cm5_grad_gfn1(g,mol%n,qq,fgb,fhb,dcm5a,lhb)
      ! solvation energy
      gbsa%gborn = gborn
      gbsa%ghb = ghb
      gsolv=gborn+gbsa%gsasa+gbsa%ghb+gshift
      etotal=etotal+gbsa%gsasa+gshift
   endif

   ! Calculating shell es gradient
   call get_gfn_coulomb_derivs(mol, basis%nshell, basis%ash, gam2sh, gfn_method, &
      &                        opt%cf, opt%lqpc, wfn%qsh, g, sigma)

   ! --- ES point charge embedding
   if (lpcem) then
      call pcem_grad_gfn1(g,pcem%grd,mol%n,pcem,mol%at,basis%nshell,mol%xyz, &
         &                basis%ash,basis%lsh,param%alphaj,wfn%qsh)
   endif

!  calculate the norm for printout
   gnorm = sqrt(sum( g**2 ))

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
!     HOMO-LUMO excitation properties if  UHF=2        
      if (wfn%nopen.eq.2) then
         call hlex(mol%n,mol%at,basis%nbf,basis%nao,wfn%ihomoa,mol%xyz,wfn%focc,S,wfn%C,wfn%emo,basis)
      endif

! ------------------------------------------------------------------------
!     LMO /xTB-IFF
      if (pr_lmo) then
         tmp=wfn%emo*evtoau
         call local(mol%n,mol%at,basis%nbf,basis%nao,wfn%ihomoa,mol%xyz,mol%z,wfn%focc,S,wfn%P,wfn%C,tmp,wfn%q,etotal,lgbsa,basis)
      endif

! ------------------------------------------------------------------------
!  exchange energy correction ala sTDA
      if (wfn%nopen.ge.2) then
         call exch(mol%n,mol%at,basis%nao,wfn%nopen,wfn%ihomoa,mol%xyz,wfn%focc,S,wfn%C,exc,basis%aoat)
         write(iunit,'(''open-shell EX :'',F16.7)') -exc
         write(iunit,'(''corrected Etot:'',F16.7, &
         &   '' (not used further except for this printout!)'')') etotal - exc
      endif

   endif printing

! ------------------------------------------------------------------------
!  get Wiberg bond orders
   call get_wiberg(mol%n,basis%nao,mol%at,mol%xyz,wfn%P,S,wfn%wbo,basis%fila2)

! ------------------------------------------------------------------------
!  dipole calculation (always done because its free)
   !call mmompop(mol%n,basis%nao,basis%aoat2,mol%xyz,wfn%p,s,dpint,qpint, &
      !&         wfn%dipm,wfn%qp)

   call calc_dipole(mol%n,mol%at,mol%xyz,mol%z,basis%nao,wfn%P,dpint,dip,dipol)

   if (profile) call timer%measure(7)

!  END OF PROPERTY & PRINTOUT BLOCK
! ========================================================================

! ========================================================================
!  SAVE FOR FINAL PRINTOUT
! ========================================================================
   eelec=etotal
!  Etot
   etotal = etotal + ep + exb
   eat = eatoms*evtoau - etotal
   etotal = etotal + ed
   res%e_elec  = eelec
   res%e_atom  = eat
   res%e_rep   = ep
   res%e_es    = ees
   res%e_aes   = eaes
   res%e_axc   = epol
   res%e_disp  = ed
   res%e_total = etotal
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

   deallocate(S,H0,tmp,X,H1,kcnsh,zsh,jab,matlist,matlist2,dpint)
   call deallocate_gbsa(gbsa)
end subroutine scf1

!> Implementation of the GFN2-xTB Hamiltonian and related parametrisations.
subroutine scf2(iunit, mol, wfn, basis, param, pcem, opt, &
      &         egap, etotal, g, sigma, res)

   use mctc_econv, only : autoev,evtoau

   ! type definitions
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data
   use tbdef_timer
   use tbdef_pcem
   use tbdef_neighbourlist

   ! global storage
   use aoparam
   use setparam

   ! interfaces
   use scc_core
   use grad_core
   use aespot
   use gbobc
   use tbmod_dftd4, only: build_wdispmat, d4dim, d4, disppot, &
      &                   mdisp, prmolc6, edisp_scc, edisp, abcappr, &
      &                   d4_gradient, dispgrad
   use ncoord, only: get_coordination_number, tb_cn_type
   use embedding
   use readin
   use lidep

   !> Unit identifier for all IO, only used for prlevel > 0.
   integer, intent(in) :: iunit
   !> Molecular structure information.
   type(tb_molecule), intent(in) :: mol
   !> Tight binding wavefunction, contains previous wavefunction for restart
   !  and is updated after this run.
   type(tb_wavefunction), intent(inout) :: wfn
   !> Atomic orbital basis set.
   type(tb_basisset), intent(in) :: basis
   !> Global method parameters.
   type(scc_parameter), intent(in) :: param
   !> External potentials, like point charge fields.
   type(tb_pcem), intent(inout) :: pcem
   !> Options for the SCF.
   type(scf_options), intent(in) :: opt
   !> HOMO-LUMO gap in eV.
   real(wp), intent(inout) :: egap
   !> Total energy in Eh.
   real(wp), intent(inout) :: etotal
   !> Gradient in Eh/a0.
   real(wp), intent(inout) :: g(:, :)
   !> Detailed results on the calculation.
   type(scc_results), intent(out) :: res
   !> Stress tensor in Eh (volume contribution not included).
   real(wp), intent(inout) :: sigma(:, :)

! ========================================================================
   type(tb_neighbourlist) :: neighlist
   type(tb_neighlist_generator) :: neighgen
   integer, allocatable :: neighs(:)
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: X(:,:)
   real(wp),allocatable :: S(:,:)
   real(wp),allocatable :: S12(:,:)
   real(wp),allocatable :: H0(:)
   real(wp),allocatable :: H1(:)
   real(wp),allocatable :: H(:,:)
   real(wp),allocatable :: jab(:,:)
   real(wp),allocatable :: ves(:) ! shell ES potential
   real(wp),allocatable :: tmp(:)
   real(wp),allocatable :: zsh(:)
   real(wp),allocatable :: qq(:)
   real(wp),allocatable :: qlmom(:,:)
   real(wp),allocatable :: cm5(:)
   real(wp),allocatable :: kcnao(:)
   real(wp),allocatable :: Xcao(:,:)
   real(wp), allocatable :: gam2sh(:)
!  AES stuff
   real(wp),allocatable  :: dpint(:,:,:),qpint(:,:,:)
   real(wp),allocatable  :: gab3(:),gab5(:)
   real(wp),allocatable  :: vs(:),vq(:,:),vd(:,:)
   real(wp),allocatable  :: gam3sh(:)
   real(wp),allocatable  :: radcn(:) ! CBNEW
   real(wp),allocatable  :: draesdcn(:)

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
   real(wp),allocatable :: dqdr(:, :, :), dqdL(:, :, :)
   integer, allocatable :: neighs3(:)
   real(wp) :: etmp

! ========================================================================
   integer,allocatable :: mdlst(:,:),mqlst(:,:)
   integer :: ndp,nqp

   integer,allocatable :: matlist (:,:)
   integer,allocatable :: matlist2(:,:)
   integer,allocatable :: xblist(:,:)
   real(wp),allocatable :: sqrab(:)
   real(wp),allocatable :: dHdcn(:)
   real(wp),allocatable :: dcndr(:,:,:)
   real(wp),allocatable :: dcndL(:,:,:)

   real(wp) :: dipol(3),dip,gsolv,eat,hlgap,efix,eelec
   real(wp) :: temp,xsum,eh1,rab,eold,dum,xx,r2
   real(wp) :: t0,t1,sum,rr,hav,alpha,ep,dx,dy,dz,dum1,r0i,r0j
   real(wp) :: efa,efb,nfoda,nfodb,hdii,hdjj,qconv,ff
   real(wp) :: x1,x2,ed,ga,gb,ehb,h0s,hmat,rab2
   real(wp) :: h0sr,scfconv,rmsq,dum2,drfdxyz(3),yy,tex,rav,tab,ljexp
   real(wp) :: ees,xa,xb,ya,yb,za,zb,repab,jmpol,kdampxb,vvxb,exb
   real(wp) :: w2,xj,gi,gj,rexp,eatoms,kexp,gnorm
   real(wp) :: w0,w1,damp,damp0,gtmp(3),exc
   real(wp) :: d3atm,esave
   real(wp) :: eaes,t6,t7,epol,kexpe
   type(scf_thresholds) :: thr

!  some parameter defaults which are not fitted
   data kexp /1.50_wp/ ! rep exp for exp(-R**kexp)
   data rexp /1.00_wp/ ! rep exp in 1/R**rexp
   data d3atm/0.00_wp/ ! ATM scal, zero in GFN1
   data ljexp/12.0_wp/ ! XB parameter for damped LJ in GFN1

   integer :: ich ! file handle
   integer :: npr,ii,jj,kk,i,j,k,m,iat,jat,mi,jter,atj,kkk,mj,mm
   integer :: ishell,jshell,np,ia,ndimv,l,nmat,nmat2
   integer :: ll,i1,i2,nn,ati,nxb,lin,startpdiag,lladr(4)
   integer :: is,js,lladr2(0:3),tmmetal
   data    lladr  /1,3,6,10/
   data    lladr2 /1,3,5,7/

   character(len=2),external :: asym
   character(len=128) :: atmp,ftmp
   logical :: ex,minpr,pr,fulldiag,xbond,fail,early3d

!  GBSA stuff
   type(tb_solvent) :: gbsa
   real(wp),allocatable :: fgb(:,:)
   real(wp),allocatable :: fhb(:)
   real(wp) :: gborn, ghb
!  for the CM5 charges
   real(wp),allocatable :: cm5a(:)
   real(wp),allocatable :: dcm5a(:,:,:)
   real(wp),allocatable :: fgba(:),dcm5(:,:)
   real(wp) :: hbpow
   real(wp),allocatable :: vborn(:)

!  point charge embedding stuff
   logical  :: lpcem
   real(wp),allocatable :: vpc(:)
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
   if (profile) call timer%new(8,.false.)
   if (profile) call timer%measure(1,"SCC setup")
   rmsq  =1.e+42_wp
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
   etotal  = 0.0_wp
   ep   = 0.0_wp
   ed   = 0.0_wp
   embd = 0.0_wp
   eat  = 0.0_wp
   egap = 0.0_wp
   molpol = 0.0_wp

   pr   = opt%prlevel.gt.1
   minpr= opt%prlevel.gt.0

   thr = scf_thresholds(opt%accuracy, 2)

   !> old/new q mixing start
   damp0=0.20
   broy = mol%n > 1
   !> when to start pseudodiag (=1 after one full diag)
   startpdiag=1000 !large number=never

!ccccccccccccccccccc
! note: H is in eV!
!ccccccccccccccccccc

   ! # atom arrays
   allocate(qq(mol%n),qlmom(3,mol%n),cm5(mol%n),sqrab(mol%n*(mol%n+1)/2), &
      &     dcndr(3,mol%n,mol%n),cn(mol%n),dcndL(3,3,mol%n),vborn(mol%n), &
      &     dhdcn(mol%n))

   ! initialize the GBSA module (GBSA works with CM5 charges)
   if(lgbsa) then
      call new_gbsa(gbsa,mol%n,mol%at)
      allocate(fgb(mol%n,mol%n),fhb(mol%n),cm5a(mol%n),dcm5a(3,mol%n,mol%n))
      gborn=0._wp
      gbsa%gsasa=0._wp
      gbsa%ghb=0._wp
      qq=0._wp
      ! initialize the neighbor list
      call update_nnlist_gbsa(gbsa,mol%xyz,.false.)
      ! compute Born radii
      call compute_brad_sasa(gbsa,mol%xyz)
      ! initialize the fgb matrix (dielectric screening of the Coulomb potential)
      call compute_fgb(gbsa,fgb,fhb)
      ! initialize the CM5 charges computation
      cm5=wfn%q
      cm5a=0.d0
      dcm5a=0.d0
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

   ! ldep J potentials (in Eh) for SCC
   allocate(gam2sh(basis%nshell), source=0.0_wp)
   do is = 1, basis%nshell
      iat=basis%ash(is)
      ati=mol%at(iat)
      gam2sh(is) = gam(ati)*(1.0_wp+lpar(basis%lsh(is),ati))
   enddo
   call get_gfn_coulomb_matrix(mol, basis%nshell, basis%ash, gam2sh, gfn_method, &
      &                        opt%cf, opt%lqpc, jab)

!  J potentials including the point charge stuff
   if(lpcem)then
      allocate( vpc(basis%nshell), source = 0.0_wp )
      call jpot_pcem_gfn2(mol%n,pcem,basis%nshell,mol%at,mol%xyz,basis%ash, &
         &                basis%lsh,vpc)
   endif

   ! set 3rd order shell gammas
   allocate(gam3sh(basis%nshell),source = 0.0_wp)
   do is=1,basis%nshell
      iat=basis%ash(is)
      ati=mol%at(iat)
      dum=param%gam3l(basis%lsh(is))  ! sp or d-pol
      if ((basis%lsh(is).eq.2).and.(tmmetal(ati).ge.1)) dum=param%gam3l(3) ! d-val
      gam3sh(is)=gam3(ati)*dum
   enddo

   if (opt%prlevel > 1) then
      write(iunit,'(/,10x,51("."))')
      write(iunit,'(10x,":",22x,a,22x,":")') "SETUP"
      write(iunit,'(10x,":",49("."),":")')
      write(iunit,intfmt) "# basis functions  ",basis%nbf
      write(iunit,intfmt) "# atomic orbitals  ",basis%nao
      write(iunit,intfmt) "# shells           ",basis%nshell
      write(iunit,intfmt) "# electrons        ",wfn%nel
      write(iunit,intfmt) "max. iterations    ",opt%maxiter
      write(iunit,chrfmt) "Hamiltonian        ","GFN2-xTB"
      write(iunit,chrfmt) "GBSA solvation     ",bool2string(lgbsa)
      write(iunit,chrfmt) "PC potential       ",bool2string(lpcem)
      if (lpcem) then
         write(iunit,intfmt) "-> # point charges ",pcem%n
         write(iunit,dblfmt) "-> sum of PC       ",sum(pcem%q),"e   "
      endif
      write(iunit,dblfmt) "electronic temp.   ",opt%etemp, "K   "
      write(iunit,dblfmt) "accuracy           ",opt%accuracy, "    "
      write(iunit,scifmt) "-> integral cutoff ",thr%intcut,  "    "
      write(iunit,scifmt) "-> integral neglect",thr%neglect, "    "
      write(iunit,scifmt) "-> SCF convergence ",thr%scfconv, "Eh  "
      write(iunit,scifmt) "-> wf. convergence ",thr%qconv,   "e   "
      write(iunit,dblfmt) "Broyden damping    ",broydamp,"    "
      write(iunit,'(10x,51("."))')
   endif

   qq    =wfn%q
   damp  =damp0

   if (profile) call timer%measure(1)
   if (profile) call timer%measure(8,"Neighbourlist")

   !> setup neighbourlist
   call neighgen%new(maxval(thr), len(mol), mol%xyz, mol%lattice, mol%npbc > 0)
   call neighlist%new(neighgen)
   allocate(neighs(len(mol)), source=0)

   if (profile) call timer%measure(8)
   if (profile) call timer%measure(2,"Dispersion")

   call neighlist%get_neighs(neighs, cutoff=thr%cutoff_cn)
   call d4dim(mol%n,mol%at,dispdim)
   allocate( gw(dispdim), &
      &      c6abns(dispdim,dispdim), &
      &      wdispmat(dispdim,dispdim), &
      &      hdisp(mol%n), &
      &      source=0.0_wp )
   call get_coordination_number(mol, neighs, neighlist, tb_cn_type%cov, &
      &                         cn, dcndr, dcndL)
   call d4(mol%n,dispdim,mol%at,param%wf,param%g_a,param%g_c,cn,gw,c6abns)
   call build_wdispmat(mol%n,dispdim,mol%at,mol%xyz,param%disp,c6abns,gw, &
      &                wdispmat)

   ! prepare aes stuff
   ! CN/dCN replaced by special smoother and faster decaying function
   call get_coordination_number(mol, neighs, neighlist, tb_cn_type%gfn, &
      &                         cn, dcndr, dcndL)

   if (profile) call timer%measure(2)
   if (profile) call timer%measure(3,"integral evaluation")

   allocate(dpint(3,basis%nao,basis%nao), &
      &     qpint(6,basis%nao,basis%nao), &
      &     source = 0.0_wp)
   call neighlist%get_neighs(neighs, cutoff=thr%cutoff_ham)
   ! compute integrals and prescreen to set up list arrays
   call build_SDQH0(mol, neighs, neighlist, basis, param, thr%intcut, cn, &
      &             S, dpint, qpint, H0)
   call count_multipole_ints(basis%nao, thr%neglect, ndp, dpint, nqp, qpint)
   call cholesky(iunit,pr,basis%nao,S,orthog) ! S is not modified
   if (orthog) then
      call prmat(6, S, basis%nao, basis%nao, 'overlap')
      call raise('E', 'System is linear dependent', 1)
   endif

   ! allocate arrays for lists and fill (to exploit sparsity)
   allocate(mdlst(2,ndp),mqlst(2,nqp))
   call setdqlist(basis%nao,ndp,nqp,thr%neglect,dpint,qpint,mdlst,mqlst)
   ! set up 1/R^n * damping function terms
   ii=mol%n*(mol%n+1)/2
   allocate(gab3(ii),gab5(ii),radcn(mol%n))
   call get_radcn(mol%n,mol%at,cn,param%cn_shift,param%cn_expo,param%cn_rmax,radcn)
   ! zero damping, xbrad=kdmp3,xbdamp=kdmp5
   call mmomgabzero(mol%n,mol%at,mol%xyz,param%xbrad,param%xbdamp,radcn,gab3,gab5)

   if (profile) call timer%measure(3)
   if (profile) call timer%measure(1)

   ! set first electrostatic potential (start with third order on-side)
   ves = wfn%qsh**2 * gam3sh * autoev
   ! add Born shifts to potential
   if (lgbsa) then
      cm5 = wfn%q + cm5a
      call get_born_shift(mol%n, cm5, fgb, fhb, vborn)
      ves = ves + vborn(basis%ash)
   else
      vborn = 0.0_wp
   endif
   ! add dispersion potential
   call disppot(mol%n,dispdim,mol%at,wfn%q,param%g_a,param%g_c,wdispmat,gw, &
      &         hdisp)
   ves = ves + hdisp(basis%ash)*autoev
   ! add external potential
   if(lpcem) then
      ves = ves + vpc
   endif
   ! get second order electrostatic potential
   call get_charge_shift(basis%nshell,wfn%qsh,jab,ves)

   ! compute intermediates for potential
   allocate(vs(mol%n),vd(3,mol%n),vq(6,mol%n), source=0.0_wp)
   call setvsdq(mol%n,mol%at,mol%xyz,wfn%q,wfn%dipm,wfn%qp,gab3,gab5,vs,vd,vq)

   ! prepare matrix indices
   nmat =0
   nmat2=0
   do ii=1,basis%nao
      iat=basis%aoat2(ii)
      do jj=1,ii-1
         jat=basis%aoat2(jj)
         if(abs(S(jj,ii)).lt.thr%neglect) then
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
      kcnao(ii)=kcnat(ishell-1,mol%at(iat))
   enddo

   if (profile) call timer%measure(1)
   if (profile) call timer%measure(6,"classical contributions")

   call neighlist%get_neighs(neighs, cutoff=thr%cutoff_rep)
   call rep_grad_gfn2(mol, neighs, neighlist, rexp, ep, g, sigma)

   if (profile) call timer%measure(6)
   if (profile) call timer%measure(4,"zeroth order Hamiltonian")

   if(pr)then
      write(iunit,'(a)')
      write(iunit,*) 'iter      E             dE          RMSdq', &
      &'      gap      omega  full diag'
   endif

   ! first order energy for given geom. and density, i.e. skip SCC and grad
   if(opt%maxiter.eq.0) then
      call qsh2qat(mol%n,mol%at,basis%nshell,wfn%qsh,wfn%q)
      call electro2(mol%n,mol%at,basis%nao,basis%nshell,jab,H0,wfn%P, &
         &          wfn%q,gam3sh,wfn%qsh,ees,etotal)
      if(lgbsa) then
         cm5=wfn%q+cm5a
         call electro_gbsa(mol%n,fgb,fhb,cm5,gborn,etotal)
      endif
      goto 9999
   endif
   if (profile) call timer%measure(4)

! ========================================================================
!  SCC iterations
! ========================================================================
   if (profile) call timer%measure(5,"iterations")
   call scc_gfn2(iunit,mol%n,wfn%nel,wfn%nopen,basis%nao,ndp,nqp,nmat,basis%nshell, &
      &          mol%at,matlist,mdlst,mqlst,basis%aoat2,basis%ao2sh,basis%ash, &
      &          wfn%q,wfn%dipm,wfn%qp,qq,qlmom,wfn%qsh,zsh, &
      &          mol%xyz,vs,vd,vq,gab3,gab5, &
      &          gbsa,fgb,fhb,cm5,cm5a,gborn,vborn, &
      &          newdisp,dispdim,param%g_a,param%g_c,gw,wdispmat,hdisp, &
      &          broy,broydamp,damp0, &
      &          lpcem,ves,vpc, &
      &          opt%etemp,wfn%focc,wfn%focca,wfn%foccb,wfn%efa,wfn%efb, &
      &          etotal,ees,eaes,epol,ed,epcem,egap, &
      &          wfn%emo,wfn%ihomo,wfn%ihomoa,wfn%ihomob, &
      &          H0,H1,wfn%C,S,dpint,qpint,X,wfn%P,jab,gam3sh, &
      &          opt%maxiter,startpdiag,thr%scfconv,thr%qconv, &
      &          minpr,pr, &
      &          fail,jter)
! ========================================================================
   ! free some memory (this stuff is not needed for gradients)
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
   !  get energy weighted density matrix
   tmp = wfn%focc * wfn%emo*evtoau
   call dmat(basis%nao,tmp,wfn%C,X)

   allocate(H(basis%nao, basis%nao), source=0.0_wp)
   !  preccalc
   do m=1,nmat2
      i=matlist2(1,m)
      j=matlist2(2,m)
      kk=j+i*(i-1)/2
      ishell = basis%ao2sh(i)
      jshell = basis%ao2sh(j)
      eh1 = -0.5_wp*(ves(ishell)+ves(jshell))
      H(j,i)=(eh1+H0(kk)/S(j,i))*evtoau*wfn%P(j,i)-X(j,i)
      H(i,j)=H(j,i)
   enddo

   ! multipole gradient stuff
   ! VS, VD, VQ-dependent potentials are changed w.r.t. SCF,
   ! since moment integrals are now computed with origin at
   ! respective atoms
   call setdvsdq(mol%n,mol%at,mol%xyz,wfn%q,wfn%dipm,wfn%qp,gab3,gab5,vs,vd,vq)
   call neighlist%get_neighs(neighs, cutoff=thr%cutoff_ham)
   call build_dSDQH0(mol, neighs, neighlist, basis, param, thr%intcut, cn, &
      &              wfn%P, H, ves, vs, vd, vq, & ! FIXME
      &              dhdcn, g, sigma)
   ! CN-level shift gradient
   call dgemv('n', 3*mol%n, mol%n,-1.0_wp, dcndr, 3*mol%n, dhdcn, 1, 1.0_wp, g, 1)
   call dgemv('n', 9, mol%n, 1.0_wp, dcndL, 9, dhdcn, 1, 1.0_wp, sigma, 1)

   allocate(draesdcn(len(mol)), source=0.0_wp)
   call dradcn(mol%n,mol%at,cn,param%cn_shift,param%cn_expo,param%cn_rmax,draesdcn)
   call aniso_grad(mol%n,mol%at,mol%xyz,wfn%q,wfn%dipm,wfn%qp,param%xbrad,param%xbdamp, &
      &          radcn,draesdcn,dcndr,gab3,gab5,g)

   ! dispersion (DFT-D type correction)
   !allocate(dqdr(3, mol%n, mol%n), dqdL(3, 3, mol%n), source=0.0_wp)
   !allocate(neighs3(mol%n), source=0)
   !call neighlist%get_neighs(neighs, cutoff=thr%cutoff_cn)
   !call get_coordination_number(mol, neighs, neighlist, tb_cn_type%cov, &
   !   &                         cn, dcndr, dcndL)
   !call neighlist%get_neighs(neighs, cutoff=thr%cutoff_disp)
   !call neighlist%get_neighs(neighs3, cutoff=thr%cutoff_atm)
   !call d4_gradient(mol, neighs, neighs3, neighlist, param%disp, param%g_a, &
   !   &             param%g_c, param%wf, cn, dcndr, dcndL, wfn%q, dqdr, dqdL, &
   !   &             etmp, g, sigma)
   !embd = etmp - ed
   call dispgrad(mol%n,dispdim,mol%at,wfn%q,mol%xyz, &
      &        param%disp,param%wf,param%g_a,param%g_c, &
      &        c6abns,mbd,g,embd)
   embd = embd-ed

   ! GBSA
   ! start GBSA gradient
   if (lgbsa) then
      call compute_gb_egrad(gbsa,wfn%q,gborn,ghb,g,minpr)
      ! solvation energy
      gbsa%gborn = gborn
      gbsa%ghb = ghb
      gsolv=gborn+gbsa%gsasa+gbsa%ghb+gshift
      etotal=etotal+gbsa%gsasa+gshift
   endif

   call get_gfn_coulomb_derivs(mol, basis%nshell, basis%ash, gam2sh, gfn_method, &
      &                        opt%cf, opt%lqpc, wfn%qsh, g, sigma)

   ! --- ES point charge embedding
   if (lpcem) then
      call pcem_grad_gfn2(g,pcem%grd,mol%n,pcem,mol%at,basis%nshell,mol%xyz, &
         &                basis%ash,basis%lsh,wfn%qsh)
   endif

!  calculate the norm for printout
   gnorm = sqrt(sum( g**2 ))
! ========================================================================
!  clear some space
   deallocate(radcn,gab3,gab5,mdlst,mqlst,vs,vd,vq,gam3sh) !CBNEW

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
!     HOMO-LUMO excitation properties if  UHF=2        
      if (wfn%nopen.eq.2) then
         call hlex(mol%n,mol%at,basis%nbf,basis%nao,wfn%ihomoa,mol%xyz,wfn%focc,S,wfn%C,wfn%emo,basis)
      endif

! ------------------------------------------------------------------------
!     LMO /xTB-IFF
      if (pr_lmo) then
         tmp=wfn%emo*evtoau
         call local(mol%n,mol%at,basis%nbf,basis%nao,wfn%ihomoa,mol%xyz,mol%z,wfn%focc,S,wfn%P,wfn%C,tmp,wfn%q,etotal,lgbsa,basis)
      endif

! ------------------------------------------------------------------------
!  exchange energy correction ala sTDA
      if (wfn%nopen.ge.2) then
         call exch(mol%n,mol%at,basis%nao,wfn%nopen,wfn%ihomoa,mol%xyz,wfn%focc,S,wfn%C,exc,basis%aoat)
         write(iunit,'(''open-shell EX :'',F16.7)') -exc
         write(iunit,'(''corrected Etot:'',F16.7, &
         &   '' (not used further except for this printout!)'')') etotal - exc
      endif

   endif printing

! ------------------------------------------------------------------------
!  get Wiberg bond orders
   call get_wiberg(mol%n,basis%nao,mol%at,mol%xyz,wfn%P,S,wfn%wbo,basis%fila2)

! ------------------------------------------------------------------------
!  dipole calculation (always done because its free)
   call calc_dipole(mol%n,mol%at,mol%xyz,mol%z,basis%nao,wfn%P,dpint,dip,dipol)

   if (profile) call timer%measure(7)

!  END OF PROPERTY & PRINTOUT BLOCK
! ========================================================================

! ========================================================================
!  SAVE FOR FINAL PRINTOUT
! ========================================================================
   eelec=etotal
!  Etot
   etotal = etotal + ep + exb + embd
   eat = eatoms*evtoau - etotal
   res%e_elec  = eelec
   res%e_atom  = eat
   res%e_rep   = ep
   res%e_es    = ees
   res%e_aes   = eaes
   res%e_axc   = epol
   res%e_disp  = ed+embd
   res%e_total = etotal
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
end subroutine scf2


subroutine count_multipole_ints(nao, thr, ndp, dpint, nqp, qpint)
   integer, intent(in) :: nao
   integer, intent(out) :: ndp
   integer, intent(out) :: nqp
   real(wp), intent(in) :: thr
   real(wp), intent(in) :: dpint(:, :, :)
   real(wp), intent(in) :: qpint(:, :, :)
   integer :: i, j, ij, k, kk
   real(wp) :: tmp1, tmp2, thr2
   thr2=(thr*1.0d-2)-thr*1.0d-12
   ndp = 0
   do i=1,nao
      do j=1,i
         ij=ij+1
         tmp1=0.0_wp
         kk=0
         do k=1,3
            tmp1=tmp1+dpint(k,j,i)*dpint(k,j,i)
         enddo
         if(tmp1.gt.thr2) ndp=ndp+1
      enddo
   enddo
   nqp = 0
   do i=1,nao
      do j=1,i
         ij=ij+1
         tmp2=0.0_wp
         kk=0
         do k=1,3
            tmp2=tmp2-qpint(k,j,i)*qpint(k,j,i)
         enddo
         do k=1,6
            tmp2=tmp2+2.0_wp*qpint(k,j,i)*qpint(k,j,i)
         enddo
         if(tmp2.gt.thr2) nqp=nqp+1
      enddo
   enddo
end subroutine count_multipole_ints

real(wp) elemental function scfthr_maxval(self) result(cutoff)
   type(scf_thresholds), intent(in) :: self
   cutoff = max(self%cutoff_disp, &
      &         self%cutoff_cn, &
      &         self%cutoff_ham, &
      &         self%cutoff_atm, &
      &         self%cutoff_rep)
end function scfthr_maxval

!> Generate a set of numerical thresholds from a given accuracy value.
type(scf_thresholds) elemental function new_thresholds(accuracy, method) result(thr)
   real(wp), intent(in) :: accuracy
   integer, intent(in) :: method
   !> primitive cut-off
   thr%intcut = max(20.0_wp, 25.0_wp - 10.0_wp*log10(accuracy))
   !> integral neglect threshold
   thr%neglect = 10.0e-9_wp * accuracy
   !> exit if E < scfconv, Hessian sensitive to this
   thr%scfconv = 1.0e-6_wp * accuracy
   !> conv threshold in dRMS q
   if (method > 1) then
      thr%qconv = 1.e-4_wp * accuracy
   else
      thr%qconv = 2.e-5_wp * accuracy
   end if
   thr%cutoff_disp = 64.0_wp - log10(accuracy)*30.0_wp
   thr%cutoff_cn = 40.0_wp - log10(accuracy)*20.0_wp
   thr%cutoff_ham = 30.0_wp - log10(accuracy)*10.0_wp
   thr%cutoff_rep = 40.0_wp - log10(accuracy)*20.0_wp
   thr%cutoff_atm = 40.0_wp - log10(accuracy)*20.0_wp
end function new_thresholds

end module scf_module
