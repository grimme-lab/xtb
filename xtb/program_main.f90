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

program XTBprog
!! ========================================================================
!  we use the official FORTRAN standard environment variables,
!  namely, the unit identifiers and the iostat specification,
!  since FORTRAN2008 there are also kindparameters available.
!  remember: output_unit = 6, input_unit = 5, error_unit = 0
   use iso_fortran_env, istdout => output_unit, istderr => error_unit, wp => real64

!! ========================================================================
!  this is part of the Grimme group's MCTC general purpose FORTRAN libs
   use mctc_timings
   use mctc_systools
   use mctc_econv
   use mctc_param

!! ========================================================================
!  class and type definitions
   use tbdef_molecule
   use tbdef_calculator
   use tbdef_wavefunction
   use tbdef_param
   use tbdef_data

!! ========================================================================
!  former common variable storage, used in the entire xtb code
   use aoparam
   use setparam
   use sphereparam
   use scanparam
   use splitparam
   use symparam
   use fixparam
   use constrain_param

!! ========================================================================
!  global modules with parametes used in the entire xtb code
   use gbobc, only: lgbsa,init_gbsa
   use shake_module, only: init_shake
   use embedding, only : init_pcem

!! ========================================================================
!  input/output modules
   use readin
   use printout
   use argparser
   use set_module
   use property_output
   use tbmod_file_utils

!! ========================================================================
!  get interfaces for methods used in this part
   use scc_core,    only : iniqshell
   use qcextern,    only : orca_chk,mopac_chk
   use single,      only : singlepoint
   use aespot,      only : get_radcn
   use iniq,        only : iniqcn
   use eeq_model
   use ncoord,      only : ncoord_gfn,dncoord_erf,dncoord_d3,ncoord_erf,ncoord_d3
   use xbasis
   use dftd4,       only : d4init,p_refq_gfn2xtb
   use axis_trafo,  only : axis3
   use qmdff,       only : ff_ini
!! ------------------------------------------------------------------------
!  interfaces from internal libraries
   use optimizer, only : wrlog
   use hessian, only : numhess
   use dynamic, only : md
   use modef,   only : modefollow
   use mdoptim, only : mdopt
   use screening, only : screen

!! ------------------------------------------------------------------------
!  periodic boundary conditions
!   use pbc, only : get_realspace_cutoff

   implicit none

   logical, parameter    :: is_release = .false.

!! ========================================================================
!  use some wrapper types to bundle information together
   type(tb_molecule) :: mol
   type(scc_results) :: res
   type(tb_calculator) :: calc
   type(freq_results) :: fres
   type(tb_wavefunction) :: wfn
   type(chrg_parameter) :: chrgeq
!  store important names and stuff like that in FORTRAN strings
   character(len=:),allocatable :: fname    ! geometry input file
   character(len=:),allocatable :: xcontrol ! instruction file
   character(len=:),allocatable :: xrc      ! global instruction file
   character(len=:),allocatable :: fnv      ! parameter file
   character(len=:),allocatable :: fnx      ! NOT USED
   character(len=:),allocatable :: tmpname  ! temporary string
   character(len=:),allocatable :: cdum     ! temporary string
   character(len=:),allocatable :: extension, basename, directory
   integer :: ftype

   integer :: nargs
   type(string),    allocatable :: argument_list(:)

!! ========================================================================
!  default names for important files in xtb
   character(len=*),parameter :: p_fname_rc = '.xtbrc'
   character(len=*),parameter :: p_fname_param_gfn0  = '.param_gfn0.xtb'
   character(len=*),parameter :: p_fname_param_gfn1  = '.param_gfn.xtb'
   character(len=*),parameter :: p_fname_param_gfn2  = '.param_gfn2.xtb'
   character(len=*),parameter :: p_fname_param_ipea  = '.param_ipea.xtb'
   character(len=*),parameter :: p_fname_param_stda1 = '.param_stda1.xtb'
   character(len=*),parameter :: p_fname_param_stda2 = '.param_stda2.xtb'

   integer :: chrg,gsolvstate
   integer :: i,j,k,l,idum
   integer :: ich,ictrl,iprop ! file handle
   real(wp) :: sigma(3,3)
   real(wp),allocatable :: cn  (:)
   real(wp),allocatable :: sat (:)
   real(wp),allocatable :: g   (:,:)
   real(wp) :: globpar(25),vec3(3)
   real(wp),allocatable :: dcn (:,:,:)
   real(wp),allocatable :: dq  (:,:,:)
   real(wp),allocatable :: dumdumdum  (:,:,:)
   real(wp),allocatable :: q  (:)
   real(wp),allocatable :: ql  (:)
   real(wp),allocatable :: qr  (:)

!! ------------------------------------------------------------------------
   character(len=2),external :: asym
   integer,external :: ncore

!! ========================================================================
!  debugging variables for numerical gradient
   logical, parameter    :: gen_param = .false.
   logical, parameter    :: debug = .false.
   type(tb_wavefunction) :: wf0
   real(wp),allocatable  :: coord(:,:),numg(:,:),gdum(:,:)
   real(wp) :: sdum(3,3)
   real(wp),parameter    :: step = 0.00001_wp, step2 = 0.5_wp/step
   real(wp) :: er,el
   logical  :: coffee ! if debugging gets really though, get a coffee

!! ------------------------------------------------------------------------
!  undocumented and unexplainable variables go here
   integer  :: rohf,err
   real(wp) :: dum5,egap,etot
   real(wp) :: zero,t0,t1,w0,w1,acc,etot2,g298
   real(wp) :: qmdff_s6,one,two
   real(wp) :: ea,ip
   real(wp) :: vomega,vfukui
   real(wp),allocatable :: f_plus(:), f_minus(:)
   type(tb_wavefunction) :: wf_p, wf_m
   parameter (zero=0.0_wp)
   parameter (one =1.0_wp)
   parameter (two =2.0_wp)
   logical :: ex,okbas
   logical :: epr,diff,murks
   logical :: exist
   logical :: lgrad,restart
   logical :: copycontrol
   logical :: newreader
   logical :: strict

!  OMP stuff
   integer :: TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
   integer :: nproc

!  start by initializing the MCTC library
   call mctc_init('xtb',10,.true.)

!  we are over this already, comment back in if you plan to break stuff
!  call raise('S','program not (yet) intended for productive runs!',1)

!! ========================================================================
!  get the XTBPATH/XTBHOME variables
   call rdvar('XTBHOME',xenv%home,iostat=err)
   if ((err.ne.0).or.(xenv%home.eq.'')) then
      call raise('S','XTBHOME is not set, using HOME instead',1)
      call rdvar('HOME',xenv%home)
      ! actually we should check if HOME is set by the OS, but ... nay
   endif
   call rdvar('XTBPATH',xenv%path,iostat=err)
   if ((err.ne.0).or.(xenv%path.eq.'')) then
      call raise('S','XTBPATH is not set, using XTBHOME instead',1)
      xenv%path = xenv%home
   endif

!! ========================================================================
!  read the command line arguments
   call rdxargs(fname,xcontrol,fnv,fnx,acc,lgrad,restart,gsolvstate,strict, &
   &            copycontrol,argument_list,nargs,coffee)

!! ========================================================================
!  read the xcontrol file
   call rdcontrol(xcontrol,copy_file=copycontrol)

!! ------------------------------------------------------------------------
!  read dot-Files before reading the rc and after reading the xcontrol
   call open_file(ich,'.CHRG','r')
   if (ich.ne.-1) then
      call getline(ich,cdum,iostat=err)
      if (err.ne.0) call raise('E','.CHRG is empty!',1)
      call set_chrg(cdum)
      call close_file(ich)
   endif

   call open_file(ich,'.UHF','r')
   if (ich.ne.-1) then
      call getline(ich,cdum,iostat=err)
      if (err.ne.0) call raise('E','.UHF is empty!',1)
      call set_spin(cdum)
      call close_file(ich)
   endif

!! ------------------------------------------------------------------------
!  read the xtbrc if you can find it (use rdpath directly instead of xfind)
   call rdpath(xenv%path,p_fname_rc,xrc,exist)
   if (exist) then
      call rdcontrol(xrc,copy_file=.false.)
   endif

!! ========================================================================
!  no user interaction up to now, time to show off!
!  print the xtb banner with version number and compilation date
!  making a fancy version of this is hard, x is difficult in ASCII art
   call xtb_header(istdout)
!  make sure you cannot blame us for destroying your computer
   call disclamer(istdout)
!  how to cite this program
   call citation(istdout)

!! ------------------------------------------------------------------------
!  lets show what we can do, because we can
!  this clutters the screen quite a bit since we started actually
!  documenting our options, so we will only do this if the user set
!  the program in verbose mode
   if (verbose) call help

!! ========================================================================
!  check if you are allowed to make a calculation today
   call prdate('S')
   if (is_release) call expire

!  C6 scaling of QMDFF dispersion to simulate solvent
   qmdff_s6=1.0_wp
   chrg=ichrg
   rohf=1 ! HS default

!! ------------------------------------------------------------------------
!  read molecule, # atoms first
   if (coffee) then ! it's coffee time
      fname = 'caffeine'
      call get_coffee(mol)
   else
      call file_generate_meta_info(fname, extension, basename, directory)
      call file_figure_out_ftype(ftype, extension, basename)
      call open_file(ich, fname, 'r')
      call mol%read(ich, format=ftype)
      call close_file(ich)
   endif

   if(mol%n.lt.1) call raise('E','no atoms!',1)

!  get some memory
   allocate(cn(mol%n),sat(mol%n),g(3,mol%n), source = 0.0_wp)

!! ------------------------------------------------------------------------
   periodic = mol%npbc > 0
   if (mol%npbc == 0) then
      if (do_cma_trafo) then
         allocate(coord(3,mol%n),source=0.0_wp)
         call axis3(1,mol%n,mol%at,mol%xyz,coord,vec3)
         mol%xyz = coord
         deallocate(coord)
      endif
   endif

!  initialize the global storage
   call init_fix(mol%n)
   call init_split(mol%n)
   call init_constr(mol%n,mol%at)
   call init_scan
   call init_walls
   call init_metadyn(mol%n,metaset%maxsave)

   do i=1,mol%n
      mol%z(i) = mol%at(i) - ncore( mol%at(i) )
      ! lanthanides without f are treated as La
      if(mol%at(i).gt.57.and.mol%at(i).lt.72) mol%z(i)=3
      atmass(i)=atomic_mass(mol%at(i)) * autoamu ! from splitparam.f90
   enddo

   ! initialize time step for MD if requested autocomplete
   if (tstep_md.lt.0.0_wp) then
      tstep_md=(minval(atmass)/(atomic_mass(1)*autoamu))**(1.0_wp/3.0_wp)
   endif

   wfn%nel = idint(sum(mol%z)) - chrg
   wfn%nopen = nalphabeta
   if(wfn%nopen.eq.0.and.mod(wfn%nel,2).ne.0) wfn%nopen=1
   call initrand

   call setup_summary(istdout,mol%n,fname,xcontrol,nargs,argument_list,wfn,xrc,exist)

   if(fit) acc=0.2 ! higher SCF accuracy during fit

!! ------------------------------------------------------------------------
!  CONSTRAINTS & SCANS
!! ------------------------------------------------------------------------
!  now we are at a point that we can check for requested constraints,
!  unfortunately, we have to do this separately from xcontrol, since
!  the xcontrol instructions are parsed before we even know that the
!  file exist (not quite true) were we could get a geometry from, which
!  is problematic in the aspect that geometry constraints tend to be
!  geometry dependent.
   call read_userdata(xcontrol,mol%n,mol%at,mol%xyz)

!  initialize metadynamics
   call load_metadynamic(metaset,mol%n,mol%at,mol%xyz)
!  restraining potential
   if (allocated(potset%xyz)) then
      if (lconstr_all_bonds)    call constrain_all_bonds(mol%n,mol%at,potset%xyz)
      if (lconstr_all_angles)   call constrain_all_angles(mol%n,mol%at,potset%xyz)
      if (lconstr_all_torsions) call constrain_all_torsions(mol%n,mol%at,potset%xyz)
      call setup_constrain_pot(mol%n,mol%at,potset%xyz)
   else
      if (lconstr_all_bonds)    call constrain_all_bonds(mol%n,mol%at,mol%xyz)
      if (lconstr_all_angles)   call constrain_all_angles(mol%n,mol%at,mol%xyz)
      if (lconstr_all_torsions) call constrain_all_torsions(mol%n,mol%at,mol%xyz)
      call setup_constrain_pot(mol%n,mol%at,mol%xyz)
   endif
   egap = 0.0_wp
!  fragmentation for CMA constrain
   if(iatf1.eq.0.and.iatf2.eq.0) then
      call ncoord_erf(mol%n,mol%at,mol%xyz,cn)
      call splitm(mol%n,mol%at,mol%xyz,cn)
   endif
   call splitprint(mol%n,mol%at,mol%xyz)

   if (verbose) then
      call fix_info(istdout,mol%n,mol%at,mol%xyz)
      call pot_info(istdout,mol%n,mol%at,mol%xyz)
   endif

!! ------------------------------------------------------------------------
!  write xcontrol
   if (copycontrol) then
      call open_set(ictrl,xcontrol)
      call write_set(ictrl)
      call close_set(ictrl)
   endif

!! ------------------------------------------------------------------------
!  if you have requested a define we stop here...
   if (define) then
      if (verbose) call main_geometry(istdout,mol)
      call eval_define(veryverbose)
   endif
!  in case you want to be productive instead of meddling around with
!  define, you wish to know if there is something gone wrong anyway
   call raise('F','Please study the warnings concerning your '// &
   &              'input carefully',1)

!  You had it coming!
   if (strict) call mctc_strict

!! ------------------------------------------------------------------------
!  PARAMETER
!! ------------------------------------------------------------------------
   if (gfn_method.eq.3) & ! lets set one thing straight:
   &  call raise('E','This is an internal error, please use gfn_method=2!',1)
   if (.not.allocated(fnv)) then
      select case(runtyp)
      case default
         call raise('E','This is an internal error, please define your runtypes!',1)
      case(p_run_nox,p_run_stda)
         fnv=xfind(p_fname_param_stda1)
         call stda_header(istdout)
      case(p_run_scc,p_run_grad,p_run_opt,p_run_hess,p_run_ohess, &
           p_run_md,p_run_omd,p_run_siman,p_run_path,p_run_screen, &
           p_run_gmd,p_run_modef,p_run_mdopt,p_run_metaopt,p_run_reactor)
         if(gfn_method.eq.0) then
            fnv=xfind(p_fname_param_gfn0)
            if(periodic)then
               call peeq_header(istdout)
            else
               call gfn0_header(istdout)
            endif
         endif
         if(gfn_method.eq.1) then
            fnv=xfind(p_fname_param_gfn1)
            call gfn1_header(istdout)
         endif
         if(gfn_method.eq.2) then
            fnv=xfind(p_fname_param_gfn2)
            call gfn2_header(istdout)
         endif
      case(p_run_vip,p_run_vea,p_run_vipea,p_run_vfukui,p_run_vomega)
         if(gfn_method.eq.0) then
            fnv=xfind(p_fname_param_gfn0)
            call gfn0_header(istdout)
         endif
         if(gfn_method.eq.1) then
            fnv=xfind(p_fname_param_ipea)
            call ipea_header(istdout)
         endif
         if(gfn_method.eq.2) then
            fnv=xfind(p_fname_param_gfn2)
            call gfn2_header(istdout)
         endif
      end select
   endif

   call open_file(ich,fnv,'r')
   exist = ich .ne. -1
   if (exist) then
      call read_gfn_param(ich,globpar,.true.)
      call close_file(ich)
   else ! no parameter file, check if we have one compiled into the code
      call use_parameterset(fnv,globpar,exist)
      if (.not.exist) call raise('E','parameter file '//fnv//' not found!',1)
   endif
   do i = 1, 86
      do j = 1, i
         if (abs(kpair(j,i)-1.0_wp).gt.1e-5_wp) &
            write(istdout,'(13x,"KAB for ",a2," - ",a2,5x,":",F22.4)') &
            asym(j),asym(i),kpair(j,i)
      enddo
   enddo


   if (gen_param) then
   !  generate a warning to keep release versions from generating huge files
      call raise('S','XTB IS DUMPING PARAMETERFILES, RESET GEN_PARAM FOR RELEASE!',1)
      call prelemparam(globpar)
   endif

   allocate(calc%param)

   if (runtyp.gt.1) then
      select case(gfn_method)
      case default
         call raise('E','Internal error, wrong GFN method passed!',1)
      case(1)
         call set_gfn1_parameter(calc%param,globpar)
         call gfn1_prparam(istdout,mol%n,mol%at,calc%param)
      case(2)
         call set_gfn2_parameter(calc%param,globpar,mol%n,mol%at)
         call gfn2_prparam(istdout,mol%n,mol%at,calc%param)
      case(0)
         call set_gfn0_parameter(calc%param,globpar,mol%n,mol%at)
         call gfn0_prparam(istdout,mol%n,mol%at,calc%param)
      end select
   else
      calc%param%kspd(1:6)=globpar(1:6)
      calc%param%wllscal  =globpar(7)
      calc%param%gscal    =globpar(8)
      calc%param%zcnf     =globpar(9)
      calc%param%tfac     =globpar(10)
      calc%param%kcn      =globpar(11)
      calc%param%fpol     =globpar(12)
      calc%param%ken1     =globpar(13)
      calc%param%lshift   =globpar(14)
      calc%param%lshifta  =globpar(15)
      calc%param%split    =globpar(16)
      calc%param%zqf      =globpar(17)
      write(istdout,'(5x,''method parameters'')')
      write(istdout,'(1x,''k(s)        :'',F8.4)') calc%param%kspd(1)
      write(istdout,'(1x,''k(p)        :'',F8.4)') calc%param%kspd(2)
      write(istdout,'(1x,''k(d)        :'',F8.4)') calc%param%kspd(3)
      write(istdout,'(1x,''k(f)        :'',F8.4)') calc%param%kspd(4)
      write(istdout,'(1x,''Tscal       :'',F8.4)') calc%param%tfac
      write(istdout,'(1x,''Gscal       :'',F8.4)') calc%param%gscal
      write(istdout,'(1x,''fpol        :'',F8.4)') calc%param%fpol
      write(istdout,'(1x,''Zcnf        :'',F8.4)') calc%param%zcnf
      write(istdout,'(1x,''Zqf         :'',F8.4)') calc%param%zqf
      write(istdout,'(1x,''kcn         :'',F8.4)') calc%param%kcn
      write(istdout,'(1x,''kEN1        :'',F8.4)') calc%param%ken1
      write(istdout,'(1x,''wllscal     :'',F8.4)') calc%param%wllscal
   endif
   write(istdout,'(a)')

!  unrestricted atomic spin constants
!  optimized scale factor for Mulliken spin densities
   call setwll_pbe (calc%param%wllscal)

!  init GBSA part
   if(lgbsa) then
      call init_gbsa(istdout,solvent,gsolvstate,temp_md,gfn_method,ngrida)
   endif
!  initialize PC embedding (set default file names and stuff)
   call init_pcem

!! ------------------------------------------------------------------------
!! ------------------------------------------------------------------------
!
!                    calcs start here
!
!! ------------------------------------------------------------------------
!! ------------------------------------------------------------------------

!! ========================================================================
   if(periodic.and.gfn_method.ne.0)then
      call raise('E', 'Periodic implementation only available at zeroth-order (GFN0).',1)
   end if

! ======================================================================
!  set up the basis set for the tb-Hamiltonian
! ======================================================================
   allocate(calc%basis)
   call xbasis0(mol%n,mol%at,calc%basis)
   select case(gfn_method)
   case default
      call raise('E','Internal error, wrong GFN method passed!',1)
   case(p_method_gfn1xtb)
      call xbasis_gfn1(mol%n,mol%at,calc%basis,okbas,diff)
   case(p_method_gfn2xtb)
      call xbasis_gfn2(mol%n,mol%at,calc%basis,okbas)
   case(p_method_gfn0xtb)
      call xbasis_gfn0(mol%n,mol%at,calc%basis,okbas,diff)
   end select
   if (.not.okbas) call raise('E','TB basis incomplete',1)
   call xbasis_cao2sao(mol%n,mol%at,calc%basis)

! ======================================================================
!  initial guess, setup wavefunction
! ======================================================================
   call wfn%allocate(mol%n,calc%basis%nshell,calc%basis%nao)

!  EN charges and CN
   if (runtyp.gt.1 .and. gfn_method.gt.0) then
      ! GFN case, CN is just for printout and will be overwritten in SCC
      if (gfn_method.lt.2) then
         call ncoord_d3(mol%n,mol%at,mol%xyz,cn)
      else
         call ncoord_gfn(mol%n,mol%at,mol%xyz,cn)
      endif
      if (guess_charges.eq.p_guess_gasteiger) then
         call iniqcn(mol%n,wfn%nel,mol%at,mol%z,mol%xyz,chrg,calc%param%ken1,wfn%q,cn,gfn_method,.true.)
      else if (guess_charges.eq.p_guess_goedecker) then
         call ncoord_erf(mol%n,mol%at,mol%xyz,cn)
         call goedecker_chrgeq(mol%n,mol%at,mol%xyz,real(chrg,wp),cn,dcn,wfn%q,dq,er,g,&
                               .false.,.false.,.false.)
      else
         call ncoord_gfn(mol%n,mol%at,mol%xyz,cn)
         wfn%q = real(chrg,wp)/real(mol%n,wp)
      endif
   else if (gfn_method.eq.0) then
      continue
   else
      ! vTB case
      call iniqcn(mol%n,wfn%nel,mol%at,mol%z,mol%xyz,chrg,calc%param%ken1,wfn%q,cn,0,.true.)
   endif
!  initialize shell charges from gasteiger charges
   call iniqshell(mol%n,mol%at,mol%z,calc%basis%nshell,wfn%q,wfn%qsh,gfn_method)

!! ========================================================================
!                         S C C  section
!! ========================================================================
   select case(mode_extrun)
   case default
      call scc_header(istdout)
   case(p_ext_eht)
      call eht_header(istdout)
   case(p_ext_qmdff)
      call qmdff_header(istdout)
      call ff_ini(mol%n,mol%at,mol%xyz,cn,qmdff_s6)
   case(p_ext_orca)
      call driver_header(istdout)
      call orca_chk
   case(p_ext_mopac)
      call driver_header(istdout)
      call mopac_chk
   case(p_ext_turbomole)
      call driver_header(istdout)
      write(*,*) 'Running Turbomole Calculation!'
   end select

   call delete_file('.sccnotconverged')

   if (restart.and.mode_extrun.eq.p_ext_xtb) then ! only in first run
      call read_restart(wfn,'xtbrestart',mol%n,mol%at,gfn_method,exist,.true.)
   endif

!  the SP energy which is always done
   call start_timing(2)
   call singlepoint &
   &       (istdout,mol,wfn,calc, &
   &        egap,etemp,maxscciter,2,exist,lgrad,acc,etot,g,sigma,res)
   call stop_timing(2)

!! ========================================================================
!  numerical gradient for debugging purposes
!! ========================================================================
   if (debug) then
!  generate a warning to keep release versions from calculating numerical gradients
   call raise('S','XTB IS CALCULATING NUMERICAL GRADIENTS, RESET DEBUG FOR RELEASE!',1)
   print'(/,"analytical gradient")'
   print *, g
   allocate( coord(3,mol%n), source = mol%xyz )
   allocate( numg(3,mol%n),gdum(3,mol%n), source = 0.0_wp )
   wf0 = wfn
   do i = 1, mol%n
      do j = 1, 3
         mol%xyz(j,i) = mol%xyz(j,i) + step
         wfn = wf0
         call singlepoint &
         &       (istdout,mol,wfn,calc, &
         &        egap,etemp,maxscciter,0,.true.,.true.,acc,er,gdum,sdum,res)
         mol%xyz(j,i) = mol%xyz(j,i) - 2*step
         wfn = wf0
         call singlepoint &
         &       (istdout,mol,wfn,calc, &
         &        egap,etemp,maxscciter,0,.true.,.true.,acc,el,gdum,sdum,res)
         mol%xyz(j,i) = mol%xyz(j,i) + step
         numg(j,i) = step2 * (er - el)
      enddo
   enddo
   print'(/,"numerical gradient")'
   print *, numg
   print'(/,"difference gradient")'
   print*,g-numg
   endif

!! ========================================================================
!  implement all further calculations here
!! ========================================================================

!! ------------------------------------------------------------------------
!  ANCopt
   if ((runtyp.eq.p_run_opt).or.(runtyp.eq.p_run_ohess).or. &
   &   (runtyp.eq.p_run_omd).or.(runtyp.eq.p_run_screen).or. &
   &   (runtyp.eq.p_run_metaopt)) then
      if (opt_engine.eq.p_engine_rf) &
      call ancopt_header(istdout,veryverbose)
      call start_timing(3)
      call geometry_optimization &
      &     (mol,wfn,calc, &
      &      egap,etemp,maxscciter,optset%maxoptcycle,etot,g,sigma,optset%optlev,.true.,.false.,murks)
      res%e_total = etot
      res%gnorm = norm2(g)
      if (nscan.gt.0) then
         call relaxed_scan(mol,wfn,calc)
      endif
      call stop_timing(3)
   endif

!! ------------------------------------------------------------------------
!  automatic VIP and VEA single point (maybe after optimization)
   if (runtyp.eq.p_run_vip.or.runtyp.eq.p_run_vipea &
      & .or.runtyp.eq.p_run_vomega) then
      call start_timing(2)
      call vip_header(istdout)
      wfn%nel = wfn%nel-1
      if (mod(wfn%nel,2).ne.0) wfn%nopen = 1
      call singlepoint &
      &       (istdout,mol,wfn,calc, &
      &        egap,etemp,maxscciter,2,.true.,.false.,acc,etot2,g,sigma,res)
      ip=etot2-etot-calc%param%ipshift
      write(istdout,'(72("-"))')
      write(istdout,'("empirical IP shift (eV):",f10.4)') &
      &                  autoev*calc%param%ipshift
      write(istdout,'("delta SCC IP (eV):",f10.4)') autoev*ip
      write(istdout,'(72("-"))')
      wfn%nel = wfn%nel+1
      call stop_timing(2)
   endif

   if (runtyp.eq.p_run_vea.or.runtyp.eq.p_run_vipea &
      & .or.runtyp.eq.p_run_vomega) then
      call start_timing(2)
      call vea_header(istdout)
      wfn%nel = wfn%nel+1
      if (mod(wfn%nel,2).ne.0) wfn%nopen = 1
      call singlepoint &
      &       (istdout,mol,wfn,calc, &
      &        egap,etemp,maxscciter,2,.true.,.false.,acc,etot2,g,sigma,res)
      ea=etot-etot2-calc%param%eashift
      write(istdout,'(72("-"))')
      write(istdout,'("empirical EA shift (eV):",f10.4)') &
      &                  autoev*calc%param%eashift
      write(istdout,'("delta SCC EA (eV):",f10.4)') autoev*ea
      write(istdout,'(72("-"))')

      wfn%nel = wfn%nel-1
      call stop_timing(2)
   endif

!  vomega (electrophilicity) index
   if (runtyp.eq.p_run_vomega) then
      write(istdout,'(a)')
      write(istdout,'(72("-"))')
      write(istdout,'(a,1x,a)') &
         "Calculation of global electrophilicity index",&
         "(IP+EA)²/(8·(IP-EA))"
      vomega=(ip+ea)**2/(8*(ip-ea))
      write(istdout,'("Global electrophilicity index (eV):",f10.4)') &
        autoev*vomega
      write(istdout,'(72("-"))')
   endif

!  Fukui Index from Mulliken population analysis
   if (runtyp.eq.p_run_vfukui) then
     allocate(f_plus(mol%n),f_minus(mol%n))
     write(istdout,'(a)')
     write(istdout,'("Fukui index Calculation")')
     wf_p=wfn
     wf_m=wfn
     wf_p%nel = wf_p%nel+1
     if (mod(wf_p%nel,2).ne.0) wf_p%nopen = 1
     call singlepoint &
     &       (istdout,mol,wf_p,calc, &
     &        egap,etemp,maxscciter,1,.true.,.false.,acc,etot2,g,sigma,res)
     f_plus=wf_p%q-wfn%q

     wf_m%nel = wf_m%nel-1
     if (mod(wf_m%nel,2).ne.0) wf_m%nopen = 1
     call singlepoint &
     &       (istdout,mol,wf_m,calc, &
     &        egap,etemp,maxscciter,1,.true.,.false.,acc,etot2,g,sigma,res)
     f_minus=wfn%q-wf_m%q
     write(istdout,'(a)')
     write(istdout, '(1x,"    #       f(+)     f(-)     f(0)")')
     do i=1,mol%n
       write(istdout,'(i6,a3,2f9.3,2f9.3,2f9.3)') i, asym(mol%at(i)), f_plus(i), f_minus(i), 0.5d0*(wf_p%q(i)-wf_m%q(i))
     enddo
     deallocate(f_plus,f_minus)
   endif

!! ------------------------------------------------------------------------
!  numerical hessian calculation
   if ((runtyp.eq.p_run_hess).or.(runtyp.eq.p_run_ohess)) then
      call numhess_header(istdout)
      if (mol%npbc > 0) then
         call raise('E',"Phonon calculations under PBC are not implemented",1)
      endif
      call start_timing(5)
      call numhess &
      &       (mol,wfn,calc, &
      &        egap,etemp,maxscciter,etot,g,sigma,fres)
      call stop_timing(5)
   endif

   ! reset the gap, since it is currently not updated in ancopt and numhess
   res%hl_gap = wfn%emo(wfn%ihomo+1)-wfn%emo(wfn%ihomo)

!! ========================================================================
!  PRINTOUT SECTION
!! ========================================================================
   if (allocated(property_file)) then
      call open_file(iprop,property_file,'w')
      if (iprop.eq.-1) then
         iprop = istdout
         deallocate(property_file)
      else
         write(istdout,'(/,a)') "Property printout bound to '"//property_file//"'"
         if (allocated(cdum)) deallocate(cdum)
         call get_command(length=l)
         allocate( character(len=l) :: cdum )
         call get_command(cdum)
         write(iprop,'("command:  ''",a,"''")') cdum
         call rdvar('HOSTNAME',cdum,err)
         if (err.eq.0) &
         write(iprop,'("hostname: ''",a,"''")') cdum
         write(iprop,'("date:     ",a)') prtimestring('S')
      endif
   else
      iprop = istdout
   endif

   call generic_header(iprop,'Property Printout',49,10)
   if (lgrad) call tmgrad(mol%n,mol%at,mol%xyz,g,etot)

   if(periodic)then
      write(*,*)'Periodic properties'
   else
      call main_property(iprop, &
           mol,wfn,calc%basis,calc%param,res,acc)
      call main_cube(verbose, &
           mol,wfn,calc%basis,calc%param,res)
   endif


   if (pr_json) then
      call open_file(ich,'xtbout.json','w')
      call main_json(ich, &
         mol,wfn,calc%basis,calc%param,res,fres)
      call close_file(ich)
   endif

   if ((runtyp.eq.p_run_opt).or.(runtyp.eq.p_run_ohess).or. &
       (runtyp.eq.p_run_omd).or.(runtyp.eq.p_run_screen).or. &
       (runtyp.eq.p_run_metaopt)) then
      call main_geometry(iprop,mol)
   endif

   if ((runtyp.eq.p_run_hess).or.(runtyp.eq.p_run_ohess)) then
      call generic_header(iprop,'Frequency Printout',49,10)
      call main_freq(iprop,mol,wfn,fres)
   endif

   if (allocated(property_file)) then
      if (iprop.ne.-1 .and. iprop.ne.istdout) then
         call write_energy(iprop,res,fres, &
            & (runtyp.eq.p_run_hess).or.(runtyp.eq.p_run_ohess))
         call close_file(iprop)
      endif
   endif

   if ((runtyp.eq.p_run_opt).or.(runtyp.eq.p_run_ohess).or. &
       (runtyp.eq.p_run_omd).or.(runtyp.eq.p_run_screen).or. &
       (runtyp.eq.p_run_metaopt)) then
      call file_generate_name(tmpname, 'xtbopt', extension, mol%ftype)
      write(istdout,'(/,a,1x,a,/)') &
         "optimized geometry written to:",tmpname
      call open_file(ich,tmpname,'w')
      call mol%write(ich, energy=res%e_total, gnorm=res%gnorm)
      call close_file(ich)
   endif

   call write_energy(istdout,res,fres, &
      & (runtyp.eq.p_run_hess).or.(runtyp.eq.p_run_ohess))

!! ------------------------------------------------------------------------
!  xtb molecular dynamics
   if ((runtyp.eq.p_run_md).or.(runtyp.eq.p_run_omd)) then
      if (metaset%maxsave .gt. 0) then
         if (mol%npbc > 0) then
            call raise('E',"Metadynamic under PBC is not implemented",1)
         endif
         call metadyn_header(istdout)
      else
         call md_header(istdout)
      endif
      fixset%n = 0 ! no fixing for MD runs
      call start_timing(6)
      idum = 0
      if (shake_md) call init_shake(mol%n,mol%at,mol%xyz,wfn%wbo)
      call md &
      &     (mol,wfn,calc, &
      &      egap,etemp,maxscciter,etot,g,sigma,0,temp_md,idum)
      call stop_timing(6)
   endif

!! ------------------------------------------------------------------------
!  metadynamics
   if (runtyp.eq.p_run_metaopt) then
      if (mol%npbc > 0) then
         call raise('S',"Metadynamic under PBC is not implemented",1)
      endif
      call metadyn_header(istdout)
      ! check if ANCOPT already convered
      if (murks) then
         call raise('E','Optimization did not converge, aborting',1)
      endif
      write(istdout,'(1x,"output written to xtbmeta.log")')
      call open_file(ich,'xtbmeta.log','w')
      call wrlog(ich,mol%n,mol%xyz,mol%at,etot,norm2(g),.true.)
      k = metaset%nstruc+1
      call start_timing(6)
      do l = k, metaset%maxsave
         metaset%nstruc = l
         metaset%xyz(:,:,metaset%nstruc) = mol%xyz
         ! randomize structure to avoid zero RMSD
         do i = 1, mol%n
            do j = 1, 3
               call random_number(er)
               mol%xyz(j,i) = mol%xyz(j,i) + 1.0e-6_wp*er
            enddo
         enddo
         call geometry_optimization &
         &     (mol,wfn,calc, &
         &      egap,etemp,maxscciter,optset%maxoptcycle,etot,g,sigma,optset%optlev,verbose,.true.,murks)
         if (.not.verbose) then
            write(istdout,'("current energy:",1x,f20.8)') etot
         endif
         if (murks) then
            call close_file(ich)
            write(istdout,'(/,3x,"***",1x,a,1x,"***",/)') &
               "FAILED TO CONVERGE GEOMETRY OPTIMIZATION"
            call touch_file('NOT_CONVERGED')
         endif
         call wrlog(ich,mol%n,mol%xyz,mol%at,etot,norm2(g),.false.)
      enddo
      call close_file(ich)
      call stop_timing(6)
   endif


!! ------------------------------------------------------------------------
!  simulated annealing
   if (runtyp.eq.p_run_siman) then
      call siman_header(istdout)
      call raise('E',"SIMAN has been deprecated in favor of CREST")
   endif

!! ------------------------------------------------------------------------
!  path finder
   if (runtyp.eq.p_run_path) then
      call rmsdpath_header(istdout)
      if (mol%npbc > 0) then
         call raise('S',"Metadynamics under PBC are not implemented",1)
      endif
      call start_timing(4)
      call bias_path(mol,wfn,calc,egap,etemp,maxscciter,etot,g,sigma)
      call stop_timing(4)
   endif

!! ------------------------------------------------------------------------
!  screen over input structures
   if (runtyp.eq.p_run_screen) then
      call start_timing(8)
      call screen(mol,wfn,calc,egap,etemp,maxscciter,etot,g,sigma)
      call stop_timing(8)
   endif

!! ------------------------------------------------------------------------
!  gmd for averaged interaction energies
   if (runtyp.eq.p_run_gmd) then
      call raise('E','GMD option is not supported anymore',1)
   endif

!! ------------------------------------------------------------------------
!  mode following for conformer search
   if (runtyp.eq.p_run_modef) then
      if (mol%npbc > 0) then
         call raise('S',"Modefollowing under PBC is not implemented",1)
      endif
      call start_timing(9)
      call modefollow(mol,wfn,calc,egap,etemp,maxscciter,etot,g,sigma)
      call stop_timing(9)
   endif

!! ------------------------------------------------------------------------
!  optimize along MD from xtb.trj for conformer searches
   if (runtyp.eq.p_run_mdopt) then
      call start_timing(10)
      call mdopt(mol,wfn,calc,egap,etemp,maxscciter,etot,g,sigma)
      call stop_timing(10)
   endif

!! ------------------------------------------------------------------------
   if (runtyp.eq.p_run_reactor) then
      call reactor_header(istdout)
      call raise('E','The nano-reactor has been moved!',1)
   endif

!! ========================================================================
!  to further speed up xtb calculations we dump our most important
!  quantities in a restart file, so we can save some precious seconds
   if (restart) then
      call write_restart(wfn,'xtbrestart',gfn_method)
   endif
   call wfn%deallocate

!! ========================================================================
!  to accelerate the fitting procedure we can print out various
!  properties when in fit-modus. So you have to return any property
!  you want to use back to the mainlevel and implement the printout here

!! ========================================================================
!  we may have generated some non-fatal errors, which have been saved,
!  so we should tell the user, (s)he may want to know what went wrong
   call raise('F','Some non-fatal runtime exceptions were caught,'// &
   &              ' please check:',1)

!! ========================================================================
!  print all files xtb interacted with while running (for debugging mainly)
   if (verbose) then
   write(istdout,'(a)')
   write(istdout,'(72("-"))')
   call print_filelist(istdout)
   endif

!! ========================================================================
!  make some post processing afterward, show some timings and stuff
   write(istdout,'(a)')
   write(istdout,'(72("-"))')
   call stop_timing_run
   call stop_timing(1)
   call prdate('E')
   write(istdout,'(72("-"))')
   call prtiming(1,'total')
   call prtiming(2,'SCF')
   if ((runtyp.eq.p_run_opt).or.(runtyp.eq.p_run_ohess).or. &
   &   (runtyp.eq.p_run_omd).or.(runtyp.eq.p_run_metaopt)) then
      call prtiming(3,'ANC optimizer')
   endif
   if ((runtyp.eq.p_run_hess).or.(runtyp.eq.p_run_ohess)) then
      call prtiming(5,'numerical hessian')
   endif
   if ((runtyp.eq.p_run_md).or.(runtyp.eq.p_run_omd).or. &
       (runtyp.eq.p_run_metaopt)) then
      call prtiming(6,'MD')
   endif
   if (runtyp.eq.p_run_siman) then
      call prtiming(7,'SIMAN')
   endif
   if (runtyp.eq.p_run_screen) then
      call prtiming(8,'screen')
   endif
   if (runtyp.eq.p_run_modef) then
      call prtiming(9,'mode following')
   endif
   if (runtyp.eq.p_run_mdopt) then
      call prtiming(10,'MD opt.')
   endif

   write(istdout,'(a)')
   if (coffee) then
      write(istderr,'(a)') "nicely done, now get yourself a"
      call get_cup(istderr,'(5x,a)')
   else
      call terminate(0)
   endif

end program XTBprog

