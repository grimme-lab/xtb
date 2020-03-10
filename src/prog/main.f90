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

program XTBprog
!! ========================================================================
!  we use the official FORTRAN standard environment variables,
!  namely, the unit identifiers and the iostat specification,
!  since FORTRAN2008 there are also kindparameters available.
!  remember: output_unit = 6, input_unit = 5, error_unit = 0
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stderr
!! ========================================================================
!  this is part of the Grimme group's MCTC general purpose FORTRAN libs
   use xtb_mctc_timings
   use xtb_mctc_systools
   use xtb_mctc_convert
   use xtb_mctc_param

!! ========================================================================
!  class and type definitions
   use xtb_type_molecule
   use xtb_type_calculator
   use xtb_type_wavefunction
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_type_latticepoint, only : TLatticePoint, init
   use xtb_type_wignerseitzcell, only : TWignerSeitzCell, init

!! ========================================================================
!  former common variable storage, used in the entire xtb code
   use xtb_aoparam
   use xtb_setparam
   use xtb_sphereparam
   use xtb_scanparam
   use xtb_splitparam
   use xtb_symparam
   use xtb_fixparam
   use xtb_constrain_param

!! ========================================================================
!  global modules with parametes used in the entire xtb code
   use xtb_solv_gbobc, only: lgbsa,init_gbsa
   use xtb_shake, only: init_shake
   use xtb_embedding, only : init_pcem

!! ========================================================================
!  input/output modules
   use xtb_io_reader, only : readMolecule
   use xtb_io_writer, only : writeMolecule
   use xtb_mctc_filetypes, only : fileType, getFileType, generateFileMetaInfo, &
      & generateFileName
   use xtb_readin
   use xtb_printout
   use xtb_argparser
   use xtb_setmod
   use xtb_propertyoutput
   use tbmod_output_writer
   use xtb_restart
   use xtb_readparam

!! ========================================================================
!  get interfaces for methods used in this part
   use xtb_scc_core, only : iniqshell
   use xtb_extern_orca, only : checkOrca
   use xtb_extern_mopac, only : checkMopac
   use xtb_single, only : singlepoint
   use xtb_aespot, only : get_radcn
   use xtb_iniq, only : iniqcn
   use xtb_eeq
   use xtb_disp_ncoord, only : ncoord_gfn, dncoord_erf, dncoord_d3, ncoord_erf, &
      & ncoord_d3
   use xtb_basis
   use xtb_axis, only : axis3
   use xtb_qmdff, only : ff_ini
!! ------------------------------------------------------------------------
!  interfaces from internal libraries
   use xtb_hessian, only : numhess
   use xtb_dynamic, only : md
   use xtb_modef, only : modefollow
   use xtb_mdoptim, only : mdopt
   use xtb_screening, only : screen

!! ------------------------------------------------------------------------
!  periodic boundary conditions
!   use xtb_pbc, only : get_realspace_cutoff

   implicit none

   !> Source of errors in the main program unit
   character(len=*), parameter :: source = "prog_main"

   logical, parameter    :: is_release = .false.

!! ========================================================================
!  use some wrapper types to bundle information together
   type(TMolecule) :: mol
   type(scc_results) :: res
   type(TCalculator) :: calc
   type(freq_results) :: fres
   type(TWavefunction) :: wfn
   type(chrg_parameter) :: chrgeq
   type(TEnvironment) :: env
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
   real(wp),allocatable :: latticePoint(:,:)

!! ------------------------------------------------------------------------
   integer,external :: ncore

!! ========================================================================
!  debugging variables for numerical gradient
   logical, parameter    :: gen_param = .false.
   logical, parameter    :: debug = .false.
   type(TWavefunction) :: wf0
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
   type(TWavefunction) :: wf_p, wf_m
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
   logical :: exitRun

!  OMP stuff
   integer :: TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
   integer :: nproc

!  start by initializing the MCTC library
   call mctc_init('xtb',10,.true.)
   call init(env)

!! ========================================================================
!  get the XTBPATH/XTBHOME variables
   call rdvar('XTBHOME',xenv%home,iostat=err)
   if ((err.ne.0).or.(xenv%home.eq.'')) then
      call env%warning('XTBHOME is not set, using HOME instead', source)
      call rdvar('HOME',xenv%home)
      ! actually we should check if HOME is set by the OS, but ... nay
   endif
   call rdvar('XTBPATH',xenv%path,iostat=err)
   if ((err.ne.0).or.(xenv%path.eq.'')) then
      call env%warning('XTBPATH is not set, using XTBHOME instead', source)
      xenv%path = xenv%home
   endif

!! ========================================================================
!  read the command line arguments
   call rdxargs(env, fname,xcontrol,fnv,fnx,acc,lgrad,restart,gsolvstate,strict, &
      & copycontrol,argument_list,nargs,coffee)

   call env%checkpoint("Command line argument parsing failed")
   

!! ========================================================================
!  read the xcontrol file
   call rdcontrol(xcontrol,env,copy_file=copycontrol)

   call env%checkpoint("Reading '"//xcontrol//"' failed")

!! ------------------------------------------------------------------------
!  read dot-Files before reading the rc and after reading the xcontrol
   call open_file(ich,'.CHRG','r')
   if (ich.ne.-1) then
      call getline(ich,cdum,iostat=err)
      if (err /= 0) then
         call env%error('.CHRG is empty!', source)
      else
         call set_chrg(env,cdum)
         call close_file(ich)
      end if
   end if

   call env%checkpoint("Reading charge from file failed")

   call open_file(ich,'.UHF','r')
   if (ich.ne.-1) then
      call getline(ich,cdum,iostat=err)
      if (err /= 0) then
         call env%error('.UHF is empty!', source)
      else
         call set_spin(env,cdum)
         call close_file(ich)
      end if
   endif

   call env%checkpoint("Reading multiplicity from file failed")

!! ------------------------------------------------------------------------
!  read the xtbrc if you can find it (use rdpath directly instead of xfind)
   call rdpath(xenv%path,p_fname_rc,xrc,exist)
   if (exist) then
      call rdcontrol(xrc,env,copy_file=.false.)

      call env%checkpoint("Reading '"//xrc//"' failed")
   endif

   ! make sure that we get a eht calculation instead of a scc for GFN0
   if(gfn_method.eq.0)  call set_exttyp('eht')

!! ========================================================================
!  no user interaction up to now, time to show off!
!  print the xtb banner with version number and compilation date
!  making a fancy version of this is hard, x is difficult in ASCII art
   call xtb_header(env%unit)
!  make sure you cannot blame us for destroying your computer
   call disclamer(env%unit)
!  how to cite this program
   call citation(env%unit)

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
      call generateFileMetaInfo(fname, directory, basename, extension)
   else                                                              
      call generateFileMetaInfo(fname, directory, basename, extension)
      ftype = getFileType(basename, extension)
      call open_file(ich, fname, 'r')
      call readMolecule(env, mol, ich, ftype)
      call close_file(ich)

      call env%checkpoint("reading geometry input '"//fname//"' failed")
   endif

   if(mol%n.lt.1) call env%terminate('System has no atoms')

!  get some memory
   allocate(cn(mol%n),sat(mol%n),g(3,mol%n), source = 0.0_wp)

!  initialize the global storage
   call init_fix(mol%n)
   call init_split(mol%n)
   call init_constr(mol%n,mol%at)
   call init_scan
   call init_walls
   call init_metadyn(mol%n,metaset%maxsave)

!! ------------------------------------------------------------------------
   atmass = atomic_mass(mol%at) * autoamu ! from splitparam.f90
   periodic = mol%npbc > 0
   if (mol%npbc == 0) then
      if (do_cma_trafo) then
         allocate(coord(3,mol%n),source=0.0_wp)
         call axis3(1,mol%n,mol%at,mol%xyz,coord,vec3)
         mol%xyz = coord
         deallocate(coord)
      endif
   endif

   do i=1,mol%n
      mol%z(i) = mol%at(i) - ncore( mol%at(i) )
      ! lanthanides without f are treated as La
      if(mol%at(i).gt.57.and.mol%at(i).lt.72) mol%z(i)=3
   enddo

   ! initialize time step for MD if requested autocomplete
   if (tstep_md.lt.0.0_wp) then
      tstep_md=(minval(atmass)/(atomic_mass(1)*autoamu))**(1.0_wp/3.0_wp)
   endif

   wfn%nel = idint(sum(mol%z)) - chrg
   wfn%nopen = nalphabeta
   if(wfn%nopen.eq.0.and.mod(wfn%nel,2).ne.0) wfn%nopen=1
   call initrand

   call setup_summary(env%unit,mol%n,fname,xcontrol,nargs,argument_list,wfn,xrc,exist)

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
   call read_userdata(xcontrol,env,mol)

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
      call fix_info(env%unit,mol%n,mol%at,mol%xyz)
      call pot_info(env%unit,mol%n,mol%at,mol%xyz)
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
      if (verbose) call main_geometry(env%unit,mol)
      call eval_define(veryverbose)
   endif
!  in case you want to be productive instead of meddling around with
!  define, you wish to know if there is something gone wrong anyway
   call env%show('Please study the warnings concerning your input carefully')
   call raise('F', 'Please study the warnings concerning your input carefully', 1)

!  You had it coming!
   if (strict) call mctc_strict

   call check_cold_fusion(mol)

   !> Finalize geometry setup
   allocate(calc%latp)
   allocate(calc%neighList)
   allocate(calc%wsCell)

   !> Setup lattice points
   write(env%unit, '(" * Setup lattice points",t40,"...")', advance='no')
   call init(calc%latp, env, mol, 60.0_wp)  ! Fixed cutoff
   call calc%latp%getLatticePoints(latticePoint, 40.0_wp)

   write(env%unit, '(i10,1x,"points")') calc%latp%nTrans

   call env%checkpoint("Generation of lattice points failed")

   !> Setup neighbour list
   write(env%unit, '(" * Setup neighbour list",t40,"...")', advance='no')
   call init(calc%neighList, len(mol))
   call calc%neighList%generate(env, mol%xyz, 40.0_wp, latticePoint, .false.)

   write(env%unit, '(i10,1x,"images")') size(calc%neighList%image)

   call env%checkpoint("Generation of neighbour list failed")

   !> Setup Wigner-Seitz cell
   write(env%unit, '(" * Setup Wigner-Seitz cell",t40,"...")', advance='no')
   call init(calc%wsCell, len(mol))
   call calc%wsCell%generate(env, mol%xyz, 40.0_wp, latticePoint, .false.)

   write(env%unit, '(i10,1x,"images")') size(calc%wsCell%image)

   call env%checkpoint("Generation of Wigner-Seitz cell failed")

!! ------------------------------------------------------------------------
!  PARAMETER
!! ------------------------------------------------------------------------
   if (gfn_method.eq.3) & ! lets set one thing straight:
   &  call env%terminate('This is an internal error, please use gfn_method=2!')
   if (.not.allocated(fnv)) then
      select case(runtyp)
      case default
         call env%terminate('This is an internal error, please define your runtypes!')
      case(p_run_nox,p_run_stda)
         fnv=xfind(p_fname_param_stda1)
         call stda_header(env%unit)
      case(p_run_scc,p_run_grad,p_run_opt,p_run_hess,p_run_ohess, &
           p_run_md,p_run_omd,p_run_siman,p_run_path,p_run_screen, &
           p_run_gmd,p_run_modef,p_run_mdopt,p_run_metaopt,p_run_reactor)
         if(gfn_method.eq.0) then
            fnv=xfind(p_fname_param_gfn0)
            if(periodic)then
               call peeq_header(env%unit)
            else
               call gfn0_header(env%unit)
            endif
         endif
         if(gfn_method.eq.1) then
            fnv=xfind(p_fname_param_gfn1)
            call gfn1_header(env%unit)
         endif
         if(gfn_method.eq.2) then
            fnv=xfind(p_fname_param_gfn2)
            call gfn2_header(env%unit)
         endif
      case(p_run_vip,p_run_vea,p_run_vipea,p_run_vfukui,p_run_vomega)
         if(gfn_method.eq.0) then
            fnv=xfind(p_fname_param_gfn0)
            call gfn0_header(env%unit)
         endif
         if(gfn_method.eq.1) then
            fnv=xfind(p_fname_param_ipea)
            call ipea_header(env%unit)
         endif
         if(gfn_method.eq.2) then
            fnv=xfind(p_fname_param_gfn2)
            call gfn2_header(env%unit)
         endif
      end select
   endif

   call open_file(ich,fnv,'r')
   exist = ich .ne. -1
   if (exist) then
      call readParam(env,ich,globpar,.true.)
      call close_file(ich)
   else ! no parameter file, check if we have one compiled into the code
      call use_parameterset(fnv,globpar,exist)
      if (.not.exist) then
         call env%error('Parameter file '//fnv//' not found!', source)
      end if
   endif

   call env%checkpoint("Could not setup parameterisation")

   do i = 1, 86
      do j = 1, i
         if (abs(kpair(j,i)-1.0_wp).gt.1e-5_wp) &
            write(env%unit,'(13x,"KAB for ",a2," - ",a2,5x,":",F22.4)') &
            toSymbol(j),toSymbol(i),kpair(j,i)
      enddo
   enddo


   if (gen_param) then
   !  generate a warning to keep release versions from generating huge files
      call env%warning('XTB IS DUMPING PARAMETERFILES, RESET GEN_PARAM FOR RELEASE!')
      call prelemparam(globpar)
   endif

   allocate(calc%param)

   if (runtyp.gt.1) then
      select case(gfn_method)
      case default
         call env%terminate('Internal error, wrong GFN method passed!')
      case(1)
         call set_gfn1_parameter(calc%param,globpar)
         call gfn1_prparam(env%unit,mol%n,mol%at,calc%param)
      case(2)
         call set_gfn2_parameter(calc%param,globpar)
         call gfn2_prparam(env%unit,mol%n,mol%at,calc%param)
      case(0)
         call set_gfn0_parameter(calc%param,globpar)
         call gfn0_prparam(env%unit,mol%n,mol%at,calc%param)
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
      write(env%unit,'(5x,''method parameters'')')
      write(env%unit,'(1x,''k(s)        :'',F8.4)') calc%param%kspd(1)
      write(env%unit,'(1x,''k(p)        :'',F8.4)') calc%param%kspd(2)
      write(env%unit,'(1x,''k(d)        :'',F8.4)') calc%param%kspd(3)
      write(env%unit,'(1x,''k(f)        :'',F8.4)') calc%param%kspd(4)
      write(env%unit,'(1x,''Tscal       :'',F8.4)') calc%param%tfac
      write(env%unit,'(1x,''Gscal       :'',F8.4)') calc%param%gscal
      write(env%unit,'(1x,''fpol        :'',F8.4)') calc%param%fpol
      write(env%unit,'(1x,''Zcnf        :'',F8.4)') calc%param%zcnf
      write(env%unit,'(1x,''Zqf         :'',F8.4)') calc%param%zqf
      write(env%unit,'(1x,''kcn         :'',F8.4)') calc%param%kcn
      write(env%unit,'(1x,''kEN1        :'',F8.4)') calc%param%ken1
      write(env%unit,'(1x,''wllscal     :'',F8.4)') calc%param%wllscal
   endif
   write(env%unit,'(a)')

!  unrestricted atomic spin constants
!  optimized scale factor for Mulliken spin densities
   call setwll_pbe (calc%param%wllscal)

!  init GBSA part
   if(lgbsa) then
      call init_gbsa(env%unit,solvent,gsolvstate,temp_md,gfn_method,ngrida)
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
      call env%terminate('Periodic implementation only available at zeroth-order.')
   end if

! ======================================================================
!  set up the basis set for the tb-Hamiltonian
! ======================================================================
   allocate(calc%basis)
   call xbasis0(mol%n,mol%at,calc%basis)
   select case(gfn_method)
   case default
      call env%terminate('Internal error, wrong GFN method passed!')
   case(p_method_gfn1xtb)
      call xbasis_gfn1(mol%n,mol%at,calc%basis,okbas,diff)
   case(p_method_gfn2xtb)
      call xbasis_gfn2(mol%n,mol%at,calc%basis,okbas)
   case(p_method_gfn0xtb)
      call xbasis_gfn0(mol%n,mol%at,calc%basis,okbas,diff)
   end select
   if (.not.okbas) call env%terminate('basis set could not be setup completely')

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
      call scc_header(env%unit)
   case(p_ext_eht)
      call eht_header(env%unit)
   case(p_ext_qmdff)
      call qmdff_header(env%unit)
      call ff_ini(mol%n,mol%at,mol%xyz,cn,qmdff_s6)
   case(p_ext_orca)
      call driver_header(env%unit)
      call checkOrca(env)
   case(p_ext_mopac)
      call driver_header(env%unit)
      call checkMopac(env)
   case(p_ext_turbomole)
      call driver_header(env%unit)
      write(*,*) 'Running Turbomole Calculation!'
   end select

   call delete_file('.sccnotconverged')

   if (restart.and.mode_extrun.eq.p_ext_xtb) then ! only in first run
      call readRestart(env,wfn,'xtbrestart',mol%n,mol%at,gfn_method,exist,.true.)
   endif

   call env%checkpoint("Setup for calculation failed")

!  the SP energy which is always done
   call start_timing(2)
   call singlepoint &
   &       (env,mol,wfn,calc, &
   &        egap,etemp,maxscciter,2,exist,lgrad,acc,etot,g,sigma,res)
   call stop_timing(2)

   call env%checkpoint("Single point calculation terminated")

!! ========================================================================
!  numerical gradient for debugging purposes
!! ========================================================================
   if (debug) then
!  generate a warning to keep release versions from calculating numerical gradients
   call env%warning('XTB IS CALCULATING NUMERICAL GRADIENTS, RESET DEBUG FOR RELEASE!')
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
         &       (env,mol,wfn,calc, &
         &        egap,etemp,maxscciter,0,.true.,.true.,acc,er,gdum,sdum,res)
         mol%xyz(j,i) = mol%xyz(j,i) - 2*step
         wfn = wf0
         call singlepoint &
         &       (env,mol,wfn,calc, &
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
      call ancopt_header(env%unit,veryverbose)
      call start_timing(3)
      call geometry_optimization &
      &     (env, mol,wfn,calc, &
      &      egap,etemp,maxscciter,optset%maxoptcycle,etot,g,sigma,optset%optlev,.true.,.false.,murks)
      res%e_total = etot
      res%gnorm = norm2(g)
      if (nscan.gt.0) then
         call relaxed_scan(env,mol,wfn,calc)
      endif
      call stop_timing(3)
   endif

!! ------------------------------------------------------------------------
!  automatic VIP and VEA single point (maybe after optimization)
   if (runtyp.eq.p_run_vip.or.runtyp.eq.p_run_vipea &
      & .or.runtyp.eq.p_run_vomega) then
      call start_timing(2)
      call vip_header(env%unit)
      wfn%nel = wfn%nel-1
      if (mod(wfn%nel,2).ne.0) wfn%nopen = 1
      call singlepoint &
      &       (env,mol,wfn,calc, &
      &        egap,etemp,maxscciter,2,.true.,.false.,acc,etot2,g,sigma,res)
      ip=etot2-etot-calc%param%ipshift
      write(env%unit,'(72("-"))')
      write(env%unit,'("empirical IP shift (eV):",f10.4)') &
      &                  autoev*calc%param%ipshift
      write(env%unit,'("delta SCC IP (eV):",f10.4)') autoev*ip
      write(env%unit,'(72("-"))')
      wfn%nel = wfn%nel+1
      call stop_timing(2)
   endif

   if (runtyp.eq.p_run_vea.or.runtyp.eq.p_run_vipea &
      & .or.runtyp.eq.p_run_vomega) then
      call start_timing(2)
      call vea_header(env%unit)
      wfn%nel = wfn%nel+1
      if (mod(wfn%nel,2).ne.0) wfn%nopen = 1
      call singlepoint &
      &       (env,mol,wfn,calc, &
      &        egap,etemp,maxscciter,2,.true.,.false.,acc,etot2,g,sigma,res)
      ea=etot-etot2-calc%param%eashift
      write(env%unit,'(72("-"))')
      write(env%unit,'("empirical EA shift (eV):",f10.4)') &
      &                  autoev*calc%param%eashift
      write(env%unit,'("delta SCC EA (eV):",f10.4)') autoev*ea
      write(env%unit,'(72("-"))')

      wfn%nel = wfn%nel-1
      call stop_timing(2)
   endif

!  vomega (electrophilicity) index
   if (runtyp.eq.p_run_vomega) then
      write(env%unit,'(a)')
      write(env%unit,'(72("-"))')
      write(env%unit,'(a,1x,a)') &
         "Calculation of global electrophilicity index",&
         "(IP+EA)²/(8·(IP-EA))"
      vomega=(ip+ea)**2/(8*(ip-ea))
      write(env%unit,'("Global electrophilicity index (eV):",f10.4)') &
        autoev*vomega
      write(env%unit,'(72("-"))')
   endif

!  Fukui Index from Mulliken population analysis
   if (runtyp.eq.p_run_vfukui) then
     allocate(f_plus(mol%n),f_minus(mol%n))
     write(env%unit,'(a)')
     write(env%unit,'("Fukui index Calculation")')
     wf_p=wfn
     wf_m=wfn
     wf_p%nel = wf_p%nel+1
     if (mod(wf_p%nel,2).ne.0) wf_p%nopen = 1
     call singlepoint &
     &       (env,mol,wf_p,calc, &
     &        egap,etemp,maxscciter,1,.true.,.false.,acc,etot2,g,sigma,res)
     f_plus=wf_p%q-wfn%q

     wf_m%nel = wf_m%nel-1
     if (mod(wf_m%nel,2).ne.0) wf_m%nopen = 1
     call singlepoint &
     &       (env,mol,wf_m,calc, &
     &        egap,etemp,maxscciter,1,.true.,.false.,acc,etot2,g,sigma,res)
     f_minus=wfn%q-wf_m%q
     write(env%unit,'(a)')
     write(env%unit, '(1x,"    #        f(+)     f(-)     f(0)")')
     do i=1,mol%n
       write(env%unit,'(i6,a4,2f9.3,2f9.3,2f9.3)') i, mol%sym(i), f_plus(i), f_minus(i), 0.5d0*(wf_p%q(i)-wf_m%q(i))
     enddo
     deallocate(f_plus,f_minus)
   endif

!! ------------------------------------------------------------------------
!  numerical hessian calculation
   if ((runtyp.eq.p_run_hess).or.(runtyp.eq.p_run_ohess)) then
      call numhess_header(env%unit)
      if (mol%npbc > 0) then
         call env%error("Phonon calculations under PBC are not implemented", source)
      endif
      call start_timing(5)
      call numhess &
      &       (env,mol,wfn,calc, &
      &        egap,etemp,maxscciter,etot,g,sigma,fres)
      call stop_timing(5)

      call env%checkpoint("Hessian calculation terminated")
   endif

   ! reset the gap, since it is currently not updated in ancopt and numhess
   res%hl_gap = wfn%emo(wfn%ihomo+1)-wfn%emo(wfn%ihomo)


   call env%checkpoint("Calculation terminated")

!! ========================================================================
!  PRINTOUT SECTION
!! ========================================================================
   if (allocated(property_file)) then
      call open_file(iprop,property_file,'w')
      if (iprop.eq.-1) then
         iprop = env%unit
         deallocate(property_file)
      else
         write(env%unit,'(/,a)') "Property printout bound to '"//property_file//"'"
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
      iprop = env%unit
   endif

   call generic_header(iprop,'Property Printout',49,10)
   if (lgrad) then
      call write_turbomole(mol, energy=etot, gradient=g, sigma=sigma)
   end if
   if (mol%ftype .eq. fileType%gaussian) then
      if (allocated(basename)) then
         cdum = basename // '.EOu'
      else
         cdum = 'xtb-gaussian.EOu'
      end if
      call open_file(ich, cdum, 'w')
      call write_gaussian_eou(ich, etot, res%dipole, g)
      call close_file(ich)
   end if

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
      if (iprop.ne.-1 .and. iprop.ne.env%unit) then
         call write_energy(iprop,res,fres, &
            & (runtyp.eq.p_run_hess).or.(runtyp.eq.p_run_ohess))
         call close_file(iprop)
      endif
   endif

   if ((runtyp.eq.p_run_opt).or.(runtyp.eq.p_run_ohess).or. &
       (runtyp.eq.p_run_omd).or.(runtyp.eq.p_run_screen).or. &
       (runtyp.eq.p_run_metaopt)) then
      call generateFileName(tmpname, 'xtbopt', extension, mol%ftype)
      write(env%unit,'(/,a,1x,a,/)') &
         "optimized geometry written to:",tmpname
      call open_file(ich,tmpname,'w')
      call writeMolecule(mol, ich, energy=res%e_total, gnorm=res%gnorm)
      call close_file(ich)
   endif

   call write_energy(env%unit,res,fres, &
      & (runtyp.eq.p_run_hess).or.(runtyp.eq.p_run_ohess))

!! ------------------------------------------------------------------------
!  xtb molecular dynamics
   if ((runtyp.eq.p_run_md).or.(runtyp.eq.p_run_omd)) then
      if (metaset%maxsave .gt. 0) then
         if (mol%npbc > 0) then
            call env%error("Metadynamic under PBC is not implemented", source)
         endif
         call metadyn_header(env%unit)
      else
         call md_header(env%unit)
      endif
      fixset%n = 0 ! no fixing for MD runs
      call start_timing(6)
      idum = 0
      if (shake_md) call init_shake(mol%n,mol%at,mol%xyz,wfn%wbo)
      call md &
      &     (env,mol,wfn,calc, &
      &      egap,etemp,maxscciter,etot,g,sigma,0,temp_md,idum)
      call stop_timing(6)
   endif

!! ------------------------------------------------------------------------
!  metadynamics
   if (runtyp.eq.p_run_metaopt) then
      if (mol%npbc > 0) then
         call env%warning("Metadynamic under PBC is not implemented", source)
      endif
      call metadyn_header(env%unit)
      ! check if ANCOPT already convered
      if (murks) then
         call env%error('Optimization did not converge, aborting', source)
      endif
      write(env%unit,'(1x,"output written to xtbmeta.log")')
      call open_file(ich,'xtbmeta.log','w')
      call writeMolecule(mol, ich, fileType%xyz, energy=etot, gnorm=norm2(g))
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
         &     (env, mol,wfn,calc, &
         &      egap,etemp,maxscciter,optset%maxoptcycle,etot,g,sigma,optset%optlev,verbose,.true.,murks)
         if (.not.verbose) then
            write(env%unit,'("current energy:",1x,f20.8)') etot
         endif
         if (murks) then
            call close_file(ich)
            write(env%unit,'(/,3x,"***",1x,a,1x,"***",/)') &
               "FAILED TO CONVERGE GEOMETRY OPTIMIZATION"
            call touch_file('NOT_CONVERGED')
         endif
         call writeMolecule(mol, ich, fileType%xyz, energy=etot, gnorm=norm2(g))
      enddo
      call close_file(ich)
      call stop_timing(6)
   endif


!! ------------------------------------------------------------------------
!  simulated annealing
   if (runtyp.eq.p_run_siman) then
      call siman_header(env%unit)
      call env%error("SIMAN has been deprecated in favor of CREST", source)
   endif

!! ------------------------------------------------------------------------
!  path finder
   if (runtyp.eq.p_run_path) then
      call rmsdpath_header(env%unit)
      if (mol%npbc > 0) then
         call env%warning("Metadynamics under PBC are not implemented", source)
      endif
      call start_timing(4)
      call bias_path(env,mol,wfn,calc,egap,etemp,maxscciter,etot,g,sigma)
      call stop_timing(4)
   endif

!! ------------------------------------------------------------------------
!  screen over input structures
   if (runtyp.eq.p_run_screen) then
      call start_timing(8)
      call screen(env,mol,wfn,calc,egap,etemp,maxscciter,etot,g,sigma)
      call stop_timing(8)
   endif

!! ------------------------------------------------------------------------
!  gmd for averaged interaction energies
   if (runtyp.eq.p_run_gmd) then
      call env%error('GMD option is not supported anymore', source)
   endif

!! ------------------------------------------------------------------------
!  mode following for conformer search
   if (runtyp.eq.p_run_modef) then
      if (mol%npbc > 0) then
         call env%warning("Modefollowing under PBC is not implemented", source)
      endif
      call start_timing(9)
      call modefollow(env,mol,wfn,calc,egap,etemp,maxscciter,etot,g,sigma)
      call stop_timing(9)
   endif

!! ------------------------------------------------------------------------
!  optimize along MD from xtb.trj for conformer searches
   if (runtyp.eq.p_run_mdopt) then
      call start_timing(10)
      call mdopt(env,mol,wfn,calc,egap,etemp,maxscciter,etot,g,sigma)
      call stop_timing(10)
   endif

!! ------------------------------------------------------------------------
   if (runtyp.eq.p_run_reactor) then
      call reactor_header(env%unit)
      call env%error('The nano-reactor has been moved!', source)
   endif

!! ========================================================================
!  to further speed up xtb calculations we dump our most important
!  quantities in a restart file, so we can save some precious seconds
   if (restart) then
      call writeRestart(env,wfn,'xtbrestart',gfn_method)
   endif
   call wfn%deallocate

!! ========================================================================
!  to accelerate the fitting procedure we can print out various
!  properties when in fit-modus. So you have to return any property
!  you want to use back to the mainlevel and implement the printout here

!! ========================================================================
!  we may have generated some non-fatal errors, which have been saved,
!  so we should tell the user, (s)he may want to know what went wrong
   call env%show("Runtime exception occurred")
   call raise('F','Some non-fatal runtime exceptions were caught,'// &
   &              ' please check:',1)

!! ========================================================================
!  print all files xtb interacted with while running (for debugging mainly)
   if (verbose) then
   write(env%unit,'(a)')
   write(env%unit,'(72("-"))')
   call print_filelist(env%unit)
   endif

!! ========================================================================
!  make some post processing afterward, show some timings and stuff
   write(env%unit,'(a)')
   write(env%unit,'(72("-"))')
   call stop_timing_run
   call stop_timing(1)
   call prdate('E')
   write(env%unit,'(72("-"))')
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

   write(env%unit,'(a)')
   call terminate(0)

contains
subroutine check_cold_fusion(mol)
   type(TMolecule), intent(in) :: mol
   integer :: iat, jat
   character(len=10) :: a10
   character(len=20) :: a20
   logical :: cold_fusion
   cold_fusion = .false.
   do iat = 1, len(mol)
      do jat = 1, iat-1
         if (mol%dist(jat, iat) < 1.0e-9_wp) then
            cold_fusion = .true.
            write(a20, '(a,i0,"-",a,i0)') &
               &  trim(mol%sym(jat)), jat, trim(mol%sym(iat)), iat
            write(a10, '(es10.3)') mol%dist(jat, iat)
            call env%error("Found *very* short distance of "//a10//" for "//&
               &           trim(a20))
         endif
      enddo
   enddo
   if (cold_fusion) then
      call env%error("XTB REFUSES TO CONTINUE WITH THIS CALCULATION!")
      call env%terminate("Some atoms in the start geometry are *very* close")
   endif
end subroutine check_cold_fusion
end program XTBprog

