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

module xtb_prog_main
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stderr
   use xtb_mctc_timings
   use xtb_mctc_systools
   use xtb_mctc_convert
   use xtb_mctc_param
   use xtb_type_molecule
   use xtb_type_calculator
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment, only : TEnvironment, init
   use xtb_prog_argparser
   use xtb_solv_state
   use xtb_setparam
   use xtb_sphereparam
   use xtb_scanparam
   use xtb_splitparam
   use xtb_fixparam
   use xtb_constrain_param, only : read_userdata
   use xtb_shake, only: init_shake
   use xtb_gfnff_shake, only: gff_init_shake => init_shake
   use xtb_embedding, only : init_pcem
   use xtb_io_reader, only : readMolecule
   use xtb_io_writer, only : writeMolecule
   use xtb_mctc_filetypes, only : fileType, getFileType, generateFileMetaInfo, &
      & generateFileName
   use xtb_readin
   use xtb_printout
   use xtb_setmod
   use xtb_propertyoutput
   use xtb_io_writer_turbomole, only : writeResultsTurbomole
   use xtb_io_writer_orca, only : writeResultsOrca
   use xtb_io_writer_gaussian, only : writeResultsGaussianExternal
   use xtb_restart
   use xtb_readparam
   use xtb_scc_core, only : iniqshell
   use xtb_aespot, only : get_radcn
   use xtb_iniq, only : iniqcn
   use xtb_eeq
   use xtb_disp_ncoord, only : ncoord_gfn, dncoord_erf, dncoord_d3, ncoord_erf, &
      & ncoord_d3
   use xtb_basis
   use xtb_axis, only : axis3
   use xtb_hessian, only : numhess
   use xtb_dynamic, only : md
   use xtb_modef, only : modefollow
   use xtb_mdoptim, only : mdopt
   use xtb_screening, only : screen
   use xtb_xtb_calculator
   use xtb_gfnff_calculator
   use xtb_iff_calculator, only : TIFFCalculator
   use xtb_paramset
   use xtb_xtb_gfn0
   use xtb_xtb_gfn1
   use xtb_xtb_gfn2
   use xtb_main_setup
   use xtb_main_defaults, only : initDefaults
   use xtb_main_json, only : main_json, write_json_gfnff_lists
   use xtb_geoopt
   use xtb_metadynamic
   use xtb_biaspath
   use xtb_coffee
   use xtb_disp_dftd3param
   use xtb_disp_dftd4
   use xtb_gfnff_param, only : gff_print
   use xtb_gfnff_topology, only : TPrintTopo
   use xtb_gfnff_convert, only : struc_convert
   use xtb_scan
   use xtb_kopt
   use xtb_iff_iffprepare, only : prepare_IFF
   use xtb_iff_data, only : TIFFData
   use xtb_oniom, only : oniom_input, TOniomCalculator
   implicit none
   private

   public :: xtbMain


contains


subroutine xtbMain(env, argParser)

   !> Source of errors in the main program unit
   character(len=*), parameter :: source = "prog_main"

   type(TEnvironment), intent(inout) :: env

   type(TArgParser), intent(inout) :: argParser

!! ========================================================================
!  use some wrapper types to bundle information together
   type(TMolecule) :: mol
   type(scc_results) :: res
   class(TCalculator), allocatable :: calc
   type(freq_results) :: fres
   type(TRestart) :: chk
   type(chrg_parameter) :: chrgeq
   type(TIFFData), allocatable :: iff_data
   type(oniom_input) :: oniom
!  store important names and stuff like that in FORTRAN strings
   character(len=:),allocatable :: fname    ! geometry input file
   character(len=:),allocatable :: xcontrol ! instruction file
   character(len=:),allocatable :: xrc      ! global instruction file
   character(len=:),allocatable :: fnv      ! parameter file
   character(len=:),allocatable :: tmpname  ! temporary string
   character(len=:),allocatable :: cdum     ! temporary string
   character(len=:),allocatable :: extension, basename, directory
   integer :: ftype

!! ========================================================================
!  default names for important files in xtb
   character(len=*),parameter :: p_fname_rc = '.xtbrc'
   character(len=*),parameter :: p_fname_param_gfn0  = 'param_gfn0-xtb.txt'
   character(len=*),parameter :: p_fname_param_gfn1  = 'param_gfn1-xtb.txt'
   character(len=*),parameter :: p_fname_param_gfn2  = 'param_gfn2-xtb.txt'
   character(len=*),parameter :: p_fname_param_gfnff = '.param_gfnff.xtb'
   character(len=*),parameter :: p_fname_param_ipea  = 'param_ipea-xtb.txt'

   integer :: gsolvstate
   integer :: i,j,k,l,idum
   integer :: ich,ictrl,iprop ! file handle
   real(wp) :: sigma(3,3)
   real(wp),allocatable :: cn  (:)
   real(wp),allocatable :: sat (:)
   real(wp),allocatable :: g   (:,:)
   real(wp) :: vec3(3)
   type(TxTBParameter) :: globpar
   real(wp),allocatable :: dcn (:,:,:)
   real(wp),allocatable :: dq  (:,:,:)
   real(wp),allocatable :: dumdumdum  (:,:,:)
   real(wp),allocatable :: q  (:)
   real(wp),allocatable :: ql  (:)
   real(wp),allocatable :: qr  (:)

!! ------------------------------------------------------------------------
   integer,external :: ncore

!! ------------------------------------------------------------------------
   logical :: struc_conversion_done = .false.

!! ========================================================================
!  debugging variables for numerical gradient
   logical, parameter    :: gen_param = .false.
   logical, parameter    :: debug = .false.
   type(TRestart) :: wf0
   real(wp),allocatable  :: coord(:,:),numg(:,:),gdum(:,:)
   real(wp) :: sdum(3,3)
   real(wp),parameter    :: step = 0.00001_wp, step2 = 0.5_wp/step
   real(wp) :: er,el
   logical  :: coffee ! if debugging gets really though, get a coffee

!! ------------------------------------------------------------------------
!  undocumented and unexplainable variables go here
   integer  :: nFiles, iFile
   integer  :: rohf,err
   real(wp) :: dum5,egap,etot,ipeashift
   real(wp) :: zero,t0,t1,w0,w1,acc,etot2,g298
   real(wp) :: one,two
   real(wp) :: ea,ip
   real(wp) :: vomega,vfukui
   real(wp),allocatable :: f_plus(:), f_minus(:)
   type(TRestart) :: wf_p, wf_m
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
   logical :: cold_fusion

!  OMP stuff
   integer :: TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
   integer :: nproc

   type(TPrintTopo) :: printTopo ! gfnff topology printout list

   xenv%home = env%xtbhome
   xenv%path = env%xtbpath


   ! ------------------------------------------------------------------------
   !> read the command line arguments
   call parseArguments(env, argParser, xcontrol, fnv, acc, lgrad, &
      & restart, gsolvstate, strict, copycontrol, coffee, printTopo, oniom)

   nFiles = argParser%countFiles()
   select case(nFiles)
   case(0)
      if (.not.coffee) then
         if(printTopo%warning) call env%error("Eventually the input file was given to wrtopo as an argument.",source)
         call env%error("No input file given, so there is nothing to do", source)
      else
         fname = 'coffee'
      end if
   case(1:)
      do iFile = 1, nFiles-1
         call argParser%nextFile(fname)
         call env%warning("Input file '"//fname//"' will be ignored", source)
      end do
      call argParser%nextFile(fname)
   end select

   if (.not.allocated(xcontrol)) then
      if (copycontrol) then
         xcontrol = 'xtb.inp'
      else
         xcontrol = fname
      end if
   end if

   call env%checkpoint("Command line argument parsing failed")


   ! ------------------------------------------------------------------------
   !> read the detailed input file
   call rdcontrol(xcontrol, env, copy_file=copycontrol)

   call env%checkpoint("Reading '"//xcontrol//"' failed")


   ! ------------------------------------------------------------------------
   !> read dot-Files before reading the rc and after reading the xcontrol
   !> Total molecular charge
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

   !> Number of unpaired electrons
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
   
   !> efield read: gfnff only
   call open_file(ich,'.EFIELD','r')
   if (ich.ne.-1) then
      call getline(ich,cdum,iostat=err)
      if (err /= 0) then
         call env%error('.EFIELD is empty!', source)
      else
         call set_efield(env,cdum)
         call close_file(ich)
      end if
   endif

   call env%checkpoint("Reading multiplicity from file failed")


   ! ------------------------------------------------------------------------
   !> read the xtbrc if you can find it (use rdpath directly instead of xfind)
   call rdpath(env%xtbpath, p_fname_rc, xrc, exist)
   if (exist) then
      call rdcontrol(xrc, env, copy_file=.false.)

      call env%checkpoint("Reading '"//xrc//"' failed")
   endif


   ! ------------------------------------------------------------------------
   !> FIXME: some settings that are still not automatic
   !> Make sure GFN0-xTB uses the correct exttyp
   if(set%gfn_method == 0)  call set_exttyp('eht')
   rohf = 1 ! HS default
   egap = 0.0_wp
   ipeashift = 0.0_wp


   ! ========================================================================
   !> no user interaction up to now, time to show off!
   !> print the xtb banner with version number and compilation date
   !> making a fancy version of this is hard, x is difficult in ASCII art
   call xtb_header(env%unit)
   !> make sure you cannot blame us for destroying your computer
   call disclamer(env%unit)
   !> how to cite this program
   call citation(env%unit)
   !> print current time
   call prdate('S')


   ! ------------------------------------------------------------------------
   !> get molecular structure
   if (coffee) then ! it's coffee time
      fname = 'caffeine'
      call get_coffee(mol)
      call generateFileMetaInfo(fname, directory, basename, extension)
   else                                                              
      call generateFileMetaInfo(fname, directory, basename, extension)
      ftype = getFileType(fname)
      call open_file(ich, fname, 'r')
      call readMolecule(env, mol, ich, ftype)
      call close_file(ich)
      if (mol%info%two_dimensional) then
         call env%warning("Two dimensional input structure detected", source)
      end if

      ! Special CT input file case
      if (mol%chrg /= 0.0_wp) then
         if (set%ichrg == 0) then
            set%ichrg = nint(mol%chrg)
         else
            call env%warning("Charge in sdf/mol input was overwritten", source)
         end if
      end if

      call env%checkpoint("reading geometry input '"//fname//"' failed")
   endif


   ! ------------------------------------------------------------------------
   !> initialize the global storage
   call init_fix(mol%n)
   call init_split(mol%n)
   call init_constr(mol%n,mol%at)
   call init_scan
   call init_walls
   call init_pcem
   if (set%runtyp.eq.p_run_bhess) then
      call init_bhess(mol%n)
   else
      call init_metadyn(mol%n,metaset%maxsave)
   end if
   call load_rmsdbias(rmsdset,mol%n,mol%at,mol%xyz)

   ! ------------------------------------------------------------------------
   !> get some memory
   allocate(cn(mol%n),sat(mol%n),g(3,mol%n), source = 0.0_wp)
   atmass = atomic_mass(mol%at) * autoamu ! from splitparam.f90
   set%periodic = mol%npbc > 0
   if (mol%npbc == 0) then
      if (set%do_cma_trafo) then
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

   !> initialize time step for MD if requested autocomplete
   if (set%tstep_md < 0.0_wp) then
      set%tstep_md = (minval(atmass)/(atomic_mass(1)*autoamu))**(1.0_wp/3.0_wp)
   endif

   mol%chrg = real(set%ichrg, wp)
   mol%uhf = set%nalphabeta
   call initrand

   call setup_summary(env%unit,mol%n,fname,xcontrol,chk%wfn,xrc)

   if(set%fit) acc=0.2 ! higher SCF accuracy during fit

   ! ------------------------------------------------------------------------
   !> 2D => 3D STRUCTURE CONVERTER
   ! ------------------------------------------------------------------------
   if (mol%info%two_dimensional) then
      call struc_convert (env,restart,mol,chk,egap,set%etemp,set%maxscciter, &
                       &  set%optset%maxoptcycle,etot,g,sigma)
      struc_conversion_done = .true.
      mol%info%two_dimensional = .false.
    end if

   ! ------------------------------------------------------------------------
   !> CONSTRAINTS & SCANS
   !> now we are at a point that we can check for requested constraints
   call read_userdata(xcontrol,env,mol)

   !> initialize metadynamics
   call load_metadynamic(metaset,mol%n,mol%at,mol%xyz)

   !> restraining potential
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
   !  fragmentation for CMA constrain
   if(iatf1.eq.0.and.iatf2.eq.0) then
      call ncoord_erf(mol%n,mol%at,mol%xyz,cn)
      call splitm(mol%n,mol%at,mol%xyz,cn)
   endif
   call splitprint(mol%n,mol%at,mol%xyz)

   if (set%verbose) then
      call fix_info(env%unit,mol%n,mol%at,mol%xyz)
      call pot_info(env%unit,mol%n,mol%at,mol%xyz)
   endif

   ! ------------------------------------------------------------------------
   !> write copy of detailed input
   if (copycontrol) then
      call open_set(ictrl,xcontrol)
      call write_set(ictrl)
      call close_set(ictrl)
   endif

   ! ------------------------------------------------------------------------
   !> if you have requested a define we stop here...
   if (set%define) then
      if (set%verbose) call main_geometry(env%unit,mol)
      call eval_define(set%veryverbose)
   endif
   call env%show('Please study the warnings concerning your input carefully')
   call raise('F', 'Please study the warnings concerning your input carefully')

   ! ========================================================================
   !> From here we switch to the method setup
   !> enable error on warnings
   if (strict) call mctc_strict
   env%strict = strict

   !> one last check on the input geometry
   call check_cold_fusion(env, mol, cold_fusion)
   if (cold_fusion) then
      call env%error("XTB REFUSES TO CONTINUE WITH THIS CALCULATION!")
      call env%terminate("Some atoms in the start geometry are *very* close")
   endif

   !> check if someone is still using GFN3...
   if (set%gfn_method.eq.3) then
      call env%terminate('This is an internal error, please use gfn_method=2!')
   end if

   ! ------------------------------------------------------------------------
   !> Print the method header and select the parameter file
   if (.not.allocated(fnv)) then
      select case(set%runtyp)
      case default
         call env%terminate('This is an internal error, please define your runtypes!')
      case(p_run_scc,p_run_grad,p_run_opt,p_run_hess,p_run_ohess,p_run_bhess, &
            p_run_md,p_run_omd,p_run_path,p_run_screen, &
            p_run_modef,p_run_mdopt,p_run_metaopt)
        if (set%mode_extrun.eq.p_ext_gfnff) then
            fnv=xfind(p_fname_param_gfnff)
        else
           if(set%gfn_method.eq.0) then
              fnv=xfind(p_fname_param_gfn0)
           endif
           if(set%gfn_method.eq.1) then
              fnv=xfind(p_fname_param_gfn1)
           endif
           if(set%gfn_method.eq.2) then
              fnv=xfind(p_fname_param_gfn2)
           endif
        end if
      case(p_run_vip,p_run_vea,p_run_vipea,p_run_vfukui,p_run_vomega)
         if(set%gfn_method.eq.0) then
            fnv=xfind(p_fname_param_gfn0)
         endif
         if(set%gfn_method.eq.1) then
            fnv=xfind(p_fname_param_ipea)
         endif
         if(set%gfn_method.eq.2) then
            fnv=xfind(p_fname_param_gfn2)
         endif
      end select
   endif

   !-------------------------------------------------------------------------
   !> Perform a precomputation of electronic properties for xTB-IFF
   if(set%mode_extrun == p_ext_iff) then
      allocate(iff_data)
      call prepare_IFF(env, mol, iff_data)
      call env%checkpoint("Could not generate electronic properties")
   end if

   ! ------------------------------------------------------------------------
   !> Obtain the parameter data
   call newCalculator(env, mol, calc, fnv, restart, acc, oniom, iff_data)
   call env%checkpoint("Could not setup single-point calculator")

   call initDefaults(env, calc, mol, gsolvstate)
   call env%checkpoint("Could not setup defaults")


   ! ------------------------------------------------------------------------
   !> initial guess, setup wavefunction
   select type(calc)
   type is(TxTBCalculator)
      call chk%wfn%allocate(mol%n,calc%basis%nshell,calc%basis%nao)

      ! Make sure number of electrons is initialized an multiplicity is consistent
      chk%wfn%nel = nint(sum(mol%z) - mol%chrg)
      chk%wfn%nopen = mol%uhf
      if(chk%wfn%nopen == 0 .and. mod(chk%wfn%nel,2) /= 0) chk%wfn%nopen=1

      !> EN charges and CN
      if (set%gfn_method.lt.2) then
         call ncoord_d3(mol%n,mol%at,mol%xyz,cn)
      else
         call ncoord_gfn(mol%n,mol%at,mol%xyz,cn)
      endif
      if (mol%npbc > 0) then
         chk%wfn%q = real(set%ichrg,wp)/real(mol%n,wp)
      else
         if (set%guess_charges.eq.p_guess_gasteiger) then
            call iniqcn(mol%n,mol%at,mol%z,mol%xyz,set%ichrg,1.0_wp,chk%wfn%q,cn,set%gfn_method,.true.)
         else if (set%guess_charges.eq.p_guess_goedecker) then
            call ncoord_erf(mol%n,mol%at,mol%xyz,cn)
            call goedecker_chrgeq(mol%n,mol%at,mol%xyz,real(set%ichrg,wp),cn,dcn,chk%wfn%q,dq,er,g,&
               .false.,.false.,.false.)
         else
            call ncoord_gfn(mol%n,mol%at,mol%xyz,cn)
            chk%wfn%q = real(set%ichrg,wp)/real(mol%n,wp)
         end if
      end if
      !> initialize shell charges from gasteiger charges
      call iniqshell(calc%xtbData,mol%n,mol%at,mol%z,calc%basis%nshell,chk%wfn%q,chk%wfn%qsh,set%gfn_method)
   end select

   ! ------------------------------------------------------------------------
   !> printout a header for the exttyp
   call calc%writeInfo(env%unit, mol)

   call delete_file('.sccnotconverged')

   call env%checkpoint("Setup for calculation failed")

   select type(calc)
   type is(TxTBCalculator)
      if (restart.and.calc%xtbData%level /= 0) then ! only in first run
         call readRestart(env,chk%wfn,'xtbrestart',mol%n,mol%at,set%gfn_method,exist,.true.)
      endif
      calc%etemp = set%etemp
      calc%maxiter = set%maxscciter
      ipeashift = calc%xtbData%ipeashift
   type is(TOniomCalculator)
      select type(xtb => calc%real_low)
      type is(TxTBCalculator)
         call chk%wfn%allocate(mol%n,xtb%basis%nshell,xtb%basis%nao)
         if (restart) then ! only in first run
            call readRestart(env,chk%wfn,'xtbrestart',mol%n,mol%at,set%gfn_method,exist,.true.)
         endif
      end select
   end select

   ! ========================================================================
   !> the SP energy which is always done
   call start_timing(2)
   call calc%singlepoint(env,mol,chk,2,exist,etot,g,sigma,egap,res)
   call stop_timing(2)
   select type(calc)
   type is(TGFFCalculator)
     gff_print=.false.
   end select
   call env%checkpoint("Single point calculation terminated")

   !> write 2d => 3d converted structure   
   if (struc_conversion_done) then
      call generateFileName(tmpname, 'gfnff_convert', extension, mol%ftype)
      write(env%unit,'(10x,a,1x,a,/)') &
         "converted geometry written to:",tmpname
      call open_file(ich,tmpname,'w')
      call writeMolecule(mol, ich, energy=res%e_total, gnorm=res%gnorm)
      call close_file(ich)
   end if
   
   ! ========================================================================
   !> determine kopt for bhess including final biased geometry optimization
   if (set%runtyp.eq.p_run_bhess) then
      call set_metadynamic(metaset,mol%n,mol%at,mol%xyz)
      call get_kopt (metaset,env,restart,mol,chk,calc,egap,set%etemp,set%maxscciter, &
         & set%optset%maxoptcycle,set%optset%optlev,etot,g,sigma,acc)
   end if

   ! ------------------------------------------------------------------------
   !> numerical gradient for debugging purposes
   if (debug) then
      !  generate a warning to keep release versions from calculating numerical gradients
      call env%warning('XTB IS CALCULATING NUMERICAL GRADIENTS, RESET DEBUG FOR RELEASE!')
      print'(/,"analytical gradient")'
      print *, g
      allocate( coord(3,mol%n), source = mol%xyz )
      allocate( numg(3,mol%n),gdum(3,mol%n), source = 0.0_wp )
      wf0 = chk
      do i = 1, mol%n
         do j = 1, 3
            mol%xyz(j,i) = mol%xyz(j,i) + step
            chk = wf0
            call calc%singlepoint(env,mol,chk,0,.true.,er,gdum,sdum,egap,res)
            mol%xyz(j,i) = mol%xyz(j,i) - 2*step
            chk = wf0
            call calc%singlepoint(env,mol,chk,0,.true.,el,gdum,sdum,egap,res)
            mol%xyz(j,i) = mol%xyz(j,i) + step
            numg(j,i) = step2 * (er - el)
         enddo
      enddo
      print'(/,"numerical gradient")'
      print *, numg
      print'(/,"difference gradient")'
      print*,g-numg
   endif


   ! ------------------------------------------------------------------------
   !  ANCopt
   if ((set%runtyp.eq.p_run_opt).or.(set%runtyp.eq.p_run_ohess).or. &
      &   (set%runtyp.eq.p_run_omd).or.(set%runtyp.eq.p_run_screen).or. &
      &   (set%runtyp.eq.p_run_metaopt)) then
      if (set%opt_engine.eq.p_engine_rf) &
         call ancopt_header(env%unit,set%veryverbose)
      call start_timing(3)
      call geometry_optimization &
         &     (env, mol,chk,calc, &
         &      egap,set%etemp,set%maxscciter,set%optset%maxoptcycle,etot,g,sigma,set%optset%optlev,.true.,.false.,murks)
      res%e_total = etot
      res%gnorm = norm2(g)
      if (nscan.gt.0) then
         call relaxed_scan(env,mol,chk,calc)
      endif
      if (murks) then
         call generateFileName(tmpname, 'xtblast', extension, mol%ftype)
         write(env%unit,'(/,a,1x,a,/)') &
            "last geometry written to:",tmpname
         call open_file(ich,tmpname,'w')
         call writeMolecule(mol, ich, energy=res%e_total, gnorm=res%gnorm)
         call close_file(ich)
         call env%terminate("Geometry optimization failed")
      end if
      call stop_timing(3)
   endif


   ! ------------------------------------------------------------------------
   !> automatic VIP and VEA single point (maybe after optimization)
   if (set%runtyp.eq.p_run_vip.or.set%runtyp.eq.p_run_vipea &
      & .or.set%runtyp.eq.p_run_vomega) then
      call start_timing(2)
      call vip_header(env%unit)
      mol%chrg = mol%chrg + 1
      chk%wfn%nel = chk%wfn%nel-1
      if (mod(chk%wfn%nel,2).ne.0) chk%wfn%nopen = 1
      call calc%singlepoint(env,mol,chk,1,exist,etot2,g,sigma,egap,res)
      ip=etot2-etot-ipeashift
      write(env%unit,'(72("-"))')
      write(env%unit,'("empirical IP shift (eV):",f10.4)') &
         &                  autoev*ipeashift
      write(env%unit,'("delta SCC IP (eV):",f10.4)') autoev*ip
      write(env%unit,'(72("-"))')
      mol%chrg = mol%chrg - 1
      chk%wfn%nel = chk%wfn%nel+1
      call stop_timing(2)
   endif

   if (set%runtyp.eq.p_run_vea.or.set%runtyp.eq.p_run_vipea &
      & .or.set%runtyp.eq.p_run_vomega) then
      call start_timing(2)
      call vea_header(env%unit)
      mol%chrg = mol%chrg - 1
      chk%wfn%nel = chk%wfn%nel+1
      if (mod(chk%wfn%nel,2).ne.0) chk%wfn%nopen = 1
      call calc%singlepoint(env,mol,chk,1,exist,etot2,g,sigma,egap,res)
      ea=etot-etot2-ipeashift
      write(env%unit,'(72("-"))')
      write(env%unit,'("empirical EA shift (eV):",f10.4)') &
         &                  autoev*ipeashift
      write(env%unit,'("delta SCC EA (eV):",f10.4)') autoev*ea
      write(env%unit,'(72("-"))')

      mol%chrg = mol%chrg + 1
      chk%wfn%nel = chk%wfn%nel-1
      call stop_timing(2)
   endif


   ! ------------------------------------------------------------------------
   !> vomega (electrophilicity) index
   if (set%runtyp.eq.p_run_vomega) then
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


   ! ------------------------------------------------------------------------
   !> Fukui Index from Mulliken population analysis
   if (set%runtyp.eq.p_run_vfukui) then
      allocate(f_plus(mol%n),f_minus(mol%n))
      write(env%unit,'(a)')
      write(env%unit,'("Fukui index Calculation")')
      wf_p%wfn=chk%wfn
      wf_m%wfn=chk%wfn
      mol%chrg = mol%chrg + 1
      wf_p%wfn%nel = wf_p%wfn%nel+1
      if (mod(wf_p%wfn%nel,2).ne.0) wf_p%wfn%nopen = 1
      call calc%singlepoint(env,mol,wf_p,1,exist,etot2,g,sigma,egap,res)
      f_plus=wf_p%wfn%q-chk%wfn%q

      mol%chrg = mol%chrg - 2
      wf_m%wfn%nel = wf_m%wfn%nel-1
      if (mod(wf_m%wfn%nel,2).ne.0) wf_m%wfn%nopen = 1
      call calc%singlepoint(env,mol,wf_m,1,exist,etot2,g,sigma,egap,res)
      f_minus=chk%wfn%q-wf_m%wfn%q
      write(env%unit,'(a)')
      write(env%unit, '(1x,"    #        f(+)     f(-)     f(0)")')
      do i=1,mol%n
         write(env%unit,'(i6,a4,2f9.3,2f9.3,2f9.3)') i, mol%sym(i), f_plus(i), f_minus(i), 0.5d0*(wf_p%wfn%q(i)-wf_m%wfn%q(i))
      enddo
      mol%chrg = mol%chrg + 1
      deallocate(f_plus,f_minus)
   endif


   ! ------------------------------------------------------------------------
   !> numerical hessian calculation
   if ((set%runtyp.eq.p_run_hess).or.(set%runtyp.eq.p_run_ohess).or.(set%runtyp.eq.p_run_bhess)) then
      if (set%runtyp.eq.p_run_bhess .and. set%mode_extrun.ne.p_ext_turbomole) then
         call generic_header(env%unit,"Biased Numerical Hessian",49,10)
      else if (set%runtyp.eq.p_run_bhess .and. set%mode_extrun.eq.p_ext_turbomole) then
         call generic_header(env%unit,"Biased Analytical TM Hessian",49,10)
      else if (set%mode_extrun.eq.p_ext_turbomole) then
         call generic_header(env%unit,"Analytical TM Hessian",49,10)
      else
         call numhess_header(env%unit)
      end if
      if (mol%npbc > 0) then
         call env%error("Phonon calculations under PBC are not implemented", source)
      endif
      call start_timing(5)
      call numhess &
         &       (env,mol,chk,calc, &
         &        egap,set%etemp,set%maxscciter,etot,g,sigma,fres)
      call stop_timing(5)

      call env%checkpoint("Hessian calculation terminated")
   endif

   ! reset the gap, since it is currently not updated in ancopt and numhess
   if (allocated(chk%wfn%emo)) then
      res%hl_gap = chk%wfn%emo(chk%wfn%ihomo+1)-chk%wfn%emo(chk%wfn%ihomo)
   end if

   call env%checkpoint("Calculation terminated")

   ! ========================================================================
   !> PRINTOUT SECTION
   if (allocated(set%property_file)) then
      call open_file(iprop,set%property_file,'w')
      if (iprop.eq.-1) then
         iprop = env%unit
         deallocate(set%property_file)
      else
         write(env%unit,'(/,a)') "Property printout bound to '"//set%property_file//"'"
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
      call writeResultsTurbomole(mol, energy=etot, gradient=g, sigma=sigma)
      if (allocated(basename)) then
         cdum = basename // '.engrad'
      else
         cdum = 'xtb-orca.engrad'
      end if
      call open_file(ich, cdum, 'w')
      call writeResultsOrca(ich, mol, etot, g)
      call close_file(ich)
   end if
   if (mol%ftype .eq. fileType%gaussian) then
      if (allocated(basename)) then
         cdum = basename // '.EOu'
      else
         cdum = 'xtb-gaussian.EOu'
      end if
      call open_file(ich, cdum, 'w')
      call writeResultsGaussianExternal(ich, etot, res%dipole, g)
      call close_file(ich)
   end if

   if(set%periodic)then
      write(*,*)'Periodic properties'
   else
      select type(calc)
      type is(TxTBCalculator)
         call main_property(iprop,env,mol,chk%wfn,calc%basis,calc%xtbData,res, &
            & calc%solvation,acc)
         call main_cube(set%verbose,mol,chk%wfn,calc%basis,res)
      end select
   endif


   if (set%pr_json) then
      select type(calc)
      type is(TxTBCalculator)
         call open_file(ich,'xtbout.json','w')
         call main_json(ich, &
            mol,chk%wfn,calc%basis,res,fres)
         call close_file(ich)
      end select
   endif
   if(printTopo%any()) then
     select type(calc)
       type is(TGFFCalculator)
         call write_json_gfnff_lists(mol%n,calc%topo,chk%nlist,printTopo)
     end select
   endif
   if ((set%runtyp.eq.p_run_opt).or.(set%runtyp.eq.p_run_ohess).or. &
      (set%runtyp.eq.p_run_omd).or.(set%runtyp.eq.p_run_screen).or. &
      (set%runtyp.eq.p_run_metaopt).or.(set%runtyp.eq.p_run_bhess)) then
      call main_geometry(iprop,mol)
   endif

   if ((set%runtyp.eq.p_run_hess).or.(set%runtyp.eq.p_run_ohess).or.(set%runtyp.eq.p_run_bhess)) then
      call generic_header(iprop,'Frequency Printout',49,10)
      call main_freq(iprop,mol,chk%wfn,fres)
   endif

   if (allocated(set%property_file)) then
      if (iprop.ne.-1 .and. iprop.ne.env%unit) then
         call write_energy(iprop,res,fres, &
            & (set%runtyp.eq.p_run_hess).or.(set%runtyp.eq.p_run_ohess).or.(set%runtyp.eq.p_run_bhess))
         call close_file(iprop)
      endif
   endif

   if ((set%runtyp.eq.p_run_opt).or.(set%runtyp.eq.p_run_ohess).or. &
      (set%runtyp.eq.p_run_omd).or.(set%runtyp.eq.p_run_screen).or. &
      (set%runtyp.eq.p_run_metaopt).or.(set%runtyp.eq.p_run_bhess)) then
      call generateFileName(tmpname, 'xtbopt', extension, mol%ftype)
      write(env%unit,'(/,a,1x,a,/)') &
         "optimized geometry written to:",tmpname
      call open_file(ich,tmpname,'w')
      call writeMolecule(mol, ich, energy=res%e_total, gnorm=res%gnorm)
      call close_file(ich)
   endif

   select type(calc)
   type is(TxTBCalculator)
      call write_energy(env%unit,res,fres, &
        & (set%runtyp.eq.p_run_hess).or.(set%runtyp.eq.p_run_ohess).or.(set%runtyp.eq.p_run_bhess))
   class default
      call write_energy_gff(env%unit,res,fres, &
        & (set%runtyp.eq.p_run_hess).or.(set%runtyp.eq.p_run_ohess).or.(set%runtyp.eq.p_run_bhess))
   end select  


   ! ------------------------------------------------------------------------
   !  xtb molecular dynamics
   if ((set%runtyp.eq.p_run_md).or.(set%runtyp.eq.p_run_omd)) then
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
      select type(calc)
      class default
         if (set%shake_md) call init_shake(mol%n,mol%at,mol%xyz,chk%wfn%wbo)
      type is(TGFFCalculator)
         if (set%shake_md) call gff_init_shake(mol%n,mol%at,mol%xyz,calc%topo)
      end select
      call md &
         &     (env,mol,chk,calc, &
         &      egap,set%etemp,set%maxscciter,etot,g,sigma,0,set%temp_md,idum)
      call stop_timing(6)
   endif


   ! ------------------------------------------------------------------------
   !  metadynamics
   if (set%runtyp.eq.p_run_metaopt) then
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
            &     (env, mol,chk,calc, &
            &      egap,set%etemp,set%maxscciter,set%optset%maxoptcycle,etot,g,sigma, &
            &      set%optset%optlev,set%verbose,.true.,murks)
         if (.not.set%verbose) then
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


   ! ------------------------------------------------------------------------
   !  path finder
   if (set%runtyp.eq.p_run_path) then
      call rmsdpath_header(env%unit)
      if (mol%npbc > 0) then
         call env%warning("Metadynamics under PBC are not implemented", source)
      endif
      call start_timing(4)
      call bias_path(env,mol,chk,calc,egap,set%etemp,set%maxscciter,etot,g,sigma)
      call stop_timing(4)
   endif


   ! ------------------------------------------------------------------------
   !  screen over input structures
   if (set%runtyp.eq.p_run_screen) then
      call start_timing(8)
      call screen(env,mol,chk,calc,egap,set%etemp,set%maxscciter,etot,g,sigma)
      call stop_timing(8)
   endif


   ! ------------------------------------------------------------------------
   !  mode following for conformer search
   if (set%runtyp.eq.p_run_modef) then
      if (mol%npbc > 0) then
         call env%warning("Modefollowing under PBC is not implemented", source)
      endif
      call start_timing(9)
      call modefollow(env,mol,chk,calc,egap,set%etemp,set%maxscciter,etot,g,sigma)
      call stop_timing(9)
   endif


   ! ------------------------------------------------------------------------
   !  optimize along MD from xtb.trj for conformer searches
   if (set%runtyp.eq.p_run_mdopt) then
      call start_timing(10)
      call mdopt(env,mol,chk,calc,egap,set%etemp,set%maxscciter,etot,g,sigma)
      call stop_timing(10)
   endif


   ! ------------------------------------------------------------------------
   !  to further speed up xtb calculations we dump our most important
   !  quantities in a restart file, so we can save some precious seconds
   select type(calc)
   type is(TxTBCalculator)
      if (restart) then
         call writeRestart(env,chk%wfn,'xtbrestart',set%gfn_method)
      endif
   end select


   ! ------------------------------------------------------------------------
   !  we may have generated some non-fatal errors, which have been saved,
   !  so we should tell the user, (s)he may want to know what went wrong
   call env%show("Runtime exception occurred")
   call raise('F','Some non-fatal runtime exceptions were caught,'// &
      &           ' please check:')

   ! ------------------------------------------------------------------------
   !  print all files xtb interacted with while running (for debugging mainly)
   if (set%verbose) then
      write(env%unit,'(a)')
      write(env%unit,'(72("-"))')
      call print_filelist(env%unit)
   endif


   ! ------------------------------------------------------------------------
   !  make some post processing afterward, show some timings and stuff
   write(env%unit,'(a)')
   write(env%unit,'(72("-"))')
   call stop_timing_run
   call stop_timing(1)
   call prdate('E')
   write(env%unit,'(72("-"))')
   call prtiming(1,'total')
   call prtiming(2,'SCF')
   if ((set%runtyp.eq.p_run_opt).or.(set%runtyp.eq.p_run_ohess).or. &
      &   (set%runtyp.eq.p_run_omd).or.(set%runtyp.eq.p_run_metaopt)) then
      call prtiming(3,'ANC optimizer')
   endif
   if (set%runtyp.eq.p_run_path) then
      call prtiming(4,'path finder')
   endif
   if (((set%runtyp.eq.p_run_hess).or.(set%runtyp.eq.p_run_ohess).or.(set%runtyp.eq.p_run_bhess))) then
      if (set%mode_extrun.ne.p_ext_turbomole) then
         call prtiming(5,'analytical hessian')
      else
         call prtiming(5,'numerical hessian')
      end if
  end if
   if ((set%runtyp.eq.p_run_md).or.(set%runtyp.eq.p_run_omd).or. &
      (set%runtyp.eq.p_run_metaopt)) then
      call prtiming(6,'MD')
   endif
   if (set%runtyp.eq.p_run_screen) then
      call prtiming(8,'screen')
   endif
   if (set%runtyp.eq.p_run_modef) then
      call prtiming(9,'mode following')
   endif
   if (set%runtyp.eq.p_run_mdopt) then
      call prtiming(10,'MD opt.')
   endif

   write(env%unit,'(a)')
   call terminate(0)


end subroutine xtbMain


!> Parse command line arguments and forward them to settings
subroutine parseArguments(env, args, inputFile, paramFile, accuracy, lgrad, &
      & restart, gsolvstate, strict, copycontrol, coffee, printTopo, oniom)
   use xtb_mctc_global, only : persistentEnv

   !> Name of error producer
   character(len=*), parameter :: source = "prog_main_parseArguments"

   !> Calculation environment
   type(TEnvironment) :: env

   !> Command line argument parser
   type(TArgParser) :: args

   !> Detailed input file name
   character(len=:),allocatable,intent(out) :: inputFile

   !> Parameter file name
   character(len=:),allocatable,intent(out) :: paramFile

   !> Accuracy number for numerical thresholds
   real(wp), intent(out) :: accuracy

   !> Reference state for solvation free energies
   integer, intent(out) :: gsolvstate

   !> Restart calculation
   logical, intent(out) :: restart

   !> Handle warnings as errors
   logical, intent(out) :: strict

   !> Debugging with a lot of caffeine
   logical, intent(out) :: coffee

   !> topology printout list
   type(TPrintTopo), intent(out) :: printTopo

   !> Print the gradient to file
   logical, intent(out) :: lgrad

   !> Copy the detailed input file
   logical, intent(out) :: copycontrol

   !> Input for ONIOM model
   type(oniom_input), intent(out) :: oniom

   !> Stuff for second argument parser
!   integer  :: narg
!   character(len=p_str_length), dimension(p_arg_length) :: argv
!   type(TAtomList) :: atl
!   integer, allocatable :: list(:)

!$ integer :: omp_get_num_threads, nproc
   integer :: nFlags
   integer :: idum, ndum
   real(wp) :: ddum
   character(len=:), allocatable :: flag, sec
   logical :: exist

   set%gfn_method = 2
   coffee = .false.
   strict = .false.
   restart = .true.
   copycontrol = .false.
   lgrad = .false.
   accuracy = 1.0_wp
   gsolvstate = solutionState%gsolv

   nFlags = args%countFlags()
   call args%nextFlag(flag)
   do while(allocated(flag))
      if (len(flag) > 2 .and. flag(1:1) == '-' .and. flag(1:2) /= '--') then
         call env%warning("the use of '"//flag//"' is discouraged, "// &
            & "please use '-"//flag//"' next time", source)
         flag = '-'//flag
      end if
      select case(flag)
      case default
         call env%warning("Unknown option '"//flag//"' provided", source)

      case('-h', '--help')
         call help(env%unit)
         call terminate(0)

      case('--citation')
         call citation(env%unit)
         call terminate(0)

      case('--license')
         call disclamer(env%unit)
         call terminate(0)

      case('--version')
         call xtb_header(env%unit)
         call terminate(0)

      case('-v','--verbose')
         set%verbose = .true.

      case('-V','--very-verbose')
         set%verbose = .true.
         set%veryverbose = .true.

      case(     '--define')
         call set_define

      case('-P','--parallel')
   !$    if (.false.) then
            call env%warning('Program compiled without threading support', source)
   !$    endif
         ! Always remove next argument to keep argument parsing consistent
         call args%nextArg(sec)
   !$    if (allocated(sec)) then
   !$    if (getValue(env,sec,idum)) then
   !$       nproc = omp_get_num_threads()
   !$       call omp_set_num_threads(idum)
#ifdef WITH_MKL
   !$       call mkl_set_num_threads(idum)
#endif
   !$    endif
   !$    endif
      case('--restart')
         restart = .true.

      case('--norestart')
         restart = .false.

      case('--copy')
         copycontrol = .true.

      case('--nocopy')
         copycontrol = .false.

      case('--strict')
         strict = .true.

      case('-I', '--input')
         call args%nextArg(inputFile)
         if (.not.allocated(inputFile)) then
            call env%error("Filename for detailed input is missing", source)
         end if

      case('--namespace')
         call args%nextArg(persistentEnv%io%namespace)
         if (.not.allocated(persistentEnv%io%namespace)) then
            call env%error("Namespace argument is missing", source)
         end if

      case('--vparam')
         call args%nextArg(paramFile)
         if (.not.allocated(paramFile)) then
            call env%error("Filename for --vparam is missing", source)
         end if

      case('--coffee')
         coffee = .true.

      case('-a', '--acc')
         call args%nextArg(sec)
         if (allocated(sec)) then
            if (getValue(env,sec,ddum)) then
               if (ddum.lt.1.e-4_wp) then
                  call env%warning("We cannot provide this level of accuracy, "//&
                     & "resetted accuracy to 0.0001", source)
                  accuracy = 1.e-4_wp
               else if (ddum.gt.1.e+3_wp) then
                  call env%warning("We cannot provide this level of accuracy, "//&
                     & "resetted accuracy to 1000", source)
                  accuracy = 1.e+3_wp
               else
                  accuracy = ddum
               endif
            end if
         else
            call env%error("Accuracy is not provided", source)
         end if

      case('-c', '--chrg', '--charge')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_chrg(env,sec)
         else
            call env%error("Molecular charge is not provided", source)
         end if

      case('-u', '--uhf')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_spin(env,sec)
         else
            call env%error("Number of unpaired electrons is not provided", source)
         end if

      case('--gfn')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_gfn(env,'method',sec)
            if (sec=='0') call set_exttyp('eht')
         else
            call env%error("No method provided for --gfn option", source)
         end if

      case('--gfn1')
         call set_gfn(env,'method','1')
         call env%warning("The use of '"//flag//"' is discouraged, " //&
            & "please use '--gfn 1' next time", source)

      case('--gfn2')
         call set_gfn(env,'method','2')
         call set_gfn(env,'d4','true')

      case('--gfn0')
         call set_gfn(env,'method','0')
         call set_exttyp('eht')
         call env%warning("The use of '"//flag//"' is discouraged, " //&
            & "please use '--gfn 0' next time", source)
      
      case('--gfnff')
         call set_exttyp('ff')
      
      case('--gff')
         call set_exttyp('ff')

      case('--iff')
         call set_exttyp('iff')

      case('--oniom')
         call set_exttyp('oniom')
         call args%nextArg(sec) 

         !> To handle no argument case
         if (.not.allocated(sec)) then
            call env%error("No inner region provided for ONIOM", source)
            return
         end if
         call move_alloc(sec, oniom%first_arg)

         call args%nextArg(sec)
         if (.not.allocated(sec)) then 
            call env%warning("No method is specified for the ONIOM calculation, default gfn2:gfnff combination will be used", source)
            call move_alloc(oniom%first_arg, sec)
            !return
         end if
         inquire(file=sec, exist=exist)
         if (exist) then
            sec = read_whole_file(sec)
         end if
         call move_alloc(sec, oniom%second_arg)

      case('--etemp')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_scc(env,'temp',sec)
         else
            call env%error("Temperature in --etemp option is missing", source)
         end if

      case('--esp')
         call set_runtyp('scc')
         call set_write(env,'esp','true')

      case('--stm')
         call set_runtyp('scc')
         call set_write(env,'stm','true')

      case('--cma')
         call set_cma

      case('--tm')
         call set_exttyp('turbomole')

      case('--enso')
         call set_enso_mode

      case('--json')
         call set_write(env,'json','true')
       
      case('--ceasefiles')
         restart = .false. 
         set%verbose=.false.
         set%ceasefiles = .true.
         call set_write(env,'wiberg','false')
         call set_write(env,'charges','false')
#ifdef _WIN32
         call set_opt(env, 'logfile', 'NUL')
#else
         call set_opt(env, 'logfile', '/dev/null')
#endif         

      case('--orca')
         call set_exttyp('orca')

      case('--driver')
         call set_exttyp('driver')
         call args%nextArg(sec)
         if (allocated(sec)) then
            set%ext_driver%executable = sec
         end if

      case('--mopac')
         call set_exttyp('mopac')

      case('--pop')
         call set_write(env,'mulliken','true')

      case('--molden')
         call set_write(env,'mos','true')

      case('--dipole')
         call set_write(env,'dipole','true')

      case('--wbo')
         call set_write(env,'wiberg','true')

      case('--lmo')
         call set_write(env,'mulliken','true')
         call set_write(env,'lmo','true')

      case('--ewin')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_siman(env,'ewin',sec)
         else
            call env%error("Real argument for --ewin is missing", source)
         end if

      case('--fod')
         call set_write(env,'fod','true')
         call set_scc(env,'temp','5000.0')

      case('--iterations')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_scc(env,'maxiterations',sec)
         else
            call env%error("Integer argument for --iterations is missing", source)
         end if

      case('--cycles')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_opt(env,'maxcycle',sec)
         else
            call env%error("Integer argument for --cycles is missing", source)
         end if

      case('-g', '--gbsa')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_gbsa(env, 'solvent', sec)
            call set_gbsa(env, 'alpb', 'false')
            call set_gbsa(env, 'kernel', 'still')
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (sec == 'reference') then
                  gsolvstate = solutionState%reference
               else if (sec == 'bar1M') then
                  gsolvstate = solutionState%mol1bar
               else
                  call env%warning("Unknown reference state '"//sec//"'", source)
               end if
            end if
         else
            call env%error("No solvent name provided for GBSA", source)
         end if

      case('--alpb')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_gbsa(env, 'solvent', sec)
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (sec == 'reference') then
                  gsolvstate = solutionState%reference
               else if (sec == 'bar1M') then
                  gsolvstate = solutionState%mol1bar
               else
                  call env%warning("Unknown reference state '"//sec//"'", source)
               end if
            end if
         else
            call env%error("No solvent name provided for GBSA", source)
         end if

      case('--cosmo')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_gbsa(env, 'solvent', sec)
            call set_gbsa(env, 'cosmo', 'true')
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (sec == 'reference') then
                  gsolvstate = 1
               else if (sec == 'bar1M') then
                  gsolvstate = 2
               else
                  call env%warning("Unknown reference state '"//sec//"'", source)
               end if
            end if
         else
            call env%error("No solvent name provided for COSMO", source)
         end if

      case('--scc', '--sp')
         call set_runtyp('scc')

      case('--vip')
         call set_gfn(env,'method','1')
         call set_runtyp('vip')

      case('--vea')
         call set_gfn(env,'method','1')
         call set_runtyp('vea')

      case('--vipea')
         call set_gfn(env,'method','1')
         call set_runtyp('vipea')

      case('--vomega')
         call set_gfn(env,'method','1')
         call set_runtyp('vomega')

      case('--vfukui')
         call set_runtyp('vfukui')

      case('--grad')
         call set_runtyp('grad')
         lgrad = .true.

      case('-o', '--opt')
         call set_runtyp('opt')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_opt(env,'optlevel',sec)
         endif

      case('--hess')
         call set_runtyp('hess')

      case('--md')
         call set_runtyp('md')

      case('--ohess')
         call set_runtyp('ohess')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_opt(env,'optlevel',sec)
         endif
      
      case('--bhess')
         call set_runtyp('bhess')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_opt(env,'optlevel',sec)
         endif

      case('--omd')
         call set_runtyp('omd')
         call set_opt(env,'optlevel','-1')

      case('--siman')
         call set_runtyp('siman')
         call set_md(env,'nvt','true')

      case('--path')
         call set_runtyp('path')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_path(env,'product',sec)
         end if

      case('--screen')
         call set_runtyp('screen')

      case('--gmd')
         call set_runtyp('gmd')
         call env%error("This feature has been deprecated, I'm sorry.", source)

      case('--modef')
         call set_runtyp('modef')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_modef(env,'mode',sec)
         end if

      case('--mdopt')
         call set_runtyp('mdopt')

      case('--metadyn')
         call set_runtyp('md')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_metadyn(env,'save',sec)
         end if
         call set_metadyn(env,'static','false')

      case('--metaopt')
         call set_runtyp('metaopt')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_opt(env,'optlevel',sec)
         end if

      case('--nat')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_natom(env,sec)
         end if

      case('--bias-input', '--gesc')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_metadyn(env, 'bias-input', sec)
         else
            call env%error("No input file for RMSD bias provided", source)
         end if

      case('--wrtopo')
         call args%nextArg(sec)
         if (allocated(sec)) then
           call setWRtopo(sec,printTopo)
           if(printTopo%warning) call env%error("A wrtopo argument has been misspelled.",source)
         else
           call env%error("The wrtopo keyword is missing an argument.",source)
         endif
      end select
      call args%nextFlag(flag)
   end do

end subroutine parseArguments

function read_whole_file(fname) result(list)
   character(len=*), intent(in) :: fname
   character(len=:), allocatable :: list
   integer :: io, stat
   character(len=:), allocatable :: line
   open(newunit=io, file=fname, iostat=stat)
   call getline(io, list, stat)
   do while(stat == 0)
      call getline(io, line, stat)
      if (stat == 0) list = list // "," // line
   end do
   close(io, iostat=stat)
end function read_whole_file

! set booleans for requested topology list printout
subroutine setWRtopo(sec,printTopo)
   ! command line argument
   character(len=*), intent(in) :: sec
   ! type holds booleans of to be printed topology lists
   type(TPrintTopo), intent(inout) :: printTopo
   ! seperator for lists is ","
   character, parameter :: sep = ","
   ! current and old position of seperator
   integer :: curr_pos, old_pos
   integer :: lenSec, i

   curr_pos = 0
   old_pos = 0
   lenSec = len(sec)
   do i=1, lenSec
     curr_pos = scan(sec(curr_pos+1:lenSec),sep)+old_pos
     if(curr_pos.ne.old_pos) then
       call selectList(sec(old_pos+1:curr_pos-1),printTopo)
     else
       call selectList(sec(old_pos+1:lenSec),printTopo)
       exit
     endif
     old_pos=curr_pos
   enddo

end subroutine setWRtopo

subroutine selectList(secSplit, printTopo)
   ! part of command line argument
   character(len=*), intent(in) :: secSplit
   ! holds booleans of to be printed topology lists
   type(TPrintTopo), intent(inout) :: printTopo

   select case(secSplit)
   case("nb")
     printTopo%nb = .true.
   case("bpair")
     printTopo%bpair = .true.
   case("alist")
     printTopo%alist = .true.
   case("blist")
     printTopo%blist = .true.
   case("tlist")
     printTopo%tlist = .true.
   case("vtors")
     printTopo%vtors = .true.
   case("vbond")
     printTopo%vbond = .true.
   case("vangl")
     printTopo%vangl = .true.
   case("hbbond")
      printTopo%hbbond = .true.
   case("eeq")
      printTopo%eeq = .true.
   case default
     printTopo%warning = .true.
   end select
end subroutine selectList

end module xtb_prog_main
