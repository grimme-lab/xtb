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
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_io, only: stderr
   use xtb_mctc_timings
   use xtb_mctc_systools
   use xtb_mctc_convert
   use xtb_mctc_param
   use xtb_type_molecule
   use xtb_type_calculator
   use xtb_type_restart
   use xtb_tblite_restart, only: loadRestart, dumpRestart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment, only: TEnvironment, init
   use xtb_prog_argparser
   use xtb_solv_state
   use xtb_setparam
   use xtb_sphereparam
   use xtb_scanparam
   use xtb_splitparam
   use xtb_fixparam
   use xtb_features, only: get_xtb_feature
   use xtb_constrain_param, only: read_userdata
   use xtb_shake, only: init_shake
   use xtb_gfnff_shake, only: gff_init_shake => init_shake
   use xtb_embedding, only: init_pcem
   use xtb_io_reader, only: readMolecule
   use xtb_io_writer, only: writeMolecule
   use xtb_mctc_filetypes, only: fileType, getFileType, generateFileMetaInfo, &
      & generateFileName
   use xtb_readin
   use xtb_printout
   use xtb_setmod
   use xtb_propertyoutput
   use xtb_io_writer_turbomole, only: writeResultsTurbomole
   use xtb_io_writer_orca, only: writeResultsOrca
   use xtb_io_writer_gaussian, only: writeResultsGaussianExternal
   use xtb_restart
   use xtb_readparam
   use xtb_scc_core, only: iniqshell
   use xtb_aespot, only: get_radcn
   use xtb_iniq, only: iniqcn
   use xtb_eeq
   use xtb_disp_ncoord, only: ncoord_gfn, dncoord_erf, dncoord_d3, ncoord_erf, &
      & ncoord_d3
   use xtb_basis
   use xtb_axis, only: axis3
   use xtb_hessian, only: numhess
   use xtb_dynamic, only: md
   use xtb_modef, only: modefollow
   use xtb_mdoptim, only: mdopt
   use xtb_screening, only: screen
   use xtb_xtb_calculator
   use xtb_gfnff_calculator
   use xtb_iff_calculator, only: TIFFCalculator
   use xtb_paramset
   use xtb_xtb_gfn0
   use xtb_xtb_gfn1
   use xtb_xtb_gfn2
   use xtb_main_setup
   use xtb_main_defaults, only: initDefaults
   use xtb_main_json, only: main_xtb_json, write_json_gfnff_lists
   use xtb_geoopt
   use xtb_metadynamic
   use xtb_biaspath
   use xtb_coffee
   use xtb_disp_dftd3param
   use xtb_disp_dftd4
   use xtb_gfnff_param, only: gff_print
   use xtb_gfnff_topology, only: TPrintTopo
   use xtb_gfnff_convert, only: struc_convert
   use xtb_scan
   use xtb_kopt
   use xtb_iff_iffprepare, only: prepare_IFF
   use xtb_iff_data, only: TIFFData
   use xtb_oniom, only: oniom_input, TOniomCalculator, calculateCharge
   use xtb_vertical, only: vfukui
   use xtb_tblite_calculator, only: TTBLiteCalculator, TTBLiteInput, &
         & newTBLiteWavefunction, get_ceh, num_grad_chrg
   use xtb_ptb_calculator, only: TPTBCalculator
   use xtb_solv_cpx, only: TCpcmx
   use xtb_dipro, only: get_jab, jab_input
   !> PTB related modules
   use xtb_main_json, only: main_ptb_json

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
      class(TCalculator), allocatable :: calc, cpxcalc
      type(freq_results) :: fres
      type(TRestart) :: chk
      type(chrg_parameter) :: chrgeq
      type(TIFFData), allocatable :: iff_data
      type(oniom_input) :: oniom
      type(jab_input) :: dipro
      type(TCpcmx) :: cpx
      type(TTBLiteInput) :: tblite
!  store important names and stuff like that in FORTRAN strings
      character(len=:), allocatable :: fname    ! geometry input file
      character(len=:), allocatable :: xcontrol ! instruction file
      character(len=:), allocatable :: xrc      ! global instruction file
      character(len=:), allocatable :: fnv      ! parameter file
      character(len=:), allocatable :: tmpname  ! temporary string
      character(len=:), allocatable :: cdum     ! temporary string
      character(len=:), allocatable :: extension, basename, directory
      integer :: ftype

!! ========================================================================
!  default names for important files in xtb
      character(len=*), parameter :: p_fname_rc = '.xtbrc'
      character(len=*), parameter :: p_fname_param_gfn0 = 'param_gfn0-xtb.txt'
      character(len=*), parameter :: p_fname_param_gfn1 = 'param_gfn1-xtb.txt'
      character(len=*), parameter :: p_fname_param_gfn2 = 'param_gfn2-xtb.txt'
      character(len=*), parameter :: p_fname_param_gfnff = '.param_gfnff.xtb'
      character(len=*), parameter :: p_fname_param_ipea = 'param_ipea-xtb.txt'
      character(len=*), parameter :: p_fname_param_ptb = 'param_ptb.txt'

      integer :: gsolvstate
      integer :: i, j, k, l, idum
      integer :: ich, ictrl, iprop ! file handle
      real(wp) :: sigma(3, 3)
      real(wp), allocatable :: cn(:)
      real(wp), allocatable :: sat(:)
      real(wp), allocatable :: g(:, :)
      real(wp), allocatable :: fukui(:, :)
      real(wp) :: vec3(3)
      type(TxTBParameter) :: globpar
      real(wp), allocatable :: dcn(:, :, :)
      real(wp), allocatable :: dq(:, :, :)
      real(wp), allocatable :: dumdumdum(:, :, :)
      real(wp), allocatable :: q(:)
      real(wp), allocatable :: ql(:)
      real(wp), allocatable :: qr(:)
!! ------------------------------------------------------------------------
      integer, external :: ncore

!! ------------------------------------------------------------------------
      logical :: struc_conversion_done = .false.
      logical :: anyopt, anyhess

!! ========================================================================
!  debugging variables for numerical gradient
      logical, parameter :: gen_param = .false.
      logical, parameter :: debug = .false.
      type(TRestart) :: wf0
      real(wp), allocatable :: coord(:, :), numg(:, :), gdum(:, :)
      real(wp) :: sdum(3, 3), nums(3, 3), eps(3, 3), latt(3, 3)

      real(wp), parameter :: step = 0.00001_wp, step2 = 0.5_wp / step ! for numerical gradient
      real(wp), parameter :: sstep = 1.0_wp * 10.0_wp**(-6), sstep2 = 0.5_wp / sstep ! for numerical sigma
      real(wp) :: er, el
      logical :: coffee ! if debugging gets really though, get a coffee

!! ------------------------------------------------------------------------
!  undocumented and unexplainable variables go here
      integer :: nFiles, iFile
      integer :: rohf, err
      real(wp) :: dum5, egap, etot, ipeashift
      real(wp) :: zero, t0, t1, w0, w1, etot2, g298
      real(wp) :: one, two
      real(wp) :: ea, ip
      real(wp) :: vomega
      real(wp) :: energy_gas
      parameter(zero=0.0_wp)
      parameter(one=1.0_wp)
      parameter(two=2.0_wp)
      logical :: ex, okbas
      logical :: epr, diff, murks
      logical :: exist
      logical :: lgrad, restart
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

      call parseArguments(env, argParser, xcontrol, fnv, lgrad, &
         & restart, gsolvstate, strict, copycontrol, coffee, printTopo, oniom, dipro, tblite)

      
      ! TEMPORARY: no solvation available for PTB and tblite !
      if (set%mode_extrun == p_ext_tblite .or. set%mode_extrun == p_ext_ptb) then
         if (allocated(set%solvInput%solvent)) then
            call env%error("Solvation is not implemented for PTB/tblite", source)
         endif
      end if

      !> Spin-polarization is only available in the tblite library
      if (set%mode_extrun /= p_ext_tblite .and. tblite%spin_polarized) then
         call env%error("Spin-polarization is only available with the tblite library! Try --tblite", source)
      end if

      !> If hessian (or ohess or bhess) is requested in combination with PTB, conduct GFN2-xTB + PTB hessian
      anyhess = (set%runtyp == p_run_hess) .or. (set%runtyp == p_run_ohess) .or. (set%runtyp == p_run_bhess)
      if (anyhess) then
         if(set%mode_extrun == p_ext_ptb) then
            set%mode_extrun = p_ext_xtb
            set%ptbsetup%ptb_in_hessian = .true.
            call set_gfn(env, 'method', '2')
            call set_gfn(env, 'd4', 'true')
            tblite%method = "gfn2"
            set%ptbsetup%hessmethod = "GFN2-xTB"
         else
            if (set%elprop == p_elprop_alpha) then
               call env%error("Raman activities are not implemented with this method",source)
            end if
         endif
      endif 

      nFiles = argParser%countFiles()
      select case (nFiles)
      case (0)
         if (.not. coffee) then
          if (printTopo%warning) call env%error("Eventually the input file was given to wrtopo as an argument.", source)
            call env%error("No input file given, so there is nothing to do", source)
         else
            fname = 'coffee'
         end if
      case (1:)
         do iFile = 1, nFiles - 1
            call argParser%nextFile(fname)
            call env%warning("Input file '"//fname//"' will be ignored", source)
         end do
         call argParser%nextFile(fname)
      end select

      if (.not. allocated(xcontrol)) then
         if (copycontrol) then
            xcontrol = 'xtb.inp'
         else
            xcontrol = fname
         end if
      end if

      anyopt = ((set%runtyp == p_run_opt) .or. (set%runtyp == p_run_ohess) .or. &
         &   (set%runtyp == p_run_omd) .or. (set%runtyp == p_run_screen) .or. &
         &   (set%runtyp == p_run_metaopt))

      if (allocated(set%solvInput%cpxsolvent).and.anyopt) then
            call env%terminate("CPCM-X not implemented for geometry optimization. &
                 & Please use another solvation model for optimization instead.")
      endif     
 
      if ((set%mode_extrun == p_ext_ptb) .and. anyopt) call env%terminate("PTB not implemented for geometry optimization. &
         &Please use another method for optimization instead.")

      call env%checkpoint("Command line argument parsing failed")

      ! ------------------------------------------------------------------------
      !> read the detailed input file
      call rdcontrol(xcontrol, env, copy_file=copycontrol)

      call env%checkpoint("Reading '"//xcontrol//"' failed")

      ! ------------------------------------------------------------------------
      !> read dot-Files before reading the rc and after reading the xcontrol
      !> Total molecular charge
      call open_file(ich, '.CHRG', 'r')
      if (ich /= -1) then
         call getline(ich, cdum, iostat=err)
         if (err /= 0) then
            call env%error('.CHRG is empty!', source)
         else
            call set_chrg(env, cdum)
            call close_file(ich)
         end if
      end if

      call env%checkpoint("Reading charge from file failed")

      !> Number of unpaired electrons
      call open_file(ich, '.UHF', 'r')
      if (ich /= -1) then
         call getline(ich, cdum, iostat=err)
         if (err /= 0) then
            call env%error('.UHF is empty!', source)
         else
            call set_spin(env, cdum)
            call close_file(ich)
         end if
      end if

      call env%checkpoint("Reading multiplicity from file failed")

      !> efield read: gfnff and PTB only
      if (set%mode_extrun == p_ext_gfnff .or. set%mode_extrun == p_ext_ptb) then
         call open_file(ich, '.EFIELD', 'r')
         if (ich /= -1) then
            call getline(ich, cdum, iostat=err)
            if (err /= 0) then
               call env%error('.EFIELD is empty!', source)
            else
               call set_efield(env, cdum)
               call close_file(ich)
            end if
         end if
      end if

      !> If EFIELD is not zero when using xtb, print a warning
      if (((set%mode_extrun /= p_ext_ptb) .and. (set%mode_extrun /= p_ext_gfnff)) &
         & .and. (sum(abs(set%efield)) /= 0.0_wp)) then
         call env%terminate("External electric field is not zero ('--efield' or file '.EFIELD'), &
            & but only supported for GFN-FF and PTB")
      end if

      ! ------------------------------------------------------------------------
      !> read the xtbrc if you can find it (use rdpath directly instead of xfind)
      call rdpath(env%xtbpath, p_fname_rc, xrc, exist)
      if (exist) then
         call rdcontrol(xrc, env, copy_file=.false.)

         call env%checkpoint("Reading '"//xrc//"' failed")
      end if

      ! ------------------------------------------------------------------------
      !> FIXME: some settings that are still not automatic
      !> Make sure GFN0-xTB uses the correct exttyp
      if (set%gfn_method == 0) call set_exttyp('eht')
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
            if (set%clichrg) then
               call env%warning("Charge in sdf/mol input was overwritten", source)
            else
               set%ichrg = nint(mol%chrg)
            end if
         end if

         call env%checkpoint("reading geometry input '"//fname//"' failed")
      end if

      ! ------------------------------------------------------------------------
      !> initialize the global storage
      call init_fix(mol%n)
      call init_split(mol%n)
      call init_constr(mol%n, mol%at)
      call init_scan
      call init_walls
      call init_pcem
      if (set%runtyp == p_run_bhess) then
         call init_bhess(mol%n)
      else
         call init_metadyn(mol%n, metaset%maxsave)
      end if
      !> Initialize the atomic masses with the physical constants
      atmass = atomic_mass(mol%at) * autoamu ! from splitparam.f90
      call load_rmsdbias(rmsdset, mol%n, mol%at, mol%xyz)
      ! ------------------------------------------------------------------------
      !> CONSTRAINTS & SCANS
      !> now we are at a point that we can check for requested constraints
      call read_userdata(xcontrol, env, mol)
      ! ------------------------------------------------------------------------
      !> get some memory
      allocate (cn(mol%n), sat(mol%n), g(3, mol%n), source=0.0_wp)
      set%periodic = mol%npbc > 0
      if (mol%npbc == 0) then
         if (set%do_cma_trafo) then
            allocate (coord(3, mol%n), source=0.0_wp)
            call axis3(1, mol%n, mol%at, mol%xyz, coord, vec3)
            mol%xyz = coord
            deallocate (coord)
         end if
      end if

      do i = 1, mol%n
         mol%z(i) = mol%at(i) - ncore(mol%at(i))
         ! lanthanides without f are treated as La
         if (mol%at(i) > 57 .and. mol%at(i) < 72) mol%z(i) = 3
      end do

      !> initialize time step for MD if requested autocomplete
      if (set%tstep_md < 0.0_wp) then
         set%tstep_md = (minval(atmass) / (atomic_mass(1) * autoamu))**(1.0_wp / 3.0_wp)
      end if

      mol%chrg = real(set%ichrg, wp)
   !! To assign charge
      mol%uhf = set%nalphabeta
      call initrand

      call setup_summary(env%unit, mol%n, fname, xcontrol, chk%wfn, xrc)

      ! ------------------------------------------------------------------------
      !> 2D => 3D STRUCTURE CONVERTER
      ! ------------------------------------------------------------------------
      if (mol%info%two_dimensional) then
         call struc_convert(env, restart, mol, chk, egap, set%etemp, set%maxscciter, &
                           &  set%optset%maxoptcycle, etot, g, sigma)
         struc_conversion_done = .true.
         mol%info%two_dimensional = .false.
      end if


      if (sum(abs(set%efield)) /= 0.0_wp) then
         write (env%unit, '(/,3x,a)') "--------------------------------------"
         write (env%unit, '(3x,a)') "--- external electric field / a.u. ---"
         write (env%unit, '(3x,3(a,f8.4))') "x = ", set%efield(1), " y = ", set%efield(2), " z = ", set%efield(3)
         write (env%unit, '(3x,a)') "--------------------------------------"
      end if

      !> initialize metadynamics
      call load_metadynamic(metaset, mol%n, mol%at, mol%xyz)

      !> restraining potential
      if (allocated(potset%xyz)) then
         if (lconstr_all_bonds) call constrain_all_bonds(mol%n, mol%at, potset%xyz)
         if (lconstr_all_angles) call constrain_all_angles(mol%n, mol%at, potset%xyz)
         if (lconstr_all_torsions) call constrain_all_torsions(mol%n, mol%at, potset%xyz)
         call setup_constrain_pot(mol%n, mol%at, potset%xyz)
      else
         if (lconstr_all_bonds) call constrain_all_bonds(mol%n, mol%at, mol%xyz)
         if (lconstr_all_angles) call constrain_all_angles(mol%n, mol%at, mol%xyz)
         if (lconstr_all_torsions) call constrain_all_torsions(mol%n, mol%at, mol%xyz)
         call setup_constrain_pot(mol%n, mol%at, mol%xyz)
      end if
      !  fragmentation for CMA constrain
      if (iatf1 == 0 .and. iatf2 == 0) then
         call ncoord_erf(mol%n, mol%at, mol%xyz, cn)
         call splitm(mol%n, mol%at, mol%xyz, cn)
      end if
      call splitprint(mol%n, mol%at, mol%xyz)

      if (set%verbose) then
         call fix_info(env%unit, mol%n, mol%at, mol%xyz)
         call pot_info(env%unit, mol%n, mol%at, mol%xyz)
      end if

      ! ------------------------------------------------------------------------
      !> write copy of detailed input
      if (copycontrol) then
         call open_set(ictrl, xcontrol)
         call write_set(ictrl)
         call close_set(ictrl)
      end if

      ! ------------------------------------------------------------------------
      !> if you have requested a define we stop here...
      if (set%define) then
         if (set%verbose) call main_geometry(env%unit, mol)
         call eval_define(set%veryverbose)
      end if
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
      end if

      !> check if someone is still using GFN3...
      if (set%gfn_method == 3) then
         call env%terminate('Wait for some months - for now, please use gfn_method=2!')
      end if

      ! ------------------------------------------------------------------------
      !> Print the method header and select the parameter file

      if (.not. allocated(fnv)) then
         select case (set%runtyp)
         case default
            call env%terminate('This is an internal error, please define your runtypes!')
         case (p_run_prescc,p_run_scc, p_run_grad, p_run_opt, p_run_hess, p_run_ohess, p_run_bhess, &
               p_run_md, p_run_omd, p_run_path, p_run_screen, &
               p_run_modef, p_run_mdopt, p_run_metaopt)
            if (set%mode_extrun == p_ext_gfnff) then
               fnv = xfind(p_fname_param_gfnff)
             elseif (set%mode_extrun .eq. p_ext_ptb) then
                  fnv = xfind(p_fname_param_ptb)
            else
               if (set%gfn_method == 0) then
                  fnv = xfind(p_fname_param_gfn0)
               end if
               if (set%gfn_method == 1) then
                  fnv = xfind(p_fname_param_gfn1)
               end if
               if (set%gfn_method == 2) then
                  fnv = xfind(p_fname_param_gfn2)
               end if
            end if
         case (p_run_vip, p_run_vea, p_run_vipea, p_run_vfukui, p_run_vomega)
            if (set%gfn_method == 0) then
               fnv = xfind(p_fname_param_gfn0)
            end if
            if (set%gfn_method == 1) then
               fnv = xfind(p_fname_param_gfn1)
            end if
            if (set%gfn_method == 2) then
               fnv = xfind(p_fname_param_gfn2)
            end if
         end select
      end if

      !-------------------------------------------------------------------------
      !> Perform a precomputation of electronic properties for xTB-IFF
      if (set%mode_extrun == p_ext_iff) then
         allocate (iff_data)
         call prepare_IFF(env, mol, iff_data)
         call env%checkpoint("Could not generate electronic properties")
      end if

      ! ------------------------------------------------------------------------
      !> Obtain the parameter data
      call newCalculator(env, mol, calc, fnv, restart, set%acc, oniom, iff_data, tblite)
      call env%checkpoint("Could not setup single-point calculator")

      call initDefaults(env, calc, mol, gsolvstate)
      call env%checkpoint("Could not setup defaults")

      ! ------------------------------------------------------------------------
      !> initial guess, setup wavefunction
      select type (calc)
      type is (TxTBCalculator)
         call chk%wfn%allocate(mol%n, calc%basis%nshell, calc%basis%nao)

         ! Make sure number of electrons is initialized an multiplicity is consistent
         chk%wfn%nel = nint(sum(mol%z) - mol%chrg)
         chk%wfn%nopen = mol%uhf
         if (chk%wfn%nopen == 0 .and. mod(chk%wfn%nel, 2) /= 0) chk%wfn%nopen = 1

         !> EN charges and CN
         if (set%gfn_method < 2) then
            call ncoord_d3(mol%n, mol%at, mol%xyz, cn)
         else
            call ncoord_gfn(mol%n, mol%at, mol%xyz, cn)
         end if
         if (mol%npbc > 0) then
            chk%wfn%q = real(set%ichrg, wp) / real(mol%n, wp)
         else
            if (set%guess_charges == p_guess_gasteiger) then
               call iniqcn(mol%n, mol%at, mol%z, mol%xyz, set%ichrg, 1.0_wp, chk%wfn%q, cn, set%gfn_method, .true.)
            else if (set%guess_charges == p_guess_goedecker) then
               call ncoord_erf(mol%n, mol%at, mol%xyz, cn)
               call goedecker_chrgeq(mol%n, mol%at, mol%xyz, real(set%ichrg, wp), cn, dcn, chk%wfn%q, dq, er, g, &
                                     .false., .false., .false.)
            else
               call ncoord_gfn(mol%n, mol%at, mol%xyz, cn)
               chk%wfn%q = real(set%ichrg, wp) / real(mol%n, wp)
            end if
         end if
         !> initialize shell charges from gasteiger charges
         call iniqshell(calc%xtbData, mol%n, mol%at, mol%z, calc%basis%nshell, chk%wfn%q, chk%wfn%qsh, set%gfn_method)
      type is (TTBLiteCalculator)
         call newTBLiteWavefunction(env, mol, calc, chk)
      end select

      ! get CEH charges !
      if (allocated(tblite%ceh)) then
         ! if numercal charges are requested !
         if (tblite%ceh%grad) then
            call num_grad_chrg(env,mol,tblite)
         ! only charges !
         else
            call get_ceh(env,mol,tblite)
         endif 

         ! stop calculation !
         call finalize_xtb(env)
      end if

      ! ------------------------------------------------------------------------
      !> printout a header for the exttyp
      call calc%writeInfo(env%unit, mol)

      call delete_file('.sccnotconverged')

      call env%checkpoint("Setup for calculation failed")

      select type (calc)
      type is (TxTBCalculator)
         if (restart .and. calc%xtbData%level /= 0) then ! only in first run
            call readRestart(env, chk%wfn, 'xtbrestart', mol%n, mol%at, set%gfn_method, exist, .true.)
         end if
         ipeashift = calc%xtbData%ipeashift
      type is (TTBLiteCalculator)
         if (restart) then
            call loadRestart(env, chk, 'xtbrestart', exist)
            if (exist) write (env%unit, "(a)") "Wavefunction read from restart file"
         end if
      type is (TOniomCalculator)
         select type (xtb => calc%real_low)
         type is (TxTBCalculator)
            call chk%wfn%allocate(mol%n, xtb%basis%nshell, xtb%basis%nao)
            call newWavefunction(env, mol, xtb, chk)
         !! assigns only partial charges q and shell charges
            if (restart) then ! only in first run
               call readRestart(env, chk%wfn, 'xtbrestart', mol%n, mol%at, set%gfn_method, exist, .true.)
            end if
         end select
         if (.not. set%oniom_settings%fixed_chrgs) then
            set%oniom_settings%innerchrg = calculateCharge(calc, env, mol, chk)
         end if

      end select
      !-------------------------------------------------------------------------
      !> DIPRO calculation of coupling integrals for dimers
      if (dipro%diprocalc) then
         call start_timing(11)
         call get_jab(env, tblite, mol, splitlist, dipro)
         call env%checkpoint("Something in your DIPRO calculation went wrong.")
         call stop_timing_run
         call stop_timing(11)
         write (*, '(A)') "----------------------------------------------------------"
         call prdate('E')
         write (*, '(A)') "----------------------------------------------------------"
         call prtiming(11, 'dipro')
         call terminate(0)
      end if

      ! ========================================================================
      !> the SP energy which is always done
      call start_timing(2)
      call calc%singlepoint(env, mol, chk, 2, exist, etot, g, sigma, egap, res)
      call stop_timing(2)
      select type (calc)
      type is (TGFFCalculator)
         gff_print = .false.
      end select
      call env%checkpoint("Single point calculation terminated")

      !> write 2d => 3d converted structure
      if (struc_conversion_done) then
         call generateFileName(tmpname, 'gfnff_convert', extension, mol%ftype)
         write (env%unit, '(10x,a,1x,a,/)') &
            "converted geometry written to:", tmpname
         call open_file(ich, tmpname, 'w')
         call writeMolecule(mol, ich, energy=res%e_total, gnorm=res%gnorm)
         call close_file(ich)
      end if

      ! ========================================================================
      !> determine kopt for bhess including final biased geometry optimization
      if (set%runtyp == p_run_bhess) then
         call set_metadynamic(metaset, mol%n, mol%at, mol%xyz)
         call get_kopt(metaset, env, restart, mol, chk, calc, egap, set%etemp, set%maxscciter, &
            & set%optset%maxoptcycle, set%optset%optlev, etot, g, sigma, set%acc)
      end if

      ! ------------------------------------------------------------------------
      !> numerical gradient for debugging purposes
      if (debug) then
         !  generate a warning to keep release versions from calculating numerical gradients
         call env%warning('XTB IS CALCULATING NUMERICAL GRADIENTS, RESET DEBUG FOR RELEASE!')
         print'(/,"analytical gradient")'
         print *, g
         allocate (coord(3, mol%n), source=mol%xyz)
         allocate (numg(3, mol%n), gdum(3, mol%n), source=0.0_wp)
         wf0 = chk
         do i = 1, mol%n
            do j = 1, 3
               mol%xyz(j, i) = mol%xyz(j, i) + step
               chk = wf0
               call calc%singlepoint(env, mol, chk, 0, .true., er, gdum, sdum, egap, res)
               mol%xyz(j, i) = mol%xyz(j, i) - 2 * step
               chk = wf0
               call calc%singlepoint(env, mol, chk, 0, .true., el, gdum, sdum, egap, res)
               mol%xyz(j, i) = mol%xyz(j, i) + step
               numg(j, i) = step2 * (er - el)
            end do
         end do
         print'(/,"numerical gradient")'
         print *, numg
         print'(/,"difference gradient")'
         print *, g - numg
         deallocate (coord)
      end if

      !> numerical sigma (=volume*stressTensor) for debugging purposes
      if (debug .and. mol%npbc == 3) then
         !  generate a warning to keep release versions from calculating numerical gradients
         call env%warning('XTB IS CALCULATING NUMERICAL STRESS, RESET DEBUG FOR RELEASE!')
         print'(/,"analytical sigma")'
         print *, sigma
         if (.not. allocated(gdum)) allocate (gdum(3, mol%n), source=0.0_wp)
         allocate (coord(3, mol%n), source=mol%xyz)
         latt = mol%lattice
         nums = 0.0_wp
         !sdum = 0.0_wp
         wf0 = chk
         do j = 1, 3
            do i = 1, j
               ! Only eps_ij=step the rest equals zero: eps_(kl.ne.ij)=0
               eps = 0.0_wp
               eps(i, j) = sstep
               ! adjust position vectors and lattice to get er
               do k = 1, 3
                  mol%xyz(k, :) = kron(k, 1) * mol%xyz(1, :) + eps(k, 1) * mol%xyz(1, :) + &
                                &   kron(k, 2) * mol%xyz(2, :) + eps(k, 2) * mol%xyz(2, :) + &
                                &   kron(k, 3) * mol%xyz(3, :) + eps(k, 3) * mol%xyz(3, :)
                  mol%lattice(k, 1) = (kron(k, 1) + eps(k, 1)) * mol%lattice(1, 1) + &
                                   & (kron(k, 2) + eps(k, 2)) * mol%lattice(2, 1) + &
                                   & (kron(k, 3) + eps(k, 3)) * mol%lattice(3, 1)
                  mol%lattice(k, 2) = (kron(k, 1) + eps(k, 1)) * mol%lattice(1, 2) + &
                                   & (kron(k, 2) + eps(k, 2)) * mol%lattice(2, 2) + &
                                   & (kron(k, 3) + eps(k, 3)) * mol%lattice(3, 2)
                  mol%lattice(k, 3) = (kron(k, 1) + eps(k, 1)) * mol%lattice(1, 3) + &
                                   & (kron(k, 2) + eps(k, 2)) * mol%lattice(2, 3) + &
                                   & (kron(k, 3) + eps(k, 3)) * mol%lattice(3, 3)
               end do
               chk = wf0
               call calc%singlepoint(env, mol, chk, 0, .true., er, gdum, sdum, egap, res)

               ! reset coordinates and lattice
               mol%xyz = coord
               mol%lattice = latt
               ! adjust position vectors and lattice to get el
               do k = 1, 3
                  mol%xyz(k, :) = kron(k, 1) * mol%xyz(1, :) - eps(k, 1) * mol%xyz(1, :) + &
                                &  kron(k, 2) * mol%xyz(2, :) - eps(k, 2) * mol%xyz(2, :) + &
                                &  kron(k, 3) * mol%xyz(3, :) - eps(k, 3) * mol%xyz(3, :)
                  mol%lattice(k, 1) = (kron(k, 1) - eps(k, 1)) * mol%lattice(1, 1) + &
                                   & (kron(k, 2) - eps(k, 2)) * mol%lattice(2, 1) + &
                                   & (kron(k, 3) - eps(k, 3)) * mol%lattice(3, 1)
                  mol%lattice(k, 2) = (kron(k, 1) - eps(k, 1)) * mol%lattice(1, 2) + &
                                   & (kron(k, 2) - eps(k, 2)) * mol%lattice(2, 2) + &
                                   & (kron(k, 3) - eps(k, 3)) * mol%lattice(3, 2)
                  mol%lattice(k, 3) = (kron(k, 1) - eps(k, 1)) * mol%lattice(1, 3) + &
                                   & (kron(k, 2) - eps(k, 2)) * mol%lattice(2, 3) + &
                                   & (kron(k, 3) - eps(k, 3)) * mol%lattice(3, 3)
               end do
               chk = wf0
               call calc%singlepoint(env, mol, chk, 0, .true., el, gdum, sdum, egap, res)

               ! numerical sigma (=volume*stressTensor)
               nums(i, j) = sstep2 * (er - el)  ! divide by 2 times step size
               nums(j, i) = nums(i, j)  ! stress tensor is symmetric
               ! reset coordinates and lattice
               mol%xyz = coord
               mol%lattice = latt
            end do
         end do

         print'(/,"numerical sigma")'
         print *, nums
         print'(/,"difference sigma")'
         print *, sigma - nums
         deallocate (coord)
      end if

!---------------------------------------------!
! Geometry optimization(ANCopt,L_ANCopt,FIRE) !
!---------------------------------------------!
      if (anyopt) then

         if (mol%npbc > 0 .and. (set%mode_extrun == p_ext_gfnff &
                 & .or. set%mode_extrun == p_ext_mcgfnff)) then  ! if(npbc)
            deallocate (set%opt_engine)
            call set_opt(env, 'engine', 'pbc_lbfgs')  ! use lbfgs
         end if

         ! start optimization timer !
         call start_timing(3)

         ! calculation !
         call geometry_optimization &
            &     (env, mol, chk, calc, &
        &      egap,set%etemp,set%maxscciter,set%optset%maxoptcycle,etot,g,sigma,set%optset%optlev,.true.,.false.,murks)

         ! save results !
         res%e_total = etot
         res%gnorm = norm2(g)

         ! constrained optimization !
         if (nscan > 0) then
            call relaxed_scan(env, mol, chk, calc)
         end if

         ! in case of failure cretae xtblast geometry !
         if (murks) then
            call generateFileName(tmpname, 'xtblast', extension, mol%ftype)
            write (env%unit, '(/,a,1x,a,/)') &
               "last geometry written to:", tmpname
            call open_file(ich, tmpname, 'w')
            call writeMolecule(mol, ich, energy=res%e_total, gnorm=res%gnorm)
            call close_file(ich)
            call env%terminate("Geometry optimization failed")
         end if

         ! stop optimization timer !
         call stop_timing(3)

      end if

      ! ------------------------------------------------------------------------
      !> automatic VIP and VEA single point (maybe after optimization)
      if (set%runtyp == p_run_vip .or. set%runtyp == p_run_vipea &
         & .or. set%runtyp == p_run_vomega) then
         call start_timing(2)
         call vip_header(env%unit)
         mol%chrg = mol%chrg + 1
         chk%wfn%nel = chk%wfn%nel - 1
         if (mod(chk%wfn%nel, 2) /= 0) chk%wfn%nopen = 1
         call calc%singlepoint(env, mol, chk, 1, exist, etot2, g, sigma, egap, res)
         ip = etot2 - etot - ipeashift
         write (env%unit, '(72("-"))')
         write (env%unit, '("empirical IP shift (eV):",f10.4)') &
            &                  autoev * ipeashift
         write (env%unit, '("delta SCC IP (eV):",f10.4)') autoev * ip
         write (env%unit, '(72("-"))')
         mol%chrg = mol%chrg - 1
         chk%wfn%nel = chk%wfn%nel + 1
         call stop_timing(2)
      end if

      if (set%runtyp == p_run_vea .or. set%runtyp == p_run_vipea &
         & .or. set%runtyp == p_run_vomega) then
         call start_timing(2)
         call vea_header(env%unit)
         mol%chrg = mol%chrg - 1
         chk%wfn%nel = chk%wfn%nel + 1
         if (mod(chk%wfn%nel, 2) /= 0) chk%wfn%nopen = 1
         call calc%singlepoint(env, mol, chk, 1, exist, etot2, g, sigma, egap, res)
         ea = etot - etot2 - ipeashift
         write (env%unit, '(72("-"))')
         write (env%unit, '("empirical EA shift (eV):",f10.4)') &
            &                  autoev * ipeashift
         write (env%unit, '("delta SCC EA (eV):",f10.4)') autoev * ea
         write (env%unit, '(72("-"))')

         mol%chrg = mol%chrg + 1
         chk%wfn%nel = chk%wfn%nel - 1
         call stop_timing(2)
      end if

      ! ------------------------------------------------------------------------
      !> vomega (electrophilicity) index
      if (set%runtyp == p_run_vomega) then
         write (env%unit, '(a)')
         write (env%unit, '(72("-"))')
         write (env%unit, '(a,1x,a)') &
            "Calculation of global electrophilicity index", &
            "(IP+EA)²/(8·(IP-EA))"
         vomega = (ip + ea)**2 / (8 * (ip - ea))
         write (env%unit, '("Global electrophilicity index (eV):",f10.4)') &
            autoev * vomega
         write (env%unit, '(72("-"))')
      end if

      ! ------------------------------------------------------------------------
      !> Fukui Index from Mulliken population analysis
      if (set%runtyp == p_run_vfukui) then
         allocate (fukui(3, mol%n))
         call vfukui(env, mol, chk, calc, fukui)
      end if

      ! ------------------------------------------------------------------------
      !> numerical hessian calculation
      if ((set%runtyp == p_run_hess) .or. (set%runtyp == p_run_ohess) .or. (set%runtyp == p_run_bhess)) then
         if (set%runtyp == p_run_bhess .and. set%mode_extrun /= p_ext_turbomole) then
            call generic_header(env%unit, "Biased Numerical Hessian", 49, 10)
         else if (set%runtyp == p_run_bhess .and. set%mode_extrun == p_ext_turbomole) then
            call generic_header(env%unit, "Biased Analytical TM Hessian", 49, 10)
         else if (set%mode_extrun == p_ext_turbomole) then
            call generic_header(env%unit, "Analytical TM Hessian", 49, 10)
         else
            call numhess_header(env%unit)
         end if
         if (mol%npbc > 0) then
            call env%error("Phonon calculations under PBC are not implemented", source)
         end if
         call start_timing(5)
         call numhess &
            &       (env, mol, chk, calc, &
            &        egap, set%etemp, set%maxscciter, etot, g, sigma, fres)
         call stop_timing(5)

         call env%checkpoint("Hessian calculation terminated")
      end if

      ! reset the gap, since it is currently not updated in ancopt and numhess
      if (allocated(chk%wfn%emo)) then
         res%hl_gap = chk%wfn%emo(chk%wfn%ihomo + 1) - chk%wfn%emo(chk%wfn%ihomo)
      end if

      !> CPCM-X post-SCF solvation
      if (allocated(calc%solvation)) then
         if (allocated(calc%solvation%cpxsolvent)) then
            select type (calc)
            type is (TxTBCalculator)
               call generic_header(env%unit, "CPCM-X post-SCF solvation evaluation", 49, 10)
               if (set%gfn_method /= 2) call env%warning("CPCM-X was parametrized for GFN2-xTB. &
                  &The results are probably inaccurate with other methods.")
               Call cpx%setup(env, calc%solvation%cpxsolvent)
               Call env%checkpoint("CPCM-X setup terminated")
               cpxcalc = calc
               deallocate (cpxcalc%solvation)
               call cpxcalc%singlepoint(env, mol, chk, 1, .false., energy_gas, g, sigma, egap, res)
          Call cpx%calc_solv(env, calc%solvation%cpxsolvent, energy_gas, 0.4_wp, 298.15_wp, 500, 0.0001_wp, res%e_total)
               Call cpx%print(set%verbose)
               Call env%checkpoint("CPCM-X post-SCF solvation evaluation terminated")
            type is (TGFFCalculator)
               call env%error("CPCM-X is not possible with a force field.", source)
            end select
         end if
      end if

      call env%checkpoint("Calculation terminated")

      ! ========================================================================
      !> PRINTOUT SECTION
      if (allocated(set%property_file)) then
         call open_file(iprop, set%property_file, 'w')
         if (iprop == -1) then
            iprop = env%unit
            deallocate (set%property_file)
         else
            write (env%unit, '(/,a)') "Property printout bound to '"//set%property_file//"'"
            if (allocated(cdum)) deallocate (cdum)
            call get_command(length=l)
            allocate (character(len=l) :: cdum)
            call get_command(cdum)
            write (iprop, '("command:  ''",a,"''")') cdum
            call rdvar('HOSTNAME', cdum, err)
            if (err == 0) &
               write (iprop, '("hostname: ''",a,"''")') cdum
            write (iprop, '("date:     ",a)') prtimestring('S')
         end if
      else
         iprop = env%unit
      end if

      call generic_header(iprop, 'Property Printout', 49, 10)
      if (lgrad) then
         call writeResultsTurbomole(mol, energy=etot, gradient=g, sigma=sigma)
         if (allocated(basename)) then
            cdum = basename//'.engrad'
         else
            cdum = 'xtb-orca.engrad'
         end if
         call open_file(ich, cdum, 'w')
         call writeResultsOrca(ich, mol, etot, g)
         call close_file(ich)
      end if
      if (mol%ftype == fileType%gaussian) then
         if (allocated(basename)) then
            cdum = basename//'.EOu'
         else
            cdum = 'xtb-gaussian.EOu'
         end if
         call open_file(ich, cdum, 'w')
         call writeResultsGaussianExternal(ich, etot, res%dipole, g)
         call close_file(ich)
      end if

      if (set%periodic) then
         write (*, *) 'Periodic properties'
      else
         select type (calc)
         type is (TxTBCalculator)
            call main_property(iprop, env, mol, chk%wfn, calc%basis, calc%xtbData, res, &
               & calc%solvation, set%acc)
            call main_cube(set%verbose, mol, chk%wfn, calc%basis, res)
         type is (TGFFCalculator)
            call gfnff_property(iprop, mol%n, mol%xyz, calc%topo, chk%nlist)
         type is (TPTBCalculator)
            call ptb_property(iprop, env, chk, calc, mol, res)
         end select
      end if

      if (set%pr_json) then
         select type (calc)
         type is (TxTBCalculator)
            call open_file(ich, 'xtbout.json', 'w')
            call main_xtb_json(ich, &
                               mol, chk%wfn, calc%basis, res, fres)
            call close_file(ich)
         type is (TPTBCalculator)
            call open_file(ich, 'xtbout.json', 'w')
            call main_ptb_json(ich, &
                               mol, chk%wfn, calc, res, fres)
            call close_file(ich)
         end select
      end if
      if (printTopo%any()) then
         select type (calc)
         type is (TGFFCalculator)
            call write_json_gfnff_lists(mol%n, res%e_total, res%gnorm, calc%topo, calc%neigh, chk%nlist, printTopo)
         end select
      end if
      if ((set%runtyp == p_run_opt) .or. (set%runtyp == p_run_ohess) .or. &
          (set%runtyp == p_run_omd) .or. (set%runtyp == p_run_screen) .or. &
          (set%runtyp == p_run_metaopt) .or. (set%runtyp == p_run_bhess)) then
         call main_geometry(iprop, mol)
      end if

      if ((set%runtyp == p_run_hess) .or. (set%runtyp == p_run_ohess) .or. (set%runtyp == p_run_bhess)) then
         call generic_header(iprop, 'Frequency Printout', 49, 10)
         call main_freq(iprop, mol, chk%wfn, fres)
      end if

      if (allocated(set%property_file)) then
         if (iprop /= -1 .and. iprop /= env%unit) then
            call write_energy(iprop, res, fres, &
               & (set%runtyp == p_run_hess) .or. (set%runtyp == p_run_ohess) .or. (set%runtyp == p_run_bhess))
            call close_file(iprop)
         end if
      end if

      if ((set%runtyp == p_run_opt) .or. (set%runtyp == p_run_ohess) .or. &
          (set%runtyp == p_run_omd) .or. (set%runtyp == p_run_screen) .or. &
          (set%runtyp == p_run_metaopt) .or. (set%runtyp == p_run_bhess)) then
         call generateFileName(tmpname, 'xtbopt', extension, mol%ftype)
         write (env%unit, '(/,a,1x,a,/)') &
            "optimized geometry written to:", tmpname
         call open_file(ich, tmpname, 'w')
         call writeMolecule(mol, ich, energy=res%e_total, gnorm=res%gnorm)
         call close_file(ich)
      end if

      select type (calc)
      type is (TxTBCalculator)
         call write_energy(env%unit, res, fres, anyhess)
      type is (TOniomCalculator)
         call write_energy_oniom(env%unit, res, fres, anyhess)
      type is (TPTBCalculator)
      class default
         call write_energy_gff(env%unit, res, fres, anyhess)
      end select

      ! ------------------------------------------------------------------------
      !  xtb molecular dynamics
      if ((set%runtyp == p_run_md) .or. (set%runtyp == p_run_omd)) then
         if (metaset%maxsave > 0) then
            if (mol%npbc > 0) then
               call env%error("Metadynamic under PBC is not implemented", source)
            end if
            call metadyn_header(env%unit)
         else
            call md_header(env%unit)
         end if
         fixset%n = 0 ! no fixing for MD runs
         call start_timing(6)
         idum = 0
         select type (calc)
         class default
            if (set%shake_md) call init_shake(mol%n, mol%at, mol%xyz, chk%wfn%wbo)
         type is (TGFFCalculator)
            if (set%shake_md) call gff_init_shake(mol%n, mol%at, mol%xyz, calc%topo)
         end select
         call md &
            &     (env, mol, chk, calc, &
            &      egap, set%etemp, set%maxscciter, etot, g, sigma, 0, set%temp_md, idum)
         call stop_timing(6)
      end if

      ! ------------------------------------------------------------------------
      !  metadynamics
      if (set%runtyp == p_run_metaopt) then
         if (mol%npbc > 0) then
            call env%warning("Metadynamic under PBC is not implemented", source)
         end if
         call metadyn_header(env%unit)
         ! check if ANCOPT already convered
         if (murks) then
            call env%error('Optimization did not converge, aborting', source)
         end if
         write (env%unit, '(1x,"output written to xtbmeta.log")')
         call open_file(ich, 'xtbmeta.log', 'w')
         call writeMolecule(mol, ich, fileType%xyz, energy=etot, gnorm=norm2(g))
         k = metaset%nstruc + 1
         call start_timing(6)
         do l = k, metaset%maxsave
            metaset%nstruc = l
            metaset%xyz(:, :, metaset%nstruc) = mol%xyz
            ! randomize structure to avoid zero RMSD
            do i = 1, mol%n
               do j = 1, 3
                  call random_number(er)
                  mol%xyz(j, i) = mol%xyz(j, i) + 1.0e-6_wp * er
               end do
            end do
            call geometry_optimization &
               &     (env, mol, chk, calc, &
               &      egap, set%etemp, set%maxscciter, set%optset%maxoptcycle, etot, g, sigma, &
               &      set%optset%optlev, set%verbose, .true., murks)
            if (.not. set%verbose) then
               write (env%unit, '("current energy:",1x,f20.8)') etot
            end if
            if (murks) then
               call close_file(ich)
               write (env%unit, '(/,3x,"***",1x,a,1x,"***",/)') &
                  "FAILED TO CONVERGE GEOMETRY OPTIMIZATION"
               call touch_file('NOT_CONVERGED')
            end if
            call writeMolecule(mol, ich, fileType%xyz, energy=etot, gnorm=norm2(g))
         end do
         call close_file(ich)
         call stop_timing(6)
      end if

      ! ------------------------------------------------------------------------
      !  path finder
      if (set%runtyp == p_run_path) then
         call rmsdpath_header(env%unit)
         if (mol%npbc > 0) then
            call env%warning("Metadynamics under PBC are not implemented", source)
         end if
         call start_timing(4)
         call bias_path(env, mol, chk, calc, egap, set%etemp, set%maxscciter, etot, g, sigma)
         call stop_timing(4)
      end if

      ! ------------------------------------------------------------------------
      !  screen over input structures
      if (set%runtyp == p_run_screen) then
         call start_timing(8)
         call screen(env, mol, chk, calc, egap, set%etemp, set%maxscciter, etot, g, sigma)
         call stop_timing(8)
      end if

      ! ------------------------------------------------------------------------
      !  mode following for conformer search
      if (set%runtyp == p_run_modef) then
         if (mol%npbc > 0) then
            call env%warning("Modefollowing under PBC is not implemented", source)
         end if
         call start_timing(9)
         call modefollow(env, mol, chk, calc, egap, set%etemp, set%maxscciter, etot, g, sigma)
         call stop_timing(9)
      end if

      ! ------------------------------------------------------------------------
      !  optimize along MD from xtb.trj for conformer searches
      if (set%runtyp == p_run_mdopt) then
         call start_timing(10)
         call mdopt(env, mol, chk, calc, egap, set%etemp, set%maxscciter, etot, g, sigma)
         call stop_timing(10)
      end if

      ! ------------------------------------------------------------------------
      !  to further speed up xtb calculations we dump our most important
      !  quantities in a restart file, so we can save some precious seconds
      select type (calc)
      type is (TxTBCalculator)
         if (restart) then
            call writeRestart(env, chk%wfn, 'xtbrestart', set%gfn_method)
         end if
      type is (TTBLiteCalculator)
         if (restart) call dumpRestart(env, chk, 'xtbrestart')
      end select

      ! ------------------------------------------------------------------------
      !  we may have generated some non-fatal errors, which have been saved,
      !  so we should tell the user, (s)he may want to know what went wrong
      call env%show("Runtime exception occurred")
      call raise('F', 'Some non-fatal runtime exceptions were caught,'// &
         &           ' please check:')

      ! ------------------------------------------------------------------------
      !  print all files xtb interacted with while running (for debugging mainly)
      if (set%verbose) then
         write (env%unit, '(a)')
         write (env%unit, '(72("-"))')
         call print_filelist(env%unit)
      end if

      ! ------------------------------------------------------------------------
      !  make some post processing afterward, show some timings and stuff
      call finalize_xtb(env) 

   end subroutine xtbMain

!> Parse command line arguments and forward them to settings
   subroutine parseArguments(env, args, inputFile, paramFile, lgrad, &
         & restart, gsolvstate, strict, copycontrol, coffee, printTopo, oniom, dipro, tblite)

      use xtb_mctc_global, only: persistentEnv

      !> Name of error producer
      character(len=*), parameter :: source = "prog_main_parseArguments"

      !> Calculation environment
      type(TEnvironment) :: env

      !> Command line argument parser
      type(TArgParser) :: args

      !> Detailed input file name
      character(len=:), allocatable, intent(out) :: inputFile

      !> Parameter file name
      character(len=:), allocatable, intent(out) :: paramFile

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

      !> Input for DIPRO
      type(jab_input), intent(inout) :: dipro

      !> Stuff for second argument parser
!   integer  :: narg
!   character(len=p_str_length), dimension(p_arg_length) :: argv
!   type(TAtomList) :: atl
!   integer, allocatable :: list(:)

      !> Input for TBLite calculator
      type(TTBLiteInput), intent(out) :: tblite

!$    integer :: omp_get_num_threads, nproc
      integer :: nFlags
      integer :: idum, ndum
      real(wp) :: ddum
      character(len=:), allocatable :: flag, sec
      logical :: exist

      set%gfn_method = 2
      dipro%diprocalc = .false.
      coffee = .false.
      strict = .false.
      restart = .true.
      copycontrol = .false.
      lgrad = .false.
      gsolvstate = solutionState%gsolv
      tblite%color = get_xtb_feature('color')

      nFlags = args%countFlags()
      call args%nextFlag(flag)
      do while (allocated(flag))
         if (len(flag) > 2 .and. flag(1:1) == '-' .and. flag(1:2) /= '--') then
            call env%warning("the use of '"//flag//"' is discouraged, "// &
               & "please use '-"//flag//"' next time", source)
            flag = '-'//flag
         end if
         select case (flag)
         case default
            call env%warning("Unknown option '"//flag//"' provided", source)

         case ('-h', '--help')
            call help(env%unit)
            call terminate(0)

         case ('--citation')
            call citation(env%unit)
            call terminate(0)

         case ('--license')
            call disclamer(env%unit)
            call terminate(0)

         case ('--version')
            call xtb_header(env%unit)
            call terminate(0)

         case ('-v', '--verbose')
            set%verbose = .true.

         case ('-V', '--very-verbose')
            set%verbose = .true.
            set%veryverbose = .true.

         case ('--define')
            call set_define

         case ('-P', '--parallel')
!$          if (.false.) then
               call env%warning('Program compiled without threading support', source)
!$          end if
            ! Always remove next argument to keep argument parsing consistent
            call args%nextArg(sec)
!$          if (allocated(sec)) then
!$             if (getValue(env, sec, idum)) then
!$                nproc = omp_get_num_threads()
!$                call omp_set_num_threads(idum)
#ifdef WITH_MKL
!$                call mkl_set_num_threads(idum)
#endif
!$             end if
!$          end if

         case ('--restart')
            restart = .true.

         case ('--norestart')
            restart = .false.

         case ('--copy')
            copycontrol = .true.

         case ('--nocopy')
            copycontrol = .false.

         case ('--strict')
            strict = .true.

         case ('-I', '--input')
            call args%nextArg(inputFile)
            if (.not. allocated(inputFile)) then
               call env%error("Filename for detailed input is missing", source)
            end if

         case ('--namespace')
            call args%nextArg(persistentEnv%io%namespace)
            if (.not. allocated(persistentEnv%io%namespace)) then
               call env%error("Namespace argument is missing", source)
            end if

         case ('--vparam')
            call args%nextArg(paramFile)
            if (.not. allocated(paramFile)) then
               call env%error("Filename for --vparam is missing", source)
            else
               tblite%param = paramFile
            end if

         case ('--coffee')
            coffee = .true.

         case ('-a', '--acc')
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (getValue(env, sec, ddum)) then
                  if (ddum < 1.e-4_wp) then
                     call env%warning("We cannot provide this level of accuracy, "//&
                        & "resetted accuracy to 0.0001", source)
                     set%acc = 1.e-4_wp
                  else if (ddum > 1.e+3_wp) then
                     call env%warning("We cannot provide this level of accuracy, "//&
                        & "resetted accuracy to 1000", source)
                     set%acc = 1.e+3_wp
                  else
                     set%acc = ddum
                  end if
               end if
               tblite%accuracy = set%acc
            else
               call env%error("Accuracy is not provided", source)
            end if

         case ('-c', '--chrg', '--charge')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_chrg(env, sec)
            else
               call env%error("Molecular charge is not provided", source)
            end if
         
         case('--ceh')
            if (get_xtb_feature('tblite')) then
               allocate(tblite%ceh)
               call set_runtyp('prescc')
               ! check if numerical gradients should be calculated !
               call args%nextArg(sec)
               if (allocated(sec)) then
                  if (sec == 'grad') then
                     tblite%ceh%grad = .true.
                     ! check if step size is provided !
                     call args%nextArg(sec)
                     if (allocated(sec)) then
                        if (getValue(env, sec, ddum)) then
                           tblite%ceh%step = ddum
                        end if
                     end if

                  else
                     call env%warning("Unknown CEH option '"//sec//"' provided", source)
                  end if
               end if 
            else
               call env%error("CEH charges are only available through tblite library", source)
               return
            end if


         case ('-u', '--uhf')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_spin(env, sec)
            else
               call env%error("Number of unpaired electrons is not provided", source)
            end if

         case ("--efield")
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_efield(env, sec)
            else
               call env%error("Electric field is not provided", source)
            end if

         case ('--gfn')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_gfn(env, 'method', sec)
               if (sec == '0') call set_exttyp('eht')
               tblite%method = "gfn"//sec
            else
               call env%error("No method provided for --gfn option", source)
            end if

         case ('--gfn1')
            call set_gfn(env, 'method', '1')
            call env%warning("The use of '"//flag//"' is discouraged, "//&
               & "please use '--gfn 1' next time", source)
            tblite%method = "gfn1"

         case ('--gfn2')
            call set_gfn(env, 'method', '2')
            call set_gfn(env, 'd4', 'true')
            tblite%method = "gfn2"

         case ('--gfn0')
            call set_gfn(env, 'method', '0')
            call set_exttyp('eht')
            call env%warning("The use of '"//flag//"' is discouraged, "//&
               & "please use '--gfn 0' next time", source)

         case ('--gfnff')
            call set_exttyp('ff')

         case ('--gff')
            call set_exttyp('ff')

         case ('--mcgfnff')
            call set_exttyp('mcff')

         case ('--iff')
            call set_exttyp('iff')

         case ('--ptb')
            call set_exttyp('ptb')
            if (.not. get_xtb_feature('tblite')) then
               call ptb_feature_not_implemented(env)
               return
            endif

         case ('--tblite')
            if (get_xtb_feature('tblite')) then
               call set_exttyp('tblite')
            else
               call env%error("Compiled without support for tblite library", source)
               return
            end if

         case ('--color')
            if (allocated(sec)) then
               select case (sec)
               case ('auto')
                  tblite%color = get_xtb_feature('color')
               case ('always')
                  tblite%color = .true.
               case ('never')
                  tblite%color = .false.
               case default
                  call env%warning("Unknown color option '"//sec//"' provided", source)
               end select
            else
               call env%error("No color scheme provided for --color option", source)
            end if

         case ('--spinpol')
            if (get_xtb_feature('tblite')) then
               tblite%spin_polarized = .true.
            else
           call env%error("Compiled without support for tblite library. This is required for spin-polarization", source)
               return
            end if

         case ('--dipro')
            if (get_xtb_feature('tblite')) then
               dipro%diprocalc = .true.
               call set_runtyp('scc')
               call args%nextArg(sec)
               if (allocated(sec)) then
                  read (sec, '(f10.3)') dipro%othr
               else
                  dipro%othr = 0.1_wp
               end if
            else
               call env%error("Compiled without support for tblite library. This is required for DIPRO", source)
               return
            end if

         case ('--oniom')
            call set_exttyp('oniom')
            call args%nextArg(sec)

            if (.not. allocated(sec)) then ! handle no argument case !
               call env%error("No inner region is  provided for ONIOM", source)
               return
            end if
            call move_alloc(sec, oniom%first_arg)

            call args%nextArg(sec)
            if (.not. allocated(sec)) then
               call env%warning("No method is specified for ONIOM," &
                     &//achar(10)//" default gfn2:gfnff combination will be used", source)
               call move_alloc(oniom%first_arg, sec)
            end if

            inquire (file=sec, exist=exist)
            if (exist) then
               sec = read_whole_file(sec)
            end if
            call move_alloc(sec, oniom%second_arg)

         case ('--cut')
            call set_cut

         case ('--etemp')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_scc(env, 'temp', sec)
            else
               call env%error("Temperature in --etemp option is missing", source)
            end if

         case ('--esp')
            call set_runtyp('scc')
            call set_write(env, 'esp', 'true')

         case ('--stm')
            call set_runtyp('scc')
            call set_write(env, 'stm', 'true')

         case ('--cma')
            call set_cma

         case ('--tm')
            call set_exttyp('turbomole')

         case ('--enso')
            call set_enso_mode

         case ('--json')
            call set_write(env, 'json', 'true')
            Call setWRtopo("json", printTopo)

         case ('--ceasefiles')
            restart = .false.
            set%verbose = .false.
            set%ceasefiles = .true.
            call set_write(env, 'wiberg', 'false')
            call set_write(env, 'charges', 'false')
#ifdef _WIN32
            call set_opt(env, 'logfile', 'NUL')
#else
            call set_opt(env, 'logfile', '/dev/null')
#endif

         case ('--orca')
            call set_exttyp('orca')

         case ('--driver')
            call set_exttyp('driver')
            call args%nextArg(sec)
            if (allocated(sec)) then
               set%ext_driver%executable = sec
            end if

         case ('--mopac')
            call set_exttyp('mopac')

         case ('--pop')
            call set_write(env, 'mulliken', 'true')

         case ('--molden')
            call set_write(env, 'mos', 'true')

         case ('--dipole')
            call set_write(env, 'dipole', 'true')

         case ('--wbo')
            call set_write(env, 'wiberg', 'true')

         case ('--lmo')
            call set_write(env, 'mulliken', 'true')
            call set_write(env, 'lmo', 'true')

         case ('--ewin')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_siman(env, 'ewin', sec)
            else
               call env%error("Real argument for --ewin is missing", source)
            end if

         case ('--fod')
            call set_write(env, 'fod', 'true')
            call set_scc(env, 'temp', '5000.0')

         case ('--iterations', '--maxiterations')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_scc(env, 'maxiterations', sec)
            else
               call env%error("Integer argument for --iterations is missing", source)
            end if

         case ('--cycles')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_opt(env, 'maxcycle', sec)
            else
               call env%error("Integer argument for --cycles is missing", source)
            end if

         case ('-g', '--gbsa')
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

         case ('--alpb')
            call args%nextArg(sec)
            call set_gbsa(env, 'alpb', 'true')
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
               call env%error("No solvent name provided for ALPB", source)
            end if

         case ('--cosmo', '--tmcosmo')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_gbsa(env, 'solvent', sec)
               call set_gbsa(env, flag(3:), 'true')
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

         case ('--cpcmx')
            if (get_xtb_feature('cpcmx')) then
               call args%nextArg(sec)
               if (allocated(sec)) then
                  call set_gbsa(env, 'solvent', 'infinity')
                  call set_gbsa(env, 'cosmo', 'true')
                  call set_gbsa(env, 'cpcmx', sec)
               else
                  call env%error("No solvent name provided for CPCM-X", source)
               end if
            else
               call env%error("The CPCM-X library was not included in this version of xTB.", source)
            end if

         case ('--scc', '--sp')
            call set_runtyp('scc')

         case ('--vip')
            call set_gfn(env, 'method', '1')
            call set_runtyp('vip')

         case ('--vea')
            call set_gfn(env, 'method', '1')
            call set_runtyp('vea')

         case ('--vipea')
            call set_gfn(env, 'method', '1')
            call set_runtyp('vipea')

         case ('--vomega')
            call set_gfn(env, 'method', '1')
            call set_runtyp('vomega')

         case ('--vfukui')
            call set_runtyp('vfukui')

         case ('--alpha')
            call set_elprop('alpha')

         case ('--raman')
            call set_elprop('raman')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_raman(env,sec)
               call args%nextArg(sec)
               if (allocated(sec)) then
                  call set_raman(env,sec)
               endif
            endif


         case ('--grad')
            call set_runtyp('grad')
            lgrad = .true.

         case ('-o', '--opt')
            call set_runtyp('opt')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_opt(env, 'optlevel', sec)
            end if

         case ('--hess')
            call set_runtyp('hess')

         case ('--md')
            call set_runtyp('md')

         case ('--ohess')
            call set_runtyp('ohess')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_opt(env, 'optlevel', sec)
            end if

         case ('--bhess')
            call set_runtyp('bhess')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_opt(env, 'optlevel', sec)
            end if

         case ('--omd')
            call set_runtyp('omd')
            call set_opt(env, 'optlevel', '-1')

         case ('--siman')
            call set_runtyp('siman')
            call set_md(env, 'nvt', 'true')

         case ('--path')
            call set_runtyp('path')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_path(env, 'product', sec)
            end if

         case ('--screen')
            call set_runtyp('screen')

         case ('--gmd')
            call set_runtyp('gmd')
            call env%error("This feature has been deprecated, I'm sorry.", source)

         case ('--modef')
            call set_runtyp('modef')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_modef(env, 'mode', sec)
            end if

         case ('--mdopt')
            call set_runtyp('mdopt')

         case ('--metadyn')
            call set_runtyp('md')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_metadyn(env, 'save', sec)
            end if
            call set_metadyn(env, 'static', 'false')

         case ('--metaopt')
            call set_runtyp('metaopt')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_opt(env, 'optlevel', sec)
            end if

         case ('--nat')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_natom(env, sec)
            end if

         case ('--bias-input', '--gesc')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_metadyn(env, 'bias-input', sec)
            else
               call env%error("No input file for RMSD bias provided", source)
            end if

         case ('--wrtopo')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call setWRtopo(sec, printTopo)
               if (printTopo%warning) call env%error("A wrtopo argument has been misspelled.", source)
            else
               call env%error("The wrtopo keyword is missing an argument.", source)
            end if

         case ('--nocellopt')
            set%optcell = .false.

         end select
         call args%nextFlag(flag)
      end do

   end subroutine parseArguments

!> kronecker delta
   function kron(i, j) result(res_kronij)
      integer, intent(in) :: i, j
      real(wp) :: res_kronij

      res_kronij = 0.0_wp
      if (i == j) res_kronij = 1.0_wp

   end function kron
   function read_whole_file(fname) result(list)
      character(len=*), intent(in) :: fname
      character(len=:), allocatable :: list
      integer :: io, stat
      character(len=:), allocatable :: line
      open (newunit=io, file=fname, iostat=stat)
      call getline(io, list, stat)
      do while (stat == 0)
         call getline(io, line, stat)
         if (stat == 0) list = list//","//line
      end do
      close (io, iostat=stat)
   end function read_whole_file

! set booleans for requested topology list printout
   subroutine setWRtopo(sec, printTopo)
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
      do i = 1, lenSec
         curr_pos = scan(sec(curr_pos + 1:lenSec), sep) + old_pos
         if (curr_pos /= old_pos) then
            call selectList(sec(old_pos + 1:curr_pos - 1), printTopo)
         else
            call selectList(sec(old_pos + 1:lenSec), printTopo)
            exit
         end if
         old_pos = curr_pos
      end do

   end subroutine setWRtopo

   subroutine selectList(secSplit, printTopo)
      ! part of command line argument
      character(len=*), intent(in) :: secSplit
      ! holds booleans of to be printed topology lists
      type(TPrintTopo), intent(inout) :: printTopo

      select case (secSplit)
      case ("etot")
         printTopo%etot = .true.
      case ("gnorm")
         printTopo%gnorm = .true.
      case ("nb")
         printTopo%nb = .true.
      case ("bpair")
         printTopo%bpair = .true.
      case ("alist")
         printTopo%alist = .true.
      case ("blist")
         printTopo%blist = .true.
      case ("tlist")
         printTopo%tlist = .true.
      case ("vtors")
         printTopo%vtors = .true.
      case ("vbond")
         printTopo%vbond = .true.
      case ("vangl")
         printTopo%vangl = .true.
      case ("hbbond")
         printTopo%hbbond = .true.
      case ("eeq")
         printTopo%eeq = .true.
      case ("json")
         printTopo%etot = .true.
         printTopo%gnorm = .true.
         printTopo%nb = .true.
         printTopo%bpair = .true.
         printTopo%alist = .true.
         printTopo%blist = .true.
         printTopo%tlist = .true.
         printTopo%vtors = .true.
         printTopo%vbond = .true.
         printTopo%vangl = .true.
         printTopo%hbbond = .true.
         printTopo%eeq = .true.
      case default
         printTopo%warning = .true.
      end select
   end subroutine selectList

   subroutine ptb_feature_not_implemented(env)
      !> Computational environment
      type(TEnvironment), intent(inout) :: env

      call env%error("PTB not available without 'tblite'. Compiled without support for 'tblite' library.")
      call env%error("Please recompile without '-Dtblite=disabled' option or change meson setup.")
   end subroutine ptb_feature_not_implemented

   !>  make some post processing afterward, show some timings and stuff
   subroutine finalize_xtb(env)
    
      !> Calculation environment
      type(TEnvironment), intent(in) :: env

      write (env%unit, '(a)')
      write (env%unit, '(72("-"))')
      call stop_timing_run
      call stop_timing(1)
      call prdate('E')
      write (env%unit, '(72("-"))')
      call prtiming(1, 'total')
      if (.not. set%runtyp == p_run_prescc) &
         call prtiming(2, 'SCF')
      if ((set%runtyp == p_run_opt) .or. (set%runtyp == p_run_ohess) .or. &
         &   (set%runtyp == p_run_omd) .or. (set%runtyp == p_run_metaopt)) then
         call prtiming(3, 'ANC optimizer')
      end if
      if (set%runtyp == p_run_path) then
         call prtiming(4, 'path finder')
      end if
      if (((set%runtyp == p_run_hess) .or. (set%runtyp == p_run_ohess) .or. (set%runtyp == p_run_bhess))) then
         if (set%mode_extrun /= p_ext_turbomole) then
            call prtiming(5, 'analytical hessian')
         else
            call prtiming(5, 'numerical hessian')
         end if
      end if
      if ((set%runtyp == p_run_md) .or. (set%runtyp == p_run_omd) .or. &
          (set%runtyp == p_run_metaopt)) then
         call prtiming(6, 'MD')
      end if
      if (set%runtyp == p_run_screen) then
         call prtiming(8, 'screen')
      end if
      if (set%runtyp == p_run_modef) then
         call prtiming(9, 'mode following')
      end if
      if (set%runtyp == p_run_mdopt) then
         call prtiming(10, 'MD opt.')
      end if

      write (env%unit, '(a)')
      call terminate(0)

   end subroutine finalize_xtb

end module xtb_prog_main
