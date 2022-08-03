! This file is part of xtb.
!
! Copyright (C) 2022 Christoph Plett
!
! SPDX-Identifier: LGPL-3.0-or-later
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

! Molecule A and B contain always the CMA transformed coords.
! Every structure is saved as array(1:6) where (1:3) are the
! translation of B starting at the origin and (4:6) are the
! rotation angles of molB.

!> Docking implementation for xtb
module xtb_prog_dock
   use xtb_type_environment, only: TEnvironment, init
   use xtb_prog_argparser, only: TArgParser
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_filetypes, only: getFileType, fileType
   use xtb_mctc_timings
   use xtb_mctc_version, only: version, author, date
   use xtb_io_reader, only: readMolecule
   use xtb_type_environment, only: TEnvironment
   use xtb_type_molecule
   use xtb_type_reader, only: TReader
   use xtb_setparam, only: initrand
   use xtb_gfnff_param, only: gff_print
   use xtb_mctc_systools
   use xtb_setmod
   use xtb_docking_set_module
   use xtb_docking_param
   use xtb_mctc_convert, only: autokcal

   implicit none

contains

!> Entry point for performing non-covalent docking
   subroutine xtbDock(env, argParser)
      use xtb_iff_iffini, only: init_iff
      use xtb_iff_iffprepare, only: precomp
      use xtb_iff_iffenergy
      use xtb_docking_search_nci, only: docking_search
      use xtb_sphereparam, only: sphere, rabc, boxr, init_walls, wpot
      use xtb_constrain_param, only: read_userdata
      use xtb_fixparam, only: init_fix
      use xtb_scanparam, only: init_constr, init_scan
      use xtb_embedding, only: init_pcem
      use xtb_splitparam, only: init_split

      !> Source of errors in the main program unit
      character(len=*), parameter :: source = "prog_dock"

      !> Calculation environment
      type(TEnvironment), intent(inout) :: env

      !> Command line arguments
      type(TArgParser), intent(inout) :: argParser

      !> Molecular structure data
      type(TMolecule) :: molA, molB

      !> Combined structure
      type(TMolecule) :: comb

      !> Storage of translational difference of molB and molA (1:3) and rotation angles of molB (4:6)
      real(wp) :: icoord(6)

      type(TReader) :: reader

      !> File names
      character(len=:), allocatable :: fnameA, fnameB, fname, fnam

      !> File types
      integer :: nFiles, ftypeA, ftypeB, ich
      !character(len=:),allocatable :: xcontrol ! instruction file

      !> global instruction file
      character(len=:), allocatable :: xrc

      logical :: exist

      !> Data of Molecules
      integer :: n, n1, n2
      integer :: nlmo1, nlmo2
      real(wp), allocatable :: xyz1(:, :)
      real(wp), allocatable :: rlmo1(:, :)
      real(wp), allocatable :: q1(:)
      real(wp), allocatable :: qdr1(:)
      real(wp), allocatable ::xyzdr1(:, :)
      real(wp), allocatable :: cn1(:)
      real(wp), allocatable :: z1(:)
      real(wp), allocatable :: alp1(:)
      real(wp), allocatable :: qct1(:, :)
      integer, allocatable :: at1(:)
      integer, allocatable :: lmo1(:)
      real(wp), allocatable :: xyz2(:, :)
      real(wp), allocatable :: rlmo2(:, :)
      real(wp), allocatable :: q2(:)
      real(wp), allocatable :: qdr2(:)
      real(wp), allocatable ::xyzdr2(:, :)
      real(wp), allocatable :: cn2(:)
      real(wp), allocatable :: z2(:)
      real(wp), allocatable :: alp2(:)
      real(wp), allocatable :: qct2(:, :)
      integer, allocatable :: at2(:)
      integer, allocatable :: lmo2(:)
      real(wp), allocatable :: c6ab(:, :)
      real(wp), allocatable :: alpab(:, :)
      real(wp), allocatable :: cprob(:)
      real(wp), allocatable :: xyz(:, :)
      real(wp), allocatable :: q(:)
      real(wp), allocatable :: cn(:)
      real(wp), allocatable :: alp0(:)
      real(wp), allocatable :: gab1(:, :)
      real(wp), allocatable :: gab2(:, :)
      real(wp), allocatable :: den1(:, :, :)
      real(wp), allocatable :: den2(:, :, :)
      real(wp), allocatable :: qcm1(:)
      real(wp), allocatable :: qcm2(:)
      integer, allocatable :: at(:)
      integer, allocatable :: neigh(:, :)
      real(wp) :: e, dum, xx(10), zeffqmdff(86) ! various
      real(wp) :: icoord0(6) ! internal coords
      real(wp) :: r(3)
      integer :: i, j, k
      real(wp) :: molA_e, molB_e

      !> Print an informative banner
      call dockingHeader(env%unit)

      !> Parse arguments
      call parseArguments(env, argParser, fname)
      if (optlvl == 'gfn0') call env%error('GFN0 is not supported for docking')
      if (optlvl == 'gfnff') set%mhset%model = p_modh_gff ! With p_modh_old some problems occur
      nFiles = argParser%countFiles()
      if (nFiles == 0) then
         call env%error("No input file given, so there is nothing to do", source)
         call dockingHelp(env%unit)
      elseif (nFiles == 1) then
         call env%error("Please provide two sets of coordinates, aborting", source)
      else
         call argParser%nextFile(fnameA)
         call argParser%nextFile(fnameB)
         if (allocated(fname)) xcontrol = fname
      end if
      call env%checkpoint("Command line argument parsing failed")

      ! ------------------------------------------------------------------------
      !> read the xtbrc if you can find it (use rdpath directly instead of xfind)
      if (allocated(xcontrol)) call rdcontrol_iff(xcontrol, env, .false.)

      call rdpath(env%xtbpath, '.xtbrc', xrc, exist)
      if (exist) then
         call rdcontrol_iff(xrc, env, .false.)
         call env%checkpoint("Reading '"//xrc//"' failed")
      end if

      call env%checkpoint("Command line argument parsing failed")
      !> Determine the file type for processing in the reader
      ftypeA = getFileType(fnameA)
      !> Generate the reader by connecting the file to the instance
      call reader%open(fnameA)
      !> Process the file to obtain the molecular structure data
      call readMolecule(env, molA, reader%unit, ftypeA)
      !> Close the file, we are done with it here
      call reader%close
      call env%checkpoint("Could not read geometry from '"//fnameA//"'")
      !> Determine the file type for processing in the reader
      ftypeB = getFileType(fnameB)
      !> Generate the reader by connecting the file to the instance
      call reader%open(fnameB)
      !> Process the file to obtain the molecular structure data
      call readMolecule(env, molB, reader%unit, ftypeB)
      !> Close the file, we are done with it here
      call reader%close
      call env%checkpoint("Could not read geometry from '"//fnameB//"'")

      !> Print current time
      call prdate('S')

      !> Random number initiation
      call initrand

      !> Printout Settings
      call dockingPrintout(env%unit, fnameA, fnameB, molA, molB)

      !> Set some parameter
      call set_iff_param
      fnam = 'xtblmoinfo'

      !> Check .CHRG, .UHF and xcontrol
      call check_for_files(env, molA, molB)

      !> Get IFF required properties with GFN2 singlepoints
      set%pr_local = .false.
      call start_timing(2)
      write(*,*)
      write (env%unit, *) 'Precomputation of electronic porperties'
      write (env%unit, *) ' For Molecule 1'
      call precomp(env, molA, molA_e, 1)
      write (env%unit, *) ' Successful'
      call stop_timing(2)
      !MolA
      call rd0(1, trim(fnam), n1, nlmo1)
      allocate (at1(n1), lmo1(10*n1), source=0)
      allocate (xyz1(3, n1), rlmo1(4, 10*n1), q1(n1),&
      &cn1(n1), alp1(n1), qct1(n1, 2), qdr1(n1), xyzdr1(3, n1),&
      &z1(n1), den1(2, 4, n1), gab1(n1, n1), qcm1(n1), cprob(n1), source=0.0_wp)
      call rd(trim(fnam), 1, n1, xyz1, at1, nlmo1, lmo1, rlmo1, q1, qct1)
      call delete_file(trim(fnam))

      !MolB
      write (env%unit, *) ' For Molecule 2'
      call precomp(env, molB, molB_e, 2)
      write (env%unit, *) ' Successful'
      call rd0(2, trim(fnam), n2, nlmo2)
      allocate (at2(n2), lmo2(10*n2), source=0)
      allocate (xyz2(3, n2), rlmo2(4, 10*n2), q2(n2),&
      & cn2(n2), alp2(n2), qct2(n2, 2), qdr2(n2), xyzdr2(3, n2),&
      & z2(n2), den2(2, 4, n2), gab2(n2, n2), qcm2(n2), source=0.0_wp)
      call rd(trim(fnam), 2, n2, xyz2, at2, nlmo2, lmo2, rlmo2, q2, qct2)
      call delete_file(trim(fnam))

      !> Special Docking CMA shift
      call cmadock(molA%n, molA%n, molA%at, molA%xyz, r)
      do i = 1, 3
         molA%xyz(i, 1:molA%n) = molA%xyz(i, 1:molA%n) - r(i)
      end do
      call cmadock(molB%n, molB%n, molB%at, molB%xyz, r)
      do i = 1, 3
         molB%xyz(i, 1:molB%n) = molB%xyz(i, 1:molB%n) - r(i)
      end do

      !>  Combined Structure
      n = n1 + n2
      allocate (at(n), neigh(0:n, n), source=0)
      allocate (xyz(3, n), q(n), c6ab(n, n), alp0(n), cn(n),&
      &        alpab(n2, n1), source=0.0_wp)

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Initialize the IFF energy stuff
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call init_iff(env, n1, n2, at1, at2, neigh, xyz1, xyz2, q1, q2, c6ab, z1, z2,&
         &          cprob, nlmo1, nlmo2, lmo1, lmo2,&
         &          qdr1, qdr2, rlmo1, rlmo2, cn1, cn2, alp1, alp2, alpab,&
         &          den1, den2, gab1, gab2, qcm1, qcm2, n, at, xyz, q, icoord, icoord0,&
         &          .false.)


      !> CONSTRAINTS & SCANS
      call init_fix(n)
      call init_split(n)
      call init_constr(n, at)
      call init_scan
      call init_walls
      call init_pcem
      !> Read the constrain
      call init(comb,at,xyz)
      !comb%n = n
      !comb%at = at
      if (allocated(xcontrol)) then
         call read_userdata(xcontrol, env, comb)
         call read_userdata_iff(xcontrol, env, comb)
      end if

      !> For directed docking, repulsive potentials for each atom other than
      !  the defined ones is setup as range dependent on the smallest distance
      !  to these atoms
      if(directedset%n > 0) then
         if(directedset%n > molA%n) call env%error(&
           &"More atoms for directed docking defined than in molecule A", source)
         if (directed_type == p_atom_pot) call get_repulsive_pot(env,xyz,comb)
         if (directed_type == p_atom_att) call get_attractive_pot(env,comb)
      end if

      deallocate (comb%at, comb%xyz)

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Single Point energy of first 'Cold fusion' structure
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call start_timing(3)
      call iff_e(env, n, n1, n2, at1, at2, neigh, xyz1, xyz2, q1, q2, c6ab, z1, z2,&
                       & nlmo1, nlmo2, lmo1, lmo2, rlmo1, rlmo2,&
                       & qdr1, qdr2, cn1, cn2, alp1, alp2, alpab, qct1, qct2,&
                       & den1, den2, gab1, gab2,&
                       & set%verbose, 0, e, icoord)
      write(*,*)
      write(*, '(''Energy of cold fusion in kcal/mol:'', F8.2)') e * autokcal
      call stop_timing(3)
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Real docking search
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call set_optlvl(env) !Sets the optimization level for search in global parameters
      call docking_search(env, molA, molB, n, n1, n2, at1, at2, neigh, xyz1,&
                    & xyz2, q1, q2, c6ab, z1, z2,&
                    &nlmo1, nlmo2, lmo1, lmo2, rlmo1, rlmo2,&
                    &qdr1, qdr2, cn1, cn2, alp1, alp2, alpab, qct1, qct2,&
                    &den1, den2, gab1, gab2, molA_e, molB_e,&
                    &cprob, e, icoord, comb)
      if(debug) call wrc('screening.coord', n1, n2, at1, at2, xyz1, xyz2, icoord)

      deallocate (at1, xyz1, rlmo1, q1, lmo1, cn1, alp1, qct1, qdr1, xyzdr1,&
     &           at2, xyz2, rlmo2, q2, lmo2, cn2, alp2, qct2, qdr2, xyzdr2,&
     &           den1, den2, gab1, gab2, at, xyz, q, c6ab, alp0, cn, neigh,&
     &           z1, z2, alpab, cprob)

      call env%checkpoint("Docking submodule")

      ! ------------------------------------------------------------------------
      !  we may have generated some non-fatal errors, which have been saved,
      !  so we should tell the user, (s)he may want to know what went wrong
      call env%show("Runtime exception occurred")
      call raise('F', 'Some non-fatal runtime exceptions were caught,'// &
         &           ' please check:')

      !  make some post processing afterward, show some timings and stuff
      write (env%unit, '(a)')
      write (env%unit, '(72("-"))')
      call stop_timing_run
      call stop_timing(1)
      call prdate('E')
      write (env%unit, '(72("-"))')
      call prtiming(1, 'total')
      call prtiming(2, 'LMO Computation')
      call prtiming(3, 'SP')
      call prtiming(4, 'Searching')
      call prtiming(5, 'ANC optimizer')

      write (env%unit, '(a)')
      call terminate(0)

   end subroutine xtbDock

   subroutine parseArguments(env, args, inputFile)
      use xtb_readin, only: getValue
      use xtb_solv_state

      !> Name of error producer
      character(len=*), parameter :: source = "prog_docking_parseArguments"

      !> Calculation environment
      type(TEnvironment), intent(inout) :: env

      !> Command line argument parser
      type(TArgParser), intent(inout) :: args

      !> Detailed input file name
      character(len=:), allocatable, intent(out) :: inputFile

      integer :: nFlags, idum
      character(len=:), allocatable :: flag, sec, trd
      real(wp) :: ddum

      !> Defaults
      optlvl = 'gfn2'

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

         case ('--help', '-h')
            call dockingHelp(env%unit)
            call terminate(0)

         case ('--dummy')
            write (*, *) 'Dummy option used'

         case ('--verbose')
            set%verbose = .true.
            debug = .true.

         case ('--very%verbose')
            set%veryverbose = .true.

         case ('--cs')
            cssym = .true.

         case ('--onlypocket')
            stack_grid = .false.
            angular_grid = .false.
            pocket_grid = .true.

         case ('-nopocket')
            call set_logicals(env, 'nopocket')

         case ('--pocket')
            call set_logicals(env, 'pocket')

         case ('--onlystack')
            pocket_grid = .false.
            stack_grid = .true.
            angular_grid = .false.

         case ('--nostack')
            call set_logicals(env, 'nostack')

         case ('--stack')
            call set_logicals(env, 'stack')

         case ('--onlyangular')
            pocket_grid = .false.
            stack_grid = .false.
            angular_grid = .true.

         case ('--noangular')
            call set_logicals(env, 'noangular')

         case ('--angular')
            call set_logicals(env, 'angular')

         case ('--org')
            call set_logicals(env, 'org')

         case ('--qcg ', '--fast')  ! fast version
            call set_logicals(env, 'fast')

         case ('--noind')
            call set_logicals(env, 'noind')

         case ('--loose')
            call set_logicals(env, 'loose')

         case ('--atm')
            call set_logicals(env, 'atm')

         case ('--stepr')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_docking(env, 'stepr', sec)
            else
               call env%error("Radial step --stepr option is missing", source)
            end if

         case ('--stepa')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_docking(env, 'stepa', sec)
            else
               call env%error("Angular step --stepa option is missing", source)
            end if

         case ('--nfinal')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_docking(env, 'nfinal', sec)
            else
               call env%error("Maxopt in --nfinal option is missing", source)
            end if

         case ('--maxgen')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_docking(env, 'maxgen', sec)
            else
               call env%error("Maxgen in --maxgen option is missing", source)
            end if

         case ('--maxparent')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_docking(env, 'maxparent', sec)
            else
            call env%error("Maxparent in --maxparent option is missing", source)
            end if

         case ('--etemp')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_scc(env, 'temp', sec)
            else
              call env%error("Temperature in --etemp option is missing", source)
            end if

         case ('--iterations')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_scc(env, 'maxiterations', sec)
            else
          call env%error("Integer argument for --iterations is missing", source)
            end if

         case ('-o', '--opt')
            call set_runtyp('opt')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_opt(env, 'optlevel', sec)
            end if

         case ('--cycles')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_opt(env, 'maxcycle', sec)
            else
              call env%error("Integer argument for --cycles is missing", source)
            end if

         case ('--optlvl')
            call args%nextArg(sec)
            select case (sec)
            case default
               call env%warning("No optlevel specified, using GFN-FF", source)
            case ('gfn')
               call args%nextArg(trd)
               if (allocated(trd)) then
                  if (trd == '0') then; optlvl = 'gfn0'; call set_gfn(env, 'method', '0'); end if
                  if (trd == '1') then; optlvl = 'gfn1'; call set_gfn(env, 'method', '1'); end if
                  if (trd == '2') then; optlvl = 'gfn2'; call set_gfn(env, 'd4', 'true'); end if
                  if (trd == 'ff') optlvl = 'gfnff'
               else
                  call env%error("No method provided for --gfn option", source)
               end if

            case ('gfn1')
               optlvl = 'gfn1'
               call set_gfn(env, 'method', '1')

            case ('gfn2')
               optlvl = 'gfn2'
               call set_gfn(env, 'd4', 'true')

            case ('gfn0')
               optlvl = 'gfn0'

            case ('gfnff')
               optlvl = 'gfnff'

            case ('gff')
               optlvl = 'gfnff'

            end select

         case ('--gfn')
            call args%nextArg(trd)
            if (allocated(trd)) then
               if (trd == '0') then; optlvl = 'gfn0'; call set_gfn(env, 'method', '0'); end if
               if (trd == '1') then; optlvl = 'gfn1'; call set_gfn(env, 'method', '1'); end if
               if (trd == '2') then; optlvl = 'gfn2'; call set_gfn(env, 'd4', 'true'); end if
               if (trd == 'ff') optlvl = 'gfnff'
            else
               call env%error("No method provided for --gfn option", source)
            end if

         case ('--gfn1')
            optlvl = 'gfn1'
            call set_gfn(env, 'method', '1')
            call env%warning("The use of '"//flag//"' is discouraged, "//&
               & "please use '--gfn 1' next time", source)

         case ('--gfn2')
            optlvl = 'gfn2'
            call set_gfn(env, 'd4', 'true')

         case ('--gfn0')
            optlvl = 'gfn0'
            call set_gfn(env, 'method', '0')
            call env%warning("The use of '"//flag//"' is discouraged, "//&
               & "please use '--gfn 0' next time", source)

         case ('--gfnff')
            optlvl = 'gfnff'

         case ('--gff')
            optlvl = 'gfnff'

         case ('-a', '--acc')
            call args%nextArg(sec)
            if (allocated(sec)) then
               if (getValue(env, sec, ddum)) then
                  if (ddum .lt. 1.e-4_wp) then
                call env%warning("We cannot provide this level of accuracy, "//&
                             & "resetted accuracy to 0.0001", source)
                     acc = 1.e-4_wp
                  else if (ddum .gt. 1.e+3_wp) then
                call env%warning("We cannot provide this level of accuracy, "//&
                             & "resetted accuracy to 1000", source)
                     acc = 1.e+3_wp
                  else
                     acc = ddum
                  end if
               end if
            else
               call env%error("Accuracy is not provided", source)
            end if

         case ('--nfrag1')
            call args%nextArg(sec)
            if (getValue(env, sec, ddum)) nfrag1 = int(ddum)

         case ('--chrg', '--chrg1')
            call args%nextArg(sec)
            if (getValue(env, sec, ddum)) chrg(1) = ddum

         case ('--chrg2')
            call args%nextArg(sec)
            if (getValue(env, sec, ddum)) chrg(2) = ddum

         case ('--uhf', '--uhf1')
            call args%nextArg(sec)
            if (getValue(env, sec, ddum)) uhf(1) = ddum

         case ('--uhf2')
            call args%nextArg(sec)
            if (getValue(env, sec, ddum)) uhf(2) = ddum

         case ('-I', '--input')
            call args%nextArg(inputFile)
            if (.not. allocated(inputFile)) then
               call env%error("Filename for detailed input is missing", source)
            end if

!Implicit solvation is missing, but shows bugs
         case ('-g', '--gbsa')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_gbsa(env, 'solvent', sec)
               call set_gbsa(env, 'alpb', 'false')
               call set_gbsa(env, 'kernel', 'still')
               call args%nextArg(sec)
               if (allocated(sec)) then
                  if (sec == 'reference') then
                     gsolvstate_iff = solutionState%reference
                  else if (sec == 'bar1M') then
                     gsolvstate_iff = solutionState%mol1bar
                  else
                 call env%warning("Unknown reference state '"//sec//"'", source)
                  end if
               end if
            else
               call env%error("No solvent name provided for GBSA", source)
            end if

         case ('--alpb')
            call args%nextArg(sec)
            if (allocated(sec)) then
               call set_gbsa(env, 'solvent', sec)
               call args%nextArg(sec)
               if (allocated(sec)) then
                  if (sec == 'reference') then
                     gsolvstate_iff = solutionState%reference
                  else if (sec == 'bar1M') then
                     gsolvstate_iff = solutionState%mol1bar
                  else
                 call env%warning("Unknown reference state '"//sec//"'", source)
                  end if
               end if
            else
               call env%error("No solvent name provided for GBSA", source)
            end if

         case ('--ensemble', '--ens')
            call set_logicals(env, 'ensemble')

         end select
         call args%nextFlag(flag)
      end do

   end subroutine parseArguments

!> Header for this submodule
   subroutine dockingHeader(unit)

      !> IO unit
      integer, intent(in) :: unit

      write (unit, '(a)') &
         "      -----------------------------------------------------------", &
         "     |                   =====================                   |", &
         "     |                          g I S S                          |", &
         "     |                   =====================                   |", &
         "     |               S. Ehlert, S. Grimme, C.Plett               |", &
         "     |          Mulliken Center for Theoretical Chemistry        |", &
         "     |                    University of Bonn                     |", &
         "      -----------------------------------------------------------", ""
      write (unit, '(3x,"*",*(1x,a))') &
         & "xtb-docking algorithm", version, "compiled by", author, "on", date
      write (unit, '(a)')

   end subroutine dockingHeader

   subroutine dockingPrintout(iunit, fnameA, fnameB, molA, molB)
      use xtb_mctc_accuracy, only: wp
      use xtb_mctc_systools, only: rdvar
      use xtb_mctc_global, only : persistentEnv
      !$ use omp_lib

      implicit none
      integer, intent(in) :: iunit
      type(TMolecule), intent(in) :: molA, molB
      character(len=:), allocatable, intent(in) :: fnameA, fnameB

      character(len=:), allocatable :: cdum
      integer :: l, err
      real(wp) :: dum5

      write (iunit, '(a)')
      call generic_header(iunit, 'Calculation Setup', 49, 10)
      write (iunit, '(a)')
      if (allocated(cdum)) deallocate (cdum)
      call get_command(length=l)
      allocate (character(len=l) :: cdum)
      call get_command(cdum)
      write (iunit, '(10x,a,":",1x,a)') 'program call               ', cdum
      call rdvar('HOSTNAME', cdum, err)
      if (err .eq. 0) &
         write (iunit, '(10x,a,":",1x,a)') 'hostname                   ', cdum
      ! Number of threads
      !$omp parallel
      !$omp master
      !$ write(iunit,'(10x,a,":",6x,i16)') 'omp threads                ',omp_get_num_threads()
      !$omp end master
      !$omp end parallel
      ! ----------------------------------------------------------------------
      !  print the home and path to check if there are set correctly
      write (iunit, '(10x,a,":",1x,a)') 'coordinate file A          ', fnameA
      write (iunit, '(10x,a,":",1x,a)') 'coordinate file B          ', fnameB
      !  print more specific calculation data
      write (iunit, '(10x,a,":",6x,i16)') 'number of atoms A          ', molA%n
      write (iunit, '(10x,a,":",6x,i16)') 'number of atoms B          ', molB%n
      write (iunit, '(10x,a,":",6x,i16)') 'charge of molecule A       ', molA%chrg
      write (iunit, '(10x,a,":",6x,i16)') 'charge of molecule B       ', molB%chrg
      write (iunit, '(10x,a,":",6x,i16)') 'spin of molecule A         ', molA%uhf
      write (iunit, '(10x,a,":",6x,i16)') 'spin of molecule B         ', molB%uhf
      call random_number(dum5)
      write (iunit, '(10x,a,":",6x,f16.14)') 'first test random number   ', dum5
      if (set%veryverbose) then
         write (iunit, '(10x,a,":",20x,a)') 'is this your card?         ', "🃓"
      end if
      write (iunit, '(a)')

   end subroutine dockingPrintout

   subroutine dockingHelp(iunit)

      implicit none
      !> IO unit
      integer, intent(in) :: iunit

      write (iunit, '(a)') &
         "Usage: xtb dock <geometry> <geometry> [-options]", &
         "", &
         "<geometry> may be provided as any valid input to xtb", &
         "", &
         "Options:", &
         "", &
         "    --I,--input <file>: Provide an control file", &
         "    --nopocket        : Perform no pocket search [Default]", &
         "    --pocket          : Perform pocket search", &
         "    --onlypocket      : Perform only pocket search", &
         "    --nostack         : Perform no stack search [Default]", &
         "    --stack           : Perform stack search", &
         "    --onlystack       : Perform only stack search", &
         "    --noangular       : Perform no angular search [Default]", &
         "    --angular         : Perform angular search", &
         "    --onlyangular     : Perform only angular search", &
         "    --atm             : Include ATM term in iff Energy", &
      "    --fast            : Faster search, sufficient for small molecules", &
         "    --chrg <real>     : Charge of first molecule", &
         "    --chrg2 <real>    : Charge of second molecule", &
     "    --uhf <real>      : Number of unpaired electrons of first molecule", &
    "    --uhf2 <real>     : Number of unpaired electrons of second molecule", &
         "    --stepr <real>    : Grid step size for rare gas prescreening", &
         "    --stepa <real>    : Grid step size for angular grid in [deg]", &
    "    --nfinal <INT>    : Number of structures that are finally optimized", &
         "    --optlvl <method> : Defines the Method of final optimizations", &
         "    --maxgen <INT>    : Number of structure generation cycles", &
         "    --maxparent <INT> : Number of parents for structure generation", &
        "    --cs              : Decreases the computational costs, if first", &
         "                        geometry is CS symmetric", &
         !   "    --sphere <REAL>   : Spherical wall potential with radius",&
         !   "    --ellips <REAL> <REAL> <REAL> : Ellipsoid wall potential given",&
         !   "                                    with three unit axis radii",&
         ""

   end subroutine dockingHelp

   subroutine check_for_files(env, molA, molB)
      use xtb_readin, only: getValue
      implicit none
      !> Calculation environment
      type(TEnvironment), intent(inout) :: env
      !> Molecular structure data
      type(TMolecule), intent(inout) :: molA, molB

      character(len=*), parameter :: source = "iff_file_read"

      integer :: ich, err
      character(len=:), allocatable :: cdum     ! temporary string
      real(wp) :: charge
      integer :: elect

      call open_file(ich, '.CHRG', 'r')
      if (ich .ne. -1) then
         !> Total Charge
         call getline(ich, cdum, iostat=err)
         if (err /= 0) then
            call env%error('.CHRG is empty!', source)
         else
            if (.not. getValue(env, cdum, charge)) then
               call env%error('.CHRG has a problem')
            end if
         end if
         !> Charge molA
         call getline(ich, cdum, iostat=err)
         if (err /= 0) then
            call env%warning('.CHRG has only one line!')
         else
            if (getValue(env, cdum, charge)) then
               molA%chrg = charge
            end if
         end if
         !> Charge molB
         call getline(ich, cdum, iostat=err)
         if (err /= 0) then
            call env%warning('.CHRG has only two lines!')
         else
            if (getValue(env, cdum, charge)) then
               molB%chrg = charge
            end if
         end if

         call close_file(ich)
      end if

      call env%checkpoint("Reading charge from file failed")

      !> Number of unpaired electrons
      call open_file(ich, '.UHF', 'r')
      if (ich .ne. -1) then
         call getline(ich, cdum, iostat=err)
         if (err /= 0) then
            call env%error('.UHF is empty!', source)
         else
            if (getValue(env, cdum, elect)) then
               molA%uhf = elect
            end if
         end if
         call getline(ich, cdum, iostat=err)
         if (err /= 0) then
            call env%warning('.UHF has only one line!')
         else
            if (getValue(env, cdum, elect)) then
               molB%uhf = elect
            end if
         end if
         call close_file(ich)
      end if

   end subroutine check_for_files

   subroutine rdcontrol_iff(fname, env, copy_file)
      use xtb_readin, only: find_new_name
      use xtb_splitparam, only: maxfrag
      use xtb_scanparam, only: maxconstr, maxscan
      use xtb_sphereparam, only: maxwalls
      use xtb_readin, only: mirror_line
      implicit none
      character(len=*), parameter :: source = 'set_rdcontrol'
      character(len=*), intent(in)  :: fname
      type(TEnvironment), intent(inout) :: env
      character(len=:), allocatable :: line
      character(len=:), allocatable :: key
      character(len=:), allocatable :: val
      character(len=:), allocatable :: newname
      logical, intent(in), optional  :: copy_file
      character, parameter :: flag = '$'
      integer :: i
      integer :: id
      integer :: ic
      integer :: ie
      integer :: ncount
      integer :: copy
      integer :: err
      logical :: exist
      logical :: do_copy
      logical :: exitRun

      if (present(copy_file)) then
         do_copy = copy_file
      else
         do_copy = .false.
      end if

      call open_file(id, fname, 'r')
      if (id .eq. -1) then
         call env%warning("could not find '"//fname//"'", source)
         return
      end if

      if (do_copy) then
         newname = find_new_name(fname)
         call open_file(copy, newname, 'w')
      else
         copy = -1 ! deactivate copy in mirror_line
      end if

!  read first line before the readloop starts, I have to do this
!  to avoid using backspace on id (dammit Turbomole format)
      call mirror_line(id, copy, line, err)
      readflags: do
         !  check if there is a $ in the *first* column
         if (index(line, flag) .eq. 1) then
            select case (line(2:))
               ! logical
            case ('fit'); call set_fit; call mirror_line(id, copy, line, err)
            case ('samerand'); call set_samerand; call mirror_line(id, copy, line, err)
            case ('cma'); call set_cma; call mirror_line(id, copy, line, err)
               ! data
            case ('cube'); call rdblock(env, set_cube, line, id, copy, err, ncount)
            case ('write'); call rdblock(env, set_write, line, id, copy, err, ncount)
            case ('gfn'); call rdblock(env, set_gfn, line, id, copy, err, ncount)
            case ('scc'); call rdblock(env, set_scc, line, id, copy, err, ncount)
            case ('opt'); call rdblock(env, set_opt, line, id, copy, err, ncount)
            case ('gbsa'); call rdblock(env, set_gbsa, line, id, copy, err, ncount)
            case ('solvation'); call rdblock(env, set_gbsa, line, id, copy, err, ncount)
            case ('thermo'); call rdblock(env, set_thermo, line, id, copy, err, ncount)
            case ('path'); call rdblock(env, set_path, line, id, copy, err, ncount)
            case ('fix'); call rdblock(env, set_fix, line, id, copy, err, ncount)
            case ('wall'); call rdblock(env, set_wall, line, id, copy, err, ncount)
               maxwalls = maxwalls + ncount
            case ('constrain'); call rdblock(env, set_constr, line, id, copy, err, ncount)
               maxconstr = maxconstr + ncount; constraint_xyz = .true.
            case ('dock'); call rdblock_docking(env, set_docking, line, id, copy, err, ncount)
            case default ! unknown keyword -> ignore, we don't raise them
!           get a new line
               call mirror_line(id, copy, line, err)
            end select
         else ! not a keyword -> ignore
            call mirror_line(id, copy, line, err)
         end if
         !  check for end of file, which I will tolerate as alternative to $end
         if (is_iostat_end(err)) exit readflags
         call env%check(exitRun)
         if (exitRun) then
            call env%error("processing of data group failed", source)
            exit
         end if
      end do readflags

      if (do_copy) call close_file(copy)
      call close_file(id)
   end subroutine rdcontrol_iff


   subroutine read_userdata_iff(fname,env,mol)
      use xtb_readin, only : find_new_name
      !use xtb_scanparam
      use xtb_type_identitymap, only : TIdentityMap,init
      implicit none
      character(len=*), parameter :: source = 'userdata_read'
      type(TEnvironment), intent(inout) :: env
      type(TMolecule), intent(inout) :: mol
      character(len=*),intent(in)  :: fname
      character(len=:),allocatable :: line
      character(len=:),allocatable :: key
      character(len=:),allocatable :: val
      character(len=:),allocatable :: newname
      type(TIdentityMap) :: idMap
      integer :: i
      integer :: id
      integer :: ic
      integer :: ie
      integer :: err
      logical :: exist
      character, parameter :: flag = '$'

      if (set%verbose) then
         write(env%unit,'(72("$"))')
         write(env%unit,'(1x,"CONSTRAINTS & SCANS: DEBUG SECTION")')
         write(env%unit,'(72("$"))')
      endif

      call open_file(id,fname,'r')
      if (id.eq.-1) then
         call env%warning("could not find '"//fname//"'",source)
         return
      endif
      rewind(id) ! not sure if this is necessary

      call init(idMap, mol)

   !  read first line before the readloop starts, I have to do this
   !  to avoid using backspace on id (dammit Turbomole format)
      call getline(id,line,err)
      readflags: do
           !  check if there is a $ in the *first* column
         if (index(line,flag).eq.1) then
            select case(line(2:))
            case('directed'      )
               if (set%verbose) write(env%unit,'(">",1x,a)') line(2:)
               call rdblock_docking2(env,set_directed,line,id,mol%n,mol%at,idMap,mol%xyz,err)
            case default ! unknown keyword -> ignore, we don't raise them
               call getline(id,line,err)
            end select
         else ! not a keyword -> ignore
            call getline(id,line,err)
         endif
      !  check for end of file, which I will tolerate as alternative to $end
         if (is_iostat_end(err)) exit readflags
   !     if (index(line,flag_end).ne.0) exit readflags ! compatibility reasons
      enddo readflags

      if (set%verbose) write(env%unit,'(72("$"))')
      call close_file(id)
   end subroutine read_userdata_iff

   subroutine get_repulsive_pot(env,xyz,comb)
      implicit none
      !> Calculation environment
      type(TEnvironment), intent(inout) :: env

      !> Combined structure
      type(TMolecule), intent(in) :: comb
      real(wp), intent(in) :: xyz(3, comb%n)

      real(wp) :: dist, min_dist, rep_pot
      integer :: i, j

      allocate(directedset%val(comb%n),directedset%expo(comb%n), source=0.0_wp)

      ! First get the smalles distance for each atom and saving it into directedset%val
      do i = 1, comb%n
         if(any(i == directedset%atoms)) cycle !Distance for atoms in defined docking region
         min_dist = 0.0_wp
         do j = 1, directedset%n
            dist = sqrt((xyz(1,i)-xyz(1,directedset%atoms(j)))**2 &
                 & +(xyz(2,i)-xyz(2,directedset%atoms(j)))**2 &
                 & +(xyz(3,i)-xyz(3,directedset%atoms(j)))**2)
            if(min_dist == 0.0_wp) then
               min_dist = dist
            elseif (dist < min_dist) then
               min_dist = dist
            end if
         end do
         directedset%val(i) = min_dist
      end do
      !> Changing the distance to a repulsive potential sitting on every atom other then
      !  the defined docking atoms. This potentail is a damped exponential increase.
      !  It is later in the energy calculation and RG screening added in sitance depdence to
      !  docked molecule via 1/r²
      do i=1, comb%n
         if(any(i == directedset%atoms)) cycle !Potential zero for atoms in defined docking region
         dist = directedset%val(i)
         !rep_pot = exp(dist - 8) * (1 / (8000 + exp(dist - 8))) !Damped to a Max of 1 Hartree
         rep_pot = 0.1*erf(0.07 * dist - 0.28) !Potential starts at distance of 4
         if(rep_pot < 0.0_wp) rep_pot = 0.0_wp
         directedset%val(i) = rep_pot !Overwrite distance with repulsive Potential
      end do
   end subroutine get_repulsive_pot

   subroutine get_attractive_pot(env,comb)
      implicit none
      !> Calculation environment
      type(TEnvironment), intent(inout) :: env

      !> Combined structure
      type(TMolecule), intent(in) :: comb
      integer :: i

      allocate(directedset%val(comb%n),directedset%expo(comb%n), source=0.0_wp)

      ! First get the smalles distance for each atom and saving it into directedset%val
      do i = 1, comb%n
         if(any(i == directedset%atoms)) then
           directedset%val(i) = attractive_pot !attractive pot is negative
         else 
           directedset%val(i) = 0.0_wp
         end if
      end do
   end subroutine get_attractive_pot

end module xtb_prog_dock
