! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> Force field topology generator
module xtb_prog_topology
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_filetypes, only : getFileType, fileType
   use xtb_mctc_systools, only : rdpath, getline
   use xtb_mctc_timings
   use xtb_mctc_version, only : version, author, date
   use xtb_io_reader, only : readMolecule, readHessian
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_reader, only : TReader
   use xtb_prog_argparser, only : TArgParser
   use xtb_gfnff_calculator, only : TGFFCalculator
   use xtb_main_setup, only : newGFFCalculator
   use xtb_setparam, only : xenv, ichrg
   use xtb_setmod, only : set_chrg
   implicit none
   private

   public :: xtbTopology


contains


!> Main program of the topology submodule
subroutine xtbTopology(env, argParser)

   !> Source of errors in the main program unit
   character(len=*), parameter :: source = "prog_topology"

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Commandline argument parser
   type(TArgParser), intent(inout) :: argParser

   !> Molecular structure data
   type(TMolecule) :: mol

   !> Force field instance
   type(TGFFCalculator) :: calc

   type(TReader) :: reader
   character(len=:), allocatable :: file, param, line
   integer :: nFiles, ftype, unit, err
   logical :: exist

   call parseArguments(env, argParser, param)

   nFiles = argParser%countFiles()
   if (nFiles == 0) then
      call env%error("No input file given, so there is nothing to do", source)
      call topologyHelp(env%unit)
   end if
   if (nFiles > 1) then
      call env%error("Multiple input files present, aborting", source)
   end if

   call env%checkpoint("Command line argument parsing failed")

   call argParser%nextFile(file)

   !> Determine the file type for processing in the reader
   ftype = getFileType(file)

   !> Generate the reader by connecting the file to the instance
   call reader%open(file)
   !> Process the file to obtain the molecular structure data
   call readMolecule(env, mol, reader%unit, ftype)
   !> Close the file, we are done with it here
   call reader%close

   call open_file(unit, '.CHRG', 'r')
   if (unit /= -1) then
      call getline(unit, line, iostat=err)
      if (err /= 0) then
         call env%error('.CHRG is empty!', source)
      else
         call set_chrg(env, line)
         mol%chrg = ichrg
         call close_file(unit)
      end if
   end if

   call env%checkpoint("Could not read geometry from '"//file//"'")

   !> Print an informative banner
   call topologyHeader(env%unit)
   !> print current time
   call prdate('S')

   if (.not.allocated(param)) then
      call rdpath(env%xtbpath, '.param_gfnff.xtb', param, exist)
      if (.not.exist) param = '.param_gfnff.xtb'
   end if

   call newGFFCalculator(env, mol, calc, param, .false.)

   call env%checkpoint("Could not setup force field topology")

   write(env%unit,'(a)')
   write(env%unit,'(72("-"))')
   call stop_timing_run
   call stop_timing(1)
   call prdate('E')
   write(env%unit,'(72("-"))')
   call prtiming(1,'total')

   write(env%unit,'(a)')
   call terminate(0)

end subroutine xtbTopology


!> Parse command line arguments
subroutine parseArguments(env, args, param)
   use xtb_mctc_global, only : persistentEnv

   !> Name of error producer
   character(len=*), parameter :: source = "prog_topology_parseArguments"

   !> Calculation environment
   type(TEnvironment) :: env

   !> Command line argument parser
   type(TArgParser) :: args

   !> Parameter file
   character(len=:), allocatable, intent(out) :: param

   integer :: nFlags
   character(len=:), allocatable :: flag, sec
   logical :: exist

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

      case('--help', '-h')
         call topologyHelp(env%unit)
         call terminate(0)

      case('--param')
         call args%nextArg(param)
         if (.not.allocated(param)) then
            call env%error("Missing argument for parameter file", source)
         else
            inquire(exist=exist, file=param)
            if (.not.exist) then
               call env%error("Parameter file does not exist", source)
            end if
         end if

      case('--namespace')
         call args%nextArg(persistentEnv%io%namespace)
         if (.not.allocated(persistentEnv%io%namespace)) then
            call env%error("Namespace argument is missing", source)
         end if

      end select
      call args%nextFlag(flag)
   end do

end subroutine parseArguments


!> Header for this submodule
subroutine topologyHeader(unit)

   !> IO unit
   integer, intent(in) :: unit

   write(unit,'(a)') &
      "      -----------------------------------------------------------",&
      "     |                   =====================                   |",&
      "     |                      x t b - t o p o                      |",&
      "     |                   =====================                   |",&
      "     |                         S. Grimme                         |",&
      "     |          Mulliken Center for Theoretical Chemistry        |",&
      "     |                    University of Bonn                     |",&
      "      -----------------------------------------------------------",""
   write(unit,'(3x,"*",*(1x,a))') &
      & "xtb-topology version", version, "compiled by", author, "on", date
   write(unit,'(a)')

end subroutine topologyHeader


subroutine topologyHelp(unit)

   !> IO unit
   integer, intent(in) :: unit

   write(unit, '(a)') &
   "Usage: xtb topology [options] <geometry>", &
   "",&
   "<geometry> may be provided as any valid input to xtb", &
   "",&
   "Options",&
   "",&
   "   --namespace STRING   Set a namespace for the calculation",&
   "",&
   "   --param FILE         Use alternative parameter file for the generator",&
   "",&
   "More information can be obtained at https://xtb-docs.readthedocs.io/",&
   ""
end subroutine topologyHelp


end module xtb_prog_topology
