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

!> The primary driver for the xtb program package
program xtb_prog_primary
   use xtb_prog_argparser
   use xtb_prog_dock, only : xtbDock
   use xtb_prog_info, only : xtbInfo
   use xtb_prog_main, only : xtbMain
   use xtb_prog_thermo, only : xtbThermo
   use xtb_prog_topology, only : xtbTopology
   use xtb_prog_ir, only: xtbIR
   use xtb_prog_submodules
   use xtb_type_environment

   implicit none

   !> Command line argument parser
   type(TArgParser) :: argParser

   !> Calculation environment
   type(TEnvironment) :: env

   !> Requested run mode
   integer :: runMode

   !> start by initializing the MCTC library
   call mctc_init('xtb',10,.true.)

   !> Initialize calculation environment
   call init(env)

   !> Initialize argument parser from command line
   call init(argParser)

   !> Get the requested run mode
   call getRunmode(argParser, runMode)

   !> Select the correct submodule to run
   select case(runMode)
   case(xtbSubmodule%main)
      !> Run the main program
      call xtbMain(env, argParser)

   case(xtbSubmodule%info)
      !> Run the info submodule
      call xtbInfo(env, argParser)

   case(xtbSubmodule%thermo)
      !> Run the thermo submodule
      call xtbThermo(env, argParser)

   case(xtbSubmodule%topo)
      !> Run the thermo submodule
      call xtbTopology(env, argParser)

   case(xtbSubmodule%dock)
      !> Run the docking submodule
      call xtbDock(env, argParser)

   case(xtbSubmodule%ir)
      !> Run the IR submodule
      call xtbIR(env, argParser)
   end select

contains

subroutine getRunmode(argParser, runMode)

   !> Command line argument parser
   type(TArgParser), intent(inout) :: argParser

   !> Requested run mode
   integer, intent(out) :: runMode

   !> First command line argument
   character(len=:), allocatable :: argument

   !> Initialize run mode with default submodule
   runMode = xtbSubmodule%main

   !> Get first argument
   call argParser%nextArg(argument)

   !> Check if we have a command line argument
   if (allocated(argument)) then

      !> Get the identifier of the submodule
      runMode = getSubmodule(argument)

      !> In case of an invalid identifier, we reset the parser and use the default
      if (runMode == xtbSubmodule%invalid) then
         call argParser%reset
         runMode = xtbSubmodule%main
      end if

   end if

end subroutine getRunmode

end program xtb_prog_primary
