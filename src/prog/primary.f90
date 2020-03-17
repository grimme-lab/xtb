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
   use xtb_prog_main, only : xtbMain
   use xtb_type_environment

   implicit none

   !> Command line argument parser
   type(TArgParser) :: argParser

   !> Calculation environment
   type(TEnvironment) :: env

   !> start by initializing the MCTC library
   call mctc_init('xtb',10,.true.)

   !> Initialize calculation environment
   call init(env)

   !> Initialize argument parser from command line
   call init(argParser)

   !> Run the main program
   call xtbMain(env, argParser)

end program xtb_prog_primary
