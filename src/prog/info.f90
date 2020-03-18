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

!> Information on the current input geometry
module xtb_prog_info
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autoamu
   use xtb_mctc_filetypes, only : getFileType
   use xtb_mctc_version, only : version, author, date
   use xtb_io_reader, only : readMolecule
   use xtb_type_environment, only : TEnvironment
   use xtb_type_identitymap, only : TIdentityMap, init
   use xtb_type_molecule, only : TMolecule, len
   use xtb_type_reader, only : TReader
   use xtb_param_atomicmass, only : getAtomicMass
   use xtb_prog_argparser, only : TArgParser
   use xtb_splitparam, only : atmass
   implicit none
   private

   public :: xtbInfo


contains


!> Main program of the info submodule
subroutine xtbInfo(env, argParser)

   !> Source of errors in the main program unit
   character(len=*), parameter :: source = "prog_info"

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Commandline argument parser
   type(TArgParser), intent(inout) :: argParser

   !> Molecular structure data
   type(TMolecule) :: mol

   !> Reader for input files
   type(TReader) :: reader

   !> Identity map for unique species in molecular structure
   type(TIdentityMap) :: idMap

   !> Read in file name
   character(len=:), allocatable :: file

   !> File type
   integer :: ftype

   !> Current run
   integer :: iRun

   !> Number of total runs
   integer :: nRun

   !> Signal from environment to terminate the program
   logical, allocatable :: exitRun(:)

   !> Check for very close distances
   logical :: coldFusion

   !> Determine the number of runs from the files in the command line
   nRun = argParser%countFiles()
   !> No input files are an error
   if (nRun == 0) then
      call env%terminate("No geometry input provided")
   end if

   !> Print an informative banner
   call infoHeader(env%unit)

   !> Save status of each run
   allocate(exitRun(nRun))

   !> Execute runs in order
   do iRun = 1, nRun

      !> Pop the first file from the parser
      call argParser%nextFile(file)
      !> Determine the file type for processing in the reader
      ftype = getFileType(file)
      !> Generate the reader by connecting the file to the instance
      call reader%open(file)
      !> Process the file to obtain the molecular structure data
      call readMolecule(env, mol, reader%unit, ftype)
      !> Close the file, we are done with it here
      call reader%close

      !> Print an informative header for this run
      call runHeader(env%unit, file, iRun, nRun, len(mol))

      call env%check(exitRun(iRun))
      if (exitRun(iRun)) then
         call env%show('Could not read geometry due to')
         cycle
      end if

      !> make sure the geometry summery does not error on splitparam globals
      atmass = getAtomicMass(mol%at) * autoamu ! from splitparam.f90

      !> Intialize and print the identity map
      call init(idMap, mol)
      call generic_header(env%unit, 'Identifiers', 49, 10)
      write(env%unit, '(a)')
      call idMap%writeInfo(env%unit)
      write(env%unit, '(a)')

      !> Check for very small distances in the input geometry
      call check_cold_fusion(env, mol, coldFusion)
      if (coldFusion) then
         call env%error("Some atoms in the geometry are very close", source)
      else
         !> Geometry summary
         call main_geometry(env%unit, mol)
      end if

      call env%check(exitRun(iRun))
      call env%show('Following warnings occurred with this input')

   end do

   if (nRun > 1) then
      write(env%unit, '(72("~"))')
      write(env%unit, '(1x, a, 1x, i0, 1x, a)') "Processed", nRun, "individual inputs"
      if (count(exitRun) > 0) then
         write(env%unit, '(1x, i0, "/", i0, 1x, a, "(", f5.1, "%)")') &
            & count(exitRun), nRun, "inputs failed info check", &
            & real(count(exitRun)*100, wp)/real(nRun, wp)
      end if
      write(env%unit, '(72("~"), /)')
   end if
   if (any(exitRun)) then
      call terminate(1)
   else
      call terminate(0)
   end if


end subroutine xtbInfo


!> Print a header for each new run
subroutine runHeader(unit, file, iRun, nRun, nAtom)

   !> IO unit
   integer, intent(in) :: unit

   !> Read in file name
   character(len=*), intent(in) :: file

   !> Current run
   integer, intent(in) :: iRun

   !> Number of total runs
   integer, intent(in) :: nRun

   !> Number of atoms in this structure
   integer, intent(in) :: nAtom

   if (nRun > 1) then
      write(unit, '(72("~"))')
      write(unit, '(1x, "Run:", t66, i0, "/", i0)') iRun, nRun
      write(unit, '(72("~"))')
   end if
   write(unit, '(1x, "Input:", 1x, a)') file
   write(unit, '(1x, "Number of atoms:", 1x, i0, /)') nAtom

end subroutine runHeader


!> Header for this submodule
subroutine infoHeader(unit)

   !> IO unit
   integer, intent(in) :: unit

   write(unit,'(a)') &
      "      -----------------------------------------------------------",&
      "     |                   =====================                   |",&
      "     |                      x t b - i n f o                      |",&
      "     |                   =====================                   |",&
      "     |                         S. Ehlert                         |",&
      "     |          Mulliken Center for Theoretical Chemistry        |",&
      "     |                    University of Bonn                     |",&
      "      -----------------------------------------------------------",""
   write(unit,'(3x,"*",*(1x,a))') &
      & "xtb version", version, "compiled by", author, "on", date
   write(unit,'(a)')

end subroutine infoHeader


end module xtb_prog_info
