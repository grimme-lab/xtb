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

!> Generic wrappers for all the readers implemented
module xtb_io_reader
   use mctc_env, only : error_type
   use mctc_io, only : structure_type, read_structure
   use xtb_io_reader_genformat, only : readHessianDFTBPlus
   use xtb_io_reader_orca, only : readHessianOrca
   use xtb_io_reader_turbomole, only : readHessianTurbomole
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_filetypes, only : fileType, hessType
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule, assignment(=)
   use xtb_type_reader, only : TReader
   implicit none
   private

   public :: readMolecule, readHessian


contains


!> Generic reader for molecules from input files
subroutine readMolecule(env, mol, unit, ftype)

   !> Name of error producer
   character(len=*), parameter :: source = "io_reader_readMolecule"

   !> Calculation environment
   class(TEnvironment), intent(inout) :: env

   !> Instance of molecular structure data
   class(TMolecule), intent(out) :: mol

   !> File reader
   !class(TReader), intent(inout) :: reader
   integer, intent(in) :: unit

   !> Idenitifier for file type
   integer, intent(in) :: ftype

   type(structure_type) :: struc
   type(error_type), allocatable :: error

   call read_structure(struc, unit, ftype, error)
   if (allocated(error)) then
      call env%error(error%message)
      return
   end if

   if (count(struc%periodic) == 1) then
      call env%error("1D periodic structures are currently unsupported", source)
      return
   end if

   if (count(struc%periodic) == 2) then
      call env%error("2D periodic structures are currently unsupported", source)
      return
   end if

   if (allocated(struc%sdf)) then
      if (any(struc%sdf%hydrogens > 0)) then
         call env%error("Hydrogen atom queries in ctfiles are currently unsupported", source)
         return
      end if
   end if

   if (allocated(struc%pdb)) then
      if (.not.any(struc%num == 1)) then
         call env%error("PDB structure without hydrogen atoms found, "//&
            &"aborting due to incomplete input geometry", source)
         return
      end if
   end if

   mol = struc
   mol%ftype = ftype

end subroutine readMolecule


subroutine readHessian(env, mol, hessian, reader, format)
   character(len=*), parameter :: source = 'io_reader_readHessian'
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol
   type(TReader), intent(inout) :: reader
   integer, intent(in) :: format
   real(wp), intent(out) :: hessian(:, :)
   character(len=:), allocatable :: message
   logical :: status

   if (any(shape(hessian) /= [3*mol%n, 3*mol%n])) then
      call env%error("Shape of hessian array does not match geometry", source)
      return
   end if

   select case(format)
   case default
      message = "Unknown hessian format"
      status = .false.

   case(hessType%tmol)
      call readHessianTurbomole(hessian, reader, mol, status, message)

   case(hessType%orca)
      call readHessianOrca(hessian, reader, mol, status, message)

   case(hessType%dftbplus)
      call readHessianDFTBPlus(hessian, reader, mol, status, message)

   end select

   if (.not.status) then
      call env%error(message, source)
      return
   end if

end subroutine readHessian


end module xtb_io_reader
