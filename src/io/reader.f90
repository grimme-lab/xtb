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
   use xtb_io_reader_ctfile, only : readMoleculeMolfile, readMoleculeSDF
   use xtb_io_reader_gaussian, only : readMoleculeGaussianExternal
   use xtb_io_reader_genformat, only : readMoleculeGenFormat, &
      & readHessianDFTBPlus
   use xtb_io_reader_pdb, only : readMoleculePDB
   use xtb_io_reader_turbomole, only : readMoleculeCoord, readHessianTurbomole
   use xtb_io_reader_vasp, only : readMoleculeVasp
   use xtb_io_reader_xyz, only : readMoleculeXYZ
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_filetypes, only : fileType
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
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

   character(len=:), allocatable :: message
   logical :: status

   select case(ftype)
   case(fileType%xyz)
      call readMoleculeXYZ(mol, unit, status, iomsg=message)
   case(fileType%tmol)
      call readMoleculeCoord(mol, unit, status, iomsg=message)
   case(fileType%molfile)
      call readMoleculeMolfile(mol, unit, status, iomsg=message)
   case(fileType%sdf)
      call readMoleculeSDF(mol, unit, status, iomsg=message)
   case(fileType%vasp)
      call readMoleculeVasp(mol, unit, status, iomsg=message)
   case(fileType%pdb)
      call readMoleculePDB(mol, unit, status, iomsg=message)
   case(fileType%gen)
      call readMoleculeGenFormat(mol, unit, status, iomsg=message)
   case(fileType%gaussian)
      call readMoleculeGaussianExternal(mol, unit, status, iomsg=message)
   case default
      status = .false.
      message = "coordinate format not known"
   end select

   if (.not.status) then
      call env%error(message, source)
      return
   end if

   call mol%update

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

   case(fileType%tmol)
      call readHessianTurbomole(hessian, reader, mol, status, message)

   case(fileType%gen)
      call readHessianDFTBPlus(hessian, reader, mol, status, message)

   end select

   if (.not.status) then
      call env%error(message, source)
      return
   end if

end subroutine readHessian


end module xtb_io_reader
