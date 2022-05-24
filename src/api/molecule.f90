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

!> API for dealing with molecular structure data
module xtb_api_molecule
   use, intrinsic :: iso_c_binding
   use xtb_mctc_accuracy, only : wp
   use xtb_api_environment
   use xtb_api_utils
   use xtb_type_molecule
   implicit none
   private

   public :: VMolecule
   public :: newMolecule_api, delMolecule_api, updateMolecule_api


   !> Void pointer to molecular structure data
   type :: VMolecule
      type(TMolecule) :: ptr
   end type VMolecule


   !> Maximum accepted atomic number, actual calculator might accept less
   integer, parameter :: maxElem = 118


contains


function newMolecule_api(venv, natoms, numbers, positions, charge, uhf, lattice, &
      & periodic) result(vmol) &
      & bind(C, name="xtb_newMolecule")
   character(len=*), parameter :: source = 'xtb_api_newMolecule'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: numbers(natoms)
   real(c_double), intent(in) :: positions(3, natoms)
   real(c_double), intent(in), optional :: charge
   integer(c_int), intent(in), optional :: uhf
   real(c_double), intent(in), optional :: lattice(3, 3)
   logical(c_bool), intent(in), optional :: periodic(3)
   type(VMolecule), pointer :: mol
   type(c_ptr) :: vmol
   integer(c_int) :: status

   vmol = c_null_ptr

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (natoms <= 0) then
         call env%ptr%error("Number of atoms must be positive", source)
         return
      end if

      if (any(numbers <= 0) .or. any(numbers > maxElem)) then
         call env%ptr%error("Invalid atomic number present", source)
         return
      end if

      allocate(mol)

      mol%ptr = TMolecule(natoms, numbers, positions, charge, uhf, &
         & lattice, periodic)

      status = verifyMolecule(mol%ptr)
      if (status /= 0_c_int) then
         deallocate(mol)
         call env%ptr%error("Could not generate molecular structure", source)
         return
      end if

      vmol = c_loc(mol)
   end if

end function newMolecule_api


subroutine delMolecule_api(vmol) &
      & bind(C, name="xtb_delMolecule")
   type(c_ptr), intent(inout) :: vmol
   type(VMolecule), pointer :: mol

   if (c_associated(vmol)) then
      call c_f_pointer(vmol, mol)
      call checkGlobalEnv

      deallocate(mol)
      vmol = c_null_ptr
   end if

end subroutine delMolecule_api


subroutine updateMolecule_api(venv, vmol, positions, lattice) &
      & bind(C, name="xtb_updateMolecule")
   character(len=*), parameter :: source = 'xtb_api_updateMolecule'
   type(c_ptr), value :: venv
   type(VEnvironment), pointer :: env
   type(c_ptr), value :: vmol
   type(VMolecule), pointer :: mol
   real(c_double), intent(in) :: positions(3, *)
   real(c_double), intent(in), optional :: lattice(3, 3)
   real(wp) :: latvecs(3, 3)
   integer(c_int) :: status

   if (c_associated(venv)) then
      call c_f_pointer(venv, env)
      call checkGlobalEnv

      if (.not.c_associated(vmol)) then
         call env%ptr%error("Molecular structure data is not allocated", source)
         return
      end if

      call c_f_pointer(vmol, mol)

      if (mol%ptr%n <= 0 .or. mol%ptr%nid <= 0 .or. .not.allocated(mol%ptr%at) &
         & .or. .not.allocated(mol%ptr%id) .or. .not.allocated(mol%ptr%xyz)) then
         call env%ptr%error("Invalid molecular structure data provided", source)
         return
      end if

      mol%ptr%xyz(:, :) = positions(1:3, 1:mol%ptr%n)
      if (present(lattice)) then
         mol%ptr%lattice(:, :) = lattice(1:3, 1:3)
      end if

      call mol%ptr%update

      status = verifyMolecule(mol%ptr)
      if (status /= 0_c_int) then
         call env%ptr%error("Could not update molecular structure", source)
         return
      end if

   end if

end subroutine updateMolecule_api


end module xtb_api_molecule
