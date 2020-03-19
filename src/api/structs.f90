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

!> Defines the interface to xtb datatypes and allows exporting them to C.
!  Mostly those will be carried around as opaque void* pointer and can
!  be manipulated with methods provided by the API.
module xtb_api_structs
   use iso_c_binding
   use xtb_mctc_accuracy, only : wp
   use xtb_api_utils
   implicit none

contains

!> Constructor for the molecular structure type, returns a nullptr in case
!  the generation fails. All values have to be provided.
type(c_ptr) function new_xtb_molecule_api &
      & (natoms, numbers, positions, charge, uhf, lattice, periodic) &
      & result(c_mol) bind(C, name="new_xTB_molecule")
   use xtb_type_molecule
   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: numbers(natoms)
   real(c_double), intent(in) :: positions(3, natoms)
   real(c_double), intent(in) :: charge
   integer(c_int), intent(in) :: uhf
   real(c_double), intent(in) :: lattice(3, 3)
   logical(c_bool), intent(in) :: periodic(3)

   type(TMolecule), pointer :: mol
   c_mol = c_null_ptr

   allocate(mol)

   mol = TMolecule(natoms, numbers, positions, charge, uhf, lattice, periodic)
   if (verify_xtb_molecule(mol) == 0) then
      c_mol = c_loc(mol)
   else
      call mol%deallocate
      deallocate(mol)
      nullify(mol)
   end if

end function new_xtb_molecule_api

!> Updates the molecular structure type. This routine checks for the allocation
!  state of the opaque pointer before copying positions and lattice.
integer(c_int) function update_xtb_molecule_api &
      & (c_mol, c_positions, c_lattice) &
      & result(status) bind(C, name="update_xTB_molecule")
   use xtb_type_molecule
   type(c_ptr), value, intent(in) :: c_mol
   type(c_ptr), value, intent(in) :: c_positions
   type(c_ptr), value, intent(in) :: c_lattice

   type(TMolecule), pointer :: mol
   real(c_double), pointer :: positions(:, :)
   real(c_double), pointer :: lattice(:, :)

   call c_f_pointer(c_mol, mol)
   if (mol%n > 0) then
      if (c_associated(c_positions)) then
         call c_f_pointer(c_positions, positions, [3, mol%n])
         mol%xyz = positions
      end if
      if (c_associated(c_lattice)) then
         call c_f_pointer(c_lattice, lattice, [3, 3])
         mol%lattice = lattice
      end if
      call mol%update
      ! check for very short distances
      status = verify_xtb_molecule(mol)
   else
      status = -1
   end if

end function update_xtb_molecule_api

!> Deconstructor for the molecular structure type.
subroutine delete_xtb_molecule_api(c_mol) &
      & bind(C, name="delete_xTB_molecule")
   use xtb_type_molecule
   type(c_ptr), value, intent(in) :: c_mol
   type(TMolecule), pointer :: mol

   if (c_associated(c_mol)) then
      call c_f_pointer(c_mol, mol)
      call mol%deallocate
      deallocate(mol)
      nullify(mol)
   end if

end subroutine delete_xtb_molecule_api

!> Constructor for the tight binding basisset.
type(c_ptr) function new_xtb_basisset_api(c_mol) &
      & result(c_basis) bind(C, name="new_xTB_basisset")
   use xtb_type_molecule
   use xtb_type_basisset
   type(c_ptr), value, intent(in) :: c_mol
   type(TBasisset), pointer :: basis
   type(TMolecule), pointer :: mol
   integer(c_int) :: status
   c_basis = c_null_ptr

   if (check_xtb_init()) then

      call c_f_pointer(c_mol, mol)

      if (mol%n > 0) then
         allocate(basis)

         call new_xtb_basisset(mol, basis, status)
         if (status == 0_c_int) then
            c_basis = c_loc(basis)
         else
            deallocate(basis)
         end if
      end if
   end if

end function new_xtb_basisset_api

!> Actual implementation of the constructor.
subroutine new_xtb_basisset(mol, basis, status)
   use xtb_setparam, only : gfn_method
   use xtb_type_molecule
   use xtb_type_basisset
   use xtb_basis
   use xtb_xtb_data
   use xtb_xtb_gfn0
   use xtb_xtb_gfn1
   use xtb_xtb_gfn2
   type(TMolecule), intent(in) :: mol
   type(TBasisset), intent(inout) :: basis
   integer(c_int), intent(out) :: status
   type(TxTBData) :: xtbData
   logical :: diff, okbas
   status = 2
   select case(gfn_method)
   case default
      okbas = .false.
   case(0)
      call initGFN0(xtbData)
      call xbasis0(xtbData, mol%n, mol%at, basis)
      call xbasis_gfn0(xtbData, mol%n, mol%at, basis, okbas, diff)
   case(1)
      call initGFN1(xtbData)
      call xbasis0(xtbData, mol%n, mol%at, basis)
      call xbasis_gfn1(xtbData, mol%n, mol%at, basis, okbas, diff)
   case(2)
      call initGFN2(xtbData)
      call xbasis0(xtbData, mol%n, mol%at, basis)
      call xbasis_gfn2(xtbData, mol%n, mol%at, basis, okbas)
   end select
   if (okbas) then
      status = 0
   else
      call basis%deallocate
   end if
end subroutine new_xtb_basisset

!> Return dimensions of this tight binding basisset.
!  For not allocated input zero is returned, every argument is optional.
subroutine dimensions_xtb_basisset_api &
      & (c_basis, nat, nbf, nao, nsh) &
      & bind(C, name="dimensions_xTB_basisset")
   use xtb_type_basisset
   type(c_ptr), value, intent(in) :: c_basis
   type(TBasisset), pointer :: basis
   integer(c_int), intent(out), optional :: nat
   integer(c_int), intent(out), optional :: nbf
   integer(c_int), intent(out), optional :: nao
   integer(c_int), intent(out), optional :: nsh

   if (c_associated(c_basis)) then
      call c_f_pointer(c_basis, basis)
      call c_return(nat, basis%n)
      call c_return(nbf, basis%nbf)
      call c_return(nao, basis%nao)
      call c_return(nsh, basis%nshell)
   else
      call c_return(nat, 0)
      call c_return(nbf, 0)
      call c_return(nao, 0)
      call c_return(nsh, 0)
   end if

end subroutine dimensions_xtb_basisset_api

!> Deconstructor for the tight binding wavefunction.
subroutine delete_xtb_basisset_api(c_basis) &
      & bind(C, name="delete_xTB_basisset")
   use xtb_type_basisset
   type(c_ptr), value, intent(in) :: c_basis
   type(TBasisset), pointer :: basis

   if (c_associated(c_basis)) then
      call c_f_pointer(c_basis, basis)
      call basis%deallocate
      deallocate(basis)
      nullify(basis)
   end if

end subroutine delete_xtb_basisset_api

!> Constructor for the tight binding wavefunction.
type(c_ptr) function new_xtb_wavefunction_api &
      & (c_mol, c_basis) &
      & result(c_wfn) bind(C, name="new_xTB_wavefunction")
   use xtb_type_molecule
   use xtb_type_basisset
   use xtb_type_wavefunction
   type(c_ptr), value, intent(in) :: c_mol
   type(c_ptr), value, intent(in) :: c_basis
   type(TMolecule), pointer :: mol
   type(TBasisset), pointer :: basis
   type(TWavefunction), pointer :: wfn
   c_wfn = c_null_ptr

   if (check_xtb_init()) then
      if (c_associated(c_basis) .and. c_associated(c_mol)) then
         call c_f_pointer(c_mol, mol)
         call c_f_pointer(c_basis, basis)
         if (mol%n > 0 .and. basis%n > 0 .and. basis%nbf > 0 .and. &
            &basis%nao > 0 .and. basis%nshell > 0) then
            allocate(wfn)

            call wfn%allocate(basis%n, basis%nshell, basis%nao)

            wfn%nel = nint(sum(mol%z) - mol%chrg)
            wfn%nopen = mol%uhf
            if (mod(wfn%nopen, 2) == 0 .and. mod(wfn%nel, 2) /= 0) wfn%nopen = 1
            if (mod(wfn%nopen, 2) /= 0 .and. mod(wfn%nel, 2) == 0) wfn%nopen = 0
            c_wfn = c_loc(wfn)
         end if
      end if
   end if

end function new_xtb_wavefunction_api

!> Return dimensions of this tight binding wavefunction.
subroutine dimensions_xtb_wavefunction_api &
      & (c_wfn, nat, nao, nsh) &
      & bind(C, name="dimensions_xTB_wavefunction")
   use xtb_type_wavefunction
   type(c_ptr), value, intent(in) :: c_wfn
   type(TWavefunction), pointer :: wfn
   integer(c_int), intent(out), optional :: nat
   integer(c_int), intent(out), optional :: nao
   integer(c_int), intent(out), optional :: nsh

   if (c_associated(c_wfn)) then
      call c_f_pointer(c_wfn, wfn)
      call c_return(nat, wfn%n)
      call c_return(nao, wfn%nao)
      call c_return(nsh, wfn%nshell)
   else
      call c_return(nat, 0)
      call c_return(nao, 0)
      call c_return(nsh, 0)
   end if

end subroutine dimensions_xtb_wavefunction_api

!> Query the wavefunction properties, pass NULL/nullptr in case you want
!  a query to be ignored.
!  In case the wavefunction is not allocated nothing is returned.
subroutine query_xtb_wavefunction_api &
      & (c_wfn, charges, dipoles, quadrupoles, bond_orders, &
      &  hl_gap, orbital_energies) &
      & bind(C, name="query_xTB_wavefunction")
   use xtb_mctc_convert, only : evtoau
   use xtb_type_wavefunction
   type(c_ptr), value, intent(in) :: c_wfn
   type(TWavefunction), pointer :: wfn
   real(c_double), intent(out), optional :: charges(*)
   real(c_double), intent(out), optional :: dipoles(3, *)
   real(c_double), intent(out), optional :: quadrupoles(6, *)
   real(c_double), intent(out), optional :: bond_orders(*)
   real(c_double), intent(out), optional :: hl_gap
   real(c_double), intent(out), optional :: orbital_energies(*)

   if (c_associated(c_wfn)) then
      call c_f_pointer(c_wfn, wfn)
      if (wfn%n > 0 .and. wfn%nao > 0 .and. wfn%nshell > 0) then
         if (present(charges)) charges(1:wfn%n) = wfn%q
         if (present(dipoles)) dipoles(1:3, 1:wfn%n) = wfn%dipm
         if (present(quadrupoles)) quadrupoles(1:6, 1:wfn%n) = wfn%qp
         if (present(bond_orders)) then
            bond_orders(1:wfn%n**2) = reshape(wfn%wbo, [wfn%n**2])
         end if
         if (present(orbital_energies)) then
            orbital_energies(1:wfn%nao) = wfn%emo * evtoau
         end if
         if (present(hl_gap) .and. wfn%ihomo+1 <= wfn%nao .and. wfn%ihomo >= 1) then
            hl_gap = (wfn%emo(wfn%ihomo+1) - wfn%emo(wfn%ihomo)) * evtoau
         end if
      end if
   end if

end subroutine query_xtb_wavefunction_api

!> Deconstructor for the tight binding basisset.
subroutine delete_xtb_wavefunction_api(c_wfn) &
      & bind(C, name="delete_xTB_wavefunction")
   use xtb_type_wavefunction
   type(c_ptr), value, intent(in) :: c_wfn
   type(TWavefunction), pointer :: wfn

   if (c_associated(c_wfn)) then
      call c_f_pointer(c_wfn, wfn)
      call wfn%deallocate
      deallocate(wfn)
      nullify(wfn)
   end if

end subroutine delete_xtb_wavefunction_api

end module xtb_api_structs
