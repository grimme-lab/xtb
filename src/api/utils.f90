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

!> Utilities to make working with C datatypes as convenient as with Fortran
!  datatypes.
module xtb_api_utils
   use iso_c_binding
   use xtb_mctc_accuracy, only : wp
   implicit none

   interface c_return
      module procedure :: c_return_int_0d
      module procedure :: c_return_double_0d
      module procedure :: c_return_double_1d
      module procedure :: c_return_double_2d
   end interface c_return

   interface c_return_alloc
      module procedure :: c_return_double_0d_alloc
      module procedure :: c_return_double_1d_alloc
      module procedure :: c_return_double_2d_alloc
   end interface c_return_alloc

   interface c_get
      module procedure :: c_get_double_0d_alloc
      module procedure :: c_get_double_1d_alloc
      module procedure :: c_get_double_2d_alloc
      module procedure :: c_get_double_0d
      module procedure :: c_get_double_1d
      module procedure :: c_get_double_2d
      module procedure :: c_get_int_0d
      module procedure :: c_get_int_1d
      module procedure :: c_get_int_2d
   end interface c_get

contains

!> Check if xTB was correctly initialized.
logical(c_bool) function check_xtb_init() &
      & result(status) bind(C, name="check_xTB_init")
   use xtb_setparam, only: gfn_method
   status = gfn_method >= 0 .and. gfn_method <= 2
end function check_xtb_init

!> Check if xTB was locked on the correct parametrisation.
logical(c_bool) function check_xtb_lock(gfn) &
      & result(status) bind(C, name="check_xTB_lock")
   use xtb_setparam, only: gfn_method
   !> Requested GFN level
   integer(c_int), intent(in) :: gfn
   status = gfn == gfn_method
end function check_xtb_lock

!> Release the lock on the xTB parametrisation.
subroutine reset_xtb_lock() bind(C, name="reset_xTB_lock")
   use xtb_setparam, only: gfn_method
   !$omp critical(xtb_load_api)
   gfn_method = -1
   !$omp end critical(xtb_load_api)
end subroutine

subroutine eeq_guess_wavefunction(mol, wfn, err)
   use xtb_setparam, only: gfn_method
   use xtb_mctc_logging
   use xtb_type_molecule
   use xtb_type_wavefunction
   use xtb_type_param
   use xtb_disp_ncoord
   use xtb_eeq
   use xtb_scc_core
   type(TMolecule), intent(in) :: mol
   type(TWavefunction), intent(inout) :: wfn
   type(mctc_error), allocatable :: err
   type(chrg_parameter) :: chrgeq
   real(wp), allocatable :: cn(:)

   ! do an EEQ guess
   allocate( cn(mol%n), source = 0.0_wp )
   call new_charge_model_2019(chrgeq,mol%n,mol%at)
   call ncoord_erf(mol%n,mol%at,mol%xyz,cn)
   call eeq_chrgeq(mol,err,chrgeq,cn,wfn%q)
   deallocate(cn)

   call iniqshell(mol%n,mol%at,mol%z,wfn%nshell,wfn%q,wfn%qsh,gfn_method)
end subroutine eeq_guess_wavefunction

!> Cold fusion check
integer function verify_xtb_molecule(mol) result(status)
   use xtb_type_molecule
   type(TMolecule), intent(in) :: mol
   integer :: iat, jat
   status = 0
   do iat = 1, mol%n
      do jat = 1, iat - 1
         if (mol%dist(jat, iat) < 1.0e-9_wp) status = status + 1
      end do
   end do
end function verify_xtb_molecule

logical function verify_xtb_basisset(mol, basis) result(status)
   use xtb_type_molecule
   use xtb_type_basisset
   use xtb_basis, only: dim_basis
   type(TMolecule), intent(in) :: mol
   type(TBasisset), intent(in) :: basis
   integer :: nshell, nao, nbf
   call dim_basis(mol%n, mol%at, nshell, nao, nbf)
   status = basis%n > 0 .and. basis%nbf > 0 .and. basis%nao > 0 &
      & .and. basis%nshell > 0 .and. mol%n == basis%n &
      & .and. nshell == basis%nshell .and. nao == basis%nao .and. nbf == basis%nbf
end function verify_xtb_basisset

logical function verify_xtb_wavefunction(basis, wfn) result(status)
   use xtb_type_basisset
   use xtb_type_wavefunction
   type(TBasisset), intent(in) :: basis
   type(TWavefunction), intent(in) :: wfn
   status = wfn%n > 0 .and. wfn%nao > 0 .and. wfn%nshell > 0 &
      & .and. wfn%n == basis%n .and. wfn%nshell == basis%nshell &
      & .and. wfn%nao == basis%nao
end function verify_xtb_wavefunction

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran value has been calculated
subroutine c_return_int_0d(c_array, f_array)
   integer(c_int), intent(out), target :: c_array
   integer, intent(in) :: f_array
   if (c_associated(c_loc(c_array))) then
      c_array = f_array
   endif
end subroutine c_return_int_0d

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran value has been calculated
subroutine c_return_double_0d(c_array, f_array)
   real(c_double), intent(out), target :: c_array
   real(wp), intent(in) :: f_array
   if (c_associated(c_loc(c_array))) then
      c_array = f_array
   endif
end subroutine c_return_double_0d

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran value has been calculated
subroutine c_return_double_0d_alloc(c_array, f_array)
   real(c_double), intent(out), target :: c_array
   real(wp), allocatable, intent(in) :: f_array
   if (c_associated(c_loc(c_array)) .and. allocated(f_array)) then
      c_array = f_array
   endif
end subroutine c_return_double_0d_alloc

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran array has been populated
subroutine c_return_double_1d(c_array, f_array)
   real(c_double), intent(out), target :: c_array(:)
   real(wp), intent(in) :: f_array(:)
   if (c_associated(c_loc(c_array))) then
      c_array = f_array
   endif
end subroutine c_return_double_1d

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran array has been populated
subroutine c_return_double_1d_alloc(c_array, f_array)
   real(c_double), intent(out), target :: c_array(:)
   real(wp), allocatable, intent(in) :: f_array(:)
   if (c_associated(c_loc(c_array)) .and. allocated(f_array)) then
      c_array = f_array
   endif
end subroutine c_return_double_1d_alloc

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran array has been populated (2D version)
subroutine c_return_double_2d(c_array, f_array)
   real(c_double), intent(out), target :: c_array(:,:)
   real(wp), intent(in) :: f_array(:,:)
   if (c_associated(c_loc(c_array))) then
      c_array = f_array
   endif
end subroutine c_return_double_2d

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran array has been populated (2D version)
subroutine c_return_double_2d_alloc(c_array, f_array)
   real(c_double), intent(out), target :: c_array(:,:)
   real(wp), allocatable, intent(in) :: f_array(:,:)
   if (c_associated(c_loc(c_array)) .and. allocated(f_array)) then
      c_array = f_array
   endif
end subroutine c_return_double_2d_alloc

!> convert vector of chars to deferred size character
subroutine c_string_convert(f_string, c_string)
   character(c_char), dimension(*), intent(in) :: c_string
   character(len=:), allocatable, intent(out) :: f_string
   integer :: i
   i = 0
   f_string = ''
   do
      i = i+1
      if (c_string(i).eq.c_null_char) exit
      f_string = f_string//c_string(i)
   enddo
end subroutine c_string_convert

subroutine c_get_double_0d_alloc(f_array, c_array)
   real(c_double), intent(in), target :: c_array
   real(wp), allocatable, intent(out) :: f_array
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   endif
end subroutine c_get_double_0d_alloc

subroutine c_get_double_1d_alloc(f_array, c_array)
   real(c_double), intent(in), target :: c_array(:)
   real(wp), allocatable, intent(out) :: f_array(:)
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   endif
end subroutine c_get_double_1d_alloc

subroutine c_get_double_2d_alloc(f_array, c_array)
   real(c_double), intent(in), target :: c_array(:,:)
   real(wp), allocatable, intent(out) :: f_array(:,:)
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   endif
end subroutine c_get_double_2d_alloc

subroutine c_get_double_0d(f_array, c_array, default)
   real(c_double), intent(in), target :: c_array
   real(wp), intent(out) :: f_array
   real(wp), intent(in) :: default
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   else
      f_array = default
   endif
end subroutine c_get_double_0d

subroutine c_get_double_1d(f_array, c_array, default)
   real(c_double), intent(in), target :: c_array(:)
   real(wp), intent(out) :: f_array(:)
   real(wp), intent(in) :: default
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   else
      f_array = default
   endif
end subroutine c_get_double_1d

subroutine c_get_double_2d(f_array, c_array, default)
   real(c_double), intent(in), target :: c_array(:,:)
   real(wp), intent(out) :: f_array(:,:)
   real(wp), intent(in) :: default
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   else
      f_array = default
   endif
end subroutine c_get_double_2d

subroutine c_get_int_0d(f_array, c_array, default)
   integer(c_int), intent(in), target :: c_array
   integer, intent(out) :: f_array
   integer, intent(in) :: default
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   else
      f_array = default
   endif
end subroutine c_get_int_0d

subroutine c_get_int_1d(f_array, c_array, default)
   integer(c_int), intent(in), target :: c_array(:)
   integer, intent(out) :: f_array(:)
   integer, intent(in) :: default
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   else
      f_array = default
   endif
end subroutine c_get_int_1d

subroutine c_get_int_2d(f_array, c_array, default)
   integer(c_int), intent(in), target :: c_array(:,:)
   integer, intent(out) :: f_array(:,:)
   integer, intent(in) :: default
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   else
      f_array = default
   endif
end subroutine c_get_int_2d

! open a unit for IO
subroutine c_open_file(unit, c_output, prlevel)
   use iso_fortran_env, only: output_unit
   integer, intent(out) :: unit
   character(kind=c_char), intent(in) :: c_output(*)
   integer, intent(in) :: prlevel
   character(len=:), allocatable :: output
   logical :: exist

   call c_string_convert(output, c_output)

   unit = output_unit
   if (output /= '-' .and. prlevel > 0) then
      inquire(file=output, exist=exist)
      if (exist) then
         open(file=output, newunit=unit)
      end if
   endif
end subroutine

end module xtb_api_utils
