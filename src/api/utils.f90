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

!> Utilities to work with data types from C
module xtb_api_utils
   use, intrinsic :: iso_c_binding
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule, only : TMolecule
   use xtb_type_environment, only : init
   implicit none
   private

   public :: c_f_character, f_c_character, verifyMolecule, checkGlobalEnv

contains


subroutine c_f_character(rhs, lhs)
   character(kind=c_char), intent(in) :: rhs(*)
   !> Resulting Fortran string
   character(len=:, kind=c_char), allocatable, intent(out) :: lhs

   integer :: ii

   do ii = 1, huge(ii) - 1
      if (rhs(ii) == c_null_char) then
         exit
      end if
   end do
   allocate(character(len=ii-1) :: lhs)
   lhs = transfer(rhs(1:ii-1), lhs)

end subroutine c_f_character


subroutine f_c_character(rhs, lhs, len)
   character(kind=c_char), intent(out) :: lhs(*)
   character(len=*), intent(in) :: rhs
   integer, intent(in) :: len
   integer :: length
   length = min(len-1, len_trim(rhs))

   lhs(1:length) = transfer(rhs(1:length), lhs(1:length)) // c_null_char

end subroutine f_c_character


!> Cold fusion check
integer function verifyMolecule(mol) result(status)
   type(TMolecule), intent(in) :: mol
   integer :: iat, jat
   status = 0
   do iat = 1, mol%n
      do jat = 1, iat - 1
         if (mol%dist(jat, iat) < 1.0e-9_wp) status = status + 1
      end do
   end do
end function verifyMolecule


subroutine checkGlobalEnv
   use xtb_mctc_global, only : persistentEnv
   if (.not.allocated(persistentEnv)) then
      allocate(persistentEnv)
      call init(persistentEnv)
   end if
end subroutine checkGlobalEnv


end module xtb_api_utils
