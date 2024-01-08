! This file is part of xtb.
!
! Copyright (C) 2023 Marcel Mueller
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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

!> Property printout of the PTB method

module xtb_ptb_property
#if WITH_TBLITE
   use mctc_env, only: wp
   use mctc_io, only: structure_type

   use tblite_basis_type, only: basis_type

   implicit none
   private

   public :: print_charges_to_screen

contains

   subroutine print_charges_to_screen(iunit, mol, bas, qat, qsh)
      implicit none
      integer, intent(in) :: iunit
      !> mctc structure type
      type(structure_type), intent(in) :: mol
      type(basis_type), intent(in) :: bas
      real(wp), intent(in) :: qat(:, :) ! Atomic partial charges
      real(wp), intent(in) :: qsh(:, :) ! Shell populations
      integer :: i, iid, ish, is
      real(wp) :: sumq

      write (iunit, '(a)')
      write (iunit, '(4x,"* Atomic partial charges (q)     Shell populations (p)")')
      write (iunit, '(4x)', advance="no")
      do i = 1,29+11*maxval(bas%nsh_at(:))
         write (iunit, '(a)', advance="no") "-"
      end do
      write (iunit, '(/)', advance="no")
      write (iunit, '(6x,"#  sym   q",11x)', advance="no")
      do i = 1, maxval(bas%nsh_at(:))
         write (iunit, '(i11)', advance="no") i
      end do
      write (iunit, '(/)', advance="no")
      write (iunit, '(4x)', advance="no")
      do i = 1,29+11*maxval(bas%nsh_at(:))
         write (iunit, '(a)', advance="no") "-"
      end do
      write (iunit, '(/)', advance="no")
      do i = 1, mol%nat
         iid = mol%id(i)
         is = bas%ish_at(i)
         write (iunit, '(i7,2x,a4,f9.5)', advance="no") &
            & i, mol%sym(iid), qat(i, 1)
         write (iunit, '(11x)', advance="no")
         do ish = 1, bas%nsh_at(i)
            write (iunit, '(2x,f9.5)', advance="no") qsh(is + ish, 1)
         end do
         write (iunit, '(/)', advance="no")
      end do
      write (iunit, '(4x)', advance="no")
      do i = 1,29+11*maxval(bas%nsh_at(:))
         write (iunit, '(a)', advance="no") "-"
      end do
      write (iunit, '(/)', advance="no")
      sumq = sum(qat(:, 1))
      if (abs(sumq) < 1.0E-7_wp) then
         sumq = 0.0_wp
      end if
      write (iunit, '(6x,a,f10.5)') &
         & "total:", sumq 

   end subroutine print_charges_to_screen

#endif
end module xtb_ptb_property
