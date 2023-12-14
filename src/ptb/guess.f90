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

!> Setup of the effective Hamiltonians for both iterations in PTB

module xtb_ptb_guess
   use mctc_io, only: structure_type
   use mctc_env, only: wp, error_type

   use tblite_wavefunction, only: wavefunction_type
   use tblite_basis_type, only: basis_type

   use xtb_ptb_vdzp, only: nshell
   implicit none

   private

   public :: guess_qsh, get_psh_from_qat, get_psh_from_qsh

contains

   subroutine guess_qsh(wfn, bas)
      !> Wavefunction of tblite type
      type(wavefunction_type), intent(inout) :: wfn
      !> Basis set and core-valence basis set data
      type(basis_type), intent(in) :: bas

      integer :: ish, iat, ii, spin

      do spin = 1, size(wfn%qat, 2)
         do iat = 1, size(wfn%qat, 1)
            ii = bas%ish_at(iat)
            do ish = 1, bas%nsh_at(iat)
               !> wfn%qat on input are actually the atomic CHARGES
               !### POP version
               ! wfn%qsh(ii + ish, spin) = wfn%n0sh(ii + ish) * &
               ! & ( wfn%n0at(iat) - wfn%qat(iat, spin) ) / wfn%n0at(iat)
               !> wfn%qsh on output are the shell POPULATIONS
               !### CHRG version
               wfn%qsh(ii+ish, spin) = (wfn%n0sh(ii+ish) / wfn%n0at(iat)) * wfn%qat(iat, spin)
               !> wfn%qsh on output are the shell CHARGES
               !##### DEV WRITE #####
               ! write(*,*) ii + ish, wfn%qsh(ii + ish, spin)
               !#####################
            end do
         end do
      end do
   end subroutine guess_qsh

   function get_psh_from_qat(wfn, bas) result(psh)
      !> Wavefunction of tblite type
      type(wavefunction_type), intent(in) :: wfn
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> Shell populations
      real(wp) :: psh(bas%nsh,wfn%nspin)

      integer :: ish, iat, ii, spin

      psh = 0.0_wp
      do spin = 1, size(wfn%qat, 2)
         do iat = 1, size(wfn%qat, 1)
            ii = bas%ish_at(iat)
            do ish = 1, bas%nsh_at(iat)
               !> wfn%qat on input are actually the atomic CHARGES
               psh(ii + ish, spin) = wfn%n0sh(ii + ish) * &
               & ( wfn%n0at(iat) - wfn%qat(iat, spin) ) / wfn%n0at(iat)
            end do
         end do
      end do
   end function get_psh_from_qat

   function get_psh_from_qsh(wfn, bas) result(psh)
      !> Wavefunction of tblite type
      type(wavefunction_type), intent(in) :: wfn
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> Shell populations
      real(wp) :: psh(bas%nsh,wfn%nspin)

      integer :: ish, iat, ii, spin

      psh = 0.0_wp
      do spin = 1, size(wfn%qat, 2)
         do iat = 1, size(wfn%qat, 1)
            ii = bas%ish_at(iat)
            do ish = 1, bas%nsh_at(iat)
               psh(ii + ish, spin) = wfn%n0sh(ii + ish) - wfn%qsh(ii + ish, spin)
            end do
         end do
      end do
   end function get_psh_from_qsh
end module xtb_ptb_guess
