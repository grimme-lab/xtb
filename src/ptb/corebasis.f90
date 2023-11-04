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

!> Core basis functions relevant for approximated effective core potential in PTB

module xtb_ptb_corebasis
   use mctc_env, only: wp
   use mctc_io, only: structure_type

   use xtb_ptb_param, only: max_elem, highest_elem, &  !> General parameters
      & max_core_shell, max_core_prim, cbas_nshell, &  !> PTB core basis parameters
      & cbas_pqn, cbas_sl_exp, cbas_angshell !> PTB core basis parameters

   use tblite_basis_type, only: cgto_type, new_basis, basis_type
   use tblite_basis_ortho, only : orthogonalize
   use tblite_basis_slater, only : slater_to_gauss

   implicit none
   private

   public :: add_PTBcore_basis

contains

   subroutine add_PTBcore_basis(mol, bas)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set type
      type(basis_type), intent(out) :: bas

      integer :: isp, izp, ish, stat, il
      integer, allocatable :: nsh_id(:)
      integer :: ang_idx(0:2), ortho(max_core_shell)
      type(cgto_type), allocatable :: cgto(:, :)

      integer :: j, k

      nsh_id = cbas_nshell(mol%num)
      allocate (cgto(maxval(nsh_id), mol%nid))
      do isp = 1, mol%nid
         ! ang_idx = 0
         ! ortho = 0
         izp = mol%num(isp)
         !##### DEV WRITE #####
         write(*,*) "number of shells: ", nsh_id(isp)
         !#####################
         do ish = 1, nsh_id(isp)
            il = cbas_angshell(ish, izp)

            ! if (ang_idx(il) > 0) then
            !    ortho(ish) = ang_idx(il)
            ! else
            !    ang_idx(il) = ish
            ! end if

            !##### DEV WRITE #####
            write(*,*) "number of primitives: ", max_core_prim
            write(*,*) "shell type: ", cbas_pqn(ish, izp), cbas_angshell(ish, izp), cbas_sl_exp(ish, izp)
            !#####################

            call slater_to_gauss(max_core_prim, cbas_pqn(ish, izp), il, &
            & cbas_sl_exp(ish, izp), cgto(ish, isp), .true., stat)

            !##### DEV WRITE #####
            write(*,*) "N_prim: ", cgto(ish,isp)%nprim
            do k = 1, cgto(ish,isp)%nprim
               write(*, *) cgto(ish,isp)%alpha(k)
            enddo
            !#####################
         end do

         ! do ish = 1, nsh_id(isp)
         !    if (ortho(ish) > 0) then
         !       call orthogonalize(cgto(ortho(ish), isp), cgto(ish, isp))
         !    end if
         ! end do
      end do

      write(*,*) "---------FINAL PTB core basis-----------"
      do isp = 1, mol%nid
         write(*,*) "Atom :", mol%num(isp)
         write(*,*) "Number of shells :", nsh_id(isp)
         do j = 1, nsh_id(isp)
            write(*,*) "N_prim: ", cgto(j,isp)%nprim
            write(*,*) cgto(j,isp)%alpha
         enddo
      end do

      call new_basis(bas, mol, nsh_id, cgto, 1.0_wp)

   end subroutine add_PTBcore_basis

end module xtb_ptb_corebasis

