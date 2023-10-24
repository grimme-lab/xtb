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

!> Two-step SCF of the PTB method

module xtb_ptb_scf
   use mctc_io, only: structure_type, new
   use mctc_env, only: wp

   use tblite_basis_type, only: basis_type, get_cutoff
   use tblite_integral_dipole, only: get_dipole_integrals
   use tblite_cutoff, only: get_lattice_points

   implicit none
   private

   public :: twostepscf

contains

   subroutine twostepscf(mol, bas)
      !> Molecular structure data
      type(structure_type), intent(inout) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      real(wp), allocatable :: lattr(:, :), overlap(:, :)
      real(wp), allocatable :: dipole(:, :, :)
      real(wp) :: cutoff

      real(wp) :: norm(bas%nao), tmp1

      integer :: i,j,ij

      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

      allocate (overlap(bas%nao, bas%nao))
      allocate (dipole(3, bas%nao, bas%nao))
      call get_dipole_integrals(mol, lattr, cutoff, bas, overlap, dipole)

      ij = 0
      do i=1,bas%nao
         norm(i)=1./sqrt(overlap(i,i))
      enddo

      ij = 0
      do i=1,bas%nao
         do j=1,i
            tmp1=norm(i)*norm(j)
            overlap(i,j)=overlap(i,j)*tmp1 
            overlap(j,i)=overlap(i,j)
         enddo
      enddo
      
      write(*,*) "Overlap:"
      do i = 1, bas%nao
         do j = 1, bas%nao
            write(*,'(f10.6)', advance="no") overlap(i, j)
         end do
         write(*,*) ""
      end do


      write(*,*) "Dipole:"
      do i = 1, bas%nao
         do j = 1, bas%nao
            write(*,'(f12.6)', advance="no") dipole(1, i, j)
         end do
         write(*,*) ""
      end do

      stop

   end subroutine twostepscf

end module xtb_ptb_scf

