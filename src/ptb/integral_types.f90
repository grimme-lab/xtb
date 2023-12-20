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

!> Declaration of an integral storage container to collect all overlap related integrals
module xtb_ptb_integral_types
   use mctc_env, only: wp
   implicit none
   private

   public :: new_aux_integral

   !> Integral container to store all overlap related integrals
   type, public :: aux_integral_type
      !> Overlap integrals
      real(wp), allocatable :: overlap_h0_1(:, :), overlap_h0_2(:, :), overlap_xc(:, :)
      !> Overlap^(x) and Overlap^(1-x) integrals
      real(wp), allocatable :: overlap_to_x(:, :), overlap_to_1_x(:, :)
      !> Norm of overlap matrix (Normalization factors)
      real(wp), allocatable :: norm(:)
   end type aux_integral_type

contains

!> Create and allocate a new integral container storage
   subroutine new_aux_integral(self, nao)
      !> Instance of the integral container
      type(aux_integral_type), intent(out) :: self
      !> Dimension of the integrals
      integer, intent(in) :: nao

      allocate (self%norm(nao), source=0.0_wp)
      allocate (self%overlap_h0_1(nao, nao), self%overlap_h0_2(nao, nao), source=0.0_wp)
      allocate (self%overlap_xc(nao, nao), source=0.0_wp)
      allocate (self%overlap_to_x(nao, nao), source=0.0_wp)
      allocate (self%overlap_to_1_x(nao, nao), source=0.0_wp)
   end subroutine new_aux_integral

end module xtb_ptb_integral_types
