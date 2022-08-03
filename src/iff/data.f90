! This file is part of xtb.
!
! Copyright (C) 2022 Christoph Plett
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

!> Topological data for force field type calculations
module xtb_iff_data
   use xtb_mctc_accuracy, only: wp
   implicit none
   private

   public :: TIFFData

   type :: TIFFData

      integer :: n !Total number atoms
      integer :: n1 !Number atoms molA
      integer :: n2 !Number atoms molB
      integer :: nlmo1 !Number lmos
      integer :: nlmo2
      real(wp), allocatable :: xyz1(:, :) !Coordinates
      real(wp), allocatable :: rlmo1(:, :) !LMO position
      real(wp), allocatable :: q1(:) !Charge
      real(wp), allocatable :: qdr1(:)
      real(wp), allocatable ::xyzdr1(:, :)
      real(wp), allocatable :: cn1(:)
      real(wp), allocatable :: z1(:)
      real(wp), allocatable :: alp1(:)
      real(wp), allocatable :: qct1(:, :)
      integer, allocatable :: at1(:)
      integer, allocatable :: lmo1(:)
      real(wp), allocatable :: xyz2(:, :)
      real(wp), allocatable :: rlmo2(:, :)
      real(wp), allocatable :: q2(:)
      real(wp), allocatable :: qdr2(:)
      real(wp), allocatable ::xyzdr2(:, :)
      real(wp), allocatable :: cn2(:)
      real(wp), allocatable :: z2(:)
      real(wp), allocatable :: alp2(:)
      real(wp), allocatable :: qct2(:, :)
      integer, allocatable :: at2(:)
      integer, allocatable :: lmo2(:)
      real(wp), allocatable :: c6ab(:, :)
      real(wp), allocatable :: alpab(:, :)
      real(wp), allocatable :: cprob(:)
      real(wp), allocatable :: xyz(:, :)
      real(wp), allocatable :: q(:)
      real(wp), allocatable :: cn(:)
      real(wp), allocatable :: alp0(:)
      real(wp), allocatable :: gab1(:, :)
      real(wp), allocatable :: gab2(:, :)
      real(wp), allocatable :: den1(:, :, :)
      real(wp), allocatable :: den2(:, :, :)
      real(wp), allocatable :: qcm1(:)
      real(wp), allocatable :: qcm2(:)
      integer, allocatable :: at(:)
      integer, allocatable :: neigh(:, :)

   contains

      procedure :: zero
      procedure :: allocateIFFData

   end type TIFFData

contains

   subroutine zero(self)
      class(TIFFData), intent(out) :: self

      if (allocated(self%xyz1)) deallocate (self%xyz1)
      if (allocated(self%rlmo1)) deallocate (self%rlmo1)
      if (allocated(self%q1)) deallocate (self%q1)
      if (allocated(self%qdr1)) deallocate (self%qdr1)
      if (allocated(self%xyzdr1)) deallocate (self%xyzdr1)
      if (allocated(self%cn1)) deallocate (self%cn1)
      if (allocated(self%z1)) deallocate (self%z1)
      if (allocated(self%alp1)) deallocate (self%alp1)
      if (allocated(self%qct1)) deallocate (self%qct1)
      if (allocated(self%at1)) deallocate (self%at1)
      if (allocated(self%lmo1)) deallocate (self%lmo1)
      if (allocated(self%xyz2)) deallocate (self%xyz2)
      if (allocated(self%rlmo2)) deallocate (self%rlmo2)
      if (allocated(self%q2)) deallocate (self%q2)
      if (allocated(self%qdr2)) deallocate (self%qdr2)
      if (allocated(self%xyzdr2)) deallocate (self%xyzdr2)
      if (allocated(self%cn2)) deallocate (self%cn2)
      if (allocated(self%z2)) deallocate (self%z2)
      if (allocated(self%alp2)) deallocate (self%alp2)
      if (allocated(self%qct2)) deallocate (self%qct2)
      if (allocated(self%at2)) deallocate (self%at2)
      if (allocated(self%lmo2)) deallocate (self%lmo2)
      if (allocated(self%c6ab)) deallocate (self%c6ab)
      if (allocated(self%alpab)) deallocate (self%alpab)
      if (allocated(self%cprob)) deallocate (self%cprob)
      if (allocated(self%xyz)) deallocate (self%xyz)
      if (allocated(self%q)) deallocate (self%q)
      if (allocated(self%cn)) deallocate (self%cn)
      if (allocated(self%alp0)) deallocate (self%alp0)
      if (allocated(self%gab1)) deallocate (self%gab1)
      if (allocated(self%gab2)) deallocate (self%gab2)
      if (allocated(self%den1)) deallocate (self%den1)
      if (allocated(self%den2)) deallocate (self%den2)
      if (allocated(self%qcm1)) deallocate (self%qcm1)
      if (allocated(self%qcm2)) deallocate (self%qcm2)
      if (allocated(self%at)) deallocate (self%at)
      if (allocated(self%neigh)) deallocate (self%neigh)

   end subroutine zero

   subroutine allocateIFFData(self, n1, n2)
      use xtb_type_molecule, only: TMolecule

      class(TIFFData), intent(out) :: self

      integer, intent(in) :: n1, n2

      integer :: n

      self%n1 = n1
      self%n2 = n2

      !Allocate Infos of molA
      allocate(self%at1(n1),self%xyz1(3,n1),self%rlmo1(4,10*n1),self%q1(n1),self%lmo1(10*n1),&
      &self%cn1(n1),self%alp1(n1),self%qct1(n1,2),self%qdr1(n1),self%xyzdr1(3,n1),&
      &self%z1(n1), self%den1(2, 4, n1), self%gab1(n1, n1), self%qcm1(n1))
      self%rlmo1 = 0.0_wp

      !Allocate Infos of molB
      allocate(self%at2(n2),self%xyz2(3,n2),self%rlmo2(4,10*n2),self%q2(n2),self%lmo2(10*n2),&
      & self%cn2(n2),self%alp2(n2),self%qct2(n2,2),self%qdr2(n2),self%xyzdr2(3,n2),&
      & self%z2(n2), self%den2(2, 4, n2), self%gab2(n2, n2), self%qcm2(n2))
      self%rlmo2 = 0.0_wp

      !Allocate combined Infos
      self%n = n1 + n2
      n = self%n
      allocate(self%at(n),self%xyz(3,n),self%q(n),self%c6ab(n,n),self%alp0(n),self%cn(n),self%neigh(0:n,n),&
      &        self%alpab(n2, n1), self%cprob(n1))

   end subroutine allocateIFFData

end module xtb_iff_data
