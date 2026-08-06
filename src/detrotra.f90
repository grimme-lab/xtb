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

module xtb_detrotra
   use xtb_mctc_accuracy, only : sp, wp
   implicit none
   private

   public :: detrotra4, detrotra8

contains

!> Determine rotational and translational modes using single precision modes
subroutine detrotra4(linear, n, xyz, h, eig)
   logical, intent(in) :: linear
   integer, intent(in) :: n
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Eigenvectors from the projected Lindh diagonalization
   real(sp), intent(in) :: h(3*n, 3*n)
   !> Eigenvalues from the projected Lindh diagonalization
   real(sp), intent(inout) :: eig(3*n)
   integer, allocatable :: rigid(:)

   call detrotra_worker(linear, n, xyz, h, real(eig, wp), rigid)
   ! Identifier for rotational and translational modes
   eig(rigid) = 0.0_sp
end subroutine detrotra4

!> Determine rotational and translational modes using double precision modes
subroutine detrotra8(linear, n, xyz, h, eig)
   logical, intent(in) :: linear
   integer, intent(in) :: n
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Eigenvectors from the projected Lindh diagonalization
   real(wp), intent(in) :: h(3*n, 3*n)
   !> Eigenvalues from the projected Lindh diagonalization
   real(wp), intent(inout) :: eig(3*n)
   integer, allocatable :: rigid(:)

   call detrotra_worker(linear, n, xyz, h, eig, rigid)
   ! Identifier for rotational and translational modes
   eig(rigid) = 0.0_wp
end subroutine detrotra8

!> Identify the modes which best preserve all interatomic distances
subroutine detrotra_worker(linear, n, xyz, h, eig, rigid)
   logical, intent(in) :: linear
   integer, intent(in) :: n
   real(wp), intent(in) :: xyz(3, n)
   class(*), intent(in) :: h(:, :)
   real(wp), intent(in) :: eig(3*n)
   integer, allocatable, intent(out) :: rigid(:)
   real(wp), parameter :: low_mode_threshold = 0.05_wp
   integer :: i, j, ii, nn, n3, nend
   integer, allocatable :: ind(:)
   real(wp), allocatable :: mode(:), tmpxyz(:, :), score(:)
   real(wp) :: a0, b0, c0
   logical :: check_all

   n3 = 3*n
   nend = min(merge(5, 6, linear), n3)
   allocate(rigid(nend))
   if (nend == 0) return

   ! Preserve the historic low-mode candidate set when it is large enough.
   ! If it is not, inspect all modes instead of indexing uninitialized entries.
   check_all = count(eig <= low_mode_threshold) < nend
   allocate(mode(n3), tmpxyz(3, n), score(n3), ind(n3))

   nn = 0
   do ii = 1, n3
      ! Check only low-lying modes unless the full fallback is required.
      if (.not.check_all .and. eig(ii) > low_mode_threshold) cycle

      select type(h)
      type is(real(sp))
         mode = real(h(:, ii), wp)
      type is(real(wp))
         mode = h(:, ii)
      class default
         error stop "Unsupported mode precision in detrotra"
      end select

      ! Distort the initial geometry along the ii-th mode.
      tmpxyz = xyz + reshape(mode, shape(tmpxyz))

      ! Compare all interatomic distances of the original and distorted
      ! geometries.
      c0 = 0.0_wp
      do i = 2, n
         do j = 1, i - 1
            a0 = sqrt(sum((xyz(:, i) - xyz(:, j))**2))
            b0 = sqrt(sum((tmpxyz(:, i) - tmpxyz(:, j))**2))
            ! Sum the squared distance differences.
            c0 = c0 + (a0 - b0)**2
         end do
      end do

      nn = nn + 1
      ! Weight the distance change by the Lindh eigenvalue.
      score(nn) = sqrt(c0/n)*abs(eig(ii))
      ind(nn) = ii
   end do

   ! Sort in ascending order.
   call qsort(score, 1, nn, ind)
   rigid = ind(:nend)
end subroutine detrotra_worker

end module xtb_detrotra
