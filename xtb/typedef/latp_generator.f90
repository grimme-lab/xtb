! This file is part of xtb.
! Based on the DFTB+ implementation: https://github.com/dftbplus/dftbplus
!
! Copyright (C) 2006-2019 DFTB+ developers group
! Copyright (C) 2017-2019 Stefan Grimme
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

!> Implements an interator over lattice points up to a certain cutoff
module tbdef_latp_generator
   use iso_fortran_env, only: wp => real64
   implicit none
   private

   public :: tb_latp_generator


   !> Lattice point iterator
   type :: tb_latp_generator
      private
      !> Square of the cutoff
      real(wp) :: cutoff2
      !> Lattice vectors
      real(wp) :: lattice(3, 3)
      !> Index ranges necessary to obtain all lattice points within given cutoff
      integer :: ranges(2, 3)
      !> Indices of the current lattice point
      integer :: curPoint(3)
      !> Increase upper limit of necessary index ranges by this amount
      integer :: pos_ext = 0
      !> Decrease lower limit of necessary index ranges by this amount
      integer :: neg_ext = 0
      !> Return all lattice vector in the index ranges, not just those shorter
      !  than cutoff
      logical :: only_inside = .false.
      !> Whether lattice vector (0, 0, 0) should be excluded when iterationg
      logical :: exclude_orig = .true.
      !> Whether to exclude points connected by inversion
      logical :: exclude_inv = .false.
      !> Whether iterator has returned last point
      logical :: finished = .false.
   contains
      procedure :: new => latpgen_new
      procedure :: get_next_point => latpgen_get_next_point
      procedure :: get_all_points => latpgen_get_all_points
   end type tb_latp_generator

contains


!> Initialises the lattice point iterator.
subroutine latpgen_new(self, lattice, inv_lat, cutoff, neg_extension, &
      &                pos_extension, only_inside, reduce_by_inversion, &
      &                exclude_origin)
   !> Instance
   class(tb_latp_generator), intent(out) :: self
   !> Lattice vectors
   real(wp), intent(in) :: lattice(:, :)
   !> Inverse lattice vectors
   real(wp), intent(in) :: inv_lat(:, :)
   !> Cutoff radius for lattice points
   real(wp), intent(in) :: cutoff
   !> Same as posExtension for negative lattice vectors
   integer, intent(in), optional :: neg_extension
   !> Extend the set along the positive lattice vectors with that many
   !  additional lattice vectors.
   integer, intent(in), optional :: pos_extension
   !> Include only lattice points being inside the cutoff radius
   logical, intent(in), optional :: only_inside
   !> whether to include time reversal symmetry when generating k-points
   logical,  intent(in), optional :: reduce_by_inversion
   !> whether to exclude the (0,0,0) point
   logical,  intent(in), optional :: exclude_origin

   integer :: ranges(2, 3)

   if (present(neg_extension)) self%neg_ext = neg_extension
   if (present(pos_extension)) self%pos_ext = pos_extension
   if (present(only_inside)) self%only_inside = only_inside
   if (present(reduce_by_inversion)) self%exclude_inv = reduce_by_inversion
   if (present(exclude_origin)) self%exclude_orig = exclude_origin

   self%lattice = lattice
   call get_ranges(cutoff, inv_lat, self%pos_ext, self%neg_ext, ranges)
   self%ranges(1, :) = minval(ranges, dim=1)
   self%ranges(2, :) = maxval(ranges, dim=1)
   self%curPoint = self%ranges(1, :)
   self%cutoff2 = cutoff**2

end subroutine latpgen_new


!> Delivers the next lattice point
subroutine latpgen_get_next_point(self, latticePoint, finished)
   !> Instance
   class(tb_latp_generator), intent(inout) :: self
   !> Next lattice point
   real(wp), intent(out) :: latticePoint(:)
   !> Whether the returned point was the last one.
   logical, intent(out) :: finished

   real(wp) :: rr(3)
   integer :: curPoint(3)

   curPoint = self%curPoint
   do
      finished = curPoint(1) > self%ranges(2, 1)
      if (finished) exit
      if (self%exclude_inv) then
         if (curPoint(1) < 0) then
            curPoint(1) = 0
         end if
         if (curPoint(2) < 0 .and. curPoint(1) == 0) then
            curPoint(2) = 0
         end if
         if (curPoint(3) < 0 .and. curPoint(1) == 0 .and. curPoint(2) == 0) then
            curPoint(3) = 0
         end if
      end if
      if (self%exclude_orig .and. all(curPoint == 0)) then
         call get_next_point(self%ranges, curPoint)
         cycle
      end if
      latticePoint = real(curPoint, wp)
      if (self%only_inside) then
         rr = latticePoint(1) * self%lattice(:,1) + &
            & latticePoint(2) * self%lattice(:,2) + &
            & latticePoint(3) * self%lattice(:,3)
         if (sum(rr**2) <= self%cutoff2) exit
         call get_next_point(self%ranges, curPoint)
      else
         exit
      end if
   end do

   if (.not. finished) then
      ! generate image for next call
      call get_next_point(self%ranges, curPoint)
   end if
   self%curPoint = curPoint

end subroutine latpgen_get_next_point


!> Returns all lattice points within a cutoff at the same time
subroutine latpgen_get_all_points(self, latticePoints)
   !> Instance
   class(tb_latp_generator), intent(inout) :: self
   !> Lattice points
   real(wp), allocatable, intent(out) :: latticePoints(:, :)

   real(wp), allocatable :: tmpLatPoints(:, :)
   integer :: maxLatPoints, iLatPoint
   logical :: finished

   maxLatPoints = product(self%ranges(2, :) - self%ranges(1, :) + 1)
   allocate(tmpLatPoints(3, maxLatPoints + 1))
   do iLatPoint = 1, maxLatPoints + 1
      call self%get_next_point(tmpLatPoints(:,iLatPoint), finished)
      if (finished) exit
   end do
   latticePoints = tmpLatPoints(:, 1:iLatPoint - 1)

end subroutine latpgen_get_all_points


!> Helper function to increase a tuple of 3 indices by one
subroutine get_next_point(ranges, point)
   !> Lower and upper bounds for the lattice point indices. Shape: (2, 3)
   integer, intent(in) :: ranges(:, :)
   !> Current lattice point, next one on exit.
   integer, intent(inout) :: point(3)

   point(3) = point(3) + 1
   if (point(3) > ranges(2, 3)) then
      point(3) = ranges(1, 3)
      point(2) = point(2) + 1
      if (point(2) > ranges(2, 2)) then
         point(2) = ranges(1, 2)
         point(1) = point(1) + 1
      end if
   end if

end subroutine get_next_point


!> Calculate the range of images of the central cell that interact
subroutine get_ranges(dist, recVec2p, pos_ext, neg_ext, ranges)
   !> distance of interaction
   real(wp), intent(in) :: dist
   !> reciprocal lattice vector
   real(wp), intent(in) :: recVec2p(:, :)
   !> Extend the set along the positive lattice vectors with that many
   !  additional lattice vectors.
   integer, intent(in) :: pos_ext
   !> Same as posExtension for negative lattice vectors
   integer, intent(in) :: neg_ext
   !> Array of the two extremal points
   integer, intent(out) :: ranges(:, :)

   integer :: ii, iTmp

   !:ASSERT(dist >= 0.0_wp)

   do ii = 1, 3
      iTmp = floor(dist * sqrt(sum(recVec2p(:, ii)**2)))
      ranges(1, ii) = -(iTmp + neg_ext)
      ranges(2, ii) = iTmp + pos_ext
   end do

end subroutine get_ranges

end module tbdef_latp_generator
