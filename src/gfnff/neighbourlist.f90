! This file is part of xtb.
!
! Copyright (C) 2021 Sebastian Ehlert
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

!> Neighbourlists for the GFN-FF
module xtb_gfnff_neighbourlist
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: TGFFNeighbourList, new


   type :: TGFFNeighbourList
      logical :: initialized = .false.
      integer :: nhb1
      integer :: nhb2
      integer :: nxb
      !> atomic charges (obtained from EEQ)
      real(wp), allocatable :: q(:)
      !> atom xyz, used to check for HB list update
      real(wp), allocatable :: hbrefgeo(:, :)
      !> HBs loose
      integer, allocatable :: hblist1(:, :)
      !> HBs bonded
      integer, allocatable :: hblist2(:, :)
      !> XBs
      integer, allocatable :: hblist3(:, :)
   end type TGFFNeighbourList

   interface new
      module procedure :: newGFFNeighbourList
   end interface


contains


subroutine newGFFNeighbourList(self, n, nhb1, nhb2, nxb)
   type(TGFFNeighbourList), intent(out) :: self
   integer, intent(in) :: n
   integer, intent(in) :: nhb1
   integer, intent(in) :: nhb2
   integer, intent(in) :: nxb
   self%initialized = .true.
   self%nhb1 = nhb1
   self%nhb2 = nhb2
   self%nxb = nxb
   allocate(self%q(n), source=0.0_wp)
   allocate(self%hbrefgeo(3, n), source=0.0_wp)
   allocate(self%hblist1(3, self%nhb1), source=0)
   allocate(self%hblist2(3, self%nhb2), source=0)
   allocate(self%hblist3(3, self%nxb), source=0)
end subroutine newGFFNeighbourList


end module xtb_gfnff_neighbourlist
