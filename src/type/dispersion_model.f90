! This file is part of xtb.
!
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

module xtb_type_dispersionmodel
   use xtb_mctc_accuracy, only : wp
   implicit none
   public :: TDispersionModel, init
   private

   integer, parameter :: max_elem = 118

   type :: TDispersionModel
      real(wp) :: g_a = 0.0_wp
      real(wp) :: g_c = 0.0_wp
      integer, allocatable :: atoms(:)
      integer, allocatable :: nref(:)
      integer, allocatable :: ncount(:, :)
      real(wp), allocatable :: cn(:, :)
      real(wp), allocatable :: q(:, :)
      real(wp), allocatable :: alpha(:, :, :)
      real(wp), allocatable :: c6(:, :, :, :)
   end type TDispersionModel

   interface init
      module procedure :: initDispersionModel
   end interface init

contains

subroutine initDispersionModel(self, maxElem, maxRef, maxFreq)
   type(TDispersionModel), intent(out) :: self
   integer, intent(in), optional :: maxElem
   integer, intent(in), optional :: maxRef
   integer, intent(in), optional :: maxFreq
   integer :: elem, ref, freq
   if (present(maxElem)) then
      elem = maxElem
   else
      elem = 118
   end if
   if (present(maxRef)) then
      ref = maxRef
   else
      ref = 7
   end if
   if (present(maxFreq)) then
      freq = maxFreq
   else
      freq = 23
   end if
   allocate(self%atoms(elem), source=0)
   allocate(self%nref(elem), source=0)
   allocate(self%ncount(ref, elem), source=0)
   allocate(self%cn(ref, elem), source=0.0_wp)
   allocate(self%q(ref, elem), source=0.0_wp)
   allocate(self%alpha(freq, ref, elem), source=0.0_wp)
   allocate(self%c6(ref, ref, elem, elem), source=0.0_wp)
end subroutine initDispersionModel

end module xtb_type_dispersionmodel
