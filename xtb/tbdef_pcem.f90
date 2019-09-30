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

module tbdef_pcem
   use iso_fortran_env, wp => real64
   implicit none

   public :: tb_pcem
   private

   type :: tb_pcem
      integer  :: n = 0
      real(wp),allocatable :: xyz(:,:)
      real(wp),allocatable :: q(:)
      real(wp),allocatable :: gam(:)
      real(wp),allocatable :: grd(:,:)
   contains
   procedure :: allocate => allocate_pcem
   procedure :: deallocate => deallocate_pcem
   end type tb_pcem

contains

subroutine allocate_pcem(self,n)
   implicit none
   class(tb_pcem),intent(inout) :: self
   integer,intent(in) :: n
   call self%deallocate
   self%n = n
   allocate(self%xyz(3,n), source = 0.0_wp )
   allocate(self%q  (n),   source = 0.0_wp )
   allocate(self%gam(n),   source = 0.0_wp )
   allocate(self%grd(3,n), source = 0.0_wp )
end subroutine allocate_pcem

subroutine deallocate_pcem(self)
   implicit none
   class(tb_pcem),intent(inout) :: self
   self%n = 0
   if (allocated(self%xyz)) deallocate(self%xyz)
   if (allocated(self%q))   deallocate(self%q)
   if (allocated(self%gam)) deallocate(self%gam)
   if (allocated(self%grd)) deallocate(self%grd)
end subroutine deallocate_pcem

end module tbdef_pcem
