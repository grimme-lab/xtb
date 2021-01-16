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

!> Proxy module to use the ELPA library, provides stubs when disabled
module xtb_mctc_elpa
   use xtb_mctc_accuracy, only : sp, dp
   use xtb_type_environment, only : TEnvironment
#ifdef WITH_ELPA
   use elpa, only : elpa_t, elpa_init, elpa_allocate, elpa_deallocate, elpa_ok
#endif
!$ use omp_lib, only : omp_get_num_threads
   implicit none
   private
   public :: TELPASolver, initELPASolver

   type :: TELPASolver
      private
#ifdef WITH_ELPA
      class(elpa_t), pointer :: elpa => null()
#endif
      integer :: n
      real(sp), allocatable :: sbmat(:, :)
      real(dp), allocatable :: dbmat(:, :)
      real(sp), allocatable :: sevec(:, :)
      real(dp), allocatable :: devec(:, :)
   contains
      generic :: solve => sgen_solve, dgen_solve
      procedure :: sgen_solve
      procedure :: dgen_solve
      final :: finalizeELPA
   end type TELPASolver

   logical, save :: elpa_initialized = .false.

contains

subroutine initELPASolver(self, env, bmat)
   character(len=*), parameter :: source = "mctc_elpa_newELPA"
   type(TELPASolver), intent(out) :: self
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(in) :: bmat(:, :)
   integer :: stat

#ifdef WITH_ELPA
   if (.not.elpa_initialized) then
      if (elpa_init(20200417) /= elpa_ok) then
         call env%error("Cannot initialize ELPA library", source)
         return
      end if
      elpa_initialized = .true.
   end if

   self%elpa => elpa_allocate(stat)
   if (stat /= elpa_ok) then
      call env%error("Cannot allocate new ELPA object", source)
      return
   end if

   self%n = size(bmat, 1)
   allocate(self%dbmat(self%n, self%n))
   allocate(self%devec(self%n, self%n))

   call self%elpa%set("na", self%n, stat)
   if (stat /= elpa_ok) then
      call env%error("ELPA set na failed", source)
      return
   end if
   call self%elpa%set("nev", self%n, stat)
   if (stat /= elpa_ok) then
      call env%error("ELPA set nev failed", source)
      return
   end if
   call self%elpa%set("local_nrows", self%n, stat)
   if (stat /= elpa_ok) then
      call env%error("ELPA set local_nrows failed", source)
      return
   end if
   call self%elpa%set("local_ncols", self%n, stat)
   if (stat /= elpa_ok) then
      call env%error("ELPA set local_ncols failed", source)
      return
   end if
   call self%elpa%set("nblk", self%n, stat)
   if (stat /= elpa_ok) then
      call env%error("ELPA set nblk failed", source)
      return
   end if
!$ call self%elpa%set("omp_threads", omp_get_num_threads(), stat)
   if (stat /= elpa_ok) then
      call env%error("ELPA set omp_threads failed", source)
      return
   end if

   stat = self%elpa%setup()
   if (stat /= elpa_ok) then
      call env%error("ELPA setup failed", source)
      return
   end if
#endif

end subroutine initELPASolver

subroutine dgen_solve(self, env, amat, bmat, eval)
   character(len=*), parameter :: source = "mctc_elpa_solve"
   class(TELPASolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(out) :: eval(:)
   integer :: stat

   self%dbmat(:, :) = bmat
#ifdef WITH_ELPA
   call self%elpa%generalized_eigenvectors(amat, self%dbmat, eval, self%devec, &
      & .false., stat)
   if (stat /= elpa_ok) then
      call env%error("ELPA eigenvalue solver failed", source)
      return
   end if
#endif
   amat(:, :) = self%devec

end subroutine dgen_solve

subroutine sgen_solve(self, env, amat, bmat, eval)
   character(len=*), parameter :: source = "mctc_elpa_solve"
   class(TELPASolver), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(out) :: eval(:)
   integer :: stat

   self%sbmat(:, :) = bmat
#ifdef WITH_ELPA
   call self%elpa%generalized_eigenvectors(amat, self%sbmat, eval, self%sevec, &
      & .false., stat)
   if (stat /= elpa_ok) then
   end if
#endif
   amat(:, :) = self%sevec

end subroutine sgen_solve

subroutine finalizeELPA(self)
   type(TELPASolver) :: self

#ifdef WITH_ELPA
   if (elpa_initialized) then
      if (associated(self%elpa)) call elpa_deallocate(self%elpa)
   end if
#endif
end subroutine finalizeELPA

end module xtb_mctc_elpa
