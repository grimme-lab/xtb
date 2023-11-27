! This file is part of ancopt.
! SPDX-Identifier: LGPL-3.0-or-later
!
! ancopt is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! ancopt is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with ancopt.  If not, see <https://www.gnu.org/licenses/>.

module xtb_pbc_optimizer_lbfgs
   use xtb_mctc_accuracy, only : wp!, error_type, fatal_error
   use xtb_type_environment, only : TEnvironment
                                      !use ancopt_blas, only : dot
   use xtb_mctc_blas, only : mctc_dot !@thoma important err? same as dot above?


   use  xtb_pbc_optimizer_filter_type, only : filter_type
   implicit none
   private

   public :: optimizer_type, lbfgs_optimizer, lbfgs_input, new_lbfgs_optimizer

   type :: lbfgs_input
      integer :: memory = 20
   end type lbfgs_input
                                                   
  type, abstract :: optimizer_type                  
  contains                                          
  !   procedure :: relax_pbc                        
     procedure(step_i), deferred :: step
  end type optimizer_type                           
                                                    
  abstract interface                                
     subroutine step_i(self, val, gcurr, glast, displ, tot_step)
        import :: optimizer_type, wp                
        class(optimizer_type), intent(inout) :: self
        real(wp), intent(in) :: val                 
        real(wp), intent(in) :: gcurr(:)            
        real(wp), intent(in) :: glast(:)            
        integer, intent(in) :: tot_step
        real(wp), intent(out) :: displ(:)     
     end subroutine step_i                           
  end interface                                     


   type, extends(optimizer_type) :: lbfgs_optimizer
      integer :: iter
      integer :: nvar
      integer :: memory
      real(wp), allocatable :: hdiag(:)
      real(wp), allocatable :: s(:, :)
      real(wp), allocatable :: y(:, :)
      real(wp), allocatable :: rho(:)
   contains
      procedure :: step
   end type lbfgs_optimizer

contains


!> Create new limited memory BFGS optimization driver
subroutine new_lbfgs_optimizer(self, env, input, filter)
   !> Instance of the optimizer
   type(lbfgs_optimizer), intent(out) :: self
  !> Calculation environment
  type(TEnvironment), intent(inout) :: env
   !> Input for the LBFGS optimizer
   type(lbfgs_input), intent(in) :: input
   !> Transformation filter
   class(filter_type), intent(in) :: filter

   self%iter = 0
   self%nvar = filter%get_dimension()
   self%memory = input%memory
   allocate(self%s(self%nvar, input%memory), source=0.0_wp)
   allocate(self%y(self%nvar, input%memory), source=0.0_wp)
   allocate(self%rho(input%memory), source=0.0_wp)
   allocate(self%hdiag(self%nvar), source=1.0_wp)

   !select type(filter)
   !class is (anc_filter)
   !   self%hdiag(:) = filter%hdiag
   !end select
end subroutine new_lbfgs_optimizer


!> updates displacement using the formula given by Nocedal, generally known
!> as limited memory BFGS algorithm
subroutine lbfgs_step(iter, memory, nvar, gradient, glast, displacement, hdiag, s, y, rho, tot_step)
   !> Current iteration step
   integer, intent(in) :: iter
   !> Memory limit for the LBFGS update
   integer, intent(in) :: memory
   !> 3*natoms
   integer, intent(in) :: nvar
   !> current gradient
   real(wp), intent(in) :: gradient(:)
   !> gradient on the last point
   real(wp), intent(in) :: glast(:)
   !> on input displacement from the last step, on exit new displacement vector
   real(wp), intent(inout) :: displacement(:)
   !> *inverse* Hessian in diagonal form
   real(wp), intent(in) :: hdiag(:)
   integer, intent(in) :: tot_step

   !> LBFGS scratch array of former displacements
   real(wp), intent(inout) :: s(:, :)
   !> LBFGS scratch array of former gradient changes
   real(wp), intent(inout) :: y(:, :)
   !> LBFGS scratch array of dot products between s and y
   real(wp), intent(inout) :: rho(:)

   real(wp), allocatable :: d(:), q(:), a(:)
   real(wp) :: b

   integer :: thisiter, lastiter, mem
   integer :: i
   real(wp) :: f_damp

   allocate(q(nvar), d(nvar), a(memory), source = 0.0_wp)

   thisiter = mod(iter-1, memory) + 1
   s(:, thisiter) = displacement
   y(:, thisiter) = gradient - glast

   b = mctc_dot(s(:, thisiter), y(:, thisiter))

   rho(thisiter) = 1.0_wp / max(abs(b), epsilon(1.0_wp))

   q = gradient
   do mem = iter, max(1, iter-memory+1), -1
      i = mod(mem-1, memory)+1
      a(i) = rho(i) * mctc_dot(s(:, i), q)
      q = q - a(i) * y(:, i)
   enddo

   d = hdiag * q
if(iter.le.memory) then
   do mem = max(1, iter-memory), iter
      i = mod(mem-1, memory)+1
      b = rho(i) * mctc_dot(y(:, i), d)
      d = d + s(:, i) * (a(i) - b)
   enddo
else
   do mem = max(1, iter-memory), iter-1
      i = mod(mem-1, memory)+1
      b = rho(i) * mctc_dot(y(:, i), d)
      d = d + s(:, i) * (a(i) - b)
   enddo
endif
   f_damp = 1.0_wp/(1.0_wp + 3000.0_wp*real(tot_step)**(-3))
   if (maxval(abs(d)).lt.0.1_wp.or.tot_step.gt.50) then
     displacement = -d
   else
     displacement = -d*f_damp
   endif
end subroutine lbfgs_step


subroutine step(self, val, gcurr, glast, displ, tot_step)
   class(lbfgs_optimizer), intent(inout) :: self
   real(wp), intent(in) :: val
   real(wp), intent(in) :: gcurr(:)
   real(wp), intent(in) :: glast(:)
   integer, intent(in) :: tot_step
   real(wp), intent(out) :: displ(:)

   self%iter = self%iter + 1

   call lbfgs_step(self%iter, self%memory, self%nvar, gcurr, glast, displ, self%hdiag, &
      & self%s, self%y, self%rho, tot_step)
end subroutine step

end module xtb_pbc_optimizer_lbfgs
