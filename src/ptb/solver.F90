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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

!> Solver for the PTB method based on the diagonalization solver from tblite
module xtb_ptb_solver
   
#if WITH_TBLITE

   use mctc_env, only : sp, dp, wp, error_type
   use tblite_blas, only : gemm
   use tblite_lapack_sygvd, only : sygvd_solver
   use tblite_wavefunction_fermi, only : get_fermi_filling
   implicit none
   private

   public :: new_ptb_solver

   !> Abstract base class for electronic solvers
   type, public, extends(sygvd_solver) :: ptb_solver_type
      !> Linear regression cofactors for the orbital energies
      real(wp) :: keps_param
      real(wp) :: keps0_param
   contains
      procedure :: get_density
      procedure :: delete
   end type ptb_solver_type

contains

   subroutine new_ptb_solver(self, nel, kt, keps_param, keps0_param)
      type(ptb_solver_type), intent(out) :: self
      real(wp), intent(in) :: nel(:)
      real(wp), intent(in) :: kt
      real(wp), intent(in), optional :: keps_param
      real(wp), intent(in), optional :: keps0_param
      self%nel = nel
      self%kt = kt

      if (present(keps_param)) then
         self%keps_param = keps_param
      else
         self%keps_param = 0.0_wp 
      end if

      if (present(keps0_param)) then
         self%keps0_param = keps0_param
      else
         self%keps0_param = 0.0_wp 
      end if 
   end subroutine

   subroutine get_density(self, hmat, smat, eval, focc, density, error)
      !> Solver for the general eigenvalue problem
      class(ptb_solver_type), intent(inout) :: self
      !> Overlap matrix
      real(wp), contiguous, intent(in) :: smat(:, :)
      !> Hamiltonian matrix, contains eigenvectors on output
      real(wp), contiguous, intent(inout) :: hmat(:, :, :)
      !> Eigenvalues
      real(wp), contiguous, intent(inout) :: eval(:, :)
      !> Occupation numbers
      real(wp), contiguous, intent(inout) :: focc(:, :)
      !> Density matrix
      real(wp), contiguous, intent(inout) :: density(:, :, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: e_fermi
      integer :: nspin, spin, homo

      nspin = min(size(eval, 2), size(hmat, 3), size(density, 3))

      select case (nspin)
      case default
         call self%solve(hmat(:, :, 1), smat, eval(:, 1), error)
         if (allocated(error)) return

         ! Linear shift of the orbital energies
         eval(:, 1) = eval(:, 1) * (1.0_wp + self%keps_param) + self%keps0_param

         focc(:, :) = 0.0_wp 
         do spin = 1, 2
            call get_fermi_filling(self%nel(spin), self%kt, eval(:, 1), &
               & homo, focc(:, spin), e_fermi)
         end do
         focc(:, 1) = focc(:, 1) + focc(:, 2)
         call get_density_matrix(focc(:, 1), hmat(:, :, 1), density(:, :, 1))
         focc(:, 2) = focc(:, 1) - focc(:, 2)
      case (2)
         hmat(:, :, :) = 2*hmat
         do spin = 1, 2
            call self%solve(hmat(:, :, spin), smat, eval(:, spin), error)
            if (allocated(error)) return

            ! Linear shift of the orbital energies
            eval(:, 1) = eval(:, 1) * (1.0_wp + self%keps_param) + self%keps0_param

            call get_fermi_filling(self%nel(spin), self%kt, eval(:, spin), &
               & homo, focc(:, spin), e_fermi)
            call get_density_matrix(focc(:, spin), hmat(:, :, spin), &
               & density(:, :, spin))
         end do
      end select
   end subroutine get_density

   !> Get the density matrix from the coefficients and occupation numbers
   subroutine get_density_matrix(focc, coeff, pmat)
      !> Occupation numbers
      real(wp), intent(in) :: focc(:)
      !> Coefficients of the wavefunction
      real(wp), contiguous, intent(in) :: coeff(:, :)
      !> Density matrix to be computed
      real(wp), contiguous, intent(out) :: pmat(:, :)

      real(wp), allocatable :: scratch(:, :)
      integer :: iao, jao

      allocate(scratch(size(pmat, 1), size(pmat, 2)))
      !$omp parallel do collapse(2) default(none) schedule(runtime) &
      !$omp shared(scratch, coeff, focc, pmat) private(iao, jao)
      do iao = 1, size(pmat, 1)
         do jao = 1, size(pmat, 2)
            scratch(jao, iao) = coeff(jao, iao) * focc(iao)
         end do
      end do
      call gemm(scratch, coeff, pmat, transb='t')
   end subroutine get_density_matrix

   !> Delete the solver instance
   subroutine delete(self)
      !> Solver for the general eigenvalue problem
      class(ptb_solver_type), intent(inout) :: self

   end subroutine delete

#endif
end module xtb_ptb_solver
