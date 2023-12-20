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

!> Response approximation within PTB

module xtb_ptb_response

   !> mctc-lib
   use mctc_env, only: error_type, wp
   use mctc_io, only: structure_type
   !> tblite-lib
   use tblite_basis_type, only: basis_type
   use tblite_context, only: context_type
   use tblite_wavefunction, only: wavefunction_type
   use tblite_integral_type, only: integral_type
   use tblite_container, only: container_cache
   use tblite_external_field, only: electric_field
   use tblite_scf_potential, only: potential_type, new_potential, add_pot_to_h1
   use tblite_scf_solver, only: solver_type
   use tblite_scf_iterator, only: get_qat_from_qsh
   !> xtb-ptb-lib
   use xtb_ptb_data, only: TPTBData
   use xtb_ptb_integral_types, only: aux_integral_type
   use xtb_ptb_param, only: ptbGlobals
   use xtb_ptb_scf, only: get_density
   use xtb_ptb_mmlpopanalysis, only: get_mml_shell_charges

   implicit none
   private

   public :: numgrad_polarizability

   !> Conversion factor from temperature to energy
   real(wp), parameter :: kt = 3.166808578545117e-06_wp

contains

   subroutine numgrad_polarizability(ctx, data, mol, bas, wfn, ints, auxints, &
         & delta, alpha, efield)
      !> Calculation context
      type(context_type), intent(inout) :: ctx
      !> Wavefunction of tblite type
      type(wavefunction_type), intent(inout) :: wfn
      type(wavefunction_type) :: wfn_tmp
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> PTB parameterization data
      type(TPTBData), intent(in) :: data
      !> Basis set and core-valence basis set data
      type(basis_type), intent(in) :: bas
      !> Integral type (last solved Hamiltonian matrix (not coeffs and not only H0!)
      !> should be on 'ints%hamiltonian')
      type(integral_type), intent(in) :: ints
      !> Auxiliary integral type
      type(aux_integral_type), intent(in) :: auxints
      !> Perturbation strength
      real(wp), intent(in) :: delta
      !> Static dipole polarizability
      real(wp), intent(out) :: alpha(3, 3)
      !> Error type
      type(error_type), allocatable :: error
      !> Potential type
      type(potential_type) :: pot
      !> Restart data for interaction containers
      type(container_cache) :: icache
      !> Electric field object
      type(electric_field) :: efield_object
      !> Electronic solver
      class(solver_type), allocatable :: solver
      !> Electronic entropy
      real(wp) :: ts
      !> Molecular dipole moment
      real(wp) :: dip_plus(3), dipminus(3)
      !> (optional) Electric field
      real(wp), intent(in), optional :: efield(:)
      real(wp) :: eff_ef(3), tmp_ef(3)
      !> Loop variables
      integer :: k, i, j

      !> Solver for the effective Hamiltonian
      call ctx%new_solver(solver, bas%nao)

      !##### DEV WRITE #####
      ! write (*, *) "Initial Hamiltonian matrix in response part ..."
      ! do i = 1, size(ints%hamiltonian, 1)
      !    do j = 1, size(ints%hamiltonian, 2)
      !       write (*, '(f11.7)', advance="no") ints%hamiltonian(i, j)
      !    end do
      !    write (*, '(/)', advance="no")
      ! end do
      !#####################

      !> Copy exact wavefunction from two step-scf to temporary wavefunction
      wfn_tmp = wfn

      alpha = 0.0_wp
      if (present(efield)) then
         eff_ef = efield
      else
         eff_ef = 0.0_wp
      end if

      cart_coord_loop: do k = 1, 3
         tmp_ef = eff_ef
         !> Add perturbation to electric field (minus sign corresponds to original PTB implementation)
         tmp_ef(k) = tmp_ef(k) - delta
         efield_object = electric_field(tmp_ef)
         !> Get potential and apply it to Hamiltonian
         call new_potential(pot, mol, bas, wfn%nspin)
         call pot%reset()
         call efield_object%update(mol, icache)
         call efield_object%get_potential(mol, icache, wfn_tmp, pot)
         call add_pot_to_h1(bas, ints, pot, wfn_tmp%coeff)
         !##### DEV WRITE #####
         ! write (*, *) "Matrix to solve ..."
         ! do i = 1, size(wfn_tmp%coeff, 1)
         !    do j = 1, size(wfn_tmp%coeff, 2)
         !       write (*, '(f11.7)', advance="no") wfn_tmp%coeff(i, j, 1)
         !    end do
         !    write (*, '(/)', advance="no")
         ! end do
         !#####################
         !> Solve effective Hamiltonian including the electric field
         call get_density(wfn_tmp, solver, ints, ts, error, ptbGlobals%geps, ptbGlobals%geps0)

         !##### DEV WRITE #####
         ! write (*, *) "Density matrix after adding field ..."
         ! do i = 1, size(wfn_tmp%density, 1)
         !    do j = 1, size(wfn_tmp%density, 2)
         !       write (*, '(f11.7)', advance="no") wfn_tmp%density(i, j, 1)
         !    end do
         !    write (*, '(/)', advance="no")
         ! end do
         !#####################

         !> Get updated atomic and shell charges
         call get_mml_shell_charges(bas, auxints%overlap_to_x, auxints%overlap_to_1_x, &
            & wfn_tmp%density, wfn_tmp%n0sh, wfn_tmp%qsh)
         call get_qat_from_qsh(bas, wfn_tmp%qsh, wfn_tmp%qat)

         write (*, *) "Atom charges after 1st iteration ..."
         do i = 1, mol%nat
            write (*, '(f8.4)') wfn_tmp%qat(i, 1)
         end do

         stop
         !> Reset Hamiltonian and enter one-step SCF routine to get updated wavefunction

         ! H = Vecp
         ! call adddsym(ndim, ffs, D(1, k), H)    ! perturb H with field only
         ! call onescf(n, ndim, nel, nopen, homo, at, rab, cns,&    ! and add 2nd iter part
         ! &              S, SS, H, Hdiag, focc, eT, scfpar, ves0, pshtmp, patmp, P1)
         ! call dipmom2(n, ndim, xyz, z, norm, P1, D, pnt, dip1)     ! get dipole moment

         ! call addsym(ndim, -ffs, Htmp, D(1, k), H)       ! other direction
         ! call solve3(ndim, nel, nopen, homo, eT, focc, H, S, P1)
         ! call mlpop2(n, ndim, P1, S1, S2, patmp, pshtmp)
         ! patmp = z - patmp
         ! H = Vecp
         ! call adddsym(ndim, -ffs, D(1, k), H)
         ! call onescf(n, ndim, nel, nopen, homo, at, rab, cns,&
         ! &              S, SS, H, Hdiag, focc, eT, scfpar, ves0, pshtmp, patmp, P1)
         ! call dipmom2(n, ndim, xyz, z, norm, P1, D, pnt, dip2)

         ! alpha(k, 1:3) = -(dip1(1:3) - dip2(1:3)) / (2_wp * ffs)                               ! numerical diff. dmu/dfield
      end do cart_coord_loop

      !> Symmetrization of polarizability tensor
      alpha(2, 1) = 0.5 * (alpha(2, 1) + alpha(1, 2))
      alpha(1, 2) = alpha(2, 1)
      alpha(3, 1) = 0.5 * (alpha(3, 1) + alpha(1, 3))
      alpha(1, 3) = alpha(3, 1)
      alpha(3, 2) = 0.5 * (alpha(3, 2) + alpha(2, 3))
      alpha(2, 3) = alpha(3, 2)

   end subroutine numgrad_polarizability

end module xtb_ptb_response
