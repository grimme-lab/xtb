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

!> Two-step SCF of the PTB method

module xtb_ptb_scf
   use mctc_io, only: structure_type
   use mctc_env, only: wp, error_type

   use tblite_basis_type, only: basis_type
   use tblite_context, only: context_type
   use tblite_scf_solver, only: solver_type

   use multicharge_model, only: mchrg_model_type

   use xtb_ptb_vdzp, only: add_vDZP_basis
   use xtb_ptb_param, only: kalphah0l, klalphaxc, &
   & nshell, max_shell, ptbGlobals, rf
   use xtb_ptb_overlaps, only: get_scaled_integrals
   use xtb_ptb_mmlpopanalysis, only: get_mml_overlaps
   use xtb_ptb_ncoord, only: ncoord_erf
   use xtb_ptb_corebasis, only: get_Vecp
   use xtb_ptb_data, only: TPTBData
   use xtb_ptb_hamiltonian, only: get_hamiltonian, get_selfenergy

   implicit none
   private

   public :: twostepscf

   real(wp), parameter :: default_cutoff = 25.0_wp

contains

   subroutine twostepscf(ctx, data, mol, bas, cbas, eeqmodel)
      !> Calculation context
      type(context_type), intent(inout) :: ctx
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> PTB parameterization data
      type(TPTBData), intent(in) :: data
      !> Basis set and core-valence basis set data
      type(basis_type), intent(in) :: bas, cbas
      !> Initialized EEQ model
      type(mchrg_model_type), intent(in) :: eeqmodel
      !> Electronic solver
      class(solver_type), allocatable :: solver
      !> Error type
      type(error_type), allocatable :: error
      !> (Scaled) overlap matrix
      real(wp), allocatable :: overlap(:, :), overlap_h0(:, :), overlap_xc(:, :)
      !> Mulliken-Loewdin overlap matrices
      real(wp) :: overlap_sx(bas%nao, bas%nao), overlap_soneminusx(bas%nao, bas%nao)
      !> Dipole integrals
      real(wp), allocatable :: dipole(:, :, :)
      !> Temporary array for exponent scaling factors specific to
      !> unique atoms in the molecule
      real(wp), allocatable :: expscal(:, :)
      !> Loop variables
      integer :: i, j, isp, izp
      !> Coordination numbers
      real(wp) :: cn_star(mol%nat), cn(mol%nat), cn_eeq(mol%nat)
      real(wp) :: radii(mol%nid)
      !> EEQ charges
      real(wp), allocatable :: q_eeq(:)
      !> Normalization factors
      real(wp), allocatable :: norm_overlap(:)
      !> Effective core potential
      real(wp), allocatable :: vecp(:, :)
      !> Effective Hamiltonian
      real(wp), allocatable :: hmat(:, :)
      !> Effective self-energies
      real(wp), allocatable :: levels(:)

      !> Solver for the effective Hamiltonian
      call ctx%new_solver(solver, bas%nao)

      allocate (expscal(max_shell, mol%nid), source=0.0_wp)
      call get_scaled_integrals(mol, overlap, dipole, norm=norm_overlap)

      !##### DEV WRITE #####
      write (*, *) "Standard overlap ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") overlap(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      ! write (*, *) "Dipole:"
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f12.6)', advance="no") dipole(1, i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################

      !> Set up a exponent scaling factors for
      !> new basis set for "H0 overlap"
      do isp = 1, mol%nid
         izp = mol%num(isp)
         expscal(:, isp) = kalphah0l(:, izp)
      end do
      call get_scaled_integrals(mol, overlap_h0, alpha_scal=expscal)
      !##### DEV WRITE #####
      write (*, *) "Overlap H0 scaled (SS) ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") overlap_h0(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################

      !> Set up a exponent scaling factors for
      !> new basis set for "V_XC overlap"
      do isp = 1, mol%nid
         izp = mol%num(isp)
         expscal(:, isp) = klalphaxc(:, izp)
      end do
      call get_scaled_integrals(mol, overlap_xc, alpha_scal=expscal)
      !##### DEV WRITE #####
      write (*, *) "Overlap XC scaled (SS) ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") overlap_xc(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################

      call get_mml_overlaps(bas, overlap, ptbGlobals%mlmix, overlap_sx, &
      & overlap_soneminusx)
      !##### DEV WRITE #####
      write (*, *) "Overlap S(1-x) ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") overlap_soneminusx(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      write (*, *) "Overlap S(x) ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") overlap_sx(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################

      !> Get first coordination number (" CN' ")
      call ncoord_erf(mol, ptbGlobals%kerfcn, default_cutoff, cn_star)

      !##### DEV WRITE #####
      write (*, *) "CN star:"
      do i = 1, mol%nat
         write (*, *) "Atom ", i, ":", cn_star(i)
      end do
      !#####################

      !> Get radii from PTB parameters for second coordination number (" CN ")
      do isp = 1, mol%nid
         izp = mol%num(isp)
         radii(isp) = rf(izp)
      end do
      !> Get second coordination number (" CN ")
      call ncoord_erf(mol, ptbGlobals%kerfcn, default_cutoff, cn, covrad=radii)

      !##### DEV WRITE #####
      write (*, *) "CN:"
      do i = 1, mol%nat
         write (*, *) "Atom ", i, ":", cn(i)
      end do
      !#####################

      !> Get third coordination number for EEQ model (" CN-EEQ ")
      call ncoord_erf(mol, ptbGlobals%kerfcn_eeq, default_cutoff, cn_eeq)

      !##### DEV WRITE #####
      write (*, *) "CN-EEQ:"
      do i = 1, mol%nat
         write (*, *) "Atom ", i, ":", cn_eeq(i)
      end do
      !#####################

      !> EEQ call
      allocate (q_eeq(mol%nat))
      call eeqmodel%solve(mol, cn_eeq, qvec=q_eeq)
      !##### DEV WRITE #####
      write (*, *) "EEQ charges:"
      do i = 1, mol%nat
         write (*, '(a,i0,a,f12.6)') "Atom ", i, ":", q_eeq(i)
      end do
      !#####################

      !> V_ECP via PTB core basis
      allocate (vecp(cbas%nao, bas%nao), source=0.0_wp)
      call get_Vecp(mol, data%corepotential, bas, cbas, norm_overlap, vecp)

      !##### DEV WRITE #####
      write (*, *) "V_ECP ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") vecp(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################

      !> Get the effective self-energies
      call get_selfenergy(mol, bas, data%hamiltonian, cn, cn_star, levels)
      !> Set up the effective Hamiltonian in the first iteration
      ! call get_hamiltonian(mol, bas, overlap, overlap_h0, overlap_xc, &
      ! & vecp, hmat)

   end subroutine twostepscf

end module xtb_ptb_scf

