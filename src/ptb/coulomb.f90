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

module xtb_ptb_coulomb
   use mctc_env, only: wp
   use mctc_io, only: structure_type

   use tblite_basis_type, only: basis_type
   use tblite_coulomb_charge_effective, only: harmonic_average
   use tblite_blas, only: symv
   use tblite_wavefunction, only: wavefunction_type
   use tblite_scf_potential, only: potential_type

   use dftd4_data_hardness, only: get_hardness

   implicit none
   private

   public :: coulomb_potential

   type :: coulomb_potential
      !> Effective Hubbard parameters
      real(wp), allocatable :: hubbard(:, :, :, :)
      !> effective gams (chemical hardness)
      real(wp), allocatable :: gam(:, :)
      !> Coulomb matrix
      real(wp), allocatable :: cmat(:, :)
      !> Ohno-Klopman contribution
      real(wp) :: cok
      !> Effective Hubbard-Derivatives for third-order contribution
      real(wp), allocatable :: hubbard_derivs(:)
   contains
      procedure :: init
      procedure :: update => get_coulomb_matrix
      procedure :: get_potential
   end type coulomb_potential

contains

   subroutine get_coulomb_matrix(self, mol, bas)
      !> Coulomb type
      class(coulomb_potential), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> tmp variables
      integer :: iat, jat, izp, jzp, ii, jj, ish, jsh
      real(wp) :: vec(3), r1, gam, tmp, r12

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, self, bas) &
      !$omp private(iat, izp, ii, ish, jat, jzp, jj, jsh, gam) &
      !$omp private(vec, r1, tmp, r12)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = bas%ish_at(iat)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            jj = bas%ish_at(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r1 = norm2(vec)
            r12 = r1**2
            do ish = 1, bas%nsh_at(iat)
               do jsh = 1, bas%nsh_at(jat)
                  gam = self%hubbard(jsh, ish, jat, iat)
                  tmp = (self%cok / sqrt(r12 + gam**(-2))) + &  !> Ohno-Klopman average
                     & ((1.0_wp - self%cok) / (r1 + gam**(-1))) !> Mataga-Nishimoto average
                  !$omp atomic
                  self%cmat(jj + jsh, ii + ish) = self%cmat(jj + jsh, ii + ish) + tmp
                  !$omp atomic
                  self%cmat(ii + ish, jj + jsh) = self%cmat(ii + ish, jj + jsh) + tmp
               end do
            end do
         end do
         do ish = 1, bas%nsh_at(iat)
            do jsh = 1, ish - 1
               gam = self%hubbard(jsh, ish, iat, iat)
               !$omp atomic
               self%cmat(ii + jsh, ii + ish) = self%cmat(ii + jsh, ii + ish) + gam
               !$omp atomic
               self%cmat(ii + ish, ii + jsh) = self%cmat(ii + ish, ii + jsh) + gam
            end do
            gam = self%hubbard(ish, ish, iat, iat)
            !$omp atomic
            self%cmat(ii + ish, ii + ish) = self%cmat(ii + ish, ii + ish) + gam 
         end do
      end do

      !##### DEV WRITE #####
      ! write (*, *) "cmat ..."
      ! do ish = 1, bas%nsh
      !    do jsh = 1, bas%nsh
      !       write (*, '(f13.9)', advance="no") self%cmat(ish, jsh)
      !    end do
      !    write (*, *)
      ! end do
      !#####################

   end subroutine get_coulomb_matrix

   subroutine init(self, mol, bas, q, gamsh, kqhubb, kok, kthirdoder)
      !> Effective Hubbard parameters
      class(coulomb_potential), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> Atomic charges
      real(wp), intent(in) :: q(:)
      !> Shell-wise gamma scaling factors
      real(wp), intent(in) :: gamsh(:, :)
      !> Scaling factor for the charge dependent part of the Hubbard parameters
      real(wp), intent(in) :: kqhubb
      !> Ohno-Klopman contribution
      real(wp), intent(in) :: kok
      !> Scaling for third-order contributions
      real(wp), intent(in) :: kthirdoder(:)

      call init_hubbard(self, mol, bas, q, gamsh, kqhubb, kok)
      call init_hubbard_derivs(self, mol, kthirdoder)

   end subroutine init

   subroutine init_hubbard(self, mol, bas, q, gamsh, kqhubb, kok)
      !> Effective Hubbard parameters
      class(coulomb_potential), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> Atomic charges
      real(wp), intent(in) :: q(:)
      !> Shell-wise gamma scaling factors
      real(wp), intent(in) :: gamsh(:, :)
      !> Scaling factor for the charge dependent part of the Hubbard parameters
      real(wp), intent(in) :: kqhubb
      !> Ohno-Klopman contribution
      real(wp), intent(in) :: kok

      integer :: izp, iat, iid, ish, jsh, jat

      !> Initialize variables
      if (.not. allocated(self%gam)) allocate (self%gam(maxval(bas%nsh_id), mol%nat), source=0.0_wp)
      if (.not. allocated(self%hubbard)) then
         allocate (self%hubbard(maxval(bas%nsh_id), maxval(bas%nsh_id), mol%nat, mol%nat), source=0.0_wp)
      end if
      if (.not. allocated(self%cmat)) then
         allocate (self%cmat(bas%nsh, bas%nsh), source=0.0_wp)
      else
         self%cmat = 0.0_wp
      end if
      self%cok = kok

      !> Atom-individual chemical hardnesses per shell; Eq. 19
      do iat = 1, mol%nat
         iid = mol%id(iat)
         izp = mol%num(iid)
         do ish = 1, bas%nsh_id(iid)
            self%gam(ish, iat) = (1.0_wp + kqhubb * q(iat)) * get_hardness(izp) * &
               & gamsh(ish, iid)
         end do
      end do

      !> Effective Hubbard parameters; Eq. 16
      do iat = 1, mol%nat
         do jat = 1, mol%nat
            self%hubbard(:, :, jat, iat) = 0.0_wp
            do ish = 1, bas%nsh_at(iat)
               do jsh = 1, bas%nsh_at(jat)
                  self%hubbard(jsh, ish, jat, iat) = &
                     & harmonic_average(self%gam(ish, iat), self%gam(jsh, jat))
               end do
            end do
         end do
      end do

   end subroutine init_hubbard

   subroutine init_hubbard_derivs(self, mol, kthirdoder)
      !> Effective Hubbard parameters
      class(coulomb_potential), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Scaling for third-order contributions
      real(wp), intent(in) :: kthirdoder(:)

      integer :: iat, iid

      !> Third-order on-center electrostatics; Eq. 17
      if (.not. allocated(self%hubbard_derivs)) then
         allocate (self%hubbard_derivs(mol%nat), source=0.0_wp)
      end if
      do iat = 1, mol%nat
         iid = mol%id(iat)
         self%hubbard_derivs(iat) = 2.0_wp * kthirdoder(iid)
      end do

   end subroutine init_hubbard_derivs

   subroutine get_potential(self, wfn, pot)
      !> Coulomb type
      class(coulomb_potential), intent(inout) :: self
      !> Wavefunction of tblite type
      type(wavefunction_type), intent(in) :: wfn
      !> Instance of the density dependent potential
      type(potential_type), intent(inout) :: pot

      call get_potential_secondorder(self, wfn, pot)
      call get_potential_thirdorder(self, wfn, pot)

   end subroutine get_potential

   subroutine get_potential_secondorder(self, wfn, pot)
      !> Coulomb type
      class(coulomb_potential), intent(inout) :: self
      !> Wavefunction of tblite type
      type(wavefunction_type), intent(in) :: wfn
      !> Instance of the density dependent potential
      type(potential_type), intent(inout) :: pot
      !> tmp variables
      integer :: i

      call symv(self%cmat, wfn%qsh(:, 1), pot%vsh(:, 1), beta=1.0_wp)

      !##### DEV WRITE #####
      ! write (*, *) "pot%vsh ..."
      ! do i = 1, size(pot%vsh, 1)
      !    write (*, *) i, pot%vsh(i,1)
      ! end do
      !#####################

   end subroutine get_potential_secondorder

   subroutine get_potential_thirdorder(self, wfn, pot)
      !> Coulomb type
      class(coulomb_potential), intent(inout) :: self
      !> Wavefunction of tblite type
      type(wavefunction_type), intent(in) :: wfn
      !> Instance of the density dependent potential
      type(potential_type), intent(inout) :: pot
      !> tmp variables
      integer :: iat

      do iat = 1, size(wfn%qat, 1)
         pot%vat(iat, 1) = pot%vat(iat, 1) + wfn%qat(iat, 1)**2 * self%hubbard_derivs(iat)
         !##### DEV WRITE #####
         ! write (*, *) "q^2 * self_hubbard_derivs: ",  wfn%qat(iat, 1)**2 * self%hubbard_derivs(iat)
         ! write (*, *) "qat: ", wfn%qat(iat, 1)
         !#####################
      end do

   end subroutine get_potential_thirdorder

end module xtb_ptb_coulomb
