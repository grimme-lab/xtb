! This file is part of xtb.
! SPDX-Identifier: LGPL-3.0-or-later
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

!> Bridge between tblite's Foster-Boys localization post-processing and
!> xtb's backend-independent LMO classification (xtb_local)
module xtb_tblite_local
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_data, only : scc_results
   use xtb_local, only : classify_lmo, get_ct_populations
#if WITH_TBLITE
   use tblite_basis_type, only : basis_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_results, only : results_type
#endif
   implicit none
   private

   public :: get_tblite_lmo

contains

!> Retrieve the Foster-Boys localized orbitals computed by tblite post-processing
#if WITH_TBLITE
subroutine get_tblite_lmo(env, mol, bas, wfn, tblite_results, etot, results)
#else
subroutine get_tblite_lmo(env, mol, etot, results)
#endif
   character(len=*), parameter :: source = 'tblite_local_get_tblite_lmo'
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol
#if WITH_TBLITE
   !> Basis set data
   type(basis_type), intent(in) :: bas
   !> Converged wavefunction
   type(wavefunction_type), intent(in) :: wfn
   !> Results container with the localization post-processing dictionary
   type(results_type), intent(in) :: tblite_results
#endif
   !> Total energy of the single point calculation
   real(wp), intent(in) :: etot
   !> Detailed results, xTB-IFF/docking part is populated here
   type(scc_results), intent(inout) :: results

#if WITH_TBLITE
   real(wp), allocatable :: cmo_lmo(:,:,:), centers(:,:,:), wbo(:,:,:)
   real(wp), allocatable :: f(:), ecent(:,:), qhl(:,:)
   real(wp) :: enhomo,enlumo,diptot
   integer :: nocc, spin, ilmo, nlmo

   if (.not.allocated(tblite_results%overlap)) then
      call env%warning("Overlap integrals not available for orbital "// &
         & "localization, skipping", source)
      return
   end if

   call tblite_results%dict%get_entry("localized-orbitals", cmo_lmo)
   call tblite_results%dict%get_entry("localized-centers", centers)
   call tblite_results%dict%get_entry("bond-orders", wbo)
   if (.not.allocated(cmo_lmo) .or. .not.allocated(centers) .or. .not.allocated(wbo)) then
      call env%warning("Could not find localized orbitals in tblite results, "// &
         & "skipping", source)
      return
   end if

   ! Charge-transfer terms of xTB-IFF are derived from the alpha-spin
   nocc = int(merge(wfn%nel(1)+1.0_wp, wfn%nel(1), mod(wfn%nel(1), 1.0_wp) > 0.5_wp))

   allocate(qhl(mol%n, 2), source=0.0_wp)
   call get_ct_populations(mol%n, bas%nao, nocc, bas%ao2at, tblite_results%overlap, &
      & wfn%coeff(:, :, 1), wfn%emo(:, 1), enhomo, enlumo, qhl)

   diptot = norm2(sum(wfn%dpat(:, :, 1), 2) + matmul(mol%xyz, wfn%qat(:, 1)))

   allocate(results%iff_results)
   call results%iff_results%allocateIFFResults(mol%n)

   ! For an open-shell (spin-unrestricted) calculation, tblite localizes
   ! each spin channel independently; classify and accumulate both into
   ! the same xTB-IFF LMO list.
   ilmo = 0
   nlmo = 0
   do spin = 1, wfn%nspin
      nocc = int(merge(wfn%nel(spin)+1.0_wp, wfn%nel(spin), &
         & mod(wfn%nel(spin), 1.0_wp) > 0.5_wp))
      if (nocc <= 0) cycle

      allocate(ecent(nocc, 4), source=0.0_wp)
      ecent(1:nocc, 1:3) = transpose(centers(1:3, 1:nocc, spin))

      ! Diagonal LMO Fock matrix element is not exposed by tblite localization
      allocate(f(nocc), source=0.0_wp)

      call classify_lmo(mol%n, mol%at, mol%xyz, wfn%qat(:, 1), bas%nao, nocc, bas%ao2at, &
         & tblite_results%overlap, cmo_lmo(:, 1:nocc, spin), f, wbo(:, :, 1), ecent, etot, &
         & diptot, enlumo, enhomo, qhl, results, ilmo0=ilmo, nlmo0=nlmo, &
         & islot_out=ilmo, nlmo_out=nlmo)

      deallocate(ecent, f)
   end do
#else
   call env%error("Compiled without support for tblite library", source)
#endif

end subroutine get_tblite_lmo

end module xtb_tblite_local
