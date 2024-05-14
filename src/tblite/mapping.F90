#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

module xtb_tblite_mapping
   
   use xtb_type_environment, only: TEnvironment
   use xtb_type_restart, only: TRestart
   use xtb_type_data, only: scc_results
   use xtb_type_molecule, only: TMolecule
   use xtb_mctc_accuracy, only: wp   
   use mctc_io_convert, only: autoev
   use xtb_type_wavefunction, only: TWavefunction

#if WITH_TBLITE
   use tblite_basis_type, only: basis_type
#endif
   implicit none
   private
   public :: convert_tblite_to_results

#if WITH_TBLITE
   public :: convert_tblite_to_wfn
   interface assignment(=)
      module procedure :: from_tblite_wfn
      module procedure :: from_tblite_basis
      module procedure :: from_struc
   end interface assignment(=)
#endif
contains

#if WITH_TBLITE
!> convert tblite wavefunction to xtb wavefunction
subroutine convert_tblite_to_wfn(env, bas, mol, chk, wbo)

   !> computational environment
   type(TEnvironment), intent(inout) :: env

   !> basis set
   type(basis_type), intent(in) :: bas
   
   !> molecular structure
   type(TMolecule), intent(in) :: mol

   !> restart data
   type(TRestart), intent(inout) :: chk

   !> wiberg bond orders
   real(wp), optional, intent(in) :: wbo(:,:,:)

   call chk%wfn%allocate(mol%n, bas%nsh, bas%nao)

   ! assignment interfaces !
   chk%wfn = chk%tblite
   chk%wfn = bas
   chk%wfn = mol
   
   if (present(wbo)) chk%wfn%wbo = wbo(:,:,1) ! only 

end subroutine convert_tblite_to_wfn
#endif

!> convert tblite results to xtb results
subroutine convert_tblite_to_results(results, mol, chk, energy, converged, gradient) 

   !> scc results
   type(scc_results), intent(inout) :: results

   !> molecular structure
   type(TMolecule), intent(in) :: mol

   !> restart data
   type(TRestart), intent(in) :: chk

   !> SCC energy
   real(wp), intent(in) :: energy

   !> convergence flag
   logical, intent(in) :: converged

   !> (analytical) gradients
   real(wp), optional, intent(in) :: gradient(:,:)

   
   results%e_total = energy
   results%converged = converged

#if WITH_TBLITE
   
   ! do not overwrite dipole moments, if calculated (PTB case) ! 
   if (all(results%dipole == 0.0_wp)) then
      results%dipole = sum(chk%tblite%dpat(:, :, 1), 2) + matmul(mol%xyz, chk%tblite%qat(: ,1))
   endif
   results%hl_gap = (chk%tblite%emo(chk%tblite%homo(1) + 1, 1) - chk%tblite%emo(chk%tblite%homo(1), 1)) * autoev

#endif

   if (present(gradient)) results%gnorm = norm2(gradient)

end subroutine convert_tblite_to_results


#if WITH_TBLITE

subroutine from_tblite_wfn(wfn, tblite)

   use tblite_wavefunction, only : wavefunction_type
   
   type(TWavefunction), intent(inout) :: wfn 
   type(wavefunction_type), intent(in) :: tblite

   wfn%dipm = tblite%dpat(:, :, 1)
   wfn%nel = nint(tblite%nocc)
   wfn%P = tblite%density(:, :, 1)
   wfn%q = tblite%qat(:, 1)
   wfn%qsh = tblite%qsh(:, 1)
   wfn%focca = tblite%focc(:, 1)
   wfn%foccb = 0.0_wp
   wfn%focc(:) = tblite%focc(:, 1)
   wfn%emo = tblite%emo(:, 1) * autoev
   wfn%C = tblite%coeff(:, :, 1)
   wfn%ihomo = tblite%homo(1)
   wfn%ihomoa = tblite%homo(1)
   wfn%ihomob = tblite%homo(2)
   wfn%qp = tblite%qpat(:, :, 1)

end subroutine from_tblite_wfn

subroutine from_tblite_basis(wfn, bas)

   use tblite_basis_type, only: basis_type
   
   type(TWavefunction), intent(inout) :: wfn
   type(basis_type), intent(in) :: bas
   
   wfn%nshell = bas%nsh
   wfn%nao = bas%nao

end subroutine from_tblite_basis

subroutine from_struc(wfn, mol)

   type(TWavefunction), intent(inout) :: wfn
   type(TMolecule), intent(in) :: mol

   wfn%n = mol%n
   wfn%nopen = mol%uhf

end subroutine from_struc
#endif

#if ! WITH_TBLITE
subroutine feature_not_implemented(env)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env

   call env%error("Compiled without support for tblite library")
end subroutine feature_not_implemented
#endif

end module xtb_tblite_mapping