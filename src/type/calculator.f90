! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

!> abstract calculator that hides implementation details from calling codes
module xtb_type_calculator
   use xtb_mctc_accuracy, only : wp
   use xtb_type_data
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_param, only : scc_parameter
   use xtb_type_wavefunction
   use xtb_setparam
   use xtb_fixparam
   use xtb_scanparam
   use xtb_sphereparam
   use xtb_solv_gbobc, only : lgbsa
   use xtb_qmdff, only : ff_eg,ff_nonb,ff_hb
   use xtb_extern_mopac, only : runMopac
   use xtb_extern_orca, only : runOrca
   use xtb_metadynamic
   use xtb_constrainpot
   implicit none

   public :: TCalculator
   private


   !> Base calculator
   type :: TCalculator

      real(wp) :: accuracy

   contains

      !> Perform single point calculation
      procedure :: singlepoint

      !> Write informative printout
      procedure :: writeInfo

   end type TCalculator

   character(len=*),private,parameter :: outfmt = &
      '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'


contains


subroutine singlepoint(self, env, mol, wfn, printlevel, restart, &
      & energy, gradient, sigma, hlgap, results)

   !> Source of the generated errors
   character(len=*), parameter :: source = 'type_calculator_singlepoint'

   !> Calculator instance
   class(TCalculator), intent(in) :: self

   !> Computational environment
   type(TEnvironment), intent(inout) :: env

   !> Molecular structure data
   type(TMolecule), intent(inout) :: mol

   !> Wavefunction data
   type(TWavefunction), intent(inout) :: wfn

   !> Print level for IO
   integer, intent(in) :: printlevel

   !> Restart from previous results
   logical, intent(in) :: restart

   !> Total energy
   real(wp), intent(out) :: energy

   !> Molecular gradient
   real(wp), intent(out) :: gradient(:, :)

   !> Strain derivatives
   real(wp), intent(out) :: sigma(:, :)

   !> HOMO-LUMO gap
   real(wp), intent(out) :: hlgap

   !> Detailed results
   type(scc_results), intent(out) :: results

   integer :: i,ich
   integer :: mode_sp_run = 1
   real(wp) :: efix
   logical :: inmol
   logical, parameter :: ccm = .true.
   logical :: exitRun

   call mol%update
   if (mol%npbc > 0) call generate_wsc(mol,mol%wsc)

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   hlgap = 0.0_wp
   efix = 0.0_wp

   ! ------------------------------------------------------------------------
   !  actual calculation
   select case(mode_extrun)
   case(p_ext_qmdff)
      call ff_eg  (mol%n,mol%at,mol%xyz,energy,gradient)
      call ff_nonb(mol%n,mol%at,mol%xyz,energy,gradient)
      call ff_hb  (mol%n,mol%at,mol%xyz,energy,gradient)

   case(p_ext_orca)
      call runOrca(env,mol,energy,gradient)

   case(p_ext_turbomole)
      call external_turbomole(mol%n,mol%at,mol%xyz,wfn%nel,wfn%nopen, &
         &                    .true.,energy,gradient,results%dipole,lgbsa)

   case(p_ext_mopac)
      call runMopac(env,mol%n,mol%at,mol%xyz,energy,gradient)

   end select

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Electronic structure method terminated", source)
      return
   end if

   ! ------------------------------------------------------------------------
   !  various external potentials
   call constrain_pot(potset,mol%n,mol%at,mol%xyz,gradient,efix)
   call constrpot   (mol%n,mol%at,mol%xyz,gradient,efix)
   call cavity_egrad(mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (metaset,mol%n,mol%at,mol%xyz,efix,gradient)

   ! ------------------------------------------------------------------------
   !  fixing of certain atoms
   !  print*,abs(efix/etot)
   energy = energy + efix
   results%e_total = energy
   results%gnorm = norm2(gradient)
   if (fixset%n.gt.0) then
      do i=1, fixset%n
         !print*,i,fixset%atoms(i)
         gradient(1:3,fixset%atoms(i))=0
      enddo
   endif

   if (printlevel.ge.2) then
      ! start with summary header
      if (.not.silent) then
         write(env%unit,'(9x,53(":"))')
         write(env%unit,'(9x,"::",21x,a,21x,"::")') "SUMMARY"
      endif
      write(env%unit,'(9x,53(":"))')
      write(env%unit,outfmt) "total energy      ", results%e_total,"Eh   "
      write(env%unit,outfmt) "gradient norm     ", results%gnorm,  "Eh/a0"
      write(env%unit,outfmt) "HOMO-LUMO gap     ", results%hl_gap, "eV   "
      write(env%unit,'(9x,53(":"))')
      write(env%unit,'(a)')
   endif

end subroutine singlepoint


subroutine print_gfn0_results(iunit,res,verbose,lgbsa)
   use xtb_type_data
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
   type(scc_results),    intent(in) :: res
   logical,intent(in) :: verbose,lgbsa
   write(iunit,outfmt) "H0 energy         ", res%e_elec, "Eh   "
   write(iunit,outfmt) "repulsion energy  ", res%e_rep,  "Eh   "
   write(iunit,outfmt) "electrostat energy", res%e_es,   "Eh   "
   write(iunit,outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
   !write(iunit,outfmt) "   -> Gborn       ", res%g_born, "Eh   " ! not saved
   write(iunit,outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
   !write(iunit,outfmt) "   -> Ghb         ", res%g_hb,   "Eh   " ! not saved
   write(iunit,outfmt) "   -> Gshift      ", res%g_shift,"Eh   "
   write(iunit,outfmt) "dispersion energy ", res%e_disp, "Eh   "
   write(iunit,outfmt) "short-range corr. ", res%e_xb,   "Eh   "
end subroutine print_gfn0_results

subroutine print_gfn1_results(iunit,res,verbose,lgbsa)
   use xtb_type_data
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
   type(scc_results),    intent(in) :: res
   logical,intent(in) :: verbose,lgbsa
   write(iunit,outfmt) "SCC energy        ", res%e_elec, "Eh   "
   write(iunit,outfmt) "-> electrostatic  ", res%e_es,   "Eh   "
   if (lgbsa) then
   write(iunit,outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
   write(iunit,outfmt) "   -> Gborn       ", res%g_born, "Eh   "
   write(iunit,outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
   write(iunit,outfmt) "   -> Ghb         ", res%g_hb,   "Eh   "
   write(iunit,outfmt) "   -> Gshift      ", res%g_shift,"Eh   "
   endif
   write(iunit,outfmt) "repulsion energy  ", res%e_rep,  "Eh   "
   write(iunit,outfmt) "dispersion energy ", res%e_disp, "Eh   "
   write(iunit,outfmt) "halogen bond corr.", res%e_xb,   "Eh   "
end subroutine print_gfn1_results

subroutine print_gfn2_results(iunit,res,verbose,lgbsa)
   use xtb_type_data
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
   type(scc_results),    intent(in) :: res
   logical,intent(in) :: verbose,lgbsa
   write(iunit,outfmt) "SCC energy        ", res%e_elec, "Eh   "
   write(iunit,outfmt) "-> isotropic ES   ", res%e_es,   "Eh   "
   write(iunit,outfmt) "-> anisotropic ES ", res%e_aes,  "Eh   "
   write(iunit,outfmt) "-> anisotropic XC ", res%e_axc,  "Eh   "
   write(iunit,outfmt) "-> dispersion     ", res%e_disp, "Eh   "
   if (lgbsa) then
   write(iunit,outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
   write(iunit,outfmt) "   -> Gborn       ", res%g_born, "Eh   "
   write(iunit,outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
   write(iunit,outfmt) "   -> Ghb         ", res%g_hb,   "Eh   "
   write(iunit,outfmt) "   -> Gshift      ", res%g_shift,"Eh   "
   endif
   write(iunit,outfmt) "repulsion energy  ", res%e_rep,  "Eh   "
end subroutine print_gfn2_results


subroutine writeInfo(self, unit, mol)

   !> Calculator instance
   class(TCalculator), intent(in) :: self

   !> Unit for I/O
   integer, intent(in) :: unit

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   select case(mode_extrun)
   case(p_ext_qmdff)
      call qmdff_header(unit)
   case(p_ext_orca, p_ext_mopac, p_ext_turbomole)
      call driver_header(unit)
   end select

end subroutine writeInfo


end module xtb_type_calculator
