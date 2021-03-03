! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> Extended tight binding calculator
module xtb_xtb_calculator
   use xtb_mctc_accuracy, only : wp
   use xtb_solv_gbsa, only : TBorn
   use xtb_solv_model, only : info, newSolvationModel, newBornModel
   use xtb_type_basisset, only : TBasisset
   use xtb_type_calculator, only : TCalculator
   use xtb_type_data
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_param, only : scc_parameter
   use xtb_type_pcem
   use xtb_type_solvation, only : TSolvation
   use xtb_type_restart, only : TRestart
   use xtb_type_wsc, only : tb_wsc
   use xtb_xtb_data, only : TxTBData
   use xtb_setparam
   use xtb_fixparam
   use xtb_scanparam
   use xtb_sphereparam
   use xtb_scf, only : scf
   use xtb_qmdff, only : ff_eg,ff_nonb,ff_hb
   use xtb_peeq, only : peeq
   use xtb_embedding, only : read_pcem
   use xtb_metadynamic
   use xtb_constrainpot
   implicit none
   interface
      subroutine generate_wsc(mol,wsc)
         import :: TMolecule, tb_wsc
         type(TMolecule), intent(inout) :: mol
         type(tb_wsc),    intent(inout) :: wsc
      end subroutine generate_wsc
   end interface

   private

   public :: TxTBCalculator


   !> Calculator interface for xTB based methods
   type, extends(TCalculator) :: TxTBCalculator

      !> Tight binding basis set
      type(TBasisset), allocatable :: basis

      !> Parametrisation data base
      type(TxTBData), allocatable :: xtbData

      !> Electronic temperature
      real(wp) :: etemp

      !> Maximum number of cycles for SCC convergence
      integer :: maxiter

      !> External potential
      type(tb_pcem) :: pcem

   contains

      !> Perform xTB single point calculation
      procedure :: singlepoint

      !> Write informative printout
      procedure :: writeInfo

   end type TxTBCalculator

   character(len=*),private,parameter :: outfmt = &
      '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'


contains


subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
      & energy, gradient, sigma, hlgap, results)

   !> Source of the generated errors
   character(len=*), parameter :: source = 'xtb_calculator_singlepoint'

   !> Calculator instance
   class(TxTBCalculator), intent(inout) :: self

   !> Computational environment
   type(TEnvironment), intent(inout) :: env

   !> Molecular structure data
   type(TMolecule), intent(inout) :: mol

   !> Wavefunction data
   type(TRestart), intent(inout) :: chk

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

   class(TSolvation), allocatable :: solvation
   type(TBorn), allocatable :: gbsa
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
   select case(self%xtbData%level)
   case(1, 2)
      if (allocated(self%solvation)) then
         call newSolvationModel(self%solvation, env, solvation, mol%at)
      end if
      call scf(env,mol,chk%wfn,self%basis,self%pcem,self%xtbData,solvation, &
         &   hlgap,self%etemp,self%maxiter,printlevel,restart,.true., &
         &   self%accuracy,energy,gradient,results)

   case(0)
      if (allocated(self%solvation)) then
         allocate(gbsa)
         call newBornModel(self%solvation, env, gbsa, mol%at)
      end if
      call peeq &
         & (env,mol,chk%wfn,self%basis,self%xtbData,gbsa,hlgap,self%etemp, &
         &  printlevel,.true.,ccm,self%accuracy,energy,gradient,sigma,results)

   end select

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Electronic structure method terminated", source)
      return
   end if

   ! ------------------------------------------------------------------------
   !  post processing of gradient and energy

   ! point charge embedding gradient file
   if (allocated(pcem_grad) .and. self%pcem%n > 0) then
      call open_file(ich,pcem_grad,'w')
      do i=1,self%pcem%n
         write(ich,'(3f12.8)')self%pcem%grd(1:3,i)
      enddo
      call close_file(ich)
   endif

   ! ------------------------------------------------------------------------
   !  various external potentials
   call constrain_pot(potset,mol%n,mol%at,mol%xyz,gradient,efix)
   call constrpot   (mol%n,mol%at,mol%xyz,gradient,efix)
   call cavity_egrad(mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (metaset,mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (rmsdset,mol%n,mol%at,mol%xyz,efix,gradient)

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
      if (.not.silent.and.allocated(self%solvation)) then
         write(env%unit,outfmt) "total w/o Gsasa/hb", &
            &  results%e_total-results%g_sasa-results%g_hb-results%g_shift, "Eh   "
      endif
      write(env%unit,outfmt) "gradient norm     ", results%gnorm,  "Eh/a0"
      write(env%unit,outfmt) "HOMO-LUMO gap     ", results%hl_gap, "eV   "
      if (.not.silent) then
         if (verbose) then
            write(env%unit,'(9x,"::",49("."),"::")')
            write(env%unit,outfmt) "HOMO orbital eigv.", chk%wfn%emo(chk%wfn%ihomo),  "eV   "
            write(env%unit,outfmt) "LUMO orbital eigv.", chk%wfn%emo(chk%wfn%ihomo+1),"eV   "
         endif
         write(env%unit,'(9x,"::",49("."),"::")')
         if (self%xtbData%level.eq.2) call print_gfn2_results(env%unit,results,verbose,allocated(self%solvation))
         if (self%xtbData%level.eq.1) call print_gfn1_results(env%unit,results,verbose,allocated(self%solvation))
         if (self%xtbData%level.eq.0) call print_gfn0_results(env%unit,results,verbose,allocated(self%solvation))
         write(env%unit,outfmt) "add. restraining  ", efix,       "Eh   "
         write(env%unit,outfmt) "total charge      ", sum(chk%wfn%q), "e    "
         if (verbose) then
            write(env%unit,'(9x,"::",49("."),"::")')
            write(env%unit,outfmt) "atomisation energy", results%e_atom, "Eh   "
         endif
      endif
      write(env%unit,'(9x,53(":"))')
      write(env%unit,'(a)')
   endif

end subroutine singlepoint


subroutine print_gfn0_results(iunit,res,verbose,lsolv)
   use xtb_type_data
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
   type(scc_results),    intent(in) :: res
   logical,intent(in) :: verbose,lsolv
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

subroutine print_gfn1_results(iunit,res,verbose,lsolv)
   use xtb_type_data
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
   type(scc_results),    intent(in) :: res
   logical,intent(in) :: verbose,lsolv
   write(iunit,outfmt) "SCC energy        ", res%e_elec, "Eh   "
   write(iunit,outfmt) "-> electrostatic  ", res%e_es,   "Eh   "
   if (lsolv) then
   write(iunit,outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
   write(iunit,outfmt) "   -> Gelec       ", res%g_born, "Eh   "
   write(iunit,outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
   write(iunit,outfmt) "   -> Ghb         ", res%g_hb,   "Eh   "
   write(iunit,outfmt) "   -> Gshift      ", res%g_shift,"Eh   "
   endif
   write(iunit,outfmt) "repulsion energy  ", res%e_rep,  "Eh   "
   write(iunit,outfmt) "dispersion energy ", res%e_disp, "Eh   "
   write(iunit,outfmt) "halogen bond corr.", res%e_xb,   "Eh   "
end subroutine print_gfn1_results

subroutine print_gfn2_results(iunit,res,verbose,lsolv)
   use xtb_type_data
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
   type(scc_results),    intent(in) :: res
   logical,intent(in) :: verbose,lsolv
   write(iunit,outfmt) "SCC energy        ", res%e_elec, "Eh   "
   write(iunit,outfmt) "-> isotropic ES   ", res%e_es,   "Eh   "
   write(iunit,outfmt) "-> anisotropic ES ", res%e_aes,  "Eh   "
   write(iunit,outfmt) "-> anisotropic XC ", res%e_axc,  "Eh   "
   write(iunit,outfmt) "-> dispersion     ", res%e_disp, "Eh   "
   if (lsolv) then
   write(iunit,outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
   write(iunit,outfmt) "   -> Gelec       ", res%g_born, "Eh   "
   write(iunit,outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
   write(iunit,outfmt) "   -> Ghb         ", res%g_hb,   "Eh   "
   write(iunit,outfmt) "   -> Gshift      ", res%g_shift,"Eh   "
   endif
   write(iunit,outfmt) "repulsion energy  ", res%e_rep,  "Eh   "
end subroutine print_gfn2_results


subroutine writeInfo(self, unit, mol)

   !> Calculator instance
   class(TxTBCalculator), intent(in) :: self

   !> Unit for I/O
   integer, intent(in) :: unit

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   call self%xtbData%writeInfo(unit, mol%at)

   if (allocated(self%solvation)) then
      call info(self%solvation, unit)
   end if

end subroutine writeInfo


end module xtb_xtb_calculator
