! This file is part of xtb.
!
! Copyright (C) 2019-2020 Stefan Grimme
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

!> Force field calculator
module xtb_gfnff_calculator
   use xtb_mctc_accuracy, only : wp
   use xtb_solv_gbsa, only : TBorn, init
   use xtb_solv_model, only : info, newBornModel
   use xtb_type_calculator, only : TCalculator
   use xtb_type_data
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart
   use xtb_setparam
   use xtb_fixparam
   use xtb_scanparam
   use xtb_sphereparam
   use xtb_metadynamic
   use xtb_constrainpot
   use xtb_gfnff_param, only : make_chrg,gff_print
   use xtb_gfnff_data, only : TGFFData
   use xtb_gfnff_topology, only : TGFFTopology
   use xtb_gfnff_neighbourlist, only : TGFFNeighbourList
   use xtb_gfnff_generator, only : TGFFGenerator
   use xtb_gfnff_eg
   implicit none
   private

   public :: TGFFCalculator, newGFFCalculator


   !> calculator interface for xTB based methods
   type, extends(TCalculator) :: TGFFCalculator

      type(TGFFData) :: param
      type(TGFFGenerator) :: gen
      type(TGFFTopology) :: topo
      logical :: update
      integer :: version

   contains

      !> Perform xTB single point calculationV
      procedure :: singlepoint

      !> Write informative printout
      procedure :: writeInfo

   end type TGFFCalculator

   character(len=*),private,parameter :: outfmt = &
      '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'


contains


subroutine newGFFCalculator(env, mol, calc, fname, restart, version)
   
   use xtb_gfnff_param
   use xtb_gfnff_setup, only : gfnff_setup
   use xtb_disp_dftd4, only : newD3Model

   character(len=*), parameter :: source = 'main_setup_newGFFCalculator'

   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(in) :: mol

   type(TGFFCalculator), intent(out) :: calc

   character(len=*), intent(in) :: fname

   logical, intent(in) :: restart

   integer, intent(in), optional :: version

   integer :: ich
   logical :: exist, okbas
   logical :: exitRun

   if (present(version)) then
      calc%version = version
   else
      calc%version = gffVersion%angewChem2020_2
   end if

   call calc%topo%zero
   calc%update = .true.

   ! global accuracy factor similar to acc in xtb used in SCF !
   calc%accuracy = 0.1_wp
   if (mol%n > 10000) then
      calc%accuracy = 2.0_wp
   end if

   ! obtain the parameter file !
   call open_file(ich, fname, 'r')
   exist = ich /= -1
   if (exist) then
      call gfnff_read_param(ich, calc%param)
      call close_file(ich) 
      
   ! no parameter file, try to load internal version !
   else
      call gfnff_load_param(calc%version, calc%param, exist)
      if (.not.exist) then
         call env%error('Parameter file '//fname//' not found!', source)
         return
      end if
   endif

   call newD3Model(calc%topo%dispm, mol%n, mol%at)

   call gfnff_setup(env, set%verbose, restart, mol, &
      & calc%gen, calc%param, calc%topo, calc%accuracy, calc%version)

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Could not create force field calculator", source)
      return
   end if

end subroutine newGFFCalculator


!> GFN-FF wrapper for the single point energy evaluation  
subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
      & energy, gradient, sigma, hlgap, results)

   !> Source of the generated errors
   character(len=*), parameter :: source = 'gfnff_calculator_singlepoint'

   !> Calculator instance
   class(TGFFCalculator), intent(inout) :: self

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

   !> HOMO-LUMO gab
   real(wp), intent(out) :: hlgap

   !> Detailed results
   type(scc_results), intent(out) :: results

   type(TBorn), allocatable :: solvation
   integer :: i,ich
   integer :: mode_sp_run = 1
   real(wp) :: efix
   logical :: inmol
   logical, parameter :: ccm = .true.
   logical :: exitRun
   logical :: pr
   
   !> GFN-FF geometry optimization
   logical :: optpr

   !-------!
   ! setup !
   !-------!

   call mol%update

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:,:) = 0.0_wp
   hlgap = 0.0_wp
   efix = 0.0_wp

   if (allocated(self%solvation)) then
      allocate(solvation)
      call newBornModel(self%solvation, env, solvation, mol%at)
   end if

   ! to distinguish optimization, final sp and sp ! 
   if ((set%runtyp.eq.p_run_opt).or.(set%runtyp.eq.p_run_ohess).or. &
      &   (set%runtyp.eq.p_run_omd).or.(set%runtyp.eq.p_run_screen).or. &
      &   (set%runtyp.eq.p_run_metaopt)) then
      optpr = printlevel < 2
   else
      optpr = .false.
   endif
   
   pr = gff_print .and. printlevel > 0
   
   !--------------------!
   ! actual calculation !
   !--------------------!

   call gfnff_eg(env,pr,mol%n,nint(mol%chrg),mol%at,mol%xyz,make_chrg, &
      & gradient,energy,results,self%param,self%topo,chk%nlist,solvation,&
      & self%update,self%version,self%accuracy,minpr=optpr)

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Force-field method terminated", source)
      return
   end if

   ! ---------------------------------------!
   ! post processing of gradient and energy !
   !----------------------------------------!
   
   ! various external potentials !
   call constrain_pot(potset,mol%n,mol%at,mol%xyz,gradient,efix)
   call constrpot   (mol%n,mol%at,mol%xyz,gradient,efix)
   call cavity_egrad(mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (metaset,mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (rmsdset,mol%n,mol%at,mol%xyz,efix,gradient)

   ! fixing of certain atoms !
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

   if (printlevel.ge.2)  then
      ! start with summary header !
      if (.not.set%silent) then
         write(env%unit,'(9x,53(":"))')
         write(env%unit,'(9x,"::",21x,a,21x,"::")') "SUMMARY"
      endif
      write(env%unit,'(9x,53(":"))')
      write(env%unit,outfmt) "total energy      ", results%e_total,"Eh   "
      if (.not.set%silent.and.allocated(solvation)) then
         write(env%unit,outfmt) "total w/o Gsolv   ", &
            &  results%e_total-results%g_solv, "Eh   "
         write(env%unit,outfmt) "total w/o Gsasa/hb", &
            &  results%e_total-results%g_sasa-results%g_hb-results%g_shift, "Eh   "
      endif
      write(env%unit,outfmt) "gradient norm     ", results%gnorm,  "Eh/a0"
      if (.not.set%silent) then
         write(env%unit,'(9x,"::",49("."),"::")')
         call print_gfnff_results(env%unit,results,set%verbose,allocated(solvation))
         write(env%unit,outfmt) "add. restraining  ", efix,       "Eh   "
         write(env%unit,outfmt) "total charge      ", sum(chk%nlist%q), "e    "
         if (set%verbose) then
            write(env%unit,'(9x,"::",49("."),"::")')
            write(env%unit,outfmt) "atomisation energy", results%e_atom, "Eh   "
         endif
      endif
      write(env%unit,'(9x,53(":"))')
      write(env%unit,'(a)')
      if (.not.set%silent) then
         call open_file(ich,'gfnff_charges','w')
         call print_charges(ich,mol%n,chk%nlist%q)
         call close_file(ich)
      end if
    
    else if (set%mode_extrun .eq. p_ext_oniom) then
      write(env%unit,outfmt) "total energy      ", results%e_total,"Eh   "
      write(env%unit,outfmt) "gradient norm     ", results%gnorm,  "Eh/a0"
    
    endif

end subroutine singlepoint

subroutine print_gfnff_results(iunit,res_gff,verbose,lsolv)
   
   use xtb_type_data
   
   !> file handle (usually output_unit=6)
   integer, intent(in) :: iunit 
   type(scc_results),    intent(in) :: res_gff
   logical,intent(in) :: verbose,lsolv
   
   write(iunit,outfmt) "bond energy       ", res_gff%e_bond, "Eh   "
   write(iunit,outfmt) "angle energy      ", res_gff%e_angl, "Eh   "
   write(iunit,outfmt) "torsion energy    ", res_gff%e_tors, "Eh   "
   write(iunit,outfmt) "repulsion energy  ", res_gff%e_rep,  "Eh   "
   write(iunit,outfmt) "electrostat energy", res_gff%e_es,   "Eh   "
   write(iunit,outfmt) "dispersion energy ", res_gff%e_disp, "Eh   "
   write(iunit,outfmt) "HB energy         ", res_gff%e_hb,   "Eh   "
   write(iunit,outfmt) "XB energy         ", res_gff%e_xb,   "Eh   "
   write(iunit,outfmt) "bonded atm energy ", res_gff%e_batm, "Eh   "
   write(iunit,outfmt) "external energy   ", res_gff%e_ext,  "Eh   "
   if (lsolv) then
      write(iunit,outfmt) "-> Gsolv          ", res_gff%g_solv, "Eh   "
      write(iunit,outfmt) "   -> Gborn       ", res_gff%g_born, "Eh   "
      write(iunit,outfmt) "   -> Gsasa       ", res_gff%g_sasa, "Eh   "
      write(iunit,outfmt) "   -> Ghb         ", res_gff%g_hb,   "Eh   "
      write(iunit,outfmt) "   -> Gshift      ", res_gff%g_shift,"Eh   "
   end if
end subroutine print_gfnff_results

subroutine writeInfo(self, unit, mol)

   !> Calculator instance
   class(TGFFCalculator), intent(in) :: self

   !> Unit for I/O
   integer, intent(in) :: unit

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   select case(set%mode_extrun)
   case(p_ext_gfnff)
     call gfnff_header(unit,self%version)
   end select

   if (allocated(self%solvation)) then
      call info(self%solvation, unit)
   end if

end subroutine writeInfo

subroutine print_charges(ifile,n,q)
   implicit none
   integer, intent(in)  :: ifile
   integer, intent(in)  :: n
   real(wp),intent(in)  :: q(n)
   integer :: i
   if (ifile.ne.-1) then
      do i = 1, n
         write(ifile,'(f14.8)') q(i)
      end do
   end if
end subroutine print_charges

end module xtb_gfnff_calculator
