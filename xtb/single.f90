! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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

!! ========================================================================
!  DRIVER FOR SINGLEPOINT CALCULATIONS
!  => call on (mini)scf to perform GFN-xTB single point calculations
!  => call on qmdff to run the QMDFF on a solvent file
!  => call on external driver library, interfaces
!     -> Turbomole, ORCA, driver (which is by itself an interface)
!  use the singlepoint function to get the appropiate call on the necessary
!  functions
!! ========================================================================
module single
   use iso_fortran_env, wp => real64
   implicit none

   character(len=*),private,parameter :: outfmt = &
      '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'

contains

subroutine singlepoint &
&                 (iunit,mol,wfn,calc, &
&                  egap,et,maxiter,prlevel,restart,lgrad,acc,etot,g,sigma,res)
   use iso_fortran_env, wp => real64
   use mctc_econv

!! ========================================================================
!  type definitions
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_calculator
   use tbdef_data
   use tbdef_pcem

!! ========================================================================
   use aoparam
   use setparam
   use fixparam
   use scanparam
   use sphereparam

!! ========================================================================
   use gbobc, only : lgbsa
   use scf_module, only : scf
   use qmdff,      only : ff_eg,ff_nonb,ff_hb
   use qcextern,   only : run_orca_egrad,run_mopac_egrad
   use peeq_module, only : peeq
   use embedding,  only : read_pcem
   implicit none

   integer, intent(in) :: iunit

!! ========================================================================
   type(tb_molecule), intent(inout) :: mol
   type(tb_wavefunction),intent(inout) :: wfn
   type(tb_calculator),intent(in) :: calc
   type(tb_pcem) :: pcem
   real(wp),intent(inout) :: egap
   real(wp),intent(in)    :: et
   integer, intent(in)    :: maxiter
   integer, intent(in)    :: prlevel
   logical, intent(in)    :: restart
   logical, intent(in)    :: lgrad
   real(wp),intent(in)    :: acc
   real(wp),intent(out)   :: etot
   real(wp),intent(out)   :: g(3,mol%n)
   type(scc_results),intent(out) :: res
   real(wp),intent(out)   :: sigma(3,3)
   integer  :: i,ich
   integer  :: mode_sp_run = 1
   real(wp) :: efix
   logical  :: inmol
   logical, parameter :: ccm = .true.
!  real(wp) :: efix1,efix2
!  real(wp),dimension(3,n) :: gfix1,gfix2

   call mol%update
   if (mol%npbc > 0) call generate_wsc(mol,mol%wsc,[1,1,1])

   etot = 0.0_wp
   efix = 0.0_wp
   g = 0.0_wp
   sigma = 0.0_wp

!! ========================================================================
!  external constrains which can be applied beforehand

   ! check for external point charge field
   call open_file(ich,pcem_file,'r')
   if (ich.ne.-1) then
      call read_pcem(ich,pcem)
      call close_file(ich)
   endif

!! ========================================================================
!  actual calculation
   select case(mode_extrun)
   case default
      call scf(iunit,mol,wfn,calc%basis,calc%param,pcem, &
         &   egap,et,maxiter,prlevel,restart,lgrad,acc,etot,g,res)

   case(p_ext_eht)
      call peeq &
         & (iunit,mol,wfn,calc%basis,calc%param,egap,et,prlevel,lgrad,ccm,acc,etot,g,sigma,res)

   case(p_ext_qmdff)
      call ff_eg  (mol%n,mol%at,mol%xyz,etot,g)
      call ff_nonb(mol%n,mol%at,mol%xyz,etot,g)
      call ff_hb  (mol%n,mol%at,mol%xyz,etot,g)

   case(p_ext_orca)
      call run_orca_egrad(mol%n,mol%at,mol%xyz,etot,g)

   case(p_ext_turbomole)
      call external_turbomole(mol%n,mol%at,mol%xyz,wfn%nel,wfn%nopen, &
         &                    lgrad,etot,g,res%dipole,lgbsa)

   case(p_ext_mopac)
      !call getmopacgrad(n,at,xyz,wfn%nopen,g,etot)
      call run_mopac_egrad(mol%n,mol%at,mol%xyz,etot,g)

   end select

!! ========================================================================
!  post processing of gradient and energy

   ! point charge embedding gradient file
   if (pcem%n > 0) then
      call open_file(ich,pcem_grad,'w')
      write(ich, '(i0)') pcem%n
      do i=1,pcem%n
         write(ich,'(3f17.12)')pcem%grd(1:3,i)
      enddo
      call close_file(ich)
   endif

!! ------------------------------------------------------------------------
!  various external potentials
   call constrain_pot(potset,mol%n,mol%at,mol%xyz,g,efix)
   call constrpot   (mol%n,mol%at,mol%xyz,g,efix)
   call cavity_egrad(mol%n,mol%at,mol%xyz,efix,g)
   call metadynamic (metaset,mol%n,mol%at,mol%xyz,efix,g)

!! ------------------------------------------------------------------------
!  fixing of certain atoms
!  print*,abs(efix/etot)
   etot = etot + efix
   res%e_total = etot
   res%gnorm = norm2(g)
   if (fixset%n.gt.0) then
      do i=1, fixset%n
         !print*,i,fixset%atoms(i)
         g(1:3,fixset%atoms(i))=0
      enddo
   endif

   if (prlevel.ge.2) then
      ! start with summary header
      if (.not.silent) then
         write(iunit,'(9x,53(":"))')
         write(iunit,'(9x,"::",21x,a,21x,"::")') "SUMMARY"
      endif
      write(iunit,'(9x,53(":"))')
      write(iunit,outfmt) "total energy      ", res%e_total,"Eh   "
      if (.not.silent.and.lgbsa) then
         write(iunit,outfmt) "total w/o Gsasa/hb", &
            &  res%e_total-res%g_sasa-res%g_hb-res%g_shift, "Eh   "
      endif
      write(iunit,outfmt) "gradient norm     ", res%gnorm,  "Eh/a0"
      write(iunit,outfmt) "HOMO-LUMO gap     ", res%hl_gap, "eV   "
      if (.not.silent) then
         if (verbose) then
            write(iunit,'(9x,"::",49("."),"::")')
            write(iunit,outfmt) "HOMO orbital eigv.", wfn%emo(wfn%ihomo),  "eV   "
            write(iunit,outfmt) "LUMO orbital eigv.", wfn%emo(wfn%ihomo+1),"eV   "
         endif
         write(iunit,'(9x,"::",49("."),"::")')
         if (gfn_method.eq.2) call print_gfn2_results(iunit,res,verbose,lgbsa)
         if (gfn_method.eq.1) call print_gfn1_results(iunit,res,verbose,lgbsa)
         if (gfn_method.eq.0) call print_gfn0_results(iunit,res,verbose,lgbsa)
         write(iunit,outfmt) "add. restraining  ", efix,       "Eh   "
         if (verbose) then
            write(iunit,'(9x,"::",49("."),"::")')
            write(iunit,outfmt) "atomisation energy", res%e_atom, "Eh   "
         endif
      endif
      write(iunit,'(9x,53(":"))')
      write(iunit,'(a)')
   endif

end subroutine singlepoint

subroutine print_gfn0_results(iunit,res,verbose,lgbsa)
   use tbdef_data
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
   use tbdef_data
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
   use tbdef_data
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


end module single
