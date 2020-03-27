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

!! ========================================================================
!  DRIVER FOR SINGLEPOINT CALCULATIONS
!  => call on (mini)scf to perform GFN-xTB single point calculations
!  => call on qmdff to run the QMDFF on a solvent file
!  => call on external driver library, interfaces
!     -> Turbomole, ORCA, driver (which is by itself an interface)
!  use the singlepoint function to get the appropiate call on the necessary
!  functions
!! ========================================================================
module xtb_single
   use xtb_mctc_accuracy, only : wp
   implicit none

   character(len=*),private,parameter :: outfmt = &
      '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'

contains

subroutine singlepoint &
&                 (env,mol,wfn,calc, &
&                  egap,et,maxiter,prlevel,restart,lgrad,acc,etot,g,sigma,res)
   use xtb_mctc_convert

!! ========================================================================
!  type definitions
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_wavefunction
   use xtb_type_calculator
   use xtb_type_data
   use xtb_type_pcem

!! ========================================================================
   use xtb_setparam
   use xtb_fixparam
   use xtb_scanparam
   use xtb_sphereparam

!! ========================================================================
   use xtb_solv_gbobc, only : lgbsa
   use xtb_scf, only : scf
   use xtb_qmdff, only : ff_eg,ff_nonb,ff_hb
   use xtb_extern_mopac, only : runMopac
   use xtb_extern_orca, only : runOrca
   use xtb_peeq, only : peeq
   use xtb_embedding, only : read_pcem
   implicit none

   character(len=*), parameter :: source = 'single'

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(inout) :: mol
   type(TWavefunction),intent(inout) :: wfn
   class(TCalculator), intent(in) :: calc
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
   logical :: exitRun
!  real(wp) :: efix1,efix2
!  real(wp),dimension(3,n) :: gfix1,gfix2

   call calc%singlepoint(env, mol, wfn, prlevel, restart, etot, g, sigma, egap, res)

end subroutine singlepoint


end module xtb_single
