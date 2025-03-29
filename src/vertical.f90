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

!-----------------------------------------------------------------------
!> Workflows employing multiple singlepoint calculations
!-----------------------------------------------------------------------
module xtb_vertical

   use xtb_mctc_accuracy, only: wp
   use xtb_type_calculator, only: TCalculator
   use xtb_type_data, only: scc_results
   use xtb_type_environment, only: TEnvironment
   use xtb_type_molecule, only: TMolecule, len
   use xtb_type_restart, only: TRestart

   implicit none
   private

   public :: vfukui

contains

!--------------------------------------------------
! Calculate Fukui indices
!--------------------------------------------------
subroutine vfukui(env, mol, chk, calc, fukui)

   implicit none
   !> Dummy-argument list
   class(TCalculator), intent(inout) :: calc
      !! instance of calculatior
   type(TEnvironment), intent(inout) :: env
      !! calculation environment
   type(TMolecule), intent(inout) :: mol
      !! molecular information
   type(TRestart), intent(inout) :: chk
   type(TRestart) :: wf_an, wf_cat

   ! fukui functions f(+), f(-), f(0)
   real(wp), intent(out) :: fukui(3,mol%n)
 
   type(scc_results) :: res
   real(wp) :: sigma(3,3)
   real(wp) :: etot2,egap, orig_chrg
   real(wp) :: g(3,mol%n)
   logical :: exist
   integer :: i, orig_uhf

   write(env%unit,'(a)')
   write(env%unit,'("Fukui index Calculation")')

   ! save original charge and uhf
   orig_chrg = mol%chrg
   orig_uhf = mol%uhf

   ! copy wavefunction
   wf_an%wfn = chk%wfn
   wf_cat%wfn = chk%wfn

   ! reduce the charge -> anion
   mol%chrg = mol%chrg - 1
   ! Create UHF entry for cation and anion -> we assume the UHF state to be equal for both species (adding or removing always towards the lowest multiplicity possible)
   if (mod(wf_an%wfn%nel, 2) == 0 .and. (mol%uhf == 0)) then
      ! if even number of electrons in the species (before charge modification), then we assume a closed-shell system
      ! -> Increase UHF to 1
      wf_an%wfn%nopen = 1
      wf_cat%wfn%nopen = 1
      ! mol object has to be modified as of today, since PTB sets its wavefunction only in the singlepoint calculation and checks for consistency with mol%uhf
      mol%uhf = 1
   elseif ((mod(wf_an%wfn%nel, 2) /= 0) .or. (mol%uhf > 0)) then
      ! open-shell system (before charge modification)
      ! -> Decrease UHF to 0
      if ((wf_an%wfn%nopen == 0) .or. (mol%uhf == 0)) then
         call env%terminate("Open shell system with closed shell wavefunction or UHF assignment.")
         ! Wavefunction does not match number of electrons (and UHF entry)
      end if
      wf_an%wfn%nopen = wf_an%wfn%nopen - 1
      wf_cat%wfn%nopen = wf_cat%wfn%nopen - 1
      ! mol object has to be modified as of today, since PTB sets its wavefunction only in the singlepoint calculation and checks for consistency with mol%uhf
      mol%uhf = mol%uhf - 1
   else
      call env%terminate("Invalid UHF and wavefunction combination.")
   end if

   ! Perform single point calculation for anion
   write(env%unit,'(a)')
   write(env%unit,'("Run single point for reduced species")')
   call calc%singlepoint(env, mol, wf_an, 1, exist, etot2, g, sigma, egap, res)

   ! increase the charge -> cation
   mol%chrg = mol%chrg + 2
   ! Perform single point calculation for cation
   write(env%unit,'(a)')
   write(env%unit,'("Run single point for oxidized species")')
   call calc%singlepoint(env, mol, wf_cat, 1, exist, etot2, g, sigma, egap, res)

   ! Calculate the fukui functions where N is the number of electrons
   ! see J. Am. Chem. Soc. 1986, 108, 19, 5708–5711 for equations
   ! keep in mind that their q are populations (p), 
   ! which are related to our q by q = Z-p 

   ! f(+) = q_(N_elec) -  q_(N_elec + 1) / neutral - anion
   fukui(1,:) = chk%wfn%q  - wf_an%wfn%q

   ! f(-) = q_(N_elec - 1) - q_(N_elec) / cation - neutral
   fukui(2,:) = wf_cat%wfn%q-chk%wfn%q

   ! f(0) = 0.5 * [q_(N_elec - 1) - q_(N_elec + 1)] / cation - anion
   fukui(3,:) = 0.5d0*(wf_cat%wfn%q-wf_an%wfn%q)

   ! Print out fukui functions
   write(env%unit,'(a)')
   write(env%unit,'("Fukui functions:")')
   write(env%unit, '(1x,"    #        f(+)     f(-)     f(0)")')
   do i=1,mol%n
      write(env%unit,'(i6,a4,2f9.3,2f9.3,2f9.3)') i, mol%sym(i), fukui(1,i), fukui(2,i), fukui(3,i)
   enddo
   write(env%unit,'(a)')
   mol%chrg = orig_chrg
   mol%uhf = orig_uhf

end subroutine vfukui

end module xtb_vertical
