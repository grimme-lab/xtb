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
! You should have received a copy of the GNU Lesser General Public Licen
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
   type(TRestart) :: wf_p, wf_m
   
   real(wp), intent(out) :: fukui(3,mol%n)
 
   type(scc_results) :: res
   real(wp) :: sigma(3,3)
   real(wp) :: etot2,egap
   real(wp) :: g(3,mol%n)
   logical :: exist
   integer :: i
   
   write(env%unit,'(a)')
   write(env%unit,'("Fukui index Calculation")')
   wf_p%wfn = chk%wfn
   wf_m%wfn = chk%wfn
   mol%chrg = mol%chrg - 1
   if (mod(wf_p%wfn%nel,2).ne.0) wf_p%wfn%nopen = 1
   call calc%singlepoint(env, mol, wf_p, 1, exist, etot2, g, sigma, egap, res)
   fukui(1,:) = wf_p%wfn%q-chk%wfn%q

   mol%chrg = mol%chrg + 2
   if (mod(wf_m%wfn%nel,2).ne.0) wf_m%wfn%nopen = 1
   call calc%singlepoint(env, mol, wf_m, 1, exist, etot2, g, sigma, egap, res)
   fukui(2,:) = chk%wfn%q-wf_m%wfn%q
   fukui(3,:) = 0.5d0*(wf_p%wfn%q-wf_m%wfn%q)
   write(env%unit,'(a)')
   write(env%unit, '(1x,"    #        f(+)     f(-)     f(0)")')
   do i=1,mol%n
      write(env%unit,'(i6,a4,2f9.3,2f9.3,2f9.3)') i, mol%sym(i), fukui(1,i), fukui(2,i), fukui(3,i)
   enddo
   mol%chrg = mol%chrg - 1

end subroutine vfukui

end module xtb_vertical
