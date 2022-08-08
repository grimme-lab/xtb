! This file is part of xtb.
!
! Copyright (C) 2022 Stefan Grimme
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

!> Intermolecular force field calculator
module xtb_iff_calculator
   use xtb_mctc_accuracy, only: wp
   use xtb_type_calculator, only: TCalculator
   use xtb_type_environment, only: TEnvironment
   use xtb_type_molecule, only: TMolecule, init
   use xtb_type_restart
   use xtb_iff_data, only: TIFFData
   use xtb_type_data, only: scc_results
   use xtb_iff_iffenergy, only: iff_e
   use xtb_docking_param
   use xtb_iff_iffini, only : init_iff

   implicit none

   private

   public :: TIFFCalculator, newIFFCalculator

   !> Calculator interface for xTB based methods
   type, extends(TCalculator) :: TIFFCalculator

      type(TIFFData) :: dat

   contains

      !> Perform xTB single point calculationV
      procedure :: singlepoint

      !> Write informative printout
      procedure :: writeInfo

   end type TIFFCalculator

   character(len=*), private, parameter :: outfmt = &
                                           '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'

contains

   subroutine newIFFCalculator(env, comb, iff_data, calc)

      character(len=*), parameter :: source = 'main_setup_newIFFCalculator'

      type(TEnvironment), intent(inout) :: env

      !> Combined structure of molA and molB (molA has to be first)
      type(TMolecule), intent(in) :: comb

      type(TIFFData), intent(in) :: iff_data

      type(TIFFCalculator), intent(out) :: calc

      integer, allocatable :: at(:)
      real(wp), allocatable :: xyz(:,:)
      logical :: exitRun
      real(wp) :: molA_e,molB_e
      character(len=:), allocatable :: fnam
      real(wp) :: icoord(6), icoord0(6)
      real(wp) :: rlmo2(4,(comb%n-natom_molA)*10)

      call set_iff_param
      call calc%dat%allocateIFFData(natom_molA,comb%n-natom_molA)
      calc%dat = iff_data

      xyz=calc%dat%xyz2
      rlmo2=calc%dat%rlmo2

      call init_iff(env,calc%dat%n1,calc%dat%n2,calc%dat%at1,calc%dat%at2,&
        &          calc%dat%neigh,calc%dat%xyz1,calc%dat%xyz2,calc%dat%q1,&
        &          calc%dat%q2,calc%dat%c6ab,calc%dat%z1,calc%dat%z2,&
        &          calc%dat%cprob,calc%dat%nlmo1,calc%dat%nlmo2,calc%dat%lmo1,calc%dat%lmo2,&
        &          calc%dat%qdr1,calc%dat%qdr2,calc%dat%rlmo1,calc%dat%rlmo2,&
        &          calc%dat%cn1,calc%dat%cn2,calc%dat%alp1,calc%dat%alp2,calc%dat%alpab,&
        &          calc%dat%den1,calc%dat%den2,calc%dat%gab1,calc%dat%gab2,calc%dat%qcm1,&
        &          calc%dat%qcm2,calc%dat%n,calc%dat%at,calc%dat%xyz,calc%dat%q,icoord,icoord0,&
        &          .false.)
      calc%dat%xyz2=xyz
      calc%dat%rlmo2=rlmo2

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not create IFF calculator", source)
         return
      end if

   end subroutine newIFFCalculator

   subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
         & energy, gradient, sigma, hlgap, results)

      !> Source of the generated errors
      character(len=*), parameter :: source = 'iff_calculator_singlepoint'

      !> Calculator instance
      class(TIFFCalculator), intent(inout) :: self

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

      call iff_e(env, self%dat%n, self%dat%n1, self%dat%n2, self%dat%at1,&
           &     self%dat%at2, self%dat%neigh, self%dat%xyz1, self%dat%xyz2,&
           &     self%dat%q1, self%dat%q2, self%dat%c6ab, self%dat%z1, self%dat%z2,&
           &     self%dat%nlmo1, self%dat%nlmo2, self%dat%lmo1, self%dat%lmo2,&
           &     self%dat%rlmo1, self%dat%rlmo2,self%dat%qdr1, self%dat%qdr2,&
           &     self%dat%cn1, self%dat%cn2, self%dat%alp1,&
           &     self%dat%alp2, self%dat%alpab, self%dat%qct1, self%dat%qct2,&
           &     self%dat%den1, self%dat%den2, self%dat%gab1, self%dat%gab2,&
           &     .true., 0, energy)

      results%e_total=energy
      gradient = 0
      sigma = 0
      hlgap = 0

   end subroutine singlepoint

   subroutine writeInfo(self, unit, mol)

      !> Calculator instance
      class(TIFFCalculator), intent(in) :: self

      !> Unit for I/O
      integer, intent(in) :: unit

      !> Molecular structure data
      type(TMolecule), intent(in) :: mol

   end subroutine writeInfo

end module xtb_iff_calculator
