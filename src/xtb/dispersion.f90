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

!> Self-consistent dispersion energy
module xtb_xtb_dispersion
   use xtb_mctc_accuracy, only : wp
   use xtb_disp_dftd4
   use xtb_disp_ncoord
   use xtb_type_molecule, only : TMolecule
   use xtb_type_dispersionmodel, only : TDispersionModel
   use xtb_xtb_data, only : TDispersionData
   implicit none
   private

   public :: TxTBDispersion, init


   type :: TxTBDispersion

      integer :: ndim

      real(wp) :: g_a

      real(wp) :: g_c

      real(wp) :: wf

      integer, allocatable :: itbl(:)

      real(wp), allocatable :: dispMat(:, :)

      real(wp), allocatable :: gweight(:)

      real(wp), allocatable :: refC6(:, :)

      type(TDispersionModel), pointer :: dispm => null()

   contains

      procedure :: addShift

      procedure :: getEnergy

   end type TxTBDispersion


   interface init
      module procedure :: initDispersion
   end interface init


contains


subroutine initDispersion(self, input, mol)

   type(TxTBDispersion), intent(out) :: self

   type(TDispersionData), intent(in), target :: input

   type(TMolecule), intent(in) :: mol

   real(wp), allocatable :: cn(:)

   integer :: k, iat
   integer, allocatable :: itbl(:)

   self%g_a = input%g_a
   self%g_c = input%g_c
   self%wf = input%wf
   self%dispm => input%dispm
   call d4dim(self%dispm, mol%n, mol%at, self%ndim)

   allocate(self%gweight(self%ndim))
   allocate(self%refC6(self%ndim, self%ndim))
   allocate(self%dispmat(self%ndim, self%ndim))
   allocate(self%itbl(mol%n))
   allocate(cn(mol%n))

   k = 0
   do iat = 1, mol%n
      self%itbl(iat) = k
      k = k + self%dispm%nref(mol%at(iat))
   enddo

   call ncoord_d4(mol%n, mol%at, mol%xyz, cn, thr=1600.0_wp)
   call d4(self%dispm, mol%n, self%ndim, mol%at, self%wf, self%g_a, &
      & self%g_c, cn, self%gweight, self%refC6)
   call build_wdispmat(self%dispm, mol%n, self%ndim, mol%at, self%itbl, mol%xyz, &
      & input%dpar, self%refC6, self%gweight, self%dispmat)

end subroutine initDispersion


subroutine addShift(self, id, qat, atomicShift)

   class(TxTBDispersion), intent(inout) :: self

   integer, intent(in) :: id(:)

   real(wp), intent(in) :: qat(:)

   real(wp), intent(inout) :: atomicShift(:)

   call disppot(self%dispm, size(id), self%ndim, id, self%itbl, qat, self%g_a, self%g_c, &
      & self%dispMat, self%gweight, atomicShift)

end subroutine addShift


subroutine getEnergy(self, id, qat, energy)

   class(TxTBDispersion), intent(inout) :: self

   integer, intent(in) :: id(:)

   real(wp), intent(in) :: qat(:)

   real(wp), intent(out) :: energy

   energy = edisp_scc(self%dispm, size(id), self%ndim, id, self%itbl, qat, self%g_a, &
      & self%g_c, self%dispMat, self%gweight)

end subroutine getEnergy


end module xtb_xtb_dispersion
