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

!> TODO
module xtb_xtb_gfn0
   use xtb_mctc_accuracy, only : wp
   use xtb_param_atomicrad, only : atomicRad
   use xtb_xtb_data
   use xtb_aoparam
   implicit none
   private

   public :: initGFN0


   interface initGFN0
      module procedure :: initData
      module procedure :: initRepulsion
      module procedure :: initCoulomb
      module procedure :: initHamiltonian
   end interface initGFN0


   !> Maximum number of elements supported by GFN0-xTB
   integer, parameter :: maxElem = 86


contains


subroutine initData(self)

   !>
   type(TxTBData), intent(out) :: self

   self%nShell = ao_n(:maxElem)

   call initGFN0(self%repulsion)
   call initGFN0(self%coulomb, self%nShell)
   call initGFN0(self%hamiltonian, self%nShell)

end subroutine initData


subroutine initRepulsion(self)

   !>
   type(TRepulsionData), intent(out) :: self

end subroutine initRepulsion


subroutine initCoulomb(self, nShell)

   !>
   type(TCoulombData), intent(out) :: self

   !>
   integer, intent(in) :: nShell(:)

end subroutine initCoulomb


subroutine initHamiltonian(self, nShell)

   !>
   type(THamiltonianData), intent(out) :: self

   !>
   integer, intent(in) :: nShell(:)

   integer :: mShell

   mShell = maxval(nShell)
   self%angShell = ao_l(:mShell, :maxElem)

   self%electronegativity = en(:maxElem)
   self%atomicRad = atomicRad(:maxElem)
   self%shellPoly = polyr(:, :maxElem)
   self%pairParam = kpair(:maxElem, :maxElem)

end subroutine initHamiltonian


end module xtb_xtb_gfn0
