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

   ! ========================================================================
   ! HAMILTONIAN DATA
   !>
   integer, parameter :: valenceShell(3, 1:maxElem) = reshape([&
      & 1, 0, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0, &
      & 1, 1, 0,  1, 1, 0,  1, 1, 1,  1, 1, 0,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 0,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 0,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 0,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 0,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 0,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1,  1, 1, 1, &
      & 1, 1, 1,  1, 1, 1,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0,  1, 1, 0, &
      & 1, 1, 1,  1, 1, 1], shape(valenceShell))

   !>
   real(wp), parameter :: referenceOcc(3, 1:maxElem) = reshape([&
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp, &
      & 2.0_wp, 3.0_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp, &
      & 2.0_wp, 6.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 2.0_wp,  2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp, &
      & 2.0_wp, 0.0_wp, 5.0_wp,  2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp, &
      & 2.0_wp, 0.0_wp, 8.0_wp,  2.0_wp, 0.0_wp, 9.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 2.0_wp,  2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp, &
      & 2.0_wp, 0.0_wp, 5.0_wp,  2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp, &
      & 2.0_wp, 0.0_wp, 8.0_wp,  2.0_wp, 0.0_wp, 9.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 2.0_wp, &
      & 2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp,  2.0_wp, 0.0_wp, 5.0_wp, &
      & 2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp,  2.0_wp, 0.0_wp, 8.0_wp, &
      & 2.0_wp, 0.0_wp, 9.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp, &
      & 2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp, &
      & 2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp], shape(referenceOcc))


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

   integer :: mShell, nPrim
   integer :: iZp, iSh

   mShell = maxval(nShell)
   self%angShell = ao_l(:mShell, :maxElem)

   self%electronegativity = en(:maxElem)
   self%atomicRad = atomicRad(:maxElem)
   self%shellPoly = polyr(:, :maxElem)
   self%pairParam = kpair(:maxElem, :maxElem)
   self%referenceOcc = referenceOcc
   self%selfEnergy = ao_lev(:mShell, :maxElem)
   self%slaterExponent = ao_exp(:mShell, :maxElem)
   self%valenceShell = valenceShell(:mShell, :maxElem)
   self%principalQuantumNumber = ao_pqn(:mShell, :maxElem)

   allocate(self%numberOfPrimitives(mShell, maxElem))
   do iZp = 1, maxElem
      do iSh = 1, nShell(iZp)
         nPrim = 0
         if (iZp <= 2) then
            select case(self%angShell(iSh, iZp))
            case(0)
               if (self%valenceShell(iSh, iZp) /= 0) then
                  nPrim = 3
               else
                  nPrim = 2
               end if
            case(1)
               nPrim = 3
            end select
         else
            select case(self%angShell(iSh, iZp))
            case(0)
               if (self%principalQuantumNumber(iSh, iZp) > 5) then
                  nPrim = 6
               else
                  nPrim = 4
               end if
            case(1)
               if (self%principalQuantumNumber(iSh, iZp) > 5) then
                  nPrim = 6
               else
                  nPrim = 3
               end if
            case(2, 3)
               nPrim = 4
            end select
         end if
         self%numberOfPrimitives(iSh, iZp) = nPrim
      end do
   end do

end subroutine initHamiltonian


end module xtb_xtb_gfn0
