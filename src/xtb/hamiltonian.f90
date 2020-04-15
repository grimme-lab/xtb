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

!> Implementation of the xTB core Hamiltonian
module xtb_xtb_hamiltonian
   use xtb_mctc_accuracy, only : wp
   use xtb_xtb_data, only : THamiltonianData
   implicit none
   private

   public :: getSelfEnergy


   interface getSelfEnergy
      module procedure :: getSelfEnergyFlat
      module procedure :: getSelfEnergy2D
   end interface getSelfEnergy


contains


subroutine getSelfEnergyFlat(hData, nShell, at, cn, qat, selfEnergy, dSEdcn, dSEdq)
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in) :: nShell(:)
   integer, intent(in) :: at(:)
   real(wp), intent(in), optional :: cn(:)
   real(wp), intent(in), optional :: qat(:)
   real(wp), intent(out) :: selfEnergy(:)
   real(wp), intent(out), optional :: dSEdcn(:)
   real(wp), intent(out), optional :: dSEdq(:)

   integer :: ind, iAt, iZp, iSh

   selfEnergy(:) = 0.0_wp
   if (present(dSEdcn)) dSEdcn(:) = 0.0_wp
   if (present(dSEdq)) dSEdq(:) = 0.0_wp
   ind = 0
   do iAt = 1, size(cn)
      iZp = at(iAt)
      do iSh = 1, nShell(iZp)
         selfEnergy(ind+iSh) = hData%selfEnergy(iSh, iZp)
      end do
      ind = ind + nShell(iZp)
   end do
   if (present(dSEdcn) .and. present(cn)) then
      ind = 0
      do iAt = 1, size(cn)
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(ind+iSh) = selfEnergy(ind+iSh) &
               & - hData%kCN(iSh, iZp) * cn(iAt)
            dSEdcn(ind+iSh) = -hData%kCN(iSh, iZp)
         end do
         ind = ind + nShell(iZp)
      end do
   end if
   if (present(dSEdq) .and. present(qat)) then
      ind = 0
      do iAt = 1, size(cn)
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(ind+iSh) = selfEnergy(ind+iSh) &
               & - hData%kQShell(iSh,iZp)*qat(iAt) - hData%kQAtom(iZp)*qat(iAt)**2
            dSEdq(ind+iSh) = -hData%kQShell(iSh,iZp) - hData%kQAtom(iZp)*2*qat(iAt)
         end do
         ind = ind + nShell(iZp)
      end do
   end if

end subroutine getSelfEnergyFlat


subroutine getSelfEnergy2D(hData, nShell, at, cn, qat, selfEnergy, dSEdcn, dSEdq)
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in) :: nShell(:)
   integer, intent(in) :: at(:)
   real(wp), intent(in), optional :: cn(:)
   real(wp), intent(in), optional :: qat(:)
   real(wp), intent(out) :: selfEnergy(:, :)
   real(wp), intent(out), optional :: dSEdcn(:, :)
   real(wp), intent(out), optional :: dSEdq(:, :)

   integer :: iAt, iZp, iSh

   selfEnergy(:, :) = 0.0_wp
   if (present(dSEdcn)) dSEdcn(:, :) = 0.0_wp
   if (present(dSEdq)) dSEdq(:, :) = 0.0_wp
   do iAt = 1, size(cn)
      iZp = at(iAt)
      do iSh = 1, nShell(iZp)
         selfEnergy(iSh, iAt) = hData%selfEnergy(iSh, iZp)
      end do
   end do
   if (present(dSEdcn) .and. present(cn)) then
      do iAt = 1, size(cn)
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(iSh, iAt) = selfEnergy(iSh, iAt) &
               & - hData%kCN(iSh, iZp) * cn(iAt)
            dSEdcn(iSh, iAt) = -hData%kCN(iSh, iZp)
         end do
      end do
   end if
   if (present(dSEdq) .and. present(qat)) then
      do iAt = 1, size(cn)
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(iSh, iAt) = selfEnergy(iSh, iAt) &
               & - hData%kQShell(iSh,iZp)*qat(iAt) - hData%kQAtom(iZp)*qat(iAt)**2
            dSEdq(iSh, iAt) = -hData%kQShell(iSh,iZp) &
               & - hData%kQAtom(iZp)*2*qat(iAt)
         end do
      end do
   end if

end subroutine getSelfEnergy2D


end module xtb_xtb_hamiltonian
