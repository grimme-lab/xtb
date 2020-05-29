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

!> Compare and generate topological information from wavefunctions
module xtb_topology
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_filetypes, only : fileType, generateFileName
   use xtb_io_writer, only : writeMolecule
   use xtb_type_molecule, only : TMolecule, len
   use xtb_type_neighbourlist, only : TNeighbourList, init, resizeNeigh
   use xtb_type_topology, only : TTopology, len
   implicit none
   private

   public :: checkTopology
   public :: makeBondTopology, compareBondTopology


contains


subroutine checkTopology(unit, mol, wbo, verbosity)

   integer, intent(in) :: unit

   type(TMolecule), intent(in) :: mol

   real(wp), intent(in) :: wbo(:, :)

   integer, intent(in) :: verbosity

   character(len=:), allocatable :: fname
   type(TMolecule) :: copy
   integer :: ich, ftype
   logical :: match

   ! CT files have some limits:
   ! bail out if we encounter a structure not suitable for this format
   if (mol%npbc /= 0 .or. mol%n > 999) return

   copy = mol
   call makeBondTopology(copy%bonds, mol, wbo)

   call compareBondTopology(unit, mol%bonds, copy%bonds, mol, verbosity, match)

   if (.not.match .and. len(copy%bonds) > 0) then
      ! prefer molfiles, except the input has been SDF
      if (mol%ftype /= fileType%sdf) then
         ftype = fileType%molfile
      else
         ftype = mol%ftype
      end if
      call generateFileName(fname, 'xtbtopo', '', ftype)
      call open_file(ich, fname, 'w')
      if (verbosity > 0 .and. len(mol%bonds) > 0) then
         write(unit, '("Writing corrected topology to", 1x, a)') fname
      else
         write(unit, '("Writing topology from bond orders to", 1x, a)') fname
      end if
      call writeMolecule(copy, ich, format=ftype)
      call close_file(ich)
   end if

end subroutine checkTopology


subroutine makeBondTopology(topo, mol, wbo)

   type(TTopology), intent(out) :: topo

   type(TMolecule), intent(in) :: mol

   real(wp), intent(in) :: wbo(:, :)

   real(wp), parameter :: wbothr = 0.3_wp
   integer :: iat, jat

   call topo%allocate(size=len(mol), order=3)
   do iat = 1, len(mol)
      do jat = 1, iat-1
         associate(w => wbo(jat, iat))
            if (any(mol%at(iat) == [6, 7, 8]) .and. any(mol%at(jat) == [6, 7, 8])) then
               if (w > wbothr .and. w <= 1.2_wp) then
                  call topo%push_back([jat, iat, 1])
               else if (w > 1.5_wp .and. w <= 2.3_wp) then
                  call topo%push_back([jat, iat, 2])
               else if (w > 2.3_wp) then
                  call topo%push_back([jat, iat, 3])
               else if (w > 1.2_wp .and. w <= 1.5_wp) then
                  call topo%push_back([jat, iat, 4])
               end if
            else
               if (w > wbothr .and. w <= 1.3_wp) then
                  call topo%push_back([jat, iat, 1])
               else if (w > 1.3_wp .and. w <= 2.3_wp) then
                  call topo%push_back([jat, iat, 2])
               else if (w > 2.3_wp) then
                  call topo%push_back([jat, iat, 3])
               end if
            end if
         end associate
      end do
   end do
end subroutine makeBondTopology


subroutine compareBondTopology(unit, topo1, topo2, mol, verbosity, match)

   integer, intent(in) :: unit

   type(TTopology), intent(in) :: topo1

   type(TTopology), intent(in) :: topo2

   type(TMolecule), intent(in) :: mol

   integer, intent(in) :: verbosity

   logical, intent(out) :: match

   type(TNeighbourList) :: neighList1, neighList2
   integer :: iat, neigh

   match = len(topo1) == len(topo2)

   if (.not.match .and. verbosity > 0) then
      write(unit, '("Topologies differ in total number of bonds")')
   end if

   if (match .or. verbosity > 1) then
      call topologyToNeighbourList(topo1, neighList1, mol)
      call topologyToNeighbourList(topo2, neighList2, mol)

      match = match .and. all(neighList1%neighs == neighList2%neighs)
      if (.not.match .and. verbosity > 0) then
         write(unit, '("Topologies differ in number of bonds")')
      end if

      if (match .or. verbosity > 1) then
         neigh = min(ubound(neighList1%iNeigh, 1), ubound(neighList2%iNeigh, 1))
         match = match .and. all(neighList1%iNeigh(:neigh, :) &
            & == neighList2%iNeigh(:neigh, :))
         if (.not.match .and. verbosity > 0) then
            write(unit, '("Topologies differ in connectivity")')
         end if

         if (match .or. verbosity > 1) then
            match = match .and. all(nint(neighList1%weight(:neigh, :)) &
               & == nint(neighList2%weight(:neigh, :)))
            if (.not.match .and. verbosity > 0) then
               write(unit, '("Topologies differ in bond orders")')
            end if
         end if
      end if
   end if

end subroutine compareBondTopology


subroutine topologyToNeighbourList(topo, neighList, mol)

   type(TTopology), intent(in) :: topo

   type(TMolecule), intent(in) :: mol

   type(TNeighbourList), intent(out) :: neighList

   integer :: iBond(3)
   integer :: mNeigh
   integer :: ii, iat, jat
   real(wp) :: r2

   call init(neighList, len(mol))

   mNeigh = ubound(neighList%iNeigh, dim=1)

   neighList%neighs(:) = 0
   neighList%iNeigh(:, :) = 0
   neighList%dist2(:, :) = 0.0_wp
   neighList%weight(:, :) = 0.0_wp
   do iat = 1, len(mol)
      neighList%iNeigh(0, iat) = iat
      neighList%weight(0, iat) = 1.0_wp
      neighList%coords(:, iat) = mol%xyz(:, iat)
      neighList%image(iat) = iat
      neighList%trans(iat) = 1
   end do

   do ii = 1, len(topo)
      call topo%get_item(ii, iBond)
      iat = min(iBond(1), iBond(2))
      jat = max(iBond(1), iBond(2))
      r2 = sum((neighList%coords(:, iat) - neighList%coords(:, jat))**2)
      neighList%neighs(iat) = neighList%neighs(iat) + 1
      if (neighList%neighs(iat) > mNeigh) then
         mNeigh = 2*mNeigh
         call resizeNeigh(mNeigh, neighList%iNeigh, neighList%dist2, &
            & neighList%weight)
      end if
      neighList%iNeigh(neighList%neighs(iat), iat) = jat
      neighList%dist2(neighList%neighs(iat), iat) = r2
      neighList%weight(neighList%neighs(iat), iat) = iBond(3)
   end do

   call neighList%sort

end subroutine topologyToNeighbourList


end module xtb_topology
