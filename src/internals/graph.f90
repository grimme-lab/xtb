! This file is part of xtb.
!
! Copyright (C) 2026 Leopold M. Seidler
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

!> For building a molecular graph, starting from a neighbour list.
!>
!> Every atom of the central cell is a node. Neighbours in periodic images
!> are folded onto the atom of the central cell and duplicate entries are
!> removed, so the adjacency is symmetric and every edge is stored once for
!> both of its atoms. The BFS then assigns each atom a parent, giving a
!> spanning tree per connected component.
module xtb_internals_graph
   use xtb_type_neighbourlist, only : TNeighbourList
   implicit none
   private

   public :: graph_type, init

   !> Molecular graph of a system
   type :: graph_type

      !> Number of atoms in the system
      integer :: n = 0

      !> Number of neighbours of each atom, the filled entries of i_neigh(:, iat)
      integer, allocatable :: neighs(:)

      !> Adjacency list, i_neigh(1:neighs(iat), iat) holds the neighbours of iat
      integer, allocatable :: i_neigh(:, :)

      !> Parent of each atom in the BFS spanning tree, zero for a component root
      integer, allocatable :: parent(:)

      !> Atoms in BFS order, order(1:n), each component starting at its root
      integer, allocatable :: order(:)

   end type graph_type

   !> Build the molecular graph from a neighbour list
   interface init
      module procedure :: new_graph
   end interface init

contains

!> Build the molecular graph from a generated neighbour list.
!>
!> The list has to be generated (and, for periodic systems, updated) by the
!> caller.  Only the atoms of the central cell are nodes and only the graph is
!> built here, the list is not stored.
pure subroutine new_graph(self, neigh_list)
   !> Graph to fill, its allocatable components are allocated here
   type(graph_type), intent(out) :: self
   !> Neighbour list of the system
   type(TNeighbourList), intent(in) :: neigh_list

   integer :: nat, iat, jat, j, k, m_neigh, head, tail, root
   ! Stamp of the last atom whose list was scanned, to drop duplicate entries
   integer, allocatable :: seen(:)
   ! Next free entry in i_neigh(:, i_at)
   integer, allocatable :: pos(:)
   logical, allocatable :: visited(:)

   nat = size(neigh_list%neighs)
   self%n = nat
   allocate(self%neighs(nat), source = 0)

   if (nat == 0) then
      allocate(self%i_neigh(1, 0), self%parent(0), self%order(0))
      return
   end if

   ! count the neighbours of every atom
   allocate(seen(nat), source = 0)
   do iat = 1, nat
      do j = 1, neigh_list%neighs(iat)
         jat = neigh_list%image(neigh_list%iNeigh(j, iat))
         if (jat <= iat) cycle
         if (seen(jat) == iat) cycle
         seen(jat) = iat
         self%neighs(iat) = self%neighs(iat) + 1
         self%neighs(jat) = self%neighs(jat) + 1
      end do
   end do

   ! fill the adjacency list
   m_neigh = max(maxval(self%neighs), 1)
   allocate(self%i_neigh(m_neigh, nat), source = 0)
   allocate(pos(nat), source = 0)
   seen(:) = 0
   do iat = 1, nat
      do j = 1, neigh_list%neighs(iat)
         jat = neigh_list%image(neigh_list%iNeigh(j, iat))
         if (jat <= iat) cycle
         if (seen(jat) == iat) cycle
         seen(jat) = iat
         pos(iat) = pos(iat) + 1
         self%i_neigh(pos(iat), iat) = jat
         pos(jat) = pos(jat) + 1
         self%i_neigh(pos(jat), jat) = iat
      end do
   end do

   ! run BFS spanning forest
   allocate(self%parent(nat), self%order(nat), source = 0)
   allocate(visited(nat), source = .false.)
   tail = 0
   do root = 1, nat
      if (visited(root)) cycle
      visited(root) = .true.
      tail = tail + 1
      self%order(tail) = root
      head = tail
      do while (head <= tail)
         iat = self%order(head)
         head = head + 1
         do k = 1, self%neighs(iat)
            jat = self%i_neigh(k, iat)
            if (visited(jat)) cycle
            visited(jat) = .true.
            self%parent(jat) = iat
            tail = tail + 1
            self%order(tail) = jat
         end do
      end do
   end do

end subroutine new_graph

end module xtb_internals_graph
