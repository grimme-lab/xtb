! This file is part of dipro.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Fragmentation based on bond orders found for a system
module xtb_dipro_fragment
   use mctc_env, only : wp
   implicit none
   private

   public :: get_wiberg_fragment

contains

!> Obtain fragment information from Wiberg/Mayer type bond orders
subroutine get_wiberg_fragment(fragment,wbo,thr)
   !> Fragment information
   integer, intent(out) :: fragment(:)
   !> Bond orders
   real(wp),intent(in) :: wbo(:, :)
   !> Threshold to count bond orders as actual bond
   real(wp),intent(in) :: thr

   real(wp),allocatable :: bond(:, :)
   logical, allocatable :: taken(:)
   integer :: i, nfrag, nat
   real(wp) :: xsum

   nat = size(fragment)
   fragment(:) = 0
   allocate(taken(nat), source=.false.)
   allocate(bond(nat, nat), source=0.0_wp)
   where(wbo > thr)
      bond = min(wbo, 1.0_wp) !filtert kleinsten Wert aus wbo und setzt den an die richtige 
!nat,nat stell in wbo, falls wbo < 1.0 wird =1 gesetzt
   elsewhere
      bond = 0.0_wp !sollte in einer symmetrischen Matrix mit 0/1/>1 resultieren
   endwhere

   nfrag = 0
   do i = 1, nat
      if(taken(i)) cycle
      nfrag=nfrag+1
      fragment(i)=nfrag
      taken(i)=.true.
      call find_neighbours(i, sum(ceiling(bond(:,:)), 1), taken, bond, fragment, nfrag)
   end do !ceiling(-63.2234)=-63, ceiling(63.2234)=63 . Aufrunden-FUnktion
end subroutine get_wiberg_fragment

!> Worker routine to transverse the graph recursively
recursive subroutine find_neighbours(i, cn, taken, bond, fragment, this_fragment)
   !> Current node in the graph
   integer, intent(in) :: i
   !> Number of edges for every node
   integer, intent(in) :: cn(:)
   !> Status of nodes
   logical, intent(inout) :: taken(:)
   !> Weights of the graph edges, destroyed while transversing the graph
   real(wp), intent(inout) :: bond(:,:)
   !> Fragments found
   integer, intent(inout) :: fragment(:)
   !> Index of the current fragmant
   integer, intent(inout) :: this_fragment

   integer :: j, k

   do k = 1, cn(i)
      j = maxloc(bond(:, i), 1)
      bond(j, i) = 0
      if (i == j .or. taken(j)) cycle
      fragment(j) = this_fragment
      taken(j) = .true.
      call find_neighbours(j, cn, taken, bond, fragment, this_fragment)
   end do
end subroutine find_neighbours

end module xtb_dipro_fragment
