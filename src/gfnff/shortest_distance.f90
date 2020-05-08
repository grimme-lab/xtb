! This file is part of xtb.
!
! Copyright (C) 2019-2020 Stefan Grimme
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calculates the shortest path between two atoms
!start and goal are integers, determining the index in xyz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*4 function shortest_distance(nspin, start, goal, neighbours, input_distances, visited, precessor)
    implicit none
    !Dummy Arguments:
    integer, intent(in) :: nspin
    integer, intent(in) :: start
    integer, intent(in) :: goal
    integer, intent(in) :: neighbours(20,nspin)
    real*4,  intent(in) :: input_distances(nspin*(nspin+1)/2)
    logical, intent(out):: visited(nspin)
    integer, intent(out):: precessor(nspin)

    integer lin

    !Stack:
    integer:: current
    integer:: neighbour
    integer:: i_neighbours
    real*4:: alt_dist
    real*4:: distance(nspin)

    !logical:: visited(nspin)

    precessor = 0
    distance = huge(distance)
    distance(start) = 0
    visited = .false.
    do while (all(visited) == .false.) !as long there are unvisited nodes
        current = minloc(distance,1, .not. visited)

        !Abort if Fragments are not connected:
        if (distance(current) == huge(distance)) then
            shortest_distance = huge(distance)
            return
        end if ! distance(current) == huge(distance)

        visited(current) = .true.
        if (current == goal) then
            shortest_distance = distance(goal)
            !route to target found
            do while (precessor(current) /= 0)
                current = precessor(current)
            end do ! End loop: while precessor(current) /= 0
            exit
        else
            !loop over all neighbours of current atom
            do i_neighbours = 1, neighbours(20,current)
                neighbour = neighbours(i_neighbours,current)
                if (visited(neighbour) == .false.) then
                    !distanzupdate
                    alt_dist = distance(current) + input_distances(lin(current,neighbour))
                    if ((alt_dist < distance(neighbour))) then
                        distance(neighbour) = alt_dist
                        precessor(neighbour) = current
                    end if ! (alt_dist < distance(neighbour))
                endif !(visited(neighbour))
            end do ! End Loop over  i_neighbours from 1 to all_Atoms(current)%n_neighbours
        endif
    end do ! End loop: while sum(visited) /= 0
        !initialize
    shortest_distance = distance(goal)
end function shortest_distance
