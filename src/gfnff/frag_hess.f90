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

module xtb_gfnff_fraghess
!$ use omp_lib
   use xtb_gfnff_mrec, only : mrecgff, mrecgffPBC
   implicit none

   interface
     function shortest_distance(nspin, start, goal, numnb, neighbours, input_distances, visited, precessor)
        use xtb_mctc_accuracy, only : wp
        implicit none
        real(wp)             :: shortest_distance
        integer, intent(in)  :: nspin
        integer, intent(in)  :: start
        integer, intent(in)  :: goal
        integer, intent(in)  :: numnb
        integer, intent(in)  :: neighbours(numnb, nspin, 1)
        real(wp), intent(in) :: input_distances(nspin, nspin)
        logical, intent(out) :: visited(nspin)
        integer, intent(out) :: precessor(nspin)
      end function shortest_distance
      subroutine eigsort4(lab,u,ew)
         use xtb_mctc_accuracy, only : sp
         implicit none
         integer,  intent(in)    :: lab
         real(sp), intent(inout) :: u(lab,lab)
         real(sp), intent(inout) :: ew(lab)
      end subroutine eigsort4
    end interface

   contains

     subroutine fragmentize(nspin, at, xyz, maxsystem, maxmagnat, jab, numnb, numctr, neigh, ispinsyst, nspinsyst, nsystem, env)
        use xtb_mctc_accuracy, only : wp, sp
        use xtb_type_environment, only: TEnvironment
        implicit none

        !Dummy Arguments:
        integer,  intent(in)  :: nspin                              ! # of atoms in the whole system
        integer,  intent(in)  :: maxsystem                          ! maximum # of fragments
        integer,  intent(in)  :: maxmagnat                          ! maximum # of atoms per fragment
        integer,  intent(in)  :: numnb, numctr
        integer,  intent(in)  :: neigh(numnb, nspin, numctr)        ! neighbour list
        integer,  intent(in)  :: at(nspin)                          ! atom nunber
        real(wp), intent(in)  :: xyz(3,nspin)                       ! xyz coordinates
        real(wp), intent(in)  :: jab(nspin*(nspin + 1)/2)           ! distances between pairs A and B
        integer, allocatable, intent(out) :: ispinsyst(:,:)         ! array with list of atoms of each fragment
        integer, allocatable, intent(out) :: nspinsyst(:)           ! array with # of atoms for each fragment
        integer,  intent(out) :: nsystem                            ! # of fragments
        real(sp)  :: rmaxab(nspin, nspin)

        !Stack
        integer  :: i, ati
        integer  :: current
        integer  :: ass
        integer  :: k
        integer  :: j, atj
        integer  :: lin
        integer  :: fragcount
        integer  :: nfrag_ini
        integer  :: nci_ini
        integer  :: eq_frag
        integer  :: maxdistatoms(2)
        integer  :: max_linkatoms(2)
        integer  :: assigned_to_frag(nspin)
        integer  :: precessor(nspin)
        integer  :: path(nspin)
        integer  :: fragvec(nspin)
        integer  :: nbox
        integer  :: nci_frag_size
        real(wp) :: grid(3,maxsystem)
        real(wp) :: magdist(nspin, nspin)
        real(wp) :: shortest_distance
        real(wp) :: maxdist
        real(wp) :: cur_dist
        real(wp) :: max_link
        real(wp) :: box_xyz(3), vec(3)
        real(wp) :: frag_cma(3,maxsystem)
        logical  :: equal(maxsystem)
        logical  :: visited(nspin)
        logical  :: assigned(nspin, nspin)

        integer, allocatable  :: ifrag_ini(:)

        type(TEnvironment), intent(inout) :: env

        allocate( ispinsyst(nspin,maxsystem), source = 0)
        allocate( nspinsyst(maxsystem), source = 0)

        if (nspin.lt.500) then
           nsystem = 1
           return
        end if

        nci_frag_size = 50
        fragcount = 0
        fragvec = 0
        eq_frag = 0
        frag_cma = 0
        grid = 0
        equal = .false.

        call mrecgffPBC(nspin,numctr, numnb, neigh, fragcount, fragvec)

        nsystem = maxval(fragvec)

        !Use NCI-Fragmentation from mrec:
        do i = 1, nspin
           nspinsyst(fragvec(i)) = nspinsyst(fragvec(i)) + 1
           ispinsyst(nspinsyst(fragvec(i)), fragvec(i)) = i
        end do ! End Loop over  i from 1 to nspin

        !Calculate COM for each NCI fragment
        do i = 1, nsystem
           nfrag_ini = nspinsyst(i)
           if (nfrag_ini.lt.nci_frag_size) then
             allocate (ifrag_ini(nfrag_ini), source = 0)
             ifrag_ini = ispinsyst(1:nfrag_ini,i)
             call com(nfrag_ini,at(ifrag_ini),xyz(:,ifrag_ini),frag_cma(:,i))
             deallocate (ifrag_ini)
           end if
        end do

        !Determine box size
        do i = 1, 3
           box_xyz(i) = maxval(frag_cma(i,:)) - minval(frag_cma(i,:))
        end do

        !Determine number of sub boxes
        nbox =  ceiling((real(nspin)/real(maxmagnat))**(1.0d0/3.0d0))

        !Determine position of nci frag on grid
        do i = 1, nsystem
           nfrag_ini = nspinsyst(i)
           if (nfrag_ini.lt.nci_frag_size) then
             do j = 1, 3
                grid(j,i) =   nbox * ( frag_cma(j,i) - minval(frag_cma(j,:)) ) / box_xyz(j)
                if (grid(j,i).eq.0) grid(j,i) = grid(j,i) + 1.0d0
                grid(j,i) = ceiling( grid(j,i) )
             end do
           end if
        end do

        !Sum over all nci fragments on same grid position
        nci_ini=nsystem
        do i = 1, nci_ini !nsystem
           eq_frag=0
           !nfrag_ini=nspinsyst(i)
           do j = i+1, nci_ini
              nfrag_ini=nspinsyst(j)
              vec(:) = grid(:,i)-grid(:,j)
              if (norm2(grid(:,i)).ne.0.and.norm2(vec).eq.0.and..not.equal(j)) then
                eq_frag = eq_frag + 1
                do k = 1, nfrag_ini
                   ispinsyst(nspinsyst(i)+k,i) = ispinsyst(k,j)
                end do
                nspinsyst(i) = nspinsyst(i) + nspinsyst(j)
                equal(j) = .true.
              end if
           end do
           nsystem = nsystem - eq_frag
        end do

        !overwrite ispinsyst/nspinsyst
        do i = 1, nsystem
           do j = 1, nci_ini
           if (.not.equal(j)) then
             nspinsyst(i) = nspinsyst(j)
             ispinsyst(:,i) = ispinsyst(:,j)
             equal(j) = .true.
             exit
           end if
           end do
        end do

        nspinsyst(nsystem+1:maxsystem) = 0
        ispinsyst(:,nsystem+1:maxsystem) = 0

        do i = 1, nsystem !loop über separierte systeme
           do j = 1, nspinsyst(i) !loop über atom1 im system
              do k = j, nspinsyst(i)
                 ati = ispinsyst(j, i)
                 atj = ispinsyst(k, i)
                 rmaxab(atj, ati) = jab(lin(atj, ati))/1.2
                 rmaxab(ati, atj) = rmaxab(atj, ati)
              end do
           end do
        end do

        magdist = rmaxab

        !Find Sub-Systems with nspinsyst > maxmagnat
        ass = maxloc(nspinsyst, 1)
        do while (nspinsyst(ass) > maxmagnat)

           !Set Spins in this list as available for algorithm
           assigned = .false.
           do i = 1, nspinsyst(ass)
              do j = 1, nspinsyst(ass)
                 assigned(ispinsyst(i, ass), ispinsyst(j, ass)) = .true.
              end do ! End Loop over  j from 1 to nspinsyst(ass)
           end do ! End Loop over  i from 1 to nspin


           !Find largest distance:
           maxdist = maxval(rmaxab, mask=assigned)
           maxdistatoms = maxloc(rmaxab, mask=assigned)

           !If a Path is found between A and B
           if (maxdist < huge(1.0_wp)) then

              !get shortest Path from A to B
              cur_dist = shortest_distance(nspin, maxdistatoms(1), maxdistatoms(2), numnb, neigh, magdist, visited, precessor)
              current = maxdistatoms(2)

              !loop while A and B are still connected
              do while (cur_dist < huge(1.0_wp))

                 !find weakest link
                 max_link = 0
                 max_linkatoms = 0
                 i = 1
                 path = 0

                 do while (precessor(current) /= 0)
                    path(i) = current
                    current = precessor(current)
                    i = i + 1
                 end do ! End loop: while precessor(current) /= 0
                 path(i) = maxdistatoms(1)
                 i = i + 1
                 max_linkatoms(1) = path(i/2)
                 max_linkatoms(2) = path((i/2) + 1)

                 !Split weakest link, set distance to infinity
                 magdist(max_linkatoms(1), max_linkatoms(2)) = huge(1.0_wp)
                 magdist(max_linkatoms(2), max_linkatoms(1)) = huge(1.0_wp)

                 !Get next-shortest Path:
                 cur_dist = shortest_distance(nspin, maxdistatoms(1), maxdistatoms(2), numnb, neigh, magdist, visited, precessor)
                 current = maxdistatoms(2)
              end do ! cur_dist < huge(1.0_wp)

              !A and B are now Seperated:
              rmaxab(maxdistatoms(1), maxdistatoms(2)) = huge(1.0_sp)
              rmaxab(maxdistatoms(2), maxdistatoms(1)) = huge(1.0_sp)
           end if ! End if: while maxdist < huge(1.0d0)

           !Split into subsystems:
           !Overwrite old spinsystem:

           cur_dist = shortest_distance(nspin, maxdistatoms(1), maxdistatoms(2), numnb, neigh, magdist, visited, precessor)
           nspinsyst(ass) = count(visited)
           ispinsyst(:, ass) = 0
           k = 1
           do i = 1, nspin
              if (visited(i)) then
                 ispinsyst(k, ass) = i
                 k = k + 1
              end if
           end do ! End Loop over  i from 1 to count(visited)

           !add new spinsystem
           cur_dist = shortest_distance(nspin, maxdistatoms(2), maxdistatoms(1), numnb, neigh, magdist, visited, precessor)
           nsystem = nsystem + 1
           nspinsyst(nsystem) = count(visited)
           k = 1
           do i = 1, nspin
              if (visited(i)) then
                 ispinsyst(k, nsystem) = i
                 k = k + 1
              end if
           end do ! End Loop over  i from 1 to count(visited)
           ass = maxloc(nspinsyst, 1)

            if (nsystem .ge. nspin) then
               call env%warning("Created more fragments than atoms. This does not look good. Turning off fragmentation.")
               nsystem=1
               return
            end if
        end do ! End loop: while nspinsyst(ass) > maxmagnat


     end subroutine fragmentize

     subroutine frag_hess_diag( nat,hess,eig_calc,ispinsyst,nspinsyst,nsystem )
     !---------------------------------------------------------------------------------------------
     ! Purpose:
     ! Subroutine performs diagonalization of fragmented hessian.
     ! fragmentize must be called beforehand, fragments are available via this module.
     !
     ! How does it work?
     ! This subroutine uses masked arrays: an array (hess_mask) consists only of logicals and has
     ! the same size as the hessian matrix of the entire system. For all entries of the hessian
     ! matrix that belong to one fragment, the corresponding position in the masked array is set
     ! true. Using the intrinsic "pack and unpack" functions in combination with mask, the hessian
     ! entries for only one fragment are written into a smaller array (mini_hess) which is faster
     ! diagonalized due to the reduzed size (# of atoms in fragment < # of atoms in the entire system).
     !
     ! Input:
     ! nat      - Number of atoms of the entire system
     ! hess     - (Lindh) Hessian of the entire system
     !
     ! Output:
     ! hess     - diagonalized hessian (eigenvectors), overwritten
     ! eig_calc - diagonalized hessian (eigenvalues)
     !---------------------------------------------------------------------------------------------
     use xtb_mctc_accuracy, only : wp, sp
        implicit none
        !Dummy Arguments
        integer,  intent(in)     :: nat                          ! # of atoms
        real(sp), intent(inout)  :: hess(3*nat,3*nat)            ! input hessian
        real(sp), intent(out)    :: eig_calc(3*nat)              ! eigenvectors od entire system
        integer,  intent(in) :: ispinsyst(:,:)         ! array with list of atoms of each fragment
        integer,  intent(in) :: nspinsyst(:)           ! array with # of atoms for each fragment
        integer,  intent(in) :: nsystem                            ! # of fragments
        !Stack
        integer                  :: isystem
        integer                  :: i,j,ii,jj,k
        integer                  :: nat3                         ! # of atoms in the system * 3
        integer                  :: nat_cur                      ! # of atoms in fragment
        integer                  :: nat3_cur                     ! nat in fragments * 3
        integer                  :: info                         ! for ssyev
        integer                  :: lwork                        ! for ssyev
        logical                  :: hess_mask(3*nat,3*nat)       ! masked hessian array
        real(sp)                 :: ev_calc(3*nat,3*nat)         ! eigenvectors of entire system
        real(sp), allocatable    :: mini_hess(:,:)               ! eigenvectors of fragment
        real(sp), allocatable    :: eig(:)                       ! eigenvalues of fragment
        real(sp), allocatable    :: aux(:)                       ! for ssyev


        ev_calc  = 0.0e0_sp
        eig_calc = 0.0e0_sp
        nat3     = 3 * nat

!$omp parallel default(none) &
!$omp private(isystem,i,ii,j,jj,nat_cur,nat3_cur,mini_hess,hess_mask,eig,lwork,aux,info) &
!$omp shared(nsystem,ev_calc,eig_calc,hess,nspinsyst,ispinsyst)
!!$omp do schedule(static)
!$omp do reduction(+:ev_calc,eig_calc)
        do isystem = 1 , nsystem
           hess_mask = .false.
           do i = 1,nspinsyst(isystem)
              do j = 1,i

                 nat_cur =  nspinsyst(isystem)
                 nat3_cur = 3 * nat_cur
                 ii = 3*ispinsyst(i,isystem)
                 jj = 3*ispinsyst(j,isystem)
                 hess_mask(ii-2:ii,jj-2:jj) = .true.
                 hess_mask(jj-2:jj,ii-2:ii) = .true.

               end do
           end do

           allocate( mini_hess(nat3_cur,nat3_cur), source = 0.0e0_sp )
           allocate( eig(nat3_cur), source = 0.0e0_sp )

           mini_hess = reshape( pack( hess, mask = hess_mask ), shape( mini_hess ) )
           lwork = 1 + 6*nat3_cur + 2*nat3_cur**2
           allocate(aux(lwork))
           call ssyev ('V','U',nat3_cur,mini_hess,nat3_cur,eig,aux,lwork,info)
           deallocate(aux)
!!$omp critical
           ev_calc  = unpack( reshape( mini_hess, [ nat3_cur*nat3_cur ]  ), mask = hess_mask, field = ev_calc )
           eig_calc = unpack( eig, mask = any(hess_mask,1), field = eig_calc )
!!$omp end critical

           deallocate( mini_hess,eig )

        end do
!$omp end do
!$omp end parallel

        do i = 1,nat3
           ev_calc(:,i) = ev_calc(:,i) - ( sum( ev_calc(:,i) ) / nat3 )
        end do

        hess = ev_calc

        call eigsort4(nat3,hess,eig_calc)

     end subroutine frag_hess_diag

end module xtb_gfnff_fraghess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calculates the shortest path between two atoms
!start and goal are integers, determining the index in xyz
!https://de.wikipedia.org/wiki/Dijkstra-Algorithmus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function shortest_distance(nspin, start, goal, numnb, neighbours, input_distances, visited, precessor)
   use xtb_mctc_accuracy, only : wp, sp
   implicit none
   !Dummy Arguments:
   real(wp)             :: shortest_distance
   integer, intent(in)  :: nspin
   integer, intent(in)  :: start
   integer, intent(in)  :: goal
   integer, intent(in)  :: numnb
   integer, intent(in)  :: neighbours(numnb, nspin, 1)
   real(wp), intent(in) :: input_distances(nspin, nspin)
   logical, intent(out) :: visited(nspin)
   integer, intent(out) :: precessor(nspin)

   !Stack:
   integer  :: current
   integer  :: neighbour
   integer  :: i_neighbours
   integer  :: bonds
   real(wp) :: alt_dist
   real(wp) :: distance(nspin)

   !logical:: visited(nspin)

   bonds = 0
   precessor = 0
   distance = huge(distance)
   distance(start) = 0
   visited = .false.
   do while (.not.all(visited)) !as long there are unvisited nodes
      current = minloc(distance, 1,.not. visited)

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
            bonds = bonds + 1
            current = precessor(current)
         end do ! End loop: while precessor(current) /= 0
         exit
      else
         !loop over all neighbours of current atom
         do i_neighbours = 1, neighbours(numnb, current, 1)
            neighbour = neighbours(i_neighbours, current, 1)
            if (.not.visited(neighbour)) then
               !distanzupdate
               alt_dist = distance(current) + input_distances(current, neighbour)
               !write( *, * ) alt_dist, distance(current),
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


subroutine eigsort4(lab,u,ew)
   use xtb_mctc_accuracy, only : wp, sp
   implicit none
   integer  :: ii,k, j, i
   real(sp) :: pp, hilf

   integer,  intent(in)    :: lab
   real(sp), intent(inout) :: u(lab,lab)
   real(sp), intent(inout) :: ew(lab)

   do ii = 2, lab
      i = ii - 1
      k = i
      pp= ew(i)
      do j = ii, lab
         if (ew(j) .gt. pp) cycle
         k = j
         pp= ew(j)
      end do
      if (k .eq. i) cycle
      ew(k) = ew(i)
      ew(i) = pp
      do j=1,lab
         hilf=u(j,i)
         u(j,i)=u(j,k)
         u(j,k)=hilf
      end do
   end do

end subroutine eigsort4

pure subroutine com(n,at,xyz,sum3)
   use xtb_mctc_accuracy, only : wp
   implicit none
   !Dummy
   integer, intent(in)   :: n
   integer, intent(in)   :: at(n)
   real(wp), intent(in)  :: xyz(3,n)
   real(wp), intent(out) :: sum3(3)
   !Stack
   integer               :: i
   real(wp)              :: sumw
   real(wp)              :: sumwx,sumwy,sumwz

   sumw=0
   sumwx=0.d0
   sumwy=0.d0
   sumwz=0.d0

   do i = 1, n
      sumw = sumw + 1
      sumwx = sumwx + xyz(1,i)
      sumwy = sumwy + xyz(2,i)
      sumwz = sumwz + xyz(3,i)
   end do

   sum3(1) = sumwx / sumw
   sum3(2) = sumwy / sumw
   sum3(3) = sumwz / sumw

end subroutine com
