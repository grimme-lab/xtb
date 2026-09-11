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

       module xtb_gfnff_shake
       use xtb_param_atomicrad, only : atomicRad
       use xtb_setparam, only: set
       implicit none
       private :: atomicRad
       integer, parameter :: ndim = 100000
       logical :: pair_exists
       integer :: ncons = 0
       integer :: nconsu = 0
       integer :: conslistu(2,ndim)
       integer, allocatable :: conslist(:,:)
       real*8,  allocatable :: distcons(:)
       real*8,  allocatable :: dro(:,:)
       real*8,  allocatable :: dr (:,:)
       integer, parameter :: maxcyc = 5000
       real*8, parameter :: tolshake = 1.d-6
       integer, private, parameter :: metal(1:86) = [&
      & 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
      & 1, 0, 0, 0, 0, 0]

       contains

       subroutine init_shake(nat,at,xyz,topo)
       use xtb_mctc_convert, only : autoaa
       use xtb_fixparam, only : shakeset
       use xtb_gfnff_topology, only : TGFFTopology
       implicit none
       type(TGFFTopology), intent(in) :: topo
       integer :: nat,at(nat)
       real*8  :: xyz(3,nat)
       real*8  :: wbo(nat,nat)
       integer :: iat,jat,i,j,jmin
       real*8  :: minrij,rij,wthr
       real*8  :: drij(3)
       
       integer list(nat*(nat+1)/2),lin,ij
       
       real*8 rco
       logical metalbond,rcut
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       nconsu=0

       ! constrain X-H only
       if(set%shake_mode.eq.1)then
          do i = 1, nat
                if(at(i).eq.1) then
                   minrij=1000.d0
                   do j = 1, nat
                      if(j.ne.i)then
                            rij=(xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2
                            if(rij.lt.minrij) then
                               minrij=rij
                               jmin=j
                            endif
                      endif
                   enddo
                   if(at(jmin).eq.1) then
                      if(jmin.gt.i) then
                            nconsu = nconsu + 1
                            conslistu(1,nconsu) = i
                            conslistu(2,nconsu) = jmin
                      endif
                   else
                      nconsu = nconsu + 1
                      conslistu(1,nconsu) = i
                      conslistu(2,nconsu) = jmin
                   endif
                endif
          enddo
       endif
       ! all bonds
       if(set%shake_mode.eq.2)then
          nconsu  = topo%nbond
          conslistu(1:2,1:nconsu)=topo%blist(1:2,1:nconsu)
       endif
       print *, "SHAKE is ", set%shake_md, " ncons = ", nconsu, "shake_mode = ", set%shake_mode
       
       if(shakeset%n.gt.0)then
             do i = 1, shakeset%n, 2
                pair_exists = .false.
                do j = 1, nconsu
                   if ((conslistu(1,j) == shakeset%atoms(i) .and. conslistu(2,j) == shakeset%atoms(i+1)) .or. &
                      (conslistu(1,j) == shakeset%atoms(i+1) .and. conslistu(2,j) == shakeset%atoms(i))) then
                      pair_exists = .true.
                      if (pair_exists) write(*,*) 'input pair',  shakeset%atoms(i), '-', shakeset%atoms(i+1), 'for shake is already considered.'
                      exit
                   endif
                enddo
                
                if (.not. pair_exists) then
                   nconsu = nconsu + 1
                   conslistu(1,nconsu) = shakeset%atoms(i)
                   conslistu(2,nconsu) = shakeset%atoms(i+1)
                   write(*,*) 'SHAKE user input constraining bond ', &
                      &      conslistu(1:2,nconsu)
                endif
             enddo
       endif
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
99     ncons  = nconsu
       if(nconsu.lt.1) then
          set%shake_md = .false.
          print *, "gfnff shake_md turned off"
          return
       endif
     print *, "SHAKE is ", set%shake_md, " ncons = ", ncons, "shake_mode = ", set%shake_mode
       allocate(conslist(2,ncons),distcons(ncons),&
       &         dro(3,ncons),dr(4,ncons))
    
       conslist(1:2,1:ncons)=conslistu(1:2,1:ncons) 
       do i = 1, ncons
          iat = conslist(1,i)
          jat = conslist(2,i)
          drij(1:3) = xyz(1:3,iat) - xyz(1:3,jat)
          distcons(i) = drij(1)**2 + drij(2)**2 + drij(3)**2
       enddo
       
       return
       end subroutine init_shake

       subroutine do_shake(nat,xyzo,xyz,velo,acc,mass,tstep)
       implicit none
       integer nat
       real*8 xyzo(3,nat),xyz(3,nat),velo(3,nat),acc(3,nat),mass(nat)
       real*8 virsh(3),tstep
       real*8 xyzt(3,nat)
       real*8 vel(3)
       integer i,icyc
       integer iat,jat
       logical conv
       real*8 maxdev
       real*8 r, dev,dist,denom
       real*8 gcons, rmi, rmj
       real*8 tau1,tau2
       integer jmaxdev

        conv = .false.
        icyc = 0

        xyzt = xyz

        tau1 = 1.d0/tstep
        tau2 = tau1*tau1

        do i = 1, ncons
         iat = conslist(1,i)
         jat = conslist(2,i)
         dro(1:3,i) = xyzo(1:3,iat) - xyzo(1:3,jat)
        enddo

100     continue

          maxdev = 0.d0

          do i = 1, ncons
           iat = conslist(1,i)
           jat = conslist(2,i)
           dr(1:3,i) = xyzt(1:3,iat) - xyzt(1:3,jat)
           dr(4,i)   = dr(1,i)**2 + dr(2,i)**2 + dr(3,i)**2
           dist = distcons(i)
           dev = abs(dr(4,i) - dist) / dist
           if(dev.gt.maxdev) then
            maxdev = dev
            jmaxdev = i
           endif
          enddo

          if(maxdev.lt.tolshake) conv = .true.

          if(.not.conv) then

           do i = 1, ncons

            iat = conslist(1,i)
            jat = conslist(2,i)
            dist = distcons(i)

            rmi = 1.d0/mass(iat)
            rmj = 1.d0/mass(jat)

            denom = 2.d0*(rmi+rmj)*(dr(1,i)*dro(1,i)+dr(2,i)*dro(2,i)+&
     &                              dr(3,i)*dro(3,i))

            gcons = (dist-dr(4,i))/denom

            xyzt(1:3,iat) = xyzt(1:3,iat) + rmi*gcons*dro(1:3,i)
            xyzt(1:3,jat) = xyzt(1:3,jat) - rmj*gcons*dro(1:3,i)

           enddo

          endif

          icyc = icyc + 1

        
        if(.not.conv.and.icyc.le.maxcyc) goto 100

        if(conv) then
         velo = velo + (xyzt-xyz)*tau1
         acc  = acc  + (xyzt-xyz)*tau2
         xyz  = xyzt
        else
         write(*,*)'SHAKE did not converge! maxdev=',maxdev
         !if(maxdev.gt.1.d-3) stop 'SHAKE error too large'
        endif

       return
       end subroutine do_shake

       end module xtb_gfnff_shake
