! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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

module shake_module
   use iso_fortran_env, only : wp => real64
   use setparam, only: xhonly,shake_mode
   implicit none
   integer, parameter :: ndim = 100000
   !integer :: shake_mode
   integer :: ncons = 0
   integer :: nconsu = 0
   integer :: conslistu(2,ndim)
   integer, allocatable :: conslist(:,:)
   real(wp),  allocatable :: distcons(:)
   real(wp),  allocatable :: dro(:,:)
   real(wp),  allocatable :: dr (:,:)
   integer, parameter :: maxcyc = 250
   real(wp), parameter :: tolshake = 1.d-7 

contains

subroutine init_shake(nat,at,xyz,wbo)
   use aoparam
   use fixparam, only : shakeset
   implicit none
   integer :: nat,at(nat)
   integer :: nbond
   real(wp)  :: xyz(3,nat)
   real(wp)  :: wbo(nat,nat)
   integer :: iat,jat,i,j,jmin
   real(wp)  :: minrij,rij,wthr
   real(wp)  :: drij(3)

   integer list(nat*(nat+1)/2),lin,ij
   real(wp) rco
   logical metalbond,rcut

   nconsu=0
   ! list of atoms interpreted as pairs: i1 j1 i2 j2 i3 j3 ....
   ! taken from fixb input
   if(shakeset%n.gt.0)then
      do i = 1, shakeset%n, 2
         nconsu = nconsu + 1
         conslistu(1,nconsu) = shakeset%atoms(i)
         conslistu(2,nconsu) = shakeset%atoms(i+1)
         write(*,*) 'SHAKE user input constraining bond ', &
            &      conslistu(1:2,nconsu)
      enddo
      !        switch of additional potential
   endif

   ! constrain X-H only
   if(shake_mode.eq.1)then
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
   if(shake_mode.eq.2)then
      list=0      
      do i = 1, nat
         do j =1, nat
            if(i.eq.j) cycle
            rij=(xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2
            rco=rad(at(j))+rad(at(i))
            ij=lin(i,j)
            rcut=0.52917726*sqrt(rij).lt.1.2*rco.and.list(ij).eq.0 ! to consider?
            metalbond=metal(at(i)).eq.1.or.metal(at(j)).eq.1       ! metal bond?
            wthr=0.5  ! WBO threshold. if WBO > thr, the bond is constrained
            if(metalbond) wthr=0.1 ! e.g. K...O in phosphates have WBO around 0.15
            !           write(*,'(2i3,3f10.4)') 
            !    .      at(i),at(j),0.52917726*sqrt(rij),rco,wbo(i,j)
            !           if(wbo(i,j).lt.wthr.and.wbo(i,j).gt.0.3.and.rcut)       ! warning
            !    .      write(*,*) 'bond slightly below threshold ', i,j, 
            !    .                 'not constrained'
            if(rcut.and.wbo(i,j).gt.wthr)then                       ! its relevant
               metalbond=metal(at(i)).eq.1.or.metal(at(j)).eq.1    ! do not constrain M bonds except Li/Be
               if(at(i).eq.3.or.at(i).eq.4.or.at(j).eq.3.or.at(j).eq.4) &
                  &         metalbond=.false.
               if(metalbond) then
                  if(i.lt.j) &
                     &          write(*,*)'init_shake: metal bond ', i,j, &
                     &                    'not constrained'
               else
                  list(ij)=1
                  nconsu = nconsu + 1
                  conslistu(1,nconsu) = i
                  conslistu(2,nconsu) = j
               endif
            endif
         enddo
      enddo

   endif

99 ncons  = nconsu
   if(nconsu.lt.1) return

   allocate(conslist(2,ncons),distcons(ncons),dro(3,ncons),dr(4,ncons))

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
   real(wp) xyzo(3,nat),xyz(3,nat),velo(3,nat),acc(3,nat),mass(nat)
   real(wp) virsh(3),tstep
   real(wp) xyzt(3,nat)
   real(wp) vel(3)
   integer i,icyc
   integer iat,jat
   logical conv
   real(wp) maxdev
   real(wp) r, dev,dist,denom
   real(wp) gcons, rmi, rmj
   real(wp) tau1,tau2
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

100 continue

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

         denom = 2.d0*(rmi+rmj)*(dr(1,i)*dro(1,i)+dr(2,i)*dro(2,i)+ &
            &                            dr(3,i)*dro(3,i))

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
      !        if(maxdev.gt.1.d-3) stop 'SHAKE error too large'        
   endif

   return
end subroutine do_shake

end module shake_module
