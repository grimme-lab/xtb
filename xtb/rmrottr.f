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


      subroutine rmrottr(natoms,mass,vel_atom,c)
      implicit none
      integer natoms
      real*8 vel_atom(3,natoms),c(3,natoms),mass(natoms)

      integer         :: i
      real*8          :: rlm(3), ram(3), omega(3)
      real*8          :: ixx,iyy,izz,ixy,ixz,iyz,dummy
      real*8          :: fixx,fiyy,fizz,fixy,fixz,fiyz,COM(3)
      real*8          :: inertia(3,3),angmom(3)
      real*8          :: tmass                  

      rlm = 0.0
      ram = 0.0
      call centerofmass(natoms,c,mass,tmass,COM)
      angmom = 0.0
      do i = 1, natoms
         c(1,i)=c(1,i)-COM(1)
         c(2,i)=c(2,i)-COM(2)
         c(3,i)=c(3,i)-COM(3)
         angmom(1) = angmom(1) + mass(i) * ( c(2,i) * vel_atom(3,i) -
     &                                      c(3,i) * vel_atom(2,i) )
         angmom(2) = angmom(2) + mass(i) * ( c(3,i) * vel_atom(1,i) -
     &                                      c(1,i) * vel_atom(3,i) )
         angmom(3) = angmom(3) + mass(i) * ( c(1,i) * vel_atom(2,i) -
     &                                      c(2,i) * vel_atom(1,i) )
      end do
      ixx = 0.0
      iyy = 0.0
      izz = 0.0
      ixy = 0.0
      ixz = 0.0
      iyz = 0.0
      do i = 1, natoms
         ixx = ixx + mass(i) * ( c(2,i) * c(2,i) + c(3,i) * c(3,i) )
         iyy = iyy + mass(i) * ( c(3,i) * c(3,i) + c(1,i) * c(1,i) )
         izz = izz + mass(i) * ( c(1,i) * c(1,i) + c(2,i) * c(2,i) )
         ixy = ixy - mass(i) *   c(1,i) * c(2,i)
         ixz = ixz - mass(i) *   c(1,i) * c(3,i)
         iyz = iyz - mass(i) *   c(2,i) * c(3,i)
      end do
      inertia(1,1) = ixx
      inertia(2,2) = iyy
      inertia(3,3) = izz
      inertia(1,2) = ixy
      inertia(2,1) = ixy
      inertia(1,3) = ixz
      inertia(3,1) = ixz
      inertia(2,3) = iyz
      inertia(3,2) = iyz
      call dmatinv(inertia,3,3,dummy)
      omega =  matmul(inertia,angmom)
      do i = 1, natoms
          rlm(1) = rlm(1) + mass(i) * vel_atom(1,i)
          rlm(2) = rlm(2) + mass(i) * vel_atom(2,i)
          rlm(3) = rlm(3) + mass(i) * vel_atom(3,i)
      end do
      do i = 1, natoms
        ram(1) =        ( omega(2) * c(3,i) - omega(3) * c(2,i) )
        ram(2) =        ( omega(3) * c(1,i) - omega(1) * c(3,i) )
        ram(3) =        ( omega(1) * c(2,i) - omega(2) * c(1,i) )

        vel_atom(1,i) = vel_atom(1,i) - rlm(1) / tmass  - ram(1)
        vel_atom(2,i) = vel_atom(2,i) - rlm(2) / tmass  - ram(2)
        vel_atom(3,i) = vel_atom(3,i) - rlm(3) / tmass  - ram(3)

        c(1,i)=c(1,i)+COM(1)
        c(2,i)=c(2,i)+COM(2)
        c(3,i)=c(3,i)+COM(3)
      end do

      end

      subroutine centerofmass(natoms,c,mass,totmass,COM)
      implicit none
      integer natoms
      real*8 c(3,natoms),totmass,mass(natoms),COM(3)

      integer         :: i,j

      COM(1)=0.0
      COM(2)=0.0
      COM(3)=0.0

      totmass=0
      do i = 1, natoms
         totmass=totmass+mass(i)
         COM(1) = COM(1) + mass(i)*c(1,i)
         COM(2) = COM(2) + mass(i)*c(2,i)
         COM(3) = COM(3) + mass(i)*c(3,i)
      enddo

      COM(1)=COM(1)/totmass
      COM(2)=COM(2)/totmass
      COM(3)=COM(3)/totmass 
c     write(9,'(3f20.9)')(COM(j),j=1,3)
      end
