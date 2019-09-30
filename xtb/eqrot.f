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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      logical function equalrot(i,j,nall,thr,rot)    
      implicit none
      integer i,j,nall
      real*8 rot(3,nall),r,r1,thr

      equalrot=.false.

      r1=rot(1,i)**2+rot(2,i)**2+rot(3,i)**2
      r=(rot(1,i)-rot(1,j))**2
     . +(rot(2,i)-rot(2,j))**2
     . +(rot(3,i)-rot(3,j))**2

      if(sqrt(r)/sqrt(r1).lt.thr) equalrot=.true.

      end

      logical function equalrot2(i,j,nall,thr,rot)    
      implicit none
      integer i,j,nall
      real*8 rot(3,0:nall),r,r1,r2,thr

      equalrot2=.false.

      r1=rot(1,i)**2+rot(2,i)**2+rot(3,i)**2
      r2=rot(1,j)**2+rot(2,j)**2+rot(3,j)**2
      r=(rot(1,i)-rot(1,j))**2
     . +(rot(2,i)-rot(2,j))**2
     . +(rot(3,i)-rot(3,j))**2

      if(2.*sqrt(r)/(sqrt(r1)+sqrt(r2)).lt.thr) equalrot2=.true.

      end

      real*8 function rotdiff(i,j,nall,rot)    
      implicit none
      integer i,j,nall
      real*8 rot(3,nall),r,r1,r2

      r1=rot(1,i)**2+rot(2,i)**2+rot(3,i)**2
      r2=rot(1,j)**2+rot(2,j)**2+rot(3,j)**2
      r=(rot(1,i)-rot(1,j))**2
     . +(rot(2,i)-rot(2,j))**2
     . +(rot(3,i)-rot(3,j))**2

      rotdiff=2.*sqrt(r)/(sqrt(r1)+sqrt(r2))

      end
