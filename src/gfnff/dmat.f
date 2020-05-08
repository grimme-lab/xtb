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

cccccccccccccccccccccccccccccccccccccccccccccc
c density matrix
c C: MO coefficient
c X: scratch
c P  dmat
cccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dmat(ndim,focc,C,P)
      implicit none
      integer ndim
      real*8 focc(*)
      real*8 C(ndim,ndim)
      real*8 P(ndim,ndim)
      integer i,m
      real*8,allocatable ::Ptmp(:,:)

      allocate(Ptmp(ndim,ndim))
      do m=1,ndim
         do i=1,ndim
            Ptmp(i,m)=C(i,m)*focc(m)
         enddo
      enddo
      call DGEMM('N','T',ndim,ndim,ndim,1.0d0,C,
     .                   ndim,Ptmp,ndim,0.0d0,P,ndim)
      deallocate(Ptmp)

      end

