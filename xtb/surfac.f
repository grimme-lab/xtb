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

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine surfac(fname,n,xyz,at)
      use iso_fortran_env, wp => real64   
      implicit none
      character(len=*),intent(in)       :: fname
      integer, intent(in)               :: n,at(n) !number of atoms,Ordnunszahlen
      real(wp), intent(in)              :: xyz(3,n)

      real(wp),allocatable              :: s(:,:) !surface 
      real(wp),allocatable              :: atom_weight(:,:)!atom and weight
      real(wp)                          :: point(3),rx,ry,rz,r,r2,rsas
      real(wp)                          :: rad(n),rad2(n),d3rad(94)

      integer, parameter                :: nangsa=86
!     integer, parameter                :: nangsa=230
      real(wp)                          :: grida(4,nangsa)
      include 'grida86.inc'
!     include 'grida230.inc'
      character(len=2), external        :: asym
      integer                           :: i,j,k,np,ip

      ! D3 radii in Bohr
      data d3rad/
     .  2.18230009,  1.73469996,  3.49559999,  3.09820008,  3.21600008,
     .  2.91030002,  2.62249994,  2.48169994,  2.29959989,  2.13739991,
     .  3.70819998,  3.48390007,  4.01060009,  3.79169989,  3.50169992,
     .  3.31069994,  3.10459995,  2.91479993,  4.24109983,  4.10349989,
     .  3.89030004,  3.76419997,  3.72110009,  3.44140005,  3.54620004,
     .  3.44210005,  3.43269992,  3.34619999,  3.30080009,  3.23090005,
     .  3.95790005,  3.86190009,  3.66249990,  3.52679992,  3.36619997,
     .  3.20959997,  4.61759996,  4.47639990,  4.21960020,  4.05970001,
     .  3.85960007,  3.75430012,  3.56900001,  3.46230006,  3.39750004,
     .  3.35249996,  3.33080006,  3.46199989,  4.26230001,  4.18739986,
     .  4.01499987,  3.89010000,  3.73799992,  3.58890009,  5.05670023,
     .  5.18139982,  4.62610006,  4.62010002,  4.57019997,  4.52710009,
     .  4.48960018,  4.45149994,  4.42339993,  4.12430000,  4.24270010,
     .  4.15409994,  4.27939987,  4.24499989,  4.22079992,  4.19859982,
     .  4.01300001,  4.24499989,  4.09800005,  3.98550010,  3.89549994,
     .  3.74900007,  3.44560003,  3.35249996,  3.25640011,  3.35990000,
     .  4.31269979,  4.27640009,  4.11749983,  4.00540018,  3.86439991,
     .  3.72160006,  5.07959986,  4.92939997,  4.70429993,  4.42519999,
     .  4.45940018,  4.39569998,  4.35389996,  4.43410015/


      do i=1,n
         rad (i) =(d3rad(at(i))*1.2)    ! scale factor adjusted to get vdW contacts right
         rad2(i) =(d3rad(at(i))*1.2)**2
      enddo

      np=n*nangsa
      allocate(s(3,np))
      allocate(atom_weight(2,np)) ! point belonging to atom and atomweight

      k=0 ! number of necessary gridpoints
      do i=1,n
        rsas=rad(i)
        grid_point: do ip=1,nangsa
!          grid point position
           rx = xyz(1,i) + rsas*grida(1,ip)
           ry = xyz(2,i) + rsas*grida(2,ip)
           rz = xyz(3,i) + rsas*grida(3,ip)
           do j=1,n
              if(i.eq.j) cycle
              r2=(xyz(1,j)-rx)**2+(xyz(2,j)-ry)**2+(xyz(3,j)-rz)**2
              if(r2.le.rad2(j)) cycle grid_point  ! closer to another atom, skip grid point
           enddo
           k=k+1
           s(1,k)=rx
           s(2,k)=ry
           s(3,k)=rz
           atom_weight(1,k) = i            ! point belonging to atom i
           atom_weight(2,k) = grida(4,ip)   ! weight belonging to point ip
         enddo grid_point
      enddo

      np=k
      write(*,*) 'generated ',np,' surface points'

      open(unit=83,file=fname)
      do i=1,np
      write(83,'(3E16.8,1x,i0,1x,f10.7)') s(1:3,i),
     & int(atom_weight(1,i)), atom_weight(2,i)
      enddo
      close(83)

! surface debug output
!     write(13,'(''$coord'')')
!     do i=1,n
!        write(13,'(3F14.6,5x,a2)') xyz(1:3,i),asym(at(i))
!     enddo
!     do i=1,np
!        write(13,'(3F14.6,5x,''XX'')') s(1:3,i)
!     enddo
!     write(13,'(''$end'')')

      deallocate(s)
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine surfac2(surface,n,xyz,at)
      use iso_fortran_env, wp => real64   
      use grid_module
      use dftd4
      implicit none
      integer, intent(in)        :: n,at(n) !number of atoms,Ordnunszahlen
      real(wp), intent(in)       :: xyz(3,n)
      type(tb_grid),intent(out)  :: surface

      real(wp),allocatable       :: s(:,:) !surface 
      real(wp),allocatable       :: atom_weight(:,:)!atom and weight
      real(wp)                   :: point(3),rx,ry,rz,r,r2,rsas
      real(wp)                   :: rad(n),rad2(n),d3rad(94)

      integer, parameter         :: nangsa=38
!     integer, parameter         :: nangsa=230
      real(wp)                   :: grida(4,nangsa)
      include 'grida38.inc'
!     include 'grida230.inc'
      character(len=2), external :: asym
      integer                    :: i,j,k,np,ip


      do i=1,n
         rad (i) =sqrt(3.0_wp)*r4r2(at(i)) ! scale factor adjusted to get vdW contacts right
         rad2(i) =3.0_wp*r4r2(at(i))**2
      enddo

      np=n*nangsa
      allocate(s(3,np))
      allocate(atom_weight(2,np)) ! point belonging to atom and atomweight

      k=0 ! number of necessary gridpoints
      do i=1,n
        rsas=rad(i)
        grid_point: do ip=1,nangsa
!          grid point position
           rx = xyz(1,i) + rsas*grida(1,ip)
           ry = xyz(2,i) + rsas*grida(2,ip)
           rz = xyz(3,i) + rsas*grida(3,ip)
           do j=1,n
              if(i.eq.j) cycle
              r2=(xyz(1,j)-rx)**2+(xyz(2,j)-ry)**2+(xyz(3,j)-rz)**2
              if(r2.le.rad2(j)) cycle grid_point  ! closer to another atom, skip grid point
           enddo
           k=k+1
           s(1,k)=rx
           s(2,k)=ry
           s(3,k)=rz
           atom_weight(1,k) = i            ! point belonging to atom i
           atom_weight(2,k) = grida(4,ip)   ! weight belonging to point ip
         enddo grid_point
      enddo

      np=k
      surface%n = np
      allocate(surface%x(3,np),surface%w(np),surface%at(np))
      write(output_unit,'("generated",1x,i0,1x,"surface points")') np

      do i=1,np
         surface%x(:,i) = s(:,i)
         surface%w(i) = atom_weight(2,i)
         surface%at(i) = nint(atom_weight(1,i))
      enddo

      deallocate(s,atom_weight)
      end subroutine surfac2
