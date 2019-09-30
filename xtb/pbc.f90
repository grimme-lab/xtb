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

module pbc
   use iso_fortran_env, wp => real64
   implicit none
!! ========================================================================
!  Helper functions for periodic boundary conditions, mostly from dftd3
!! ------------------------------------------------------------------------

contains

!! ------------------------------------------------------------------------
!  generate a supercell based on a realspace cutoff, this subroutine
!  doesn't know anything about the convergence behaviour of the
!  associated property. rthr is assumed to be the *quadratic* threshold!
pure subroutine get_realspace_cutoff(lat,rthr,tau_max)
   implicit none
   real(wp),intent(in)  :: rthr
   real(wp),intent(in)  :: lat(3,3)
   integer, intent(out) :: tau_max(3)
   real(wp) :: r_cutoff
   real(wp) :: normx(3),normy(3),normz(3)
   real(wp) :: cos10,cos21,cos32

   r_cutoff=sqrt(rthr)
!  write(*,*) 'lat',lat
   ! find normal to the plane...
   call crossproduct(lat(:,2),lat(:,3),normx)
   call crossproduct(lat(:,3),lat(:,1),normy)
   call crossproduct(lat(:,1),lat(:,2),normz)
!  write(*,*) 'normy',normy
   ! ...normalize it...
   normx=normx/norm2(normx)
   normy=normy/norm2(normy)
   normz=normz/norm2(normz)
!  write(*,*) 'normy_',normy
   ! cos angles between normals and lattice vectors
   cos10=sum(normx*lat(:,1))
   cos21=sum(normy*lat(:,2))
   cos32=sum(normz*lat(:,3))
!  write(*,*) 'cos32',cos32
!  tau_max(1)=abs(2*r_cutoff/cos10)
!  tau_max(2)=abs(2*r_cutoff/cos21)
!  tau_max(3)=abs(2*r_cutoff/cos32)
   tau_max(1)=ceiling(abs(r_cutoff/cos10))
   tau_max(2)=ceiling(abs(r_cutoff/cos21))
   tau_max(3)=ceiling(abs(r_cutoff/cos32))

contains

pure subroutine crossproduct(a,b,c)
   implicit none

   real(wp),intent(in)  :: a(3),b(3)
   real(wp),intent(out) :: c(3)
   real(wp) :: x,y,z

   x=a(2)*b(3)-b(2)*a(3)
   y=a(3)*b(1)-b(3)*a(1)
   z=a(1)*b(2)-b(1)*a(2)
   c=(/x,y,z/)
end subroutine crossproduct

end subroutine get_realspace_cutoff

real(wp) function volume(lat)
   implicit none
   real(wp),intent(in) ::lat(3,3)
   real(wp) :: det

   det=lat(1,1)*lat(2,2)*lat(3,3)+lat(1,2)*lat(2,3)*lat(3,1) + &
   &   lat(1,3)*lat(2,1)*lat(3,2)-lat(1,3)*lat(2,2)*lat(3,1) - &
   &   lat(1,2)*lat(2,1)*lat(3,3)-lat(1,1)*lat(2,3)*lat(3,2)
   volume=abs(det)
end function volume

subroutine xyz_to_abc(xyz,abc,lat,n)
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: lat(3,3)
   real(wp),intent(out) :: abc(3,n)

   real(wp) :: lat_1(3,3)
   integer  :: i,j,k

   call inv_cell(lat,lat_1)

   abc(:,:n)=0.0d0
   do i=1,n
      do j=1,3
         do k=1,3
            abc(j,i)=abc(j,i)+lat_1(j,k)*xyz(k,i)
         enddo !k
         abc(j,i)=mod(abc(j,i),1.0_wp)
      enddo !j
   enddo !i

end subroutine xyz_to_abc

pure subroutine inv_cell(x,a) !x is normal lat, a is lat^(-1)
   implicit none
   real(wp),intent(in)  :: x(3,3) !unitcell vectors in direct space
   real(wp),intent(out) :: a(3,3) !unitcell vectors in reciprocal space
   real(wp) :: det

   a=0.0_wp
   det=x(1,1)*x(2,2)*x(3,3)+x(1,2)*x(2,3)*x(3,1)+x(1,3)*x(2,1)* &
   &   x(3,2)-x(1,3)*x(2,2)*x(3,1)-x(1,2)*x(2,1)*x(3,3)-x(1,1)* &
   &   x(2,3)*x(3,2)
!  write(*,*)'Det:',det
   a(1,1)=x(2,2)*x(3,3)-x(2,3)*x(3,2)
   a(2,1)=x(2,3)*x(3,1)-x(2,1)*x(3,3)
   a(3,1)=x(2,1)*x(3,2)-x(2,2)*x(3,1)
   a(1,2)=x(1,3)*x(3,2)-x(1,2)*x(3,3)
   a(2,2)=x(1,1)*x(3,3)-x(1,3)*x(3,1)
   a(3,2)=x(1,2)*x(3,1)-x(1,1)*x(3,2)
   a(1,3)=x(1,2)*x(2,3)-x(1,3)*x(2,2)
   a(2,3)=x(1,3)*x(2,1)-x(1,1)*x(2,3)
   a(3,3)=x(1,1)*x(2,2)-x(1,2)*x(2,1)
   a=a/det
end subroutine inv_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine abc_to_xyz(abc,xyz,lat,n)
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: abc(3,n)
   real(wp),intent(in)  :: lat(3,3)
   real(wp),intent(out) :: xyz(3,n)

   integer i,j,k

   xyz(:,:n)=0.0_wp
   do i=1,n
      do j=1,3
         do k=1,3
            xyz(j,i)=xyz(j,i)+lat(j,k)*abc(k,i)
         enddo !k
      enddo !j
   enddo !i

end subroutine abc_to_xyz

pure subroutine get_translation(t1,t2,t3,xyz,abc,out)
  implicit none
  integer,  intent(in)  :: t1,t2,t3
  real(wp), intent(in)  :: xyz(3)
  real(wp), intent(in)  :: abc(3,3)
  real(wp), intent(inout) :: out(3)
  out=0.0_wp
  out(1)=(xyz(1)+t1*abc(1,1)+t2*abc(1,2)+t3*abc(1,3))
  out(2)=(xyz(2)+t2*abc(2,1)+t2*abc(2,2)+t3*abc(2,3))
  out(3)=(xyz(3)+t3*abc(3,1)+t2*abc(3,2)+t3*abc(3,3))
end subroutine get_translation

pure function mdet3 (a) result (det)
  !! determinant of 3x3 matrix
  implicit none
  real*8, intent(in)  :: a(3,3)
  real*8 :: det

  det =   a(1,1)*a(2,2)*a(3,3)  &
       - a(1,1)*a(2,3)*a(3,2)  &
       - a(1,2)*a(2,1)*a(3,3)  &
       + a(1,2)*a(2,3)*a(3,1)  &
       + a(1,3)*a(2,1)*a(3,2)  &
       - a(1,3)*a(2,2)*a(3,1)
end function mdet3

pure function minv3(a) result(b)
  !! http://fortranwiki.org/fortran/show/Matrix+inversion
  !! Performs a direct calculation of the inverse of a 3×3 matrix.
  real*8, intent(in) :: a(3,3)   !! Matrix
  real*8             :: b(3,3)   !! Inverse matrix
  real*8             :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1/(a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2)&
       - a(1,2)*a(2,1)*a(3,3) + a(1,2)*a(2,3)*a(3,1)&
       + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1))

  ! Calculate the inverse of the matrix
  b(1,1) = +detinv * (a(2,2)*a(3,3) - a(2,3)*a(3,2))
  b(2,1) = -detinv * (a(2,1)*a(3,3) - a(2,3)*a(3,1))
  b(3,1) = +detinv * (a(2,1)*a(3,2) - a(2,2)*a(3,1))
  b(1,2) = -detinv * (a(1,2)*a(3,3) - a(1,3)*a(3,2))
  b(2,2) = +detinv * (a(1,1)*a(3,3) - a(1,3)*a(3,1))
  b(3,2) = -detinv * (a(1,1)*a(3,2) - a(1,2)*a(3,1))
  b(1,3) = +detinv * (a(1,2)*a(2,3) - a(1,3)*a(2,2))
  b(2,3) = -detinv * (a(1,1)*a(2,3) - a(1,3)*a(2,1))
  b(3,3) = +detinv * (a(1,1)*a(2,2) - a(1,2)*a(2,1))
end function minv3

pure logical function sametest(v)
  implicit none
  real*8, intent(in) ::v(3)
  real*8 ::l
  real*8, parameter :: s=1.0E-9

  sametest=.false.
  l=sqrt(dot_product(v,v))
  if (l<=s) sametest=.true.
end function sametest
end module pbc
