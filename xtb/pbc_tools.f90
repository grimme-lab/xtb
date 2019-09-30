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

!> @brief helper tools for periodic boundary conditions
module pbc_tools
   use iso_fortran_env, only : wp => real64
   use mctc_constants,  only : pi
! ---------------------------> IMPORTANT <--------------------------- !
!  DO NOT INCLUDE ANY DERIVED TYPES FROM THE TBDEFS INTO THIS MODULE  !
!   features that need to work on the tbdefs should go into pbc.f90   !
! ------------------------------------------------------------------- !
   implicit none
   private :: wp, pi
   public

   real(wp),parameter :: tpi = 2.0_wp*pi

   interface coord_trafo
      module procedure :: coord_trafo_12
      module procedure :: coord_trafo_inplace
   end interface coord_trafo

contains

!> @brief convert cell parameters to direct lattice
pure subroutine cell_to_dlat(cellpar,lattice)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp),intent(in)  :: cellpar(6)   !< cell parameters
   real(wp),intent(out) :: lattice(3,3) !< direct lattice
   real(wp)             :: dvol

   dvol = cell_to_dvol(cellpar)

   associate( alen => cellpar(1), blen => cellpar(2), clen => cellpar(3), &
      &       alp  => cellpar(4), bet  => cellpar(5), gam  => cellpar(6) )

   lattice(1,1) = alen
   lattice(2,1) = 0.0_wp
   lattice(3,1) = 0.0_wp
   lattice(3,2) = 0.0_wp
   lattice(1,2) = blen*cos(gam)
   lattice(2,2) = blen*sin(gam)
   lattice(1,3) = clen*cos(bet)
   lattice(2,3) = clen*(cos(alp) - cos(bet)*cos(gam))/sin(gam);
   lattice(3,3) = dvol/(alen*blen*sin(gam))

   end associate

end subroutine cell_to_dlat

!> @brief convert cell parameters to reciprocal lattice
pure subroutine cell_to_rlat(cellpar,rec_lat)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp),intent(in)  :: cellpar(6)   !< cell parameters
   real(wp),intent(out) :: rec_lat(3,3) !< reciprocal lattice
   real(wp)             :: dvol

   dvol = cell_to_dvol(cellpar)

   associate( alen => cellpar(1), blen => cellpar(2), clen => cellpar(3), &
      &       alp  => cellpar(4), bet  => cellpar(5), gam  => cellpar(6) )

   rec_lat(1,1) = tpi/alen
   rec_lat(1,2) = 0.0_wp
   rec_lat(1,3) = 0.0_wp
   rec_lat(2,1) =-tpi*cos(gam)/(alen*sin(gam))
   rec_lat(2,2) = tpi/(blen*sin(gam))
   rec_lat(2,3) = 0.0_wp
   rec_lat(3,1) = tpi*(cos(gam)*cos(alp) - cos(bet))*blen*clen/(dvol*sin(gam))
   rec_lat(3,2) = tpi*(cos(gam)*cos(bet) - cos(alp))*alen*clen/(dvol*sin(gam))
   rec_lat(3,3) = tpi*sin(gam)*alen*blen/dvol

   end associate

end subroutine cell_to_rlat

!> @brief convert direct lattice to cell parameters
pure subroutine dlat_to_cell(lattice,cellpar)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp),intent(in)  :: lattice(3,3) !< direct lattice
   real(wp),intent(out) :: cellpar(6)   !< cell parameters

   associate( alen => cellpar(1), blen => cellpar(2), clen => cellpar(3), &
      &       alp  => cellpar(4), bet  => cellpar(5), gam  => cellpar(6) )

   alen = norm2(lattice(:,1))
   blen = norm2(lattice(:,2))
   clen = norm2(lattice(:,3))

   alp = acos(dot_product(lattice(:,2),lattice(:,3))/(blen*clen))
   bet = acos(dot_product(lattice(:,1),lattice(:,3))/(alen*clen))
   gam = acos(dot_product(lattice(:,1),lattice(:,2))/(alen*blen))

   end associate

end subroutine dlat_to_cell

!> @brief convert direct lattice to cell parameters
pure subroutine dlat_to_rlat(lattice,rec_lat)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp),intent(in)  :: lattice(3,3) !< direct lattice
   real(wp),intent(out) :: rec_lat(3,3) !< reciprocal lattice
!  real(wp)             :: dvol
   rec_lat = tpi*transpose(mat_inv_3x3(lattice))
!  alternative implementation
!  dvol = dlat_volume(lattice)
!  associate( a => lattice(:,1), b => lattice(:,2), c => lattice(:,3) )
!     rec_lat(:,1) = tpi*cross(b,c)/dvol
!     rec_lat(:,2) = tpi*cross(c,a)/dvol
!     rec_lat(:,3) = tpi*cross(a,b)/dvol
!  end associate
end subroutine dlat_to_rlat

!> @brief implements the cross/vector product between two 3D vectors
pure function cross(a,b) result(c)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp),intent(in) :: a(3)
   real(wp),intent(in) :: b(3)
   real(wp)            :: c(3)
   c(1)=a(2)*b(3)-b(2)*a(3)
   c(2)=a(3)*b(1)-b(3)*a(1)
   c(3)=a(1)*b(2)-b(1)*a(2)
end function cross

!> @brief determinat of 3×3 matrix
pure function mat_det_3x3(a) result (det)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp),intent(in) :: a(3,3)
   real(wp)            :: det

   det =  a(1,1)*a(2,2)*a(3,3)  &
      & - a(1,1)*a(2,3)*a(3,2)  &
      & - a(1,2)*a(2,1)*a(3,3)  &
      & + a(1,2)*a(2,3)*a(3,1)  &
      & + a(1,3)*a(2,1)*a(3,2)  &
      & - a(1,3)*a(2,2)*a(3,1)

end function mat_det_3x3

!> Performs a direct calculation of the inverse of a 3×3 matrix.
!
!  reference: http://fortranwiki.org/fortran/show/Matrix+inversion
pure function mat_inv_3x3(a) result(b)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp),intent(in) :: a(3,3)   !< Matrix
   real(wp)            :: b(3,3)   !< Inverse matrix
   real(wp)            :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1.0_wp/mat_det_3x3(a)

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
end function mat_inv_3x3

!> @brief calculate the cell volume from the cell parameters
pure function cell_to_dvol(cellpar) result(dvol)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp),intent(in) :: cellpar(6) !< cell parameters
   real(wp)            :: vol2
   real(wp)            :: dvol

   associate( alen => cellpar(1), blen => cellpar(2), clen => cellpar(3), &
      &       alp  => cellpar(4), bet  => cellpar(5), gam  => cellpar(6) )

   vol2 = 1.0_wp - cos(alp)**2 - cos(bet)**2 - cos(gam)**2 &
      & + 2.0_wp*cos(alp)*cos(bet)*cos(gam)

   dvol = sqrt(abs(vol2))*alen*blen*clen
   ! return negative volume instead of imaginary one (means bad cell parameters)
   if (vol2 < 0.0_wp) dvol = -dvol ! this should not happen, but who knows...

   end associate
end function cell_to_dvol

!> @brief calculate the cell volume from the direct lattice
pure function dlat_to_dvol(lattice) result(dvol)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp),intent(in) :: lattice(3,3) !< direct lattice
   real(wp)            :: dvol

   dvol = abs(mat_det_3x3(lattice))

end function dlat_to_dvol

!> @brief inverts volume of direct unit cell
pure function dvol_to_rvol(dvol) result(rvol)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp),intent(in) :: dvol !< direct volume
   real(wp)            :: rvol !< reciprocal volume
   rvol = tpi**3/dvol
end function dvol_to_rvol

!> @brief transform from fractional coordinates into cartesian coordinates
pure subroutine abc_to_xyz(n,dlat,abc,xyz)
   use iso_fortran_env, wp => real64
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: dlat(3,3)
   real(wp),intent(in)  :: abc(3,n)
   real(wp),intent(out) :: xyz(3,n)

   call coord_trafo(n,dlat,abc,xyz)

end subroutine abc_to_xyz

!> @brief transform from cartesian coordinates into fractional coordinates
pure subroutine xyz_to_abc(n,dlat,xyz,abc,pbc)
   use iso_fortran_env, wp => real64
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: dlat(3,3)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(out) :: abc(3,n)
   logical, intent(in)  :: pbc(3)
   real(wp)             :: invlat(3,3)
   real(wp)             :: tmp(3)
   integer              :: iat,idir

   ! note that this is not the reciprocal lattice
   invlat = mat_inv_3x3(dlat)
   do iat = 1, n
      tmp = matmul(invlat,xyz(:,iat))
      do idir = 1, 3
         if (pbc(idir)) &
            tmp(idir) = shift_back_abc(tmp(idir))
      enddo
      abc(:,iat) = tmp
   enddo

end subroutine xyz_to_abc

!> @brief shift back fractional coordinates into the range [0,1)
pure elemental function shift_back_abc(in) result(out)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp),intent(in) :: in   !< fractional coordinate in (-∞,+∞)
   real(wp)            :: out  !< fractional coordinate in [0,1)
   real(wp),parameter  :: p_pbc_eps = 1.0e-14_wp
   out = in
   if(in < (0.0_wp - p_pbc_eps)) &
      out = in + real(ceiling(-in),wp)
   if(in > (1.0_wp + p_pbc_eps)) &
      out = in - real(floor  ( in),wp)
   if (abs(in - 1.0_wp) < p_pbc_eps) &
      out = in - 1.0_wp
end function shift_back_abc

!> @brief transform from one coordinate system into another
pure subroutine coord_trafo_12(dim,trafo,coord1,coord2)
   use iso_fortran_env, wp => real64
   implicit none
   integer, intent(in)  :: dim
   real(wp),intent(in)  :: trafo(3,3)
   real(wp),intent(in)  :: coord1(3,dim)
   real(wp),intent(out) :: coord2(3,dim)

   integer i

   coord2 = 0.0_wp
   do concurrent (i = 1:dim)
      coord2(:,i) = matmul(trafo,coord1(:,i))
   enddo

end subroutine coord_trafo_12

!> @brief transform from one coordinate system into another
pure subroutine coord_trafo_inplace(dim,trafo,coord)
   use iso_fortran_env, wp => real64
   implicit none
   integer, intent(in)    :: dim
   real(wp),intent(in)    :: trafo(3,3)
   real(wp),intent(inout) :: coord(3,dim)
   real(wp) :: tmp(3)

   integer i

   do concurrent (i = 1:dim)
      tmp = coord(:,i)
      coord(:,i) = matmul(trafo,tmp)
   enddo

end subroutine coord_trafo_inplace

!> @brief calculate distance between to atoms under minimum image convention
pure function minimum_image_distance(lsame,fi,fj,dlat,lpbc) result(dist)
   use iso_fortran_env, wp => real64
   implicit none
   logical, intent(in) :: lsame
   real(wp),intent(in) :: fi(3)
   real(wp),intent(in) :: fj(3)
   real(wp),intent(in) :: dlat(3,3)
   logical, intent(in) :: lpbc(3)
   real(wp) :: dist
   real(wp) :: fij(3),tfij(3),rij(3),tmp
   integer  :: rep(3),tx,ty,tz,idir
   logical  :: first

   fij = fi - fj

   if (lsame) then
      do idir = 1, 3
         if (lpbc(idir)) then
            rep(idir) = 1
         else
            rep(idir) = 0
         endif
      enddo
      first = .true.
      dist  = 0.0_wp
      tmp   = 0.0_wp

      do tx = -rep(1),rep(1)
         do ty = -rep(2),rep(2)
            do tz = -rep(3),rep(3)
               if (tx==0 .and. ty==0 .and. tz==0) cycle
               tfij = fij + [tx, ty, tz]

               rij = matmul(dlat,tfij)
               tmp = norm2(rij)

               if (first) then
                  dist = tmp
                  first = .false.
               else
                  dist = min(dist,tmp)
               endif

            enddo
         enddo
      enddo

   else
      do idir = 1, 3
         if (lpbc(idir)) fij(idir) = fij(idir) - idnint(fij(idir))
      enddo

      rij = matmul(dlat,fij)

      dist = norm2(rij)
   endif

end function minimum_image_distance

!> @brief calculate center of unit cell
pure function get_center_dlat(dlat) result(center)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp),intent(in) :: dlat(3,3)
   real(wp) :: center(3)
   real(wp) :: half(3)

   half = 0.5_wp

   center = matmul(dlat,half)

end function get_center_dlat

pure function outer_prod_3x3(a,b) result(c)
   real(wp),intent(in) :: a(3),b(3)
   real(wp) :: c(3,3)
   c(1,1) = a(1)*b(1); c(2,1) = a(2)*b(1); c(3,1) = a(3)*b(1)
   c(1,2) = a(1)*b(2); c(2,2) = a(2)*b(2); c(3,2) = a(3)*b(2)
   c(1,3) = a(1)*b(3); c(2,3) = a(2)*b(3); c(3,3) = a(3)*b(3)
end function outer_prod_3x3

pure subroutine latgrad_to_sigma(latgrad,lattice,sigma)
   implicit none
   real(wp), intent(in) :: latgrad(3,3)
   real(wp), intent(in) :: lattice(3,3)
   real(wp), intent(out) :: sigma(3,3)

   integer :: i,j,k

   sigma = 0.0_wp

   do i = 1, 3
      do j = 1, 3
         do k = 1,3
            sigma(i,j) = sigma(i,j) + latgrad(i,k)*lattice(j,k)
         enddo
      enddo
   enddo

end subroutine latgrad_to_sigma
pure subroutine sigma_to_latgrad(sigma,inv_lat,latgrad)
   implicit none
   real(wp), intent(in) :: sigma(3,3)
   real(wp), intent(in) :: inv_lat(3,3)
   real(wp), intent(out) :: latgrad(3,3)

   integer :: i,j,k

   latgrad(1,1) = sigma(1,1)*inv_lat(1,1) + sigma(1,2)*inv_lat(1,2) + sigma(1,3)*inv_lat(1,3)
   latgrad(1,2) = sigma(1,1)*inv_lat(2,1) + sigma(1,2)*inv_lat(2,2) + sigma(1,3)*inv_lat(2,3)
   latgrad(1,3) = sigma(1,1)*inv_lat(3,1) + sigma(1,2)*inv_lat(3,2) + sigma(1,3)*inv_lat(3,3)
   latgrad(2,1) = sigma(2,1)*inv_lat(1,1) + sigma(2,2)*inv_lat(1,2) + sigma(2,3)*inv_lat(1,3)
   latgrad(2,2) = sigma(2,1)*inv_lat(2,1) + sigma(2,2)*inv_lat(2,2) + sigma(2,3)*inv_lat(2,3)
   latgrad(2,3) = sigma(2,1)*inv_lat(3,1) + sigma(2,2)*inv_lat(3,2) + sigma(2,3)*inv_lat(3,3)
   latgrad(3,1) = sigma(3,1)*inv_lat(1,1) + sigma(3,2)*inv_lat(1,2) + sigma(3,3)*inv_lat(1,3)
   latgrad(3,2) = sigma(3,1)*inv_lat(2,1) + sigma(3,2)*inv_lat(2,2) + sigma(3,3)*inv_lat(2,3)
   latgrad(3,3) = sigma(3,1)*inv_lat(3,1) + sigma(3,2)*inv_lat(3,2) + sigma(3,3)*inv_lat(3,3)

end subroutine sigma_to_latgrad

end module pbc_tools
