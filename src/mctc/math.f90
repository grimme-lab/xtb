! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> Simple algebraic functions
module xtb_mctc_math
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_lapack, only : lapack_syev
   use xtb_mctc_constants, only : twothirdpi
   implicit none
   private

   public :: matDet3x3, matInv3x3, crossProd, derivDet3x3, eigval3x3, eigvec3x3


contains


!> Performs a direct calculation of the inverse of a 3×3 matrix.
!
!  reference: http://fortranwiki.org/fortran/show/Matrix+inversion
pure function matInv3x3(a) result(b)

   !> Matrix
   real(wp), intent(in) :: a(3, 3)

   !> Inverse matrix
   real(wp) :: b(3, 3)

   real(wp) :: detinv

   ! Calculate the inverse determinant of the matrix
   detinv = 1.0_wp/matDet3x3(a)

   ! Calculate the inverse of the matrix
   b(1, 1) = +detinv * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2))
   b(2, 1) = -detinv * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1))
   b(3, 1) = +detinv * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
   b(1, 2) = -detinv * (a(1, 2) * a(3, 3) - a(1, 3) * a(3, 2))
   b(2, 2) = +detinv * (a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1))
   b(3, 2) = -detinv * (a(1, 1) * a(3, 2) - a(1, 2) * a(3, 1))
   b(1, 3) = +detinv * (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2))
   b(2, 3) = -detinv * (a(1, 1) * a(2, 3) - a(1, 3) * a(2, 1))
   b(3, 3) = +detinv * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1))

end function matInv3x3


!> Determinat of 3×3 matrix
pure function matDet3x3(a) result (det)

   !> Matrix
   real(wp), intent(in) :: a(3, 3)

   !> Determinant
   real(wp) :: det

   det =  a(1, 1) * a(2, 2) * a(3, 3)  &
      & - a(1, 1) * a(2, 3) * a(3, 2)  &
      & - a(1, 2) * a(2, 1) * a(3, 3)  &
      & + a(1, 2) * a(2, 3) * a(3, 1)  &
      & + a(1, 3) * a(2, 1) * a(3, 2)  &
      & - a(1, 3) * a(2, 2) * a(3, 1)

end function matDet3x3


!> Derivative of determinant of a 3x3 matrix
subroutine derivDet3x3(deriv, mat)

   !> derivative of the determinant
   real(wp), intent(out) :: deriv(3, 3)

   !> The matrix from which to calculate the determinant.
   real(wp), intent(in) :: mat(3, 3)

   deriv(1, 1) = mat(2, 2) * mat(3, 3) - mat(3, 2) * mat(2, 3)
   deriv(1, 2) = mat(2, 3) * mat(3, 1) - mat(3, 3) * mat(2, 1)
   deriv(1, 3) = mat(2, 1) * mat(3, 2) - mat(3, 1) * mat(2, 2)
   deriv(2, 1) = mat(1, 3) * mat(3, 2) - mat(1, 2) * mat(3, 3)
   deriv(2, 2) = mat(1, 1) * mat(3, 3) - mat(1, 3) * mat(3, 1)
   deriv(2, 3) = mat(1, 2) * mat(3, 1) - mat(1, 1) * mat(3, 2)
   deriv(3, 1) = mat(1, 2) * mat(2, 3) - mat(1, 3) * mat(2, 2)
   deriv(3, 2) = mat(1, 3) * mat(2, 1) - mat(1, 1) * mat(2, 3)
   deriv(3, 3) = mat(1, 1) * mat(2, 2) - mat(1, 2) * mat(2, 1)

   deriv = deriv * sign(1.0_wp, matDet3x3(mat))

end subroutine derivDet3x3


!> Implements the cross/vector product between two 3D vectors
pure function crossProd(a,b) result(c)

   !> First vector
   real(wp), intent(in) :: a(3)

   !> Second vector
   real(wp), intent(in) :: b(3)

   !> Orthogonal vector
   real(wp) :: c(3)

   c(1) = a(2) * b(3) - b(2) * a(3)
   c(2) = a(3) * b(1) - b(3) * a(1)
   c(3) = a(1) * b(2) - b(1) * a(2)

end function crossProd


!> Calculates eigenvalues based on the trigonometric solution of A = pB + qI
pure subroutine eigval3x3(a, w)

   !> The symmetric input matrix
   real(wp), intent(in) :: a(3, 3)

   !> Contains eigenvalues on exit
   real(wp), intent(out) :: w(3)

   real(wp) :: q, p, r

   r = a(1, 2) * a(1, 2) + a(1, 3) * a(1, 3) + a(2, 3) * a(2, 3)
   q = (a(1, 1) + a(2, 2) + a(3, 3)) / 3.0_wp
   w(1) = a(1, 1) - q
   w(2) = a(2, 2) - q
   w(3) = a(3, 3) - q
   p = sqrt((w(1) * w(1) + w(2) * w(2) + w(3) * w(3) + 2*r) / 6.0_wp)
   r = (w(1) * (w(2) * w(3) - a(2, 3) * a(2, 3)) &
      & - a(1, 2) * (a(1, 2) * w(3) - a(2, 3) * a(1, 3)) &
      & + a(1, 3) * (a(1, 2) * a(2, 3) - w(2) * a(1, 3))) / (p*p*p) * 0.5_wp

   if (r <= -1.0_wp) then
      r = 0.5_wp * twothirdpi
   else if (r >= 1.0_wp) then
      r = 0.0_wp
   else
      r = acos(r) / 3.0_wp
   end if

   w(3) = q + 2 * p * cos(r)
   w(1) = q + 2 * p * cos(r + twothirdpi)
   w(2) = 3 * q - w(1) - w(3)

end subroutine eigval3x3


!> Calculates eigenvector using an analytical method based on vector cross
!  products.
pure subroutine eigvec3x3(a, w, q)
   ! .. Arguments ..
   real(wp), intent(inout) :: a(3,3)
   real(wp), intent(out) :: w(3)
   real(wp), intent(out) :: q(3,3)

   real(wp), parameter :: eps = epsilon(1.0_wp)
   ! .. Local Variables ..
   real(wp) norm, n1, n2, n3, precon
   integer   i, j

   w(1) = max(abs(a(1, 1)), abs(a(1, 2)))
   w(2) = max(abs(a(1, 3)), abs(a(2, 2)))
   w(3) = max(abs(a(2, 3)), abs(a(3, 3)))
   precon = max(w(1), max(w(2), w(3)))

   ! null matrix
   if (precon < eps) then
      w(1) = 0.0_wp
      w(2) = 0.0_wp
      w(3) = 0.0_wp
      q(1, 1) = 1.0_wp
      q(2, 2) = 1.0_wp
      q(3, 3) = 1.0_wp
      q(1, 2) = 0.0_wp
      q(1, 3) = 0.0_wp
      q(2, 3) = 0.0_wp
      q(2, 1) = 0.0_wp
      q(3, 1) = 0.0_wp
      q(3, 2) = 0.0_wp
      return
   end if

   norm = 1.0_wp / precon

   a(1, 1) = a(1, 1) * norm
   a(1, 2) = a(1, 2) * norm
   a(2, 2) = a(2, 2) * norm
   a(1, 3) = a(1, 3) * norm
   a(2, 3) = a(2, 3) * norm
   a(3, 3) = a(3, 3) * norm

   ! Calculate eigenvalues
   call eigval3x3(a, w)

   ! Compute first eigenvector
   a(1, 1) = a(1, 1) - w(1)
   a(2, 2) = a(2, 2) - w(1)
   a(3, 3) = a(3, 3) - w(1)

   q(1, 1) = a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2)
   q(2, 1) = a(1, 3) * a(1, 2) - a(1, 1) * a(2, 3)
   q(3, 1) = a(1, 1) * a(2, 2) - a(1, 2) * a(1, 2)
   q(1, 2) = a(1, 2) * a(3, 3) - a(1, 3) * a(2, 3)
   q(2, 2) = a(1, 3) * a(1, 3) - a(1, 1) * a(3, 3)
   q(3, 2) = a(1, 1) * a(2, 3) - a(1, 2) * a(1, 3)
   q(1, 3) = a(2, 2) * a(3, 3) - a(2, 3) * a(2, 3)
   q(2, 3) = a(2, 3) * a(1, 3) - a(1, 2) * a(3, 3)
   q(3, 3) = a(1, 2) * a(2, 3) - a(2, 2) * a(1, 3)
   n1 = q(1, 1) * q(1, 1) + q(2, 1) * q(2, 1) + q(3, 1) * q(3, 1)
   n2 = q(1, 2) * q(1, 2) + q(2, 2) * q(2, 2) + q(3, 2) * q(3, 2)
   n3 = q(1, 3) * q(1, 3) + q(2, 3) * q(2, 3) + q(3, 3) * q(3, 3)

   norm = n1
   i = 1
   if (n2 > norm) then
      i = 2
      norm = n1
   end if
   if (n3 > norm) then
      i = 3
   end if

   if (i == 1) then
      norm = sqrt(1.0_wp / n1)
      q(1, 1) = q(1, 1) * norm
      q(2, 1) = q(2, 1) * norm
      q(3, 1) = q(3, 1) * norm
   else if (i == 2) then
      norm = sqrt(1.0_wp / n2)
      q(1, 1) = q(1, 2) * norm
      q(2, 1) = q(2, 2) * norm
      q(3, 1) = q(3, 2) * norm
   else
      norm = sqrt(1.0_wp / n3)
      q(1, 1) = q(1, 3) * norm
      q(2, 1) = q(2, 3) * norm
      q(3, 1) = q(3, 3) * norm
   end if

   ! Robustly compute a right-hand orthonormal set (ev1, u, v)
   if (abs(q(1, 1)) > abs(q(2, 1))) then
      norm = sqrt(1.0_wp / (q(1, 1) * q(1, 1) + q(3, 1) * q(3, 1)))
      q(1, 2) = -q(3, 1) * norm
      q(2, 2) = 0.0_wp
      q(3, 2) = +q(1, 1) * norm
   else
      norm = sqrt(1.0_wp / (q(2, 1) * q(2, 1) + q(3, 1) * q(3, 1)))
      q(1, 2) = 0.0_wp
      q(2, 2) = +q(3, 1) * norm
      q(3, 2) = -q(2, 1) * norm
   end if
   q(1, 3) = q(2, 1) * q(3, 2) - q(3, 1) * q(2, 2)
   q(2, 3) = q(3, 1) * q(1, 2) - q(1, 1) * q(3, 2)
   q(3, 3) = q(1, 1) * q(2, 2) - q(2, 1) * q(1, 2)

   ! Reset A
   a(1, 1) = a(1, 1) + w(1)
   a(2, 2) = a(2, 2) + w(1)
   a(3, 3) = a(3, 3) + w(1)

   ! A*U
   n1 = a(1, 1) * q(1, 2) + a(1, 2) * q(2, 2) + a(1, 3) * q(3, 2)
   n2 = a(1, 2) * q(1, 2) + a(2, 2) * q(2, 2) + a(2, 3) * q(3, 2)
   n3 = a(1, 3) * q(1, 2) + a(2, 3) * q(2, 2) + a(3, 3) * q(3, 2)

   ! A*V, note out of order computation
   a(3, 3) = a(1, 3) * q(1, 3) + a(2, 3) * q(2, 3) + a(3, 3) * q(3, 3)
   a(1, 3) = a(1, 1) * q(1, 3) + a(1, 2) * q(2, 3) + a(1, 3) * q(3, 3)
   a(2, 3) = a(1, 2) * q(1, 3) + a(2, 2) * q(2, 3) + a(2, 3) * q(3, 3)

   ! UT*(A*U) - l2*E
   n1 = q(1, 2) * n1 + q(2, 2) * n2 + q(3, 2) * n3 - w(2)
   ! UT*(A*V)
   n2 = q(1, 2) * a(1, 3) + q(2, 2) * a(2, 3) + q(3, 2) * a(3, 3)
   ! VT*(A*V) - l2*E
   n3 = q(1, 3) * a(1, 3) + q(2, 3) * a(2, 3) + q(3, 3) * a(3, 3) - w(2)

   if (abs(n1) >= abs(n3)) then
      norm = max(abs(n1), abs(n2))
      if (norm > eps) then
         if (abs(n1) >= abs(n2)) then
            n2 = n2 / n1
            n1 = sqrt(1.0_wp / (1.0_wp + n2 * n2))
            n2 = n2 * n1
         else
            n1 = n1 / n2
            n2 = sqrt(1.0_wp / (1.0_wp + n1 * n1))
            n1 = n1 * n2
         end if
         q(1, 2) = n2 * q(1, 2) - n1 * q(1, 3)
         q(2, 2) = n2 * q(2, 2) - n1 * q(2, 3)
         q(3, 2) = n2 * q(3, 2) - n1 * q(3, 3)
      end if
   else
      norm = max(abs(n3), abs(n2))
      if (norm > eps) then
         if (abs(n3) >= abs(n2)) then
            n2 = n2 / n3
            n3 = sqrt(1.0_wp / (1.0_wp + n2 * n2))
            n2 = n2 * n3
         else
            n3 = n3 / n2
            n2 = sqrt(1.0_wp / (1.0_wp + n3 * n3))
            n3 = n3 * n2
         end if
         q(1, 2) = n3 * q(1, 2) - n2 * q(1, 3)
         q(2, 2) = n3 * q(2, 2) - n2 * q(2, 3)
         q(3, 2) = n3 * q(3, 2) - n2 * q(3, 3)
      end if
   end if

   ! Calculate third eigenvector from cross product
   q(1, 3) = q(2, 1) * q(3, 2) - q(3, 1) * q(2, 2)
   q(2, 3) = q(3, 1) * q(1, 2) - q(1, 1) * q(3, 2)
   q(3, 3) = q(1, 1) * q(2, 2) - q(2, 1) * q(1, 2)

   w(1) = w(1) * precon
   w(2) = w(2) * precon
   w(3) = w(3) * precon

end subroutine eigvec3x3


end module xtb_mctc_math
