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
   implicit none
   private

   public :: matDet3x3, matInv3x3, crossProd, derivDet3x3


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


end module xtb_mctc_math
