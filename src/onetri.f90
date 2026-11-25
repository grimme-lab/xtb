! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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
module xtb_onetri
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas, only : mctc_symm, mctc_gemm
   use xtb_blowsy
   implicit none
   private

   public :: onetri

contains

subroutine onetri(ity,s,s1,array,n,ival)
   !     ******designed for abelian groups only******
   !
   !     calling sequence:
   !     ity       =1 sym ao
   !               =-1 anti sym ao
   !     s         input property matrix over ao's
   !     s1        transformed integrals on output
   !     s2        scratch arrays large enough to hold a square matrix
   !     array     mo matrix over so's
   !     n         linear dimension of arrays
   !     bernd hess, university of bonn, january 1991
   integer, intent(in) :: n
   integer, intent(in) :: ity
   integer, intent(in) :: ival
   real(wp), intent(in) :: s(*)
   real(wp), intent(inout), target :: s1(*)
   real(wp), intent(in) :: array(n,ival)
   real(wp) :: s2(n,n)
   real(wp), pointer :: s1_ao(:, :), s1_mo(:, :)

   s1_ao(1:n, 1:n) => s1(1:n*n)
   s1_mo(1:ival, 1:ival) => s1(1:ival*ival)

   !
   !     determine if we have an antisymmetric integral
   if (ity /= -1) then
      !
      !     blow up symmetric matrix s
      call blowsy(ity,s,s1,n)

      !
      !     transformation of s
      call mctc_symm(s1_ao, array, s2, side='l', uplo='l', m=n, n=ival)
      call mctc_gemm(array, s2, s1_mo, transa='t', m=ival, n=ival, k=n)
   else
      !
      !     blow up anti-symmetric matrix s
      call blowsy(ity,s,s1,n)
      !
      !     transformation of s
      call mctc_gemm(s1_ao, array, s2, m=n, n=ival, k=n)
      call mctc_gemm(array, s2, s1_mo, transa='t', m=ival, n=ival, k=n)
   end if
end subroutine

end module xtb_onetri
