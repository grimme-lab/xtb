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
! You should have received a copy of the GNU Lesser General Public Licen
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.
module xtb_onetri
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas, only : blas_symm, blas_gemm
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
   real(wp), intent(inout) :: s1(*)
   real(wp), intent(in) :: array(n,ival)
   real(wp) :: s2(n,n)
   external :: dsymm, dgemm ! can't use wrappers due to change of leading dim

   !
   !     determine if we have an antisymmetric integral
   if (ity /= -1) then
      !
      !     blow up symmetric matrix s
      call blowsy(ity,s,s1,n)

      !
      !     transformation of s
      call dsymm('l','l',n,ival,1.d0,s1,n,array,n,0.d0,s2,n)
      call dgemm('t','n',ival,ival,n,1.d0,array,n,s2,n,0.d0,s1,ival)
   else
      !
      !     blow up anti-symmetric matrix s
      call blowsy(ity,s,s1,n)
      !
      !     transformation of s
      call dgemm('n','n',n,ival,n,1.d0,s1,n,array,n,0.d0,s2,n)
      call dgemm('t','n',ival,ival,n,1.d0,array,n,s2,n,0.d0,s1,ival)
   end if
end subroutine

end module xtb_onetri
