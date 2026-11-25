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

!**********************************************************************
!                                                                     *
!  make the loewdin orthogonalization matrix x = u' s1/2 u            *
!  and the Loewdin orbitals cl from canonical ca
!                                                                     *
!**********************************************************************

subroutine makel(nao, s, can, clo)
   use xtb_mctc_lapack, only: mctc_syev
   use xtb_mctc_blas, only: mctc_gemm
   implicit integer (i - n)
   implicit real * 8(a - h, o - z)
   dimension s(nao, nao)
   dimension can(nao, nao)
   dimension clo(nao, nao)
   real*8, allocatable :: vecs(:, :), e(:), cc(:, :), x(:, :)

   allocate (vecs(nao, nao), e(nao), cc(nao, nao))
   allocate (x(nao, nao))

   vecs = s

   call mctc_syev(vecs, e, jobz='V', uplo='U', info=info)

   do i = 1, nao
      if (e(i) < 0) stop 'sorry, must stop in S^1/2!'
      e(i) = dsqrt(e(i))
   end do

   do m = 1, nao
      do i = 1, nao
         x(i, m) = vecs(i, m)
         cc(i, m) = e(m) * vecs(i, m)
      end do
   end do

   call mctc_gemm(x, cc, vecs, transb='T')

   x = vecs
   deallocate (e, cc, vecs)

   call mctc_gemm(x, can, clo)

   deallocate (x)
   return
end subroutine makel

! unrestricted version
subroutine umakel(nao, s, cana, canb, cloa, clob)
   use xtb_mctc_lapack, only: mctc_syev
   use xtb_mctc_blas, only: mctc_gemm
   implicit integer (i-n)
   implicit real * 8(a - h, o - z)
   dimension s(nao, nao)
   dimension cana(nao, nao)
   dimension cloa(nao, nao)
   dimension canb(nao, nao)
   dimension clob(nao, nao)
   real*8, allocatable :: vecs(:, :), e(:), cc(:, :), x(:, :)

   allocate (vecs(nao, nao), e(nao), cc(nao, nao))
   allocate (x(nao, nao))

   vecs = s

   call mctc_syev(vecs, e, jobz='V', uplo='U', info=info)

   do i = 1, nao
      if (e(i) < 0) stop 'sorry, must stop in S^1/2!'
      e(i) = dsqrt(e(i))
   end do

   do m = 1, nao
      do i = 1, nao
         x(i, m) = vecs(i, m)
         cc(i, m) = e(m) * vecs(i, m)
      end do
   end do

   call mctc_gemm(x, cc, vecs, transb='T')

   x = vecs
   deallocate (e, cc, vecs)

   call mctc_gemm(x, cana, cloa)
   call mctc_gemm(x, canb, clob)

   deallocate (x)
   return
end subroutine umakel
