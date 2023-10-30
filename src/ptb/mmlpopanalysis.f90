! This file is part of xtb.
!
! Copyright (C) 2023 Marcel Mueller
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

!> Module for the mixed Mulliken-Loewdin population analysis

module xtb_ptb_mmlpopanalysis
   use mctc_env, only: wp

   use xtb_mctc_lapack, only: lapack_syev

   use tblite_basis_type, only: basis_type
   use tblite_blas, only: gemm
   implicit none
   private

   public :: get_mml_overlaps

contains

   subroutine get_mml_overlaps(bas, overlap, ratio, sx, soneminusx)
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> Overlap matrix
      real(wp), intent(in) :: overlap(:, :)
      !> Ratio of Mulliken to Loewdin population
      real(wp), intent(in) :: ratio
      !> Overlap matrix for Mulliken-Lowdin population analysis
      real(wp), intent(out) :: sx(bas%nao, bas%nao), soneminusx(bas%nao, bas%nao)
      !> Variables for ML-pop:
      !> Eigenvalues of overlap matrix (exponentiated)
      real(wp), allocatable :: seig(:), seig1(:), seig2(:)
      !> Auxiliary variables for diagonalization
      integer :: lwork, info, i, j
      real(wp), allocatable :: aux(:)
      !> Temporary overlap matrix as a result of the diagonalization
      real(wp), allocatable :: tmps(:, :)
      !> Temporary overlap matrix for the matrix multiplication
      real(wp) :: tmp2s(bas%nao, bas%nao)

      allocate (tmps(bas%nao, bas%nao), seig(bas%nao), &
      & seig1(bas%nao), seig2(bas%nao))
      tmps = overlap

      allocate (aux(1))
      lwork = -1
      call lapack_syev('V', 'U', bas%nao, tmps, bas%nao, seig, aux, lwork, info)
      lwork = int(aux(1))
      deallocate (aux)
      allocate (aux(lwork))
      call lapack_syev('V', 'U', bas%nao, tmps, bas%nao, seig, aux, lwork, info)
      ! write (*, *) "Eigenvalues:"
      do i = 1, bas%nao
         !write (*, '(a,f10.6)') "Eigenvalues: ", seig(i)
         seig1(i) = seig(i)**(1.0_wp - ratio)
         seig2(i) = seig(i)**ratio
         ! write (*, '(a,f10.6)') "Eigenvalues (1-x): ", seig1(i)
         ! write (*, '(a,f10.6)') "Eigenvalues (x): ", seig2(i)
      end do

      !##### DEV WRITE #####
      ! write (*, *) "Overlap tmpS:"
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.5)', advance="no") tmps(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################
      do i = 1, bas%nao
         do j = 1, bas%nao
            tmp2s(j, i) = tmps(j, i)*seig1(i)
         end do
      end do
      !##### DEV WRITE #####
      ! write (*, *) "Overlap tmp2S:"
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") tmp2s(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################
      call gemm(tmps, tmp2s, soneminusx, 'N', 'T')

      do i = 1, bas%nao
         do j = 1, bas%nao
            tmp2s(j, i) = tmps(j, i)*seig2(i)
         end do
      end do
      call gemm(tmps, tmp2s, sx, 'N', 'T')

   end subroutine get_mml_overlaps

end module xtb_ptb_mmlpopanalysis
