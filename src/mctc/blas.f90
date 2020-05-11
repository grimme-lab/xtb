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

!> Interfaces to BLAS
module xtb_mctc_blas
   use xtb_mctc_blas_level1, only : blas_asum, blas_axpy, blas_copy, blas_dot, &
      & blas_dotc, blas_dotu, blas_nrm2, blas_rot, blas_scal, blas_swap, &
      & blas_iamax, blas_cabs1
   use xtb_mctc_blas_level2, only : blas_gbmv, blas_gemv, blas_ger, blas_gerc, &
      & blas_geru, blas_hbmv, blas_hemv, blas_her, blas_her2, blas_hpmv, &
      & blas_hpr, blas_hpr2, blas_sbmv, blas_spmv, blas_spr, blas_spr2, &
      & blas_symv, blas_syr, blas_syr2, blas_tbmv, blas_tbsv, blas_tpmv, &
      & blas_tpsv, blas_trmv, blas_trsv
   use xtb_mctc_blas_level3, only : blas_gemm, blas_hemm, blas_herk, blas_her2k, &
      & blas_symm, blas_syrk, blas_syr2k, blas_trmm, blas_trsm
   use xtb_mctc_blas_wrap1
   use xtb_mctc_blas_wrap2
   use xtb_mctc_blas_wrap3
   implicit none
   public

end module xtb_mctc_blas
