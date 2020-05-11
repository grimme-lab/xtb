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

!> Interfaces to LAPACK
module xtb_mctc_lapack
   use xtb_mctc_lapack_geneigval, only : lapack_sygv, lapack_sygvd, lapack_sygvx, &
      & lapack_spgv, lapack_spgvd, lapack_spgvx,lapack_hegv, lapack_hegvd, &
      & lapack_hegvx, lapack_hpgv, lapack_hpgvd,  lapack_hpgvx
   use xtb_mctc_lapack_gst, only : lapack_sygst, lapack_spgst, lapack_hegst, &
      & lapack_hpgst
   use xtb_mctc_lapack_stdeigval, only : lapack_syev, lapack_syevd, lapack_syevx, &
      & lapack_syevr, lapack_spev, lapack_spevd, lapack_spevx, lapack_heev, &
      & lapack_heevd, lapack_heevx, lapack_heevr, lapack_hpev, lapack_hpevd, &
      & lapack_hpevx
   use xtb_mctc_lapack_trf, only : lapack_getrf, lapack_sytrf, lapack_sptrf, &
      & lapack_potrf, lapack_pptrf
   use xtb_mctc_lapack_tri, only : lapack_getri, lapack_sytri, lapack_sptri, &
      & lapack_potri, lapack_pptri
   use xtb_mctc_lapack_trs, only : lapack_getrs, lapack_sytrs, lapack_sptrs, &
      & lapack_potrs, lapack_pptrs
   use xtb_mctc_lapack_wrap
   implicit none
   public


end module xtb_mctc_lapack
