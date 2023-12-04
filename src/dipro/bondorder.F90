! This file is part of dipro.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

!> Calculation of Wiberg/Mayer type bond orders
module xtb_dipro_bondorder
   use mctc_env, only : wp
#if WITH_TBLITE
   use tblite_basis_type, only : basis_type
   use tblite_blas, only : gemm
   implicit none

   public :: get_wiberg_bondorder

contains

!> Calculate Wiberg/Mayer type bond orders from density and overlap matrices
subroutine get_wiberg_bondorder(bas, smat, pmat, wbo)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap matrix
   real(wp), intent(in) :: smat(:, :)
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :)
   !> Wiberg/Mayer type bond orders
   real(wp), intent(out) :: wbo(:, :)

   real(wp), allocatable :: spmat(:, :)
   integer :: iao, jao, iat, jat

   allocate(spmat(bas%nao,bas%nao)) !XXXX  (mold=pmat)
   call gemm(pmat, smat, spmat)

   wbo(:, :) = 0.0_wp
   do iao = 1, bas%nao
      iat = bas%ao2at(iao)
      do jao = 1, bas%nao
         jat = bas%ao2at(jao)
         wbo(jat, iat) = wbo(jat, iat) + spmat(iao, jao) * spmat(jao, iao)
      end do
   end do

   do iat = 1, size(wbo, 2)
      wbo(iat, iat) = 0.0_wp
   end do

end subroutine get_wiberg_bondorder
#endif

end module xtb_dipro_bondorder
