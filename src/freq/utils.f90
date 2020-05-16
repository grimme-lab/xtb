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

!> Utililty functions for dealing with hessian and frequency calculations
module xtb_freq_utils
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_lapack, only : lapack_syev
   use xtb_type_environment, only : TEnvironment
   implicit none
   private

   public :: massWeightHessian, diagHessian


contains


subroutine massWeightHessian(hessian, atmass)
   real(wp), intent(inout) :: hessian(:, :)
   real(wp), intent(in) :: atmass(:)

   integer :: ndim, nat
   integer :: iat, ii, jj, ic
   real(wp), allocatable :: isqm(:)

   ndim = size(hessian, 1)
   nat = size(atmass)

   allocate(isqm(ndim))

   do iat = 1, nat
      ic = (iat - 1)*3
      isqm(ic+1:ic+3) = sqrt(1.0_wp/atmass(iat))
   end do

   do ii = 1, ndim
      do jj = 1, ndim
         hessian(jj, ii) = hessian(jj, ii) * isqm(ii) * isqm(jj)
      end do
   end do

end subroutine massWeightHessian


subroutine diagHessian(env, hessian, freq)
   character(len=*), parameter :: source = 'freq_utils_diagHessian'
   type(TEnvironment), intent(inout) :: env
   real(wp), intent(inout) :: hessian(:, :)
   real(wp), intent(out) :: freq(:)
   real(wp), allocatable :: work(:)
   integer :: ndim, lwork, info, ii, nn, zero(6)
   ndim = size(hessian, 1)
   lwork = 1 + 6*ndim + 2*ndim**2
   allocate(work(lwork))
   call lapack_syev('v', 'u', ndim, hessian, ndim, freq, work, lwork, info)

   if (info /= 0) then
      call env%error("Could not diagonalize Hessian", source)
      return
   end if

   do ii = 1, ndim
      freq(ii) = sign(sqrt(abs(freq(ii))), freq(ii))
   end do

   nn = 0
   do ii = 1, ndim
      if (abs(freq(ii)) < 1.0e-2_wp .and. nn < size(zero)) then
         nn = nn + 1
         zero(nn) = ii
      end if
   end do

end subroutine diagHessian


end module xtb_freq_utils
