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

!> Generalized Born interaction kernels
module xtb_solv_kernel
   implicit none
   public :: gbKernel
   private


   !> Possible kernels for the generalized Born model
   type :: TGBKernelEnum

      !> Classical Still kernel
      integer :: still = 1

      !> P16 kernel by Lange (JCTC 2012, 8, 1999-2011)
      integer :: p16 = 2

   end type TGBKernelEnum

   !> Actual enumerator for the generalized Born kernels
   type(TGBKernelEnum), parameter :: gbKernel = TGBKernelEnum()


contains


end module xtb_solv_kernel
