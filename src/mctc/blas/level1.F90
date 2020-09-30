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
module xtb_mctc_blas_level1
   use xtb_mctc_accuracy, only : sp, dp
#ifdef USE_CUBLAS
   use cublas, only : cublasSasum, cublasDasum, cublasSaxpy, cublasDaxpy, &
      & cublasScopy, cublasDcopy, cublasSdot, cublasDdot, &
      & cublasSnrm2, cublasDnrm2, cublasSrot, cublasDrot, &
      & cublasSscal, cublasDscal, cublasSswap, cublasDswap, &
      & cublasSiamax, cublasDiamax
#endif
   implicit none
   private

   public :: mctc_asum, mctc_axpy, mctc_copy, mctc_dot, mctc_nrm2, mctc_rot
   public :: mctc_scal, mctc_swap, mctc_iamax

   public :: blas_asum, blas_axpy, blas_copy, blas_dot, blas_dotc, blas_dotu
   public :: blas_nrm2, blas_rot, blas_scal, blas_swap, blas_iamax, blas_cabs1


   !> Takes the sum of the absolute values.
   interface mctc_asum
      module procedure :: mctc_sasum
      module procedure :: mctc_dasum
   end interface mctc_asum

   !> Constant times a vector plus a vector.
   interface mctc_axpy
      module procedure :: mctc_saxpy
      module procedure :: mctc_daxpy
   end interface mctc_axpy

   !> Copies a vector, x, to a vector, y.
   interface mctc_copy
      module procedure :: mctc_scopy
      module procedure :: mctc_dcopy
   end interface mctc_copy

   !> Forms the dot product of two vectors.
   interface mctc_dot
      module procedure :: mctc_sdot
      module procedure :: mctc_ddot
   end interface mctc_dot

   !> Returns the euclidean norm of a vector via the function
   !  name, so that
   !
   !    NRM2 := sqrt( x'*x )
   interface mctc_nrm2
      module procedure :: mctc_snrm2
      module procedure :: mctc_dnrm2
   end interface mctc_nrm2

   !> Applies a plane rotation.
   interface mctc_rot
      module procedure :: mctc_srot
      module procedure :: mctc_drot
   end interface mctc_rot

   !> Scales a vector by a constant.
   interface mctc_scal
      module procedure :: mctc_sscal
      module procedure :: mctc_dscal
   end interface mctc_scal

   !> Interchanges two vectors.
   interface mctc_swap
      module procedure :: mctc_sswap
      module procedure :: mctc_dswap
   end interface mctc_swap

   !> Finds the index of the first element having maximum absolute value.
   interface mctc_iamax
      module procedure :: mctc_isamax
      module procedure :: mctc_idamax
   end interface mctc_iamax


   !> Takes the sum of the absolute values.
   interface blas_asum
      pure function sasum(n, x, incx)
         import :: sp
         real(sp) :: sasum
         real(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function sasum
      pure function scasum(n, x, incx)
         import :: sp
         real(sp) :: scasum
         complex(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function scasum
      pure function dasum(n, x, incx)
         import :: dp
         real(dp) :: dasum
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function dasum
      pure function dzasum(n, x, incx)
         import :: dp
         real(dp) :: dzasum
         complex(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function dzasum
   end interface blas_asum

   !> Constant times a vector plus a vector.
   !  Uses unrolled loops for increments equal to one.
   interface blas_axpy
      pure subroutine saxpy(n, a, x, incx, y, incy)
         import :: sp
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         real(sp), intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine saxpy
      pure subroutine daxpy(n, a, x, incx, y, incy)
         import :: dp
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         real(dp), intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine daxpy
      pure subroutine caxpy(n, a, x, incx, y, incy)
         import :: sp
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(inout) :: y(*)
         complex(sp), intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine caxpy
      pure subroutine zaxpy(n, a, x, incx, y, incy)
         import :: dp
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(inout) :: y(*)
         complex(dp), intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine zaxpy
   end interface blas_axpy

   !> Copies a vector, x, to a vector, y.
   !  Uses unrolled loops for increments equal to 1.
   interface blas_copy
      pure subroutine scopy(n, x, incx, y, incy)
         import :: sp
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine scopy
      pure subroutine dcopy(n, x, incx, y, incy)
         import :: dp
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine dcopy
      pure subroutine ccopy(n, x, incx, y, incy)
         import :: sp
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine ccopy
      pure subroutine zcopy(n, x, incx, y, incy)
         import :: dp
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine zcopy
   end interface blas_copy

   !> Forms the dot product of two vectors.
   !  Uses unrolled loops for increments equal to one.
   interface blas_dot
      pure function sdot(n, x, incx, y, incy)
         import :: sp
         real(sp) :: sdot
         real(sp), intent(in) :: x(*)
         real(sp), intent(in) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end function sdot
      pure function ddot(n, x, incx, y, incy)
         import :: dp
         real(dp) :: ddot
         real(dp), intent(in) :: x(*)
         real(dp), intent(in) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end function ddot
   end interface blas_dot

   !> Forms the dot product of two complex vectors
   !    DOTC = X^H * Y
   interface blas_dotc
      pure function cdotc(n, x, incx, y, incy)
         import :: sp
         complex(sp) :: cdotc
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(in) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end function cdotc
      pure function zdotc(n, x, incx, y, incy)
         import :: dp
         complex(dp) :: zdotc
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(in) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end function zdotc
   end interface blas_dotc

   !> Forms the dot product of two complex vectors
   !    DOTU = X^T * Y
   interface blas_dotu
      pure function cdotu(n, x, incx, y, incy)
         import :: sp
         complex(sp) :: cdotu
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(in) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end function cdotu
      pure function zdotu(n, x, incx, y, incy)
         import :: dp
         complex(dp) :: zdotu
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(in) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end function zdotu
   end interface blas_dotu

   !> Returns the euclidean norm of a vector via the function
   !  name, so that
   !
   !    NRM2 := sqrt( x'*x )
   interface blas_nrm2
      pure function snrm2(n, x, incx)
         import :: sp
         real(sp) :: snrm2
         real(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function snrm2
      pure function dnrm2(n, x, incx)
         import :: dp
         real(dp) :: dnrm2
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function dnrm2
      pure function scnrm2(n, x, incx)
         import :: sp
         real(sp) :: scnrm2
         complex(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function scnrm2
      pure function dznrm2(n, x, incx)
         import :: dp
         real(dp) :: dznrm2
         complex(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function dznrm2
   end interface blas_nrm2

   !> Applies a plane rotation.
   interface blas_rot
      pure subroutine srot(n, x, incx, y, incy, c, s)
         import :: sp
         real(sp), intent(inout) :: x(*)
         real(sp), intent(inout) :: y(*)
         real(sp), intent(in) :: c
         real(sp), intent(in) :: s
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine srot
      pure subroutine drot(n, x, incx, y, incy, c, s)
         import :: dp
         real(dp), intent(inout) :: x(*)
         real(dp), intent(inout) :: y(*)
         real(dp), intent(in) :: c
         real(dp), intent(in) :: s
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine drot
      pure subroutine csrot(n, x, incx, y, incy, c, s)
         import :: sp
         complex(sp), intent(inout) :: x(*)
         complex(sp), intent(inout) :: y(*)
         real(sp), intent(in) :: c
         real(sp), intent(in) :: s
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine csrot
      pure subroutine zdrot(n, x, incx, y, incy, c, s)
         import :: dp
         complex(dp), intent(inout) :: x(*)
         complex(dp), intent(inout) :: y(*)
         real(dp), intent(in) :: c
         real(dp), intent(in) :: s
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine zdrot
   end interface blas_rot

   !> Scales a vector by a constant.
   !  Uses unrolled loops for increment equal to 1.
   interface blas_scal
      pure subroutine sscal(n, a, x, incx)
         import :: sp
         real(sp), intent(inout) :: x(*)
         real(sp), intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine sscal
      pure subroutine dscal(n, a, x, incx)
         import :: dp
         real(dp), intent(inout) :: x(*)
         real(dp), intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine dscal
      pure subroutine cscal(n, a, x, incx)
         import :: sp
         complex(sp), intent(inout) :: x(*)
         complex(sp), intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine cscal
      pure subroutine zscal(n, a, x, incx)
         import :: dp
         complex(dp), intent(inout) :: x(*)
         complex(dp), intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine zscal
      pure subroutine csscal(n, a, x, incx)
         import :: sp
         complex(sp), intent(inout) :: x(*)
         real(sp), intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine csscal
      pure subroutine zdscal(n, a, x, incx)
         import :: dp
         complex(dp), intent(inout) :: x(*)
         real(dp), intent(in) :: a
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end subroutine zdscal
   end interface blas_scal

   !> Interchanges two vectors.
   !  Uses unrolled loops for increments equal to 1.
   interface blas_swap
      pure subroutine sswap(n, x, incx, y, incy)
         import :: sp
         real(sp), intent(inout) :: x(*)
         real(sp), intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine sswap
      pure subroutine dswap(n, x, incx, y, incy)
         import :: dp
         real(dp), intent(inout) :: x(*)
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine dswap
      pure subroutine cswap(n, x, incx, y, incy)
         import :: sp
         complex(sp), intent(inout) :: x(*)
         complex(sp), intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine cswap
      pure subroutine zswap(n, x, incx, y, incy)
         import :: dp
         complex(dp), intent(inout) :: x(*)
         complex(dp), intent(inout) :: y(*)
         integer, intent(in) :: incx
         integer, intent(in) :: incy
         integer, intent(in) :: n
      end subroutine zswap
   end interface blas_swap

   !> Finds the index of the first element having maximum absolute value.
   interface blas_iamax
      pure function isamax(n, x, incx)
         import :: sp
         integer :: isamax
         real(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function isamax
      pure function idamax(n, x, incx)
         import :: dp
         integer :: idamax
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function idamax
      pure function icamax(n, x, incx)
         import :: sp
         integer :: icamax
         complex(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function icamax
      pure function izamax(n, x, incx)
         import :: dp
         integer :: izamax
         complex(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         integer, intent(in) :: n
      end function izamax
   end interface blas_iamax

   !> Computes |Re(.)| + |Im(.)| of a complex number
   interface blas_cabs1
      pure function scabs1(c)
         import :: sp
         real(sp) :: scabs1
         complex(sp), intent(in) :: c
      end function scabs1
      pure function dcabs1(z)
         import :: dp
         real(dp) :: dcabs1
         complex(dp), intent(in) :: z
      end function dcabs1
   end interface blas_cabs1


contains


function mctc_sasum(xvec) result(asum)
   real(sp) :: asum
   real(sp), intent(in) :: xvec(:)
   integer :: incx, n
   incx = 1
   n = size(xvec)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec) create(asum)
   !$acc host_data use_device(xvec, asum)
   asum = cublasSasum(n, xvec, incx)
   !$acc end host_data
   !$acc exit data copyout(asum) delete(xvec)
#else
   asum = blas_asum(n, xvec, incx)
#endif
end function mctc_sasum


function mctc_dasum(xvec) result(asum)
   real(dp) :: asum
   real(dp), intent(in) :: xvec(:)
   integer :: incx, n
   incx = 1
   n = size(xvec)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec) create(asum)
   !$acc host_data use_device(xvec, asum)
   asum = cublasDasum(n, xvec, incx)
   !$acc end host_data
   !$acc exit data copyout(asum) delete(xvec)
#else
   asum = blas_asum(n, xvec, incx)
#endif
end function mctc_dasum


subroutine mctc_saxpy(xvec, yvec, alpha)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp) :: a
   integer :: incx, incy, n
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   incx = 1
   incy = 1
   n = size(xvec)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec)
   !$acc host_data use_device(xvec, yvec)
   call cublasSaxpy(n, a, xvec, incx, yvec, incy)
   !$acc end host_data
   !$acc exit data copyout(yvec) delete(xvec)
#else
   call blas_axpy(n, a, xvec, incx, yvec, incy)
#endif
end subroutine mctc_saxpy


subroutine mctc_daxpy(xvec, yvec, alpha)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp) :: a
   integer :: incx, incy, n
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   incx = 1
   incy = 1
   n = size(xvec)
#ifdef USE_CUBLAS
   !$acc enter data copyin(xvec, yvec)
   !$acc host_data use_device(xvec, yvec)
   call cublasDaxpy(n, a, xvec, incx, yvec, incy)
   !$acc end host_data
   !$acc exit data copyout(yvec) delete(xvec)
#else
   call blas_axpy(n, a, xvec, incx, yvec, incy)
#endif
end subroutine mctc_daxpy


subroutine mctc_scopy(xvec, yvec)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   integer :: incx, incy, n
   incx = 1
   incy = 1
   n = size(xvec)
   call blas_copy(n, xvec, incx, yvec, incy)
end subroutine mctc_scopy


subroutine mctc_dcopy(xvec, yvec)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   integer :: incx, incy, n
   incx = 1
   incy = 1
   n = size(xvec)
   call blas_copy(n, xvec, incx, yvec, incy)
end subroutine mctc_dcopy


function mctc_sdot(xvec, yvec) result(dot)
   real(sp) :: dot
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(in) :: yvec(:)
   integer :: incx, incy, n
   incx = 1
   incy = 1
   n = size(xvec)
   dot = blas_dot(n, xvec, incx, yvec, incy)
end function mctc_sdot


function mctc_ddot(xvec, yvec) result(dot)
   real(dp) :: dot
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(in) :: yvec(:)
   integer :: incx, incy, n
   incx = 1
   incy = 1
   n = size(xvec)
   dot = blas_dot(n, xvec, incx, yvec, incy)
end function mctc_ddot


function mctc_snrm2(xvec) result(nrm2)
   real(sp) :: nrm2
   real(sp), intent(in) :: xvec(:)
   integer :: incx, n
   incx = 1
   n = size(xvec)
   nrm2 = blas_nrm2(n, xvec, incx)
end function mctc_snrm2


function mctc_dnrm2(xvec) result(nrm2)
   real(dp) :: nrm2
   real(dp), intent(in) :: xvec(:)
   integer :: incx, n
   incx = 1
   n = size(xvec)
   nrm2 = blas_nrm2(n, xvec, incx)
end function mctc_dnrm2


subroutine mctc_srot(xvec, yvec, c, s)
   real(sp), intent(in) :: c
   real(sp), intent(in) :: s
   real(sp), intent(inout) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   integer :: incx, incy, n
   incx = 1
   incy = 1
   n = size(xvec)
   call blas_rot(n, xvec, incx, yvec, incy, c, s)
end subroutine mctc_srot


subroutine mctc_drot(xvec, yvec, c, s)
   real(dp), intent(in) :: c
   real(dp), intent(in) :: s
   real(dp), intent(inout) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   integer :: incx, incy, n
   incx = 1
   incy = 1
   n = size(xvec)
   call blas_rot(n, xvec, incx, yvec, incy, c, s)
end subroutine mctc_drot


subroutine mctc_sscal(xvec, alpha)
   real(sp), intent(in) :: alpha
   real(sp), intent(inout) :: xvec(:)
   integer :: incx, n
   incx = 1
   n = size(xvec)
   call blas_scal(n, alpha, xvec, incx)
end subroutine mctc_sscal


subroutine mctc_dscal(xvec, alpha)
   real(dp), intent(in) :: alpha
   real(dp), intent(inout) :: xvec(:)
   integer :: incx, n
   incx = 1
   n = size(xvec)
   call blas_scal(n, alpha, xvec, incx)
end subroutine mctc_dscal


subroutine mctc_sswap(xvec, yvec)
   real(sp), intent(inout) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   integer :: incx, incy, n
   incx = 1
   incy = 1
   n = size(xvec)
   call blas_swap(n, xvec, incx, yvec, incy)
end subroutine mctc_sswap


subroutine mctc_dswap(xvec, yvec)
   real(dp), intent(inout) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   integer :: incx, incy, n
   incx = 1
   incy = 1
   n = size(xvec)
   call blas_swap(n, xvec, incx, yvec, incy)
end subroutine mctc_dswap


function mctc_isamax(xvec) result(iamax)
   integer :: iamax
   real(sp), intent(in) :: xvec(:)
   integer :: incx, n
   incx = 1
   n = size(xvec)
   iamax = blas_iamax(n, xvec, incx)
end function mctc_isamax


function mctc_idamax(xvec) result(iamax)
   integer :: iamax
   real(dp), intent(in) :: xvec(:)
   integer :: incx, n
   incx = 1
   n = size(xvec)
   iamax = blas_iamax(n, xvec, incx)
end function mctc_idamax


end module xtb_mctc_blas_level1
