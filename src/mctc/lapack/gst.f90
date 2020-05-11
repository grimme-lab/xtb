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

!> Helper routines to transform generalized eigenvalue problems in standard form.
module xtb_mctc_lapack_gst
   use xtb_mctc_accuracy, only : sp, dp
   use xtb_type_environment, only : TEnvironment
   implicit none
   private

   public :: mctc_sygst, mctc_spgst

   public :: lapack_sygst, lapack_spgst, lapack_hegst, lapack_hpgst


   !> Reduces a real symmetric-definite generalized eigenproblem
   !  to standard form.
   !
   !  If ITYPE = 1, the problem is A*x = lambda*B*x,
   !  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
   !
   !  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
   !  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
   !
   !  B must have been previously factorized as U**T*U or L*L**T by DPOTRF.
   interface mctc_sygst
      module procedure :: mctc_ssygst
      module procedure :: mctc_dsygst
   end interface mctc_sygst

   !> reduces a real symmetric-definite generalized eigenproblem
   !  to standard form, using packed storage.
   !
   !  If ITYPE = 1, the problem is A*x = lambda*B*x,
   !  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
   !
   !  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
   !  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
   !
   !  B must have been previously factorized as U**T*U or L*L**T by DPPTRF.
   interface mctc_spgst
      module procedure :: mctc_sspgst
      module procedure :: mctc_dspgst
   end interface mctc_spgst


   !> DSYGST reduces a real symmetric-definite generalized eigenproblem
   !  to standard form.
   !
   !  If ITYPE = 1, the problem is A*x = lambda*B*x,
   !  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
   !
   !  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
   !  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
   !
   !  B must have been previously factorized as U**T*U or L*L**T by DPOTRF.
   interface lapack_sygst
      pure subroutine ssygst(itype, uplo, n, a, lda, b, ldb, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(in) :: b(ldb, *)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine ssygst
      pure subroutine dsygst(itype, uplo, n, a, lda, b, ldb, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(in) :: b(ldb, *)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine dsygst
   end interface lapack_sygst

   !> reduces a complex Hermitian-definite generalized
   !  eigenproblem to standard form.
   !
   !  If ITYPE = 1, the problem is A*x = lambda*B*x,
   !  and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
   !
   !  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
   !  B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
   !
   !  B must have been previously factorized as U**H*U or L*L**H by CPOTRF.
   interface lapack_hegst
      pure subroutine chegst(itype, uplo, n, a, lda, b, ldb, info)
         import :: sp
         complex(sp), intent(inout) :: a(lda, *)
         complex(sp), intent(in) :: b(ldb, *)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine chegst
      pure subroutine zhegst(itype, uplo, n, a, lda, b, ldb, info)
         import :: dp
         complex(dp), intent(inout) :: a(lda, *)
         complex(dp), intent(in) :: b(ldb, *)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine zhegst
   end interface lapack_hegst

   !> reduces a real symmetric-definite generalized eigenproblem
   !  to standard form, using packed storage.
   !
   !  If ITYPE = 1, the problem is A*x = lambda*B*x,
   !  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
   !
   !  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
   !  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
   !
   !  B must have been previously factorized as U**T*U or L*L**T by DPPTRF.
   interface lapack_spgst
      pure subroutine sspgst(itype, uplo, n, ap, bp, info)
         import :: sp
         real(sp), intent(inout) :: ap(*)
         real(sp), intent(in) :: bp(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine sspgst
      pure subroutine dspgst(itype, uplo, n, ap, bp, info)
         import :: dp
         real(dp), intent(inout) :: ap(*)
         real(dp), intent(in) :: bp(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine dspgst
   end interface lapack_spgst

   !> reduces a complex Hermitian-definite generalized
   !  eigenproblem to standard form, using packed storage.
   !
   !  If ITYPE = 1, the problem is A*x = lambda*B*x,
   !  and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
   !
   !  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
   !  B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
   !
   !  B must have been previously factorized as U**H*U or L*L**H by CPPTRF.
   interface lapack_hpgst
      pure subroutine chpgst(itype, uplo, n, ap, bp, info)
         import :: sp
         complex(sp), intent(inout) :: ap(*)
         complex(sp), intent(in) :: bp(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine chpgst
      pure subroutine zhpgst(itype, uplo, n, ap, bp, info)
         import :: dp
         complex(dp), intent(inout) :: ap(*)
         complex(dp), intent(in) :: bp(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
      end subroutine zhpgst
   end interface lapack_hpgst


contains


subroutine mctc_ssygst(env, amat, bmat, itype, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sygst'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(in) :: bmat(:, :)
   integer, intent(in), optional :: itype
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: ity, info, n, lda, ldb
   if(present(itype)) then
      ity = itype
   else
      ity = 1
   endif
   if(present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   endif
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   call lapack_sygst(ity, ula, n, amat, lda, bmat, ldb, info)
   if (info /= 0) then
      call env%error("Transformation of the eigenvalue problem failed", source)
   end if
end subroutine mctc_ssygst


subroutine mctc_dsygst(env, amat, bmat, itype, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_sygst'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(in) :: bmat(:, :)
   integer, intent(in), optional :: itype
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: ity, info, n, lda, ldb
   if(present(itype)) then
      ity = itype
   else
      ity = 1
   endif
   if(present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   endif
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   call lapack_sygst(ity, ula, n, amat, lda, bmat, ldb, info)
   if (info /= 0) then
      call env%error("Transformation of the eigenvalue problem failed", source)
   end if
end subroutine mctc_dsygst


subroutine mctc_sspgst(env, amat, bmat, itype, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_spgst'
   type(TEnvironment), intent(inout) :: env
   real(sp), intent(inout) :: amat(:)
   real(sp), intent(in) :: bmat(:)
   integer, intent(in), optional :: itype
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: ity, info, n, nn
   if(present(itype)) then
      ity = itype
   else
      ity = 1
   endif
   if(present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   endif
   nn = size(amat)
   if(nn <= 0) then
      n = nn
   else
      n = int((-1+sqrt(1+8*real(nn, sp))))/2
   endif
   call lapack_spgst(ity, ula, n, amat, bmat, info)
   if (info /= 0) then
      call env%error("Transformation of the eigenvalue problem failed", source)
   end if
end subroutine mctc_sspgst


subroutine mctc_dspgst(env, amat, bmat, itype, uplo)
   character(len=*), parameter :: source = 'mctc_lapack_spgst'
   type(TEnvironment), intent(inout) :: env
   real(dp), intent(inout) :: amat(:)
   real(dp), intent(in) :: bmat(:)
   integer, intent(in), optional :: itype
   character(len=1), intent(in), optional :: uplo
   character(len=1) :: ula
   integer :: ity, info, n, nn
   if(present(itype)) then
      ity = itype
   else
      ity = 1
   endif
   if(present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   endif
   nn = size(amat)
   if(nn <= 0) then
      n = nn
   else
      n = int((-1+sqrt(1+8*real(nn, dp))))/2
   endif
   call lapack_spgst(ity, ula, n, amat, bmat, info)
   if (info /= 0) then
      call env%error("Transformation of the eigenvalue problem failed", source)
   end if
end subroutine mctc_dspgst


end module xtb_mctc_lapack_gst
