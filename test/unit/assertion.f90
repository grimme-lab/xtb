module assertion
   use, intrinsic :: iso_fortran_env, only: int16, int32, int64, real32, real64
   use xtb_mctc_io, only : stderr
   implicit none
   private :: stderr

   interface assert_eq
      module procedure :: assert_eq_char
      module procedure :: assert_eq_int16
      module procedure :: assert_eq_int32
      module procedure :: assert_eq_int64
      module procedure :: assert_eq_int16_array
      module procedure :: assert_eq_int32_array
      module procedure :: assert_eq_int64_array
   end interface assert_eq

   interface assert_close
      module procedure :: assert_close_int16
      module procedure :: assert_close_int32
      module procedure :: assert_close_int64
      module procedure :: assert_close_real32
      module procedure :: assert_close_real64
   end interface assert_close

   integer, public :: afail = 0

contains

subroutine assert(bool)
   logical,intent(in) :: bool

   if (.not.bool) then
      write(stderr,'("assertion FAILED")')
      afail = afail+1
   endif
end subroutine assert

subroutine assert_eq_char(val1,val2)
   character(len=*),intent(in) :: val1,val2

   if (val1 /= val2) then
      write(stderr,'("assertion:",1x,a," == ",a,1x,"FAILED")') &
         val1,val2
      afail = afail+1
   endif
end subroutine assert_eq_char

subroutine assert_eq_int16(val1,val2)
   integer(int16),intent(in) :: val1,val2

   if (val1 /= val2) then
      write(stderr,'("assertion:",1x,g21.14," == ",g21.14,1x,"FAILED")') &
         val1,val2
      afail = afail+1
   endif
end subroutine assert_eq_int16

subroutine assert_eq_int32(val1,val2)
   integer(int32),intent(in) :: val1,val2

   if (val1 /= val2) then
      write(stderr,'("assertion:",1x,g21.14," == ",g21.14,1x,"FAILED")') &
         val1,val2
      afail = afail+1
   endif
end subroutine assert_eq_int32

subroutine assert_eq_int64(val1,val2)
   integer(int64),intent(in) :: val1,val2

   if (val1 /= val2) then
      write(stderr,'("assertion:",1x,g21.14," == ",g21.14,1x,"FAILED")') &
         val1,val2
      afail = afail+1
   endif
end subroutine assert_eq_int64

subroutine assert_eq_int16_array(val1,val2)
   integer(int16),intent(in) :: val1(:),val2(:)
   integer :: i

   if (size(val1) .ne. size(val2)) then
      write(stderr,'("shape missmatch:",1x,i0," == ",i0,1x,"FAILED")') &
         val1,val2
      afail = afail+1
   endif
   if (any(val1 /= val2)) then
      do i = 1, size(val1)
         if (val1(i) /= val2(i)) &
         write(stderr,'("assertion:",1x,g21.14," == ",g21.14,1x,"FAILED at ",i0)')&
            val1(i),val2(i),i
      enddo
      afail = afail+1
   endif
end subroutine assert_eq_int16_array

subroutine assert_eq_int32_array(val1,val2)
   integer(int32),intent(in) :: val1(:),val2(:)
   integer :: i

   if (size(val1) .ne. size(val2)) then
      write(stderr,'("shape missmatch:",1x,i0," == ",i0,1x,"FAILED")') &
         val1,val2
      afail = afail+1
   endif
   if (any(val1 /= val2)) then
      do i = 1, size(val1)
         if (val1(i) /= val2(i)) &
         write(stderr,'("assertion:",1x,g21.14," == ",g21.14,1x,"FAILED at ",i0)')&
            val1(i),val2(i),i
      enddo
      afail = afail+1
   endif
end subroutine assert_eq_int32_array

subroutine assert_eq_int64_array(val1,val2)
   integer(int64),intent(in) :: val1(:),val2(:)
   integer :: i

   if (size(val1) .ne. size(val2)) then
      write(stderr,'("shape missmatch:",1x,i0," == ",i0,1x,"FAILED")') &
         val1,val2
      afail = afail+1
   endif
   if (any(val1 /= val2)) then
      do i = 1, size(val1)
         if (val1(i) /= val2(i)) &
         write(stderr,'("assertion:",1x,g21.14," == ",g21.14,1x,"FAILED at ",i0)')&
            val1(i),val2(i),i
      enddo
      afail = afail+1
   endif
end subroutine assert_eq_int64_array

subroutine assert_close_real64(val1,val2,thr)
   real(real64),intent(in) :: val1,val2,thr
   real(real64) :: diff

   diff = val1 - val2
   if (abs(diff) > thr) then
      write(stderr,'("assertion:",1x,g21.14," == ",g21.14,1x,"FAILED")') &
         val1,val2
      afail = afail+1
   endif
end subroutine assert_close_real64

subroutine assert_close_real32(val1,val2,thr)
   real(real32),intent(in) :: val1,val2,thr
   real(real32) :: diff

   diff = val1 - val2
   if (abs(diff) > thr) then
      write(stderr,'("assertion:",1x,g21.14," == ",g21.14,1x,"FAILED")') &
         val1,val2
      afail = afail+1
   endif
end subroutine assert_close_real32

subroutine assert_close_int16(val1,val2,thr)
   integer(int16),intent(in) :: val1,val2,thr
   integer(int16) :: diff

   diff = val1 - val2
   if (abs(diff) > thr) then
      write(stderr,'("assertion:",1x,g21.14," == ",g21.14,1x,"FAILED")') &
         val1,val2
      afail = afail+1
   endif
end subroutine assert_close_int16

subroutine assert_close_int32(val1,val2,thr)
   integer(int32),intent(in) :: val1,val2,thr
   integer(int32) :: diff

   diff = val1 - val2
   if (abs(diff) > thr) then
      write(stderr,'("assertion:",1x,g21.14," == ",g21.14,1x,"FAILED")') &
         val1,val2
      afail = afail+1
   endif
end subroutine assert_close_int32

subroutine assert_close_int64(val1,val2,thr)
   integer(int64),intent(in) :: val1,val2,thr
   integer(int64) :: diff

   diff = val1 - val2
   if (abs(diff) > thr) then
      write(stderr,'("assertion:",1x,g21.14," == ",g21.14,1x,"FAILED")') &
         val1,val2
      afail = afail+1
   endif
end subroutine assert_close_int64

end module assertion
