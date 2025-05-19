! This file is part of xtb.
!
! Copyright (C) 2025 Igor S. Gerasimov
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

!> Helpers to profile xtb with Tracy
module xtb_tracying
#ifdef WITH_TRACY
  use tracy
#endif
  use, intrinsic :: iso_c_binding, only: c_char, c_int32_t, c_int64_t
  implicit none
  private

  !>
  !> @brief Context manager for profiling zones in xtb.
  !>
  type :: xtb_zone_context
#ifdef WITH_TRACY
    type(tracy_zone_context) :: ctx
    logical                  :: inited = .false.
#endif
  contains
    procedure :: start => zone_start
    procedure :: end   => zone_end
#ifdef WITH_TRACY
    final     :: zone_final
#endif
  end type

  public :: xtb_zone_context
contains
  !> @brief Starts a new profiling zone
  !>
  !> @param[inout] zone           profiling zone context to start
  !> @param[in]    source         name of the source file
  !> @param[in]    function_name  name of the function where the zone is started
  !> @param[in]    line           line number where the zone is started
  !> @param[in]    zone_name      (Optional) name for the profiling zone
  !> @param[in]    color          (Optional) color for profiling zone
  !>
  subroutine zone_start(zone, source, function_name, line, zone_name, color)
    class(xtb_zone_context),                         intent(inout) :: zone
    character(kind=c_char, len=*),           target, intent(in)    :: source, function_name
    integer,                                         intent(in)    :: line
    character(kind=c_char, len=*), optional, target, intent(in)    :: zone_name
    integer(c_int32_t),            optional,         intent(in)    :: color
    ! local variables
    integer(c_int64_t) :: srcloc_id
    integer(c_int32_t) :: line_
#ifdef WITH_TRACY
    if (zone%inited) then
      call zone%end()
    end if
    line_ = int(line, kind=c_int32_t)
    srcloc_id = tracy_alloc_srcloc(line_, source, function_name, zone_name, color)
    zone%inited = .true.
    zone%ctx = tracy_zone_begin(srcloc_id)
#endif
  end subroutine zone_start
  !>
  !> @brief Ends profiling zone
  !>
  !> @param[inout] zone           profiling zone context to stop
  !>
  subroutine zone_end(zone)
    class(xtb_zone_context), intent(inout) :: zone
#ifdef WITH_TRACY
    if (zone%inited) then
      call tracy_zone_end(zone%ctx)
    end if
    zone%inited = .false.
#endif
  end subroutine zone_end
  !>
  !> @brief Ends profiling zone
  !>
  !> @details Using for automatic zone ends at the end of scope (as compiler defines it)
  !>
  !> @param[inout] zone           profiling zone context to stop
  !>
  subroutine zone_final(zone)
    type(xtb_zone_context), intent(inout) :: zone
    call zone%end()
  end subroutine zone_final
end module xtb_tracying
