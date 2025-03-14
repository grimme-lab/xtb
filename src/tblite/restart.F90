! This file is part of xtb.
! SPDX-Identifier: LGPL-3.0-or-later
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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

!> Restart file for calculations with tblite library
module xtb_tblite_restart
   use xtb_mctc_accuracy, only : wp, i8
   use xtb_type_environment, only : TEnvironment
   use xtb_type_restart, only : TRestart
   implicit none
   private

   public :: loadRestart, dumpRestart

contains

!> Load restart file
subroutine loadRestart(env, chk, fname, success)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Restart file
   type(TRestart), intent(inout) :: chk
   !> File name
   character(len=*), intent(in) :: fname
   !> Successfully loaded restart
   logical, intent(out) :: success

#if WITH_TBLITE
   external :: open_binary, close_file
   integer :: io, info
   integer(i8) :: version, dummy, nat, nsh, nel, nuhf

   ! call env%io%readBinary(io, fname)
   call open_binary(io, fname, "r")
   if (io == -1) then
      return
   end if

   read(io, iostat=info) &
      version, dummy, nat, nsh, nel, nuhf

   if (info == 0 &
      .and. nat == size(chk%tblite%n0at, kind=kind(nat)) &
      .and. nsh == size(chk%tblite%n0sh, kind=kind(nsh))) then
      read(io, iostat=info) chk%tblite%qsh
      if (version >= 2) then
         if (info == 0) &
            read(io, iostat=info) chk%tblite%dpat
         if (info == 0) &
            read(io, iostat=info) chk%tblite%qpat
      end if
   end if

   success = info == 0
   ! call env%io%closeFile(io)
   call close_file(io)
#else
   call feature_not_implemented(env)
#endif
end subroutine loadRestart

!> Dump restart file
subroutine dumpRestart(env, chk, fname)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Restart file
   type(TRestart), intent(inout) :: chk
   !> File name
   character(len=*), intent(in) :: fname

#if WITH_TBLITE
   external :: open_binary, close_file
   integer :: io

   ! call env%io%writeBinary(io, fname)
   call open_binary(io, fname, "w")
   if (io == -1) then
      return
   end if

   write(io) &
      2_i8, &
      0_i8, &
      size(chk%tblite%n0at, kind=i8), &
      size(chk%tblite%n0sh, kind=i8), &
      nint(chk%tblite%nocc, i8), &
      nint(chk%tblite%nuhf, i8)

   write(io) chk%tblite%qsh

   ! Always write multipole moments
   write(io) chk%tblite%dpat
   write(io) chk%tblite%qpat

   ! call env%io%closeFile(io)
   call close_file(io)
#else
   call feature_not_implemented(env)
#endif
end subroutine dumpRestart


#if ! WITH_TBLITE
subroutine feature_not_implemented(env)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env

   call env%error("Compiled without support for tblite library")
end subroutine feature_not_implemented
#endif

end module xtb_tblite_restart
