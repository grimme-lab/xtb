! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

subroutine mctc_init(progname,ntimer,verbose)
   use, intrinsic :: iso_c_binding, only : c_int, c_funptr, c_funloc
   use xtb_mctc_global
   use xtb_mctc_timings
   use xtb_mctc_filetools
   use xtb_type_environment, only : init
   implicit none
   character(len=*),intent(in) :: progname
   integer,         intent(in) :: ntimer
   logical,         intent(in) :: verbose
#ifdef USE_CUSOLVER
   integer :: err
#endif

!! ========================================================================
!  signal processing
   interface
      subroutine xtb_wsigint(signal) bind(C)
         import :: c_int
         integer(c_int), value :: signal
      end subroutine xtb_wsigint
      subroutine xtb_wsigterm(signal) bind(C)
         import :: c_int
         integer(c_int), value :: signal
      end subroutine xtb_wsigterm
      subroutine xtb_signal_handler(signal, handler) bind(C)
         import :: c_int, c_funptr
         integer(c_int), value :: signal
         type(c_funptr), value :: handler
      end subroutine xtb_signal_handler
   end interface
!  here two important signal handlers are installed, it seems that
!  FORTRAN by itself cannot handle signals in the way I expected it
!  to do, but this will force it to die when I hit CTRL^C.
   call xtb_signal_handler(2_c_int, c_funloc(xtb_wsigint))
   call xtb_signal_handler(15_c_int, c_funloc(xtb_wsigterm))

!  initialize the timing system
   call start_timing_run
   call init_timing(ntimer,verb=verbose) ! verbosity allows printing of cputime
   call start_timing(1)

!  initialize the messagebuffer for the error handler
   call init_errorbuffer

!  initialize the filelist
   call init_filelist(20)

!  set this for xtb_mctc_global
   name = progname

   allocate(persistentEnv)
   call init(persistentEnv)

#ifdef USE_CUSOLVER
   err = cusolverDnCreate(cusolverDnH)
   if (err /= 0) then
      call persistentEnv%error("failed to create cusolver handle", "mctc_init")
   end if
#endif
end subroutine mctc_init

subroutine mctc_sanity(sane)
   use xtb_mctc_global
   logical,intent(out) :: sane
   sane = good
end subroutine mctc_sanity

subroutine mctc_strict
   use xtb_mctc_global
   strict = .true.
end subroutine mctc_strict

subroutine mctc_mute
   use xtb_mctc_global
   mute = .true.
end subroutine mctc_mute
