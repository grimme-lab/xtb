! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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
   use mctc_global
   use mctc_timings
   use mctc_filetools
   implicit none
   character(len=*),intent(in) :: progname
   integer,         intent(in) :: ntimer
   logical,         intent(in) :: verbose

!! ========================================================================
!  signal processing
   external :: wSIGTERM
   external :: wSIGINT
!  here two important signal handlers are installed, it seems that
!  FORTRAN by itself cannot handle signals in the way I expected it
!  to do, but this will force it to die when I hit CTRL^C.
!  UPDATE: Since the FORTRAN standard does not specify which signals
!          are handled by the language, some compilers may handle no
!          signals at all. This is *non-standard* so it will break
!          your program on some platforms with some compilers!
   call signal(2,wSIGINT)
   call signal(15,wSIGTERM)

!  initialize the timing system
   call start_timing_run
   call init_timing(ntimer,verb=verbose) ! verbosity allows printing of cputime
   call start_timing(1)

!  initialize the messagebuffer for the error handler
   call init_errorbuffer

!  initialize the filelist
   call init_filelist(20)

!  set this for mctc_global
   name = progname

end subroutine mctc_init

subroutine mctc_sanity(sane)
   use mctc_global
   logical,intent(out) :: sane
   sane = good
end subroutine mctc_sanity

subroutine mctc_strict
   use mctc_global
   strict = .true.
end subroutine mctc_strict

subroutine mctc_mute
   use mctc_global
   mute = .true.
end subroutine mctc_mute
