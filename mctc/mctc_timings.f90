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

module mctc_timings
   use iso_fortran_env, wp => real64
   implicit none

   public  :: init_timing,start_timing_run,stop_timing_run
   public  :: start_timing,stop_timing
   public  :: prdate,prtiming

   private :: timing_wall,timing_cpu,timing_max
   private :: start_date,start_time,start_zone,start_values
   private :: stop_date,stop_time,stop_zone,stop_values
   private :: timing
   private :: verbose

   real(wp),allocatable :: timing_wall(:)
   real(wp),allocatable :: timing_cpu(:)
   integer :: timing_max
   logical  :: verbose = .false.
   
   character(len=8)  :: start_date,stop_date
   character(len=10) :: start_time,stop_time
   character(len=5)  :: start_zone,stop_zone
   integer(int64) :: start_values(8),stop_values(8)

   intrinsic :: kind,date_and_time,system_clock,cpu_time,present

contains

subroutine prtiming(i,inmsg)
   implicit none
   integer,intent(in) :: i
   character(len=*),intent(in),optional :: inmsg
   character(len=:),allocatable :: msg
   real(wp) :: cputime,walltime
   integer(int64) ::  cpudays, cpuhours, cpumins
   integer(int64) :: walldays,wallhours,wallmins
!                      <7-------> <5---------> <7-----------> <10------------>
!  '(1x,a,1x,"time:",1x,5i,1x,"d",1x,2i,1x,"h",1x,2i,1x,"min",1x,f5.2,1x,"sec")'
!                     '(5i,1x,"d")'
!                             '(a,1x,2i,1x,"h")'
!                                            '(1x,2i,1x,"min")'                 
!                                                           '(1x,f5.2,1x,"sec")'
!  '(1x,a,1x,"time:",1x,a)'
   if (present(inmsg)) then
      msg = inmsg
   else
      msg = '*'
   endif
   !           DAYS   HOURS   MINUTES   SECONDS
   ! DAYS        1     1/24    1/1440   1/86400
   ! HOURS      24      1       1/60     1/3600
   ! MINUTES   1440    60        1        1/60
   ! SECONDS  86400   3600      60         1
   cputime = timing_cpu (i)
!  if (cputime.ge.86400._wp) then
      cpudays = int(cputime/86400._wp)
      cputime = cputime - cpudays*86400._wp
!  endif
!  if (cputime.ge.3600._wp) then
      cpuhours = int(cputime/3600._wp)
      cputime = cputime - cpuhours*3600._wp
!  endif
!  if (cputime.ge.60._wp) then
      cpumins = int(cputime/60._wp)
      cputime = cputime - cpumins*60._wp
!  endif
   walltime = timing_wall(i)
!  if (walltime.ge.86400._wp) then
      walldays = int(walltime/86400._wp)
      walltime = walltime - walldays*86400._wp
!  endif
!  if (walltime.ge.3600._wp) then
      wallhours = int(walltime/3600._wp)
      walltime = walltime - wallhours*3600._wp
!  endif
!  if (walltime.ge.60._wp) then
      wallmins = int(walltime/60._wp)
      walltime = walltime - wallmins*60._wp
!  endif
   
   if (verbose) then
      write(output_unit,'(1x,a,":")') msg
      write(output_unit,'(" * wall-time: ",i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') walldays,wallhours,wallmins,walltime
      write(output_unit,'(" *  cpu-time: ",i5," d, ",i2," h, ",i2," min, ",f6.3," sec")')  cpudays, cpuhours, cpumins, cputime
      write(output_unit,'(1x,"*",1x,"ratio c/w:",1x,f9.3,1x,"speedup")') timing_cpu (i)/timing_wall(i)
   else
      write(output_unit,'(1x,a,1x,"time: ",i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') msg,walldays,wallhours,wallmins,walltime
   endif
end subroutine prtiming

function prtimestring(mode) result(timestring)
   character,intent(in) :: mode
   character(len=*),parameter :: timefmt = &
   &  '(a,''/'',a,''/'',a,1x,''at'',1x,a,'':'',a,'':'',a)'
   character(len=8)  :: date
   character(len=10) :: time
   character(len=5)  :: zone
   character(len=31) :: timestring
   select case(mode)
   case('S','s')
   write(timestring,timefmt) &
   &   start_date(:4),start_date(5:6),start_date(7:), &
   &   start_time(:2),start_time(3:4),start_time(5:)
   case('E','e')
   write(timestring,timefmt) &
   &   stop_date(:4),stop_date(5:6),stop_date(7:), &
   &   stop_time(:2),stop_time(3:4),stop_time(5:)
   case default
   call date_and_time(date,time,zone)
   write(timestring,timefmt) &
   &   date(:4),date(5:6),date(7:), &
   &   time(:2),time(3:4),time(5:)
   end select
end function prtimestring

subroutine prdate(mode)
   implicit none
   character,intent(in) :: mode
   character(len=*),parameter :: outfmt = '('' * '',a,1x,a)'
   character(len=8)  :: date
   character(len=10) :: time
   character(len=5)  :: zone
   character(len=31) :: timestring
   integer(int64) :: values(8)
   select case(mode)
   case('S','s')
   write(output_unit,outfmt) 'started run on', prtimestring(mode)
   case('E','e')
   write(output_unit,outfmt) 'finished run on', prtimestring(mode)
   case default
   write(output_unit,outfmt) 'current time:', prtimestring(mode)
   end select
end subroutine prdate

subroutine init_timing(i,verb)
   implicit none
   integer,intent(in) :: i
   logical,intent(in),optional :: verb
   if (allocated(timing_wall)) deallocate(timing_wall)
   if (allocated(timing_cpu))  deallocate(timing_cpu)
   if (present(verb)) verbose = verb
   timing_max = i
   allocate( timing_wall(i),  &
   &         timing_cpu (i),  &
   &         source=0.0_wp )
end subroutine init_timing

subroutine start_timing_run
   call date_and_time(start_date,start_time,start_zone,start_values)
end subroutine start_timing_run

subroutine stop_timing_run
   call date_and_time(stop_date,stop_time,stop_zone,stop_values)
end subroutine stop_timing_run

subroutine start_timing(i)
   implicit none
   integer,intent(in) :: i
   real(wp) :: time_cpu
   real(wp) :: time_wall
   call timing(time_cpu,time_wall)
   timing_cpu (i) = timing_cpu (i) - time_cpu
   timing_wall(i) = timing_wall(i) - time_wall
end subroutine start_timing

subroutine stop_timing(i)
   implicit none
   integer,intent(in) :: i
   real(wp) :: time_cpu
   real(wp) :: time_wall
   call timing(time_cpu,time_wall)
   timing_cpu (i) = timing_cpu (i) + time_cpu
   timing_wall(i) = timing_wall(i) + time_wall
end subroutine stop_timing

subroutine timing(time_cpu,time_wall)
   implicit none
   real(wp),intent(out) :: time_cpu
   real(wp),intent(out) :: time_wall
   integer(int64) :: time_count,time_rate,time_max
   call system_clock(time_count,time_rate,time_max)
   call cpu_time(time_cpu)
   time_wall = real(time_count,wp)/real(time_rate,wp)
end subroutine timing

end module mctc_timings
