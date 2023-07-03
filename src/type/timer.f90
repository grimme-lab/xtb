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

module xtb_type_timer
   use xtb_mctc_accuracy, only : wp, int64 => i8
   implicit none

   public :: tb_timer
   private

   type :: tb_timer
      integer, private :: n = 0
         !! number of timers
      logical, private :: verbose = .false.
         !! printlevel 
      real(wp),private :: totwall = 0.0_wp
      real(wp),private :: totcpu = 0.0_wp
      logical, private,allocatable :: running(:)
      real(wp),private,allocatable :: twall(:)
      real(wp),private,allocatable :: tcpu(:)
      character(len=40),private,allocatable :: tag(:)
   contains
   procedure :: new => allocate_timer
   procedure :: allocate => allocate_timer
   procedure :: deallocate => deallocate_timer
   procedure :: measure => timer
   procedure :: write_timing
   procedure :: write => write_all_timings
   procedure :: get => get_timer
   procedure,private :: start_timing
   procedure,private :: stop_timing
   end type tb_timer

contains

!> To initialize timer
subroutine allocate_timer(self,n,verbose)
   
   implicit none
   class(tb_timer),intent(inout) :: self
      !! instance of timer
   integer, intent(in)           :: n
      !! number of timers
   logical, intent(in), optional :: verbose
      !! if verbose

   real(wp) :: time_cpu
   real(wp) :: time_wall

   call self%deallocate
   if (n < 1) return
      !! capture negative values
   
   self%n = n
   if (present(verbose)) self%verbose = verbose
   allocate( self%twall(0:n), source = 0.0_wp )
   allocate( self%tcpu(0:n),  source = 0.0_wp )
   allocate( self%running(n), source =.false. )
   allocate( self%tag(n) ); self%tag = ' '

   call self%start_timing(0)
      !! launch timer

end subroutine allocate_timer

!--------------------------------------------
!> To deallocate memory
!--------------------------------------------
subroutine deallocate_timer(self)
   
   implicit none
   !> dummy arg list
   class(tb_timer),intent(inout) :: self
      !! instance of timer
   
   self%n = 0
   self%totwall = 0
   self%totcpu  = 0
   self%verbose = .false.
   if (allocated(self%twall))   deallocate(self%twall)
   if (allocated(self%tcpu))    deallocate(self%tcpu)
   if (allocated(self%running)) deallocate(self%running)

end subroutine deallocate_timer

!--------------------------------------------
!> To obtain current elapsed time
!--------------------------------------------
function get_timer(self,i) result(time)

   !> dummy arg list
   class(tb_timer),intent(inout) :: self
      !! instance of timer
   integer,intent(in),optional :: i
      !! specific timer

   !> local vars
   integer  :: it
   real(wp) :: tcpu,twall
   real(wp) :: time
   logical  :: running
   
   !> if i is not given, calculate overall elapsed time
   if (present(i)) then
      it = i
   else
      it = 0
   endif
   
   if (it > 0) then
      running = self%running(it)
   else
      running = .true.
   endif
   
   if (running) then
      call timing(tcpu,twall)
      time = self%twall(it) + twall
   else
      time = self%twall(it)
   endif

end function get_timer

!--------------------------------------------
!> To write timing for specific timer 
!--------------------------------------------
subroutine write_timing(self,iunit,i,inmsg,verbose)

   implicit none
   !> dummy arg list
   class(tb_timer),intent(inout) :: self
      !! instance of timer
   integer,intent(in) :: iunit
      !! I/O unit
   integer,intent(in) :: i
      !! index 
   character(len=*),intent(in),optional :: inmsg
      !! raw message text 
   logical,intent(in),optional :: verbose
      !! if verbose
   
   !> local vars
   character(len=26) :: msg
   real(wp) :: cputime,walltime
   integer(int64) ::  cpudays, cpuhours, cpumins
   integer(int64) :: walldays,wallhours,wallmins
   logical :: lverbose

!  '(1x,a,1x,"time:",1x,a)'
   !> check if a tag should be added
   if (present(inmsg)) then
      msg = inmsg
   else
      msg = self%tag(i)
   endif
   
   !> verbosity settings
   if (present(verbose)) then
      lverbose = verbose
   else
      lverbose = self%verbose
   endif
   
   !           DAYS   HOURS   MINUTES   SECONDS
   ! DAYS        1     1/24    1/1440   1/86400
   ! HOURS      24      1       1/60     1/3600
   ! MINUTES   1440    60        1        1/60
   ! SECONDS  86400   3600      60         1
   
   !> convert elapsed CPU time into days, hours, minutes
   cputime = self%tcpu (i)
   cpudays = int(cputime/86400._wp)
   cputime = cputime - cpudays*86400._wp
   cpuhours = int(cputime/3600._wp)
   cputime = cputime - cpuhours*3600._wp
   cpumins = int(cputime/60._wp)
   cputime = cputime - cpumins*60._wp

   !> convert elapsed wall time into days, hours, minutes
   walltime = self%twall(i)
   walldays = int(walltime/86400._wp)
   walltime = walltime - walldays*86400._wp
   wallhours = int(walltime/3600._wp)
   walltime = walltime - wallhours*3600._wp
   wallmins = int(walltime/60._wp)
   walltime = walltime - wallmins*60._wp
   
   !> printout
   if (lverbose) then
      write(iunit,'(1x,a)') msg
      write(iunit,'(" * wall-time: ",i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
         walldays,wallhours,wallmins,walltime
      write(iunit,'(" *  cpu-time: ",i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
         cpudays, cpuhours, cpumins, cputime
      write(iunit,'(1x,"*",1x,"ratio c/w:",1x,f9.3,1x,"speedup")') self%tcpu (i)/self%twall(i)
   else
      write(iunit,'(1x,a30,1x,"...",i9," min, ",f6.3," sec")') &
         msg, wallmins, walltime
   endif

end subroutine write_timing

!--------------------------------------------
!> To write timing for all timers
!--------------------------------------------
subroutine write_all_timings(self,iunit,inmsg)
   
   implicit none
   !> dummy args list
   class(tb_timer),intent(inout) :: self
      !! instance of timer
   integer,intent(in) :: iunit
      !! I/O unit
   character(len=*),intent(in),optional :: inmsg
      !! message 

   !> local vars
   character(len=26) :: msg
   real(wp) :: cputime,walltime
   integer  :: i
   integer(int64) ::  cpudays, cpuhours, cpumins
   integer(int64) :: walldays,wallhours,wallmins

   call self%stop_timing(0)

!  '(1x,a,1x,"time:",1x,a)'
   !> check if an external message should be added
   if (present(inmsg)) then
      msg = inmsg // " (total)"
   else
      msg = "total time"
   endif
   !           DAYS   HOURS   MINUTES   SECONDS
   ! DAYS        1     1/24    1/1440   1/86400
   ! HOURS      24      1       1/60     1/3600
   ! MINUTES   1440    60        1        1/60
   ! SECONDS  86400   3600      60         1
   
   !> convert overall elapsed CPU time into days, hours, minutes
   cputime = self%tcpu (0)
   cpudays = int(cputime/86400._wp)
   cputime = cputime - cpudays*86400._wp
   cpuhours = int(cputime/3600._wp)
   cputime = cputime - cpuhours*3600._wp
   cpumins = int(cputime/60._wp)
   cputime = cputime - cpumins*60._wp

   !> convert overall elapsed wall time into days, hours, minutes
   walltime = self%twall(0)
   walldays = int(walltime/86400._wp)
   walltime = walltime - walldays*86400._wp
   wallhours = int(walltime/3600._wp)
   walltime = walltime - wallhours*3600._wp
   wallmins = int(walltime/60._wp)
   walltime = walltime - wallmins*60._wp
   
   !> overall printout
   write(iunit,'(a)')
   if (self%verbose) then
      write(iunit,'(1x,a,":")') msg
      write(iunit,'(" * wall-time: ",i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
         walldays,wallhours,wallmins,walltime
      write(iunit,'(" *  cpu-time: ",i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
         cpudays, cpuhours, cpumins, cputime
      write(iunit,'(1x,"*",1x,"ratio c/w:",1x,f9.3,1x,"speedup")') self%tcpu (0)/self%twall(0)
   else
      write(iunit,'(1x,a26,i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
         msg,walldays,wallhours,wallmins,walltime
   endif

   !> printout every timer and corresponding speedup
   do i = 1, self%n
      walltime = self%twall(i)
      wallmins = int(walltime/60._wp)
      walltime = walltime - wallmins*60._wp
      write(iunit,'(1x,a30,1x,"...",i9," min, ",f6.3," sec (",f7.3,"%)")') &
         self%tag(i), wallmins, walltime, 100*self%twall(i)/self%twall(0)
   enddo
   write(iunit,'(a)')
   
end subroutine write_all_timings

!--------------------------------------------
!> start/stop button
!--------------------------------------------
subroutine timer(self,i,inmsg)
   
   implicit none
   !> dummy arg list
   class(tb_timer),intent(inout) :: self
      !! instance of timer
   integer,intent(in) :: i
      !! index
   character(len=*),intent(in),optional :: inmsg
      !! raw message text

   if (i > self%n .or. i < 1) return
      !! check if appropriate index is given
   
   !> switcher between start/stop status
   if (self%running(i)) then
      call self%stop_timing(i)
   else
      call self%start_timing(i)
   endif
   self%running(i) = .not.self%running(i)
      !! update status
   
   if (present(inmsg)) self%tag(i) = trim(inmsg)
      !! assign tag to specific timer 

end subroutine timer

!--------------------------------------------
!> To start counting
!--------------------------------------------
subroutine start_timing(self,i)
   
   implicit none
   !> dummy arg list
   class(tb_timer),intent(inout) :: self
      !! instance of timer
   integer,intent(in) :: i
      !! specific timer index 
   
   !> local vars
   real(wp) :: time_cpu
   real(wp) :: time_wall
   
   call timing(time_cpu,time_wall)
   self%tcpu (i) = self%tcpu (i) - time_cpu
   self%twall(i) = self%twall(i) - time_wall

end subroutine start_timing

!--------------------------------------------
!> To stop counting
!--------------------------------------------
subroutine stop_timing(self,i)
   
   implicit none
   !> dummy arg list
   class(tb_timer),intent(inout) :: self
      !! instance of timer
   integer,intent(in) :: i
      !! specific timer index
   
   !> local vars
   real(wp) :: time_cpu
   real(wp) :: time_wall
   
   call timing(time_cpu,time_wall)
   self%tcpu (i) = self%tcpu (i) + time_cpu
   self%twall(i) = self%twall(i) + time_wall

end subroutine stop_timing

!--------------------------------------------
!> To retrieve the current CPU and wall time
!--------------------------------------------
subroutine timing(time_cpu,time_wall)
   
   implicit none
   !> dummy arg list
   real(wp),intent(out) :: time_cpu
   real(wp),intent(out) :: time_wall
   
   !> local vars
   integer(int64) :: time_count
      !! current value of system clock (time passed from arbitary point)
   integer(int64) :: time_rate
      !! number of clock ticks per second (conversion factor b/n ticks and seconds)
   integer(int64) :: time_max
   
   call system_clock(time_count,time_rate,time_max)
   call cpu_time(time_cpu)

   time_wall = real(time_count,wp)/real(time_rate,wp)
      !! elapsed time in seconds

end subroutine timing

end module xtb_type_timer
