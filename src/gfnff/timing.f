ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      subroutine timing(t,w)
      use iso_fortran_env, only : wp => real64
         real(wp) :: t,w
         integer :: time_count,time_rate,time_max
         call system_clock(time_count,time_rate,time_max)
         call cpu_time(t)
         w = real(time_count,wp)/real(time_rate,wp)
      end subroutine timing


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine prtime(io,tt,ww,string)
      integer io
      real*8 ww,tt,t,tsec,wsec
      integer tday,thour,tmin
      integer wday,whour,wmin
      character*(*) string

      t=tt
      tday=idint(t/86400)
      t=t-tday*86400
      thour=idint(t/3600)
      t=t-thour*3600
      tmin=idint(t/60)
      t=t-60*tmin
      tsec=t

      t=ww
      wday=idint(t/86400)
      t=t-wday*86400
      whour=idint(t/3600)
      t=t-whour*3600
      wmin=idint(t/60)
      t=t-60*wmin
      wsec=t

      if(tday.ne.0)then
         write(io,2 )trim(string),tday,thour,tmin,tsec
         write(io,20)trim(string),wday,whour,wmin,wsec
         return
      endif
      if(thour.ne.0)then
         write(io,3 )trim(string),thour,tmin,tsec
         write(io,30)trim(string),whour,wmin,wsec
         return
      endif
      if(tmin .ne.0)then
         write(io,4 )trim(string),tmin,tsec
         write(io,40)trim(string),wmin,wsec
         return
      endif
      write(io,5 )trim(string),tsec
      write(io,50)trim(string),wsec
      return

 2    format('cpu  time for ',a,2x,i3,' d',i3,' h',i3,' m',f5.1,' s')
 3    format('cpu  time for ',a,2x,i3,' h',i3,' m',f5.1,' s')
 4    format('cpu  time for ',a,2x,i3,' m',f5.1,' s')
 5    format('cpu  time for ',a,2x,f6.2,' s')
 20   format('wall time for ',a,2x,i3,' d',i3,' h',i3,' m',f5.1,' s')
 30   format('wall time for ',a,2x,i3,' h',i3,' m',f5.1,' s')
 40   format('wall time for ',a,2x,i3,' m',f5.1,' s')
 50   format('wall time for ',a,2x,f6.2,' s')

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
