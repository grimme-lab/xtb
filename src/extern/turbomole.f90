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

subroutine external_turbomole(n,at,xyz,nel,nopen,grd,eel,g,dip,lsolv)
   use xtb_mctc_accuracy, only : wp
   use xtb_setparam
   implicit none
   integer n, at(n), nel, nopen
   logical grd,lsolv
   real(wp) xyz(3,n)
   real(wp) xyz_cached(3,n)
   real(wp) g  (3,n)
   real(wp) eel
   real(wp) dip(3)
   character(len=255) atmp
   logical :: cache, exist

   cache = .false.
   dip=0


   ! TM (RI)
   if(extcode.eq.1)then
      !$omp critical (turbo_lock)
      inquire(file='gradient', exist=exist)
      if (exist .and. grd) then
         call rdtm(n,grd,eel,g,xyz_cached)
         cache = all(abs(xyz_cached - xyz) < 1.e-10_wp)
      end if
      if (.not.cache) then
         call wrtm(n,at,xyz)
         if(extmode.eq.1)then
            call execute_command_line('exec ridft  >  job.last 2>> /dev/null')
            if(grd)call execute_command_line('exec rdgrad >> job.last 2>> /dev/null')
         endif
         call extcodeok(extcode)
         call rdtm(n,grd,eel,g,xyz_cached)
      end if
      !$omp end critical (turbo_lock)
      return
   endif

   ! TM+d3+gcp
   if(extcode.eq.2)then
      !$omp critical (turbo_lock)
      inquire(file='gradient', exist=exist)
      if (exist .and. grd) then
         call rdtm(n,grd,eel,g,xyz_cached)
         cache = all(abs(xyz_cached - xyz) < 1.e-10_wp)
      end if
      if (.not.cache) then
         call wrtm(n,at,xyz)
         if(extmode.le.2)then
            call execute_command_line('exec ridft  >  job.last 2>> /dev/null')
            call execute_command_line('exec rdgrad >> job.last 2>> /dev/null')
            call execute_command_line('exec dftd3 coord -grad >> job.last 2>> /dev/null')
            call execute_command_line('exec gcp coord -file -grad >>job.last 2>>/dev/null')
         endif
         call extcodeok(extcode)
         call rdtm(n,.true.,eel,g,xyz_cached)
      end if
      !$omp end critical (turbo_lock)
      return
   endif

   ! TM (NORI)
   if(extcode.eq.3)then
      !$omp critical (turbo_lock)
      inquire(file='gradient', exist=exist)
      if (exist .and. grd) then
         call rdtm(n,grd,eel,g,xyz_cached)
         cache = all(abs(xyz_cached - xyz) < 1.e-10_wp)
      end if
      if (.not.cache) then
         call wrtm(n,at,xyz)
         if(extmode.eq.1)then
            call execute_command_line('exec dscf  > job.last 2>> /dev/null')
            if(grd)call execute_command_line('exec grad >> job.last 2>> /dev/null')
         endif
         call extcodeok(extcode)
         call rdtm(n,grd,eel,g,xyz_cached)
      end if
      !$omp end critical (turbo_lock)
      return
   endif


   call raise('E','This external code is not implemented')

end subroutine external_turbomole

!ccccccccccccccccccccccccccccccccc
! TM
!ccccccccccccccccccccccccccccccccc

subroutine wrtm(n,at,xyz)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_symbols, only : toSymbol
   implicit none
   integer n, at(n), iunit, i
   real(wp) xyz(3,n)

   iunit=33
   open(newunit=iunit,file='coord')

   write(iunit,'(a)')'$coord'
   do i=1,n
      write(iunit,'(3F24.14,6x,a2)') &
      &   xyz(1,i),xyz(2,i),xyz(3,i),toSymbol(at(i))
   enddo
   write(iunit,'(a)')'$end'

   close(iunit)

end subroutine wrtm

subroutine rdtm(n,grd,e,g,xyz)
   use xtb_mctc_accuracy, only : wp
   implicit none
   integer n, iunit, i, nl, j, nn
   logical grd
   real(wp) g(3,n), e, xx(10), x, y, z, xyz(3,n)
   logical ex
   character(len=128) a1

   iunit=33

   if(.not.grd)then
      open(newunit=iunit,file='energy')
      101   read(iunit,'(a)',end=102)a1
      call readl(a1,xx,nn)
      if(nn.ge.4) e=xx(2)
      goto 101
      102   close(iunit)
      return
   endif

   inquire(file='gradient',exist=ex)
   if(.not.ex) then
      call raise('E','no gradient file found!')
   endif

   j=0
   open(newunit=iunit,file='gradient')
201 read(iunit,'(a)',end=301)a1
   j=j+1
   if(index(a1,'cycle').ne.0)nl=j
   goto 201
   301   continue

   if(nl.lt.2)then
      call raise('E','illegal gradient file!')
   endif

   rewind iunit
   do i=1,nl
      read(iunit,'(a)')a1
   enddo
   call readl(a1,xx,nn)
   e=xx(2)
   do i=1,n
      read(iunit,*)xyz(1,i),xyz(2,i),xyz(3,i)
   enddo
   do i=1,n
      read(iunit,*)g(1,i),g(2,i),g(3,i)
   enddo

   close(iunit)

end subroutine rdtm

subroutine rijcheck(extcode)
   implicit none
   integer, intent (inout) :: extcode
   integer icheck,ich
   icheck=1
   call execute_command_line('exec sdg rij | wc -l > TmPfIlE')
   open(newunit=ich,file='TmPfIlE',status='old')
   read(ich,*)icheck
   close(ich,status='delete')
   if(icheck.lt.1) extcode=3
   return
end subroutine rijcheck

subroutine extcodeok(extcode)
   implicit none
   integer, intent (in) :: extcode
   integer :: ich
   character(len=80) atmp
   ! TM
   if(extcode.le.3)then
      call execute_command_line('exec grep "actual step" control > TmPfIlE')
      open(newunit=ich,file='TmPfIlE',status='old')
      read(ich,'(a)',end=100)atmp
 100  close(ich,status='delete')
      if(index(atmp,'actual').ne.0) call raise('E','external code error: '//&
      &                                            trim(atmp))
   endif

   return
end subroutine extcodeok
