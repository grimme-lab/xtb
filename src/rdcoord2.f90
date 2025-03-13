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

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! read coordinates in au or Ang. if its a xmol file
! redone by S.E. to avoid some input errors. Looks for $coord, ang,
! bohr or number (xmol) in the first line
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine rdcoord(fname,n,xyz,iat)
   use xtb_mctc_global, only : persistentEnv
   use xtb_mctc_io, only : stdout
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   use xtb_mctc_systools
   use xtb_setmod
   implicit none

   character(len=*),intent(in) :: fname
   integer, intent(inout) :: n
   real(wp),intent(out) :: xyz(3,n)
   integer, intent(out) :: iat(n)

   real(wp) :: floats(3),f
   character(len=128) :: line
   character(len=80)  :: strings(3)
   integer  :: j,k,ich,cs,cf,ncheck

   if(index(fname,'.sdf').ne.0)then
      call set_geopref('sdf')
      call rdsdf(fname,n,xyz,iat) ! sdf is not allowed to overwrite charge
      return
   endif

   f=0.0_wp
   call open_file(ich,fname,'r')
   ncheck=0
   rewind(ich)
   do
      read(ich,'(a)',end=200)line
      if(line.ne."") exit
   enddo

   call readline(line,floats,strings,cs,cf)
   if(cf.eq.1.and.floats(1).gt.0) then
      call set_geopref('xmol')
      f=aatoau
      call getline(ich,commentline)
   else if (index(line,'$coord').ne.0) then
      call set_geopref('coord')
      f=1.0_wp
   else if (index(line,'ang').ne.0) then
      f=aatoau
   else if (index(line,'bohr').ne.0) then
      f=1.0_wp
   endif
   if(f.lt.1.0_wp) then
      call persistentEnv%error('Coordinate format not recognized!')
      return
   endif
   do
      read(ich,'(a)',end=200)line
!     print*,trim(line)
      if(index(line,'$').ne.0) exit
      !        if(index(line,'$redu').ne.0) exit
      !        if(index(line,'$user').ne.0) exit
      !        if(index(line,'$end' ).ne.0) exit
      !        if(index(line,'$set' ).ne.0) exit
      call readline(line,floats,strings,cs,cf)
      if(cf.ne.3) cycle
      !        call readl(line,floats,k)
      call elem(strings(1),j)
      if(j.le.0) cycle !ignores dummies and unknown elements
      ncheck=ncheck+1
      xyz(1,ncheck)=floats(1)*f
      xyz(2,ncheck)=floats(2)*f
      xyz(3,ncheck)=floats(3)*f
      iat(ncheck)=j
   enddo

   200 continue

   if (n.ne.ncheck) then
      write(stdout,'(i0,1x,''/='',1x,i0)') n,ncheck
      call persistentEnv%error('reading coord file failed')
   endif
   call close_file(ich)

   !  321 FORMAT(F20.10,F20.10,F20.10,1X,A3,1X,A3,1X,A3,I3,L) !debug output
end subroutine rdcoord

subroutine rdatomnumber(fname,n)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   implicit none
   character(len=*) :: fname
   integer,intent(out) :: n

   real(wp) :: floats(3),f
   character(len=128) :: line
   character(len=128) :: strings(3)
   integer :: j,ich,cs,cf

   f=0.0_wp
   call open_file(ich,fname,'r')
   if(index(fname,'.sdf').ne.0)then
      read(ich,'(a)')
      read(ich,'(a)')
      read(ich,'(a)')
      read(ich,'(i3)') n
      return
   endif
   n=0
   do
      read(ich,'(a)',end=200)line
      if(line.ne."") exit
   enddo
   call readline(line,floats,strings,cs,cf)
   if(cf.eq.1.and.floats(1).gt.0.and.cs.eq.0) then
      f=aatoau
      !         write(*,*)floats(1)
      n=int(floats(1))
      call close_file(ich)
      return
   else if (index(line,'$coord').ne.0) then
      f=1.0_wp
   else if (index(line,'ang').ne.0) then
      f=aatoau
   else if (index(line,'bohr').ne.0) then
      f=1.0_wp
   endif
   if(f.lt.1.0_wp) then
      write(*,*) f
      call raise('E','Coordinate format not recognized!')
   endif
   DO
      read(ich,'(a)',end=200)line
      if(index(line,'$').ne.0) exit
      !        if(index(line,'$redu').ne.0) exit
      !        if(index(line,'$user').ne.0) exit
      !        if(index(line,'$end' ).ne.0) exit
      call readline(line,floats,strings,cs,cf)
      if(cf.ne.3) exit
      call elem(strings(1),j)
      if(j.eq.0) cycle
      n=n+1
   ENDDO

   200 continue

   call close_file(ich)

   !  321 FORMAT(F20.10,F20.10,F20.10,1X,A3,1X,A3,1X,A3,I3,L) !debug output
end subroutine rdatomnumber


!reads a line cuts the at blanks and tabstops and returns all floats and strings in order of occurence
subroutine readline(line,floats,strings,cs,cf)
   use xtb_mctc_io, only : stdout
   use xtb_mctc_accuracy, only : wp
   implicit none
   real(wp) :: floats(3)
   character(len=*),intent(in) :: line
   character(len=80)  :: strings(3)

   real(wp) :: num
   character(len=80) :: stmp,str
   character(len=1)  :: digit
   integer  :: i,ty,cs,cf

   stmp=''
   cs=1
   cf=1
   strings(:)=''
   do i=1,len(trim(line))
      digit=line(i:i)
      if(digit.ne.' '.and.digit.ne.char(9)) then  !should exclude tabstops and blanks, 9 is ascii code for tab
         stmp=trim(stmp)//trim(digit)
      elseif(stmp.ne.'')then
         call checktype(stmp,num,str,ty)      !get type of string, 0=number, 1=character
         if(ty.eq.0) then
            floats(cf)=num
            cf=cf+1
         elseif(ty.eq.1) then
            strings(cs)=str
            cs=cs+1
         else
            write(stdout, &
               & '(''readline: problem in checktype, must abort'')')
            exit
         endif
         stmp=''
      endif
      if(i.eq.len(trim(line))) then  !special case: end of line
         call checktype(stmp,num,str,ty)
         if(ty.eq.0) then
            floats(cf)=num
            cf=cf+1
         elseif(ty.eq.1) then
            strings(cs)=str
            cs=cs+1
         else
            write(stdout, &
               & '(''readline: problem in checktype, must abort'')')
            exit
         endif
         stmp=''
      endif
   enddo
   cs=cs-1
   cf=cf-1
end subroutine readline

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! read coordinates, xyz only
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine rdxyz(fname,n,xyz)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   implicit integer (i-n)
   implicit real(wp) (a-h,o-z)
   integer, intent(in)  :: n
   real(wp),intent(out) :: xyz(3,n)
   character(len=128)   :: line
   character(len=*),intent(in) :: fname
   real(wp) :: xx(10)

   call open_file(ich,fname,'r')

   if(index(fname,'xyz').ne.0)then
      read(ich,'(a)')line
      read(ich,'(a)')line
      do i=1,n
         read(ich,'(a)')line
         call readl(line,xx,nn)
         xyz(1:3,i)=xx(1:3)
      enddo
      xyz = xyz * aatoau
   else
      read(ich,'(a)')line
      do i=1,n
         read(ich,'(a)')line
         call readl(line,xx,nn)
         xyz(1:3,i)=xx(1:3)
      enddo
   endif

   call close_file(ich)
end subroutine rdxyz

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! read sdf coordinate file
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine rdsdf(fname,n,xyz,iat)
   use xtb_mctc_global, only : persistentEnv
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   use xtb_mctc_systools
   use xtb_setmod
   implicit none
   integer, intent(inout) :: n
   real(wp),intent(out) :: xyz(3,n)
   integer, intent(out) :: iat(n)
   character(len=128)   :: line
   character(len=*),intent(in)  :: fname
   character(len=:),allocatable :: wfile

   real(wp) :: xx(20)
   character(len=20)   :: chdum
   integer  :: nn,i,j,ich,idum
   logical  :: exist

   idum=0
   call open_file(ich,fname,'r')

   call getline(ich,molnameline)
   read(ich,'(a)')!skip
   call getline(ich,commentline)
   read(ich,'(i3)') n

   do i=1,n
      read(ich,'(a)')line
      call readl(line,xx,nn)
      call elem(line,j)
      iat(i) = j
      xyz(1:3,i)=xx(1:3)
   enddo

   xyz = xyz * aatoau

   100 read(ich,'(a)',end=200)line
   if(index(line,'$$$$').ne.0) goto 200
   if(index(line,'M  CHG').ne.0) then
      call readl(line,xx,nn)
      idum=idum+idint(xx(nn))
   endif
   goto 100
   200 continue
   if (idum.ne.0) then
      write(chdum,'(i0)') idum
      call set_chrg(persistentEnv, trim(chdum))
      if (idum.ne.set%ichrg) then
         call persistentEnv%warning('sdf input attempted to set charge, '//&
            &'but variables is already locked')
      endif
   endif

   call close_file(ich)

end subroutine rdsdf
