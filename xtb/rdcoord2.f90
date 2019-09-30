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

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! read coordinates in au or Ang. if its a xmol file
! redone by S.E. to avoid some input errors. Looks for $coord, ang,
! bohr or number (xmol) in the first line
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine rdcoord(fname,n,xyz,iat)
   use iso_fortran_env, only: wp => real64,output_unit
   use mctc_econv
   use mctc_systools
   use set_module
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
      call raise('E','Coordinate format not recognized!',1)
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
      write(output_unit,'(i0,1x,''/='',1x,i0)') n,ncheck
      call raise('E','reading coord file failed',1)
   endif
   call close_file(ich)

   !  321 FORMAT(F20.10,F20.10,F20.10,1X,A3,1X,A3,1X,A3,I3,L) !debug output
end subroutine rdcoord

subroutine rdatomnumber(fname,n)
   use iso_fortran_env, only : wp => real64
   use mctc_econv, only : aatoau
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
      call raise('E','Coordinate format not recognized!',1)
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
   use iso_fortran_env, only: wp => real64,output_unit
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
            write(output_unit, &
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
            write(output_unit, &
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
   use iso_fortran_env, only : wp => real64
   use mctc_econv, only : aatoau
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
   use iso_fortran_env, only : wp => real64
   use mctc_econv
   use mctc_systools
   use set_module
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
      call set_chrg(trim(chdum))
      if (idum.ne.ichrg) then
         call raise('S','sdf input attempted to set charge, but variables is already locked',1)
      endif
   endif

   call close_file(ich)

end subroutine rdsdf

subroutine pbcrdcoord(fname,lattice,n,xyz,iat)
  use mctc_strings
  use mctc_systools
  use set_module
  implicit none             

  real*8                :: xyz(3,*)
  real*8, INTENT(OUT)   ::lattice(3,3)
  integer, INTENT(out)               :: iat(*) 
  integer, INTENT(in)               :: n 
  character*(*), INTENT(IN)          :: fname
  logical              :: selective=.FALSE. ! Selective dynamics
  logical              :: cartesian=.TRUE.  ! Cartesian or direct
  real*8, parameter :: autoang = 0.52917726d0

  real*8 xx(10),scalar
  character*200 line
  character*80 args(90),args2(90)

  integer i,j,ich,nn,ntype,ntype2,atnum,i_dummy1,i_dummy2,ncheck

  call set_geopref('poscar')

  lattice=0

  ich=142
  call open_file(ich,fname,'r')
  rewind(ich)
  ncheck=0
  ntype=0
  read(ich,'(a)',end=200)line !first line must contain Element Info
  call parse(line,' ',args,ntype)
  read(ich,'(a)',end=200)line !second line contains global scaling factor
  call readl(line,xx,nn)
  scalar=xx(1)/autoang        !the Ang->au conversion is included in the scaling factor
  !c      write(*,'(F8.6)')scalar
  DO i=1,3            ! reading the lattice constants
     read(ich,'(a)',end=200)line
     call readl(line,xx,nn)
     IF (nn < 3) call raise('E', 'Error reading unit cell vectors' ,1)
     lattice(1,i)=xx(1)*scalar
     lattice(2,i)=xx(2)*scalar
     lattice(3,i)=xx(3)*scalar
     !  write(*,'(3F6.2)')lattice(1,i),lattice(2,i),lattice(3,i)
  ENDDO  
  read(ich,'(a)',end=200)line !Ether here are the numbers of each element, or (>vasp.5.1) here are the element symbols
  line=adjustl(line)
  call readl(line,xx,nn)
  IF (nn.eq.0) then      ! CONTCAR files have additional Element line here since vasp.5.1
     call parse(line,' ',args,ntype)
     read(ich,'(a)',end=200)line
     line=adjustl(line)
     call readl(line,xx,nn)
  ENDIF
  !       call elem(args(1),i_dummy2)
  !       IF (i_dummy2<1 .OR. i_dummy2>94) THEN
  !          args=args2
  !       ENDIF
  IF (nn.NE.ntype ) THEN
     call raise('E', 'Error reading number of atomtypes',1)
  ENDIF
  ncheck=0
  DO i=1,nn
     i_dummy1=INT(xx(i))
     call elem(args(i),i_dummy2)
     IF (i_dummy2<1 .OR. i_dummy2>94) call raise('E', 'Error: unknown element.',1)
     DO j=1,i_dummy1
        ncheck=ncheck+1
        iat(ncheck)=i_dummy2
     ENDDO
  ENDDO
  if (n.ne.ncheck) call raise('E','Error reading Number of Atoms',1)

  read(ich,'(a)',end=200)line
  line=adjustl(line)
  IF (line(:1).EQ.'s' .OR. line(:1).EQ.'S') THEN
     selective=.TRUE.
     read(ich,'(a)',end=200)line
     line=adjustl(line)
  ENDIF

  !c      write(*,*)line(:1)
  cartesian=(line(:1).EQ.'c' .OR. line(:1).EQ.'C' .OR. &
       &line(:1).EQ.'k' .OR. line(:1).EQ.'K')
  DO i=1,n
     read(ich,'(a)',end=200)line
     call readl(line,xx,nn)
     IF (nn.NE.3) call raise('E', 'Error reading coordinates.',1)

     IF (cartesian) THEN
        xyz(1,i)=xx(1)*scalar
        xyz(2,i)=xx(2)*scalar
        xyz(3,i)=xx(3)*scalar
     ELSE
        xyz(1,i)=lattice(1,1)*xx(1)+lattice(1,2)*&
             &    xx(2)+lattice(1,3)*xx(3)
        xyz(2,i)=lattice(2,1)*xx(1)+lattice(2,2)*xx(2)+lattice(2,3)*&
             &    xx(3)
        xyz(3,i)=lattice(3,1)*xx(1)+lattice(3,2)*xx(2)+lattice(3,3)*&
             &    xx(3)
     ENDIF

     !c      write(*,'(3F20.10,1X,I3)')xyz(:,i),iat(i)   !debug printout

  ENDDO


200 continue

  call close_file(ich)
end subroutine pbcrdcoord

subroutine pbcrdatomnumber(fname,n)
  use mctc_strings
  use mctc_systools
  implicit none

  integer, INTENT(out)               :: n 
  character*(*), INTENT(IN)          :: fname
  logical              :: selective=.FALSE. ! Selective dynamics
  logical              :: cartesian=.TRUE.  ! Cartesian or direct

  real*8 xx(10),scalar,fdum
  character*80 line,args(90),args2(90)

  integer i,j,ich,nn,ntype,ntype2,atnum,i_dummy1,i_dummy2

  ich=142
  call open_file(ich,fname,'r')
  n=0
  ntype=0
  read(ich,'(a)',end=200)line !first line must contain Element Info
  call parse(line,' ',args,ntype)
  read(ich,'(a)',end=200)line !second line contains global scaling factor
  call readl(line,xx,nn)
  !c      write(*,'(F8.6)')scalar
  DO i=1,3            ! reading the lattice constants
     read(ich,'(a)',end=200)line
     call readl(line,xx,nn)
     IF (nn < 3) call raise('E', 'Error reading unit cell vectors' ,1)
     !  write(*,'(3F6.2)')lattice(1,i),lattice(2,i),lattice(3,i)
  ENDDO
  read(ich,'(a)',end=200)line !Ether here are the numbers of each element, or (>vasp.5.1) here are the element symbols
  line=adjustl(line)
  call readl(line,xx,nn)
  IF (nn.eq.0) then      ! CONTCAR files have additional Element line here since vasp.5.1
     call parse(line,' ',args,ntype)
     read(ich,'(a)',end=200)line
     line=adjustl(line)
     call readl(line,xx,nn)
  ENDIF
  !       call elem(args(1),i_dummy2)
  !       IF (i_dummy2<1 .OR. i_dummy2>94) THEN
  !          args=args2
  !       ENDIF
  IF (nn.NE.ntype ) THEN
     !         IF(nn.NE.ntype2) THEN
     call raise('E', 'Error reading number of atomtypes',1)
     !         ELSE
     !           ntype=ntype2
     !         ENDIF
  ENDIF
  n=0
  DO i=1,nn
     i_dummy1=INT(xx(i))
     n=n+i_dummy1
  ENDDO

200 continue

  call close_file(ich)
end subroutine pbcrdatomnumber

