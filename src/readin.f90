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

module xtb_readin
   
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_strings, only : value
   use xtb_type_environment, only : TEnvironment
   implicit none

   character,private,parameter :: flag = '$'
   character,private,parameter :: space = ' '
   character,private,parameter :: equal = '='
   character,private,parameter :: hash = '#'
   character,private,parameter :: dot = '.'
   character,private,parameter :: comma = ','
   character,private,parameter :: minus = '-'
   character(len=*),private,parameter :: flag_end = '$end'

!! ------------------------------------------------------------------[SAW]-
!  this function returns a logical and is always evaluated for its
!  side effect (parsing the given string for its real/int/bool value)
   interface getValue
      module procedure getIntValue
      module procedure getIntArray
      module procedure getRealValue
      module procedure getRealArray
      module procedure getBoolValue
   end interface getValue

contains

!! ------------------------------------------------------------------[SAW]-
!  I could use rdpath directly, but this would require access to xpath,
!  so I use xfind as a wrapper with access to the xpath variable to
!  circumvent this. Also as a gimmick, I do not return a logical, but
!  some file name. xfind succeeds if fname.ne.name, but if you inquire
!  for fname in case of failure, you might hit a local file, which you
!  than can read. This is intended as a feature (just saying).
function xfind(name) result(fname)
   use xtb_mctc_systools, only : rdpath
   use xtb_setparam, only : xenv
   character(len=*),intent(in)  :: name
   character(len=:),allocatable :: fname
   character(len=:),allocatable :: dum
   logical :: exist

   call rdpath(xenv%path,name,dum,exist)
   if (exist) then
      fname = dum
   else
      fname = name
   endif

end function xfind

!! ------------------------------------------------------------------[SAW]-
!>  wrapper around getline from the MCTC lib that strips comments
!>  automatically und removes all leading and trailing whitespace
subroutine strip_line(in,line,err)
   
   use xtb_mctc_systools, only : getline
   implicit none
   integer,intent(in)  :: in
      !! input file unit
   character(len=:),allocatable,intent(out) :: line
      !! output line
   integer,intent(out) :: err
   integer :: ic

   call getline(in,line,iostat=err)
   if (err.ne.0) return
   
   !> check for comment characters
   ic = index(line,hash)
   if (ic.eq.1) then
      line = ''
      return
   else if (ic.gt.1) then
      line = line(:ic-1)
   endif
   !> to remove all leading and trailing whitespaces
   line = trim(adjustl(line))

end subroutine strip_line

!! ------------------------------------------------------------------[SAW]-
!  same as strip_line, but has the additional function of copying to
!  one unit while reading from another, which is helpful for backing up
!  files you plan to replace in the next step of you program.
!  Funnily this subroutine exist way before strip_line...
subroutine mirror_line(in,out,line,err)
   
   use xtb_mctc_systools, only : getline
   implicit none
   
   integer,intent(in)  :: in
   integer,intent(in)  :: out
   character(len=:),allocatable,intent(out) :: line
   integer,intent(out) :: err
   
   integer :: ic

   call getline(in,line,iostat=err)
      !! returns the line from input file 
   if (err.ne.0) return
      !! iostat=-1 will return err=0
!  now write the line to the copy, we write what we read, not what we see
!  if (out.ne.-1) write(out,'(a)') trim(line)
   
   !> Check for comment characters
   ic = index(line,hash)
   if (ic.eq.1) then
      line = ''
      return
   else if (ic.gt.1) then
      line = line(:ic-1)
   endif
!  now write the line to the copy, we write what we see, not what we read
   if (out.ne.-1) write(out,'(a)') trim(line)
!  strip line from space, but after printing
   line = trim(adjustl(line))

end subroutine mirror_line

function find_new_name(fname) result(newname)
   character(len=*),intent(in)  :: fname
   character(len=:),allocatable :: newname
   character(len=5) :: dum ! five digit number must be enough
   character,parameter :: hash = '#'
   character,parameter :: dot = '.'
   integer :: i,idot
   logical :: exist

   idot = index(fname,dot)
   i = 0
   do
      i = i+1
      write(dum,'(i0)') i
      if (idot == 1) then
         newname = dot//hash//trim(dum)//fname
      else
         newname = hash//trim(dum)//dot//fname
      endif
      inquire(file=newname,exist=exist)
      if (.not.exist) return
   enddo

end function find_new_name

function getIntValue(env,val,dum) result(status)
   implicit none
   character(len=*), parameter :: source = 'readin_getIntValue'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: val
   integer,intent(out) :: dum
   integer :: err
   logical :: status
   
!  call value(val,dum,ios=err)
   read(val,*,iostat=err) dum
   if (err.eq.0) then
      status = .true.
   else
      call env%warning('could not parse '''//val//'''',source)
      status = .false.
   endif
end function getIntValue

function getRealValue(env,val,dum) result(status)
   implicit none
   character(len=*), parameter :: source = 'readin_getRealValue'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: val
   real(wp),intent(out) :: dum
   integer :: err
   logical :: status
   
!  call value(val,dum,ios=err)
   read(val,*,iostat=err) dum
   if (err.eq.0) then
      status = .true.
   else
      call env%warning('could not parse '''//val//'''',source)
      status = .false.
   endif
end function getRealValue

function getBoolValue(env,val,dum) result(status)
   implicit none
   character(len=*), parameter :: source = 'readin_getBoolValue'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: val
   logical,intent(out) :: dum
   logical :: status
   
   select case(val)
   case('Y','y','Yes','yes','T','t','true','True','1')
      status = .true.
      dum = .true.
   case('N','n','No','no','F','f','false','False','0')
      status = .true.
      dum = .false.
   case default
      call env%warning('could not parse '''//val//'''',source)
      status = .false.
   end select

end function getBoolValue

function getIntArray(env,val,dum) result(status)
   implicit none
   character(len=*), parameter :: source = 'readin_getIntArray'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: val
   integer,intent(out) :: dum(:)
   integer :: i,err
   logical :: status
  
!  call value(val,dum,ios=err)
   read(val,*,iostat=err) (dum(i),i=1,size(dum,1))
   if (err.eq.0) then
      status = .true.
   else
      call env%warning('could not parse '''//val//'''',source)
      status = .false.
   endif

end function getIntArray

function getRealArray(env,val,dum) result(status)
   implicit none
   character(len=*), parameter :: source = 'readin_getRealArray'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: val
   real(wp),intent(out) :: dum(:)
   integer :: i,err
   logical :: status
  
!  call value(val,dum,ios=err)
   read(val,*,iostat=err) (dum(i),i=1,size(dum,1))
   if (err.eq.0) then
      status = .true.
   else
      call env%warning('could not parse '''//val//'''',source)
      status = .false.
   endif

end function getRealArray

pure elemental function bool2int(bool) result(int)
   logical,intent(in) :: bool
   integer :: int
   if (bool) then
      int = 1
   else
      int = 0
   endif
end function bool2int

pure function bool2string(bool) result(string)
   logical,intent(in) :: bool
   character(len=:),allocatable :: string
   if (bool) then
      string = 'true'
   else
      string = 'false'
   endif
end function bool2string

function getListValue(env,val,dum,n) result(status)
   implicit none
   character(len=*), parameter :: source = 'readin_getListValue'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: val
   integer,intent(out) :: dum(:)
   integer,intent(out) :: n
   integer :: i,j,k,l,err
   logical :: status
  
   i = index(val,minus)
   if (i.eq.0) then
      read(val,*,iostat=err) dum(1)
      if (err.ne.0) then
         call env%warning('could not parse '''//val//'''',source)
         status = .false.
         return
      endif
      n = 1
      status = .true.
   else
      read(val(:i-1),*,iostat=err) j
      if (err.ne.0) then
         call env%warning('could not parse '''//val(:i-1)//''' in '''//val//'''',source)
         status = .false.
         return
      endif
      read(val(i+1:),*,iostat=err) k
      if (err.ne.0) then
         call env%warning('could not parse '''//val(i+1:)//''' in '''//val//'''',source)
         status = .false.
         return
      endif
      if (k.lt.j) then
         call env%warning('end is lower than start in list '''//val//'''',source)
         status = .false.
         return
      endif
      if ((k-j).gt.size(dum,1)) then
         call env%warning('too many list items in '''//val//'''',source)
         status = .false.
         return
      endif
      n = 0
      do i = j, k
         n = n+1
         dum(n) = i
      enddo
      status = .true.
   endif
end function getListValue

subroutine readlog(fname,nat,at,xyz,nstruc)
   use xtb_mctc_convert
   use xtb_mctc_systools
   implicit none
   character(len=*),intent(in) :: fname
   integer, intent(in)    :: nat
   integer, intent(in)    :: at(nat)
   real(wp),intent(out)   :: xyz(3,nat,nstruc)
   integer, intent(inout) :: nstruc

   character(len=:),allocatable :: line
   real(wp) :: floats(3),f
   character(len=80)  :: strings(3)
   integer  :: j,k,cs,cf,ncheck
   integer  :: ich ! file handle
   integer  :: idum,i,err

   call open_file(ich,fname,'r')
   if (ich.eq.-1) then
      call raise('E',"Could not find '"//fname//"'!")
   endif

   j = 0
   read_struc: do
      read(ich,*,iostat=err) idum
      if (err /= 0) exit read_struc
      if (idum /= nat) then
         call raise('E',"Atom number missmatch in '"//fname//"'!")
      endif
      if (j .ge. nstruc) exit read_struc
      read(ich,'(a)')
      if (err /= 0) exit read_struc
      do i = 1, nat
         call getline(ich,line,err)
         if (err /= 0) exit read_struc
         call readline(line,floats,strings,cs,cf)
         xyz(:,i,j+1)=floats*aatoau
      enddo
      j = j+1
   enddo read_struc
   nstruc = j

   call close_file(ich)

end subroutine readlog


end module xtb_readin
