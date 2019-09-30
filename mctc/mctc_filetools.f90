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

!! ------------------------------------------------------------------------
!  From the Fortran 2003 Standard
!! ------------------------------------------------------------------------
!  OPEN(UNIT=integer,           -- unit specifier
!       ACCESS=character,       -- Specifies access mode for the connection
!                                  values: 'direct', 'sequential' (default),
!                                          'stream'
!                                  if ACCESS is 'direct', RECL is required
!       ACTION=character,       -- Specifies the actions allowed for the connection
!                                  values: 'read', 'write', 'readwrite'
!       ASYNCHRONOUS=character, -- Specifiers whether or not asynchronous
!                                  input/output is allowed
!                                  values: 'yes', 'no' (default)
!       BLANK=character,        -- values: 'null' (default), 'zero'
!       DECIMAL=character,      -- values: 'comma', 'point' (default)
!       DELIM=character,        -- values: 'apostrophe', 'quote', 'none' (default)
!       ENCODING=character,     -- Specifies the character encoding method
!                                  for formatted files
!                                  values: 'utf-8', 'default' (default)
!       ERR=label,              -- 
!       FILE=character,         -- Specifies the name of the file being connected
!                                  rules:
!                                  1. Trailing blanks in the name are ignored.
!                                  2. The name must be a file name allowed by the 
!                                     processor.
!                                  3. The interpretation of case is processor 
!                                     dependent; for example, the processor
!                                     may distinguish file names by case or
!                                     it may interpret the name all in uppercase
!                                     or lowercase letters.
!                                  4. The FILE specifier must appear if the STATUS
!                                     value is 'new' or 'replace'.
!                                  5. If the STATUS value is 'old', the FILE
!                                     specifier must appear unless the unit is
!                                     already connected to a unit that exists.
!                                  6. If the STATUS specifier has the value
!                                     'scratch', the FILE specifier must not
!                                     appear.
!       FORM=character,         -- Specifies formatted or unformatted input/output
!                                  values: 'formated', 'unformatted'
!                                  rules:
!                                  1. The default value is 'unformatted', if the
!                                     file is connected for direct or stream access.
!                                  2. The default value is 'formatted', if the file 
!                                     is connected for sequential access.
!       IOMSG=character,        -- 
!       IOSTAT=integer,         -- 
!       PAD=character,          -- values: 'yes' (default), 'no'
!       POSITION=character,     -- Specifies the initial position of the file
!                                  values: 'asis' (default), 'rewind', 'append'
!                                  rules:
!                                  1. The default value is 'asis', permitting an
!                                     OPEN statement to change other connection
!                                     properties of a file that is already
!                                     connected without changing its position.
!                                  2. The file must be connected for sequential or
!                                     stream access.
!                                  3. If the file is new, it is positioned at its
!                                     initial point.
!       RECL=integer,           -- Specifies record length for direct- or
!                                  sequential-access files. The expression
!                                  specifies the length of each direct-access
!                                  record or the maximum record length if the
!                                  access method is sequential.
!       ROUND=character,        -- values: 'up', 'down', 'zero', 'nearest',
!                                          'compatible', 'processor_defined'
!       SIGN=character,         -- values: 'plus', 'suppress', 'processor_defined'
!       STATUS=character)       -- Specifies the file existence before the
!                                  connection is made
!                                  values:
!                                 'old'     -- file must exist
!                                 'new'     -- file must not exist, will be created
!                                 'unknown' -- default
!                                 'replace' -- if file does not exist, it is
!                                              created and given STATUS 'old';
!                                              if file does exist, it is deleted,
!                                              and a new file is created with the
!                                              same name and STATUS 'old'.
!                                 'scratch' -- file is created and connected to the
!                                              specified unit; it is deleted when
!                                              the program terminates normally
!                                              or a CLOSE statement is executed on
!                                              that unit (file must be unnamed)
!! ------------------------------------------------------------------------
module mctc_filetools
implicit none

type :: filelog
   character(len=:),allocatable :: name
   character(len=1) :: status
   integer :: unit
   logical :: open
end type filelog

type(filelog),allocatable :: filelist(:)
integer :: nfiles = 0

contains

subroutine init_filelist(maxfiles)
   integer,intent(in) :: maxfiles
   nfiles = 0
   if (allocated(filelist)) deallocate(filelist)
   allocate( filelist(maxfiles) )
end subroutine init_filelist

subroutine reallocate_filelist
   integer :: current_size,new_size
   type(filelog),allocatable :: tmplist(:)
   current_size = size(filelist,1)
   new_size = current_size+current_size/2+1
   allocate( tmplist(new_size) )
   tmplist(:current_size) = filelist
   deallocate( filelist )
   call move_alloc(tmplist,filelist)
end subroutine reallocate_filelist

subroutine push_file(unit,name,status)
   implicit none
   character(len=*),intent(in)  :: name
   character(len=1),intent(in)  :: status
   integer,         intent(in)  :: unit
   character(len=:),allocatable :: tmp

   nfiles = nfiles + 1
   if (nfiles .gt. size(filelist,1)) call reallocate_filelist

   filelist(nfiles)%name = name
   filelist(nfiles)%status = status
   filelist(nfiles)%unit = unit
   if (status.eq.'d' .or. status.eq.'t') then
      filelist(nfiles)%open = .false.
   else
      filelist(nfiles)%open = .true.
   endif
end subroutine push_file

subroutine pop_file(unit,status)
   implicit none
   integer,intent(in) :: unit
   character(len=1),intent(in),optional :: status
   logical :: found
   integer :: i
   found = .false.

   do i = 1, nfiles
      found = unit.eq.filelist(i)%unit .and. filelist(i)%open
      if (found) then
         if (present(status)) filelist(i)%status = status
         filelist(i)%open = .false.
         exit
      endif
   enddo
end subroutine pop_file

end module mctc_filetools
