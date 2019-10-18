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

!! ========================================================================
!  * WELCOME TO THE   C O N S T R A I N S   &   S C A N S   MODULE IN XTB *
!! ------------------------------------------------------------------------
!  The syntax is:
!> $constrain
!> ## this part is read by setparam.f90
!>    force constant=<real>
!>    all bonds=<bool>
!>    all angles=<bool>
!>    all torsions=<bool>
!> ## this part is read in this module
!>    distance: <i>,<j>,        auto
!>    distance: <i>,<j>,        <real>
!>    angle:    <i>,<j>,<k>,    auto
!>    angle:    <i>,<j>,<k>,    <real>
!>    dihedral: <i>,<j>,<k>,<l>,auto
!>    dihedral: <i>,<j>,<k>,<l>,<real>
!>    center:   <int>,<real>
!> $fix
!>    force constant=<real>
!>    spring exponent=<int>
!>    atoms: <list>
!> $split
!>    fragment1: <list>
!>    fragment2: <list>
!> $wall
!>    sphere: auto,all
!>    sphere: auto,<list>
!>    sphere: <real>,all
!>    sphere: <real>,<list>
!>    ellipsoid: auto,all
!>    ellipsoid: auto,<list>
!>    ellipsoid: <real>,<real>,<real>,all
!>    ellipsoid: <real>,<real>,<real>,<list>
!> $scan
!>    ...
!> $end
!! ========================================================================
module constrain_param
   use iso_fortran_env, wp => real64

   use mctc_strings, only : parse
   use readin, only : getline => strip_line,get_value,get_list_value
   use setparam, only : verbose

   implicit none

   private :: getline,get_value

   character,private,parameter :: flag = '$'
   character,private,parameter :: colon = ':'
   character,private,parameter :: space = ' '
   character,private,parameter :: equal = '='
   character,private,parameter :: hash = '#'
   character,private,parameter :: comma = ','
   character,private,parameter :: semicolon = ';'
   character,private,parameter :: dot = '.'
   character(len=*),private,parameter :: flag_end = flag//'end'

!  Using allocatable arrays of dynamic length strings is only possible
!  with a lot of hacks, so we use good'ol fixed size stack arrays.
!  Let's choose something different from 42 that is not dividable by 10... ;)
!  Happy debugging!
   integer,private,parameter :: p_str_length = 48
   integer,private,parameter :: p_arg_length = 24

   public

contains

subroutine read_userdata(fname,nat,at,xyz)
   use iso_fortran_env, only : output_unit,iostat_end
   use readin, only : find_new_name
   use scanparam
   implicit none
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   character(len=*),intent(in)  :: fname
   character(len=:),allocatable :: line
   character(len=:),allocatable :: key
   character(len=:),allocatable :: val
   character(len=:),allocatable :: newname
   integer :: i
   integer :: id
   integer :: ic
   integer :: ie
   integer :: err
   logical :: exist

   if (verbose) then
      write(output_unit,'(72("$"))')
      write(output_unit,'(1x,"CONSTRAINTS & SCANS: DEBUG SECTION")')
      write(output_unit,'(72("$"))')
   endif

   call open_file(id,fname,'r')
   if (id.eq.-1) then
      call raise('S',"could not find '"//fname//"'",1)
      return
   endif
   rewind(id) ! not sure if this is necessary

!  read first line before the readloop starts, I have to do this
!  to avoid using backspace on id (dammit Turbomole format)
   call getline(id,line,err)
   readflags: do
      !  check if there is a $ in the *first* column
      if (index(line,flag).eq.1) then
         select case(line(2:))
         case('fix'      )
            if (verbose) write(output_unit,'(">",1x,a)') line(2:)
            call rdblock(set_fix,    line,id,nat,at,xyz,err)
         case('split'    )
            if (verbose) write(output_unit,'(">",1x,a)') line(2:)
            call rdblock(set_split,  line,id,nat,at,xyz,err)
         case('constrain')
            if (verbose) write(output_unit,'(">",1x,a)') line(2:)
            if (allocated(potset%xyz)) then
               call rdblock(set_constr, line,id,nat,at,potset%xyz,err)
            else
               call rdblock(set_constr, line,id,nat,at,xyz,err)
            endif
         case('scan'     )
            if (verbose) write(output_unit,'(">",1x,a)') line(2:)
            call rdblock(set_scan,   line,id,nat,at,xyz,err)
         case('wall'     )
            if (verbose) write(output_unit,'(">",1x,a)') line(2:)
            call rdblock(set_wall,   line,id,nat,at,xyz,err)
         case('metadyn'  )
            if (verbose) write(output_unit,'(">",1x,a)') line(2:)
            call rdblock(set_metadyn,line,id,nat,at,xyz,err)
         case('hess'     )
            if (verbose) write(output_unit,'(">",1x,a)') line(2:)
            call rdblock(set_hess,   line,id,nat,at,xyz,err)
         case('path'     )
            if (verbose) write(output_unit,'(">",1x,a)') line(2:)
            call rdblock(set_path,   line,id,nat,at,xyz,err)
         case('reactor'  )
            if (verbose) write(output_unit,'(">",1x,a)') line(2:)
            call rdblock(set_reactor,line,id,nat,at,xyz,err)
         case('set'      ); call rdsetbl(set_legacy,line,id,nat,at,xyz,err)
         case default ! unknown keyword -> ignore, we don't raise them
            call getline(id,line,err)
         end select
      else ! not a keyword -> ignore
         call getline(id,line,err)
      endif
   !  check for end of file, which I will tolerate as alternative to $end
      if (err.eq.iostat_end) exit readflags
!     if (index(line,flag_end).ne.0) exit readflags ! compatibility reasons
   enddo readflags

   if (verbose) write(output_unit,'(72("$"))')
   call close_file(id)
end subroutine read_userdata

subroutine rdsetbl(handler,line,id,nat,at,xyz,err)
   use iso_fortran_env, only : iostat_end
   implicit none
   interface
      subroutine handler(key,val,nat,at,xyz)
      use iso_fortran_env, only : wp => real64
      character(len=*),intent(in) :: key
      character(len=*),intent(in) :: val
      integer, intent(in) :: nat
      integer, intent(in) :: at(nat)
      real(wp),intent(in) :: xyz(3,nat)
      end subroutine handler
   end interface
   integer,intent(in) :: id
   integer,intent(in) :: nat
   integer,intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   external :: handler
   integer,intent(out) :: err
   character(len=:),allocatable,intent(out) :: line
   character(len=:),allocatable :: key
   character(len=:),allocatable :: val
   integer :: ie
   do
      call getline(id,line,err)
      if (err.eq.iostat_end) exit
      if (index(line,flag).ne.0) exit
      if (verbose) write(output_unit,'("->",1x,a)') line

      ! find the first colon
      ie = index(line,space)
      if ((line.eq.'').or.(ie.eq.0)) cycle
      key = trim(line(:ie-1))
      val = trim(adjustl(line(ie+1:)))
      call handler(key,val,nat,at,xyz)
   enddo

end subroutine rdsetbl

subroutine rdblock(handler,line,id,nat,at,xyz,err)
   use iso_fortran_env, only : iostat_end
   implicit none
   interface
      subroutine handler(key,val,nat,at,xyz)
      use iso_fortran_env, only : wp => real64
      character(len=*),intent(in) :: key
      character(len=*),intent(in) :: val
      integer, intent(in) :: nat
      integer, intent(in) :: at(nat)
      real(wp),intent(in) :: xyz(3,nat)
      end subroutine handler
   end interface
   integer,intent(in) :: id
   integer,intent(in) :: nat
   integer,intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   external :: handler
   integer,intent(out) :: err
   character(len=:),allocatable,intent(out) :: line
   character(len=:),allocatable :: key
   character(len=:),allocatable :: val
   integer :: ie
   do
      call getline(id,line,err)
      if (err.eq.iostat_end) exit
      if (index(line,flag).ne.0) exit
      if (verbose) write(output_unit,'("->",1x,a)') line

      ! find the first colon
      ie = index(line,colon)
      if ((line.eq.'').or.(ie.eq.0)) cycle
      key = trim(line(:ie-1))
      val = trim(adjustl(line(ie+1:)))
      call handler(key,val,nat,at,xyz)
   enddo

end subroutine rdblock

subroutine set_fix(key,val,nat,at,xyz)
   use tbdef_atomlist
   use fixparam
   use setparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   type(tb_atomlist) :: atl

   integer  :: i
   integer  :: iat
   integer  :: idum
   integer  :: nlist
   integer, allocatable :: list(:)
   real(wp) :: ddum
   logical  :: ldum

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   call atl%resize(nat)

   call parse(val,comma,argv,narg)
!  some debug printout
   if (verbose) then
      do idum = 1, narg
         write(output_unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('elements')
      call atl%new
      do idum = 1, narg
         ! get element by symbol
         call elem(argv(idum),iat)
         ! alternatively, try ordinal number
         if (iat.eq.0) then
            ldum = get_value(trim(argv(idum)),iat)
            if (.not.ldum) cycle ! skip garbage input
         endif
         ! check for unreasonable input
         if (iat.le.0) then
            call raise('S',"Unknown element: '"//trim(argv(idum))//"'",1)
            cycle
         endif
         ! now find the elements in the geometry
         call atl%add(at.eq.iat)
      enddo
      if (fixset%n > 0) call atl%add(fixset%atoms(:fixset%n))
      call atl%to_list(list)
      fixset%atoms = list
      fixset%n = size(list)
   case('atoms')
      call atl%new(val)
      if (atl%get_error()) then
         call raise('S','something is wrong in the fixing list',1)
         return
      endif
      if (fixset%n > 0) call atl%add(fixset%atoms(:fixset%n))
      call atl%to_list(list)
      fixset%atoms = list
      fixset%n = size(list)
   case('freeze')
      call atl%new(val)
      if (atl%get_error()) then
         call raise('S','something is wrong in the freezing list',1)
         return
      endif
      if (freezeset%n > 0) call atl%add(freezeset%atoms(:freezeset%n))
      call atl%to_list(list)
      freezeset%atoms = list
      freezeset%n = size(list)
   case('shake')
      allocate(list(nat*(nat+1)/2), source=0)
      if (mod(narg,2).ne.0) then
         call raise('S',"could not read input for user defined shake!",1)
         return
      endif
         if (narg+shakeset%n > nat*(nat+1)/2) then
         call raise('S',"too many SHAKE constraints!",1)
         return
      endif
      if (.not.shake_md) shake_md = .true.
      do idum = 1, narg
         if (get_value(trim(argv(idum)),iat)) then
            if (iat.gt.nat) then
               call raise('S','Attempted constrain atom not present in molecule.',1)
               cycle
            endif
            shakeset%n = shakeset%n+1
            shakeset%atoms(shakeset%n) = iat
         else
            call raise('S',"Something went wrong in set_fix_ 'shake'.",1)
            return ! you screwed it, let's get out of here
         endif
      enddo
   end select
 
end subroutine set_fix

subroutine set_constr(key,val,nat,at,xyz)
   use mctc_constants
   use mctc_econv
   use tbdef_atomlist
   use scanparam
   use splitparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   type(tb_atomlist) :: atl

   integer  :: iat
   integer  :: ioffset
   integer  :: idum
   real(wp) :: ddum
   integer  :: nlist
   integer, allocatable :: list(:)
   logical  :: ldum
   integer  :: i,j,k,l
   real(wp) :: phi,dist,ra(3),rb(3)
   real(wp),external :: valijkl

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   call atl%resize(nat)

   call parse(val,comma,argv,narg)
   if (verbose) then
      do idum = 1, narg
         write(output_unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them

   case('elements')
      call atl%new
      do idum = 1, narg
         ! get element by symbol
         call elem(argv(idum),iat)
         ! alternatively, try ordinal number
         if (iat.eq.0) then
            ldum = get_value(trim(argv(idum)),iat)
            if (.not.ldum) cycle ! skip garbage input
         endif
         ! check for unreasonable input
         if (iat.le.0) then
            call raise('S',"Unknown element: '"//trim(argv(idum))//"'",1)
            cycle
         endif
         ! now find the elements in the geometry
         call atl%add(at.eq.iat)
      enddo
      if (potset%pos%n > 0) call atl%add(potset%pos%atoms(:potset%pos%n))
      call atl%to_list(list)
      potset%pos%atoms = list
      potset%pos%n = size(list)
   case('atoms')
      call atl%new(val)
      if (atl%get_error()) then
         call raise('S','something is wrong in the fixing list',1)
         return
      endif
      if (potset%pos%n > 0) call atl%add(potset%pos%atoms(:potset%pos%n))
      call atl%to_list(list)
      potset%pos%atoms = list
      potset%pos%n = size(list)

   case('DISTANCE')
      if (narg.ne.3) then
         call raise('E','not enough arguments to constrain a distance',1)
         return
      endif
      ioffset = 2*potset%dist%n
      potset%dist%n = potset%dist%n+1
   !  part 1: get the constrained atoms
      do i = 1, narg-1
         if (get_value(trim(argv(i)),idum)) then
            potset%dist%atoms(ioffset+i) = idum
         else
            call raise('S',"Something went wrong in set_constr_ 'distance'. (1)",1)
            potset%dist%n = potset%dist%n-1 ! remove invalid contrain
            return
         endif
      enddo
   !  part 2: get the distance between those atoms
      i = potset%dist%atoms(ioffset+1) 
      j = potset%dist%atoms(ioffset+2)
      dist = norm2(xyz(:,i)-xyz(:,j))
      if (trim(argv(narg)).eq.'auto') then
         potset%dist%val(potset%dist%n) = dist
      else
         if (get_value(trim(argv(narg)),ddum)) then
            potset%dist%val(potset%dist%n) = ddum * aatoau
         else
            call raise('S',"Something went wrong in set_constr_ 'distance'. (2)",1)
            potset%dist%n = potset%dist%n-1 ! remove invalid contrain
            return
         endif
      endif

   case('ANGLE')
      if (narg.ne.4) then
         call raise('E','not enough arguments to constrain an angle',1)
         return
      endif
      ioffset = 3*potset%angle%n
      potset%angle%n = potset%angle%n+1
   !  part 1: get the constrained atoms
      do i = 1, narg-1
         if (get_value(trim(argv(i)),idum)) then
            potset%angle%atoms(ioffset+i) = idum
         else
            call raise('S',"Something went wrong in set_constr_ 'angle'. (1)",1)
            potset%angle%n = potset%angle%n-1 ! remove invalid contrain
            return
         endif
      enddo
   !  part 2: get the angle between the constrained atoms
      i = potset%angle%atoms(ioffset+1)
      j = potset%angle%atoms(ioffset+2)
      k = potset%angle%atoms(ioffset+3)
      call bangl(xyz,i,j,k,phi)
      if (trim(argv(narg)).eq.'auto') then
         potset%angle%val(potset%angle%n) = phi
      else
         if (get_value(trim(argv(narg)),ddum)) then
            potset%angle%val(potset%angle%n) = pi/180.0_wp * ddum
         else
            call raise('S',"Something went wrong in set_constr_ 'angle'. (2)",1)
            potset%angle%n = potset%angle%n-1 ! remove invalid contrain
            return
         endif
      endif

   case('DIHEDRAL')
      if (narg.ne.5) then
         call raise('E','not enough arguments to constrain a dihedral',1)
         return
      endif
      ioffset = 4*potset%dihedral%n
      potset%dihedral%n = potset%dihedral%n+1
      if (nconstr.gt.maxconstr) & ! double check this
      &  call raise('E','This should never happen! Let somebody check set_constr',1)
   !  part 1: get the constrained atoms
      do i = 1, narg-1
         if (get_value(trim(argv(i)),idum)) then
            potset%dihedral%atoms(ioffset+i) = idum
         else
            call raise('S',"Something went wrong in set_constr_ 'dihedral'. (1)",1)
            potset%dihedral%n = potset%dihedral%n-1 ! remove invalid contrain
            return
         endif
      enddo
   !  part 2: get the angle between the constrained atoms
      i = potset%dihedral%atoms(ioffset+1)
      j = potset%dihedral%atoms(ioffset+2)
      k = potset%dihedral%atoms(ioffset+3)
      l = potset%dihedral%atoms(ioffset+4)
      phi=valijkl(nat,xyz,i,j,k,l)
      if (trim(argv(narg)).eq.'auto') then
         potset%dihedral%val(potset%dihedral%n) = phi
      else
         if (get_value(trim(argv(narg)),ddum)) then
            potset%dihedral%val(potset%dihedral%n) = pi/180.0_wp * ddum
         else
            call raise('S',"Something went wrong in set_constr_ 'dihedral'. (2)",1)
            potset%dihedral%n = potset%dihedral%n-1 ! remove invalid contrain
            return
         endif
      endif

   case('distance')
      if (narg.ne.3) then
         call raise('E','not enough arguments to constrain a distance',1)
         return
      endif
      nconstr = nconstr+1
      if (nconstr.gt.maxconstr) & ! double check this
      &  call raise('E','This should never happen! Let somebody check set_constr',1)
   !  part 1: get the constrained atoms
      do i = 1, narg-1
         if (get_value(trim(argv(i)),idum)) then
            atconstr(i,nconstr) = idum
         else
            call raise('S',"Something went wrong in set_constr_ 'distance'. (1)",1)
            nconstr = nconstr-1 ! remove invalid contrain
            return
         endif
      enddo
   !  part 2: get the distance between those atoms
      i = atconstr(1,nconstr) 
      j = atconstr(2,nconstr)
      dist = norm2(xyz(:,i)-xyz(:,j))
      if (trim(argv(narg)).eq.'auto') then
         valconstr(nconstr) = dist
      else
         if (get_value(trim(argv(narg)),ddum)) then
            valconstr(nconstr) = ddum * aatoau
         else
            call raise('S',"Something went wrong in set_constr_ 'distance'. (2)",1)
            nconstr = nconstr-1 ! remove invalid contrain
            return
         endif
      endif
      write(output_unit,'("constraining bond",2(1x,i0),1x,"to",'//&
         '1x,f12.7,1x,"Å, actual value:",1x,f12.7,1x,"Å")') &
         i,j, valconstr(nconstr)*autoaa, dist*autoaa

   case('angle')
      if (narg.ne.4) then
         call raise('E','not enough arguments to constrain an angle',1)
         return
      endif
      nconstr = nconstr+1
      if (nconstr.gt.maxconstr) & ! double check this
      &  call raise('E','This should never happen! Let somebody check set_constr',1)
   !  part 1: get the constrained atoms
      do i = 1, narg-1
         if (get_value(trim(argv(i)),idum)) then
            atconstr(i,nconstr) = idum
         else
            call raise('S',"Something went wrong in set_constr_ 'angle'. (1)",1)
            nconstr = nconstr-1 ! remove invalid contrain
            return
         endif
      enddo
   !  part 2: get the angle between the constrained atoms
      i = atconstr(1,nconstr)
      j = atconstr(2,nconstr)
      k = atconstr(3,nconstr)
      call bangl(xyz,i,j,k,phi)
      if (trim(argv(narg)).eq.'auto') then
         valconstr(nconstr) = phi
      else
         if (get_value(trim(argv(narg)),ddum)) then
            valconstr(nconstr) = pi/180.0_wp * ddum
         else
            call raise('S',"Something went wrong in set_constr_ 'angle'. (2)",1)
            nconstr = nconstr-1 ! remove invalid contrain
            return
         endif
      endif
      write(output_unit,'("constraining angle",3(1x,i0),1x,"to",'//&
         '1x,f12.7,"°, actual value:",1x,f12.7,"°")') &
         i,j,k,180.0_wp/pi * valconstr(nconstr),180.0_wp/pi * phi

   case('dihedral')
      if (narg.ne.5) then
         call raise('E','not enough arguments to constrain a dihedral',1)
         return
      endif
      nconstr = nconstr+1
      if (nconstr.gt.maxconstr) & ! double check this
      &  call raise('E','This should never happen! Let somebody check set_constr',1)
   !  part 1: get the constrained atoms
      do i = 1, narg-1
         if (get_value(trim(argv(i)),idum)) then
            atconstr(i,nconstr) = idum
         else
            call raise('S',"Something went wrong in set_constr_ 'dihedral'. (1)",1)
            nconstr = nconstr-1 ! remove invalid contrain
            return
         endif
      enddo
   !  part 2: get the angle between the constrained atoms
      i = atconstr(1,nconstr)
      j = atconstr(2,nconstr)
      k = atconstr(3,nconstr)
      l = atconstr(4,nconstr)
      phi=valijkl(nat,xyz,i,j,k,l)
      if (trim(argv(narg)).eq.'auto') then
         valconstr(nconstr) = phi
      else
         if (get_value(trim(argv(narg)),ddum)) then
            valconstr(nconstr) = pi/180.0_wp * ddum
         else
            call raise('S',"Something went wrong in set_constr_ 'dihedral'. (2)",1)
            nconstr = nconstr-1 ! remove invalid contrain
            return
         endif
      endif
      write(output_unit,'("constraining angle",4(1x,i0),1x,"to",'//&
         '1x,f12.7,"°, actual value:",1x,f12.7,"°")') &
         i,j,k,l,180.0_wp/pi * valconstr(nconstr),180.0_wp/pi * phi

   case('center')
      ! copied from molbld.f without modification (except for error handling)
      if (narg.ne.2) then
         call raise('E','not enough argument to constrain center',1)
         return
      endif
      nconstr = nconstr+1
      atconstr(1,nconstr) = -2
      if (get_value(trim(argv(2)),ddum)) then
         valconstr(nconstr) = ddum
      else
         call raise('S',"Something went wrong in set_constr_ 'center'. (1)",1)
         nconstr = nconstr-1
         return
      endif
      if (get_value(trim(argv(1)),idum)) then
         iatf1 = idum
      else
         call raise('S',"Something went wrong in set_constr_ 'center'. (2)",1)
         nconstr = nconstr-1
         return
      endif
      massf1 = 0.0_wp
      do i = 1, iatf1
         massf1 = massf1 + atmass(i)
      enddo

   case('cma','cma interface')
      ! copied from molbld.f without modification (except for error handling)
      if (narg.ne.1) then
         call raise('E','not enough argument to constrain cma',1)
         return
      endif
      if (key.eq.'cma interface') call cmaiface(nat,at,xyz)
      nconstr = nconstr+1
      atconstr(1,nconstr) = 0
      if (trim(argv(1)).eq.'auto') then
         massf1 = 0.0_wp
         massf2 = 0.0_wp
         do i = 1, nat
            if (splitlist(i).eq.1) then
               massf1 = massf1 + atmass(i)
            else
               massf2 = massf2 + atmass(i)
            endif
         enddo
         call cmafrag(nat,at,xyz,ra,rb)
         valconstr(nconstr) = rcma
         write(output_unit,'("constraining fragment CMA to initial R(Ang.)=",f8.3)')&
            valconstr(nconstr)*autoaa
      else
         if (get_value(trim(argv(1)),ddum)) then
            valconstr(nconstr) = ddum*aatoau
         else
            call raise('S',"Something went wrong in set_constr_ 'cma'.",1)
            nconstr = nconstr-1
            return
         endif
      endif

!   case('x','y','z')
!      idum = index('xyz',key)
   case('z')
      ! copied from molbld.f without modification (except for error handling)
      if (narg.ne.1) then
         call raise('E','not enough argument to constrain z coordinate',1)
         return
      endif
      zconstr = 1
      nconstr = nconstr+1
      atconstr(1,nconstr) = -1
      if (get_value(trim(argv(1)),ddum)) then
         valconstr(nconstr) = ddum*aatoau
      else
         call raise('S',"Something went wrong in set_constr_ 'z'.",1)
         nconstr = nconstr-1
         return
      endif
      if (sum(abs(xyz(3,1:iatf1))).gt.1.d-3) then
         call raise('S','z-coordinates of fragment 1 must be 0!',1)
         nconstr = nconstr-1
         return
      endif
      write(output_unit,'("constraining fragment 2 to Z=0 plane at (Ang.)=",f8.3)') &
         valconstr(nconstr)*autoaa

   end select
end subroutine set_constr

!! --------------------------------------------------------------[SAW1809]-
!  this is the new version of the scan routine exploiting all features
subroutine set_scan(key,val,nat,at,xyz)
   use scanparam
   use mctc_econv, only : aatoau
   use mctc_constants, only : pi
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum

   integer  :: i,ie
   real(wp) :: start_value,end_value
   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv
   character(len=:),allocatable :: temp

   nscan = nscan+1 ! make a new scan

   ie = index(val,semicolon)
   if (ie.ne.0) then
      temp = val(:ie-1)
      idum = nconstr
      call set_constr(key,temp,nat,at,xyz) ! generate a new constraint
      scan_list(nscan)%iconstr = nconstr ! new generated
      if (idum.eq.nconstr) call raise('E','Failed to generate constraint',1)
      temp = val(1+ie:)
   else
      if (get_value(key,idum)) then
         if (idum.gt.nconstr) then
            call raise('E','Constraint '''//key//''' is not defined',1)
         endif
         scan_list(nscan)%iconstr = idum
      else
         call raise('E','Constraint '''//key//''' is invalid in this context',1)
      endif
      temp = val
   endif

   call parse(temp,comma,argv,narg)
   if (verbose) then
      do idum = 1, narg
         write(output_unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif

   ! get back the constraint number which should be scanned
   i = scan_list(nscan)%iconstr
   ! instead of simply saving the kind of constraint while reading,
   ! we used some really obvious integer code system depending on
   ! the number of elements in the atconstr-Array.
   ! Expected to break at some point, add a FIXME and wait for complaints...
   idum = atconstr(1,i)
   if (atconstr(2,i).gt.0) idum = 1
   if (atconstr(3,i).gt.0) idum = 2
   if (atconstr(4,i).gt.0) idum = 3

   if (get_value(trim(argv(1)),ddum)) then
      if (idum.le.1) then
         start_value = ddum * aatoau
      else
         start_value = ddum * pi/180.0_wp
      endif
   else
      call raise('E','Invalid start value for scan',1)
   endif
   if (get_value(trim(argv(2)),ddum)) then
      if (idum.le.1) then
         end_value = ddum * aatoau
      else
         end_value = ddum * pi/180.0_wp
      endif
   else
      call raise('E','Invalid end value for scan',1)
   endif
   if (get_value(trim(argv(3)),idum)) then
      scan_list(nscan)%nscan = idum
      allocate( scan_list(nscan)%valscan(idum), source = 0.0_wp )
   else
      call raise('E','Invalid step number for scan',1)
   endif

   do i = 1, scan_list(nscan)%nscan
      scan_list(nscan)%valscan(i) = start_value + (end_value-start_value) &
         &                       * real(i-1,wp)/real(scan_list(nscan)%nscan-1,wp)
      if (verbose) write(output_unit,'(i5,1x,f12.8)') i,scan_list(nscan)%valscan(i)
   enddo

end subroutine set_scan

subroutine set_wall(key,val,nat,at,xyz)
   use mctc_econv, only : autoaa
   use sphereparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   integer  :: idum
   real(wp) :: ddum,darray(3)
   logical  :: ldum
   integer  :: list(nat),nlist
   integer  :: tlist(nat),ntlist
   integer  :: iarg,iaxis
   real(wp) :: radius,radii(3),center(3)

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   nlist = 0
   ntlist = 0
   list = 0
   tlist = 0

   radius = 0.0_wp
   radii  = 0.0_wp
   center = 0.0_wp

   call parse(val,comma,argv,narg)
   if (verbose) then
      do idum = 1, narg
         write(output_unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('sphere')
      if (narg.lt.2) then
         call raise('S',"Not enough arguments to set up a spherical wall",1)
         return
      endif
   !  part 1: get the sphere radius
      if (trim(argv(1)).eq.'auto') then
         call get_sphere_radius(nat,at,xyz,center,radius,do_trafo=.false.)
      else
         if (get_value(trim(argv(1)),ddum)) then
            radius = ddum
            center = 0.0_wp
         else
            ! warning already generated by get_value
            return ! something went wrong
         endif
      endif
   !  part 2: get atoms
      if (trim(argv(2)).eq.'all') then
         call set_sphere_radius(radius,center)
      else
         do iarg = 2, narg
            if (get_list_value(trim(argv(iarg)),tlist,ntlist)) then
               if (nlist+ntlist.gt.nat) then
                  call raise('S',"Too many atoms in list for spherical wall.",1)
                  return ! something went wrong
               endif
               if (maxval(tlist(:ntlist)).gt.nat) then
                  call raise('S',"Attempted to wall in a non-existing atom",1)
                  cycle ! skip crappy input
               endif
               list(nlist+1:nlist+ntlist) = tlist
               nlist = nlist + ntlist
            else
               ! warning already generated by get_list_value
               return ! something went wrong
            endif
         enddo
         call set_sphere_radius(radius,center,nlist,list)
      endif
      write(output_unit,'("spherical wallpotenial with radius",'//&
         '1x,f12.7,1x,"Å")') radius*autoaa

   case('ellipsoid')
      if (narg.lt.4) then
         call raise('S',"Not enough arguments to set up an ellipsoidal wall",1)
         return
      endif
   !  part 1: get ellipsoid axis
      if ((trim(argv(1)).eq.'auto').or.(trim(argv(2)).eq.'auto').or. &
          (trim(argv(3)).eq.'auto')) then
         ! use get_ellipsoid_radius if you manage to implement it
         call get_sphere_radius(nat,at,xyz,center,radius,do_trafo=.false.)
      endif
      do iaxis = 1, 3
         if (trim(argv(iaxis)).eq.'auto') then
            radii(iaxis) = radius
         else
            if (get_value(trim(argv(iaxis)),ddum)) then
               radii(iaxis) = ddum
               center(iaxis) = 0.0_wp
            else
               ! warning already generated by get_value
               return ! something went wrong
            endif
         endif
      enddo
   !  part 2: get atoms
      if (trim(argv(4)).eq.'all') then
         call set_sphere_radius(radii,center)
      else
         do iarg = 4, narg
            if (get_list_value(trim(argv(iarg)),tlist,ntlist)) then
               if (nlist+ntlist.gt.nat) then
                  call raise('S',"Too many atoms in list for spherical wall.",1)
                  return ! something went wrong
               endif
               if (maxval(tlist(:ntlist)).gt.nat) then
                  call raise('S',"Attempted to wall in a non-existing atom",1)
                  cycle ! skip crappy input
               endif
               list(nlist+1:nlist+ntlist) = tlist
               nlist = nlist + ntlist
            else
               ! warning already generated by get_list_value
               return ! something went wrong
            endif
         enddo
         call set_sphere_radius(radii,center,nlist,list)
      endif
      write(output_unit,'("ellipsoidal wallpotenial with radii",'//&
         '3(1x,f12.7,1x,"Å"))') radii*autoaa

   end select
 
end subroutine set_wall

subroutine set_split(key,val,nat,at,xyz)
   use splitparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   integer  :: list(nat),nlist
   integer  :: i,j

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   call parse(val,comma,argv,narg)
   if (verbose) then
      do idum = 1, narg
         write(output_unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('fragment1')
      do i = 1, narg
         if (get_list_value(trim(argv(i)),list,nlist)) then
            do j = 1, nlist
               splitlist(list(j)) = 1
            enddo
            iatf1 = iatf1 + nlist ! compatibility stuff, as usual
            iatf(1) = iatf(1) + nlist
         else
            return ! you screwed it, let's get out of here
         endif
      enddo
   case('fragment2')
      do i = 1, narg
         if (get_list_value(trim(argv(i)),list,nlist)) then
            do j = 1, nlist
               splitlist(list(j)) = 2
            enddo
            iatf2 = iatf2 + nlist ! compatibility stuff, as usual
            iatf(2) = iatf(2) + nlist
         else
            return ! you screwed it, let's get out of here
         endif
      enddo
   case('fragment')
      if (get_value(trim(argv(1)),idum)) then
         if (idum.gt.nat) then
            call raise('S',"rejecting fragment number greater than number of atoms",1)
            return ! doesn't really make sense, sorry, your problem
         endif
         do i = 2, narg
            if (get_list_value(trim(argv(i)),list,nlist)) then
               do j = 1, nlist
                  splitlist(list(j)) = idum
               enddo
               iatf(idum) = iatf(idum) + nlist
               if (idum.eq.1) then ! legacy, don't remove or something breaks
                  iatf1 = iatf1 + nlist
               else if (idum.eq.2) then
                  iatf2 = iatf2 + nlist
               endif
            else
               return ! you screwed it, let's get out of here
            endif
         enddo
      endif
   end select
end subroutine set_split

subroutine set_hess(key,val,nat,at,xyz)
   use splitparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   integer  :: i,j

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   call parse(val,comma,argv,narg)
   if (verbose) then
      do idum = 1, narg
         write(output_unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('element mass')
      if (mod(narg,2).ne.0) then
         call raise('S',"Something went wrong in set_mass_ 'element'.",1)
      endif
      do i = 1, narg, 2
         j = i+1
         if (get_value(trim(argv(i)),idum).and.&
            &get_value(trim(argv(j)),ddum)) then
            where(at.eq.idum) atmass = ddum
            write(output_unit,'(a,1x,i0,1x,a,1x,g0)') &
               'mass of elements ',idum,' changed to', ddum
         endif
      enddo
   case('modify mass','isotope')
      if (mod(narg,2).ne.0) then
         call raise('S',"Something went wrong in set_mass_ 'modify'.",1)
      endif
      do i = 1, narg, 2
         j = i+1
         if (get_value(trim(argv(i)),idum).and.&
             get_value(trim(argv(j)),ddum)) then
            if (idum.gt.nat) then
               call raise('S','Attempted setting atom mass not present in system.',1)
               cycle
            endif
            atmass(idum) = ddum
            write(output_unit,'(a,1x,i0,1x,a,1x,g0)') &
               'mass of atom ',idum,' changed to',atmass(idum)
         endif
      enddo
   case('scale mass')
      if (mod(narg,2).ne.0) then
         call raise('S',"Something went wrong in set_mass_ 'scale'.",1)
      endif
      do i = 1, narg, 2
         j = i+1
         if (get_value(trim(argv(i)),idum).and.&
             get_value(trim(argv(j)),ddum)) then
            if (idum.gt.nat) then
               call raise('S','Attempted scaling atom not present in system.',1)
               cycle
            endif
            atmass(idum) = atmass(idum)*ddum
            write(output_unit,'(a,1x,i0,1x,a,1x,g0)') &
               'mass of atom ',idum,' changed to',atmass(idum)
         endif
      enddo
   end select
 
end subroutine set_hess

subroutine set_reactor(key,val,nat,at,xyz)
   use tbdef_atomlist
   use setparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   type(tb_atomlist) :: atl

   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   integer  :: nlist
   integer, allocatable :: list(:)
   integer  :: i,j

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   call parse(val,comma,argv,narg)
   if (verbose) then
      do idum = 1, narg
         write(output_unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('atoms')
      call atl%new(val)
      if (atl%get_error()) then
         call raise('S','something is wrong in the reactor atom list',1)
         return
      endif
      if (reactset%nat > 0) call atl%add(reactset%atoms(:reactset%nat))
      call atl%to_list(list)
      reactset%atoms = list
      reactset%nat = size(list)
   end select
 
end subroutine set_reactor

subroutine set_path(key,val,nat,at,xyz)
   use tbdef_atomlist
   use setparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   type(tb_atomlist) :: atl

   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   integer  :: nlist
   integer, allocatable :: list(:)
   integer  :: i,j

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   call parse(val,comma,argv,narg)
   if (verbose) then
      do idum = 1, narg
         write(output_unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('atoms')
      call atl%new(val)
      if (atl%get_error()) then
         call raise('S','something is wrong in the bias path atom list',1)
         return
      endif
      if (pathset%nat > 0) call atl%add(pathset%atoms(:pathset%nat))
      call atl%to_list(list)
      pathset%atoms = list
      pathset%nat = size(list)
   end select
 
end subroutine set_path

subroutine set_metadyn(key,val,nat,at,xyz)
   use tbdef_atomlist
   use fixparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   type(tb_atomlist) :: atl

   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   integer  :: nlist
   integer, allocatable :: list(:)
   integer  :: i,j

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   call parse(val,comma,argv,narg)
   if (verbose) then
      do idum = 1, narg
         write(output_unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('atoms')
      call atl%new(val)
      if (atl%get_error()) then
         call raise('S','something is wrong in the metadynamic atom list',1)
         return
      endif
      if (metaset%nat > 0) call atl%add(metaset%atoms(:metaset%nat))
      call atl%to_list(list)
      metaset%atoms = list
      metaset%nat = size(list)
   case('modify factor')
      if (mod(narg,2).ne.0) then
         call raise('S',"Something went wrong in set_metadyn_ 'modify'.",1)
      endif
      do i = 1, narg, 2
         j = i+1
         if (get_value(trim(argv(i)),idum).and.&
             get_value(trim(argv(j)),ddum)) then
            if (idum.gt.metaset%maxsave) then
               call raise('S','Attempted using factor not present in system.',1)
               cycle
            endif
            metaset%factor(idum) = ddum
         endif
      enddo
   case('scale factor')
      if (mod(narg,2).ne.0) then
         call raise('S',"Something went wrong in set_metadyn_ 'scale'.",1)
      endif
      do i = 1, narg, 2
         j = i+1
         if (get_value(trim(argv(i)),idum).and.&
             get_value(trim(argv(j)),ddum)) then
            if (idum.gt.metaset%maxsave) then
               call raise('S','Attempted scaling factor not present in system.',1)
               cycle
            endif
            metaset%factor(idum) = metaset%factor(idum)*ddum
         endif
      enddo
   end select
 
end subroutine set_metadyn

subroutine set_freeze(key,val,nat,at,xyz)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)

   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   integer  :: list(nat),nlist
   integer  :: i,j

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   call parse(val,comma,argv,narg)
   if (verbose) then
      do idum = 1, narg
         write(output_unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('hessf')
   case('hessa')
   end select
 
end subroutine set_freeze

subroutine set_legacy(key,val,nat,at,xyz)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   select case(key)
   case default ! complaining about unknown keywords should already be done
      continue  ! so we do nothing here
   case('hessf');       call set_fix('freeze',val,nat,at,xyz)
!   case('hessa');       call set_frozh('hessa',val)
   case('fragment1'); call set_split('fragment1',val,nat,at,xyz)
   case('fragment2'); call set_split('fragment1',val,nat,at,xyz)
   case('constrxyz'); call set_fix('atoms',val,nat,at,xyz)
!   case('constrainel')
!   case('constrain')
!   case('scan')
   case('ellips'); call set_wall('ellipsoid',val,nat,at,xyz)
   case('sphere'); call set_wall('sphere',val,nat,at,xyz)
   case('fix'); call set_fix('atoms',val,nat,at,xyz)
   case('atomlist+'); call set_metadyn('atoms',val,nat,at,xyz)
   end select

end subroutine set_legacy

end module constrain_param
