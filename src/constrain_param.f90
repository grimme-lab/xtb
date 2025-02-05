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
!>    sandwich: auto,all
!>    sandwich: auto,<list>
!>    sandwich: <real>,all
!>    sandwich: <real>,<list>
!> $scan
!>    ...
!> $end
!! ========================================================================
module xtb_constrain_param
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_strings, only : parse
   use xtb_readin, only : getline => strip_line,getValue,getListValue
   use xtb_setparam, only : set
   use xtb_type_environment, only : TEnvironment
   use xtb_type_identitymap, only : TIdentityMap, init
   use xtb_type_molecule, only : TMolecule

   implicit none

   private :: getline,getValue, handlerInterface

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
   integer,private,parameter :: p_arg_length = 96

   public

   abstract interface
      subroutine handlerInterface(env,key,val,nat,at,idMap,xyz)
         import :: wp, TEnvironment, TIdentityMap
         type(TEnvironment), intent(inout) :: env
         character(len=*),intent(in) :: key
         character(len=*),intent(in) :: val
         integer, intent(in) :: nat
         type(TIdentityMap), intent(in) :: idMap
         integer, intent(in) :: at(nat)
         real(wp),intent(in) :: xyz(3,nat)
      end subroutine handlerInterface
   end interface

contains

subroutine read_userdata(fname,env,mol)
   use xtb_readin, only : find_new_name
   use xtb_scanparam
   implicit none
   character(len=*), parameter :: source = 'userdata_read'
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(inout) :: mol
   character(len=*),intent(in)  :: fname
   character(len=:),allocatable :: line
   character(len=:),allocatable :: key
   character(len=:),allocatable :: val
   character(len=:),allocatable :: newname
   type(TIdentityMap) :: idMap
   integer :: i
   integer :: id
   integer :: ic
   integer :: ie
   integer :: err
   logical :: exist

   if (set%verbose) then
      write(env%unit,'(72("$"))')
      write(env%unit,'(1x,"CONSTRAINTS & SCANS: DEBUG SECTION")')
      write(env%unit,'(72("$"))')
   endif

   call open_file(id,fname,'r')
   if (id.eq.-1) then
      call env%warning("could not find '"//fname//"'",source)
      return
   endif
   rewind(id) ! not sure if this is necessary

   call init(idMap, mol)

   call idMap%writeInfo(env%unit)

!  read first line before the readloop starts, I have to do this
!  to avoid using backspace on id (dammit Turbomole format)
   call getline(id,line,err)
   readflags: do
      !  check if there is a $ in the *first* column
      if (index(line,flag).eq.1) then
         select case(line(2:))
         case('fix'      )
            if (set%verbose) write(env%unit,'(">",1x,a)') line(2:)
            call rdblock(env,set_fix,    line,id,mol%n,mol%at,idMap,mol%xyz,err)
         case('split'    )
            if (set%verbose) write(env%unit,'(">",1x,a)') line(2:)
            call rdblock(env,set_split,  line,id,mol%n,mol%at,idMap,mol%xyz,err)
         case('constrain')
            if (set%verbose) write(env%unit,'(">",1x,a)') line(2:)
            if (allocated(potset%xyz)) then
               call rdblock(env,set_constr, line,id,mol%n,mol%at,idMap,potset%xyz,err)
            else
               call rdblock(env,set_constr, line,id,mol%n,mol%at,idMap,mol%xyz,err)
            endif
         case('scan'     )
            if (set%verbose) write(env%unit,'(">",1x,a)') line(2:)
            call rdblock(env,set_scan,   line,id,mol%n,mol%at,idMap,mol%xyz,err)
         case('wall'     )
            if (set%verbose) write(env%unit,'(">",1x,a)') line(2:)
            call rdblock(env,set_wall,   line,id,mol%n,mol%at,idMap,mol%xyz,err)
         case('metadyn'  )
            if (set%verbose) write(env%unit,'(">",1x,a)') line(2:)
            call rdblock(env,set_metadyn,line,id,mol%n,mol%at,idMap,mol%xyz,err)
         case('hess'     )
            if (set%verbose) write(env%unit,'(">",1x,a)') line(2:)
            call rdblock(env,set_hess,   line,id,mol%n,mol%at,idMap,mol%xyz,err)
         case('path'     )
            if (set%verbose) write(env%unit,'(">",1x,a)') line(2:)
            call rdblock(env,set_path,   line,id,mol%n,mol%at,idMap,mol%xyz,err)
         case('reactor'  )
            if (set%verbose) write(env%unit,'(">",1x,a)') line(2:)
            call rdblock(env,set_reactor,line,id,mol%n,mol%at,idMap,mol%xyz,err)
         case('set'      ); call rdsetbl(env,set_legacy,line,id,mol%n,mol%at,idMap,mol%xyz,err)
         case default ! unknown keyword -> ignore, we don't raise them
            call getline(id,line,err)
         end select
      else ! not a keyword -> ignore
         call getline(id,line,err)
      endif
   !  check for end of file, which I will tolerate as alternative to $end
      if (is_iostat_end(err)) exit readflags
!     if (index(line,flag_end).ne.0) exit readflags ! compatibility reasons
   enddo readflags

   if (set%verbose) write(env%unit,'(72("$"))')
   call close_file(id)
end subroutine read_userdata

subroutine rdsetbl(env,handler,line,id,nat,at,idMap,xyz,err)
   implicit none
   character(len=*), parameter :: source = 'userdata_rdsetbl'
   type(TEnvironment), intent(inout) :: env
   integer,intent(in) :: id
   procedure(handlerInterface) :: handler
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   type(TIdentityMap), intent(in) :: idMap
   real(wp),intent(in) :: xyz(3,nat)
   integer,intent(out) :: err
   character(len=:),allocatable,intent(out) :: line
   character(len=:),allocatable :: key
   character(len=:),allocatable :: val
   integer :: ie
   logical :: exitRun
   do
      call getline(id,line,err)
      if (is_iostat_end(err)) exit
      if (index(line,flag).ne.0) exit
      if (set%verbose) write(env%unit,'("->",1x,a)') line

      ! find the first colon
      ie = index(line,space)
      if ((line.eq.'').or.(ie.eq.0)) cycle
      key = trim(line(:ie-1))
      val = trim(adjustl(line(ie+1:)))
      call handler(env,key,val,nat,at,idMap,xyz)
      call env%check(exitRun)
      if (exitRun) then
         call env%error("handler could not process input", source)
         return
      end if
   enddo

end subroutine rdsetbl

subroutine rdblock(env,handler,line,id,nat,at,idMap,xyz,err)
   implicit none
   character(len=*), parameter :: source = 'userdata_rdblock'
   type(TEnvironment), intent(inout) :: env
   integer,intent(in) :: id
   procedure(handlerInterface) :: handler
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   type(TIdentityMap), intent(in) :: idMap
   real(wp),intent(in) :: xyz(3,nat)
   integer,intent(out) :: err
   character(len=:),allocatable,intent(out) :: line
   character(len=:),allocatable :: key
   character(len=:),allocatable :: val
   integer :: ie
   logical :: exitRun
   do
      call getline(id,line,err)
      if (is_iostat_end(err)) exit
      if (index(line,flag).ne.0) exit
      if (set%verbose) write(env%unit,'("->",1x,a)') line

      ! find the first colon
      ie = index(line,colon)
      if ((line.eq.'').or.(ie.eq.0)) cycle
      key = trim(line(:ie-1))
      val = trim(adjustl(line(ie+1:)))
      call handler(env,key,val,nat,at,idMap,xyz)
      call env%check(exitRun)
      if (exitRun) then
         call env%error("handler could not process input", source)
         return
      end if
   enddo

end subroutine rdblock

subroutine set_fix(env,key,val,nat,at,idMap,xyz)
   use xtb_type_atomlist
   use xtb_fixparam
   use xtb_setparam
   implicit none
   character(len=*), parameter :: source = 'userdata_fix'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   type(TIdentityMap), intent(in) :: idMap
   real(wp),intent(in) :: xyz(3,nat)

   type(TAtomList) :: atl

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
!  some debug xtb_printout
   if (set%verbose) then
      do idum = 1, narg
         write(env%unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('elements')
      call atl%new
      do idum = 1, narg
         ! get element by symbol
         if (idMap%has(argv(idum))) then
            call idMap%get(list, argv(idum))
            if (allocated(list)) then
               call atl%add(list)
            else
               call env%warning("Unknown element: '"//trim(argv(idum))//"'",source)
               cycle
            end if
         else
            ldum = getValue(env,trim(argv(idum)),iat)
            if (.not.ldum) cycle ! skip garbage input
            ! check for unreasonable input
            if (iat > 0) then
               call atl%add(at.eq.iat)
            else
               call env%warning("Unknown element: '"//trim(argv(idum))//"'",source)
               cycle
            endif
         endif
      enddo
      if (fixset%n > 0) call atl%add(fixset%atoms(:fixset%n))
      call atl%to_list(list)
      fixset%atoms = list
      fixset%n = size(list)
   case('atoms')
      call atl%new(val)
      if (atl%get_error()) then
         call env%warning('something is wrong in the fixing list',source)
         return
      endif
      if (fixset%n > 0) call atl%add(fixset%atoms(:fixset%n))
      call atl%to_list(list)
      fixset%atoms = list
      fixset%n = size(list)
   case('freeze')
      call atl%new(val)
      if (atl%get_error()) then
         call env%warning('something is wrong in the freezing list',source)
         return
      endif
      if (freezeset%n > 0) call atl%add(freezeset%atoms(:freezeset%n))
      call atl%to_list(list)
      freezeset%atoms = list
      freezeset%n = size(list)
   case('shake')
      allocate(list(nat*(nat+1)/2), source=0)
      if (mod(narg,2).ne.0) then
         call env%warning("could not read input for user defined shake!",source)
         return
      endif
         if (narg+shakeset%n > nat*(nat+1)/2) then
         call env%warning("too many SHAKE constraints!",source)
         return
      endif
      if (.not.set%shake_md) set%shake_md = .true.
      do idum = 1, narg
         if (getValue(env,trim(argv(idum)),iat)) then
            if (iat.gt.nat) then
               call env%warning('Attempted constrain atom not present in molecule.',source)
               cycle
            endif
            shakeset%n = shakeset%n+1
            shakeset%atoms(shakeset%n) = iat
         else
            call env%warning("Something went wrong in set_fix_ 'shake'.",source)
            return ! you screwed it, let's get out of here
         endif
      enddo
   end select

end subroutine set_fix

subroutine set_constr(env,key,val,nat,at,idMap,xyz)
   use xtb_mctc_constants
   use xtb_mctc_convert
   use xtb_type_atomlist
   use xtb_scanparam
   use xtb_splitparam
   implicit none
   character(len=*), parameter :: source = 'userdata_constr'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   type(TIdentityMap), intent(in) :: idMap
   real(wp),intent(in) :: xyz(3,nat)

   type(TAtomList) :: atl

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
   if (set%verbose) then
      do idum = 1, narg
         write(env%unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them

   case('elements')
      call atl%new
      do idum = 1, narg
         ! get element by symbol
         if (idMap%has(argv(idum))) then
            call idMap%get(list, argv(idum))
            if (allocated(list)) then
               call atl%add(list)
            else
               call env%warning("Unknown element: '"//trim(argv(idum))//"'",source)
               cycle
            end if
         else
            ldum = getValue(env,trim(argv(idum)),iat)
            if (.not.ldum) cycle ! skip garbage input
            ! check for unreasonable input
            if (iat > 0) then
               call atl%add(at.eq.iat)
            else
               call env%warning("Unknown element: '"//trim(argv(idum))//"'",source)
               cycle
            endif
         endif
      enddo
      if (potset%pos%n > 0) call atl%add(potset%pos%atoms(:potset%pos%n))
      call atl%to_list(list)
      potset%pos%atoms = list
      potset%pos%n = size(list)
   case('atoms')
      call atl%new(val)
      if (val.eq."all") then
         allocate(list(nat))
         do i=1,nat
            list(i)=i
         enddo 
      else
         if (atl%get_error()) then
            call env%warning('something is wrong in the fixing list',source)
            return
         endif
         if (potset%pos%n > 0) call atl%add(potset%pos%atoms(:potset%pos%n))
         call atl%to_list(list)
      endif
      potset%pos%atoms = list
      potset%pos%n = size(list)

   case('DISTANCE')
      if (narg.lt.3 .or. narg.gt.4) then
         call env%error('not enough arguments to constrain a distance',source)
         return
      endif
      ioffset = 2*potset%dist%n
      potset%dist%n = potset%dist%n+1
   !  part 1: get the constrained atoms
      do i = 1, 2
         if (getValue(env,trim(argv(i)),idum)) then
            potset%dist%atoms(ioffset+i) = idum
         else
            call env%warning("Something went wrong in set_constr_ 'distance'. (1)",source)
            potset%dist%n = potset%dist%n-1 ! remove invalid contrain
            return
         endif
      enddo
   !  part 2: get the distance between those atoms
      i = potset%dist%atoms(ioffset+1)
      j = potset%dist%atoms(ioffset+2)
      if (any([i, j] > nat .or. [i, j] < 1)) then
         call env%warning("Atomic index is out of bounds for distance constraint",source)
         potset%dist%n = potset%dist%n-1 ! remove invalid contrain
         return
      end if
      dist = sqrt(sum((xyz(:,i)-xyz(:,j))**2))
      if (trim(argv(3)).eq.'auto') then
         potset%dist%val(potset%dist%n) = dist
      else
         if (getValue(env,trim(argv(3)),ddum)) then
            potset%dist%val(potset%dist%n) = ddum * aatoau
         else
            call env%warning("Something went wrong in set_constr_ 'distance'. (2)",source)
            potset%dist%n = potset%dist%n-1 ! remove invalid contrain
            return
         endif
      endif
      if (narg.eq.4) then
         if (getValue(env,trim(argv(4)),idum)) then
            if (idum < 2 .or. mod(idum, 2).ne.0) then
               call env%warning("Invalid spring exponent given", source)
               potset%dist%n = potset%dist%n-1 ! remove invalid contrain
               return
            end if
            potset%dist%expo(potset%dist%n) = real(idum, wp)
         else
            call env%warning("Something went wrong in set_constr_ 'distance'. (3)",source)
            potset%dist%n = potset%dist%n-1 ! remove invalid contrain
            return
         end if
         write(env%unit,'("constraining bond",2(1x,i0),1x,"to",'//&
            '1x,f12.7,1x,"Å, actual value:",1x,f12.7,1x,"Å",1x,"with expo",1x,i0)') &
            i,j, potset%dist%val(potset%dist%n)*autoaa, dist*autoaa, &
            nint(potset%dist%expo(potset%dist%n))
      else
         write(env%unit,'("constraining bond",2(1x,i0),1x,"to",'//&
            '1x,f12.7,1x,"Å, actual value:",1x,f12.7,1x,"Å")') &
            i,j, potset%dist%val(potset%dist%n)*autoaa, dist*autoaa
      end if

   case('ANGLE')
      if (narg.ne.4) then
         call env%error('not enough arguments to constrain an angle',source)
         return
      endif
      ioffset = 3*potset%angle%n
      potset%angle%n = potset%angle%n+1
   !  part 1: get the constrained atoms
      do i = 1, narg-1
         if (getValue(env,trim(argv(i)),idum)) then
            potset%angle%atoms(ioffset+i) = idum
         else
            call env%warning("Something went wrong in set_constr_ 'angle'. (1)",source)
            potset%angle%n = potset%angle%n-1 ! remove invalid contrain
            return
         endif
      enddo
   !  part 2: get the angle between the constrained atoms
      i = potset%angle%atoms(ioffset+1)
      j = potset%angle%atoms(ioffset+2)
      k = potset%angle%atoms(ioffset+3)
      if (any([i, j, k] > nat .or. [i, j, k] < 1)) then
         call env%warning("Atomic index is out of bounds for angle constraint",source)
         potset%angle%n = potset%angle%n-1 ! remove invalid contrain
         return
      end if
      call bangl(xyz,i,j,k,phi)
      if (trim(argv(narg)).eq.'auto') then
         potset%angle%val(potset%angle%n) = phi
      else
         if (getValue(env,trim(argv(narg)),ddum)) then
            potset%angle%val(potset%angle%n) = pi/180.0_wp * ddum
         else
            call env%warning("Something went wrong in set_constr_ 'angle'. (2)",source)
            potset%angle%n = potset%angle%n-1 ! remove invalid contrain
            return
         endif
      endif

   case('DIHEDRAL')
      if (narg.ne.5) then
         call env%error('not enough arguments to constrain a dihedral',source)
         return
      endif
      ioffset = 4*potset%dihedral%n
      potset%dihedral%n = potset%dihedral%n+1
      if (nconstr.gt.maxconstr) then ! double check this
         call env%error('This should never happen! Let somebody check set_constr',source)
         return
      end if
   !  part 1: get the constrained atoms
      do i = 1, narg-1
         if (getValue(env,trim(argv(i)),idum)) then
            potset%dihedral%atoms(ioffset+i) = idum
         else
            call env%warning("Something went wrong in set_constr_ 'dihedral'. (1)",source)
            potset%dihedral%n = potset%dihedral%n-1 ! remove invalid contrain
            return
         endif
      enddo
   !  part 2: get the angle between the constrained atoms
      i = potset%dihedral%atoms(ioffset+1)
      j = potset%dihedral%atoms(ioffset+2)
      k = potset%dihedral%atoms(ioffset+3)
      l = potset%dihedral%atoms(ioffset+4)
      if (any([i, j, k, l] > nat .or. [i, j, k, l] < 1)) then
         call env%warning("Atomic index is out of bounds for dihedral constraint",source)
         potset%dihedral%n = potset%dihedral%n-1 ! remove invalid contrain
         return
      end if
      phi=valijkl(nat,xyz,i,j,k,l)
      if (trim(argv(narg)).eq.'auto') then
         potset%dihedral%val(potset%dihedral%n) = phi
      else
         if (getValue(env,trim(argv(narg)),ddum)) then
            potset%dihedral%val(potset%dihedral%n) = pi/180.0_wp * ddum
         else
            call env%warning("Something went wrong in set_constr_ 'dihedral'. (2)",source)
            potset%dihedral%n = potset%dihedral%n-1 ! remove invalid contrain
            return
         endif
      endif

   case('distance')
      if (narg.ne.3) then
         call env%error('not enough arguments to constrain a distance',source)
         return
      endif
      nconstr = nconstr+1
      if (nconstr.gt.maxconstr) then ! double check this
         call env%error('This should never happen! Let somebody check set_constr',source)
         return
      end if
   !  part 1: get the constrained atoms
      do i = 1, narg-1
         if (getValue(env,trim(argv(i)),idum)) then
            atconstr(i,nconstr) = idum
         else
            call env%warning("Something went wrong in set_constr_ 'distance'. (1)",source)
            nconstr = nconstr-1 ! remove invalid contrain
            return
         endif
      enddo
   !  part 2: get the distance between those atoms
      i = atconstr(1,nconstr)
      j = atconstr(2,nconstr)
      if (any([i, j] > nat .or. [i, j] < 1)) then
         call env%warning("Atomic index is out of bounds for distance constraint",source)
         nconstr = nconstr-1 ! remove invalid contrain
         return
      end if
      dist = sqrt(sum((xyz(:,i)-xyz(:,j))**2))
      if (trim(argv(narg)).eq.'auto') then
         valconstr(nconstr) = dist
      else
         if (getValue(env,trim(argv(narg)),ddum)) then
            valconstr(nconstr) = ddum * aatoau
         else
            call env%warning("Something went wrong in set_constr_ 'distance'. (2)",source)
            nconstr = nconstr-1 ! remove invalid contrain
            return
         endif
      endif
      write(env%unit,'("constraining bond",2(1x,i0),1x,"to",'//&
         '1x,f12.7,1x,"Å, actual value:",1x,f12.7,1x,"Å")') &
         i,j, valconstr(nconstr)*autoaa, dist*autoaa

   case('angle')
      if (narg.ne.4) then
         call env%error('not enough arguments to constrain an angle',source)
         return
      endif
      nconstr = nconstr+1
      if (nconstr.gt.maxconstr) then ! double check this
         call env%error('This should never happen! Let somebody check set_constr',source)
         return
      end if
   !  part 1: get the constrained atoms
      do i = 1, narg-1
         if (getValue(env,trim(argv(i)),idum)) then
            atconstr(i,nconstr) = idum
         else
            call env%warning("Something went wrong in set_constr_ 'angle'. (1)",source)
            nconstr = nconstr-1 ! remove invalid contrain
            return
         endif
      enddo
   !  part 2: get the angle between the constrained atoms
      i = atconstr(1,nconstr)
      j = atconstr(2,nconstr)
      k = atconstr(3,nconstr)
      if (any([i, j, k] > nat .or. [i, j, k] < 1)) then
         call env%warning("Atomic index is out of bounds for angle constraint",source)
         nconstr = nconstr-1 ! remove invalid contrain
         return
      end if
      call bangl(xyz,i,j,k,phi)
      if (trim(argv(narg)).eq.'auto') then
         valconstr(nconstr) = phi
      else
         if (getValue(env,trim(argv(narg)),ddum)) then
            valconstr(nconstr) = pi/180.0_wp * ddum
         else
            call env%warning("Something went wrong in set_constr_ 'angle'. (2)",source)
            nconstr = nconstr-1 ! remove invalid contrain
            return
         endif
      endif
      write(env%unit,'("constraining angle",3(1x,i0),1x,"to",'//&
         '1x,f12.7,"°, actual value:",1x,f12.7,"°")') &
         i,j,k,180.0_wp/pi * valconstr(nconstr),180.0_wp/pi * phi

   case('dihedral')
      if (narg.ne.5) then
         call env%error('not enough arguments to constrain a dihedral',source)
         return
      endif
      nconstr = nconstr+1
      if (nconstr.gt.maxconstr) then ! double check this
         call env%error('This should never happen! Let somebody check set_constr',source)
         return
      end if
   !  part 1: get the constrained atoms
      do i = 1, narg-1
         if (getValue(env,trim(argv(i)),idum)) then
            atconstr(i,nconstr) = idum
         else
            call env%warning("Something went wrong in set_constr_ 'dihedral'. (1)",source)
            nconstr = nconstr-1 ! remove invalid contrain
            return
         endif
      enddo
   !  part 2: get the angle between the constrained atoms
      i = atconstr(1,nconstr)
      j = atconstr(2,nconstr)
      k = atconstr(3,nconstr)
      l = atconstr(4,nconstr)
      if (any([i, j, k, l] > nat .or. [i, j, k, l] < 1)) then
         call env%warning("Atomic index is out of bounds for dihedral constraint",source)
         nconstr = nconstr-1 ! remove invalid contrain
         return
      end if
      phi=valijkl(nat,xyz,i,j,k,l)
      if (trim(argv(narg)).eq.'auto') then
         valconstr(nconstr) = phi
      else
         if (getValue(env,trim(argv(narg)),ddum)) then
            valconstr(nconstr) = pi/180.0_wp * ddum
         else
            call env%warning("Something went wrong in set_constr_ 'dihedral'. (2)",source)
            nconstr = nconstr-1 ! remove invalid contrain
            return
         endif
      endif
      write(env%unit,'("constraining angle",4(1x,i0),1x,"to",'//&
         '1x,f12.7,"°, actual value:",1x,f12.7,"°")') &
         i,j,k,l,180.0_wp/pi * valconstr(nconstr),180.0_wp/pi * phi

   case('center')
      ! copied from molbld.f without modification (except for error handling)
      if (narg.ne.2) then
         call env%error('not enough argument to constrain center',source)
         return
      endif
      nconstr = nconstr+1
      atconstr(1,nconstr) = -2
      if (getValue(env,trim(argv(2)),ddum)) then
         valconstr(nconstr) = ddum
      else
         call env%warning("Something went wrong in set_constr_ 'center'. (1)",source)
         nconstr = nconstr-1
         return
      endif
      if (getValue(env,trim(argv(1)),idum)) then
         iatf1 = idum
      else
         call env%warning("Something went wrong in set_constr_ 'center'. (2)",source)
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
         call env%error('not enough argument to constrain cma',source)
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
         write(env%unit,'("constraining fragment CMA to initial R(Ang.)=",f8.3)')&
            valconstr(nconstr)*autoaa
      else
         if (getValue(env,trim(argv(1)),ddum)) then
            valconstr(nconstr) = ddum*aatoau
         else
            call env%warning("Something went wrong in set_constr_ 'cma'.",source)
            nconstr = nconstr-1
            return
         endif
      endif

!   case('x','y','z')
!      idum = index('xyz',key)
   case('z')
      ! copied from molbld.f without modification (except for error handling)
      if (narg.ne.1) then
         call env%error('not enough argument to constrain z coordinate',source)
         return
      endif
      zconstr = 1
      nconstr = nconstr+1
      atconstr(1,nconstr) = -1
      if (getValue(env,trim(argv(1)),ddum)) then
         valconstr(nconstr) = ddum*aatoau
      else
         call env%warning("Something went wrong in set_constr_ 'z'.",source)
         nconstr = nconstr-1
         return
      endif
      if (sum(abs(xyz(3,1:iatf1))).gt.1.d-3) then
         call env%warning('z-coordinates of fragment 1 must be 0!',source)
         nconstr = nconstr-1
         return
      endif
      write(env%unit,'("constraining fragment 2 to Z=0 plane at (Ang.)=",f8.3)') &
         valconstr(nconstr)*autoaa

   end select
end subroutine set_constr

!! --------------------------------------------------------------[SAW1809]-
!  this is the new version of the scan routine exploiting all features
subroutine set_scan(env,key,val,nat,at,idMap,xyz)
   use xtb_scanparam
   use xtb_mctc_convert, only : aatoau
   use xtb_mctc_constants, only : pi
   implicit none
   character(len=*), parameter :: source = 'userdata_scan'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   type(TIdentityMap), intent(in) :: idMap
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
      call set_constr(env,key,temp,nat,at,idMap,xyz) ! generate a new constraint
      scan_list(nscan)%iconstr = nconstr ! new generated
      if (idum.eq.nconstr) then
         call env%error('Failed to generate constraint',source)
         return
      end if
      temp = val(1+ie:)
   else
      if (getValue(env,key,idum)) then
         if (idum.gt.nconstr) then
            call env%error('Constraint '''//key//''' is not defined',source)
            return
         endif
         scan_list(nscan)%iconstr = idum
      else
         call env%error('Constraint '''//key//''' is invalid in this context',source)
         return
      endif
      temp = val
   endif

   call parse(temp,comma,argv,narg)
   if (set%verbose) then
      do idum = 1, narg
         write(env%unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
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

   if (getValue(env,trim(argv(1)),ddum)) then
      if (idum.le.1) then
         start_value = ddum * aatoau
      else
         start_value = ddum * pi/180.0_wp
      endif
   else
      call env%error('Invalid start value for scan',source)
      return
   endif
   if (getValue(env,trim(argv(2)),ddum)) then
      if (idum.le.1) then
         end_value = ddum * aatoau
      else
         end_value = ddum * pi/180.0_wp
      endif
   else
      call env%error('Invalid end value for scan',source)
      return
   endif
   if (getValue(env,trim(argv(3)),idum)) then
      scan_list(nscan)%nscan = idum
      allocate( scan_list(nscan)%valscan(idum), source = 0.0_wp )
   else
      call env%error('Invalid step number for scan',source)
      return
   endif

   do i = 1, scan_list(nscan)%nscan
      scan_list(nscan)%valscan(i) = start_value + (end_value-start_value) &
         &                       * real(i-1,wp)/real(scan_list(nscan)%nscan-1,wp)
      if (set%verbose) write(env%unit,'(i5,1x,f12.8)') i,scan_list(nscan)%valscan(i)
   enddo

end subroutine set_scan

subroutine set_wall(env,key,val,nat,at,idMap,xyz)
   use xtb_mctc_convert, only : autoaa
   use xtb_sphereparam
   implicit none
   character(len=*), parameter :: source = 'userdata_wall'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   type(TIdentityMap), intent(in) :: idMap
   real(wp),intent(in) :: xyz(3,nat)

   integer  :: idum,i
   real(wp) :: ddum,darray(3),min_z,max_z
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
   if (set%verbose) then
      do idum = 1, narg
         write(env%unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('sphere')
      if (narg.lt.2) then
         call env%error("Not enough arguments to set up a spherical wall",source)
         return
      endif
   !  part 1: get the sphere radius
      if (trim(argv(1)).eq.'auto') then
         call get_sphere_radius(nat,at,xyz,center,radius,do_trafo=.false.)
      else
         if (getValue(env,trim(argv(1)),ddum)) then
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
            if (getListValue(env,trim(argv(iarg)),tlist,ntlist)) then
               if (nlist+ntlist.gt.nat) then
                  call env%error("Too many atoms in list for spherical wall.",source)
                  return ! something went wrong
               endif
               if (maxval(tlist(:ntlist)).gt.nat) then
                  call env%error("Attempted to wall in a non-existing atom",source)
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
      write(env%unit,'("spherical wallpotential with radius",'//&
         '1x,f12.7,1x,"Å")') radius*autoaa

   case('ellipsoid')
      if (narg.lt.4) then
         call env%error("Not enough arguments to set up an ellipsoidal wall",source)
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
            if (getValue(env,trim(argv(iaxis)),ddum)) then
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
            if (getListValue(env,trim(argv(iarg)),tlist,ntlist)) then
               if (nlist+ntlist.gt.nat) then
                  call env%error("Too many atoms in list for spherical wall.",source)
                  return ! something went wrong
               endif
               if (maxval(tlist(:ntlist)).gt.nat) then
                  call env%error("Attempted to wall in a non-existing atom",source)
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
      write(env%unit,'("ellipsoidal wallpotenial with radii",'//&
         '3(1x,f12.7,1x,"Å"))') radii*autoaa

   case('sandwich')
      if (narg.lt.2) then
         call env%error("Not enough arguments to set up sandwich walls",source)
         return
      endif
   !  part 1: get the sandwich distance
      wpot%sandwich = .true.
      set%do_cma_trafo = .true.
      center = 0.0_wp
!      number_walls=1
      call get_sphere_radius(nat,at,xyz,center,radius,do_trafo=.true.)
      if (trim(argv(1)).eq.'auto') then
            radius=(maxval(xyz(3,:))-minval(xyz(3,:)))/2.0_wp ! need to set $cma in xcontrol, not done automatically
            wpot(1)%radius=radius
      else
            if (getValue(env,trim(argv(1)),ddum)) then
                radius = ddum !in Bohr!!!
                wpot(1)%radius=radius
            else
                call env%error("Undefined arguments for sandwich: ... in your xcontrol file!",source)
                return ! something went wrong
            endif
      endif
   
   !  part 2: get atoms 
      if (trim(argv(2)).eq.'all') then
         call set_sphere_radius(radius,center)
      else
         do iarg = 2, narg
            if (getListValue(env,trim(argv(iarg)),tlist,ntlist)) then
               if (nlist+ntlist.gt.nat) then
                  call env%error("Too many atoms in list for spherical wall.",source)
                  return ! something went wrong
               endif
               if (maxval(tlist(:ntlist)).gt.nat) then
                  call env%error("Attempted to wall in a non-existing atom.",source)
                  cycle ! skip crappy input
               endif
               list(nlist+1:nlist+ntlist) = tlist
               nlist = nlist + ntlist

               !get auto sandwich distance for list of atoms
               max_z = 0.0_wp
               min_z = 0.0_wp
               do i = 1, nat
                  if (any(list == i)) then
                     max_z = max(max_z, xyz(3,i))
                     min_z = min(min_z, xyz(3,i))
                  end if
               end do
               radius=(max_z - min_z)/2.0_wp
               wpot(1)%radius=radius

            else
               ! warning already generated by get_list_value
               return ! something went wrong
            endif
         enddo
         call set_sphere_radius(radius,center,nlist,list)
      endif

      write(env%unit,'("sandwich wallpotential with radius in A (diameter=2*radius+2*4A safety buffer) ",'//&
         '1x,f12.7,1x,"Å")') radius*autoaa

   end select

end subroutine set_wall

subroutine set_split(env,key,val,nat,at,idMap,xyz)
   use xtb_splitparam
   implicit none
   character(len=*), parameter :: source = 'userdata_split'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   type(TIdentityMap), intent(in) :: idMap
   real(wp),intent(in) :: xyz(3,nat)

   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   integer  :: list(nat),nlist
   integer  :: i,j

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   call parse(val,comma,argv,narg)
   if (set%verbose) then
      do idum = 1, narg
         write(env%unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('fragment1')
      do i = 1, narg
         if (getListValue(env,trim(argv(i)),list,nlist)) then
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
         if (getListValue(env,trim(argv(i)),list,nlist)) then
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
      if (getValue(env,trim(argv(1)),idum)) then
         if (idum.gt.nat) then
            call env%warning("rejecting fragment number greater than number of atoms",source)
            return ! doesn't really make sense, sorry, your problem
         endif
         do i = 2, narg
            if (getListValue(env,trim(argv(i)),list,nlist)) then
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

subroutine set_hess(env,key,val,nat,at,idMap,xyz)
   use xtb_type_atomlist, only : TAtomList
   use xtb_splitparam
   implicit none
   character(len=*), parameter :: source = 'userdata_hess'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   type(TIdentityMap), intent(in) :: idMap
   real(wp),intent(in) :: xyz(3,nat)

   type(TAtomList) :: atl
   integer, allocatable :: list(:)
   real(wp) :: ddum
   integer  :: i,j,idum,iat,narg
   character(len=p_str_length),dimension(p_arg_length) :: argv
   character(len=256) :: warningstring

   call parse(val,comma,argv,narg)
   if (set%verbose) then
      do idum = 1, narg
         write(env%unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('element mass')
      if (mod(narg,2).ne.0) then
         call env%warning("Something went wrong in set_mass_ 'element'.",source)
      endif
      do i = 1, narg, 2
         j = i+1
         if (getValue(env,trim(argv(j)),ddum)) then
            if (idMap%has(argv(i))) then
               call idMap%set(atmass, argv(i), ddum)
               write(env%unit,'(a,a,a,1x,g0)') &
                  "mass of elements '",trim(argv(i)),"' changed to", ddum
            else if (getValue(env,trim(argv(i)),idum)) then
               where(at.eq.idum) atmass = ddum
               write(env%unit,'(a,i0,a,1x,g0)') &
                  "mass of elements with Z=",idum," changed to", ddum
            end if
         end if
      enddo
   case('modify mass','isotope')
      if (mod(narg,2).ne.0) then
         call env%warning("Something went wrong in set_mass_ 'modify'.",source)
      endif
      do i = 1, narg, 2
         j = i+1
         if (getValue(env,trim(argv(j)),ddum)) then
            call atl%new(argv(i))
            if (atl%get_error()) then
               call env%warning('something is wrong in the mass list',source)
               cycle
            endif
            call atl%to_list(list)
            do idum = 1, size(list)
               iat = list(idum)
               if (iat.gt.nat) then
                  write(warningstring, '(a, i0, a)') 'Attempted setting atom mass for atom ', &
                     & iat, ' that is not present in system.'
                  call env%warning(trim(warningstring), source)
                  cycle
               endif
               atmass(iat) = ddum
               write(env%unit,'(a,1x,i0,1x,a,1x,g0)') &
                  & 'mass of atom ',iat,' changed to',atmass(iat)
            enddo
            call atl%destroy()
         endif
      enddo
   case('scale mass')
      if (mod(narg,2).ne.0) then
         call env%warning("Something went wrong in set_mass_ 'scale'.",source)
      endif
      do i = 1, narg, 2
         j = i+1
         if (getValue(env,trim(argv(j)),ddum)) then
            call atl%new(argv(i))
            if (atl%get_error()) then
               call env%warning('something is wrong in the mass list',source)
               cycle
            endif
            call atl%to_list(list)
            do idum = 1, size(list)
               iat = list(idum)
               if (iat.gt.nat) then
                  write(warningstring, '(a, i0, a)') 'Attempted setting atom mass for atom ', &
                  & iat, ' that is not present in system.'
                  call env%warning(trim(warningstring), source)
                  cycle
               endif
               atmass(iat) = atmass(iat)*ddum
               write(env%unit,'(a,1x,i0,1x,a,1x,g0)') &
                  'mass of atom ',iat,' changed to',atmass(iat)
            enddo
            call atl%destroy()
         endif
      enddo
   end select

end subroutine set_hess

subroutine set_reactor(env,key,val,nat,at,idMap,xyz)
   use xtb_type_atomlist
   use xtb_setparam
   implicit none
   character(len=*), parameter :: source = 'userdata_reactor'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   type(TIdentityMap), intent(in) :: idMap
   real(wp),intent(in) :: xyz(3,nat)

   type(TAtomList) :: atl

   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   integer  :: nlist
   integer, allocatable :: list(:)
   integer  :: i,j

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   call parse(val,comma,argv,narg)
   if (set%verbose) then
      do idum = 1, narg
         write(env%unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('atoms')
      call atl%new(val)
      if (atl%get_error()) then
         call env%warning('something is wrong in the reactor atom list',source)
         return
      endif
      if (set%reactset%nat > 0) call atl%add(set%reactset%atoms(:set%reactset%nat))
      call atl%to_list(list)
      set%reactset%atoms = list
      set%reactset%nat = size(list)
   end select

end subroutine set_reactor

subroutine set_path(env,key,val,nat,at,idMap,xyz)
   use xtb_type_atomlist
   use xtb_setparam
   implicit none
   character(len=*), parameter :: source = 'userdata_path'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   type(TIdentityMap), intent(in) :: idMap
   real(wp),intent(in) :: xyz(3,nat)

   type(TAtomList) :: atl

   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   integer  :: nlist
   integer, allocatable :: list(:)
   integer  :: i,j

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   call parse(val,comma,argv,narg)
   if (set%verbose) then
      do idum = 1, narg
         write(env%unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('atoms')
      call atl%new(val)
      if (atl%get_error()) then
         call env%warning('something is wrong in the bias path atom list',source)
         return
      endif
      if (set%pathset%nat > 0) call atl%add(set%pathset%atoms(:set%pathset%nat))
      call atl%to_list(list)
      set%pathset%atoms = list
      set%pathset%nat = size(list)
   end select

end subroutine set_path

subroutine set_metadyn(env,key,val,nat,at,idMap,xyz)
   use xtb_type_atomlist
   use xtb_fixparam
   implicit none
   character(len=*), parameter :: source = 'userdata_metadyn'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   type(TIdentityMap), intent(in) :: idMap
   real(wp),intent(in) :: xyz(3,nat)

   type(TAtomList) :: atl

   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   integer  :: nlist
   integer, allocatable :: list(:)
   integer  :: i,j,iat

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   call parse(val,comma,argv,narg)
   if (set%verbose) then
      do idum = 1, narg
         write(env%unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('bias atoms')
      call atl%new(val)
      if (atl%get_error()) then
         call env%warning('something is wrong in the metadynamic atom list',source)
         return
      endif
      if (rmsdset%nat > 0) call atl%add(rmsdset%atoms(:rmsdset%nat))
      call atl%to_list(list)
      rmsdset%atoms = list
      rmsdset%nat = size(list)
   case('bias elements')
      call atl%new
      do idum = 1, narg
         ! get element by symbol
         if (idMap%has(argv(idum))) then
            call idMap%get(list, argv(idum))
            if (allocated(list)) then
               call atl%add(list)
            else
               call env%warning("Unknown element: '"//trim(argv(idum))//"'",source)
               cycle
            end if
         else
            ldum = getValue(env,trim(argv(idum)),iat)
            if (.not.ldum) cycle ! skip garbage input
            ! check for unreasonable input
            if (iat > 0) then
               call atl%add(at.eq.iat)
            else
               call env%warning("Unknown element: '"//trim(argv(idum))//"'",source)
               cycle
            endif
         endif
      enddo
      if (rmsdset%nat > 0) call atl%add(rmsdset%atoms(:rmsdset%nat))
      call atl%to_list(list)
      rmsdset%atoms = list
      rmsdset%nat = size(list)
   case('atoms')
      call atl%new(val)
      if (atl%get_error()) then
         call env%warning('something is wrong in the metadynamic atom list',source)
         return
      endif
      if (metaset%nat > 0) call atl%add(metaset%atoms(:metaset%nat))
      call atl%to_list(list)
      metaset%atoms = list
      metaset%nat = size(list)
   case('modify factor')
      if (mod(narg,2).ne.0) then
         call env%warning("Something went wrong in set_metadyn_ 'modify'.",source)
      endif
      do i = 1, narg, 2
         j = i+1
         if (getValue(env,trim(argv(i)),idum).and.&
             getValue(env,trim(argv(j)),ddum)) then
            if (idum.gt.metaset%maxsave) then
               call env%warning('Attempted using factor not present in system.',source)
               cycle
            endif
            metaset%factor(idum) = ddum
         endif
      enddo
   case('scale factor')
      if (mod(narg,2).ne.0) then
         call env%warning("Something went wrong in set_metadyn_ 'scale'.",source)
      endif
      do i = 1, narg, 2
         j = i+1
         if (getValue(env,trim(argv(i)),idum).and.&
             getValue(env,trim(argv(j)),ddum)) then
            if (idum.gt.metaset%maxsave) then
               call env%warning('Attempted scaling factor not present in system.',source)
               cycle
            endif
            metaset%factor(idum) = metaset%factor(idum)*ddum
         endif
      enddo
   end select

end subroutine set_metadyn

subroutine set_freeze(env,key,val,nat,at,idMap,xyz)
   implicit none
   character(len=*), parameter :: source = 'userdata_freeze'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   type(TIdentityMap), intent(in) :: idMap
   real(wp),intent(in) :: xyz(3,nat)

   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   integer  :: list(nat),nlist
   integer  :: i,j

   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   call parse(val,comma,argv,narg)
   if (set%verbose) then
      do idum = 1, narg
         write(env%unit,'("-->",1x,i0,":",1x,a)') idum, trim(argv(idum))
      enddo
   endif
   select case(key)
   case default ! ignore, don't even think about raising them
   case('hessf')
   case('hessa')
   end select

end subroutine set_freeze

subroutine set_legacy(env,key,val,nat,at,idMap,xyz)
   implicit none
   character(len=*), parameter :: source = 'userdata_legacy'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   type(TIdentityMap), intent(in) :: idMap
   real(wp),intent(in) :: xyz(3,nat)
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   select case(key)
   case default ! complaining about unknown keywords should already be done
      continue  ! so we do nothing here
   case('hessf');       call set_fix(env,'freeze',val,nat,at,idMap,xyz)
!   case('hessa');       call set_frozh('hessa',val)
   case('fragment1'); call set_split(env,'fragment1',val,nat,at,idMap,xyz)
   case('fragment2'); call set_split(env,'fragment1',val,nat,at,idMap,xyz)
   case('constrxyz'); call set_fix(env,'atoms',val,nat,at,idMap,xyz)
!   case('constrainel')
!   case('constrain')
!   case('scan')
   case('ellips'); call set_wall(env,'ellipsoid',val,nat,at,idMap,xyz)
   case('sphere'); call set_wall(env,'sphere',val,nat,at,idMap,xyz)
   case('sandwich'); call set_wall(env,'sandwich',val,nat,at,idMap,xyz)
   case('fix'); call set_fix(env,'atoms',val,nat,at,idMap,xyz)
   case('atomlist+'); call set_metadyn(env,'atoms',val,nat,at,idMap,xyz)
   end select

end subroutine set_legacy

end module xtb_constrain_param
