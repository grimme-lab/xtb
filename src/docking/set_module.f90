! This file is part of xtb.
!
! Copyright (C) 2022 Christoph Plett
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
!  * WELCOME TO THE   C O N S T R A I N S   &   S C A N S   MODULE OF
!  THE DOCKING ALGORITHM *
!! ------------------------------------------------------------------------
module xtb_docking_set_module
   use xtb_mctc_accuracy, only: wp
   use xtb_readin, only: mirror_line, getValue
   use xtb_docking_param
   use xtb_type_environment, only: TEnvironment
   use xtb_type_identitymap, only : TIdentityMap
   use xtb_type_atomlist
   use xtb_mctc_strings, only : parse
   use xtb_setmod

   implicit none

   private :: wp, mirror_line, getValue

   character, private, parameter :: flag = '$'
   character, private, parameter :: colon = ':'
   character, private, parameter :: space = ' '
   character, private, parameter :: equal = '='
   character, private, parameter :: hash = '#'
   character, private, parameter :: dot = '.'
   character, private, parameter :: comma = ','
   character(len=*), private, parameter :: flag_end = flag//'end'

!  Using allocatable arrays of dynamic length strings is only possible
!  with a lot of hacks, so we use good'ol fixed size stack arrays.
!  Let's choose something different from 42 that is not dividable by 10... ;)
!  Happy debugging!
   integer, private, parameter :: p_str_length = 48
   integer, private, parameter :: p_arg_length = 24

   public

      !abstract interface
      !   subroutine handlerInterface(env, key, val)
      !      import TEnvironment
      !      type(TEnvironment), intent(inout) :: env
      !      character(len=*), intent(in) :: key
      !      character(len=*), intent(in) :: val
      !   end subroutine handlerInterface
      !end interface

      abstract interface
         subroutine handlerInterface2(env,key,val,nat,at,idMap,xyz)
            import :: wp, TEnvironment, TIdentityMap
            type(TEnvironment), intent(inout) :: env
            character(len=*),intent(in) :: key
            character(len=*),intent(in) :: val
            integer, intent(in) :: nat
            type(TIdentityMap), intent(in) :: idMap
            integer, intent(in) :: at(nat)
            real(wp),intent(in) :: xyz(3,nat)
         end subroutine handlerInterface2
      end interface

contains

   subroutine rdblock_docking(env, handler, line, id, copy, err, ncount)
      !use xtb_setmod, only : handlerInterface
      character(len=*), parameter :: source = 'set_rdblock'
      type(TEnvironment), intent(inout) :: env
      procedure(handlerInterface) :: handler
      integer, intent(in) :: id
      integer, intent(in) :: copy
      integer, intent(out) :: err
      integer, intent(out) :: ncount
      character(len=:), allocatable, intent(out) :: line
      character(len=:), allocatable :: key
      character(len=:), allocatable :: val
      integer :: ie
      logical :: exitRun
      ncount = 0
      do
         call mirror_line(id, copy, line, err)
         if (is_iostat_end(err)) return
         if (index(line, flag) .ne. 0) return
         ! find the equal sign
         ie = index(line, equal)
         if (line .eq. '') cycle ! skip empty lines
         ncount = ncount + 1   ! but count non-empty lines first
         if (ie .eq. 0) then! cycle
            call set_logicals(env, line)
         else
            key = trim(line(:ie - 1))
            val = trim(adjustl(line(ie + 1:)))
            call handler(env, key, val)
            call env%check(exitRun)
            if (exitRun) then
               call env%error("handler could not process input", source)
               return
            end if
         end if
      end do

   end subroutine rdblock_docking

   subroutine rdblock_docking2(env,handler,line,id,nat,at,idMap,xyz,err)
      !use xtb_constrain_param, only : handlerInterface
      use xtb_readin, only : getline => strip_line,getValue,getListValue
      character(len=*), parameter :: source = 'userdata_rdblock'
      type(TEnvironment), intent(inout) :: env
      integer,intent(in) :: id
      procedure(handlerInterface2) :: handler
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
         if (debug) write(env%unit,'("->",1x,a)') line
   
         ! find the first colon
         ie = index(line,colon)
         if ((line.eq.'')) cycle
         if (ie .eq. 0) then! cycle
            call set_logicals(env, line)
         else
            key = trim(line(:ie-1))
            val = trim(adjustl(line(ie+1:)))
            call handler(env,key,val,nat,at,idMap,xyz)
            call env%check(exitRun)
            if (exitRun) then
               call env%error("handler could not process input", source)
               return
            end if
         end if
      enddo

   end subroutine rdblock_docking2

   subroutine set_docking(env, key, val)
      character(len=*), parameter :: source = 'set_docking'
      type(TEnvironment), intent(inout) :: env
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: val
      integer  :: err
      integer  :: idum
      real(wp) :: ddum
      logical  :: ldum
      logical, save :: set1 = .true.
      logical, save :: set2 = .true.
      logical, save :: set3 = .true.
      logical, save :: set4 = .true.
      logical, save :: set5 = .true.
      logical, save :: set6 = .true.
      select case (key)
      case default ! do nothing
        call env%warning("the key '"//key//"' is not recognized by scc", source)
      case ('stepr')
         if (getValue(env, val, ddum) .and. set1) stepr = ddum
         if (stepr < 1.0) then
            stepr = 1.0
            call env%warning("Too small step for radial grid, taking 1.0.")
         end if
         set1 = .false.
      case ('stepa')
         if (getValue(env, val, ddum) .and. set2) stepa = ddum
         if (stepa < 1.0) then
            stepa = 1.0
         call env%warning("Too small step for angular grid, taking 1.0. This is really expensive!")
         end if
         set2 = .false.
      case ('nfinal')
         if (getValue(env, val, idum) .and. set3) n_opt = idum
         set3 = .false.
      case ('maxgen')
         if (getValue(env, val, idum) .and. set4) maxgen = idum
         set4 = .false.
      case ('maxparent')
         if (getValue(env, val, idum) .and. set5) maxparent = idum
         set5 = .false.
      case ('nstack')
         if (getValue(env, val, idum) .and. set6) mxcma = idum
         set6 = .false.
      end select
   end subroutine set_docking

   subroutine write_set_directed(ictrl)
      use xtb_type_atomlist
      integer, intent(in) :: ictrl
      type(TAtomList) :: atl
      character(len=:), allocatable :: string
      integer :: i
      if (directedset%n.eq.0) return
   
      write(ictrl,'(a,"directed")') flag
      if (directedset%n > 0) then
         call atl%new(directedset%atoms(:directedset%n))
         call atl%to_string(string)
         write(ictrl,'(3x,"atoms:",1x,a)') string
      endif
   
   end subroutine write_set_directed

   subroutine set_directed(env,key,val,nat,at,idMap,xyz)
      character(len=*), parameter :: source = 'userdata_directed'
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
      if (debug) then
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
         if (directedset%n > 0) call atl%add(directedset%atoms(:directedset%n))
         call atl%to_list(list)
         directedset%atoms = list
         directedset%n = size(list)
      case('atoms')
         call atl%new(val)
         if (atl%get_error()) then
            call env%warning('something is wrong in the fixing list',source)
            return
         endif
         if (directedset%n > 0) call atl%add(directedset%atoms(:directedset%n))
         call atl%to_list(list)
         directedset%atoms = list
         directedset%n = size(list)
      end select
      call write_set_directed(env%unit)
   
   end subroutine set_directed

   subroutine set_logicals(env, key)
      character(len=*), parameter :: source = 'set_logicals'
      type(TEnvironment), intent(inout) :: env
      character(len=*), intent(in) :: key
      integer  :: err
      integer  :: idum
      real(wp) :: ddum
      logical  :: ldum
      logical, save :: set1 = .true.
      logical, save :: set2 = .true.
      logical, save :: set3 = .true.
      logical, save :: set4 = .true.
      logical, save :: set5 = .true.
      logical, save :: set6 = .true.
      logical, save :: set7 = .true.
      logical, save :: set8 = .true.
      logical, save :: set9 = .true.
      logical, save :: set10 = .true.
      logical, save :: set11 = .true.
      logical, save :: set12 = .true.
      logical, save :: set13 = .true.
      logical, save :: set14 = .true.
      logical, save :: set15 = .true.
      select case (key)
      case default ! do nothing
        call env%warning("the key '"//key//"' is not recognized by scc", source)
      case ('nopocket')
         if (set1) pocket_grid = .false.
         set1 = .false.
      case ('pocket')
         if (set2) pocket_grid = .true.
         set2 = .false.
      case ('nostack')
         if (set3) stack_grid = .false.
         set3 = .false.
      case ('stack')
         if (set4) stack_grid = .true.
         set4 = .false.
      case ('noangular')
         if (set5) angular_grid = .false.
         set5 = .false.
      case ('angular')
         if (set6) angular_grid = .true.
         set6 = .false.
      case ('atm')
         if (set7) fulle = .true.
         set7 = .false.
      case ('fast')
         if (set8) then
            maxparent = 30      ! # of parents in gene pool 100
            maxgen = 7          ! # of generations 10
            mxcma = 250         ! R points in CMA search 1000
            stepr = 4.0         ! R grid step in Bohr 2.5
            stepa = 60          ! angular grid size in deg. 45
            n_opt = 4          ! # of final grad opts 15
         end if
         set8 = .false.
      case('qcg')
         if (set8) then
            maxparent = 50      ! # of parents in gene pool 100
            maxgen = 7          ! # of generations 10
!            mxcma = 250         ! R points in CMA search 1000
!            stepr = 4.0         ! R grid step in Bohr 2.5
!            stepa = 60          ! angular grid size in deg. 45
            n_opt = 5          ! # of final grad opts 15
            qcg = .true.
         end if
         set8 = .false.
      case ('noind')
         if (set9) mode = 1
         set9 = .false.
      case ('loose')
         if (set10) mode = 2
         set10 = .false.
      case ('cs')
         if (set11) cssym = .true.
         set11 = .false.
      case ('org')
         if (set12) incl_org = .true.
         set12 = .false.
      case ('ensemble')
         if (set13) docking_ens = .true.
         set13 = .false.
      case('attractive')
         if (set14) directed_type = p_atom_att
         set14 = .false.
      case('repulsive')
         if (set15) directed_type = p_atom_pot
         set15 = .false.
      end select
   end subroutine set_logicals

   subroutine set_optlvl(env)

      !> Calculation environment
      type(TEnvironment), intent(inout) :: env

      if (optlvl == 'gfn2') then
         call set_gfn(env, 'method', '2')
         call set_gfn(env, 'd4', 'true')
      end if

      if (optlvl == 'gfn1') then
         call set_gfn(env, 'method', '1')
      end if

      if (optlvl == 'gfn0') then
         call set_gfn(env, 'method', '0')
         call set_exttyp('eht')
      end if

      if (optlvl == 'gfnff') then
         call set_exttyp('ff')
      end if

   end subroutine set_optlvl

end module xtb_docking_set_module
