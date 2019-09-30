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

module argparser
   use iso_fortran_env, wp => real64, istdout => output_unit
   use mctc_systools
   implicit none

   type :: string
      character(len=:),allocatable :: val
   endtype string

contains

subroutine rdxargs(fname,xcontrol,fnv,fnx,acc,lgrad,restart,gsolvstate,strict,  &
           &       copycontrol,argument_list,nargs,coffee)
   use set_module
   use readin, only : xfind,get_value
   use gbobc,  only : lgbsa
   implicit none
   character(len=:),allocatable,intent(out) :: fname
   character(len=:),allocatable,intent(out) :: xcontrol
   character(len=:),allocatable,intent(out) :: fnv
   character(len=:),allocatable,intent(out) :: fnx
   type(string),    allocatable,intent(out) :: argument_list(:)
   real(wp),intent(out) :: acc
   integer, intent(out) :: gsolvstate
   logical, intent(out) :: restart
   logical, intent(out) :: strict
   logical, intent(out) :: coffee
   logical, intent(out) :: lgrad
   logical, intent(out) :: copycontrol
   integer, intent(out) :: nargs

!$ integer  :: TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
!$ integer  :: nproc

   real(wp) :: ddum
   integer  :: iarg,skip,icount
   integer  :: id
   integer  :: idum
   integer  :: err
   integer  :: argc
   character(len=:),allocatable :: arg
   character(len=:),allocatable :: sec
   logical  :: exist

   logical  :: getopts

   coffee = .false.
   strict = .false.
   restart = .true.
   copycontrol = .false.
   lgrad = .false.
   acc = 1.0_wp
   gsolvstate = 0

   argc = command_argument_count()
   allocate( argument_list(argc) )

   getopts = .true.
   icount = 0
   skip = 0
   do iarg = 1, argc
      call rdarg(iarg,arg)
      exist = .false.
      if (skip.gt.0) then
         skip = skip-1
         argument_list(icount)%val = argument_list(icount)%val//" "//arg
         cycle
      endif
      if (arg.eq.'--') then
         getopts=.false.
         cycle
      endif
      if (getopts) then
!        write(output_unit,'(i0,'':'',x,a)') iarg,arg ! debugging stuff

         if ((len(arg).gt.2).and. &
         &  (index(arg,'--').eq.0).and.(index(arg,'-').eq.1)) then
            call raise('S',"the use of '"//arg//"' is discouraged, "// &
            &              "please use '-"//arg//"' next time",1)
            arg = '-'//arg
         endif
         select case(arg)
!! ========================================================================
!  check the help, version and citation flag, exit if found
!! ========================================================================
         case('-h','--help')
            call help!(istdout)
            call terminate(0)
         case(     '--citation')
            call citation!(istdout)
            call terminate(0)
         case(     '--license')
            call disclamer(istdout)
            call terminate(0)
         case(     '--version')
            call xtb_header(istdout)
            call terminate(0)

!! ========================================================================
!  now for true options, like verbosity or printouts
!! ========================================================================
         case('-v','--verbose')
            verbose = .true.
         case('-V','--very-verbose')
            verbose = .true.
            veryverbose = .true.
      !* this option is special, no calculation will be done,
      !  but all information about the run will be dumped and
      !  then the program exits with non-zero status if *any*
      !  warning was issued.
         case(     '--define')
            call set_define
      !$ case('-P','--parallel')
      !$    skip = 1
      !$    call rdarg(iarg+1,sec)
      !$    if (get_value(sec,idum)) then
      !$       nproc = omp_get_num_threads()
      !$       if (idum.gt.nproc) &
      !$       & call raise('S','Process number higher than OMP_NUM_THREADS, '//&
      !$       &                'I hope you know what you are doing.',1)
      !$       call omp_set_num_threads(idum)
      !$       call mkl_set_num_threads(idum)
      !$    endif
         case(     '--restart')
            restart = .true.
         case(     '--norestart')
            restart = .false.
         case(     '--copy')
            copycontrol = .true.
         case(     '--nocopy')
            copycontrol = .false.
         case(     '--strict')
            strict = .true.
         case('-I','--input','--by')
            skip = 1
            call rdarg(iarg+1,xcontrol)
         case(     '--namespace')
            skip = 1
            call rdarg(iarg+1,xenv%namespace)

         case(     '--vparam')
            skip = 1
            call rdarg(iarg+1,fnv)
         case(     '--xparam')
            skip = 1
            call rdarg(iarg+1,fnx)

         case(     '--coffee')
            coffee = .true.

!! ========================================================================
!  now we check for specifications concerning the version number,
!  runtyp and accuracy provided by the user via command line.
!  Be aware that this options interfere with the xcontrol and xtbrc.
!  If it ever happens that you have to implement a flag below here,
!  make sure you also make it work with the xcontrol-system.
!! ========================================================================

         case('-a','--acc')
            skip = 1
            call rdarg(iarg+1,sec)
            if (get_value(sec,ddum)) then
               if (ddum.lt.1.e-4_wp) then
                  call raise('S',"We cannot provide this level of accuracy, "//&
                                 "resetted accuracy to 0.0001",1)
                  acc = 1.e-4_wp
               else if (ddum.gt.1.e+3_wp) then
                  call raise('S',"We cannot provide this level of accuracy, "//&
                                 "resetted accuracy to 1000",1)
                  acc = 1.e+3_wp
               else
                  acc = ddum
               endif
            endif

         case('-c','--chrg')
            skip = 1
            call rdarg(iarg+1,sec)
            call set_chrg(sec)

         case('-u','--uhf')
            skip = 1
            call rdarg(iarg+1,sec)
            call set_spin(sec)

         case(     '--gfn')
            skip = 1
            call rdarg(iarg+1,sec)
            call set_gfn('method',sec)
         !  legacy stuff
            if(sec=='0')call set_exttyp('eht')    !ppracht 10/2018 GFN0
         case(     '--gfn1')
            call set_gfn('method','1')
            call raise('S',"The use of '"//arg//"' is discouraged, " //&
                           "please use '--gfn 1' next time",1)
         case(     '--gfn2')
            call set_gfn('method','2')
            call set_gfn('d4','true')

         case(     '--gfn0')
            call set_gfn('method','0')
            call set_exttyp('eht')    !ppracht 10/2018 GFN0
            call raise('S',"The use of '"//arg//"' is discouraged, " //&
                           "please use '--gfn 0' next time",1)

         case(     '--etemp')
            skip = 1
            call rdarg(iarg+1,sec)
            call set_scc('temp',sec)

         case(     '--esp')
            call set_runtyp('scc')
            call set_write('esp','true')

         case(     '--stm')
            call set_runtyp('scc')
            call set_write('stm','true')

         case(     '--periodic','--pbc')
            call set_gfn('periodic','true')
            call set_gfn('method','0')
         case(     '--cma')
            call set_cma

         case(     '--qmdff')
            call set_exttyp('qmdff')

         case(     '--tm')
            call set_exttyp('turbomole')

         case(     '--enso')
            call set_enso_mode

         case(     '--orca')
            call set_exttyp('orca')

         case(     '--mopac')
            call set_exttyp('mopac')

         case(     '--pop')
            call set_write('mulliken','true')

         case(     '--molden')
            call set_write('mos','true')

         case(     '--dipole')
            call set_write('dipole','true')

         case(     '--wbo')
            call set_write('wiberg','true')

         case(     '--lmo')
            call set_write('mulliken','true')
            call set_write('lmo','true')

         case(     '--ewin')
            call rdarg(iarg+1,sec)
            if (index(sec,'-').ne.1) then
                 skip = 1
                 call set_siman('ewin',sec)
            endif

         case(     '--fod')
            call set_write('fod','true')
            call set_scc('etemp','12500.0')

         case(     '--iterations')
            skip = 1
            call rdarg(iarg+1,sec)
            call set_scc('maxiterations',sec)

         case(     '--cycles')
            skip = 1
            call rdarg(iarg+1,sec)
            call set_opt('maxcycle',sec)

         case('-g','--gbsa')
            lgbsa = .true.
            if (iarg.lt.argc) then
               call rdarg(iarg+1,sec)
               if (index(sec,'-').ne.1) then
                  skip = 1
                  call set_gbsa('solvent',sec)

                  if (iarg+skip.lt.argc) then
                     call rdarg(iarg+skip+1,sec)
                     if (index(sec,'-').ne.1) then
                        if (sec == 'reference') then
                           skip = skip + 1
                           gsolvstate = 1
                        else if (sec == 'bar1M') then
                           skip = skip + 1
                           gsolvstate = 2
                        endif
                     endif
                  endif

                  if (iarg+skip.lt.argc) then
                     call rdarg(iarg+skip+1,sec)
                     if (index(sec,'-').ne.1) then
                        if (sec == 'normal'.OR.sec == 'tight'.OR.&
                           &sec == 'verytight'.OR.sec == 'extreme') then
                           skip = skip + 1
                           call set_gbsa('gbsagrid',sec)
                        endif
                     endif
                  endif

               endif
            endif

!! ------------------------------------------------------------------------
!  check for the runtyp, this will be determining what is done
!  in the main program
!! ------------------------------------------------------------------------
!        case(     '--nox','--nodiff')
!           call set_runtyp('nox')

!        case(     '--stda')
!           call set_runtyp('stda')

         case(     '--scc','--sp')
            call set_runtyp('scc')

         case(     '--vip')
            call set_gfn('method','1')
            call set_runtyp('vip')

         case(     '--vea')
            call set_gfn('method','1')
            call set_runtyp('vea')

         case(     '--vipea')
            call set_gfn('method','1')
            call set_runtyp('vipea')

          case(     '--vomega')
            call set_gfn('method','1')
            call set_runtyp('vomega')
        
          case(     '--vfukui')
            call set_runtyp('vfukui')
         
          case(     '--grad')
            call set_runtyp('grad')
            lgrad = .true.

         case('-o','--opt')
            call set_runtyp('opt')
            if (iarg.lt.argc) then
               call rdarg(iarg+1,sec)
               if (index(sec,'-').ne.1) then
                  skip = 1
                  call set_opt('optlevel',sec)
               endif
            endif

         case(     '--hess')
            call set_runtyp('hess')

         case(     '--md')
            call set_runtyp('md')


         case(     '--ohess')
            call set_runtyp('ohess')
            if (iarg.lt.argc) then
               call rdarg(iarg+1,sec)
               if (index(sec,'-').ne.1) then
                  skip = 1
                  call set_opt('optlevel',sec)
               endif
            endif

         case(     '--omd')
            call set_runtyp('omd')
            call set_opt('optlevel','-1')

         case(     '--siman')
            call set_runtyp('siman')
            call set_md('nvt','true')

         case(     '--path')
            call set_runtyp('path')
            if (iarg.lt.argc) then
               call rdarg(iarg+1,sec)
               if (index(sec,'-').ne.1) then
                  skip = 1
                  call set_path('product',sec)
               endif
            endif

        case(     '--screen')
           call set_runtyp('screen')

         case(     '--gmd')
            call set_runtyp('gmd')
            call raise('E',"This feature has been deprecated, I'm sorry.",1)

         case(     '--modef')
            call set_runtyp('modef')
            call rdarg(iarg+1,sec)
            if (index(sec,'-').ne.1) then
               skip = 1
               call set_modef('mode',sec)
            endif

         case(     '--mdopt')
            call set_runtyp('mdopt')

         case(     '--metadyn')
            call set_runtyp('md')
            if (iarg.lt.argc) then
               call rdarg(iarg+1,sec)
               if (index(sec,'-').ne.1) then
                  skip = 1
                  call set_metadyn('save',sec)
               endif
            endif
         case(     '--metaopt')
            call set_runtyp('metaopt')
            if (iarg.lt.argc) then
               call rdarg(iarg+1,sec)
               if (index(sec,'-').ne.1) then
                  skip = 1
                  call set_opt('optlevel',sec)
               endif
            endif

         case(     '--reactor')
            call raise('E','The nano-reactor has been moved to CREST',1)
            call set_runtyp('reactor')

!! ------------------------------------------------------------------------
!  no match => take as file name
!! ------------------------------------------------------------------------
         case default
            inquire(file=arg,exist=exist)
            if (exist) then
               if (allocated(fname)) call raise('S',  &
               &  "There are multiple files provided. '"//fname//  &
               &  "' will be ignored in this run.",1)
               fname = arg
            else
               if (index(arg,'-').eq.1) then
                  call raise('S',"Unfortunately, '"//arg// &
                  &    "' is not supported in this program. Check with --help.",1)
               else
                  call raise('E',"You don't have a file named '"//arg// &
                  &    "' here",1)
               endif
            endif

         end select

      else ! getopts?
         inquire(file=arg,exist=exist)
         if (exist) then
            if (allocated(fname)) call raise('S',  &
            &  "There are multiple files provided. '"//fname//  &
            &  "' will be ignored in this run.",1)
            fname = arg
         else
            call raise('E',"You don't have a file named '"//arg//"' here.",1)
         endif
      endif ! getopts?
      icount = icount+1
      if (exist) then
         argument_list(icount)%val = "file '"//arg//"'"
      else
         argument_list(icount)%val = "flag "//arg
      endif
   enddo

   nargs = icount

   if (.not.(allocated(fname).or.coffee)) then
      if (periodic) then
         fname = 'POSCAR'
      else
         call raise('E',"No geometry given, so there is nothing to do.",1)
      endif
   endif

!  make sure that we get a eht calculation instead of a scc for GFN0
   if(gfn_method.eq.0)  call set_exttyp('eht')

   if (.not.allocated(xcontrol)) then
      if (.not.copycontrol) then
         xcontrol = fname
      else
         xcontrol = 'xcontrol'
      endif
   endif


end subroutine rdxargs

end module argparser
