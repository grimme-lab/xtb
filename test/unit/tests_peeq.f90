program peeq_tester
   use xtb_mctc_accuracy, only : wp
!$ use omp_lib
! ------------------------------------------------------------------------
!  general purpose library
! ------------------------------------------------------------------------
   use xtb_mctc_global
   use xtb_mctc_timings
   use xtb_mctc_systools
   use xtb_mctc_convert

   implicit none

! ------------------------------------------------------------------------
!  local variables
! ------------------------------------------------------------------------
   integer  :: testid
   integer  :: idum,nargs
   character(len=:),allocatable :: arg,sec

!! ------------------------------------------------------------------------
!!  signal processing
!! ------------------------------------------------------------------------
!   external :: wSIGTERM
!   external :: wSIGINT
!!  here two important signal handlers are installed, it seems that
!!  FORTRAN by itself cannot handle signals in the way I expected it
!!  to do, but this will force it to die when I hit CTRL^C.
!   call signal(2,wSIGINT)
!   call signal(15,wSIGTERM)

   call mctc_init('test',10,.true.)

   nargs = command_argument_count()
   if (nargs.lt.2) then
      call raise('E',"Please give the tester a test to run!")
   endif

   call rdarg(1,arg)
   call rdarg(2,sec)

   select case(arg)
   case('repulsion')
      select case(sec)
      case('cluster'); call test_repulsion_cluster
      case('pbc3d'); call test_repulsion_pbc3d
      end select
   case('gfn2')
      select case(sec)
      case('basic'); call test_gfn2_mindless_basic
      case('solvation'); call test_gfn2_mindless_solvation
      case('scc'); call test_gfn2_scc
      case('api'); call test_gfn2_api
      case('gbsa'); call test_gfn2gbsa_api
      case('salt'); call test_gfn2salt_api
      case('pcem'); call test_gfn2_pcem_api
      case('d-metal'); call test_gfn2_dmetal
      case('cosmo'); call test_gfn2_mindless_cosmo
      end select
   case('gfn1')
      select case(sec)
      case('basic'); call test_gfn1_mindless_basic
      case('solvation'); call test_gfn1_mindless_solvation
      case('scc'); call test_gfn1_scc
      case('api'); call test_gfn1_api
      case('gbsa'); call test_gfn1gbsa_api
      case('pcem'); call test_gfn1_pcem_api
      case('xb'); call test_gfn1_xb
      case('pbc3d'); call test_gfn1_pbc3d
      case('ipea'); call test_ipea_indole
      case('cosmo'); call test_gfn1_mindless_cosmo
      end select
   case('gfn0')
      select case(sec)
      case('basic'); call test_gfn0_mindless_basic
      case('solvation'); call test_gfn0_mindless_solvation
      case('sp');  call test_gfn0_sp
      case('api'); call test_gfn0_api
      case('srb'); call test_gfn0_api_srb
      end select
   case('gfnff')
      select case(sec)
      case('basic'); call test_gfnff_mindless_basic
      case('solvation'); call test_gfnff_mindless_solvation
      case('scaleup'); call test_gfnff_scaleup
      case('pdb'); call test_gfnff_pdb
      case('sdf'); call test_gfnff_sdf
      case('sp');  call test_gfnff_sp
      case('hb');  call test_gfnff_hb
      case('gbsa');call test_gfnff_gbsa
      end select
   case('peeq')
      select case(sec)
      case('sp');  call test_peeq_sp
      case('api'); call test_peeq_api
      case('srb'); call test_peeq_api_srb
      end select
   case('pbc_tools')
      select case(sec)
      case('convert'); call test_pbc_tools_convert
      case('cutoff');  call test_pbc_tools_cutoff
      end select
   case('xtb_type_atomlist')
      select case(sec)
      case('list'); call test_atomlist
      end select
   case('latticepoint')
      select case(sec)
      case('pbc3d'); call test_latticepoint_pbc3d
      end select
   end select

   ! falling through the tester is always an error
   call terminate(1)

end program peeq_tester
