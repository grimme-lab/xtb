program peeq_tester
   use iso_fortran_env, istdout => output_unit, kdp => real64
!$ use omp_lib
! ------------------------------------------------------------------------
!  general purpose library
! ------------------------------------------------------------------------
   use mctc_global
   use mctc_timings
   use mctc_systools
   use mctc_econv

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
   case('eeq_model')
      select case(sec)
      case('water'); call test_eeq_model_water
      case('ewald'); call test_eeq_model_ewald
      case('gbsa');  call test_eeq_model_gbsa
      case('hbond'); call test_eeq_model_hbond
      case('salt');  call test_eeq_model_salt
      end select
   case('gfn2')
      select case(sec)
      case('scc'); call test_gfn2_scc
      case('api'); call test_gfn2_api
      case('gbsa'); call test_gfn2gbsa_api
      case('salt'); call test_gfn2salt_api
      case('pcem'); call test_gfn2_pcem_api
      end select
   case('gfn1')
      select case(sec)
      case('scc'); call test_gfn1_scc
      case('api'); call test_gfn1_api
      case('gbsa'); call test_gfn1gbsa_api
      case('pcem'); call test_gfn1_pcem_api
      end select
   case('gfn0')
      select case(sec)
      case('sp');  call test_gfn0_sp
      case('api'); call test_gfn0_api
      end select
   case('peeq')
      select case(sec)
      case('sp');  call test_peeq_sp
      case('api'); call test_peeq_api
      end select
   case('dftd4')
      select case(sec)
      case('properties'); call test_dftd4_properties
      case('energies');   call test_dftd4_energies
      case('pbc_disp');   call test_dftd4_pbc_energies
      case('cell_gradient'); call test_dftd4_cell_gradient
      case('api');           call test_dftd4_api
      case('pbc_api');       call test_dftd4_pbc_api
      end select
   case('geometry_reader')
      select case(sec)
      case('coord_3d_a'); call test_geometry_reader_file_coord_CaMgCO_3d
      case('coord_3d_b'); call test_geometry_reader_file_coord_CaF2_3d
      case('coord_2d');   call test_geometry_reader_file_coord_C_2d
      case('coord_1d');   call test_geometry_reader_file_coord_C_1d
      case('coord_0d');   call test_geometry_reader_file_coord_general_0d
      case('xmol_0d');    call test_geometry_reader_file_xmol_water_0d
      case('poscar_3d');  call test_geometry_reader_file_poscar_sio2_3d
      case('molfile'); call test_geometry_reader_molfile_benzen
      case('molfile_flat'); call test_geometry_reader_molfile_benzen_flat
      case('sdfile'); call test_geometry_reader_file_sdf_h2o
      case('sdfile_flat'); call test_geometry_reader_file_sdf_h2o_flat
      case('sdfile_noh'); call test_geometry_reader_file_sdf_benzen_hquery
      case('pdb'); call test_geometry_reader_file_pdb_4qxx
      case('pdb_noh'); call test_geometry_reader_file_pdb_4qxx_noh
      end select
   case('pbc_tools')
      select case(sec)
      case('convert'); call test_pbc_tools_convert
      case('cutoff');  call test_pbc_tools_cutoff
      end select
   case('tbdef_molecule')
      select case(sec)
      case('mic');  call test_class_molecule_mic_distances
      case('axis'); call test_class_molecule_axis_trafo
      end select
   case('tbdef_wsc')
      select case(sec)
      case('0d'); call test_wigner_seitz_0d
      case('3d'); call test_wigner_seitz_3d
      end select
   case('tbdef_atomlist')
      select case(sec)
      case('list'); call test_atomlist
      end select
   case('symmetry')
      select case(sec)
      case('water');  call test_symmetry_water
      case('li8'); call test_symmetry_li8
      case('pcl3'); call test_symmetry_pcl3
      case('c20'); call test_symmetry_c20
      end select
   case('thermo')
      select case(sec)
      case('axis'); call test_axis
      case('calc'); call test_thermo_calc
      case('print'); call test_print_thermo
      end select
   end select

   ! falling through the tester is always an error
   call terminate(1)

end program peeq_tester
