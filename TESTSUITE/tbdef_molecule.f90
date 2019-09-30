subroutine test_class_molecule_mic_distances
   use iso_fortran_env, wp => real64, istdout => output_unit
   use assertion
   use tbdef_molecule
   use tbdef_param
   use eeq_model
   use ncoord
   use pbc_tools
   implicit none
   real(wp),parameter :: thr = 1.0e-10_wp
   ! SiO2 (random, no symmetry)
   integer, parameter :: nat = 6
   integer, parameter :: at(nat) = [14,14,8,8,8,8]
   real(wp),parameter :: abc(3,nat) = reshape(&
      &[.095985472469032_wp,     .049722204206931_wp,     0.10160624337938_wp, &
      & 0.54722204206931_wp,     0.52863628207623_wp,     0.38664208660311_wp, &
      & 0.29843937068984_wp,     0.39572194413818_wp,     0.20321248675876_wp, &
      & 0.23364982659922_wp,     0.85647058758674_wp,     0.31884968761485_wp, &
      & 0.72250232459952_wp,     0.65548544066844_wp,     .056207709103487_wp, &
      & 0.70514214000043_wp,     0.28321754549582_wp,     0.36424822189074_wp],&
      & shape(abc))
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[ 8.7413053236641_wp,      0.0000000000000_wp,      0.0000000000000_wp,   &
      &  0.0000000000000_wp,      8.7413053236641_wp,      0.0000000000000_wp,   &
      &  0.0000000000000_wp,      0.0000000000000_wp,      8.7413053236641_wp],  &
      & shape(lattice))
   integer, parameter :: wsc_rep(3) = [1,1,1]

   type(tb_molecule)       :: mol

   call mol%allocate(nat)
   mol%at   = at
   mol%abc  = abc
   mol%npbc = 3
   mol%pbc  = .true.
   mol%lattice = lattice
   mol%volume = dlat_to_dvol(lattice)
   call dlat_to_cell(lattice,mol%cellpar)
   call dlat_to_rlat(lattice,mol%rec_lat)
   call coord_trafo(nat,lattice,abc,mol%xyz)
   call mol%wrap_back
   call mol%calculate_distances

   call assert_close(mol%dist(1,1),8.7413053236641_wp,thr)
   call assert_close(mol%dist(3,6),3.9480992656526_wp,thr)
   call assert_close(mol%dist(2,3),2.9411012169549_wp,thr)
   call assert_close(mol%dist(1,4),2.8120981086302_wp,thr)

   call mol%deallocate

   call terminate(afail)

end subroutine test_class_molecule_mic_distances

subroutine test_class_molecule_axis_trafo
   use iso_fortran_env, wp => real64
   use assertion
   use tbdef_molecule
   use tbdef_param
   use ncoord
   use eeq_model
   implicit none
   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 7
   integer, parameter :: at(nat) = [6,8,1,6,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.236995341_wp,     0.794935672_wp,     0.000000000_wp,  &
      & 2.159203072_wp,    -0.439230068_wp,     0.000000000_wp,  &
      & 0.329518229_wp,     2.884317091_wp,     0.000000000_wp,  &
      &-2.377225608_wp,    -0.313030650_wp,     0.000000000_wp,  &
      &-2.288077922_wp,    -2.365932814_wp,     0.000000000_wp,  &
      &-3.405901173_wp,     0.355579334_wp,     1.660010899_wp,  &
      &-3.405901173_wp,     0.355579335_wp,    -1.660010899_wp   &
      & ],shape(xyz))

   type(tb_molecule)       :: mol
   real(wp) :: molmass
   real(wp) :: center(3)
   real(wp) :: moments(3)

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp
   call mol%set_atomic_masses

   center = mol%center_of_geometry()
   call assert_close(center(1),-1.2501984620000_wp,thr)
   call assert_close(center(2), .18174541428571_wp,thr)
   call assert_close(center(3), 0.0000000000000_wp,thr)

   call mol%shift_to_center_of_geometry
   call assert_close(mol%xyz(2,1),.61319025771429_wp,thr)
   call assert_close(mol%xyz(1,3),1.5797166910000_wp,thr)
   call assert_close(mol%xyz(3,6),1.6600108990000_wp,thr)

   molmass = sum(mol%atmass)
   call assert_close(molmass,80303.049694083_wp,1.0e-5_wp) ! in au
   center = mol%center_of_mass()
   call assert_close(center(1),1.2502034210661_wp,thr)

   call mol%shift_to_center_of_mass
   call assert_close(mol%xyz(2,4),-0.31302884934293_wp,thr)
   call assert_close(mol%xyz(1,1), 0.23699038193392_wp,thr)

   moments = mol%moments_of_inertia()
   call assert_close(moments(1),57768.744315301_wp,1.0e-5_wp)
   call assert_close(moments(3),361019.51690412_wp,1.0e-5_wp)

!   call mol%align_to_principal_axes(break_symmetry = .false.)
!   call assert_close(mol%xyz(2,1),-.79493732581159_wp,thr)
!   call assert_close(mol%xyz(1,7), 3.4059059117382_wp,thr)
!
!   call mol%align_to_principal_axes(break_symmetry = .true.)
!   call assert_close(mol%xyz(2,1),-.79493703211973_wp,thr)
!   call assert_close(mol%xyz(1,7), 3.4059054710786_wp,thr)

   call mol%deallocate

   call terminate(afail)

end subroutine test_class_molecule_axis_trafo
