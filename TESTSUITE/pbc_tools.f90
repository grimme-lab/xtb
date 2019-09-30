subroutine test_pbc_tools_convert
   use iso_fortran_env, wp => real64
   use assertion
   use pbc_tools
   implicit none
   real(wp),parameter :: thr = 1.0e-10_wp
   real(wp),parameter :: cellpar(6) = &
      [9.09903131_wp, 9.09903131_wp, 30.46049560_wp, 90.0_wp, 90.0_wp, 120.0_wp]
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[5.9598811567890_wp,      2.1071361905157_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      6.3214085715472_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      0.0000000000000_wp,      7.2993338808807_wp],   &
      & shape(lattice))
   real(wp) :: dlat(3,3),rlat(3,3),volume,cpar(6),center(3)

   call cell_to_dlat(cellpar,dlat)
   call assert_close(cellpar(1),dlat(1,1),thr)
   call assert_close(dlat(1,3),-13.648544412579_wp,thr)
   call assert_close(dlat(3,3), 26.878966813429_wp,thr)

   call cell_to_rlat(cellpar,rlat)
   call assert_close(rlat(1,1), 0.69053343077018_wp,thr)
   call assert_close(rlat(2,1),-0.96832302603374_wp,thr)
   call assert_close(rlat(3,2), 0.19327597283376_wp,thr)

   volume = cell_to_dvol(cellpar)
   call assert_close(volume,1292.0766773144_wp,thr)

   call dlat_to_cell(lattice,cpar)
   call assert_close(cpar(1),cpar(2),thr)
   call assert_close(cpar(3),7.2993338808807_wp,thr)
   call assert_close(cpar(4),1.0471975511966_wp,thr)

   call dlat_to_rlat(lattice,rlat)
   call assert_close(rlat(1,1), 1.0542467445047_wp,thr)
   call assert_close(rlat(1,3),-.35141558150155_wp,thr)
   call assert_close(rlat(2,2), .99395336277746_wp,thr)

   volume = dlat_to_dvol(lattice)
   call assert_close(volume, 275.00126402469_wp,thr)

   center = get_center_dlat(lattice)
   call assert_close(center(1),2.9799405783945_wp,thr)
   call assert_close(center(2),4.2142723810314_wp,thr)

   call terminate(afail)

end subroutine test_pbc_tools_convert

subroutine test_pbc_tools_cutoff
   use iso_fortran_env, wp => real64
   use assertion
   use pbc_tools
   use pbc
   implicit none
   real(wp),parameter :: thr = 1.0e-10_wp
   real(wp),parameter :: lattice_1(3,3) = reshape(&
      &[9.0990313100000_wp,      0.0000000000000_wp,      0.0000000000000_wp,    &
      & 7.4082581428274_wp,      5.2829993440840_wp,      0.0000000000000_wp,    &
      &-13.648544412579_wp,     -4.3680854682659_wp,      26.878966813429_wp],   &
      & shape(lattice_1))
   real(wp),parameter :: lattice_2(3,3) = reshape(&
      &[5.9598811567890_wp,      2.1071361905157_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      6.3214085715472_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      0.0000000000000_wp,      7.2993338808807_wp],   &
      & shape(lattice_2))
   real(wp),parameter :: rthr_1 =  40.0_wp ** 2
   real(wp),parameter :: rthr_2 =  70.0_wp ** 2
   real(wp),parameter :: rthr_3 = 100.0_wp ** 2
   integer :: rep(3)

   call get_realspace_cutoff(lattice_1,rthr_1,rep)
   call assert_eq(rep(1), 8)
   call assert_eq(rep(3), 2)

   call get_realspace_cutoff(lattice_1,rthr_2,rep)
   call assert_eq(rep(2),14)
   call assert_eq(rep(3), 3)

   call get_realspace_cutoff(lattice_1,rthr_3,rep)
   call assert_eq(rep(1),20)
   call assert_eq(rep(2),20)
   call assert_eq(rep(3), 4)

   call get_realspace_cutoff(lattice_2,rthr_1,rep)
   call assert_eq(rep(1),rep(2))
   call assert_eq(rep(3), 7)

   call get_realspace_cutoff(lattice_2,rthr_2,rep)
   call assert_eq(rep(1),rep(3))
   call assert_eq(rep(2),12)

   call get_realspace_cutoff(lattice_2,rthr_3,rep)
   call assert_eq(rep(3),rep(2))
   call assert_eq(rep(1),17)

   call terminate(afail)

end subroutine test_pbc_tools_cutoff
