subroutine test_axis
   use iso_fortran_env, wp => real64, istdout => output_unit
   use assertion

   use mctc_econv
   use tbdef_molecule
   use splitparam

   use axis_trafo

   implicit none

   real(wp), parameter :: thr = 1.0e-8_wp
   real(wp), parameter :: thr2 = 1.0e-5_wp
   integer,  parameter :: nat = 6
   integer,  parameter :: at(nat) = [8,8,1,1,1,1]
   real(wp), parameter :: ams(nat) = [15.99940492_wp,15.99940492_wp, &
      &  1.00794075_wp, 1.00794075_wp, 1.00794075_wp, 1.00794075_wp]
   real(wp), parameter :: iso(nat) = [16.99913170_wp,17.9991610_wp, &
      &  2.01107778_wp, 3.01042777_wp, 1.00794075_wp, 1.00794075_wp]
   real(wp), parameter :: xyz(3,nat) = reshape([&
      &-2.98956539380840_wp,-0.24213847149862_wp,-0.00000088744317_wp, &
      & 2.35000755760470_wp, 0.22356837963350_wp,-0.00000838160401_wp, &
      &-1.19731697901940_wp, 0.13238632078926_wp,-0.00000485980301_wp, &
      &-3.84082093208368_wp, 1.35616304830930_wp, 0.00002536003281_wp, &
      & 2.83883868349296_wp,-0.73499856747260_wp,-1.46075352362473_wp, &
      & 2.83885706381384_wp,-0.73498070976087_wp, 1.46074229244211_wp],&
      & shape(xyz))
   type(tb_molecule) :: mol
   real(wp) :: moments(3), rot(3)
   real(wp) :: rot1(3)
   real(wp) :: rot2(3), avmom2, mass2
   real(wp) :: eval3(3), coord3(3,nat)
   real(wp) :: rot4(3), evec4(3,3)

   call mol%allocate(nat)
   mol%at = at
   mol%xyz = xyz
   mol%atmass = ams * amutoau
   call mol%update

   moments = mol%moments_of_inertia()

   call assert_close(moments(1), 15649.5560990126_wp, thr2)
   call assert_close(moments(2), 484991.786527778_wp, thr2)
   call assert_close(moments(3), 485024.382510520_wp, thr2)

   rot = 0.5_wp*autorcm/[moments(3), moments(2), moments(1)]

   call assert_close(rot(1), .226251131473004_wp, thr)
   call assert_close(rot(2), .226266337664493_wp, thr)
   call assert_close(rot(3), 7.01216792608729_wp, thr)

   atmass = ams

   call axis(mol%n, mol%at, mol%xyz, rot1(1), rot1(2), rot1(3))

   call assert_close(rot1(1), .226251131473004_wp, thr2)
   call assert_close(rot1(2), .226266337664493_wp, thr2)
   call assert_close(rot1(3), 7.01216792608729_wp, thr2)

   call axis2(mol%n, mol%at, mol%xyz, rot2(1), rot2(2), rot2(3), avmom2, mass2)

   call assert_close(rot2(1), .226251131473004_wp, thr2)
   call assert_close(rot2(2), .226266337664493_wp, thr2)
   call assert_close(rot2(3), 7.01216792608729_wp, thr2)

   call axisvec(mol%n, mol%at, mol%xyz, rot4(1), rot4(2), rot4(3), evec4)

   call assert_close(rot4(1), .226251131473004_wp, thr2)
   call assert_close(rot4(2), .226266337664493_wp, thr2)
   call assert_close(rot4(3), 7.01216792608729_wp, thr2)

   mol%atmass = iso * amutoau

   moments = mol%moments_of_inertia()

   call assert_close(moments(1), 23062.951395017_wp, thr2)
   call assert_close(moments(2), 569661.08474644_wp, thr2)
   call assert_close(moments(3), 577041.88473780_wp, thr2)

   rot = 0.5_wp*autorcm/[moments(3), moments(2), moments(1)]

   call assert_close(rot(1), .19017218374861_wp, thr)
   call assert_close(rot(2), .19263614502269_wp, thr)
   call assert_close(rot(3), 4.7581644454539_wp, thr)

   atmass = iso

   call axis(mol%n, mol%at, mol%xyz, rot1(1), rot1(2), rot1(3))

   call assert_close(rot1(1), .19017218374861_wp, thr2)
   call assert_close(rot1(2), .19263614502269_wp, thr2)
   call assert_close(rot1(3), 4.7581644454539_wp, thr2)

   call axis2(mol%n, mol%at, mol%xyz, rot2(1), rot2(2), rot2(3), avmom2, mass2)

   call assert_close(rot2(1), .19017218374861_wp, thr2)
   call assert_close(rot2(2), .19263614502269_wp, thr2)
   call assert_close(rot2(3), 4.7581644454539_wp, thr2)

   call axisvec(mol%n, mol%at, mol%xyz, rot4(1), rot4(2), rot4(3), evec4)

   call assert_close(rot4(1), .19017218374861_wp, thr2)
   call assert_close(rot4(2), .19263614502269_wp, thr2)
   call assert_close(rot4(3), 4.7581644454539_wp, thr2)

   call terminate(afail)

end subroutine test_axis
