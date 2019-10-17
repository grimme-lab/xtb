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

subroutine test_thermo_calc
   use iso_fortran_env, wp => real64, istdout => output_unit
   use assertion

   use mctc_econv
   use thermo

   implicit none

   real(wp), parameter :: thr = 1.0e-8_wp
   integer,  parameter :: nat = 6
   integer,  parameter :: nvibs = 12
   real(wp), parameter :: ams(nat) = [15.99940492_wp,15.99940492_wp, &
      &  1.00794075_wp, 1.00794075_wp, 1.00794075_wp, 1.00794075_wp]
   real(wp), parameter :: moments(3) = &
      &[15649.5560990126_wp, 484991.786527778_wp, 485024.382510520_wp]
   real(wp), parameter :: vibs(nvibs) = &
      &[117.29518854_wp,      161.91929949_wp,      164.05816261_wp, &
      & 218.07506352_wp,      402.96868229_wp,      558.86387616_wp, &
      &1523.62725135_wp,     1561.31279119_wp,     3460.43559680_wp, &
      &3633.65239963_wp,     3636.81476002_wp,     3665.34271369_wp]*rcmtoau
   real(wp), parameter :: vibs_iso(nvibs) = &
      &[ 86.38105508090_wp,  141.77084264327_wp,  155.89354123588_wp, &
      & 182.95872682340_wp,  333.60535802793_wp,  421.06198802042_wp, &
      &1060.49921429103_wp, 1518.72673328712_wp, 2224.56212998692_wp, &
      &2512.05098623212_wp, 3621.54401864530_wp, 3625.88416620323_wp]*rcmtoau
   real(wp), parameter :: rotational_constants(3) = &
      & 0.5_wp*autorcm/[moments(3), moments(2), moments(1)]
   real(wp), parameter :: averaged_moment = &
      & sum(moments)/3 /(kgtome*aatoau**2*1.0e+20_wp)
   real(wp), parameter :: molecular_mass = sum(ams)
   logical, parameter :: atom = .false.
   logical, parameter :: linear = .false.
   real(wp), parameter :: symmetry_number = 1
   real(wp), parameter :: energy = 0.0_wp
   real(wp), parameter :: temperature = 298.15_wp
   real(wp), parameter :: rotor_cutoff = 100.0_wp
   logical, parameter :: pr = .true.
   real(wp), parameter :: zero_point_energy = 0.5_wp * sum(vibs)
   real(wp) :: et, ht, g, ts

   call thermodyn(output_unit,rotational_constants(1),rotational_constants(2), &
      &           rotational_constants(3),averaged_moment,linear,atom, &
      &           symmetry_number,molecular_mass,vibs,nvibs,energy,temperature, &
      &           rotor_cutoff,et,ht,g,ts,zero_point_energy,pr)

   call assert_close(et, 0.50275771916811E-01_wp, thr)
   call assert_close(ht, 0.67528241233247E-02_wp, thr)
   call assert_close(g,  0.17417501635698E-01_wp, thr)
   call assert_close(ts, 0.32858270281112E-01_wp, thr)

   call terminate(afail)
end subroutine test_thermo_calc

subroutine test_print_thermo
   use iso_fortran_env, wp => real64, istdout => output_unit
   use assertion

   use mctc_econv
   use splitparam

   use axis_trafo
   use property_output

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
   integer,  parameter :: nvibs = nat*3
   real(wp), parameter :: vibs(nvibs) = &
      &[  0.00000000_wp,        0.00000000_wp,        0.00000000_wp, &
      &   0.00000000_wp,        0.00000000_wp,        0.00000000_wp, &
      & 117.29518854_wp,      161.91929949_wp,      164.05816261_wp, &
      & 218.07506352_wp,      402.96868229_wp,      558.86387616_wp, &
      &1523.62725135_wp,     1561.31279119_wp,     3460.43559680_wp, &
      &3633.65239963_wp,     3636.81476002_wp,     3665.34271369_wp]
   real(wp), parameter :: vibs_iso(nvibs) = &
      &[  0.00000000000_wp,    0.00000000000_wp,    0.00000000000_wp, &
      &   0.00000000000_wp,    0.00000000000_wp,    0.00000000000_wp, &
      &  86.38105508090_wp,  141.77084264327_wp,  155.89354123588_wp, &
      & 182.95872682340_wp,  333.60535802793_wp,  421.06198802042_wp, &
      &1060.49921429103_wp, 1518.72673328712_wp, 2224.56212998692_wp, &
      &2512.05098623212_wp, 3621.54401864530_wp, 3625.88416620323_wp]
   real(wp), parameter :: energy = 0.0_wp
   real(wp) :: etot, htot, gtot
   integer  :: nimag

   atmass = ams

   call print_thermo(output_unit,nat,nvibs,at,xyz,vibs,energy, &
      &              htot,gtot,nimag,.true.)

   call assert_eq(nimag, 0)
   call assert_close(htot, 0.50275771916811E-01_wp, thr)
   call assert_close(gtot, 0.17337250373280E-01_wp, thr)

   atmass = iso

   call print_thermo(output_unit,nat,nvibs,at,xyz,vibs_iso,energy, &
      &              htot,gtot,nimag,.true.)

   call assert_eq(nimag, 0)
   call assert_close(htot, 0.43308512037251E-01_wp, thr)
   call assert_close(gtot, 0.88976699718681E-02_wp, thr)

   call terminate(afail)

end subroutine test_print_thermo
