! This file is part of xtb.
! SPDX-Identifier: LGPL-3.0-or-later
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

module test_coulomb
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
   implicit none
   private

   public :: collect_coulomb

contains

!> Collect all exported unit tests
subroutine collect_coulomb(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("point-cluster", test_coulomb_point_cluster), &
      new_unittest("point-pbc3d", test_coulomb_point_pbc3d), &
      new_unittest("gfn1-cluster", test_coulomb_gfn1_cluster), &
      new_unittest("gfn1-pbc3d", test_coulomb_gfn1_pbc3d), &
      new_unittest("gfn2-cluster", test_coulomb_gfn2_cluster), &
      new_unittest("gfn2-pbc3d", test_coulomb_gfn2_pbc3d), &
      new_unittest("gaussian-cluster", test_coulomb_gaussian_cluster), &
      new_unittest("gaussian-pbc3d", test_coulomb_gaussian_pbc3d) &
      ]

end subroutine collect_coulomb


subroutine test_coulomb_point_cluster(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_la, only : contract
   use xtb_type_coulomb
   use xtb_type_environment
   use xtb_type_molecule
   type(error_type), allocatable, intent(out) :: error
   real(wp), parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 24
   integer, parameter :: at(nat) = [6,7,6,7,6,6,6,8,7,6,8,7,6,6, &
      &                             1,1,1,1,1,1,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[ 2.02799738646442_wp,  0.09231312124713_wp, -0.14310895950963_wp, &
      &  4.75011007621000_wp,  0.02373496014051_wp, -0.14324124033844_wp, &
      &  6.33434307654413_wp,  2.07098865582721_wp, -0.14235306905930_wp, &
      &  8.72860718071825_wp,  1.38002919517619_wp, -0.14265542523943_wp, &
      &  8.65318821103610_wp, -1.19324866489847_wp, -0.14231527453678_wp, &
      &  6.23857175648671_wp, -2.08353643730276_wp, -0.14218299370797_wp, &
      &  5.63266886875962_wp, -4.69950321056008_wp, -0.13940509630299_wp, &
      &  3.44931709749015_wp, -5.48092386085491_wp, -0.14318454855466_wp, &
      &  7.77508917214346_wp, -6.24427872938674_wp, -0.13107140408805_wp, &
      & 10.30229550927022_wp, -5.39739796609292_wp, -0.13672168520430_wp, &
      & 12.07410272485492_wp, -6.91573621641911_wp, -0.13666499342053_wp, &
      & 10.70038521493902_wp, -2.79078533715849_wp, -0.14148379504141_wp, &
      & 13.24597858727017_wp, -1.76969072232377_wp, -0.14218299370797_wp, &
      &  7.40891694074004_wp, -8.95905928176407_wp, -0.11636933482904_wp, &
      &  1.38702118184179_wp,  2.05575746325296_wp, -0.14178615122154_wp, &
      &  1.34622199478497_wp, -0.86356704498496_wp,  1.55590600570783_wp, &
      &  1.34624089204623_wp, -0.86133716815647_wp, -1.84340893849267_wp, &
      &  5.65596919189118_wp,  4.00172183859480_wp, -0.14131371969009_wp, &
      & 14.67430918222276_wp, -3.26230980007732_wp, -0.14344911021228_wp, &
      & 13.50897177220290_wp, -0.60815166181684_wp,  1.54898960808727_wp, &
      & 13.50780014200488_wp, -0.60614855212345_wp, -1.83214617078268_wp, &
      &  5.41408424778406_wp, -9.49239668625902_wp, -0.11022772492007_wp, &
      &  8.31919801555568_wp, -9.74947502841788_wp,  1.56539243085954_wp, &
      &  8.31511620712388_wp, -9.76854236502758_wp, -1.79108242206824_wp],&
      &  shape(xyz))
   real(wp), parameter :: charges(nat) = [&
      &-0.05476196_wp, 0.00248420_wp, 0.08809027_wp,-0.28313223_wp, 0.12260575_wp, &
      &-0.04041664_wp, 0.26374283_wp,-0.43589031_wp,-0.11704515_wp, 0.31376962_wp, &
      &-0.44064154_wp,-0.07607650_wp,-0.05087676_wp,-0.04160170_wp, 0.06582772_wp, &
      & 0.08200852_wp, 0.08188349_wp, 0.05475583_wp, 0.09364446_wp, 0.07149007_wp, &
      & 0.07154523_wp, 0.09011030_wp, 0.06924195_wp, 0.06924253_wp]

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TCoulomb) :: coulomb

   integer :: ii, jj
   logical :: fail
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp), allocatable :: gradient(:, :)
   real(wp), allocatable :: jmat(:, :)
   real(wp), allocatable :: shift(:)
   real(wp), allocatable :: djdr(:, :, :)
   real(wp), allocatable :: djdtr(:, :)
   real(wp), allocatable :: djdL(:, :, :)

   call init(env)
   call init(mol, at, xyz)
   call init(coulomb, env, mol)

   allocate(shift(nat))
   allocate(jMat(nat, nat))
   allocate(gradient(3, nat))
   allocate(djdr(3, nat, nat))
   allocate(djdtr(3, nat))
   allocate(djdL(3, 3, nat))

   call coulomb%getCoulombMatrix(mol, jmat)

   do ii = 1, nat
      call check(error, jmat(ii,ii), 0.0_wp, thr=thr)
      do jj = 1, ii-1
         call check(error, jmat(jj,ii), jmat(ii,jj), thr=thr)
      end do
   end do

   call check(error, jmat( 1, 3), 0.21100723251445_wp, thr=thr)
   call check(error, jmat( 2, 3), 0.38630364669842_wp, thr=thr)
   call check(error, jmat( 4, 6), 0.23442574831659_wp, thr=thr)
   call check(error, jmat( 3,10), 0.11824472153399_wp, thr=thr)
   call check(error, jmat( 6, 3), 0.24063746843517_wp, thr=thr)
   call check(error, jmat( 7,18), 0.11492593111318_wp, thr=thr)
   call check(error, jmat(12,20), 0.25392041646803_wp, thr=thr)

   shift(:) = matmul(jmat, charges)
   energy = 0.5_wp*dot_product(charges, shift)
   call check(error, energy, -0.16130778864155_wp, thr=thr)

   call coulomb%getCoulombDerivs(mol, charges, djdr, djdtr, djdL)
   call contract(djdr, charges, gradient)
   call contract(djdL, charges, sigma)

   ! check numerical gradient
   do ii = 1, nat
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call check(error, gradient(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   ! check numerical strain derivatives
   eps = unity
   do ii = 1, 3
      do jj = 1, ii
         eps(jj, ii) = eps(jj, ii) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) + step

         call check(error, sigma(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

end subroutine test_coulomb_point_cluster

subroutine test_coulomb_point_pbc3d(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi, sqrtpi
   use xtb_mctc_convert
   use xtb_mctc_la, only : contract
   use xtb_type_coulomb
   use xtb_type_environment
   use xtb_type_molecule
   type(error_type), allocatable, intent(out) :: error

   real(wp),parameter :: thr = 1.0e-9_wp
   integer, parameter :: nat = 32
   integer, parameter :: at(nat) = [8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[4.5853168464880421_wp,  4.9392326929575878_wp,  4.1894081210748118_wp,  &
      & 5.8862267152491423_wp,  1.1425258978245871_wp,  6.5015058204768126_wp,  &
      & 1.4284279616220412_wp,  4.6017511285540875_wp,  3.1465884436348119_wp,  &
      & 2.3323404704521411_wp,  1.4471801154820869_wp,  0.7121932185858125_wp,  &
      & 4.7333543155561415_wp,  2.7747291872305868_wp,  3.1951352976178122_wp,  &
      & 5.6617754419101418_wp,  2.2133191164485870_wp,  2.1235404838618126_wp,  &
      & 5.3107598618381422_wp,  2.6902988056185868_wp,  0.7384466319968125_wp,  &
      & 4.4947071071761426_wp,  3.9530790692635867_wp,  0.6801776747598124_wp,  &
      & 4.8005171923760424_wp,  4.9185102874975870_wp,  1.8186363449528122_wp,  &
      & 4.6951362687070421_wp,  4.2781752812835867_wp,  3.1816411821728123_wp,  &
      & 1.3838574419160412_wp,  5.2817805008910863_wp,  4.1482702947948136_wp,  &
      & 1.0268974195990415_wp,  4.7234752637800881_wp,  5.4989995400388123_wp,  &
      & 2.0852659694760409_wp,  5.0956317453800875_wp,  6.5351699846458127_wp,  &
      & 2.3344644666691412_wp, -0.0736561690909131_wp,  6.4245628001158135_wp,  &
      & 2.4894017448231409_wp,  0.6213510313930869_wp,  5.0967297417158131_wp,  &
      & 1.5745272273791413_wp,  0.1243470825760870_wp,  3.9731040773988129_wp,  &
      & 5.8221065925130420_wp,  5.3013563342055878_wp,  1.7264876737078123_wp,  &
      & 3.4487807319551416_wp,  3.6355832152975864_wp,  0.7429568016758125_wp,  &
      & 4.8499393376520423_wp,  3.4713855169305874_wp,  6.4691872586348129_wp,  &
      & 0.2495364434351412_wp,  2.4795455690160870_wp,  2.1043557230378123_wp,  &
      & 5.6691068338331423_wp,  1.1234174220755870_wp,  2.1414388326468128_wp,  &
      & 3.7072009289431418_wp,  2.4357632918535872_wp,  3.0094700999208119_wp,  &
      & 4.1414520030430415_wp,  5.7877262477775879_wp,  1.7803680119358125_wp,  &
      & 5.0142851411171421_wp,  2.4165926460955873_wp,  4.1857610486448129_wp,  &
      & 3.0280930003030413_wp,  4.6201081184690871_wp,  6.2533190952188136_wp,  &
      & 0.5863628696651412_wp,  0.5757236365910867_wp,  4.1021714214668128_wp,  &
      & 2.3776130524831411_wp,  1.6969724987740866_wp,  5.2327688986668139_wp,  &
      & 1.9486148363011413_wp,  0.4390675147070869_wp,  2.9999022491838123_wp,  &
      & 3.5312997625581413_wp,  0.4467415528495868_wp,  4.8114121395028135_wp,  &
      & 6.5089895990100421_wp,  5.2100409408535882_wp,  6.0066553789008132_wp,  &
      & 0.9001165013630412_wp,  3.6420787128610868_wp,  5.4413106648508132_wp,  &
      & 1.6012116650460413_wp,  5.6845471271780879_wp,  0.7675566847298124_wp], &
      & shape(xyz)) * aatoau
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[6.4411018522600_wp,    0.0492571261505_wp,    0.2192046129910_wp,  &
      & 0.0462076831749_wp,    6.6435057067500_wp,    0.1670513770770_wp,  &
      & 0.2262248220170_wp,   -0.9573234940220_wp,    6.7608039126200_wp], &
      & shape(lattice)) * aatoau
   real(wp),parameter :: charges(nat) = [&
      &-0.40074825616499E+0_wp,-0.39675916994874E+0_wp,-0.39566141463345E+0_wp, &
      &-0.40205336237269E+0_wp,-7.70569833745062E-2_wp,-6.93201005060526E-2_wp, &
      & 0.16161438132489E+0_wp,-8.53298056844963E-2_wp,-7.28764536006634E-2_wp, &
      & 0.17055907439056E+0_wp, 0.17070980395417E+0_wp,-7.04067235322931E-2_wp, &
      &-7.59964390504423E-2_wp, 0.16964839407562E+0_wp,-6.20585040406314E-2_wp, &
      &-7.29250708812517E-2_wp, 8.64174740477892E-2_wp, 0.10406729706244E+0_wp, &
      & 8.62223301936871E-2_wp, 0.11318668876957E+0_wp, 8.39374413133319E-2_wp, &
      & 9.08063716407997E-2_wp, 9.33750029719196E-2_wp, 9.02414021711847E-2_wp, &
      & 8.52471706343666E-2_wp, 9.46559327232836E-2_wp, 8.25241730550529E-2_wp, &
      & 0.10241788528707E+0_wp, 0.10484272561566E+0_wp, 0.10417532504838E+0_wp, &
      & 8.96531455310284E-2_wp, 9.68902639794006E-2_wp]

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TCoulomb) :: coulomb

   integer :: ii, jj
   logical :: fail
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp), allocatable :: gradient(:, :)
   real(wp), allocatable :: jmat(:, :)
   real(wp), allocatable :: shift(:)
   real(wp), allocatable :: djdr(:, :, :)
   real(wp), allocatable :: djdtr(:, :)
   real(wp), allocatable :: djdL(:, :, :)

   call init(env)
   call init(mol, at, xyz, lattice=lattice)
   call init(coulomb, env, mol, tolerance=1.0e-8_wp)

   call check(error, coulomb%rCutoff, 18.11939328_wp, thr=thr)
   call check(error, coulomb%gCutoff, 1.4680064_wp, thr=thr)
   call check(error, coulomb%alpha, 0.2097152_wp, thr=thr)

   allocate(shift(nat))
   allocate(jMat(nat, nat))
   allocate(gradient(3, nat))
   allocate(djdr(3, nat, nat))
   allocate(djdtr(3, nat))
   allocate(djdL(3, 3, nat))

   call coulomb%getCoulombMatrix(mol, jmat)

   do ii = 1, nat
      call check(error, jmat(ii,ii), -0.22693841283719_wp, thr=thr)
      do jj = 1, ii-1
         call check(error, jmat(jj,ii), jmat(ii,jj), thr=thr)
      end do
   end do
   call check(error, jmat( 1, 3),-0.11073168961581E-01_wp, thr=thr)
   call check(error, jmat( 2, 3),-0.59052463388022E-01_wp, thr=thr)
   call check(error, jmat( 4, 6),-0.21264008600309E-01_wp, thr=thr)
   call check(error, jmat( 3,10),-0.24000444326437E-02_wp, thr=thr)
   call check(error, jmat( 6, 3),-0.36037976935086E-01_wp, thr=thr)
   call check(error, jmat( 7,18), 0.45116733777831E-01_wp, thr=thr)
   call check(error, jmat(12,20),-0.47663492130087E-01_wp, thr=thr)

   shift(:) = matmul(jmat, charges)
   energy = 0.5_wp*dot_product(charges, shift)
   call check(error, energy,-0.19871381077095_wp, thr=thr)

   call coulomb%getCoulombDerivs(mol, charges, djdr, djdtr, djdL)
   call contract(djdr, charges, gradient)
   call contract(djdL, charges, sigma)

   ! check numerical gradient
   do ii = 1, nat, 5
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call check(error, gradient(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   ! check numerical strain derivatives
   eps = unity
   do ii = 1, 3
      do jj = 1, ii
         eps(jj, ii) = eps(jj, ii) + step
         mol%lattice(:, :) = matmul(eps, lattice)
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%lattice(:, :) = matmul(eps, lattice)
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) + step

         call check(error, sigma(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

end subroutine test_coulomb_point_pbc3d


subroutine test_coulomb_gfn1_cluster(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_la, only : contract
   use xtb_type_environment
   use xtb_coulomb_klopmanohno
   use xtb_type_molecule
   type(error_type), allocatable, intent(out) :: error
   real(wp), parameter :: thr = 1.0e-10_wp
   real(wp), parameter :: thr2 = 1.0e-7_wp
   integer, parameter :: nat = 24
   integer, parameter :: at(nat) = [6,7,6,7,6,6,6,8,7,6,8,7,6,6, &
      &                             1,1,1,1,1,1,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[ 2.02799738646442_wp,  0.09231312124713_wp, -0.14310895950963_wp, &
      &  4.75011007621000_wp,  0.02373496014051_wp, -0.14324124033844_wp, &
      &  6.33434307654413_wp,  2.07098865582721_wp, -0.14235306905930_wp, &
      &  8.72860718071825_wp,  1.38002919517619_wp, -0.14265542523943_wp, &
      &  8.65318821103610_wp, -1.19324866489847_wp, -0.14231527453678_wp, &
      &  6.23857175648671_wp, -2.08353643730276_wp, -0.14218299370797_wp, &
      &  5.63266886875962_wp, -4.69950321056008_wp, -0.13940509630299_wp, &
      &  3.44931709749015_wp, -5.48092386085491_wp, -0.14318454855466_wp, &
      &  7.77508917214346_wp, -6.24427872938674_wp, -0.13107140408805_wp, &
      & 10.30229550927022_wp, -5.39739796609292_wp, -0.13672168520430_wp, &
      & 12.07410272485492_wp, -6.91573621641911_wp, -0.13666499342053_wp, &
      & 10.70038521493902_wp, -2.79078533715849_wp, -0.14148379504141_wp, &
      & 13.24597858727017_wp, -1.76969072232377_wp, -0.14218299370797_wp, &
      &  7.40891694074004_wp, -8.95905928176407_wp, -0.11636933482904_wp, &
      &  1.38702118184179_wp,  2.05575746325296_wp, -0.14178615122154_wp, &
      &  1.34622199478497_wp, -0.86356704498496_wp,  1.55590600570783_wp, &
      &  1.34624089204623_wp, -0.86133716815647_wp, -1.84340893849267_wp, &
      &  5.65596919189118_wp,  4.00172183859480_wp, -0.14131371969009_wp, &
      & 14.67430918222276_wp, -3.26230980007732_wp, -0.14344911021228_wp, &
      & 13.50897177220290_wp, -0.60815166181684_wp,  1.54898960808727_wp, &
      & 13.50780014200488_wp, -0.60614855212345_wp, -1.83214617078268_wp, &
      &  5.41408424778406_wp, -9.49239668625902_wp, -0.11022772492007_wp, &
      &  8.31919801555568_wp, -9.74947502841788_wp,  1.56539243085954_wp, &
      &  8.31511620712388_wp, -9.76854236502758_wp, -1.79108242206824_wp],&
      &  shape(xyz))
   integer, parameter :: nsh = 2*nat
   integer, parameter :: nshell(4) = [2, 2, 2, 2]
   real(wp), parameter :: charges(nat) = [&
      &-0.05476196_wp, 0.00248420_wp, 0.08809027_wp,-0.28313223_wp, 0.12260575_wp, &
      &-0.04041664_wp, 0.26374283_wp,-0.43589031_wp,-0.11704515_wp, 0.31376962_wp, &
      &-0.44064154_wp,-0.07607650_wp,-0.05087676_wp,-0.04160170_wp, 0.06582772_wp, &
      & 0.08200852_wp, 0.08188349_wp, 0.05475583_wp, 0.09364446_wp, 0.07149007_wp, &
      & 0.07154523_wp, 0.09011030_wp, 0.06924195_wp, 0.06924253_wp]
   real(wp), parameter :: shellCharges(nsh) = [&
      & 1.01616655_wp,-0.95681690_wp, 0.55966943_wp,-0.72093054_wp, &
      & 1.00614202_wp,-0.77968123_wp, 0.44540376_wp,-0.84062833_wp, &
      & 0.99973101_wp,-0.75581436_wp, 0.97432321_wp,-0.97083994_wp, &
      & 1.02254933_wp,-0.56987953_wp, 0.27277623_wp,-0.83701638_wp, &
      & 0.55069736_wp,-0.81347267_wp, 1.04846089_wp,-0.50836208_wp, &
      & 0.27378768_wp,-0.83862010_wp, 0.56409485_wp,-0.80234455_wp, &
      & 1.02001741_wp,-0.94938645_wp, 1.02076323_wp,-0.94757830_wp, &
      & 0.05218447_wp,-0.01315042_wp, 0.07996410_wp,-0.01800416_wp, &
      & 0.07997803_wp,-0.01800646_wp, 0.05094088_wp,-0.01121930_wp, &
      & 0.08384962_wp,-0.01890929_wp, 0.06080732_wp,-0.01556872_wp, &
      & 0.06083555_wp,-0.01557668_wp, 0.07824990_wp,-0.01849406_wp, &
      & 0.06637552_wp,-0.01700968_wp, 0.06654954_wp,-0.01700776_wp]
   real(wp), parameter :: atomicHardness(1, 4) = reshape([&
      & 0.479988_wp, 0.476106_wp, 0.583349_wp, 0.470099_wp], [1, 4])
   real(wp), parameter :: shellHardness(2, 4) = reshape([&
      & 0.479988_wp, 0.457372_wp, 0.476106_wp, 0.491108_wp, &
      & 0.583349_wp, 0.605202_wp, 0.470099_wp, 0.470099_wp], [2, 4])

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TKlopmanOhno) :: coulomb

   integer :: ii, jj
   logical :: fail
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp), allocatable :: gradient(:, :)
   real(wp), allocatable :: jmat(:, :)
   real(wp), allocatable :: shift(:)
   real(wp), allocatable :: djdr(:, :, :)
   real(wp), allocatable :: djdtr(:, :)
   real(wp), allocatable :: djdL(:, :, :)

   call init(env)
   call init(mol, at, xyz)
   call init(coulomb, env, mol, gamAverage%harmonic, atomicHardness, 2.0_wp)

   allocate(shift(nat))
   allocate(jMat(nat, nat))
   allocate(gradient(3, nat))
   allocate(djdr(3, nat, nat))
   allocate(djdtr(3, nat))
   allocate(djdL(3, 3, nat))

   call coulomb%getCoulombMatrix(mol, jmat)

   do ii = 1, nat
      call check(error, jmat(ii,ii), atomicHardness(1, mol%id(ii)), thr=thr)
      do jj = 1, ii-1
         call check(error, jmat(jj,ii), jmat(ii,jj), thr=thr)
      end do
   end do

   call check(error, jmat( 1, 3), 0.19316589570303_wp, thr=thr)
   call check(error, jmat( 2, 3), 0.30046155646117_wp, thr=thr)
   call check(error, jmat( 4, 6), 0.21047957782049_wp, thr=thr)
   call check(error, jmat( 3,10), 0.11481217301814_wp, thr=thr)
   call check(error, jmat( 6, 3), 0.21511721746087_wp, thr=thr)
   call check(error, jmat( 7,18), 0.11170281004947_wp, thr=thr)
   call check(error, jmat(12,20), 0.22373063680923_wp, thr=thr)

   shift(:) = matmul(jmat, charges)
   energy = 0.5_wp*dot_product(charges, shift)
   call check(error, energy, 0.74487442026461E-01_wp, thr=thr)

   call coulomb%getCoulombDerivs(mol, charges, djdr, djdtr, djdL)
   call contract(djdr, charges, gradient)
   call contract(djdL, charges, sigma)

   ! check numerical gradient
   do ii = 1, nat
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call check(error, gradient(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   ! check numerical strain derivatives
   ! Excluded trace due to large values
   eps = unity
   do ii = 1, 3
      do jj = 1, ii-1
         eps(jj, ii) = eps(jj, ii) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) + step

         call check(error, sigma(jj, ii), (er - el)*step2, thr=thr2)
      end do
   end do

   ! Reset
   deallocate(shift, jMat, djdr, djdtr, djdL)
   allocate(shift(nsh))
   allocate(jMat(nsh, nsh))
   allocate(djdr(3, nat, nsh))
   allocate(djdtr(3, nsh))
   allocate(djdL(3, 3, nsh))

   ! Second run now for shellwise electrostatics
   call init(mol, at, xyz)
   call init(coulomb, env, mol, gamAverage%harmonic, shellHardness, 2.0_wp, &
      & nshell=nshell)

   call coulomb%getCoulombMatrix(mol, jmat)

   do ii = 1, nsh, 2
      jj = ii/2+1
      call check(error, jmat(ii,ii), shellHardness(1, mol%id(jj)), thr=thr)
      call check(error, jmat(ii+1,ii+1), shellHardness(2, mol%id(jj)), thr=thr)
   end do
   do ii = 1, nsh
      do jj = 1, ii-1
         call check(error, jmat(jj,ii), jmat(ii,jj), thr=thr)
      end do
   end do

   call check(error, jmat( 1, 3), 0.29122779054758_wp, thr=thr)
   call check(error, jmat( 2, 3), 0.28857037810257_wp, thr=thr)
   call check(error, jmat( 4, 6), 0.29936008553464_wp, thr=thr)
   call check(error, jmat( 3,10), 0.21662900538106_wp, thr=thr)
   call check(error, jmat( 6, 3), 0.29754583720422_wp, thr=thr)
   call check(error, jmat( 7,18), 0.12567214529183_wp, thr=thr)
   call check(error, jmat(12,20), 0.17601945603238_wp, thr=thr)

   shift(:) = matmul(jmat, shellCharges)
   energy = 0.5_wp*dot_product(shellCharges, shift)
   call check(error, energy, 0.15670711466457_wp, thr=thr)

   call coulomb%getCoulombDerivs(mol, shellCharges, djdr, djdtr, djdL)
   call contract(djdr, shellCharges, gradient)
   call contract(djdL, shellCharges, sigma)

   ! check numerical gradient
   do ii = 1, nat
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, shellCharges)
         er = 0.5_wp*dot_product(shellCharges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, shellCharges)
         el = 0.5_wp*dot_product(shellCharges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call check(error, gradient(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   ! check numerical strain derivatives
   eps = unity
   do ii = 1, 3
      do jj = 1, ii
         eps(jj, ii) = eps(jj, ii) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, shellCharges)
         er = 0.5_wp*dot_product(shellCharges, shift)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, shellCharges)
         el = 0.5_wp*dot_product(shellCharges, shift)

         eps(jj, ii) = eps(jj, ii) + step

         call check(error, sigma(jj, ii), (er - el)*step2, thr=thr2)
      end do
   end do

end subroutine test_coulomb_gfn1_cluster

subroutine test_coulomb_gfn1_pbc3d(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi, sqrtpi
   use xtb_mctc_convert
   use xtb_mctc_la, only : contract
   use xtb_type_environment
   use xtb_coulomb_klopmanohno
   use xtb_type_molecule
   type(error_type), allocatable, intent(out) :: error

   real(wp),parameter :: thr = 1.0e-9_wp
   integer, parameter :: nat = 32
   integer, parameter :: at(nat) = [8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[4.5853168464880421_wp,  4.9392326929575878_wp,  4.1894081210748118_wp,  &
      & 5.8862267152491423_wp,  1.1425258978245871_wp,  6.5015058204768126_wp,  &
      & 1.4284279616220412_wp,  4.6017511285540875_wp,  3.1465884436348119_wp,  &
      & 2.3323404704521411_wp,  1.4471801154820869_wp,  0.7121932185858125_wp,  &
      & 4.7333543155561415_wp,  2.7747291872305868_wp,  3.1951352976178122_wp,  &
      & 5.6617754419101418_wp,  2.2133191164485870_wp,  2.1235404838618126_wp,  &
      & 5.3107598618381422_wp,  2.6902988056185868_wp,  0.7384466319968125_wp,  &
      & 4.4947071071761426_wp,  3.9530790692635867_wp,  0.6801776747598124_wp,  &
      & 4.8005171923760424_wp,  4.9185102874975870_wp,  1.8186363449528122_wp,  &
      & 4.6951362687070421_wp,  4.2781752812835867_wp,  3.1816411821728123_wp,  &
      & 1.3838574419160412_wp,  5.2817805008910863_wp,  4.1482702947948136_wp,  &
      & 1.0268974195990415_wp,  4.7234752637800881_wp,  5.4989995400388123_wp,  &
      & 2.0852659694760409_wp,  5.0956317453800875_wp,  6.5351699846458127_wp,  &
      & 2.3344644666691412_wp, -0.0736561690909131_wp,  6.4245628001158135_wp,  &
      & 2.4894017448231409_wp,  0.6213510313930869_wp,  5.0967297417158131_wp,  &
      & 1.5745272273791413_wp,  0.1243470825760870_wp,  3.9731040773988129_wp,  &
      & 5.8221065925130420_wp,  5.3013563342055878_wp,  1.7264876737078123_wp,  &
      & 3.4487807319551416_wp,  3.6355832152975864_wp,  0.7429568016758125_wp,  &
      & 4.8499393376520423_wp,  3.4713855169305874_wp,  6.4691872586348129_wp,  &
      & 0.2495364434351412_wp,  2.4795455690160870_wp,  2.1043557230378123_wp,  &
      & 5.6691068338331423_wp,  1.1234174220755870_wp,  2.1414388326468128_wp,  &
      & 3.7072009289431418_wp,  2.4357632918535872_wp,  3.0094700999208119_wp,  &
      & 4.1414520030430415_wp,  5.7877262477775879_wp,  1.7803680119358125_wp,  &
      & 5.0142851411171421_wp,  2.4165926460955873_wp,  4.1857610486448129_wp,  &
      & 3.0280930003030413_wp,  4.6201081184690871_wp,  6.2533190952188136_wp,  &
      & 0.5863628696651412_wp,  0.5757236365910867_wp,  4.1021714214668128_wp,  &
      & 2.3776130524831411_wp,  1.6969724987740866_wp,  5.2327688986668139_wp,  &
      & 1.9486148363011413_wp,  0.4390675147070869_wp,  2.9999022491838123_wp,  &
      & 3.5312997625581413_wp,  0.4467415528495868_wp,  4.8114121395028135_wp,  &
      & 6.5089895990100421_wp,  5.2100409408535882_wp,  6.0066553789008132_wp,  &
      & 0.9001165013630412_wp,  3.6420787128610868_wp,  5.4413106648508132_wp,  &
      & 1.6012116650460413_wp,  5.6845471271780879_wp,  0.7675566847298124_wp], &
      & shape(xyz)) * aatoau
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[6.4411018522600_wp,    0.0492571261505_wp,    0.2192046129910_wp,  &
      & 0.0462076831749_wp,    6.6435057067500_wp,    0.1670513770770_wp,  &
      & 0.2262248220170_wp,   -0.9573234940220_wp,    6.7608039126200_wp], &
      & shape(lattice)) * aatoau
   real(wp),parameter :: charges(nat) = [&
      &-0.40074825616499E+0_wp,-0.39675916994874E+0_wp,-0.39566141463345E+0_wp, &
      &-0.40205336237269E+0_wp,-7.70569833745062E-2_wp,-6.93201005060526E-2_wp, &
      & 0.16161438132489E+0_wp,-8.53298056844963E-2_wp,-7.28764536006634E-2_wp, &
      & 0.17055907439056E+0_wp, 0.17070980395417E+0_wp,-7.04067235322931E-2_wp, &
      &-7.59964390504423E-2_wp, 0.16964839407562E+0_wp,-6.20585040406314E-2_wp, &
      &-7.29250708812517E-2_wp, 8.64174740477892E-2_wp, 0.10406729706244E+0_wp, &
      & 8.62223301936871E-2_wp, 0.11318668876957E+0_wp, 8.39374413133319E-2_wp, &
      & 9.08063716407997E-2_wp, 9.33750029719196E-2_wp, 9.02414021711847E-2_wp, &
      & 8.52471706343666E-2_wp, 9.46559327232836E-2_wp, 8.25241730550529E-2_wp, &
      & 0.10241788528707E+0_wp, 0.10484272561566E+0_wp, 0.10417532504838E+0_wp, &
      & 8.96531455310284E-2_wp, 9.68902639794006E-2_wp]
   real(wp),parameter :: atomicHardness(1, 3) = reshape([&
      & 0.583349_wp, 0.479988_wp, 0.470099_wp], [1, 3])

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TKlopmanOhno) :: coulomb

   integer :: ii, jj
   logical :: fail
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp), allocatable :: gradient(:, :)
   real(wp), allocatable :: jmat(:, :)
   real(wp), allocatable :: shift(:)
   real(wp), allocatable :: djdr(:, :, :)
   real(wp), allocatable :: djdtr(:, :)
   real(wp), allocatable :: djdL(:, :, :)

   call init(env)
   call init(mol, at, xyz, lattice=lattice)
   call init(coulomb, env, mol, gamAverage%harmonic, atomicHardness, 2.0_wp, &
      & tolerance=1.0e-8_wp)

   call check(error, coulomb%rCutoff, 18.11939328_wp, thr=thr)
   call check(error, coulomb%gCutoff, 1.4680064_wp, thr=thr)
   call check(error, coulomb%alpha, 0.2097152_wp, thr=thr)

   allocate(shift(nat))
   allocate(jMat(nat, nat))
   allocate(gradient(3, nat))
   allocate(djdr(3, nat, nat))
   allocate(djdtr(3, nat))
   allocate(djdL(3, 3, nat))

   call coulomb%getCoulombMatrix(mol, jmat)

   do ii = 1, nat
      do jj = 1, ii-1
         call check(error, jmat(jj,ii), jmat(ii,jj), thr=thr)
      end do
   end do
   call check(error, jmat( 1, 3),-0.26692685831907E-01_wp, thr=thr)
   call check(error, jmat( 2, 3),-0.68773950642218E-01_wp, thr=thr)
   call check(error, jmat( 4, 6),-0.38988468002079E-01_wp, thr=thr)
   call check(error, jmat( 3,10),-0.23000529589385E-01_wp, thr=thr)
   call check(error, jmat( 6, 3),-0.51594544630953E-01_wp, thr=thr)
   call check(error, jmat( 7,18), 0.41437913557408E-02_wp, thr=thr)
   call check(error, jmat(12,20),-0.64306568194767E-01_wp, thr=thr)

   shift(:) = matmul(jmat, charges)
   energy = 0.5_wp*dot_product(charges, shift)
   call check(error, energy, 0.10334568151034_wp, thr=thr)

   call coulomb%getCoulombDerivs(mol, charges, djdr, djdtr, djdL)
   call contract(djdr, charges, gradient)
   call contract(djdL, charges, sigma)

   if (allocated(error)) return

   ! check numerical gradient
   do ii = 1, nat, 5
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call check(error, gradient(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   if (allocated(error)) return

   ! check numerical strain derivatives
   eps = unity
   do ii = 1, 3
      do jj = 1, ii
         eps(jj, ii) = eps(jj, ii) + step
         mol%lattice(:, :) = matmul(eps, lattice)
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%lattice(:, :) = matmul(eps, lattice)
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) + step

         call check(error, sigma(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   ! reset and try another shape
   call init(coulomb, env, mol, gamAverage%geometric, atomicHardness, 3.0_wp, &
      & tolerance=1.0e-8_wp)

   call coulomb%getCoulombMatrix(mol, jmat)

   do ii = 1, nat
      do jj = 1, ii-1
         call check(error, jmat(jj,ii), jmat(ii,jj), thr=thr)
      end do
   end do
   call check(error, jmat( 1, 3),-0.13464382457757E-01_wp, thr=thr)
   call check(error, jmat( 2, 3),-0.60078828891465E-01_wp, thr=thr)
   call check(error, jmat( 4, 6),-0.24046964967552E-01_wp, thr=thr)
   call check(error, jmat( 3,10),-0.60407457734382E-02_wp, thr=thr)
   call check(error, jmat( 6, 3),-0.38251294880785E-01_wp, thr=thr)
   call check(error, jmat( 7,18), 0.32150302911807E-01_wp, thr=thr)
   call check(error, jmat(12,20),-0.50061175269104E-01_wp, thr=thr)

   shift(:) = matmul(jmat, charges)
   energy = 0.5_wp*dot_product(charges, shift)
   call check(error, energy, 0.86046710402599E-01_wp, thr=thr)

end subroutine test_coulomb_gfn1_pbc3d


subroutine test_coulomb_gfn2_cluster(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_la, only : contract
   use xtb_type_environment
   use xtb_coulomb_klopmanohno
   use xtb_type_molecule
   type(error_type), allocatable, intent(out) :: error
   real(wp), parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 24
   integer, parameter :: at(nat) = [6,7,6,7,6,6,6,8,7,6,8,7,6,6, &
      &                             1,1,1,1,1,1,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[ 2.02799738646442_wp,  0.09231312124713_wp, -0.14310895950963_wp, &
      &  4.75011007621000_wp,  0.02373496014051_wp, -0.14324124033844_wp, &
      &  6.33434307654413_wp,  2.07098865582721_wp, -0.14235306905930_wp, &
      &  8.72860718071825_wp,  1.38002919517619_wp, -0.14265542523943_wp, &
      &  8.65318821103610_wp, -1.19324866489847_wp, -0.14231527453678_wp, &
      &  6.23857175648671_wp, -2.08353643730276_wp, -0.14218299370797_wp, &
      &  5.63266886875962_wp, -4.69950321056008_wp, -0.13940509630299_wp, &
      &  3.44931709749015_wp, -5.48092386085491_wp, -0.14318454855466_wp, &
      &  7.77508917214346_wp, -6.24427872938674_wp, -0.13107140408805_wp, &
      & 10.30229550927022_wp, -5.39739796609292_wp, -0.13672168520430_wp, &
      & 12.07410272485492_wp, -6.91573621641911_wp, -0.13666499342053_wp, &
      & 10.70038521493902_wp, -2.79078533715849_wp, -0.14148379504141_wp, &
      & 13.24597858727017_wp, -1.76969072232377_wp, -0.14218299370797_wp, &
      &  7.40891694074004_wp, -8.95905928176407_wp, -0.11636933482904_wp, &
      &  1.38702118184179_wp,  2.05575746325296_wp, -0.14178615122154_wp, &
      &  1.34622199478497_wp, -0.86356704498496_wp,  1.55590600570783_wp, &
      &  1.34624089204623_wp, -0.86133716815647_wp, -1.84340893849267_wp, &
      &  5.65596919189118_wp,  4.00172183859480_wp, -0.14131371969009_wp, &
      & 14.67430918222276_wp, -3.26230980007732_wp, -0.14344911021228_wp, &
      & 13.50897177220290_wp, -0.60815166181684_wp,  1.54898960808727_wp, &
      & 13.50780014200488_wp, -0.60614855212345_wp, -1.83214617078268_wp, &
      &  5.41408424778406_wp, -9.49239668625902_wp, -0.11022772492007_wp, &
      &  8.31919801555568_wp, -9.74947502841788_wp,  1.56539243085954_wp, &
      &  8.31511620712388_wp, -9.76854236502758_wp, -1.79108242206824_wp],&
      &  shape(xyz))
   integer, parameter :: nsh = 38
   integer, parameter :: nshell(4) = [2, 2, 2, 1]
   real(wp), parameter :: charges(nat) = [&
      &-0.05476196_wp, 0.00248420_wp, 0.08809027_wp,-0.28313223_wp, 0.12260575_wp, &
      &-0.04041664_wp, 0.26374283_wp,-0.43589031_wp,-0.11704515_wp, 0.31376962_wp, &
      &-0.44064154_wp,-0.07607650_wp,-0.05087676_wp,-0.04160170_wp, 0.06582772_wp, &
      & 0.08200852_wp, 0.08188349_wp, 0.05475583_wp, 0.09364446_wp, 0.07149007_wp, &
      & 0.07154523_wp, 0.09011030_wp, 0.06924195_wp, 0.06924253_wp]
   real(wp), parameter :: shellCharges(nsh) = [&
      &-0.03587416_wp,-0.01858174_wp, 0.27383743_wp,-0.27841269_wp, &
      &-0.07642190_wp, 0.16034079_wp, 0.02101744_wp,-0.29972495_wp, &
      &-0.04569512_wp, 0.16484435_wp,-0.04778340_wp, 0.02157296_wp, &
      &-0.02744706_wp, 0.28860665_wp, 0.25700597_wp,-0.69772420_wp, &
      & 0.24467004_wp,-0.35271751_wp,-0.02181401_wp, 0.32593100_wp, &
      & 0.25819509_wp,-0.69903269_wp, 0.27169266_wp,-0.34626972_wp, &
      &-0.03452070_wp,-0.01338789_wp,-0.02989810_wp,-0.00748430_wp, &
      & 0.06457802_wp, 0.08293905_wp, 0.08296802_wp, 0.05698136_wp, &
      & 0.09025556_wp, 0.07152988_wp, 0.07159003_wp, 0.08590674_wp, &
      & 0.06906357_wp, 0.06926350_wp]
   real(wp), parameter :: atomicHardness(1, 4) = reshape([&
      & 0.538015_wp, 0.461493_wp, 0.451896_wp, 0.405771_wp], [1, 4])
   real(wp), parameter :: shellHardness(2, 4) = reshape([&
      & 0.538015_wp, 0.5948486449370_wp, 0.461493_wp, 0.5152519503756_wp, &
      & 0.451896_wp, 0.5195457349920_wp, 0.405771_wp, 0.0_wp], [2, 4])

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TKlopmanOhno) :: coulomb

   integer :: ii, jj
   logical :: fail
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp), allocatable :: gradient(:, :)
   real(wp), allocatable :: jmat(:, :)
   real(wp), allocatable :: shift(:)
   real(wp), allocatable :: djdr(:, :, :)
   real(wp), allocatable :: djdtr(:, :)
   real(wp), allocatable :: djdL(:, :, :)

   call init(env)
   call init(mol, at, xyz)
   call init(coulomb, env, mol, gamAverage%geometric, atomicHardness, 2.0_wp)

   allocate(shift(nat))
   allocate(jMat(nat, nat))
   allocate(gradient(3, nat))
   allocate(djdr(3, nat, nat))
   allocate(djdtr(3, nat))
   allocate(djdL(3, 3, nat))

   call coulomb%getCoulombMatrix(mol, jmat)

   do ii = 1, nat
      call check(error, jmat(ii,ii), atomicHardness(1, mol%id(ii)), thr=thr)
      do jj = 1, ii-1
         call check(error, jmat(jj,ii), jmat(ii,jj), thr=thr)
      end do
   end do

   call check(error, jmat( 1, 3), 0.19643947717683_wp, thr=thr)
   call check(error, jmat( 2, 3), 0.30530133574756_wp, thr=thr)
   call check(error, jmat( 4, 6), 0.21212302738038_wp, thr=thr)
   call check(error, jmat( 3,10), 0.11548839173240_wp, thr=thr)
   call check(error, jmat( 6, 3), 0.21966640779889_wp, thr=thr)
   call check(error, jmat( 7,18), 0.11159957571378_wp, thr=thr)
   call check(error, jmat(12,20), 0.21900194564411_wp, thr=thr)

   shift(:) = matmul(jmat, charges)
   energy = 0.5_wp*dot_product(charges, shift)
   call check(error, energy, 0.54562180505117E-01_wp, thr=thr)

   call coulomb%getCoulombDerivs(mol, charges, djdr, djdtr, djdL)
   call contract(djdr, charges, gradient)
   call contract(djdL, charges, sigma)

   ! check numerical gradient
   do ii = 1, nat
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call check(error, gradient(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   ! check numerical strain derivatives
   eps = unity
   do ii = 1, 3
      do jj = 1, ii
         eps(jj, ii) = eps(jj, ii) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) + step

         call check(error, sigma(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   ! Reset
   deallocate(shift, jMat, djdr, djdtr, djdL)
   allocate(shift(nsh))
   allocate(jMat(nsh, nsh))
   allocate(djdr(3, nat, nsh))
   allocate(djdtr(3, nsh))
   allocate(djdL(3, 3, nsh))

   ! Second run now for shellwise electrostatics
   call init(mol, at, xyz)
   call init(coulomb, env, mol, gamAverage%arithmetic, shellHardness, 2.0_wp, &
      & nshell=nshell)

   call coulomb%getCoulombMatrix(mol, jmat)

   do ii = 1, 14, 2
      jj = ii/2+1
      call check(error, jmat(ii,ii), shellHardness(1, mol%id(jj)), thr=thr)
      call check(error, jmat(ii+1,ii+1), shellHardness(2, mol%id(jj)), thr=thr)
   end do
   do ii = 1, nsh
      do jj = 1, ii-1
         call check(error, jmat(jj,ii), jmat(ii,jj), thr=thr)
      end do
   end do

   call check(error, jmat( 1, 3), 0.29593390363812_wp, thr=thr)
   call check(error, jmat( 2, 3), 0.30152111619404_wp, thr=thr)
   call check(error, jmat( 4, 6), 0.31707000700913_wp, thr=thr)
   call check(error, jmat( 3,10), 0.22194971747742_wp, thr=thr)
   call check(error, jmat( 6, 3), 0.31180455546450_wp, thr=thr)
   call check(error, jmat( 7,18), 0.12575667968267_wp, thr=thr)
   call check(error, jmat(12,20), 0.18160321973606_wp, thr=thr)

   shift(:) = matmul(jmat, shellCharges)
   energy = 0.5_wp*dot_product(shellCharges, shift)
   call check(error, energy, 0.76930095102192E-01_wp, thr=thr)

   call coulomb%getCoulombDerivs(mol, shellCharges, djdr, djdtr, djdL)
   call contract(djdr, shellCharges, gradient)
   call contract(djdL, shellCharges, sigma)

   ! check numerical gradient
   do ii = 1, nat
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, shellCharges)
         er = 0.5_wp*dot_product(shellCharges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, shellCharges)
         el = 0.5_wp*dot_product(shellCharges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call check(error, gradient(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   ! check numerical strain derivatives
   eps = unity
   do ii = 1, 3
      do jj = 1, ii
         eps(jj, ii) = eps(jj, ii) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, shellCharges)
         er = 0.5_wp*dot_product(shellCharges, shift)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, shellCharges)
         el = 0.5_wp*dot_product(shellCharges, shift)

         eps(jj, ii) = eps(jj, ii) + step

         call check(error, sigma(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

end subroutine test_coulomb_gfn2_cluster

subroutine test_coulomb_gfn2_pbc3d(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi, sqrtpi
   use xtb_mctc_convert
   use xtb_mctc_la, only : contract
   use xtb_type_environment
   use xtb_coulomb_klopmanohno
   use xtb_type_molecule
   type(error_type), allocatable, intent(out) :: error

   real(wp),parameter :: thr = 1.0e-9_wp
   integer, parameter :: nat = 32
   integer, parameter :: nsh = 48
   integer, parameter :: at(nat) = [8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[4.5853168464880421_wp,  4.9392326929575878_wp,  4.1894081210748118_wp,  &
      & 5.8862267152491423_wp,  1.1425258978245871_wp,  6.5015058204768126_wp,  &
      & 1.4284279616220412_wp,  4.6017511285540875_wp,  3.1465884436348119_wp,  &
      & 2.3323404704521411_wp,  1.4471801154820869_wp,  0.7121932185858125_wp,  &
      & 4.7333543155561415_wp,  2.7747291872305868_wp,  3.1951352976178122_wp,  &
      & 5.6617754419101418_wp,  2.2133191164485870_wp,  2.1235404838618126_wp,  &
      & 5.3107598618381422_wp,  2.6902988056185868_wp,  0.7384466319968125_wp,  &
      & 4.4947071071761426_wp,  3.9530790692635867_wp,  0.6801776747598124_wp,  &
      & 4.8005171923760424_wp,  4.9185102874975870_wp,  1.8186363449528122_wp,  &
      & 4.6951362687070421_wp,  4.2781752812835867_wp,  3.1816411821728123_wp,  &
      & 1.3838574419160412_wp,  5.2817805008910863_wp,  4.1482702947948136_wp,  &
      & 1.0268974195990415_wp,  4.7234752637800881_wp,  5.4989995400388123_wp,  &
      & 2.0852659694760409_wp,  5.0956317453800875_wp,  6.5351699846458127_wp,  &
      & 2.3344644666691412_wp, -0.0736561690909131_wp,  6.4245628001158135_wp,  &
      & 2.4894017448231409_wp,  0.6213510313930869_wp,  5.0967297417158131_wp,  &
      & 1.5745272273791413_wp,  0.1243470825760870_wp,  3.9731040773988129_wp,  &
      & 5.8221065925130420_wp,  5.3013563342055878_wp,  1.7264876737078123_wp,  &
      & 3.4487807319551416_wp,  3.6355832152975864_wp,  0.7429568016758125_wp,  &
      & 4.8499393376520423_wp,  3.4713855169305874_wp,  6.4691872586348129_wp,  &
      & 0.2495364434351412_wp,  2.4795455690160870_wp,  2.1043557230378123_wp,  &
      & 5.6691068338331423_wp,  1.1234174220755870_wp,  2.1414388326468128_wp,  &
      & 3.7072009289431418_wp,  2.4357632918535872_wp,  3.0094700999208119_wp,  &
      & 4.1414520030430415_wp,  5.7877262477775879_wp,  1.7803680119358125_wp,  &
      & 5.0142851411171421_wp,  2.4165926460955873_wp,  4.1857610486448129_wp,  &
      & 3.0280930003030413_wp,  4.6201081184690871_wp,  6.2533190952188136_wp,  &
      & 0.5863628696651412_wp,  0.5757236365910867_wp,  4.1021714214668128_wp,  &
      & 2.3776130524831411_wp,  1.6969724987740866_wp,  5.2327688986668139_wp,  &
      & 1.9486148363011413_wp,  0.4390675147070869_wp,  2.9999022491838123_wp,  &
      & 3.5312997625581413_wp,  0.4467415528495868_wp,  4.8114121395028135_wp,  &
      & 6.5089895990100421_wp,  5.2100409408535882_wp,  6.0066553789008132_wp,  &
      & 0.9001165013630412_wp,  3.6420787128610868_wp,  5.4413106648508132_wp,  &
      & 1.6012116650460413_wp,  5.6845471271780879_wp,  0.7675566847298124_wp], &
      & shape(xyz)) * aatoau
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[6.4411018522600_wp,    0.0492571261505_wp,    0.2192046129910_wp,  &
      & 0.0462076831749_wp,    6.6435057067500_wp,    0.1670513770770_wp,  &
      & 0.2262248220170_wp,   -0.9573234940220_wp,    6.7608039126200_wp], &
      & shape(lattice)) * aatoau
   real(wp),parameter :: charges(nat) = [&
      &-0.40074825616499E+0_wp,-0.39675916994874E+0_wp,-0.39566141463345E+0_wp, &
      &-0.40205336237269E+0_wp,-7.70569833745062E-2_wp,-6.93201005060526E-2_wp, &
      & 0.16161438132489E+0_wp,-8.53298056844963E-2_wp,-7.28764536006634E-2_wp, &
      & 0.17055907439056E+0_wp, 0.17070980395417E+0_wp,-7.04067235322931E-2_wp, &
      &-7.59964390504423E-2_wp, 0.16964839407562E+0_wp,-6.20585040406314E-2_wp, &
      &-7.29250708812517E-2_wp, 8.64174740477892E-2_wp, 0.10406729706244E+0_wp, &
      & 8.62223301936871E-2_wp, 0.11318668876957E+0_wp, 8.39374413133319E-2_wp, &
      & 9.08063716407997E-2_wp, 9.33750029719196E-2_wp, 9.02414021711847E-2_wp, &
      & 8.52471706343666E-2_wp, 9.46559327232836E-2_wp, 8.25241730550529E-2_wp, &
      & 0.10241788528707E+0_wp, 0.10484272561566E+0_wp, 0.10417532504838E+0_wp, &
      & 8.96531455310284E-2_wp, 9.68902639794006E-2_wp]
   real(wp),parameter :: shellCharges(nsh) = [&
      & 2.548824270679E-01_wp, -6.417655149522E-01_wp,  2.550635758392E-01_wp, &
      &-6.396203663170E-01_wp,  2.532654296214E-01_wp, -6.361388572143E-01_wp, &
      & 2.540628366400E-01_wp, -6.399929765001E-01_wp, -4.278481458258E-02_wp, &
      &-4.516926912533E-02_wp, -4.246973895241E-02_wp, -4.510205724522E-02_wp, &
      &-7.814112584798E-02_wp,  3.246268964607E-01_wp, -4.508032230865E-02_wp, &
      &-4.909932323115E-02_wp, -4.266711021723E-02_wp, -4.507074501219E-02_wp, &
      &-7.801055759726E-02_wp,  3.277862846551E-01_wp, -7.706736816894E-02_wp, &
      & 3.274667529347E-01_wp, -4.274056057957E-02_wp, -4.198972178466E-02_wp, &
      &-4.402591663970E-02_wp, -4.284660977700E-02_wp, -7.796687803440E-02_wp, &
      & 3.284463358790E-01_wp, -4.319146370456E-02_wp, -3.898718552660E-02_wp, &
      &-4.414934521005E-02_wp, -4.299559035202E-02_wp,  8.243712818772E-02_wp, &
      & 9.197499978461E-02_wp,  6.661156570939E-02_wp,  9.654423384175E-02_wp, &
      & 7.038786665539E-02_wp,  8.271657181661E-02_wp,  7.609575142603E-02_wp, &
      & 6.807174316298E-02_wp,  8.294849862734E-02_wp,  8.523505177342E-02_wp, &
      & 6.304638754319E-02_wp,  6.764681513028E-02_wp,  9.295478478980E-02_wp, &
      & 8.550915995666E-02_wp,  6.157630928522E-02_wp,  6.771601209250E-02_wp]
   integer, parameter :: nshell(3) = [2, 2, 1]
   real(wp), parameter :: shellHardness(2, 3) = reshape([&
      & 0.451896_wp, 0.5195457349920_wp, 0.538015_wp, 0.5948486449370_wp, &
      & 0.405771_wp, 0.0_wp], [2, 3])

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TKlopmanOhno) :: coulomb

   integer :: ii, jj
   logical :: fail
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp), allocatable :: gradient(:, :)
   real(wp), allocatable :: jmat(:, :)
   real(wp), allocatable :: shift(:)
   real(wp), allocatable :: djdr(:, :, :)
   real(wp), allocatable :: djdtr(:, :)
   real(wp), allocatable :: djdL(:, :, :)

   call init(env)
   call init(mol, at, xyz, lattice=lattice)
   call init(coulomb, env, mol, gamAverage%arithmetic, shellHardness, 2.0_wp, &
      & nshell=nshell)

   call check(error, coulomb%rCutoff, 18.11939328_wp, thr=thr)
   call check(error, coulomb%gCutoff, 1.4680064_wp, thr=thr)
   call check(error, coulomb%alpha, 0.2097152_wp, thr=thr)

   allocate(shift(nsh))
   allocate(jMat(nsh, nsh))
   allocate(gradient(3, nat))
   allocate(djdr(3, nat, nsh))
   allocate(djdtr(3, nsh))
   allocate(djdL(3, 3, nsh))

   call coulomb%getCoulombMatrix(mol, jmat)

   do ii = 1, nsh
      do jj = 1, ii-1
         call check(error, jmat(jj,ii), jmat(ii,jj), thr=thr)
      end do
   end do

   call check(error, jmat( 1, 3), -0.65948223431173E-01_wp, thr=thr)
   call check(error, jmat( 2, 3), -0.63531253417773E-01_wp, thr=thr)
   call check(error, jmat( 4, 6), -0.71249446497053E-01_wp, thr=thr)
   call check(error, jmat( 3,10), -0.47826859938370E-01_wp, thr=thr)
   call check(error, jmat( 6, 3), -0.72960918626276E-01_wp, thr=thr)
   call check(error, jmat( 7,18), -0.57502477782259E-01_wp, thr=thr)
   call check(error, jmat(12,20), -0.13191548797843E-01_wp, thr=thr)

   shift(:) = matmul(jmat, shellCharges)
   energy = 0.5_wp*dot_product(shellCharges, shift)
   call check(error, energy, 0.93929971020624E-01_wp, thr=thr)

   call coulomb%getCoulombDerivs(mol, shellCharges, djdr, djdtr, djdL)
   call contract(djdr, shellCharges, gradient)
   call contract(djdL, shellCharges, sigma)

   if (allocated(error)) return

   ! check numerical gradient
   do ii = 1, nat, 5
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, shellCharges)
         er = 0.5_wp*dot_product(shellCharges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, shellCharges)
         el = 0.5_wp*dot_product(shellCharges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call check(error, gradient(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   if (allocated(error)) return

   ! check numerical strain derivatives
   eps = unity
   do ii = 1, 3
      do jj = 1, ii
         eps(jj, ii) = eps(jj, ii) + step
         mol%lattice(:, :) = matmul(eps, lattice)
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, shellCharges)
         er = 0.5_wp*dot_product(shellCharges, shift)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%lattice(:, :) = matmul(eps, lattice)
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, shellCharges)
         el = 0.5_wp*dot_product(shellCharges, shift)

         eps(jj, ii) = eps(jj, ii) + step

         call check(error, sigma(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

end subroutine test_coulomb_gfn2_pbc3d


subroutine test_coulomb_gaussian_cluster(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : sqrtpi
   use xtb_mctc_la, only : contract
   use xtb_coulomb_gaussian
   use xtb_type_environment
   use xtb_type_molecule
   type(error_type), allocatable, intent(out) :: error
   real(wp), parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 24
   integer, parameter :: at(nat) = [6,7,6,7,6,6,6,8,7,6,8,7,6,6, &
      &                             1,1,1,1,1,1,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[ 2.02799738646442_wp,  0.09231312124713_wp, -0.14310895950963_wp, &
      &  4.75011007621000_wp,  0.02373496014051_wp, -0.14324124033844_wp, &
      &  6.33434307654413_wp,  2.07098865582721_wp, -0.14235306905930_wp, &
      &  8.72860718071825_wp,  1.38002919517619_wp, -0.14265542523943_wp, &
      &  8.65318821103610_wp, -1.19324866489847_wp, -0.14231527453678_wp, &
      &  6.23857175648671_wp, -2.08353643730276_wp, -0.14218299370797_wp, &
      &  5.63266886875962_wp, -4.69950321056008_wp, -0.13940509630299_wp, &
      &  3.44931709749015_wp, -5.48092386085491_wp, -0.14318454855466_wp, &
      &  7.77508917214346_wp, -6.24427872938674_wp, -0.13107140408805_wp, &
      & 10.30229550927022_wp, -5.39739796609292_wp, -0.13672168520430_wp, &
      & 12.07410272485492_wp, -6.91573621641911_wp, -0.13666499342053_wp, &
      & 10.70038521493902_wp, -2.79078533715849_wp, -0.14148379504141_wp, &
      & 13.24597858727017_wp, -1.76969072232377_wp, -0.14218299370797_wp, &
      &  7.40891694074004_wp, -8.95905928176407_wp, -0.11636933482904_wp, &
      &  1.38702118184179_wp,  2.05575746325296_wp, -0.14178615122154_wp, &
      &  1.34622199478497_wp, -0.86356704498496_wp,  1.55590600570783_wp, &
      &  1.34624089204623_wp, -0.86133716815647_wp, -1.84340893849267_wp, &
      &  5.65596919189118_wp,  4.00172183859480_wp, -0.14131371969009_wp, &
      & 14.67430918222276_wp, -3.26230980007732_wp, -0.14344911021228_wp, &
      & 13.50897177220290_wp, -0.60815166181684_wp,  1.54898960808727_wp, &
      & 13.50780014200488_wp, -0.60614855212345_wp, -1.83214617078268_wp, &
      &  5.41408424778406_wp, -9.49239668625902_wp, -0.11022772492007_wp, &
      &  8.31919801555568_wp, -9.74947502841788_wp,  1.56539243085954_wp, &
      &  8.31511620712388_wp, -9.76854236502758_wp, -1.79108242206824_wp],&
      &  shape(xyz))
   real(wp), parameter :: charges(nat) = [&
      &-0.05476196_wp, 0.00248420_wp, 0.08809027_wp,-0.28313223_wp, 0.12260575_wp, &
      &-0.04041664_wp, 0.26374283_wp,-0.43589031_wp,-0.11704515_wp, 0.31376962_wp, &
      &-0.44064154_wp,-0.07607650_wp,-0.05087676_wp,-0.04160170_wp, 0.06582772_wp, &
      & 0.08200852_wp, 0.08188349_wp, 0.05475583_wp, 0.09364446_wp, 0.07149007_wp, &
      & 0.07154523_wp, 0.09011030_wp, 0.06924195_wp, 0.06924253_wp]
   real(wp), parameter :: rad(1, 4) = reshape([&
      & 1.88862966_wp, 1.32250290_wp, 1.23166285_wp, 0.55159092_wp], [1, 4])

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TGaussianSmeared) :: coulomb

   integer :: ii, jj
   logical :: fail
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp), allocatable :: gradient(:, :)
   real(wp), allocatable :: jmat(:, :)
   real(wp), allocatable :: shift(:)
   real(wp), allocatable :: djdr(:, :, :)
   real(wp), allocatable :: djdtr(:, :)
   real(wp), allocatable :: djdL(:, :, :)

   call init(env)
   call init(mol, at, xyz)
   call init(coulomb, env, mol, rad)

   allocate(shift(nat))
   allocate(jMat(nat, nat))
   allocate(gradient(3, nat))
   allocate(djdr(3, nat, nat))
   allocate(djdtr(3, nat))
   allocate(djdL(3, 3, nat))

   call coulomb%getCoulombMatrix(mol, jmat)

   do ii = 1, nat
      call check(error, jmat(ii,ii), sqrt(2.0_wp)/sqrtpi/rad(1, mol%id(ii)), thr=thr)
      do jj = 1, ii-1
         call check(error, jmat(jj,ii), jmat(ii,jj), thr=thr)
      end do
   end do

   call check(error, jmat( 1, 3), 0.20845479968556_wp, thr=thr)
   call check(error, jmat( 2, 3), 0.34290971840222_wp, thr=thr)
   call check(error, jmat( 4, 6), 0.23234316136565_wp, thr=thr)
   call check(error, jmat( 3,10), 0.11824383006732_wp, thr=thr)
   call check(error, jmat( 6, 3), 0.23395190326154_wp, thr=thr)
   call check(error, jmat( 7,18), 0.11492593106728_wp, thr=thr)
   call check(error, jmat(12,20), 0.25389462706068_wp, thr=thr)

   shift(:) = matmul(jmat, charges)
   energy = 0.5_wp*dot_product(charges, shift)
   call check(error, energy, 0.10159878036650_wp, thr=thr)

   call coulomb%getCoulombDerivs(mol, charges, djdr, djdtr, djdL)
   call contract(djdr, charges, gradient)
   call contract(djdL, charges, sigma)

   ! check numerical gradient
   do ii = 1, nat
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call check(error, gradient(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   ! check numerical strain derivatives
   eps = unity
   do ii = 1, 3
      do jj = 1, ii
         eps(jj, ii) = eps(jj, ii) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) + step

         call check(error, sigma(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

end subroutine test_coulomb_gaussian_cluster


subroutine test_coulomb_gaussian_pbc3d(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi, sqrtpi
   use xtb_mctc_convert
   use xtb_mctc_la, only : contract
   use xtb_coulomb_gaussian
   use xtb_type_environment
   use xtb_type_molecule
   type(error_type), allocatable, intent(out) :: error

   real(wp),parameter :: thr = 1.0e-9_wp
   integer, parameter :: nat = 32
   integer, parameter :: at(nat) = [8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[4.5853168464880421_wp,  4.9392326929575878_wp,  4.1894081210748118_wp,  &
      & 5.8862267152491423_wp,  1.1425258978245871_wp,  6.5015058204768126_wp,  &
      & 1.4284279616220412_wp,  4.6017511285540875_wp,  3.1465884436348119_wp,  &
      & 2.3323404704521411_wp,  1.4471801154820869_wp,  0.7121932185858125_wp,  &
      & 4.7333543155561415_wp,  2.7747291872305868_wp,  3.1951352976178122_wp,  &
      & 5.6617754419101418_wp,  2.2133191164485870_wp,  2.1235404838618126_wp,  &
      & 5.3107598618381422_wp,  2.6902988056185868_wp,  0.7384466319968125_wp,  &
      & 4.4947071071761426_wp,  3.9530790692635867_wp,  0.6801776747598124_wp,  &
      & 4.8005171923760424_wp,  4.9185102874975870_wp,  1.8186363449528122_wp,  &
      & 4.6951362687070421_wp,  4.2781752812835867_wp,  3.1816411821728123_wp,  &
      & 1.3838574419160412_wp,  5.2817805008910863_wp,  4.1482702947948136_wp,  &
      & 1.0268974195990415_wp,  4.7234752637800881_wp,  5.4989995400388123_wp,  &
      & 2.0852659694760409_wp,  5.0956317453800875_wp,  6.5351699846458127_wp,  &
      & 2.3344644666691412_wp, -0.0736561690909131_wp,  6.4245628001158135_wp,  &
      & 2.4894017448231409_wp,  0.6213510313930869_wp,  5.0967297417158131_wp,  &
      & 1.5745272273791413_wp,  0.1243470825760870_wp,  3.9731040773988129_wp,  &
      & 5.8221065925130420_wp,  5.3013563342055878_wp,  1.7264876737078123_wp,  &
      & 3.4487807319551416_wp,  3.6355832152975864_wp,  0.7429568016758125_wp,  &
      & 4.8499393376520423_wp,  3.4713855169305874_wp,  6.4691872586348129_wp,  &
      & 0.2495364434351412_wp,  2.4795455690160870_wp,  2.1043557230378123_wp,  &
      & 5.6691068338331423_wp,  1.1234174220755870_wp,  2.1414388326468128_wp,  &
      & 3.7072009289431418_wp,  2.4357632918535872_wp,  3.0094700999208119_wp,  &
      & 4.1414520030430415_wp,  5.7877262477775879_wp,  1.7803680119358125_wp,  &
      & 5.0142851411171421_wp,  2.4165926460955873_wp,  4.1857610486448129_wp,  &
      & 3.0280930003030413_wp,  4.6201081184690871_wp,  6.2533190952188136_wp,  &
      & 0.5863628696651412_wp,  0.5757236365910867_wp,  4.1021714214668128_wp,  &
      & 2.3776130524831411_wp,  1.6969724987740866_wp,  5.2327688986668139_wp,  &
      & 1.9486148363011413_wp,  0.4390675147070869_wp,  2.9999022491838123_wp,  &
      & 3.5312997625581413_wp,  0.4467415528495868_wp,  4.8114121395028135_wp,  &
      & 6.5089895990100421_wp,  5.2100409408535882_wp,  6.0066553789008132_wp,  &
      & 0.9001165013630412_wp,  3.6420787128610868_wp,  5.4413106648508132_wp,  &
      & 1.6012116650460413_wp,  5.6845471271780879_wp,  0.7675566847298124_wp], &
      & shape(xyz)) * aatoau
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[6.4411018522600_wp,    0.0492571261505_wp,    0.2192046129910_wp,  &
      & 0.0462076831749_wp,    6.6435057067500_wp,    0.1670513770770_wp,  &
      & 0.2262248220170_wp,   -0.9573234940220_wp,    6.7608039126200_wp], &
      & shape(lattice)) * aatoau
   real(wp),parameter :: charges(nat) = [&
      &-0.40074825616499E+0_wp,-0.39675916994874E+0_wp,-0.39566141463345E+0_wp, &
      &-0.40205336237269E+0_wp,-7.70569833745062E-2_wp,-6.93201005060526E-2_wp, &
      & 0.16161438132489E+0_wp,-8.53298056844963E-2_wp,-7.28764536006634E-2_wp, &
      & 0.17055907439056E+0_wp, 0.17070980395417E+0_wp,-7.04067235322931E-2_wp, &
      &-7.59964390504423E-2_wp, 0.16964839407562E+0_wp,-6.20585040406314E-2_wp, &
      &-7.29250708812517E-2_wp, 8.64174740477892E-2_wp, 0.10406729706244E+0_wp, &
      & 8.62223301936871E-2_wp, 0.11318668876957E+0_wp, 8.39374413133319E-2_wp, &
      & 9.08063716407997E-2_wp, 9.33750029719196E-2_wp, 9.02414021711847E-2_wp, &
      & 8.52471706343666E-2_wp, 9.46559327232836E-2_wp, 8.25241730550529E-2_wp, &
      & 0.10241788528707E+0_wp, 0.10484272561566E+0_wp, 0.10417532504838E+0_wp, &
      & 8.96531455310284E-2_wp, 9.68902639794006E-2_wp]
   real(wp), parameter :: rad(1, 3) = reshape([&
      & 1.23166285_wp, 1.88862966_wp, 0.55159092_wp], [1, 3])

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TGaussianSmeared) :: coulomb

   integer :: ii, jj
   logical :: fail
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp), allocatable :: gradient(:, :)
   real(wp), allocatable :: jmat(:, :)
   real(wp), allocatable :: shift(:)
   real(wp), allocatable :: djdr(:, :, :)
   real(wp), allocatable :: djdtr(:, :)
   real(wp), allocatable :: djdL(:, :, :)

   call init(env)
   call init(mol, at, xyz, lattice=lattice)
   call init(coulomb, env, mol, rad, tolerance=1.0e-8_wp)

   call check(error, coulomb%rCutoff, 18.11939328_wp, thr=thr)
   call check(error, coulomb%gCutoff, 1.4680064_wp, thr=thr)
   call check(error, coulomb%alpha, 0.2097152_wp, thr=thr)

   allocate(shift(nat))
   allocate(jMat(nat, nat))
   allocate(gradient(3, nat))
   allocate(djdr(3, nat, nat))
   allocate(djdtr(3, nat))
   allocate(djdL(3, 3, nat))

   call coulomb%getCoulombMatrix(mol, jmat)

   do ii = 1, nat
      do jj = 1, ii-1
         call check(error, jmat(jj,ii), jmat(ii,jj), thr=thr)
      end do
   end do
   call check(error, jmat( 1, 3),-0.11073244345567E-01_wp, thr=thr)
   call check(error, jmat( 2, 3),-0.59052463388029E-01_wp, thr=thr)
   call check(error, jmat( 4, 6),-0.21273998758336E-01_wp, thr=thr)
   call check(error, jmat( 3,10),-0.24405713033154E-02_wp, thr=thr)
   call check(error, jmat( 6, 3),-0.36042624892326E-01_wp, thr=thr)
   call check(error, jmat( 7,18), 0.43960420991048E-01_wp, thr=thr)
   call check(error, jmat(12,20),-0.47663494495435E-01_wp, thr=thr)

   shift(:) = matmul(jmat, charges)
   energy = 0.5_wp*dot_product(charges, shift)
   call check(error, energy, 0.17261986036367_wp, thr=thr)

   call coulomb%getCoulombDerivs(mol, charges, djdr, djdtr, djdL)
   call contract(djdr, charges, gradient)
   call contract(djdL, charges, sigma)

   ! check numerical gradient
   do ii = 1, nat, 5
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call check(error, gradient(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   ! check numerical strain derivatives
   eps = unity
   do ii = 1, 3
      do jj = 1, ii
         eps(jj, ii) = eps(jj, ii) + step
         mol%lattice(:, :) = matmul(eps, lattice)
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         er = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%lattice(:, :) = matmul(eps, lattice)
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call coulomb%update(env, mol)
         call coulomb%getCoulombMatrix(mol, jmat)
         shift(:) = matmul(jmat, charges)
         el = 0.5_wp*dot_product(charges, shift)

         eps(jj, ii) = eps(jj, ii) + step

         call check(error, sigma(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

end subroutine test_coulomb_gaussian_pbc3d

end module test_coulomb
