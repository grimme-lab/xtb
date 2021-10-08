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

module test_repulsion
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
   implicit none
   private

   public :: collect_repulsion

contains

!> Collect all exported unit tests
subroutine collect_repulsion(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("cluster", test_repulsion_cluster), &
      new_unittest("pbc3d", test_repulsion_pbc3d) &
      ]

end subroutine collect_repulsion


subroutine test_repulsion_cluster(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : sqrtpi
   use xtb_mctc_la, only : contract
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_xtb_data
   use xtb_xtb_gfn2
   use xtb_xtb_repulsion
   type(error_type), allocatable, intent(out) :: error
   real(wp), parameter :: thr = 1.0e-9_wp
   real(wp), parameter :: thr2 = 1.0e-8_wp
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

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRepulsionData) :: rep

   integer :: ii, jj
   logical :: fail
   real(wp) :: energy, sigma(3, 3), sdum(3, 3)
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: gradient(:, :)
   real(wp), allocatable :: gdum(:, :)

   call init(env)
   call init(mol, at, xyz)
   call initGFN2(rep)

   allocate(gradient(3, nat))
   allocate(gdum(3, nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call repulsionEnGrad(mol, rep, trans, 40.0_wp, energy, gradient, sigma)

   call check(error, energy, 0.49222837261241_wp, thr=thr)

   ! check numerical gradient
   do ii = 1, nat
      do jj = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call repulsionEnGrad(mol, rep, trans, 40.0_wp, er, gdum, sdum)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call repulsionEnGrad(mol, rep, trans, 40.0_wp, el, gdum, sdum)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call check(error, gradient(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   ! check numerical strain derivatives
   eps = unity
   do ii = 1, 3
      do jj = 1, ii
         er = 0.0_wp
         el = 0.0_wp
         eps(jj, ii) = eps(jj, ii) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call repulsionEnGrad(mol, rep, trans, 40.0_wp, er, gdum, sdum)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call repulsionEnGrad(mol, rep, trans, 40.0_wp, el, gdum, sdum)

         eps(jj, ii) = eps(jj, ii) + step
         call check(error, sigma(jj, ii), (er - el)*step2, thr=thr2)
      end do
   end do

end subroutine test_repulsion_cluster

subroutine test_repulsion_pbc3d(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi, sqrtpi
   use xtb_mctc_convert
   use xtb_mctc_la, only : contract
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_latticepoint, only : TLatticePoint, init
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_xtb_data
   use xtb_xtb_gfn1
   use xtb_xtb_repulsion
   type(error_type), allocatable, intent(out) :: error

   real(wp),parameter :: thr = 1.0e-10_wp
   real(wp),parameter :: thr2 = 1.0e-8_wp
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

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRepulsionData) :: rep
   type(TLatticePoint) :: latp
   type(TNeighbourlist) :: neighList

   integer :: ii, jj
   logical :: fail
   real(wp) :: energy, sigma(3, 3), sdum(3, 3)
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   integer, allocatable :: neighs(:)
   real(wp), allocatable :: gradient(:, :), gdum(:, :)
   real(wp), allocatable :: trans(:, :), trans2(:, :)

   call init(env)
   call init(mol, at, xyz, lattice=lattice)
   call initGFN1(rep)
   call init(latp, env, mol, 40.0_wp)  ! Fixed cutoff
   call init(neighList, len(mol))

   allocate(neighs(nat))
   allocate(gradient(3, nat))
   allocate(gdum(3, nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call latp%getLatticePoints(trans, 40.0_wp)
   call neighList%generate(env, mol%xyz, 40.0_wp, trans, .false.)
   call neighlist%getNeighs(neighs, 40.0_wp)

   call repulsionEnGrad(mol, rep, trans, 40.0_wp, energy, gradient, sigma)

   call check(error, energy, 0.58472111816325_wp, thr=thr)

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call repulsionEnGrad(mol, rep, neighs, neighList, energy, gradient, sigma)

   call check(error, energy, 0.58472111816325_wp, thr=thr)

   ! check numerical gradient
   do ii = 1, nat, 5
      do jj = modulo(ii, 2)+1, 3, 2
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call repulsionEnGrad(mol, rep, trans, 40.0_wp, er, gdum, sdum)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call repulsionEnGrad(mol, rep, trans, 40.0_wp, el, gdum, sdum)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call check(error, gradient(jj, ii), (er - el)*step2, thr=thr)
      end do
   end do

   ! check numerical strain derivatives
   trans2 = trans
   eps = unity
   do ii = 1, 3
      do jj = 1, ii
         er = 0.0_wp
         el = 0.0_wp
         eps(jj, ii) = eps(jj, ii) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         trans(:, :) = matmul(eps, trans2)
         call mol%update
         call repulsionEnGrad(mol, rep, trans, 40.0_wp, er, gdum, sdum)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         trans(:, :) = matmul(eps, trans2)
         call mol%update
         call repulsionEnGrad(mol, rep, trans, 40.0_wp, el, gdum, sdum)

         eps(jj, ii) = eps(jj, ii) + step
         call check(error, sigma(jj, ii), (er - el)*step2, thr=thr2)
      end do
   end do

end subroutine test_repulsion_pbc3d

end module test_repulsion
