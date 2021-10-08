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

module test_coordinationnumber
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
   implicit none
   private

   public :: collect_coordinationnumber

contains

!> Collect all exported unit tests
subroutine collect_coordinationnumber(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("lp-pbc3d", test_ncoord_pbc3d_latticepoints), &
      new_unittest("lp-pbc3d", test_ncoord_pbc3d_neighbourlist) &
      ]

end subroutine collect_coordinationnumber


subroutine test_ncoord_pbc3d_latticepoints(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_disp_coordinationnumber, only : getCoordinationNumber, cnType
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init, len
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_type_latticepoint, only : TLatticePoint, init
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-10_wp
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

   type(TMolecule) :: mol
   type(TEnvironment) :: env
   type(TLatticePoint) :: latp

   real(wp), allocatable :: trans(:, :)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)

   call init(env)

   call init(mol, at, xyz, lattice=lattice)

   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat))

   call init(latp, env, mol, 40.0_wp)  ! Fixed cutoff
   call latp%getLatticePoints(trans, 40.0_wp)

   call getCoordinationNumber(mol, trans, 40.0_wp, cnType%exp, &
      & cn, dcndr, dcndL)

   call check(error, cn( 1), 1.0644889257925_wp, thr=thr)
   call check(error, cn( 2), 1.0631767469669_wp, thr=thr)
   call check(error, cn(11), 3.1006970695816_wp, thr=thr)
   call check(error, cn(18), 1.0075732171735_wp, thr=thr)
   call check(error, cn(19), 1.0059732047433_wp, thr=thr)
   call check(error, cn(26), 1.0076129647756_wp, thr=thr)
   call check(error, norm2(dcndr), 1.0111771304250_wp, thr=thr)

   call getCoordinationNumber(mol, trans, 40.0_wp, cnType%erf, &
      & cn, dcndr, dcndL)

   call check(error, cn( 4), 1.0017771035280_wp, thr=thr)
   call check(error, cn( 8), 3.9823122257845_wp, thr=thr)
   call check(error, cn(12), 3.9825267700108_wp, thr=thr)
   call check(error, cn(14), 2.9944763199487_wp, thr=thr)
   call check(error, cn(15), 3.9804643916712_wp, thr=thr)
   call check(error, cn(16), 3.9809919526339_wp, thr=thr)
   call check(error, cn(31), 0.9937992456497_wp, thr=thr)
   call check(error, norm2(dcndr), 0.56981134655575_wp, thr=thr)

   call getCoordinationNumber(mol, trans, 40.0_wp, cnType%cov, &
      & cn, dcndr, dcndL)

   call check(error, cn( 9), 3.8048129507609_wp, thr=thr)
   call check(error, cn(13), 3.8069332740701_wp, thr=thr)
   call check(error, cn(20), 0.9239504636072_wp, thr=thr)
   call check(error, cn(21), 0.9245161490576_wp, thr=thr)
   call check(error, cn(22), 0.9236911432582_wp, thr=thr)
   call check(error, cn(30), 0.9241833212868_wp, thr=thr)
   call check(error, norm2(dcndr), 0.53663107191832_wp, thr=thr)

   call getCoordinationNumber(mol, trans, 40.0_wp, cnType%gfn, &
      & cn, dcndr, dcndL)

   call check(error, cn( 3), 1.1819234421527_wp, thr=thr)
   call check(error, cn( 5), 4.2453004956734_wp, thr=thr)
   call check(error, cn(10), 3.3142686980438_wp, thr=thr)
   call check(error, cn(23), 1.0257176118018_wp, thr=thr)
   call check(error, cn(24), 1.0254266219049_wp, thr=thr)
   call check(error, cn(32), 1.0249293241088_wp, thr=thr)
   call check(error, norm2(dcndr), 2.9659913317090_wp, thr=thr)

end subroutine test_ncoord_pbc3d_latticepoints


subroutine test_ncoord_pbc3d_neighbourlist(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_disp_coordinationnumber, only : getCoordinationNumber, cnType
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init, len
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_type_latticepoint, only : TLatticePoint, init
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-10_wp
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

   type(TMolecule) :: mol
   type(TEnvironment) :: env
   type(TLatticePoint) :: latp
   type(TNeighbourlist) :: neighList

   integer, allocatable :: neighs(:)
   real(wp), allocatable :: trans(:, :)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)

   call init(env)

   call init(mol, at, xyz, lattice=lattice)

   allocate(neighs(nat))
   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat))

   call init(latp, env, mol, 40.0_wp)  ! Fixed cutoff
   call latp%getLatticePoints(trans, 40.0_wp)
   call init(neighList, len(mol))
   call neighList%generate(env, mol%xyz, 40.0_wp, trans, .false.)

   call neighlist%getNeighs(neighs, 40.0_wp)
   call getCoordinationNumber(mol, neighs, neighList, cnType%exp, &
      & cn, dcndr, dcndL)

   call check(error, cn( 1), 1.0644889257925_wp, thr=thr)
   call check(error, cn( 2), 1.0631767469669_wp, thr=thr)
   call check(error, cn(11), 3.1006970695816_wp, thr=thr)
   call check(error, cn(18), 1.0075732171735_wp, thr=thr)
   call check(error, cn(19), 1.0059732047433_wp, thr=thr)
   call check(error, cn(26), 1.0076129647756_wp, thr=thr)
   call check(error, norm2(dcndr), 1.0111771304250_wp, thr=thr)

   call getCoordinationNumber(mol, neighs, neighList, cnType%erf, &
      & cn, dcndr, dcndL)

   call check(error, cn( 4), 1.0017771035280_wp, thr=thr)
   call check(error, cn( 8), 3.9823122257845_wp, thr=thr)
   call check(error, cn(12), 3.9825267700108_wp, thr=thr)
   call check(error, cn(14), 2.9944763199487_wp, thr=thr)
   call check(error, cn(15), 3.9804643916712_wp, thr=thr)
   call check(error, cn(16), 3.9809919526339_wp, thr=thr)
   call check(error, cn(31), 0.9937992456497_wp, thr=thr)
   call check(error, norm2(dcndr), 0.56981134655575_wp, thr=thr)

   call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
      & cn, dcndr, dcndL)

   call check(error, cn( 9), 3.8048129507609_wp, thr=thr)
   call check(error, cn(13), 3.8069332740701_wp, thr=thr)
   call check(error, cn(20), 0.9239504636072_wp, thr=thr)
   call check(error, cn(21), 0.9245161490576_wp, thr=thr)
   call check(error, cn(22), 0.9236911432582_wp, thr=thr)
   call check(error, cn(30), 0.9241833212868_wp, thr=thr)
   call check(error, norm2(dcndr), 0.53663107191832_wp, thr=thr)

   call getCoordinationNumber(mol, neighs, neighList, cnType%gfn, &
      & cn, dcndr, dcndL)

   call check(error, cn( 3), 1.1819234421527_wp, thr=thr)
   call check(error, cn( 5), 4.2453004956734_wp, thr=thr)
   call check(error, cn(10), 3.3142686980438_wp, thr=thr)
   call check(error, cn(23), 1.0257176118018_wp, thr=thr)
   call check(error, cn(24), 1.0254266219049_wp, thr=thr)
   call check(error, cn(32), 1.0249293241088_wp, thr=thr)
   call check(error, norm2(dcndr), 2.9659913317090_wp, thr=thr)

end subroutine test_ncoord_pbc3d_neighbourlist

end module test_coordinationnumber
