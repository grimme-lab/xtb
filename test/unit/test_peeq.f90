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

module test_peeq
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check, test_failed
   implicit none
   private

   public :: collect_peeq

contains

!> Collect all exported unit tests
subroutine collect_peeq(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("sp", test_peeq_sp), &
      new_unittest("api", test_peeq_api), &
      new_unittest("srb", test_peeq_api_srb) &
      ]

end subroutine collect_peeq


subroutine test_peeq_sp(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_mctc_systools

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_wavefunction
   use xtb_type_basisset
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment
   use xtb_solv_gbsa

   use xtb_setparam, only : set

   use xtb_pbc_tools
   use xtb_basis
   use xtb_peeq
   use xtb_readparam
   use xtb_paramset

   use xtb_xtb_data
   use xtb_xtb_gfn0

   type(error_type), allocatable, intent(out) :: error

   real(wp),parameter :: thr = 1.0e-7_wp
   ! SiO2
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,14,8]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[4.51379029_wp,    8.29324256_wp,   -0.60000000_wp, &
      & 0.56666667_wp,    0.56666667_wp,   -0.60000000_wp, &
      & 0.96666667_wp,    0.96666667_wp,    1.20000001_wp],&
      & shape(xyz))
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[6.04712363_wp,    0.00000000_wp,    0.00000000_wp,    &
      & 0.00000000_wp,    9.82657590_wp,    0.00000000_wp,    &
      & 0.00000000_wp,    0.00000000_wp,   13.22808290_wp],   &
      & shape(lattice))
   real(wp),parameter :: et = 300.0_wp
   integer, parameter :: prlevel = 2
   logical, parameter :: lgrad = .true.
   logical, parameter :: ccm = .true.
   logical, parameter :: restart = .false.
   real(wp),parameter :: acc = 1.0_wp

   type(TMolecule)     :: mol
   type(TEnvironment)  :: env
   type(TWavefunction) :: wfn
   type(TBasisset)     :: basis
   type(scc_results)     :: res
   type(TxTBData) :: xtbData
   type(TBorn), allocatable :: gbsa

   real(wp)              :: energy
   real(wp)              :: hl_gap
   real(wp), allocatable :: gradient(:,:)
   real(wp)              :: sigma(3,3)

   integer, parameter    :: wsc_rep(3) = [1,1,1] ! FIXME

   character(len=*), parameter   :: p_fnv_gfn0 = 'param_gfn0-xtb.txt'
   character(len=:), allocatable :: fnv
   type(TxTBParameter) :: globpar
   integer  :: ipar
   logical  :: exist

   logical  :: okbas,diff

   set%gfn_method = 0
   call init(env)

   call init(mol, at, xyz, lattice=lattice)

   allocate( gradient(3,mol%n) )
   energy = 0.0_wp
   gradient = 0.0_wp
   sigma = 0.0_wp

   ! give an optional summary on the geometry used
   call print_pbcsum(stdout,mol)

   ! we will try an internal parameter file first to avoid IO
   call use_parameterset(p_fnv_gfn0,globpar,xtbData,exist)
   ! no luck, we have to fire up some IO to get our parameters
   if (.not.exist) then
      ! let's check if we can find the parameter file
      call rdpath(env%xtbpath,p_fnv_gfn0,fnv,exist)
      ! maybe the user provides a local parameter file, this was always
      ! an option in `xtb', so we will give it a try
      if (.not.exist) fnv = p_fnv_gfn0
      call open_file(ipar,fnv,'r')
      if (ipar.eq.-1) then
         call test_failed(error, "Parameter file '"//fnv//"' not found!")
         return
      endif
      call readParam(env,ipar,globpar,xtbData,.true.)
      call close_file(ipar)
   endif

   call newBasisset(xtbData,mol%n,mol%at,basis,okbas)

   call wfn%allocate(mol%n,basis%nshell,basis%nao)
   wfn%nel = idint(sum(mol%z)) - mol%chrg
   wfn%nopen = 0
   wfn%q = mol%chrg / real(mol%n,wp)

   call mctc_mute

   call peeq(env,mol,wfn,basis,xtbData,gbsa,hl_gap,et,prlevel,lgrad,.true.,acc, &
      &      energy,gradient,sigma,res)

   call check_(error, energy,-7.3570001012578_wp, thr=thr)
   call check_(error, hl_gap, 2.1170422163611_wp, thr=1.0e-4_wp)
   call check_(error, norm2(gradient),4.6531382100157_wp, thr=thr)
   call check_(error, gradient(1,1),0.0171831109909_wp, thr=thr)
   call check_(error, gradient(3,2),3.1002099423654_wp, thr=thr)

   call check_(error, sigma(1,1),-0.48669769112351_wp, thr=thr)
   call check_(error, sigma(2,1),-0.33613693973170_wp, thr=thr)
   call check_(error, sigma(2,3), -1.3504590911322_wp, thr=thr)

   call check_(error, res%e_elec,-8.303090012925_wp, thr=thr)
   call check_(error, res%e_es,  -0.126523327000_wp, thr=thr)
   call check_(error, res%e_disp,-0.004101334765_wp, thr=thr)
   call check_(error, res%e_rep,  1.076714573432_wp, thr=thr)
   call check_(error, res%e_xb,   0.000000000000_wp, thr=thr)

   ! reset for reevaluatuation without CCM
   energy = 0.0_wp
   gradient = 0.0_wp
   sigma = 0.0_wp

   call peeq(env,mol,wfn,basis,xtbData,gbsa,hl_gap,et,prlevel,lgrad,.false.,acc, &
      &      energy,gradient,sigma,res)

   call check_(error, energy,-7.3514275392244_wp, thr=thr)
   call check_(error, hl_gap, 2.2269146636198_wp, thr=1.0e-4_wp)
   call check_(error, norm2(gradient),4.6440176418778_wp, thr=thr)
   call check_(error, gradient(1,1),0.0078432890861_wp, thr=thr)
   call check_(error, gradient(3,2),3.0969842115048_wp, thr=thr)

   call check_(error, sigma(1,1),-0.51138347441108_wp, thr=thr)
   call check_(error, sigma(2,1),-0.31805950090099_wp, thr=thr)
   call check_(error, sigma(2,3), -1.3396127658627_wp, thr=thr)

   call check_(error, res%e_elec,-8.297517450892_wp, thr=thr)
   call check_(error, res%e_es,  -0.126523327000_wp, thr=thr)
   call check_(error, res%e_disp,-0.004101334765_wp, thr=thr)
   call check_(error, res%e_rep,  1.076714573432_wp, thr=thr)
   call check_(error, res%e_xb,   0.000000000000_wp, thr=thr)

end subroutine test_peeq_sp

subroutine test_peeq_api(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_data
   use xtb_type_param
   use xtb_type_environment

   use xtb_pbc_tools

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   type(error_type), allocatable, intent(out) :: error

   real(wp),parameter :: thr = 1.0e-10_wp
   ! CaF2
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [9,9,20]
   real(wp),parameter :: abc(3,nat) = reshape(&
      &[0.25_wp, 0.25_wp, 0.25_wp, &
      & 0.75_wp, 0.75_wp, 0.75_wp, &
      & 0.00_wp, 0.00_wp, 0.00_wp], shape(abc))
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[5.9598811567890_wp,      2.1071361905157_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      6.3214085715472_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      0.0000000000000_wp,      7.2993338808807_wp],   &
      & shape(lattice))
   type(peeq_options),parameter :: opt = peeq_options( &
      &  prlevel = 2, ccm = .true., acc = 1.0_wp, etemp = 300.0_wp, grad = .true. )

   type(TMolecule)    :: mol
   type(TEnvironment) :: env
   type(TRestart) :: chk
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   logical :: failed
   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp) :: gradlatt(3,3)
   real(wp) :: sigma(3,3),inv_lat(3,3)
   real(wp),allocatable :: gradient(:,:)
   real(wp),allocatable :: xyz(:,:)

   ! setup the environment variables
   call init(env)

   allocate(xyz(3, nat))
   call coord_trafo(nat,lattice,abc,xyz)
   call init(mol, at, xyz, lattice=lattice)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp
   gradlatt = 0.0_wp

   call mctc_mute

   call newXTBCalculator(env, mol, calc, method=0)
   call env%check(failed)
   if (failed) then
      call test_failed(error, "Setup failed")
      return
   end if
   call newWavefunction(env, mol, calc, chk)

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)
   inv_lat = mat_inv_3x3(mol%lattice)
   call sigma_to_latgrad(sigma,inv_lat,gradlatt)

   call check_(error, hl_gap, 4.9685235017906_wp, thr=thr)
   call check_(error, energy,-8.4863996084661_wp, thr=thr)
   call check_(error, norm2(gradient),0.33507483384363E-04_wp, thr=thr)
   call check_(error, norm2(gradlatt),0.33064163041261E-02_wp, thr=thr)

end subroutine test_peeq_api

subroutine test_peeq_api_srb(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_mctc_convert
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_data
   use xtb_type_param
   use xtb_type_environment

   use xtb_pbc_tools

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   type(error_type), allocatable, intent(out) :: error

   real(wp),parameter :: thr = 1.0e-9_wp
   ! CaF2
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
   type(peeq_options),parameter :: opt = peeq_options( &
      &  prlevel = 2, ccm = .true., acc = 1.0_wp, etemp = 300.0_wp, grad = .true. )

   type(TMolecule)    :: mol
   type(TEnvironment) :: env
   type(TRestart) :: chk
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   logical :: failed
   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp) :: gradlatt(3,3)
   real(wp) :: sigma(3,3),inv_lat(3,3)
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call init(env)

   call init(mol, at, xyz, lattice=lattice)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp
   gradlatt = 0.0_wp

   call mctc_mute

   call newXTBCalculator(env, mol, calc, method=0)
   call env%check(failed)
   if (failed) then
      call test_failed(error, "Setup failed")
      return
   end if
   call newWavefunction(env, mol, calc, chk)

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)
   inv_lat = mat_inv_3x3(mol%lattice)
   call sigma_to_latgrad(sigma,inv_lat,gradlatt)

   call check_(error, hl_gap, 3.2452476555284_wp, thr=thr)
   call check_(error, energy,-47.338099017467_wp, thr=thr)
   call check_(error, norm2(gradient),0.60674690405096E-01_wp, thr=thr)
   call check_(error, norm2(gradlatt),0.25787608030965E-01_wp, thr=thr)

   call check_(error, gradient(1, 3), 0.15527599551578E-02_wp, thr=thr)
   call check_(error, gradient(2,11), 0.14689883866846E-01_wp, thr=thr)
   call check_(error, gradient(1, 6),-0.18911924607680E-02_wp, thr=thr)
   call check_(error, gradient(3, 5), 0.15935221260425E-02_wp, thr=thr)
   call check_(error, gradient(1, 8), 0.25516596273654E-02_wp, thr=thr)

   call check_(error, gradlatt(1,3), 0.24324041066017E-03_wp, thr=thr)
   call check_(error, gradlatt(2,2), 0.11757572269151E-01_wp, thr=thr)
   call check_(error, gradlatt(1,2), 0.60490544331116E-03_wp, thr=thr)

end subroutine test_peeq_api_srb

end module test_peeq
