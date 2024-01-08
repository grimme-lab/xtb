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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

module test_ptb

   use mctc_env, only: wp
   use mctc_io, only: structure_type, new

   use xtb_type_molecule, only: TMolecule, assignment(=)
   use xtb_test_molstock, only: getMolecule
   use xtb_features, only : get_xtb_feature

   use testdrive, only: new_unittest, unittest_type, error_type, check_ => check, test_failed, skip_test

   implicit none
   private

   real(wp), parameter :: thr = 1.0e-7_wp
   real(wp), parameter :: thr2 = 1.0e-5_wp
   real(wp), parameter :: thr3 = 1.0e-4_wp
   real(wp), parameter :: thr4 = 1.0e-3_wp
   !> Reference implementation of PTB used partially single precision
   !> Numerical derivatives enhance small errors, which is why we need a
   !> higher threshold
   real(wp), parameter :: thr_alpha = 1.0e-1_wp

   public :: collect_ptb

contains

!> Collect all exported unit tests
   subroutine collect_ptb(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
#if WITH_TBLITE
                  new_unittest("basis", test_ptb_basis), &
                  new_unittest("eeq", test_ptb_eeq), &
                  new_unittest("overlap", test_ptb_overlap), &
                  new_unittest("overlap_h0", test_ptb_overlap_h0), &
                  new_unittest("overlap_sx", test_ptb_overlap_SX), &
                  new_unittest("v_ecp", test_ptb_V_ECP), &
                  new_unittest("selfenergies", test_ptb_selfenergies), &
                  new_unittest("hamiltonian_h0", test_ptb_hamiltonian_h0), &
                  new_unittest("v_xc", test_ptb_V_XC), &
                  new_unittest("hubbard", test_ptb_hubbard), &
                  new_unittest("coulomb_pot", test_ptb_coulomb_potential), &
                  new_unittest("plus_U_pot", test_ptb_plus_U_potential), &
                  new_unittest("mb16-43-01", test_ptb_mb16_43_01), &
                  new_unittest("mb16-43-01_charged", test_ptb_mb16_43_01_charged), &
                  new_unittest("mb16-43-01_efield", test_ptb_mb16_43_01_efield), &
                  new_unittest("dipole_moment", test_ptb_dipmom_caffeine), &
                  new_unittest("polarizability", test_ptb_polarizability) &
#else
                  new_unittest("ptb_not_present", test_ptb_not_present) &
#endif
                  ]

   end subroutine collect_ptb

#if WITH_TBLITE
   subroutine test_ptb_basis(error)
      use xtb_ptb_vdzp, only: add_vDZP_basis
      use tblite_basis_type, only: basis_type

      type(error_type), allocatable, intent(out) :: error
      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> Basis set type
      type(basis_type) :: bas

      integer :: nsh_exp = 11
      integer :: nao_exp = 23
      integer :: maxl_exp = 2

      call getMolecule(struc, "h2o")
      mol = struc

      call add_vDZP_basis(mol, bas)

      call check_(error, bas%nsh, nsh_exp, "Number of shells not matching to expected value.")
      call check_(error, bas%nao, nao_exp, "Number of AOs not matching to expected value.")
      call check_(error, bas%maxl, maxl_exp, "Maximum angular momentum not matching to expected value.")
      call check_(error, bas%cgto(1, 1)%alpha(1), 3.54364182058582_wp, thr=thr)
      call check_(error, bas%cgto(1, 2)%alpha(1), 81.8867808750392_wp, thr=thr)

   end subroutine test_ptb_basis

   subroutine test_ptb_eeq(error)
      use xtb_ptb_param, only: ptbGlobals, initPTB
      use xtb_ptb_ncoord, only: ncoord_erf
      use multicharge_model, only: new_mchrg_model, mchrg_model_type
      use xtb_ptb_data, only: TPTBData

      type(error_type), allocatable, intent(out) :: error
      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData
      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> EEQ Model
      type(mchrg_model_type) :: eeqmodel
      real(wp), allocatable :: cn_eeq(:)
      integer :: i
      !> EEQ charges
      real(wp), allocatable :: q_eeq(:)

      real(wp), parameter :: q_exp(16) = [ &
      &  0.191209985_wp, &
      &  0.093681828_wp, &
      &  0.013831373_wp, &
      &  0.086328769_wp, &
      &  0.023833233_wp, &
      & -0.024182298_wp, &
      &  0.127341995_wp, &
      & -0.033180774_wp, &
      &  0.195514665_wp, &
      &  0.033261333_wp, &
      &  0.004616166_wp, &
      &  0.478869532_wp, &
      & -0.573638388_wp, &
      & -0.379500359_wp, &
      &  0.060803901_wp, &
      & -0.298790960_wp]

      call getMolecule(struc, "mindless01")
      mol = struc
      allocate (ptbData)
      call initPTB(ptbData, mol%num)
      allocate (cn_eeq(mol%nat))
      call ncoord_erf(mol, ptbGlobals%kerfcn_eeq, 25.0_wp, cn_eeq)
      call new_mchrg_model(eeqmodel, chi=ptbData%eeq%chi, &
      & rad=ptbData%eeq%alp, eta=ptbData%eeq%gam, kcn=ptbData%eeq%cnf)
      allocate (q_eeq(mol%nat))
      call eeqmodel%solve(mol, cn_eeq, qvec=q_eeq)

      do i = 1, mol%nat
         call check_(error, q_eeq(i), q_exp(i), thr=thr, &
         & message="EEQ charge for an atom not matching to expected value.")
      end do

   end subroutine test_ptb_eeq

   subroutine test_ptb_overlap(error)
      use xtb_ptb_integrals, only: get_integrals
      use xtb_ptb_integral_types, only: aux_integral_type, new_aux_integral
      use xtb_ptb_vdzp, only: add_vDZP_basis
      use tblite_basis_type, only: basis_type, get_cutoff
      use tblite_integral_type, only: integral_type, new_integral
      use tblite_adjlist, only: adjacency_list, new_adjacency_list
      use tblite_cutoff, only: get_lattice_points

      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> PTB vDZP basis set
      type(basis_type) :: bas
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> Integral type
      type(integral_type) :: ints
      !> Adjacency list
      type(adjacency_list) :: list
      !> Error type
      type(error_type), allocatable, intent(out) :: error
      !> (Scaled) overlap matrix
      character(len=:), allocatable :: message
      real(wp), parameter :: overlap_exp(6) = [ &
      & 0.93209460_wp, & ! 1,2
      & 0.35489609_wp, & ! 1,3
      & 0.65682608_wp, & ! 2,3
      & 0.05627743_wp, & ! 1,15
      & -0.14217162_wp, &  ! 1,24; diffferent because of tblite ordering
      & 0.41844087_wp] ! 14,23; diffferent because of tblite ordering
      real(wp), allocatable :: lattr(:, :)
      real(wp) :: cutoff

      message = "Overlap matrix element not matching to expected value."

      call getMolecule(struc, "mgh2")
      mol = struc
      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, bas)

      !> Get the cutoff for the lattice points
      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      !> Get the adjacency list for iteration through the Hamiltonian
      call new_adjacency_list(list, mol, lattr, cutoff)

      call new_integral(ints, bas%nao)
      call get_integrals(mol, bas, lattr, list, ints%overlap)
      call check_(error, ints%overlap(1, 2), overlap_exp(1), thr=thr)
      call check_(error, ints%overlap(1, 3), overlap_exp(2), thr=thr)
      call check_(error, ints%overlap(2, 3), overlap_exp(3), thr=thr)
      call check_(error, ints%overlap(1, 15), overlap_exp(4), thr=thr)
      call check_(error, ints%overlap(1, 23), overlap_exp(5), thr=thr)
      call check_(error, ints%overlap(12, 22), overlap_exp(6), thr=thr)

   end subroutine test_ptb_overlap

   subroutine test_ptb_overlap_h0(error)
      use xtb_ptb_integrals, only: get_integrals
      use xtb_ptb_param, only: kalphah0l
      use xtb_ptb_vdzp, only: add_vDZP_basis
      use xtb_ptb_integral_types, only: aux_integral_type, new_aux_integral
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_param, only: initPTB
      use tblite_basis_type, only: basis_type, get_cutoff
      use tblite_integral_type, only: integral_type, new_integral
      use tblite_adjlist, only: adjacency_list, new_adjacency_list
      use tblite_cutoff, only: get_lattice_points

      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> PTB vDZP basis set
      type(basis_type) :: bas
      !> Integral type
      type(integral_type) :: ints
      !> Auxiliary integral type
      type(aux_integral_type) :: auxints
      !> Adjacency list
      type(adjacency_list) :: list
      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData
      !> Error type
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: message
      real(wp), parameter :: overlap_exp(6) = [ &
      & 0.95689468_wp, & ! 1,2
      & 0.39195790_wp, & ! 1,3
      & 0.62961212_wp, & ! 2,3
      & 0.03782850_wp, & ! 1,15
      &-0.13826216_wp, &  ! 1,24; diffferent because of tblite ordering
      & 0.43334922_wp] ! 14,23; diffferent because of tblite ordering
      real(wp), allocatable :: lattr(:, :)
      real(wp) :: cutoff

      message = "Scaled overlap matrix element not matching to expected value."

      call getMolecule(struc, "mgh2")
      mol = struc
      allocate (ptbData)
      call initPTB(ptbData, mol%num)
      call add_vDZP_basis(mol, id_to_atom(mol, ptbData%hamiltonian%kalphah0l), bas)

      !> Get the cutoff for the lattice points
      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      !> Get the adjacency list for iteration through the Hamiltonian
      call new_adjacency_list(list, mol, lattr, cutoff)

      !> set up the basis set for the PTB-Hamiltonian
      call new_integral(ints, bas%nao)
      call new_aux_integral(auxints, bas%nao)
      call get_integrals(mol, bas, lattr, list, auxints%overlap_h0_1)
      call check_(error, auxints%overlap_h0_1(1, 2), overlap_exp(1), thr=thr, &
      & message=message)
      call check_(error, auxints%overlap_h0_1(1, 3), overlap_exp(2), thr=thr, &
      & message=message)
      call check_(error, auxints%overlap_h0_1(2, 3), overlap_exp(3), thr=thr, &
      & message=message)
      call check_(error, auxints%overlap_h0_1(1, 15), overlap_exp(4), thr=thr, &
      & message=message)
      call check_(error, auxints%overlap_h0_1(1, 23), overlap_exp(5), thr=thr, &
      & message=message)
      call check_(error, auxints%overlap_h0_1(12, 22), overlap_exp(6), thr=thr, &
      & message=message)

   end subroutine test_ptb_overlap_h0

   subroutine test_ptb_overlap_SX(error)
      !> PTB dependencies
      use xtb_ptb_integrals, only: get_integrals
      use xtb_ptb_mmlpopanalysis, only: get_mml_overlaps
      use xtb_ptb_param, only: ptbGlobals
      use xtb_ptb_vdzp, only: add_vDZP_basis
      !> tblite dependencies
      use tblite_basis_type, only: basis_type, get_cutoff
      use tblite_integral_type, only: integral_type, new_integral
      use tblite_adjlist, only: adjacency_list, new_adjacency_list
      use tblite_cutoff, only: get_lattice_points

      !> PTB vDZP basis set
      type(basis_type) :: bas
      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> Integral type
      type(integral_type) :: ints
      !> Adjacency list
      type(adjacency_list) :: list
      !> Error type
      type(error_type), allocatable, intent(out) :: error
      !> Mulliken-Loewdin overlap matrices
      real(wp), allocatable :: overlap_sx(:, :), overlap_oneminusx(:, :)
      character(len=:), allocatable :: message
      real(wp), parameter :: overlap_oneminusx_exp(6) = [ &
      & 0.70788_wp, & ! 1,2
      & 0.16203_wp, & ! 1,3
      & 0.41532_wp, & ! 2,3
      & 0.01449_wp, & ! 1,15
      &-0.07203_wp, &  ! 1,24; diffferent because of tblite ordering
      & 0.28751_wp] ! 14,23; diffferent because of tblite ordering
      real(wp), allocatable :: lattr(:, :)
      real(wp) :: cutoff

      message = "Solved overlap matrix (^(1-x)) element not matching to expected value."

      call getMolecule(struc, "mgh2")
      mol = struc
      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, bas)
      !> Get the cutoff for the lattice points
      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      !> Get the adjacency list for iteration through the Hamiltonian
      call new_adjacency_list(list, mol, lattr, cutoff)

      call new_integral(ints, bas%nao)
      call get_integrals(mol, bas, lattr, list, ints%overlap)
      allocate (overlap_sx(bas%nao, bas%nao), overlap_oneminusx(bas%nao, bas%nao))
      call get_mml_overlaps(bas, ints%overlap, ptbGlobals%mlmix, overlap_sx, &
         & overlap_oneminusx)

      call check_(error, overlap_oneminusx(1, 2), overlap_oneminusx_exp(1), thr=thr2, &
      & message=message)
      call check_(error, overlap_oneminusx(1, 3), overlap_oneminusx_exp(2), thr=thr2, &
      & message=message)
      call check_(error, overlap_oneminusx(2, 3), overlap_oneminusx_exp(3), thr=thr2, &
      & message=message)
      call check_(error, overlap_oneminusx(1, 15), overlap_oneminusx_exp(4), thr=thr2, &
      & message=message)
      call check_(error, overlap_oneminusx(1, 23), overlap_oneminusx_exp(5), thr=thr2, &
      & message=message)
      call check_(error, overlap_oneminusx(12, 22), overlap_oneminusx_exp(6), thr=thr2, &
      & message=message)

   end subroutine test_ptb_overlap_SX

   subroutine test_ptb_V_ECP(error)
      !> tblite dependencies
      use tblite_basis_type, only: basis_type, get_cutoff
      use tblite_integral_type, only: integral_type, new_integral
      use tblite_adjlist, only: adjacency_list, new_adjacency_list
      use tblite_cutoff, only: get_lattice_points
      !> PTB dependencies
      use xtb_ptb_vdzp, only: add_vDZP_basis
      use xtb_ptb_corebasis, only: add_core_basis, get_Vecp
      use xtb_ptb_integrals, only: get_integrals
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_param, only: initPTB
      use xtb_ptb_integral_types, only: aux_integral_type, new_aux_integral

      !> PTB vDZP basis set and core basis set
      type(basis_type) :: bas, cbas
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> Integral type
      type(integral_type) :: ints
      !> Adjacency list
      type(adjacency_list) :: list
      !> Auxiliary integral type
      type(aux_integral_type) :: auxints
      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData
      !> Effective core potential
      real(wp), allocatable :: vecp(:, :)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: message
      real(wp), parameter :: vecp_ref(4) = [ &
      &  0.077719_wp, & ! 1,1 ; diffferent because of tblite ordering
      & -0.059122_wp, & ! 1,3 ; diffferent because of tblite ordering
      &  0.052775_wp, & ! 3,5 ; diffferent because of tblite ordering
      &  0.117176_wp]   ! 9,9 ; diffferent because of tblite ordering
      real(wp), parameter :: xyz(3, 2) = reshape([ &
      & 2.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp], [3, 2])
      integer, parameter :: nat = 2
      integer, parameter :: at(nat) = [5, 17]
      real(wp), allocatable :: lattr(:, :)
      real(wp) :: cutoff

      call new(mol, at, xyz)
      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, bas)

      allocate (ptbData)
      call initPTB(ptbData, mol%num)

      !> Get the cutoff for the lattice points
      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      !> Get the adjacency list for iteration through the Hamiltonian
      call new_adjacency_list(list, mol, lattr, cutoff)

      !> New integrals
      call new_integral(ints, bas%nao)
      call new_aux_integral(auxints, bas%nao)
      !> Add the core basis set to 'cbas' basis set type
      call add_core_basis(mol, ptbData%corepotential, cbas)
      !> -> for normalization factors
      call get_integrals(mol, bas, lattr, list, ints%overlap, norm=auxints%norm)
      !> V_ECP via PTB core basis
      call get_Vecp(mol, ptbData%corepotential, bas, cbas, auxints%norm, vecp)

      message = "V_ecp matrix element not matching to expected value."
      call check_(error, vecp(1, 1), vecp_ref(1), thr=thr2, &
      & message=message)
      call check_(error, vecp(1, 5), vecp_ref(2), thr=thr2, &
      & message=message)
      call check_(error, vecp(5, 8), vecp_ref(3), thr=thr2, &
      & message=message)
      call check_(error, vecp(12, 12), vecp_ref(4), thr=thr2, &
      & message=message)
   end subroutine test_ptb_V_ECP

   subroutine test_ptb_selfenergies(error)
      use tblite_basis_type, only: basis_type
      use xtb_ptb_param, only: initPTB
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_hamiltonian, only: get_selfenergy
      use xtb_ptb_vdzp, only: add_vDZP_basis

      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> Basis set type
      type(basis_type) :: bas
      !> Error type
      type(error_type), allocatable, intent(out) :: error
      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData
      !> Self-energies
      real(wp), allocatable :: selfenergies(:)
      !> Loop indices
      integer :: i
      real(wp), parameter :: cn_normal(16) = [ &
      &     0.029861202_wp, &
      &     0.000001336_wp, &
      &     0.130180690_wp, &
      &     0.000001491_wp, &
      &     0.045177881_wp, &
      &     0.033255113_wp, &
      &     0.000001968_wp, &
      &     0.092410992_wp, &
      &     0.218995650_wp, &
      &     0.043474970_wp, &
      &     0.032456272_wp, &
      &     0.113876443_wp, &
      &     0.054029961_wp, &
      &     0.186931056_wp, &
      &     0.222153392_wp, &
      &     0.154564877_wp &
      & ]

      real(wp), parameter :: cn_star(16) = [ &
      &     4.180004755_wp, &
      &     0.803964926_wp, &
      &     1.830423819_wp, &
      &     1.291529751_wp, &
      &     1.254518733_wp, &
      &     0.976661010_wp, &
      &     1.524506244_wp, &
      &     1.645124982_wp, &
      &     3.435549964_wp, &
      &     1.012204252_wp, &
      &     1.002889652_wp, &
      &     2.106643916_wp, &
      &     3.669906182_wp, &
      &     3.358227639_wp, &
      &     3.430865101_wp, &
      &     5.328552015_wp &
      & ]

      real(wp), parameter :: selfenergies_exp(69) = [ &
      &    -0.979211416_wp, &
      &    -0.322921123_wp, &
      &    -0.697551918_wp, &
      &    -1.023303502_wp, &
      &    -0.527093446_wp, &
      &    -0.514468061_wp, &
      &    -0.663960241_wp, &
      &    -0.281743767_wp, &
      &    -0.225076174_wp, &
      &    -1.383725200_wp, &
      &    -0.320247843_wp, &
      &    -0.673454332_wp, &
      &    -0.380112291_wp, &
      &    -0.359748390_wp, &
      &    -0.691229491_wp, &
      &    -0.286971656_wp, &
      &    -0.231757459_wp, &
      &    -1.205066048_wp, &
      &    -0.457780964_wp, &
      &    -0.711328979_wp, &
      &    -0.424791274_wp, &
      &    -0.315612564_wp, &
      &    -0.670965339_wp, &
      &    -0.283512863_wp, &
      &    -0.227190530_wp, &
      &    -0.704259715_wp, &
      &    -0.289469733_wp, &
      &    -0.234950020_wp, &
      &    -1.369436961_wp, &
      &    -0.319173701_wp, &
      &    -0.666351872_wp, &
      &    -0.375372836_wp, &
      &    -0.357968581_wp, &
      &    -1.200385286_wp, &
      &    -0.309750483_wp, &
      &    -0.624221535_wp, &
      &    -0.327315414_wp, &
      &    -0.350288899_wp, &
      &    -0.672137693_wp, &
      &    -0.283868580_wp, &
      &    -0.227600095_wp, &
      &    -0.672496043_wp, &
      &    -0.283796084_wp, &
      &    -0.227556008_wp, &
      &    -1.025588207_wp, &
      &    -0.388417755_wp, &
      &    -0.718646241_wp, &
      &    -0.447066846_wp, &
      &    -0.571777637_wp, &
      &    -0.794691589_wp, &
      &    -0.268730643_wp, &
      &    -0.574208389_wp, &
      &    -0.289554304_wp, &
      &    -0.194234001_wp, &
      &    -0.795098088_wp, &
      &    -0.268832436_wp, &
      &    -0.580445264_wp, &
      &    -0.290173217_wp, &
      &    -0.197401138_wp, &
      &    -1.200410535_wp, &
      &    -0.309694149_wp, &
      &    -0.624487718_wp, &
      &    -0.327397806_wp, &
      &    -0.350251414_wp, &
      &    -0.737447660_wp, &
      &    -0.274291608_wp, &
      &    -0.488598371_wp, &
      &    -0.279699675_wp, &
      &    -0.270934337_wp &
      & ]

      call getMolecule(struc, "mindless01")
      mol = struc
      allocate (ptbData)
      call initPTB(ptbData, mol%num)
      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, bas)

      !> Get self-energies
      call get_selfenergy(mol, bas, ptbData%hamiltonian, cn_normal, cn_star, selfenergies)

      do i = 1, size(selfenergies_exp)
         call check_(error, selfenergies(i), selfenergies_exp(i), thr=thr, &
         & message="Self-energies not matching to expected value.")
      end do

   end subroutine test_ptb_selfenergies

   subroutine test_ptb_hamiltonian_h0(error)
      !> tblite basis set type
      use tblite_basis_type, only: basis_type, get_cutoff
      use tblite_integral_type, only: new_integral, integral_type
      use tblite_cutoff, only: get_lattice_points
      use tblite_adjlist, only: adjacency_list, new_adjacency_list
      !> PTB core basis set generation
      use xtb_ptb_vdzp, only: add_vDZP_basis
      use xtb_ptb_integrals, only: get_integrals
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_param, only: initPTB, ptbGlobals
      use xtb_ptb_integral_types, only: aux_integral_type, new_aux_integral
      use xtb_ptb_hamiltonian, only: get_hamiltonian

      !> PTB vDZP basis set and core basis set
      type(basis_type) :: bas
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData
      !> Integral type
      type(integral_type) :: ints
      !> Auxiliary integral type
      type(aux_integral_type) :: auxints
      !> Adjacency list
      type(adjacency_list) :: list
      !> Error type
      type(error_type), allocatable, intent(out) :: error
      !> Cutoff for lattice points
      real(wp) :: cutoff
      !> Lattice points
      real(wp), allocatable :: lattr(:, :)
      character(len=:), allocatable :: message
      !> Structure
      real(wp), parameter :: xyz(3, 2) = reshape([ &
      & 2.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp], [3, 2])
      integer, parameter :: nat = 2
      integer, parameter :: at(nat) = [5, 17]
      !> Dummy vecp
      real(wp), allocatable :: vecp(:, :)

      real(wp), parameter :: h0_ref(6) = [ &
      &  -1.59330281_wp, & ! 1,1
      &  -2.24996207_wp, & ! 1,2
      &   0.34974782_wp, & ! 1,23 ; diffferent because of tblite ordering
      &   0.0_wp, & ! 7,11 ; different because of tblite ordering
      &  -1.17757007_wp, & ! 3,6 ; different because of tblite ordering
      &   0.48301561_wp]   ! 11,24 ; diffferent because of tblite ordering
      real(wp), parameter :: levels(10) = [ &
      &    -0.796651404_wp, &
      &    -0.269771638_wp, &
      &    -0.593749262_wp, &
      &    -0.292154638_wp, &
      &    -0.204518309_wp, &
      &    -1.022956257_wp, &
      &    -0.356719004_wp, &
      &    -0.736092857_wp, &
      &    -0.464644474_wp, &
      &    -0.572539094_wp]

      call new(mol, at, xyz)
      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, bas)

      allocate (ptbData)
      call initPTB(ptbData, mol%num)

      !> Get the cutoff for the lattice points
      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      !> Get the adjacency list for iteration through the Hamiltonian
      call new_adjacency_list(list, mol, lattr, cutoff)
      call new_integral(ints, bas%nao)
      call new_aux_integral(auxints, bas%nao)
      call get_integrals(mol, lattr, list, auxints%overlap_h0_1, &
         & alpha_scal=id_to_atom(mol, ptbData%hamiltonian%kalphah0l))
      allocate (vecp(bas%nao, bas%nao), source=0.0_wp)

      call get_hamiltonian(mol, list, bas, ptbData%hamiltonian, ptbData%hamiltonian%kla, auxints%overlap_h0_1, &
      & levels, ints%hamiltonian, ptbGlobals%kpol, ptbGlobals%kitr, ptbGlobals%kitocod)
      message = "H0 matrix element not matching to expected value."
      call check_(error, ints%hamiltonian(1, 1), h0_ref(1), thr=thr)
      call check_(error, ints%hamiltonian(1, 2), h0_ref(2), thr=thr)
      call check_(error, ints%hamiltonian(1, 22), h0_ref(3), thr=thr)
      call check_(error, ints%hamiltonian(13, 6), h0_ref(4), thr=thr)
      call check_(error, ints%hamiltonian(8, 5), h0_ref(5), thr=thr)
      call check_(error, ints%hamiltonian(13, 26), h0_ref(6), thr=thr)
   end subroutine test_ptb_hamiltonian_h0

   subroutine test_ptb_V_XC(error)
      !> tblite basis set type
      use tblite_basis_type, only: basis_type, get_cutoff
      use tblite_integral_type, only: integral_type, new_integral
      use tblite_adjlist, only: adjacency_list, new_adjacency_list
      use tblite_cutoff, only: get_lattice_points
      !> PTB core basis set generation
      use xtb_ptb_vdzp, only: add_vDZP_basis
      use xtb_ptb_corebasis, only: add_core_basis, get_Vecp
      !> PTB overlap matrix calculation
      use xtb_ptb_integrals, only: get_integrals
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_param, only: initPTB
      use xtb_ptb_paulixc, only: calc_Vxc_pauli
      use xtb_ptb_integral_types, only: aux_integral_type, new_aux_integral

      !> PTB vDZP basis set
      type(basis_type) :: bas
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> Integral type
      type(integral_type) :: ints
      !> Auxiliary integral type
      type(aux_integral_type) :: auxints
      !> Adjacency list
      type(adjacency_list) :: list
      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData
      !> Error type
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: message
      real(wp), parameter :: Vxc_ref(4) = [ &
      & -0.92793357_wp, & ! 1,1
      & -0.85981333_wp, & ! 1,2
      &  0.06632750_wp, & ! 1,23 ; diffferent because of tblite ordering
      &  0.00151880_wp]   ! 11,24 ; diffferent because of tblite ordering
      real(wp), parameter :: xyz(3, 2) = reshape([ &
      & 2.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp], [3, 2])
      integer, parameter :: nat = 2
      integer, parameter :: at(nat) = [5, 17]
      !> (Scaled) overlap matrix
      real(wp), allocatable :: Vxc(:, :)
      real(wp), parameter :: shellpops(10) = [ &
      &     0.562668128_wp, &
      &     0.182476968_wp, &
      &     1.714291065_wp, &
      &     0.447798837_wp, &
      &     0.583667232_wp, &
      &     1.180050163_wp, &
      &     0.466096258_wp, &
      &     3.688775422_wp, &
      &     1.071618511_wp, &
      &     0.102557478_wp]
      real(wp), parameter :: levels(10) = [ &
      &    -0.796651404_wp, &
      &    -0.269771638_wp, &
      &    -0.593749262_wp, &
      &    -0.292154638_wp, &
      &    -0.204518309_wp, &
      &    -1.022956257_wp, &
      &    -0.356719004_wp, &
      &    -0.736092857_wp, &
      &    -0.464644474_wp, &
      &    -0.572539094_wp]
      real(wp), allocatable :: lattr(:, :)
      real(wp) :: cutoff

      call new(mol, at, xyz)

      allocate (ptbData)
      call initPTB(ptbData, mol%num)

      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, id_to_atom(mol, ptbData%pauli%klalphaxc), bas)

      !> Get the cutoff for the lattice points
      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      !> Get the adjacency list for iteration through the Hamiltonian
      call new_adjacency_list(list, mol, lattr, cutoff)
      !> New integrals
      call new_integral(ints, bas%nao)
      call new_aux_integral(auxints, bas%nao)
      call get_integrals(mol, bas, lattr, list, auxints%overlap_xc)
      allocate (Vxc(bas%nao, bas%nao), source=0.0_wp)
      call calc_Vxc_pauli(mol, bas, shellpops, auxints%overlap_xc, levels, ptbData%pauli%kxc1, Vxc)

      message = "V_XC matrix element not matching to expected value."
      call check_(error, Vxc(1, 1), Vxc_ref(1), thr=thr, &
      & message=message)
      call check_(error, Vxc(1, 2), Vxc_ref(2), thr=thr, &
      & message=message)
      call check_(error, Vxc(1, 22), Vxc_ref(3), thr=thr, &
      & message=message)
      call check_(error, Vxc(13, 26), Vxc_ref(4), thr=thr, &
      & message=message)

   end subroutine test_ptb_V_XC

   subroutine test_ptb_hubbard(error)
      !> PTB overlap matrix calculation
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_param, only: initPTB
      use xtb_ptb_coulomb, only: coulomb_potential
      use xtb_ptb_vdzp, only: add_vDZP_basis
      !> tblite basis set type
      use tblite_basis_type, only: basis_type

      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> Coulomb potential
      type(coulomb_potential) :: coulomb
      !> PTB vDZP basis set
      type(basis_type) :: bas
      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData
      !> Error type
      type(error_type), allocatable, intent(out) :: error
      !> Structure type (xtb)
      type(TMolecule) :: struc
      integer :: iat, ish, ii

      real(wp), parameter :: q(16) = [ &
      &     0.191209985_wp, &
      &     0.093681828_wp, &
      &     0.013831373_wp, &
      &     0.086328769_wp, &
      &     0.023833233_wp, &
      &    -0.024182298_wp, &
      &     0.127341995_wp, &
      &    -0.033180774_wp, &
      &     0.195514665_wp, &
      &     0.033261333_wp, &
      &     0.004616166_wp, &
      &     0.478869532_wp, &
      &    -0.573638388_wp, &
      &    -0.379500359_wp, &
      &     0.060803901_wp, &
      &    -0.298790960_wp]

      real(wp), parameter :: gam_ref(69) = [ &
      &     0.572444634_wp, &
      &     0.174965574_wp, &
      &     0.278145596_wp, &
      &     0.398307229_wp, &
      &     0.140308229_wp, &
      &     0.071711282_wp, &
      &     0.699475104_wp, &
      &     0.013446193_wp, &
      &     1.442753158_wp, &
      &     0.288724569_wp, &
      &     0.503017753_wp, &
      &     0.756666283_wp, &
      &     0.576968584_wp, &
      &     1.661483681_wp, &
      &     0.698218709_wp, &
      &     0.013422041_wp, &
      &     1.440161689_wp, &
      &     0.227163646_wp, &
      &     3.918675421_wp, &
      &     0.645967076_wp, &
      &     0.871429742_wp, &
      &     1.585273418_wp, &
      &     0.679336015_wp, &
      &     0.013059054_wp, &
      &     1.401213818_wp, &
      &     0.705226516_wp, &
      &     0.013556754_wp, &
      &     1.454616149_wp, &
      &     0.285342872_wp, &
      &     0.497126139_wp, &
      &     0.747803801_wp, &
      &     0.570210818_wp, &
      &     1.642023493_wp, &
      &     0.208808961_wp, &
      &     0.276840832_wp, &
      &     0.652246415_wp, &
      &     0.299293868_wp, &
      &     0.853512542_wp, &
      &     0.689151236_wp, &
      &     0.013247734_wp, &
      &     1.421458914_wp, &
      &     0.684256722_wp, &
      &     0.013153646_wp, &
      &     1.411363380_wp, &
      &    11.195497673_wp, &
      &     0.011092823_wp, &
      &     0.414254987_wp, &
      &     0.324375967_wp, &
      &     0.072787435_wp, &
      &     0.332654683_wp, &
      &     0.007146331_wp, &
      &     0.772424282_wp, &
      &     0.025798653_wp, &
      &     1.373674042_wp, &
      &     0.351502934_wp, &
      &     0.007551243_wp, &
      &     0.816189927_wp, &
      &     0.027260407_wp, &
      &     1.451506566_wp, &
      &     0.202104465_wp, &
      &     0.267951949_wp, &
      &     0.631303903_wp, &
      &     0.289684057_wp, &
      &     0.826107721_wp, &
      &     3.596643423_wp, &
      &     0.007302541_wp, &
      &     0.215477130_wp, &
      &     0.082115458_wp, &
      &     1.009534102_wp]

      call getMolecule(struc, "mindless01")
      mol = struc
      allocate (ptbData)
      call initPTB(ptbData, mol%num)

      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, bas)
      call coulomb%init(mol, bas, q, ptbData%coulomb%shellHardnessFirstIter, &
      & ptbData%coulomb%kQHubbard, ptbData%coulomb%kOK1, ptbData%coulomb%kTO)

      !> Check the shell gammas
      ii = 0
      do iat = 1, mol%nat
         do ish = 1, bas%nsh_at(iat)
            ii = ii + 1
            call check_(error, coulomb%gam(ish, iat), gam_ref(ii), thr=thr2)
         end do
      end do

      !> Check the Hubbard matrix
      call check_(error, coulomb%hubbard(4, 5, 12, 16), 0.490990522465014_wp, thr=thr2)
      call check_(error, coulomb%hubbard(1, 2, 6, 9), 0.393374820246109_wp, thr=thr2)

   end subroutine test_ptb_hubbard

   subroutine test_ptb_coulomb_potential(error)
      !> PTB overlap matrix calculation
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_param, only: initPTB
      use xtb_ptb_coulomb, only: coulomb_potential
      use xtb_ptb_vdzp, only: add_vDZP_basis
      use xtb_ptb_guess, only: guess_qsh
      use xtb_ptb_hamiltonian, only: get_occupation
      use xtb_ptb_integrals, only: get_integrals
      !> tblite basis set type
      use tblite_basis_type, only: basis_type, get_cutoff
      use tblite_wavefunction, only: wavefunction_type, new_wavefunction
      use tblite_scf_potential, only: potential_type, new_potential, add_pot_to_h1
      use tblite_integral_type, only: integral_type, new_integral
      use tblite_adjlist, only: adjacency_list, new_adjacency_list
      use tblite_cutoff, only: get_lattice_points

      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> Coulomb potential
      type(coulomb_potential) :: coulomb
      !> PTB vDZP basis set
      type(basis_type) :: bas
      !> Wavefunction data
      type(wavefunction_type) :: wfn
      !> Potential type
      type(potential_type) :: pot
      !> Integral type
      type(integral_type) :: ints
      !> Adjacency list
      type(adjacency_list) :: list
      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData
      !> Error type
      type(error_type), allocatable, intent(out) :: error
      !> Structure type (xtb)
      type(TMolecule) :: struc
      real(wp), parameter :: q(3) = [ &
      &    0.42739211_wp, &
      &    -0.21369606_wp, &
      &    -0.21369606_wp]
      !> Conversion factor from temperature to energy
      real(wp), parameter :: kt = 3.166808578545117e-06_wp
      real(wp), parameter :: coulomb_pot_ref(4) = [ &
      &  -0.05693153_wp, & ! 1,1
      &  -0.33917531_wp, & ! 1,2
      &  -0.00539212_wp, & ! 1,21 ; diffferent because of tblite ordering
      &   0.01305793_wp]   ! 6,24 ; diffferent because of tblite ordering
      real(wp), allocatable :: lattr(:, :)
      real(wp) :: cutoff

      call getMolecule(struc, "mgh2")
      mol = struc
      allocate (ptbData)
      call initPTB(ptbData, mol%num)

      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, bas)
      !> Get the cutoff for the lattice points
      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      !> Get the adjacency list for iteration through the Hamiltonian
      call new_adjacency_list(list, mol, lattr, cutoff)
      !> New integrals
      call new_integral(ints, bas%nao)
      call get_integrals(mol, bas, lattr, list, ints%overlap)

      call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, &
         & nspin=1, kt=300.0_wp * kt)
      call new_potential(pot, mol, bas, wfn%nspin)
      call pot%reset()

      wfn%qat(:, 1) = q
      !> Project reference occupation on wavefunction and use EEQ charges as guess
      call get_occupation(mol, bas, ptbData%hamiltonian%refocc, wfn%nocc, wfn%n0at, wfn%n0sh)
      call guess_qsh(wfn, bas)

      call coulomb%init(mol, bas, wfn%qat(:, 1), ptbData%coulomb%shellHardnessFirstIter, &
      & ptbData%coulomb%kQHubbard, ptbData%coulomb%kOK1, ptbData%coulomb%kTO)
      call coulomb%update(mol, bas)
      call coulomb%get_potential(wfn, pot)
      ints%hamiltonian = 0.0_wp
      call add_pot_to_h1(bas, ints, pot, wfn%coeff)

      call check_(error, wfn%coeff(1, 1, 1), coulomb_pot_ref(1), thr=thr)
      call check_(error, wfn%coeff(1, 2, 1), coulomb_pot_ref(2), thr=thr)
      call check_(error, wfn%coeff(1, 21, 1), coulomb_pot_ref(3), thr=thr)
      call check_(error, wfn%coeff(5, 23, 1), coulomb_pot_ref(4), thr=thr)

   end subroutine test_ptb_coulomb_potential

   subroutine test_ptb_plus_U_potential(error)

      !> xtb-ptb lib
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_param, only: initPTB
      use xtb_ptb_plusu, only: plusu_potential_type
      use xtb_ptb_vdzp, only: add_vDZP_basis
      !> tblite lib
      use tblite_basis_type, only: basis_type
      use tblite_wavefunction, only: wavefunction_type, new_wavefunction
      use tblite_scf_potential, only: potential_type, new_potential, add_pot_to_h1

      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> PTB vDZP basis set
      type(basis_type) :: bas
      !> Wavefunction data
      type(wavefunction_type) :: wfn
      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData
      !> Error type
      type(error_type), allocatable, intent(out) :: error
      !> +U potential
      type(plusu_potential_type) :: plusu
      real(wp), parameter :: q(2) = [ &
         &   -1.1113011_wp, &
         &    1.1113086_wp &
         & ]
      real(wp), parameter :: cn(2) = [ &
         &    0.9397642_wp, &
         &    0.9397642_wp &
         & ]
      !> Conversion factor from temperature to energy
      real(wp), parameter :: kt = 3.166808578545117e-06_wp
      real(wp), parameter :: xyz(3, 2) = reshape([ &
      & 2.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp], [3, 2])
      integer, parameter :: nat = 2
      integer, parameter :: at(nat) = [5, 17]
      real(wp), parameter :: plusU_pot_ref(4) = [ &
      &  -0.0023185_wp, & ! 1,1
      &  -0.0018289_wp, & ! 1,2
      &  -0.5266562_wp, & ! 1,21 ; diffferent because of tblite ordering
      &  -1.6745659_wp]   ! 6,24 ; diffferent because of tblite ordering

      call new(mol, at, xyz)
      allocate (ptbData)
      call initPTB(ptbData, mol%num)

      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, bas)
      call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, &
         & nspin=1, kt=300.0_wp * kt)

      wfn%density = 2.0_wp
      wfn%coeff = 0.0_wp
      wfn%qat(:, 1) = q
      call plusu%init(ptbData%plusU, mol, bas, wfn%qat(:, 1), cn)
      call plusu%get_potential(mol, bas, wfn%density(:, :, 1), wfn%coeff(:, :, 1))

      call check_(error, wfn%coeff(1, 1, 1), plusU_pot_ref(1), thr=thr)
      call check_(error, wfn%coeff(1, 2, 1), plusU_pot_ref(2), thr=thr)
      call check_(error, wfn%coeff(1, 20, 1), plusU_pot_ref(3), thr=thr)
      call check_(error, wfn%coeff(8, 26, 1), plusU_pot_ref(4), thr=thr)

   end subroutine test_ptb_plus_U_potential

   subroutine test_ptb_mb16_43_01(error)
      !> PTB overlap matrix calculation
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_calculator, only: TPTBCalculator, newPTBcalculator
      use xtb_ptb_guess, only: get_psh_from_qsh
      !> tblite lib
      use xtb_type_environment, only: TEnvironment, init
      use xtb_type_calculator, only: TCalculator
      use xtb_type_restart, only: TRestart
      use xtb_type_data, only: scc_results

      !> Error type
      type(error_type), allocatable, intent(out) :: error
      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> Generic calculator class
      class(TCalculator), allocatable :: calc
      !> PTB calculator class
      type(TPTBCalculator), allocatable :: ptb, ptb_save
      !> Environment type
      type(TEnvironment) :: env
      !> Restart type
      type(TRestart) :: chk
      !> Calculation results
      type(scc_results) :: res
      !> variables that are needed for calc%singlepoint but useless for PTB
      real(wp) :: energy, gap
      real(wp), allocatable :: gradient(:, :)
      real(wp) :: sigma(3, 3)
      integer :: iat, ish
      !> PTB shell populations
      real(wp), allocatable :: psh(:, :)
      !> PTB reference atomic charges of MB16-43-01
      real(wp), parameter :: q_ref(16) = [ &
         &     0.68170445_wp, &
         &    -0.02179514_wp, &
         &    -0.39271375_wp, &
         &    -0.04860209_wp, &
         &    -0.16998965_wp, &
         &     0.05243959_wp, &
         &    -0.02162296_wp, &
         &    -0.27122167_wp, &
         &    -0.10233830_wp, &
         &     0.05425879_wp, &
         &     0.06345147_wp, &
         &     0.18291012_wp, &
         &    -0.23545020_wp, &
         &    -0.17286158_wp, &
         &     0.01710696_wp, &
         &     0.38467781_wp &
         & ]
      !> PTB reference shell populations of Na atom in MB16-43-01
      real(wp), parameter :: psh_ref(6) = [ &
         &     1.42836761_wp, &
         &     0.55988306_wp, &
         &     0.07385892_wp, &
         &     5.16059363_wp, &
         &     0.86661109_wp, &
         &     0.22898123_wp &
         & ]
      !> PTB reference WBOs of MB16-43-01
      real(wp), parameter :: wbo_ref(2) = [ &
         & 0.616933710808355_wp, &
         & 0.608892999962882_wp &
         & ]

      !> Initialize calculation environment
      call init(env)

      call getMolecule(struc, "mindless01")

      allocate (ptb, ptb_save)
      call newPTBCalculator(env, struc, ptb)
      ptb_save = ptb
      call move_alloc(ptb, calc)

      gap = 0.0_wp
      allocate (gradient(3, struc%n), source=0.0_wp)
      call calc%singlepoint(env, struc, chk, 2, .false., energy, gradient, sigma, &
         & gap, res)
      allocate (psh(ptb_save%bas%nsh, chk%tblite%nspin), source=0.0_wp)
      psh = get_psh_from_qsh(chk%tblite, ptb_save%bas)

      do iat = 1, struc%n
         call check_(error, chk%wfn%q(iat), q_ref(iat), thr=thr2)
      end do
      do ish = 1, ptb_save%bas%nsh_at(1)
         call check_(error, psh(ish, 1), psh_ref(ish), thr=thr2)
      end do
      call check_(error, chk%wfn%wbo(16, 12), wbo_ref(1), thr=thr2)
      call check_(error, chk%wfn%wbo(13, 12), wbo_ref(2), thr=thr2)

   end subroutine test_ptb_mb16_43_01

   subroutine test_ptb_mb16_43_01_charged(error)
      !> PTB overlap matrix calculation
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_calculator, only: TPTBCalculator, newPTBcalculator
      use xtb_ptb_guess, only: get_psh_from_qsh
      !> tblite lib
      use xtb_type_environment, only: TEnvironment, init
      use xtb_type_calculator, only: TCalculator
      use xtb_type_restart, only: TRestart
      use xtb_type_data, only: scc_results

      !> Error type
      type(error_type), allocatable, intent(out) :: error
      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> Generic calculator class
      class(TCalculator), allocatable :: calc
      !> PTB calculator class
      type(TPTBCalculator), allocatable :: ptb, ptb_save
      !> Environment type
      type(TEnvironment) :: env
      !> Restart type
      type(TRestart) :: chk
      !> Calculation results
      type(scc_results) :: res
      !> variables that are needed for calc%singlepoint but useless for PTB
      real(wp) :: energy, gap
      real(wp), allocatable :: gradient(:, :)
      real(wp) :: sigma(3, 3)
      integer :: iat, ish
      !> PTB shell populations
      real(wp), allocatable :: psh(:, :)
      !> PTB reference atomic charges of MB16-43-01
      real(wp), parameter :: q_ref(16) = [ &
         &     0.77426386_wp, &
         &     0.26649239_wp, &
         &     0.30837190_wp, &
         &     0.03983377_wp, &
         &    -0.11302422_wp, &
         &     0.11155364_wp, &
         &    -0.00620840_wp, &
         &    -0.21730028_wp, &
         &    -0.44619143_wp, &
         &     0.10306089_wp, &
         &     0.12243405_wp, &
         &     0.39382297_wp, &
         &     0.21994776_wp, &
         &    -0.01642522_wp, &
         &     0.07058748_wp, &
         &     0.38873424_wp &
         & ]
      !> PTB reference shell populations of Na atom in MB16-43-01
      real(wp), parameter :: psh_ref(6) = [ &
         &     1.42766702_wp, &
         &     0.56155699_wp, &
         &     0.03990374_wp, &
         &     5.16610789_wp, &
         &     0.85170540_wp, &
         &     0.17879509_wp &
         & ]
      !> PTB reference WBOs of MB16-43-01
      real(wp), parameter :: wbo_ref(2) = [ &
         & 0.772617695042871_wp, &
         & 0.526583157894007_wp &
         & ]

      !> Initialize calculation environment
      call init(env)

      call getMolecule(struc, "mindless01")
      struc%chrg = 2.0_wp

      allocate (ptb, ptb_save)
      call newPTBCalculator(env, struc, ptb)
      ptb_save = ptb
      call move_alloc(ptb, calc)

      gap = 0.0_wp
      allocate (gradient(3, struc%n), source=0.0_wp)
      call calc%singlepoint(env, struc, chk, 2, .false., energy, gradient, sigma, &
         & gap, res)
      allocate (psh(ptb_save%bas%nsh, chk%tblite%nspin), source=0.0_wp)
      psh = get_psh_from_qsh(chk%tblite, ptb_save%bas)

      do iat = 1, struc%n
         call check_(error, chk%wfn%q(iat), q_ref(iat), thr=thr3)
      end do
      do ish = 1, ptb_save%bas%nsh_at(1)
         call check_(error, psh(ish, 1), psh_ref(ish), thr=thr3)
      end do
      call check_(error, chk%wfn%wbo(16, 12), wbo_ref(1), thr=thr3)
      call check_(error, chk%wfn%wbo(13, 12), wbo_ref(2), thr=thr3)

   end subroutine test_ptb_mb16_43_01_charged

   subroutine test_ptb_mb16_43_01_efield(error)
      !> PTB overlap matrix calculation
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_calculator, only: TPTBCalculator, newPTBcalculator
      use xtb_ptb_guess, only: get_psh_from_qsh
      !> xtb lib
      use xtb_type_environment, only: TEnvironment, init
      use xtb_type_calculator, only: TCalculator
      use xtb_type_restart, only: TRestart
      use xtb_type_data, only: scc_results
      use xtb_setparam

      !> Error type
      type(error_type), allocatable, intent(out) :: error
      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> Generic calculator class
      class(TCalculator), allocatable :: calc
      !> PTB calculator class
      type(TPTBCalculator), allocatable :: ptb, ptb_save
      !> Environment type
      type(TEnvironment) :: env
      !> Restart type
      type(TRestart) :: chk
      !> Calculation results
      type(scc_results) :: res
      !> variables that are needed for calc%singlepoint but useless for PTB
      real(wp) :: energy, gap
      real(wp), allocatable :: gradient(:, :)
      real(wp) :: sigma(3, 3)
      integer :: iat, ish
      !> PTB shell populations
      real(wp), allocatable :: psh(:, :)
      !> PTB reference atomic charges of MB16-43-01
      real(wp), parameter :: q_ref(16) = [ &
         &    -1.87083162_wp, &
         &    -3.61266413_wp, &
         &     5.08098387_wp, &
         &    -3.30344743_wp, &
         &     6.89679486_wp, &
         &     0.00416323_wp, &
         &    -0.24406779_wp, &
         &     3.84145983_wp, &
         &     2.70491922_wp, &
         &     0.41943598_wp, &
         &     0.78384251_wp, &
         &   -11.00332421_wp, &
         &    -4.96965395_wp, &
         &     2.25944041_wp, &
         &     2.46367706_wp, &
         &     0.54923654_wp &
         & ]
      !> PTB reference shell populations of Na atom in MB16-43-01
      real(wp), parameter :: psh_ref(6) = [ &
         &     1.45512176_wp, &
         &     0.57954907_wp, &
         &     0.65686184_wp, &
         &     5.17897618_wp, &
         &     0.99362525_wp, &
         &     2.00669752_wp &
         & ]
      !> PTB reference WBOs of MB16-43-01
      real(wp), parameter :: wbo_ref(2) = [ &
         & 0.839744207773282_wp, &
         & 1.31529170060079_wp &
         & ]
      set%efield(1) = 0.0_wp
      set%efield(2) = 0.1_wp
      set%efield(3) = 0.3_wp

      !> Initialize calculation environment
      call init(env)

      call getMolecule(struc, "mindless01")

      allocate (ptb, ptb_save)
      call newPTBCalculator(env, struc, ptb)
      ptb_save = ptb
      call move_alloc(ptb, calc)

      gap = 0.0_wp
      allocate (gradient(3, struc%n), source=0.0_wp)
      call calc%singlepoint(env, struc, chk, 2, .false., energy, gradient, sigma, &
         & gap, res)
      allocate (psh(ptb_save%bas%nsh, chk%tblite%nspin), source=0.0_wp)
      psh = get_psh_from_qsh(chk%tblite, ptb_save%bas)

      do iat = 1, struc%n
         call check_(error, chk%wfn%q(iat), q_ref(iat), thr=thr3)
      end do
      do ish = 1, ptb_save%bas%nsh_at(1)
         call check_(error, psh(ish, 1), psh_ref(ish), thr=thr3)
      end do
      call check_(error, chk%wfn%wbo(16, 12), wbo_ref(1), thr=thr3)
      call check_(error, chk%wfn%wbo(13, 12), wbo_ref(2), thr=thr3)
      set%efield = 0.0_wp

   end subroutine test_ptb_mb16_43_01_efield

   subroutine test_ptb_dipmom_caffeine(error)
      !> PTB overlap matrix calculation
      use xtb_ptb_calculator, only: TPTBCalculator, newPTBcalculator
      !> xtb lib
      use xtb_type_environment, only: TEnvironment, init
      use xtb_type_calculator, only: TCalculator
      use xtb_type_restart, only: TRestart
      use xtb_type_data, only: scc_results
      use xtb_setparam

      !> Error type
      type(error_type), allocatable, intent(out) :: error
      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> Generic calculator class
      class(TCalculator), allocatable :: calc
      !> PTB calculator class
      type(TPTBCalculator), allocatable :: ptb
      !> Environment type
      type(TEnvironment) :: env
      !> Restart type
      type(TRestart) :: chk
      !> Calculation results
      type(scc_results) :: res
      !> variables that are needed for calc%singlepoint but useless for PTB
      real(wp) :: energy, gap
      real(wp), allocatable :: gradient(:, :)
      real(wp) :: sigma(3, 3)
      real(wp), parameter :: dipmom_ref(3) = [ &
         & -1.029811_wp, &
         & 1.487256_wp, &
         & 0.005426_wp &
         & ]
      real(wp), parameter :: dipmom_anion_ref(3) = [ &
         & -14.460494, &
         & 6.045541, &
         & 0.284856 &
         & ]
      real(wp), parameter :: dipmom_efield_ref(3) = [ &
         & -6.754478, &
         & 80.047721, &
         & 26.646085 &
         & ]
      integer :: i

      !> Initialize calculation environment
      call init(env)

      call getMolecule(struc, "caffeine")

      !> Neutral
      allocate (ptb)
      call newPTBCalculator(env, struc, ptb)
      call move_alloc(ptb, calc)

      gap = 0.0_wp
      allocate (gradient(3, struc%n), source=0.0_wp)
      call calc%singlepoint(env, struc, chk, 2, .false., energy, gradient, sigma, &
         & gap, res)
      do i = 1, 3
         call check_(error, res%dipole(i), dipmom_ref(i), thr=thr2)
      end do
      deallocate (calc)

      !> Now with anion
      struc%chrg = -2.0_wp
      allocate (ptb)
      call newPTBCalculator(env, struc, ptb)
      call move_alloc(ptb, calc)

      gap = 0.0_wp
      gradient = 0.0_wp
      call calc%singlepoint(env, struc, chk, 2, .false., energy, gradient, sigma, &
         & gap, res)
      do i = 1, 3
         call check_(error, res%dipole(i), dipmom_anion_ref(i), thr=thr2)
      end do
      deallocate (calc)
      struc%chrg = 0.0_wp

      !> Now with electric field
      set%efield(1) = 0.0_wp
      set%efield(2) = 0.1_wp
      set%efield(3) = 0.3_wp
      allocate (ptb)
      call newPTBCalculator(env, struc, ptb)
      call move_alloc(ptb, calc)

      gap = 0.0_wp
      gradient = 0.0_wp
      call calc%singlepoint(env, struc, chk, 2, .false., energy, gradient, sigma, &
         & gap, res)
      do i = 1, 3
         call check_(error, res%dipole(i), dipmom_efield_ref(i), thr=thr4)
      end do
      deallocate (calc)
      set%efield = 0.0_wp

   end subroutine test_ptb_dipmom_caffeine

   subroutine test_ptb_polarizability(error)
      !> PTB overlap matrix calculation
      use xtb_ptb_calculator, only: TPTBCalculator, newPTBcalculator
      !> xtb lib
      use xtb_type_environment, only: TEnvironment, init
      use xtb_type_calculator, only: TCalculator
      use xtb_type_restart, only: TRestart
      use xtb_type_data, only: scc_results
      use xtb_setparam

      !> Error type
      type(error_type), allocatable, intent(out) :: error
      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> Generic calculator class
      class(TCalculator), allocatable :: calc
      !> PTB calculator class
      type(TPTBCalculator), allocatable :: ptb
      !> Environment type
      type(TEnvironment) :: env
      !> Restart type
      type(TRestart) :: chk
      !> Calculation results
      type(scc_results) :: res
      !> variables that are needed for calc%singlepoint but useless for PTB
      real(wp) :: energy, gap
      real(wp), allocatable :: gradient(:, :)
      real(wp) :: sigma(3, 3)
      real(wp), parameter :: alpha_ref(6) = [ &
      & 88.369115_wp, &
      &  0.006262_wp, &
      & 89.779851_wp, &
      & -0.002677_wp, &
      &  0.000768_wp, &
      & 99.524132_wp &
      & ]

      set%elprop = p_elprop_alpha
      !> Initialize calculation environment
      call init(env)

      call getMolecule(struc, "feco5")

      allocate (ptb)
      call newPTBCalculator(env, struc, ptb)
      call move_alloc(ptb, calc)

      gap = 0.0_wp
      allocate (gradient(3, struc%n), source=0.0_wp)
      call calc%singlepoint(env, struc, chk, 2, .false., energy, gradient, sigma, &
         & gap, res)

      call check_(error, res%alpha(1, 1), alpha_ref(1), thr=thr_alpha)
      call check_(error, res%alpha(1, 2), alpha_ref(2), thr=thr_alpha)
      call check_(error, res%alpha(2, 2), alpha_ref(3), thr=thr_alpha)
      call check_(error, res%alpha(1, 3), alpha_ref(4), thr=thr_alpha)
      call check_(error, res%alpha(2, 3), alpha_ref(5), thr=thr_alpha)
      call check_(error, res%alpha(3, 3), alpha_ref(6), thr=thr_alpha)
      set%elprop = p_elprop_dipole

   end subroutine test_ptb_polarizability

   pure function id_to_atom(mol, idparam) result(atomparam)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> PTB parameterization data for each species
      real(wp), intent(in) :: idparam(:, :)
      !> PTB parameterization data for each atom
      real(wp), allocatable :: atomparam(:, :)

      integer :: iat, iid

      allocate (atomparam(size(idparam, 1), mol%nat), source=0.0_wp)
      do iat = 1, mol%nat
         iid = mol%id(iat)
         atomparam(:, iat) = idparam(:, iid)
      end do

   end function id_to_atom
#else
   subroutine test_ptb_not_present(error)
      !> Error type
      type(error_type), allocatable, intent(out) :: error
      if (.not. get_xtb_feature('tblite')) then
         call skip_test(error, "xtb not compiled with tblite support")
         return
      end if
   end subroutine test_ptb_not_present

#endif
end module test_ptb
