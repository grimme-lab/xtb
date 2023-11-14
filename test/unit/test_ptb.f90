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

module test_ptb

   use mctc_env, only: wp

   use xtb_type_molecule, only: TMolecule, assignment(=)
   use xtb_test_molstock, only: getMolecule

   use testdrive, only: new_unittest, unittest_type, error_type, check_ => check, test_failed

   use mctc_io, only: structure_type, new

   implicit none
   private

   real(wp), parameter :: thr = 1.0e-7_wp
   real(wp), parameter :: thr2 = 1.0e-5_wp

   public :: collect_ptb

contains

!> Collect all exported unit tests
   subroutine collect_ptb(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("basis", test_ptb_basis), &
                  new_unittest("eeq", test_ptb_eeq), &
                  new_unittest("overlap", test_ptb_overlap), &
                  new_unittest("overlap_h0", test_ptb_overlap_h0), &
                  new_unittest("overlap_sx", test_ptb_overlap_SX), &
                  new_unittest("v_ecp", test_ptb_V_ECP), &
                  new_unittest("selfenergies", test_ptb_selfenergies), &
                  new_unittest("hamiltonian_h0", test_ptb_hamiltonian_h0), &
                  new_unittest("v_xc", test_ptb_V_XC) &
                  ]

   end subroutine collect_ptb

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
      use xtb_ptb_overlaps, only: get_scaled_integrals

      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      type(error_type), allocatable, intent(out) :: error
      !> (Scaled) overlap matrix
      real(wp), allocatable :: overlap(:, :)
      character(len=:), allocatable :: message
      real(wp), parameter :: overlap_exp(6) = [ &
      & 0.93209460_wp, & ! 1,2
      & 0.35489609_wp, & ! 1,3
      & 0.65682608_wp, & ! 2,3
      & 0.05627743_wp, & ! 1,15
      & -0.14217162_wp, &  ! 1,24; diffferent because of tblite ordering
      & 0.41844087_wp] ! 14,23; diffferent because of tblite ordering

      message = "Overlap matrix element not matching to expected value."

      call getMolecule(struc, "mgh2")
      mol = struc

      call get_scaled_integrals(mol, overlap)
      call check_(error, overlap(1, 2), overlap_exp(1), thr=thr, &
      & message=message)
      call check_(error, overlap(1, 3), overlap_exp(2), thr=thr, &
      & message=message)
      call check_(error, overlap(2, 3), overlap_exp(3), thr=thr, &
      & message=message)
      call check_(error, overlap(1, 15), overlap_exp(4), thr=thr, &
      & message=message)
      call check_(error, overlap(1, 23), overlap_exp(5), thr=thr, &
      & message=message)
      call check_(error, overlap(12, 22), overlap_exp(6), thr=thr, &
      & message=message)

   end subroutine test_ptb_overlap

   subroutine test_ptb_overlap_h0(error)
      use xtb_ptb_overlaps, only: get_scaled_integrals
      use xtb_ptb_param, only: kalphah0l
      use xtb_ptb_vdzp, only: max_shell

      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      type(error_type), allocatable, intent(out) :: error
      !> (Scaled) overlap matrix
      real(wp), allocatable :: overlap(:, :)
      !> Temporary array for exponent scaling factors specific to
      !> unique atoms in the molecule
      real(wp), allocatable :: expscal(:, :)
      !> Loop variables
      integer :: isp, izp
      character(len=:), allocatable :: message
      real(wp), parameter :: overlap_exp(6) = [ &
      & 0.95689468_wp, & ! 1,2
      & 0.39195790_wp, & ! 1,3
      & 0.62961212_wp, & ! 2,3
      & 0.03782850_wp, & ! 1,15
      &-0.13826216_wp, &  ! 1,24; diffferent because of tblite ordering
      & 0.43334922_wp] ! 14,23; diffferent because of tblite ordering

      message = "Scaled overlap matrix element not matching to expected value."

      call getMolecule(struc, "mgh2")
      mol = struc

      allocate (expscal(max_shell, mol%nid), source=0.0_wp)
      do isp = 1, mol%nid
         izp = mol%num(isp)
         expscal(:, isp) = kalphah0l(:, izp)
      end do
      call get_scaled_integrals(mol, overlap, expscal)
      call check_(error, overlap(1, 2), overlap_exp(1), thr=thr, &
      & message=message)
      call check_(error, overlap(1, 3), overlap_exp(2), thr=thr, &
      & message=message)
      call check_(error, overlap(2, 3), overlap_exp(3), thr=thr, &
      & message=message)
      call check_(error, overlap(1, 15), overlap_exp(4), thr=thr, &
      & message=message)
      call check_(error, overlap(1, 23), overlap_exp(5), thr=thr, &
      & message=message)
      call check_(error, overlap(12, 22), overlap_exp(6), thr=thr, &
      & message=message)

   end subroutine test_ptb_overlap_h0

   subroutine test_ptb_overlap_SX(error)
      use xtb_ptb_overlaps, only: get_scaled_integrals
      use xtb_ptb_mmlpopanalysis, only: get_mml_overlaps
      use xtb_ptb_param, only: ptbGlobals
      use tblite_basis_type, only: basis_type
      use xtb_ptb_vdzp, only: add_vDZP_basis

      !> PTB vDZP basis set
      type(basis_type) :: bas
      !> Structure type (xtb)
      type(TMolecule) :: struc
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      type(error_type), allocatable, intent(out) :: error
      !> (Scaled) overlap matrix
      real(wp), allocatable :: overlap(:, :)
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

      message = "Solved overlap matrix (^(1-x)) element not matching to expected value."

      call getMolecule(struc, "mgh2")
      mol = struc

      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, bas)
      call get_scaled_integrals(mol, overlap)
      allocate (overlap_sx(bas%nao, bas%nao), overlap_oneminusx(bas%nao, bas%nao))
      call get_mml_overlaps(bas, overlap, ptbGlobals%mlmix, overlap_sx, &
         & overlap_oneminusx)

      write (*, *) overlap_oneminusx(1, 2), overlap_oneminusx_exp(1)
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
      !> tblite basis set type
      use tblite_basis_type, only: basis_type
      !> PTB core basis set generation
      use xtb_ptb_vdzp, only: add_vDZP_basis
      use xtb_ptb_corebasis, only: add_core_basis, get_Vecp
      !> PTB overlap matrix calculation
      use xtb_ptb_overlaps, only: get_scaled_integrals
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_param, only: initPTB

      !> PTB vDZP basis set and core basis set
      type(basis_type) :: bas, cbas
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData
      !> (Scaled) overlap matrix
      real(wp), allocatable :: overlap(:, :)
      !> Normalization factors
      real(wp), allocatable :: norm_overlap(:)
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

      call new(mol, at, xyz)
      allocate (ptbData)
      call initPTB(ptbData, mol%num)

      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, bas)
      !> Add the core basis set to 'cbas' basis set type
      call add_core_basis(mol, ptbData%corepotential, cbas)
      !> -> for normalization factors
      call get_scaled_integrals(mol, overlap, norm=norm_overlap)
      !> V_ECP via PTB core basis
      call get_Vecp(mol, ptbData%corepotential, bas, cbas, norm_overlap, vecp)

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
      use tblite_basis_type, only: basis_type
      !> PTB core basis set generation
      use xtb_ptb_vdzp, only: add_vDZP_basis
      use xtb_ptb_overlaps, only: get_scaled_integrals
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_param, only: initPTB

      !> PTB vDZP basis set and core basis set
      type(basis_type) :: bas
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData
      !> (Scaled) overlap matrix
      real(wp), allocatable :: overlap_h0(:, :)
      !> Error type
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: message
      !> Structure
      real(wp), parameter :: xyz(3, 2) = reshape([ &
      & 2.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp], [3, 2])
      integer, parameter :: nat = 2
      integer, parameter :: at(nat) = [5, 17]

      real(wp), parameter :: h0_ref(4) = [ &
      &  0.0_wp, & ! 1,1 ; diffferent because of tblite ordering
      &  0.0_wp, & ! 1,3 ; diffferent because of tblite ordering
      &  0.0_wp, & ! 3,5 ; diffferent because of tblite ordering
      &  0.0_wp]   ! 9,9 ; diffferent because of tblite ordering

      !> Hamiltonian matrix
      real(wp), allocatable :: hamiltonian(:, :)

      call new(mol, at, xyz)
      allocate (ptbData)
      call initPTB(ptbData, mol%num)
      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, bas)

      allocate (hamiltonian(bas%nao, bas%nao), source=0.0_wp)

      message = "H0 matrix element not matching to expected value."
      call check_(error, hamiltonian(1, 1), h0_ref(1), thr=thr2, &
      & message=message)
      call check_(error, hamiltonian(1, 5), h0_ref(2), thr=thr2, &
      & message=message)
      call check_(error, hamiltonian(5, 8), h0_ref(3), thr=thr2, &
      & message=message)
      call check_(error, hamiltonian(12, 12), h0_ref(4), thr=thr2, &
      & message=message)
   end subroutine test_ptb_hamiltonian_h0

   subroutine test_ptb_V_XC(error)
      !> tblite basis set type
      use tblite_basis_type, only: basis_type
      !> PTB core basis set generation
      use xtb_ptb_vdzp, only: add_vDZP_basis
      use xtb_ptb_corebasis, only: add_core_basis, get_Vecp
      !> PTB overlap matrix calculation
      use xtb_ptb_overlaps, only: get_scaled_integrals
      use xtb_ptb_data, only: TPTBData
      use xtb_ptb_param, only: initPTB
      use xtb_ptb_paulixc, only: calc_Vxc_pauli

      !> PTB vDZP basis set
      type(basis_type) :: bas
      !> Structure type (mctc-lib)
      type(structure_type) :: mol
      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData
      !> Error type
      type(error_type), allocatable, intent(out) :: error
      !> (Scaled) overlap matrix
      real(wp), allocatable :: overlap(:, :)
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
      real(wp), allocatable :: overlap_xc(:, :), Vxc(:, :)
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

      call new(mol, at, xyz)
      allocate (ptbData)
      call initPTB(ptbData, mol%num)

      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, bas)
      call get_scaled_integrals(mol, overlap_xc, alpha_scal=ptbData%pauli%klalphaxc)
      call calc_Vxc_pauli(mol, bas, shellpops, overlap_xc, levels, ptbData%pauli%kxc1, Vxc)

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

end module test_ptb
