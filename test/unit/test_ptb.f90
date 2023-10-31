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

   use xtb_type_molecule, only: TMolecule
   use xtb_test_molstock, only: getMolecule

   use testdrive, only: new_unittest, unittest_type, error_type, check_ => check, test_failed

   use mctc_io, only: structure_type, new

   implicit none
   private

   real(wp), parameter :: thr = 1.0e-7_wp

   public :: collect_ptb

contains

!> Collect all exported unit tests
   subroutine collect_ptb(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("basis", test_ptb_basis) &
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

      call new(mol, struc%at, struc%xyz, struc%chrg, struc%uhf, struc%lattice)
      if (allocated(mol%pdb)) mol%pdb = struc%pdb
      if (allocated(mol%sdf)) mol%sdf = struc%sdf
      mol%periodic = .false.

      call add_vDZP_basis(mol, bas)

      call check_(error, bas%nsh, nsh_exp, "Number of shells not matching to expected value.")
      call check_(error, bas%nao, nao_exp, "Number of AOs not matching to expected value.")
      call check_(error, bas%maxl, maxl_exp, "Maximum angular momentum not matching to expected value.")
      call check_(error, bas%cgto(1,1)%alpha(1), 3.54364182058582_wp, thr=thr)
      call check_(error, bas%cgto(1,2)%alpha(1), 81.8867808750392_wp, thr=thr)

   end subroutine test_ptb_basis

end module test_ptb
