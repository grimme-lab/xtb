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

module test_molecule
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
   implicit none
   private

   public :: collect_molecule

contains

!> Collect all exported unit tests
subroutine collect_molecule(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("mic-distances", test_class_molecule_mic_distances), &
      new_unittest("axis-trafo", test_class_molecule_axis_trafo), &
      new_unittest("xyz-writer-nonpbc", test_xyz_writer_nonpbc), &
      new_unittest("xyz-writer-pbc", test_xyz_writer_pbc), &
      new_unittest("update-nowrap-pbc", test_update_nowrap_pbc) &
      ]

end subroutine collect_molecule


subroutine test_class_molecule_mic_distances(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule
   use xtb_type_param
   use xtb_eeq
   use xtb_disp_ncoord
   use xtb_pbc_tools
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-10_wp
   ! SiO2 (random, no symmetry)
   integer, parameter :: nat = 6
   integer, parameter :: at(nat) = [14,14,8,8,8,8]
   real(wp),parameter :: abc(3,nat) = reshape(&
      &[.095985472469032_wp,     .049722204206931_wp,     0.10160624337938_wp, &
      & 0.54722204206931_wp,     0.52863628207623_wp,     0.38664208660311_wp, &
      & 0.29843937068984_wp,     0.39572194413818_wp,     0.20321248675876_wp, &
      & 0.23364982659922_wp,     0.85647058758674_wp,     0.31884968761485_wp, &
      & 0.72250232459952_wp,     0.65548544066844_wp,     .056207709103487_wp, &
      & 0.70514214000043_wp,     0.28321754549582_wp,     0.36424822189074_wp],&
      & shape(abc))
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[ 8.7413053236641_wp,      0.0000000000000_wp,      0.0000000000000_wp,   &
      &  0.0000000000000_wp,      8.7413053236641_wp,      0.0000000000000_wp,   &
      &  0.0000000000000_wp,      0.0000000000000_wp,      8.7413053236641_wp],  &
      & shape(lattice))
   integer, parameter :: wsc_rep(3) = [1,1,1]
   real(wp),allocatable :: xyz(:,:)

   type(TMolecule)       :: mol

   allocate(xyz(3, nat))
   call coord_trafo(nat,lattice,abc,xyz)
   call init(mol, at, xyz, lattice=lattice)

   call check(error, mol%dist(1,1),8.7413053236641_wp, thr=thr)
   call check(error, mol%dist(3,6),3.9480992656526_wp, thr=thr)
   call check(error, mol%dist(2,3),2.9411012169549_wp, thr=thr)
   call check(error, mol%dist(1,4),2.8120981086302_wp, thr=thr)

   call mol%deallocate

end subroutine test_class_molecule_mic_distances

subroutine test_class_molecule_axis_trafo(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule
   use xtb_type_param
   use xtb_disp_ncoord
   use xtb_eeq
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 7
   integer, parameter :: at(nat) = [6,8,1,6,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.236995341_wp,     0.794935672_wp,     0.000000000_wp,  &
      & 2.159203072_wp,    -0.439230068_wp,     0.000000000_wp,  &
      & 0.329518229_wp,     2.884317091_wp,     0.000000000_wp,  &
      &-2.377225608_wp,    -0.313030650_wp,     0.000000000_wp,  &
      &-2.288077922_wp,    -2.365932814_wp,     0.000000000_wp,  &
      &-3.405901173_wp,     0.355579334_wp,     1.660010899_wp,  &
      &-3.405901173_wp,     0.355579335_wp,    -1.660010899_wp   &
      & ],shape(xyz))

   type(TMolecule)       :: mol
   real(wp) :: molmass
   real(wp) :: center(3)
   real(wp) :: moments(3)

   call init(mol, at, xyz)

   center = mol%center_of_geometry()
   call check(error, center(1),-1.2501984620000_wp, thr=thr)
   call check(error, center(2), .18174541428571_wp, thr=thr)
   call check(error, center(3), 0.0000000000000_wp, thr=thr)

   call mol%shift_to_center_of_geometry
   call check(error, mol%xyz(2,1),.61319025771429_wp, thr=thr)
   call check(error, mol%xyz(1,3),1.5797166910000_wp, thr=thr)
   call check(error, mol%xyz(3,6),1.6600108990000_wp, thr=thr)

   molmass = sum(mol%atmass)
   call check(error, molmass,80303.049694083_wp, thr=1.0e-5_wp) ! in au
   center = mol%center_of_mass()
   call check(error, center(1),1.2502034210661_wp, thr=thr)

   call mol%shift_to_center_of_mass
   call check(error, mol%xyz(2,4),-0.31302884934293_wp, thr=thr)
   call check(error, mol%xyz(1,1), 0.23699038193392_wp, thr=thr)

   moments = mol%moments_of_inertia()
   call check(error, moments(1),57768.744315301_wp, thr=1.0e-5_wp)
   call check(error, moments(3),361019.51690412_wp, thr=1.0e-5_wp)

!   call mol%align_to_principal_axes(break_symmetry = .false.)
!   call check(error, mol%xyz(2,1),-.79493732581159_wp, thr=thr)
!   call check(error, mol%xyz(1,7), 3.4059059117382_wp, thr=thr)
!
!   call mol%align_to_principal_axes(break_symmetry = .true.)
!   call check(error, mol%xyz(2,1),-.79493703211973_wp, thr=thr)
!   call check(error, mol%xyz(1,7), 3.4059054710786_wp, thr=thr)

   call mol%deallocate

end subroutine test_class_molecule_axis_trafo

subroutine test_xyz_writer_nonpbc(error)
   use mctc_io, only : structure_type
   use xtb_io_writer, only : writeMolecule
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_filetypes, only : fileType
   use xtb_type_molecule, only : TMolecule, init, assignment(=)
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8, 1, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      &  0.0000000_wp,  0.0000000_wp,  0.2216670_wp, &
      &  0.0000000_wp,  1.4309280_wp, -0.8866680_wp, &
      &  0.0000000_wp, -1.4309280_wp, -0.8866680_wp], [3, nat])

   type(TMolecule) :: mol
   type(structure_type) :: struc
   integer :: unit
   character(len=256) :: line1, line2

   call init(mol, at, xyz)
   struc = mol

   call check(error, size(struc%lattice), 0)
   if (allocated(error)) return
   call check(error, any(struc%periodic), .false.)
   if (allocated(error)) return

   open(newunit=unit, status='scratch')
   call writeMolecule(mol, unit, format=fileType%xyz)
   rewind(unit)
   read(unit, '(a)') line1
   read(unit, '(a)') line2
   close(unit)

   call check(error, index(line2, 'Lattice='), 0)
   if (allocated(error)) return
   call check(error, index(line2, 'pbc='), 0)
   if (allocated(error)) return
   call check(error, index(line2, 'Properties='), 0)
   if (allocated(error)) return

   call mol%deallocate

end subroutine test_xyz_writer_nonpbc

subroutine test_xyz_writer_pbc(error)
   use mctc_io, only : structure_type
   use xtb_io_writer, only : writeMolecule
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_filetypes, only : fileType
   use xtb_type_molecule, only : TMolecule, init, assignment(=)
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 2
   integer, parameter :: at(nat) = [1, 8]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp], [3, nat])
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      & 10.0_wp,  0.0_wp,  0.0_wp, &
      &  0.0_wp, 10.0_wp,  0.0_wp, &
      &  0.0_wp,  0.0_wp, 10.0_wp], [3, 3])

   type(TMolecule) :: mol
   type(structure_type) :: struc
   integer :: unit
   character(len=512) :: line1, line2

   call init(mol, at, xyz, lattice=lattice)
   struc = mol

   call check(error, size(struc%lattice), 9)
   if (allocated(error)) return
   call check(error, all(struc%periodic), .true.)
   if (allocated(error)) return

   open(newunit=unit, status='scratch')
   call writeMolecule(mol, unit, format=fileType%xyz)
   rewind(unit)
   read(unit, '(a)') line1
   read(unit, '(a)') line2
   close(unit)

   call check(error, index(line2, 'Lattice=') > 0, .true.)
   if (allocated(error)) return
   call check(error, index(line2, 'pbc="T T T"') > 0, .true.)
   if (allocated(error)) return

   call mol%deallocate

end subroutine test_xyz_writer_pbc

subroutine test_update_nowrap_pbc(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule, only : TMolecule, init
   type(error_type), allocatable, intent(out) :: error
   real(wp), parameter :: thr = 1.0e-12_wp
   integer, parameter :: nat = 2
   integer, parameter :: at(nat) = [1, 8]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 1.0_wp, 2.0_wp, 3.0_wp, &
      & 4.0_wp, 5.0_wp, 6.0_wp], [3, nat])
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      & 10.0_wp,  0.0_wp,  0.0_wp, &
      &  0.0_wp, 10.0_wp,  0.0_wp, &
      &  0.0_wp,  0.0_wp, 10.0_wp], [3, 3])
   type(TMolecule) :: mol

   call init(mol, at, xyz, lattice=lattice)
   mol%xyz(1, 1) = -1.5_wp
   mol%lattice = 12.0_wp * reshape([ &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 1.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])

   call mol%update(wrap=.false.)
   call check(error, mol%xyz(1, 1), -1.5_wp, thr=thr)
   if (allocated(error)) return
   call check(error, mol%abc(1, 1), 0.875_wp, thr=thr)
   if (allocated(error)) return
   call check(error, mol%cellpar(1), 12.0_wp, thr=thr)
   if (allocated(error)) return
   call check(error, mol%volume, 1728.0_wp, thr=thr)
   if (allocated(error)) return

   call mol%update
   call check(error, mol%xyz(1, 1), 10.5_wp, thr=thr)

   call mol%deallocate
end subroutine test_update_nowrap_pbc

end module test_molecule
