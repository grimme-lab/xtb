! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

module xtb_io_writer_gaussian
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule, only : TMolecule, len
   implicit none
   private

   public :: writeMoleculeGaussianExternal
   public :: writeResultsGaussianExternal


contains


subroutine writeMoleculeGaussianExternal(mol, unit)
   type(TMolecule), intent(in) :: mol
   integer, intent(in) :: unit
   integer :: iat

   write(unit, '(4i10)') len(mol), 1, nint(mol%chrg), mol%uhf
   do iat = 1, len(mol)
      write(unit, '(i10,4f20.12)') mol%at(iat), mol%xyz(:, iat), 0.0_wp
   end do

end subroutine writeMoleculeGaussianExternal


subroutine writeResultsGaussianExternal(unit, energy, dipole, gradient, hessian, dipgrad)
   integer, intent(in)  :: unit
   real(wp), intent(in) :: energy
   real(wp), intent(in) :: dipole(:)
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: hessian(:, :)
   real(wp), intent(in), optional :: dipgrad(:)

   integer :: i,j,ij,nat

   if (present(gradient)) then
      nat = size(gradient, 2)
   else if (present(hessian)) then
      nat = size(hessian, 2) / 3
   else if (present(dipgrad)) then
      nat = size(dipgrad) / 9
   else
      nat = 0
   end if

   write(unit, '(4D20.12)') energy, dipole

   if (present(gradient)) then
      write(unit, '(3D20.12)') gradient
   else
      ! In case we have no gradient, but either hessian or dipole gradient
      ! the gradient entry is padded with zeros
      if (nat > 0) write(unit, '(3D20.12)') spread(0.0_wp, 1, 3*nat)
   endif

   ! First, we fake the polarizability and the dipole derivatives, which the
   ! Gaussian Manual says should be of this form,
   !
   ! Polar(I), I=1,6        3D20.12   2 rows of 3 x 0.0
   write(unit, '(3D20.12)') spread(0.0_wp, 1, 6)

   ! DDip(I), I=1,9*NAtoms  3D20.12   9 x Natoms or 3 natoms rows of 3 x .0.0
   if (present(dipgrad)) then
      write(unit, '(3D20.12)') dipgrad
   else
      if (nat > 0) write(unit, '(3D20.12)') spread(0.0_wp, 1, 9*nat)
   end if

   if (present(hessian)) then
      ! Now we are iterating over Hessian matrix elements. We will
      ! append those to the output file in the correct format, given in
      ! the Gaussian manual as
      !
      ! FFX(I), I=1,(3*NAtoms*(3*NAtoms+1))/2      3D20.12

      ij = 0
      do i = 1, size(hessian, 2)
         do j = 1, i
            ij = ij + 1
            write(unit, '(D20.12)', advance='no') hessian(j, i)
            if (mod(ij, 3) == 0) write(unit, '(a)')
         end do
      end do
   end if

end subroutine writeResultsGaussianExternal


end module xtb_io_writer_gaussian
