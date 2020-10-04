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


subroutine writeResultsGaussianExternal(mol, unit, energy, dipole, gradient, hess)
   implicit none
   type(TMolecule), intent(in) :: mol
   integer, intent(in)  :: unit
   real(wp), intent(in) :: energy
   real(wp), intent(in) :: dipole(:)
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: hess(:, :)

   integer :: i,j,count
   real(wp) :: buf(3)

   write(unit, '(4D20.12)') energy, dipole

   if (present(gradient)) then
      write(unit, '(3D20.12)') gradient
   elseif (present(hess)) then
      ! In case the gradient is not computed but the Hessian is, we output
      ! gradients of zero. I don't know what would cause that to be honest...
      write(unit, '(3D20.12)') (/ 0.0, 0.0, 0.0 /)
   endif

   if (present(hess)) then
      ! First, we fake the polarizability and the dipole derivatives, which the
      ! Gaussian Manual says should be of this form,
      !
      ! Polar(I), I=1,6        3D20.12   2 rows of 3 x 0.0
      ! DDip(I), I=1,9*NAtoms  3D20.12   9 x Natoms or 3 natoms rows of 3 x .0.0
      do i=1,3*mol%n + 2
         write(unit, '(3D20.12)') (/ 0.0, 0.0, 0.0 /)
      enddo

      ! Now we are iterating over Hessian matrix elements. We will
      ! append those to the output file in the correct format, given in
      ! the Gaussian manual as
      !
      ! FFX(I), I=1,(3*NAtoms*(3*NAtoms+1))/2      3D20.12

      ! that is, we need to print in group of three numbers, which is why the
      ! whole buffer / count thing is present
      count = 1
      do i=1,3*mol%n
         do j=1,i
            buf(count) = hess(i,j)

            if (count .eq. 3) then
               count = 1
               write(unit, '(3D20.12)') buf
            else
               count = count + 1
            end if

         enddo
      enddo
   endif

end subroutine writeResultsGaussianExternal


end module xtb_io_writer_gaussian
