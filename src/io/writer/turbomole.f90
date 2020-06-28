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

module xtb_io_writer_turbomole
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule, only : TMolecule, len
   use xtb_pbc_tools
   implicit none
   private

   public :: writeMoleculeCoord
   public :: writeResultsTurbomole
   public :: writeNormalModesTurbomole


contains


subroutine writeMoleculeCoord(mol, unit)
   class(TMolecule), intent(in) :: mol
   integer, intent(in) :: unit
   integer :: i

   write(unit,'(a)') "$coord"
   do i = 1, len(mol)
      write(unit,'(3f20.14,6x,a)') mol%xyz(:,i), trim(mol%sym(i))
   enddo
   write(unit,'(a,1x,i0)') "$periodic", mol%npbc
   if (mol%npbc > 0) then
      write(unit,'(a)') "$lattice bohr"
      write(unit,'(3f20.14)') mol%lattice
   endif
   write(unit,'(a)') "$end"

end subroutine writeMoleculeCoord


subroutine writeResultsTurbomole(mol, unit, energy, gradient, sigma)
   type(TMolecule), intent(in) :: mol
   integer, intent(in), optional :: unit
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gradient(:, :)
   real(wp), intent(in), optional :: sigma(:, :)
   integer :: ien, igr, igl
   real(wp) :: en

   if (present(unit)) then
      ien = unit
      igr = unit
      igl = unit
   else
      if (present(energy)) then
         call open_file(ien, "energy", "w")
      end if
      if (present(gradient)) then
         call open_file(igr, "gradient", "w")
      end if
      if (present(sigma) .and. mol%npbc > 0) then
         call open_file(igl, "gradlatt", "w")
      end if
   end if

   if (present(energy)) then
      en = energy
      call writeEnergyTurbomole(ien, en)
   else
      en = 0.0_wp
   end if

   if (present(gradient)) then
      call writeGradientTurbomole(igr, mol%xyz, mol%sym, en, gradient)
   end if

   if (present(sigma) .and. mol%npbc > 0) then
      call writeGradLattTurbomole(igl, mol%lattice, en, sigma)
   end if

   if (present(unit)) then
      if (present(energy) .or. present(gradient) &
         & .or. (present(sigma) .and. mol%npbc > 0)) then
         write(unit, '("$end")')
      end if
   else
      if (present(energy)) then
         write(ien, '("$end")')
         call close_file(ien)
      end if
      if (present(gradient)) then
         write(igr, '("$end")')
         call close_file(igr)
      end if
      if (present(sigma) .and. mol%npbc > 0) then
         write(igl, '("$end")')
         call close_file(igl)
      end if
   end if

end subroutine writeResultsTurbomole


subroutine writeEnergyTurbomole(unit, energy)
   integer, intent(in) :: unit
   real(wp), intent(in) :: energy
   write(unit, '("$energy")')
   write(unit, '(i6,3F18.11)') 1, spread(energy, 1, 3)
end subroutine writeEnergyTurbomole


subroutine writeGradientTurbomole(unit, xyz, sym, energy, gradient)
   integer, intent(in) :: unit
   real(wp), intent(in) :: xyz(:, :)
   character(len=*), intent(in) :: sym(:)
   real(wp), intent(in) :: energy
   real(wp), intent(in) :: gradient(:, :)
   integer :: i

   write(unit, '("$grad")')
   write(unit, '(2x,"cycle =",1x,i6,4x,"SCF energy =",f18.11,3x,'//&
      &        '"|dE/dxyz| =",f10.6)') 1, energy, norm2(gradient)
   do i = 1, size(xyz, dim=2)
      write(unit, '(3(F20.14,2x),4x,a2)') xyz(1,i), xyz(2,i), xyz(3,i), sym(i)
   end do
   do i = 1, size(gradient, dim=2)
      write(unit, '(3ES22.13)') gradient(1,i), gradient(2,i), gradient(3,i)
   end do

end subroutine writeGradientTurbomole


subroutine writeGradLattTurbomole(unit, lattice, energy, sigma)
   integer, intent(in) :: unit
   real(wp), intent(in) :: lattice(:, :)
   real(wp), intent(in) :: energy
   real(wp), intent(in) :: sigma(:, :)
   real(wp) :: gradlatt(3, 3), inv_lat(3, 3)
   integer :: i

   inv_lat = mat_inv_3x3(lattice)
   call sigma_to_latgrad(sigma, inv_lat, gradlatt)

   write(unit, '("$gradlatt")')
   write(unit, '(2x,"cycle =",1x,i6,4x,"SCF energy =",f18.11,3x,'//&
      &        '"|dE/dlatt| =",f10.6)') 1, energy, norm2(gradlatt)
   do i = 1, size(lattice, dim=2)
      write(unit, '(3(F20.14,2x))') lattice(1, i), lattice(2, i), lattice(3, i)
   end do
   do i = 1, size(gradlatt, dim=2)
      write(unit, '(3ES22.13)') gradlatt(1, i), gradlatt(2, i), gradlatt(3, i)
   end do
end subroutine writeGradLattTurbomole


subroutine writeNormalModesTurbomole(unit, atmass, normalModes)
   implicit none
   integer, intent(in) :: unit
   real(wp), intent(in) :: atmass(:) !< atom mass
   real(wp), intent(in) :: normalModes(:, :) !< NM
   real(wp), allocatable :: h(:, :)

   integer :: i,j,k,mincol,maxcol,ic,n,nat3
   character(len=:), allocatable :: frmt

   h = normalModes
   n = size(atmass)
   nat3 = size(normalModes, 1)

   do j=1,nat3
      do i=1,n
         do k=1,3
            ic = 3*(i-1) + k
            h(ic,j) = h(ic,j) / sqrt(atmass(i))  ! remove mass (needed for vibesTM compatibility)
         enddo
      enddo
   enddo

   frmt = '(i2,i3,5f15.10)'
   if (nat3.gt.99) frmt = '(i3,i2,5f15.10)'

   write(unit, '(a)') '$vibrational normal modes'

   do i = 1, nat3
      ic = 0
      maxcol = 0
      do while(maxcol.lt.nat3)
         mincol = maxcol+1
         maxcol = min(maxcol+5,nat3)
         ic=ic+1
         write(unit,frmt) i,ic,(h(i,j),j=mincol,maxcol) ! right TM order
!        write(unit,frmt) i,ic,(h(j,i),j=mincol,maxcol) ! natural, state order
      end do
   end do

   write(unit, '(a)') '$end'

end subroutine writeNormalModesTurbomole


end module xtb_io_writer_turbomole
