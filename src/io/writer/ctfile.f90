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

module xtb_io_writer_ctfile
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_symbols, only : toNumber, toSymbol
   use xtb_mctc_version, only : version
   use xtb_type_buffer, only : tb_buffer
   use xtb_type_molecule, only : TMolecule, len
   implicit none
   private

   public :: writeMoleculeMolfile, writeMoleculeSDF


contains


subroutine writeMoleculeSDF(mol, unit, energy, gnorm)
   class(TMolecule), intent(in) :: mol
   integer, intent(in) :: unit
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gnorm
   type(tb_buffer) :: sd_values
   character(len=:), allocatable :: line
   character(len=*), parameter :: sd_format = &
      & '("> <",a,">",/,f20.12,/)'

   call writeMoleculeMolfile(mol, unit, "xtb: "//version)

   sd_values = mol%info
   call sd_values%reset
   do while(sd_values%next())
      call sd_values%getline(line)
      write(unit, '(a)') line
   enddo

   if (present(energy)) then
      write(unit, sd_format) "total energy / Eh", energy
   endif

   if (present(gnorm)) then
      write(unit, sd_format) "gradient norm / Eh/a0", gnorm
   endif

   write(unit, '("$$$$")')

end subroutine writeMoleculeSDF


subroutine writeMoleculeMolfile(mol, unit, comment_line)
   use xtb_mctc_convert
   class(TMolecule), intent(in) :: mol
   integer, intent(in) :: unit
   character(len=*), intent(in) :: comment_line
   integer, parameter :: list4(4) = 0
   integer :: iatom, ibond, iatoms(3), list12(12)
   logical :: has_sdf_data
   integer, parameter :: charge_to_ccc(-3:3) = [7, 6, 5, 0, 3, 2, 1]
   character(len=8)  :: date
   character(len=10) :: time

   call date_and_time(date, time)

   write(unit, '(a)') mol%name
   write(unit, '(2x,"xtb",5x,3a2,a4,"3D")') &
      &  date(5:6), date(7:8), date(3:4), time(:4)
   write(unit, '(a)') comment_line
   write(unit, '(3i3,3x,2i3,12x,i3,1x,a5)') &
      &  len(mol), len(mol%bonds), 0, 0, 0, 999, 'V2000'

   has_sdf_data = allocated(mol%sdf)

   do iatom = 1, mol%n
      if (has_sdf_data) then
         list12 = [mol%sdf(iatom)%isotope, 0, 0, 0, 0, mol%sdf(iatom)%valence, &
            & 0, 0, 0, 0, 0, 0]
      else
         list12 = 0
      endif
      write(unit, '(3f10.4,1x,a3,i2,11i3)') &
         & mol%xyz(:, iatom)*autoaa, mol%sym(iatom), list12
   enddo

   do ibond = 1, len(mol%bonds)
      call mol%bonds%get_item(ibond, iatoms)
      write(unit, '(7i3)') &
         & iatoms, list4
   enddo

   if (has_sdf_data) then
      if (sum(mol%sdf%charge) /= nint(mol%chrg)) then
         write(unit, '(a,*(i3,1x,i3,1x,i3))') "M  CHG", 1, 1, nint(mol%chrg)
      else
         do iatom = 1, mol%n
            if (mol%sdf(iatom)%charge /= 0) then
               write(unit, '(a,*(i3,1x,i3,1x,i3))') &
                  & "M  CHG", 1, iatom, mol%sdf(iatom)%charge
            end if
         end do
      end if
   else
      if (nint(mol%chrg) /= 0) then
         write(unit, '(a,*(i3,1x,i3,1x,i3))') "M  CHG", 1, 1, nint(mol%chrg)
      end if
   end if

   write(unit, '(a)') "M  END"

end subroutine writeMoleculeMolfile


end module xtb_io_writer_ctfile
