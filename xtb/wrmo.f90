! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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

subroutine write_tm_mos(iunit,n,at,basis,wfn)
   use iso_fortran_env, wp => real64, istdout => output_unit
   use mctc_econv
   use tbdef_wavefunction
   use tbdef_basisset
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   type(tb_basisset),    intent(in) :: basis
   type(tb_wavefunction),intent(in) :: wfn

   character(len=*),parameter :: keyword = "$scfmo   expanded   format(4d20.14)"
   character(len=*),parameter :: format = "(4d20.14)"
   integer :: i,j

   write(iunit,'(a)') keyword
   write(iunit,'(a)') &
      '# molecular orbitals from extended tight binding calculation'
   do i = 1, wfn%nao
      write(iunit,'(i6,2x,"a",6x,"eigenvalue=",d20.14,3x,"nsaos=",i0)') &
         i, wfn%emo(i)*evtoau, wfn%nao
      write(iunit,format) (wfn%C(j,i),j=1,wfn%nao)
   enddo
   write(iunit,'(a)') "$end"
end subroutine write_tm_mos
