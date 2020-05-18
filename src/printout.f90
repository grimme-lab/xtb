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

module xtb_printout
   use xtb_mctc_accuracy, only : wp, sp
   implicit none

contains

!! ---------------------------------------------------------------[FB1808]-
subroutine writecosmofile(np,pa,espe,fname,nat,at,xyz,atom_weight)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autoaa
   use xtb_mctc_symbols, only : toLcSymbol, toSymbol
   implicit none
   integer, intent(in)             :: np
   real(wp), intent(in)            :: pa(3,np)
   real(wp), intent(in)            :: espe(np)
   integer, intent(in)             :: at(nat)
   real(wp), intent(in)            :: xyz(3,nat)
   integer, intent(in)             :: nat
   character(len=*),intent(in)     :: fname
   real(wp),intent(in)             :: atom_weight(2,np)
   logical                         :: exist
   integer                         :: id, i

   call open_file(id,fname,'w')
   write(id,'(a)') '$coord_car'
   write(id,'(a,/,a)') '!BIOSYM archive 3','coordinates from COSMO calculation'
   do i=1,nat
      write(id,'("X1",1x,3f22.14,1x,"COSM 1",1x,a,1x,a,1x,"0.000")')&
         xyz(:,i)*autoaa,toLcSymbol(at(i)),toSymbol(at(i))
   enddo
   write(id,'(a)') 'end'
   write(id,'(a)') '$segment_information'
   do i=1,np
      write(id,'(2x,i5,2x,i0,4f22.14,1x,f22.14,1x,f22.14,1x,"0.000")')&
         i,int(atom_weight(1,i)), pa(:,i)*autoaa,espe(i)/10, &
         atom_weight(2,i)*100,espe(i)/10
   enddo
   call close_file(id)

end subroutine writecosmofile

subroutine setup_summary(iunit,n,fname,xcontrol,wfx,xrc,exist)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_global, only : persistentEnv
   use xtb_mctc_systools
   use xtb_type_wavefunction
   use xtb_setparam
   !$ use omp_lib
   implicit none
   type(TWavefunction),intent(in) :: wfx
   integer, intent(in) :: iunit
   character(len=*),intent(in) :: xrc
   character(len=*),intent(in) :: xcontrol
   character(len=*),intent(in) :: fname
   integer :: i,l,err
   integer,intent(in) :: n
   logical,intent(in) :: exist
   real(wp) :: dum5
   character(len=:),allocatable :: cdum
   write(iunit,'(a)')
   call generic_header(iunit,'Calculation Setup',49,10)
   write(iunit,'(a)')
   if (allocated(cdum)) deallocate(cdum)
   call get_command(length=l)
   allocate( character(len=l) :: cdum )
   call get_command(cdum)
   write(iunit,'(10x,a,":",1x,a)') 'program call               ',cdum
   call rdvar('HOSTNAME',cdum,err)
   if (err.eq.0) &
      write(iunit,'(10x,a,":",1x,a)') 'hostname                   ',cdum
   if (allocated(persistentEnv%io%namespace)) &
      write(iunit,'(10x,a,":",1x,a)') 'calculation namespace      ',persistentEnv%io%namespace
   ! ----------------------------------------------------------------------
   !  print the home and path to check if there are set correctly
   write(iunit,'(10x,a,":",1x,a)') 'coordinate file            ',fname
   if (verbose) then
      write(iunit,'(10x,a,":",1x,a)') 'xtbhome directory          ',xenv%home
      write(iunit,'(10x,a,":",1x,a)') 'path for xtb               ',xenv%path
      write(iunit,'(10x,a,":",1x,a)') 'xcontrol input file        ',xcontrol
      if (exist) &
         write(iunit,'(10x,a,":",1x,a)') 'global configurations file ',xrc
   endif
   ! ----------------------------------------------------------------------
   !  technical data
   !$omp parallel
   !$omp master
   !$ write(iunit,'(10x,a,":",6x,i16)') 'omp threads                ',omp_get_num_threads()
   !$omp end master
   !$omp end parallel
   ! ----------------------------------------------------------------------
   !  print more specific calculation data
   write(iunit,'(10x,a,":",6x,i16)')   'number of atoms            ',n
   write(iunit,'(10x,a,":",6x,i16)')   'number of electrons        ',wfx%nel
   write(iunit,'(10x,a,":",6x,i16)')   'charge                     ',ichrg
   write(iunit,'(10x,a,":",6x,f16.1)') 'spin                       ',0.5_wp*wfx%nopen
   call random_number(dum5)
   write(iunit,'(10x,a,":",6x,f16.14)') 'first test random number   ',dum5
   if (verbose) then
      write(iunit,'(10x,a,":",12x,"0x",z8)') 'a pointer address          ',real(dum5,sp)
      write(iunit,'(10x,a,":",6x,z16)') 'random memory content      ',dum5
   endif
   if (veryverbose) then
      write(iunit,'(10x,a,":",20x,a)') 'is this your card?         ',"🃓"
      write(iunit,'(10x,a,":",15x,"so true")')'this was released?         '
   endif
   write(iunit,'(a)')

end subroutine

end module xtb_printout
