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

module xtb_extern_mopac
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_symbols, only : toSymbol
   use xtb_type_environment, only : TEnvironment
   use xtb_mctc_systools
   use xtb_setparam
   use xtb_readin
   use xtb_mctc_convert
   use xtb_mctc_systools
   use xtb_setparam
   implicit none
   private

   public :: checkMopac, runMopac


contains


subroutine checkMopac(env)
   character(len=*), parameter :: source = 'extern_mopac_checkMopac'
   type(TEnvironment), intent(inout) :: env
   character(len=:),allocatable :: homedir,syspath
   character(len=5) :: chdum
   integer :: i
   logical :: exist

   ! check for MOPAC executable
   if (allocated(ext_mopac%executable)) then ! user input
      if (ext_mopac%executable(1:1).eq.'~') then
         ! this is relative to the users home, expand it
         call rdvar('HOME',homedir)
         ext_mopac%executable = homedir // ext_mopac%executable(2:)
         if (verbose) then
            write(stdout,'(a,1x,a)') &
               "user home directory        :",homedir
         endif
      endif
      inquire(file=ext_mopac%executable,exist=exist)
      if (.not.exist) then
         call env%error("'"//ext_mopac%executable//"' was not found, please check",source)
         return
      endif
   else ! no executable provided, lets find it
      call rdvar('PATH',syspath)
      if (verbose) then
         write(stdout,'(a,1x,a)') &
            "system path                :",syspath
      endif
      call rdpath(syspath,'mopac',ext_mopac%executable,exist)
      if (.not.exist) then
         call env%error('Could not locate mopac executable',source)
         return
      endif
   endif
   if (verbose) then
      write(stdout,'(a,1x,a)') &
         "mopac executable           :",ext_mopac%executable
   endif

   ! see if there is a preference for an input file
   if (allocated(ext_mopac%input_file)) then
      inquire(file=ext_mopac%input_file,exist=exist)
      ext_mopac%exist = exist
   else
      ext_mopac%input_file = 'GCOORD'
      inquire(file=ext_mopac%input_file,exist=exist)
   endif
   if (verbose) then
      write(stdout,'(a,1x,a)') &
         "mopac input file           :",ext_mopac%input_file,&
         "mopac input present        :",bool2string(exist)
         !"mopac input override       :",bool2string(.not.ext_mopac%exist)
   endif

   ! check for the input line
   if (allocated(ext_mopac%input_string)) then
      if (index(ext_mopac%input_string,'grad') == 0) then
         call env%warning('added grad keyword to mopac input', source)
         ext_mopac%input_string = ext_mopac%input_string //' grad'
      end if
      if (index(ext_mopac%input_string,'charge=') == 0) then
         write(chdum,'(i5)') ichrg
         ! add total charge
         ext_mopac%input_string = ext_mopac%input_string //' charge='//trim(adjustl(chdum))
      endif
      if (nalphabeta > 0) then
         if (index(ext_mopac%input_string,'uhf') == 0) then
            ! write spin state if necessary
            select case(nalphabeta)
            case default ! skip
            case(1); ext_mopac%input_string = ext_mopac%input_string // ' uhf doublet'
            case(2); ext_mopac%input_string = ext_mopac%input_string // ' uhf triplet'
            case(3); ext_mopac%input_string = ext_mopac%input_string // ' uhf quartet'
            case(4); ext_mopac%input_string = ext_mopac%input_string // ' uhf quintet'
            case(5); ext_mopac%input_string = ext_mopac%input_string // ' uhf sextet'
            case(6); ext_mopac%input_string = ext_mopac%input_string // ' uhf septet'
            case(7); ext_mopac%input_string = ext_mopac%input_string // ' uhf octet'
            end select
         endif
      endif
   else
      ! general input
      ext_mopac%input_string = '1scf pm6-d3h4 aux(42,PRECISION=12,MOS=-99999,COMP)'
      ! write spin state if necessary
      select case(nalphabeta)
      case default ! skip
      case(1); ext_mopac%input_string = ext_mopac%input_string // ' uhf doublet'
      case(2); ext_mopac%input_string = ext_mopac%input_string // ' uhf triplet'
      case(3); ext_mopac%input_string = ext_mopac%input_string // ' uhf quartet'
      case(4); ext_mopac%input_string = ext_mopac%input_string // ' uhf quintet'
      case(5); ext_mopac%input_string = ext_mopac%input_string // ' uhf sextet'
      case(6); ext_mopac%input_string = ext_mopac%input_string // ' uhf septet'
      case(7); ext_mopac%input_string = ext_mopac%input_string // ' uhf octet'
      end select
      ! convergence criterium in kcal/mol
      ext_mopac%input_string = ext_mopac%input_string // &
         ' scfcrt=6.d-5 geo-ok mmok grad xyz charge='
      write(chdum,'(i5)') ichrg
      ! add total charge
      ext_mopac%input_string = ext_mopac%input_string // trim(adjustl(chdum))
   endif
   if (verbose) then
      write(stdout,'(a,1x,a)') &
         "mopac input line           :",ext_mopac%input_string
   endif

end subroutine checkMopac


subroutine runMopac(env,nat,at,xyz,energy,gradient,dipole)
   character(len=*), parameter :: source = 'extern_mopac_checkMopac'
   type(TEnvironment), intent(inout) :: env
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: energy
   real(wp),intent(out) :: gradient(3,nat)
   real(wp),intent(out) :: dipole(3)

   integer :: i,j,err
   integer :: imopac ! file handle
   logical :: exist
   character(len=:),allocatable :: line
   integer :: num
   real(wp) :: dum(10),edum

   !$omp critical(mopac_lock)
   call open_file(imopac,ext_mopac%input_file,'w')
   write(imopac,'(a,/,/)') ext_mopac%input_string
   do i = 1, nat
      write(imopac,'(3x,a2,3(f20.14,i5))') &
         toSymbol(at(i)),autoaa*xyz(1,i),1,autoaa*xyz(2,i),1,autoaa*xyz(3,i),1
   enddo
   call close_file(imopac)

   write(stdout,'(72("="))')
   write(stdout,'(1x,"*",1x,a)') &
      "handing control over to mopac..."
   call execute_command_line('exec 2>&1 '//ext_mopac%executable//' '// &
                             ext_mopac%input_file,exitstat=err)
   if (err.ne.0) then
      call env%error('mopac returned with non-zero exit status, following this',source)
   else
      if (verbose) then
         inquire(file=ext_mopac%input_file//'.arc',exist=exist)
         if (exist) then
            call open_file(imopac,ext_mopac%input_file//'.arc','r')
            print_mopac_output: do
               call getline(imopac,line,iostat=err)
               if (is_iostat_end(err)) exit print_mopac_output
            enddo print_mopac_output
            call close_file(imopac)
         endif
      endif

      write(stdout,'(1x,"*",1x,a)') &
         "regaining control after successful mopac run..."
   endif
   write(stdout,'(72("="))')

   call open_file(imopac,ext_mopac%input_file//'.aux','r')
   if (imopac.eq.-1) then
      call env%error("Could not find '"//ext_mopac%input_file//".aux'",source)
   else
      read_mopac_output: do
         call getline(imopac,line,iostat=err)
         if (is_iostat_end(err)) then
            call env%error('Could not find gradient in mopac output',source)
         endif
         if (index(line,'HEAT_OF_FORMATION:KCAL/MOL') > 0) then
            call readl(line,dum,num)
            energy = dum(num)*kcaltoau
            cycle read_mopac_output
         endif
         if (index(line,'DIP_VEC:DEBYE') > 0)then
            call readl(line,dum,num)
            dipole(1:3) = dum(2:4)*dtoau
            cycle read_mopac_output
         endif
         if (index(line,'GRADIENTS:KCAL/MOL/ANGSTROM') > 0) then
            read(imopac,*)((gradient(j,i),j=1,3),i=1,nat)
            exit read_mopac_output
         endif
      enddo read_mopac_output
      call close_file(imopac)
      gradient = gradient * kcaltoau / aatoau
   endif
   !$omp end critical (mopac_lock)

end subroutine runMopac


end module xtb_extern_mopac
