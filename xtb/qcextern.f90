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

module qcextern
   use iso_fortran_env, only : wp => real64, istdout => output_unit
   implicit none

contains

subroutine mopac_chk()
   use mctc_systools
   use setparam
   use readin
   implicit none
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
            write(istdout,'(a,1x,a)') &
               "user home directory        :",homedir
         endif
      endif
      inquire(file=ext_mopac%executable,exist=exist)
      if (.not.exist) then
         call raise('E',"'"//ext_mopac%executable//"' was not found, please check",1)
      endif
   else ! no executable provided, lets find it
      call rdvar('PATH',syspath)
      if (verbose) then
         write(istdout,'(a,1x,a)') &
            "system path                :",syspath
      endif
      call rdpath(syspath,'mopac',ext_mopac%executable,exist)
      if (.not.exist) then
         call raise('E','Could not locate mopac executable',1)
      endif
   endif
   if (verbose) then
      write(istdout,'(a,1x,a)') &
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
      write(istdout,'(a,1x,a)') &
         "mopac input file           :",ext_mopac%input_file,&
         "mopac input present        :",bool2string(exist)
         !"mopac input override       :",bool2string(.not.ext_mopac%exist)
   endif

   ! check for the input line
   if (allocated(ext_mopac%input_string)) then
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
      write(istdout,'(a,1x,a)') &
         "mopac input line           :",ext_mopac%input_string
   endif

end subroutine mopac_chk

subroutine run_mopac_egrad(nat,at,xyz,energy,gradient)
   use mctc_econv
   use mctc_systools
   use setparam
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: energy
   real(wp),intent(out) :: gradient(3,nat)

   integer :: i,j,err
   integer :: imopac ! file handle
   logical :: exist
   character(len=:),allocatable :: line
   character(len=2),external :: asym
   integer :: num
   real(wp) :: dum(10),edum

   !$omp critical(mopac_lock)
   call open_file(imopac,ext_mopac%input_file,'w')
   write(imopac,'(a,/,/)') ext_mopac%input_string
   do i = 1, nat
      write(imopac,'(3x,a2,3(f20.14,i5))') &
         asym(at(i)),autoaa*xyz(1,i),1,autoaa*xyz(2,i),1,autoaa*xyz(3,i),1
   enddo
   call close_file(imopac)

   write(istdout,'(72("="))')
   write(istdout,'(1x,"*",1x,a)') &
      "handing control over to mopac..."
   call execute_command_line('exec 2>&1 '//ext_mopac%executable//' '// &
                             ext_mopac%input_file,exitstat=err)
   if (err.ne.0) then
      call raise('E','mopac returned with non-zero exit status, following this',1)
   endif
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

   write(istdout,'(1x,"*",1x,a)') &
      "regaining control after successful mopac run..."
   write(istdout,'(72("="))')

   call open_file(imopac,ext_mopac%input_file//'.aux','r')
   if (imopac.eq.-1) then
      call raise('E',"Could not find '"//ext_mopac%input_file//".aux'",1)
   endif
   read_mopac_output: do
      call getline(imopac,line,iostat=err)
      if (is_iostat_end(err)) then
         call raise('E','Could not find gradient in mopac output',1)
      endif
      if (index(line,'TOTAL_ENERGY:EV') > 0) then
         call readl(line,dum,num)
         energy = dum(num)*evtoau
         cycle read_mopac_output
      endif
      if (index(line,'GRADIENTS:KCAL/MOL/ANGSTROM') > 0) then
         read(imopac,*)((gradient(j,i),j=1,3),i=1,nat)
         exit read_mopac_output
      endif
   enddo read_mopac_output
   call close_file(imopac)
   !$omp end critical (mopac_lock)
   gradient = gradient * kcaltoau / aatoau
end subroutine run_mopac_egrad

subroutine orca_chk()
   use mctc_systools
   use mctc_strings
   use setparam
   use readin
   implicit none
   character(len=:),allocatable :: homedir,syspath
   character(len=:),allocatable :: line
   character(len=5) :: chdum
   logical :: exist
   logical :: chk_xyzfile
   logical :: chk_engrad
   integer :: iorca ! file handle
   integer :: err
   integer :: idx_bang,idx_star,idx_engr,idx_file,idx_xmol,idx_spce
!$ integer,external :: omp_get_num_threads

   ! check for ORCA executable
   if (allocated(ext_orca%executable)) then ! user input
      if (ext_orca%executable(1:1).eq.'~') then
         ! this is relative to the users home, expand it
         call rdvar('HOME',homedir)
         ext_orca%executable = homedir // ext_orca%executable(2:)
         if (verbose) then
            write(istdout,'(a,1x,a)') &
               "user home directory        :",homedir
         endif
      endif
      inquire(file=ext_orca%executable,exist=exist)
      if (.not.exist) then
         call raise('E',"'"//ext_orca%executable//"' was not found, please check",1)
      endif
   else ! no executable provided, lets find it
      call rdvar('PATH',syspath)
      if (verbose) then
         write(istdout,'(a,1x,a)') &
            "system path                :",syspath
      endif
      call rdpath(syspath,'orca',ext_orca%executable,exist)
      if (.not.exist) then
         call raise('E','Could not locate orca executable',1)
      endif
   endif
   if (verbose) then
      write(istdout,'(a,1x,a)') &
         "orca executable           :",ext_orca%executable
      if (index(ext_orca%executable,'/usr') == 1) then
         write(istdout,'(a)') &
            "are you attempting to perform a calculation with the GNOME screen reader?"
      endif
   endif

   ! see if there is a preference for an input file
   if (allocated(ext_orca%input_file)) then
      inquire(file=ext_orca%input_file,exist=exist)
      ext_orca%exist = exist
   else
      ext_orca%input_file = 'orca.inp'
         inquire(file=ext_orca%input_file,exist=exist)
   endif
   if (verbose) then
      write(istdout,'(a,1x,a)') &
         "orca input file           :",ext_orca%input_file,&
         "orca input present        :",bool2string(exist),&
         "orca input override       :",bool2string(.not.ext_orca%exist)
   endif
   ! sanity check
   if (ext_orca%exist) then
      call open_file(iorca,ext_orca%input_file,'r')
      if (iorca.eq.-1) then
         call raise('E',"ORCA input file '"//ext_orca%input_file//"' just vanished!",1)
      endif
      chk_engrad = .false.
      chk_xyzfile = .false.
      do
         call strip_line(iorca,line,err)
         if (err.ne.0) exit
         idx_bang = index(line,'!')
         idx_star = index(line,'*')
         idx_engr = index(lowercase(line),'engrad')
         idx_file = index(lowercase(line),'xyzfile')
         idx_xmol = index(line,'.xyz')
         idx_spce = index(trim(line),' ',back=.true.)
         if (idx_bang.gt.0 .and. idx_engr.gt.0 .and. idx_engr.gt.idx_bang) &
            chk_engrad = .true.
         if (idx_star.gt.0 .and. idx_file.gt.0 .and. idx_xmol.gt.0 .and. &
            &idx_file.gt.idx_star .and. idx_xmol.gt.idx_file) then
            chk_xyzfile = .true.
            ext_orca%input_string = trim(adjustl(line(idx_spce:)))
         endif
      enddo
      call close_file(iorca)
      if (.not.(chk_engrad.and.chk_xyzfile)) then
         call raise('E',"Please add '! ENGRAD' and/or '* xyzfile' to '"//&
         & ext_orca%input_file //"'!",1)
      endif

   else
      ! check for the input line
      if (allocated(ext_orca%input_string)) then
         if (index(ext_orca%input_string,'!') == 1) then
            ext_orca%input_string = ext_orca%input_string(2:)
         endif
      else
         ! general input
         ext_orca%input_string = 'b97-3c'
      endif
   endif
   if (verbose) then
      write(istdout,'(a,1x,a)') &
      !$ "orca parallel             :",bool2string(omp_get_num_threads() > 1),&
      !$ "orca number of threads    :",omp_get_num_threads(),&
         "orca input line           :",ext_orca%input_string
   endif

end subroutine orca_chk

subroutine run_orca_egrad(nat,at,xyz,energy,gradient)
   use mctc_econv
   use setparam
   use write_geometry
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: energy
   real(wp),intent(out) :: gradient(3,nat)

   integer :: i,j,err
   integer :: iorca ! file handle
   logical :: exist
   character(len=:),allocatable :: line
   character(len=:),allocatable :: outfile
   character(len=2),external :: asym
!$ integer,external :: omp_get_num_threads

   !$omp critical (orca_lock)
   if (ext_orca%exist) then
      ! we dump the name of the external xyz file to input_string... not cool
      call open_file(iorca,ext_orca%input_string,'w')
      call write_xyz(iorca,nat,at,xyz)
      call close_file(iorca)
   else
      call open_file(iorca,ext_orca%input_file,'w')
      write(iorca,'("#",1x,a)') &
         "orca input generated by xtb, this is not the right way of using orca!"
      !$ write(iorca,'("%",a,/,3x,a,1x,i0,/,a)') &
      !$    "pal","nprocs",omp_get_num_threads(),"end"
      write(iorca,'("%",a,1x,i0)') &
         "MaxCore",3000 ! hard coded, might be replaced at some point
      write(iorca,'("!",1x,a)') &
      ext_orca%input_string
!      if(index(solv,'h2o').ne.0)   solv='water'
!      if(index(solv,'chcl3').ne.0) solv='chloroform'
!      if(index(solv,'ether').ne.0) solv='diethyl ether'
!      if(smd)then
!         write(ich,'(''! cpcm('',a,'')'')')trim(solv)
!         write(ich,'(''%cpcm'')')
!         write(ich,'('' smd     true'')')
!         write(ich,'('' solvent '',''"'',a,''"'')')trim(solv)
!         write(ich,'(''end'')')
!      endif
 
      write(iorca,'("%",a,/,3x,a,1x,a,/,a)') &
         "method","runtyp","gradient","end"
      write(iorca,'("*",1x,a,1x,i0,1x,i0)') &
         "xyz",ichrg,nalphabeta+1
      do i = 1, nat
         write(iorca,'(3x,a2,3(2x,F20.14))') &
            asym(at(i)),xyz(1,i)*autoaa,xyz(2,i)*autoaa,xyz(3,i)*autoaa
      enddo
      write(iorca,'("*",/)')
      call close_file(iorca)
   endif

   write(istdout,'(72("="))')
   write(istdout,'(1x,"*",1x,a)') &
      "letting orca take over the control..."
   call execute_command_line('exec 2>&1 '//ext_orca%executable//' '// &
                             ext_orca%input_file,exitstat=err)
   if (err.ne.0) then
      call raise('E','orca returned with non-zero exit status, doing the same',1)
   endif
   write(istdout,'(1x,"*",1x,a)') &
      "successful orca run, taking over control again..."
   write(istdout,'(72("="))')

   i = index(ext_orca%input_file,'.inp')
   if (i > 0) then
      outfile = ext_orca%input_file(:i-1)//'.engrad'
   else
      outfile = ext_orca%input_file//'.engrad'
   endif
   inquire(file=outfile,exist=exist)
   if (.not.exist) then
      call raise('E',"Could not find '"//outfile//"', aborting driver run",1)
   endif
   call open_file(iorca,outfile,'r')
   read(iorca,'(a)')
   read(iorca,'(a)')
   read(iorca,'(a)')
   read(iorca,*) i
   read(iorca,'(a)')
   read(iorca,'(a)')
   read(iorca,'(a)')
   read(iorca,*) energy
   read(iorca,'(a)')
   read(iorca,'(a)')
   read(iorca,'(a)')
   do j=1,nat
      read(iorca,*)gradient(1,j)
      read(iorca,*)gradient(2,j)
      read(iorca,*)gradient(3,j)
   enddo
   call close_file(iorca)
   !$omp end critical (orca_lock)

end subroutine run_orca_egrad

end module qcextern
