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

module xtb_extern_orca
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_filetypes, only : fileType
   use xtb_mctc_symbols, only : toSymbol
   use xtb_type_calculator, only : TCalculator
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule, len
   use xtb_type_param, only : scc_parameter
   use xtb_type_restart, only : TRestart
   use xtb_type_wsc, only : tb_wsc
   use xtb_io_writer, only : writeMolecule
   use xtb_mctc_systools
   use xtb_mctc_strings
   use xtb_setparam
   use xtb_readin
   use xtb_mctc_convert
   implicit none
   private

   public :: checkOrca, runOrca


   type, extends(TCalculator) :: TOrcaCalculator
   contains
      !> Perform single point calculation
      procedure :: singlepoint

      !> Write informative printout
      procedure :: writeInfo
   end type TOrcaCalculator


contains


subroutine checkOrca(env)
   character(len=*), parameter :: source = 'extern_orca_checkOrca'
   type(TEnvironment), intent(inout) :: env
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
   if (allocated(set%ext_orca%executable)) then ! user input
      if (set%ext_orca%executable(1:1).eq.'~') then
         ! this is relative to the users home, expand it
         call rdvar('HOME',homedir)
         set%ext_orca%executable = homedir // set%ext_orca%executable(2:)
         if (set%verbose) then
            write(stdout,'(a,1x,a)') &
               "user home directory        :",homedir
         endif
      endif
      inquire(file=set%ext_orca%executable,exist=exist)
      if (.not.exist) then
         call env%error("'"//set%ext_orca%executable//"' was not found, please check",source)
         return
      endif
   else ! no executable provided, lets find it
      call rdvar('PATH',syspath)
      if (set%verbose) then
         write(stdout,'(a,1x,a)') &
            "system path                :",syspath
      endif
      call rdpath(syspath,'orca',set%ext_orca%executable,exist)
      if (.not.exist) then
         call env%error('Could not locate orca executable',source)
         return
      endif
   endif
   if (set%verbose) then
      write(stdout,'(a,1x,a)') &
         "orca executable           :",set%ext_orca%executable
      if (index(set%ext_orca%executable,'/usr') == 1) then
         write(stdout,'(a)') &
            "are you attempting to perform a calculation with the GNOME screen reader?"
      endif
   endif

   ! see if there is a preference for an input file
   if (allocated(set%ext_orca%input_file)) then
      inquire(file=set%ext_orca%input_file,exist=exist)
      set%ext_orca%exist = exist
   else
      set%ext_orca%input_file = 'orca.inp'
         inquire(file=set%ext_orca%input_file,exist=exist)
   endif
   if (set%verbose) then
      write(stdout,'(a,1x,a)') &
         "orca input file           :",set%ext_orca%input_file,&
         "orca input present        :",bool2string(exist),&
         "orca input override       :",bool2string(.not.set%ext_orca%exist)
   endif
   ! sanity check
   if (set%ext_orca%exist) then
      call open_file(iorca,set%ext_orca%input_file,'r')
      if (iorca.eq.-1) then
         call env%error("ORCA input file '"//set%ext_orca%input_file//"' just vanished!",source)
         return
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
            set%ext_orca%input_string = trim(adjustl(line(idx_spce:)))
         endif
      enddo
      call close_file(iorca)
      if (.not.(chk_engrad.and.chk_xyzfile)) then
         call env%error("Please add '! ENGRAD' and/or '* xyzfile' to '"//&
         & set%ext_orca%input_file //"'!",source)
         return
      endif

   else
      ! check for the input line
      if (allocated(set%ext_orca%input_string)) then
         if (index(set%ext_orca%input_string,'!') == 1) then
            set%ext_orca%input_string = set%ext_orca%input_string(2:)
         endif
      else
         ! general input
         set%ext_orca%input_string = 'b97-3c'
      endif
   endif
   if (set%verbose) then
      write(stdout,'(a,1x,a)') &
      !$ "orca parallel             :",bool2string(omp_get_num_threads() > 1),&
      !$ "orca number of threads    :",omp_get_num_threads(),&
         "orca input line           :",set%ext_orca%input_string
   endif

end subroutine checkOrca


subroutine runOrca(env,mol,energy,gradient)
   character(len=*), parameter :: source = 'extern_orca_runOrca'
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol
   real(wp),intent(out) :: energy
   real(wp),intent(out) :: gradient(:, :)

   integer :: i,j,err
   integer :: iorca ! file handle
   logical :: exist
   character(len=:),allocatable :: line
   character(len=:),allocatable :: outfile
!$ integer,external :: omp_get_num_threads

   !$omp critical (orca_lock)
   if (set%ext_orca%exist) then
      ! we dump the name of the external xyz file to input_string... not cool
      call open_file(iorca,set%ext_orca%input_string,'w')
      call writeMolecule(mol, iorca, format=fileType%xyz)
      call close_file(iorca)
   else
      call open_file(iorca,set%ext_orca%input_file,'w')
      write(iorca,'("#",1x,a)') &
         "orca input generated by xtb, this is not the right way of using orca!"
      !$ write(iorca,'("%",a,/,3x,a,1x,i0,/,a)') &
      !$    "pal","nprocs",omp_get_num_threads(),"end"
      write(iorca,'("%",a,1x,i0)') &
         "MaxCore",3000 ! hard coded, might be replaced at some point
      write(iorca,'("!",1x,a)') &
      set%ext_orca%input_string
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
         "xyz",set%ichrg,set%nalphabeta+1
      do i = 1, len(mol)
         write(iorca,'(3x,a2,3(2x,F20.14))') toSymbol(mol%at(i)), &
            & mol%xyz(1,i)*autoaa,mol%xyz(2,i)*autoaa,mol%xyz(3,i)*autoaa
      enddo
      write(iorca,'("*",/)')
      call close_file(iorca)
   endif

   write(stdout,'(72("="))')
   write(stdout,'(1x,"*",1x,a)') &
      "letting orca take over the control..."
   call execute_command_line('exec 2>&1 '//set%ext_orca%executable//' '// &
                             set%ext_orca%input_file,exitstat=err)
   if (err.ne.0) then
      call env%error('orca returned with non-zero exit status, doing the same',source)
   else
      write(stdout,'(1x,"*",1x,a)') &
         "successful orca run, taking over control again..."
   endif
   write(stdout,'(72("="))')

   i = index(set%ext_orca%input_file,'.inp')
   if (i > 0) then
      outfile = set%ext_orca%input_file(:i-1)//'.engrad'
   else
      outfile = set%ext_orca%input_file//'.engrad'
   endif
   inquire(file=outfile,exist=exist)
   if (.not.exist) then
      call env%error("Could not find '"//outfile//"', aborting driver run",source)
   else
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
      do j=1,len(mol)
         read(iorca,*)gradient(1,j)
         read(iorca,*)gradient(2,j)
         read(iorca,*)gradient(3,j)
      enddo
      call close_file(iorca)
   endif
   !$omp end critical (orca_lock)

end subroutine runOrca


subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
      & energy, gradient, sigma, hlgap, results)

   !> Source of the generated errors
   character(len=*), parameter :: source = 'type_calculator_singlepoint'

   !> Calculator instance
   class(TOrcaCalculator), intent(inout) :: self

   !> Computational environment
   type(TEnvironment), intent(inout) :: env

   !> Molecular structure data
   type(TMolecule), intent(inout) :: mol

   !> Wavefunction data
   type(TRestart), intent(inout) :: chk

   !> Print level for IO
   integer, intent(in) :: printlevel

   !> Restart from previous results
   logical, intent(in) :: restart

   !> Total energy
   real(wp), intent(out) :: energy

   !> Molecular gradient
   real(wp), intent(out) :: gradient(:, :)

   !> Strain derivatives
   real(wp), intent(out) :: sigma(:, :)

   !> HOMO-LUMO gap
   real(wp), intent(out) :: hlgap

   !> Detailed results
   type(scc_results), intent(out) :: results

end subroutine singlepoint

subroutine writeInfo(self, unit, mol)

   !> Calculator instance
   class(TOrcaCalculator), intent(in) :: self

   !> Unit for I/O
   integer, intent(in) :: unit

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

end subroutine writeInfo

end module xtb_extern_orca
