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
!$ use omp_lib
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_filetypes, only : fileType
   use xtb_mctc_symbols, only : toSymbol
   use xtb_type_calculator, only : TCalculator
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule, len
   use xtb_type_param, only : scc_parameter
   use xtb_type_restart, only : TRestart
   use xtb_io_writer, only : writeMolecule
   use xtb_mctc_systools
   use xtb_mctc_strings
   use xtb_setparam
   use xtb_readin
   use xtb_mctc_convert
   use xtb_fixparam
   use xtb_scanparam
   use xtb_sphereparam
   use xtb_metadynamic
   use xtb_constrainpot
   implicit none
   private

   public :: TOrcaCalculator, newOrcaCalculator


   type, extends(TCalculator) :: TOrcaCalculator
      type(qm_external) :: ext
   contains
      !> Perform single point calculation
      procedure :: singlepoint
      !> Calculate hessian
      procedure :: hessian
      !> Write informative printout
      procedure :: writeInfo
   end type TOrcaCalculator


contains

!-----------------------------------------------------
! Check the ORCA execution settings 
!-----------------------------------------------------
subroutine checkOrca(env, ext)
   character(len=*), parameter :: source = 'extern_orca_checkOrca'
      !! the name of the error producer routine
   
   type(TEnvironment), intent(inout) :: env
      !! Calculation environment to handle I/O stream and error log
   type(qm_external), intent(inout) :: ext
      !! Settings for the Orca calculation

   character(len=:),allocatable :: homedir,syspath
   character(len=:),allocatable :: line
   character(len=5) :: chdum
   logical :: exist
   logical :: chk_xyzfile
   logical :: chk_engrad
   integer :: iorca 
      !! file handle
   integer :: err
   integer :: idx_bang,idx_star,idx_engr,idx_file,idx_xmol,idx_spce
!$ integer,external :: omp_get_num_threads

   !> check if the ORCA executable exist
   if (allocated(ext%executable)) then 
      !! user input
      
      if (ext%executable(1:1).eq.'~') then
         !! this is relative to the users home, expand it
         call rdvar('HOME',homedir) 
            !! read environment variable  
         ext%executable = homedir // ext%executable(2:)
      endif
      inquire(file=ext%executable,exist=exist)
      if (.not.exist) then
         call env%error("'"//ext%executable//"' was not found, please check",source)
         return
      endif
   
   else 
      !! no executable provided, lets find it 
      call rdvar('PATH',syspath)
      call rdpath(syspath,'orca',ext%executable,exist)
         !! find ORCA excutable in PATH variable
      if (.not.exist) then
         call env%error('Could not locate orca executable',source)
         return
      endif
   endif
   
   !> check for wrong executable with the same name
   if (index(ext%executable,'/usr') == 1) then
      inquire(file=ext%executable//'_scf',exist=exist)
      if (.not. exist) then
         call env%error("Executable '"//ext%executable//"' is in the system path, "//&
            & "is this the GNOME screen reader?")
         return
      endif
   endif

   !> see if there is a preference for an input file
   if (allocated(ext%input_file)) then
      inquire(file=ext%input_file,exist=exist)
      ext%exist = exist
   else
      ext%input_file = "orca-"//get_random_name()//'.inp'
         !! Create random file for the ORCA input file
      inquire(file=ext%input_file,exist=exist)
   endif

   !> sanity check
   if (ext%exist) then
      
      call open_file(iorca,ext%input_file,'r')
      !> if exists, but empty
      if (iorca.eq.-1) then
         call env%error("ORCA input file '"//ext%input_file//"' just vanished!",source)
         return
      endif

      chk_engrad = .false.
      chk_xyzfile = .false.
      do
         call strip_line(iorca,line,err)
            !! to read line from iorca unit
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
            ext%input_string = trim(adjustl(line(idx_spce:)))
               !! calculation method
         endif
      enddo
      call close_file(iorca)
      if (.not.(chk_engrad.and.chk_xyzfile)) then
         call env%error("Please add '! ENGRAD' and/or '* xyzfile' to '"//&
         & ext%input_file //"'!",source)
         return
      endif

   else
      !! if not exist, check for the input line
      if (allocated(ext%input_string)) then
         if (index(ext%input_string,'!') == 1) then
            ext%input_string = ext%input_string(2:)
         endif
      else
         ! general input
         ext%input_string = 'b97-3c'
      endif
   endif

end subroutine checkOrca

!-----------------------------------------------------
! Create a new calculator for driving the Orca program
!-----------------------------------------------------
subroutine newOrcaCalculator(self, env, ext,oniom)
   
   implicit none
   type(TOrcaCalculator), intent(out) :: self
      !! Instance of the Orca calculator
   type(TEnvironment), intent(inout) :: env
      !! Calculation environment
   type(qm_external), intent(in) :: ext
      !! Settings for the Orca calculator
   logical, intent(in),optional :: oniom

   self%ext = ext
      !! to save external settings from TSet%ext_orca to the new caluclator
   self%threadsafe = .false.
      !! not to call orca in parallel
   if (present(oniom)) self%ext%oniom=oniom
      !! if oniom calc
   call checkOrca(env, self%ext)


end subroutine newOrcaCalculator

!-----------------------------------------------------
! Create *.inp for the ORCA run
!-----------------------------------------------------
subroutine writeOrcaInp(io,mol,input,mode)
   integer, intent(in) :: io
   type(TMolecule), intent(in) :: mol
   character(len=*), intent(in) :: input
   character(len=*), intent(in) :: mode
   integer :: i, nel, mult

   write(io,'("#",1x,a)') &
      "orca input generated by xtb, this is not the right way of using orca!"
   !$omp parallel
   !$omp master
   !$ write(io,'("%",a,/,3x,a,1x,i0,/,a)') &
   !$    "pal","nprocs",omp_get_num_threads()-1,"end"
   !$omp end master
   !$omp end parallel
   write(io,'("%",a,1x,i0)') &
      "MaxCore",3000 ! hard coded, might be replaced at some point
   write(io,'("!",1x,a)') &
      "nopop", "miniprint", input, mode
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

   nel = sum(mol%at) - nint(mol%chrg)
   mult = mod(nel, 2) + 1
   if (mod(nel, 2) == mod(mol%uhf, 2)) mult = mol%uhf + 1

   write(io,'("*",1x,a,1x,i0,1x,i0)') &
      "xyz",nint(mol%chrg),mult
   do i = 1, mol%n
      write(io,'(3x,a2,3(2x,F20.14))') toSymbol(mol%at(i)), mol%xyz(:,i)*autoaa
   end do
   write(io,'("*",/)')

end subroutine writeOrcaInp


subroutine readOrcaEngrad(io, energy, gradient)
   integer, intent(in) :: io
   real(wp), intent(out) :: energy
   real(wp), intent(out) :: gradient(:, :)

   integer :: iat, nat

   read(io,'(a)')
   read(io,'(a)')
   read(io,'(a)')
   read(io,*) nat
   read(io,'(a)')
   read(io,'(a)')
   read(io,'(a)')
   read(io,*) energy
   read(io,'(a)')
   read(io,'(a)')
   read(io,'(a)')
   do iat = 1, nat
      read(io,*)gradient(1,iat)
      read(io,*)gradient(2,iat)
      read(io,*)gradient(3,iat)
   end do
end subroutine readOrcaEngrad

subroutine readOrcaHess(io, hess, dipgrad)
   integer, intent(in) :: io
   real(wp), intent(out) :: hess(:, :)
   real(wp), intent(out) :: dipgrad(:, :)

   integer :: ndim, irow, ii, stat
   character(len=:), allocatable :: line
   real(wp), allocatable :: buffer(:)
   integer, allocatable :: icol(:)

   stat = 0
   call getline(io, line, stat)
   do while(stat == 0)
      if (line == "$hessian") then
         call getline(io, line, stat)
         read(line, *) ndim
         allocate(icol(ndim), buffer(ndim))
         call getline(io, line, stat)
         do while(len(line) > 0 .and. stat == 0)
            icol(:) = -1
            read(line, *, iostat=stat) icol
            icol = pack(icol, icol >= 0) + 1
            do irow = 1, ndim
               read(io, *) ii, buffer(:size(icol))
               hess(icol, irow) = buffer(:size(icol))
            end do
            call getline(io, line, stat)
         end do
      end if

      if (line == "$dipole_derivatives") then
        call getline(io, line, stat)
         read(line, *) ndim
         do ii = 1, ndim
            read(io, *) dipgrad(:, ii)
         end do
      end if
      call getline(io, line, stat)
   end do
end subroutine readOrcaHess

!-----------------------------------------------------
! To create or modify input file and execute the ORCA
!-----------------------------------------------------
subroutine runOrca(env,ext,mol,energy,gradient)

   character(len=*), parameter :: source = 'extern_orca_runOrca'
      !! the name of the error producer routine
   type(TEnvironment), intent(inout) :: env
      !! Calculation environment to handle I/O stream and error log
   type(qm_external), intent(in) :: ext
      !! Settings for the Orca calculation
   type(TMolecule), intent(in) :: mol
      !! Molecular structure data 
   real(wp),intent(out) :: energy
   real(wp),intent(out) :: gradient(:, :)

   integer :: i,j,err
   integer :: iorca 
      !! file IO unit
   logical :: exist
   character(len=:),allocatable :: line
   character(len=:),allocatable :: outfile
   character(len=:),allocatable :: tmpfile

   !$omp critical (orca_lock)
   
   !> To decide whether to create or modify the input_file
   if (ext%exist) then
      !! we dump the name of the external xyz file to input_string... not cool
      call open_file(iorca,ext%input_string,'w')
      call writeMolecule(mol, iorca, format=fileType%xyz)
      call close_file(iorca)
   else
      !! create new.inp file 
      call open_file(iorca,ext%input_file,'w')
      call writeOrcaInp(iorca,mol,ext%input_string, "engrad")
      call close_file(iorca)
   endif
   !> Actual orca cml run 
   write(env%unit,'(72("="))')
   write(env%unit,'(1x,"*",1x,a)') &
      "letting orca take over the control..."
   if (set%oniom_settings%silent) then 
      call execute_command_line('exec 2>&1 '//ext%executable//' '// &
                             ext%input_file//'>orca.out',exitstat=err)
   else
      call execute_command_line('exec 2>&1 '//ext%executable//' '// &
                             ext%input_file,exitstat=err)
   endif

   if (err.ne.0) then
      call env%error('orca returned with non-zero exit status, doing the same',source)
   else
      write(env%unit,'(1x,"*",1x,a)') &
         "successful orca run, taking over control again..."
   endif
   write(env%unit,'(72("="))')
   
   !> find output .engrad file
   i = index(ext%input_file,'.inp')
   
   if (i > 0) then
      outfile = ext%input_file(:i-1)//'.engrad'
      tmpfile = ext%input_file(:i-1)
   else
      outfile = ext%input_file//'.engrad'
      tmpfile = ext%input_file
   endif
   inquire(file=outfile,exist=exist)
   if (.not.exist) then
      call env%error("Could not find '"//outfile//"', aborting driver run",source)
   else
      call open_file(iorca,outfile,'r')
      call readOrcaEngrad(iorca, energy, gradient)
      if (set%ceasefiles) then
         call remove_file(iorca)
         call env%io%deleteFile(tmpfile//".gbw")
         call env%io%deleteFile(tmpfile//".densities")
         call env%io%deleteFile(tmpfile//".rr")
         call env%io%deleteFile(tmpfile//"_property.txt")
         call env%io%deleteFile("orca.out")
      else
         call close_file(iorca)
      endif
   endif

   !$omp end critical (orca_lock)

end subroutine runOrca


subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
      & energy, gradient, sigma, hlgap, results)
   !> Source of the generated errors
   character(len=*), parameter :: source = 'extern_orca_singlepoint'
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

   integer :: i,ich
   integer :: mode_sp_run = 1
   real(wp) :: efix
   real(wp) :: dipole(3)
   logical :: inmol
   logical, parameter :: ccm = .true.
   logical :: exitRun
   character(len=*),parameter :: outfmt = &
      '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'

   call mol%update

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   hlgap = 0.0_wp
   efix = 0.0_wp
   dipole(:) = 0.0_wp

   call runOrca(env,self%ext,mol,energy,gradient)

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Electronic structure method terminated", source)
      return
   end if

   ! ------------------------------------------------------------------------
   !  various external potentials
   call constrain_pot(potset,mol%n,mol%at,mol%xyz,gradient,efix)
   call constrpot   (mol%n,mol%at,mol%xyz,gradient,efix)
   call cavity_egrad(mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (metaset,mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (rmsdset,mol%n,mol%at,mol%xyz,efix,gradient)

   ! ------------------------------------------------------------------------
   !  fixing of certain atoms
   !  print*,abs(efix/etot)
   energy = energy + efix
   results%e_total = energy
   results%gnorm = norm2(gradient)
   results%dipole = dipole
   if (fixset%n.gt.0) then
      do i=1, fixset%n
         !print*,i,fixset%atoms(i)
         gradient(1:3,fixset%atoms(i))=0
      enddo
   endif

   if (printlevel.ge.2) then
      ! start with summary header
      if (.not.set%silent) then
         write(env%unit,'(9x,53(":"))')
         write(env%unit,'(9x,"::",21x,a,21x,"::")') "SUMMARY"
      endif
      write(env%unit,'(9x,53(":"))')
      write(env%unit,outfmt) "total energy      ", results%e_total,"Eh   "
      write(env%unit,outfmt) "gradient norm     ", results%gnorm,  "Eh/a0"
      write(env%unit,outfmt) "HOMO-LUMO gap     ", results%hl_gap, "eV   "
      write(env%unit,'(9x,53(":"))')
      write(env%unit,'(a)')
   endif
end subroutine singlepoint

!> Evaluate hessian by finite difference for all atoms
subroutine hessian(self, env, mol0, chk0, list, step, hess, dipgrad)
   character(len=*), parameter :: source = "extern_turbomole_hessian"
   !> Single point calculator
   class(TOrcaCalculator), intent(inout) :: self
   !> Computation environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol0
   !> Restart data
   type(TRestart), intent(in) :: chk0
   !> List of atoms to displace
   integer, intent(in) :: list(:)
   !> Step size for numerical differentiation
   real(wp), intent(in) :: step
   !> Array to add Hessian to
   real(wp), intent(inout) :: hess(:, :)
   !> Array to add dipole gradient to
   real(wp), intent(inout) :: dipgrad(:, :)

   integer :: i,j,err
   integer :: iorca ! file handle
   logical :: exist
   character(len=:),allocatable :: line
   character(len=:),allocatable :: outfile
!$ integer,external :: omp_get_num_threads

   !$omp critical (orca_lock)
   if (self%ext%exist) then
      ! we dump the name of the external xyz file to input_string... not cool
      call open_file(iorca,self%ext%input_string,'w')
      call writeMolecule(mol0, iorca, format=fileType%xyz)
      call close_file(iorca)
   else
      call open_file(iorca,self%ext%input_file,'w')
      call writeOrcaInp(iorca,mol0,self%ext%input_string, "anfreq")
      call close_file(iorca)
   endif

   write(env%unit,'(72("="))')
   write(env%unit,'(1x,"*",1x,a)') &
      "letting orca take over the control..."
   call execute_command_line('exec 2>&1 '//self%ext%executable//' '// &
                             self%ext%input_file,exitstat=err)
   if (err.ne.0) then
      call env%error('orca returned with non-zero exit status, doing the same',source)
   else
      write(env%unit,'(1x,"*",1x,a)') &
         "successful orca run, taking over control again..."
   endif
   write(env%unit,'(72("="))')

   i = index(self%ext%input_file,'.inp')
   if (i > 0) then
      outfile = self%ext%input_file(:i-1)//'.hess'
   else
      outfile = self%ext%input_file//'.hess'
   endif
   inquire(file=outfile,exist=exist)
   if (.not.exist) then
      call env%error("Could not find '"//outfile//"', aborting driver run",source)
   else
      call open_file(iorca,outfile,'r')
      call readOrcaHess(iorca, hess, dipgrad)
      call close_file(iorca)
   endif
   !$omp end critical (orca_lock)

end subroutine hessian

subroutine writeInfo(self, unit, mol)

   !> Calculator instance
   class(TOrcaCalculator), intent(in) :: self

   !> Unit for I/O
   integer, intent(in) :: unit

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   call generic_header(unit,"Orca driver",49,10)

   write(unit,'(a,1x,a)') &
      "orca executable           :",self%ext%executable, &
      "orca input file           :",self%ext%input_file,&
      "orca input line           :",self%ext%input_string
   !$omp parallel
   !$omp master
   !$ write(unit,'(a,1x,g0)') &
   !$ "orca parallel             :",bool2string(omp_get_num_threads() > 1),&
   !$ "orca number of threads    :",omp_get_num_threads()
   !$omp end master
   !$omp end parallel
end subroutine writeInfo

function get_random_name() result(str)
   character(len=:), allocatable :: str
   real :: rnd
   integer :: irnd

   call random_number(rnd)
   irnd = transfer(rnd, irnd)
   str = to_hex(abs(irnd))
end function

pure function to_hex(val, width) result(string)
   integer, intent(in) :: val
   integer, intent(in), optional :: width
   character(len=:), allocatable :: string
   integer, parameter :: buffer_len = range(val)+2
   character(len=buffer_len) :: buffer
   integer :: pos
   integer :: n
   integer, parameter :: base = 16
   character(len=1), parameter :: numbers(0:base-1) = &
      ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f"]

   if (val == 0) then
      string = numbers(0)
      return
   end if

   n = abs(val)
   buffer = ""

   pos = buffer_len + 1
   do while (n > 0)
      pos = pos - 1
      buffer(pos:pos) = numbers(mod(n, base))
      n = n/base
   end do
   if (val < 0) then
      pos = pos - 1
      buffer(pos:pos) = '-'
   end if

   string = buffer(pos:)
end function to_hex

end module xtb_extern_orca
