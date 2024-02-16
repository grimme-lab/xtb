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
      
      !> perform single-point calculation
      procedure :: singlepoint
      
      !> calculate hessian
      procedure :: hessian
      
      !> write informative printout
      procedure :: writeInfo
   
   end type TOrcaCalculator


contains

!> create a new calculator for driving ORCA
subroutine newOrcaCalculator(self, env, ext, oniom)
   
   !> ORCA calculator instance
   type(TOrcaCalculator), intent(out) :: self
   
   !> calculation environment
   type(TEnvironment), intent(inout) :: env
   
   !> ORCA settings
   type(qm_external), intent(in) :: ext
   
   !> if ONIOM calculation 
   logical, intent(in), optional :: oniom

   ! settings !
   self%ext = ext
   if (present(oniom)) self%ext%oniom=oniom
   self%threadsafe = .false. ! parallelization !

   ! setup ORCA !
   call checkOrca(env, self%ext)

end subroutine newOrcaCalculator

!> setup environment for ORCA run 
subroutine checkOrca(env, ext)

   use xtb_extern_driver, only : checkExe

   !> traceback 
   character(len=*), parameter :: source = 'extern_orca_checkOrca'
   
   !> calculation environment
   type(TEnvironment), intent(inout) :: env
   
   !> ORCA settings
   type(qm_external), intent(inout) :: ext


   !> file existence 
   logical :: exist
   
   !> XMol file  
   logical :: chk_xyzfile
   
   !> ENGRAD keyword
   logical :: chk_engrad
  
   !> ORCA input file unit
   integer :: iorca 
   
   ! find ORCA !
   call checkExe(env,ext%executable,"orca")
   
   ! check for wrong executable [GNOME screen reader] !
   if (index(ext%executable,'/usr') == 1) then
      inquire(file=ext%executable//'_scf',exist=exist)
      if (.not. exist) then
         call env%error("Executable '"//ext%executable//"' is in the system path, "//&
            & "is this the GNOME screen reader?", source)
         return
      endif
   endif
   
   ! preference for an input file !
   if (allocated(ext%input_file)) then
      inquire(file=ext%input_file,exist=exist)
      ext%exist = exist

   else ! assign input file name !
      ext%input_file = "orca-"//get_random_name()//'.inp'
      inquire(file=ext%input_file,exist=exist)  
      ext%exist = exist
   endif

   ! check existing input !
   if (ext%exist) then
      
      call open_file(iorca,ext%input_file,'r')
      
      ! sanity check !
      if (iorca.eq.-1) then 
         call env%error("ORCA input file '"//ext%input_file//"' just vanished!",source)
         return
      endif

      ! find engrad & xyzfile keys !
      call chkInput(iorca,chk_engrad,chk_xyzfile,ext%str_file)
      
      call close_file(iorca)

      if (.not.(chk_engrad.and.chk_xyzfile)) then
         call env%error("Please add '! ENGRAD' and/or '* xyzfile' to '"//&
         & ext%input_file //"'!",source)
         return
      endif
   
   else   ! input string: find || create !

      if (allocated(ext%input_string)) then
         if (index(ext%input_string,'!') == 1) then
            ext%input_string = ext%input_string(2:)
         endif
      else
         ext%input_string = 'b97-3c'
      endif
   
   endif

end subroutine checkOrca

!> check ORCA input file for certain keywords
subroutine chkInput(iunit,grad,xyz,geo)

   !> ORCA file unit
   integer, intent(in) :: iunit 
   
   !> if certain keywords are present
   logical, intent(out) :: grad, xyz

   !> buffer string
   character(len=:), allocatable :: line
   
   !> molecular structure file name
   character(len=:), allocatable :: geo

   !> error handling
   integer :: err

   !> position holders
   integer :: startPos, endPos
   
   ! Intilialization !
   grad = .false.
   xyz = .false.
   
   ! read file !
   do
      
      call strip_line(iunit,line,err)
      if (err.ne.0) exit
      
      line = lowercase(line)
      
      if (firstArgComesFirst(line,'!','engrad')) &
            grad = .true.

      if(firstArgComesFirst(line,'*','xyzfile')) then
         if (firstArgComesFirst(line,'xyzfile','.xyz')) then
            
            ! get the name of geometry file ! 
            endPos = index(line,'.xyz')
            startPos = index(line(:endPos),' ',back=.true.)
            geo = line(startPos+1:endPos+3)
            
            xyz = .true.
         endif
      endif

   enddo

end subroutine chkInput

!> check argument positions wrt. each other
function firstArgComesFirst(str,arg1,arg2)

   !> raw line
   character(len=*), intent(in) :: str

   !> arguments
   character(len=*), intent(in) :: arg1, arg2

   !> return value
   logical :: firstArgComesFirst
   
   !> arg positions
   integer :: pos1, pos2

   ! find arg positions !
   pos1 = index(str,arg1)
   pos2 = index(str,arg2)

   if (pos1 > 0 .and. pos2 > 0) then
      firstArgComesFirst = pos1 < pos2
   else 
      firstArgComesFirst = .false.
   endif 

end function

!> create *.inp for ORCA run
subroutine writeOrcaInp(env,io,mol,input,mode)
   
   !> calculation enironment
   type(TEnvironment), intent(inout) :: env

   !> ORCA input file unit
   integer, intent(in) :: io

   !> molecular structure data
   type(TMolecule), intent(in) :: mol

   !> input string 
   character(len=*), intent(in) :: input

   !> calculation type
   character(len=*), intent(in) :: mode
   
   integer :: i, nel, mult, etc
   integer :: num_threads

   write(io,'("#",1x,a)') &
      "ORCA input is generated automatically. Not correct way."
   
   ! number of cores !
   !$omp parallel
   !$omp master
      num_threads = omp_get_num_threads()
   !$omp end master
   !$omp end parallel 
   
   write(io,'("%",a,/,3x,a,1x,i0,/,a)') &
      "pal","nprocs", num_threads ,"end"
   
   ! memory usage !
   write(io,'("%",a,1x,i0)') &
      "MaxCore", get_mem(env,num_threads)
   
   ! orca keywords !
   write(io,'("!",1x,a)') &
      "nopop", "miniprint", input, mode

   ! calculate number of electrons and multiplicity !
   nel = sum(mol%at) - nint(mol%chrg)
   mult = mod(nel, 2) + 1
   if (mod(nel, 2) == mod(mol%uhf, 2)) mult = mol%uhf + 1

   ! print molecular geometry !
   write(io,'("*",1x,a,1x,i0,1x,i0)') &
      "xyz",nint(mol%chrg),mult

   do i = 1, mol%n
      write(io,'(3x,a2,3(2x,F20.14))') toSymbol(mol%at(i)), mol%xyz(:,i)*autoaa
   end do

   write(io,'("*",/)')

end subroutine writeOrcaInp

!> get available memory from /proc/meminfo
function get_mem(env,nThread)

   !> calculation environment
   type(TEnvironment), intent(inout) :: env 

   !> number of threads
   integer, intent(in) :: nThread

   !> available memory
   integer :: get_mem
   
   !> IO for 
   integer :: meminfo

   !> buffer string 
   character(len=:), allocatable :: line

   !> error handling
   integer :: err, mem

   !> array for readl
   real(wp) :: xx 

   call open_file(meminfo, '/proc/meminfo','r')
   
   ! default case !
   get_mem = 3000
   
   ! read from /proc/meminfo ! 
   if (meminfo.ne.-1) then

      call strip_line(meminfo,line,err) 
      call strip_line(meminfo,line,err)
      call strip_line(meminfo,line,err)
      if (index(line,"MemAvailable").ne.0) then
         call readl(line,xx,mem)
         get_mem = (int(xx)/1024)/nThread ! memory per core !
         get_mem = get_mem - 100 ! 100 MB that maybe needed for add jobs !
      endif

   endif

endfunction get_mem

!> read energy and gradients from ORCA calculation
subroutine readOrcaEngrad(io, energy, gradient)
   
   !> orca unit number
   integer, intent(in) :: io

   !> electronic energy
   real(wp), intent(out) :: energy

   !> gradients 
   real(wp), intent(out) :: gradient(:, :)

   !> number of atoms
   integer :: nat

   !> loop index
   integer :: iat
   
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

!>  wrapper for ORCA system call 
subroutine runOrca(env,ext,mol,energy,gradient)

   !> traceback
   character(len=*), parameter :: source = 'extern_orca_runOrca'
   
   !> calculation environment 
   type(TEnvironment), intent(inout) :: env
   
   !> ORCA settings 
   type(qm_external), intent(in) :: ext
   
   !> molecular structure data 
   type(TMolecule), intent(in) :: mol

   !> ORCA SP energy
   real(wp),intent(out) :: energy
   
   !> ORCA gradients
   real(wp),intent(out) :: gradient(:, :)

   integer :: i,j,err, err2
   !> file IO unit
   integer :: iorca, eunit
   logical :: exist
   character(len=:),allocatable :: line
   character(len=:),allocatable :: outfile
   character(len=:),allocatable :: tmpfile

   character(Len=:), allocatable :: message

   !$omp critical (orca_lock)
   
   ! input file handling !
   if (ext%exist) then
      
      ! rewrite *.xyz file from ORCA input !
      call open_file(iorca,ext%str_file,'w')
      call writeMolecule(mol, iorca, format=fileType%xyz)
      call close_file(iorca)

   else
      
      ! create ORCA input !
      call open_file(iorca,ext%input_file,'w')
      call writeOrcaInp(env,iorca,mol,ext%input_string, "engrad")
      call close_file(iorca)
   
   endif
   
   ! shift responsibility from xtb !  
   write(env%unit,'(25x,"*",1x,a)') &
      "ORCA takes the control..."
   write(env%unit,'(2x,72("-"))')

   ! ORCA system call !
   if (set%oniom_settings%silent) then ! silent mode ! 
      call execute_command_line('exec '//ext%executable//' '// &
                             ext%input_file//'>orca.out 2> err',exitstat=err)
   else
      call execute_command_line('exec '//ext%executable//' '// &
                             ext%input_file//' 2> err',exitstat=err)
   endif
   
   ! ORCA dumps mpirun errors into stderr !
   call open_file(eunit,'err','r')
   call getline(eunit,message,err2)

   ! shift responsibility back !
   write(env%unit,'(2x,72("-"))')
   if (err.ne.0 .or. .not. is_iostat_end(err2)) then
      call raise('E','ORCA returned with non-zero exit status, doing the same')
   else
      write(env%unit,'(10x,"*",1x,a)') &
         "successful ORCA run, taking over control again..."
      call env%io%deleteFile("err")
   endif
   
   ! find output .engrad file !
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
      call raise('E',"Could not find '"//outfile//"', aborting driver run")
  
   else

      ! read .engrad !
      call open_file(iorca,outfile,'r')
      call readOrcaEngrad(iorca, energy, gradient)

      ! clean calculation directory !
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

!> wrapper for ORCA single-point run 
subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
      & energy, gradient, sigma, hlgap, results)
   
   !> traceback
   character(len=*), parameter :: source = 'extern_orca_singlepoint'
   
   !> calculator instance
   class(TOrcaCalculator), intent(inout) :: self
   
   !> computational environment
   type(TEnvironment), intent(inout) :: env
   
   !> molecular structure data
   type(TMolecule), intent(inout) :: mol
   
   !> wavefunction data
   type(TRestart), intent(inout) :: chk
   
   !> print level for I/O
   integer, intent(in) :: printlevel
   
   !> restart from previous results
   logical, intent(in) :: restart
   
   !> total energy
   real(wp), intent(out) :: energy
   
   !> molecular gradient
   real(wp), intent(out) :: gradient(:, :)
   
   !> strain derivatives
   real(wp), intent(out) :: sigma(:, :)
   
   !> HOMO-LUMO gap
   real(wp), intent(out) :: hlgap
   
   !> detailed results
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

   ! Intialization !
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   hlgap = 0.0_wp
   efix = 0.0_wp
   dipole(:) = 0.0_wp
   call mol%update 

   ! system call wrapper !
   call runOrca(env,self%ext,mol,energy,gradient)

   ! error handling !
   call env%check(exitRun)
   if (exitRun) then
      call env%error("Electronic structure method terminated", source)
      return
   end if

   ! various external potentials !
   call constrain_pot(potset,mol%n,mol%at,mol%xyz,gradient,efix)
   call constrpot   (mol%n,mol%at,mol%xyz,gradient,efix)
   call cavity_egrad(mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (metaset,mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (rmsdset,mol%n,mol%at,mol%xyz,efix,gradient)
   energy = energy + efix

   ! save results !
   results%e_total = energy
   results%gnorm = norm2(gradient)
   results%dipole = dipole

   ! exact fixing !
   if (fixset%n.gt.0) then
      do i=1, fixset%n
         gradient(1:3,fixset%atoms(i))=0
      enddo
   endif

   ! print calculation results !  
   if (printlevel.ge.2) then

      ! start with summary header !
      if (.not.set%silent) then
         write(env%unit,'(/,9x,53(":"))')
         write(env%unit,'(9x,"::",18x,a,19x,"::")') "ORCA SUMMARY"
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
subroutine hessian(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad)
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
   !> Array to add polarizability gradient to
   real(wp), intent(inout), optional :: polgrad(:, :)

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
      call writeOrcaInp(env,iorca,mol0,self%ext%input_string, "anfreq")
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

!> random name generator
function get_random_name() result(str)

   !> random name
   character(len=:), allocatable :: str
   
   !> random number
   real :: rnd

   !> reinterpretation of real number as an integer 
   integer :: irnd

   call random_number(rnd)
   irnd = transfer(rnd, irnd)
   str = to_hex(abs(irnd))

end function

!> transform integer value into hexadecimal
pure function to_hex(val) result(string)

   !> integer input
   integer, intent(in) :: val
   
   !> resulting name
   character(len=:), allocatable :: string

   !> length of the buffer
   integer, parameter :: buffer_len = range(val) + 2

   !> intermediate value buffer
   character(len=buffer_len) :: buffer
   
   !> position in the buffer
   integer :: pos

   !> current integer value
   integer :: n

   !> hexadecimal base
   integer, parameter :: base = 16
   
   !> hexadecimal digits mapping 
   character(len=1), parameter :: numbers(0:base-1) = &
      ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f"]


   ! initialization !
   n = val
   buffer = ""
   pos = buffer_len + 1 
   
   ! return 0 !
   if (val == 0) then
      string = numbers(0)
      return
   end if
   
   ! decimal -> hexadecimal !
   do while (n > 0)
      pos = pos - 1
      buffer(pos:pos) = numbers(mod(n, base))
      n = n/base
   end do

   string = buffer(pos:)

end function to_hex

end module xtb_extern_orca
