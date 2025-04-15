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
   use xtb_type_calculator, only : TCalculator
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule, len
   use xtb_type_param, only : scc_parameter
   use xtb_type_restart, only : TRestart
   use xtb_io_writer, only : writeMolecule
   use xtb_mctc_systools
   use xtb_setparam
   use xtb_readin
   use xtb_mctc_convert
   use xtb_mctc_systools
   use xtb_setparam
   use xtb_fixparam
   use xtb_scanparam
   use xtb_sphereparam
   use xtb_metadynamic
   use xtb_constrainpot
   implicit none
   private

   public :: TMopacCalculator, newMopacCalculator


   type, extends(TCalculator) :: TMopacCalculator
      type(qm_external) :: ext
   contains
      !> Perform single point calculation
      procedure :: singlepoint
      !> Write informative printout
      procedure :: writeInfo
   end type TMopacCalculator


contains


subroutine checkMopac(env,ext)
   character(len=*), parameter :: source = 'extern_mopac_checkMopac'
   type(TEnvironment), intent(inout) :: env
   type(qm_external), intent(inout) :: ext
   character(len=:),allocatable :: homedir,syspath
   character(len=5) :: chdum
   integer :: i
   logical :: exist

   ! check for MOPAC executable
   if (allocated(ext%executable)) then ! user input
      if (ext%executable(1:1).eq.'~') then
         ! this is relative to the users home, expand it
         call rdvar('HOME',homedir)
         ext%executable = homedir // ext%executable(2:)
         if (set%verbose) then
            write(stdout,'(a,1x,a)') &
               "user home directory        :",homedir
         endif
      endif
      inquire(file=ext%executable,exist=exist)
      if (.not.exist) then
         call env%error("'"//ext%executable//"' was not found, please check",source)
         return
      endif
   else ! no executable provided, lets find it
      call rdvar('PATH',syspath)
      if (set%verbose) then
         write(stdout,'(a,1x,a)') &
            "system path                :",syspath
      endif
      call rdpath(syspath,'mopac',ext%executable,exist)
      if (.not.exist) then
         call env%error('Could not locate mopac executable',source)
         return
      endif
   endif
   if (set%verbose) then
      write(stdout,'(a,1x,a)') &
         "mopac executable           :",ext%executable
   endif

   ! see if there is a preference for an input file
   if (allocated(ext%input_file)) then
      inquire(file=ext%input_file,exist=exist)
      ext%exist = exist
   else
      ext%input_file = 'GCOORD'
      inquire(file=ext%input_file,exist=exist)
   endif
   if (set%verbose) then
      write(stdout,'(a,1x,a)') &
         "mopac input file           :",ext%input_file,&
         "mopac input present        :",bool2string(exist)
         !"mopac input override       :",bool2string(.not.ext%exist)
   endif

   ! check for the input line
   if (allocated(ext%input_string)) then
      if (index(ext%input_string,'grad') == 0) then
         call env%warning('added grad keyword to mopac input', source)
         ext%input_string = ext%input_string //' grad'
      end if
      if (index(ext%input_string,'charge=') == 0) then
         write(chdum,'(i5)') set%ichrg
         ! add total charge
         ext%input_string = ext%input_string //' charge='//trim(adjustl(chdum))
      endif
      if (set%nalphabeta > 0) then
         if (index(ext%input_string,'uhf') == 0) then
            ! write spin state if necessary
            select case(set%nalphabeta)
            case default ! skip
            case(1); ext%input_string = ext%input_string // ' uhf doublet'
            case(2); ext%input_string = ext%input_string // ' uhf triplet'
            case(3); ext%input_string = ext%input_string // ' uhf quartet'
            case(4); ext%input_string = ext%input_string // ' uhf quintet'
            case(5); ext%input_string = ext%input_string // ' uhf sextet'
            case(6); ext%input_string = ext%input_string // ' uhf septet'
            case(7); ext%input_string = ext%input_string // ' uhf octet'
            end select
         endif
      endif
   else
      ! general input
      ext%input_string = '1scf pm6-d3h4 aux(42,PRECISION=12,MOS=-99999,COMP)'
      ! write spin state if necessary
      select case(set%nalphabeta)
      case default ! skip
      case(1); ext%input_string = ext%input_string // ' uhf doublet'
      case(2); ext%input_string = ext%input_string // ' uhf triplet'
      case(3); ext%input_string = ext%input_string // ' uhf quartet'
      case(4); ext%input_string = ext%input_string // ' uhf quintet'
      case(5); ext%input_string = ext%input_string // ' uhf sextet'
      case(6); ext%input_string = ext%input_string // ' uhf septet'
      case(7); ext%input_string = ext%input_string // ' uhf octet'
      end select
      ! convergence criterium in kcal/mol
      ext%input_string = ext%input_string // &
         ' scfcrt=6.d-5 geo-ok mmok grad xyz charge='
      write(chdum,'(i5)') set%ichrg
      ! add total charge
      ext%input_string = ext%input_string // trim(adjustl(chdum))
   endif
   if (set%verbose) then
      write(stdout,'(a,1x,a)') &
         "mopac input line           :",ext%input_string
   endif

end subroutine checkMopac


!> Create a new calculator for driving the Mopac program
subroutine newMopacCalculator(self, env, ext)
   !> Instance of the Mopac calculator
   type(TMopacCalculator), intent(out) :: self
   !> Calculation environment
   type(TEnvironment), intent(inout) :: env
   !> Settings for the Mopac calculator
   type(qm_external), intent(in) :: ext

   self%ext = ext
   self%threadsafe = .false.
   call checkMopac(env, self%ext)
end subroutine newMopacCalculator


subroutine runMopac(env,ext,nat,at,xyz,energy,gradient,dipole)
   character(len=*), parameter :: source = 'extern_mopac_checkMopac'
   type(TEnvironment), intent(inout) :: env
   type(qm_external), intent(in) :: ext
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
   call open_file(imopac,ext%input_file,'w')
   write(imopac,'(a,/,/)') ext%input_string
   do i = 1, nat
      write(imopac,'(3x,a2,3(f20.14,i5))') &
         toSymbol(at(i)),autoaa*xyz(1,i),1,autoaa*xyz(2,i),1,autoaa*xyz(3,i),1
   enddo
   call close_file(imopac)

   write(stdout,'(72("="))')
   write(stdout,'(1x,"*",1x,a)') &
      "handing control over to mopac..."
   call execute_command_line('exec 2>&1 '//ext%executable//' '// &
                             ext%input_file,exitstat=err)
   if (err.ne.0) then
      call env%error('mopac returned with non-zero exit status, following this',source)
   else
      if (set%verbose) then
         inquire(file=ext%input_file//'.arc',exist=exist)
         if (exist) then
            call open_file(imopac,ext%input_file//'.arc','r')
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

   call open_file(imopac,ext%input_file//'.aux','r')
   if (imopac.eq.-1) then
      call env%error("Could not find '"//ext%input_file//".aux'",source)
   else
      read_mopac_output: do
         call getline(imopac, line, iostat=err)
         if (is_iostat_end(err)) then
            if (nat > 1) then ! Gradient is only contained if more than one atom
               call env%error('Could not find gradient in mopac output', source)
            else
               gradient = 0.0_wp
            endif
            exit read_mopac_output  ! Exit the loop on end-of-file
         endif
         if (index(line, 'HEAT_OF_FORMATION:KCAL/MOL') > 0) then
            call readl(line, dum, num)
            energy = dum(num) * kcaltoau
            cycle read_mopac_output
         endif
         if (index(line, 'DIP_VEC:DEBYE') > 0) then
            call readl(line, dum, num)
            dipole(1:3) = dum(2:4) * dtoau
            cycle read_mopac_output
         endif
         if (index(line, 'GRADIENTS:KCAL/MOL/ANGSTROM') > 0) then ! Gradient is only contained if more than one atom
            read(imopac, *) ((gradient(j, i), j = 1, 3), i = 1, nat)
            exit read_mopac_output
         endif
         ! Add a fallback exit condition to avoid infinite loops
         if (err /= 0) then
            call env%error('Unexpected error while reading mopac output', source)
            exit read_mopac_output
         endif
      enddo read_mopac_output
      call close_file(imopac)
      gradient = gradient * kcaltoau / aatoau
   endif
   !$omp end critical (mopac_lock)

end subroutine runMopac


subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
      & energy, gradient, sigma, hlgap, results)
   !> Source of the generated errors
   character(len=*), parameter :: source = 'extern_mopac_singlepoint'
   !> Calculator instance
   class(TMopacCalculator), intent(inout) :: self
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

   call runMopac(env,self%ext,mol%n,mol%at,mol%xyz,energy,gradient,dipole)

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

subroutine writeInfo(self, unit, mol)

   !> Calculator instance
   class(TMopacCalculator), intent(in) :: self

   !> Unit for I/O
   integer, intent(in) :: unit

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   call generic_header(unit,"Mopac driver",49,10)
end subroutine writeInfo


end module xtb_extern_mopac
