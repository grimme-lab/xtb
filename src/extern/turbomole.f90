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

module xtb_extern_turbomole
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
   use xtb_io_writer, only : writeMolecule
   use xtb_io_reader, only : readHessian
   use xtb_type_reader, only : TReader
   use xtb_mctc_filetypes, only : fileType
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

   public :: wrtm, rdtm
   public :: TTMCalculator, newTMCalculator

   type, extends(TCalculator) :: TTMCalculator
      type(qm_external) :: ext
      integer :: extcode
      integer :: extmode
   contains
      !> Perform single point calculation
      procedure :: singlepoint
      !> Calculate hessian
      procedure :: hessian
      !> Write informative printout
      procedure :: writeInfo
   end type TTMCalculator

contains

!----------------------------------------------------------
! Create a new calculator for driving the Turbomole program
!----------------------------------------------------------
subroutine newTMCalculator(self,extcode,extmode)
   
   implicit none
   type(TTMCalculator), intent(out) :: self
      !! Instance of the Turbomole calculator
   integer, intent(in) :: extcode
      !! RI, RI+d3+gCP,NoRI
   integer, intent(in) :: extmode
      !! ridft+rdgrad, ridft+rdgrad+dftd3+gcp

   self%threadsafe = .false.
      !! No parallelization
   self%extcode = extcode
   self%extmode = extmode

end subroutine newTMCalculator


subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
      & energy, gradient, sigma, hlgap, results)
   !> Source of the generated errors
   character(len=*), parameter :: source = 'extern_turbomole_singlepoint'
   !> Calculator instance
   class(TTMCalculator), intent(inout) :: self
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

   call external_turbomole(env,mol%n,mol%at,mol%xyz,chk%wfn%nel,chk%wfn%nopen, &
      & self%extcode,self%extmode,.true.,energy,gradient,results%dipole,self%lSolv,mol%chrg,mol%uhf)

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

!---------------------------------------------------------
! To run TURBOMOLE, and cefine if needed 
!---------------------------------------------------------
subroutine external_turbomole(env,n,at,xyz,nel,nopen,extcode,extmode,grd,eel,g,dip,lsolv,chrg,uhf)
   
   use xtb_mctc_accuracy, only : wp
   use xtb_setparam
   
   implicit none
   character(len=*),parameter :: source = 'external_turbomole'
   type(TEnvironment), intent(inout) :: env
      !! Calculation environment to handle I/O stream and error log
   integer :: n, at(n), nel, nopen
      !! structural information from mol(Tmolecule) and chk(TRestart) 
   logical ::  grd
      !! if gradient calcultion is needed
   logical :: lsolv
      !! if solvated
   integer, intent(in) :: extcode, extmode
      !! turbomole settings
   real(wp) :: xyz(3,n)
      !! coordinates
   real(wp) :: xyz_cached(3,n)
      !! to save xyz values
   real(wp) :: g (3,n)
      !! gradients
   real(wp) :: eel
      !! electronic energy
   real(wp) :: dip(3)
      !! dipole moments
   real(wp), intent(in) :: chrg
      !! charge
   integer, intent(in) :: uhf
      !! multiplicity, number of unpaired electrons

   character(len=255) atmp
   character(len=:), allocatable :: syspath, cefine
   logical :: cache
      !! to check if molecule moves b/n optimization runs
   logical :: exist
   integer :: err

   cache = .false.
   dip=0
   
   !> TM input
   inquire(file="control", exist=exist)
      !! if control file present in the calc dir
   if (.not.exist) then
      
      !> cefine 
      if(len(set%ext_turbo%input_string).ne.0) then      
         call rdvar("PATH", syspath)
         call rdpath(syspath, "cefine", cefine, exist)
         if (exist) then
            call wrtm(n,at,xyz)
            if (allocated(set%ext_turbo%input_string)) then 
               call execute_command_line("exec "//cefine//" "//set%ext_turbo%input_string)
            else
               call execute_command_line("exec "//cefine//" --func b97-3c")
            endif
         else
            call env%error("No cefine binary is found", source)
            return
         end if

      else
         call writeControl(chrg,uhf)
            !! default control file
      endif

   end if

   inquire(file="control", exist=exist)
   if (.not.exist) then
      call env%error("No 'control' file in current directory")
      return
   end if
   
  
   ! TM (RI)
   if(extcode.eq.1)then
      !$omp critical (turbo_lock)
      inquire(file='gradient', exist=exist)
      if (exist .and. grd) then
         call rdtm(env,n,grd,eel,g,xyz_cached)
         cache = all(abs(xyz_cached - xyz) < 1.e-7_wp)
      end if
      if (.not.cache) then
         call wrtm(n,at,xyz)
         if(extmode.eq.1)then
            
            write(env%unit,'(72("="))')
            write(env%unit,'(1x,"*",1x,a)') &
               "letting TURBOMOLE take over the control..."
            
            if (set%oniom_settings%silent) then
               call execute_command_line('exec ridft >  job.last 2>> /dev/null',exitstat=err)
               if(grd)call execute_command_line('exec rdgrad >> job.last 2>> /dev/null ',exitstat=err)
            else
               call execute_command_line('exec ridft  | tee  job.last 2>> /dev/null',exitstat=err)
               if(grd)call execute_command_line('exec rdgrad |tee -a job.last 2>> /dev/null ',exitstat=err)
            endif

            if (err.ne.0) then
               call env%error('TURBOMOLE returned with non-zero exit status, doing the same',source)
            else
               write(env%unit,'(1x,"*",1x,a)') &
               "successful TURBOMOLE run, taking over control again..."
            endif
               write(env%unit,'(72("="))')
       endif
         call extcodeok(extcode)
         call rdtm(env,n,grd,eel,g,xyz_cached)
         if (set%ceasefiles) then
            call env%io%deleteFile('job.last') 
            call env%io%deleteFile('mos') 
            call env%io%deleteFile('statistics') 
            call env%io%deleteFile('basis') 
            call env%io%deleteFile('auxbasis') 
            call env%io%deleteFile('coord') 
         endif
      end if
      !$omp end critical (turbo_lock)
      return
   endif

   ! TM+d3+gcp
   if(extcode.eq.2)then
      !$omp critical (turbo_lock)
      inquire(file='gradient', exist=exist)
      if (exist .and. grd) then
         call rdtm(env,n,grd,eel,g,xyz_cached)
         cache = all(abs(xyz_cached - xyz) < 1.e-7_wp)
      end if
      if (.not.cache) then
         call wrtm(n,at,xyz)
         if(extmode.le.2)then
            call execute_command_line('exec ridft  >  job.last 2>> /dev/null')
            call execute_command_line('exec rdgrad >> job.last 2>> /dev/null')
            call execute_command_line('exec dftd3 coord -grad >> job.last 2>> /dev/null')
            call execute_command_line('exec gcp coord -file -grad >>job.last 2>>/dev/null')
         endif
         call extcodeok(extcode)
         call rdtm(env,n,.true.,eel,g,xyz_cached)
         
         if (set%ceasefiles) then
            call env%io%deleteFile('job.last') 
            call env%io%deleteFile('mos') 
            call env%io%deleteFile('statistics') 
            call env%io%deleteFile('basis') 
            call env%io%deleteFile('auxbasis') 
            call env%io%deleteFile('coord') 
         endif
      end if
      !$omp end critical (turbo_lock)
      return
   endif

   ! TM (NORI)
   if(extcode.eq.3)then
      !$omp critical (turbo_lock)
      inquire(file='gradient', exist=exist)
      if (exist .and. grd) then
         call rdtm(env,n,grd,eel,g,xyz_cached)
         cache = all(abs(xyz_cached - xyz) < 1.e-7_wp)
      end if
      if (.not.cache) then
         call wrtm(n,at,xyz)
         if(extmode.eq.1)then
            call execute_command_line('exec dscf  > job.last 2>> /dev/null')
            if(grd)call execute_command_line('exec grad >> job.last 2>> /dev/null')
         endif
         call extcodeok(extcode)
         call rdtm(env,n,grd,eel,g,xyz_cached)
         if (set%ceasefiles) then
            call env%io%deleteFile('job.last') 
            call env%io%deleteFile('mos') 
            call env%io%deleteFile('statistics') 
            call env%io%deleteFile('basis') 
            call env%io%deleteFile('auxbasis') 
            call env%io%deleteFile('coord') 
         endif
      end if
      !$omp end critical (turbo_lock)
      return
   endif


   call raise('E','This external code is not implemented')

end subroutine external_turbomole

!-------------------------------------------------------
! create default TM control -> future TMprep
!-------------------------------------------------------
subroutine writeControl(chrg,uhf)

   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_symbols, only : toSymbol
   implicit none
   real(wp), intent(in) :: chrg
      !! charge
   integer, intent(in) :: uhf
      !! number of unpaired electrons
   integer :: iunit
      !! number of electrins

   iunit=39
   open(newunit=iunit,file='control')

   write(iunit,'(a)')'$coord file=coord'
   write(iunit,'(a,i0,a,i0)')'$eht charge=',nint(chrg),' unpaired=',uhf
   write(iunit,'(a)')'$symmetry c1'
   write(iunit,'(a)')'$atoms'
   write(iunit,'(3x,a)')'basis =def2-mTZVP'
   write(iunit,'(a)')'$dft'
   write(iunit,'(3x,a)')'functional b97-3c'
   write(iunit,'(3x,a)')'gridsize m4'
   write(iunit,'(a)')'$rij'
   write(iunit,'(a)')'$energy file=energy'
   write(iunit,'(a)')'$grad file=gradient'
   write(iunit,'(a)')'$disp3 -bj -abc'
   write(iunit,'(a)')'$end'

   close(iunit)

end subroutine writeControl

subroutine wrtm(n,at,xyz)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_symbols, only : toSymbol
   implicit none
   integer n, at(n), iunit, i
   real(wp) xyz(3,n)

   iunit=33
   open(newunit=iunit,file='coord')

   write(iunit,'(a)')'$coord'
   do i=1,n
      write(iunit,'(3F24.14,6x,a2)') &
      &   xyz(1,i),xyz(2,i),xyz(3,i),toSymbol(at(i))
   enddo
   write(iunit,'(a)')'$end'

   close(iunit)

end subroutine wrtm

!----------------------------------------
! read TM energy and gradient files 
!----------------------------------------
!subroutine readTM(n,ifgrd,energy,gradient,xyz)
   
 !  use xtb_mctc_accuracy, only : wp
   !use xtb_filetools, only : open_file, close_file, remove_file
 !  use xtb_readin , only : strip_line,getValue

!   implicit none
!   integer, intent(in) :: n
      !! number of atoms
!   logical, intent(in) :: ifgrd
      !! if gradient fileshould be read
!   real(wp), intent(out) :: energy, gradient(3,n), xyz(3,n)
      !! TM output
   
!   integer :: grd, e
      !! file units for 'gradient' and 'energy' files
!   character(len=:), allocatable :: line
      !! tmp line
!   logical :: exist


!   if (.not.ifgrd) then
!      call open_file(e,'energy','r')
!      do
!         call strip_line(e,line,exist)
!         if (exist)
!            call readline()
!      enddo
!   endif
!
!
!end subroutine readTM
subroutine rdtm(env,n,grd,e,g,xyz)
   
   use xtb_mctc_accuracy, only : wp
   implicit none
   type(TEnvironment), intent(inout) :: env
   integer n, iunit, i, nl, j, nn
   logical grd
   real(wp) g(3,n), e, xx(10), x, y, z, xyz(3,n)
   logical ex
   character(len=128) a1

   iunit=33

   if(.not.grd)then
      open(newunit=iunit,file='energy')
      101   read(iunit,'(a)',end=102)a1
      call readl(a1,xx,nn)
      if(nn.ge.4) e=xx(2)
         !! to assign the last entry 
      goto 101
      102 continue
      if (set%ceasefiles) then
         call env%io%deleteFile('energy') 
      else
         close(iunit)
      endif
      return
   endif

   inquire(file='gradient',exist=ex)
   if(.not.ex) then
      call raise('E','no gradient file found!')
   endif

   j=0
   open(newunit=iunit,file='gradient')
201 read(iunit,'(a)',end=301)a1
   j=j+1
   if(index(a1,'cycle').ne.0)nl=j
   goto 201
   301   continue

   if(nl.lt.2)then
      call raise('E','illegal gradient file!')
   endif

   rewind iunit
   do i=1,nl
      read(iunit,'(a)')a1
   enddo
   call readl(a1,xx,nn)
   e=xx(2)
   do i=1,n
      read(iunit,*)xyz(1,i),xyz(2,i),xyz(3,i)
   enddo
   do i=1,n
      read(iunit,*)g(1,i),g(2,i),g(3,i)
   enddo
   if (set%ceasefiles) then
      call env%io%deleteFile('energy') 
      call env%io%deleteFile('gradient') 
   else
      close(iunit)
   endif

end subroutine rdtm

subroutine extcodeok(extcode)
   implicit none
   integer, intent (in) :: extcode
   integer :: ich
   character(len=80) atmp
   ! TM
   if(extcode.le.3)then
      call execute_command_line('exec grep "actual step" control > TmPfIlE')
      open(newunit=ich,file='TmPfIlE',status='old')
      read(ich,'(a)',end=100)atmp
 100  close(ich,status='delete')
      if(index(atmp,'actual').ne.0) call raise('E','external code error: '//&
      &                                            trim(atmp))
   endif

   return
end subroutine extcodeok

subroutine writeInfo(self, unit, mol)

   !> Calculator instance
   class(TTMCalculator), intent(in) :: self

   !> Unit for I/O
   integer, intent(in) :: unit

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   call generic_header(unit,"Orca driver",49,10)
end subroutine writeInfo

!> Evaluate hessian by finite difference for all atoms
subroutine hessian(self, env, mol0, chk0, list, step, hess, dipgrad)
   character(len=*), parameter :: source = "extern_turbomole_hessian"
   !> Single point calculator
   class(TTMCalculator), intent(inout) :: self
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

   integer :: idipd, stat
   type(TReader) :: reader

   call wrtm(mol0%n,mol0%at,mol0%xyz) !Overwrite coord with RAM-xyz file

   call execute_command_line('exec aoforce > job.last2>> /dev/null')
   call reader%open('hessian')
   call readHessian(env, mol0, hess, reader, fileType%tmol)
   call reader%close

   call open_file(idipd,'dipgrad','r')
   if(idipd == -1) then
      call env%error("No dipolegradient found", source)
      return
   end if

   call read_dipgrad(idipd, mol0%n, dipgrad, stat)

   if(stat /=0 ) then
      call env%error('An error occurred while reading the dipolegradient', source)
      return
   end if

   call close_file(idipd)
end subroutine hessian

subroutine read_dipgrad(idipd, n, dipd, error)
   use xtb_mctc_systools

   implicit none

   integer, intent(out) :: error
   character(len=:), allocatable :: line
   integer, intent(in) :: n
   real(wp), intent(inout) :: dipd(:,:)
   integer, intent(in) :: idipd
   integer :: i,j

   error = 0

   do while(error == 0)
      call getline(idipd, line, error)
      if (index(line, '$dipgrad          cartesian dipole gradients') == 1) then
         do i=1, 3*n
            read(idipd,*) (dipd(j,i),j=1,3)
         end do
         exit
      end if
   enddo
end subroutine read_dipgrad

end module xtb_extern_turbomole
