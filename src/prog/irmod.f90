! This file is part of xtb.
!
! Copyright (C) 2022 Julia Kohn, Sebastian Ehlert
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

!> IR / Raman intensities from DFTB+ Hessians
module xtb_prog_ir
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autorcm, autoamu
   use xtb_mctc_filetypes, only : getFileType, hessType
   use xtb_mctc_timings
   use xtb_mctc_version, only : version, author, date
   use xtb_io_reader, only : readMolecule, readHessian
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_reader, only : TReader
   use xtb_prog_argparser, only : TArgParser
   use xtb_freq_project, only : projectHessian
   use xtb_freq_utils, only : massWeightHessian, diagHessian
   use xtb_propertyoutput, only : print_thermo
   use xtb_splitparam, only : atmass
   use xtb_setmod, only : set_thermo, set_symmetry
   implicit none
   private

   public :: xtbIR


contains


!> Main program of the IR submodule
subroutine xtbIR(env, argParser)

    use xtb_freq_io, only: write_tm_vibspectrum, wrhess, writeHessianOut, g98fake2

   !> Source of errors in the main program unit
   character(len=*), parameter :: source = "prog_ir"

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Commandline argument parser
   type(TArgParser), intent(inout) :: argParser

   !> Molecular structure data
   type(TMolecule) :: mol

   type(TReader) :: reader
   character(len=:), allocatable :: file,born
   real(wp), allocatable :: freq(:), hessian(:, :), dipd(:, :), dipt(:),atmasslong(:),rmass(:)
   real(wp) :: htot, gtot, zp, sum2, trdip(3), autokmmol=1.7770969e+6_wp,xsum
   integer :: nFiles, ftype, hFormat, nimag, i,j,k, ifile, ia, ic,ios
   logical :: massWeighted

   !> dipd are 1. mixed 2nd derivative of energy wrt atomic positions and electric field
   !           2. 1st derivative of dipole moments wrt to atomic positions
   !           3. 1st derivative of forces wrt electric field
   ! 1. == 2.== 3. == Born charges

   call parseArguments(env, argParser, hFormat, massWeighted, born)

   nFiles = argParser%countFiles()
   if (nFiles == 0) then
      call env%error("No input file given, so there is nothing to do", source)
      call IRHelp(env%unit)
   end if
   if (nFiles > 2) then
      call env%error("Multiple input files present, aborting", source)
   end if

   call env%checkpoint("Command line argument parsing failed")

   call argParser%nextFile(file)

   !> Determine the file type for processing in the reader
   ftype = getFileType(file)

   !> Generate the reader by connecting the file to the instance
   call reader%open(file)
   !> Process the file to obtain the molecular structure data
   call readMolecule(env, mol, reader%unit, ftype)
   !> Close the file, we are done with it here
   call reader%close

   call env%checkpoint("Could not read geometry from '"//file//"'")

   !> Print an informative banner
   call IRHeader(env%unit)
   !> print current time
   call prdate('S')

   !> Also store a copy of the atomic masses in amu in global storage
   atmass = mol%atmass * autoamu

   !>read in dipole derivatives from dftb+ custom file (born.out)
   allocate(dipd(3,3*mol%n))

   call reader%open(born)
   do i=1,3*mol%n
      read(reader%unit,*,iostat=ios) dipd(1:3,i)  
      if (ios < 0 ) then
          !>abort run if born charges can't be read
          write(env%unit,*) "'"//born//"' file corrupted ... aborting" 
          stop 
      end if
   end do
   call reader%close

   write(env%unit,*) "  "
   write(env%unit,*) "Dipole derivatives (DFTB+) read in from file '"//born//"' "
   write(env%unit,*) "  "

   !>read in hessian from dftb+ custom file (hessian.out)
   call argParser%nextFile(file)
   allocate(freq(3*mol%n),rmass(3*mol%n))

   if (allocated(file)) then
      allocate(hessian(3*mol%n, 3*mol%n))
      call reader%open(file)
      call readHessian(env, mol, hessian, reader, format=hFormat)
      call reader%close
      
      write(env%unit,*) "Hessian (DFTB+) read in from file '"//file//"' "
      write(env%unit,*) "  "

      if (.not.massWeighted) then
         call massWeightHessian(hessian, mol%atmass)
      end if
      if (hFormat == hessType%dftbplus) then
         call projectHessian(hessian, mol, .true., .true.)
      end if

      call diagHessian(env, hessian, freq)

      call preigf(env%unit, freq, autorcm)

      call env%checkpoint("Could not read hessian from '"//file//"'")
   else
      freq(:) = -1.0_wp
   end if
   freq(:) = freq * autorcm

   !>expansion from atmass(n) to atmasslong(3*n)
   allocate(atmasslong(3*mol%n))
   do j=1,mol%n
      do k=1,3
         atmasslong((j-1)*3+k)=mol%atmass(j)
      end do
   end do

   !>IR intensities copied from src/hessian.F90
   allocate(dipt(3*mol%n))
    do i = 1, 3*mol%n
      do k = 1, 3
         sum2 = 0.0_wp
         do j = 1, 3*mol%n
            sum2 = sum2 + dipd(k,j)*(hessian(j,i)/sqrt(atmasslong(j)))              
         end do
         trdip(k) = sum2
      end do
      dipt(i) = autokmmol*(trdip(1)**2+trdip(2)**2+trdip(3)**2)
   end do

   !>write IR vib_spectrum
   call open_file(ifile,'vib_spectrum','w')
   call write_tm_vibspectrum(ifile,3*mol%n,freq,dipt)
   call close_file(ifile)

   k=0
   do i=1,3*mol%n
      xsum=0
      k=k+1
      do ia=1,mol%n
         do ic=1,3
            j = (ia-1)*3+ic
            xsum=xsum+hessian(j,i)**2/atmasslong(ia)
         enddo
      enddo
      rmass(i)=xsum
   enddo

   !>write fake molden input for mode visualization
   call g98fake2('g98.out',mol%n,mol%at,mol%xyz,freq,rmass,dipt,hessian)

   call print_thermo(env%unit,mol%n, size(freq), mol%at, mol%xyz, freq, 0.0_wp, &
      & htot, gtot, nimag, .true., zp)

   write(env%unit,*) "  "
   write(env%unit,*) "IR spectrum can be found in file vib_spectrum "
   write(env%unit,*) "  "
   write(env%unit,*) "molden input for visualization of spectrum and vibrations can be found in g98.out"

   write(env%unit,'(a)')
   write(env%unit,'(72("-"))')
   call stop_timing_run
   call stop_timing(1)
   call prdate('E')
   write(env%unit,'(72("-"))')
   call prtiming(1,'total')

   write(env%unit,'(a)')
   call terminate(0)

end subroutine xtbIR


subroutine preigf(unit, freq, scale)
   real(wp), intent(in) :: freq(:)
   real(wp), intent(in) :: scale
   integer, intent(in) :: unit
   integer, parameter :: nRows = 6
   integer :: ntimes, nrest, i, j, n2, k
   ntimes = size(freq)/nRows
   nrest = mod(size(freq), nRows)
   if (ntimes.eq.0) nrest = size(freq)
   j = 1
   n2 = nRows
   do k = 1, ntimes
      write(unit,'("eigval :",3x,10f9.2)') (freq(i)*scale,i=j,n2)
      j = j + nRows
      n2 = n2 + nRows
   end do
   if (nrest.gt.0.or.ntimes.eq.0) then
      write(unit,'("eigval :",3x,10f9.2)') (freq(i)*scale,i=j,j+nrest-1)
   end if

end subroutine preigf


!> Parse command line arguments
subroutine parseArguments(env, args, htype, massWeighted, born)

   !> Name of error producer
   character(len=*), parameter :: source = "prog_ir_parseArguments"

   !> Calculation environment
   type(TEnvironment) :: env

   !> Command line argument parser
   type(TArgParser) :: args

   !> File type of the hessian
   integer, intent(out) :: htype

   !> Hessian is already mass weighted
   logical, intent(out) :: massWeighted

   integer :: nFlags
   character(len=:), allocatable :: flag, sec
   character(len=:), allocatable, intent(out) :: born

   htype = hessType%tmol
   massWeighted = .false.

   nFlags = args%countFlags()
   call args%nextFlag(flag)
   do while(allocated(flag))
      if (len(flag) > 2 .and. flag(1:1) == '-' .and. flag(1:2) /= '--') then
         call env%warning("the use of '"//flag//"' is discouraged, "// &
            & "please use '-"//flag//"' next time", source)
         flag = '-'//flag
      end if
      select case(flag)
      case default
         call env%warning("Unknown option '"//flag//"' provided", source)

      case('--help', '-h')
         call IRHelp(env%unit)
         call terminate(0)

      case('--dftb+','--dftbplus') !> '+' is not accepted by every command-line parser
         htype = hessType%dftbplus

      case('--born')   
         call args%nextArg(sec)
         if (allocated(sec)) then
            born=sec
         else
            call env%error("dipole derivatives are missing", source)
         end if

      case('--scale')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_thermo(env, 'scale', sec)
         else
            call env%error("Freq. scaling factor missing", source)
         end if

      case('--sthr')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_thermo(env, 'sthr', sec)
         else
            call env%error("Rotor cutoff is missing", source)
         end if
     
      case('--ithr')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_thermo(env, 'imagthr','-'//sec)
         else
            call env%error("Imaginary cutoff is missing", source)
         end if
      
      case('--desy')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_symmetry(env,'desy',sec)
         else
            call env%error("Threshold for symmetry recognition missing", source)
         end if

      case('--temp')
         call args%nextArg(sec)
         if (allocated(sec)) then
            call set_thermo(env, 'temp', sec)
         else
            call env%error("Temperature in K is missing", source)
         end if

      end select
      call args%nextFlag(flag)
   end do

end subroutine parseArguments


!> Header for this submodule
subroutine IRHeader(unit)

   !> IO unit
   integer, intent(in) :: unit

   write(unit,'(a)') &
      "      -----------------------------------------------------------",&
      "     |                   =====================                   |",&
      "     |                            I R                            |",&
      "     |                   =====================                   |",&
      "     |                    S. Ehlert , J. Kohn                    |",&
      "     |          Mulliken Center for Theoretical Chemistry        |",&
      "     |                    University of Bonn                     |",&
      "      -----------------------------------------------------------",""
   write(unit,'(3x,"*",*(1x,a))') &
      & "xtb-IR version", version, "compiled by", author, "on", date
   write(unit,'(a)')

end subroutine IRHeader


subroutine IRHelp(unit)

   !> IO unit
   integer, intent(in) :: unit

   write(unit, '(a)') &
   "Usage: xtb ir [options] <geometry> --dftbplus <hessian> --born <born>", &
   "",&
   "<geometry> may be provided as any valid input to xtb", &
   "the [hessian] file is read and processed depending on the selected options,", &
   "the default format is a non-massweighted Turbomole hessian file", &
   "",&
   "Options",&
   "",&
   "   --sthr REAL         Rotor cutoff for RRHO partition function in rcm",&
   "",&
   "   --ithr REAL         Imag. freq. cutoff for RRHO in rcm",&
   "                       enter a positive value, the sign will be inverted",&
   "",&
   "   --scale REAL        Frequency scaling factor for RRHO",&
   "",&
   "   --desy REAL         Threshold for symmetry recognition",&
   "",&
   "   --temp REAL[,...]   Temperature for thermodynamic functions in K,",&
   "                       takes a comma separated list of temperatures",&
   "",&
   "   --dftb+             Read a DFTB+ hessian file, implies projection",&
   "   --dftbplus ",&
   "",&
   "   --born              Read a DFTB+ dipole derivatives / born charges file",&
   "",&
   "Output Conventions:",&
   "",&
   "total energies are given in atomic units (Eh), frequencies in reciprocal cm",&
   "",&
   "More information can be obtained at https://xtb-docs.readthedocs.io/",&
   ""
end subroutine IRHelp


end module xtb_prog_ir
