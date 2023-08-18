! This file is part of xtb.
! SPDX-Identifier: LGPL-3.0-or-later
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

module xtb_extern_driver
  use xtb_mctc_accuracy, only: wp
  use xtb_mctc_io, only: stdout
  use xtb_mctc_filetypes, only: fileType, generateFileName
  use xtb_mctc_symbols, only: toSymbol
  use xtb_type_calculator, only: TCalculator
  use xtb_type_data, only: scc_results
  use xtb_type_environment, only: TEnvironment
  use xtb_type_molecule, only: TMolecule, len
  use xtb_type_param, only: scc_parameter
  use xtb_type_restart, only: TRestart
  use xtb_io_writer, only: writeMolecule
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
  use xtb_extern_turbomole, only: wrtm, rdtm
  implicit none
  private

  public :: TDriverCalculator, newDriverCalculator

  type, extends(TCalculator) :: TDriverCalculator
    type(qm_external) :: ext
  contains
    !> Perform single point calculation
    procedure :: singlepoint
    !> Write informative printout
    procedure :: writeInfo
  end type TDriverCalculator

contains

!> Create a new calculator for driving the external driver program
  subroutine newDriverCalculator(self, env, ext)
    !> Instance of the external driver calculator
    type(TDriverCalculator), intent(out) :: self
    !> Calculation environment
    type(TEnvironment), intent(inout) :: env
    !> Settings for the external driver calculator
    type(qm_external), intent(in) :: ext

    self%threadsafe = .false.
    self%ext = ext
  end subroutine newDriverCalculator

  subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
        & energy, gradient, sigma, hlgap, results)
    !> Source of the generated errors
    character(len=*), parameter :: source = 'extern_driver_singlepoint'
    !> Calculator instance
    class(TDriverCalculator), intent(inout) :: self
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

    integer :: i, ich
    integer :: mode_sp_run = 1
    real(wp) :: efix
    real(wp) :: dipole(3)
    logical :: cache, exist
    logical, parameter :: ccm = .true.
    logical :: exitRun
    character(len=*), parameter :: outfmt = &
                                   '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'
    real(wp) :: xyz_cached(3, mol%n)
    integer :: err
    character(len=:),allocatable :: extension
    character(len=:),allocatable :: tmpname

    call mol%update

    cache = .false.
    energy = 0.0_wp
    gradient(:, :) = 0.0_wp
    sigma(:, :) = 0.0_wp
    hlgap = 0.0_wp
    efix = 0.0_wp
    dipole(:) = 0.0_wp

    !$omp critical (turbo_lock)
    inquire (file='gradient', exist=exist)
    if (exist) then
      ! ### only TM output is supported for now ###
      call rdtm(env,mol%n, .true., energy, gradient, xyz_cached)
      cache = all(abs(xyz_cached - mol%xyz) < 1.e-10_wp)
      if (cache .and. printlevel > 0) then 
         write(env%unit,'(/,a,/)') &
         & "Geometry is equivalent to the one in &
         & 'gradient'. Reading gradient from: 'gradient'."
      endif
    end if
    if (.not. cache) then
      call generateFileName(tmpname, 'xtbdriver', extension, mol%ftype)
      if (printlevel > 0) then
        write(env%unit,'(/,a,1x,a,/)') &
           "updated geometry written to:",tmpname
      endif
      call open_file(ich,tmpname,'w')
      if (exist) then
        call writeMolecule(mol, ich, format=mol%ftype, energy=energy, &
              & gnorm=norm2(gradient))
      else
        call writeMolecule(mol, ich, format=mol%ftype)
      end if
      call close_file(ich)

      write (env%unit, '(72("="))')
      write (env%unit, '(1x,"*",1x,a)') &
        "letting driver take over the control..."
      call execute_command_line('exec 2>&1 '//self%ext%executable, exitstat=err)
      if (err /= 0) then
        call env%error('driver returned with non-zero exit status, doing the same', source)
      else
        write (env%unit, '(1x,"*",1x,a)') &
          "successful driver run, taking over control again..."
      end if
      write (env%unit, '(72("="))')

      ! ### only TM output is supported for now ###
      call rdtm(env,mol%n, .true., energy, gradient, xyz_cached)
    end if
    !$omp end critical (turbo_lock)

    call env%check(exitRun)
    if (exitRun) then
      call env%error("Electronic structure method terminated", source)
      return
    end if

    ! ------------------------------------------------------------------------
    !  various external potentials
    call constrain_pot(potset, mol%n, mol%at, mol%xyz, gradient, efix)
    call constrpot(mol%n, mol%at, mol%xyz, gradient, efix)
    call cavity_egrad(mol%n, mol%at, mol%xyz, efix, gradient)
    call metadynamic(metaset, mol%n, mol%at, mol%xyz, efix, gradient)
    call metadynamic(rmsdset, mol%n, mol%at, mol%xyz, efix, gradient)

    ! ------------------------------------------------------------------------
    !  fixing of certain atoms
    !  print*,abs(efix/etot)
    energy = energy + efix
    results%e_total = energy
    results%gnorm = norm2(gradient)
    results%dipole = dipole
    if (fixset%n .gt. 0) then
      do i = 1, fixset%n
        !print*,i,fixset%atoms(i)
        gradient(1:3, fixset%atoms(i)) = 0
      end do
    end if

    if (printlevel .ge. 2) then
      ! start with summary header
      if (.not. set%silent) then
        write (env%unit, '(9x,53(":"))')
        write (env%unit, '(9x,"::",21x,a,21x,"::")') "SUMMARY"
      end if
      write (env%unit, '(9x,53(":"))')
      write (env%unit, outfmt) "total energy      ", results%e_total, "Eh   "
      write (env%unit, outfmt) "gradient norm     ", results%gnorm, "Eh/a0"
      write (env%unit, outfmt) "HOMO-LUMO gap     ", results%hl_gap, "eV   "
      write (env%unit, '(9x,53(":"))')
      write (env%unit, '(a)')
    end if
  end subroutine singlepoint

  subroutine writeInfo(self, unit, mol)

    !> Calculator instance
    class(TDriverCalculator), intent(in) :: self

    !> Unit for I/O
    integer, intent(in) :: unit

    !> Molecular structure data
    type(TMolecule), intent(in) :: mol

    call generic_header(unit, "Driver", 49, 10)
  end subroutine writeInfo

end module xtb_extern_driver
