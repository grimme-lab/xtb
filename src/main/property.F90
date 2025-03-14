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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

module xtb_propertyoutput
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_io, only: stdout
   use xtb_mctc_symbols, only: toSymbol
   use xtb_mctc_convert, only: evtoau, autod, autoaa
   use xtb_solv_cm5
   use xtb_cube
   use xtb_topology
#if WITH_TBLITE
   use tblite_basis_type, only: basis_type
   use tblite_wavefunction_type, only: wavefunction_type
#endif
   contains

   subroutine write_energy(iunit, sccres, frqres, hess)
      use xtb_type_data
      implicit none
      integer, intent(in) :: iunit ! file handle (usually output_unit=6)
      logical, intent(in) :: hess
      type(scc_results), intent(in) :: sccres
      type(freq_results), intent(in) :: frqres
      character(len=*), parameter :: outfmt = '(10x,"|",1x,a,f24.12,1x,a,1x,"|")'
      write (iunit, '(a)')
      write (iunit, '(11x,49("-"))')
      if (hess) then
         write (iunit, outfmt) "TOTAL ENERGY      ", frqres%etot, "Eh  "
         write (iunit, outfmt) "TOTAL ENTHALPY    ", frqres%etot + frqres%htot, "Eh  "
         write (iunit, outfmt) "TOTAL FREE ENERGY ", frqres%etot + frqres%gtot, "Eh  "
         write (iunit, outfmt) "GRADIENT NORM     ", frqres%gnorm, "Eh/α"
      else
         write (iunit, outfmt) "TOTAL ENERGY      ", sccres%e_total, "Eh  "
         write (iunit, outfmt) "GRADIENT NORM     ", sccres%gnorm, "Eh/α"
      end if
      write (iunit, outfmt) "HOMO-LUMO GAP     ", sccres%hl_gap, "eV  "
      write (iunit, '(11x,49("-"))')
   end subroutine write_energy

   subroutine write_energy_gff(iunit, sccres, frqres, hess)
      use xtb_type_data
      implicit none
      integer, intent(in) :: iunit ! file handle (usually output_unit=6)
      logical, intent(in) :: hess
      type(scc_results), intent(in) :: sccres
      type(freq_results), intent(in) :: frqres
      character(len=*), parameter :: outfmt = '(10x,"|",1x,a,f24.12,1x,a,1x,"|")'
      write (iunit, '(a)')
      write (iunit, '(11x,49("-"))')
      if (hess) then
         write (iunit, outfmt) "TOTAL ENERGY      ", frqres%etot, "Eh  "
         write (iunit, outfmt) "TOTAL ENTHALPY    ", frqres%etot + frqres%htot, "Eh  "
         write (iunit, outfmt) "TOTAL FREE ENERGY ", frqres%etot + frqres%gtot, "Eh  "
         write (iunit, outfmt) "GRADIENT NORM     ", frqres%gnorm, "Eh/α"
      else
         write (iunit, outfmt) "TOTAL ENERGY      ", sccres%e_total, "Eh  "
         write (iunit, outfmt) "GRADIENT NORM     ", sccres%gnorm, "Eh/α"
      end if
      write (iunit, '(11x,49("-"))')
   end subroutine write_energy_gff

   subroutine write_energy_oniom(iunit, sccres, frqres, hess)
      use xtb_type_data
      implicit none
      integer, intent(in) :: iunit ! file handle (usually output_unit=6)
      logical, intent(in) :: hess
      type(scc_results), intent(in) :: sccres
      type(freq_results), intent(in) :: frqres
      character(len=*), parameter :: outfmt = '(10x,"|",1x,a,f18.12,1x,a,1x,"|")'

      write (iunit, '(a)')
      write (iunit, '(11x,49("-"))')
      if (hess) then
         write (iunit, outfmt) "ONIOM TOTAL ENERGY      ", frqres%etot, "Eh  "
         write (iunit, outfmt) "ONIOM TOTAL ENTHALPY    ", frqres%etot + frqres%htot, "Eh  "
         write (iunit, outfmt) "ONIOM TOTAL FREE ENERGY ", frqres%etot + frqres%gtot, "Eh  "
         write (iunit, outfmt) "ONIOM GRADIENT NORM     ", frqres%gnorm, "Eh/α"
      else
         write (iunit, outfmt) "ONIOM TOTAL ENERGY      ", sccres%e_total, "Eh  "
         write (iunit, outfmt) "ONIOM GRADIENT NORM     ", sccres%gnorm, "Eh/α"
      end if
      write (iunit, '(11x,49("-"))')
   end subroutine write_energy_oniom

   subroutine main_property &
      (iunit, env, mol, wfx, basis, xtbData, res, solvModel, acc)

      use xtb_mctc_convert

!! ========================================================================
!  load class definitions
      use xtb_type_molecule
      use xtb_type_wavefunction
      use xtb_type_environment
      use xtb_type_basisset
      use xtb_type_data
      use xtb_type_param
      use xtb_solv_model
      use xtb_solv_gbsa, only: TBorn
      use xtb_xtb_data
      use xtb_intgrad

!! ========================================================================
!  global storage of options, parameters and basis set
      use xtb_setparam

!! ------------------------------------------------------------------------
      use xtb_aespot
      use xtb_dtrafo

      implicit none

!! ========================================================================
      integer, intent(in) :: iunit ! file handle (usually output_unit=6)
!  molecule data
      type(TMolecule), intent(in) :: mol
      type(TEnvironment), intent(inout) :: env
      type(TxTBData), intent(in) :: xtbData
      real(wp), intent(in) :: acc      ! accuracy of integral calculation
      type(TWavefunction), intent(inout) :: wfx
      type(TBasisset), intent(in) :: basis
      type(scc_results), intent(in) :: res
      type(TSolvModel), allocatable, intent(in) :: solvModel

      real(wp), allocatable :: S(:, :)     ! overlap integrals
      real(wp), allocatable :: dpint(:, :, :) ! dipole integrals
      real(wp), allocatable :: qpint(:, :, :) ! quadrupole integrals
      real(wp), allocatable :: C(:, :)     ! molecular orbitals
      real(wp), allocatable :: emo(:)     ! orbital energies
      real(wp), allocatable :: focc(:)    ! fractional occupation numbers
      integer :: ifile
      integer :: ndim, ndp, nqp
      real(wp) :: dip, dipol(3)
      real(wp) :: intcut, neglect
      real(wp), parameter :: trans(3, 1) = 0.0_wp

      type(TBorn) :: gbsa

!  primitive cut-off
      intcut = 25.0_wp - 10.0 * log10(acc)
      intcut = max(20.0_wp, intcut)
!  integral neglect threshold
      neglect = 10.0d-9 * acc
      ndim = basis%nao * (basis%nao + 1) / 2
      allocate (S(basis%nao, basis%nao), dpint(3, basis%nao, basis%nao), &
         & qpint(6, basis%nao, basis%nao), source=0.0_wp)
#ifdef XTB_GPU
      call sdqint_gpu(xtbData%nShell, xtbData%hamiltonian%angShell, mol%n, mol%at, &
         &        basis%nbf, basis%nao, mol%xyz, trans, intcut, &
         &        basis%caoshell, basis%saoshell, basis%nprim, basis%primcount, &
         &        basis%alp, basis%cont, S, dpint, qpint)
#else
      call sdqint(xtbData%nShell, xtbData%hamiltonian%angShell, mol%n, mol%at, &
         &        basis%nbf, basis%nao, mol%xyz, intcut, &
         &        basis%caoshell, basis%saoshell, basis%nprim, basis%primcount, &
         &        basis%alp, basis%cont, S, dpint, qpint)
#endif

!! orbital energies and occupation
      if (set%pr_eig) then
         write (iunit, '(/,4x,"*",1x,a)') "Orbital Energies and Occupations"
         call print_orbital_eigenvalues(iunit, wfx, 11)
      end if

!! Mulliken and CM5 charges
      if (set%pr_mulliken .and. set%gfn_method == 1) then
         call print_mulliken(iunit, mol%n, mol%at, mol%sym, mol%xyz, mol%z, &
            & basis%nao, S, wfx%P, basis%aoat2, basis%lao2)
      end if
      if (set%pr_charges) then
         call open_file(ifile, 'charges', 'w')
         call print_charges(ifile, mol%n, wfx%q)
         call close_file(ifile)
      end if

      ! GBSA information
      if (allocated(solvModel) .and. set%pr_gbsa) then
         call newBornModel(solvModel, env, gbsa, mol%at)
         call gbsa%update(env, mol%at, mol%xyz)
         call print_gbsa_info(iunit, mol%sym, gbsa)
      end if

!! D4 molecular dispersion printout
      if ((set%newdisp .and. set%gfn_method == 2) .and. set%pr_mulliken) then
         call print_molpol(iunit, mol%n, mol%at, mol%sym, mol%xyz, wfx%q, &
            & xtbData%dispersion%wf, xtbData%dispersion%g_a, xtbData%dispersion%g_c, &
            & xtbData%dispersion%dispm)
      end if
      if (set%gfn_method == 0 .and. set%pr_mulliken) then
         call print_molpol(iunit, mol%n, mol%at, mol%sym, mol%xyz, wfx%q, &
            & xtbData%dispersion%wf, xtbData%dispersion%g_a, xtbData%dispersion%g_c,&
            & xtbData%dispersion%dispm)
      end if

!! Spin population
      if (set%pr_spin_population .and. wfx%nopen /= 0) then
         call print_spin_population(iunit, mol%n, mol%at, mol%sym, basis%nao, wfx%focca,&
            & wfx%foccb, S, wfx%C, basis%aoat2, basis%lao2)
      end if

      if (set%pr_fod_pop) then
         call open_file(ifile, 'fod', 'w')
         call print_fod_population(iunit, ifile, mol%n, mol%at, mol%sym, basis%nao, S, &
            & wfx%C, set%etemp, wfx%emo, wfx%ihomoa, wfx%ihomob, basis%aoat2, basis%lao2)
         call close_file(ifile)
      end if

!! wiberg bond orders
      if (set%pr_wiberg) then
         call open_file(ifile, 'wbo', 'w')
         call print_wbofile(ifile, mol%n, wfx%wbo, 0.1_wp)
         call close_file(ifile)
         call print_wiberg(iunit, mol%n, mol%at, mol%sym, wfx%wbo, 0.1_wp)

         call checkTopology(iunit, mol, wfx%wbo, 1)
      end if

      if (set%pr_wbofrag) &
         call print_wbo_fragment(iunit, mol%n, mol%at, wfx%wbo, 0.1_wp)

!! molden file
      if (set%pr_molden_input) then
         allocate (C(basis%nbf, basis%nao), focc(basis%nao), emo(basis%nao), source=0.0_wp)
         if (basis%nbf == basis%nao) then
            C = wfx%C
         else
            call sao2cao(basis%nao, wfx%C, basis%nbf, C, basis)
         end if
         emo = wfx%emo * evtoau
         focc = wfx%focca + wfx%foccb
         call printmold(mol%n, basis%nao, basis%nbf, mol%xyz, mol%at, C, emo, focc, 2.0_wp, basis)
         write (iunit, '(/,"MOs/occ written to file <molden.input>",/)')
         deallocate (C, focc, emo)
      end if

      if (set%pr_gbw) &
         call wrgbw(xtbData, mol%n, mol%at, mol%xyz, mol%z, basis, wfx)

      if (set%pr_tmbas .or. set%pr_tmmos) then
         call open_file(ifile, 'basis', 'w')
         call write_tm_basis(ifile, xtbData, mol%n, mol%at, basis, wfx)
         call close_file(ifile)
      end if

      if (set%pr_tmmos) then
         call open_file(ifile, 'mos', 'w')
         call write_tm_mos(ifile, mol%n, mol%at, basis, wfx)
         call close_file(ifile)
      end if

!! multipole moment prinout
      if (set%pr_dipole) then
         if (set%gfn_method > 1) then
            ! print overall multipole moment
            call molmom(iunit, mol%n, mol%xyz, wfx%q, wfx%dipm, wfx%qp, dip, dipol)
            write (iunit, '(a)')
         else
            call print_dipole(iunit, mol%n, mol%at, mol%xyz, mol%z, wfx%nao, wfx%P, dpint)
         end if
      end if

   end subroutine main_property

   !> wrapper for tblite-PTB property output
   subroutine ptb_property&
                     (iunit, env, chk, calc, mol, res)
      
      use xtb_type_molecule, only: TMolecule
      use xtb_type_restart, only: TRestart
      use xtb_type_environment,  only: TEnvironment
      use xtb_type_data, only: scc_results
      use xtb_type_calculator, only: TCalculator
      use xtb_ptb_calculator, only: TPTBCalculator

      integer, intent(in) :: iunit
      type(TMolecule), intent(in) :: mol
      type(TEnvironment), intent(inout) :: env
      type(TRestart),  intent(inout) :: chk
      type(TPTBCalculator), intent(in) :: calc
      type(scc_results), intent(in) :: res
      
#if WITH_TBLITE
   call tblite_ptb_property(iunit, env, chk%tblite, calc%bas, mol, chk%wfn, res)
#endif

   end subroutine ptb_property

#if WITH_TBLITE
   subroutine tblite_ptb_property &
      (iunit, env, wfn, bas, struc, wfx, res)

      use xtb_mctc_convert
      use xtb_type_molecule
      use xtb_type_wavefunction
      use xtb_type_environment
      use xtb_type_basisset
      use xtb_type_data

      !========================================================================
      !> global storage of options, parameters and basis set
      use xtb_setparam

      !========================================================================
      !> PTB specific property output
      use xtb_ptb_property, only: print_charges_to_screen
      use xtb_ptb_guess, only: get_psh_from_qsh

      use mctc_io_structure, only: structure_type

      implicit none

      !========================================================================
      integer, intent(in) :: iunit ! file handle (usually output_unit=6)
      !> tblite data formats
      type(structure_type) :: mol
      type(wavefunction_type), intent(in) :: wfn
      type(basis_type), intent(in) :: bas
      !> molecule data
      type(TMolecule), intent(in) :: struc
      type(TEnvironment), intent(inout) :: env
      type(TWavefunction), intent(inout) :: wfx
      type(scc_results), intent(in) :: res
      integer :: ifile, i
      real(wp), allocatable :: psh(:, :)
      real(wp) :: dip, isotropic_alpha

      mol = struc

      !> orbital energies and occupation
      if (set%pr_eig) then
         write (iunit, '(/,4x,"*",1x,a)') "Orbital Energies and Occupations"
         call print_orbital_eigenvalues(iunit, wfx, 11)
      end if

      !> Mixed Mulliken-Loewdin atomic charges and shell populations
      allocate (psh(bas%nsh, wfn%nspin), source=0.0_wp)
      psh = get_psh_from_qsh(wfn, bas)
      call print_charges_to_screen(iunit, mol, bas, wfn%qat, psh)
      if (set%pr_charges) then
         call open_file(ifile, 'charges', 'w')
         call print_charges(ifile, struc%n, wfx%q)
         call close_file(ifile)
      end if

      !> Spin population
      ! if (set%pr_spin_population .and. wfx%nopen .ne. 0) then
      !    call print_spin_population(iunit, mol%n, mol%at, mol%sym, basis%nao, wfx%focca,&
      !       & wfx%foccb, S, wfx%C, basis%aoat2, basis%lao2)
      ! end if

!! wiberg bond orders
      if (set%pr_wiberg) then
         call open_file(ifile, 'wbo', 'w')
         call print_wbofile(ifile, struc%n, wfx%wbo, 0.1_wp)
         call close_file(ifile)
         call print_wiberg(iunit, struc%n, struc%at, struc%sym, wfx%wbo, 0.1_wp)

         call checkTopology(iunit, struc, wfx%wbo, 1)
      end if

      if (set%pr_wbofrag) &
         call print_wbo_fragment(iunit, struc%n, struc%at, wfx%wbo, 0.1_wp)

      ! if (set%pr_tmmos) then
      !    call open_file(ifile, 'mos', 'w')
      !    call write_tm_mos(ifile, struc%n, struc%at, basis, wfx)
      !    call close_file(ifile)
      ! end if

      dip = norm2(res%dipole)

      write (iunit, '(a)')
      write (iunit, '(1x)', advance="no")
      do i = 1, 38
         write (iunit, '(a)', advance="no") "-"
      end do
      write (iunit, '(/)', advance="no")
      write (iunit, '(4x,"Molecular dipole moment (a.u.)")')
      write (iunit, '(4x,"X        Y        Z")')
      write (iunit, '(1x)', advance="no")
      do i = 1, 38
         write (iunit, '(a)', advance="no") "-"
      end do
      write (iunit, '(/)', advance="no")
      write (iunit, '(1x,3f9.4)') &
            & res%dipole(1), res%dipole(2), res%dipole(3)
      write (iunit, '(1x)', advance="no")
      do i = 1, 38
         write (iunit, '(a)', advance="no") "-"
      end do
      write (iunit, '(/)', advance="no")
      write (iunit, '(4x,"Total dipole moment (a.u. / Debye):",/,1x,2f9.4)') &
            & dip, dip * autod

      write (iunit, '(a)')
      write (iunit, '(1x)', advance="no")
      do i = 1, 38
         write (iunit, '(a)', advance="no") "-"
      end do
      write (iunit, '(/)', advance="no")
      write (iunit, '(4x,"Molecular quadrupole tensor: (a.u.)")')
      write (iunit, '(9x,"X         Y         Z")')
      write (iunit, '(4x,a,f10.4)') "X", res%quadrupole(1)
      write (iunit, '(4x,a,2f10.4)') "Y", res%quadrupole(2:3)
      write (iunit, '(4x,a,3f10.4)') "Z", res%quadrupole(4:6)
      write (iunit, '(1x)', advance="no")

      if (set%elprop == p_elprop_alpha) then
         isotropic_alpha = (res%alpha(1, 1) + res%alpha(2, 2) + res%alpha(3, 3)) / 3.0_wp
         write (iunit, '(a)')
         write (iunit, '(1x)', advance="no")
         do i = 1, 38
            write (iunit, '(a)', advance="no") "-"
         end do
         write (iunit, '(/)', advance="no")
         write (iunit, '(4x,"Numerical polarizability tensor: (a.u.)")')
         write (iunit, '(9x,"X         Y         Z")')
         write (iunit, '(4x,a,3f10.4)') "X", res%alpha(1, 1:3)
         write (iunit, '(4x,a,3f10.4)') "Y", res%alpha(2, 1:3)
         write (iunit, '(4x,a,3f10.4)') "Z", res%alpha(3, 1:3)
         write (iunit, '(1x)', advance="no")
         do i = 1, 38
            write (iunit, '(a)', advance="no") "-"
         end do
         write (iunit, '(/)', advance="no")
         write (iunit, '(4x,"Total isotropic dipole polarizability (a.u. / Å³):",/,1x,2f9.4)') &
               & isotropic_alpha, isotropic_alpha * (autoaa**3)
      end if

   end subroutine tblite_ptb_property
#endif

   subroutine gfnff_property(iunit, n, xyz, topo, nlist)
      use xtb_gfnff_topology, only: TGFFTopology
      use xtb_gfnff_neighbourlist, only: TGFFNeighbourList
      use xtb_aespot, only: molqdip
!! ========================================================================
      !  global storage of options, parameters and basis set
      use xtb_setparam
      integer, intent(in) :: iunit, n
      real(wp), intent(in) :: xyz(3, n)
      type(TGFFTopology), intent(in) :: topo
      type(TGFFNeighbourList), intent(in) :: nlist

      ! dipole moment from charge
      if (set%pr_dipole) then
         call molqdip(iunit, n, xyz, nlist%q)
      end if

   end subroutine gfnff_property

   subroutine main_cube &
      (lverbose, mol, wfx, basis, res)

      use xtb_mctc_convert

!! ========================================================================
!  load class definitions
      use xtb_type_molecule
      use xtb_type_wavefunction
      use xtb_type_basisset
      use xtb_type_data
      use xtb_type_param

!! ========================================================================
!  global storage of options, parameters and basis set
      use xtb_setparam

!! ------------------------------------------------------------------------
      use xtb_aespot
      use xtb_scc_core
      use esp
      use stm
      use xtb_dtrafo

      implicit none

!! ========================================================================
      logical, intent(in) :: lverbose
!  molecule data
      type(TMolecule), intent(in) :: mol
      type(TWavefunction), intent(in) :: wfx
      type(TBasisset), intent(in) :: basis
      type(scc_results), intent(in) :: res

      real(wp), allocatable :: C(:, :)     ! molecular orbitals
      real(wp), allocatable :: emo(:)     ! orbital energies
      real(wp), allocatable :: focc(:)    ! fractional occupation numbers
      real(wp), allocatable :: focca(:)   ! fractional occupation numbers (alpha)
      real(wp), allocatable :: foccb(:)   ! fractional occupation numbers (beta)
      integer :: ndim, ndp, nqp
      real(wp) :: dip, dipol(3)
      real(wp) :: acc, intcut, neglect
      real(wp) :: efa, efb, ga, gb, nfoda, nfodb

!! ------------------------------------------------------------------------
!  FOD
      if (set%pr_fod) then
         allocate (C(basis%nbf, basis%nao), focca(basis%nao), foccb(basis%nao), focc(basis%nao), emo(basis%nao), &
                   source=0.0_wp)
         if (wfx%ihomoa + 1 <= basis%nao) &
            call fermismear(.false., basis%nao, wfx%ihomoa, set%etemp, wfx%emo, focca, nfoda, efa, ga)
         if (wfx%ihomob + 1 <= basis%nao) &
            call fermismear(.false., basis%nao, wfx%ihomob, set%etemp, wfx%emo, foccb, nfodb, efb, gb)
         emo = wfx%emo * evtoau
         call fodenmak(.true., basis%nao, emo, focca, efa)
         call fodenmak(.true., basis%nao, emo, foccb, efb)
         focc = focca + foccb
         if (basis%nbf == basis%nao) then
            C = wfx%C
         else
            call sao2cao(basis%nao, wfx%C, basis%nbf, C, basis)
         end if
         if (lverbose) &
            write (stdout, '(/,"FOD written to file: ''fod.cub''",/)')
         call cube(mol%n, basis%nao, basis%nbf, mol%xyz, mol%at, C, emo, focc, 'fod.cub', basis)
         deallocate (C, focca, foccb, focc, emo)
      end if

!! ------------------------------------------------------------------------
!  print spin density to cube file
      if (set%pr_spin_density .and. wfx%nopen /= 0) then
         allocate (C(basis%nbf, basis%nao), focc(basis%nao), emo(basis%nao), source=0.0_wp)
         if (basis%nbf == basis%nao) then
            C = wfx%C
         else
            call sao2cao(basis%nao, wfx%C, basis%nbf, C, basis)
         end if
         if (lverbose) &
            write (stdout, '(/,"(R)spin-density written to file: ''spindensity.cub''",/)')
         emo = wfx%emo * evtoau
         focc = wfx%focca - wfx%foccb
         call cube(mol%n, basis%nao, basis%nbf, mol%xyz, mol%at, C, emo, focc, 'spindensity.cub', basis)
         deallocate (C, focc, emo)
      end if

!! ------------------------------------------------------------------------
!  print density to cube file
      if (set%pr_density) then
         allocate (C(basis%nbf, basis%nao), emo(basis%nao), source=0.0_wp)
         if (basis%nbf == basis%nao) then
            C = wfx%C
         else
            call sao2cao(basis%nao, wfx%C, basis%nbf, C, basis)
         end if
         if (lverbose) &
            write (stdout, '(/,"density written to file: ''density.cub''",/)')
         emo = wfx%emo * evtoau
         call cube(mol%n, basis%nao, basis%nbf, mol%xyz, mol%at, C, emo, wfx%focc, 'density.cub', basis)
         deallocate (C, emo)
      end if

!! ------------------------------------------------------------------------
!  make an ESP plot
      if (set%pr_esp) then
         allocate (C(basis%nbf, basis%nao), source=0.0_wp)
         if (basis%nbf == basis%nao) then
            C = wfx%C
         else
            call sao2cao(basis%nao, wfx%C, basis%nbf, C, basis)
         end if
         call espplot(mol%n, basis%nao, basis%nbf, mol%at, mol%xyz, mol%z, wfx%focc, C, basis)
         deallocate (C)
      end if

!! ------------------------------------------------------------------------
!  make a STM image
      if (set%pr_stm) then
         allocate (C(basis%nbf, basis%nao), focc(basis%nao), source=0.0_wp)
         if (basis%nbf == basis%nao) then
            C = wfx%C
         else
            call sao2cao(basis%nao, wfx%C, basis%nbf, C, basis)
         end if
         if (wfx%ihomoa + 1 <= wfx%nao) &
            call fermismear(.false., basis%nao, wfx%ihomoa, set%etemp, wfx%emo, focc, nfoda, efa, ga)
         if (wfx%ihomob + 1 <= wfx%nao) &
            call fermismear(.false., basis%nao, wfx%ihomob, set%etemp, wfx%emo, focc, nfodb, efb, gb)
         call stmpic(mol%n, basis%nao, basis%nbf, mol%at, mol%xyz, C, 0.5_wp * (efa + efb), wfx%emo, basis)
         deallocate (C, focc)
      end if

   end subroutine main_cube

   subroutine main_freq &
      (iunit, mol, wfx, res)

      use xtb_mctc_convert

!! ========================================================================
!  load class definitions
      use xtb_type_molecule
      use xtb_type_wavefunction
      use xtb_type_basisset
      use xtb_type_data
      use xtb_type_param

!! ========================================================================
!  global storage of options, parameters and basis set
      use xtb_setparam
      use xtb_splitparam, only: atmass

!! ------------------------------------------------------------------------
      use xtb_hessian
      use xtb_disp_ncoord
      use xtb_io_writer_turbomole, only: writeNormalModesTurbomole

      implicit none

!! ========================================================================
      integer, intent(in) :: iunit
!  molecule data
      type(TMolecule), intent(inout) :: mol
      type(TWavefunction), intent(in) :: wfx
      type(freq_results), intent(inout) :: res

      integer :: ifile
      integer :: i, ii, j, jj, k, l
      character(len=:), allocatable :: hname
      real(wp), allocatable :: bond(:, :)
      integer, allocatable :: molvec(:)
      real(wp), allocatable :: cn(:)
      real(wp), allocatable :: xyz0(:, :)
      real(wp), allocatable :: h(:, :)
      real(wp) :: etot, h298, dum
      integer :: lowmode

      allocate (molvec(mol%n), source=0)
      allocate (xyz0(3, mol%n), h(3 * mol%n, 3 * mol%n), bond(mol%n, mol%n), cn(mol%n), source=0.0_wp)

      if (res%linear) then
         write (iunit, '(1x,a)') 'vibrational frequencies (cm⁻¹)'
      else
         write (iunit, '(1x,a)') 'projected vibrational frequencies (cm⁻¹)'
      end if
      call PREIGF(iunit, res%freq, res%n3true)

      write (iunit, '(1x,a)') 'reduced masses (amu)'
      write (iunit, '(8(i4,'':'',f6.2))') (i, res%rmass(i), i=1, res%n3)
      write (iunit, '(1x,a)') 'IR intensities (km·mol⁻¹)'
      write (iunit, '(8(i4,'':'',f6.2))') (i, res%dipt(i), i=1, res%n3)
      write (iunit, '(1x,a)') 'Raman intensities (Ä⁴*amu⁻¹)'
      write (iunit, '(8(i4,'':'',f6.2))') (i, res%polt(i), i=1, res%n3)

      call open_file(ifile, 'vibspectrum', 'w')
      if (set%elprop == p_elprop_alpha) then
         call write_tm_vibspectrum(ifile, res%n3, res%freq, res%dipt, res%polt,&
                                           set%ptbsetup%raman_temp, set%ptbsetup%raman_lambda)
      else
         call write_tm_vibspectrum(ifile, res%n3, res%freq, res%dipt, res%polt)
      end if
      call close_file(ifile)

      write (iunit, '(1x,a)') 'output can be read by thermo (or use thermo option).'
      write (iunit, '(1x,a)') 'writing <g98.out> molden fake output.'
      write (iunit, '(1x,a)') &
         & 'recommended (thermochemical) frequency scaling factor: 1.0'
      call g98fake2('g98.out', mol%n, mol%at, mol%xyz, res%freq, res%rmass, res%dipt, res%hess)

      if (set%pr_nmtm) then
         call open_file(ifile, "vib_normal_modes", 'w')
         if (ifile /= -1) then
            call writeNormalModesTurbomole(ifile, atmass, res%hess)
            call close_file(ifile)
         end if
      end if

      call generic_header(iunit, "Thermodynamic Functions", 49, 10)
      call print_thermo(iunit, mol%n, res%n3true, mol%at, mol%xyz, res%freq, res%etot, res%htot, res%gtot, &
                        res%nimag, .true., res%zp)
      res%pg = trim(set%pgroup)
      res%temp = set%thermotemp(set%nthermo)
      if (set%enso_mode) then
         call open_file(ifile, "xtb_enso.json", 'w')
         if (ifile /= -1) then
            call enso_printout(ifile, res)
            call close_file(ifile)
         end if
      end if

      ! distort along imags if present
      call distort(mol, res%freq, res%hess)

      if (set%pr_modef .and. (mol%n > 3)) then

         ! do analysis and write mode following file
         call wrmodef(0, mol%n, mol%at, mol%xyz, wfx%wbo, res%rmass, res%freq, res%hess, h, set%mode_vthr, res%linear)

         ! localize the modes
         if (set%mode_vthr > 1.d-6) then
            ! determine molecular fragments
            call ncoord_erf(mol%n, mol%at, mol%xyz, cn)
            call cutcov(mol%n, mol%at, mol%xyz, cn, wfx%wbo, bond)
            call mrec(i, mol%xyz, cn, bond, mol%n, mol%at, molvec)
            call locmode(mol%n, res%n3, mol%at, mol%xyz, set%mode_vthr, res%freq, res%rmass, res%hess, &
                         i, molvec)
            call PREIGF0(iunit, res%freq, res%n3true)
            write (iunit, '("written to xtb_localmodes and g98l.out")')
            call wrmodef(1, mol%n, mol%at, mol%xyz, wfx%wbo, res%rmass, res%freq, res%hess, &
                         h, set%mode_vthr + 200.0_wp, res%linear)
         end if

         call open_file(ifile, '.tmpxtbmodef', 'w')
         write (ifile, *) res%lowmode, res%lowmode
         write (ifile, *) res%etot    ! energy for comparison
         call close_file(ifile)

      end if

   end subroutine main_freq

   subroutine print_charges(ifile, n, q)
      implicit none
      integer, intent(in) :: ifile
      integer, intent(in) :: n
      real(wp), intent(in) :: q(n)
      integer :: i
      if (ifile /= -1) then
         do i = 1, n
            write (ifile, '(f14.8)') q(i)
         end do
      end if
   end subroutine print_charges

   subroutine print_mulliken(iunit, n, at, sym, xyz, z, nao, S, P, aoat2, lao2)
      use xtb_scc_core, only: mpop
      implicit none
      integer, intent(in) :: iunit
      integer, intent(in) :: n
      integer, intent(in) :: at(n)
      character(len=*), intent(in) :: sym(n)
      real(wp), intent(in) :: xyz(3, n)
      real(wp), intent(in) :: z(n)
      integer, intent(in) :: nao
      real(wp), intent(in) :: S(nao, nao)
      real(wp), intent(in) :: P(nao, nao)
      integer, intent(in) :: aoat2(nao)
      integer, intent(in) :: lao2(nao)
      real(wp), allocatable :: q(:)       ! Mulliken partial charges
      real(wp), allocatable :: qlmom(:, :) ! population per shell
      real(wp), allocatable :: cm5(:)     ! CM5 partial charges
      real(wp), allocatable :: cm5a(:)     ! CM5 partial charges
      real(wp), allocatable :: dcm5a(:, :, :)! CM5 partial charges
      integer :: i

      allocate (cm5(n), q(n), qlmom(3, n), cm5a(n), dcm5a(3, n, n), source=0.0_wp)
      call mpop(n, nao, aoat2, lao2, S, P, q, qlmom)
      q = z - q
      call calc_cm5(n, at, xyz, cm5a, dcm5a)
      cm5 = q + cm5a
      write (iunit, '(a)')
      write (iunit, '(2x,"Mulliken/CM5 charges         n(s)   n(p)   n(d)")')
      do i = 1, n
         write (iunit, '(i6,a4,2f9.5,1x,4f7.3)') &
            i, sym(i), q(i), cm5(i), qlmom(1, i), qlmom(2, i), qlmom(3, i)
      end do
   end subroutine print_mulliken

   subroutine print_wbofile(iunit, n, wbo, thr)
      implicit none
      integer, intent(in) :: iunit
      integer, intent(in) :: n
      real(wp), intent(in) :: wbo(n, n)
      real(wp), intent(in) :: thr
      integer :: i, j
      do i = 1, n
         do j = 1, i - 1
            if (wbo(j, i) > thr) write (iunit, *) j, i, wbo(j, i)
         end do
      end do
   end subroutine print_wbofile

   subroutine print_wiberg(iunit, n, at, sym, wbo, thr)
      implicit none
      integer, intent(in) :: iunit
      integer, intent(in) :: n
      integer, intent(in) :: at(n)
      character(len=*), intent(in) :: sym(n)
      real(wp), intent(in) :: wbo(n, n)
      real(wp), intent(in) :: thr

      real(wp), allocatable :: wbr(:, :)
      integer, allocatable :: imem(:)
      integer :: i, j, k, ibmax
      real(wp) :: xsum

      allocate (wbr(n, n), source=wbo)
      allocate (imem(n), source=0)

      write (iunit, '(a)')
      write (iunit, '("Wiberg/Mayer (AO) data.")')
      write (iunit, '("largest (>",f4.2,") Wiberg bond orders for each atom")') thr
      write (iunit, '(a)')
      write (iunit, '(1x,75("-"))')
      write (iunit, '(5x,"#",3x,"Z",1x,"sym",2x,"total",t25,3(5x,"#",1x,"sym",2x,"WBO",2x))')
      write (iunit, '(1x,75("-"))')
      do i = 1, n
         do j = 1, n
            imem(j) = j
         end do
         call wibsort(n, i, imem, wbr)
         ibmax = 0
         xsum = 0.0_wp
         do j = 1, n
            if (wbr(i, j) > thr) ibmax = j
            xsum = xsum + wbr(i, j)
         end do
         if (ibmax > 0) then
            write (iunit, '(i6,1x,i3,1x,a4,f6.3,1x,"--")', advance='no') &
               & i, at(i), sym(i), xsum
         else
            write (iunit, '(i6,1x,i3,1x,a4,f6.3)') &
               & i, at(i), sym(i), xsum
         end if
         do j = 1, ibmax, 3
            if (j > 1) then
               write (iunit, '(t25)', advance='no')
            end if
            do k = j, min(ibmax, j + 2)
               write (iunit, '(i6,1x,a4,f6.3)', advance='no') &
                  & imem(k), sym(imem(k)), wbr(i, k)
            end do
            write (iunit, '(a)')
         end do
      end do
      write (iunit, '(1x,75("-"))')
      write (iunit, '(a)')

      deallocate (wbr, imem)

   contains

      SUBROUTINE wibsort(ncent, imo, imem, qmo)
         implicit none
         integer :: ncent
         integer :: imo
         real(wp) :: qmo(ncent, ncent)
         integer :: imem(ncent)
         integer :: ii, i, j, k, ihilf
         real(wp) :: pp

         do ii = 2, ncent
            i = ii - 1
            k = i
            pp = qmo(imo, i)
            do j = ii, ncent
               if (qmo(imo, j) < pp) cycle
               k = j
               pp = qmo(imo, j)
            end do
            if (k == i) cycle
            qmo(imo, k) = qmo(imo, i)
            qmo(imo, i) = pp

            ihilf = imem(i)
            imem(i) = imem(k)
            imem(k) = ihilf
         end do

      end SUBROUTINE wibsort

   end subroutine print_wiberg

   subroutine print_wbo_fragment(iunit, n, at, wbo, thr)
      use xtb_type_atomlist
      implicit none
      integer, intent(in) :: iunit
      integer, intent(in) :: n
      integer, intent(in) :: at(n)
      real(wp), intent(in) :: wbo(n, n)
      real(wp), intent(in) :: thr

      type(TAtomList) :: atl

      real(wp), allocatable :: bond(:, :)
      integer, allocatable :: cn(:)
      integer, allocatable :: fragment(:)
      integer, allocatable :: list(:)
      character(len=:), allocatable :: string
      integer :: i, j, k, nfrag
      real(wp) :: xsum

      allocate (fragment(n), cn(n), list(n), source=0)
      allocate (bond(n, n), source=0.0_wp)
      where (wbo > thr)
         bond = min(wbo, 1.0_wp)
      elsewhere
         bond = 0.0_wp
      end where
      forall (i=1:n) cn(i) = sum(ceiling(bond(:, i)))

      call mrec(nfrag, cn, bond, n, at, fragment)

      write (iunit, '(a)')
      if (nfrag > 1) then
         write (iunit, '(1x,"*",1x,i0,1x,a)', advance='no') &
            nfrag, "fragments found"
      else
         write (iunit, '(1x,"*",1x,a)', advance='no') &
            "no fragments found"
      end if
      write (iunit, '(1x,"(WBO >",f5.2,")")') thr
      write (iunit, '(a)')
      do i = 1, nfrag
         call atl%new
         call atl%add(fragment == i)
         call atl%to_string(string)
         write (iunit, '(3x,a,"(",i0,"):",1x,a)') "fragment", i, string
      end do

   contains
      subroutine mrec(molcount, cn, bond, n, at, molvec)
         ! molcount: number of total fragments (increased during search)
         ! xyz: overall Cart. coordinates
         ! n: overall number of atoms
         ! at: atomic number array
         ! molvec: assignment vector of atom to fragment
         implicit none
         integer, intent(in) :: cn(n)
         integer, intent(in) :: n, at(n)
         integer, intent(inout) :: molvec(n), molcount
         real(wp), intent(inout) :: bond(n, n)
         logical, allocatable :: taken(:)
         integer :: i
         allocate (taken(n))
         molvec = 0
         molcount = 1
         taken = .false.
         do i = 1, n
            if (.not. taken(i)) then
               molvec(i) = molcount
               taken(i) = .true.
               call neighbours(i, cn, at, taken, n, bond, molvec, molcount)
               molcount = molcount + 1
            end if
         end do
         molcount = molcount - 1
      end subroutine mrec

      recursive subroutine neighbours(i, cn, at, taken, n, bond, &
         &                                molvec, molcnt)
         implicit none
         integer, intent(in) :: cn(n)
         real(wp), intent(inout) :: bond(n, n)
         integer, intent(in) :: i, n, at(n)
         integer, intent(inout) :: molcnt, molvec(n)
         logical, intent(inout) :: taken(n)
         integer :: j, icn, k

         icn = cn(i)
         do k = 1, icn
            j = maxloc(bond(:, i), 1)
            bond(j, i) = 0
            if (i == j) cycle
            if (.not. taken(j)) then
               molvec(j) = molcnt
               taken(j) = .true.
               call neighbours(j, cn, at, taken, n, bond, molvec, molcnt)
            end if
         end do
      end subroutine neighbours

   end subroutine print_wbo_fragment

   subroutine print_molpol(iunit, n, at, sym, xyz, q, wf, g_a, g_c, dispm)
      use xtb_disp_dftd4
      use xtb_disp_ncoord
      use xtb_eeq
      use xtb_type_dispersionmodel
      implicit none
      integer, intent(in) :: iunit
      integer, intent(in) :: n
      integer, intent(in) :: at(n)
      character(len=*), intent(in) :: sym(n)
      real(wp), intent(in) :: xyz(3, n)
      real(wp), intent(in) :: q(n)
      real(wp), intent(in) :: wf
      real(wp), intent(in) :: g_a
      real(wp), intent(in) :: g_c
      type(TDispersionModel), intent(in) :: dispm

      integer :: i
      integer :: dispdim
      real(wp) :: molpol, molc6, molc8
      real(wp), allocatable :: covcn(:)   ! covalent coordination number
      real(wp), allocatable :: gw(:)      ! gaussian weights for references
      real(wp), allocatable :: c6ref(:, :) ! unscaled reference C6
      real(wp), allocatable :: aw(:, :)    ! frequency dependent polarizibilities
      real(wp), allocatable :: c6ab(:, :)  ! actual C6 coeffients

      call d4dim(dispm, n, at, dispdim)
      allocate (covcn(n), aw(23, n), c6ab(n, n), gw(dispdim), &
                c6ref(dispdim, dispdim), source=0.0_wp)

      call ncoord_d4(n, at, xyz, covcn, thr=1600.0_wp)
      call d4(dispm, n, dispdim, at, wf, g_a, g_c, covcn, gw, c6ref)
      call mdisp(dispm, n, dispdim, at, q, xyz, g_a, g_c, gw, c6ref, &
                 molc6, molc8, molpol, aout=aw, cout=c6ab)

      write (iunit, '(a)')
      write (iunit, '("     #   Z     ")', advance='no')
      write (iunit, '("     covCN")', advance='no')
      write (iunit, '("         q")', advance='no')
      write (iunit, '("      C6AA")', advance='no')
      write (iunit, '("      α(0)")', advance='no')
      write (iunit, '(a)')
      do i = 1, n
         write (iunit, '(i6,1x,i3,1x,a4)', advance='no') &
         &     i, at(i), sym(i)
         write (iunit, '(f10.3)', advance='no') covcn(i)
         write (iunit, '(f10.3)', advance='no') q(i)
         write (iunit, '(f10.3)', advance='no') c6ab(i, i)
         write (iunit, '(f10.3)', advance='no') aw(1, i)
         write (iunit, '(a)')
      end do
      write (iunit, '(/,1x,"Mol. C6AA /au·bohr⁶  :",f18.6,'// &
      &            '/,1x,"Mol. C8AA /au·bohr⁸  :",f18.6,'// &
      &            '/,1x,"Mol. α(0) /au        :",f18.6,/)') &
      &             molc6, molc8, molpol

   end subroutine print_molpol

   subroutine print_dipole(iunit, n, at, xyz, z, nao, P, dpint)
      use xtb_mctc_convert
      implicit none
      integer, intent(in) :: iunit
      integer, intent(in) :: n
      integer, intent(in) :: at(n)
      real(wp), intent(in) :: xyz(3, n)
      real(wp), intent(in) :: z(n)
      integer, intent(in) :: nao
      real(wp), intent(in) :: P(nao, nao)
      real(wp), intent(in) :: dpint(3, nao, nao)

      integer :: i, j, k
      real(wp) :: d(3), dip

      ! core part
      d = 0.0_wp
      do i = 1, n
         d = d + xyz(:, i) * z(i)
      end do

      ! contraction with P
      k = 0
      do i = 1, nao
         do j = 1, i - 1
            k = k + 1
            d = d - 2.0_wp * P(j, i) * dpint(:, i, j)
         end do
         k = k + 1
         d = d - P(i, i) * dpint(:, i, i)
      end do

      dip = norm2(d)

      write (iunit, '(a)')
      write (iunit, '(1x,"dipole moment from electron density (au)")')
      write (iunit, '(1x,"    X       Y       Z   ")')
      write (iunit, '(3f9.4,"  total (Debye): ",f8.3)') &
            & d(1), d(2), d(3), dip * autod
      write (iunit, '(a)')

   end subroutine print_dipole

   subroutine print_spin_population(iunit, n, at, sym, nao, focca, foccb, S, C, aoat2, lao2)
      use xtb_scc_core, only: dmat, mpop
      implicit none
      integer, intent(in) :: iunit       ! STDOUT
      integer, intent(in) :: n           ! number of atoms
      integer, intent(in) :: at(n)       ! atom types
      character(len=*), intent(in) :: sym(n) ! atom symbols
      integer, intent(in) :: nao         ! number of spherical atomic orbitals
      real(wp), intent(in) :: focca(nao)  ! fractional occupation numbers (alpha)
      real(wp), intent(in) :: foccb(nao)  ! fractional occupation numbers (beta)
      real(wp), intent(in) :: S(nao, nao)  ! overlap matrix
      real(wp), intent(in) :: C(nao, nao)  ! eigenvector/orbitals
      integer, intent(in) :: aoat2(nao)
      integer, intent(in) :: lao2(nao)

      integer :: i
      real(wp), allocatable :: tmp(:)
      real(wp), allocatable :: q(:)
      real(wp), allocatable :: qlmom(:, :)
      real(wp), allocatable :: X(:, :)

      allocate (tmp(nao), q(n), qlmom(3, n), X(nao, nao), source=0.0_wp)

      write (iunit, '("(R)spin-density population")')
      tmp = focca - foccb
      call dmat(nao, tmp, C, X) ! X is scratch
      call mpop(n, nao, aoat2, lao2, S, X, q, qlmom)
      write (iunit, '(a)')
      write (iunit, '(1x,"Mulliken population  n(s)   n(p)   n(d)")')
      do i = 1, n
         write (iunit, '(i6,a4,1f8.4,1x,4f7.3)') &
            &  i, sym(i), q(i), qlmom(1, i), qlmom(2, i), qlmom(3, i)
      end do

   end subroutine print_spin_population

   subroutine print_fod_population(iunit, ifile, n, at, sym, nao, S, C, etemp, emo, ihomoa, &
         & ihomob, aoat2, lao2)
      use xtb_mctc_convert
      use xtb_scc_core
      implicit none
      integer, intent(in) :: iunit       ! STDOUT
      integer, intent(in) :: ifile       ! file handle for printout of FOD population
      integer, intent(in) :: n           ! number of atoms
      integer, intent(in) :: at(n)       ! atom types
      character(len=*), intent(in) :: sym(n)  ! atom symbols
      integer, intent(in) :: nao         ! number of spherical atomic orbitals
      real(wp), intent(in) :: S(nao, nao)  ! overlap matrix
      real(wp), intent(in) :: C(nao, nao)  ! eigenvector/orbitals
      real(wp), intent(in) :: etemp       ! electronic temperature
      real(wp), intent(in) :: emo(nao)    ! orbital energies
      integer, intent(in) :: ihomoa      ! position of HOMO in alpha space
      integer, intent(in) :: ihomob      ! position of HOMO in beta space
      integer, intent(in) :: aoat2(nao)
      integer, intent(in) :: lao2(nao)

      integer :: i
      real(wp), allocatable :: focc(:)    ! fractional occupation numbers
      real(wp), allocatable :: q(:)       ! FOD populations
      real(wp), allocatable :: qlmom(:, :) ! FOD populations per shell
      real(wp), allocatable :: X(:, :)     ! Loewdin orthonormalizer
      real(wp), allocatable :: focca(:)   ! fractional occupation numbers (alpha)
      real(wp), allocatable :: foccb(:)   ! fractional occupation numbers (beta)
      real(wp) :: efa, efb, ga, gb, nfoda, nfodb

      allocate (q(n), qlmom(3, n), X(nao, nao), focca(nao), foccb(nao), focc(nao), &
                source=0.0_wp)

      call makel(nao, S, C, X)
      if (ihomoa + 1 <= nao) &
         call fermismear(.false., nao, ihomoa, etemp, emo, focca, nfoda, efa, ga)
      if (ihomob + 1 <= nao) &
         call fermismear(.false., nao, ihomob, etemp, emo, foccb, nfodb, efb, gb)
      call fodenmak(.true., nao, emo * evtoau, focca, efa)
      call fodenmak(.true., nao, emo * evtoau, foccb, efb)

      focc = focca + foccb
      write (iunit, '(/,"NFOD :",1x,F10.4)') sum(focc)
      q = 0
      qlmom = 0
      call lpop(n, nao, aoat2, lao2, focc, X, 1.0d0, q, qlmom)
      write (iunit, '(a)')
      write (iunit, '(" Loewdin FODpop      n(s)   n(p)   n(d)")')
      do i = 1, n
         write (iunit, '(i6,a4,f8.4,1x,4f7.3)') &
            i, sym(i), q(i), qlmom(1, i), qlmom(2, i), qlmom(3, i)
      end do
      if (ifile /= -1) then
         do i = 1, n
            write (ifile, '(F14.8)') q(i)
         end do
      end if

   end subroutine print_fod_population

   subroutine print_thermo(iunit, nat, nvib_in, at, xyz, freq, etot, htot, gtot, nimag, pr, zp)
      use xtb_mctc_convert
      use xtb_readin
      use xtb_setparam
      use xtb_axis, only: axis2
      use xtb_thermo
      implicit none
      integer, intent(in) :: iunit
      logical, intent(in) :: pr
      integer, intent(in) :: nat
      integer, intent(in) :: at(nat)
      integer, intent(in) :: nvib_in
      real(wp), intent(in) :: freq(3 * nat)
      real(wp), intent(in) :: xyz(3, nat)
      real(wp), intent(in) :: etot
      real(wp), intent(out) :: gtot
      real(wp), intent(out) :: htot
      real(wp), intent(out) :: zp

      real(wp) xx(10), sthr, temp, scale_factor
      real(wp) aa, bb, cc, vibthr, ithr
      real(wp) escf, symnum, wt, avmom, diff
      real(wp) :: omega, maxfreq, fswitch, lnq_r, lnq_v
      real(wp), allocatable :: et(:), ht(:), gt(:), ts(:)
      integer nn, nvib, i, j, k, n, nvib_theo, isthr
      integer, intent(out) :: nimag
      real(wp), allocatable :: vibs(:), tmp(:)
      character(len=*), parameter :: outfmt = &
                                     '(9x,"::",1x,a,f24.12,1x,a,1x,"::")'
      character(len=*), parameter :: dblfmt = &
                                     '(10x,":",2x,a,f24.7,1x,a,1x,":")'
      character(len=*), parameter :: intfmt = &
                                     '(10x,":",2x,a,i24,       6x,":")'
      character(len=*), parameter :: chrfmt = &
                                     '(10x,":",2x,a,a24,       6x,":")'

      logical linear, atom, da

      allocate (et(set%nthermo), ht(set%nthermo), gt(set%nthermo), ts(set%nthermo), &
         &      vibs(3 * nat), tmp(3 * nat), source=0.0_wp)

      ! frequencies read in are considered
      ! as being real if .gt. this value in cm-1
      ! this threshold requires projected freqs.!
      vibthr = 1.0
      ithr = set%thermo_ithr

      atom = .false.
      linear = .false.
      sthr = set%thermo_sthr
      if (abs(set%thermo_fscal - 1.0_wp) > 1.0e-8_wp) then
         scale_factor = set%thermo_fscal
      else
         if (set%mode_extrun == p_ext_gfnff) then
            scale_factor = 1.03_wp
         else
            scale_factor = 1.0_wp
         end if
      end if
      nvib = 0
      nimag = 0

   call axis2(nat,xyz,aa,bb,cc,avmom,wt)

      nvib_theo = 3 * nat - 6
      if (cc < 1.d-10) linear = .true.
      if (linear) nvib_theo = 3 * nat - 5

      if (aa + bb + cc < 1.d-6) then
         atom = .true.
         nvib = 0
         nvib_theo = 0
      end if

      ! the rotational number
      call getsymmetry(pr, iunit, nat, at, xyz, set%desy, set%maxatdesy, set%pgroup)
      call getsymnum(set%pgroup, linear, symnum)

      vibs = 0
      do i = 1, 3 * nat
         if (abs(freq(i)) > vibthr .and. i <= nvib_in) then
            nvib = nvib + 1
            vibs(nvib) = freq(i)
         end if
      end do

      ! scale
      vibs(1:nvib) = vibs(1:nvib) * scale_factor

      do i = 1, nvib
         ! artifacts
         if (vibs(i) < 0 .and. vibs(i) > ithr) then
            vibs(i) = -vibs(i)
            if (pr) write (iunit, *) 'inverting freq ', i, vibs(i)
         end if
      end do
      tmp = vibs

      k = nvib
      nvib = 0
      j = 0
      diff = abs(maxval(vibs) - set%thermo_sthr)
      do i = 1, k
         if (tmp(i) > 0) then
            nvib = nvib + 1
            if (abs(tmp(i) - set%thermo_sthr) < diff) then
               diff = abs(tmp(i) - set%thermo_sthr)
               isthr = nvib
            end if
            vibs(nvib) = tmp(i) * rcmtoau ! work in atomic units, seriously
         else
            j = j + 1
         end if
      end do
      nimag = j

      if (pr) then
         write (iunit, '(a)')
         write (iunit, '(10x,51("."))')
         write (iunit, '(10x,":",22x,a,22x,":")') "SETUP"
         write (iunit, '(10x,":",49("."),":")')
         write (iunit, intfmt) "# frequencies    ", nvib
         write (iunit, intfmt) "# imaginary freq.", nimag
         write (iunit, chrfmt) "linear?          ", bool2string(linear)
         write (iunit, chrfmt) "only rotor calc. ", bool2string(nvib == 0)
         write (iunit, chrfmt) "symmetry         ", trim(set%pgroup)
         write (iunit, intfmt) "rotational number", int(symnum)
         write (iunit, dblfmt) "scaling factor   ", scale_factor, "    "
         write (iunit, dblfmt) "rotor cutoff     ", set%thermo_sthr, "cm⁻¹"
         write (iunit, dblfmt) "imag. cutoff     ", ithr, "cm⁻¹"
         write (iunit, '(10x,":",49("."),":")')
      end if

      call print_thermo_sthr_ts(iunit, nvib, vibs, avmom, set%thermo_sthr, set%thermotemp(set%nthermo))

      ! do calc.
      zp = 0.5_wp * sum(vibs(1:nvib))
      do i = 1, set%nthermo
         temp = set%thermotemp(i)
         call thermodyn(iunit, aa, bb, cc, avmom, linear, atom, symnum, wt, vibs, nvib, escf, &
            & temp, sthr, et(i), ht(i), gt(i), ts(i), zp, pr)
         !call oldthermo(aa,bb,cc,avmom,linear,atom,symnum,wt,vibs,nvib,escf, &
         !   & temp,sthr,et(i),ht(i),gt(i),ts(i),zp,pr)
      end do

      write (iunit, '(a)')
      write (iunit, '(a10)', advance='no') "T/K"
      write (iunit, '(a16)', advance='no') "H(0)-H(T)+PV"
      write (iunit, '(a16)', advance='no') "H(T)/Eh"
      write (iunit, '(a16)', advance='no') "T*S/Eh"
      write (iunit, '(a16)', advance='no') "G(T)/Eh"
      write (iunit, '(a)')
      write (iunit, '(3x,72("-"))')
      do i = 1, set%nthermo
         write (iunit, '(3f10.2)', advance='no') set%thermotemp(i)
         write (iunit, '(3e16.6)', advance='no') ht(i)
         write (iunit, '(3e16.6)', advance='no') et(i)
         write (iunit, '(3e16.6)', advance='no') ts(i)
         write (iunit, '(3e16.6)', advance='no') gt(i)
         if (i == set%nthermo .and. set%nthermo > 1) then
            write (iunit, '(1x,"(used)")')
         else
            write (iunit, '(a)')
         end if
      end do
      write (iunit, '(3x,72("-"))')

      gtot = gt(set%nthermo)
      htot = et(set%nthermo)

      write (iunit, '(a)')
      write (iunit, '(9x,53(":"))')
      write (iunit, '(9x,"::",18x,a,18x,"::")') "THERMODYNAMIC"
      write (iunit, '(9x,53(":"))')
      write (iunit, outfmt) "total free energy ", gtot + etot, "Eh  "
      write (iunit, '(9x,"::",49("."),"::")')
      write (iunit, outfmt) "total energy      ", etot, "Eh  "
      write (iunit, outfmt) "zero point energy ", zp, "Eh  "
      write (iunit, outfmt) "G(RRHO) w/o ZPVE  ", gtot - zp, "Eh  "
      write (iunit, outfmt) "G(RRHO) contrib.  ", gtot, "Eh  "
      write (iunit, '(9x,53(":"))')

   end subroutine print_thermo

   subroutine print_thermo_sthr_lnq(iunit, nvib, vibs, avmom, sthr, temp)
      use xtb_mctc_convert
      use xtb_thermo
      implicit none
      integer, intent(in) :: iunit
      integer, intent(in) :: nvib
      real(wp), intent(in) :: vibs(nvib)
      real(wp), intent(in) :: avmom
      real(wp), intent(in) :: sthr
      real(wp), intent(in) :: temp

      integer :: i
      real(wp) :: maxfreq, omega, lnq_r, lnq_v, fswitch

      write (iunit, '(a)')
      maxfreq = max(300.0_wp, chg_inverted(0.99_wp, sthr))
      write (iunit, '(a8,a14,a12,10x,a12,10x,a12)') &
         "mode", "ω/cm⁻¹", "ln{qvib}", "ln{qrot}", "ln{qtot}"
      write (iunit, '(3x,72("-"))')
      do i = 1, nvib
         omega = vibs(i) * autorcm
         lnq_r = lnqvib(temp, omega)
         lnq_v = lnqrot(temp, omega, avmom)
         fswitch = 1.0_wp - chg_switching(omega, sthr)
         if (omega > maxfreq) exit
         write (iunit, '(i8,f10.2,2(f12.5,1x,"(",f6.2,"%)"),f12.5)') &
            i, omega, lnq_v, (1.0_wp - fswitch) * 100, &
            lnq_r, fswitch * 100, (1.0_wp - fswitch) * lnq_v + fswitch * lnq_r
      end do
      write (iunit, '(3x,72("-"))')

   end subroutine print_thermo_sthr_lnq

   subroutine print_thermo_sthr_ts(iunit, nvib, vibs, avmom_si, sthr_rcm, temp)
      use xtb_mctc_constants
      use xtb_mctc_convert
      use xtb_thermo
      implicit none

      integer, intent(in) :: iunit      !< output unit, usually STDOUT
      integer, intent(in) :: nvib       !< number of frequencies
      real(wp), intent(in) :: vibs(nvib) !< frequencies in Eh
      real(wp), intent(in) :: avmom_si   !< average moment
      real(wp), intent(in) :: sthr_rcm   !< rotor cutoff
      real(wp), intent(in) :: temp       !< temperature

      integer :: i
      real(wp) :: maxfreq, omega, s_r, s_v, fswitch
      real(wp) :: beta, xxmom, e, ewj, mu, RT, sthr, avmom
      beta = 1.0_wp / kB / temp ! beta in 1/Eh
      sthr = sthr_rcm * rcmtoau ! sthr in Eh
      RT = kb * temp * autokcal ! RT in kcal/mol for printout
      avmom = avmom_si * kgtome * aatoau**2 * 1.0e+20_wp ! in me·α²

      write (iunit, '(a)')
      maxfreq = max(300.0_wp, chg_inverted(0.99_wp, sthr_rcm))
      write (iunit, '(a8,a14,1x,a27,a27,a12)') &
         "mode", "ω/cm⁻¹", "T·S(HO)/kcal·mol⁻¹", "T·S(FR)/kcal·mol⁻¹", "T·S(vib)"
      write (iunit, '(3x,72("-"))')
      do i = 1, nvib
         ! frequency is Eh
         omega = vibs(i)
         ! omega in Eh, beta in 1/Eh
         ewj = exp(-omega * beta)
         ! moment of intertia corresponding to the rotor with frequency omega
         ! mu is in me·α² (au)
         mu = 0.5_wp / (omega + 1.0e-14_wp)
         ! this reduced moment limits the rotational moment of inertia for
         ! this vibration to that of the total molecule rotation/3
         ! avmom and mu are in au
         mu = mu * avmom / (mu + avmom)
         !              free rotor entropy
         ! Cramer, page 328 for one degree of freedom or
         ! http://cccbdb.nist.gov/thermo.asp, eq. 35, sigma=1
         !              harm. osc. entropy
         if (omega > 0) then
            ! this is S/R which is dimensionless
            s_v = omega * beta * ewj / (1.0_wp - ewj) - log(1.0_wp - ewj)
            s_r = 0.5_wp + log(sqrt(pi / beta * 2.0_wp * mu))
         else
            s_v = 0.0_wp
            s_r = 0.0_wp
         end if
         ! Head-Gordon weighting
         fswitch = 1.0_wp - chg_switching(omega, sthr)
         if (omega > maxfreq * rcmtoau) exit
         write (iunit, '(i8,f10.2,2(f12.5,1x,"(",f6.2,"%)"),f12.5)') &
            i, omega * autorcm, -RT * s_v, (1.0_wp - fswitch) * 100, &
            -RT * s_r, fswitch * 100, -RT * ((1.0_wp - fswitch) * s_v + fswitch * s_r)
      end do
      write (iunit, '(3x,72("-"))')

   end subroutine print_thermo_sthr_ts

   subroutine print_gbsa_info(iunit, sym, gbsa)
      use xtb_mctc_constants
      use xtb_mctc_convert
      use xtb_solv_gbsa, only: TBorn
      implicit none
      integer, intent(in) :: iunit
      character(len=*), intent(in) :: sym(:)
      type(TBorn), intent(in) :: gbsa

      integer :: i

      write (iunit, '(a)')
      write (iunit, '(1x,"*",1x,a)') &
         &  "generalized Born model for continuum solvation"
      write (iunit, '(a)')
      if (gbsa%lhb) then
         write (iunit, '(2x,2a4,5x,3a)') "#", "Z", "Born rad/Å", "   SASA/Å²", "    H-bond"
         do i = 1, size(sym)
            write (iunit, '(i6,1x,i3,1x,a4,3f10.3)') &
               &  i, gbsa%at(i), sym(i), &
               &  gbsa%brad(i) * autoaa, gbsa%sasa(i) * fourpi * autoaa**2, &
               &  gbsa%hbw(i)
         end do
      else
         write (iunit, '(2x,2a4,5x,2a)') "#", "Z", "Born rad/Å", "   SASA/Å²"
         do i = 1, size(sym)
            write (iunit, '(i6,1x,i3,1x,a4,2f10.3)') &
               &  i, gbsa%at(i), sym(i), &
               &  gbsa%brad(i) * autoaa, gbsa%sasa(i) * fourpi * autoaa**2
         end do
      end if
      write (iunit, '(/,1x,"total SASA / Å² :",f13.3)') &
         &  sum(gbsa%sasa) * fourpi * autoaa**2

   end subroutine print_gbsa_info

end module xtb_propertyoutput

subroutine print_orbital_eigenvalues(iunit, wfn, range)
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert
   use xtb_type_wavefunction
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: range
   type(TWavefunction), intent(in) :: wfn
   character(len=*), parameter :: hlfmt = '(    a24,f21.7,1x,"Eh",f18.4,1x,"eV")'
   integer :: maxorb, minorb, iorb
   real(wp) :: gap

   minorb = max(wfn%ihomoa - (range + 1), 1)
   maxorb = min(wfn%ihomoa + range, wfn%nao)
   gap = wfn%emo(wfn%ihomoa + 1) - wfn%emo(wfn%ihomoa)

   write (iunit, '(a)')
   write (iunit, '(a10,a14,a21,a21)') "#", "Occupation", "Energy/Eh", "Energy/eV"
   write (iunit, '(6x,61("-"))')
   if (minorb > 1) then
      call write_line(1, wfn%focc, wfn%emo, wfn%ihomo)
      if (minorb > 2) &
         write (iunit, '(a10,a14,a21,a21)') "...", "...", "...", "..."
   end if
   do iorb = minorb, maxorb
      call write_line(iorb, wfn%focc, wfn%emo, wfn%ihomo)
   end do
   if (maxorb < wfn%nao) then
      if (maxorb < wfn%nao - 1) then
         if (wfn%focc(maxorb) > 1.0e-7_wp) then
            write (iunit, '(a10,a14,a21,a21)') "...", "...", "...", "..."
         else
            write (iunit, '(a10,a14,a21,a21)') "...", "", "...", "..."
         end if
      end if
      call write_line(wfn%nao, wfn%focc, wfn%emo, wfn%ihomo)
   end if
   write (iunit, '(6x,61("-"))')
   write (iunit, hlfmt) "HL-Gap", gap * evtoau, gap
   write (iunit, hlfmt) "Fermi-level", (wfn%efa + wfn%efb) / 2 * evtoau, (wfn%efa + wfn%efb) / 2
contains
   subroutine write_line(iorb, focc, emo, ihomo)
      integer, intent(in) :: iorb
      integer, intent(in) :: ihomo
      real(wp), intent(in) :: focc(:)
      real(wp), intent(in) :: emo(:)
      character(len=*), parameter :: mofmt = '(i10,f14.4,f21.7,f21.4)'
      character(len=*), parameter :: vofmt = '(i10,14x,  f21.7,f21.4)'
      if (focc(iorb) < 1.0e-7_wp) then
         write (iunit, vofmt, advance='no') iorb, emo(iorb) * evtoau, emo(iorb)
      else
         write (iunit, mofmt, advance='no') iorb, focc(iorb), emo(iorb) * evtoau, emo(iorb)
      end if
      if (iorb == ihomo) then
         write (iunit, '(1x,"(HOMO)")')
      elseif (iorb == ihomo + 1) then
         write (iunit, '(1x,"(LUMO)")')
      else
         write (iunit, '(a)')
      end if
   end subroutine write_line
end subroutine print_orbital_eigenvalues
