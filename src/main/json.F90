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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

! % cat xtbout.json
! {
! # calculation information
!    "coordinate file": coord,
!    "method": "GFN2-xTB",
!    "parameter file": "/home/awvwgk/.parameter/.param_gfn2.xtb"
! # molecule information
!    "number of atoms": 52,
!    "number of electron": 160,
!    "total charge": 0,
!    "spin": 1,
! # basis set information
!    "number of basisfunctions": 157,
!    "number of atomic orbitals": 155,
!    "number of shells": 85,
! # wavefunction information
!    "partial charges":
!    [
!       <charges>
!    ],
!    "atomic dipole moments":
!    [
!       [<dipx>,<dipy>,<dipz>]
!    ],
!    "atomic quadrupole moments":
!    [
!       [<quadxx>,<quadyy>,<qyadzz>,<quadxy>,<quadxz>,<quadyz>]
!    ],
!    "orbital energies":
!    [
!       <emo>
!    ],
!    "occupation numbers":
!    [
!       <focc>
!    ],
! # program information
!    "version": 6.1
! }
module xtb_main_json
   implicit none

   private

   public :: main_xtb_json, write_json_gfnff_lists
   public :: main_ptb_json

contains

   subroutine main_xtb_json &
      (ijson, mol, wfx, xbas, sccres, freqres)
      use xtb_mctc_accuracy, only: wp

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

      implicit none

!! ========================================================================
      integer, intent(in) :: ijson ! file handle (usually json-file)
!  molecule data
      type(TMolecule), intent(in) :: mol
      type(TWavefunction), intent(in) :: wfx
      type(TBasisset), intent(in) :: xbas
      type(scc_results), intent(in) :: sccres
      type(freq_results), intent(in) :: freqres
      logical :: alpha
      alpha = set%elprop == p_elprop_alpha

      call write_json_header(ijson)
      call write_json_scc_results(ijson, sccres)
      if (freqres%gtot > 0.0_wp) then
         call write_json_thermo(ijson, freqres)
      end if
      call write_json_charges(ijson, wfx)
      if (set%gfn_method == 2) then
         call write_json_dipole_moments(ijson, wfx)
         call write_json_quadrupole_moments(ijson, wfx)
      end if
      call write_json_wavefunction(ijson, wfx)
      if (freqres%n3true > 0) then
         call write_json_frequencies(ijson, freqres)
         call write_json_reduced_masses(ijson, freqres)
         call write_json_intensities(ijson, freqres, alpha)
      end if
      call write_json_footer(ijson)

   end subroutine main_xtb_json

   subroutine write_json_header(ijson)
      integer, intent(in) :: ijson
      write (ijson, '("{")')
   end subroutine write_json_header

   subroutine write_json_footer(ijson)
      use xtb_setparam
      include 'xtb_version.fh'
      integer, intent(in) :: ijson
      character(len=:), allocatable :: cmdline
      integer :: l
      call get_command(length=l)
      allocate (character(len=l) :: cmdline)
      call get_command(cmdline)
      write (ijson, '(3x,''"program call":'',1x,''"'',a,''",'')') cmdline
      write (ijson, '(3x,''"method": "GFN'',i0,''-xTB",'')') set%gfn_method
      write (ijson, '(3x,a)') '"xtb version": "'//version//'"'
      write (ijson, '("}")')
   end subroutine write_json_footer

   subroutine write_json_ptb_footer(ijson)
      use xtb_setparam
      include 'xtb_version.fh'
      integer, intent(in) :: ijson
      character(len=:), allocatable :: cmdline
      integer :: l
      call get_command(length=l)
      allocate (character(len=l) :: cmdline)
      call get_command(cmdline)
      write (ijson, '(3x,''"program call":'',1x,''"'',a,''",'')') cmdline
      write (ijson, '(3x,''"method": "PTB",'')')
      write (ijson, '(3x,a)') '"xtb version": "'//version//'"'
      write (ijson, '("}")')
   end subroutine write_json_ptb_footer

   subroutine write_json_scc_results(ijson, sccres)
      use xtb_type_data
      integer, intent(in) :: ijson
      type(scc_results), intent(in) :: sccres
      character(len=*), parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
      write (ijson, jfmtf) 'total energy', sccres%e_total
      write (ijson, jfmtf) 'HOMO-LUMO gap / eV', sccres%hl_gap
      write (ijson, jfmtf) 'electronic energy', sccres%e_elec
      write (ijson, '(3x,''"'',a,''":'',1x,"[",2(f15.8,","),f15.8,"],")') &
         'dipole / a.u.', sccres%dipole
      !write(ijson,jfmtf) 'classical repulsion energy',sccres%e_rep
      !write(ijson,jfmtf) 'isotropic electrostatic energy',sccres%e_es
      !write(ijson,jfmtf) 'anisotropic electrostatic energy',sccres%e_aes
      !write(ijson,jfmtf) 'anisotropic XC energy',sccres%e_axc
      !write(ijson,jfmtf) 'classical halogen bound energy',sccres%e_xb
      !write(ijson,jfmtf) 'Generalized Born free energy',sccres%g_born
      !write(ijson,jfmtf) 'SASA free energy',sccres%g_born
      !write(ijson,jfmtf) 'Hydrogen bound free energy',sccres%g_born
   end subroutine write_json_scc_results

   subroutine write_json_charges(ijson, wfn)
      use xtb_type_wavefunction
      integer, intent(in) :: ijson
      type(TWavefunction), intent(in) :: wfn
      character(len=*), parameter :: jfmta = '(3x,''"'',a,''": ['')'
      character(len=*), parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
      integer :: i
      write (ijson, jfmta) 'partial charges'
      write (ijson, '(3x,f15.8,",")') (wfn%q(i), i=1, wfn%n - 1)
      write (ijson, '(3x,f15.8,"],")') wfn%q(wfn%n)
   end subroutine write_json_charges

   subroutine write_json_bondorder(ijson, mol, wfn)
      use xtb_type_molecule, only: TMolecule
      use xtb_type_wavefunction, only: TWavefunction
      integer, intent(in) :: ijson
      type(TMolecule), intent(in) :: mol
      type(TWavefunction), intent(in) :: wfn
      character(len=*), parameter :: jfmta = '(3x,''"'',a,''": ['')'
      integer :: i, j
      logical :: first
      write (ijson, jfmta) 'bond orders'
      do i = 1, mol%n - 1
         do j = i, mol%n
            if (wfn%wbo(j, i) > 0.01) then
               if (first) then
                  write (ijson, '(a)') ','
               end if
               first = .true.
               write (ijson, '(3x,"[ ",i5,",",i5,",",f8.4,"]")', advance='no') i, j, wfn%wbo(j, i)
            end if
         end do
      end do
      write (ijson, '(a,/)', advance='no') '],'
   end subroutine write_json_bondorder

   subroutine write_json_dipole_moments(ijson, wfn)
      use xtb_type_wavefunction
      integer, intent(in) :: ijson
      type(TWavefunction), intent(in) :: wfn
      character(len=*), parameter :: jfmta = '(3x,''"'',a,''": ['')'
      character(len=*), parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
      integer :: i, j
      write (ijson, jfmta) 'atomic dipole moments'
      do i = 1, wfn%n - 1
         write (ijson, '(3x,"[",2(f15.8,","),f15.8,"],")') (wfn%dipm(j, i), j=1, 3)
      end do
      write (ijson, '(3x,"[",2(f15.8,","),f15.8,"]],")') (wfn%dipm(j, wfn%n), j=1, 3)
   end subroutine write_json_dipole_moments

   subroutine write_json_quadrupole_moments(ijson, wfn)
      use xtb_type_wavefunction
      integer, intent(in) :: ijson
      type(TWavefunction), intent(in) :: wfn
      character(len=*), parameter :: jfmta = '(3x,''"'',a,''": ['')'
      character(len=*), parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
      integer :: i, j
      write (ijson, jfmta) 'atomic quadrupole moments'
      do i = 1, wfn%n - 1
         write (ijson, '(3x,"[",5(f15.8,","),f15.8,"],")') (wfn%qp(j, i), j=1, 6)
      end do
      write (ijson, '(3x,"[",5(f15.8,","),f15.8,"]],")') (wfn%qp(j, wfn%n), j=1, 6)
   end subroutine write_json_quadrupole_moments

   subroutine write_json_wavefunction(ijson, wfn)
      use xtb_type_wavefunction
      integer, intent(in) :: ijson
      type(TWavefunction), intent(in) :: wfn
      character(len=*), parameter :: jfmta = '(3x,''"'',a,''": ['')'
      character(len=*), parameter :: jfmti = '(3x,''"'',a,''":'',1x,i0,",")'
      character(len=*), parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
      integer :: i
      write (ijson, jfmti) 'number of molecular orbitals', wfn%nao
      write (ijson, jfmti) 'number of electrons', wfn%nel
      write (ijson, jfmti) 'number of unpaired electrons', wfn%nopen
      write (ijson, jfmta) 'orbital energies / eV'
      write (ijson, '(3x,f15.8,",")') (wfn%emo(i), i=1, wfn%nao - 1)
      write (ijson, '(3x,f15.8,"],")') wfn%emo(wfn%nao)
      write (ijson, jfmta) 'fractional occupation'
      write (ijson, '(3x,f15.8,",")') (wfn%focc(i), i=1, wfn%nao - 1)
      write (ijson, '(3x,f15.8,"],")') wfn%focc(wfn%nao)
   end subroutine write_json_wavefunction

   subroutine write_json_ptb_wavefunction(ijson, wfn)
      use xtb_type_wavefunction
      integer, intent(in) :: ijson
      type(TWavefunction), intent(in) :: wfn
      character(len=*), parameter :: jfmta = '(3x,''"'',a,''": ['')'
      character(len=*), parameter :: jfmti = '(3x,''"'',a,''":'',1x,i0,",")'
      character(len=*), parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
      integer :: i, max_print
      max_print = min(wfn%nao, wfn%ihomo + 8)
      write (ijson, jfmti) 'number of molecular orbitals', wfn%nao
      write (ijson, jfmti) 'number of electrons', wfn%nel
      write (ijson, jfmti) 'number of unpaired electrons', wfn%nopen
      write (ijson, jfmta) 'orbital energies / eV'
      write (ijson, '(3x,f15.8,",")') (wfn%emo(i), i=1, max_print - 1)
      write (ijson, '(3x,f15.8,"],")') wfn%emo(max_print)
      write (ijson, jfmta) 'fractional occupation'
      write (ijson, '(3x,f15.8,",")') (wfn%focc(i), i=1, max_print - 1)
      write (ijson, '(3x,f15.8,"],")') wfn%focc(max_print)
   end subroutine write_json_ptb_wavefunction

   subroutine write_json_thermo(ijson, freqres)
      use xtb_type_data
      integer, intent(in) :: ijson
      type(freq_results), intent(in) :: freqres
      character(len=*), parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
      write (ijson, jfmtf) 'total enthalpy', freqres%htot
      write (ijson, jfmtf) 'total free energy', freqres%gtot
   end subroutine write_json_thermo

   subroutine write_json_frequencies(ijson, freqres)
      use xtb_type_data
      integer, intent(in) :: ijson
      type(freq_results), intent(in) :: freqres
      character(len=*), parameter :: jfmta = '(3x,''"'',a,''": ['')'
      integer :: i
      write (ijson, jfmta) 'vibrational frequencies / rcm'
      write (ijson, '(3x,f15.8,",")') (freqres%freq(i), i=1, freqres%n3true - 1)
      write (ijson, '(3x,f15.8,"],")') freqres%freq(freqres%n3true)
   end subroutine write_json_frequencies

   subroutine write_json_intensities(ijson, freqres, printalpha)
      use xtb_type_data
      use xtb_mctc_accuracy, only: wp
      integer, intent(in) :: ijson
      type(freq_results), intent(in) :: freqres
      logical, intent(in) :: printalpha
      character(len=*), parameter :: jfmta = '(3x,''"'',a,''": ['')'
      integer :: i
      write (ijson, jfmta) 'IR intensities / km/mol'
      do i = 1, freqres%n3true - 1
         if (abs(freqres%freq(i)) < 1.0e-2_wp) then
            write (ijson, '(3x,f15.8,",")') 0.0_wp
         else
            write (ijson, '(3x,f15.8,",")') freqres%dipt(i)
         end if
      end do
      write (ijson, '(3x,f15.8,"],")') freqres%dipt(freqres%n3true)
      if (printalpha) then
         write (ijson, jfmta) 'Raman activities / A^4/amu'
         do i = 1, freqres%n3true - 1
            if (abs(freqres%freq(i)) < 1.0e-2_wp) then
               write (ijson, '(3x,f15.8,",")') 0.0_wp
            else
               write (ijson, '(3x,f15.8,",")') freqres%polt(i)
            end if
         end do
         write (ijson, '(3x,f15.8,"],")') freqres%polt(freqres%n3true)
      end if
   end subroutine write_json_intensities

   subroutine write_json_reduced_masses(ijson, freqres)
      use xtb_type_data
      integer, intent(in) :: ijson
      type(freq_results), intent(in) :: freqres
      character(len=*), parameter :: jfmta = '(3x,''"'',a,''": ['')'
      integer :: i
      write (ijson, jfmta) 'reduced masses'
      write (ijson, '(3x,f15.8,",")') (freqres%rmass(i), i=1, freqres%n3true - 1)
      write (ijson, '(3x,f15.8,"],")') freqres%rmass(freqres%n3true)
   end subroutine write_json_reduced_masses

   subroutine write_json_gfnff_lists(n, etot, gnorm, topo, neigh, nlist, printTopo)
      use xtb_gfnff_topology, only: TGFFTopology
      use xtb_gfnff_neighbourlist, only: TGFFNeighbourList
      use xtb_gfnff_topology, only: TPrintTopo
      use xtb_mctc_accuracy, only: wp
      use xtb_gfnff_neighbor
      include 'xtb_version.fh'
      !> gfnff topology lists
      type(TGFFTopology), intent(in) :: topo
      !> gfnff neighbourlist
      type(TNeigh) :: neigh
      !> gfnff neighbourlist
      type(TGFFNeighbourList), intent(in) :: nlist
      !> topology printout booleans
      type(TPrintTopo), intent(in) :: printTopo
      !> total energy and gradient norm
      real(wp), intent(in) :: etot, gnorm
      character(len=:), allocatable :: cmdline
      integer :: iunit, i, j, n, l

      call open_file(iunit, 'gfnff_lists.json', 'w')
      ! header
      write (iunit, '("{")')
      ! lists printout
      if (printTopo%etot) then ! total energy is scalar
         write (iunit, '(3x,''"total energy":'',f25.15,",")') etot
      end if
      if (printTopo%gnorm) then ! gradient norm is scalar
         write (iunit, '(3x,''"gradient norm":'',f25.15,",")') gnorm
      end if
      if (printTopo%nb) then ! nb(numnb, n, numctr)
         write (iunit, '(3x,''"nb":'',"[")') ! open nb
         if (neigh%numctr == 1) then
            do j = 1, n - 1
               write (iunit, '(3x,"[",*(i7,:,","))', advance='no') neigh%nb(:, j, 1) ! open nb entry
               write (iunit, '("],")') ! close nb entry
            end do
            write (iunit, '(3x,"[",*(i7,:,","),"]",/)', advance='no') neigh%nb(:, n, 1)
            write (iunit, '("]")')
         else ! periodic boundary conditions
            do i = 1, neigh%numctr - 1 ! iterate over all cells
               write (iunit, '(3x,"[")') ! open cell
               do j = 1, n - 1
                  write (iunit, '(3x,"[",*(i7,:,","))', advance='no') neigh%nb(:, j, i)
                  write (iunit, '("],")')
               end do
               write (iunit, '(3x,"[",*(i7,:,","),"]",/)', advance='no') neigh%nb(:, n, i)
               write (iunit, '("]")')
               write (iunit, '(3x,"],")') ! close cell
            end do
            write (iunit, '(3x,"[")') ! open last cell
            do j = 1, n - 1
               write (iunit, '(3x,"[",*(i7,:,","))', advance='no') neigh%nb(:, j, neigh%numctr)
               write (iunit, '("],")')
            end do
            write (iunit, '(3x,"[",*(i7,:,","),"]",/)', advance='no') neigh%nb(:, n, neigh%numctr)
            write (iunit, '("]")')
            write (iunit, '(3x,"]")') ! close last cell
         end if
         write (iunit, '(3x,"],")') ! close nb
      end if
      ! bpair(j,i,iTr) number bonds between i and j when j is translated by iTr
      if (printTopo%bpair) then
         write (iunit, '(3x,''"bpair":'',"[")') ! open bpair
         if (neigh%numctr == 1) then
            do i = 1, n - 1
               write (iunit, '(3x,"[",*(i7,:,","))', advance='no') neigh%bpair(:, i, 1) ! open entry
               write (iunit, '("],")') ! close entry
            end do
            write (iunit, '(3x,"[",*(i7,:,","),"]",/)', advance='no') neigh%bpair(:, n, 1)
            write (iunit, '("]")')
         else ! periodic boundary conditions
            do i = 1, neigh%numctr - 1 ! iterate over all cells
               write (iunit, '(3x,"[")') ! open cell
               do j = 1, n - 1
                  write (iunit, '(3x,"[",*(i7,:,","))', advance='no') neigh%bpair(:, j, i)
                  write (iunit, '("],")')
               end do
               write (iunit, '(3x,"[",*(i7,:,","),"]",/)', advance='no') neigh%bpair(:, n, i)
               write (iunit, '("]")')
               write (iunit, '(3x,"],")') ! close cell
            end do
            write (iunit, '(3x,"[")') ! open last cell
            do j = 1, n - 1
               write (iunit, '(3x,"[",*(i7,:,","))', advance='no') neigh%bpair(:, j, neigh%numctr)
               write (iunit, '("],")')
            end do
            write (iunit, '(3x,"[",*(i7,:,","),"]",/)', advance='no') neigh%bpair(:, n, neigh%numctr)
            write (iunit, '("]")')
            write (iunit, '(3x,"]")') ! close last cell
         end if
         write (iunit, '(3x,"],")') ! close bpair
      end if
      if (printTopo%alist) then ! alist(3,nangl)
         write (iunit, '(3x,''"alist":'',"[")')
         do j = 1, topo%nangl - 1
            write (iunit, '(3x,"[",*(i8,:,","))', advance='no') topo%alist(:, j)
            write (iunit, '("],")')
         end do
         write (iunit, '(3x,"[",*(i8,:,","),"]",/)', advance='no') topo%alist(:, topo%nangl)
         write (iunit, '("]")')
         write (iunit, '(3x,"],")')
      end if
      if (printTopo%blist) then ! blist(3,nbond)
         write (iunit, '(3x,''"blist":'',"[")')
         do j = 1, neigh%nbond - 1
            write (iunit, '(3x,"[",*(i8,:,","))', advance='no') neigh%blist(:, j)
            write (iunit, '("],")')
         end do
         write (iunit, '(3x,"[",*(i8,:,","),"]",/)', advance='no') neigh%blist(:, neigh%nbond)
         write (iunit, '("]")')
         write (iunit, '(3x,"],")')
      end if
      if (printTopo%tlist) then ! tlist(5,ntors)
         write (iunit, '(3x,''"tlist":'',"[")')
         do j = 1, topo%ntors - 1
            write (iunit, '(3x,"[",*(i8,:,","))', advance='no') topo%tlist(:, j)
            write (iunit, '("],")')
         end do
         write (iunit, '(3x,"[",*(i8,:,","),"]",/)', advance='no') topo%tlist(:, topo%ntors)
         write (iunit, '("]")')
         write (iunit, '(3x,"],")')
      end if
      if (printTopo%vtors) then ! vtors(2,ntors)
         write (iunit, '(3x,''"vtors":'',"[")')
         do j = 1, topo%ntors - 1
            write (iunit, '(3x,"[",*(f25.15,:,","))', advance='no') topo%vtors(:, j)
            write (iunit, '("],")')
         end do
         write (iunit, '(3x,"[",*(f25.15,:,","),"]",/)', advance='no') topo%vtors(:, topo%ntors)
         write (iunit, '("]")')
         write (iunit, '(3x,"],")')
      end if
      if (printTopo%vbond) then ! vbond(3,nbond)
         write (iunit, '(3x,''"vbond":'',"[")')
         do j = 1, neigh%nbond - 1
            write (iunit, '(3x,"[",*(f25.15,:,","))', advance='no') topo%vbond(:, j)
            write (iunit, '("],")')
         end do
         write (iunit, '(3x,"[",*(f25.15,:,","),"]",/)', advance='no') topo%vbond(:, neigh%nbond)
         write (iunit, '("]")')
         write (iunit, '(3x,"],")')
      end if
      if (printTopo%vangl) then ! vangl(2,nangl)
         write (iunit, '(3x,''"vangl":'',"[")')
         do j = 1, topo%nangl - 1
            write (iunit, '(3x,"[",*(f25.15,:,","))', advance='no') topo%vangl(:, j)
            write (iunit, '("],")')
         end do
         write (iunit, '(3x,"[",*(f25.15,:,","),"]",/)', advance='no') topo%vangl(:, topo%nangl)
         write (iunit, '("]")')
         write (iunit, '(3x,"],")')
      end if
      if (printTopo%hbbond) then ! hbbond: 3x(3,nhb) energies: 3x(1,nhb)
         write (iunit, '(3x,''"hbl":'',"[")') !> HBs loose
         if (nlist%nhb1 >= 1) then
            do j = 1, nlist%nhb1 - 1
               write (iunit, '(3x,"[",*(i7,:,","))', advance='no') nlist%hblist1(:, j)
               write (iunit, '("],")')
            end do
            write (iunit, '(3x,"[",*(i7,:,","),"]",/)', advance='no') nlist%hblist1(:, nlist%nhb1)
            write (iunit, '("]")')
            write (iunit, '(3x,"],")')
         else
            write (iunit, '(3x,"[",*(i7,:,""))', advance='no') 0
            write (iunit, '("]")')
            write (iunit, '(3x,"],")')
         end if

         write (iunit, '(3x,''"hbb":'',"[")') !> HBs bonded
         if (nlist%nhb2 >= 1) then
            do j = 1, nlist%nhb2 - 1
               write (iunit, '(3x,"[",*(i7,:,","))', advance='no') nlist%hblist2(:, j)
               write (iunit, '("],")')
            end do
            write (iunit, '(3x,"[",*(i7,:,","),"]",/)', advance='no') nlist%hblist2(:, nlist%nhb2)
            write (iunit, '("]")')
            write (iunit, '(3x,"],")')
         else
            write (iunit, '(3x,"[",*(i7,:,""))', advance='no') 0
            write (iunit, '("]")')
            write (iunit, '(3x,"],")')
         end if

         write (iunit, '(3x,''"xb":'',"[")') !> XBs
         if (nlist%nxb >= 1) then
            do j = 1, nlist%nxb - 1
               write (iunit, '(3x,"[",*(i7,:,","))', advance='no') nlist%hblist3(:, j)
               write (iunit, '("],")')
            end do
            write (iunit, '(3x,"[",*(i7,:,","),"]",/)', advance='no') nlist%hblist3(:, nlist%nxb)
            write (iunit, '("]")')
            write (iunit, '(3x,"],")')
         else
            write (iunit, '(3x,"[",*(i7,:,""))', advance='no') 0
            write (iunit, '("]")')
            write (iunit, '(3x,"],")')
         end if

         ! energies
         write (iunit, '(3x,''"hbl_e":'',"[")')
         do j = 1, nlist%nhb1 - 1
            write (iunit, '(3x,"[",*(f25.15,:,","))', advance='no') nlist%hbe1(j)
            write (iunit, '("],")')
         end do
         write (iunit, '(3x,"[",*(f25.15,:,","),"]",/)', advance='no') nlist%hbe1(nlist%nhb1)
         write (iunit, '("]")')
         write (iunit, '(3x,"],")')

         write (iunit, '(3x,''"hbb_e":'',"[")')
         do j = 1, nlist%nhb2 - 1
            write (iunit, '(3x,"[",*(f25.15,:,","))', advance='no') nlist%hbe2(j)
            write (iunit, '("],")')
         end do
         write (iunit, '(3x,"[",*(f25.15,:,","),"]",/)', advance='no') nlist%hbe2(nlist%nhb2)
         write (iunit, '("]")')
         write (iunit, '(3x,"],")')

         write (iunit, '(3x,''"xb_e":'',"[")')
         do j = 1, nlist%nxb - 1
            write (iunit, '(3x,"[",*(f25.15,:,","))', advance='no') nlist%hbe3(j)
            write (iunit, '("],")')
         end do
         write (iunit, '(3x,"[",*(f25.15,:,","),"]",/)', advance='no') nlist%hbe3(nlist%nxb)
         write (iunit, '("]")')
         write (iunit, '(3x,"],")')
      end if
      if (printTopo%eeq) then ! eeq(3,n)
         write (iunit, '(3x,''"eeq":'',"[")') !> EEQ charges
         do j = 1, size(nlist%q) - 1
            write (iunit, '(3x,"[",*(f25.15,:,","))', advance='no') nlist%q(j)
            write (iunit, '("],")')
         end do
         write (iunit, '(3x,"[",*(f25.15,:,","),"]",/)', advance='no') nlist%q(size(nlist%q))
         write (iunit, '("]")')
         write (iunit, '(3x,"],")')
      end if

      ! footer
      call get_command(length=l)
      allocate (character(len=l) :: cmdline)
      call get_command(cmdline)
      write (iunit, '(3x,''"program call":'',1x,''"'',a,''",'')') cmdline
      write (iunit, '(3x,''"method": "GFN-FF"'',",")')
      write (iunit, '(3x,a)') '"xtb version": "'//version//'"'
      write (iunit, '("}")')
      call close_file(iunit)

   end subroutine write_json_gfnff_lists
   
   !> wrapper for tblite-PTB JSON output
   subroutine main_ptb_json &
      (ijson, mol, wfx, calc, sccres, freqres)

      use xtb_type_molecule, only: TMolecule
      use xtb_type_wavefunction, only: TWavefunction
      use xtb_type_data, only: scc_results, freq_results
      use xtb_ptb_calculator, only: TPTBCalculator

      integer, intent(in) :: ijson 
      type(TMolecule), intent(in) :: mol
      type(TWavefunction), intent(in) :: wfx
      type(TPTBCalculator), intent(in) :: calc
      type(scc_results), intent(in) :: sccres
      type(freq_results), intent(in) :: freqres

#if WITH_TBLITE
      call tblite_ptb_json(ijson, mol, wfx, calc%bas, sccres, freqres)
#endif

   end subroutine main_ptb_json

#if WITH_TBLITE
   subroutine tblite_ptb_json &
      (ijson, mol, wfx, bas, sccres, freqres)
      use mctc_env, only: wp
   !! ========================================================================
      !  load class definitions
      use xtb_type_molecule, only: TMolecule
      use xtb_type_wavefunction, only: TWavefunction
      use xtb_type_data, only: scc_results, freq_results
      !> tblite-specific types
      use tblite_basis_type, only: basis_type
   !! ========================================================================
      !  global storage of options, parameters and basis set
      use xtb_setparam
      implicit none

   !! ========================================================================
      integer, intent(in) :: ijson ! file handle (usually json-file)
      !  molecule data
      type(TMolecule), intent(in) :: mol
      type(TWavefunction), intent(in) :: wfx
      type(basis_type), intent(in) :: bas
      type(scc_results), intent(in) :: sccres
      type(freq_results), intent(in) :: freqres
      logical :: alpha
      alpha = set%elprop == p_elprop_alpha

      call write_json_header(ijson)
      call write_json_scc_results(ijson, sccres)
      if (freqres%gtot > 0.0_wp) then
         call write_json_thermo(ijson, freqres)
      end if
      call write_json_charges(ijson, wfx)
      call write_json_ptb_shell_charges(ijson, mol, bas, wfx)
      call write_json_bondorder(ijson, mol, wfx)
      call write_json_dipole_moments(ijson, wfx)
      call write_json_quadrupole_moments(ijson, wfx)
      call write_json_ptb_wavefunction(ijson, wfx)
      if (freqres%n3true > 0) then
         call write_json_frequencies(ijson, freqres)
         call write_json_reduced_masses(ijson, freqres)
         call write_json_intensities(ijson, freqres, alpha)
      end if
      call write_json_ptb_footer(ijson)

   end subroutine tblite_ptb_json

   subroutine write_json_ptb_shell_charges(ijson, mol, bas, wfn)
      use xtb_type_wavefunction, only: TWavefunction
      use xtb_type_molecule, only: TMolecule
      use xtb_ptb_vdzp, only: max_shell
      use tblite_basis_type, only: basis_type
      integer, intent(in) :: ijson
      type(TMolecule), intent(in) :: mol
      type(basis_type), intent(in) :: bas
      type(TWavefunction), intent(in) :: wfn
      character(len=*), parameter :: jfmta = '(3x,''"'',a,''": ['')'
      character(len=*), parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
      integer :: iat, ish, ii
      write (ijson, jfmta) 'shell charges'
      do iat = 1, mol%n
         ii = bas%ish_at(iat)
         write (ijson, '(3x,a)', advance='no') "["
         do ish = 1, bas%nsh_at(iat) - 1
            write (ijson, '(f15.8,",")', advance="no") wfn%qsh(ii + ish)
         end do
         if (iat == mol%n) then
            write (ijson, '(f15.8,"]],",/)', advance="no") wfn%qsh(ii + bas%nsh_at(iat))
         else
            write (ijson, '(f15.8,"],",/)', advance="no") wfn%qsh(ii + bas%nsh_at(iat))
         end if
      end do
   end subroutine write_json_ptb_shell_charges
#endif

end module xtb_main_json
