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
subroutine main_json &
      (ijson,mol,wfx,xbas,sccres,freqres)
   use xtb_mctc_accuracy, only : wp

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
   type(TWavefunction),intent(in) :: wfx
   type(TBasisset),    intent(in) :: xbas
   type(scc_results),    intent(in) :: sccres
   type(freq_results),   intent(in) :: freqres

   call write_json_header(ijson)
   call write_json_scc_results(ijson,sccres)
   if (freqres%gtot.gt.0.0_wp) then
   call write_json_thermo(ijson,freqres)
   endif
   call write_json_charges(ijson,wfx)
   if (gfn_method.eq.2) then
   call write_json_dipole_moments(ijson,wfx)
   call write_json_quadrupole_moments(ijson,wfx)
   endif
   call write_json_wavefunction(ijson,wfx)
   if (freqres%n3true.gt.0) then
   call write_json_frequencies(ijson,freqres)
   call write_json_reduced_masses(ijson,freqres)
   call write_json_intensities(ijson,freqres)
   endif
   call write_json_footer(ijson)

end subroutine main_json

subroutine write_json_header(ijson)
   integer,intent(in) :: ijson
   write(ijson,'("{")')
end subroutine write_json_header

subroutine write_json_footer(ijson)
   use xtb_setparam
   include 'xtb_version.fh'
   integer,intent(in) :: ijson
   character(len=:),allocatable :: cmdline
   integer :: l
   call get_command(length=l)
   allocate( character(len=l) :: cmdline )
   call get_command(cmdline)
   write(ijson,'(3x,''"program call":'',1x,''"'',a,''",'')') cmdline
   write(ijson,'(3x,''"method": "GFN'',i0,''-xTB",'')') gfn_method
   write(ijson,'(3x,a)') '"xtb version": "'//version//'"'
   write(ijson,'("}")')
end subroutine write_json_footer

subroutine write_json_scc_results(ijson,sccres)
   use xtb_type_data
   integer,intent(in) :: ijson
   type(scc_results),intent(in) :: sccres
   character(len=*),parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
   write(ijson,jfmtf) 'total energy',sccres%e_total
   write(ijson,jfmtf) 'HOMO-LUMO gap/eV',sccres%hl_gap
   write(ijson,jfmtf) 'electronic energy',sccres%e_elec
   write(ijson,'(3x,''"'',a,''":'',1x,"[",2(f15.8,","),f15.8,"],")') &
      'dipole',sccres%dipole
   !write(ijson,jfmtf) 'classical repulsion energy',sccres%e_rep
   !write(ijson,jfmtf) 'isotropic electrostatic energy',sccres%e_es
   !write(ijson,jfmtf) 'anisotropic electrostatic energy',sccres%e_aes
   !write(ijson,jfmtf) 'anisotropic XC energy',sccres%e_axc
   !write(ijson,jfmtf) 'classical halogen bound energy',sccres%e_xb
   !write(ijson,jfmtf) 'Generalized Born free energy',sccres%g_born
   !write(ijson,jfmtf) 'SASA free energy',sccres%g_born
   !write(ijson,jfmtf) 'Hydrogen bound free energy',sccres%g_born
end subroutine write_json_scc_results

subroutine write_json_charges(ijson,wfn)
   use xtb_type_wavefunction
   integer,intent(in) :: ijson
   type(TWavefunction),intent(in) :: wfn
   character(len=*),parameter :: jfmta = '(3x,''"'',a,''": ['')'
   character(len=*),parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
   integer :: i
   write(ijson,jfmta) 'partial charges'
   write(ijson,'(3x,f15.8,",")') (wfn%q(i),i=1,wfn%n-1)
   write(ijson,'(3x,f15.8,"],")')  wfn%q(wfn%n)
end subroutine write_json_charges

subroutine write_json_dipole_moments(ijson,wfn)
   use xtb_type_wavefunction
   integer,intent(in) :: ijson
   type(TWavefunction),intent(in) :: wfn
   character(len=*),parameter :: jfmta = '(3x,''"'',a,''": ['')'
   character(len=*),parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
   integer :: i,j
   write(ijson,jfmta) 'atomic dipole moments'
   do i = 1, wfn%n-1
   write(ijson,'(3x,"[",2(f15.8,","),f15.8,"],")') (wfn%dipm(j,i),j=1,3)
   enddo
   write(ijson,'(3x,"[",2(f15.8,","),f15.8,"]],")') (wfn%dipm(j,wfn%n),j=1,3)
end subroutine write_json_dipole_moments

subroutine write_json_quadrupole_moments(ijson,wfn)
   use xtb_type_wavefunction
   integer,intent(in) :: ijson
   type(TWavefunction),intent(in) :: wfn
   character(len=*),parameter :: jfmta = '(3x,''"'',a,''": ['')'
   character(len=*),parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
   integer :: i,j
   write(ijson,jfmta) 'atomic quadrupole moments'
   do i = 1, wfn%n-1
   write(ijson,'(3x,"[",5(f15.8,","),f15.8,"],")') (wfn%qp(j,i),j=1,6)
   enddo
   write(ijson,'(3x,"[",5(f15.8,","),f15.8,"]],")') (wfn%qp(j,wfn%n),j=1,6)
end subroutine write_json_quadrupole_moments

subroutine write_json_wavefunction(ijson,wfn)
   use xtb_type_wavefunction
   integer,intent(in) :: ijson
   type(TWavefunction),intent(in) :: wfn
   character(len=*),parameter :: jfmta = '(3x,''"'',a,''": ['')'
   character(len=*),parameter :: jfmti = '(3x,''"'',a,''":'',1x,i0,",")'
   character(len=*),parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
   integer :: i
   write(ijson,jfmti) 'number of molecular orbitals',wfn%nao
   write(ijson,jfmti) 'number of electrons',wfn%nel
   write(ijson,jfmti) 'number of unpaired electrons',wfn%nopen
   write(ijson,jfmta) 'orbital energies/eV'
   write(ijson,'(3x,f15.8,",")') (wfn%emo(i),i=1,wfn%nao-1)
   write(ijson,'(3x,f15.8,"],")')  wfn%emo(wfn%nao)
   write(ijson,jfmta) 'fractional occupation'
   write(ijson,'(3x,f15.8,",")') (wfn%focc(i),i=1,wfn%nao-1)
   write(ijson,'(3x,f15.8,"],")')  wfn%focc(wfn%nao)
end subroutine write_json_wavefunction

subroutine write_json_thermo(ijson,freqres)
   use xtb_type_data
   integer,intent(in) :: ijson
   type(freq_results),intent(in) :: freqres
   character(len=*),parameter :: jfmtf = '(3x,''"'',a,''":'',1x,f20.8,",")'
   write(ijson,jfmtf) 'total enthalpy',freqres%htot
   write(ijson,jfmtf) 'total free energy',freqres%gtot
end subroutine write_json_thermo

subroutine write_json_frequencies(ijson,freqres)
   use xtb_type_data
   integer,intent(in) :: ijson
   type(freq_results),intent(in) :: freqres
   character(len=*),parameter :: jfmta = '(3x,''"'',a,''": ['')'
   write(ijson,jfmta) 'vibrational frequencies/rcm'
   write(ijson,'(3x,f15.8,",")') (freqres%freq(i),i=1,freqres%n3true-1)
   write(ijson,'(3x,f15.8,"],")') freqres%freq(freqres%n3true)
end subroutine write_json_frequencies

subroutine write_json_intensities(ijson,freqres)
   use xtb_type_data
   integer,intent(in) :: ijson
   type(freq_results),intent(in) :: freqres
   character(len=*),parameter :: jfmta = '(3x,''"'',a,''": ['')'
   write(ijson,jfmta) 'IR intensities/amu'
   write(ijson,'(3x,f15.8,",")') (freqres%dipt(i),i=1,freqres%n3true-1)
   write(ijson,'(3x,f15.8,"],")') freqres%dipt(freqres%n3true)
end subroutine write_json_intensities

subroutine write_json_reduced_masses(ijson,freqres)
   use xtb_type_data
   integer,intent(in) :: ijson
   type(freq_results),intent(in) :: freqres
   character(len=*),parameter :: jfmta = '(3x,''"'',a,''": ['')'
   write(ijson,jfmta) 'reduced masses'
   write(ijson,'(3x,f15.8,",")') (freqres%rmass(i),i=1,freqres%n3true-1)
   write(ijson,'(3x,f15.8,"],")') freqres%rmass(freqres%n3true)
end subroutine write_json_reduced_masses

subroutine write_json_gfnff_lists(n, topo, nlist, printTopo)
  use xtb_gfnff_topology, only : TGFFTopology
  use xtb_gfnff_neighbourlist, only : TGFFNeighbourList
  use xtb_gfnff_topology, only : TPrintTopo
  include 'xtb_version.fh'
  !> gfnff topology lists
  type(TGFFTopology), intent(in) :: topo
  !> gfnff neighbourlist
  type(TGFFNeighbourList), intent(in) :: nlist
  !> topology printout booleans
  type(TPrintTopo), intent(in) :: printTopo
  character(len=:),allocatable :: cmdline
  integer :: iunit, j, n

  call open_file(iunit,'gfnff_lists.json','w')
  ! header
  write(iunit,'("{")')
  ! lists printout
  if(printTopo%nb)then ! nb(20,n)
    write(iunit,'(3x,''"nb":'',"[")')
    do j=1, n-1
      write(iunit,'(3x,"[",*(i7,:,","))',advance='no') topo%nb(:,j)
      write(iunit,'("],")')
    enddo
    write(iunit,'(3x,"[",*(i7,:,","),"]",/)',advance='no') topo%nb(:,n)
    write(iunit,'("]")')
    write(iunit,'(3x,"],")')
  endif
  if(printTopo%bpair)then ! bpair(n*(n+1)/2) packed symmetric matrix
    write(iunit,'(3x,''"bpair":'',"[")')
    write(iunit,'(3x,*(i7,:,","))',advance='no') topo%bpair
    write(iunit,'(3x,"],")')
  endif
  if(printTopo%alist)then ! alist(3,nangl)
    write(iunit,'(3x,''"alist":'',"[")')
    do j=1, topo%nangl-1
      write(iunit,'(3x,"[",*(i8,:,","))',advance='no') topo%alist(:,j)
      write(iunit,'("],")')
    enddo
    write(iunit,'(3x,"[",*(i8,:,","),"]",/)',advance='no') topo%alist(:,topo%nangl)
    write(iunit,'("]")')
    write(iunit,'(3x,"],")')
  endif
  if(printTopo%blist)then ! blist(2,nbond)
    write(iunit,'(3x,''"blist":'',"[")')
    do j=1, topo%nbond-1
      write(iunit,'(3x,"[",*(i8,:,","))',advance='no') topo%blist(:,j)
      write(iunit,'("],")')
    enddo
    write(iunit,'(3x,"[",*(i8,:,","),"]",/)',advance='no') topo%blist(:,topo%nbond)
    write(iunit,'("]")')
    write(iunit,'(3x,"],")')
  endif
  if(printTopo%tlist)then ! tlist(5,ntors)
    write(iunit,'(3x,''"tlist":'',"[")')
    do j=1, topo%ntors-1
      write(iunit,'(3x,"[",*(i8,:,","))',advance='no') topo%tlist(:,j)
      write(iunit,'("],")')
    enddo
    write(iunit,'(3x,"[",*(i8,:,","),"]",/)',advance='no') topo%tlist(:,topo%ntors)
    write(iunit,'("]")')
    write(iunit,'(3x,"],")')
  endif
  if(printTopo%vtors)then ! vtors(2,ntors)
    write(iunit,'(3x,''"vtors":'',"[")')
    do j=1, topo%ntors-1
      write(iunit,'(3x,"[",*(f25.15,:,","))',advance='no') topo%vtors(:,j)
      write(iunit,'("],")')
    enddo
    write(iunit,'(3x,"[",*(f25.15,:,","),"]",/)',advance='no') topo%vtors(:,topo%ntors)
    write(iunit,'("]")')
    write(iunit,'(3x,"],")')
  endif
  if(printTopo%vbond)then ! vbond(3,nbond)
    write(iunit,'(3x,''"vbond":'',"[")')
    do j=1, topo%nbond-1
      write(iunit,'(3x,"[",*(f25.15,:,","))',advance='no') topo%vbond(:,j)
      write(iunit,'("],")')
    enddo
    write(iunit,'(3x,"[",*(f25.15,:,","),"]",/)',advance='no') topo%vbond(:,topo%nbond)
    write(iunit,'("]")')
    write(iunit,'(3x,"],")')
  endif
  if(printTopo%vangl)then ! vangl(2,nangl)
    write(iunit,'(3x,''"vangl":'',"[")')
    do j=1, topo%nangl-1
      write(iunit,'(3x,"[",*(f25.15,:,","))',advance='no') topo%vangl(:,j)
      write(iunit,'("],")')
    enddo
    write(iunit,'(3x,"[",*(f25.15,:,","),"]",/)',advance='no') topo%vangl(:,topo%nangl)
    write(iunit,'("]")')
    write(iunit,'(3x,"],")')
  endif
  if(printTopo%hbbond)then ! hbbond
   write(iunit,'(3x,''"hbl":'',"[")') !> HBs loose
   if (nlist%nhb1.ge.1) then
    do j=1, nlist%nhb1-1
      write(iunit,'(3x,"[",*(i7,:,","))',advance='no') nlist%hblist1(:,j)
      write(iunit,'("],")')
    enddo
    write(iunit,'(3x,"[",*(i7,:,","),"]",/)',advance='no') nlist%hblist1(:,nlist%nhb1)
    write(iunit,'("]")')
    write(iunit,'(3x,"],")')
   else
     write(iunit,'(3x,"[",*(i7,:,""))',advance='no') 0
     write(iunit,'("]")')
     write(iunit,'(3x,"],")')
   end if

   write(iunit,'(3x,''"hbb":'',"[")') !> HBs bonded
   if (nlist%nhb2.ge.1) then
    do j=1, nlist%nhb2-1
      write(iunit,'(3x,"[",*(i7,:,","))',advance='no') nlist%hblist2(:,j)
      write(iunit,'("],")')
    enddo
    write(iunit,'(3x,"[",*(i7,:,","),"]",/)',advance='no') nlist%hblist2(:,nlist%nhb2)
    write(iunit,'("]")')
    write(iunit,'(3x,"],")')
   else
     write(iunit,'(3x,"[",*(i7,:,""))',advance='no') 0
     write(iunit,'("]")')
     write(iunit,'(3x,"],")')
   end if

   write(iunit,'(3x,''"xb":'',"[")') !> XBs
   if (nlist%nxb.ge.1) then
    do j=1, nlist%nxb-1
      write(iunit,'(3x,"[",*(i7,:,","))',advance='no') nlist%hblist3(:,j)
      write(iunit,'("],")')
    enddo
    write(iunit,'(3x,"[",*(i7,:,","),"]",/)',advance='no') nlist%hblist3(:,nlist%nxb)
    write(iunit,'("]")')
    write(iunit,'(3x,"],")')
   else
     write(iunit,'(3x,"[",*(i7,:,""))',advance='no') 0
     write(iunit,'("]")')
     write(iunit,'(3x,"],")')
   end if

   ! energies
   write(iunit,'(3x,''"hbl_e":'',"[")') 
   do j=1, nlist%nhb1-1
     write(iunit,'(3x,"[",*(f25.15,:,","))',advance='no') nlist%hbe1(j)
     write(iunit,'("],")')
   enddo
   write(iunit,'(3x,"[",*(f25.15,:,","),"]",/)',advance='no') nlist%hbe1(nlist%nhb1)
   write(iunit,'("]")')
   write(iunit,'(3x,"],")')

   write(iunit,'(3x,''"hbb_e":'',"[")') 
   do j=1, nlist%nhb2-1
     write(iunit,'(3x,"[",*(f25.15,:,","))',advance='no') nlist%hbe2(j)
     write(iunit,'("],")')
   enddo
   write(iunit,'(3x,"[",*(f25.15,:,","),"]",/)',advance='no') nlist%hbe2(nlist%nhb2)
   write(iunit,'("]")')
   write(iunit,'(3x,"],")')

   write(iunit,'(3x,''"xb_e":'',"[")') 
   do j=1, nlist%nxb-1
     write(iunit,'(3x,"[",*(f25.15,:,","))',advance='no') nlist%hbe3(j)
     write(iunit,'("],")')
   enddo
   write(iunit,'(3x,"[",*(f25.15,:,","),"]",/)',advance='no') nlist%hbe3(nlist%nxb)
   write(iunit,'("]")')
   write(iunit,'(3x,"],")')
  endif
  if(printTopo%eeq)then ! eeq charges
   write(iunit,'(3x,''"eeq":'',"[")') !> EEQ charges
   do j=1, size(nlist%q)-1
     write(iunit,'(3x,"[",*(f25.15,:,","))',advance='no') nlist%q(j)
     write(iunit,'("],")')
   enddo
   write(iunit,'(3x,"[",*(f25.15,:,","),"]",/)',advance='no')nlist%q(size(nlist%q))
   write(iunit,'("]")')
   write(iunit,'(3x,"],")')
  endif
  ! footer
  call get_command(length=l)
  allocate( character(len=l) :: cmdline )
  call get_command(cmdline)
  write(iunit,'(3x,''"program call":'',1x,''"'',a,''",'')') cmdline
  write(iunit,'(3x,''"method": "GFN-FF"'',",")')
  write(iunit,'(3x,a)') '"xtb version": "'//version//'"'
  write(iunit,'("}")')
  call close_file(iunit)

end subroutine write_json_gfnff_lists
