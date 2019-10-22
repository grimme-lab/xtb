! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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

module property_output

contains

subroutine write_energy(iunit,sccres,frqres,hess)
   use iso_fortran_env, wp => real64
   use tbdef_data
   implicit none
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
   logical, intent(in) :: hess
   type(scc_results), intent(in) :: sccres
   type(freq_results),intent(in) :: frqres
   character(len=*),parameter :: outfmt = '(10x,"|",1x,a,f24.12,1x,a,1x,"|")'
   write(iunit,'(a)')
   write(iunit,'(11x,49("-"))')
   if (hess) then
      write(iunit,outfmt) "TOTAL ENERGY      ", frqres%etot,            "Eh  "
      write(iunit,outfmt) "TOTAL ENTHALPY    ", frqres%etot+frqres%htot,"Eh  "
      write(iunit,outfmt) "TOTAL FREE ENERGY ", frqres%etot+frqres%gtot,"Eh  "
      write(iunit,outfmt) "GRADIENT NORM     ", frqres%gnorm,           "Eh/α"
   else
      write(iunit,outfmt) "TOTAL ENERGY      ", sccres%e_total,"Eh  "
      write(iunit,outfmt) "GRADIENT NORM     ", sccres%gnorm,  "Eh/α"
   endif
   write(iunit,outfmt) "HOMO-LUMO GAP     ", sccres%hl_gap, "eV  "
   write(iunit,'(11x,49("-"))')
end subroutine write_energy

subroutine main_property &
      (iunit,mol,wfx,basis,xpar,res,acc)
   use iso_fortran_env, wp => real64

   use mctc_econv

!! ========================================================================
!  load class definitions
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_data
   use tbdef_param

!! ========================================================================
!  global storage of options, parameters and basis set
   use setparam
   use aoparam
   use gbobc

!! ------------------------------------------------------------------------
   use scc_core, only : wiberg
   use aespot
   use dtrafo

   implicit none

!! ========================================================================
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
!  molecule data
   type(tb_molecule), intent(in) :: mol
   real(wp),intent(in) :: acc      ! accuracy of integral calculation
   type(tb_wavefunction),intent(inout) :: wfx
   type(tb_basisset),    intent(in) :: basis
   type(scc_parameter),  intent(in) :: xpar
   type(scc_results),    intent(in) :: res

   real(wp),allocatable :: S(:,:)     ! overlap integrals
   real(wp),allocatable :: dpint(:,:) ! dipole integrals
   real(wp),allocatable :: qpint(:,:) ! quadrupole integrals
   real(wp),allocatable :: C(:,:)     ! molecular orbitals
   real(wp),allocatable :: emo(:)     ! orbital energies
   real(wp),allocatable :: focc(:)    ! fractional occupation numbers
   integer  :: ifile
   integer  :: ndim,ndp,nqp
   real(wp) :: dip,dipol(3)
   real(wp) :: intcut,neglect

   type(tb_solvent) :: gbsa

!  primitive cut-off
   intcut=25.0_wp-10.0*log10(acc)
   intcut=max(20.0_wp,intcut)
!  integral neglect threshold
   neglect =10.0d-9*acc
   ndim = basis%nao*(basis%nao+1)/2
   allocate(S(basis%nao,basis%nao), dpint(3,ndim), qpint(6,ndim), source = 0.0_wp )
   call sdqint(mol%n,mol%at,basis%nbf,basis%nao,mol%xyz,neglect,ndp,nqp,intcut, &
      &        basis%caoshell,basis%saoshell,basis%nprim,basis%primcount, &
      &        basis%alp,basis%cont,S,dpint,qpint)

!! orbital energies and occupation
    if (pr_eig) then
       write(iunit,'(/,4x,"*",1x,a)') "Orbital Energies and Occupations"
       call print_orbital_eigenvalues(iunit,wfx,11)
    endif

!! Mulliken and CM5 charges
   if (pr_mulliken.and.gfn_method.eq.1) then
      call open_file(ifile,'charges','w')
      call print_mulliken(iunit,ifile,mol%n,mol%at,mol%xyz,mol%z,basis%nao,S,wfx%P,basis%aoat2,basis%lao2)
      call close_file(ifile)
   else if (pr_charges) then
      call open_file(ifile,'charges','w')
      call print_charges(ifile,mol%n,wfx%q)
      call close_file(ifile)
   endif

   ! GBSA information
   if (lgbsa.and.pr_gbsa) then
      call new_gbsa(gbsa,mol%n,mol%at)
      call update_nnlist_gbsa(gbsa,mol%xyz,.false.)
      call compute_brad_sasa(gbsa,mol%xyz)
      call print_gbsa_info(iunit,gbsa)
   endif

!! D4 molecular dispersion printout
   if ((newdisp.and.gfn_method.eq.2).and.pr_mulliken) &
   call print_molpol(iunit,mol%n,mol%at,mol%xyz,wfx%q,xpar%wf,xpar%g_a,xpar%g_c)
   if (gfn_method.eq.0.and.pr_mulliken) &
   call print_molpol(iunit,mol%n,mol%at,mol%xyz,wfx%q,xpar%wf,xpar%g_a,xpar%g_c)

!! Spin population
   if (pr_spin_population .and. wfx%nopen.ne.0) &
   call print_spin_population(iunit,mol%n,mol%at,basis%nao,wfx%focca,wfx%foccb,S,wfx%C, &
   &                          basis%aoat2,basis%lao2)

   if (pr_fod_pop) then
      call open_file(ifile,'fod','w')
      call print_fod_population(iunit,ifile,mol%n,mol%at,basis%nao,S,wfx%C,etemp,wfx%emo, &
                                wfx%ihomoa,wfx%ihomob,basis%aoat2,basis%lao2)
      call close_file(ifile)
   endif


!! wiberg bond orders
   if (pr_wiberg.and.gfn_method.eq.0) &
   call wiberg(mol%n,basis%nao,mol%at,mol%xyz,wfx%P,S,wfx%wbo,.false.,.false.,basis%fila2)
   if (pr_wiberg) &
   call print_wiberg(iunit,mol%n,mol%at,wfx%wbo,0.1_wp)

   if (pr_wbofrag) &
   call print_wbo_fragment(iunit,mol%n,mol%at,wfx%wbo,0.1_wp)

!! molden file
   if (pr_molden_input) then
      allocate(C(basis%nbf,basis%nao),focc(basis%nao),emo(basis%nao), source = 0.0_wp)
      if (basis%nbf.eq.basis%nao) then
         C = wfx%C
      else
         call sao2cao(basis%nao,wfx%C,basis%nbf,C,basis)
      endif
      emo  = wfx%emo * evtoau
      focc = wfx%focca + wfx%foccb
      call printmold(mol%n,basis%nao,basis%nbf,mol%xyz,mol%at,C,emo,focc,2.0_wp,basis)
      write(iunit,'(/,"MOs/occ written to file <molden.input>",/)')
      deallocate(C,focc,emo)
   endif

   if (pr_gbw) &
   call wrgbw(mol%n,mol%at,mol%xyz,mol%z,basis,wfx)

   if (pr_tmbas .or. pr_tmmos) then
      call open_file(ifile,'basis','w')
      call write_tm_basis(ifile,mol%n,mol%at,basis,wfx)
      close(ifile)
   endif

   if (pr_tmmos) then
      call open_file(ifile,'mos','w')
      call write_tm_mos(ifile,mol%n,mol%at,basis,wfx)
      close(ifile)
   endif

!! multipole moment prinout
   if (pr_dipole) then
      if (gfn_method.gt.1) then
         ! print overall multipole moment
         call molmom(iunit,mol%n,mol%xyz,wfx%q,wfx%dipm,wfx%qp,dip,dipol)
         write(iunit,'(a)')
      else
         call print_dipole(iunit,mol%n,mol%at,mol%xyz,mol%z,wfx%nao,wfx%P,dpint)
      endif
   endif

end subroutine main_property

subroutine main_cube &
      (lverbose,mol,wfx,basis,xpar,res)
   use iso_fortran_env, wp => real64, istdout => output_unit

   use mctc_econv

!! ========================================================================
!  load class definitions
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_data
   use tbdef_param

!! ========================================================================
!  global storage of options, parameters and basis set
   use setparam
   use aoparam

!! ------------------------------------------------------------------------
   use aespot
   use scc_core
   use esp
   use stm
   use dtrafo

   implicit none

!! ========================================================================
   logical, intent(in) :: lverbose
!  molecule data
   type(tb_molecule), intent(in) :: mol
   type(tb_wavefunction),intent(in) :: wfx
   type(tb_basisset),    intent(in) :: basis
   type(scc_parameter),  intent(in) :: xpar
   type(scc_results),    intent(in) :: res

   real(wp),allocatable :: C(:,:)     ! molecular orbitals
   real(wp),allocatable :: emo(:)     ! orbital energies
   real(wp),allocatable :: focc(:)    ! fractional occupation numbers
   real(wp),allocatable :: focca(:)   ! fractional occupation numbers (alpha)
   real(wp),allocatable :: foccb(:)   ! fractional occupation numbers (beta)
   integer  :: ndim,ndp,nqp
   real(wp) :: dip,dipol(3)
   real(wp) :: acc,intcut,neglect
   real(wp) :: efa,efb,ga,gb,nfoda,nfodb

!! ------------------------------------------------------------------------
!  FOD
   if (pr_fod) then
      allocate( C(basis%nbf,basis%nao), focca(basis%nao), foccb(basis%nao), focc(basis%nao), emo(basis%nao), &
                source = 0.0_wp )
      if(wfx%ihomoa+1.le.basis%nao) &
         call fermismear(.false.,basis%nao,wfx%ihomoa,etemp,wfx%emo,focca,nfoda,efa,ga)
      if(wfx%ihomob+1.le.basis%nao) &
         call fermismear(.false.,basis%nao,wfx%ihomob,etemp,wfx%emo,foccb,nfodb,efb,gb)
      emo = wfx%emo * evtoau
      call fodenmak(.true.,basis%nao,emo,focca,efa)
      call fodenmak(.true.,basis%nao,emo,foccb,efb)
      focc = focca+foccb
      if(basis%nbf.eq.basis%nao) then
         C = wfx%C
      else
         call sao2cao(basis%nao,wfx%C,basis%nbf,C,basis)
      endif
      if (lverbose) &
      write(istdout,'(/,"FOD written to file: ''fod.cub''",/)')
      call cube(mol%n,basis%nao,basis%nbf,mol%xyz,mol%at,C,emo,focc,'fod.cub',basis)
      deallocate(C, focca, foccb, focc, emo)
   endif

!! ------------------------------------------------------------------------
!  print spin density to cube file
   if (pr_spin_density.and.wfx%nopen.ne.0) then
      allocate( C(basis%nbf,basis%nao), focc(basis%nao), emo(basis%nao), source = 0.0_wp )
      if(basis%nbf.eq.basis%nao) then
         C = wfx%C
      else
         call sao2cao(basis%nao,wfx%C,basis%nbf,C,basis)
      endif
      if (lverbose) &
      write(istdout,'(/,"(R)spin-density written to file: ''spindensity.cub''",/)')
      emo = wfx%emo * evtoau
      focc = wfx%focca - wfx%foccb
      call cube(mol%n,basis%nao,basis%nbf,mol%xyz,mol%at,C,emo,focc,'spindensity.cub',basis)
      deallocate(C, focc, emo)
   endif

!! ------------------------------------------------------------------------
!  print density to cube file
   if (pr_density) then
      allocate( C(basis%nbf,basis%nao), emo(basis%nao), source = 0.0_wp )
      if(basis%nbf.eq.basis%nao) then
         C = wfx%C
      else
         call sao2cao(basis%nao,wfx%C,basis%nbf,C,basis)
      endif
      if (lverbose) &
      write(istdout,'(/,"density written to file: ''density.cub''",/)')
      emo = wfx%emo * evtoau
      call cube(mol%n,basis%nao,basis%nbf,mol%xyz,mol%at,C,emo,wfx%focc,'density.cub',basis)
      deallocate(C, emo)
   endif

!! ------------------------------------------------------------------------
!  make an ESP plot
   if (pr_esp) then
      allocate( C(basis%nbf,basis%nao), source = 0.0_wp )
      if(basis%nbf.eq.basis%nao) then
         C = wfx%C
      else
         call sao2cao(basis%nao,wfx%C,basis%nbf,C,basis)
      endif
      call espplot(mol%n,basis%nao,basis%nbf,mol%at,mol%xyz,mol%z,wfx%focc,C,basis)
      deallocate(C)
   endif

!! ------------------------------------------------------------------------
!  make a STM image
   if (pr_stm) then
      allocate( C(basis%nbf,basis%nao), focc(basis%nao), source = 0.0_wp )
      if(basis%nbf.eq.basis%nao) then
         C = wfx%C
      else
         call sao2cao(basis%nao,wfx%C,basis%nbf,C,basis)
      endif
      if(wfx%ihomoa+1.le.wfx%nao) &
         call fermismear(.false.,basis%nao,wfx%ihomoa,etemp,wfx%emo,focc,nfoda,efa,ga)
      if(wfx%ihomob+1.le.wfx%nao) &
         call fermismear(.false.,basis%nao,wfx%ihomob,etemp,wfx%emo,focc,nfodb,efb,gb)
      call stmpic(mol%n,basis%nao,basis%nbf,mol%at,mol%xyz,C,0.5_wp*(efa+efb),wfx%emo,basis)
      deallocate(C, focc)
   endif


end subroutine main_cube

subroutine main_freq &
      (iunit,mol,wfx,res)
   use iso_fortran_env, wp => real64

   use mctc_econv

!! ========================================================================
!  load class definitions
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_data
   use tbdef_param

!! ========================================================================
!  global storage of options, parameters and basis set
   use setparam
   use aoparam

!! ------------------------------------------------------------------------
   use hessian
   use ncoord

   implicit none

!! ========================================================================
   integer, intent(in) :: iunit
!  molecule data
   type(tb_molecule), intent(in) :: mol
   type(tb_wavefunction),intent(in) :: wfx
   type(freq_results),   intent(inout) :: res

   integer  :: ifile
   integer  :: i,ii,j,jj,k,l
   character(len=:),allocatable  :: hname
   integer, allocatable :: bond(:,:)
   integer, allocatable :: molvec(:)
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: xyz0(:,:)
   real(wp),allocatable :: h(:,:)
   real(wp) :: etot,h298,dum
   integer  :: lowmode

   allocate( molvec(mol%n), bond(mol%n,mol%n), source = 0 )
   allocate( xyz0(3,mol%n), h(3*mol%n,3*mol%n), cn(mol%n), source = 0.0_wp )

   if(res%linear)then
      write(iunit,'(1x,a)') 'vibrational frequencies (cm-1)'
   else
      write(iunit,'(1x,a)') 'projected vibrational frequencies (cm-1)'
   endif
   call PREIGF(iunit,res%freq,res%n3true)

   write(iunit,'(1x,a)') 'reduced masses (amu)'
   write(iunit,'(8(i4,'':'',f6.2))') (i,res%rmass(i),i=1,res%n3)
   write(iunit,'(1x,a)') 'IR intensities (amu)'
   write(iunit,'(8(i4,'':'',f6.2))') (i,res%dipt(i),i=1,res%n3)
   write(iunit,'(1x,a)') 'Raman intensities (amu)'
   write(iunit,'(8(i4,'':'',f6.2))') (i,res%polt(i),i=1,res%n3)

   call open_file(ifile,'vibspectrum','w')
   call write_tm_vibspectrum(ifile,res%n3,res%freq,res%dipt)
   call close_file(ifile)

   write(iunit,'(1x,a)') 'output can be read by thermo (or use thermo option).'
   write(iunit,'(1x,a)') 'writing <g98.out> molden fake output.'
   write(iunit,'(1x,a)') &
      & 'recommended (thermochemical) frequency scaling factor: 1.0'
   call g98fake2('g98.out',mol%n,mol%at,mol%xyz,res%freq,res%rmass,res%dipt,res%hess)

   call generic_header(iunit,"Thermodynamic Functions",49,10)
   call print_thermo(iunit,mol%n,res%n3true,mol%at,mol%xyz,res%freq,res%etot,res%htot,res%gtot,res%nimag,.true.)
   res%pg = trim(pgroup)
   res%temp = thermotemp(nthermo)
   if (enso_mode) then
      call open_file(ifile,"xtb_enso.json",'w')
      if (ifile .ne. -1) then
         call enso_printout(ifile,res)
         call close_file(ifile)
      endif
   endif

   ! distort along imags if present
   xyz0 = mol%xyz
   call distort(mol%n,mol%at,xyz0,res%freq,res%hess)

   if(pr_modef .and. (mol%n.gt.3)) then

      ! do analysis and write mode following file
      call wrmodef(0,mol%n,mol%at,mol%xyz,wfx%wbo,res%rmass,res%freq,res%hess,h,mode_vthr,res%linear)

      ! localize the modes
      if(mode_vthr.gt.1.d-6)then
         ! determine molecular fragments
         call ncoord_erf(mol%n,mol%at,mol%xyz,cn)
         call cutcov(mol%n,mol%at,mol%xyz,cn,wfx%wbo,bond)
         call mrec(i,mol%xyz,cn,bond,mol%n,mol%at,molvec)
         call locmode(mol%n,res%n3,mol%at,mol%xyz,mode_vthr,res%freq,res%rmass,res%hess, &
                      i,molvec)
         call PREIGF0(iunit,res%freq,res%n3true)
         write(iunit,'("written to xtb_localmodes and g98l.out")')
         call wrmodef(1,mol%n,mol%at,mol%xyz,wfx%wbo,res%rmass,res%freq,res%hess, &
                      h,mode_vthr+200.0_wp,res%linear)
      endif

      call open_file(ifile,'.tmpxtbmodef','w')
      write(ifile,*) res%lowmode,res%lowmode
      write(ifile,*) res%etot    ! energy for comparison
      call close_file(ifile)

   endif

end subroutine main_freq

subroutine print_charges(ifile,n,q)
   use iso_fortran_env, wp => real64
   implicit none
   integer, intent(in)  :: ifile
   integer, intent(in)  :: n
   real(wp),intent(in)  :: q(n)
   integer :: i
   if (ifile.ne.-1) then
      do i = 1, n
         write(ifile,'(f14.8)') q(i)
      enddo
   endif
end subroutine print_charges

subroutine print_mulliken(iunit,ifile,n,at,xyz,z,nao,S,P,aoat2,lao2)
   use iso_fortran_env, wp => real64
   use scc_core, only : mpop
   implicit none
   integer, intent(in)  :: iunit
   integer, intent(in)  :: ifile
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: z(n)
   integer, intent(in)  :: nao
   real(wp),intent(in)  :: S(nao,nao)
   real(wp),intent(in)  :: P(nao,nao)
   integer, intent(in)  :: aoat2(nao)
   integer, intent(in)  :: lao2(nao)
   real(wp),allocatable :: q(:)       ! Mulliken partial charges
   real(wp),allocatable :: qlmom(:,:) ! population per shell
   real(wp),allocatable :: cm5(:)     ! CM5 partial charges
   character(len=2),external :: asym
   integer :: i

   allocate( cm5(n), q(n), qlmom(3,n), source = 0.0_wp )
   call mpop(n,nao,aoat2,lao2,S,P,q,qlmom)
   q = q - z
   call docm5(n,at,.false.,xyz,q,cm5)
   write(iunit,'(a)')
   write(iunit,'(2x,"Mulliken/CM5 charges        n(s)   n(p)   n(d)")')
   do i=1,n
      write(iunit,'(i6,a3,2f9.5,1x,4f7.3)') &
         i,asym(at(i)),q(i),cm5(i),qlmom(1,i),qlmom(2,i),qlmom(3,i)
   enddo
   if (ifile.ne.-1) then
      do i = 1, n
         write(ifile,'(f14.8)') q(i)
      enddo
   endif
end subroutine print_mulliken

subroutine print_wiberg(iunit,n,at,wbo,thr)
   use iso_fortran_env, wp => real64
   use scc_core, only : wibsort
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: wbo(n,n)
   real(wp),intent(in) :: thr

   real(wp),allocatable :: wbr(:,:)
   integer, allocatable :: imem(:)
   character(len=2),external :: asym
   integer  :: i,j,ibmax
   real(wp) :: xsum

   allocate( wbr(n,n), source = wbo )
   allocate( imem(n),  source = 0 )

   write(iunit,'(a)')
   write(iunit,'("Wiberg/Mayer (AO) data.")')
   write(iunit,'("largest (>",f4.2,") Wiberg bond orders for each atom")') thr
   write(iunit,'("          total WBO             WBO to atom ...")')
   do i=1,n
      do j=1,n
         imem(j)=j
      enddo
      call wibsort(n,i,imem,wbr)
      ibmax=0
      xsum =0.0_wp
      do j=1,n
         if (wbr(i,j).gt.thr) ibmax=j
         xsum=xsum+wbr(i,j)
      enddo
      write(iunit,'(i6,a4,1x,f6.3,9(4x,a2,i4,f6.3))',advance='no') &
         & i,asym(at(i)),xsum
      do j = 1, ibmax
         write(iunit,'(4x,a2,i4,f6.3)',advance='no') &
         & asym(at(imem(j))),imem(j),wbr(i,j)
      enddo
      write(iunit,'(a)')
   enddo

   deallocate(wbr,imem)

end subroutine print_wiberg

subroutine print_wbo_fragment(iunit,n,at,wbo,thr)
   use iso_fortran_env, wp => real64
   use tbdef_atomlist
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: wbo(n,n)
   real(wp),intent(in) :: thr

   type(tb_atomlist) :: atl

   real(wp),allocatable :: bond(:,:)
   integer, allocatable :: cn(:)
   integer, allocatable :: fragment(:)
   integer, allocatable :: list(:)
   character(len=2),external :: asym
   character(len=:),allocatable :: string
   integer  :: i,j,k,nfrag
   real(wp) :: xsum

   allocate( fragment(n),cn(n), list(n), source = 0 )
   allocate( bond(n,n), source = 0.0_wp )
   where( wbo > thr )
      bond = min(wbo,1.0_wp)
   elsewhere
      bond = 0.0_wp
   endwhere
   forall(i = 1:n) cn(i) = sum(ceiling(bond(:,i)))

   call mrec(nfrag,cn,bond,n,at,fragment)

   write(iunit,'(a)')
   if (nfrag > 1) then
      write(iunit,'(1x,"*",1x,i0,1x,a)',advance='no') &
         nfrag, "fragments found"
   else
      write(iunit,'(1x,"*",1x,a)',advance='no') &
         "no fragments found"
   endif
   write(iunit,'(1x,"(WBO >",f5.2,")")') thr
   write(iunit,'(a)')
   do i = 1, nfrag
      call atl%new
      call atl%add(fragment.eq.i)
      call atl%to_string(string)
      write(iunit,'(3x,a,"(",i0,"):",1x,a)') "fragment", i, string
   enddo

contains
   subroutine mrec(molcount,cn,bond,n,at,molvec)
      ! molcount: number of total fragments (increased during search)
      ! xyz: overall Cart. coordinates
      ! n: overall number of atoms
      ! at: atomic number array
      ! molvec: assignment vector of atom to fragment
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in)    :: cn(n)
      integer, intent(in)    :: n,at(n)
      integer, intent(inout) :: molvec(n),molcount
      real(wp),intent(inout) :: bond(n,n)
      logical, allocatable   :: taken(:)
      integer :: i
      allocate( taken(n) )
      molvec=0
      molcount=1
      taken=.false.
      do i=1,n
       if(.not.taken(i)) then
         molvec(i)=molcount
         taken(i)=.true.
         call neighbours(i,cn,at,taken,n,bond,molvec,molcount)
         molcount=molcount+1
      endif
      enddo
      molcount=molcount-1
   end subroutine mrec

   recursive subroutine neighbours(i,cn,at,taken,n,bond, &
      &                                molvec,molcnt)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in)    :: cn(n)
      real(wp),intent(inout) :: bond(n,n)
      integer, intent(in)    :: i,n,at(n)
      integer, intent(inout) :: molcnt,molvec(n)
      logical, intent(inout) :: taken(n)
      integer :: j,icn,k

      icn=cn(i)
      do k=1,icn
         j=maxloc(bond(:,i),1)
         bond(j,i)=0
         if (i .eq. j) cycle
         if (.not.taken(j)) then
            molvec(j)=molcnt
            taken(j)=.true.
            call neighbours(j,cn,at,taken,n,bond,molvec,molcnt)
         endif
      enddo
   end subroutine neighbours

end subroutine print_wbo_fragment

subroutine print_molpol(iunit,n,at,xyz,q,wf,g_a,g_c)
   use iso_fortran_env, wp => real64
   use dftd4
   use eeq_model
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: xyz(3,n)
   real(wp),intent(in) :: q(n)
   real(wp),intent(in) :: wf
   real(wp),intent(in) :: g_a
   real(wp),intent(in) :: g_c

   integer  :: i
   integer  :: dispdim
   real(wp) :: molpol,molc6,molc8
   real(wp),allocatable :: covcn(:)   ! covalent coordination number
   real(wp),allocatable :: gw(:)      ! gaussian weights for references
   real(wp),allocatable :: c6ref(:,:) ! unscaled reference C6
   real(wp),allocatable :: aw(:,:)    ! frequency dependent polarizibilities
   real(wp),allocatable :: c6ab(:,:)  ! actual C6 coeffients
   character(len=2),external :: asym

   call d4dim(n,at,dispdim)
   allocate( covcn(n), aw(23,n), c6ab(n,n), gw(dispdim), &
             c6ref(dispdim,dispdim), source = 0.0_wp )

   call covncoord(n,at,xyz,covcn,thr=1600.0_wp)
   call d4(n,dispdim,at,wf,g_a,g_c,covcn,gw,c6ref)
   call mdisp(n,dispdim,at,q,xyz,g_a,g_c,gw,c6ref, &
              molc6,molc8,molpol,aout=aw,cout=c6ab)

   write(iunit,'(a)')
   write(iunit,'("     #   Z   ")',advance='no')
   write(iunit,'("     covCN")',advance='no')
   write(iunit,'("         q")',advance='no')
   write(iunit,'("      C6AA")',advance='no')
   write(iunit,'("      α(0)")',advance='no')
   write(iunit,'(a)')
   do i=1,n
      write(iunit,'(i6,1x,i3,1x,a2)',advance='no') &
      &     i,at(i),asym(at(i))
      write(iunit,'(f10.3)',advance='no')covcn(i)
      write(iunit,'(f10.3)',advance='no')q(i)
      write(iunit,'(f10.3)',advance='no')c6ab(i,i)
      write(iunit,'(f10.3)',advance='no')aw(1,i)
      write(iunit,'(a)')
   enddo
   write(iunit,'(/,1x,"Mol. C6AA /au·bohr⁶  :",f18.6,'// &
   &            '/,1x,"Mol. C8AA /au·bohr⁸  :",f18.6,'// &
   &            '/,1x,"Mol. α(0) /au        :",f18.6,/)') &
   &             molc6,molc8,molpol

end subroutine print_molpol

subroutine print_dipole(iunit,n,at,xyz,z,nao,P,dpint)
   use iso_fortran_env, wp => real64
   use mctc_econv
  implicit none
  integer, intent(in) :: iunit
  integer, intent(in) :: n
  integer, intent(in) :: at(n)
  real(wp),intent(in) :: xyz(3,n)
  real(wp),intent(in) :: z(n)
  integer, intent(in) :: nao
  real(wp),intent(in) :: P(nao,nao)
  real(wp),intent(in) :: dpint(3,nao*(1+nao)/2)

  integer  :: i,j,k
  real(wp) :: d(3),dip

  ! core part
  d = 0.0_wp
  do i = 1, n
     d = d + xyz(:,i)*z(i)
  enddo

  ! contraction with P
  k = 0
  do i = 1, nao
     do j = 1, i-1
        k = k+1
        d = d - 2.0_wp*P(j,i)*dpint(:,k)
     enddo
     k = k+1
     d = d - P(i,i)*dpint(:,k)
  enddo

  dip = norm2(d)

  write(iunit,'(a)')
  write(iunit,'(1x,"dipole moment from electron density (au)")')
  write(iunit,'(1x,"    X       Y       Z   ")')
  write(iunit,'(3f9.4,"  total (Debye): ",f8.3)') &
       & d(1),   d(2),   d(3), dip*autod
  write(iunit,'(a)')

end subroutine print_dipole

subroutine print_spin_population(iunit,n,at,nao,focca,foccb,S,C,aoat2,lao2)
   use iso_fortran_env, wp => real64
   use scc_core, only : dmat, mpop
   implicit none
   integer, intent(in) :: iunit       ! STDOUT
   integer, intent(in) :: n           ! number of atoms
   integer, intent(in) :: at(n)       ! atom types
   integer, intent(in) :: nao         ! number of spherical atomic orbitals
   real(wp),intent(in) :: focca(nao)  ! fractional occupation numbers (alpha)
   real(wp),intent(in) :: foccb(nao)  ! fractional occupation numbers (beta)
   real(wp),intent(in) :: S(nao,nao)  ! overlap matrix
   real(wp),intent(in) :: C(nao,nao)  ! eigenvector/orbitals
   integer, intent(in)  :: aoat2(nao)
   integer, intent(in)  :: lao2(nao)

   integer :: i
   character(len=2),external :: asym
   real(wp),allocatable :: tmp(:)
   real(wp),allocatable :: q(:)
   real(wp),allocatable :: qlmom(:,:)
   real(wp),allocatable :: X(:,:)

   allocate( tmp(nao),q(n),qlmom(3,n),X(nao,nao), source = 0.0_wp )

   write(iunit,'("(R)spin-density population")')
   tmp = focca - foccb
   call dmat(nao,tmp,C,X) ! X is scratch
   call mpop(n,nao,aoat2,lao2,S,X,q,qlmom)
   write(iunit,'(a)')
   write(iunit,'(1x,"Mulliken population n(s)   n(p)   n(d)")')
   do i=1,n
      write(iunit,'(i6,a3,1f8.4,1x,4f7.3)') &
         &  i,asym(at(i)),q(i),qlmom(1,i),qlmom(2,i),qlmom(3,i)
   enddo

end subroutine print_spin_population

subroutine print_fod_population(iunit,ifile,n,at,nao,S,C,etemp,emo,ihomoa,ihomob, &
      &                         aoat2,lao2)
   use iso_fortran_env, wp => real64
   use mctc_econv
   use scc_core
   implicit none
   integer, intent(in) :: iunit       ! STDOUT
   integer, intent(in) :: ifile       ! file handle for printout of FOD population
   integer, intent(in) :: n           ! number of atoms
   integer, intent(in) :: at(n)       ! atom types
   integer, intent(in) :: nao         ! number of spherical atomic orbitals
   real(wp),intent(in) :: S(nao,nao)  ! overlap matrix
   real(wp),intent(in) :: C(nao,nao)  ! eigenvector/orbitals
   real(wp),intent(in) :: etemp       ! electronic temperature
   real(wp),intent(in) :: emo(nao)    ! orbital energies
   integer, intent(in) :: ihomoa      ! position of HOMO in alpha space
   integer, intent(in) :: ihomob      ! position of HOMO in beta space
   integer, intent(in)  :: aoat2(nao)
   integer, intent(in)  :: lao2(nao)

   integer :: i
   character(len=2),external :: asym
   real(wp),allocatable :: focc(:)    ! fractional occupation numbers
   real(wp),allocatable :: q(:)       ! FOD populations
   real(wp),allocatable :: qlmom(:,:) ! FOD populations per shell
   real(wp),allocatable :: X(:,:)     ! Loewdin orthonormalizer
   real(wp),allocatable :: focca(:)   ! fractional occupation numbers (alpha)
   real(wp),allocatable :: foccb(:)   ! fractional occupation numbers (beta)
   real(wp) :: efa,efb,ga,gb,nfoda,nfodb

   allocate( q(n), qlmom(3,n), X(nao,nao), focca(nao), foccb(nao), focc(nao), &
             source = 0.0_wp )

   call makel(nao, S, C, X)
   if(ihomoa+1.le.nao) &
      call fermismear(.false.,nao,ihomoa,etemp,emo,focca,nfoda,efa,ga)
   if(ihomob+1.le.nao) &
      call fermismear(.false.,nao,ihomob,etemp,emo,foccb,nfodb,efb,gb)
   call fodenmak(.true.,nao,emo * evtoau,focca,efa)
   call fodenmak(.true.,nao,emo * evtoau,foccb,efb)

   focc = focca+foccb
   write(iunit,'(/,"NFOD :",1x,F10.4)') sum(focc)
   q=0
   qlmom=0
   call lpop(n,nao,aoat2,lao2,focc,X,1.0d0,q,qlmom)
   write(iunit,'(a)')
   write(iunit,'(" Loewdin FODpop     n(s)   n(p)   n(d)")')
   do i = 1, n
      write(iunit,'(i6,a3,f8.4,1x,4f7.3)') &
         i,asym(at(i)),q(i),qlmom(1,i),qlmom(2,i),qlmom(3,i)
   enddo
   if (ifile.ne.-1) then
      do i = 1, n
         write(ifile,'(F14.8)') q(i)
      enddo
   endif

end subroutine print_fod_population

subroutine print_thermo(iunit,nat,nvib_in,at,xyz,freq,etot,htot,gtot,nimag,pr)
   use iso_fortran_env, only : wp => real64
   use mctc_econv
   use readin
   use setparam
   use axis_trafo, only : axis2
   use thermo
   implicit none
   integer, intent(in) :: iunit
   logical, intent(in) :: pr
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   integer, intent(in) :: nvib_in
   real(wp),intent(in) :: freq(3*nat)
   real(wp),intent(in) :: xyz(3,nat)
   real(wp),intent(in) :: etot
   real(wp),intent(out) :: gtot
   real(wp),intent(out) :: htot

   real(wp) xx(10),sthr,temp,scale_factor
   real(wp) aa,bb,cc,vibthr,ithr
   real(wp) escf,symnum,wt,avmom,zp,diff
   real(wp) :: omega,maxfreq,fswitch,lnq_r,lnq_v
   real(wp),allocatable :: et(:),ht(:),gt(:),ts(:)
   integer nn,nvib,i,j,k,n,nvib_theo,isthr
   integer, intent(out) :: nimag
   real(wp),allocatable :: vibs(:),tmp(:)
   character(len=*),parameter :: outfmt = &
      '(9x,"::",1x,a,f24.12,1x,a,1x,"::")'
   character(len=*),parameter :: dblfmt = &
      '(10x,":",2x,a,f24.7,1x,a,1x,":")'
   character(len=*),parameter :: intfmt = &
      '(10x,":",2x,a,i24,       6x,":")'
   character(len=*),parameter :: chrfmt = &
      '(10x,":",2x,a,a24,       6x,":")'

   logical linear,atom,da

   allocate( et(nthermo), ht(nthermo), gt(nthermo), ts(nthermo), &
      &      vibs(3*nat), tmp(3*nat), source = 0.0_wp )

   ! frequencies read in are considered
   ! as being real if .gt. this value in cm-1
   ! this threshold requires projected freqs.!
   vibthr=1.0
   ithr=-20

   atom=.false.
   linear=.false.
   sthr=thermo_sthr
   scale_factor=1.0
   nvib=0
   nimag=0

   call axis2(nat,at,xyz,aa,bb,cc,avmom,wt)

   nvib_theo=3*nat-6
   if(cc.lt.1.d-10) linear=.true.
   if(linear) nvib_theo=3*nat-5

   if(aa+bb+cc.lt.1.d-6)then
      atom=.true.
      nvib=0
      nvib_theo=0
   endif

   ! the rotational number
   call getsymmetry(pr,iunit,nat,at,xyz,desy,maxatdesy,pgroup)
   call getsymnum(pgroup,linear,symnum)

   vibs=0
   do i=1,3*nat
      if(abs(freq(i)).gt.vibthr.and.i.le.nvib_in)then
         nvib=nvib+1
         vibs(nvib)=freq(i)
      endif
   enddo

   ! scale
   vibs(1:nvib)=vibs(1:nvib)*scale_factor

   do i=1,nvib
      ! artifacts
      if(vibs(i).lt.0.and.vibs(i).gt.ithr) then
         vibs(i)=-vibs(i)
         if(pr)write(iunit,*)'inverting freq ',i,vibs(i)
      endif
   enddo
   tmp=vibs

   k=nvib
   nvib=0
   j=0
   diff = abs(maxval(vibs) - thermo_sthr)
   do i=1,k
      if(tmp(i).gt.0) then
         nvib=nvib+1
         if (abs(tmp(i) - thermo_sthr) < diff) then
            diff = abs(tmp(i) - thermo_sthr)
            isthr = nvib
         endif
         vibs(nvib)=tmp(i)*rcmtoau ! work in atomic units, seriously
      else
         j=j+1
      endif
   enddo
   nimag=j

   if(pr)then
      write(iunit,'(a)')
      write(iunit,'(10x,51("."))')
      write(iunit,'(10x,":",22x,a,22x,":")') "SETUP"
      write(iunit,'(10x,":",49("."),":")')
      write(iunit,intfmt) "# frequencies    ",nvib
      write(iunit,intfmt) "# imaginary freq.",nimag
      write(iunit,chrfmt) "linear?          ",bool2string(linear)
      write(iunit,chrfmt) "only rotor calc. ",bool2string(nvib.eq.0)
      write(iunit,chrfmt) "symmetry         ",trim(pgroup)
      write(iunit,intfmt) "rotational number",int(symnum)
      write(iunit,dblfmt) "scaling factor   ",scale_factor,"    "
      write(iunit,dblfmt) "rotor cutoff     ",thermo_sthr, "cm⁻¹"
      write(iunit,dblfmt) "imag. cutoff     ",ithr, "cm⁻¹"
      write(iunit,'(10x,":",49("."),":")')
   endif

   call print_thermo_sthr_ts(iunit,nvib,vibs,avmom,thermo_sthr,thermotemp(nthermo))

   ! do calc.
   zp = 0.5_wp * sum(vibs(1:nvib))
   do i = 1, nthermo
      temp=thermotemp(i)
      call thermodyn(iunit,aa,bb,cc,avmom,linear,atom,symnum,wt,vibs,nvib,escf, &
         & temp,sthr,et(i),ht(i),gt(i),ts(i),zp,pr)
      !call oldthermo(aa,bb,cc,avmom,linear,atom,symnum,wt,vibs,nvib,escf, &
      !   & temp,sthr,et(i),ht(i),gt(i),ts(i),zp,pr)
   enddo

   write(iunit,'(a)')
   write(iunit,'(a10)',advance='no') "T/K"
   write(iunit,'(a16)',advance='no') "H(0)-H(T)+PV"
   write(iunit,'(a16)',advance='no') "H(T)/Eh"
   write(iunit,'(a16)',advance='no') "T*S/Eh"
   write(iunit,'(a16)',advance='no') "G(T)/Eh"
   write(iunit,'(a)')
   write(iunit,'(3x,72("-"))')
   do i = 1, nthermo
      write(iunit,'(3f10.2)',advance='no') thermotemp(i)
      write(iunit,'(3e16.6)',advance='no') ht(i)
      write(iunit,'(3e16.6)',advance='no') et(i)
      write(iunit,'(3e16.6)',advance='no') ts(i)
      write(iunit,'(3e16.6)',advance='no') gt(i)
      if (i.eq.nthermo .and. nthermo.gt.1) then
         write(iunit,'(1x,"(used)")')
      else
         write(iunit,'(a)')
      endif
   enddo
   write(iunit,'(3x,72("-"))')

   gtot = gt(nthermo)
   htot = et(nthermo)

   write(iunit,'(a)')
   write(iunit,'(9x,53(":"))')
   write(iunit,'(9x,"::",18x,a,18x,"::")') "THERMODYNAMIC"
   write(iunit,'(9x,53(":"))')
   write(iunit,outfmt) "total free energy ", gtot+etot,"Eh  "
   write(iunit,'(9x,"::",49("."),"::")')
   write(iunit,outfmt) "total energy      ",      etot,"Eh  "
   write(iunit,outfmt) "zero point energy ",        zp,"Eh  "
   write(iunit,outfmt) "G(RRHO) w/o ZPVE  ",   gtot-zp,"Eh  "
   write(iunit,outfmt) "G(RRHO) contrib.  ",      gtot,"Eh  "
   write(iunit,'(9x,53(":"))')

end subroutine print_thermo

subroutine print_thermo_sthr_lnq(iunit,nvib,vibs,avmom,sthr,temp)
   use iso_fortran_env, wp => real64
   use mctc_econv
   use thermo
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: nvib
   real(wp),intent(in) :: vibs(nvib)
   real(wp),intent(in) :: avmom
   real(wp),intent(in) :: sthr
   real(wp),intent(in) :: temp

   integer  :: i
   real(wp) :: maxfreq,omega,lnq_r,lnq_v,fswitch

   write(iunit,'(a)')
   maxfreq = max(300.0_wp,chg_inverted(0.99_wp,sthr))
   write(iunit,'(a8,a14,a12,10x,a12,10x,a12)') &
      "mode","ω/cm⁻¹","ln{qvib}","ln{qrot}","ln{qtot}"
   write(iunit,'(3x,72("-"))')
   do i = 1, nvib
      omega = vibs(i)*autorcm
      lnq_r = lnqvib(temp,omega)
      lnq_v = lnqrot(temp,omega,avmom)
      fswitch = 1.0_wp - chg_switching(omega,sthr)
      if (omega > maxfreq) exit
      write(iunit,'(i8,f10.2,2(f12.5,1x,"(",f6.2,"%)"),f12.5)') &
         i,omega,lnq_v,(1.0_wp-fswitch)*100, &
         lnq_r,fswitch*100,(1.0_wp-fswitch) * lnq_v + fswitch * lnq_r
   enddo
   write(iunit,'(3x,72("-"))')

end subroutine print_thermo_sthr_lnq

subroutine print_thermo_sthr_ts(iunit,nvib,vibs,avmom_si,sthr_rcm,temp)
   use iso_fortran_env, wp => real64
   use mctc_constants
   use mctc_econv
   use thermo
   implicit none
 
   integer, intent(in) :: iunit      !< output unit, usually STDOUT
   integer, intent(in) :: nvib       !< number of frequencies
   real(wp),intent(in) :: vibs(nvib) !< frequencies in Eh
   real(wp),intent(in) :: avmom_si   !< average moment
   real(wp),intent(in) :: sthr_rcm   !< rotor cutoff
   real(wp),intent(in) :: temp       !< temperature

   integer  :: i
   real(wp) :: maxfreq,omega,s_r,s_v,fswitch
   real(wp) :: beta,xxmom,e,ewj,mu,RT,sthr,avmom
   beta = 1.0_wp/kB/temp ! beta in 1/Eh
   sthr = sthr_rcm * rcmtoau ! sthr in Eh
   RT = kb*temp*autokcal ! RT in kcal/mol for printout
   avmom = avmom_si*kgtome*aatoau**2*1.0e+20_wp ! in me·α²

   write(iunit,'(a)')
   maxfreq = max(300.0_wp,chg_inverted(0.99_wp,sthr_rcm))
   write(iunit,'(a8,a14,1x,a27,a27,a12)') &
      "mode","ω/cm⁻¹","T·S(HO)/kcal·mol⁻¹","T·S(FR)/kcal·mol⁻¹","T·S(vib)"
   write(iunit,'(3x,72("-"))')
   do i = 1, nvib
      ! frequency is Eh
      omega=vibs(i)
      ! omega in Eh, beta in 1/Eh
      ewj=exp(-omega*beta)
      ! moment of intertia corresponding to the rotor with frequency omega
      ! mu is in me·α² (au)
      mu = 0.5_wp / (omega+1.0e-14_wp)
      ! this reduced moment limits the rotational moment of inertia for
      ! this vibration to that of the total molecule rotation/3
      ! avmom and mu are in au
      mu=mu*avmom/(mu+avmom)
      !              free rotor entropy
      ! Cramer, page 328 for one degree of freedom or
      ! http://cccbdb.nist.gov/thermo.asp, eq. 35, sigma=1
      !              harm. osc. entropy
      if(omega.gt.0)then
         ! this is S/R which is dimensionless
         s_v = omega*beta*ewj/(1.0_wp-ewj) - log(1.0_wp-ewj)
         s_r = 0.5_wp + log(sqrt(pi/beta*2.0_wp*mu))
      else
         s_v = 0.0_wp
         s_r = 0.0_wp
      endif
      ! Head-Gordon weighting
      fswitch=1.0_wp-chg_switching(omega,sthr)
      if (omega > maxfreq*rcmtoau) exit
      write(iunit,'(i8,f10.2,2(f12.5,1x,"(",f6.2,"%)"),f12.5)') &
         i,omega*autorcm,-RT*s_v,(1.0_wp-fswitch)*100, &
         -RT*s_r,fswitch*100,-RT*((1.0_wp-fswitch) * s_v + fswitch * s_r)
   enddo
   write(iunit,'(3x,72("-"))')

end subroutine print_thermo_sthr_ts

subroutine print_gbsa_info(iunit,gbsa)
   use iso_fortran_env, wp => real64
   use mctc_constants
   use mctc_econv
   use gbobc
   implicit none
   integer, intent(in) :: iunit
   type(tb_solvent), intent(in) :: gbsa

   integer :: i
   character(len=2),external :: asym

   write(iunit,'(a)')
   write(iunit,'(1x,"*",1x,a)') &
      &  "generalized Born model for continuum solvation"
   write(iunit,'(a)')
   write(iunit,'(2x,2a4,3x,3a)') "#","Z","Born rad/Å","   SASA/Å²","    H-bond"
   do i = 1, gbsa%nat
      write(iunit,'(i6,1x,i3,1x,a2,3f10.3)') &
         &  i,gbsa%at(i),asym(gbsa%at(i)), &
         &  gbsa%brad(i)*autoaa,gbsa%sasa(i)*fourpi/gbsa%gamsasa(i)*autoaa**2, &
         &  gbsa%hbw(i)
   enddo
   write(iunit,'(/,1x,"total SASA / Å² :",f13.3)') &
      &  sum(gbsa%sasa/gbsa%gamsasa)*fourpi*autoaa**2
      

end subroutine print_gbsa_info

end module property_output

subroutine print_orbital_eigenvalues(iunit,wfn,range)
   use iso_fortran_env, wp => real64
   use mctc_econv
   use tbdef_wavefunction
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: range
   type(tb_wavefunction),intent(in) :: wfn
   character(len=*),parameter :: hlfmt = '(    a24,f21.7,1x,"Eh",f18.4,1x,"eV")'
   integer :: maxorb,minorb,iorb
   real(wp) :: gap

   minorb = max(wfn%ihomoa - (range+1), 1)
   maxorb = min(wfn%ihomoa +  range, wfn%nao)
   gap = wfn%emo(wfn%ihomoa+1) - wfn%emo(wfn%ihomoa)

   write(iunit,'(a)')
   write(iunit,'(a10,a14,a21,a21)') "#","Occupation","Energy/Eh","Energy/eV"
   write(iunit,'(6x,61("-"))')
   if (minorb .gt. 1) then
      call write_line(1,wfn%focc,wfn%emo,wfn%ihomo)
      if (minorb .gt. 2) &
         write(iunit,'(a10,a14,a21,a21)') "...","...","...","..."
   endif
   do iorb = minorb,maxorb
      call write_line(iorb,wfn%focc,wfn%emo,wfn%ihomo)
   enddo
   if (maxorb .lt. wfn%nao) then
      if (maxorb .lt. wfn%nao-1) then
         if (wfn%focc(maxorb) > 1.0e-7_wp) then
            write(iunit,'(a10,a14,a21,a21)') "...","...","...","..."
         else
            write(iunit,'(a10,a14,a21,a21)') "...",   "","...","..."
         endif
      endif
      call write_line(wfn%nao,wfn%focc,wfn%emo,wfn%ihomo)
   endif
   write(iunit,'(6x,61("-"))')
   write(iunit,hlfmt) "HL-Gap",gap*evtoau,gap
   write(iunit,hlfmt) "Fermi-level",(wfn%efa+wfn%efb)/2*evtoau,(wfn%efa+wfn%efb)/2
contains
subroutine write_line(iorb,focc,emo,ihomo)
   integer, intent(in) :: iorb
   integer, intent(in) :: ihomo
   real(wp),intent(in) :: focc(:)
   real(wp),intent(in) :: emo (:)
   character(len=*),parameter :: mofmt = '(i10,f14.4,f21.7,f21.4)'
   character(len=*),parameter :: vofmt = '(i10,14x,  f21.7,f21.4)'
   if (focc(iorb) < 1.0e-7_wp) then
      write(iunit,vofmt,advance='no') iorb,             emo(iorb)*evtoau, emo(iorb)
   else
      write(iunit,mofmt,advance='no') iorb, focc(iorb), emo(iorb)*evtoau, emo(iorb)
   endif
   if (iorb == ihomo) then
      write(iunit,'(1x,"(HOMO)")')
   elseif (iorb == ihomo+1) then
      write(iunit,'(1x,"(LUMO)")')
   else
      write(iunit,'(a)')
   endif
end subroutine write_line
end subroutine print_orbital_eigenvalues
