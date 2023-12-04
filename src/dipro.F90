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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

!> Implementation of the dimer projection method for extended tight binding methods.
module xtb_dipro
   use mctc_env, only : error_type, fatal_error, get_argument, wp
   use mctc_io, only : structure_type, new
   use mctc_io_convert, only : autoev
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule, assignment(=)
   use xtb_type_calculator, only : TCalculator
   use xtb_tblite_calculator, only : TTBLiteCalculator, TTBLiteInput, newTBLiteCalculator
   use xtb_setparam
#if WITH_TBLITE
   use xtb_dipro_bondorder, only : get_wiberg_bondorder
   use xtb_dipro_fragment, only : get_wiberg_fragment
   use xtb_dipro_output, only : format_list, to_string
   use xtb_dipro_xtb, only : get_calculator
   use tblite_basis_type, only : get_cutoff, basis_type
   use tblite_blas, only : dot, gemv, gemm
   use tblite_context_type, only : context_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_overlap, only : get_overlap
   use tblite_output_format, only : format_string
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction, &
      & get_density_matrix
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
#endif
   implicit none
   private
   public :: get_jab, jab_input

   !> Configuration data for calculation
   type :: jab_input
      !> Flag for evoking DIPRO
      logical :: diprocalc
      !>Orbital degeneracy threshold
      real(wp) :: othr
      !> Name of the requested tight binding method
      character(len=:), allocatable :: method
      !> List of fragments, generated if not given here
      integer, allocatable :: fragment(:)
      !> Threshold for generating fragments from bond orders
      real(wp) :: thr = 0.1_wp
      !> Calculation accuracy
      real(wp) :: accuracy = 1.0_wp
      !> Output verbosity
      integer :: verbosity = 2
      !> total j_ab,eff in eV for test module
      real(wp) :: totjab(3)
   end type jab_input

   character(len=*), parameter :: source = 'xtb_dipro'
   !> Conversion factor from temperature to energy (Boltzmann's constant in atomic units)
   real(wp), parameter :: ktoau = 3.166808578545117e-06_wp

contains

!> Entry point for calculation of dimer projection (DIPRO) related properties
subroutine get_jab(env, tblite, mol, fragment, dipro)
   use, intrinsic :: iso_fortran_env, only : output_unit
   implicit none
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol  
   !requested gfn method for calculations Input
   type(TTBLiteInput), intent(inout) :: tblite

   type(jab_input),intent(inout) :: dipro

   type(TEnvironment),intent(inout) :: env

   integer, allocatable, intent(inout) :: fragment(:)
#if WITH_TBLITE
   !> structure_type /= molecular structure
   type(structure_type) :: struc
   type(error_type),allocatable :: error

   type(context_type) :: ctx
   type(basis_type) :: bas
   !> fcalc is =xcalc just for fragments
   type(xtb_calculator) :: xcalc
   type(xtb_calculator), allocatable :: fcalc(:)
   !> mfrag is =struc just for fragments
   type(structure_type), allocatable :: mfrag(:)
   type(wavefunction_type) :: wfn
   !> wfx is =wfn just for fragments
   type(wavefunction_type), allocatable :: wfx(:)

   !> Molecular gradient, strain derivatives
   real(wp), allocatable :: gradient(:, :), sigma(:,:), nel(:)
   real(wp), allocatable :: overlap(:, :), trans(:, :), wbo(:, :), chrg(:), p2mat(:,:), coeff2(:,:),loc(:,:)
   real(wp), allocatable :: orbital(:, :, :), scmat(:, :), fdim(:, :), scratch(:), efrag(:),y(:,:),Edim(:,:)
   integer, allocatable :: spinfrag(:), start_index(:),end_index(:),orbprint(:)

   integer :: charge, stat, unit, ifr, nfrag, nao, i, j, k

   logical :: exist

   character(3) :: adv='NO '

   real(wp) :: energy, cutoff, jab, sab, jeff, Vtot(3)

!======================================================================================

   struc=mol

   if ( tblite%method == '' ) then 
      tblite%method = 'gfn2'
      call env%warning("No method provided, falling back to default GFN2-xTB.", source)
   end if

!=========================set up calculator===========================================   

   call get_calculator(xcalc, struc, tblite%method, error)  
   call new_wavefunction(wfn, struc%nat, xcalc%bas%nsh, xcalc%bas%nao, & 
      & 1, set%etemp * ktoau)
   wfn%nspin=1

!=========================print Header===============================================

   call generic_header(6,'D I P R O',49,10)

!=========================calculation for dimer======================================             

   write(*,'(A)') "Calculation for dimer "
   write(*,'(A)') "--------------------- "
   write(*,'(A)') "  "

   write(*,'(A,F4.0)') "charge of dimer : ", mol%chrg
   write(*,'(A,I2)') "unpaired e- of dimer : ", set%nalphabeta

   call xtb_singlepoint(ctx, struc, xcalc, wfn, tblite%accuracy, energy,gradient,sigma,2)
   if (ctx%failed()) then
       call env%error("Single point calculation for dimer failed.", source)
      return
   end if

   allocate(overlap(xcalc%bas%nao, xcalc%bas%nao),y(xcalc%bas%nao,2))
   cutoff = get_cutoff(xcalc%bas)
   call get_lattice_points(struc%periodic, struc%lattice, cutoff, trans) 
   call get_overlap(struc, trans, cutoff, xcalc%bas, overlap)  

!==================set up fragments if not given by xcontrol=========================

   if (allocated(fragment)) then
   else
      !> wfn%density is [nao,nao,spin], pmat in wiberg bondorder is [nao,nao], thus pmat2 is introduced here     
      allocate(wbo(struc%nat, struc%nat),p2mat(xcalc%bas%nao,xcalc%bas%nao))
      do i=1,xcalc%bas%nao
         do j=1,xcalc%bas%nao
            p2mat(i,j)=wfn%density(i,j,1)
         end do
      end do
      call get_wiberg_bondorder(xcalc%bas, overlap, p2mat, wbo)

      allocate(fragment(struc%nat))
      call get_wiberg_fragment(fragment, wbo, 0.1_wp)
   end if

   nfrag = maxval(fragment)
   select case(nfrag)
   case(:1)
      call ctx%message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      call ctx%message("Found no fragments in the input structure.")
      call ctx%message("Aborting...")
      call env%error("Found no fragments in input structure.", source)
      return
   case(2:)
      call ctx%message("Found "//to_string(nfrag)//" fragments!")
      do ifr = 1, nfrag
         call ctx%message("Fragment "//to_string(ifr)//":  "//format_list(fragment == ifr))
      end do
      call ctx%message("")
      if (nfrag > 2) then
         call ctx%message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         call ctx%message("Found "//to_string(nfrag)// &
            & " fragments in the input structure, too many fragments.")
         call ctx%message("Aborting...")
         call env%error("Found too many fragments in input structure.", source)
         return
      end if
   end select

   nao = size(wfn%emo)
   allocate(orbital(nao, nfrag,nao), efrag(nfrag), scmat(nao, nao), fdim(nao, nao), &
      & scratch(nao), mfrag(nfrag), wfx(nfrag), chrg(nfrag), spinfrag(nfrag), &
      & start_index(nfrag),end_index(nfrag),orbprint(nfrag),nel(nfrag))

!==================================external files CHRG & UHF read-in====================================

   inquire(file='.UHFfrag', exist=exist)
   if (exist) then
      open(file='.UHFfrag', newunit=unit)
      write(output_unit, '(a)') "[Info] Fragment spin read from .UHFfrag"
      read(unit,*,iostat=stat) spinfrag(1),spinfrag(2)
      close(unit)
   else if (set%nalphabeta .ne. 0) then
      call ctx%message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      call ctx%message("Total spin specified but no fragment spins. Spins are not determined &
              &automatically by xtb. Please set up .UHFfrag or use total spin =0." )
      call ctx%message("Aborting...")
      call env%error("Total spin specified but no fragment spins.", source)
      return
   else
      spinfrag=0
   end if

   inquire(file='.CHRGfrag', exist=exist)
   if (exist) then
      open(file='.CHRGfrag', newunit=unit)
      write(output_unit, '(a)') "[Info] Fragment charge read from .CHRGfrag"
      read(unit,*,iostat=stat) chrg(1),chrg(2)
      close(unit)
   else
      chrg=0
   end if

!=================================fragment calculations=============================================   

   allocate(fcalc(nfrag))
   do ifr = 1, nfrag 
      call ctx%message("Calculation for fragment "//to_string(ifr))
      write(*,*) "------------------------------"
      write(*,*) "  "
      call get_structure_fragment(mfrag(ifr), struc, fragment == ifr)

      !> summation of fragment charges stored in chrg(nfrag)
      if ( all(chrg.eq.0) ) then
         chrg(ifr)=0
         !> fragment mask generated on the fly 
         chrg(ifr)=sum(pack(wfn%qat(:,1), fragment == ifr)) !> wfn%qat is [nat,nspin=1]
         mfrag(ifr)%charge=nint(chrg(ifr))
      else
         mfrag(ifr)%charge=nint(chrg(ifr))
      end if
      write(*,'(A,F4.0)') "charge of fragment : ", mfrag(ifr)%charge

      !> uhf fragments spins
      mfrag(ifr)%uhf = spinfrag(ifr)
      write(*,'(A,I2)') "unpaired e- of fragment : ", mfrag(ifr)%uhf

      call get_calculator(fcalc(ifr), mfrag(ifr), tblite%method, error)
      !> mol%charge is updated automatically from wfn by tblite library 
      call new_wavefunction(wfx(ifr), mfrag(ifr)%nat, fcalc(ifr)%bas%nsh, fcalc(ifr)%bas%nao, &
         & 1, set%etemp * ktoau)

     !> mol%type (dimer) == mfrag%type (fragments), wfn (dimer) == wfx (fragments), calc (dimer)==fcalc(fragments)
      wfx%nspin=1
      call xtb_singlepoint(ctx, mfrag(ifr), fcalc(ifr), wfx(ifr), tblite%accuracy, energy)
      if (ctx%failed()) then
         call env%error("Single point calculation for fragment failed.", source)
         return
      end if

      nel(ifr)=wfx(ifr)%nel(1)+wfx(ifr)%nel(2)

!==================================DIPRO==================================================

      do j = 1, fcalc(ifr)%bas%nao
         !> coeff is [nao,nao,spin=1]
         call unpack_coeff(xcalc%bas, fcalc(ifr)%bas, orbital(:, ifr, j), &
         & wfx(ifr)%coeff(:, j,1), fragment == ifr)
      end do 
   end do

   allocate(coeff2(nao,nao),Edim(nao,nao))
   Edim=0
   do i=1,nao
      Edim(i,i)=wfn%emo(i,1)
      do j=1,nao
         coeff2(i,j)=wfn%coeff(i,j,1) 
      end do
   end do

   start_index = -1
   end_index = -1

   !> find out which orbitals are within [HOMO-othr,LUMO+othr] and should be considered for DIPRO
   write(*,*) "  "
   write(*,*) "  ::::::::::::::::::::::::::::::::::::::::::::::"
   write(*,*) "  ::     D I P R O      C O U P L I N G S     ::"
   write(*,*) "  ::::::::::::::::::::::::::::::::::::::::::::::"
   write(*,*) "  "

   call ctx%message("energy threshhold for near-degenerate orbitals near HOMO and LUMO &
           &considered for DIPRO: "//format_string(dipro%othr, '(f6.3)')//" eV")

   do ifr=1,nfrag
      do j = 1, fcalc(ifr)%bas%nao
         if (wfx(ifr)%emo(j,1) .ge. (wfx(ifr)%emo(wfx(ifr)%homo(2),1) - dipro%othr/autoev) .and.&
           & wfx(ifr)%emo(j,1) .le. (wfx(ifr)%emo(wfx(ifr)%homo(2)+1,1) + dipro%othr/autoev)) then
           if (start_index(ifr).eq.-1) then 
                   start_index(ifr) = j
              end if
              end_index(ifr) = j
         endif
      end do
      write(*,'(A, I2, I10)') "no. of orbitals within energy window of frag", ifr, end_index(ifr)-start_index(ifr)+1
   end do
   write(*,*) "  "

!========================================DIPRO equations===========================================
!> equations after B. Baumeier, J. Kirkpatrick, D. Andrienko, PCCP 2010, 12, 11103.
!> y_A^i = C_A^i * S_AB * C_AB                  orbital projection of monomer A onto dimer
!> y_B^j = C_B^j * S_AB * C_AB                  orbital projection of monomer B onto dimer
!> S_ab^ij = y_A^i * y_B^j                      projected overlap
!> E_A^i = y_A^i * E_AB * y_A^i                 site energy monomeCALL execute_command_line('wc -l < file.txt > wc.txt' )r A
!> E_B^j = y_B^j * E_AB * y_B^j                 site energy monomer B
!> J_AB^ij = y_A^i * E_AB * y_B^j               coupling integral J_AB
!> J_AB^ij effective = (J_AB^ij - (E_A^i + E_B^j) / 2 * S_ab^ij) / (1 - S_ab^ij * S_ab^ij)            
!> A,B: monomers   AB: dimer
!> i,j: orbitals of the monomers
!> C: orbital coefficients     S: overlap matrix      E: orbital energy of the dimer

   Vtot=0
   !> gemm(amat,bmat,cmat,transa,transb,a1,a2): X=a1*Amat*Bmat+a2*Cmat
   !> scmat=S_dim*C_dim
   call gemm(overlap, coeff2, scmat)
   do j = start_index(1), end_index(1)
      orbprint(1)=wfx(1)%homo(max(2,1))-j

      y(:,1)=0
      !> gemv(amat, xvec,yvec,a1,a2,transa): X=a1*Amat*xvec+a2*yvec
      !> y_mon1(ifr)=C_mon1(j)*S_dim*C_dim
      call gemv( scmat, orbital(:, 1, j), y(:,1), trans="t" )
      !> scratch=E_dim*y1(j)
      call gemv( Edim, y(:,1), scratch )
      !> E_mon=y1(j)*E_dim*y1(j)
      efrag(1)=dot( y(:,1), scratch)

      do k = start_index(2), end_index(2)
         orbprint(2)=wfx(2)%homo(max(2,1))-k

         y(:,2)=0
         !> y_mon2(ifr)=C_mon2(k)*S_dim*C_dim
         call gemv( scmat, orbital(:, 2, k), y(:,2), trans="t" )
         !> scratch=E_dim*y2(k)
         call gemv( Edim, y(:,2), scratch )
         !> E_mon=y(ifr)*E_dim*y2(k)
         efrag(2)=dot( y(:,2), scratch)

         !> sab=y1(j)*y2(k)
         sab=dot( y(:,1), y(:,2) )
         !> jab=y1(j)*E_dim*y2(k)
         jab=dot( y(:,1), scratch )

         jeff = (jab - sum(efrag) / nfrag * sab) / (1.0_wp - sab**2)

!=======================================Printout============================================         

         do ifr=1,2

            !> both coupling monomer orbitals are printed in the same line
            if (ifr.eq.2) then 
               adv='YES'
            else
               adv='NO '
            end if

            !> if the number of electrons is odd, the SOMO is the new HOMO
            if (mod(nint(nel(ifr)),2).eq.1) then
               if (orbprint(ifr).gt.0) then
                    write(*, '( "  Fragment ", I1, " SOMO-", I1 )', ADVANCE=adv ) ifr, orbprint(ifr)
                  else if (orbprint(ifr).eq.0) then
                    write(*, '( "  Fragment ", I1, " SOMO ")', ADVANCE=adv ) ifr
                  else if (orbprint(ifr).eq.-1) then
                    write(*, '( "  Fragment ", I1, " LUMO ")', ADVANCE=adv ) ifr
                  else if (orbprint(ifr).lt.-1) then
                    write(*, '( "  Fragment ", I1, " LUMO+", I1 )', ADVANCE=adv ) ifr, abs(orbprint(ifr))-1
               end if
            else if (mod(nint(nel(ifr)),2).eq.0) then
               if (orbprint(ifr).gt.0) then
                    write(*, '( "  Fragment ", I1, " HOMO-", I1 )', ADVANCE=adv ) ifr, orbprint(ifr)
                  else if (orbprint(ifr).eq.0) then
                    write(*, '( "  Fragment ", I1, " HOMO ")', ADVANCE=adv ) ifr
                  else if (orbprint(ifr).eq.-1) then
                    write(*, '( "  Fragment ", I1, " LUMO ")', ADVANCE=adv ) ifr
                  else if (orbprint(ifr).lt.-1) then
                    write(*, '( "  Fragment ", I1, " LUMO+", I1 )', ADVANCE=adv ) ifr, abs(orbprint(ifr))-1
               end if
            end if
         end do
         write(*,'(A)') "--------------------------------------"

         call ctx%message("E_mon(orb) frag1 frag2"//format_string(efrag(1)*autoev, '(f20.3)')// &
           &format_string(efrag(2)*autoev, '(f10.3)')//" eV")
         call ctx%message("J(AB): "//format_string(jab*autoev, '(f20.3)')//" eV")    
         call ctx%message("S(AB): "//format_string(sab, '(f22.8)')//" Eh")
         call ctx%message("|J(AB,eff)|: "//format_string(abs(jeff)*autoev, '(f16.3)')//" eV")
         write(*,*) "  "

         !> Vtotal=sqrt(sum(jeff^2)) for 1. occ-->occ; 2. virt-->virt; 3. occ-->virt/virt-->occ transitions
         if(orbprint(1).ge.0.and.orbprint(2).ge.0) then
            Vtot(1)=Vtot(1)+jeff**2
         else if (orbprint(1).le.-1.and.orbprint(2).le.-1) then
            Vtot(2)=Vtot(2)+jeff**2
         else 
            Vtot(3)=Vtot(3)+jeff**2
         end if

      end do
   end do

   dipro%totjab(1)=sqrt(Vtot(1))*autoev
   dipro%totjab(2)=sqrt(Vtot(2))*autoev
   dipro%totjab(3)=sqrt(Vtot(3))*autoev
   
   write(*,'(A)') ".............................................................................."
   call ctx%message(":  total |J(AB,eff)| for hole transport (occ. MOs) :"//format_string(sqrt(Vtot(1))*autoev, '(f20.3)')//" eV  :")
   call ctx%message(":  total |J(AB,eff)| for charge transport (unocc. MOs) :"//format_string(sqrt(Vtot(2))*autoev,& 
          &'(f16.3)')//" eV  :")
   call ctx%message(":  total |J(AB,eff)| for charge transfer (CT) :"//format_string(sqrt(Vtot(3))*autoev, '(f25.3)')//" eV  :")
   write(*,'(A)') ".............................................................................."

   write(*,*) " "
   write(*,*) "Please remember, DIPRO is not available for restart!"
   write(*,*) "  "
   write(*,'(A)') "normal termination of dipro"
#else
   Call env%error("DIPRO is not available without TBLite.", source)
#endif

end subroutine get_jab

#if WITH_TBLITE
!> Unpack coefficients from a fragment orbital space in the full orbital space
subroutine unpack_coeff(full_bas, frag_bas, full, frag, mask)
   !> Basis set information for full system
   type(basis_type), intent(in) :: full_bas
   !> Basis set information for fragment
   type(basis_type), intent(in) :: frag_bas
   !> Quantity in full orbital space
   real(wp), intent(out) :: full(:)
   !> Quantity in fragment orbital space
   real(wp), intent(in) :: frag(:)
   !> Atom resolved mask for this fragment
   logical, intent(in) :: mask(:)

   integer :: iat, jat, ish, ii, jj, iao, jao, nao

   jat = 0
   do iat = 1, size(mask)
      if (mask(iat)) then
         jat = jat + 1
         ii = full_bas%ish_at(iat)
         jj = frag_bas%ish_at(jat)
         do ish = 1, full_bas%nsh_at(iat)
            iao = full_bas%iao_sh(ii+ish)
            jao = frag_bas%iao_sh(jj+ish)
            nao = full_bas%nao_sh(ii+ish)
            full(iao+1:iao+nao) = frag(jao+1:jao+nao)
         end do
      else
         ii = full_bas%ish_at(iat)
         do ish = 1, full_bas%nsh_at(iat)
            iao = full_bas%iao_sh(ii+ish)
            nao = full_bas%nao_sh(ii+ish)
            full(iao+1:iao+nao) = 0.0_wp
         end do
      end if
   end do
end subroutine unpack_coeff


!> Extract a the fragment structure from the full structure
subroutine get_structure_fragment(frag, struc, mask) !mol
   !> Molecular structure data of the fragment
   type(structure_type), intent(out) :: frag
   !> Molecular structure data of the full system
   type(structure_type), intent(in) :: struc  !mol
   !> Atom resolved mask for this fragment
   logical, intent(in) :: mask(:)

   integer :: nat
   integer, allocatable :: num(:)
   real(wp), allocatable :: xyz(:, :)

   nat = count(mask)
   num = pack(struc%num(struc%id), mask)
   xyz = reshape(pack(struc%xyz, spread(mask, 1, 3)), [3, nat])
   call new(frag, num, xyz)
end subroutine get_structure_fragment
#endif

end module xtb_dipro
