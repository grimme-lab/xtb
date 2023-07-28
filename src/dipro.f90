! This file is part of dipro.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

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
   use dipro_bondorder, only : get_wiberg_bondorder
   use dipro_fragment, only : get_wiberg_fragment
   use dipro_output, only : format_list, to_string
   use dipro_xtb, only : get_calculator
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
   implicit none
   private

   public :: get_jab, jab_input

   !> Configuration data for calculation
   type :: jab_input
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
      !> Electronic temperature in Kelvin
      real(wp) :: etemp = 300.0_wp
   end type jab_input

   !> Conversion factor from temperature to energy (Boltzmann's constant in atomic units)
   real(wp), parameter :: ktoau = 3.166808578545117e-06_wp

contains

!> Entry point for calculation of dipole projection related properties
subroutine get_jab(set, tblite, mol, fragment, error)
   use, intrinsic :: iso_fortran_env, only : output_unit
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol  !structure_type
   type(structure_type) :: struc

   !> Acc, Etemp, guess, chrg Input
   type(TSet), intent(in) :: set
   !requested gfn method for calculations Input
   type(TTBLiteInput), intent(inout) :: tblite
   !> Error handling
   type(error_type), allocatable, intent(out) :: error


   integer :: spin, charge, stat, unit, ifr, nfrag, nao, i, j, no, start_index,end_index, orbprint
   logical :: exist
   real(wp) :: energy, cutoff, jab, sab, jeff, ethresh
   real(wp), allocatable :: loc(:,:)
   type(context_type) :: ctx
   type(basis_type) :: bas
   type(xtb_calculator) :: xcalc, fcalc
   type(structure_type), allocatable :: mfrag(:)
   type(wavefunction_type) :: wfn
   type(wavefunction_type), allocatable :: wfx(:)
   real(wp), allocatable :: overlap(:, :), trans(:, :), wbo(:, :), chrg(:), p2mat(:,:), coeff2(:,:)
   real(wp), allocatable :: orbital(:, :, :), scmat(:, :), fdim(:, :), scratch(:), efrag(:),y(:,:),Edim(:,:)
   integer, allocatable :: fragment(:), spinfrag(:)
   !> Molecular gradient
   real(wp), allocatable :: gradient(:, :)
   !> Strain derivatives
   real(wp), allocatable :: sigma(:, :)

   struc=mol

   if ( tblite%method == '' ) then 
      tblite%method = 'gfn2'
      write(*,*) "==> No method provided, falling back to default GFN2-xTB."
   end if

   call get_calculator(xcalc, struc, tblite%method, error)  !mol
!   if (allocated(error)) return
   call new_wavefunction(wfn, struc%nat, xcalc%bas%nsh, xcalc%bas%nao, &   !mol%nat
      & 1, set%etemp * ktoau)
   wfn%nspin=1 !XXXX  das ist number of spins, nicht spin S, warum hat das überhaupt eine dimension, xtb kann doch gar nicht mehrere spinkanäle
               !parallel speichern und rechnen wie zb singlet triplet dublet

   call ctx%message("Calculation for dimer ")
   call ctx%message("charge of dimer : "//format_string(mol%chrg, '(f7.0)'))
   write(*,*) "unpaired e- of dimer : ", set%nalphabeta

   call xtb_singlepoint(ctx, struc, xcalc, wfn, tblite%accuracy, energy,gradient,sigma,2) !, &  !mol
 !     & verbosity-1) !input%verbosity-1
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if

   allocate(overlap(xcalc%bas%nao, xcalc%bas%nao),y(xcalc%bas%nao,2))
   cutoff = get_cutoff(xcalc%bas)
   call get_lattice_points(struc%periodic, struc%lattice, cutoff, trans)  !mol,mol
   call get_overlap(struc, trans, cutoff, xcalc%bas, overlap)  !mol

   if (allocated(fragment)) then
   else
      allocate(wbo(struc%nat, struc%nat),p2mat(xcalc%bas%nao,xcalc%bas%nao)) !mol,mol
      do i=1,xcalc%bas%nao
         do j=1,xcalc%bas%nao
            p2mat(i,j)=wfn%density(i,j,1)
         end do
      end do
      call get_wiberg_bondorder(xcalc%bas, overlap, p2mat, wbo) !wfn%density hat [nao,nao,spin], pmat in wiberg bondorder nur [nao,nao]

      allocate(fragment(struc%nat)) !mol
      call get_wiberg_fragment(fragment, wbo, 0.1_wp)
   end if

   nfrag = maxval(fragment)
   select case(nfrag)
   case(:1)
      call fatal_error(error, "Found no fragments in the input structure")
      return
   case(2:)
      call ctx%message("Found "//to_string(nfrag)//" fragments:")
      do ifr = 1, nfrag
         call ctx%message("Fragment "//to_string(ifr)//":  "//format_list(fragment == ifr))
      end do
      call ctx%message("")
      if (nfrag > 2) then
         call fatal_error(error, "Found "//to_string(nfrag)// &
            & " fragments in the input structure, too many fragments.")
         return
      end if
   end select

   nao = size(wfn%emo)
   allocate(orbital(nao, nfrag,nao), efrag(nfrag), scmat(nao, nao), fdim(nao, nao), &
      & scratch(nao), mfrag(nfrag), wfx(nfrag), chrg(nfrag), spinfrag(nfrag))

   !--------------------------------------------------------------------------------------

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
      stop
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

   do ifr = 1, nfrag 
      call ctx%message("Calculation for fragment "//to_string(ifr))
      call get_structure_fragment(mfrag(ifr), struc, fragment == ifr) !mol

      !------------------summation of fragment charges stored in chrg(nfrag)-------
      if ( all(chrg.eq.0) ) then
         chrg(ifr)=0
         !> fragment mask generated on the fly 
         chrg(ifr)=sum(pack(wfn%qat(:,1), fragment == ifr)) !wfn%qat ist [nat,nspin=1]
         mfrag(ifr)%charge=nint(chrg(ifr))
      else
         mfrag(ifr)%charge=nint(chrg(ifr))
      end if
      call ctx%message("charge of fragment : "//format_string(mfrag(ifr)%charge , '(f7.0)'))
      !---------------------uhf fragments spins------------------------------------
      mfrag(ifr)%uhf = spinfrag(ifr)
!      call ctx%message("unpaired e- of fragment : "//format_string(mfrag(ifr)%uhf , '(f7.0)'))
      write(*,*) "unpaired e- of fragment : ", mfrag(ifr)%uhf
      !----------------------------------------------------------------------------

      call get_calculator(fcalc, mfrag(ifr), tblite%method, error)
      !> mol%charge is updated automatically from wfn by tblite library 
      if (allocated(error)) return
         call new_wavefunction(wfx(ifr), mfrag(ifr)%nat, fcalc%bas%nsh, fcalc%bas%nao, &
         & 1, set%etemp * ktoau)

         !> mol%type (dimer) == mfrag%type (fragments), wfn (dimer) == wfx (fragments), calc (dimer)==fcalc(fragments)
         wfx%nspin=1 !XXXX
         call xtb_singlepoint(ctx, mfrag(ifr), fcalc, wfx(ifr), tblite%accuracy, energy) !, &
!         & verbosity-1)  !=input%verbosity-1)
         if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if

      no = wfx(ifr)%homo(max(2,1))  !homo hat dimension [nao]
!---------------------------------------------------------
!write(*,*) "emo von wfx homo", wfx(ifr)%emo(no,1)  !geht
!---------------------------------------------------------

      do j = 1, nao
         !> for homo,homo-1,lumo,lumo+
!         no = wfx(ifr)%homo(max(2,1))+(j-2)
         !> coeff(homo) stored in orbital(:,ifr,orb)
         call unpack_coeff(xcalc%bas, fcalc%bas, orbital(:, ifr, j), &
         & wfx(ifr)%coeff(:, j,1), fragment == ifr)  !coeff hat dimension [nao,nao,spin]
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

   ethresh=0.1d0
   call ctx%message("energy threshhold for near-degenerate orbitals near HOMO and LUMO &
           &considered for DIPRO: "//format_string(ethresh, '(f20.3)')//" eV")
   do ifr=1,nfrag
      do j = 1, nao
         if (wfx(ifr)%emo(j,1) .ge. (wfx(ifr)%emo(wfx(ifr)%homo(max(2,1)),1) - ethresh/autoev) .and.&
           & wfx(ifr)%emo(j,1) .le. (wfx(ifr)%emo(wfx(ifr)%homo(max(2,1))+1,1) + ethresh/autoev)) then
              if (start_index.eq.-1) then 
                 start_index = j
              end if
              end_index = j
         endif
      end do
      call ctx%message("fragment no. "//format_string(ifr, '(i0)'))
      call ctx%message("no. of orbitals considered within energy window: "//format_string(end_index-start_index+1, '(i4)'))
   end do

   !> gemm(amat,bmat,cmat,transa,transb,a1,a2): X=a1*Amat*Bmat+a2*Cmat
   call gemm(overlap, coeff2, scmat)  !scmat=S_dim*C_dim
   do j = start_index, end_index !1, nao
   y=0
      do ifr = 1, nfrag
      !> gemv(amat, xvec,yvec,a1,a2,transa): X=a1*Amat*xvec+a2*yvec
         call gemv( scmat, orbital(:, ifr, j), y(:,ifr), trans="t" ) !y_mon(ifr)=C_mon(j)*S_dim*C_dim
         call gemv( Edim, y(:,ifr), scratch ) !scratch=E_dim*y(ifr)
         efrag(ifr)=dot( y(:,ifr), scratch) !E_mon=y(ifr)*E_dim*y(ifr)
      end do

      sab=dot( y(:,1), y(:,2) )
      call gemv( Edim, y(:,2), scratch ) !scratch=E_dim*y(2)
      jab=dot( y(:,1), scratch ) !jab=y(1)*E_dim*y(2)

   jeff = (jab - sum(efrag) / nfrag * sab) / (1.0_wp - sab**2)

   orbprint=wfx(1)%homo(max(2,1))-j

   if (orbprint.gt.0) then
      call ctx%message("HOMO - "//format_string(orbprint, '(i0)')//" orbital") 
      else if (orbprint.eq.0) then
      call ctx%message("HOMO orbital")
      else if (orbprint.eq.-1) then
      call ctx%message("LUMO orbital")
      else if (orbprint.lt.-1) then
      call ctx%message("LUMO + "//format_string(abs(orbprint)-1, '(i0)')//" orbital")
   end if

   call ctx%message("E_mon(orb) frag1 frag2"//format_string(efrag(1)*autoev, '(f20.3)')// &
           &format_string(efrag(2)*autoev, '(f20.3)')//" eV")
   !> abs = absolute values
   call ctx%message("|J(AB)|: "//format_string(abs(jab)*autoev, '(f20.3)')//" eV")    
!   call ctx%message("S(AB): "//format_string(sab, '(f20.8)'))
   call ctx%message("|J(AB,eff)|: "//format_string(abs(jeff)*autoev, '(f16.3)')//" eV")
   end do

end subroutine get_jab


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
!>
!> Todo: This routine can currently only create neutral fragments,
!>       charge and spin information is not partition or transferred to fragments
!> Done: Fragment charges added in get_jab routine after get_structure_fragment call
!>       Fragment spin separate read-in in get_jab routine after fragment charges
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
   num = pack(struc%num(struc%id), mask)  !mol,mol
   xyz = reshape(pack(struc%xyz, spread(mask, 1, 3)), [3, nat]) !mol
   call new(frag, num, xyz)
end subroutine get_structure_fragment


end module xtb_dipro
