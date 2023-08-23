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
! You should have received a copy of the GNU Lesser General Public Licen
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> Implementaion of the ONIOM method
!> publication: https://doi.org/10.1039/D3CP02178E (further reference) 
module xtb_oniom

   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_convert, only: aatoau
   use xtb_type_atomlist, only: TAtomlist, len
   use xtb_type_calculator, only: TCalculator
   use xtb_type_data, only: scc_results
   use xtb_type_environment, only: TEnvironment
   use xtb_type_molecule, only: TMolecule, len
   use xtb_type_restart, only: TRestart
   use xtb_type_neighbourlist, only: TNeighbourlist
   use xtb_gfnff_calculator, only: TGFFCalculator, newGFFCalculator
   use xtb_type_topology, only: TTopology, len
   use xtb_xtb_calculator, only: TxTBCalculator, newXTBCalculator, newWavefunction
   use xtb_extern_orca, only: TOrcaCalculator, newOrcaCalculator
   use xtb_extern_turbomole, only: TTMCalculator, newTMCalculator
   use xtb_setparam, only: set

   implicit none
   private

   public :: TOniomCalculator, newOniomCalculator, oniom_input, calculateCharge
   
   !> To handle CML arguments
   type :: oniom_input
      character(len=:), allocatable :: first_arg
      character(len=:), allocatable :: second_arg
   
      end type oniom_input
   
   !> ONIOM calculator 
   type, extends(TCalculator) :: TOniomCalculator 

      !> methods specified by the user, or default: gfnff, gfn2
      integer :: method_low, method_high

      !> charge for inner region
      integer :: chrg_model
         ! rewrite

      !> inner region list
      type(TAtomList) :: list 

      !> whole region, low-level 
      class(TCalculator), allocatable :: real_low

      !> inner region, low-level
      class(TCalculator), allocatable :: model_low

      !> inner region, high-level
      class(TCalculator), allocatable :: model_high

      !> index list of the inner region
      integer, allocatable :: idx(:)
      
      !> topology used for cutting bonds
      type(TTopology), allocatable :: topo

      !> wavefunctions for inner region calculations
      type(TRestart) :: chk_low, chk_high

   contains

      procedure :: singlepoint
      procedure :: hessian
      procedure :: writeInfo
      procedure :: cutbond

   end type TOniomCalculator

   interface coord
      module procedure :: newcoord
   end interface

   interface resize
      module procedure :: new_atom
      module procedure :: new_coordinates
   end interface


contains

!> create ONIOM Calcultor
subroutine newOniomCalculator(self, env, mol, input)
   
   implicit none
   
   !> new ONIOM calculator
   type(TOniomCalculator), intent(out) :: self
   
   !> calculation environment
   type(TEnvironment), intent(inout) :: env
   
   !> molecular structure data
   type(TMolecule), intent(in) :: mol
   
   !> cml input
   type(oniom_input), intent(in) :: input
   
   type(TxTBCalculator), allocatable :: xtb
   type(TOrcaCalculator), allocatable :: orca
   type(TGFFCalculator), allocatable :: gff
   integer :: icol
   integer :: i
   
   !--------------!
   ! method setup !
   !--------------!

   ! methods explicitly specified !
   if (allocated(input%first_arg)) then
      
      icol = index(input%first_arg, ':')
      if (icol == 0) then
         call env%error("Invalid method '"//input%first_arg//"' provided")
         return
      end if

      self%method_high = string_to_id(input%first_arg(:icol - 1))
      self%method_low = string_to_id(input%first_arg(icol + 1:))
      
      if (self%method_high < 0 .or. self%method_high > 5) then
         call env%error("Invalid high-level method")
         return
      end if
      
      if (self%method_low < 0 .or. self%method_low > 3) then
         call env%error("Invalid low-level method")
         return
      end if
   
   ! default, gfn2:gfnff !
   else

      self%method_high = 2
      self%method_low = 3
   
   endif
   
   ! write user-defined inner region list as raw string into array !
   self%list = TAtomList(list=input%second_arg)
   call self%list%to_list(self%idx)

   if (len(self%list) .eq. 0 .or. self%list%error) then
      call env%error("Invalid inner region: '"//input%second_arg//"'")
      return
   end if

   ! check inner region !
   if (any(self%idx < 1) .or. any(self%idx > mol%n)) then
      call env%error('Out of bound inner region')
      return
   end if

   ! whole region calculator ! 
   select case (self%method_low)
   
   ! gfn1/2 !
   case default
      allocate (xtb)
      call newXTBCalculator(env, mol, xtb, method=self%method_low)
      call move_alloc(xtb, self%real_low)

   ! gfnff !
   case (3)
      allocate (gff)
      call newGFFCalculator(env, mol, gff, "", .false.)
      call move_alloc(gff, self%real_low)

   end select

end subroutine newOniomCalculator

!> 3 singlepoint energy calculations
subroutine singlepoint(self, env, mol, chk, printlevel, restart, energy, gradient, sigma, hlgap, results)

   use xtb_io_writer
   use xtb_mctc_filetypes, only : fileType
   use xtb_xtb_calculator, only: TxTBCalculator

   implicit none
   
   !> instance of TOniomCalculator
   class(TOniomCalculator), intent(inout) :: self
   
   !> calculation environment      
   type(TEnvironment), intent(inout) :: env
   
   !> molecular structure data
   type(TMolecule), intent(inout) :: mol
   
   !> restart data wrapper
   type(TRestart), intent(inout) :: chk
   
   !> print level for IO
   integer, intent(in) :: printlevel 
   
   !> if restarted
   logical, intent(in) :: restart
   
   
   !> ONIOM total energy
   real(wp), intent(out) :: energy
   
   !> ONIOM gradients
   real(wp), intent(out) :: gradient(:, :)
   
   !> strain derivatives
   real(wp), intent(out) :: sigma(:, :) 
   
   !> HOMO-LUMO gap
   real(wp), intent(out) :: hlgap
   
   !> final output stream
   type(scc_results), intent(out) :: results
   
   
   ! temporary storages !
   type(TxTBCalculator), allocatable :: tmp 
   type(TGFFCalculator), allocatable :: gff
   type(TOrcaCalculator), allocatable :: orca
   type(TTMCalculator), allocatable :: turbo
   
   ! for inner region high- and low-level calculation !
   type(scc_results) :: results_low, results_high
   type(TMolecule) :: inner_mol
   real(wp), allocatable :: gradient_low(:, :), gradient_high(:, :)
   real(wp) :: energy_model_low, energy_model_high, hlgap_low, hlgap_high
   real(wp) :: sigma_low(3, 3), sigma_high(3, 3)

   !> arrays for matmul of gradient matrix and Jacobian  
   real(wp),allocatable :: arr_gh(:), arr_gl(:)
   
   !> Jacobian matrix
   real(wp), allocatable :: jacobian(:,:)
   
   integer,allocatable :: idx2(:)
   integer :: i, coord_unit
   logical :: exitRun

   ! check whether the calculator is initialized !
   if (.not. allocated(self%real_low)) then
      call env%error("Outer region calculator not provided")
      return
   end if

   ! forward solvation to outer region !
   if (allocated(self%solvation)) then
      call move_alloc(self%solvation, self%real_low%solvation)
   end if
   
   !-------------------------!
   ! whole system, low-level !
   !-------------------------!

   if (.not.set%oniom_settings%cut_inner) then
      
      if (printlevel > 0) then
         write(env%unit,'(a)')
         write(env%unit,'(2x,72("-"))')
         write (env%unit, '(/6x,a/)') "Singlepoint calculation of whole system with low-level method"
         write(env%unit,'(2x,72("-"))')
         write(env%unit,'(a)')
      endif

      call self%real_low%singlepoint(env, mol, chk, printlevel, restart, &
         & energy, gradient, sigma, hlgap, results)
   
   endif
   
   ! creating Linked Atoms !
   call self%cutbond(env, mol, chk, self%topo, inner_mol,jacobian,idx2)
   
   ! define inner region charge !
   inner_mol%chrg = real(set%oniom_settings%innerchrg)

   ! --cut flag termination !
   if (set%oniom_settings%cut_inner) then
      
      write(env%unit,'(a)')
      write(env%unit,'(2x,72("-"))')
      write(env%unit,'(2x,"|",24x,a,1x,i0,22x,"|")') "INNER REGION CHARGE = ", nint(inner_mol%chrg)
      write(env%unit,'(2x,72("-"))')
      write(env%unit,'(a)')
      call terminate(0)
   
   endif
    
   call env%check(exitRun)
   if (exitRun) then
      call env%error("Could not create the linked atoms")
      return
   end if

   !-------------------------!
   ! inner region, low-level !
   !-------------------------!
   
   if (.not. allocated(self%model_low)) then

      select case (self%method_low)
      case default
         call env%error("Invalid low-level inner method")
         return

      ! gfn1/2 !
      case (1, 2)
         allocate (tmp)
         call newXTBCalculator(env, inner_mol, tmp, method=self%method_low)
         call move_alloc(tmp, self%model_low)

      ! gfnff !
      case (3)
         allocate (gff)
         call newGFFCalculator(env, inner_mol, gff, "", .false.)
         call move_alloc(gff, self%model_low)
      
      end select
      
      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not setup low-level method")
         return
      end if
   
   end if
  
   !--------------------------!
   ! inner region, high-level !
   !--------------------------!

   if (.not. allocated(self%model_high)) then

      select case (self%method_high)
      case default
         call env%error("Invalid high-level inner method")
         return
      
      ! gfn1/2 !
      case (1, 2)
         allocate (tmp)
         call newXTBCalculator(env, inner_mol, tmp, method=self%method_high)
         call move_alloc(tmp, self%model_high)

      ! gfnff !
      case (3)
         allocate (gff)
         call newGFFCalculator(env, inner_mol, gff, "", .false.)
         call move_alloc(gff, self%model_high)

      ! orca !
      case (4)
         allocate (orca)
         call newOrcaCalculator(orca, env, set%ext_orca,oniom=.true.)
         call move_alloc(orca, self%model_high)

      ! turbomole !
      case (5)  
         allocate (turbo)
         call newTMCalculator(turbo, 1, 1)
         call move_alloc(turbo, self%model_high)
         
         ! to copy coord file into origin.coord !
         call protectCoord(env)

      end select
      
      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not setup high-level method")
         return
      end if

   end if
 
   ! setup partial and shell charges for gnf1/2 !
   if (.not.allocated(self%chk_low%wfn%qsh)) then
      select type (calc => self%model_low)
      type is (TxTBCalculator)
         call newWavefunction(env, inner_mol, calc, self%chk_low)
      end select
   end if
   if (.not.allocated(self%chk_high%wfn%qsh)) then
      select type (calc => self%model_high)
      type is (TxTBCalculator)
         call newWavefunction(env, inner_mol, calc, self%chk_high)
      end select
   end if
   

   ! SP: inner region low-level !
   if (printlevel > 0) then
      write(env%unit,'(a)')
      write(env%unit,'(2x,72("-"))')
      write (env%unit, '(/6x,a/)') "Singlepoint calculation of inner region with low-level method"
      write(env%unit,'(2x,72("-"))')
      write(env%unit,'(a)')
   endif

   allocate (gradient_low(3, inner_mol%n))
   energy_model_low = 0.0_wp
   gradient_low = 0.0_wp

   call self%model_low%singlepoint(env, inner_mol, self%chk_low, printlevel, restart, &
       & energy_model_low, gradient_low, sigma_low, hlgap_low, results_low)


   ! SP: inner region high-level !
   if (printlevel > 0) then
      write(env%unit,'(a)')
      write(env%unit,'(2x,72("-"))')
      write (env%unit, '(/6x,a/)') "Singlepoint calculation of inner region with high-level method"
      write(env%unit,'(2x,72("-"))')
      write(env%unit,'(a)')
   endif

   allocate (gradient_high(3, inner_mol%n))
   energy_model_high = 0.0_wp
   gradient_high = 0.0_wp

   call self%model_high%singlepoint(env, inner_mol, self%chk_high, printlevel, restart, &
       & energy_model_high, gradient_high, sigma_high, hlgap_high, results_high)
   

   ! write opt logs for inner region !
   if (set%oniom_settings%logs)then
      call writeMolecule(inner_mol, set%oniom_settings%ilog1, format=filetype%xyz,energy=energy_model_low)
      call writeMolecule(inner_mol, set%oniom_settings%ilog2, format=filetype%xyz,energy=energy_model_high)
   endif


   results%dipole = results%dipole + results_high%dipole - results_low%dipole 
   sigma = sigma - sigma_low + sigma_high
   results%hl_gap = hlgap - hlgap_low + hlgap_high

   !----------------!
   ! Postprocessing !
   !----------------!

   ! ONIOM energy / reference formula (1,2) !
   energy = energy + energy_model_high - energy_model_low
   results%e_total = energy
   
   ! [gradient*Jacobian] with forward and backward transformation !
   call matrix_to_array(gradient_high,arr_gh)
   call matrix_to_array(gradient_low,arr_gl)
   
   arr_gh=matmul(jacobian,arr_gh)
   arr_gl=matmul(jacobian,arr_gl)
   
   call array_to_matrix(arr_gh,gradient_high)
   call array_to_matrix(arr_gl,gradient_low)
   

   ! ONIOM gradients !
   do i = 1, size(idx2)
      gradient(:, idx2(i)) = gradient(:, idx2(i)) + gradient_high(:, i) - gradient_low(:, i) 
   end do 
   results%gnorm=norm2(gradient)   
   deallocate(gradient_high)
   deallocate(gradient_low)

end subroutine singlepoint

!> Evaluate hessian by finite difference for all atoms
subroutine hessian(self, env, mol0, chk0, list, step, hess, dipgrad)

   character(len=*), parameter :: source = "extern_oniom_hessian"
   
   !> single point calculator
   class(TOniomCalculator), intent(inout) :: self
   
   !> computation environment
   type(TEnvironment), intent(inout) :: env
   
   !> molecular structure data
   type(TMolecule), intent(in) :: mol0
   
   !> restart data
   type(TRestart), intent(in) :: chk0
   
   !> list of atoms to displace
   integer, intent(in) :: list(:)
   
   !> step size for numerical differentiation
   real(wp), intent(in) :: step
   
   !> array to add Hessian to
   real(wp), intent(inout) :: hess(:, :)
   
   !> array to add dipole gradient to
   real(wp), intent(inout) :: dipgrad(:, :)

   real(wp), allocatable :: jacobian(:,:)
   integer,allocatable :: idx2(:)

   integer :: ii, jj, ic, jc, im, jm
   real(wp), allocatable :: hess_model(:, :), dipgrad_model(:, :)
   type(TMolecule) :: mol_model
   integer, allocatable :: list_model(:)

   ! compute complete hessian for outer region !
   call self%real_low%hessian(env, mol0, chk0, list, step, &
       & hess, dipgrad)

   !call chk0%wfn%allocate(mol0%n,chk0%basis%nshell,chk0%real_low%basis%nao)
   
   ! create linked atoms !
   call self%cutbond(env, mol0, chk0, self%topo, mol_model,jacobian,idx2)
   mol_model%chrg = float(self%chrg_model)
   list_model = [(ii, ii = 1, size(self%idx))]

   allocate(hess_model(3*mol_model%n, 3*mol_model%n))
   allocate(dipgrad_model(3, 3*mol_model%n))

   hess_model(:, :) = 0.0_wp
   dipgrad_model(:, :) = 0.0_wp

   call self%model_low%hessian(env, mol_model, self%chk_low, list_model, step, &
       & hess_model, dipgrad_model)

   do ii = 1, size(self%idx)
      im = 3 * (ii - 1)
      ic = 3 * (self%idx(ii) - 1)
      do jj = 1, size(self%idx)
         jm = 3 * (jj - 1)
         jc = 3 * (self%idx(jj) - 1)
         hess(jc+1:jc+3, ic+1:ic+3) = hess(jc+1:jc+3, ic+1:ic+3) &
            & - hess_model(jm+1:jm+3, im+1:im+3)
      end do
      dipgrad(:, ic+1:ic+3) = dipgrad(:, ic+1:ic+3) - dipgrad_model(:, im+1:im+3)
   end do

   hess_model(:, :) = 0.0_wp
   dipgrad_model(:, :) = 0.0_wp
   call self%model_high%hessian(env, mol_model, self%chk_high, list_model, step, &
       & hess_model, dipgrad_model)

   do ii = 1, size(self%idx)
      im = 3 * (ii - 1)
      ic = 3 * (self%idx(ii) - 1)
      do jj = 1, size(self%idx)
         jm = 3 * (jj - 1)
         jc = 3 * (self%idx(jj) - 1)
         hess(jc+1:jc+3, ic+1:ic+3) = hess(jc+1:jc+3, ic+1:ic+3) &
            & + hess_model(jm+1:jm+3, im+1:im+3)
      end do
      dipgrad(:, ic+1:ic+3) = dipgrad(:, ic+1:ic+3) + dipgrad_model(:, im+1:im+3)
   end do

end subroutine hessian

subroutine writeInfo(self, unit, mol)

   class(TOniomCalculator), intent(in) :: self

   !> unit for I/O
   integer, intent(in) :: unit
   
   !> molecular structural data
   type(Tmolecule), intent(in) :: mol

end subroutine writeInfo

!> create inner region
subroutine cutbond(self, env, mol, chk, topo, inner_mol, jacobian, idx2)

   use xtb_type_molecule, only: init
   use xtb_topology, only: makeBondTopology, topologyToNeighbourList
   use xtb_io_writer
   use xtb_mctc_filetypes, only : fileType
   
   implicit none
   character(len=*), parameter :: source = "xtb_oniom_cutbond" 

   !> polymorphic calculator
   class(TOniomCalculator), intent(in) :: self
   
   !> calculation environment
   type(TEnvironment), intent(inout) :: env
   
   !> wavefunction wrapper
   type(TRestart), intent(in) :: chk
   
   !> molecular structure data 
   type(TMolecule), intent(in) :: mol
   
   !> inner region mol str data
   type(TMolecule), intent(inout) :: inner_mol
   
   !> topology info
   type(TTopology), allocatable, intent(inout) :: topo
   
   !> jacobian matrix 
   real(wp), allocatable, intent(inout) :: jacobian(:,:)
   
   !> list of inner region atom indices + host atoms
   integer, allocatable, intent(out) :: idx2(:)
  
   
   !> outer region geometry
   type(TMolecule), allocatable :: outer_mol

   !> neighbour list
   type(TNeighbourList) :: neighList
      
   !> pair-wise indices for cutbound  
   integer, allocatable :: brokenBondPairs(:,:)
   
   integer, allocatable :: at(:), at_out(:)
   integer, allocatable :: bonded(:, :)
   real(wp), allocatable :: xyz(:, :), xyz_out(:, :)
   character(len=:),allocatable :: fname_inner
   
   !> if the both bonded atoms are inside the inner region
   logical :: inside
   
   !> control output
   logical :: set1 = .true.
   logical :: set2 = .true.

   !> number of LAs
   integer :: nla

   integer :: i, j, k, pre_last, pre_last_out, iterator
   integer :: io

   !> initial no. atoms in inner & outer regions
   integer :: in, out

   !> loop indices
   integer :: in_itr, out_itr

   !-------!
   ! SETUP !
   !-------!

   inside = .FALSE.
   
   ! initial number of atoms in inner/outer region without LAs !
   in = len(self%list)
   out = mol%n - in 
   in_itr = 1
   out_itr = 1
   nla = 0

   ! allocate accordingly the basic molecular data !
   allocate (at(in))
   allocate (xyz(3, in))
   allocate (at_out(out))
   allocate (xyz_out(3, out))

   ! save inner region list in the matching array !
   idx2=self%idx
   
   ! divide initial mol into inner and outer regions !
   do i = 1, mol%n
      if (i==self%idx(in_itr)) then
         
         if (in_itr > in) then
            call env%error("The internal error, inconsistent molecular dimensionality", source=source) 
            return
         endif
         
         at(in_itr) = mol%at(i)
         xyz(:, in_itr) = mol%xyz(:, i)
         in_itr = in_itr + 1
      
      else
         
         if (out_itr > out) then 
            call env%error("The internal error, inconsistent molecular dimensionality", source=source) 
            return
         endif
         
         at_out(out_itr) = mol%at(i)
         xyz_out(:, out_itr) = mol%xyz(:, i)
         out_itr = out_itr + 1
      
      endif
   enddo 
 
   ! initialiaze jacobian as identity !
   call create_jacobian(jacobian,at)
   
   ! identify bonded atoms and save them into an array + assign iterator !
   select type (calc => self%real_low)
   class default
      call env%error("Topology information could not be derived from the given calculator",source)
      return

   ! gfnff !
   type is (TGFFCalculator)
      
      ! bonded atom list !
      bonded = calc%topo%blist 
      
      ! number of bonds !
      iterator = size(bonded,2)

   ! gfn1/2 !
   type is (TxTBCalculator)
      if (.not. allocated(topo)) then
         allocate (topo)

         ! return assigned topo%list !
         call makeBondTopology(topo, mol, chk%wfn%wbo) 
         
         ! return neighList !
         call topologyToNeighbourList(topo, neighList, mol) 
      
      end if

      allocate (bonded(2, len(topo)))
      do i = 1, len(topo)
         bonded(:, i) = topo%list(1:2, i)
      end do
      
      iterator = size(bonded,2)
   
   end select

   !-----------!
   ! ALGORITHM !
   !-----------!

   ! iterate for all atoms in the user-provided list !
   do i = 1, size(self%idx)

      ! iterate through pair of atoms that are bonded !
      do j = 1, iterator

         ! if atom in the list is bonded !
         if (bonded(1, j) == self%idx(i)) then

            ! iterate again through the list !
            do k = 1, size(self%idx)
               
               ! inside inner region !
               if (self%idx(k) == bonded(2, j)) then
                  inside = .TRUE.
               end if

            end do

            ! bond is broken !
            if (.not. inside) then
               
               if (.not.set%oniom_settings%ignore_topo) then

                  ! check if single bond is broken  !
                  select type (calc => self%real_low)
                  class default
                     call checkfororder(env, mol, self%idx(i), bonded(2, j), bond=topo%list(3, j))
                  type is (TGFFCalculator)
                     call checkfororder(env, mol, self%idx(i), bonded(2, j), hybrid=calc%topo%hyb)
                  end select
               
               endif

               ! adjust ordinal numbers !
               call resize(at)
               call resize(idx2)
               call resize(at_out)
         
               ! save index of the host atom -> jacobian ! 
               idx2(size(idx2)) = bonded(2,j)

               ! assign new atom as H !
               at(size(at)) = 1
               at_out(size(at_out)) = 1
               
               ! number of LAs !
               nla = nla + 1
               
               ! adjust coordinate matrix !
               call resize(xyz)
               call resize(xyz_out)
      
               ! increment Jacobian ! 
               call resize_jacobian(jacobian)   
               
               ! determine new position of added H atom !
               call coord(env,mol,xyz,xyz_out,self%idx(i),bonded(2,j),jacobian,i)
               
            end if
            
            inside = .FALSE.
         
         ! if atom in the list is bonded !
         else if (bonded(2, j) == self%idx(i)) then
            
            ! iterate again through the list !
            do k = 1, size(self%idx)
                  
               ! inside inner region !
               if (self%idx(k) == bonded(1, j)) then
                  inside = .TRUE.
               end if

            end do

            ! bond is broken !
            if (.not. inside) then
               
               if (.not.set%oniom_settings%ignore_topo) then
                  
                  ! check if single bond ! 
                  select type (calc => self%real_low)   
                  class default
                     call checkfororder(env, mol, self%idx(i), bonded(1, j), bond=topo%list(3, j))
                  type is (TGFFCalculator)
                     call checkfororder(env, mol, self%idx(i), bonded(1, j), hybrid=calc%topo%hyb)
                  end select
               
               endif

               ! adjust ordinal numbers !
               call resize(at)
               call resize(at_out)
               call resize(idx2)
               
               ! save index of the host atom -> jacobian ! 
               idx2(size(idx2)) = bonded(1,j)
               
               ! assign new atom as H !
               at(size(at)) = 1
               at_out(size(at_out)) = 1
               
               ! number of LAs !
               nla = nla + 1
               
               ! adjust coordinate matrix !
               call resize(xyz)
               call resize(xyz_out)
               
               ! increment Jacobian matrix !
               call resize_jacobian(jacobian)   
               
               ! determine new position of added H atom !
               call coord(env,mol, xyz, xyz_out, self%idx(i), bonded(1, j),jacobian,i)
            
            end if

            inside = .FALSE.
         
         end if
      end do
   end do

   !----------------!
   ! POSTPROCESSING !
   !----------------!
   
   ! outer region saturation !
   if (set%oniom_settings%outer) then
      
      ! initialize !
      allocate(outer_mol)
      call init(outer_mol, at_out, xyz_out)

      ! create Xmol file !
      call open_file(io, "outer_region.xyz", "w")
      call writeMolecule(outer_mol, io, filetype%xyz)
      call close_file(io)
   
      ! check LAs postions !
      if (set2) call check_dist(outer_mol,env,nla)
   
   end if
   
   ! initialize inner region mol !
   call init(inner_mol, at, xyz)
   
   ! check LAs postions !
   if (set1) call check_dist(inner_mol,env,nla)
  
   ! to distinguish cases with/without SP !
   if (set%oniom_settings%cut_inner) then
      fname_inner = "inner_region_without_h.xyz" 
   else
      fname_inner = "inner_region.xyz"
   endif

   ! create Xmol file for inner region !
   call open_file(io, fname_inner, "w")
   call writeMolecule(inner_mol, io, filetype%xyz)
   call close_file(io) 
   
   call env%checkpoint("ONIOM is terminated")

   set1 = .false.
   set2 = .false.  

end subroutine cutbond

subroutine check_dist(mol, env, nla)

   use xtb_mctc_convert, only : aatoau

   !> source
   character(len=*), parameter :: source = "check_dist"

   !> molecular structure
   type(TMolecule), intent(in) :: mol 

   !> calculation environment 
   type(TEnvironment), intent(inout) :: env

   !> LA's  indices
   integer, intent(in) :: nla


   !> warning message
   character(len=:), allocatable :: warn
   
   !> ordinal numbers as characters 
   character(len=3) :: dummy1, dummy2

   !> minimal allowed distance
   real(wp) :: min_dist

   !> loop indices
   integer :: i,j

   !> inner/outer region with LA's
   integer :: regsize_saturated

   !> inner/outer region without LA's 
   integer :: regsize_raw

   ! calculate the start and end iterator for loop !
   regsize_raw = mol%n - nla + 1
   regsize_saturated = mol%n

   ! H2 bond length !
   min_dist = 0.74 * aatoau   

   ! check all LA's !
   do i = regsize_raw, regsize_saturated
      do j = 1, regsize_saturated
         if (i.ne.j) then
            if (mol%dist(j,i) < min_dist) then
               write (dummy1, '(I3)') i
               write (dummy2, '(I3)') j
               warn = "The distance b/n atoms "//dummy1//" (LA) and " &
                     & //dummy2//" is less then min," &
                     & //achar(10) //"please examine carefully your cutout region"  
               call env%warning(warn,source)
            endif
         endif
      enddo
   enddo

end subroutine check_dist

!> increase atomic number array by 1 
subroutine new_atom(at)
   
   implicit none
   integer, allocatable, intent(inout) :: at(:)
   integer, allocatable :: tmp2(:)

   if (.not. allocated(at)) then
      allocate(at(1))
   else
      allocate (tmp2(size(at) + 1))
      tmp2(:size(at)) = at(:size(at))
      deallocate (at)
      call move_alloc(tmp2, at)
   endif

end subroutine new_atom

!> increase coordinate matrix  by 1 
subroutine new_coordinates(xyz)
   
   implicit none
   real, allocatable :: xyz(:,:)
   real, allocatable :: tmp1(:,:)
   integer :: atom_num

   atom_num = size(xyz, 2)
   allocate (tmp1(size(xyz,1), atom_num + 1))
   tmp1(:, :size(xyz,2)) = xyz(:, :size(xyz,2))
   deallocate (xyz)
   call move_alloc(tmp1, xyz)

end subroutine new_coordinates

!> create identity matrix 
subroutine create_jacobian(matrix,at)
   
   implicit none
   real(wp),allocatable :: matrix(:,:)
   integer,intent(in) :: at(:)
   integer :: i, j

   allocate (matrix(size(at)*3,size(at)*3))
   do i=1,size(at)*3
      do j=1,size(at)*3
         if (i==j) then
            matrix(i,j)=1.0_wp
         else
            matrix(i,j)=0.0_wp
         endif
      enddo
   enddo
   
end subroutine create_jacobian

!> increase matrix dimensionality
! (3 new entries in diagonal and subsequent increase) !
subroutine resize_jacobian(matrix)
   
   implicit none
   real(wp), allocatable :: matrix(:,:)
   
   !> temporary storage
   real(wp), allocatable :: tmp(:,:)
   
   !> current number of coordinates
   integer :: coord_num

   coord_num = size(matrix,1)

   allocate(tmp(coord_num+3,coord_num+3))
   tmp(:coord_num,:coord_num) = matrix(:coord_num,:coord_num)
   deallocate(matrix)
   call move_alloc(tmp, matrix)

end subroutine resize_jacobian

!> calculate new postion for LA and corresponding J
subroutine newcoord(env,mol,xyz,xyz_out,idx1,idx2,jacobian,connectorPosition)
   
   implicit none
   
   !> name of error producer routine
   character(len=*), parameter :: source = "oniom_newcoord"
   
   !> calculation environment 
   type(TEnvironment), intent(inout) :: env
   
   !> molecular structure data
   type(TMolecule), intent(in) :: mol
   
   !> inreg coordinate matrix 
   real(wp), intent(inout) :: xyz(:, :)
   
   !> outreg coordinate matrix
   real(wp), intent(inout) :: xyz_out(:, :)
   
   !> connector (one that stays in the inner region)
   integer, intent(in) :: idx1
   
   !> host (one that is substitued) 
   integer, intent(in) :: idx2
   
   !> jacobian maxtrix
   real(wp), intent(inout) :: jacobian(:,:)
   
   !> ordinal number of connector atom in mol%at
   integer, intent(in) :: connectorPosition
   
   !> standard bond length between LAC(linked atom connector)-H
   real(wp) :: dist 
   
   !> standard bond length between LAC(linked atom connector)-LAH(linked atom host)
   real(wp) :: dist2
   
   !> scaling parameter
   real(wp) :: prefactor
   
   !> squared vector length
   real(wp) :: dist_12
   
   !> cleaved atom pair coordinates
   real(wp) :: xyz1(3), xyz2(3)

   !> message for default parameter
   character(len=:), allocatable :: warning
   
   !> message for auto-switching between derived and fixed modes 
   character(len=:), allocatable :: warning2
   
   !> default is used
   logical :: def
   
   !> control warnings
   logical, save :: rep=.false.
   
   !> ordinal numbers as characters 
   character(len=3) :: dummy1, dummy2
   
   !> difference between new and old coordinates
   real(wp) :: xyz_difference(3)
   
   integer :: i,j,k
   
   !-------!
   ! SETUP !
   !-------!
   
   write (dummy1, '(I3)') mol%at(idx1)
   write (dummy2, '(I3)') mol%at(idx2)
   
   def = .false.
   warning = "Atoms "//dummy1//" and "//dummy2//" are not accounted in the parameter suite(S,P,N,C,O), the default distance values will be used."
   warning2 = "The distance between atoms "//dummy1//" and "//dummy2//" is almost the same, switching to fixed regime"
   
   ! coordinates of connector atom !
   xyz1 = mol%xyz(:, idx1)

   ! coordinates of host atom !
   xyz2 = mol%xyz(:, idx2) 
   
   dist_12=sum((xyz1 - xyz2)**2)

   ! identify average bond distances !
   ! b/n connector-H and connector-host !
   select case (mol%at(idx1))
   case default
      def = .true.
      
      ! C-H, def !
      dist = 1.084*aatoau
      ! C-C, def !
      dist2 = 1.528*aatoau
   
   ! H !
   case (1) 
      
      ! H-H !    
      dist = 0.740*aatoau
      
      select case (mol%at(idx2))
      case default
         dist2 = 1.084*aatoau
         def = .true.
      case(1)
         ! H-H !
         dist2 = 0.740*aatoau
      case (6)
         ! H-C !
         dist2 = 1.084*aatoau
      case (8)
         ! H-O !
         dist2 = 0.964*aatoau
      case (7)
         ! H-N !
         dist2 = 1.024*aatoau
      case (15)
         ! H-P !
         dist2 = 1.414*aatoau
      case (16)
         ! H-S !
         dist2 = 1.389*aatoau
      end select
   
   ! C !
   case (6)
      
      ! C-H !
      dist = 1.084*aatoau
      
      select case (mol%at(idx2))
      case default
         dist2 = 1.528*aatoau
         def = .true.
      case(1)
         ! C-H !
         dist2 = 1.084*aatoau
      case (6)
         ! C-C !
         dist2 = 1.528*aatoau
      case (8)
         ! C-O !
         dist2 = 1.430*aatoau
      case (7)
         ! C-N !
         dist2 = 1.475*aatoau
      case (15)
         ! C-P !
         dist2 = 1.860*aatoau
      case (16)
         ! C-S !
         dist2 = 1.750*aatoau
      end select
    
   ! N !
    case (7)
      
      ! N-H !
      dist = 1.024*aatoau
      
      select case (mol%at(idx2))
      case default
         dist2 = 1.470*aatoau
         def = .true.
      case(1)
         ! N-H !
         dist2 = 1.024*aatoau
      case (6)
         ! N-C !
         dist2 = 1.475*aatoau
      case (8)
         ! N-O !
         dist2 = 1.360*aatoau
      case (7)
         ! N-N !
         dist2 = 1.470*aatoau
      case (15)
         ! N-P !
         dist2 = 1.770*aatoau
      case (16)
         ! N-S !
         dist2 = 1.650*aatoau
      end select

  
   ! O !
   case (8)
      
      ! O-H !
      dist = 0.964*aatoau

      select case (mol%at(idx2))
      case default
         dist2 = 1.450*aatoau
         def = .true.
      case(1)
         ! O-H !
         dist2 = 0.964*aatoau
      case (6)
         ! O-C !
         dist2 = 1.430*aatoau
      case (8)
         ! O-O !
         dist2 = 1.450*aatoau
      case (7)
         ! O-N !
         dist2 = 1.360*aatoau
      case (15)
         ! O-P !
         dist2 = 1.750*aatoau
      case (16)
         ! O-S !
         dist2 = 1.500*aatoau
      end select

   ! P !
   case (15)
      
      ! P-H !
      dist = 1.414*aatoau
      
      select case (mol%at(idx2))
      case default
         dist2 = 1.770*aatoau
         def = .true.
      case(1)
         ! P-H !
         dist2 = 1.414*aatoau
      case (6)
         ! P-C !
         dist2 = 1.860*aatoau
      case (8)
         ! P-O !
         dist2 = 1.750*aatoau
      case (7)
         ! P-N !
         dist2 = 1.770*aatoau
      endselect

   ! S !
   case (16)
      
      ! S-H !
      dist = 1.389*aatoau
      
      select case (mol%at(idx2))
      case default
         dist2 = 1.650*aatoau
         def = .true.
      case(1)
         ! S-H !
         dist2 = 1.389*aatoau
      case (6)
         ! S-C !
         dist2 = 1.750*aatoau
      case (8)
         ! S-O !
         dist2 = 1.500*aatoau
      case (7)
         ! S-N !
         dist2 = 1.650*aatoau
      case (16)
         ! S-S !
         dist2 = 2.040*aatoau
      endselect

   end select 
   
   ! derived mode !
   if (set%oniom_settings%derived) then
      
      prefactor = dist/sqrt(dist_12)
   
   ! fixed mode !
   else
      
      prefactor = dist/dist2 
      
      ! default values !
      if(def.and. .not.rep) then
         rep=.true.
         call env%warning(warning,source)
      endif
   endif

   ! LA coordinates / reference formula (3) !
   xyz(:, size(xyz, 2)) = xyz1 + (xyz2 - xyz1) * prefactor
   xyz_out(:,size(xyz_out,2)) = xyz2 + (xyz1 - xyz2) * prefactor 
      
   ! determine the difference between LA and LAH cooordinates !
   ! (yes -> change from derived to fixed) ! 
   xyz_difference=xyz2-xyz(:,size(xyz,2))
   if (all(xyz_difference<1.0E-5).and.set%oniom_settings%derived) then 
      set%oniom_settings%derived=.false.
      call env%warning(warning2,source)
      prefactor = dist/dist2 
   endif

   ! take derivatives !
   call derivative(jacobian,connectorPosition,size(xyz,2),prefactor,mol%xyz,idx1,idx2,dist_12,set%oniom_settings%derived)

end subroutine newcoord

!> increment Jacobian matrix for newly added atoms
subroutine derivative(jacobian,con,link,prefactor,xyz,idx1,idx2,dist_12,derived)
   
   implicit none
   
   !> jacobian matrix
   real(wp), intent(inout) :: jacobian(:,:)
   
   !> position of connector atom 
   integer,intent(in) :: con
   
   !> position of linked atom 
   integer, intent(in) :: link

   !> scaling factor 
   real(wp),intent(in) :: prefactor
   
   !> square of distance vector
   real(wp),intent(in) :: dist_12

   !> coordinates of whole molecule
   real(wp),intent(in) :: xyz(:,:)
   
   !> connector (one that stays in inner region)
   integer, intent(in) :: idx1
   
   !> host (one that is substitued)
   integer, intent(in) :: idx2
   
   !> fix value of prefactor variable
   logical, intent(in) :: derived
   
   !> account for all 3 coordinates in J matrix
   integer :: con3, link3
   
   !> save the positions of changed matrix elements
   integer :: counter1(3), counter2(3)
   
   integer :: i, j
   
   !-------!
   ! SETUP !
   !-------!

   con3=con*3
   link3=link*3
   
   do i=1,3
      counter1(i)=con3-i+1
      counter2(i)=link3-i+1
   enddo
   
   ! nullify all new matrix elements !
   jacobian(counter2(3):,:) = 0.0_wp
   jacobian(:,counter2(3):) = 0.0_wp

   !-------------------------!
   ! JACOBIAN INCREMENTATION !
   !-------------------------!
   do i=1,3
      
      ! fixed mode / reference ESI formula (3,4) !
      if(.not.derived) then
         jacobian(counter1(i),counter2(i)) = 1 - prefactor
         jacobian(counter2(i),counter2(i)) = prefactor
      
      ! derived mode / reference ESI formula (6,7) !
      else

         ! x coordinate !
         if (i==1) then
            jacobian(counter2(i),counter2(i)) = prefactor * ( (xyz(2,idx1)-xyz(2,idx2)**2) & 
                                                            & + ((xyz(3,idx1)-xyz(3,idx2))**2)  ) / dist_12
         
         ! y coordinate !
         else if (i==2) then
            jacobian(counter2(i),counter2(i)) = prefactor * ( (xyz(1,idx1)-xyz(1,idx2)**2) & 
                                                            & + ((xyz(3,idx1)-xyz(3,idx2))**2)  ) / dist_12
         
         ! z coordinate !
         else
            jacobian(counter2(i),counter2(i)) = prefactor * ( (xyz(1,idx1)-xyz(1,idx2)**2) & 
                                                            & + ((xyz(2,idx1)-xyz(2,idx2))**2)  ) / dist_12
         
         endif
         
         jacobian(counter1(i),counter2(i)) = 1.0_wp - prefactor * (1.0_wp - (((xyz(i,idx2)-xyz(i,idx1))**2) / dist_12 ))
      
      endif
   
   enddo
   
end subroutine derivative

!> assign methods
function string_to_id(string) result(id)

   implicit none
   character(len=*), intent(in) :: string
   integer :: id

   select case (string)

   case default
      id = -1
   case ('gfn2')
      id = 2
   case ('gfn1')
      id = 1
   case ('gfnff')
      id = 3
   case ('orca')
      id = 4
   case ('turbomole')
      id = 5
   end select

end function string_to_id

!> check bond order
subroutine checkfororder(env, mol, idx1, idx2, bond, hybrid)
   
   implicit none
   
   !> name of error producer routine
   character(len=*), parameter :: source = 'xtb_oniom_checkfororder' 
   
   !> hybridization info from GFN-FF; topo%hyb
   integer, intent(in), optional :: hybrid(:)
   
   !> wiberg bond orders 
   integer, intent(in), optional :: bond
   
   !> calculation environment 
   type(TEnvironment), intent(inout) :: env
   
   !> molecular structure data
   type(TMolecule), intent(in) :: mol
   
   !> connector (one that stays in inner region)
   integer, intent(in) :: idx1

   !> host (one that is substitued)
   integer, intent(in) :: idx2
   
   character(len=:), allocatable :: warning
   integer :: b
   character(len=5) :: dummy1, dummy2

   write (dummy1, '(I5)') idx1
   write (dummy2, '(I5)') idx2
   warning = "You are cutting a bond with the order higher than 1 between "//dummy1//" and "//dummy2
   
   if (present(bond)) then
      if (bond > 1) then
         call env%error(warning, source)
         return
      end if
   else
      if (hybrid(idx1) /= 0 .and. hybrid(idx1) < 3) then
         if (hybrid(idx2) /= 0 .and. hybrid(idx2) < 3) then
            call env%error(warning, source)
            return
         end if
      end if
   end if

end subroutine checkfororder

!> automatic inner region charge determination
function calculateCharge(self, env, mol, chk) result(chrg_model)
   
   implicit none
   
   !> name of error producer routine
   character(len=*), parameter :: source = 'xtb_oniom_calculateCharge'
   
   !> polymorhic calculator
   class(TOniomCalculator), intent(inout) :: self
   
   !> calculation environment
   type(TEnvironment), intent(inout) :: env
   
   !> molecular structure data
   type(TMolecule), intent(in) :: mol
   
   !> wavefuntion wrapper
   type(TRestart), intent(in) :: chk
   
   !> inner region charge
   real(wp) :: charge
   
   integer :: i, j, n, k, pre_last
   integer :: chrg_model
   integer, allocatable :: at(:)

   charge = 0.0_wp

   select type (calc => self%real_low)
   
   ! gfnff !
   type is (TGFFCalculator)
      do i = 1, size(self%idx)
         charge = charge + calc%topo%qa(self%idx(i))
      end do

   ! GFN1/2 !
   type is (TxTBCalculator)
      do i = 1, size(self%idx)
         charge = charge + chk%wfn%q(self%idx(i))         
      end do

   class default
      call env%error("Not possible to calculate with external methods for real region", source)
      return
   end select
      
   chrg_model = nint(charge)

end function calculateCharge

!> transform matrix into 1-dim array
subroutine matrix_to_array(mtrx,arr)
   
   implicit none
   real(wp), intent(in) :: mtrx (:,:)
   real(wp), allocatable, intent(out) :: arr(:)
   integer :: i, j, k

   allocate(arr(size(mtrx,1)*size(mtrx,2)))
   
   k=1
   do i=1, size(mtrx,2)
      do j=1, size(mtrx,1)
         arr(k)=mtrx(j,i)
         k=k+1
      enddo
   enddo

end subroutine matrix_to_array

!> transform 1-dim array to matrix
subroutine array_to_matrix(arr,mtrx)
   
   implicit none
   real(wp), intent(out) :: mtrx (:,:)
   real(wp), allocatable, intent(inout) :: arr(:)
   integer :: i, j, k
   
   k=1
   do i=1, size(mtrx,2)
      do j=1, size(mtrx,1)
         mtrx(j,i)=arr(k)
         k=k+1
      enddo
   enddo

   deallocate(arr)

end subroutine array_to_matrix
   
!> create origin.coord (if coord exist)
subroutine protectCoord(env)
   
   use xtb_readin, only : mirror_line
   implicit none
   character(len=*),parameter :: source = "xtb_oniom_protectCoord"
   type(TEnvironment),intent(inout) :: env
   integer :: cunit, new_cunit, err
   character(len=:),allocatable :: line
   logical :: exist

   inquire(file='coord',exist=exist)
   if(exist) then   
      call open_file(cunit,'coord','r')
      call open_file(new_cunit,'origin.coord','w')
      do
         call mirror_line(cunit,new_cunit,line,err)
         if(is_iostat_end(err)) exit
      enddo
      call env%warning("coord file will be renamed as unopt.coord",source)

   endif

end subroutine protectCoord

end module xtb_oniom
