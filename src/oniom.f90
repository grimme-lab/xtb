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

   public :: TOniomCalculator, newOniomCalculator, oniom_input

   type :: oniom_input
      character(len=:), allocatable :: first_arg
      character(len=:), allocatable :: second_arg
      character(len=:), allocatable :: chrg
      logical :: g
   end type oniom_input

   type, extends(TCalculator) :: TOniomCalculator !> new calculator type

      !> The method specified by the usernan
      integer :: method_low, method_high

      !> Charge, if present
      integer :: chrg_real, chrg_model

      type(TAtomList) :: list !> list of atoms in inner region

      !> QM 1, low level, whole region
      class(TCalculator), allocatable :: real_low

      !> QM 1, low level, inner region
      class(TCalculator), allocatable :: model_low

      !> QM 2, high level, inner region
      class(TCalculator), allocatable :: model_high

      !> Index list of the inner region
      integer, allocatable :: idx(:)
      
      !> Bond topology
      type(TTopology), allocatable :: topo

      type(TRestart) :: chk_low, chk_high
      
      logical :: fixed  

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

!> Create new calculator
subroutine newOniomCalculator(self, env, mol, input)
   !> Calculator instance
   type(TOniomCalculator), intent(out) :: self
   !> Computation environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol
   !> Restart checkpoint
   !type(TRestart), intent(inout) :: chk
   !> Input for Oniom calculation
   type(oniom_input), intent(in) :: input

   type(TxTBCalculator), allocatable :: xtb
   type(TOrcaCalculator), allocatable :: orca
   type(TGFFCalculator), allocatable :: gff
   integer :: icol
   integer :: i
   
   
   !> Method identification
   if (allocated(input%first_arg)) then
      !! if methods specified
      icol = index(input%first_arg, ':')
      if (icol == 0) then
         call env%error("Invalid method '"//input%first_arg//"' provided")
         return
      end if
      self%method_high = string_to_id(input%first_arg(:icol - 1))
      self%method_low = string_to_id(input%first_arg(icol + 1:))
      if (self%method_high < 0 .or. self%method_high > 5) then
         call env%error("Invalid inner method")
         return
      end if
      if (self%method_low < 0 .or. self%method_low > 3) then
         call env%error("Invalid outer method")
         return
      end if
   else
      !! default, gfn2:gfnff
      self%method_high = 2
      self%method_low = 3
   end if
   self%fixed = input%g

   self%list = TAtomList(list=input%second_arg)
   call self%list%to_list(self%idx)

   if (len(self%list) == 0) then
      call env%error("No atoms in inner region '"//input%second_arg//"'")
      return
   end if

   ! Check if the user-defined inner region is valid
   if (any(self%idx < 1) .or. any(self%idx > mol%n)) then
      call env%error('The inner region specification is not correct')
      return
   end if

   ! The allocation of new calculator
   select case (self%method_low)
   case default
      allocate (xtb)
      call newXTBCalculator(env, mol, xtb, method=self%method_low)
      call move_alloc(xtb, self%real_low)

   case (3)
      allocate (gff)
      call newGFFCalculator(env, mol, gff, "", .false.)
      call move_alloc(gff, self%real_low)

   end select

end subroutine newOniomCalculator

subroutine singlepoint(self, env, mol, chk, printlevel, restart, energy, gradient, sigma, hlgap, results)

   !? Is it right? What is inside init function? What interface/module procedure do?
   use xtb_xtb_calculator, only: TxTBCalculator
   use xtb_setmod, only: set_chrg

   class(TOniomCalculator), intent(inout) :: self

   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(inout) :: mol

   type(TRestart), intent(inout) :: chk

   integer, intent(in) :: printlevel !> Print level for IO

   logical, intent(in) :: restart

   real(wp), intent(out) :: energy

   real(wp), intent(out) :: gradient(:, :)

   real(wp), intent(out) :: sigma(:, :) !> Strain derivatives

   real(wp), intent(out) :: hlgap
   
   real(wp), allocatable :: jacobian(:,:)

   type(scc_results), intent(out) :: results

   type(TxTBCalculator), allocatable :: tmp !> Creating dummy object
   type(scc_results) :: results_low, results_high
   type(TGFFCalculator), allocatable :: gff
   type(TOrcaCalculator), allocatable :: orca
   type(TTMCalculator), allocatable :: turbo
   type(TMolecule) :: inner_mol
   logical :: exitRun
   real(wp), allocatable :: gradient_low(:, :), gradient_high(:, :)
   real(wp) :: energy_model_low, energy_model_high, hlgap_low, hlgap_high
   real(wp) :: sigma_low(3, 3), sigma_high(3, 3)
   integer :: i
   integer,allocatable :: idx2(:)

   real(wp),allocatable :: arr_gh(:), arr_gl(:)
      !! 1-dim arrays for matmul of gradient matrix and jacobian  
   

   !> Check whether the calculator is initialized
   if (.not. allocated(self%real_low)) then
      call env%error("Outer region calculator not provided")
      return
   end if

   !> Forward solvation to outer region
   if (allocated(self%solvation)) then
      call move_alloc(self%solvation, self%real_low%solvation)
   end if

   !> Check whether the low-level calculator needs a wavefunction
   if (.not.allocated(chk%wfn%qsh)) then
   !   print *, 'overwritten' - for restart
      select type(xtb => self%real_low)
      type is(TxTBCalculator)
         call newWavefunction(env, mol, xtb, chk)
      end select
   end if


   ! Perform calculation on outer region
   call self%real_low%singlepoint(env, mol, chk, printlevel, restart, &
       & energy, gradient, sigma, hlgap, results)
   !> check for charges
!   call calculateCharge(self, env, mol, chk)



   !> Creating Linked atoms
   call self%cutbond(env, mol, chk, self%topo, inner_mol,jacobian,idx2)
   call env%check(exitRun)
   if (exitRun) then
      call env%error("Could not create linked atoms")
      return
   end if
   inner_mol%chrg = float(self%chrg_model)

   !> Inner region, low method
   if (.not. allocated(self%model_low)) then

      if (printlevel > 0) &
         write (env%unit, '(/a/)') "Initializing low level method"

      select case (self%method_low)
      case default
         call env%error("Invalid low-level inner method")

      case (1, 2)
         allocate (tmp)
         call newXTBCalculator(env, inner_mol, tmp, method=self%method_low)
         call move_alloc(tmp, self%model_low)

      case (3)
         allocate (gff)
         call newGFFCalculator(env, inner_mol, gff, "", .false.)
         call move_alloc(gff, self%model_low)
      
      end select
      
      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not setup low level method")
         return
      end if
   
   end if

   if (.not.allocated(self%chk_low%wfn%qsh)) then
      select type (calc => self%model_low)
      type is (TxTBCalculator)
         call newWavefunction(env, inner_mol, calc, self%chk_low)
      end select
   end if

   if (allocated(self%real_low%solvation)) then
      !allocate(self%model_low%solvation)
      self%model_low%solvation = self%real_low%solvation
   end if
   
   ! Inner region, high method
   if (.not. allocated(self%model_high)) then

      if (printlevel > 0) &
         write (env%unit, '(/a/)') "Initializing high level method"

      select case (self%method_high)
      case (1, 2)
         allocate (tmp)
         call newXTBCalculator(env, inner_mol, tmp, method=self%method_high)
         call move_alloc(tmp, self%model_high)

      case (3)

         allocate (gff)
         call newGFFCalculator(env, inner_mol, gff, "", .false.)
         call move_alloc(gff, self%model_high)

      case (4)

         allocate (orca)
         call newOrcaCalculator(orca, env, set%ext_orca)
         call move_alloc(orca, self%model_high)

      case (5)

         allocate (turbo)
         call newTMCalculator(turbo, 1, 1)
         call move_alloc(turbo, self%model_high)

      end select
      
      if (allocated(self%real_low%solvation)) then
         
         allocate(self%model_high%solvation)
         self%model_high%solvation = self%real_low%solvation
      end if
      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not setup high level method")
         return
      end if

   end if

   if (.not.allocated(self%chk_high%wfn%qsh)) then
      select type (calc => self%model_high)
      type is (TxTBCalculator)
         call newWavefunction(env, inner_mol, calc, self%chk_high)
      end select
   end if

   !> Inner region low method, actual calculatoin

   if (printlevel > 0) &
      write (env%unit, '(/a/)') "SP of inner region with low method"

   allocate (gradient_low(3, inner_mol%n))
   energy_model_low = 0.0_wp
   gradient_low = 0.0_wp

   call self%model_low%singlepoint(env, inner_mol, self%chk_low, printlevel, restart, &
       & energy_model_low, gradient_low, sigma_low, hlgap_low, results_low)

   !> Inner region high method, actual calculatoin

   if (printlevel > 0) &
      write (env%unit, '(/a/)') "SP of inner region with high method"

   allocate (gradient_high(3, inner_mol%n))
   energy_model_high = 0.0_wp
   gradient_high = 0.0_wp

   call self%model_high%singlepoint(env, inner_mol, self%chk_high, printlevel, restart, &
       & energy_model_high, gradient_high, sigma_high, hlgap_high, results_high)

   results%dipole = results%dipole + results_high%dipole - results_low%dipole
   energy = energy + energy_model_high - energy_model_low !> The actual Oniom energy
   
   results%hl_gap=hlgap - hlgap_low + hlgap_high
   sigma=sigma- sigma_low + sigma_high

   call matrix_to_array(gradient_high,arr_gh)
   call matrix_to_array(gradient_low,arr_gl)
  
   arr_gh=matmul(jacobian,arr_gh)
   arr_gl=matmul(jacobian,arr_gl)

   
   call array_to_matrix(arr_gh,gradient_high)
   call array_to_matrix(arr_gl,gradient_low)



   !> ONIOM gradients
   do i = 1, size(idx2)
      gradient(:, idx2(i)) = gradient(:, idx2(i)) + gradient_high(:, i) - gradient_low(:, i)
   end do
   
   results%e_total = energy

   deallocate(gradient_high)
   deallocate(gradient_low)
end subroutine singlepoint

!> Evaluate hessian by finite difference for all atoms
subroutine hessian(self, env, mol0, chk0, list, step, hess, dipgrad)
   character(len=*), parameter :: source = "extern_oniom_hessian"
   !> Single point calculator
   class(TOniomCalculator), intent(inout) :: self
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

   real(wp), allocatable :: jacobian(:,:)
   integer,allocatable :: idx2(:)

   integer :: ii, jj, ic, jc, im, jm
   real(wp), allocatable :: hess_model(:, :), dipgrad_model(:, :)
   type(TMolecule) :: mol_model
   integer, allocatable :: list_model(:)

   ! Compute complete hessian for outer region
   call self%real_low%hessian(env, mol0, chk0, list, step, &
       & hess, dipgrad)

   !call chk0%wfn%allocate(mol0%n,chk0%basis%nshell,chk0%real_low%basis%nao)
   ! Creating Linked atoms
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

   integer, intent(in) :: unit !> unit for I/O

   type(Tmolecule), intent(in) :: mol !> molecular structural data

end subroutine writeInfo

subroutine cutbond(self, env, mol, chk, topo, inner_mol,jacobian,idx2)

   !> to initialize mol object
   use xtb_type_molecule, only: init

   !> To get topology info from wfn
   use xtb_topology, only: makeBondTopology, topologyToNeighbourList
   

   !> actual calculator
   class(TOniomCalculator), intent(in) :: self

   !> Environment
   type(TEnvironment), intent(inout) :: env

   !> Wavefunction
   type(TRestart), intent(in) :: chk

   !>  Molecular data storage
   type(TMolecule), intent(in) :: mol

   !> Inner region geometry
   type(TMolecule), intent(inout) :: inner_mol

   !> Some buffer geometry
   type(TMolecule) :: cash

   !> topology info
   type(TTopology), allocatable, intent(inout) :: topo

   !> neighbour lists
   type(TNeighbourList) :: neighList
   
   !>
   real(wp),allocatable, intent(inout) :: jacobian(:,:)
   
   integer,allocatable, intent(out) :: idx2(:)

   integer :: i, j, n, k, pre_last,iterator

   !> Arrays for inner region atoms and their atomic numbers
   integer, allocatable :: at(:), at2(:)
   real(wp), allocatable :: xyz(:, :), xyz2(:, :)
   logical :: inside
   integer, allocatable :: bonded(:, :)
   
   
   !> The pair-wise indices for cutbound  
   integer, allocatable :: brokenBondPairs(:,:)
   
   
   inside = .FALSE.

   n = len(self%list)
   allocate (at2(size(self%idx)))
   allocate (xyz2(3, size(self%idx)))
   allocate (at(n))
   allocate (xyz(3, n))
   allocate (idx2(size(self%idx)))
   
   !> to save inner atom list in the matching array
   idx2=self%idx

   !> Assignment of initial inner region
   do i = 1, size(self%idx)
      at(i) = mol%at(self%idx(i))
      at2(i) = mol%at(self%idx(i))
      xyz(:, i) = mol%xyz(:, self%idx(i))
      xyz2(:, i) = mol%xyz(:, self%idx(i))
   end do
 
   call create_jacobian(jacobian,at)
   
   !> To identify bonded atoms and save them into an array + iterator
   select type (calc => self%real_low)
   type is (TGFFCalculator)
      bonded = calc%topo%blist !> automatic allocation
      iterator = size(bonded,2)

   type is (TxTBCalculator)
      if (.not. allocated(topo)) then
         allocate (topo)
         call makeBondTopology(topo, mol, chk%wfn%wbo) !> return assigned topo%list

         call topologyToNeighbourList(topo, neighList, mol) !> return neighList
      end if
      allocate (bonded(2, len(topo)))
      do i = 1, len(topo)
         bonded(:, i) = topo%list(1:2, i)
      end do
      iterator = size(bonded,2)
   end select
   !> Actual bond cutting and creating linked atom
   do i = 1, size(self%idx)
      do j = 1, iterator
         if (bonded(1, j) == self%idx(i)) then
            do k = 1, size(self%idx)
               if (self%idx(k) == bonded(2, j)) then
                  inside = .TRUE.
               end if
            end do
            if (.not. inside) then

               select type (calc => self%real_low)
               class default
                  call checkfororder(env, mol, self%idx(i), bonded(2, j), bond=topo%list(3, j))
               type is (TGFFCalculator)
                  call checkfororder(env, mol, self%idx(i), bonded(2, j), hybrid=calc%topo%hyb)
               end select

               pre_last = size(at)
               call resize(at)
                  !! Increase no of atoms(size of array) by 1
               call resize(idx2)
               idx2(pre_last+1) = bonded(2,j)
               at(pre_last + 1) = 1
                  !! Transform the newly created atom to hydrogen atom(Link Atom)
               call resize(xyz)
                  !! Adjust accordingly the coordinate matrix 
               call resize_jacobian(jacobian)   
                  !! Readjust jacobian
               call coord(mol, xyz, at(pre_last), self%idx(i), bonded(2, j),jacobian,i)
                  !! Determine the coordinates or added H atom

               
            end if
            inside = .FALSE.
         else if (bonded(2, j) == self%idx(i)) then
            do k = 1, size(self%idx)
               if (self%idx(k) == bonded(1, j)) then
                  inside = .TRUE.
               end if
            end do
            if (.not. inside) then
               select type (calc => self%real_low)
               class default
                  call checkfororder(env, mol, self%idx(i), bonded(1, j), bond=topo%list(3, j))
               type is (TGFFCalculator)
                  call checkfororder(env, mol, self%idx(i), bonded(1, j), hybrid=calc%topo%hyb)
               end select
               pre_last = size(at)
               call resize(at)
               call resize(idx2)
               idx2(pre_last+1) = bonded(1,j)
               at(pre_last + 1) = 1
               call resize(xyz)
               call resize_jacobian(jacobian)   
                  !! Readjust jacobian
               call coord(mol, xyz, at(pre_last), self%idx(i), bonded(1, j),jacobian,i)
            end if
            inside = .FALSE.
         end if
      end do
   end do
   
   
   call init(inner_mol, at, xyz)
   
   
   !call init(cash,at2,xyz2)
  ! block
   !   use xtb_io_writer
   !   use xtb_mctc_filetypes, only : fileType
   !   integer :: io
   !   call open_file(io, "inner-region_with_h.xyz", "w")
   !   call writeMolecule(inner_mol, io, filetype%xyz)
   !   call close_file(io)
   !end block


end subroutine cutbond

subroutine new_atom(at)
   
   implicit none
   integer, allocatable, intent(inout) :: at(:)
   integer, allocatable :: tmp2(:)
   allocate (tmp2(size(at) + 1))
   tmp2(:size(at)) = at(:size(at))
   deallocate (at)
   call move_alloc(tmp2, at)

end subroutine new_atom

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
!
end subroutine new_coordinates

!> Create jacobian matrix
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


!> To increase matrix(jacobian) for new atom 
!> (3 new entries in diagonal and subsequent increase)
subroutine resize_jacobian(matrix)
   
   implicit none
   real(wp), allocatable :: matrix(:,:)
   real(wp), allocatable :: tmp(:,:)
      !! To store old values while reallocation
   integer :: coord_num
      !! The current number of coordinates

   coord_num = size(matrix,1)

   allocate(tmp(coord_num+3,coord_num+3))
   tmp(:coord_num,:coord_num) = matrix(:coord_num,:coord_num)
   deallocate(matrix)
   call move_alloc(tmp, matrix)

end subroutine resize_jacobian

!call coord(mol, xyz, pre_last, at(self%idx(i)), self%idx(i), bonded(1, j))
subroutine newcoord(mol, xyz,  at, idx1, idx2, jacobian, connectorPosition)
   
   use xtb_setparam, only : set
  
   implicit none
   type(TMolecule), intent(in) :: mol
   real(wp), intent(inout) :: xyz(:, :)
   !integer, intent(in) :: pre_last
   integer, intent(in) :: at, idx1, idx2
   real(wp), intent(inout) :: jacobian(:,:)
   integer, intent(in) :: connectorPosition

   real(wp), allocatable :: tmp_mtrx(:)
   real(wp) :: dist, prefactor,dist_13
   real(wp) :: xyz1(3), xyz2(3)
   integer :: i,j,k


   select case (at)
   case default
      dist = 1.084*aatoau
      !> Carbon
   case (6)
      dist = 1.09*aatoau
      !> Oxygen
   case (8)
      dist = 0.964*aatoau
      !> Nitrogen
   case (7)
      dist = 1.024*aatoau
      !> Phosphorus
   case (15)
      dist = 1.414*aatoau
      !> Sulfur
   case (16)
      dist = 1.389*aatoau
   end select
   

   xyz1 = mol%xyz(:, idx1)
   xyz2 = mol%xyz(:, idx2)
   prefactor = dist/sqrt(sum((xyz1 - xyz2)**2))
  
   dist_13=sum((xyz1 - xyz2)**2)
      !! dist_13 = dist**2
   if (set%g_fixed) prefactor=0.409*aatoau 

   !> new H coordinates
   xyz(:, size(xyz, 2)) = xyz1 + (xyz2 - xyz1) * prefactor
   
   call derivative(jacobian,connectorPosition,size(xyz,2),prefactor,xyz,idx1,idx2,dist_13,set%g_fixed)

end subroutine newcoord

!> To take derivative of model system coordinates wrt real system coordinates
!> see https://xtb-docs.readthedocs.io/ for more info
subroutine derivative(jacobian,con,link,prefactor,xyz,idx1,idx2,dist_13,fixed)
   
   implicit none
   real(wp), intent(inout) :: jacobian(:,:)
   integer,intent(in) :: con
      !! position of connector atom in the model system atom list
   integer, intent(in) :: link
      !! position of linked atom in the model system atom list
   real(wp),intent(in) :: prefactor,dist_13
   real(wp),intent(in) :: xyz(:,:)
   integer, intent(in) :: idx1, idx2
   logical, intent(in) :: fixed
      !! To fix value of prefactor variable
   integer :: i, j
   integer :: con3, link3
      !! to account for all 3 coordinates
   integer :: counter1(3), counter2(3)
      !! to save the positions of changed matrix elements
   
   con3=con*3
   link3=link*3

   do i=1,3
      counter1(i)=con3-i+1
      counter2(i)=link3-i+1
   enddo
   
   !> Assign all new matrix elements to 0
   jacobian(counter2(3):,:) = 0.0_wp
   jacobian(:,counter2(3):) = 0.0_wp

   !> Algorithm to find non-zero elements of the matrix
   do i=1,3
      
      if(fixed) then
         !> Fixed Jacobian
         jacobian(counter1(i),counter2(i))=1-prefactor
         jacobian(counter2(i),counter2(i))=prefactor
      else
         !> Derived Jacobian
         jacobian(counter1(i),counter2(i))=1 + (prefactor*((xyz(i,idx2)-xyz(i,idx1))**2))/(dist_13**(3.0_wp/2.0_wp))-(prefactor/(dist_13**(1.0_wp/2.0_wp))) 
         jacobian(counter2(i),counter2(i))=(prefactor*(sqrt(dist_13)-xyz(i,idx2)+xyz(i,idx1))*(sqrt(dist_13)+xyz(i,idx2)-xyz(i,idx1)))/(dist_13**(3.0_wp/2.0_wp))
      endif
      
   
   enddo
   
end subroutine derivative

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

subroutine checkfororder(env, mol, idx1, idx2, bond, hybrid)

   character(len=*), parameter :: source = 'oniom_cutbond'

   integer, intent(in), optional :: hybrid(:)

   integer, intent(in), optional :: bond

   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(in) :: mol

   integer, intent(in) :: idx1, idx2

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

subroutine calculateCharge(self, env, mol, chk)

   character(len=*), parameter :: source = 'charge for inner region'

   class(TOniomCalculator), intent(inout) :: self

   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(in) :: mol

   type(TRestart), intent(in) :: chk

   integer :: i, j, n, k, pre_last
   integer, allocatable :: at(:)
   real(wp) :: charge

   charge = 0.0_wp

   select type (calc => self%real_low)
   type is (TGFFCalculator)

      do i = 1, size(self%idx)
         charge = charge + calc%topo%qa(self%idx(i))

      end do
     type is (TxTBCalculator)
     do i = 1, size(self%idx)
         charge = charge + chk%wfn%q(self%idx(i))
     end do

    class default
        call env%error("Not possible to calculate with external methods for real region", source)
        return
    end select

    self%chrg_model = nint(charge)

end subroutine calculateCharge

!> To transform matrix into 1 dimensional array
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
 
!> To transform 1 dimensional array to matrix
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
   
  
end module xtb_oniom
