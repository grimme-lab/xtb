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

!> Suspend and restore the global geometry terms, i.e. fixed atoms,
!> constraints, walls and bias structures, which are numbered for the whole
!> system, while calculations on a subsystem with its own numbering are done
module xtb_geometry_terms
   use xtb_fixparam, only : fixset, metaset, rmsdset
   use xtb_scanparam, only : potset, nconstr
   use xtb_sphereparam, only : number_walls, wpot
   implicit none
   private

   public :: geometry_state_type, suspend_geometry_terms, restore_geometry_terms


   !> Saved counts of global geometry terms suspended for a subsystem
   type :: geometry_state_type
      !> Number of fixed atoms
      integer :: nfix = 0
      !> Number of constraining potentials
      integer :: npot = 0
      !> Number of constrained atom positions
      integer :: npos = 0
      !> Number of distance constraints
      integer :: ndist = 0
      !> Number of angle constraints
      integer :: nangle = 0
      !> Number of dihedral constraints
      integer :: ndihedral = 0
      !> Number of legacy constraints
      integer :: nconstr = 0
      !> Number of wall potentials
      integer :: nwalls = 0
      !> Number of metadynamics bias structures
      integer :: nmeta = 0
      !> Number of RMSD bias structures
      integer :: nrmsd = 0
   end type geometry_state_type

contains


!> Suspend global geometry terms during calculations on a subsystem
subroutine suspend_geometry_terms(state)

   !> Saved counts for fixed atoms, restraints, walls and bias structures
   type(geometry_state_type), intent(out) :: state

   state%nfix = fixset%n
   state%npot = potset%n
   state%npos = potset%pos%n
   state%ndist = potset%dist%n
   state%nangle = potset%angle%n
   state%ndihedral = potset%dihedral%n
   state%nconstr = nconstr
   state%nwalls = 0
   if (allocated(wpot)) state%nwalls = number_walls
   state%nmeta = metaset%nstruc
   state%nrmsd = rmsdset%nstruc

   fixset%n = 0
   potset%n = 0
   potset%pos%n = 0
   potset%dist%n = 0
   potset%angle%n = 0
   potset%dihedral%n = 0
   nconstr = 0
   number_walls = 0
   metaset%nstruc = 0
   rmsdset%nstruc = 0
end subroutine suspend_geometry_terms


!> Restore global geometry terms after inner-region calculations
subroutine restore_geometry_terms(state)

   !> Counts saved by suspend_geometry_terms
   type(geometry_state_type), intent(in) :: state

   fixset%n = state%nfix
   potset%n = state%npot
   potset%pos%n = state%npos
   potset%dist%n = state%ndist
   potset%angle%n = state%nangle
   potset%dihedral%n = state%ndihedral
   nconstr = state%nconstr
   number_walls = state%nwalls
   metaset%nstruc = state%nmeta
   rmsdset%nstruc = state%nrmsd
end subroutine restore_geometry_terms


end module xtb_geometry_terms
