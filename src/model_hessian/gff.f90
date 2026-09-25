! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
! Copyright (C) 2026 Leopold M. Seidler
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

!> GFN-FF-specific model Hessian implementation
module xtb_modelhessian_gff
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_math, only : crossProd
   use xtb_bmatrix, only : bmat_bond, bmat_angle, bmat_linbend, linbend_frame, &
      & bmat_torsion, bmat_accum_packed, bmat_accum_pairblock_packed
   use xtb_gfnff_data, only : TGFFData
   use xtb_gfnff_neighbor, only : TNeigh
   use xtb_gfnff_topology, only : TGFFTopology
   use xtb_modelhessian_shared, only : itabrow
   use xtb_param_model_hessian, only : rav => legacy_rav, &
      & aav => legacy_aav, c6 => gff_c6, gff_stretch_constant, &
      & gff_bend_constant, gff_torsion_constant, &
      & distance_threshold => gff_distance_threshold
   use xtb_modelhessian_type, only : TModelHessian
   use xtb_type_environment, only : TEnvironment
   use xtb_type_setvar, only : modhess_setvar
   implicit none(type, external)
   private

   !> GFN-FF model Hessian using calculator-owned topology data
   type, public, extends(TModelHessian) :: TGFFModelHessian
      private
      type(TGFFData), pointer :: param => null()
      type(TGFFTopology), pointer :: topo => null()
      type(TNeigh), pointer :: neigh => null()
      real(wp) :: s6
      logical :: do_bend, do_torsion
   contains
      procedure, public :: compute_packed
   end type TGFFModelHessian

   public :: newGFFModelHessian

contains

!> Create a GFN-FF model Hessian bound to calculator-owned data
function newGFFModelHessian(param, topo, neigh, modh) result(model_hessian)
   !> GFN-FF parameters
   type(TGFFData), intent(in), target :: param
   !> Molecular GFN-FF topology
   type(TGFFTopology), intent(in), target :: topo
   !> GFN-FF neighbor list
   type(TNeigh), intent(in), target :: neigh
   !> Model Hessian configuration
   type(modhess_setvar), intent(in) :: modh

   type(TGFFModelHessian) :: model_hessian

   model_hessian%param => param
   model_hessian%topo => topo
   model_hessian%neigh => neigh
   model_hessian%s6 = modh%s6
   model_hessian%do_bend = modh%kf /= 0.0_wp
   model_hessian%do_torsion = modh%kt /= 0.0_wp
end function newGFFModelHessian

!> Compute a packed GFN-FF model Hessian
subroutine compute_packed(self, env, xyz, n, hess, at)
   !> GFN-FF model Hessian holding parameters, topology, and neighbour list
   class(TGFFModelHessian), intent(in) :: self
   !> Calculation environment
   type(TEnvironment), intent(inout) :: env
   !> Number of atoms
   integer, intent(in) :: n
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Packed lower-triangle Hessian
   real(wp), intent(out) :: hess((3*n)*(3*n + 1)/2)
   !> Atomic numbers
   integer, intent(in) :: at(n)

   hess = 0.0_wp
   call stretch(self, xyz, n, hess, at)
   if (self%do_bend) call bend(self, xyz, n, hess, at)
   if (self%do_torsion) call torsion(self, xyz, n, hess, at)
end subroutine compute_packed


!> Add GFN-FF bond, Coulomb, and dispersion contributions
subroutine stretch(self, xyz, n, hess, at)
   !> GFN-FF model Hessian holding parameters, topology, and neighbour list
   class(TGFFModelHessian), intent(in) :: self
   !> Number of atoms
   integer, intent(in) :: n
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Packed lower-triangle Hessian, updated in place
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
   !> Atomic numbers
   integer, intent(in) :: at(n)

   integer :: ibond, i, j, ir, jr
   integer :: mapped_at(n)
   real(wp) :: vec(3), r2, gmm, cdisp, qq, r0_squared
   real(wp) :: pair_hessian(3, 3)

   mapped_at = gff_atomic_number(at)

   do ibond = 1, self%neigh%nbond
      i = self%neigh%blist(1, ibond)
      j = self%neigh%blist(2, ibond)
      ir = itabrow(mapped_at(i))
      jr = itabrow(mapped_at(j))
      vec = xyz(:, i) - xyz(:, j)
      r2 = dot_product(vec, vec)
      gmm = gff_stretch_constant * exp(aav(ir, jr)*(rav(ir, jr)**2 - r2))
      call bmat_accum_packed(n, hess, [i, j], bmat_bond(vec), gmm)
   end do

   do i = 1, n
      do j = 1, i - 1
         vec = xyz(:, i) - xyz(:, j)
         r2 = dot_product(vec, vec)
         if (r2 > 1600.0_wp) cycle
         cdisp = -self%s6 * sqrt(c6(mapped_at(i))*c6(mapped_at(j)))
         qq = 2.0_wp * self%topo%qa(i) * self%topo%qa(j)
         r0_squared = self%param%d3r0(pair_index(at(i), at(j)))
         call get_pair_hessian(vec, qq, cdisp, r0_squared, pair_hessian)
         call bmat_accum_pairblock_packed(n, hess, i, j, -pair_hessian)
      end do
   end do
end subroutine stretch

!> Add GFN-FF angle-bending contributions
subroutine bend(self, xyz, n, hess, at)
   !> GFN-FF model Hessian holding parameters, topology, and neighbour list
   class(TGFFModelHessian), intent(in) :: self
   !> Number of atoms
   integer, intent(in) :: n
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Packed lower-triangle Hessian, updated in place
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
   !> Atomic numbers
   integer, intent(in) :: at(n)

   integer :: iangl, i, j, m, ir, jr, mr
   integer :: mapped_at(n)
   real(wp) :: vec_mi(3), vec_mj(3), vec_ij(3)
   real(wp) :: rmi2, rmj2, rmi, rmj, rij
   real(wp) :: gij, cross_vec(3), sinphi
   real(wp) :: bmat9(9), evec1(3), evec2(3)

   mapped_at = gff_atomic_number(at)
   do iangl = 1, self%topo%nangl
      m = self%topo%alist(1, iangl)
      i = self%topo%alist(2, iangl)
      j = self%topo%alist(3, iangl)
      mr = itabrow(mapped_at(m))
      ir = itabrow(mapped_at(i))
      jr = itabrow(mapped_at(j))

      vec_mi = xyz(:, i) - xyz(:, m)
      vec_mj = xyz(:, j) - xyz(:, m)
      vec_ij = xyz(:, j) - xyz(:, i)
      rmi2 = dot_product(vec_mi, vec_mi)
      rmj2 = dot_product(vec_mj, vec_mj)
      rmi = sqrt(rmi2)
      rmj = sqrt(rmj2)
      rij = norm2(vec_ij)
      if (rmi <= distance_threshold .or. rmj <= distance_threshold &
            & .or. rij <= distance_threshold) cycle

      gij = gff_bend_constant * exp( &
         aav(mr, ir) * rav(mr, ir)**2 + aav(mr, jr) * rav(mr, jr)**2 &
         - aav(mr, ir) * rmi2 - aav(mr, jr) * rmj2)
      cross_vec = crossProd(vec_mi, vec_mj)
      sinphi = norm2(cross_vec) / (rmi*rmj)
      if (sinphi > distance_threshold) then
         bmat9 = bmat_angle(vec_mi, vec_mj)
         call bmat_accum_packed(n, hess, [i, m, j], bmat9, gij)
      else
         ! linear centre: the two Decius bends along a fixed frame
         ! perpendicular to the axis, k B^T B summed over both
         call linbend_frame(vec_mi/rmi, evec1, evec2)
         call bmat_accum_packed(n, hess, [i, m, j], &
            & bmat_linbend(vec_mi, vec_mj, evec1), gij)
         call bmat_accum_packed(n, hess, [i, m, j], &
            & bmat_linbend(vec_mi, vec_mj, evec2), gij)
      end if
   end do
end subroutine bend

!> Add GFN-FF torsional contributions
subroutine torsion(self, xyz, n, hess, at)
   !> GFN-FF model Hessian holding parameters, topology, and neighbour list
   class(TGFFModelHessian), intent(in) :: self
   !> Number of atoms
   integer, intent(in) :: n
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Packed lower-triangle Hessian, updated in place
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
   !> Atomic numbers
   integer, intent(in) :: at(n)

   integer :: itors, i, j, k, l, ir, jr, kr, lr
   integer :: mapped_at(n)
   real(wp) :: torsion_xyz(3, 4), bmat(3, 4), brow12(12)
   real(wp) :: rij(3), rjk(3), rkl(3), tij

   mapped_at = gff_atomic_number(at)
   do itors = 1, self%topo%ntors
      i = self%topo%tlist(3, itors)
      j = self%topo%tlist(1, itors)
      k = self%topo%tlist(2, itors)
      l = self%topo%tlist(4, itors)
      ir = itabrow(mapped_at(i))
      jr = itabrow(mapped_at(j))
      kr = itabrow(mapped_at(k))
      lr = itabrow(mapped_at(l))

      torsion_xyz = xyz(:, [i, j, k, l])
      rij = xyz(:, i) - xyz(:, j)
      rjk = xyz(:, j) - xyz(:, k)
      rkl = xyz(:, k) - xyz(:, l)
      tij = gff_torsion_constant * exp( &
         aav(ir, jr) * (rav(ir, jr)**2 - dot_product(rij, rij)) &
         + aav(jr, kr) * (rav(jr, kr)**2 - dot_product(rjk, rjk)) &
         + aav(kr, lr) * (rav(kr, lr)**2 - dot_product(rkl, rkl)))
      bmat = bmat_torsion(torsion_xyz)
      brow12 = [bmat(:, 1), bmat(:, 2), bmat(:, 3), bmat(:, 4)]
      call bmat_accum_packed(n, hess, [i, j, k, l], brow12, tij)
   end do
end subroutine torsion

!> Map heavy elements to the legacy GFN-FF model-Hessian C6 table
pure elemental function gff_atomic_number(at) result(mapped_at)
   !> Atomic number
   integer, intent(in) :: at

   integer :: mapped_at

   mapped_at = at
   if (at > 54) mapped_at = at - 18
   if (at > 72) mapped_at = at - 32
   if (at > 56 .and. at < 72) mapped_at = 39
   if (at > 86) mapped_at = 55
end function gff_atomic_number

!> Return the packed pair index for two atomic numbers
pure elemental function pair_index(i, j) result(index)
   !> Atomic numbers of the pair
   integer, intent(in) :: i, j

   integer :: index

   index = min(i, j) + max(i, j) * (max(i, j) - 1) / 2
end function pair_index

!> Evaluate the damped Coulomb and dispersion Cartesian Hessian block
pure subroutine get_pair_hessian(vec, qq, cdisp, r0_squared, hessian)
   !> Cartesian pair displacement
   real(wp), intent(in) :: vec(3)
   !> Coulomb pair factor
   real(wp), intent(in) :: qq
   !> Dispersion pair factor
   real(wp), intent(in) :: cdisp
   !> Squared D3 reference distance
   real(wp), intent(in) :: r0_squared
   !> Symmetric Cartesian pair-Hessian block
   real(wp), intent(out) :: hessian(3, 3)

   real(wp) :: r2, r, r3, damped_r

   r2 = dot_product(vec, vec)
   r = sqrt(r2)
   r3 = r * r2
   damped_r = r + sqrt(r0_squared)
   call getqqxx(vec(1), qq, cdisp, r, r2, r3, damped_r, r0_squared, hessian(1, 1))
   call getqqxy(vec(1), vec(2), qq, cdisp, r, r2, r3, damped_r, r0_squared, hessian(1, 2))
   call getqqxy(vec(1), vec(3), qq, cdisp, r, r2, r3, damped_r, r0_squared, hessian(1, 3))
   call getqqxx(vec(2), qq, cdisp, r, r2, r3, damped_r, r0_squared, hessian(2, 2))
   call getqqxy(vec(2), vec(3), qq, cdisp, r, r2, r3, damped_r, r0_squared, hessian(2, 3))
   call getqqxx(vec(3), qq, cdisp, r, r2, r3, damped_r, r0_squared, hessian(3, 3))
   hessian(2, 1) = hessian(1, 2)
   hessian(3, 1) = hessian(1, 3)
   hessian(3, 2) = hessian(2, 3)
end subroutine get_pair_hessian

!> Evaluate a diagonal pair-Hessian element
pure subroutine getqqxx(dx, qq, cdisp, r, r2, r3, damped_r, r0_squared, d2)
   !> Cartesian separation component
   real(wp), intent(in) :: dx
   !> Coulomb pair factor
   real(wp), intent(in) :: qq
   !> Dispersion pair factor
   real(wp), intent(in) :: cdisp
   !> Pair distance
   real(wp), intent(in) :: r
   !> Squared pair distance
   real(wp), intent(in) :: r2
   !> Cubed pair distance
   real(wp), intent(in) :: r3
   !> Damped pair distance
   real(wp), intent(in) :: damped_r
   !> Squared D3 reference distance
   real(wp), intent(in) :: r0_squared
   !> Diagonal second derivative
   real(wp), intent(out) :: d2

   real(wp) :: damped_r2, dx2, denominator, r6, r8

   damped_r2 = damped_r**2
   dx2 = dx**2
   d2 = qq * (2.0_wp*dx2/(r2*damped_r*damped_r2) &
      + dx2 / (r3*damped_r2) - 1.0_wp / (r*damped_r2))
   r6 = r3 * r3
   r8 = r6 * r2
   denominator = r0_squared**3 + r6
   d2 = d2 + cdisp * (dx2*72.0_wp*r8/denominator**3 &
      - dx2 * 24.0_wp * r2 / denominator**2 - 6.0_wp * r2 * r2 / denominator**2)
end subroutine getqqxx

!> Evaluate a mixed pair-Hessian element
!>
!> `r` is kept for a shared argument list with `getqqxx` and is unused here.
pure subroutine getqqxy(dx, dy, qq, cdisp, r, r2, r3, damped_r, r0_squared, d2)
   !> First Cartesian separation component
   real(wp), intent(in) :: dx
   !> Second Cartesian separation component
   real(wp), intent(in) :: dy
   !> Coulomb pair factor
   real(wp), intent(in) :: qq
   !> Dispersion pair factor
   real(wp), intent(in) :: cdisp
   !> Pair distance, unused
   real(wp), intent(in) :: r
   !> Squared pair distance
   real(wp), intent(in) :: r2
   !> Cubed pair distance
   real(wp), intent(in) :: r3
   !> Damped pair distance
   real(wp), intent(in) :: damped_r
   !> Squared D3 reference distance
   real(wp), intent(in) :: r0_squared
   !> Mixed second derivative
   real(wp), intent(out) :: d2

   real(wp) :: damped_r2, denominator, r6, r8

   damped_r2 = damped_r**2
   d2 = qq * (2.0_wp*dx*dy/(r2*damped_r*damped_r2) + dx*dy/(r3*damped_r2))
   r6 = r3 * r3
   r8 = r6 * r2
   denominator = r0_squared**3 + r6
   d2 = d2 + cdisp * (dx*dy*72.0_wp*r8/denominator**3 &
      - dx * dy * 24.0_wp * r2 / denominator**2)
end subroutine getqqxy

end module xtb_modelhessian_gff
