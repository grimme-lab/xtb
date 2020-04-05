! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> Implementation of DFT-D3 with Becke-Johnson damping
module xtb_disp_dftd3
   use xtb_mctc_accuracy, only : wp
   use xtb_disp_dftd3param
   use xtb_type_molecule
   use xtb_type_neighbourlist
   use xtb_type_param
   implicit none
   private

   public :: d3_gradient


   interface d3_gradient
      module procedure :: d3_gradient_neigh
      module procedure :: d3_gradient_latp
   end interface d3_gradient


contains


!> Calculate the weights of the reference system and the derivatives w.r.t.
!  coordination number for later use.
subroutine weight_references(nat, atoms, wf, cn, gwvec, gwdcn)
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> Atomic numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> Exponent for the Gaussian weighting.
   real(wp), intent(in) :: wf
   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)
   !> weighting for the atomic reference systems
   real(wp), intent(out) :: gwvec(:, :)
   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(out) :: gwdcn(:, :)

   integer :: iat, ati, iref, icount
   real(wp) :: norm, dnorm, gw, expw, expd, gwk, dgwk

   gwvec = 0.0_wp
   gwdcn = 0.0_wp

   do iat = 1, nat
      ati = atoms(iat)
      norm = 0.0_wp
      dnorm = 0.0_wp
      do iref = 1, number_of_references(ati)
         gw = weight_cn(wf, cn(iat), reference_cn(iref, ati))
         norm = norm + gw
         dnorm = dnorm + 2*wf*(reference_cn(iref, ati) - cn(iat)) * gw
      end do
      norm = 1.0_wp / norm
      do iref = 1, number_of_references(ati)
         expw = weight_cn(wf, cn(iat), reference_cn(iref, ati))
         expd = 2*wf*(reference_cn(iref, ati) - cn(iat)) * expw

         gwk = expw * norm
         if (gwk /= gwk) then
            if (maxval(reference_cn(:number_of_references(ati), ati)) &
               & == reference_cn(iref, ati)) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         gwvec(iref, iat) = gwk

         dgwk = expd*norm-expw*dnorm*norm**2
         if (dgwk /= dgwk) then
            dgwk = 0.0_wp
         endif
         gwdcn(iref, iat) = dgwk

      end do
   end do

end subroutine weight_references


!> calculate atomic dispersion coefficients and their derivatives w.r.t.
!  the coordination number.
subroutine get_atomic_c6(nat, atoms, gwvec, gwdcn, c6, dc6dcn)
   !> Nr. of atoms (without periodic images)
   integer, intent(in) :: nat
   !> numbers of every atom.
   integer, intent(in) :: atoms(:)
   !> weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)
   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in) :: gwdcn(:, :)
   !> C6 coefficients for all atom pairs.
   real(wp), intent(out) :: c6(:, :)
   !> derivative of the C6 w.r.t. the coordination number
   real(wp), intent(out) :: dc6dcn(:, :)

   integer :: iat, jat, ati, atj, iref, jref
   real(wp) :: refc6, dc6, dc6dcni, dc6dcnj

   c6 = 0.0_wp
   dc6dcn = 0.0_wp

   do iat = 1, nat
      ati = atoms(iat)
      do jat = 1, iat
         atj = atoms(jat)
         dc6 = 0.0_wp
         dc6dcni = 0.0_wp
         dc6dcnj = 0.0_wp
         do iref = 1, number_of_references(ati)
            do jref = 1, number_of_references(atj)
               refc6 = get_c6(iref, jref, ati, atj)
               dc6 = dc6 + gwvec(iref, iat) * gwvec(jref, jat) * refc6
               dc6dcni = dc6dcni + gwdcn(iref, iat) * gwvec(jref, jat) * refc6
               dc6dcnj = dc6dcnj + gwvec(iref, iat) * gwdcn(jref, jat) * refc6
            end do
         end do
         c6(iat, jat) = dc6
         c6(jat, iat) = dc6
         dc6dcn(iat, jat) = dc6dcni
         dc6dcn(jat, iat) = dc6dcnj
      end do
   end do
end subroutine get_atomic_c6


subroutine d3_gradient_latp &
      & (mol, trans, par, weighting_factor, r4r2, cutoff, &
      &  cn, dcndr, dcndL, energy, gradient, sigma)

   type(TMolecule), intent(in) :: mol
   type(dftd_parameter), intent(in) :: par

   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: weighting_factor
   real(wp), intent(in) :: r4r2(:)
   real(wp), intent(in) :: cutoff
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)

   real(wp), intent(inout) :: energy
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)

   integer :: nat, max_ref
   integer :: iat, jat, ati, atj, itr

   real(wp) :: cutoff2
   real(wp) :: r4r2ij, r0, rij(3), r2, t6, t8, t10, d6, d8, d10
   real(wp) :: dE, dG(3), dS(3, 3), disp, ddisp

   real(wp), allocatable :: gw(:, :), dgwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: energies(:), dEdcn(:)

   nat = len(mol)
   max_ref = maxval(number_of_references(mol%at))
   cutoff2 = cutoff**2
   allocate(gw(max_ref, nat), dgwdcn(max_ref, nat), c6(nat, nat), &
      &     dc6dcn(nat, nat), energies(nat), dEdcn(nat), source=0.0_wp)

   call weight_references(nat, mol%at, weighting_factor, cn, gw, dgwdcn)

   call get_atomic_c6(nat, mol%at, gw, dgwdcn, c6, dc6dcn)

   !$omp parallel do default(none) &
   !$omp reduction(+:energies, gradient, sigma, dEdcn) &
   !$omp shared(mol, trans, cutoff2, par, r4r2, c6, dc6dcn) &
   !$omp private(iat, jat, itr, ati, atj, r2, rij, r4r2ij, r0, t6, t8, t10, &
   !$omp&        d6, d8, d10, disp, ddisp, dE, dG, dS)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      do jat = 1, iat
         atj = mol%at(jat)
         r4r2ij = 3*r4r2(ati)*r4r2(atj)
         r0 = par%a1*sqrt(r4r2ij) + par%a2
         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-10_wp) cycle


            t6 = 1._wp/(r2**3+r0**6)
            t8 = 1._wp/(r2**4+r0**8)
            t10 = 1._wp/(r2**5+r0**10)

            d6 = -6*r2**2*t6**2
            d8 = -8*r2**3*t8**2
            d10 = -10*r2**4*t10**2

            disp = par%s6*t6 + par%s8*r4r2ij*t8 &
               &  + par%s10*49.0_wp/40.0_wp*r4r2ij**2*t10
            ddisp= par%s6*d6 + par%s8*r4r2ij*d8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*d10

            dE = -c6(iat, jat)*disp * 0.5_wp
            dG = -c6(iat, jat)*ddisp*rij
            dS = spread(dG, 1, 3) * spread(rij, 2, 3) * 0.5_wp

            energies(iat) = energies(iat) + dE
            dEdcn(iat) = dEdcn(iat) - dc6dcn(iat, jat) * disp
            sigma = sigma + dS
            if (iat /= jat) then
               energies(jat) = energies(jat) + dE
               dEdcn(jat) = dEdcn(jat) - dc6dcn(jat, iat) * disp
               gradient(:, iat) = gradient(:, iat) + dG
               gradient(:, jat) = gradient(:, jat) - dG
               sigma = sigma + dS
            end if

         end do
      end do
   end do
   !$omp end parallel do

   call dgemv('n', 3*nat, nat, 1.0_wp, dcndr, 3*nat, dEdcn, 1, 1.0_wp, gradient, 1)
   call dgemv('n', 9, nat, 1.0_wp, dcndL, 9, dEdcn, 1, 1.0_wp, sigma, 1)

   energy = sum(energies)

end subroutine d3_gradient_latp


subroutine d3_gradient_neigh &
      & (mol, neighs, neighlist, par, weighting_factor, r4r2, &
      &  cn, dcndr, dcndL, energy, gradient, sigma)

   type(TMolecule), intent(in) :: mol
   type(TNeighbourList), intent(in) :: neighlist
   type(dftd_parameter), intent(in) :: par

   integer, intent(in) :: neighs(:)
   real(wp), intent(in) :: weighting_factor
   real(wp), intent(in) :: r4r2(:)
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)

   real(wp), intent(inout) :: energy
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)

   integer :: nat, max_ref
   integer :: iat, jat, ati, atj, ij, img

   real(wp) :: r4r2ij, r0, rij(3), r2, t6, t8, t10, d6, d8, d10
   real(wp) :: dE, dG(3), dS(3, 3), disp, ddisp

   real(wp), allocatable :: gw(:, :), dgwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: energies(:), dEdcn(:)

   nat = len(mol)
   max_ref = maxval(number_of_references(mol%at))
   allocate(gw(max_ref, nat), dgwdcn(max_ref, nat), c6(nat, nat), &
      &     dc6dcn(nat, nat), energies(nat), dEdcn(nat), source=0.0_wp)

   call weight_references(nat, mol%at, weighting_factor, cn, gw, dgwdcn)

   call get_atomic_c6(nat, mol%at, gw, dgwdcn, c6, dc6dcn)

   !$omp parallel do default(none) &
   !$omp reduction(+:energies, gradient, sigma, dEdcn) &
   !$omp shared(mol, neighs, neighlist, par, r4r2, c6, dc6dcn) &
   !$omp private(ij, img, jat, ati, atj, r2, rij, r4r2ij, r0, t6, t8, t10, &
   !$omp&        d6, d8, d10, disp, ddisp, dE, dG, dS)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      do ij = 1, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dist2(ij, iat)
         rij = mol%xyz(:, iat) - neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)

         r4r2ij = 3*r4r2(ati)*r4r2(atj)
         r0 = par%a1*sqrt(r4r2ij) + par%a2

         t6 = 1._wp/(r2**3+r0**6)
         t8 = 1._wp/(r2**4+r0**8)
         t10 = 1._wp/(r2**5+r0**10)

         d6 = -6*r2**2*t6**2
         d8 = -8*r2**3*t8**2
         d10 = -10*r2**4*t10**2

         disp = par%s6*t6 + par%s8*r4r2ij*t8 &
            &  + par%s10*49.0_wp/40.0_wp*r4r2ij**2*t10
         ddisp= par%s6*d6 + par%s8*r4r2ij*d8 &
            & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*d10

         dE = -c6(iat, jat)*disp * 0.5_wp
         dG = -c6(iat, jat)*ddisp*rij
         dS = spread(dG, 1, 3) * spread(rij, 2, 3) * 0.5_wp

         energies(iat) = energies(iat) + dE
         dEdcn(iat) = dEdcn(iat) - dc6dcn(iat, jat) * disp
         sigma = sigma + dS
         if (iat /= jat) then
            energies(jat) = energies(jat) + dE
            dEdcn(jat) = dEdcn(jat) - dc6dcn(jat, iat) * disp
            gradient(:, iat) = gradient(:, iat) + dG
            gradient(:, jat) = gradient(:, jat) - dG
            sigma = sigma + dS
         endif

      enddo
   enddo
   !$omp end parallel do

   call dgemv('n', 3*nat, nat, 1.0_wp, dcndr, 3*nat, dEdcn, 1, 1.0_wp, gradient, 1)
   call dgemv('n', 9, nat, 1.0_wp, dcndL, 9, dEdcn, 1, 1.0_wp, sigma, 1)

   energy = sum(energies)

end subroutine d3_gradient_neigh


real(wp) pure elemental function weight_cn(wf,cn,cnref) result(cngw)
   real(wp),intent(in) :: wf, cn, cnref
   intrinsic :: exp
   cngw = exp ( -wf * ( cn - cnref )**2 )
end function weight_cn


real(wp) pure elemental function pair_scale(iat, jat) result(scale)
   integer, intent(in) :: iat, jat
   if (iat == jat) then
      scale = 0.5_wp
   else
      scale = 1.0_wp
   endif
end function pair_scale


end module xtb_disp_dftd3
