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
   use xtb_mctc_la, only : contract
   use xtb_disp_dftd3param
   use xtb_param_sqrtzr4r2, only : sqrtZr4r2
   use xtb_type_molecule
   use xtb_type_neighbourlist
   use xtb_type_param
   implicit none
   private

   public :: d3_gradient, d3_atm_gradient


   interface d3_gradient
      module procedure :: d3_full_gradient_neigh
      module procedure :: d3_full_gradient_latp
      module procedure :: d3_gradient_neigh
      module procedure :: d3_gradient_latp
   end interface d3_gradient


   interface d3_atm_gradient
      module procedure :: d3_atm_gradient_neigh
      module procedure :: d3_atm_gradient_latp
   end interface d3_atm_gradient


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
      if (norm > 1e-80_wp) then
         norm = 1.0_wp / norm
      else
         norm = 0.0_wp
      end if
      do iref = 1, number_of_references(ati)
         expw = weight_cn(wf, cn(iat), reference_cn(iref, ati))
         expd = 2*wf*(reference_cn(iref, ati) - cn(iat)) * expw

         gwk = expw * norm
         if (norm == 0.0_wp) then
            if (abs(maxval(reference_cn(:number_of_references(ati), ati)) &
               & - reference_cn(iref, ati)) < 1e-12_wp) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         gwvec(iref, iat) = gwk

         dgwk = expd*norm-expw*dnorm*norm**2
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


subroutine d3_full_gradient_latp &
      & (mol, trans, par, weighting_factor, cutoff, cutoff3, &
      &  cn, dcndr, dcndL, energy, gradient, sigma, e2, e3)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Damping parameters
   type(dftd_parameter), intent(in) :: par

   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: weighting_factor
   real(wp), intent(in) :: cutoff
   real(wp), intent(in) :: cutoff3
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)

   real(wp), intent(inout) :: energy
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(out), optional :: e2
   real(wp), intent(out), optional :: e3

   integer :: nat, max_ref
   real(wp), allocatable :: gw(:, :), dgwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: energies(:), energies3(:), dEdcn(:)

   nat = len(mol)
   max_ref = maxval(number_of_references(mol%at))
   allocate(gw(max_ref, nat), dgwdcn(max_ref, nat), c6(nat, nat), &
      & dc6dcn(nat, nat), energies(nat), energies3(nat), &
      & dEdcn(nat), source=0.0_wp)

   call weight_references(nat, mol%at, weighting_factor, cn, gw, dgwdcn)

   call get_atomic_c6(nat, mol%at, gw, dgwdcn, c6, dc6dcn)

   call disp_gradient_latp(mol, trans, cutoff, par, sqrtZr4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)

   if (present(e2)) e2 = sum(energies)
   if (par%s9 /= 0.0_wp) then
      call atm_gradient_latp(mol, trans, cutoff3, par, sqrtZr4r2, c6, dc6dcn, &
         &  energies3, gradient, sigma, dEdcn)
   end if
   if (present(e3)) e3 = sum(energies3)

   call contract(dcndr, dEdcn, gradient, beta=1.0_wp)
   call contract(dcndL, dEdcn, sigma, beta=1.0_wp)

   energy = energy + sum(energies) + sum(energies3)

end subroutine d3_full_gradient_latp


subroutine d3_gradient_latp &
      & (mol, trans, par, weighting_factor, cutoff, &
      &  cn, dcndr, dcndL, energy, gradient, sigma)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Damping parameters
   type(dftd_parameter), intent(in) :: par

   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: weighting_factor
   real(wp), intent(in) :: cutoff
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)

   real(wp), intent(inout) :: energy
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)

   integer :: nat, max_ref
   real(wp), allocatable :: gw(:, :), dgwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: energies(:), dEdcn(:)

   nat = len(mol)
   max_ref = maxval(number_of_references(mol%at))
   allocate(gw(max_ref, nat), dgwdcn(max_ref, nat), c6(nat, nat), &
      &     dc6dcn(nat, nat), energies(nat), dEdcn(nat), source=0.0_wp)

   call weight_references(nat, mol%at, weighting_factor, cn, gw, dgwdcn)

   call get_atomic_c6(nat, mol%at, gw, dgwdcn, c6, dc6dcn)

   call disp_gradient_latp(mol, trans, cutoff, par, sqrtZr4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)

   call contract(dcndr, dEdcn, gradient, beta=1.0_wp)
   call contract(dcndL, dEdcn, sigma, beta=1.0_wp)

   energy = energy + sum(energies)

end subroutine d3_gradient_latp


!> Implementation of the pairwise dispersion energy
subroutine disp_gradient_latp &
      & (mol, trans, cutoff, par, r4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Damping parameters
   type(dftd_parameter), intent(in) :: par

   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: cutoff
   real(wp), intent(in) :: r4r2(:)
   real(wp), intent(in) :: c6(:, :)
   real(wp), intent(in) :: dc6dcn(:, :)

   real(wp), intent(inout) :: energies(:)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(inout) :: dEdcn(:)

   integer :: iat, jat, ati, atj, itr

   real(wp) :: cutoff2
   real(wp) :: r4r2ij, r0, rij(3), r2, t6, t8, t10, d6, d8, d10
   real(wp) :: dE, dG(3), dS(3, 3), disp, ddisp

   cutoff2 = cutoff**2
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

end subroutine disp_gradient_latp


subroutine d3_atm_gradient_latp &
      & (mol, trans, par, weighting_factor, cutoff, &
      &  cn, dcndr, dcndL, energy, gradient, sigma)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Damping parameters
   type(dftd_parameter), intent(in) :: par

   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: weighting_factor
   real(wp), intent(in) :: cutoff
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)

   real(wp), intent(inout) :: energy
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)

   integer :: nat, max_ref
   real(wp), allocatable :: gw(:, :), dgwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: energies(:), dEdcn(:)

   nat = len(mol)
   max_ref = maxval(number_of_references(mol%at))
   allocate(gw(max_ref, nat), dgwdcn(max_ref, nat), c6(nat, nat), &
      &     dc6dcn(nat, nat), energies(nat), dEdcn(nat), source=0.0_wp)

   call weight_references(nat, mol%at, weighting_factor, cn, gw, dgwdcn)

   call get_atomic_c6(nat, mol%at, gw, dgwdcn, c6, dc6dcn)

   call atm_gradient_latp(mol, trans, cutoff, par, sqrtZr4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)

   call contract(dcndr, dEdcn, gradient, beta=1.0_wp)
   call contract(dcndL, dEdcn, sigma, beta=1.0_wp)

   energy = energy + sum(energies)

end subroutine d3_atm_gradient_latp


subroutine atm_gradient_latp &
      & (mol, trans, cutoff, par, r4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Damping parameters
   type(dftd_parameter), intent(in) :: par

   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: r4r2(:)
   real(wp), intent(in) :: cutoff
   real(wp), intent(in) :: c6(:, :)
   real(wp), intent(in) :: dc6dcn(:, :)

   real(wp), intent(inout) :: energies(:)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(inout) :: dEdcn(:)

   integer :: iat, jat, kat, ati, atj, atk, jtr, ktr
   real(wp) :: cutoff2
   real(wp) :: rij(3), rjk(3), rik(3), r2ij, r2jk, r2ik
   real(wp) :: c6ij, c6jk, c6ik, cij, cjk, cik, scale
   real(wp) :: dE, dG(3, 3), dS(3, 3), dCN(3)
   real(wp), parameter :: sr = 4.0_wp/3.0_wp

   cutoff2 = cutoff**2

   do iat = 1, len(mol)
      ati = mol%at(iat)
      do jat = 1, iat
         atj = mol%at(jat)

         c6ij = c6(jat,iat)
         !cij = par%a1*sqrt(3.0_wp*r4r2(ati)*r4r2(atj))+par%a2
         cij = sr*get_vdwrad(ati, atj)

         do kat = 1, jat
            atk = mol%at(kat)

            c6ik = c6(kat,iat)
            c6jk = c6(kat,jat)

            !cik = par%a1*sqrt(3.0_wp*r4r2(ati)*r4r2(atk))+par%a2
            cik = sr*get_vdwrad(ati, atk)
            !cjk = par%a1*sqrt(3.0_wp*r4r2(atj)*r4r2(atk))+par%a2
            cjk = sr*get_vdwrad(atj, atk)

            do jtr = 1, size(trans, dim=2)
               rij = mol%xyz(:, jat) - mol%xyz(:, iat) + trans(:, jtr)
               r2ij = sum(rij**2)
               if (r2ij > cutoff2 .or. r2ij < 1.0e-14_wp) cycle
               do ktr = 1, size(trans, dim=2)
                  if (jat == kat .and. jtr == ktr) cycle
                  rik = mol%xyz(:, kat) - mol%xyz(:, iat) + trans(:, ktr)
                  r2ik = sum(rik**2)
                  if (r2ik > cutoff2 .or. r2ik < 1.0e-14_wp) cycle
                  rjk = mol%xyz(:, kat) - mol%xyz(:, jat) + trans(:, ktr) &
                     & - trans(:, jtr)
                  r2jk = sum(rjk**2)
                  if (r2jk > cutoff2 .or. r2jk < 1.0e-14_wp) cycle

                  call deriv_atm_triple(c6ij, c6ik, c6jk, cij, cjk, cik, &
                     & r2ij, r2jk, r2ik, dc6dcn(iat,jat), dc6dcn(jat,iat), &
                     & dc6dcn(jat,kat), dc6dcn(kat,jat), dc6dcn(iat,kat), &
                     & dc6dcn(kat,iat), rij, rjk, rik, par%alp, dE, dG, dS, dCN)

                  scale = par%s9 * triple_scale(iat, jat, kat)
                  energies(iat) = energies(iat) + dE * scale/3
                  energies(jat) = energies(jat) + dE * scale/3
                  energies(kat) = energies(kat) + dE * scale/3
                  gradient(:, iat) = gradient(:, iat) + dG(:, 1) * scale
                  gradient(:, jat) = gradient(:, jat) + dG(:, 2) * scale
                  gradient(:, kat) = gradient(:, kat) + dG(:, 3) * scale
                  sigma(:, :) = sigma + dS * scale
                  dEdcn(iat) = dEdcn(iat) + dCN(1) * scale
                  dEdcn(jat) = dEdcn(jat) + dCN(2) * scale
                  dEdcn(kat) = dEdcn(kat) + dCN(3) * scale

               end do
            end do

         end do
      end do
   end do

end subroutine atm_gradient_latp


subroutine d3_full_gradient_neigh &
      & (mol, neighs, neighs3, neighlist, par, weighting_factor, &
      &  cn, dcndr, dcndL, energy, gradient, sigma, e2, e3)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Neighbour list
   type(TNeighbourList), intent(in) :: neighlist

   !> Damping parameters
   type(dftd_parameter), intent(in) :: par

   integer, intent(in) :: neighs(:)
   integer, intent(in) :: neighs3(:)
   real(wp), intent(in) :: weighting_factor
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)

   real(wp), intent(inout) :: energy
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(out), optional :: e2
   real(wp), intent(out), optional :: e3

   integer :: nat, max_ref
   real(wp), allocatable :: gw(:, :), dgwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: energies(:), energies3(:), dEdcn(:)

   nat = len(mol)
   max_ref = maxval(number_of_references(mol%at))
   allocate(gw(max_ref, nat), dgwdcn(max_ref, nat), c6(nat, nat), &
      & dc6dcn(nat, nat), energies(nat), energies3(nat), &
      & dEdcn(nat), source=0.0_wp)

   call weight_references(nat, mol%at, weighting_factor, cn, gw, dgwdcn)

   call get_atomic_c6(nat, mol%at, gw, dgwdcn, c6, dc6dcn)

   call disp_gradient_neigh(mol, neighs, neighlist, par, sqrtZr4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)

   if (present(e2)) e2 = sum(energies)
   if (par%s9 /= 0.0_wp) then
      call atm_gradient_neigh(mol, neighs, neighlist, par, sqrtZr4r2, c6, dc6dcn, &
         & energies3, gradient, sigma, dEdcn)
   end if
   if (present(e3)) e3 = sum(energies3)

   call contract(dcndr, dEdcn, gradient, beta=1.0_wp)
   call contract(dcndL, dEdcn, sigma, beta=1.0_wp)

   energy = energy + sum(energies) + sum(energies3)

end subroutine d3_full_gradient_neigh


subroutine d3_gradient_neigh &
      & (mol, neighs, neighlist, par, weighting_factor, &
      &  cn, dcndr, dcndL, energy, gradient, sigma)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Neighbour list
   type(TNeighbourList), intent(in) :: neighlist

   !> Damping parameters
   type(dftd_parameter), intent(in) :: par

   integer, intent(in) :: neighs(:)
   real(wp), intent(in) :: weighting_factor
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)

   real(wp), intent(inout) :: energy
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)

   integer :: nat, max_ref
   real(wp), allocatable :: gw(:, :), dgwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: energies(:), dEdcn(:)

   nat = len(mol)
   max_ref = maxval(number_of_references(mol%at))
   allocate(gw(max_ref, nat), dgwdcn(max_ref, nat), c6(nat, nat), &
      &     dc6dcn(nat, nat), energies(nat), dEdcn(nat), source=0.0_wp)

   call weight_references(nat, mol%at, weighting_factor, cn, gw, dgwdcn)

   call get_atomic_c6(nat, mol%at, gw, dgwdcn, c6, dc6dcn)

   call disp_gradient_neigh(mol, neighs, neighlist, par, sqrtZr4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)

   call contract(dcndr, dEdcn, gradient, beta=1.0_wp)
   call contract(dcndL, dEdcn, sigma, beta=1.0_wp)

   energy = energy + sum(energies)

end subroutine d3_gradient_neigh


subroutine disp_gradient_neigh &
      & (mol, neighs, neighlist, par, r4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Neighbour list
   type(TNeighbourList), intent(in) :: neighlist

   !> Damping parameters
   type(dftd_parameter), intent(in) :: par

   integer, intent(in) :: neighs(:)
   real(wp), intent(in) :: r4r2(:)
   real(wp), intent(in) :: c6(:, :)
   real(wp), intent(in) :: dc6dcn(:, :)

   real(wp), intent(inout) :: energies(:)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(inout) :: dEdcn(:)

   integer :: iat, jat, ati, atj, ij, img
   real(wp) :: r4r2ij, r0, rij(3), r2, t6, t8, t10, d6, d8, d10
   real(wp) :: dE, dG(3), dS(3, 3), disp, ddisp

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

end subroutine disp_gradient_neigh


subroutine d3_atm_gradient_neigh &
      & (mol, neighs, neighlist, par, weighting_factor, &
      &  cn, dcndr, dcndL, energy, gradient, sigma)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Neighbour list
   type(TNeighbourList), intent(in) :: neighlist

   !> Damping parameters
   type(dftd_parameter), intent(in) :: par

   integer, intent(in) :: neighs(:)
   real(wp), intent(in) :: weighting_factor
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

   call atm_gradient_neigh(mol, neighs, neighlist, par, sqrtZr4r2, c6, dc6dcn, &
      & energies, gradient, sigma, dEdcn)

   call contract(dcndr, dEdcn, gradient, beta=1.0_wp)
   call contract(dcndL, dEdcn, sigma, beta=1.0_wp)

   energy = energy + sum(energies)

end subroutine d3_atm_gradient_neigh


subroutine atm_gradient_neigh &
      & (mol, neighs, neighlist, par, r4r2, c6, dc6dcn, &
      &  energies, gradient, sigma, dEdcn)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Neighbour list
   type(TNeighbourList), intent(in) :: neighlist

   !> Damping parameters
   type(dftd_parameter), intent(in) :: par

   integer, intent(in) :: neighs(:)
   real(wp), intent(in) :: r4r2(:)
   real(wp), intent(in) :: c6(:, :)
   real(wp), intent(in) :: dc6dcn(:, :)

   real(wp), intent(inout) :: energies(:)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(inout) :: dEdcn(:)

   integer :: iat, jat, kat, ati, atj, atk, jtr, ktr, ij, jk, ik
   real(wp) :: rij(3), rjk(3), rik(3), r2ij, r2jk, r2ik
   real(wp) :: c6ij, c6jk, c6ik, cij, cjk, cik, scale
   real(wp) :: dE, dG(3, 3), dS(3, 3), dCN(3)
   real(wp), parameter :: sr = 4.0_wp/3.0_wp

   do iat = 1, len(mol)
      ati = mol%at(iat)
      do ij = 1, neighs(iat)
         jtr = neighlist%ineigh(ij, iat)
         r2ij = neighlist%dist2(ij, iat)
         rij = neighlist%coords(:, jtr) - neighlist%coords(:, iat)
         jat = neighlist%image(jtr)
         atj = mol%at(jat)

         c6ij = c6(jat,iat)
         !cij = par%a1*sqrt(3.0_wp*r4r2(ati)*r4r2(atj))+par%a2
         cij = sr*get_vdwrad(ati, atj)

         do ik = 1, ij-1
            ktr = neighlist%ineigh(ik, iat)
            rik = neighlist%coords(:, ktr) - neighlist%coords(:, iat)
            r2ik = neighlist%dist2(ik, iat)
            rjk = neighlist%coords(:, ktr) - neighlist%coords(:, jtr)
            r2jk = sum(rjk**2)
            kat = neighlist%image(ktr)
            atk = mol%at(kat)

            c6ik = c6(kat,iat)
            c6jk = c6(kat,jat)

            !cik = par%a1*sqrt(3.0_wp*r4r2(ati)*r4r2(atk))+par%a2
            cik = sr*get_vdwrad(ati, atk)
            !cjk = par%a1*sqrt(3.0_wp*r4r2(atj)*r4r2(atk))+par%a2
            cjk = sr*get_vdwrad(atj, atk)

            call deriv_atm_triple(c6ij, c6ik, c6jk, cij, cjk, cik, &
               & r2ij, r2jk, r2ik, dc6dcn(iat,jat), dc6dcn(jat,iat), &
               & dc6dcn(jat,kat), dc6dcn(kat,jat), dc6dcn(iat,kat), &
               & dc6dcn(kat,iat), rij, rjk, rik, par%alp, dE, dG, dS, dCN)

            scale = par%s9 * triple_scale(iat, jat, kat)
            energies(iat) = energies(iat) + dE * scale/3
            energies(jat) = energies(jat) + dE * scale/3
            energies(kat) = energies(kat) + dE * scale/3
            gradient(:, iat) = gradient(:, iat) + dG(:, 1) * scale
            gradient(:, jat) = gradient(:, jat) + dG(:, 2) * scale
            gradient(:, kat) = gradient(:, kat) + dG(:, 3) * scale
            sigma(:, :) = sigma + dS * scale
            dEdcn(iat) = dEdcn(iat) + dCN(1) * scale
            dEdcn(jat) = dEdcn(jat) + dCN(2) * scale
            dEdcn(kat) = dEdcn(kat) + dCN(3) * scale

         end do
      end do
   end do

end subroutine atm_gradient_neigh


pure subroutine deriv_atm_triple(c6ij, c6ik, c6jk, cij, cjk, cik, &
      & r2ij, r2jk, r2ik, dc6ij, dc6ji, dc6jk, dc6kj, dc6ik, dc6ki, &
      & rij, rjk, rik, alp, dE, dG, dS, dCN)

   real(wp), intent(in) :: c6ij, c6ik, c6jk
   real(wp), intent(in) :: cij, cjk, cik
   real(wp), intent(in) :: r2ij, r2jk, r2ik
   real(wp), intent(in) :: dc6ij, dc6ji, dc6jk, dc6kj, dc6ik, dc6ki
   real(wp), intent(in) :: rij(3), rjk(3), rik(3)
   integer, intent(in) :: alp
   real(wp), intent(out) :: dE, dG(3, 3), dS(3, 3), dCN(3)

   real(wp) :: c9, dc9, ccc1, rrr1, rrr2, rrr3, ang, dang, fdmp, dfdmp, dGr, cr

   c9 = -sqrt(c6ij*c6ik*c6jk)

   ccc1 = cij*cjk*cik

   rrr2 = r2ij*r2jk*r2ik
   rrr1 = sqrt(rrr2)
   rrr3 = rrr1*rrr2

   ang = 0.375_wp * (r2ij+r2jk-r2ik)*(r2ij-r2jk+r2ik)*(-r2ij+r2jk+r2ik) &
      & / (rrr3*rrr2) + 1.0_wp/(rrr3)

   cr = (ccc1/rrr1)**(1.0_wp/3.0_wp)
   fdmp = 1.0_wp/(1.0_wp + 6.0_wp*cr**alp)
   dfdmp = -(2.0_wp*alp*cr**alp) * fdmp**2

   ! Energy contribution
   dE = -fdmp*ang*c9

   ! Derivative w.r.t. i-j distance
   dang = -0.375_wp*(r2ij**3+r2ij**2*(r2jk+r2ik) &
      & +r2ij*(3.0_wp*r2jk**2+2.0_wp*r2jk*r2ik+3.0_wp*r2ik**2) &
      & -5.0_wp*(r2jk-r2ik)**2*(r2jk+r2ik)) / (rrr3*rrr2)
   dGr = (-dang*c9*fdmp + dfdmp*c9*ang)/r2ij
   dG(:, 1) = -dGr * rij
   dG(:, 2) = +dGr * rij 
   dS(:, :) = 0.5_wp * dGr * spread(rij, 1, 3) * spread(rij, 2, 3)

   ! Derivative w.r.t. i-k distance
   dang = -0.375_wp*(r2ik**3+r2ik**2*(r2jk+r2ij) &
      & +r2ik*(3.0_wp*r2jk**2+2.0*r2jk*r2ij+3.0_wp*r2ij**2) &
      & -5.0_wp*(r2jk-r2ij)**2*(r2jk+r2ij)) / (rrr3*rrr2)
   dGr = (-dang*c9*fdmp + dfdmp*c9*ang)/r2ik
   dG(:, 1) = -dGr * rik + dG(:, 1)
   dG(:, 3) = +dGr * rik 
   dS(:, :) = 0.5_wp * dGr * spread(rik, 1, 3) * spread(rik, 2, 3) + dS

   ! Derivative w.r.t. j-k distance
   dang=-0.375_wp*(r2jk**3+r2jk**2*(r2ik+r2ij) &
      & +r2jk*(3.0_wp*r2ik**2+2.0_wp*r2ik*r2ij+3.0_wp*r2ij**2) &
      & -5.0_wp*(r2ik-r2ij)**2*(r2ik+r2ij)) / (rrr3*rrr2)
   dGr = (-dang*c9*fdmp + dfdmp*c9*ang)/r2jk
   dG(:, 2) = -dGr * rjk + dG(:, 2)
   dG(:, 3) = +dGr * rjk + dG(:, 3)
   dS(:, :) = 0.5_wp * dGr * spread(rjk, 1, 3) * spread(rjk, 2, 3) + dS

   ! CN derivative
   dc9 = 0.5_wp*c9*(dc6ij/c6ij+dc6ik/c6ik)
   dCN(1) = -ang*fdmp*dc9
   dc9 = 0.5_wp*c9*(dc6ji/c6ij+dc6jk/c6jk)
   dCN(2) = -ang*fdmp*dc9
   dc9 = 0.5_wp*c9*(dc6ki/c6ik+dc6kj/c6jk)
   dCN(3) = -ang*fdmp*dc9

end subroutine deriv_atm_triple


elemental function weight_cn(wf,cn,cnref) result(cngw)
   real(wp),intent(in) :: wf, cn, cnref
   real(wp) :: cngw
   real(wp) :: val
   intrinsic :: exp

   val = -wf * ( cn - cnref )**2
   if (val < -200.0_wp) then ! technically, exp(-200) -> 1.383897e-87
     cngw = 0.0_wp
   else
     cngw = exp ( val )
   end if

end function weight_cn


!> Correct counting of energy contributions for atomic pairs.
elemental function pair_scale(iat, jat) result(scale)

   !> Atom indices
   integer, intent(in) :: iat, jat

   !> Fraction of energy
   real(wp) :: scale

   if (iat == jat) then
      scale = 0.5_wp
   else
      scale = 1.0_wp
   endif

end function pair_scale


!> Logic exercise to distribute a triple energy to atomwise energies.
elemental function triple_scale(ii, jj, kk) result(scale)

   !> Atom indices
   integer, intent(in) :: ii, jj, kk

   !> Fraction of energy
   real(wp) :: scale

   if (ii == jj) then
      if (ii == kk) then
         ! ii'i" -> 1/6
         scale = 1.0_wp/6.0_wp
      else
         ! ii'j -> 1/2
         scale = 0.5_wp
      end if
   else
      if (ii /= kk .and. jj /= kk) then
         ! ijk -> 1 (full)
         scale = 1.0_wp
      else
         ! ijj' and iji' -> 1/2
         scale = 0.5_wp
      end if
   end if

end function triple_scale


end module xtb_disp_dftd3
