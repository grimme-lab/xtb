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

!> Implements the integral evaluation for the SCF
submodule(scf_module) aes_integrals
   use mctc_constants, only: pi
   use mctc_econv, only: evtoau

   use tbdef_molecule
   use tbdef_neighbourlist
   use tbdef_basisset
   use tbdef_param

   use aoparam

   use scc_core
   use grad_core
   use intgrad

   implicit none

contains

!> Computes the overlap integrals and the overlap dependent core Hamiltonian.
!
!  Also calculates dipole integrals for evalulation of the dipole moment later.
module subroutine build_SDH0(mol, neighs, neighlist, basis, param, intcut, &
      &                      cn,  kcnsh, sint, dpint, H0)
   !> Molecular structure information.
   type(tb_molecule), intent(in) :: mol
   !> Static neighbourlist.
   type(tb_neighbourlist), intent(in) :: neighlist
   !> Tight-binding basis set.
   type(tb_basisset), intent(in) :: basis
   !> Global parameters.
   type(scc_parameter), intent(in) :: param
   !> Number of neighbours for each atom.
   integer, intent(in) :: neighs(:)
   !> Integral cutoff according to prefactor from Gaussian product theorem.
   real(wp), intent(in) :: intcut
   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)
   !> Shellwise scaling parameter for the atomic levels (in eV).
   real(wp), intent(in) :: kcnsh(:)
   !> Core Hamiltonian (packed).
   real(wp), intent(out) :: H0(:)
   !> Overlap integrals.
   real(wp), intent(out) :: sint(:,:)
   !> Dipole integrals.
   real(wp), intent(out) :: dpint(:,:,:)
   !> Temporary unpacked core Hamiltonian.
   real(wp), allocatable :: h(:,:)

   integer :: i, j, k, l, ii, jj, ij, ijao, img, io, jo
   integer :: iat, jat, ati, atj, ish, jsh, icao, jcao, iao, jao, jshmax
   integer :: il, jl, iptyp, jptyp, naoi, naoj
   real(wp) :: r1, r2, ri(3), rj(3)
   real(wp) :: dtmp(3), qtmp(6), ss(6,6), dd(3,6,6), qq(6,6,6), tmp(6,6)
   real(wp) :: hii, hjj, km, hav
   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)
   integer, parameter :: itt(0:3) = (/0,1,4,10/)
   real(wp), parameter :: point(3) = 0.0_wp

   allocate(h(basis%nao, basis%nao), source=0.0_wp)

   !call timing(t1,t3)
   sint=0.0_wp
   dpint=0.0_wp
   H0=0.0_wp
   H = 0.0_wp

   ! diagonal elements
   do iao = 1, basis%nao
      iat = basis%aoat2(iao)
      ati = mol%at(iat)
      ish = basis%ao2sh(iao)
      il  = basis%lsh(ish)+1

      ! calculate environment dependent shift
      hii = basis%level(ish) + kcnsh(ish)*cn(iat)
      H(iao,iao) = hii
   end do

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:sint, dpint, H) shared(mol, neighs, neighlist, intcut) &
   !$omp shared(basis, param, ao_n, ao_l, cn, kcnsh) &
   !$omp private(ri, rj, r2, jat, ish, jsh, ati, atj, icao, jcao, ij, img, &
   !$omp&        naoi, naoj, jshmax, il, jl, iptyp, jptyp, iao, jao, &
   !$omp&        io, jo, hii, hjj, km, hav, ss, dd, qq, tmp)
   do iat=1, len(mol)
      ati = mol%at(iat)
      ri = mol%xyz(:, iat)
      io = basis%shells(1, iat)-1
      do ij = 0, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dists2(ij, iat)
         rj = neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         if (ij > 0 .and. iat == jat) cycle
         jo = basis%shells(1, jat)-1
         ishells: do ish = 1, ao_n(ati)
            il = basis%lsh(ish+io)+1
            icao = basis%caoshell(ish, iat)
            naoi = llao(il-1)
            iptyp = itt(il-1)
            jshmax = ao_n(atj)
            if(iat == jat) jshmax = ish
            jshells: do jsh = 1, jshmax
               jl = basis%lsh(jsh+jo)+1
               jcao = basis%caoshell(jsh, jat)
               naoj = llao(jl-1)
               jptyp = itt(jl-1)

               hii = basis%level(ish+io) + kcnsh(ish+io)*cn(iat)
               hjj = basis%level(jsh+jo) + kcnsh(jsh+jo)*cn(jat)
               km = h0scal(il, jl, ati, atj, basis%valsh(ish+io)/=0, &
                  &        basis%valsh(jsh+jo)/=0, param%kspd, param%kmagic, &
                  &        param%kenscal)

               hav = 0.5_wp * (hii + hjj) * km * rfactor(il, jl, ati, atj, ri, rj)

               ss = 0.0_wp
               dd = 0.0_wp
               qq = 0.0_wp
               call get_multiints(icao, jcao, naoi, naoj, iptyp, jptyp, ri, rj, &
                  &               point, intcut, basis%nprim, basis%primcount, &
                  &               basis%alp, basis%cont, ss, dd, qq)
               ! transform from CAO to SAO
               call dtrf2(ss, il-1, jl-1)
               do k = 1, 3
                  tmp(1:6, 1:6) = dd(k, 1:6, 1:6)
                  call dtrf2(tmp, il-1, jl-1)
                  dd(k, 1:6, 1:6) = tmp(1:6, 1:6)
               enddo
               do ii = 1, llao2(il-1)
                  iao = ii + basis%saoshell(ish, iat)
                  do jj = 1, llao2(jl-1)
                     jao = jj + basis%saoshell(jsh, jat)
                     sint(jao, iao) = sint(jao, iao) + ss(jj, ii)
                     dpint(:, jao, iao) = dpint(:, jao, iao) + dd(:, jj, ii)
                     if (ij > 0) H(jao, iao) = H(jao, iao) + hav*ss(jj, ii)
                     if(jao >= iao .and. iat == jat) cycle
                     sint(iao, jao) = sint(iao, jao) + ss(jj, ii)
                     dpint(:, iao, jao) = dpint(:, iao, jao) + dd(:, jj, ii)
                     if (ij > 0) H(iao, jao) = H(iao, jao) + hav*ss(jj, ii)
                  enddo
               enddo
            enddo jshells
         enddo ishells
      enddo
   enddo
   !$omp end parallel do

   do iao = 1, basis%nao
      do jao = 1, iao
         H0(jao+iao*(iao-1)/2) = H(jao,iao)
      enddo
   enddo

end subroutine build_SDH0

!> Computes the overlap integrals and the overlap dependent core Hamiltonian.
!
!  Also calculates dipole integrals for evalulation of the dipole moment later.
module subroutine build_SDQH0(mol, neighs, neighlist, basis, param, intcut, cn, &
      &                       sint, dpint, qpint, H0)
   !> Molecular structure information.
   type(tb_molecule), intent(in) :: mol
   !> Static neighbourlist.
   type(tb_neighbourlist), intent(in) :: neighlist
   !> Tight-binding basis set.
   type(tb_basisset), intent(in) :: basis
   !> Global parameters.
   type(scc_parameter), intent(in) :: param
   !> Number of neighbours for each atom.
   integer, intent(in) :: neighs(:)
   !> Integral cutoff according to prefactor from Gaussian product theorem.
   real(wp), intent(in) :: intcut
   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)
   !> Core Hamiltonian (packed).
   real(wp), intent(out) :: H0(:)
   !> Overlap integrals.
   real(wp), intent(out) :: sint(:,:)
   !> Dipole integrals.
   real(wp), intent(out) :: dpint(:,:,:)
   !> Quadrupole integrals.
   real(wp), intent(out) :: qpint(:,:,:)
   !> Temporary unpacked core Hamiltonian.
   real(wp), allocatable :: h(:,:)

   integer :: i, j, k, l, ii, jj, ij, ijao, img, io, jo
   integer :: iat, jat, ati, atj, ish, jsh, icao, jcao, iao, jao, jshmax
   integer :: il, jl, iptyp, jptyp, naoi, naoj
   real(wp) :: r1, r2, ri(3), rj(3), zi, zj
   real(wp) :: dtmp(3), qtmp(6), ss(6,6), dd(3,6,6), qq(6,6,6), tmp(6,6)
   real(wp) :: hii, hjj, km, hav
   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)
   integer, parameter :: itt(0:3) = (/0,1,4,10/)
   real(wp), parameter :: point(3) = 0.0_wp

   allocate(h(basis%nao, basis%nao), source=0.0_wp)

   !call timing(t1,t3)
   sint=0.0_wp
   dpint=0.0_wp
   H0=0.0_wp
   H = 0.0_wp

   ! diagonal elements
   do iao = 1, basis%nao
      iat = basis%aoat2(iao)
      ati = mol%at(iat)
      ish = basis%ao2sh(iao)
      il  = basis%lsh(ish)+1

      ! calculate environment dependent shift
      hii = basis%level(ish) - kcnat(il-1, ati)*cn(iat)
      H(iao,iao) = hii
   end do

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:sint, dpint, qpint, H) shared(mol, neighs, neighlist) &
   !$omp shared(basis, param, ao_n, ao_l, cn, kcnat, intcut) &
   !$omp private(ri, rj, r2, jat, ish, jsh, ati, atj, icao, jcao, ij, img, &
   !$omp&        naoi, naoj, jshmax, il, jl, iptyp, jptyp, iao, jao, &
   !$omp&        io, jo, zi, zj, hii, hjj, km, hav, ss, dd, qq, tmp)
   do iat=1, len(mol)
      ati = mol%at(iat)
      ri = mol%xyz(:, iat)
      io = basis%shells(1, iat)-1
      do ij = 0, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dists2(ij, iat)
         rj = neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         jo = basis%shells(1, jat)-1
         ishells: do ish = 1, ao_n(ati)
            il = basis%lsh(ish+io)+1
            icao = basis%caoshell(ish, iat)
            naoi = llao(il-1)
            iptyp = itt(il-1)
            jshmax = ao_n(atj)
            if(iat == jat) jshmax = ish
            jshells: do jsh = 1, jshmax
               jl = basis%lsh(jsh+jo)+1
               jcao = basis%caoshell(jsh, jat)
               naoj = llao(jl-1)
               jptyp = itt(jl-1)

               hii = basis%level(ish+io) - kcnat(il-1, ati)*cn(iat)
               hjj = basis%level(jsh+jo) - kcnat(jl-1, atj)*cn(jat)

               zi = basis%zeta(ish+io)
               zj = basis%zeta(jsh+jo)
               km = h0scal(il, jl, ati, atj, basis%valsh(ish+io)/=0, &
                  &        basis%valsh(jsh+jo)/=0, param%kspd, param%kmagic, &
                  &        param%kenscal) &
                  & * (0.5_wp*(zi+zj)/sqrt(zi*zj))**(-0.5_wp)

               hav = 0.5_wp * (hii + hjj) * km * rfactor(il, jl, ati, atj, ri, rj)

               ss = 0.0_wp
               dd = 0.0_wp
               qq = 0.0_wp
               call get_multiints(icao, jcao, naoi, naoj, iptyp, jptyp, ri, rj, &
                  &               point, intcut, basis%nprim, basis%primcount, &
                  &               basis%alp, basis%cont, ss, dd, qq)
               ! transform from CAO to SAO
               call dtrf2(ss, il-1, jl-1)
               do k = 1, 3
                  tmp(1:6, 1:6) = dd(k, 1:6, 1:6)
                  call dtrf2(tmp, il-1, jl-1)
                  dd(k, 1:6, 1:6) = tmp(1:6, 1:6)
               enddo
               do k = 1, 6
                  tmp(1:6, 1:6) = qq(k, 1:6, 1:6)
                  call dtrf2(tmp, il-1, jl-1)
                  qq(k, 1:6, 1:6) = tmp(1:6, 1:6)
               enddo
               do ii = 1, llao2(il-1)
                  iao = ii + basis%saoshell(ish, iat)
                  do jj = 1, llao2(jl-1)
                     jao = jj + basis%saoshell(jsh, jat)
                     sint(jao, iao) = sint(jao, iao) + ss(jj, ii)
                     dpint(:, jao, iao) = dpint(:, jao, iao) + dd(:, jj, ii)
                     qpint(:, jao, iao) = qpint(:, jao, iao) + qq(:, jj, ii)
                     if (ij > 0) H(jao, iao) = H(jao, iao) + hav*ss(jj, ii)
                     if(jao >= iao .and. iat == jat) cycle
                     sint(iao, jao) = sint(iao, jao) + ss(jj, ii)
                     dpint(:, iao, jao) = dpint(:, iao, jao) + dd(:, jj, ii)
                     qpint(:, iao, jao) = qpint(:, iao, jao) + qq(:, jj, ii)
                     if (ij > 0) H(iao, jao) = H(iao, jao) + hav*ss(jj, ii)
                  enddo
               enddo
            enddo jshells
         enddo ishells
      enddo
   enddo
   !$omp end parallel do

   do iao = 1, basis%nao
      do jao = 1, iao
         H0(jao+iao*(iao-1)/2) = H(jao,iao)
      enddo
   enddo

end subroutine build_SDQH0

!> Computes the gradient of the overlap integral contributions to the Hamiltonian.
module subroutine build_dSDH0(mol, neighs, neighlist, basis, param, intcut, &
      &                       cn, kcnsh, P, Pew, ves, dhdcn, gradient, sigma)
   !> Molecular structure information.
   type(tb_molecule), intent(in) :: mol
   !> Static neighbourlist.
   type(tb_neighbourlist), intent(in) :: neighlist
   !> Tight-binding basis set.
   type(tb_basisset), intent(in) :: basis
   !> Global parameters.
   type(scc_parameter), intent(in) :: param
   !> Number of neighbours for each atom.
   integer, intent(in) :: neighs(:)
   !> Integral cutoff according to prefactor from Gaussian product theorem.
   real(wp), intent(in) :: intcut
   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)
   !> Shellwise scaling parameter for the atomic levels (in eV).
   real(wp), intent(in) :: kcnsh(:)
   !> Density matrix.
   real(wp), intent(in) :: P(:, :)
   !> Energy weighted density matrix.
   real(wp), intent(in) :: Pew(:, :)
   !> Potential on each shell from the charge field (in eV).
   real(wp), intent(in) :: ves(:)
   !> Derivative of the Hamiltonian w.r.t. the coordination number.
   real(wp), intent(out) :: dhdcn(:)
   !> Molecular gradient.
   real(wp), intent(inout) :: gradient(:,:)
   !> Stress tensor (not scaled by volume).
   real(wp), intent(inout) :: sigma(:,:)

   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)
   integer, parameter :: itt(0:3)=(/0,1,4,10/)
   integer :: iat, jat, ati, atj, ish, jsh, icao, jcao, iao, jao, ixyz
   integer :: ii, jj, ij, jshmax, img, io, jo, il, jl
   integer :: ishtyp, jshtyp, iptyp, jptyp, naoi, naoj
   real(wp) :: r2, rij(3), ri(3), rj(3)
   real(wp) :: tmp(6,6), dumdum(3), dum, sdq(6,6), sdqg(3,6,6)
   real(wp) :: dG(3), dS(3,3)
   real(wp) :: hii, hjj, km, hav, Pij, HPij, vij, rfc, dfc(3)

   dhdcn = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:gradient, sigma, dhdcn) &
   !$omp shared(mol, neighs, neighlist, basis, param, cn, kcnsh, P, Pew, ves, &
   !$omp&       ao_n, ao_l, intcut) &
   !$omp private(jat, ij, img, ati, atj, ishtyp, jshtyp, icao, jcao, naoi, naoj, &
   !$omp&        jshmax, iptyp, jptyp, iao, jao, ii, jj, io, jo, il, jl, &
   !$omp&        ixyz, ri, rj, r2, rij, hii, hjj, km, hav, vij, rfc, dfc, &
   !$omp&        Pij, HPij, sdq, sdqg, tmp, dG, dS)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      ri = mol%xyz(:, iat)
      io = basis%shells(1, iat)-1
      do ij = 1, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dists2(ij, iat)
         rj = neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         rij = ri - rj
         if (iat == jat) cycle
         jo = basis%shells(1, jat)-1
         do ish = 1, ao_n(ati)
            il = basis%lsh(ish+io)+1
            icao = basis%caoshell(ish,iat)
            naoi = llao(il-1)
            iptyp = itt(il-1)
            jshmax = ao_n(atj)
            if(iat.eq.jat) jshmax=ish
            do jsh = 1, jshmax
               jl = basis%lsh(jsh+jo)+1
               jcao = basis%caoshell(jsh,jat)
               naoj = llao(jl-1)
               jptyp = itt(jl-1)

               hii = basis%level(ish+io) + kcnsh(ish+io)*cn(iat)
               hjj = basis%level(jsh+jo) + kcnsh(jsh+jo)*cn(jat)
               km = h0scal(il, jl, ati, atj, basis%valsh(ish+io)/=0, &
                  &        basis%valsh(jsh+jo)/=0, param%kspd, param%kmagic, &
                  &        param%kenscal)

               call drfactor(il, jl, 0, ati, atj, r2, ri, rj, rfc, dfc)
               hav = 0.5_wp * (hii + hjj) * km * rfc

               vij = -0.5_wp*(ves(ish+io) + ves(jsh+jo))

               ! we go through the primitives
               ! (because the screening is the same for all of them)
               sdqg = 0.0_wp
               sdq = 0.0_wp
               call get_grad_overlap(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj, &
                  &                  rj,intcut,basis%nprim,basis%primcount, &
                  &                  basis%alp,basis%cont,sdq,sdqg)
               call dtrf2(sdq,il-1,jl-1)
               do ixyz = 1, 3
                  ! transform from CAO to SAO, only transform the gradient
                  tmp(1:6, 1:6)=sdqg(ixyz, 1:6, 1:6)
                  call dtrf2(tmp, il-1, jl-1)
                  sdqg(ixyz, 1:6, 1:6)=tmp(1:6, 1:6)
               enddo
               do ii = 1, llao2(il-1)
                  iao = ii + basis%saoshell(ish,iat)
                  do jj = 1, llao2(jl-1)
                     jao = jj + basis%saoshell(jsh,jat)
                     Pij = P(jao, iao) * evtoau
                     HPij = (hav + vij)*Pij - Pew(jao, iao)
                     dG = sdqg(:, jj, ii) * 2*HPij
                     ! add polynomial gradient
                     dG = dG + 2*hav*dfc/rfc * Pij * sdq(jj,ii)
                     ! generate strain contribution
                     dS = spread(dG, 1, 3) * spread(rij, 2, 3)
                     ! reset for CN derivative
                     HPij = km*rfc*Pij*sdq(jj,ii)
                     if (iat /= jat) then
                        sigma = sigma + dS
                        gradient(:, iat) = gradient(:, iat) + dG
                        gradient(:, jat) = gradient(:, jat) - dG
                        dhdcn(iat) = dhdcn(iat) + HPij*kcnsh(ish+io)
                        dhdcn(jat) = dhdcn(jat) + HPij*kcnsh(jsh+jo)
                     else
                        dhdcn(iat) = dhdcn(iat) + HPij*(kcnsh(ish+io)+kcnsh(jsh+jo))
                        sigma = sigma + 0.5_wp * dS
                     endif
                  enddo
               enddo
            enddo ! jsh : loop over shells on jat
         enddo ! ish : loop over shells on iat
      enddo ! jat
   enddo ! iat
   !$omp end parallel do

   do iao = 1, basis%nao
      iat = basis%aoat2(iao)
      ish = basis%ao2sh(iao)
      dhdcn(iat) = dhdcn(iat)+P(iao,iao)*kcnsh(ish)*evtoau ! diagonal contribution
   enddo

end subroutine build_dSDH0

!> Computes the overlap and multipole integral contributions to the gradient
!  of the Hamiltonian.
module subroutine build_dSDQH0(mol, neighs, neighlist, basis, param, intcut, cn, &
      &                        P, Pew, ves, vs, vd, vq, dhdcn, gradient, sigma)
   !> Molecular structure information.
   type(tb_molecule), intent(in) :: mol
   !> Static neighbourlist.
   type(tb_neighbourlist), intent(in) :: neighlist
   !> Tight-binding basis set.
   type(tb_basisset), intent(in) :: basis
   !> Global parameters.
   type(scc_parameter), intent(in) :: param
   !> Number of neighbours for each atom.
   integer, intent(in) :: neighs(:)
   !> Integral cutoff according to prefactor from Gaussian product theorem.
   real(wp), intent(in) :: intcut
   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)
   !> Potential on each shell from the charge field (in eV).
   real(wp), intent(in) :: ves(:)
   !> Potential on each atom from the multipole field (in Eh).
   real(wp), intent(in) :: vs(:)
   !> Potential on each atom from the multipole field (in Eh).
   real(wp), intent(in) :: vd(:,:)
   !> Potential on each atom from the multipole field (in Eh).
   real(wp), intent(in) :: vq(:,:)
   !> Density matrix.
   real(wp), intent(in) :: P(:,:)
   !> Energy weighted density matrix.
   real(wp), intent(in) :: Pew(:,:)
   !> Derivative of the Hamiltonian w.r.t. the coordination number.
   real(wp), intent(out) :: dhdcn(:)
   !> Molecular gradient.
   real(wp), intent(inout) :: gradient(:,:)
   !> Stress tensor (not scaled by volume).
   real(wp), intent(inout) :: sigma(:,:)

   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)
   integer, parameter :: itt(0:3)=(/0,1,4,10/)
   integer :: iat, jat, ati, atj, ish, jsh, icao, jcao, iao, jao, ixyz, k
   integer :: ii, jj, ij, jshmax, img, iptyp, jptyp, naoi, naoj, il, jl, io, jo
   real(wp) :: r2, rij(3), ri(3), rj(3), rfc, dfc(3)
   real(wp) :: tmp(6,6), sdq(10,6,6), sdqg(3,19,6,6), hii, hjj, zi, zj, km, hav
   real(wp) :: vspd(19), vij, Pij, HPij, dG(3), dS(3,3)

   dhdcn = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:gradient, sigma, dhdcn) &
   !$omp shared(mol, neighs, neighlist, basis, param, P, Pew, ves, vs, vd, vq, &
   !$omp&       ao_n, ao_l, intcut, kcnat, cn) &
   !$omp private(jat, ij, img, ati, atj, il, jl, icao, jcao, naoi, naoj, &
   !$omp&        jshmax, iptyp, jptyp, iao, jao, ii, jj, io, jo, k, ixyz, &
   !$omp&        ri, rj, r2, rij, sdq, sdqg, tmp, zi, zj, hii, hjj, hav, &
   !$omp&        km, rfc, dfc, vspd, vij, Pij, HPij, dG, dS)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      ri = mol%xyz(:, iat)
      io = basis%shells(1, iat)-1
      do ij = 1, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dists2(ij, iat)
         rj = neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         rij = ri - rj
         jo = basis%shells(1, jat)-1

         ! generate atomic potentials from multipole field
         vspd = [vs(iat)+vs(jat), &
            &    vd(1, jat), vd(2, jat), vd(3, jat), &
            &    vq(1, jat), vq(2, jat), vq(3, jat), &
            &    vq(4, jat), vq(5, jat), vq(6, jat), &
            &    vd(1, iat), vd(2, iat), vd(3, iat), &
            &    vq(1, iat), vq(2, iat), vq(3, iat), &
            &    vq(4, iat), vq(5, iat), vq(6, iat)]

         do ish = 1, ao_n(ati)
            il = basis%lsh(ish+io)+1
            icao = basis%caoshell(ish, iat)
            naoi = llao(il-1)
            iptyp = itt(il-1)
            jshmax = ao_n(atj)
            if(iat.eq.jat) jshmax = ish
            do jsh = 1, jshmax ! jshells
               jl = basis%lsh(jsh+jo)+1
               jcao = basis%caoshell(jsh, jat)
               naoj = llao(jl-1)
               jptyp = itt(jl-1)

               hii = basis%level(ish+io) - kcnat(il-1, ati)*cn(iat)
               hjj = basis%level(jsh+jo) - kcnat(jl-1, atj)*cn(jat)

               zi = basis%zeta(ish+io)
               zj = basis%zeta(jsh+jo)
               km = h0scal(il, jl, ati, atj, basis%valsh(ish+io)/=0, &
                  &        basis%valsh(jsh+jo)/=0, param%kspd, param%kmagic, &
                  &        param%kenscal) &
                  & * (0.5_wp*(zi+zj)/sqrt(zi*zj))**(-0.5_wp)

               call drfactor(il, jl, 0, ati, atj, r2, ri, rj, rfc, dfc)
               hav = 0.5_wp * (hii + hjj) * km * rfc * evtoau

               vij = -0.5_wp*(ves(ish+io) + ves(jsh+jo)) * evtoau

               sdqg = 0.0_wp
               sdq = 0.0_wp
               call get_grad_multiint(icao, jcao, naoi, naoj, iptyp, jptyp, &
                  &                   ri, rj, intcut, basis%nprim, &
                  &                   basis%primcount, basis%alp, basis%cont, &
                  &                   sdq, sdqg)
               tmp(1:6, 1:6) = sdq(1, 1:6, 1:6)
               call dtrf2(tmp, il-1, jl-1)
               sdq(1, 1:6, 1:6) = tmp(1:6, 1:6)
               do k = 1, 19 ! 1 S, 2-4 D, 5-10 Q, 11-13 D, 14-19 Q
                  do ixyz = 1, 3
                     tmp(1:6, 1:6) = sdqg(ixyz, k, 1:6, 1:6)
                     call dtrf2(tmp, il-1, jl-1)
                     sdqg(ixyz, k, 1:6, 1:6) = tmp(1:6, 1:6)
                  enddo
               enddo
               do ii = 1, llao2(il-1)
                  iao = ii + basis%saoshell(ish, iat)
                  do jj = 1, llao2(jl-1)
                     jao = jj + basis%saoshell(jsh, jat)
                     Pij = P(jao, iao)
                     HPij = (hav + vij)*Pij - Pew(jao, iao)
                     dG = sdqg(:, 1, jj, ii)*2*Pew(jao, iao) & ! FIXME
                       & + Pij*sum(spread(vspd, 1, 3) * sdqg(:, :, jj, ii), dim=2)
                     ! add polynomial gradient
                     dG = dG + 2*hav*dfc/rfc * Pij * sdq(1, jj, ii)
                     ! generate strain contribution
                     dS = spread(dG, 1, 3) * spread(rij, 2, 3)
                     ! reset for CN derivative
                     HPij = km*rfc*Pij*sdq(1,jj,ii)*evtoau
                     if (iat /= jat) then
                        sigma = sigma + dS
                        gradient(:, iat) = gradient(:, iat) + dG
                        gradient(:, jat) = gradient(:, jat) - dG
                        dhdcn(iat) = dhdcn(iat) - HPij*kcnat(il-1, ati)
                        dhdcn(jat) = dhdcn(jat) - HPij*kcnat(jl-1, atj)
                     else
                        sigma = sigma + 0.5_wp * dS
                     endif
                  enddo
               enddo
            enddo ! jsh : loop over shells on jat
         enddo ! ish : loop over shells on iat
      enddo ! jat
   enddo ! iat
   !$omp end parallel do

   do iao = 1, basis%nao
      iat = basis%aoat2(iao)
      ati = mol%at(iat)
      ish = basis%ao2sh(iao)
      il = basis%lsh(ish)+1
      ! diagonal contribution
      dhdcn(iat) = dhdcn(iat)-P(iao,iao)*kcnat(il-1,ati)*evtoau
   enddo

end subroutine build_dSDQH0

end submodule aes_integrals
