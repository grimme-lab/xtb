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

!> Implementation of the xTB core Hamiltonian
module xtb_xtb_hamiltonian
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_xtb_data, only : THamiltonianData
   use xtb_intgrad
   implicit none
   private

   public :: getSelfEnergy, build_SDQH0
   public :: count_dpint, count_qpint


   interface getSelfEnergy
      module procedure :: getSelfEnergyFlat
      module procedure :: getSelfEnergy2D
   end interface getSelfEnergy


   integer,private, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer,private, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

contains


subroutine getSelfEnergyFlat(hData, nShell, at, cn, qat, selfEnergy, dSEdcn, dSEdq)
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in) :: nShell(:)
   integer, intent(in) :: at(:)
   real(wp), intent(in), optional :: cn(:)
   real(wp), intent(in), optional :: qat(:)
   real(wp), intent(out) :: selfEnergy(:)
   real(wp), intent(out), optional :: dSEdcn(:)
   real(wp), intent(out), optional :: dSEdq(:)

   integer :: ind, iAt, iZp, iSh

   selfEnergy(:) = 0.0_wp
   if (present(dSEdcn)) dSEdcn(:) = 0.0_wp
   if (present(dSEdq)) dSEdq(:) = 0.0_wp
   ind = 0
   do iAt = 1, size(cn)
      iZp = at(iAt)
      do iSh = 1, nShell(iZp)
         selfEnergy(ind+iSh) = hData%selfEnergy(iSh, iZp)
      end do
      ind = ind + nShell(iZp)
   end do
   if (present(dSEdcn) .and. present(cn)) then
      ind = 0
      do iAt = 1, size(cn)
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(ind+iSh) = selfEnergy(ind+iSh) &
               & - hData%kCN(iSh, iZp) * cn(iAt)
            dSEdcn(ind+iSh) = -hData%kCN(iSh, iZp)
         end do
         ind = ind + nShell(iZp)
      end do
   end if
   if (present(dSEdq) .and. present(qat)) then
      ind = 0
      do iAt = 1, size(cn)
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(ind+iSh) = selfEnergy(ind+iSh) &
               & - hData%kQShell(iSh,iZp)*qat(iAt) - hData%kQAtom(iZp)*qat(iAt)**2
            dSEdq(ind+iSh) = -hData%kQShell(iSh,iZp) - hData%kQAtom(iZp)*2*qat(iAt)
         end do
         ind = ind + nShell(iZp)
      end do
   end if

end subroutine getSelfEnergyFlat


subroutine getSelfEnergy2D(hData, nShell, at, cn, qat, selfEnergy, dSEdcn, dSEdq)
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in) :: nShell(:)
   integer, intent(in) :: at(:)
   real(wp), intent(in), optional :: cn(:)
   real(wp), intent(in), optional :: qat(:)
   real(wp), intent(out) :: selfEnergy(:, :)
   real(wp), intent(out), optional :: dSEdcn(:, :)
   real(wp), intent(out), optional :: dSEdq(:, :)

   integer :: iAt, iZp, iSh

   selfEnergy(:, :) = 0.0_wp
   if (present(dSEdcn)) dSEdcn(:, :) = 0.0_wp
   if (present(dSEdq)) dSEdq(:, :) = 0.0_wp
   do iAt = 1, size(cn)
      iZp = at(iAt)
      do iSh = 1, nShell(iZp)
         selfEnergy(iSh, iAt) = hData%selfEnergy(iSh, iZp)
      end do
   end do
   if (present(dSEdcn) .and. present(cn)) then
      do iAt = 1, size(cn)
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(iSh, iAt) = selfEnergy(iSh, iAt) &
               & - hData%kCN(iSh, iZp) * cn(iAt)
            dSEdcn(iSh, iAt) = -hData%kCN(iSh, iZp)
         end do
      end do
   end if
   if (present(dSEdq) .and. present(qat)) then
      do iAt = 1, size(cn)
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(iSh, iAt) = selfEnergy(iSh, iAt) &
               & - hData%kQShell(iSh,iZp)*qat(iAt) - hData%kQAtom(iZp)*qat(iAt)**2
            dSEdq(iSh, iAt) = -hData%kQShell(iSh,iZp) &
               & - hData%kQAtom(iZp)*2*qat(iAt)
         end do
      end do
   end if

end subroutine getSelfEnergy2D

!> Computes the dipole and quadrupole integrals and performs screening to
!  determine, which contribute to potential
subroutine build_SDQH0(nShell, hData, nat, at, nbf, nao, xyz, selfEnergy, &
      & intcut, caoshell, saoshell, nprim, primcount, alp, cont, &
      & sint, dpint, qpint, H0)
   implicit none
   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
   !> # of atoms
   integer, intent(in)  :: nat
   !> # of spherical AOs (SAOs)
   integer, intent(in)  :: nao
   !> # of Cartesian AOs (CAOs)
   integer, intent(in)  :: nbf
   integer, intent(in)  :: at(nat)
   !> Cartesian coordinates
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp), intent(in) :: selfEnergy(:)
   !> Integral cutoff according to prefactor from Gaussian product theorem
   real(wp),intent(in)  :: intcut
   !> Map shell of atom to index in CAO space (lowest Cart. component is taken)
   integer, intent(in)  :: caoshell(:,:)
   !> Map shell of atom to index in SAO space (lowest m_l component is taken)
   integer, intent(in)  :: saoshell(:,:)
   integer, intent(in)  :: nprim(:)
   !> Index of first primitive (over entire system) of given CAO
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)
   !> Overlap integral matrix
   real(wp),intent(out) :: sint(nao,nao)
   !> Dipole integral matrix
   real(wp),intent(out) :: dpint(3,nao,nao)
   !> Quadrupole integral matrix
   real(wp),intent(out) :: qpint(6,nao,nao)
   !> Core Hamiltonian
   real(wp),intent(out) :: H0(:)


   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,km,mi,mj,ij
   real(wp) tmp1,tmp2,tmp3,tmp4,step,step2
   real(wp) dx,dy,dz,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cc,cj,alpi,rab2,ab,est

   real(wp)  ra(3),rb(3),f1,f2,point(3)
   real(wp) dtmp(3),qtmp(6),ss(6,6),dd(3,6,6),qq(6,6,6),tmp(6,6)
   integer ip,jp,iat,jat,izp,jzp,ish,jsh,icao,jcao,iao,jao,jshmax
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer itt(0:3)
   parameter(itt  =(/0,1,4,10/))
   real(wp) :: saw(10)


   ! integrals
   sint = 0.0_wp
   dpint = 0.0_wp
   qpint = 0.0_wp
   ! --- Aufpunkt for moment operator
   point = 0.0_wp

   !$OMP PARALLEL &
   !$omp PRIVATE (iat,jat,izp,cc,ci,ra,rb,saw, &
   !$omp& rab2,jzp,ish,ishtyp,icao,naoi,iptyp, &
   !$omp& jsh,jshmax,jshtyp,jcao,naoj,jptyp,ss,dd,qq, &
   !$omp& est,alpi,alpj,ab,iprim,jprim,ip,jp, &
   !$omp& mli,mlj,tmp,tmp1,tmp2,iao,jao,ii,jj,k,ij) &
   !$omp reduction(+:sint,dpint,qpint)
   !$OMP DO schedule(dynamic)
   do iat = 1,nat
      ra(1:3) = xyz(1:3,iat)
      izp = at(iat)
      do jat = 1,iat-1
         rb(1:3) = xyz(1:3,jat)
         jzp = at(jat)
         rab2 = sum( (rb-ra)**2 )
         !           ints < 1.d-9 for RAB > 40 Bohr
         if(rab2.gt.2000) cycle
         do ish = 1,nShell(izp)
            ishtyp = hData%angShell(ish,izp)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            do jsh = 1,nShell(jzp)
               jshtyp = hData%angShell(jsh,jzp)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

!               ! get indices
!               i = 1+basis%saoshell(ish, iat)
!               j = 1+basis%saoshell(jsh, jat)
!
!               call h0scal(hData, ishtyp, jshtyp, izp, jzp, valao2(i).ne.0, &
!                  & valao2(j).ne.0, km)
!
!               ! we scale the two shells depending on their exponent
!               zi = basis%aoexp(i)
!               zj = basis%aoexp(j)
!               zetaij = (2 * sqrt(zi*zj)/(zi+zj))**hData%wExp

               ss = 0.0_wp
               dd = 0.0_wp
               qq = 0.0_wp
               call get_multiints(icao,jcao,naoi,naoj,iptyp,jptyp,ra,rb,point, &
                  &               intcut,nprim,primcount,alp,cont,ss,dd,qq)
               !transform from CAO to SAO
               call dtrf2(ss,ishtyp,jshtyp)
               do k = 1,3
                  tmp(1:6,1:6) = dd(k,1:6,1:6)
                  call dtrf2(tmp,ishtyp,jshtyp)
                  dd(k,1:6,1:6) = tmp(1:6,1:6)
               enddo
               do k = 1,6
                  tmp(1:6,1:6) = qq(k,1:6,1:6)
                  call dtrf2(tmp,ishtyp,jshtyp)
                  qq(k,1:6,1:6) = tmp(1:6,1:6)
               enddo
               do ii = 1,llao2(ishtyp)
                  iao = ii+saoshell(ish,iat)
                  do jj = 1,llao2(jshtyp)
                     jao = jj+saoshell(jsh,jat)
                     sint(iao, jao) = sint(iao, jao) + ss(jj, ii)
                     sint(jao, iao) = sint(jao, iao) + ss(jj, ii)
                     dpint(1:3, iao, jao) = dpint(1:3, iao, jao) + dd(1:3, jj, ii)
                     dpint(1:3, jao, iao) = dpint(1:3, jao, iao) + dd(1:3, jj, ii)
                     qpint(1:6, iao, jao) = qpint(1:6, iao, jao) + qq(1:6, jj, ii)
                     qpint(1:6, jao, iao) = qpint(1:6, jao, iao) + qq(1:6, jj, ii)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$OMP END DO
   !$OMP END PARALLEL

   ! diagonal elements
   do iat = 1, nat
      ra = xyz(:, iat)
      izp = at(iat)
      do ish = 1, nShell(izp)
         ishtyp = hData%angShell(ish,izp)
         do iao = 1, llao2(ishtyp)
            i = iao+saoshell(ish,iat)
            sint(i,i) = 1.0_wp + sint(i,i)
         end do

         icao = caoshell(ish,iat)
         naoi = llao(ishtyp)
         iptyp = itt(ishtyp)
         do jsh = 1, ish
            jshtyp = hData%angShell(jsh,izp)
            jcao = caoshell(jsh,iat)
            naoj = llao(jshtyp)
            jptyp = itt(jshtyp)
            ss = 0.0_wp
            dd = 0.0_wp
            qq = 0.0_wp
            call get_multiints(icao,jcao,naoi,naoj,iptyp,jptyp,ra,ra,point, &
               &               intcut,nprim,primcount,alp,cont,ss,dd,qq)
            !transform from CAO to SAO
            !call dtrf2(ss,ishtyp,jshtyp)
            do k = 1,3
               tmp(1:6, 1:6) = dd(k,1:6, 1:6)
               call dtrf2(tmp, ishtyp, jshtyp)
               dd(k, 1:6, 1:6) = tmp(1:6, 1:6)
            enddo
            do k = 1,6
               tmp(1:6, 1:6) = qq(k, 1:6, 1:6)
               call dtrf2(tmp, ishtyp, jshtyp)
               qq(k, 1:6, 1:6) = tmp(1:6, 1:6)
            enddo
            do ii = 1, llao2(ishtyp)
               iao = ii + saoshell(ish,iat)
               do jj = 1, llao2(jshtyp)
                  jao = jj + saoshell(jsh,iat)
                  if (jao > iao .and. ish ==  jsh) cycle
                  dpint(1:3, iao, jao) = dpint(1:3, iao, jao) + dd(1:3, jj, ii)
                  if (iao /= jao) then
                     dpint(1:3, jao, iao) = dpint(1:3, jao, iao) + dd(1:3, jj, ii)
                  end if
                  qpint(1:6, iao, jao) = qq(1:6, jj, ii)
                  if (jao /= iao) then
                     qpint(1:6, jao, iao) = qq(1:6, jj, ii)
                  end if
               end do
            end do
         end do
      end do
   end do

end subroutine build_SDQH0


!> Count number of significant dipole integrals
subroutine count_dpint(ndp, dpint, thr)

   !> Number of significant dipole integrals
   integer, intent(out) :: ndp

   !> Dipole integrals
   real(wp), intent(in) :: dpint(:, :, :)

   !> Neglect threshold for dipole integrals
   real(wp), intent(in) :: thr

   integer :: i, j
   real(wp) :: tmp1, thr2

   ndp = 0
   thr2 = (thr*1.0e-2_wp)-thr*1.0e-12_wp

   do i = 1, size(dpint, dim=3)
      do j = 1, i
         tmp1 = sum(dpint(1:3, j, i)**2)
         if (tmp1 > thr2) ndp = ndp + 1
      enddo
   enddo

end subroutine count_dpint


!> Count number of significant quadrupole integrals
subroutine count_qpint(nqp, qpint, thr)

   !> Number of significant quadrupole integrals
   integer, intent(out) :: nqp

   !> Quadrupole integrals
   real(wp), intent(in) :: qpint(:, :, :)

   !> Neglect threshold for quadrupole integrals
   real(wp), intent(in) :: thr

   integer :: i, j
   real(wp) :: tmp2, thr2

   nqp = 0
   thr2 = (thr*1.0e-2_wp)-thr*1.0e-12_wp

   do i = 1, size(qpint, dim=3)
      do j = 1, i
         tmp2 = sum(qpint(1:3, j, i)**2) + 2*sum(qpint(4:6, j, i)**2)
         if (tmp2 > thr2) nqp = nqp + 1
      enddo
   enddo

end subroutine count_qpint


end module xtb_xtb_hamiltonian
