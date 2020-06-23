! This file is part of xtb.
!
! Copyright (C) 2020 Sebastian Ehlert
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

!> Implementation of the Born radii integrator
module xtb_solv_born
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas, only : mctc_scal
   implicit none
   private

   public :: compute_bornr

   !> van der Waals to Lee-Richard's surface correction (GBOBCII parameter)
   real(wp), parameter :: alp = 1._wp
   real(wp), parameter :: bet = 0.8_wp
   real(wp), parameter :: gam = 4.85_wp

contains


subroutine compute_bornr(nat, nnrad, nnlistr, ddpair, vdwr, rho, svdw, c1, &
      & brad, brdr)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Number of neighbours to consider for radii
   integer, intent(in) :: nnrad

   !> Neighbourlist
   integer, intent(in) :: nnlistr(:, :)

   !> Position and distances
   real(wp), intent(in) :: ddpair(:, :)

   !> Van-der-Waals radii
   real(wp), intent(in) :: vdwr(:)

   !> Descreened van-der-Waals radii
   real(wp), intent(in) :: rho(:)

   !> van-der-Waals radii with offset
   real(wp), intent(in) :: svdw(:)

   !> Scaling factor for the Born radii
   real(wp), intent(in) :: c1

   !> Born radii
   real(wp), intent(out) :: brad(:)

   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), intent(out) :: brdr(:, :, :)

   integer :: iat
   real(wp) :: br, dpsi, svdwi, vdwri, s1, v1, s2, arg, arg2
   real(wp) :: th, ch, alpi, beti, gami

   call compute_psi(nat, nnrad, nnlistr, ddpair, vdwr, rho, brad, brdr)

   do iat = 1, nat

      br = brad(iat)

      svdwi = svdw(iat)
      vdwri = vdwr(iat)
      s1 = 1.0_wp/svdwi
      v1 = 1.0_wp/vdwri
      s2 = 0.5_wp*svdwi

      br = br*s2

      arg2 = br*(gam*br-bet)
      arg = br*(alp+arg2)
      arg2 = 2.0_wp*arg2+alp+gam*br*br

      th = tanh(arg)
      ch = cosh(arg)

      br = 1.0_wp/(s1-v1*th)
      ! Include GBMV2-like scaling
      br = c1*br

      dpsi = ch*(s1-v1*th)
      dpsi = s2*v1*arg2/(dpsi*dpsi)
      dpsi = c1*dpsi

      brad(iat) = br
      call mctc_scal(brdr(:, :, iat), dpsi)

   end do

end subroutine compute_bornr


pure subroutine compute_psi(nat, nnrad, nnlistr, ddpair, vdwr, rho, psi, dpsidr)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Number of neighbours to consider for radii
   integer, intent(in) :: nnrad

   !> Neighbourlist
   integer, intent(in) :: nnlistr(:, :)

   !> Position and distances
   real(wp), intent(in) :: ddpair(:, :)

   !> Van-der-Waals radii
   real(wp), intent(in) :: vdwr(:)

   !> Descreened van-der-Waals radii
   real(wp), intent(in) :: rho(:)

   !> Integrated value of Psi
   real(wp), intent(out) :: psi(:)

   !> Derivative of Psi w.r.t. cartesian coordinates
   real(wp), intent(out) :: dpsidr(:, :, :)

   real(wp), allocatable :: dpsitr(:, :)
   integer  :: kk
   integer  :: i, ii, jj, nn
   real(wp) :: dr(3), r, rhoi, rhoj
   real(wp) :: gi, gj, ap, am, lnab, rhab, ab, dgi, dgj
   real(wp) :: drjj(3)
   real(wp) :: rh1, rhr1, r24, rh2, r1, aprh1, r12
   real(wp) :: rvdwi, rvdwj
   integer  :: ovij, ovji, ov

   allocate(dpsitr(3, nat))
   psi(:) = 0.0_wp
   dpsidr(:, :, :) = 0.0_wp
   dpsitr(:, :) = 0.0_wp

   do kk = 1, nnrad

      ii = nnlistr(1, kk)
      jj = nnlistr(2, kk)
      nn = nnlistr(3, kk)

      r = ddpair(1, nn)
      dr(:) = ddpair(2:4, nn)

      rhoi = rho(ii)
      rhoj = rho(jj)
      rvdwi = vdwr(ii)
      rvdwj = vdwr(jj)

      ovij = 1
      ovji = 1
      if(r.ge.(rvdwi+rhoj)) ovij = 0
      if(r.ge.(rhoi+rvdwj)) ovji = 0
      ov = ovij+10*ovji

      select case(ov)
      case(0) ! ij do not overlap; ji do not overlap
         ! nonoverlaping spheres
         if(abs(rhoi-rhoj).lt.1.d-8) then
            ! equal reduced radii
            r1 = 1.0_wp/r
            ap = r+rhoj
            am = r-rhoj
            ab = ap*am
            rhab = rhoj/ab
            lnab = 0.5_wp*log(am/ap)*r1
            gi = rhab+lnab
            dgi = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
            ! accumulate psi
            psi(ii) = psi(ii)+gi
            psi(jj) = psi(jj)+gi
            ! accumulate psi gradient
            drjj(:) = dgi*dr(:)
            dpsitr(:, ii) = dpsitr(:, ii)+drjj(:)
            dpsidr(:, jj, ii) = dpsidr(:, jj, ii)-drjj(:)
            dpsitr(:, jj) = dpsitr(:, jj)-drjj(:)
            dpsidr(:, ii, jj) = dpsidr(:, ii, jj)+drjj(:)
         else
            ! unequal reduced radii
            ! ij contribution
            r1 = 1.0_wp/r
            ap = r+rhoj
            am = r-rhoj
            ab = ap*am
            rhab = rhoj/ab
            lnab = 0.5_wp*log(am/ap)*r1
            gi = rhab+lnab
            dgi = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
            ! ji contribution
            ap = r+rhoi
            am = r-rhoi
            ab = ap*am
            rhab = rhoi/ab
            lnab = 0.5_wp*log(am/ap)*r1
            gj = rhab+lnab
            dgj = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
            ! accumulate psi
            psi(ii) = psi(ii)+gi
            psi(jj) = psi(jj)+gj
            ! accumulate psi gradient
            drjj(:) = dgi*dr(:)
            dpsitr(:, ii) = dpsitr(:, ii)+drjj(:)
            dpsidr(:, jj, ii) = dpsidr(:, jj, ii)-drjj(:)

            drjj(:) = dgj*dr(:)
            dpsitr(:, jj) = dpsitr(:, jj)-drjj(:)
            dpsidr(:, ii, jj) = dpsidr(:, ii, jj)+drjj(:)
         endif

      case(10) ! ij do not overlap; ji overlap

         ! ij contribution
         r1 = 1.0_wp/r
         ap = r+rhoj
         am = r-rhoj
         ab = ap*am
         rhab = rhoj/ab
         lnab = 0.5_wp*log(am/ap)*r1
         gi = rhab+lnab
         dgi = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
         ! accumulate psi
         psi(ii) = psi(ii)+gi
         ! accumulate psi gradient
         drjj(:) = dgi*dr(:)
         dpsitr(:, ii) = dpsitr(:, ii)+drjj(:)
         dpsidr(:, jj, ii) = dpsidr(:, jj, ii)-drjj(:)

         if((r+rhoi).gt.rvdwj) then
            ! ji contribution
            r1 = 1.0_wp/r
            r12 = 0.5_wp*r1
            r24 = r12*r12

            ap = r+rhoi
            am = r-rhoi
            rh1 = 1.0_wp/rvdwj
            rhr1 = 1.0_wp/ap
            aprh1 = ap*rh1
            lnab = log(aprh1)

            gj = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

            dgj = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
               &         rhoi*r24*(rhr1-rh1*aprh1)+ &
               &         r12*(r1*lnab-rhr1)
            dgj = dgj*r1
            ! accumulate psi
            psi(jj) = psi(jj)+gj
            ! accumulate psi gradient
            drjj(:) = dgj*dr(:)
            dpsitr(:, jj) = dpsitr(:, jj)-drjj(:)
            dpsidr(:, ii, jj) = dpsidr(:, ii, jj)+drjj(:)
         endif

      case(1) ! ij overlap; ji do not overlap

         if((r+rhoj).gt.rvdwi) then
            ! ij contribution
            r1 = 1.0_wp/r
            r12 = 0.5_wp*r1
            r24 = r12*r12

            ap = r+rhoj
            am = r-rhoj
            rh1 = 1.0_wp/rvdwi
            rhr1 = 1.0_wp/ap
            aprh1 = ap*rh1
            lnab = log(aprh1)

            gi = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

            dgi = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
               &         rhoj*r24*(rhr1-rh1*aprh1)+ &
               &         r12*(r1*lnab-rhr1)
            dgi = dgi*r1
            ! accumulate psi
            psi(ii) = psi(ii)+gi
            ! accumulate psi gradient
            drjj(:) = dgi*dr(:)
            dpsitr(:, ii) = dpsitr(:, ii)+drjj(:)
            dpsidr(:, jj, ii) = dpsidr(:, jj, ii)-drjj(:)
         endif

         ! ji contribution
         ap = r+rhoi
         am = r-rhoi
         ab = ap*am
         rhab = rhoi/ab
         lnab = 0.5_wp*log(am/ap)*r1
         gj = rhab+lnab
         dgj = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
         ! accumulate psi
         psi(jj) = psi(jj)+gj
         ! accumulate psi gradient
         drjj(:) = dgj*dr(:)
         dpsitr(:, jj) = dpsitr(:, jj)-drjj(:)
         dpsidr(:, ii, jj) = dpsidr(:, ii, jj)+drjj(:)

      case(11) ! ij and ji overlap
         ! overlaping spheres
         if((r+rhoj).gt.rvdwi) then
            ! ij contribution
            r1 = 1.0_wp/r
            r12 = 0.5_wp*r1
            r24 = r12*r12

            ap = r+rhoj
            am = r-rhoj
            rh1 = 1.0_wp/rvdwi
            rhr1 = 1.0_wp/ap
            aprh1 = ap*rh1
            lnab = log(aprh1)

            gi = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

            dgi = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
               &         rhoj*r24*(rhr1-rh1*aprh1)+ &
               &         r12*(r1*lnab-rhr1)
            dgi = dgi*r1
            ! accumulate psi
            psi(ii) = psi(ii)+gi
            ! accumulate psi gradient
            drjj(:) = dgi*dr(:)
            dpsitr(:, ii) = dpsitr(:, ii)+drjj(:)
            dpsidr(:, jj, ii) = dpsidr(:, jj, ii)-drjj(:)
         endif

         if((r+rhoi).gt.rvdwj) then
            ! ji contribution
            r1 = 1.0_wp/r
            r12 = 0.5_wp*r1
            r24 = r12*r12

            ap = r+rhoi
            am = r-rhoi
            rh1 = 1.0_wp/rvdwj
            rhr1 = 1.0_wp/ap
            aprh1 = ap*rh1
            lnab = log(aprh1)

            gj = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

            dgj = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
               &         rhoi*r24*(rhr1-rh1*aprh1)+ &
               &         r12*(r1*lnab-rhr1)
            dgj = dgj*r1
            ! accumulate psi
            psi(jj) = psi(jj)+gj
            ! accumulate psi gradient
            drjj(:) = dgj*dr(:)
            dpsitr(:, jj) = dpsitr(:, jj)-drjj(:)
            dpsidr(:, ii, jj) = dpsidr(:, ii, jj)+drjj(:)
         endif

      end select

   enddo

   ! save one-center terms
   do i = 1, nat
      dpsidr(:, i, i) = dpsitr(:, i)
   enddo

end subroutine compute_psi

end module xtb_solv_born
