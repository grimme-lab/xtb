! This file is part of xtb.
!
! Copyright (C) 2022 Christoph Plett
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

module xtb_iff_iffini
   use xtb_mctc_accuracy, only: wp
   use xtb_docking_param
   use xtb_iff_ifflmo
   use xtb_sphereparam
   use xtb_disp_ncoord
   use xtb_type_environment, only: TEnvironment
   use xtb_setmod
   use xtb_type_param
   use xtb_disp_dftd4param
   use xtb_disp_dftd4, only : newD4Model, zeta, weight_references, get_atomic_c6, d4, d4dim
   use xtb_type_dispersionmodel, only: TDispersionModel
   use xtb_mctc_param, only: gamma => chemical_hardness
   use xtb_type_molecule, only: TMolecule, init
   implicit none

   private
   public :: init_iff

contains
   subroutine init_iff(env,n1,n2,at1,at2,neigh,xyz1,xyz2,q1,q2,c6ab,z1,z2,&
   &                   cprob, nlmo1, nlmo2, lmo1, lmo2,&
   &                   qdr1, qdr2, rlmo1, rlmo2,&
   &                   cn1, cn2, alp1, alp2, alpab,&
   &                   den1, den2, gab1, gab2, qcm1, qcm2,&
   &                   n, at, xyz, q, icoord, icoord0,&
   &                   pr)

      type(TEnvironment), intent(inout) :: env
      integer, intent(in) :: n, n1, n2
      integer, intent(inout) :: nlmo1, nlmo2
      integer, intent(out) :: neigh(0:n, n)
      integer, intent(out) :: at(n)
      integer, intent(in) :: at1(n1)
      integer, intent(in) :: at2(n2)
      integer, intent(inout) :: lmo1(n1*10)
      integer, intent(inout) :: lmo2(n2*10)
      real(wp), intent(out) :: xyz(3, n)
      real(wp), intent(out) :: q(n)
      real(wp), intent(inout) :: xyz1(3, n1)
      real(wp), intent(inout) :: xyz2(3, n2)
      real(wp), intent(inout) :: rlmo1(4, n1*10)
      real(wp), intent(inout) :: rlmo2(4, n2*10)
      real(wp), intent(in) :: q1(n1)
      real(wp), intent(out) :: qdr1(n1)
      real(wp), intent(in) :: q2(n2)
      real(wp), intent(out) :: qdr2(n2)
      real(wp), intent(out) :: cn1(n1)
      real(wp), intent(out) :: cn2(n2)
      real(wp), intent(out) :: alp1(n1)
      real(wp), intent(out) :: alp2(n2)
      real(wp), intent(out) :: c6ab(n, n)
      real(wp), intent(out) :: z1(n1)
      real(wp), intent(out) :: z2(n2)
      real(wp), intent(in) :: qcm1(n1)
      real(wp), intent(in) :: qcm2(n2)
      real(wp), intent(out) :: alpab(n2, n1)
      real(wp), intent(out) :: cprob(n1)
      real(wp), intent(inout) :: icoord(6), icoord0(6)
      real(wp), intent(out) :: den1(2,4,n1),den2(2,4,n2),gab1(n1,n1),gab2(n2,n2)
      logical, intent(in) :: pr

      real(wp) :: dum_array(3, 1)
      integer, allocatable :: dum_at(:)
      real(wp) :: distance, max_distanceA, max_distanceB
      integer :: i, j, k, nfreq, np
      real(wp) :: edum, r(3), rr, rrr, r42, c6d4, tmp, thr, a1, a2, s8, ai, aj, p1, fc
      real(wp) :: aa, bb, cc, avmom, wt, h298, g298, ts298, zp, symn 
      real(wp) :: freq1(3*n1), freq2(3*n2), dum(n), dum1(n1), dum2(n2)
      real(wp) :: cn(n), alp0(n), c6ab1(n1, n1), c6ab2(n2, n2), gsolv_inf
      real(wp) :: f,fall,fhyd,c6(94),qt1(1),qt2(1),r0,e1(4),c1(4),c6aad4,boxoffset
      !> D3 CAA(CNmax) coefficients
      data c6/& 
      &  3.03, 1.56, 85.32, 55.14, 28.03, 18.21, 15.58, 10.37,&
      &  7.13, 6.29, 186.11, 175.56, 153.59, 149.77, 151.69, 125.81,&
      & 90.40, 64.65, 338.02, 343.33, 313.98, 321.59, 269.33, 153.34,&
      &113.85, 109.50, 118.29, 141.88, 175.00, 153.45, 172.32, 180.29,&
      &224.29, 210.52, 169.24, 130.40, 476.26, 491.74, 451.08, 458.69,&
      &428.90, 337.11, 265.87, 239.01, 244.66, 265.87, 268.61, 263.45,&
      &310.83, 336.10, 417.38, 418.58, 358.50, 290.22, 649.09, 716.55,&
      &766.02, 688.04, 647.02, 620.52, 581.74, 587.79, 549.03, 595.36,&
      &524.92, 544.44, 472.80, 446.23, 440.52, 487.08, 411.57, 420.98,&
      &458.77, 398.10, 332.14, 297.83, 305.53, 336.87, 317.16, 298.86,&
      &389.61, 455.29, 533.88, 558.97, 495.45, 412.83, 822.05, 867.56,&
      &906.04, 1224.59, 980.38, 753.86, 698.16, 1051.25/
      real(wp), parameter ::au = 627.509541d0
      character*4 pgroup
      logical :: linear

      alp0 = 0.0_wp
      c6ab = 0.0_wp
      qt1 = 0.0_wp
      qt2 = 0.0_wp

! copy the whole system for CM5, GBSA, D4
      do i = 1, n1
         xyz(1:3, i) = xyz1(1:3, i)
         at(i) = at1(i)
         q(i) = q1(i)
      end do
      do i = 1, n2
         xyz(1:3, i + n1) = xyz2(1:3, i)
         at(i + n1) = at2(i)
         q(i + n1) = q2(i)
      end do
! q1 & q2 eingelesen
!     used for xtb recalc
      chrg(1) = sum(q1) + 1.d-8
      chrg(2) = sum(q2) + 1.d-8

!     ccccccccccccccccccccccccccccccccccccccccccccccccc
!               preCOMPUTE D3,D4 and REP stuff
!     ccccccccccccccccccccccccccccccccccccccccccccccccc

      if (sphere .eq. 1) write (*, '(''spherical box radius     :'',F10.3)') boxr
      if (sphere .eq. 2) write (*, '(''ellipsoid axis lengths   :'',3F10.3)') rabc

      !> Determine a factor to shift molB away from molA so that the Disp coefficients are those of seperated mols
      max_distanceA = 0.0_wp
      distance = 0.0_wp
      do i = 1, n1
         do j = 1, i - 1
            distance = sqrt(sum((xyz1(:, j) - xyz1(:, i))**2))
            max_distanceA = max(max_distanceA, distance)
         end do
      end do

      max_distanceB = 0.0_wp
      distance = 0.0_wp
      do i = 1, n2
         do j = 1, i - 1
            distance = sqrt(sum((xyz2(:, j) - xyz2(:, i))**2))
            max_distanceB = max(max_distanceB, distance)
         end do
      end do

      shift_geo = max_distanceA + max_distanceB + 1.d+6
!      shift_geo = 1.d+8

      do i = 1, n2  ! shift (later done reversely) to get seperated disp coefficients
         xyz(1:3, i + n1) = xyz2(1:3, i) + shift_geo
      end do

      a1 = par_d3_a1
      a2 = par_d3_a2
      s8 = par_d3_s8
      do i = 1, 86 ! precompute BJ damping terms
         dum_array(1:3, 1) = 0.0_wp
         call c6_alp0(1, [i], dum_array, [0.0d0], [0.d0], qt1, qt2) !qt1=a(0), qt2=c6ab
         r0_atom(i) = qt1(1)**(1./3.)   ! free atom values alpha(0) = radius | Saved in dockingparam
         do j = 1, i
            r42 = r2r4(i)*r2r4(j)
            rrab(j, i) = r42*3.0d0*s8  !include s8
            rrab(i, j) = r42*3.0d0*s8
            r0ab6(j, i) = (a1*dsqrt(3.0d0*r42) + a2)**6
            r0ab8(j, i) = (a1*dsqrt(3.0d0*r42) + a2)**8
            r0ab6(i, j) = r0ab6(j, i)
            r0ab8(i, j) = r0ab8(j, i)
         end do
      end do

      !> All C6 and alp0 coefficients for molA and molB
      call ncoord_d4(n, at, xyz, cn, 500.d0)
      call c6_alp0(n, at, xyz, cn, q, alp0, c6ab)
      c6aad4 = 0
      do i = 1, n
      do j = 1, n
         c6aad4 = c6aad4 + c6ab(j, i)
      end do
      end do

! at this point we have
! q1/q2    : fragment atomic Mulliken charges
! q        : full system atomic CM5 charges

      c6d4 = 0
      do i = 1, n1
         cprob(i) = sqrt(c6ab(i, i)*c6(probe_atom_type)) ! probe RG atom C6AB
         ai = alp0(i)**(1./3)
         do j = 1, n2
            aj = alp0(j + n1)**(1./3)
            alpab(j, i) = 2.6*(ai + aj)*0.5               ! the factor just converts the values to the
            ! r0ab scale so that the element factors are about one
            c6d4 = c6d4 + c6ab(j + n1, i)
         end do
      end do

      if(set%verbose) write (env%unit, '(''molecular C6AA D4 /au    :'',2F12.1)') c6aad4
      if(set%verbose) write (env%unit, '(''molecular C6AB D4 /au    :'',2F12.1)') c6d4

! setup ATM and other stuff
      thr = 25.0 ! intra mol cut-off to avoid too many pairs right from the start
      do i = 1, n1
         alp1(i) = alp0(i)
         cn1(i) = cn(i)
         z1(i) = val_e(at1(i))
         k = 0
         do j = 1, n1
            rr = (xyz1(1, j) - xyz1(1, i))**2 + (xyz1(2, j) - xyz1(2, i))**2 +&
              &(xyz1(3, j) - xyz1(3, i))**2
            r0 = sqrt(rr)
            gab1(j, i) = 1.0d0/(r0 + 1./sqrt(gam(at1(i))*gam(at1(j))))
            if (j .eq. i) cycle
            ai = alp0(i)**(1./3)
            aj = alp0(j)**(1./3)
            rrr = 1.13*2.6*(ai + aj)*0.5/r0
            tmp = 1.0d0/(1.d0 + 6.d0*rrr**14)
            if (tmp .gt. 0.90 .and. r0 .lt. thr) then
               k = k + 1
               neigh(k, i) = j
            end if
         end do
         neigh(0, i) = k
      end do
!     molB
      do i = 1, n2
         alp2(i) = alp0(i + n1)
         cn2(i) = cn(i + n1)
         z2(i) = val_e(at2(i))
         k = 0
         do j = 1, n2
            rr = (xyz2(1, j) - xyz2(1, i))**2 + (xyz2(2, j) - xyz2(2, i))**2 +&
              &(xyz2(3, j) - xyz2(3, i))**2
            r0 = sqrt(rr)
            gab2(j, i) = 1.0d0/(r0 + 1./sqrt(gam(at2(i))*gam(at2(j))))
            if (j .eq. i) cycle
            ai = alp0(i + n1)**(1./3)
            aj = alp0(j + n1)**(1./3)
            rrr = 1.13*2.6*(ai + aj)*0.5/r0
            tmp = 1.0d0/(1.d0 + 6.d0*rrr**14)
            if (tmp .gt. 0.90 .and. r0 .lt. thr) then
               k = k + 1
               neigh(k, i + n1) = j
            end if
         end do
         neigh(0, i + n1) = k
      end do

!     define LP centers
      call setlmo(1, n1, nlmo1, at1, xyz1, alp1, q1, cn1, rlmo1, lmo1)
      call setlmo(2, n2, nlmo2, at2, xyz2, alp2, q2, cn2, rlmo2, lmo2)

      if (debug) then ! write actual centers to file
         call wrc('xtbiff.coord', n1, n2, at1, at2, xyz1, xyz2, icoord)
         call wrlmocoord('lmocent2.coord', n1, n2, xyz1, xyz2,&
                        &at1, at2, nlmo1, nlmo2,&
                        &lmo1, lmo2, rlmo1, rlmo2)
      end if

      icoord = 0
!     shift to CMA and xyz1 is not changed anymore
      if (nfrag1 .eq. 0) then
         call cmadock(n1, n1, at1, xyz1, r) !cma kommt als r(3) zurück
      else
         call cmadock(nfrag1, n1, at1, xyz1, r)
      end if
      icoord(1:3) = -r(1:3)
      do i = 1, 3
         xyz1(i, 1:n1) = xyz1(i, 1:n1) - r(i)
         rlmo1(i, 1:nlmo1) = rlmo1(i, 1:nlmo1) - r(i)
      end do
      call cmadock(n2, n2, at2, xyz2, r)
      icoord(1:3) = icoord(1:3) + r(1:3)
      do i = 1, 3
         xyz2(i, 1:n2) = xyz2(i, 1:n2) - r(i)
         rlmo2(i, 1:nlmo2) = rlmo2(i, 1:nlmo2) - r(i)
      end do
! save orig coord
      icoord0 = icoord

!     Drude and Gauss ini
      np = 4 ! # of prims
      fall = 0.50d0
      fhyd = 3.00d0
      p1 = 4.00   ! just for convenience so that r0scal is around 1
      do i = 1, n1
         f = fall
         if (at1(i) .le. 2) f = fhyd
         fc = par_drude_fc*gam(at1(i))
         qdr1(i) = -f*(fc*alp1(i))**0.50d0
         alp1(i) = alp1(i)**(1./3.)           ! make radius
         if (at1(i) .eq. 1) then
            aa = p1*r0scal(at1(i))/alp1(i)    ! Gaussian density
            aa = aa*(1.+cn1(i))**(1./2.)
         else
            aa = p1*r0scal(at1(i))/alp1(i)    ! Gaussian density
            aa = aa*(1.+cn1(i))**(-1./3.)
         end if
         dum1(i) = aa
         call setsto4(np, 1, 1, aa, e1, c1)
         den1(1, 1:np, i) = e1(1:np)            ! exp
         den1(2, 1:np, i) = c1(1:np)            ! contract
      end do
      do i = 1, n2
         f = fall
         if (at2(i) .le. 2) f = fhyd
         fc = par_drude_fc*gam(at2(i))
         qdr2(i) = -f*(fc*alp2(i))**0.50d0
         alp2(i) = alp2(i)**(1./3.)           ! make radius
         if (at2(i) .eq. 1) then
            aa = p1*r0scal(at2(i))/alp2(i)    ! Gaussian density
            aa = aa*(1.+cn2(i))**(1./2.)
         else
            aa = p1*r0scal(at2(i))/alp2(i)    ! Gaussian density
            aa = aa*(1.+cn2(i))**(-1./3.)
         end if
         dum2(i) = aa
         call setsto4(np, 1, 1, aa, e1, c1)
         den2(1, 1:np, i) = e1(1:np)            ! exp
         den2(2, 1:np, i) = c1(1:np)            ! contract
      end do

!     ccccccccccccccccccccccccccccccccccccccccccccccccc
!                       some output
!     ccccccccccccccccccccccccccccccccccccccccccccccccc

      if (debug) then
      write(*,*) '  #   Z           coordinates (CMA)   CN(D4)    qM/CM5','   alp   dexpo'
         do i = 1, n1
            write (*, '(2i4,3f10.5,3f6.2,3f7.3)')&
          &i, at1(i), xyz1(1:3, i), cn1(i), q1(i), q(i), alp1(i), dum1(i)
         end do
         write (*, *) 'LMOcent on atom     coordinates    charge'
         do i = 1, nlmo1
            write (*, '(2i7,3f10.5,5x,f7.3)') i, lmo1(i), rlmo1(1:4, i)
         end do
         write (*, *) &
        &'  #   Z           coordinates (CMA)   CN(D4)     qM/CM5',&
        &'   alp   dexpo'
         do i = 1, n2
            write (*, '(2i4,3f10.5,3f6.2,3f7.3)')&
           &i, at2(i), xyz2(1:3, i), cn2(i), q2(i), q(i + n1), alp2(i), dum2(i)
         end do
         write (*, *) 'LMOcent on atom     coordinates    charge'
         do i = 1, nlmo2
            write (*, '(2i7,3f10.5,5x,f7.3)') i, lmo2(i), rlmo2(1:4, i)
         end do
      end if

      write (*, '(''  LUMO energy 1 (read)           :'',F8.3)')&
      &elumo(1)*27.212
      write (*, '(''  HOMO energy 1 (read)           :'',F8.3)')&
      &ehomo(1)*27.212

      write (*, '(''  LUMO energy 2 (read)           :'',F8.3)')&
      &elumo(2)*27.212
      write (*, '(''  HOMO energy 2 (read)           :'',F8.3)')&
      &ehomo(2)*27.212

   end subroutine

   subroutine setsto4(nprim, n, l, zeta, expo, cont)
      integer, intent(in) :: n, l
      integer, intent(inout) :: nprim
      real(wp), intent(in) :: zeta
      real(wp), intent(out) :: expo(nprim), cont(nprim)

      real(wp) :: ALLC(4, 5, 4), ALLZ(4, 5, 4)
      integer :: iam
      real(wp) :: DEX(-1:96)
      integer :: I, J, DEX2
      real(wp) :: XNORM, PI

      nprim = 4

      PI = 4.D0*ATAN(1.D0)
      DO I = -1, 10
         if (I .lt. 2) then
            DEX2 = 1
         else
            DEX2 = 1
            do j = 1, I, 2
               DEX2 = DEX2*j
            end do
         end if
         DEX(I) = DEX2
      END DO
!                      1S
      ALLZ(1, 1, 1) = 5.216844534
      ALLZ(2, 1, 1) = 9.546182760D-1
      ALLZ(3, 1, 1) = 2.652034102D-1
      ALLZ(4, 1, 1) = 8.801862774D-2
      ALLC(1, 1, 1) = 5.675242080D-2
      ALLC(2, 1, 1) = 2.601413550D-1
      ALLC(3, 1, 1) = 5.328461143D-1
      ALLC(4, 1, 1) = 2.916254405D-1
!                      2S
      ALLZ(1, 2, 1) = 1.161525551D1
      ALLZ(2, 2, 1) = 2.000243111
      ALLZ(3, 2, 1) = 1.607280687D-1
      ALLZ(4, 2, 1) = 6.125744532D-2
      ALLC(1, 2, 1) = -1.198411747D-2
      ALLC(2, 2, 1) = -5.472052539D-2
      ALLC(3, 2, 1) = 5.805587176D-1
      ALLC(4, 2, 1) = 4.770079976D-1
!                     2P
      ALLZ(1, 2, 2) = 1.798260992
      ALLZ(2, 2, 2) = 4.662622228D-1
      ALLZ(3, 2, 2) = 1.643718620D-1
      ALLZ(4, 2, 2) = 6.543927065D-2
      ALLC(1, 2, 2) = 5.713170255D-2
      ALLC(2, 2, 2) = 2.857455515D-1
      ALLC(3, 2, 2) = 5.517873105D-1
      ALLC(4, 2, 2) = 2.632314924D-1
!                      3S
      ALLZ(1, 3, 1) = 1.513265591
      ALLZ(2, 3, 1) = 4.262497508D-1
      ALLZ(3, 3, 1) = 7.643320863D-2
      ALLZ(4, 3, 1) = 3.760545063D-2
      ALLC(1, 3, 1) = -3.295496352D-2
      ALLC(2, 3, 1) = -1.724516959D-1
      ALLC(3, 3, 1) = 7.518511194D-1
      ALLC(4, 3, 1) = 3.589627317D-1
!                     3P
      ALLZ(1, 3, 2) = 1.853180239
      ALLZ(2, 3, 2) = 1.915075719D-1
      ALLZ(3, 3, 2) = 8.655487938D-2
      ALLZ(4, 3, 2) = 4.184253862D-2
      ALLC(1, 3, 2) = -1.434249391D-2
      ALLC(2, 3, 2) = 2.755177580D-1
      ALLC(3, 3, 2) = 5.846750879D-1
      ALLC(4, 3, 2) = 2.144986514D-1
!                      4S
      ALLZ(1, 4, 1) = 3.242212833D-1
      ALLZ(2, 4, 1) = 1.663217177D-1
      ALLZ(3, 4, 1) = 5.081097451D-2
      ALLZ(4, 4, 1) = 2.829066600D-2
      ALLC(1, 4, 1) = -1.120682822D-1
      ALLC(2, 4, 1) = -2.845426863D-1
      ALLC(3, 4, 1) = 8.909873788D-1
      ALLC(4, 4, 1) = 3.517811205D-1
!                     4P
      ALLZ(1, 4, 2) = 1.492607880
      ALLZ(2, 4, 2) = 4.327619272D-1
      ALLZ(3, 4, 2) = 7.553156064D-2
      ALLZ(4, 4, 2) = 3.706272183D-2
      ALLC(1, 4, 2) = -6.035216774D-3
      ALLC(2, 4, 2) = -6.013310874D-2
      ALLC(3, 4, 2) = 6.451518200D-1
      ALLC(4, 4, 2) = 4.117923820D-1
!                      5S
      ALLZ(1, 5, 1) = 8.602284252D-1
      ALLZ(2, 5, 1) = 1.189050200D-1
      ALLZ(3, 5, 1) = 3.446076176D-2
      ALLZ(4, 5, 1) = 1.974798796D-2
      ALLC(1, 5, 1) = 1.103657561D-2
      ALLC(2, 5, 1) = -5.606519023D-1
      ALLC(3, 5, 1) = 1.179429987
      ALLC(4, 5, 1) = 1.734974376D-1
!                     5P
      ALLZ(1, 5, 2) = 3.962838833D-1
      ALLZ(2, 5, 2) = 1.838858552D-1
      ALLZ(3, 5, 2) = 4.943555157D-2
      ALLZ(4, 5, 2) = 2.750222273D-2
      ALLC(1, 5, 2) = -1.801459207D-2
      ALLC(2, 5, 2) = -1.360777372D-1
      ALLC(3, 5, 2) = 7.533973719D-1
      ALLC(4, 5, 2) = 3.409304859D-1
! 3d
      ALLZ(1, 3, 3) = 9.185846715D-1
      ALLZ(2, 3, 3) = 2.920461109D-1
      ALLZ(3, 3, 3) = 1.187568890D-1
      ALLZ(4, 3, 3) = 5.286755896D-2
      ALLC(1, 3, 3) = 5.799057705D-2
      ALLC(2, 3, 3) = 3.045581349D-1
      ALLC(3, 3, 3) = 5.601358038D-1
      ALLC(4, 3, 3) = 2.432423313D-1
! 4d
      ALLZ(1, 4, 3) = 1.995825422
      ALLZ(2, 4, 3) = 1.823461280D-1
      ALLZ(3, 4, 3) = 8.197240896D-2
      ALLZ(4, 4, 3) = 4.000634951D-2
      ALLC(1, 4, 3) = -2.816702620D-3
      ALLC(2, 4, 3) = 2.177095871D-1
      ALLC(3, 4, 3) = 6.058047348D-1
      ALLC(4, 4, 3) = 2.717811257D-1
! 5d
      ALLZ(1, 5, 3) = 4.230617826D-1
      ALLZ(2, 5, 3) = 8.293863702D-2
      ALLZ(3, 5, 3) = 4.590326388D-2
      ALLZ(4, 5, 3) = 2.628744797D-2
      ALLC(1, 5, 3) = -2.421626009D-2
      ALLC(2, 5, 3) = 3.937644956D-1
      ALLC(3, 5, 3) = 5.489520286D-1
      ALLC(4, 5, 3) = 1.190436963D-1

! 4f (the add. factors are required because the normal. function
!     below is not correct for f)
      ALLZ(1, 4, 4) = 0.05691670217
      ALLZ(2, 4, 4) = 0.2074585819
      ALLZ(3, 4, 4) = 0.09298346885
      ALLZ(4, 4, 4) = 0.04473508853
      ALLC(1, 4, 4) = 0.05902730589/1.063832358490576**0.5
      ALLC(2, 4, 4) = 0.3191828952/1.063832358490576**0.5
      ALLC(3, 4, 4) = 0.5639423893/1.063832358490576**0.5
      ALLC(4, 4, 4) = 0.2284796537/1.063832358490576**0.5

      iam = l - 1
      DO J = 1, nprim
         cont(J) = ALLC(J, N, L)
         expo(J) = ALLZ(J, N, L)*zeta**2
         XNORM = (2.D0*EXPO(J)/PI)**0.75D0*(4.D0*EXPO(J))**(IAM/2.D0)/&
               &SQRT(DEX(2*IAM - 1))
         cont(j) = cont(j)*xnorm
      END DO

   end subroutine setsto4

   subroutine c6_alp0(nat, at, xyz, cn, q, azero, c6ab)

      integer, intent(in) :: nat, at(nat)
      real(wp), intent(in) :: xyz(3, nat)
      !> Coordination number of every atom and charge.
      real(wp), intent(in) :: cn(nat), q(nat)
      !> C6AB Coefficients
      real(wp), intent(out) :: c6ab(nat, nat)
      !> Polarizabilities
      real(wp), intent(out) :: azero(nat)

      !> Molecular Structure information.
      type(TMolecule) :: mol

      type(TDispersionModel) :: dispm

      !> Charge scaling height.
      real(wp) :: g_a = 3.0_wp

      !> Charge scaling steepness.
      real(wp) :: g_c = 2.0_wp

      !> Exponent for the Gaussian weighting.
      real(wp) :: wf = 6.0_wp

      integer :: max_ref, i, ii, ia, k, iz

      real(wp) :: aw(23, nat), covcn(nat)

      real(wp), allocatable :: zetavec(:, :), zetadcn(:, :), zetadq(:, :)
      real(wp), allocatable :: zerovec(:, :), zerodcn(:, :), zerodq(:, :)
      real(wp), allocatable :: dc6dcn(:, :), dc6dq(:, :), gw(:), c6ref(:,:), zetvec(:)

      integer :: dispdim

      aw = 0.0_wp
      call init(mol, at, xyz)

      call newD4Model(dispm, g_a, g_c, p_refq_gfn2xtb)
      call d4dim(dispm, nat, at, dispdim) !To get disp dimension

      max_ref = maxval(dispm%nref(at))
      allocate (zetavec(max_ref, nat), zetadcn(max_ref, nat), zetadq(max_ref, nat), &
      &         zerovec(max_ref, nat), zerodcn(max_ref, nat), zerodq(max_ref, nat), &
      &         dc6dcn(nat, nat), dc6dq(nat, nat), gw(dispdim), c6ref(dispdim, dispdim), &
      &         zetvec(dispdim), source=0.0_wp)

      call weight_references(dispm, nat, at, g_a, g_c, wf, q, cn, zeff, gamma, &
         & zetavec, zerovec, zetadcn, zerodcn, zetadq) !generated zetavecotrs
      call get_atomic_c6(dispm, nat, mol%at, zetavec, zetadcn, zetadq, &
         & c6ab, dc6dcn, dc6dq) !c6ab are generated here
      call d4(dispm, nat, dispdim, at, wf, g_a, g_c, cn, gw, c6ref) !To get gw

      k = 0
      do i = 1, nat
         ia = at(i)
         iz = zeff(ia)
         do ii = 1, dispm%nref(ia)
            k = k + 1
            zetvec(k) = gw(k)*zeta(g_a, gamma(ia)*g_c, dispm%q(ii, ia) + iz, q(i) + iz)
            aw(:, i) = aw(:, i) + zetvec(k)*dispm%alpha(:, ii, ia)
         end do
      end do

      azero(:) = aw(1, :)

      deallocate(dispm%atoms,dispm%nref,dispm%ncount,dispm%cn,dispm%q,dispm%alpha,dispm%c6)
      dispm%g_a = 0.0_wp
      dispm%g_c = 0.0_wp

      deallocate(zetavec,zetadcn,zetadq,zerovec,zerodcn,zerodq,dc6dcn,dc6dq,gw,c6ref,zetvec)

   end subroutine c6_alp0

end module xtb_iff_iffini
