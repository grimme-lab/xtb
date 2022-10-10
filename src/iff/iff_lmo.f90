! This file is part of xtb.
!
! Copyright (C) 2022 Stefan Grimme, Christoph Plett
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

module xtb_iff_ifflmo
   use xtb_docking_param
   use xtb_axis, only: axisvec
   use xtb_splitparam, only: atmass
   use mctcpar_atomic_masses
   use xtb_setmod
   implicit none

   private
   public :: setlmo

contains

   subroutine trflmo(mol, n, nl, c, cl)
      integer, intent(in) :: n
      integer, intent(in) :: nl
      integer, intent(in) :: mol
      real(wp), intent(in) :: c(3, n)
      real(wp), intent(out) :: cl(4, n*10)

      integer :: i, iat, jat, kat
      real(wp) :: r3(3), r, r1, r2

      do i = 1, nl
!        LMO on atom
         if (lmoatom(1, i, mol) .eq. 1) then
            iat = lmoatom(2, i, mol)
            cl(1:3, i) = c(1:3, iat)
         end if
!        LP or pi
         if (lmoatom(1, i, mol) .eq. 4) then
            iat = lmoatom(2, i, mol)
            jat = lmoatom(3, i, mol)
            kat = lmoatom(4, i, mol)
            call trflp(n, c, iat, jat, kat,&
     &           lmoint(1, i, mol), lmoint(2, i, mol), lmoint(3, i, mol), r3)
            cl(1:3, i) = r3(1:3)
         end if
!        sigma
         if (lmoatom(1, i, mol) .eq. 2) then
            iat = lmoatom(2, i, mol)
            jat = lmoatom(3, i, mol)
            r3(1:3) = c(1:3, iat) - c(1:3, jat)
            r = dsqrt(r3(1)**2 + r3(2)**2 + r3(3)**2)
            r1 = lmoint(1, i, mol)
            r2 = lmoint(2, i, mol)
            if (r2 .gt. r) then
               cl(1:3, i) = c(1:3, iat) + r1*r3(1:3)/r
            else
               cl(1:3, i) = c(1:3, iat) - r1*r3(1:3)/r
            end if
         end if
      end do

   end subroutine trflmo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! important LMO prep routine for anisotropic ES
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setlmo(mol, n, nl, at, xyz, alp0, q, cn, c, lmo)

      integer, intent(in) :: mol, n
      real(wp), intent(in) :: xyz(3, n)
      real(wp), intent(in) :: alp0(n)
      real(wp), intent(in) :: q(n)
      real(wp), intent(in) :: cn(n)
      integer, intent(in) :: at(n)
      real(wp), intent(inout) :: c(4, 10*n)
      integer, intent(inout) :: lmo(10*n), nl

      integer :: i, j, ii, ia1, imin, ind(n), ltmp(10*n), nlp, npi, nsig, ia2, ia3
      real(wp) :: ex, ey, ez, rmin, f, norm(n), lpcharge, dip(3), cnfac, di, ff, rr1
      real(wp) :: d2,alp,p6pi,p6sig,r1,r2,qq,dex,c1,c2,ex1,ey1,ez1,rr,rr2
      real(wp) :: r0_free, ctmp(4, 10*n), prxs(3, 3), aa, bb, cc, vtmp(3), vnew(3)
      real(wp) :: s1, s2
      logical :: sigbond, pibond

      lmoatom(1, 1:maxlmo, mol) = 0
      ltmp = lmo
      ctmp = c
      c = 0

      if (n .eq. 1) then ! if just an atom, exit here (LPs break symmetry)
         nl = 0
         return
      end if

      npi = 0
      nlp = 0
      nsig = 0
      i = 0
      do ii = 1, nl
         if (n .eq. 1) ltmp(ii) = 2  ! a single atom always has only LP
         ex = ctmp(1, ii)
         ey = ctmp(2, ii)
         ez = ctmp(3, ii)
         rmin = 1.d+42
         do ia1 = 1, n
            norm(ia1) = sqrt((xyz(1, ia1) - ex)**2 + (xyz(2, ia1) - ey)**2&
                         &+ (xyz(3, ia1) - ez)**2)
            ind(ia1) = ia1
         end do
         call Qsort(norm, 1, n, ind)
         ia3 = 0
         ia1 = ind(1)
         ia2 = ind(1)
         if (n .gt. 1) ia2 = ind(2)
         if (n .gt. 2) ia3 = ind(3)
         if (ltmp(ii) .eq. 2 .and. lpatom(at(ia1)) .eq. 1) then
!           put it on the line along atom...LP center
            i = i + 1
            if (i .gt. maxlmo) stop 'too many LMO'
            f = -par_pos1_lp
            ex1 = xyz(1, ia1) + f*(xyz(1, ia1) - ex)
            ey1 = xyz(2, ia1) + f*(xyz(2, ia1) - ey)
            ez1 = xyz(3, ia1) + f*(xyz(3, ia1) - ez)
            c(1, i) = ex1
            c(2, i) = ey1
            c(3, i) = ez1
            c(4, i) = -0.010*par_chrg_lp*alp0(ia1)/r0_atom(at(ia1))
            !  speical hack for heavy elements
            if (at(ia1) .gt. 10)&
           &c(4, i) = -0.010*par_chrg_lp*alp0(ia1)/r0_atom(at(ia1))/1.8 
            lmo(i) = ia1 ! important overwrite of LMO type info by atom number to which LMO i is attached
            lmoatom(1, i, mol) = 4
            lmoatom(2, i, mol) = ia1
            lmoatom(3, i, mol) = ia2
            lmoatom(4, i, mol) = ia3
            call internal4(n, xyz, c(1, i), ia1, ia2, ia3,&
                &lmoint(1, i, mol), lmoint(2, i, mol), lmoint(3, i, mol))
            ! to distinguish LPs, set lmoatom(3:4,...) = 0
            lmoatom(3, i, mol) = 0
            lmoatom(4, i, mol) = 0
            !!!
            if (i .gt. maxlmo) stop 'too many LMO'
            i = i + 1
            f = -par_pos2_lp                             
            ! along the line but close to nuc (or at nuc if par_pos2_lp is 0)
            ! same side a + charge
            ex1 = xyz(1, ia1) + f*(xyz(1, ia1) - ex)
            ey1 = xyz(2, ia1) + f*(xyz(2, ia1) - ey)
            ez1 = xyz(3, ia1) + f*(xyz(3, ia1) - ez)
            c(1, i) = ex1
            c(2, i) = ey1
            c(3, i) = ez1
            c(4, i) = -c(4, i - 1)
            lmo(i) = ia1
            lmoatom(1, i, mol) = 1
            lmoatom(2, i, mol) = ia1
            nlp = nlp + 1
         end if

         sigbond = (sigatom(at(ia1)) .eq. 1 .or. sigatom(at(ia2)) .eq. 1)&   
                   &.and. (ltmp(ii) .eq. 1) .and. (n .ne. 1) .and. at(ia1) .ne. 1&
                   &.and. at(ia2) .ne. 1
         pibond = ltmp(ii) .ge. 3
         if (sigbond .or. pibond) then       ! sigma and pi part -
            i = i + 1
            if (i .gt. maxlmo) stop 'too many LMO'
            alp = (alp0(ia1)*alp0(ia2))**0.5 ! average of two closest atom values
            r0_free = (r0_atom(at(ia1))*r0_atom(at(ia2)))**0.5
            c(1:3, i) = ctmp(1:3, ii)         ! the - charge point at LMO center
            lmo(i) = ia1
            ! rest is + charge position
            r1 = sqrt((ex - xyz(1, ia1))**2 + (ey - xyz(2, ia1))**2&
                  &+ (ez - xyz(3, ia1))**2)
            r2 = sqrt((ex - xyz(1, ia2))**2 + (ey - xyz(2, ia2))**2&
                  &+ (ez - xyz(3, ia2))**2)
            c1 = r1/(r1 + r2)
            c2 = r2/(r1 + r2)
            if (sigbond) then               ! sigma
               rr = sqrt((xyz(1, ia1) - xyz(1, ia2))**2 +&
                      &(xyz(2, ia1) - xyz(2, ia2))**2 +&
                      &(xyz(3, ia1) - xyz(3, ia2))**2)
               ex1 = xyz(1, ia2) - xyz(1, ia1)
               ey1 = xyz(2, ia2) - xyz(2, ia1)
               ez1 = xyz(3, ia2) - xyz(3, ia1)
               qq = 0.01*par_chrg_sig*alp/r0_free
               c(4, i) = -qq
               lmoatom(1, i, mol) = 2
               lmoatom(2, i, mol) = ia1
               lmoatom(3, i, mol) = ia2
               lmoint(1, i, mol) = r1
               lmoint(2, i, mol) = r2
               i = i + 1
               if (i .gt. maxlmo) stop 'too many LMO'
               c(1, i) = ex - par_pos_sig*ex1/rr
               c(2, i) = ey - par_pos_sig*ey1/rr
               c(3, i) = ez - par_pos_sig*ez1/rr
               c(4, i) = qq*c1
               lmo(i) = ia1
               rr1 = sqrt((xyz(1, ia1) - c(1, i))**2 +&
                       &(xyz(2, ia1) - c(2, i))**2 +&
                       &(xyz(3, ia1) - c(3, i))**2)
               rr2 = sqrt((xyz(1, ia2) - c(1, i))**2 +&
                       &(xyz(2, ia2) - c(2, i))**2 +&
                       &(xyz(3, ia2) - c(3, i))**2)
               lmoatom(1, i, mol) = 2
               lmoatom(2, i, mol) = ia1
               lmoatom(3, i, mol) = ia2
               lmoint(1, i, mol) = rr1
               lmoint(2, i, mol) = rr2
               i = i + 1
               if (i .gt. maxlmo) stop 'too many LMO'
               c(1, i) = ex + par_pos_sig*ex1/rr
               c(2, i) = ey + par_pos_sig*ey1/rr
               c(3, i) = ez + par_pos_sig*ez1/rr
               c(4, i) = qq*c2
               lmo(i) = ia2
               rr2 = sqrt((xyz(1, ia2) - c(1, i))**2 +&
                       &(xyz(2, ia2) - c(2, i))**2 +&
                       &(xyz(3, ia2) - c(3, i))**2)
               rr1 = sqrt((xyz(1, ia1) - c(1, i))**2 +&
                       &(xyz(2, ia1) - c(2, i))**2 +&
                       &(xyz(3, ia1) - c(3, i))**2)
               lmoatom(1, i, mol) = 2
               lmoatom(2, i, mol) = ia2
               lmoatom(3, i, mol) = ia1
               lmoint(1, i, mol) = rr2
               lmoint(2, i, mol) = rr1
               nsig = nsig + 1
            else                          ! pi
               ff = 1.0
               if (ltmp(ii) .eq. 4) ff = 3.0   ! deloc pi
               qq = ff*0.01*par_chrg_pi*alp/r0_free
               c(4, i) = -qq
               lmoatom(1, i, mol) = 4
               lmoatom(2, i, mol) = ia1
               lmoatom(3, i, mol) = ia2
               lmoatom(4, i, mol) = ia3
               call internal4(n, xyz, c(1, i), ia1, ia2, ia3,&
                   &lmoint(1, i, mol), lmoint(2, i, mol), lmoint(3, i, mol))
               i = i + 1
               if (i .gt. maxlmo) stop 'too many LMO'
               c(1:3, i) = xyz(1:3, ia1)      ! + on atoms
               c(4, i) = qq*c1
               lmo(i) = ia1
               lmoatom(1, i, mol) = 1
               lmoatom(2, i, mol) = ia1
               i = i + 1
               if (i .gt. maxlmo) stop 'too many LMO'
               c(1:3, i) = xyz(1:3, ia2)
               c(4, i) = qq*c2
               lmo(i) = ia2
               lmoatom(1, i, mol) = 1
               lmoatom(2, i, mol) = ia2
               npi = npi + 1
            end if
         end if
!     next LMO
      end do

! special treatment for linear molecules to avoid the symmetry breaking
! this is done by projecting the off-center charges onto the molecular axis, thus
! done after positioning of the charge centers (most importantly the LPs)
      aa = 0.0d0
      bb = 0.0d0
      cc = 0.0d0
      prxs = 0.0d0
      atmass = atomic_mass(at)
      !Both yield the same result (delete axisvec2 here)
      call axisvec(n, at, xyz, aa, bb, cc, prxs)

      if (cc .lt. 1.d-8) then  ! check, whether molecule is linear
         do j = 1, i
            if (lmoatom(1, j, mol) .le. 1) cycle
            vtmp = 0.0d0
            vnew = 0.0d0
            ia1 = lmoatom(2, j, mol)
            vtmp(1:3) = c(1:3, j) - xyz(1:3, ia1)
            s2 = vtmp(1)**2
            s2 = s2 + vtmp(2)**2
            s2 = s2 + vtmp(3)**2
            s2 = sqrt(s2)
            ! project off-center charge onto molecular axis
            s1 = vtmp(1)*prxs(1, 1)
            s1 = vtmp(2)*prxs(2, 1)
            s1 = vtmp(3)*prxs(3, 1)
            vnew(1) = s1*prxs(1, 1)
            vnew(2) = s1*prxs(2, 1)
            vnew(3) = s1*prxs(3, 1)
            c(1:3, j) = xyz(1:3, ia1) + vnew(1:3)

            ! rescale charges
            s1 = vnew(1)**2
            s1 = s1 + vnew(2)**2
            s1 = s1 + vnew(3)**2
            s1 = sqrt(s1)
            aa = vtmp(1)*vnew(1)
            aa = aa + vtmp(2)*vnew(2)
            aa = aa + vtmp(3)*vnew(3)
            aa = aa/(s1*s2) ! cosine
            aa = 2.0d0 - aa**2 ! 1+ sine**2
            ia2 = lmoatom(2, j + 1, mol)
            if (ia1 .eq. ia2) then
               if (lmoatom(3, j, mol) .eq. 0) then   !LP case
                  c(4, j) = c(4, j)*aa
                  c(4, j + 1) = c(4, j + 1)*aa
                  cycle
               end if
               if (lmoatom(3, j + 1, mol) .eq. 0) then ! deloc pi case
                  c(4, j) = c(4, j)*aa
                  c(4, j + 1) = c(4, j + 1)*aa
                  c(4, j + 2) = c(4, j + 2)*aa
                  cycle
               end if
               ! sigma or pi
               c(4, j) = c(4, j)*aa
               c(4, j + 1) = c(4, j + 1)*aa
               c(4, j + 2) = c(4, j + 2)*aa
            end if
         end do
      end if

      nl = i
      if(set%verbose) write (*, '(''# of LP/pi/sigma LMOs    :'',3i5)') nlp, npi, nsig
      if(set%verbose) write (*, '(''# of off-center charges  :'', i5)') nl

      dip = 0
      do i = 1, nl
         dip(1:3) = dip(1:3) + c(1:3, i)*c(4, i)
      end do
      do i = 1, n
         dip(1:3) = dip(1:3) + xyz(1:3, i)*q(i)
      end do
      di = sqrt(dip(1)**2 + dip(2)**2 + dip(3)**2)
      if(set%verbose) write (*, '(''dipole moment (read/calc):'',2F8.3)') dipol(mol), di

   end subroutine setlmo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! internal coords (R,angle,dihedral) for LP and pi LMO center
! as defined wrt next neighbors, i.e., ijkl
! c(3) is the LMO center
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine internal4(n, xyz, c, j, k, l, r, alp, bet)
      integer, intent(in) :: n, j, k, l
      real(wp), intent(in) :: xyz(3, n), c(3)
      real(wp), intent(out) :: r, alp, bet

      real(wp) :: xyz2(3, 4)

      alp = 10.0 ! nonsense value
      bet = 10.0 ! nonsense value
      r = sqrt((c(1) - xyz(1, j))**2 +&
              &(c(2) - xyz(2, j))**2 +&
              &(c(3) - xyz(3, j))**2)

      if (k .gt. 0) then
         xyz2(1:3, 1) = c(1:3)
         xyz2(1:3, 2) = xyz(1:3, j)
         xyz2(1:3, 3) = xyz(1:3, k)
         call xbangle(XYZ2, 1, 2, 3, alp)
         if (l .gt. 0) then
            xyz2(1:3, 4) = xyz(1:3, l)
            call xDIHED(XYZ2, 1, 2, 3, 4, bet)
         end if
      end if

   end subroutine internal4

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine trflp(n, coord, nai, nbi, nci, geo1i, geo2i, geo3i, ci)
      integer, intent(in) :: n, nai, nbi, nci
      real(wp), intent(in) :: coord(3, n)
      real(wp), intent(in) :: geo1i, geo2i, geo3i
      real(wp), intent(out) :: ci(3)

      real(wp) :: xb, yb, zb, xpa, xpb, xa, ya, za, xyb, xrd, yqd, xqd, zqd
      real(wp) :: costh, sinth, cosph, sinph, cosd, sind, xd, yd, zd, xpd
      real(wp) :: ypd, zpd, rbc, cosa, ypa, coskh, sinkh, sina, yza, zqa, xqa
      integer :: ma, mb, mc, k

      COSA = COS(GEO2I)
      MB = NBI
      MC = NAI
      XB = COORD(1, MB) - COORD(1, MC)
      YB = COORD(2, MB) - COORD(2, MC)
      ZB = COORD(3, MB) - COORD(3, MC)
      RBC = 1.0D00/DSQRT(XB*XB + YB*YB + ZB*ZB)
!
!     THE ATOMS ARE NOT COLLINEAR
!
      MA = NCI
      XA = COORD(1, MA) - COORD(1, MC)
      YA = COORD(2, MA) - COORD(2, MC)
      ZA = COORD(3, MA) - COORD(3, MC)
!
!     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
!     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
!
      XYB = DSQRT(XB*XB + YB*YB)
      K = -1
      IF (XYB .LE. 0.1D00) then
         XPA = ZA
         ZA = -XA
         XA = XPA
         XPB = ZB
         ZB = -XB
         XB = XPB
         XYB = DSQRT(XB*XB + YB*YB)
         K = 1
      end if
!
!     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
!
      COSTH = XB/XYB
      SINTH = YB/XYB
      XPA = XA*COSTH + YA*SINTH
      YPA = YA*COSTH - XA*SINTH
      SINPH = ZB*RBC
      COSPH = DSQRT(ABS(1.D00 - SINPH*SINPH))
      XQA = XPA*COSPH + ZA*SINPH
      ZQA = ZA*COSPH - XPA*SINPH
!
!     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
!
      YZA = DSQRT(YPA*YPA + ZQA*ZQA)
      IF (YZA .LT. 1.D-4) stop 'lp coord gen error'
      COSKH = YPA/YZA
      SINKH = ZQA/YZA
!
!     COORDINATES :-   A=(XQA,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
!     NONE ARE NEGATIVE.
!     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
!
      SINA = SIN(GEO2I)
      SIND = -SIN(GEO3I)
      COSD = COS(GEO3I)
      XD = GEO1I*COSA
      YD = GEO1I*SINA*COSD
      ZD = GEO1I*SINA*SIND
!
!     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
!
      YPD = YD*COSKH - ZD*SINKH
      ZPD = ZD*COSKH + YD*SINKH
      XPD = XD*COSPH - ZPD*SINPH
      ZQD = ZPD*COSPH + XD*SINPH
      XQD = XPD*COSTH - YPD*SINTH
      YQD = YPD*COSTH + XPD*SINTH
      IF (K .GE. 1) then
         XRD = -ZQD
         ZQD = XQD
         XQD = XRD
      end if
      CI(1) = XQD + COORD(1, MC)
      CI(2) = YQD + COORD(2, MC)
      CI(3) = ZQD + COORD(3, MC)

   end subroutine trflp
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   subroutine xBANGLE(XYZ, I, J, K, ANGLE)
      integer, intent(in) :: I, J, K
      real(wp), intent(in) :: XYZ(3, *)
      real(wp), intent(out) :: ANGLE

      real(wp) :: D2IJ, D2JK, D2IK, XY, TEMP
! BANGLE CALCULATES THE ANGLE BETWEEN ATOMS I,J, AND K. THE
!        CARTESIAN COORDINATES ARE IN XYZ.
      D2IJ = (XYZ(1, I) - XYZ(1, J))**2 +&
            &(XYZ(2, I) - XYZ(2, J))**2 +&
            &(XYZ(3, I) - XYZ(3, J))**2
      D2JK = (XYZ(1, J) - XYZ(1, K))**2 +&
            &(XYZ(2, J) - XYZ(2, K))**2 +&
            &(XYZ(3, J) - XYZ(3, K))**2
      D2IK = (XYZ(1, I) - XYZ(1, K))**2 +&
            &(XYZ(2, I) - XYZ(2, K))**2 +&
            &(XYZ(3, I) - XYZ(3, K))**2
      XY = SQRT(D2IJ*D2JK + 1.d-14)
      TEMP = 0.5D0*(D2IJ + D2JK - D2IK)/XY
      IF (TEMP .GT. 1.0D0) TEMP = 1.0D0
      IF (TEMP .LT. -1.0D0) TEMP = -1.0D0
      ANGLE = ACOS(TEMP)
      RETURN
   end subroutine xBANGLE

   subroutine xDIHED(XYZ, I, J, K, L, ANGLE)
      real(wp), intent(in) :: XYZ(3, 4)
      integer, intent(in) :: I, J, K, L
      real(wp), intent(out) :: ANGLE

      real(wp) :: XI1, XJ1, XL1, YI1, YJ1, YL1, ZI1, ZJ1, ZL1
      real(wp) :: XI2, XJ2, XL2, YI2, YJ2, YL2, ZI2, ZJ2, ZL2
      real(wp) :: YI3, YL3
      real(wp) :: DIST, COSA, DDD, YXDIST, COSTH, SINTH, COSPH, SINPH

      real(wp), parameter :: PI = 3.14159265358979_wp
      XI1 = XYZ(1, I) - XYZ(1, K)
      XJ1 = XYZ(1, J) - XYZ(1, K)
      XL1 = XYZ(1, L) - XYZ(1, K)
      YI1 = XYZ(2, I) - XYZ(2, K)
      YJ1 = XYZ(2, J) - XYZ(2, K)
      YL1 = XYZ(2, L) - XYZ(2, K)
      ZI1 = XYZ(3, I) - XYZ(3, K)
      ZJ1 = XYZ(3, J) - XYZ(3, K)
      ZL1 = XYZ(3, L) - XYZ(3, K)
!      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
      DIST = SQRT(XJ1**2 + YJ1**2 + ZJ1**2 + 1.d-14)
      COSA = ZJ1/DIST
      IF (COSA .GT. 1.0D0) COSA = 1.0D0
      IF (COSA .LT. -1.0D0) COSA = -1.0D0
      DDD = 1.0D0 - COSA**2
      IF (DDD .GT. 0.0) then
         YXDIST = DIST*SQRT(DDD)
      end if
      IF (YXDIST .LE. 1.0D-9 .AND. DDD .GT. 0.0) then
         XI2 = XI1
         XL2 = XL1
         YI2 = YI1
         YL2 = YL1
         COSTH = COSA
         SINTH = 0.D0
      end if
      IF (YXDIST .GT. 1.0D-9 .OR. DDD .LE. 0.0) then
         COSPH = YJ1/YXDIST
         SINPH = XJ1/YXDIST
         XI2 = XI1*COSPH - YI1*SINPH
         XJ2 = XJ1*COSPH - YJ1*SINPH
         XL2 = XL1*COSPH - YL1*SINPH
         YI2 = XI1*SINPH + YI1*COSPH
         YJ2 = XJ1*SINPH + YJ1*COSPH
         YL2 = XL1*SINPH + YL1*COSPH
!      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
         COSTH = COSA
         SINTH = YJ2/DIST
      end if
      YI3 = YI2*COSTH - ZI1*SINTH
      YL3 = YL2*COSTH - ZL1*SINTH
      CALL xDANG(XL2, YL3, XI2, YI3, ANGLE)
      IF (ANGLE .LT. 0.) ANGLE = 2.0D0*PI + ANGLE
      IF (ANGLE .GE. 2.0d0*PI) ANGLE = 0.D0
      RETURN
   end subroutine xDIHED

   subroutine xDANG(A1, A2, B1, B2, RCOS)
      real(wp), intent(inout) :: A1, A2, B1, B2
      real(wp), intent(out) :: RCOS

      real(wp) :: ZERO, SINTH, COTH, ANORM, BNORM, COSTH
      real(wp), parameter :: PI = 3.14159265358979_wp
      ZERO = 1.0D-6
      IF (ABS(A1) .LT. ZERO .AND. ABS(A2) .LT. ZERO) then
         RCOS = 0.0D0
         return
      end if
      IF (ABS(B1) .LT. ZERO .AND. ABS(B2) .LT. ZERO) then
         RCOS = 0.0D0
         return
      end if
      ANORM = 1.0D0/SQRT(A1**2 + A2**2)
      BNORM = 1.0D0/SQRT(B1**2 + B2**2)
      A1 = A1*ANORM
      A2 = A2*ANORM
      B1 = B1*BNORM
      B2 = B2*BNORM
      SINTH = (A1*B2) - (A2*B1)
      COSTH = A1*B1 + A2*B2
      IF (COSTH .GT. 1.0D0) COSTH = 1.0D0
      IF (COSTH .LT. -1.0D0) COSTH = -1.0D0
      RCOS = ACOS(COSTH)
      IF (ABS(RCOS) .GE. 4.0D-4) then
         IF (SINTH .GT. 0.D0) RCOS = 2.0D0*PI - RCOS
         RCOS = -RCOS
         RETURN
      END IF
      RCOS = 0.0D0
      RETURN
   end subroutine xDANG

end module xtb_iff_ifflmo
