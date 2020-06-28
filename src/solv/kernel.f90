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

!> Generalized Born interaction kernels
module xtb_solv_kernel
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas, only : mctc_dot, mctc_gemv, mctc_symv
   use xtb_mctc_constants, only : fourpi, pi
   implicit none
   private

   public :: gbKernel
   public :: addBornMatSaltStill, addBornMatStill, addBornMatP16
   public :: addGradientSaltStill, addGradientStill, addGradientP16
   public :: addBornDerivSaltStill, addBornDerivStill


   !> Possible kernels for the generalized Born model
   type :: TGBKernelEnum

      !> Classical Still kernel
      integer :: still = 1

      !> P16 kernel by Lange (JCTC 2012, 8, 1999-2011)
      integer :: p16 = 2

   end type TGBKernelEnum

   !> Actual enumerator for the generalized Born kernels
   type(TGBKernelEnum), parameter :: gbKernel = TGBKernelEnum()


   !> P16 zeta parameter
   real(wp), parameter :: zetaP16 = 1.028_wp

   !> P16 zeta parameter over 16
   real(wp), parameter :: zetaP16o16 = zetaP16 / 16.0_wp


contains


subroutine addGradientSaltStill(nat, ntpair, ppind, ddpair, qat, kappa, &
      & brad, brdr, ionscr, discr, energy, gradient)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Number of all interacting pairs
   integer, intent(in) :: ntpair

   !> Index list
   integer, intent(in) :: ppind(:, :)

   !> Distances of all pairs
   real(wp), intent(in) :: ddpair(:, :)

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Debye screening length
   real(wp), intent(in) :: kappa

   !> Born radii
   real(wp), intent(in) :: brad(:)

   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), intent(in) :: brdr(:, :, :)

   !> Ion screening
   real(wp), intent(in) :: ionscr(:)

   !> Derivative of ion screening w.r.t. Born radii
   real(wp), intent(in) :: discr(:)

   !> Total Born solvation energy
   real(wp), intent(out) :: energy

   !> Deriatives of Born solvation energy
   real(wp), intent(inout) :: gradient(:, :)

   integer :: i, j, k, nnj
   integer :: kk
   real(wp), parameter :: a13=1._wp/3._wp
   real(wp), parameter :: a4=0.25_wp
   real(wp) :: aa, r2, fgb, fgb2, br3
   real(wp) :: qq, dd, expd, dfgb, dfgb2, dfgb3, egb, ap, bp, qfg
   real(wp) :: gg, expa
   real(wp) :: r0vdw, r01, r02, ar02
   real(wp) :: grddbi, grddbj
   real(wp) :: dr(3), r
   real(wp), allocatable :: grddb(:)

   allocate(grddb(nat), source = 0.0_wp )

   egb = 0._wp
   grddb(:) = 0._wp

   ! GB-SE energy and gradient

   ! compute energy and fgb direct and radii derivatives
   do kk = 1, ntpair
      r = ddpair(1, kk)
      r2 = r*r

      i = ppind(1, kk)
      j = ppind(2, kk)

      qq = qat(i)*qat(j)
      aa = brad(i)*brad(j)
      dd = a4*r2/aa
      expd = exp(-dd)
      dfgb = r2+aa*expd
      fgb = sqrt(dfgb)
      aa = kappa*fgb
      expa = exp(-aa)
      gg = (ionscr(i)+ionscr(j))*expa
      qfg = qq/fgb

      egb = egb + qfg*(gg-1.0_wp)

      dfgb3 = qfg*(gg*(1._wp+aa)-1.0_wp)/dfgb

      ap = (1._wp-a4*expd)*dfgb3
      dr = ap*ddpair(2:4, kk)
      gradient(:, i) = gradient(:, i) - dr
      gradient(:, j) = gradient(:, j) + dr

      qfg = qfg*expa
      bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
      grddbi = (brad(j)*bp+qfg*discr(i))
      grddbj = (brad(i)*bp+qfg*discr(j))
      grddb(i) = grddb(i) + grddbi
      grddb(j) = grddb(j) + grddbj

   enddo

   ! self-energy part
   do i = 1, nat
      gg = exp(-kappa*brad(i))
      aa = 2._wp*ionscr(i)*gg-1.0_wp
      qq = qat(i)/brad(i)
      egb = egb + 0.5_wp*qq*qat(i)*aa
      ap = aa-brad(i)*2._wp*(discr(i)+ionscr(i)*kappa)*gg
      grddbi = -0.5_wp*qq*qq*ap
      grddb(i) = grddb(i) + grddbi
   enddo

   ! contract with the Born radii derivatives
   call mctc_gemv(brdr, grddb, gradient, beta=1.0_wp)

   energy = egb

end subroutine addGradientSaltStill


subroutine addGradientStill(nat, ntpair, ppind, ddpair, qat, keps, &
      & brad, brdr, energy, gradient)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Number of all interacting pairs
   integer, intent(in) :: ntpair

   !> Index list
   integer, intent(in) :: ppind(:, :)

   !> Distances of all pairs
   real(wp), intent(in) :: ddpair(:, :)

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Dielectric screening
   real(wp), intent(in) :: keps

   !> Born radii
   real(wp), intent(in) :: brad(:)

   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), intent(in) :: brdr(:, :, :)

   !> Total Born solvation energy
   real(wp), intent(out) :: energy

   !> Deriatives of Born solvation energy
   real(wp), intent(inout) :: gradient(:, :)

   integer :: i, j, k, nnj
   integer :: kk
   real(wp), parameter :: a13=1._wp/3._wp
   real(wp), parameter :: a4=0.25_wp
   real(wp) :: aa, r2, fgb, fgb2, br3
   real(wp) :: qq, dd, expd, dfgb, dfgb2, dfgb3, egb, ap, bp, qfg
   real(wp) :: gg, expa
   real(wp) :: r0vdw, r01, r02, ar02
   real(wp) :: grddbi, grddbj
   real(wp) :: dr(3), r
   real(wp), allocatable :: grddb(:)

   allocate(grddb(nat), source = 0.0_wp )

   egb = 0._wp
   grddb(:) = 0._wp

   ! GB energy and gradient

   ! compute energy and fgb direct and radii derivatives
   do kk = 1, ntpair
      r = ddpair(1, kk)
      r2 = r*r

      i = ppind(1, kk)
      j = ppind(2, kk)

      ! dielectric scaling of the charges
      qq = qat(i)*qat(j)
      aa = brad(i)*brad(j)
      dd = a4*r2/aa
      expd = exp(-dd)
      fgb2 = r2+aa*expd
      dfgb2 = 1._wp/fgb2
      dfgb = sqrt(dfgb2)
      dfgb3 = dfgb2*dfgb*keps

      egb = egb + qq*keps*dfgb

      ap = (1._wp-a4*expd)*dfgb3
      dr = ap*ddpair(2:4, kk)
      gradient(:, i) = gradient(:, i) - dr*qq
      gradient(:, j) = gradient(:, j) + dr*qq

      bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
      grddbi = brad(j)*bp
      grddbj = brad(i)*bp
      grddb(i) = grddb(i) + grddbi*qq
      !gradient = gradient + brdr(:, :, i) * grddbi*qq
      grddb(j) = grddb(j) + grddbj*qq
      !gradient = gradient + brdr(:, :, j) * grddbj*qq

   enddo

   ! self-energy part
   do i = 1, nat
      bp = 1._wp/brad(i)
      qq = qat(i)*bp
      egb = egb + 0.5_wp*qat(i)*qq*keps
      grddbi = -0.5_wp*keps*qq*bp
      grddb(i) = grddb(i) + grddbi*qat(i)
      !gradient = gradient + brdr(:, :, i) * grddbi*qat(i)
   enddo

   ! contract with the Born radii derivatives
   call mctc_gemv(brdr, grddb, gradient, beta=1.0_wp)

   energy = egb

end subroutine addGradientStill


subroutine addGradientP16(nat, ntpair, ppind, ddpair, qat, keps, &
      & brad, brdr, energy, gradient)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Number of all interacting pairs
   integer, intent(in) :: ntpair

   !> Index list
   integer, intent(in) :: ppind(:, :)

   !> Distances of all pairs
   real(wp), intent(in) :: ddpair(:, :)

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Dielectric screening
   real(wp), intent(in) :: keps

   !> Born radii
   real(wp), intent(in) :: brad(:)

   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), intent(in) :: brdr(:, :, :)

   !> Total Born solvation energy
   real(wp), intent(out) :: energy

   !> Deriatives of Born solvation energy
   real(wp), intent(inout) :: gradient(:, :)

   integer :: iat, jat, kk
   real(wp) :: vec(3), r2, r1, ab, arg1, arg16, qq, fgb, fgb2, dfgb, dfgb2, egb
   real(wp) :: dEdbri, dEdbrj, dG(3), ap, bp, dS(3, 3)
   real(wp), allocatable :: dEdbr(:)

   allocate(dEdbr(nat), source = 0.0_wp )

   egb = 0._wp
   dEdbr(:) = 0._wp

   ! GB energy and gradient
   ! omp parallel do default(none) reduction(+:egb, gradient, dEdbr) &
   ! omp private(iat, jat, vec, r1, r2, ab, arg1, arg16, fgb, dfgb, dfgb2, ap, &
   ! omp& bp, qq, dEdbri, dEdbrj, dG, dS) &
   ! omp shared(keps, qat, ntpair, ddpair, ppind, brad)
   do kk = 1, ntpair
      vec(:) = ddpair(2:4, kk)
      r1 = ddpair(1, kk)
      r2 = r1*r1

      iat = ppind(1, kk)
      jat = ppind(2, kk)
      qq = qat(iat)*qat(jat)

      ab = sqrt(brad(iat) * brad(jat))
      arg1 = ab / (ab + zetaP16o16*r1) ! 1 / (1 + ζR/(16·ab))
      arg16 = arg1 * arg1 ! 1 / (1 + ζR/(16·ab))²
      arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))⁴
      arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))⁸
      arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))¹⁶

      fgb = r1 + ab*arg16
      dfgb = 1.0_wp / fgb
      dfgb2 = dfgb * dfgb

      egb = egb + qq*keps*dfgb

      ! (1 - ζ/(1 + Rζ/(16 ab))^17)/(R + ab/(1 + Rζ/(16 ab))¹⁶)²
      ap = (1.0_wp - zetaP16 * arg1 * arg16) * dfgb2
      dG(:) = ap * vec * keps / r1 * qq
      gradient(:, iat) = gradient(:, iat) - dG
      gradient(:, jat) = gradient(:, jat) + dG

      ! -(Rζ/(2·ab²·(1 + Rζ/(16·ab))¹⁷) + 1/(2·ab·(1 + Rζ/(16·ab))¹⁶))/(R + ab/(1 + Rζ/(16·ab))¹⁶)²
      bp = -0.5_wp*(r1 * zetaP16 / ab * arg1 + 1.0_wp) / ab * arg16 * dfgb2
      dEdbri = brad(jat) * bp * keps * qq
      dEdbrj = brad(iat) * bp * keps * qq
      dEdbr(iat) = dEdbr(iat) + dEdbri
      dEdbr(jat) = dEdbr(jat) + dEdbrj

   end do

   ! self-energy part
   do iat = 1, nat
      bp = 1._wp/brad(iat)
      qq = qat(iat)*bp
      egb = egb + 0.5_wp*qat(iat)*qq*keps
      dEdbri = -0.5_wp*keps*qq*bp
      dEdbr(iat) = dEdbr(iat) + dEdbri*qat(iat)
      !gradient = gradient + brdr(:, :, i) * dEdbri*qat(i)
   enddo

   ! contract with the Born radii derivatives
   call mctc_gemv(brdr, dEdbr, gradient, beta=1.0_wp)

   energy = egb

end subroutine addGradientP16


pure subroutine addBornDerivSaltStill(nat, ntpair, ppind, ddpair, qat, kappa, &
      & brad, brdr, ionscr, discr, gborn, dAmatdr, Afac)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Number of all interacting pairs
   integer, intent(in) :: ntpair

   !> Index list
   integer, intent(in) :: ppind(:, :)

   !> Distances of all pairs
   real(wp), intent(in) :: ddpair(:, :)

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Debye screening length
   real(wp), intent(in) :: kappa

   !> Born radii
   real(wp), intent(in) :: brad(:)

   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), intent(in) :: brdr(:, :, :)

   !> Ion screening
   real(wp), intent(in) :: ionscr(:)

   !> Derivative of ion screening w.r.t. Born radii
   real(wp), intent(in) :: discr(:)

   real(wp), intent(inout) :: dAmatdr(:, :, :)
   real(wp), intent(inout) :: Afac(:, :)
   real(wp), intent(out)   :: gborn

   integer :: i,j,k,nnj
   integer :: kk
   real(wp), parameter :: a13=1._wp/3._wp
   real(wp), parameter :: a4=0.25_wp
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp) :: aa,r2,fgb,fgb2,br3
   real(wp) :: qq,dd,expd,dfgb,dfgb2,dfgb3,ap,bp,qfg
   real(wp) :: gg,expa,aii,egb
   real(wp) :: r0vdw,r01,r02,ar02
   real(wp) :: grddbi,grddbj
   real(wp) :: dr(3),r

   egb = 0._wp

   ! GB-SE energy and dAmatdr

   ! compute energy and fgb direct and radii derivatives
   do kk = 1, ntpair
      r = ddpair(1,kk)
      r2 = r*r

      i = ppind(1,kk)
      j = ppind(2,kk)

      qq = qat(i)*qat(j)
      aa = brad(i)*brad(j)
      dd = a4*r2/aa
      expd = exp(-dd)
      fgb2 = r2+aa*expd
      fgb = sqrt(fgb2)
      dfgb2 = 1._wp/fgb2
      dfgb = sqrt(dfgb2)
      aa = kappa*fgb
      expa = exp(-aa)
      gg = (ionscr(i)+ionscr(j))*expa

      egb = egb + qq*dfgb*(gg-1.0_wp)

      dfgb3 = (gg*(1._wp+aa)-1.0_wp)*dfgb*dfgb2

      ap = (1._wp-a4*expd)*dfgb3
      dr = ap*ddpair(2:4,kk)
      dAmatdr(:,i,j) = dAmatdr(:,i,j) - dr*qat(i)
      dAmatdr(:,j,i) = dAmatdr(:,j,i) + dr*qat(j)
      Afac(:,i) = Afac(:,i) - dr*qat(j)
      Afac(:,j) = Afac(:,j) + dr*qat(i)

      qfg = dfgb*expa
      bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
      grddbi = (brad(j)*bp+qfg*discr(i))*qat(j)
      grddbj = (brad(i)*bp+qfg*discr(j))*qat(i)

      dAmatdr(:,:,i) = dAmatdr(:,:,i) + brdr(:,:,i)*grddbi
      dAmatdr(:,:,j) = dAmatdr(:,:,j) + brdr(:,:,j)*grddbj

   enddo

   ! self-energy part
   do i = 1, nat
      gg = exp(-kappa*brad(i))
      aa = 2._wp*ionscr(i)*gg-1.0_wp
      qq = qat(i)/brad(i)
      egb = egb + 0.5_wp*qq*qat(i)*aa
      ap = aa-brad(i)*2._wp*(discr(i)+ionscr(i)*kappa)*gg
      grddbi = -0.5_wp*qq*ap/brad(i)
      dAmatdr(:,:,i) = dAmatdr(:,:,i) + brdr(:,:,i)*grddbi
   enddo

   gborn = egb

end subroutine addBornDerivSaltStill


pure subroutine addBornDerivStill(nat, ntpair, ppind, ddpair, qat, keps, &
      & brad, brdr, gborn, dAmatdr, Afac)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Number of all interacting pairs
   integer, intent(in) :: ntpair

   !> Index list
   integer, intent(in) :: ppind(:, :)

   !> Distances of all pairs
   real(wp), intent(in) :: ddpair(:, :)

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Dielectric screening
   real(wp), intent(in) :: keps

   !> Born radii
   real(wp), intent(in) :: brad(:)

   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), intent(in) :: brdr(:, :, :)

   real(wp), intent(inout) :: dAmatdr(:, :, :)
   real(wp), intent(inout) :: Afac(:, :)
   real(wp), intent(out)   :: gborn

   integer :: i,j,k,nnj
   integer :: kk
   real(wp), parameter :: a13=1._wp/3._wp
   real(wp), parameter :: a4=0.25_wp
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp) :: aa,r2,fgb,fgb2,br3
   real(wp) :: qq,dd,expd,dfgb,dfgb2,dfgb3,ap,bp,qfg
   real(wp) :: gg,expa,aii,egb
   real(wp) :: r0vdw,r01,r02,ar02
   real(wp) :: grddbi,grddbj
   real(wp) :: dr(3),r

   egb = 0._wp

   ! GB energy and gradient

   ! compute energy and fgb direct and radii derivatives
   do kk = 1, ntpair
      r = ddpair(1,kk)
      r2 = r*r

      i = ppind(1,kk)
      j = ppind(2,kk)

      ! dielectric scaling of the charges
      qq = qat(i)*qat(j)*keps
      aa = brad(i)*brad(j)
      dd = a4*r2/aa
      expd = exp(-dd)
      fgb2 = r2+aa*expd
      dfgb2 = 1._wp/fgb2
      dfgb = sqrt(dfgb2)
      dfgb3 = dfgb2*dfgb*keps

      egb = egb + qq*dfgb

      ap = (1._wp-a4*expd)*dfgb3
      dr = ap*ddpair(2:4,kk)
      dAmatdr(:,i,j) = dAmatdr(:,i,j) - dr*qat(i)
      dAmatdr(:,j,i) = dAmatdr(:,j,i) + dr*qat(j)
      Afac(:,i) = Afac(:,i) - dr*qat(j)
      Afac(:,j) = Afac(:,j) + dr*qat(i)

      bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
      grddbi = brad(j)*bp
      grddbj = brad(i)*bp

      dAmatdr(:,:,j) = dAmatdr(:,:,j) + brdr(:,:,j)*grddbj*qat(i)
      dAmatdr(:,:,i) = dAmatdr(:,:,i) + brdr(:,:,i)*grddbi*qat(j)

   enddo

   ! self-energy part
   do i = 1, nat
      bp = 1._wp/brad(i)
      qq = qat(i)*bp
      egb = egb + 0.5_wp*qat(i)*qq*keps
      grddbi = -keps*bp*bp*0.5_wp
      dAmatdr(:,:,i) = dAmatdr(:,:,i) + brdr(:,:,i)*grddbi*qat(i)
   enddo

   gborn = egb

end subroutine addBornDerivStill


pure subroutine addBornMatSaltStill(nat, ntpair, ppind, ddpair, kappa, brad, ionscr, &
      & Amat)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Number of all interacting pairs
   integer, intent(in) :: ntpair

   !> Index list
   integer, intent(in) :: ppind(:, :)

   !> Distances of all pairs
   real(wp), intent(in) :: ddpair(:, :)

   !> Debye screening length
   real(wp), intent(in) :: kappa

   !> Born radii
   real(wp), intent(in) :: brad(:)

   !> Ion screening
   real(wp), intent(in) :: ionscr(:)

   !> Interaction matrix
   real(wp), intent(inout) :: Amat(:, :)

   integer  :: i, j, nnj
   integer  :: kk
   real(wp), parameter :: a13=1.0_wp/3.0_wp
   real(wp), parameter :: a4=0.25_wp
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp) :: aa, r2, gg, arg, bp
   real(wp) :: dd, expd, fgb, fgb2, dfgb

   ! compute energy and Amat direct and radii derivatives
   do kk = 1, ntpair
      r2=ddpair(1, kk)
      r2=r2*r2

      i=ppind(1, kk)
      j=ppind(2, kk)

      aa=brad(i)*brad(j)
      dd=a4*r2/aa
      expd=exp(-dd)
      fgb2=r2+aa*expd
      fgb=sqrt(r2+aa*expd)
      dfgb=1._wp/fgb
      gg=ionscr(i)+ionscr(j)
      Amat(i, j)=(exp(-kappa*fgb)*gg-1.0_wp)*dfgb + Amat(i, j)
      Amat(j, i)=(exp(-kappa*fgb)*gg-1.0_wp)*dfgb + Amat(j, i)
   enddo

   ! self-energy part
   do i = 1, nat
      gg=ionscr(i)*2.0_wp
      Amat(i, i)= Amat(i, i) + (exp(-kappa*brad(i))*gg-1.0_wp)/brad(i)
   enddo

end subroutine addBornMatSaltStill


pure subroutine addBornMatStill(nat, ntpair, ppind, ddpair, keps, brad, Amat)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Number of all interacting pairs
   integer, intent(in) :: ntpair

   !> Index list
   integer, intent(in) :: ppind(:, :)

   !> Distances of all pairs
   real(wp), intent(in) :: ddpair(:, :)

   !> Dielectric screening
   real(wp), intent(in) :: keps

   !> Born radii
   real(wp), intent(in) :: brad(:)

   !> Interaction matrix
   real(wp), intent(inout) :: Amat(:, :)

   integer  :: i, j, nnj
   integer  :: kk
   real(wp), parameter :: a13=1.0_wp/3.0_wp
   real(wp), parameter :: a4=0.25_wp
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp) :: aa, r2, gg, arg, bp
   real(wp) :: dd, expd, fgb, fgb2, dfgb

   ! compute energy and Amat direct and radii derivatives
   do kk = 1, ntpair
      r2 = ddpair(1, kk)
      r2 = r2*r2

      i = ppind(1, kk)
      j = ppind(2, kk)

      aa = brad(i)*brad(j)
      dd = a4*r2/aa
      expd = exp(-dd)
      fgb2 = r2+aa*expd
      dfgb = 1.0_wp/sqrt(fgb2)
      Amat(i, j) = keps*dfgb + Amat(i, j)
      Amat(j, i) = keps*dfgb + Amat(j, i)
   enddo

   ! self-energy part
   do i = 1, nat
      bp = 1._wp/brad(i)
      Amat(i, i) = Amat(i, i) + keps*bp
   enddo

end subroutine addBornMatStill


subroutine addBornMatP16(nat, ntpair, ppind, ddpair, keps, brad, Amat)

   !> Number of atoms
   integer, intent(in) :: nat

   !> Number of all interacting pairs
   integer, intent(in) :: ntpair

   !> Index list
   integer, intent(in) :: ppind(:, :)

   !> Distances of all pairs
   real(wp), intent(in) :: ddpair(:, :)

   !> Dielectric screening
   real(wp), intent(in) :: keps

   !> Born radii
   real(wp), intent(in) :: brad(:)

   !> Interaction matrix
   real(wp), intent(inout) :: Amat(:, :)

   integer :: kk
   integer :: iat, jat
   real(wp) :: r1, ab, arg, eab, fgb, dfgb, bp

   ! omp parallel do default(none) shared(Amat, ntpair, ppind, ddpair, brad, keps) &
   ! omp private(kk, iat, jat, r1, ab, arg, fgb, dfgb)
   do kk = 1, ntpair
      r1 = ddpair(1, kk)

      iat = ppind(1, kk)
      jat = ppind(2, kk)

      ab = sqrt(brad(iat) * brad(jat))
      arg = ab / (ab + zetaP16o16*r1) ! ab / (1 + ζR/(16·ab))
      arg = arg * arg ! ab / (1 + ζR/(16·ab))²
      arg = arg * arg ! ab / (1 + ζR/(16·ab))⁴
      arg = arg * arg ! ab / (1 + ζR/(16·ab))⁸
      arg = arg * arg ! ab / (1 + ζR/(16·ab))¹⁶
      fgb = r1 + ab*arg
      dfgb = 1.0_wp / fgb

      Amat(iat, jat) = keps*dfgb + Amat(iat, jat)
      Amat(jat, iat) = keps*dfgb + Amat(jat, iat)
   enddo

   ! self-energy part
   do iat = 1, nat
      bp = 1.0_wp/brad(iat)
      Amat(iat, iat) = Amat(iat, iat) + keps*bp
   enddo

end subroutine addBornMatP16


end module xtb_solv_kernel
