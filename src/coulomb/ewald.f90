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

!> Implementation of Ewald summation specific tasks.
!
!  Part of this code originates from the DFTB+ codebase.
module xtb_coulomb_ewald
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_type_environment, only : TEnvironment
   implicit none
   private

   public :: ewaldMatPBC3D, ewaldDerivPBC3D, ewaldDerivPBC3D_alp
   public :: getOptimalAlpha, getMaxG, getMaxR


   !> Nr. of max. bisection steps
   integer, parameter :: nSearchIter = 30

   abstract interface
   !> Returns the max. value of a term in the reciprocal space part of the Ewald
   !  summation for a given vector length.
   pure function getGTermGen(gg, alpha, vol) result(gTerm)
      import :: wp

      !> Length of the reciprocal space vector
      real(wp), intent(in) :: gg

      !> Parameter of the Ewald summation
      real(wp), intent(in) :: alpha

      !> Volume of the real space unit cell
      real(wp), intent(in) :: vol

      !> Reciprocal term
      real(wp) :: gTerm

   end function getGTermGen
   end interface

contains


pure function ewaldMatPBC3D(vec, gTrans, qpc, volume, alpha, scale) result(Amat)

   !> Distance from i to WSC atom
   real(wp), intent(in) :: vec(:)

   !> Reciprocal lattice
   real(wp), intent(in) :: gTrans(:, :)

   !> Pseudo-quadrupole charge
   real(wp), intent(in) :: qpc

   !> Direct cell volume
   real(wp), intent(in) :: volume

   !> Convergence factor
   real(wp), intent(in) :: alpha

   !> Additional scaling factor
   real(wp), intent(in) :: scale

   !> Element of interaction matrix
   real(wp) :: Amat

   integer :: iG
   real(wp) :: rik2, rik(3), expterm

   Amat = 0.0_wp
   do iG = 1, size(gTrans, dim=2)
      rik(:) = gTrans(:, iG)
      rik2 = dot_product(rik, rik)
      expterm = exp(-rik2/(4.0_wp*alpha**2))/rik2
      Amat = Amat + cos(dot_product(rik, vec)) * expterm * (1.0_wp + 2*rik2*qpc**2)
   end do
   Amat = Amat * 4.0_wp*pi/volume * scale

end function ewaldMatPBC3D


pure subroutine ewaldDerivPBC3D_alp(vec, gTrans, qpc, volume, alpha, scale, &
      & dAmat, sigma)

   !> Distance from i to WSC atom
   real(wp),intent(in) :: vec(:)

   !> Reciprocal lattice
   real(wp),intent(in) :: gTrans(:, :)

   !> Pseudo-quadrupole charge
   real(wp), intent(in) :: qpc

   !> Direct cell volume
   real(wp),intent(in) :: volume

   !> Convergence factor
   real(wp),intent(in) :: alpha

   !> Additional scaling factor
   real(wp), intent(in) :: scale

   !> Derivative of interaction matrix
   real(wp),intent(out) :: dAmat(:)

   !> Strain of interaction matrix
   real(wp),intent(out) :: sigma(:, :)

   integer  :: iG
   real(wp) :: rik2, rik(3), dtmp, expterm, arg
   real(wp) :: fqpc, falp, dS(3, 3)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   dAmat = 0.0_wp
   sigma = 0.0_wp
   fqpc = 2*qpc**2
   falp = 2.0_wp/3.0_wp/(4.0_wp*alpha**2)
   do iG = 1, size(gTrans, dim=2)
      rik = gTrans(:, iG)
      rik2 = dot_product(rik,rik)
      expterm = exp(-rik2/(4.0_wp*alpha**2))/rik2
      arg = dot_product(rik,vec)
      dtmp = -sin(arg) * expterm
      dAmat = dAmat + rik*dtmp
      dS = spread(rik,1,3)*spread(rik,2,3)
      sigma = sigma + expterm * cos(arg) * ( &
         & - unity * (1.0_wp + rik2*falp + rik2*fqpc) &
         & + (2.0_wp/rik2 + 0.5_wp/alpha**2 + 0.5_wp*fqpc) * dS)
   end do
   dAmat = dAmat * 4.0_wp*pi/volume * scale
   sigma = sigma * 4.0_wp*pi/volume * scale

end subroutine ewaldDerivPBC3D_alp


pure subroutine ewaldDerivPBC3D(vec, gTrans, qpc, volume, alpha, scale, &
      & dAmat, sigma)

   !> Distance from i to WSC atom
   real(wp),intent(in) :: vec(:)

   !> Reciprocal lattice
   real(wp),intent(in) :: gTrans(:, :)

   !> Pseudo-quadrupole charge
   real(wp), intent(in) :: qpc

   !> Direct cell volume
   real(wp),intent(in) :: volume

   !> Convergence factor
   real(wp),intent(in) :: alpha

   !> Additional scaling factor
   real(wp), intent(in) :: scale

   !> Derivative of interaction matrix
   real(wp),intent(out) :: dAmat(:)

   !> Strain of interaction matrix
   real(wp),intent(out) :: sigma(:, :)

   integer  :: iG
   real(wp) :: rik2, rik(3), dtmp, expterm, arg
   real(wp) :: fqpc, dS(3, 3)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   dAmat = 0.0_wp
   sigma = 0.0_wp
   fqpc = 2*qpc**2
   do iG = 1, size(gTrans, dim=2)
      rik = gTrans(:, iG)
      rik2 = dot_product(rik,rik)
      expterm = exp(-rik2/(4.0_wp*alpha**2))/rik2
      arg = dot_product(rik,vec)
      dtmp = -sin(arg) * expterm
      dAmat = dAmat + rik*dtmp
      dS = spread(rik,1,3)*spread(rik,2,3)
      sigma = sigma + 0.5_wp * expterm * cos(arg) * ( &
         & - unity * (1.0_wp + rik2*fqpc) &
         & + (2.0_wp/rik2 + 0.5_wp/alpha**2 + 0.5_wp*fqpc) * dS)
   end do
   dAmat = dAmat * 4.0_wp*pi/volume * scale
   sigma = sigma * 4.0_wp*pi/volume * scale

end subroutine ewaldDerivPBC3D


!> Get optimal alpha-parameter for the Ewald summation by finding alpha, where
!  decline of real and reciprocal part of Ewald are equal.
subroutine getOptimalAlpha(env, latVec, recVec, volume, tolerance, alpha)

   !> Source of the generated error
   character(len=*), parameter :: source = 'coulomb_ewald_getOptimalAlpha'

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Lattice vectors
   real(wp), intent(in) :: latVec(:,:)

   !> Reciprocal vectors
   real(wp), intent(in) :: recVec(:,:)

   !> Volume of the unit cell
   real(wp), intent(in) :: volume

   !> Tolerance for difference in real and rec. part
   real(wp), intent(in) :: tolerance

   !> Optimal alpha
   real(wp) :: alpha

   real(wp) :: alphaLeft, alphaRight
   real(wp), parameter :: alphaInit = 1.0e-8_wp

   real(wp) :: minG, minR, diff
   integer :: iIter
   integer :: iError
   character(len=100) :: errorString

   minG = sqrt(minval(sum(recVec(:,:)**2, dim=1)))
   minR = sqrt(minval(sum(latVec(:,:)**2, dim=1)))

   iError = 0
   alpha = alphaInit
   diff = diffRecReal(alpha, getGTermPBC3D, minG, minR, volume)
   do while (diff < -tolerance .and. alpha <= huge(1.0_wp))
      alpha = 2.0_wp * alpha
      diff = diffRecReal(alpha, getGTermPBC3D, minG, minR, volume)
   end do
   if (alpha > huge(1.0_wp)) then
      iError = 1
   elseif (alpha == alphaInit) then
      iError = 2
   end if

   if (iError == 0) then
      alphaLeft = 0.5_wp * alpha
      do while (diff < tolerance .and. alpha <= huge(1.0_wp))
         alpha = 2.0_wp * alpha
         diff = diffRecReal(alpha, getGTermPBC3D, minG, minR, volume)
      end do
      if (alpha > huge(1.0_wp)) then
         iError = 3
      end if
   end if

   if (iError == 0) then
      alphaRight = alpha
      alpha = (alphaLeft + alphaRight) / 2.0
      iIter = 0
      diff = diffRecReal(alpha, getGTermPBC3D, minG, minR, volume)
      do while (abs(diff) > tolerance .and. iIter <= nSearchIter)
         if (diff < 0) then
            alphaLeft = alpha
         else
            alphaRight = alpha
         end if
         alpha = (alphaLeft + alphaRight) / 2.0
         diff = diffRecReal(alpha, getGTermPBC3D, minG, minR, volume)
         iIter = iIter + 1
      end do
      if (iIter > nSearchIter) then
         iError = 4
      end if
   end if

   if (iError /= 0) then
      !alpha = exp(-0.310104 * log(volume) + 0.786382) / 2.0
      call env%error('Failed to determine the optimal alpha for Ewald sum', source)
   end if

end subroutine getOptimalAlpha


!> Returns the longest reciprocal vector which gives a bigger contribution to
!  the Ewald sum than a certain tolerance.
subroutine getMaxG(env, alpha, volume, minValue, maxG)

   !> Source of the generated error
   character(len=*), parameter :: source = 'coulomb_ewald_getMaxG'

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Parameter of the ewald summation
   real(wp), intent(in) :: alpha

   !> Volume of the unit cell
   real(wp), intent(in) :: volume

   !> Tolerance value
   real(wp), intent(in) :: minValue

   !> Magnitude of reciprocal vector
   real(wp) :: maxG

   real(wp), parameter :: gInit = 1.0e-8_wp
   real(wp) :: xLeft, xRight, yLeft, yRight, yy
   integer :: iError, iIter
   character(len=100) :: errorString

   iError = 0
   maxG = gInit
   yy = getGTermPBC3D(maxG, alpha, volume)
   do while (yy > minValue .and. maxG <= huge(1.0_wp))
      maxG = 2.0_wp * maxG
      yy = getGTermPBC3D(maxG, alpha, volume)
   end do
   if (maxG > huge(1.0_wp)) then
      iError = 1
   elseif (maxG == gInit) then
      iError = 2
   end if

   if (iError == 0) then
      xLeft = 0.5_wp * maxG
      yLeft = getGTermPBC3D(xLeft, alpha, volume)
      xRight = maxG
      yRight = yy

      iIter = 1
      do while (yLeft - yRight > minValue .and. iIter <= nSearchIter)
         maxG = 0.5_wp * (xLeft + xRight)
         yy = getGTermPBC3D(maxG, alpha, volume)
         if (yy >= minValue) then
            xLeft = maxG
            yLeft = yy
         else
            xRight = maxG
            yRight = yy
         end if
         iIter = iIter + 1
      end do
      if (iIter > nSearchIter) then
         iError = 3
      end if
   end if

   if (iError /= 0) then
      call env%error('Failed to determine max. reciprocal lattice vector', source)
   end if

end subroutine getMaxG


!> Returns the longest real space vector which gives a bigger contribution to
!  the Ewald sum than a certain tolerance.
subroutine getMaxR(env, alpha, minValue, maxR)

   !> Source of the generated error
   character(len=*), parameter :: source = 'coulomb_ewald_getMaxR'

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Parameter of the ewald summation
   real(wp), intent(in) :: alpha

   !> Tolerance value
   real(wp), intent(in) :: minValue

   !> Magnitude of real space vector
   real(wp) :: maxR

   real(wp), parameter :: rInit = 1.0e-8_wp
   real(wp) :: xLeft, xRight, yLeft, yRight, yy
   integer :: iError, iIter
   character(len=100) :: errorString

   iError = 0
   maxR = rInit
   yy = getRTerm(maxR, alpha)
   do while (yy > minValue .and. maxR <= huge(1.0_wp))
      maxR = 2.0_wp * maxR
      yy = getRTerm(maxR, alpha)
   end do
   if (maxR > huge(1.0_wp)) then
      iError = 1
   elseif (maxR == rInit) then
      iError = 2
   end if

   if (iError == 0) then
      xLeft = 0.5_wp * maxR
      yLeft = getRTerm(xLeft, alpha)
      xRight = maxR
      yRight = yy

      iIter = 1
      do while (yLeft - yRight > minValue .and. iIter <= nSearchIter)
         maxR = 0.5_wp * (xLeft + xRight)
         yy = getRTerm(maxR, alpha)
         if (yy >= minValue) then
            xLeft = maxR
            yLeft = yy
         else
            xRight = maxR
            yRight = yy
         end if
         iIter = iIter + 1
      end do
      if (iIter > nSearchIter) then
         iError = 3
      end if
   end if

   if (iError /= 0) then
      call env%error('Failed to determine max. real lattice vector', source)
   end if

end subroutine getMaxR


!> Returns the difference in the decrease of the real and reciprocal parts of the
!  Ewald sum. In order to make the real space part shorter than the reciprocal
!  space part, the values are taken at different distances for the real and the
!  reciprocal space parts.
pure function diffRecReal(alpha, getGTerm, minG, minR, volume) result(diff)

   !> Parameter for the Ewald summation
   real(wp), intent(in) :: alpha

   !> Procedure pointer to reciprocal routine
   procedure(getGTermGen) :: getGTerm

   !> Length of the shortest reciprocal space vector in the sum
   real(wp), intent(in) :: minG

   !> Length of the shortest real space vector in the sum
   real(wp), intent(in) :: minR

   !> Volume of the real space unit cell
   real(wp), intent(in) :: volume

   !> Difference between changes in the two terms
   real(wp) :: diff

   diff = ((getGTerm(4.0_wp*minG, alpha, volume) &
      & - getGTerm(5.0_wp*minG, alpha, volume))) &
      & - (getRTerm(2.0_wp*minR, alpha) - getRTerm(3.0_wp*minR, alpha))

end function diffRecReal


!> Returns the max. value of a term in the reciprocal space part of the Ewald
!  summation for a given vector length.
pure function getGTermPBC3D(gg, alpha, vol) result(gTerm)

   !> Length of the reciprocal space vector
   real(wp), intent(in) :: gg

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol

   !> Reciprocal term
   real(wp) :: gTerm

   gTerm = 4.0_wp*pi*(exp(-0.25_wp*gg**2/(alpha**2))/(vol*gg**2))

end function getGTermPBC3D


!> Returns the max. value of a term in the real space part of the Ewald summation
!  for a given vector length.
pure function getRTerm(rr, alpha) result(rTerm)

   !> Length of the real space vector
   real(wp), intent(in) :: rr

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Real space term
   real(wp) :: rTerm

   rTerm = erfc(alpha*rr)/rr

end function getRTerm


end module xtb_coulomb_ewald
