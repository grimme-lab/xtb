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

!> TODO
module xtb_disp_coordinationnumber
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_param_covalentradd3, only : covalentRadD3
   use xtb_param_paulingen, only : paulingEN
   use xtb_type_molecule, only : TMolecule, len
   use xtb_type_neighbourlist, only : TNeighbourList
   implicit none
   private

   public :: cnType, getCoordinationNumber, cutCoordinationNumber


   !> Possible counting functions for calculating coordination numbers
   type :: TCNTypeEnum

      !> Counting function not specified
      integer :: invalid = 0

      !> Original DFT-D3 coordination number
      integer :: exp = 1

      !> Faster decaying error function CN, better for dense systems
      integer :: erf = 2

      !> Error function CN with covalency correction
      integer :: cov = 3

      !> Particular long-ranged version of the DFT-D3 coordination number
      integer :: gfn = 4

   end type TCNTypeEnum

   !> Enumerator for different coordination number types
   type(TCNTypeEnum), parameter :: cnType = TCNTypeEnum()


   abstract interface
   !> Abstract interface for the counting function (and its derivative)
   pure function countingFunction(k, r, r0)
      import :: wp

      !> Constant for counting function
      real(wp), intent(in) :: k

      !> Actual distance
      real(wp), intent(in) :: r

      !> Critical distance
      real(wp), intent(in) :: r0

      !> Value of the counting function in the range of [0,1]
      real(wp) :: countingFunction

   end function countingFunction
   end interface


   !> Parameter for electronegativity scaling
   real(wp),parameter :: k4=4.10451_wp

   !> Parameter for electronegativity scaling
   real(wp),parameter :: k5=19.08857_wp

   !> Parameter for electronegativity scaling
   real(wp),parameter :: k6=2*11.28174_wp**2


contains


!> Geometric fractional coordination number, supports both error function
!  and exponential counting functions.
subroutine getCoordinationNumber(mol, neighs, neighlist, cf, cn, dcndr, dcndL)

   !> Molecular structure information.
   type(TMolecule), intent(in) :: mol

   !> Number of interacting neighbours.
   integer, intent(in) :: neighs(:)

   !> Neighbourlist.
   type(TNeighbourList), intent(in) :: neighlist

   !> Coordination number type (by counting function).
   integer, intent(in) :: cf

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out) :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out) :: dcndL(:, :, :)

   real(wp), parameter :: kcn_exp = 16.0_wp
   real(wp), parameter :: kcn_erf = 7.5_wp
   real(wp), parameter :: kcn_gfn = 10.0_wp

   select case(cf)
   case(cnType%exp)
      call ncoord_impl(mol, neighs, neighlist, kcn_exp, expCount, dexpCount, &
         &             .false., covalentRadD3, paulingEN, cn, dcndr, dcndL)
   case(cnType%erf)
      call ncoord_impl(mol, neighs, neighlist, kcn_erf, erfCount, derfCount, &
         &             .false., covalentRadD3, paulingEN, cn, dcndr, dcndL)
   case(cnType%cov)
      call ncoord_impl(mol, neighs, neighlist, kcn_erf, erfCount, derfCount, &
         &             .true., covalentRadD3, paulingEN, cn, dcndr, dcndL)
   case(cnType%gfn)
      call ncoord_impl(mol, neighs, neighlist, kcn_gfn, gfnCount, dgfnCount, &
         &             .false., covalentRadD3, paulingEN, cn, dcndr, dcndL)
   end select

end subroutine getCoordinationNumber

!> Actual implementation of the coordination number, takes a generic counting
!  function to return the respective CN.
subroutine ncoord_impl(mol, neighs, neighlist, kcn, cfunc, dfunc, enscale, &
      &                rcov, en, cn, dcndr, dcndL)

   !> Molecular structure information.
   type(TMolecule), intent(in) :: mol

   !> Number of interacting neighbours.
   integer, intent(in) :: neighs(:)

   !> Neighbourlist.
   type(TNeighbourList), target, intent(in) :: neighlist

   !> Function implementing the counting function.
   procedure(countingFunction) :: cfunc

   !> Function implementing the derivative of counting function w.r.t. distance.
   procedure(countingFunction) :: dfunc

   !> Use a covalency criterium by Pauling EN's.
   logical, intent(in) :: enscale

   !> Steepness of counting function
   real(wp), intent(in) :: kcn

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Electronegativity
   real(wp), intent(in) :: en(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out) :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out) :: dcndL(:, :, :)

   integer :: iat, jat, ati, atj, ij, img
   real(wp) :: r2, r1, rc, rij(3), countf, countd(3), stress(3, 3), den

   cn = 0.0_wp
   dcndr = 0.0_wp
   dcndL = 0.0_wp

   !$omp parallel do default(none) private(den) shared(enscale, rcov, en)&
   !$omp reduction(+:cn, dcndr, dcndL) shared(mol, kcn, neighlist, neighs) &
   !$omp private(ij, img, jat, ati, atj, r2, rij, r1, rc, countf, countd, stress)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      do ij = 1, neighs(iat)
         img = neighlist%iNeigh(ij, iat)
         r2 = neighlist%dist2(ij, iat)
         rij = -neighlist%coords(:, iat) + neighlist%coords(:, img) ! FIXME
         jat = neighlist%image(img)
         atj = mol%at(jat)
         r1 = sqrt(r2)

         rc = rcov(ati) + rcov(atj)

         if (enscale) then
            den = k4*exp(-(abs(en(ati)-en(atj)) + k5)**2/k6)
         else
            den = 1.0_wp
         endif

         countf = den * cfunc(kcn, r1, rc)
         countd = den * dfunc(kcn, r1, rc) * rij/r1

         cn(iat) = cn(iat) + countf
         if (iat /= jat) then
            cn(jat) = cn(jat) + countf
         endif

         dcndr(:, iat, iat) = dcndr(:, iat, iat) + countd
         dcndr(:, jat, jat) = dcndr(:, jat, jat) - countd
         dcndr(:, iat, jat) = dcndr(:, iat, jat) + countd
         dcndr(:, jat, iat) = dcndr(:, jat, iat) - countd

         stress = spread(countd, 1, 3) * spread(rij, 2, 3)

         dcndL(:, :, iat) = dcndL(:, :, iat) + stress
         if (iat /= jat) then
            dcndL(:, :, jat) = dcndL(:, :, jat) + stress
         endif

      enddo
   enddo
   !$omp end parallel do

end subroutine ncoord_impl


!> Error function counting function for coordination number contributions.
pure function erfCount(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count

   count = 0.5_wp * (1.0_wp + erf(-k*(r-r0)/r0))

end function erfCount


!> Derivative of the counting function w.r.t. the distance.
pure function derfCount(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp), parameter :: sqrtpi = sqrt(pi)

   real(wp) :: count

   count = -k/sqrtpi/r0*exp(-k**2*(r-r0)**2/r0**2)

end function derfCount


!> Exponential counting function for coordination number contributions.
pure function expCount(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count

   count =1.0_wp/(1.0_wp+exp(-k*(r0/r-1.0_wp)))

end function expCount


!> Derivative of the counting function w.r.t. the distance.
pure function dexpCount(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count
   real(wp) :: expterm

   expterm = exp(-k*(r0/r-1._wp))

   count = (-k*r0*expterm)/(r**2*((expterm+1._wp)**2))

end function dexpCount


!> Exponential counting function for coordination number contributions.
pure function gfnCount(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count

   count = expCount(k, r, r0) * expCount(2*k, r, r0+2)

end function gfnCount


!> Derivative of the counting function w.r.t. the distance.
pure function dgfnCount(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count

   count = dexpCount(k, r, r0) * expCount(2*k, r, r0+2) &
      &  + expCount(k, r, r0) * dexpCount(2*k, r, r0+2)

end function dgfnCount


!> Cutoff function for large coordination numbers
pure subroutine cutCoordinationNumber(nAtom, cn, dcndr, dcndL, maxCN)

   !> number of atoms
   integer, intent(in) :: nAtom

   !> on input coordination number, on output modified CN
   real(wp), intent(inout) :: cn(:)

   !> on input derivative of CN w.r.t. cartesian coordinates,
   !> on output derivative of modified CN
   real(wp), intent(inout), optional :: dcndr(:, :, :)

   !> on input derivative of CN w.r.t. strain deformation,
   !> on output derivative of modified CN
   real(wp), intent(inout), optional :: dcndL(:, :, :)

   !> maximum CN (not strictly obeyed)
   real(wp), intent(in), optional :: maxCN

   real(wp) :: cnmax
   integer :: iAt

   if (present(maxCN)) then
      cnmax = maxCN
   else
      cnmax = 4.5_wp
   end if

   if (cnmax <= 0.0_wp) return

   if (present(dcndL)) then
      do iAt = 1, nAtom
         dcndL(:, :, iAt) = dcndL(:, :, iAt) * dCutCN(cn(iAt), cnmax)
      end do
   end if

   if (present(dcndr)) then
      do iAt = 1, nAtom
         dcndr(:, :, iAt) = dcndr(:, :, iAt) * dCutCN(cn(iAt), cnmax)
      end do
   end if

   do iAt = 1, nAtom
      cn(iAt) = cutCN(cn(iAt), cnmax)
   end do

end subroutine cutCoordinationNumber


!> Cutting function for the coordination number.
elemental function cutCN(cn, cut) result(cnp)

   !> Current coordination number.
   real(wp), intent(in) :: cn

   !> Cutoff for the CN, this is not the maximum value.
   real(wp), intent(in) :: cut

   !> Cuting function vlaue
   real(wp) :: cnp

   cnp = log(1.0_wp + exp(cut)) - log(1.0_wp + exp(cut - cn))

end function cutCN


!> Derivative of the cutting function w.r.t. coordination number
elemental function dCutCN(cn, cut) result(dcnpdcn)

   !> Current coordination number.
   real(wp), intent(in) :: cn

   !> Cutoff for the CN, this is not the maximum value.
   real(wp), intent(in) :: cut

   !> Derivative of the cutting function
   real(wp) :: dcnpdcn

   dcnpdcn = exp(cut)/(exp(cut) + exp(cn))

end function dCutCn


end module xtb_disp_coordinationnumber
