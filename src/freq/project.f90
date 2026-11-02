! This file is part of xtb.
!
! Copyright (C) 2019-2021 Sebastian Ehlert
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

!> Implementation of Hessian projections
module xtb_freq_project
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas, only : mctc_syrk, mctc_symm
   use xtb_mctc_la, only : blckmgs, syprj
   use xtb_mctc_math, only : eigvec3x3, crossProd
   use xtb_axis, only : is_linear
   use xtb_type_molecule, only : TMolecule
   implicit none
   private

   public :: projectHessian, projectFixedAtoms, fixedAtomModes, trproj


contains


!> Project translations and rotations from the Hessian matrix
subroutine projectHessian(hessian, mol, removeTrans, removeRot, fixedAtoms, &
      & nModes)

   !> Mass weighted Hessian matrix
   real(wp), intent(inout) :: hessian(:,:)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Translation should be projected
   logical, intent(in) :: removeTrans

   !> Rotations should be projected, will be automatically disabled for PBC input
   logical, intent(in) :: removeRot

   !> Atoms whose three Cartesian coordinates are fixed. The structure can then
   !> only rotate about axes running through all fixed atoms, and it cannot
   !> translate, so removeTrans should be false.
   integer, intent(in), optional :: fixedAtoms(:)

   !> Number of modes removed by the projection
   integer, intent(out), optional :: nModes

   real(wp), allocatable :: nullvec(:, :), projector(:, :), scratch(:, :)
   logical, allocatable :: fixed(:)
   integer :: ndim, ii, jj, iat, ic

   call getNullSpace(nullvec, fixed, mol, removeTrans, removeRot, .true., &
      & fixedAtoms, nModes)

   if (size(nullvec, 2) == 0 .and. .not.any(fixed)) then
      return
   end if

   ndim = 3*mol%n

   ! symmetrize
   do ii = 1, ndim
      do jj = 1, ii - 1
         hessian(jj, ii) = 0.5_wp*(hessian(jj, ii) + hessian(ii, jj))
         hessian(ii, jj) = hessian(jj, ii)
      end do
   end do

   ! fixed coordinates are removed by deleting their rows and columns
   do iat = 1, mol%n
      if (.not.fixed(iat)) cycle
      ic = 3*(iat - 1)
      hessian(ic+1:ic+3, :) = 0.0_wp
      hessian(:, ic+1:ic+3) = 0.0_wp
   end do

   if (size(nullvec, 2) == 0) then
      return
   end if

   allocate(projector(ndim, ndim))
   projector(:, :) = 0.0_wp

   do jj = 1, ndim
      projector(jj, jj) = 1.0_wp
   end do

   call mctc_syrk(nullvec, projector, alpha=-1.0_wp, beta=1.0_wp)

   deallocate(nullvec)
   allocate(scratch(ndim, ndim))

   call mctc_symm(projector, hessian, scratch, side='r')
   call mctc_symm(projector, scratch, hessian, side='l')

end subroutine projectHessian


!> Project the fixed atoms and the rotations still allowed by them from a
!> packed Cartesian Hessian, as used by the geometry optimizers
subroutine projectFixedAtoms(hessp, mol, fixedAtoms)

   !> Packed Cartesian Hessian matrix (no mass weighting)
   real(wp), intent(inout) :: hessp(:)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Atoms whose three Cartesian coordinates are fixed
   integer, intent(in) :: fixedAtoms(:)

   real(wp), allocatable :: nullvec(:, :)
   logical, allocatable :: fixed(:)
   integer :: ndim, ii, jj, ij

   call getNullSpace(nullvec, fixed, mol, .false., .true., .false., fixedAtoms)

   ndim = 3*mol%n

   ! fixed coordinates are removed by deleting their rows and columns
   ij = 0
   do ii = 1, ndim
      do jj = 1, ii
         ij = ij + 1
         if (fixed((ii - 1)/3 + 1) .or. fixed((jj - 1)/3 + 1)) then
            hessp(ij) = 0.0_wp
         end if
      end do
   end do

   if (size(nullvec, 2) > 0) then
      call syprj(ndim, size(nullvec, 2), nullvec, ndim, hessp)
   end if

end subroutine projectFixedAtoms


!> Number of Cartesian degrees of freedom removed by fixing atoms, i.e. the
!> fixed coordinates and the rotations still allowed by the fixed atoms
function fixedAtomModes(mol, fixedAtoms) result(nModes)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Atoms whose three Cartesian coordinates are fixed
   integer, intent(in) :: fixedAtoms(:)

   !> Number of removed degrees of freedom
   integer :: nModes

   real(wp), allocatable :: nullvec(:, :)
   logical, allocatable :: fixed(:)

   call getNullSpace(nullvec, fixed, mol, .false., .true., .false., &
      & fixedAtoms, nModes)

end function fixedAtomModes


!> Orthonormal basis of the translations and rotations to be projected from a
!> Hessian, linearly dependent vectors are set to zero. Fixed coordinates are
!> only marked, since they are removed by deleting their rows and columns,
!> the basis vectors are orthogonal to them. At most six vectors are built,
!> independent of the number of fixed atoms.
subroutine getNullSpace(nullvec, fixed, mol, removeTrans, removeRot, &
      & massWeight, fixedAtoms, nModes)

   !> Basis vectors of the projected rigid-body motions
   real(wp), allocatable, intent(out) :: nullvec(:, :)

   !> Atoms whose Cartesian coordinates are fixed
   logical, allocatable, intent(out) :: fixed(:)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Translation should be projected
   logical, intent(in) :: removeTrans

   !> Rotations should be projected, will be automatically disabled for PBC input
   logical, intent(in) :: removeRot

   !> Displacement vectors are mass weighted, otherwise unit masses are used.
   !> Rotations without fixed atoms are always about the center of mass, which
   !> spans the same space as about the centroid if translations are projected.
   logical, intent(in) :: massWeight

   !> Atoms whose three Cartesian coordinates are fixed. The structure can then
   !> only rotate about axes running through all fixed atoms, and it cannot
   !> translate, so removeTrans should be false.
   integer, intent(in), optional :: fixedAtoms(:)

   !> Number of removed degrees of freedom, fixed coordinates and linearly
   !> independent basis vectors
   integer, intent(out), optional :: nModes

   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp) :: center(3), vec(3), r2, inertia(3,3), moments(3), molMass, axes(3,3)
   integer :: ndim, nproj, nrot, nfixed, rotOffset, linearAxis
   integer :: ii, iat, ic
   logical :: isLinear, freeAxis(3)

   if (present(nModes)) nModes = 0

   ndim = 3*mol%n

   isLinear = mol%linear

   ! duplicated entries are only counted once, assign element-wise since a
   ! vector subscript with repeated indices must not be defined
   allocate(fixed(mol%n), source=.false.)
   if (present(fixedAtoms)) then
      do ii = 1, size(fixedAtoms)
         fixed(fixedAtoms(ii)) = .true.
      end do
   end if
   nfixed = count(fixed)

   nproj = 0
   if (removeTrans) then
      nproj = nproj + 3
   end if
   nrot = 0
   if (removeRot .and. mol%npbc == 0) then
      nrot = merge(2, 3, isLinear)
      nproj = nproj + nrot
   end if

   allocate(nullvec(ndim, nproj))
   nullvec(:, :) = 0.0_wp

   if (nproj == 0) then
      if (present(nModes)) nModes = 3*nfixed
      return
   end if

   if (removeTrans) then
      do iat = 1, mol%n
         ic = (iat - 1) * 3
         do ii = 1, 3
            nullvec(ic+ii, ii) = 1.0_wp
         end do
      end do
   end if

   if (removeRot .and. mol%npbc == 0) then
      center(:) = 0.0_wp
      molMass = 0.0_wp
      do iat = 1, mol%n
         center(:) = center + mol%xyz(:, iat) * mol%atmass(iat)
         molMass = molMass + mol%atmass(iat)
      end do
      center(:) = center / molMass

      inertia(:, :) = 0.0_wp
      do iat = 1, mol%n
         vec(:) = mol%xyz(:, iat) - center
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         inertia(:, :) = inertia + mol%atmass(iat) * unity * r2
         inertia(:, 1) = inertia(:, 1) - mol%atmass(iat) * vec(1) * vec
         inertia(:, 2) = inertia(:, 2) - mol%atmass(iat) * vec(2) * vec
         inertia(:, 3) = inertia(:, 3) - mol%atmass(iat) * vec(3) * vec
      end do

      call eigvec3x3(inertia, moments, axes)

      ! axis to project with respect to
      linearAxis = 0
      if (isLinear) linearAxis = minloc(moments, dim=1)
      do ii = 1, 3
         freeAxis(ii) = ii /= linearAxis &
            & .and. (isLinear .or. moments(ii) >= epsilon(0.0_wp))
      end do

      if (nfixed > 0) then
         ! Rotations about the center of mass would move the fixed atoms. The
         ! remaining rotations are those about axes through all fixed atoms,
         ! i.e. axes through their centroid with vanishing moment of inertia
         ! of the fixed atoms (unit masses).
         center(:) = 0.0_wp
         do iat = 1, mol%n
            if (fixed(iat)) center(:) = center + mol%xyz(:, iat)
         end do
         center(:) = center / real(nfixed, wp)

         ! A single fixed atom can rotate about any axis, keep the principal
         ! axes of the molecule to treat linear molecules as above.
         if (nfixed > 1) then
            inertia(:, :) = 0.0_wp
            do iat = 1, mol%n
               if (.not.fixed(iat)) cycle
               vec(:) = mol%xyz(:, iat) - center
               r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
               inertia(:, :) = inertia + unity * r2
               inertia(:, 1) = inertia(:, 1) - vec(1) * vec
               inertia(:, 2) = inertia(:, 2) - vec(2) * vec
               inertia(:, 3) = inertia(:, 3) - vec(3) * vec
            end do

            call eigvec3x3(inertia, moments, axes)

            ! Only collinear fixed atoms leave a rotation, about their common
            ! axis. For a linear molecule this is the molecular axis, rotation
            ! about which does not displace any atom.
            freeAxis(:) = .false.
            if (.not.isLinear .and. is_linear(mol%xyz(:, fixedAtoms))) then
               freeAxis(minloc(moments, dim=1)) = .true.
            end if
         end if
      end if

      rotOffset = merge(3, 0, removeTrans)
      do ii = 1, 3
         if (.not.freeAxis(ii)) cycle
         rotOffset = rotOffset + 1
         do iat = 1, mol%n
            ic = (iat - 1) * 3
            vec(:) = mol%xyz(:,iat) - center
            vec(:) = crossProd(axes(:, ii), vec)
            nullvec(ic+1:ic+3, rotOffset) = vec
         end do
      end do
   end if

   ! mass weight the displacement vectors
   if (massWeight) then
      do iat = 1, mol%n
         ic = (iat - 1) * 3
         nullvec(ic+1:ic+3, :) = nullvec(ic+1:ic+3, :) * sqrt(mol%atmass(iat))
      end do
   end if

   ! orthogonalize against the fixed coordinates
   do iat = 1, mol%n
      if (.not.fixed(iat)) cycle
      ic = (iat - 1) * 3
      nullvec(ic+1:ic+3, :) = 0.0_wp
   end do

   ! Build an orthonormal basis for the span of all displacement vectors,
   ! linearly dependent vectors are set to zero.
   call blckmgs(ndim, nproj, ndim, nullvec)

   if (present(nModes)) then
      nModes = 3*nfixed + count(sum(nullvec**2, dim=1) > 0.5_wp)
   end if

end subroutine getNullSpace

subroutine trproj(natoms,nat3,xyz,hess,ldebug,nmode,mode,ndim)

   !----------------------------------------------------------------------
   !  subroutine trproj drives projection of hessian out of the
   !  space of translational and rotational motions:
   !  first get xyz c.m.; second get transl. and rot. projection matrix
   !
   !  get center of mass coordinates with unit mass
   !
   ! Input
   !   natoms  = number of atoms
   !   nat3    = 3*natoms
   !   xyz     = cartesian coordinates
   !       ldebug  = debug flag = .true. for debugging
   !
   ! Ouput fo gtrprojm.f
   !   xyzucm  = temporary c.m. coordinates
   !
   !       hess        = projected hessian out of space of transl. and rot.
   !                 motion
   !
   !----------------------------------------------------------------------

   implicit none

   ! Input
   logical, intent(in) :: ldebug
   integer, intent(in) :: natoms,nat3,nmode,ndim
   real(wp), dimension(3,natoms) :: xyz
   real(wp), dimension(nat3,ndim):: mode
   ! Ouput
   real(wp), dimension(nat3*(nat3+1)/2) :: hess
   ! Local
   integer :: i
   real(wp) :: xm,ym,zm
   real(wp), dimension(3,natoms) ::xyzucm

   xyzucm(:,:) = xyz(:,:)

   xm = 0.0_wp
   ym = 0.0_wp
   zm = 0.0_wp

   do i=1,natoms
      xm = xm + xyzucm(1,i)
      ym = ym + xyzucm(2,i)
      zm = zm + xyzucm(3,i)
   end do

   xm = xm/natoms
   ym = ym/natoms
   zm = zm/natoms


   do i=1,natoms
      xyzucm(1,i) = xyzucm(1,i) - xm
      xyzucm(2,i) = xyzucm(2,i) - ym
      xyzucm(3,i) = xyzucm(3,i) - zm
   end do

   ! get translational and rotational projection matrix

   call gtrprojm(natoms,nat3,xyzucm,hess,ldebug,nmode,mode,ndim)

end subroutine trproj

!----------------------------------------------------------------------
subroutine gtrprojm(natoms,nat3,xyzucm,hess,ldebug,nmode,mode,ndim)
   !----------------------------------------------------------------------
   ! calculating the translational-rotational projection matrix
   !
   ! Input
   !   natoms  = number of atoms
   !   nat3    = 3*natoms
   !   xyzucm  = coords c.m. from gxyzucm.f
   !           hess    = hessian
   !       ldebug      = debug flag = .true. for debugging
   !
   ! Ouput
   !   fmat    = F-matrix with translational and rotational vectors
   !       pmat        = projection matrix P = (1-FFt)
   !   hess    = projected hessian
   !----------------------------------------------------------------------
   implicit none

   ! Input
   logical, intent(in) :: ldebug
   integer, intent(in) :: natoms,nat3,nmode,ndim
   real(wp), dimension(3,natoms) :: xyzucm
   real(wp), dimension(nat3,ndim):: mode
   ! Ouput
   real(wp), dimension(nat3*(nat3+1)/2) :: hess

   ! Local
   integer :: i,ii
   real(wp), allocatable :: fmat(:,:)
   integer :: nprj

   nprj=6
   if(nmode.gt.0) nprj=nprj+nmode
   allocate(fmat(nat3,nprj))
   fmat(:,:) = 0.0_wp

   do i=1,natoms
      do ii=1,3
         !        translation vectors
         fmat(3*(i-1)+ii,ii) = 1.0_wp
      end do
      !        rotational vectors
      fmat(3*(i-1)+1,4) =  0.0_wp
      fmat(3*(i-1)+2,4) = -xyzucm(3,i)
      fmat(3*(i-1)+3,4) =  xyzucm(2,i)
      fmat(3*(i-1)+1,5) =  xyzucm(3,i)
      fmat(3*(i-1)+2,5) =  0.0_wp
      fmat(3*(i-1)+3,5) = -xyzucm(1,i)
      fmat(3*(i-1)+1,6) = -xyzucm(2,i)
      fmat(3*(i-1)+2,6) =  xyzucm(1,i)
      fmat(3*(i-1)+3,6) =  0.0_wp
   end do

   if(nmode.gt.0) then  ! NMF
      do i=1,nmode
         fmat(1:nat3,6+i)=mode(1:nat3,i)
      enddo
   endif

   if(ldebug) then
      write(*,'(a)')
      write(*,'(a)') ' Basis vectors before orthonormalization'
      write(*,'(3e22.14)') fmat
   end if

   ! do orthogonalization

   call  blckmgs(nat3,nprj,nat3,fmat)

   ! do projection

   call syprj(nat3,nprj,fmat,nat3,hess)

   deallocate(fmat)

end subroutine gtrprojm

end module xtb_freq_project
