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
   use xtb_mctc_blas, only : mctc_syrk, mctc_nrm2, mctc_symm
   use xtb_mctc_math, only : eigvec3x3, crossProd
   use xtb_type_molecule, only : TMolecule
   implicit none
   private

   public :: projectHessian, trproj


contains


!> Project translations and rotations from the Hessian matrix
subroutine projectHessian(hessian, mol, removeTrans, removeRot)

   !> Mass weighted Hessian matrix
   real(wp), intent(inout) :: hessian(:,:)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Translation should be projected
   logical, intent(in) :: removeTrans

   !> Rotations should be projected, will be automatically disabled for PBC input
   logical, intent(in) :: removeRot

   real(wp), allocatable :: nullvec(:, :), projector(:, :), scratch(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp) :: center(3), vec(3), r2, inertia(3,3), moments(3), molMass, axes(3,3)
   integer :: ndim, nproj, ii, jj, iat, ic

   ndim = 3*mol%n

   nproj = 0
   if (removeTrans) then
      nproj = nproj + 3
   end if
   if (removeRot .and. mol%npbc == 0) then
      nproj = nproj + 3
   end if

   if (nproj == 0) then
      return
   end if

   allocate(nullvec(ndim, nproj))
   allocate(projector(ndim, ndim))
   projector(:, :) = 0.0_wp
   nullvec(:, :) = 0.0_wp

   ! symmetrize
   do ii = 1, ndim
      do jj = 1, ii - 1
         hessian(jj, ii) = 0.5_wp*(hessian(jj, ii) + hessian(ii, jj))
         hessian(ii, jj) = hessian(jj, ii)
      end do
   end do

   do jj = 1, ndim
      projector(jj, jj) = 1.0_wp
   end do

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
         inertia(:, :) = inertia + mol%atmass(iat) &
            & * (unity*r2 - spread(vec, 1, 3)*spread(vec, 2, 3))
      end do

      call eigvec3x3(inertia, moments, axes)

      ! axis to project with respect to
      do ii = 1, 3
         if (moments(ii) < epsilon(0.0_wp)) cycle
         do iat = 1, mol%n
            ic = (iat - 1) * 3
            vec(:) = mol%xyz(:,iat) - center
            vec(:) = crossProd(axes(:, ii), vec)
            nullvec(ic+1:ic+3, nproj-ii+1) = vec
         end do
      end do
   end if

   ! mass weight the displacement vectors
   do iat = 1, mol%n
      ic = (iat - 1) * 3
      nullvec(ic+1:ic+3, :) = nullvec(ic+1:ic+3, :) * sqrt(mol%atmass(iat))
   end do

   ! renormalise
   do ii = 1, nproj
      if (sum(nullvec(:,ii)**2) > epsilon(1.0_wp)) then
         nullvec(:,ii) = nullvec(:,ii) / mctc_nrm2(nullvec(:, ii))
      end if
   end do

   call mctc_syrk(nullvec, projector, alpha=-1.0_wp, beta=1.0_wp)

   deallocate(nullvec)
   allocate(scratch(ndim, ndim))

   call mctc_symm(projector, hessian, scratch, side='r')
   call mctc_symm(projector, scratch, hessian, side='l')

end subroutine projectHessian

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
   use xtb_fixparam
   use xtb_mctc_la, only : blckmgs,syprj

   implicit none

   ! Input
   logical, intent(in) :: ldebug
   integer, intent(in) :: natoms,nat3,nmode,ndim
   real(wp), dimension(3,natoms) :: xyzucm
   real(wp), dimension(nat3,ndim):: mode
   ! Ouput
   real(wp), dimension(nat3*(nat3+1)/2) :: hess

   ! Local
   integer :: i,ii,iii
   real(wp), allocatable :: fmat(:,:)
   integer :: nprj

   nprj=6
   if(nmode.gt.0) nprj=nprj+nmode
   if(nmode.lt.0) nprj=nprj+fixset%n*3
   allocate(fmat(nat3,nprj))
   fmat(:,:) = 0.0_wp

   if(nmode.ge.0) then
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
   endif

   if(nmode.gt.0) then  ! NMF
      do i=1,nmode
         fmat(1:nat3,6+i)=mode(1:nat3,i)
      enddo
   endif

   if(nmode.lt.0) then ! exact fixing
      do i=1,natoms
         !        rotational vectors
         fmat(3*(i-1)+1,1) =  0.0_wp
         fmat(3*(i-1)+2,1) = -xyzucm(3,i)
         fmat(3*(i-1)+3,1) =  xyzucm(2,i)
         fmat(3*(i-1)+1,2) =  xyzucm(3,i)
         fmat(3*(i-1)+2,2) =  0.0_wp
         fmat(3*(i-1)+3,2) = -xyzucm(1,i)
         fmat(3*(i-1)+1,3) = -xyzucm(2,i)
         fmat(3*(i-1)+2,3) =  xyzucm(1,i)
         fmat(3*(i-1)+3,3) =  0.0_wp
      enddo
      do i=1,fixset%n
         iii=fixset%atoms(i)
         do ii=1,3
            fmat(3*(iii-1)+ii,3+(i-1)*3+ii) = 1.0_wp
         end do
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
