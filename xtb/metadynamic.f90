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

subroutine metadynamic(metavar,nat,at,xyz,ebias,g)
   use iso_fortran_env, wp => real64
   use tbdef_setvar
   use ls_rmsd
   implicit none
   type(metadyn_setvar),intent(in) :: metavar
   integer, intent(in)    :: nat
   integer, intent(in)    :: at(nat)
   real(wp),intent(in)    :: xyz(3,nat)
   real(wp),intent(inout) :: ebias
   real(wp),intent(inout) :: g(3,nat)

   real(wp),allocatable   :: xyzref(:,:),grad(:,:),xyzdup(:,:)
   real(wp) :: U(3,3), x_center(3), y_center(3)
   real(wp) :: etmp,rmsdval,e
   integer  :: i,j,k,iref,iat

   if (metavar%nat == 0) then
      allocate( xyzref(3,nat), grad(3,nat),source = 0.0_wp )
      !$omp parallel default(none) &
      !$omp shared(metavar,nat,xyz) &
      !$omp private(grad,xyzref,U,x_center,y_center,rmsdval,e,etmp) &
      !$omp reduction(+:ebias,g)
      !$omp do schedule(dynamic)
      do iref = 1, metavar%nstruc
         grad = 0.0_wp
         xyzref = metavar%xyz(:,:,iref)
         call rmsd(nat,xyz,xyzref,1,U,x_center,y_center,rmsdval, &
                   .true.,grad)
         e = metavar%factor(iref) * exp(-metavar%width * rmsdval**2)
         ebias = ebias + e
         etmp = -2.0_wp * metavar%width * e * rmsdval
         g = g + etmp*grad
      enddo
      !$omp enddo
      !$omp end parallel
   else
      allocate( xyzref(3,metavar%nat), xyzdup(3,metavar%nat), grad(3,metavar%nat) )
      !$omp parallel default(none) &
      !$omp shared(metavar,nat,xyz) &
      !$omp private(grad,xyzref,xyzdup,U,x_center,y_center,rmsdval,e,etmp,i,iat) &
      !$omp reduction(+:ebias,g)
      !$omp do schedule(dynamic)
      do iref = 1, metavar%nstruc
         grad = 0.0_wp
         do i = 1, metavar%nat
            iat = metavar%atoms(i)
            xyzref(:,i) = metavar%xyz(:,iat,iref)
            xyzdup(:,i) = xyz(:,iat)
         enddo
         call rmsd(metavar%nat,xyzdup,xyzref,1,U,x_center,y_center,rmsdval, &
                   .true.,grad)
         e = metavar%factor(iref) * exp(-metavar%width * rmsdval**2)
         ebias = ebias + e
         etmp = -2.0_wp * metavar%width * e * rmsdval
         do i = 1, metavar%nat
            iat = metavar%atoms(i)
            g(:,iat) = g(:,iat) + etmp*grad(:,i)
         enddo
      enddo
      !$omp enddo
      !$omp end parallel
      deallocate( xyzdup )
   endif
   deallocate( xyzref, grad )

end subroutine metadynamic

subroutine load_metadynamic(metavar,nat,at,xyz)
   use iso_fortran_env, wp => real64, istdout => output_unit
   use fixparam
   use readin
   implicit none
   type(metadyn_setvar),intent(inout) :: metavar
   integer, intent(in)    :: nat
   integer, intent(in)    :: at(nat)
   real(wp),intent(in)    :: xyz(3,nat)

   integer  :: nstruc

   if (.not.allocated(metavar%fname)) return

   nstruc = metavar%maxsave
   call readlog(metavar%fname,nat,at,metavar%xyz,nstruc)
   metavar%nstruc = nstruc
   write(istdout,'(a,1x,i0,1x,a)') &
      "metadynamics with", nstruc, "initial structures loaded"

end subroutine load_metadynamic

