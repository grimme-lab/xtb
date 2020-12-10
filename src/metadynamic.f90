! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

module xtb_metadynamic
contains

subroutine metadynamic(metavar,nat,at,xyz,ebias,g)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_setvar
   use xtb_lsrmsd
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

   if(metavar%nstruc < 1 ) return

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
         e = metavar%factor(iref) * exp(-metavar%width(iref) * rmsdval**2)
         ebias = ebias + e
         etmp = -2.0_wp * metavar%width(iref) * e * rmsdval
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
         e = metavar%factor(iref) * exp(-metavar%width(iref) * rmsdval**2)
         ebias = ebias + e
         etmp = -2.0_wp * metavar%width(iref) * e * rmsdval
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
   use xtb_mctc_io, only : stdout
   use xtb_mctc_accuracy, only : wp
   use xtb_fixparam
   use xtb_readin
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
   write(stdout,'(a,1x,i0,1x,a)') &
      "metadynamics with", nstruc, "initial structures loaded"
end subroutine load_metadynamic

subroutine load_rmsdbias(metavar,nat,at,xyz)
   use xtb_mctc_io, only : stdout
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_mctc_systools, only : getline
   use xtb_fixparam
   use xtb_readin
   implicit none
   type(metadyn_setvar),intent(inout) :: metavar
   integer, intent(in)    :: nat
   integer, intent(in)    :: at(nat)
   real(wp),intent(in)    :: xyz(3,nat)

   integer :: ii, nn, mm, cc, stat, unit
   integer, parameter :: initial_size = 16
   real(wp) :: kk, aa, xx, yy, zz
   character(len=:), allocatable :: line
   character(len=4) :: sym
   real(wp), allocatable :: tmp_xyz(:, :, :), tmp_realloc3(:, :, :)
   real(wp), allocatable :: tmp_par(:, :), tmp_realloc2(:, :)

   if (.not.allocated(metavar%fname)) return

   write(stdout, '("#", *(1x, g0))') &
      "Reading bias information from", metavar%fname

   allocate(tmp_xyz(3, nat, initial_size), tmp_par(2, initial_size))

   open(unit, file=metavar%fname, iostat=stat)
   cc = 0
   rdxyz: do while(stat == 0)
      call getline(unit, line, stat)
      if (stat /= 0) exit rdxyz
      read(line, *, iostat=stat) nn
      if (stat /= 0) exit rdxyz
      if (nn /= nat) then
         stat = 1
         exit rdxyz
      end if

      call getline(unit, line, stat)
      read(line, *, iostat=stat) kk, aa
      if (stat /= 0) exit rdxyz

      nn = size(tmp_xyz, 3)
      if (cc >= nn) then
         mm = nn + nn/2 + 1
         nn = min(nn, mm)
         call move_alloc(tmp_xyz, tmp_realloc3)
         call move_alloc(tmp_par, tmp_realloc2)
         allocate(tmp_xyz(3, nat, mm), tmp_par(2, mm))
         tmp_xyz(:, :, :nn) = tmp_realloc3(:, :, :nn)
         tmp_par(:, :nn) = tmp_realloc2(:, :nn)
         deallocate(tmp_realloc3, tmp_realloc2)
      end if
      cc = cc + 1

      tmp_par(:, cc) = [kk, aa]

      do ii = 1, nat
         call getline(unit, line, stat)
         if (stat /= 0) exit rdxyz
         read(line, *, iostat=stat) sym, xx, yy, zz
         if (stat /= 0) exit rdxyz
         tmp_xyz(:, ii, cc) = [xx, yy, zz] * aatoau
      end do
   end do rdxyz
   if (is_iostat_end(stat)) stat = 0
   if (stat /= 0) return
   close(unit, iostat=stat)

   write(stdout, '("#", *(1x, g0))') &
      "Read bias potential for", cc, "structures"

   call metavar%allocate(nat, cc)
   metavar%nstruc = cc
   metavar%xyz(:, :, :) = tmp_xyz(:, :, :cc)
   metavar%factor(:) = tmp_par(1, :cc)
   metavar%width(:) = tmp_par(2, :cc)

end subroutine load_rmsdbias

subroutine set_metadynamic(metavar,nat,at,xyz)
   use xtb_mctc_io, only : stdout
   use xtb_mctc_accuracy, only : wp
   use xtb_fixparam
   use xtb_readin
   implicit none
   type(metadyn_setvar),intent(inout) :: metavar
   integer, intent(in)    :: nat
   integer, intent(in)    :: at(nat)
   real(wp),intent(in)    :: xyz(3,nat)

   integer  :: nstruc

   nstruc = metavar%maxsave
   metavar%nstruc = nstruc
   metavar%xyz(:,:,nstruc) = xyz
   metavar%width=1.0
   write(stdout,'(a,1x,i0,1x,a)') &
      "metadynamics with", nstruc, "initial structures loaded"
end subroutine set_metadynamic

end module xtb_metadynamic
