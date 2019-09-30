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

subroutine bfgs(nat3,gnorm,grad,grado,dx,hess)
   use iso_fortran_env, wp => real64
   !-------------------------------------------------------------------
   ! Purpose:
   ! Performs BFGS update of Hessian matrix
   !
   ! Input:
   ! nat3	= dimension parameter as declared in the calling routine
   ! grad	= actual gradient
   ! grado = gradient one cycle before
   ! dx    = displ = displacement = coords(k) - coord(k-1) ; k=cycle
   ! hess	= hessian matrix and in Output updated hessian
   !--------------------------------------------------------------------
   implicit none

   ! Input:    
   integer, intent(in) :: nat3
   real(wp),intent(in) :: grad(nat3)
   real(wp),intent(in) :: grado(nat3)
   real(wp),intent(in) :: dx(nat3)
   real(wp),intent(in) :: gnorm
   ! Output:
   real(wp),intent(inout) :: hess(nat3*(nat3+1)/2)
   ! Local:
   integer  :: i,j,ij,ii
   real(wp),allocatable :: svec(:),tvec(:)
   real(wp) :: ddtd, dds, temp
   real(wp) :: ddot, thrs, scal, damp, dampO,dampD,thr
   real(wp) :: ooddtd, oodds, sdds, tddtd
   !---------------------------------------------------------------------  
   allocate( svec(nat3),tvec(nat3), source = 0.0_wp )

   ! damping of H update
!  call hdamp(gnorm,dampO,dampD)

   thrs=1.d-12
   ! calculate dg = grad(k+1) - grad(k)
   svec(1:nat3) = grad(1:nat3) - grado(1:nat3) 

   ! calculate tvec = h*dx
   call dspmv('u',nat3,1.0_wp,hess,dx,1,0.0_wp,tvec,1)

   ! calculate scalar dxdx and jtdx
   ddtd = ddot(nat3,tvec,1,dx,1)
   dds  = ddot(nat3,svec,1,dx,1)
   ooddtd = 1.0_wp / ddtd
   oodds  = 1.0_wp / dds

   if(dds > thrs .and. ddtd > thrs) then
   !$omp parallel default(none) &
   !$omp shared(nat3,oodds,ooddtd,svec,tvec) &
   !$omp private(i,j,ii,ij,sdds,tddtd,temp) &
   !$omp shared(hess)
   !$omp do
      do i=1,nat3
         ii = i*(i-1)/2
         sdds  = svec(i)*oodds
         tddtd = tvec(i)*ooddtd
         do j=1,i
            ij = ii + j
!           scal=dampd
!           if(i.ne.j)scal=dampo
            !temp= (svec(i)*svec(j))/dds - (tvec(i)*tvec(j))/ddtd
            temp = svec(j)*sdds - tvec(j)*tddtd
            hess(ij) = hess(ij) + temp 
         end do
      end do
   !$omp end do
   !$omp end parallel
!  else
!     write(*,'(a)') ' ******* Hesse update not performed ******* '
!     write(*,*    ) dds,ddtd,thrs                                   
   endif

   ! limit diagonal to (0.01 slightly better than 0.001)
   thr=1.d-2
   ij=0
   do i=1,nat3
      ij=ij+i
      if(abs(hess(ij)).lt.thr)hess(ij)=thr   
   enddo

end subroutine bfgs

subroutine powell(nat3,gnorm,grad,grado,dx,hess)
   use iso_fortran_env, wp => real64

   !-------------------------------------------------------------------
   ! Purpose:
   ! Performs Powell update of Hessian matrix
   !
   ! Input:
   ! nat3	= dimension parameter as declared in the calling routine=
   !         3*natoms
   ! grad	= actual gradient
   ! grado = gradient one cycle before
   ! dx    = displ = displacement = coords(k) - coord(k-1) ; k=cycle
   ! hess	= hessian matrix and in Output updated hessian
   !--------------------------------------------------------------------
   implicit none

   ! Input:    
   integer, intent(in) :: nat3
   real(wp),intent(in) :: grad(nat3),grado(nat3),dx(nat3),gnorm
   ! Output:
   real(wp),intent(inout) :: hess(nat3*(nat3+1)/2)
   ! Local:
   integer :: i,j,ij
   real(wp) :: dds,ddtd,temp
   real(wp), dimension(nat3) :: tvec
   real(wp) :: ddot, thrs, scal, damp, dampO,dampD
   !---------------------------------------------------------------------  

   ! damping of H update
!  call hdamp(gnorm,dampD,dampO)

   thrs=1.d-14

   call dspmv('u',nat3,1.0d0,hess,dx,1,0.0d0,tvec,1)

   tvec(1:nat3) = grad(1:nat3) - grado(1:nat3) - tvec(1:nat3)

   ! calculate scalar dxdx and jtdx
   dds  = ddot(nat3,dx,1,dx,1)

   if(dds > thrs) then
      ddtd = ddot(nat3,tvec,1,dx,1)/dds
      do i=1,nat3
         do j=1,i
            ij = i*(i-1)/2 + j
!           scal=dampD
!           if(i.ne.j)scal=dampo
            temp=tvec(i)*dx(j) + dx(i)*tvec(j) - dx(i)*ddtd*dx(j)
!           hess(ij) = hess(ij) + temp*scal/dds
            hess(ij) = hess(ij) + temp/dds
         end do
      end do
   else
!     write(*,'(a)') ' ******* Hesse update not performed ******* '
   endif

end subroutine powell

subroutine hdamp(gnorm,dampO,dampD)
   use iso_fortran_env, wp => real64
   implicit none
   real(wp) gnorm,dampO,dampD
   dampD = 1.0
   dampO = 1.0
   return
   ! damping of H update
   if(gnorm > 0.5)then
      dampD = 1.0
      dampO = 0.0
   endif
   if(gnorm < 0.5 .and. gnorm > 0.2)then
      dampD = 1.0 
      dampO = 0.2 
   endif
   if(gnorm < 0.2 .and. gnorm > 0.05) then
      dampD = 1.0
      dampO = 0.5
   endif
end subroutine hdamp
