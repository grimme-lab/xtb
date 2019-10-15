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

module mdoptim
   use iso_fortran_env, wp => real64
contains

subroutine mdopt(mol,wfx,calc,egap,et,maxiter,epot,grd,sigma)
   use tbdef_molecule
   use tbdef_calculator
   use tbdef_wavefunction
   use tbdef_data

   use setparam
   use splitparam

   implicit none
   type(tb_molecule), intent(inout) :: mol
   type(tb_wavefunction),intent(inout) :: wfx
   type(tb_calculator),intent(in) :: calc
   integer  :: icall,maxiter
   real(wp) :: epot,et,egap
   real(wp), intent(inout) :: grd(3,mol%n)
   real(wp), intent(inout) :: sigma(3,3)
   type(scc_results) :: res

   real(wp),allocatable :: xyznew(:,:,:), eread(:)
   real(wp) :: e,angtoau
   parameter (angtoau=1.0d0/0.529177260d0)
   integer i,j,iz1,iz2,nall,ol
   character(80) atmp
   logical ex, fail
   integer ich

   write(*,*)
   write(*,'(7x,''======================================='')')
   write(*,'(7x,''|              M D O P T              |'')')
   write(*,'(7x,''|    optimizations along trajectory   |'')')
   write(*,'(7x,''======================================='')')
   write(*,*)
   write(*,*)'skipping interval            :',skip_md
   write(*,*)'opt level                    :',optset%optlev

   atmp='xtb.trj'
   call cqpath_read_pathfile_parameter(atmp,iz1,iz2,nall)
   if(iz2.ne.mol%n) stop 'read error in mdopt'
   allocate(xyznew(3,mol%n,nall),eread(nall))
   call cqpath_read_pathfile(atmp,iz1,iz2,nall,xyznew,mol%at,eread)
   write(*,*)'total number of points on trj:',nall
   write(*,*)'total number of optimizations:',nall/skip_md

   call open_file(ich,'xtb_ensemble.xyz','w')

   j   =0
   do i=1,nall

      j = j +1

      if(mod(j,skip_md).eq.0)then
         mol%xyz=xyznew(:,:,i)*angtoau

         call geometry_optimization &
            &      (mol,wfx,calc, &
            &       egap,etemp,maxscciter,optset%maxoptcycle,epot,grd,sigma, &
            &       optset%optlev,.false.,.true.,fail)


         if(.not.fail)then
            call wrxyz(ich,mol%n,mol%at,mol%xyz,epot,0.0d0)    ! write ensemble file
         endif
      endif

   enddo

!     include the input (may be lower than those on trj)
!     call ancopt(n,at,xyz,q,qsh,nshell,z,cn,nel,nopen,nbf,nao,P,wb, &
!    &        kspd,gscal,kcn,xbrad,xbdamp,alphaj,d3a1,d3a2,d3s8, &
!    &        d3atm,egap,etemp,maxscciter,maxoptcycle,epot,grd,ol,   &
!    &        .false.,fail)
!     if(.not.fail)then
!        call wrxyz(ich,n,at,xyz,epot,0.0d0)    ! write ensemble file
!     endif

   call close_file(ich)

end subroutine mdopt

subroutine wrxyz(iu,n,at,xyz,e,e2)
   implicit none
   integer n,at(n),iu
   real*8 xyz(3,n),e,e2
   integer i,j
   character*2 asym
   real*8, parameter ::autoang=0.52917726d0

   write(iu,*) n
   write(iu,'('' SCF done '',2F16.8)') e,e2
   do j=1,n
      write(iu,'(a2,3F24.10)')asym(at(j)),xyz(1:3,j)*autoang
   enddo

end subroutine wrxyz


end module mdoptim
