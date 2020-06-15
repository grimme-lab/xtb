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

module xtb_mdoptim
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_filetypes, only : fileType
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_calculator
   use xtb_type_restart
   use xtb_type_data
   use xtb_io_writer, only : writeMolecule
   use xtb_setparam
   use xtb_splitparam
   use xtb_geoopt
   use xtb_cqpath
   implicit none
   private

   public :: mdopt


contains


subroutine mdopt(env, mol, chk, calc, egap, et, maxiter, epot, grd, sigma)

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(inout) :: mol
   type(TRestart),intent(inout) :: chk
   class(TCalculator), intent(inout) :: calc
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
   if(iz2.ne.mol%n) call env%terminate('read error in mdopt')
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
            &      (env, mol,chk,calc, &
            &       egap,etemp,maxscciter,optset%maxoptcycle,epot,grd,sigma, &
            &       optset%optlev,.false.,.true.,fail)

         if(.not.fail)then
            call writeMolecule(mol, ich, fileType%xyz, energy=epot)
         endif
      endif

   enddo

   call close_file(ich)

end subroutine mdopt


end module xtb_mdoptim
