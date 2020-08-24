!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
module xtb_kopt
  use xtb_mctc_accuracy, only : wp
  use xtb_mctc_convert, only : autoaa
  use xtb_type_environment
  use xtb_type_molecule
  use xtb_type_restart
  use xtb_type_calculator
  use xtb_type_data
  use xtb_type_setvar
  use xtb_restart
  use xtb_setmod
  use xtb_setparam
  use xtb_geoopt
  use xtb_lsrmsd
  implicit none
  private :: wp
  public

contains
  
subroutine get_kopt( &
         & metaset,env,restart,mol,chk,calc,egap,et,maxiter,maxcycle,optlev,&
         & etot,g,sigma,acc)
  implicit none
! Dummy -----------------------------------------------------------------------
  class(TCalculator),intent(inout)            :: calc
  type(metadyn_setvar),intent(inout)          :: metaset
  type(TEnvironment),intent(inout)            :: env
  type(TMolecule),intent(inout)               :: mol
  type(TRestart),intent(inout)                :: chk
  integer,intent(in)                          :: maxiter
  integer,intent(in)                          :: maxcycle
  integer,intent(in)                          :: optlev
  real(wp),intent(inout)                      :: etot
  real(wp),intent(in)                         :: et
  real(wp),intent(in)                         :: acc
  real(wp),intent(inout)                      :: egap
  real(wp),intent(inout)                      :: g(3,mol%n)
  real(wp),intent(inout)                      :: sigma(3,3)
  logical,intent(in)                          :: restart
! Stack -----------------------------------------------------------------------
  type(TMolecule)                             :: tmol
  integer                                     :: i
  integer                                     :: idum
  real(wp)                                    :: current_rmsd
  real(wp)                                    :: kopt
  real(wp)                                    :: kbias
  real(wp)                                    :: ax,cx
  logical                                     :: fail
!! ========================================================================
  call generic_header(env%unit,"Optimal kpush determination",49,10)
!------------------------------------------------------------------------------
! copy input coordinates
  tmol=mol

! set kopt boundaries
  ax         =0.0_wp
  cx         =0.005_wp*(dble(mol%n)/(target_rmsd+0.01_wp))

  kbias = ax
  metaset%factor = -kbias
  call get_rmsd( &
                 & calc,env,restart,mol,chk,egap,et,maxiter,maxcycle,optlev,&
                 & etot,g,sigma,current_rmsd)
  write(env%unit,'("target rmsd / Å        ",f9.6)') target_rmsd  
  write(env%unit,'("unbiased initial rmsd  ",f9.6)') current_rmsd  
  write(env%unit,*)
  write(env%unit,'("iter. min.        max.        rmsd       kpush")') 
  if ( current_rmsd.gt.target_rmsd ) then
     do i = 1,10
        mol=tmol
        kbias = 0.5_wp * (cx +ax)
        metaset%factor = -kbias
        call get_rmsd( &
                       & calc,env,restart,mol,chk,egap,et,maxiter,maxcycle,optlev,&
                       & etot,g,sigma,current_rmsd)
        write(*,'(i2,4f12.6)') i,&
        &ax,cx,current_rmsd,metaset%factor(metaset%nstruc)
        if ( current_rmsd.gt.target_rmsd ) then
           ax = kbias
        else
           cx = kbias
        end if
        if ( abs(current_rmsd-target_rmsd) .lt. 1.0d-3 ) exit
     end do
  end if
  kopt = kbias
  metaset%factor = -kopt
  write(env%unit,'("final kpush: ",f9.6)') -kopt   

end subroutine get_kopt

subroutine get_rmsd( &
         & calc,env,restart,mol,chk,egap,et,maxiter,maxcycle,optlev,&
         & etot,g,sigma,rmsdval)
  implicit none
! Dummy -----------------------------------------------------------------------
  class(TCalculator),intent(inout)            :: calc
  type(TEnvironment),intent(inout)            :: env
  type(TMolecule),intent(inout)               :: mol
  type(TRestart),intent(inout)                :: chk
  integer,intent(in)                          :: maxiter
  integer,intent(in)                          :: maxcycle
  integer,intent(in)                          :: optlev
  real(wp),intent(inout)                      :: etot
  real(wp),intent(in)                         :: et
  real(wp),intent(out)                        :: rmsdval
  real(wp),intent(inout)                      :: egap
  real(wp),intent(inout)                      :: g(3,mol%n)
  real(wp),intent(inout)                      :: sigma(3,3)
  logical,intent(in)                          :: restart
! Stack -----------------------------------------------------------------------
  type(TMolecule)                             :: tmol
  integer                                     :: idum
  real(wp)                                    :: U(3,3)
  real(wp)                                    :: x_center(3),y_center(3)
  real(wp),allocatable                        :: grad(:,:)
  logical                                     :: fail
!------------------------------------------------------------------------------
! save input geometry
  tmol=mol
!------------------------------------------------------------------------------
! geometry optimization
  call geometry_optimization &
      &     (env,mol,chk,calc,   &
      &      egap,etemp,maxiter,maxcycle,etot,g,sigma,optlev,.false.,.true.,fail)

  allocate( grad(3,mol%n), source = 0.0_wp )
  call rmsd(mol%n,mol%xyz,tmol%xyz,1,U,x_center,y_center,rmsdval,.true.,grad)
  rmsdval=rmsdval*autoaa

end subroutine get_rmsd

end module xtb_kopt
