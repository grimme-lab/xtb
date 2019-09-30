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

module fixparam
   use iso_fortran_env, only : wp => real64
   use tbdef_setvar
   implicit none
   private :: wp
   public

!! ------------------------------------------------------------------------
!  exact fixing and freezing for Hessian
!! ------------------------------------------------------------------------
   type(fix_setvar) :: fixset
   type(fix_setvar) :: shakeset
   type(fix_setvar) :: freezeset

!! ------------------------------------------------------------------------
!  RMSD based meta dynamic feature
!! ------------------------------------------------------------------------
   type(metadyn_setvar) :: metaset

contains

subroutine init_fix(nfix)
   implicit none
   integer, intent(in) :: nfix
   call clear_fix
   call fixset%allocate(nfix)
   call shakeset%allocate(nfix*(nfix+1)/2)
   call freezeset%allocate(nfix)
end subroutine init_fix

subroutine clear_fix
   call fixset%deallocate
   call shakeset%deallocate
   call freezeset%deallocate
end subroutine clear_fix

subroutine init_metadyn(nat,nstruc)
   implicit none
   integer, intent(in) :: nat
   integer, intent(in) :: nstruc
   call metaset%allocate(nat,nstruc)
   metaset%factor = metaset%global_factor
end subroutine init_metadyn

subroutine clear_metadyn
   call metaset%deallocate
end subroutine clear_metadyn

subroutine fix_info(iunit,n,at,xyz)
   use mctc_constants
   use mctc_econv
   implicit none
   integer, intent(in)  :: iunit
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: xyz(3,n)

   character(len=2),external :: asym

   integer  :: i,ii,j,k,l,m,mm
   real(wp) :: val

   if (fixset%n.gt.0 .or. freezeset%n.gt.0 .or. shakeset%n.gt.0) then
      call generic_header(iunit,"Fixed Atoms",49,10)
      write(iunit,'(a)')
   endif
   if (fixset%n.gt.0) then
      write(iunit,'(1x,"*",1x,i0,1x,a)') fixset%n, &
         "fixed atom positions, i.e. in gradient"
      write(iunit,'(a)')
      write(iunit,'(5x,"#",3x,"Z",3x,32x,"position/Å")')
      do i = 1, fixset%n
         ii = fixset%atoms(i)
         write(iunit,'(i6,1x,i3,1x,a2,3f14.7)') &
            ii,at(ii),asym(at(ii)),xyz(:,ii)*autoaa
      enddo
      write(iunit,'(a)')
   endif
   if (freezeset%n.gt.0) then
      write(iunit,'(1x,"*",1x,i0,1x,a)') freezeset%n, &
         "frozen atom positions, i.e. in hessian"
      write(iunit,'(a)')
      write(iunit,'(5x,"#",3x,"Z",3x,32x,"position/Å")')
      do i = 1, freezeset%n
         ii = freezeset%atoms(i)
         write(iunit,'(i6,1x,i3,1x,a2,3f14.7)') &
            ii,at(ii),asym(at(ii)),xyz(:,ii)*autoaa
      enddo
      write(iunit,'(a)')
   endif
   if (shakeset%n.gt.0) then
      write(iunit,'(1x,"*",1x,i0,1x,a)') shakeset%n/2, &
         "constrained distances for dynamics (SHAKE)"
      write(iunit,'(a)')
      write(iunit,'(2(5x,"#",3x,"Z",3x),8x,"value/Å")')
      do m = 1, shakeset%n, 2
         i = shakeset%atoms(m)
         j = shakeset%atoms(m+1)
         val = norm2(xyz(:,i)-xyz(:,j))
         write(iunit,'(2(i6,1x,i3,1x,a2),1x,f14.7)') &
            i,at(i),asym(at(i)), &
            j,at(j),asym(at(j)), &
            val*autoaa
      enddo
      write(iunit,'(a)')

   endif

end subroutine fix_info

end module fixparam
