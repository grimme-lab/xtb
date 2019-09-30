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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine getname1(i,atmp)
      use setparam, only : get_namespace
      integer,intent(in) :: i
      character(len=*),intent(out) :: atmp

      write(atmp,'(''scoord.'',i0)')i          
      !atmp = get_namespace(trim(atmp))

      end subroutine getname1

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!      subroutine getname3(i,atmp)
!      integer,intent(in) :: i
!      character(len=*),intent(out) :: atmp
!
!      write(atmp,'(''tmpcoord_siman.'',i0)')i          
!      atmp = get_namespace(trim(atmp))
!
!      end subroutine getname3
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!      subroutine getname2(i,j,a)
!      integer,intent(in) :: i,j
!      character(len=*),intent(out) :: a
!      write(a,'(''rmsd scoord.'',i0,'' scoord.'',i0,'' &>/dev/null'')')
!     . i,j
!      end subroutine getname2
!
