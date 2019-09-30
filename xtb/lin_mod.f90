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

!---------------------------------------------------------------------
module lin_mod
!---------------------------------------------------------------------


contains

  !***********************************************************************
  !* address in packed array
  !***********************************************************************

  pure elemental integer function lin(i1,i2)
    integer,intent(in) :: i1,i2
    integer :: idum1,idum2
    idum1=max(i1,i2)
    idum2=min(i1,i2)
    lin=idum2+idum1*(idum1-1)/2        
    return
  end function lin

  pure elemental integer function lina(i1,i2)
    integer,intent(in) :: i1,i2
    integer :: idum1,idum2
    idum1=max(i1,i2)
    idum2=min(i1,i2)
    lina=idum2+idum1*(idum1-1)/2        
    return
  end function lina


end module lin_mod
