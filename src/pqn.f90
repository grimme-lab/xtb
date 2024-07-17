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


! principal quantum number of valence shell
pure elemental integer function pqn(at)
  integer,intent(in) :: at
  if(at.le.2)then
     pqn=1
  elseif(at.le.10)then
     pqn=2
  elseif(at.le.18)then
     pqn=3
  elseif(at.le.36)then
     pqn=4
  elseif(at.le.54)then
     pqn=5
  else
     pqn=6
  endif
end function pqn

pure elemental integer function ncore(at)
  integer,intent(in) :: at
  if(at.le.2)then
     ncore=0
  elseif(at.le.10)then
     ncore=2
  elseif(at.le.18)then
     ncore=10
  elseif(at.le.29)then   !zn
     ncore=18
  elseif(at.le.36)then
     ncore=28
  elseif(at.le.47)then
     ncore=36
  elseif(at.le.54)then
     ncore=46
  elseif(at.le.71)then
     ncore=54
  elseif(at.le.79)then
     ncore=68
  elseif(at.le.86)then
     ncore=78
  elseif(at.le.103)then
     ncore=86
  endif
end function ncore
