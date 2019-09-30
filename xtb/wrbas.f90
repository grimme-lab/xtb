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

subroutine write_tm_basis(iunit,nat,at)
   use iso_fortran_env, wp => real64
   use aoparam
   use ehtparam
   implicit none
   character(2),external :: esym
   integer,intent(in)  :: iunit
   integer,intent(in)  :: nat
   integer,intent(in)  :: at(nat)
   character(1) :: lnam(0:3)
   integer :: nn(94)
   integer :: iat,iatyp,ish,icao,ip,iprim,ishtyp

   lnam(0)='s'
   lnam(1)='p'
   lnam(2)='d'
   lnam(3)='f'

   write(iunit,'(a)') '$basis'
   nn=0
   do iat = 1, nat
      nn(at(iat)) = iat
   enddo
   write(iunit,'(a)') '*'
   do iatyp = 1, 94
      iat = nn(iatyp)
      if (iat.eq.0) cycle
      write(iunit,'(a,1x,a)') trim(esym(iatyp)),'tbbas'
      write(iunit,'(a)') '*'
      do ish = 1, ao_n(iatyp)
         ishtyp = ao_l(ish,iatyp)
         icao = caoshell(ish,iat)
         write(iunit,'(1x,i3,2x,a1,25x,a)') &
            nprim(icao+1),lnam(ishtyp)
         do ip = 1, nprim(icao+1)
            iprim = ip + primcount(icao+1)
            if (cont(iprim) < 0) then
               write(iunit,'(2x,g16.11,1x,g17.11))') alp(iprim), cont(iprim)
            else
               write(iunit,'(2(2x,g16.11))') alp(iprim), cont(iprim)
            endif
         enddo
      enddo
      write(iunit,'(a)') '*'
   enddo
   write(iunit,'(a)') '$end'
end subroutine write_tm_basis

