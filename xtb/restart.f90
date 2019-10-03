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

subroutine read_restart(wfx,fname,n,at,gfn_method,success,verbose)
   use iso_fortran_env, wp => real64, istdout => output_unit
   use tbdef_wavefunction
   use scc_core, only : qsh2qat
   implicit none
   type(tb_wavefunction),intent(inout) :: wfx
   character(len=*),intent(in) :: fname
   integer,intent(in)  :: n
   integer,intent(in)  :: at(n)
   integer,intent(in)  :: gfn_method
   logical,intent(out) :: success
   logical,intent(in)  :: verbose

   integer(int64) :: iver8,idum8,n8,nshell8,nel8,nopen8
   integer :: ich ! file handle
   integer :: err
   logical :: exist

   success = .false.
   call open_binary(ich,fname,'r')
   if (ich.ne.-1) then
      ! read the first 48 byte, which identify the calculation specs
      read(ich,iostat=err) iver8,idum8,n8,nshell8,nel8,nopen8
      if(err.eq.0) then
         if (iver8.ne.int(gfn_method,int64).and.verbose) &
            &  call raise('S','Version number missmatch in restart file.',1)
         if (nel8.ne.int(wfx%nel,int64).and.verbose) &
            &  call raise('S','Number of electron missmatch in restart file.',1)
         if (nopen8.ne.int(wfx%nopen,int64).and.verbose) &
            &  call raise('S','Multiplicity missmatch in restart file.',1)
         if ((n8.eq.int(wfx%n,int64)).and.(nshell8.eq.int(wfx%nshell,int64))) then
            success = .true.
            read(ich) wfx%qsh
            if (verbose) &
            write(istdout,'("q/qsh data taken from xtbrestart")')
            call qsh2qat(n,at,wfx%nshell,wfx%qsh,wfx%q)
            if ((gfn_method.gt.1).and.(iver8.gt.1)) then
!              read dipole and qpole CAMM
               read(ich) wfx%dipm
               read(ich) wfx%qp
               if (verbose) &
               write(istdout,'("CAMM data taken from xtbrestart")')
            endif
         else
            if (verbose) &
            call raise('S','Dimension missmatch in restart file.',1)
            success = .false.
         endif
      else
         if (verbose) &
         call raise('S',"Dimension missmatch in restart file.",1)
         success = .false.
      endif
      call close_file(ich)
   endif

end subroutine read_restart

subroutine write_restart(wfx,fname,gfn_method)
   use iso_fortran_env, wp => real64, istdout => output_unit
   use tbdef_wavefunction
   implicit none
   type(tb_wavefunction),intent(inout) :: wfx
   character(len=*),intent(in) :: fname
   integer,intent(in)  :: gfn_method
   integer :: ich ! file handle

   call open_binary(ich,fname,'w')
   write(ich) int(gfn_method,int64),int(gfn_method,int64), &
              int(wfx%n,int64),int(wfx%nshell,int64), &
              int(wfx%nel,int64),int(wfx%nopen,int64)
   write(ich) wfx%qsh
   if (gfn_method.gt.1) then
      write(ich) wfx%dipm
      write(ich) wfx%qp
   endif
   call close_file(ich)

end subroutine write_restart
