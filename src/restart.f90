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

module xtb_restart
   use xtb_mctc_accuracy, only : wp, i8
   use xtb_mctc_io, only : stdout
   use xtb_type_environment, only : TEnvironment
   use xtb_type_wavefunction, only : TWavefunction
   use xtb_scc_core, only : qsh2qat
   implicit none

   public :: readRestart, writeRestart


contains


subroutine readRestart(env,wfx,fname,n,at,gfn_method,success,verbose)
   character(len=*), parameter :: source = 'restart_readRestart'
   type(TEnvironment), intent(inout) :: env
   type(TWavefunction),intent(inout) :: wfx
   character(len=*),intent(in) :: fname
   integer,intent(in)  :: n
   integer,intent(in)  :: at(n)
   integer,intent(in)  :: gfn_method
   logical,intent(out) :: success
   logical,intent(in)  :: verbose

   integer(i8) :: iver8,idum8,n8,nshell8,nel8,nopen8
   integer :: ich ! file handle
   integer :: err
   logical :: exist

   success = .false.
   call open_binary(ich,fname,'r')
   if (ich.ne.-1) then
      ! read the first 48 byte, which identify the calculation specs
      read(ich,iostat=err) iver8,idum8,n8,nshell8,nel8,nopen8
      if(err.eq.0) then
         if (iver8.ne.int(gfn_method,i8).and.verbose) &
            &  call env%warning('Version number missmatch in restart file.', source)
         if (nel8.ne.int(wfx%nel,i8).and.verbose) &
            &  call env%warning('Number of electron missmatch in restart file.', source)
         if (nopen8.ne.int(wfx%nopen,i8).and.verbose) &
            &  call env%warning('Multiplicity missmatch in restart file.', source)
         if ((n8.eq.int(wfx%n,i8)).and.(nshell8.eq.int(wfx%nshell,i8))) then
            success = .true.
            read(ich) wfx%qsh
            if (verbose) &
            write(stdout,'("q/qsh data taken from xtbrestart")')
            call qsh2qat(n,at,wfx%nshell,wfx%qsh,wfx%q)
            if ((gfn_method.gt.1).and.(iver8.gt.1)) then
!              read dipole and qpole CAMM
               read(ich) wfx%dipm
               read(ich) wfx%qp
               if (verbose) &
               write(stdout,'("CAMM data taken from xtbrestart")')
            endif
         else
            if (verbose) &
            call env%warning('Dimension missmatch in restart file.', source)
            success = .false.
         endif
      else
         if (verbose) &
         call env%warning("Dimension missmatch in restart file.", source)
         success = .false.
      endif
      call close_file(ich)
   endif

end subroutine readRestart


subroutine writeRestart(env,wfx,fname,gfn_method)
   type(TEnvironment), intent(inout) :: env
   type(TWavefunction),intent(inout) :: wfx
   character(len=*),intent(in) :: fname
   integer,intent(in)  :: gfn_method
   integer :: ich ! file handle

   call open_binary(ich,fname,'w')
   write(ich) int(gfn_method,i8),int(gfn_method,i8), &
              int(wfx%n,i8),int(wfx%nshell,i8), &
              int(wfx%nel,i8),int(wfx%nopen,i8)
   write(ich) wfx%qsh
   if (gfn_method.gt.1) then
      write(ich) wfx%dipm
      write(ich) wfx%qp
   endif
   call close_file(ich)

end subroutine writeRestart


end module xtb_restart
