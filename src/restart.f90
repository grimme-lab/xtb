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
   implicit none

   public :: readRestart, writeRestart
   public :: read_restart_gff, write_restart_gff


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
      end if
   end if

end subroutine readRestart


subroutine read_restart_gff(fname,n,p_ext_gfnff,success,verbose)
   use iso_fortran_env, wp => real64, istdout => output_unit
   use xtb_gfnff_param
   use xtb_gfnff_fraghess, only : nsystem,ispinsyst,nspinsyst
   implicit none
   character(len=*),intent(in) :: fname
   integer,intent(in)  :: n
   integer,intent(in)  :: p_ext_gfnff
   logical,intent(out) :: success
   logical,intent(in)  :: verbose

   integer(int64) :: iver8,nat8

   integer :: ich ! file handle
   integer :: err
   logical :: exist

   success = .false.
   call open_binary(ich,fname,'r')
   if (ich.ne.-1) then
      !read the first byte, which identify the calculation specs
      read(ich,iostat=err) iver8,nat8
      if(err.eq.0) then
         if (iver8.ne.int(p_ext_gfnff,int64).and.verbose) &
            &  call raise('W','Version number missmatch in restart file.',1)
         if (nat8.ne.n.and.verbose) then
            call raise('W','Atom number missmatch in restart file.',1)
            success=.false.
            return
         else if (iver8.eq.int(p_ext_gfnff)) then
            success = .true.
            read(ich) ffTopo%nbond,ffTopo%nangl,ffTopo%ntors,ffTopo%nhb1,ffTopo%nhb2,ffTopo%nxb,ffTopo%nathbH,ffTopo%nathbAB,  &
                    & ffTopo%natxbAB,ffTopo%nbatm,ffTopo%nfrag,nsystem,ffTopo%maxsystem        
            read(ich) ffTopo%nbond_blist,ffTopo%nbond_vbond,ffTopo%nangl_alloc,ffTopo%ntors_alloc,ffTopo%bond_hb_nr,ffTopo%b_max
            call gfnff_param_alloc(n)
            if (.not.allocated(ispinsyst)) allocate( ispinsyst(n,ffTopo%maxsystem), source = 0 )
            if (.not.allocated(nspinsyst)) allocate( nspinsyst(ffTopo%maxsystem), source = 0 )
            read(ich) ffTopo%nb,ffTopo%bpair,ffTopo%blist,ffTopo%alist, &
               & ffTopo%tlist,ffTopo%b3list,ffTopo%hblist1,ffTopo%hblist2, &
               & ffTopo%hblist3,ffTopo%fraglist,ffTopo%hbatHl,ffTopo%hbatABl, &
               & ffTopo%xbatABl,ispinsyst,nspinsyst,ffTopo%bond_hb_AH, &
               & ffTopo%bond_hb_B,ffTopo%bond_hb_Bn,ffTopo%nr_hb
            read(ich) ffTopo%vbond,ffTopo%vangl,ffTopo%vtors,ffTopo%chieeq, &
               & ffTopo%gameeq,ffTopo%alpeeq,ffTopo%alphanb,ffTopo%qa, &
               & ffTopo%q,ffTopo%xyze0,ffTopo%zetac6,&
               & ffTopo%qfrag,ffTopo%hbbas
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

end subroutine read_restart_gff


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


subroutine write_restart_gff(fname,nat,p_ext_gfnff)
   use iso_fortran_env, wp => real64, istdout => output_unit
   use xtb_gfnff_param
   use xtb_gfnff_fraghess, only : nsystem,ispinsyst,nspinsyst
   implicit none
   character(len=*),intent(in) :: fname
   integer,intent(in)  :: nat
   integer,intent(in)  :: p_ext_gfnff
   integer :: ich ! file handle

   call open_binary(ich,fname,'w')
   !Dimensions
   write(ich) int(p_ext_gfnff,int64),int(nat,int64)
   write(ich) ffTopo%nbond,ffTopo%nangl,ffTopo%ntors,   &
            & ffTopo%nhb1,ffTopo%nhb2,ffTopo%nxb,ffTopo%nathbH,ffTopo%nathbAB,ffTopo%natxbAB,ffTopo%nbatm,ffTopo%nfrag,nsystem,  &
            & ffTopo%maxsystem
   write(ich) ffTopo%nbond_blist,ffTopo%nbond_vbond,ffTopo%nangl_alloc,ffTopo%ntors_alloc,ffTopo%bond_hb_nr,ffTopo%b_max
   !Arrays Integers
   write(ich) ffTopo%nb,ffTopo%bpair,ffTopo%blist,ffTopo%alist,ffTopo%tlist,ffTopo%b3list,ffTopo%hblist1,ffTopo%hblist2,ffTopo%hblist3, &
      & ffTopo%fraglist,ffTopo%hbatHl,ffTopo%hbatABl,ffTopo%xbatABl,ispinsyst,nspinsyst,             &
      & ffTopo%bond_hb_AH,ffTopo%bond_hb_B,ffTopo%bond_hb_Bn,ffTopo%nr_hb
   !Arrays Reals
   write(ich) ffTopo%vbond,ffTopo%vangl,ffTopo%vtors,ffTopo%chieeq,ffTopo%gameeq,ffTopo%alpeeq,ffTopo%alphanb,ffTopo%qa,ffTopo%q,       &
      & ffTopo%xyze0,ffTopo%zetac6,ffTopo%qfrag,ffTopo%hbbas
   call close_file(ich)
end subroutine write_restart_gff


end module xtb_restart
