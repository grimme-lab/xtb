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
   use xtb_gfnff_neighbor
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
      call close_file(ich)
   end if

end subroutine readRestart


subroutine read_restart_gff(env,fname,n,version,success,verbose,topo,neigh)
   use xtb_gfnff_param
   use xtb_gfnff_topology, only : TGFFTopology
   implicit none
   character(len=*), parameter :: source = 'restart_read_restart_gff'
   type(TEnvironment), intent(inout) :: env
   type(TGFFTopology), intent(inout) :: topo
   type(TNeigh), intent(inout) :: neigh
   character(len=*),intent(in) :: fname
   integer,intent(in)  :: n
   integer,intent(in)  :: version
   logical,intent(out) :: success
   logical,intent(in)  :: verbose

   integer(i8) :: iver8,nat8

   integer :: ich ! file handle
   integer :: err
   logical :: exist

   success = .false.
   call open_binary(ich,fname,'r')
   if (ich.ne.-1) then
      !read the first byte, which identify the calculation specs
      read(ich,iostat=err) iver8,nat8
      if(err.eq.0) then
         if (iver8.ne.int(version,i8).and.verbose) &
            &  call env%warning('Version number missmatch in restart file.',source)
         if (nat8.ne.n.and.verbose) then
            call env%warning('Atom number missmatch in restart file.',source)
            success=.false.
            call close_file(ich)
            return
         else if (iver8.eq.int(version)) then
            success = .true.
            read(ich) topo%nangl,topo%ntors,topo%nstors,topo%nathbH,topo%nathbAB,  &
                    & topo%natxbAB,topo%nbatm,topo%nfrag,topo%nsystem,topo%maxsystem
            read(ich) topo%nbond_blist,topo%nbond_vbond,topo%nangl_alloc,topo%ntors_alloc,topo%bond_hb_nr,topo%b_max
            read(ich) neigh%numnb, neigh%numctr, neigh%nbond, neigh%iTrDim
            read(ich) topo%hb_mapNAB, topo%hb_mapNH
            topo%nbond = neigh%nbond
            call gfnff_param_alloc(topo,neigh, n)
            if (.not.allocated(topo%ispinsyst)) allocate( topo%ispinsyst(n,topo%maxsystem), source = 0 )
            if (.not.allocated(topo%nspinsyst)) allocate( topo%nspinsyst(topo%maxsystem), source = 0 )
            if (.not.allocated(neigh%nr_hb)) allocate(neigh%nr_hb(neigh%nbond))
            if (.not.allocated(topo%hb_mapABH)) allocate(topo%hb_mapABH(n))
            if (.not.allocated(topo%isABH)) allocate(topo%isABH(n))
            read(ich) neigh%blist
            read(ich) topo%alist
            read(ich) topo%tlist
            read(ich) topo%b3list
            read(ich) topo%sTorsl
            read(ich) topo%fraglist
            read(ich) topo%hbatHl
            read(ich) topo%hbatABl
            read(ich) topo%xbatABl
            read(ich) topo%ispinsyst
            read(ich) topo%nspinsyst
            read(ich) topo%bond_hb_AH
            read(ich) topo%bond_hb_B
            read(ich) topo%bond_hb_Bn
            read(ich) neigh%nr_hb
            read(ich) topo%hb_mapABH
            read(ich) topo%isABH
            read(ich) topo%vangl,topo%vtors,topo%chieeq, &
               & topo%gameeq,topo%alpeeq,topo%alphanb,topo%qa, &
               & topo%xyze0,topo%zetac6,&
               & topo%qfrag,topo%hbbas,topo%hbaci
            if (.not.allocated(neigh%vbond)) allocate (neigh%vbond(3,neigh%nbond))
            if (.not.allocated(neigh%nb)) allocate (neigh%nb(neigh%numnb,n,neigh%numctr))
            if (.not.allocated(neigh%blist)) allocate (neigh%blist(3,neigh%nbond))
            if (.not.allocated(neigh%vbond)) allocate(neigh%vbond(3,neigh%nbond))
            if (.not.allocated(neigh%bpair)) allocate(neigh%bpair(n,n,neigh%numctr))
            if (.not.allocated(neigh%iTrSum)) allocate(neigh%iTrSum(neigh%iTrDim*(neigh%iTrDim+1)/2))
            if (.not.allocated(neigh%iTrNeg)) allocate(neigh%iTrNeg(neigh%iTrDim))
            read(ich) neigh%vbond, neigh%nb, neigh%blist, neigh%nr_hb, neigh%vbond, &
                    & neigh%bpair, neigh%iTrSum, neigh%iTrNeg
         else
            if (verbose) &
               call env%warning("Dimension missmatch in restart file.",source)
            success = .false.
         endif
      else
         if (verbose) &
            call env%warning("Dimension missmatch in restart file.",source)
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


subroutine write_restart_gff(env,fname,nat,version,topo,neigh)
   use xtb_gfnff_param
   use xtb_gfnff_topology, only : TGFFTopology
   implicit none
   type(TEnvironment), intent(inout) :: env
   type(TGFFTopology), intent(in) :: topo
   type(TNeigh), intent(inout) :: neigh
   character(len=*),intent(in) :: fname
   integer,intent(in)  :: nat
   integer,intent(in)  :: version
   integer :: ich ! file handle

   call open_binary(ich,fname,'w')
   !Dimensions
   write(ich) int(version,i8),int(nat,i8)
   write(ich) topo%nangl,topo%ntors,topo%nstors, topo%nathbH,topo%nathbAB, &
            & topo%natxbAB,topo%nbatm,topo%nfrag,topo%nsystem, topo%maxsystem
   write(ich) topo%nbond_blist,topo%nbond_vbond,topo%nangl_alloc,topo%ntors_alloc,topo%bond_hb_nr,topo%b_max
   write(ich) neigh%numnb, neigh%numctr, neigh%nbond, neigh%iTrDim
   write(ich) topo%hb_mapNAB, topo%hb_mapNH
   !Arrays Integers
   write(ich) neigh%blist
   write(ich) topo%alist
   write(ich) topo%tlist
   write(ich) topo%b3list
   write(ich) topo%sTorsl
   write(ich) topo%fraglist
   write(ich) topo%hbatHl
   write(ich) topo%hbatABl
   write(ich) topo%xbatABl
   write(ich) topo%ispinsyst
   write(ich) topo%nspinsyst
   write(ich) topo%bond_hb_AH
   write(ich) topo%bond_hb_B
   write(ich) topo%bond_hb_Bn
   write(ich) neigh%nr_hb
   write(ich) topo%hb_mapABH
   write(ich) topo%isABH
   !Arrays Reals
   write(ich) topo%vangl,topo%vtors,topo%chieeq,&
           & topo%gameeq,topo%alpeeq,topo%alphanb,topo%qa, &
           & topo%xyze0,topo%zetac6, &
           & topo%qfrag,topo%hbbas,topo%hbaci
   write(ich) neigh%vbond, neigh%nb, neigh%blist, neigh%nr_hb, neigh%vbond, &
           & neigh%bpair, neigh%iTrSum, neigh%iTrNeg
   call close_file(ich)
end subroutine write_restart_gff


end module xtb_restart
