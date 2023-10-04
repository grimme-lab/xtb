! This file is part of xtb.
!
! Copyright (C) 2019-2020 Stefan Grimme
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

module xtb_gfnff_setup
  use xtb_gfnff_ini, only : gfnff_ini
  use xtb_gfnff_data, only : TGFFData
  use xtb_gfnff_topology, only : TGFFTopology
  use xtb_gfnff_generator, only : TGFFGenerator
  implicit none
  character(len=*), parameter :: source = 'gfnff_setup'
  private
  public :: gfnff_setup, gfnff_input

contains

subroutine gfnff_setup(env,verbose,restart,mol,gen,param,topo,accuracy,version)
  use xtb_restart
  use xtb_type_environment, only : TEnvironment
  use xtb_type_molecule, only : TMolecule
  use xtb_gfnff_param, only : ini, gfnff_set_param
  use xtb_setparam, only : set
  implicit none
! Dummy
  !integer,intent(in) :: ich
  type(TGFFTopology), intent(inout) :: topo
  type(TGFFGenerator), intent(inout) :: gen
  type(TGFFData), intent(inout) :: param
  integer,intent(in) :: version
  logical,intent(in) :: restart
  logical,intent(in) :: verbose
  real(wp),intent(in) :: accuracy
  type(TMolecule)  :: mol
  type(TEnvironment), intent(inout) :: env
! Stack
  logical            :: ex
  logical            :: success
  logical :: exitRun

  call gfnff_input(env, mol, topo)
  call env%check(exitRun)
  if (exitRun) then
     call env%error("Failed to prepare topology from geometry input", source)
     return
  end if

  call gfnff_set_param(mol%n, gen, param)
  param%dispscale = set%dispscale
  if (restart) then
     call read_restart_gff(env,'gfnff_topo',mol%n,version,success,.true.,topo)
     if (success) then
        write(env%unit,'(10x,"GFN-FF topology read from file successfully!")')
        return
     else
        call env%warning("Could not read topology file.", source)
        call env%check(exitRun)
        if (exitRun) return
     end if
  end if

  call gfnff_ini(env,verbose,ini,mol,gen,param,topo,accuracy)

  call env%check(exitRun)
  if (exitRun) then
     call env%error("Failed to generate topology", source)
     return
  end if

  if (.not.mol%info%two_dimensional) then
     call write_restart_gff(env,'gfnff_topo',mol%n,version,topo)
     call write_gfnff_adjacency('gfnff_adjacency',topo)
  end if

end subroutine gfnff_setup

subroutine gfnff_input(env, mol, topo)
  use xtb_mctc_accuracy, only : wp
  use xtb_type_environment, only : TEnvironment
  use xtb_type_molecule
  use xtb_mctc_filetypes, only : fileType
  use xtb_gfnff_param
  implicit none
  ! Dummy
  type(TMolecule),intent(in) :: mol
  type(TEnvironment), intent(inout) :: env
  type(TGFFTopology), intent(inout) :: topo
  ! Stack
  integer           :: i,j,k
  integer           :: ni
  integer           :: ns
  integer           :: nf
  integer           :: ich
  integer           :: iatom
  integer           :: iresidue
  integer           :: ifrag
  integer           :: ibond
  integer           :: bond_ij(3)
  real(wp)          :: r
  real(wp)          :: dum1
  real(wp)          :: floats(10)
  logical           :: ex
  character(len=80) :: atmp, atmp_0
  character(len=80) :: s(10)
  integer, allocatable :: rn(:)
  ! IO Error
  integer :: err

  if (.not.allocated(topo%nb))       allocate( topo%nb(20,mol%n), source = 0 )
  if (.not.allocated(topo%qfrag))    allocate( topo%qfrag(mol%n), source = 0.0d0 )
  if (.not.allocated(topo%fraglist)) allocate( topo%fraglist(mol%n), source = 0 )

  select case(mol%ftype)
  !--------------------------------------------------------------------
  ! PDB case
  case(fileType%pdb)
    ini = .true.
    ifrag=0
    allocate(rn(mol%n))
    rn(:) = mol%pdb%residue_number
    do iresidue = minval(rn),maxval(rn)
      if (any(iresidue .eq. rn)) then
        ifrag=ifrag+1
        where(iresidue .eq. rn) topo%fraglist = ifrag
      end if
    end do
    deallocate(rn)
    topo%nfrag = maxval(topo%fraglist)
    if (.not.allocated(topo%qpdb)) allocate(topo%qpdb(mol%n))
    do iatom=1,mol%n
      topo%qfrag(topo%fraglist(iatom)) = topo%qfrag(topo%fraglist(iatom)) &
        & + dble(mol%pdb(iatom)%charge)
      topo%qpdb(iatom) = mol%pdb(iatom)%charge
    end do
    if (abs(mol%chrg - sum(topo%qfrag(1:topo%nfrag))) < sqrt(epsilon(1.0_wp))) then
       write(env%unit,'(10x,"charge from pdb residues: ",i0)') &
          & nint(sum(topo%qfrag(1:topo%nfrag)))
    else
       ! ignore fragment charges if they are not consistent with the total charge
       call env%warning("Fragment charges from PDB file are not consistent with the total charge (is this a manual override?)", source)
       if (allocated(topo%qpdb)) deallocate(topo%qpdb)
       topo%qfrag(1) = mol%chrg
       topo%qfrag(2:topo%nfrag) = 0.0_wp
       topo%nfrag = 0
    end if
  !--------------------------------------------------------------------
  ! SDF case
  case(fileType%sdf,fileType%molfile)
    ini = .false.
    topo%nb=0
    topo%nfrag=0
    do ibond = 1, len(mol%bonds)
      call mol%bonds%get_item(ibond,bond_ij)
      i = bond_ij(1)
      j = bond_ij(2)
      ni=topo%nb(20,i)
      ex=.false.
      do k=1,ni
        if(topo%nb(k,i).eq.j) then
          ex=.true.
          exit
        endif
      enddo
      if(.not.ex)then
        topo%nb(20,i)=topo%nb(20,i)+1
        topo%nb(topo%nb(20,i),i)=j
        topo%nb(20,j)=topo%nb(20,j)+1
        topo%nb(topo%nb(20,j),j)=i
      endif
    end do
    do i=1,mol%n
      if(topo%nb(20,i).eq.0)then
        dum1=1.d+42
        k = 0
        do j=1,i
          r=sqrt(sum((mol%xyz(:,i)-mol%xyz(:,j))**2))
          if(r.lt.dum1.and.r.gt.0.001)then
            dum1=r
            k=j
          endif
        enddo
        if (k > 0) then
          topo%nb(20,i)=1
          topo%nb(1,i)=k
        end if
      endif
    end do
    ! initialize qfrag as in the default case
    topo%qfrag(1)=mol%chrg
    topo%qfrag(2:mol%n)=0
  !--------------------------------------------------------------------
  ! General case: input = xyz or coord
  case default
    if (mol%npbc > 0) then
      call env%error("Input file format not suitable for GFN-FF!")
      return
    end if
    ini = .true.
    call open_file(ich,'.CHRG','r')
    if (ich.ne.-1) then
      read(ich,'(a)') atmp_0 ! first line contains total charge
      read(ich,'(a)',iostat=err)atmp  ! second line contains fragment charges
      if (err .ne. 0) atmp=atmp_0
      call close_file(ich)
      call readline(atmp,floats,s,ns,nf)
      topo%qfrag(1:nf)=floats(1:nf)
      if (abs(mol%chrg - sum(topo%qfrag(1:nf))) < sqrt(epsilon(1.0_wp))) then
         write(env%unit,'(10x,"charge from .CHRG file: ",i0)') &
            & nint(sum(topo%qfrag(1:nf)))
      else
         ! ignore fragment charges if they are not consistent with the total charge
         topo%qfrag(1:nf) = 0.0_wp
      end if
    else
      topo%qfrag(1)=mol%chrg
      topo%qfrag(2:mol%n)=0
    end if
  end select

  !-------------------------------------------------------------------

end subroutine gfnff_input


subroutine write_gfnff_adjacency(fname, topo)
   implicit none
   character(len=*),intent(in) :: fname
   integer :: ifile ! file handle
   type(TGFFTopology) :: topo
   integer :: i, j

   call open_file(ifile,fname,'w') 
   ! looping over topology neighboring list
   if (ifile.ne.-1) then
      write(ifile, '(a)') '# indices of neighbouring atoms (max seven)'  
      do i = 1, size(topo%nb, 2)
        write(ifile, '(*(i0:, 1x))') (topo%nb(j, i), j = 1, topo%nb(size(topo%nb, 1), i))
      end do
   end if
   call close_file(ifile)

end subroutine write_gfnff_adjacency


end module xtb_gfnff_setup
