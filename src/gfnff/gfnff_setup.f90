subroutine gfnff_setup(verbose,restart,mol,p_ext_gfnff)
  use iso_fortran_env
  use re_start
  use tbdef_molecule
  use gff_param
  use setparam, only : ichrg
  implicit none
! Dummy
  !integer,intent(in) :: ich
  integer,intent(in) :: p_ext_gfnff
  logical,intent(in) :: restart
  logical,intent(in) :: verbose
  type(tb_molecule)  :: mol
! Stack
  logical            :: ex
  logical            :: success

  call gfnff_input(mol)
  call gfnff_set_param(mol%n)
  if (restart) then
     inquire(file='gfnff_topo', exist=ex)
     if (ex) then
       call read_restart('gfnff_topo',mol%n,p_ext_gfnff,success,.true.)
       if (success) write(*,'(10x,"GFN-FF topology read from file successfully!")')
       !hbrefgeo is usually set within gfnff_ini2/gfnff_hbset0 equal to initial xyz
       hbrefgeo=mol%xyz
       if (.not.success) then
          write(*,'(10x,"GFN-FF topology read in did not work!")')    
          write(*,'(10x,"Generating new topology file!")')
          call gfnff_ini(verbose,ini,mol%n,ichrg,mol%at,mol%xyz)
          call write_restart('gfnff_topo',mol%n,p_ext_gfnff)
       end if
     else
       call gfnff_ini(verbose,ini,mol%n,ichrg,mol%at,mol%xyz)
       if (.not.mol%struc%two_dimensional) then
          call write_restart('gfnff_topo',mol%n,p_ext_gfnff)
       end if   
     end if
  else if (.not.restart) then
     call gfnff_ini(verbose,ini,mol%n,ichrg,mol%at,mol%xyz)
     call write_restart('gfnff_topo',mol%n,p_ext_gfnff)
  end if

end subroutine gfnff_setup  

subroutine gfnff_input(mol)
  use iso_fortran_env, only : wp => real64
  use tbdef_molecule
  use gff_param
  use setparam, only : ichrg
  implicit none
! Dummy  
  type(tb_molecule),intent(in) :: mol
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
  integer           :: bond_ij(2)
  real(wp)          :: r
  real(wp)          :: dum1
  real(wp)          :: floats(10)
  logical           :: ex
  character(len=80) :: atmp
  character(len=80) :: s(10)

  if (.not.allocated(nb))       allocate( nb(20,mol%n), source = 0 )
  if (.not.allocated(qfrag))    allocate( qfrag(mol%n), source = 0.0d0 )
  if (.not.allocated(fraglist)) allocate( fraglist(mol%n), source = 0 )
  if (.not.allocated(q))        allocate( q(mol%n), source = 0.0d0 )

  if (allocated(mol%pdb)) then
    read_file_type = 2
    ini = .true.
  else if (allocated(mol%sdf)) then
    read_file_type = 1
    ini = .false.
  else
    read_file_type = 0
    ini = .true.
  end if

  select case(read_file_type)
  !--------------------------------------------------------------------
  ! PDB case
    case(2)
        ifrag=0
        associate(rn => mol%pdb%residue_number, qatom => mol%pdb%charge)
          do iresidue = minval(rn),maxval(rn)
            if (any(iresidue .eq. rn)) then
              ifrag=ifrag+1
              where(iresidue .eq. rn) fraglist = ifrag
            end if
          end do
          nfrag = maxval(fraglist)
          do iatom=1,mol%n
             qfrag(fraglist(iatom)) = qfrag(fraglist(iatom)) + dble(qatom(iatom))
          end do
          qpdb = qatom
        end associate
        ichrg=idint(sum(qfrag(1:nfrag)))
        write(*,'(10x,"charge from pdb residues: ",i0)') ichrg
  !--------------------------------------------------------------------
  ! SDF case
    case(1)
       nb=0
       nfrag=0
       do ibond = 1, len(mol%bonds)
          call mol%bonds%get_item(ibond,bond_ij)
          i = bond_ij(1)
          j = bond_ij(2)
          ni=nb(20,i)
          ex=.false.
          do k=1,ni
             if(nb(k,i).eq.j) then 
                ex=.true.
                exit
             endif
          enddo
          if(.not.ex)then
             nb(20,i)=nb(20,i)+1
             nb(nb(20,i),i)=j
             nb(20,j)=nb(20,j)+1
             nb(nb(20,j),j)=i
          endif
       end do
       associate(xyz => mol%xyz)
         do i=1,mol%n
            if(nb(20,i).eq.0)then
              dum1=1.d+42
              do j=1,i
                 r=sqrt((xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2)
                 if(r.lt.dum1.and.r.gt.0.001)then
                    dum1=r
                    k=j
                 endif
              enddo
              nb(20,i)=1
              nb(1,i)=k
           endif
         end do
       end associate
  !--------------------------------------------------------------------
  ! General case: input = xyz or coord
    case(0)
    !ex=.false.  
    !inquire(file='.CHRG',exist=ex)
       ! if (ichrg.ne.0) then
       !    qfrag(1) = ichrg
       !    qfrag(2:mol%n) = 9999
       !if(ex)then
       call open_file(ich,'.CHRG','r')
       if (ich.ne.-1) then
           !open(unit=1,file='.CHRG')
           read(ich,'(a)')atmp
           call close_file(ich)
           call readline(atmp,floats,s,ns,nf)
           qfrag(1:nf)=floats(1:nf)
           ichrg=int(sum(qfrag(1:nf)))
           qfrag(nf+1:mol%n)=9999
        else
           qfrag=0
        end if
  !-------------------------------------------------------------------
  ! Default
    case default
        write(*,'(10x,"Input file format not recognized!")')
  end select

 !-------------------------------------------------------------------
  efield = 0

end subroutine gfnff_input
