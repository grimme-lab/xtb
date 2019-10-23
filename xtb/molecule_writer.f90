submodule(tbdef_molecule) molecule_writer
   use tbdef_molecule
   implicit none

contains

module subroutine write_molecule_generic(self, unit, format, energy, gnorm)
   use tbmod_file_utils
   include 'xtb_version.fh'
   class(tb_molecule), intent(in) :: self
   integer, intent(in) :: unit
   integer, intent(in), optional :: format
   real(wp), intent(in), optional :: energy
   real(wp), intent(in), optional :: gnorm
   character(len=:), allocatable :: comment_line
   character(len=20) :: energy_line
   character(len=20) :: gnorm_line
   integer :: ftype
   if (present(format)) then
      ftype = format
   else
      ftype = self%ftype
   endif

   comment_line = ''
   if (present(energy)) then
      write(energy_line, '(f20.12)') energy
      comment_line = comment_line // " energy: " // trim(adjustl(energy_line))
   endif
   if (present(gnorm)) then
      write(gnorm_line, '(f20.12)') gnorm
      comment_line = comment_line // " gnorm: " // trim(adjustl(gnorm_line))
   endif
   comment_line = comment_line // " xtb: " // version

   select case(ftype)
   case(p_ftype%xyz)
      call write_xyz(self, unit, trim(comment_line))
   case(p_ftype%tmol)
      call write_tmol(self, unit, trim(comment_line))
   case(p_ftype%sdf)
      call write_sdf(self, unit, trim(comment_line))
   case(p_ftype%vasp)
      call write_vasp(self, unit, trim(comment_line))
   end select

end subroutine write_molecule_generic

subroutine write_xyz(mol, unit, comment_line)
   use mctc_econv
   class(tb_molecule), intent(in) :: mol
   integer, intent(in) :: unit
   character(len=*), intent(in) :: comment_line
   integer :: i

   write(unit, '(i0)') mol%n
   write(unit, '(a)') comment_line
   do i = 1, mol%n
      write(unit, '(a2,6x,3f20.14)') mol%sym(i), mol%xyz(:,i)*autoaa
   enddo

end subroutine write_xyz

subroutine write_tmol(mol, unit, comment_line)
   class(tb_molecule), intent(in) :: mol
   integer, intent(in) :: unit
   character(len=*), intent(in) :: comment_line
   integer :: i

   write(unit,'(a,1x,"#",1x,a)') "$coord bohr", comment_line
   do i = 1, mol%n
      write(unit,'(3f20.14,6x,a2)') mol%xyz(:,i), mol%sym(i)
   enddo
   write(unit,'(a,1x,i0)') "$periodic", mol%npbc
   if (mol%npbc > 0) then
      write(unit,'(a)') "$lattice bohr"
      write(unit,'(3f20.14)') mol%lattice
   endif
   write(unit,'(a)') "$end"

end subroutine write_tmol

subroutine write_sdf(mol, unit, comment_line)
   use mctc_econv
   class(tb_molecule), intent(in) :: mol
   integer, intent(in) :: unit
   character(len=*), intent(in) :: comment_line
   integer, parameter :: list4(4) = 0, list12(12) = 0
   integer :: iatom

   write(unit, '(a)')
   write(unit, '(a)') comment_line
   write(unit, '(a)')
   write(unit, '(3i3,3x,2i3,12x,i3,1x,a5)') &
      & mol%n, 0, 0, 0, 0, 999, 'V2000'

   do iatom = 1, mol%n
      write(unit, '(3f10.4,1x,a3,i2,11i3)') &
         & mol%xyz(:, iatom)*autoaa, mol%sym(iatom), list12
   enddo

   write(unit, '(a,*(1x,i3))') "M  CHG", 1, 1, nint(mol%chrg)
   write(unit, '(a)') "M  END"
   write(unit, '(a)') "$$$$"

end subroutine write_sdf

subroutine write_vasp(mol, unit, comment_line)
   use mctc_econv
   class(tb_molecule), intent(in) :: mol
   integer, intent(in) :: unit
   character(len=*), intent(in) :: comment_line
   integer :: i,j,iat
   integer,allocatable :: kinds(:)
   character(len=2),external :: asym

   allocate(kinds(mol%n), source = 1)

   ! use vasp 5.x format
   write(unit,'(a)') comment_line

   ! scaling factor for lattice parameters is always one
   write(unit,'(f20.14)') 1.0_wp
   ! write the lattice parameters
   do i = 1, 3
      write(unit,'(3f20.14)') mol%lattice(:,i)*autoaa
   enddo

   j = 0
   iat = 0
   do i = 1, mol%n
      if (iat.eq.mol%at(i)) then
         kinds(j) = kinds(j)+1
      else
         j = j+1
         iat = mol%at(i)
         write(unit,'(1x,a)',advance='no') mol%sym(i)
      endif
   enddo
   write(unit,'(a)')

   ! write the count of the consequtive atom types
   do i = 1, j
      write(unit,'(1x,i0)',advance='no') kinds(i)
   enddo
   write(unit,'(a)')
   deallocate(kinds)

   ! we write cartesian coordinates
   write(unit,'("Cartesian")')

   ! now write the cartesian coordinates
   do i = 1, mol%n
      write(unit,'(3f20.14)') mol%xyz(:,i)*autoaa
   enddo

end subroutine write_vasp

end submodule molecule_writer
