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

module xtb_io_reader_turbomole
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_boundaryconditions, only : boundaryCondition
   use xtb_mctc_constants
   use xtb_mctc_convert
   use xtb_mctc_resize
   use xtb_mctc_symbols, only : toNumber, toSymbol
   use xtb_pbc_tools
   use xtb_readin, getline => strip_line
   use xtb_type_molecule
   use xtb_type_vendordata, only : turbo_info
   implicit none
   private

   public :: readMoleculeCoord


   logical, parameter :: debug = .false.


contains


subroutine readMoleculeCoord(mol, unit, status, iomsg)
   type(TMolecule),intent(inout) :: mol
   integer,intent(in) :: unit !< file handle
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   character(len=:), allocatable :: line
   character(len=:), allocatable :: cell_string, lattice_string
   character(len=1), parameter :: flag = '$'
   character(len=4), allocatable :: sym(:)
   real(wp), allocatable :: coord(:, :)
   real(wp) :: latvec(9), conv
   integer :: error
   integer :: iatom, i, j
   integer :: periodic, cell_vectors
   integer, parameter :: p_initial_size = 100
   integer, parameter :: p_nlv(3) = [1, 4, 9]
   integer, parameter :: p_ncp(3) = [1, 3, 6]
   logical :: has_coord, has_periodic, has_lattice, has_cell
   logical :: cartesian, coord_in_bohr, lattice_in_bohr

   status = .false.

   allocate(sym(p_initial_size), source='    ')
   allocate(coord(3, p_initial_size), source=0.0_wp)

   iatom = 0
   periodic = 0
   has_coord = .false.
   has_periodic = .false.
   has_lattice = .false.
   has_cell = .false.
   cartesian = .true.
   coord_in_bohr = .true.
   lattice_in_bohr = .true.

   error = 0
   do while(error == 0)
      call getline(unit, line, error)
      if (index(line, flag) == 1) then
         if (index(line, 'end') == 2) then
            exit

         elseif (.not.has_coord .and. index(line, 'coord') == 2) then
            has_coord = .true.
            ! $coord angs / $coord bohr / $coord frac
            call select_unit(line, coord_in_bohr, cartesian)
            coord_group: do while(error == 0)
               call getline(unit, line, error)
               if (index(line, flag) == 1) then
                  backspace(unit)
                  exit coord_group
               endif
               if (iatom >= size(coord, 2)) call resize(coord)
               if (iatom >= size(sym)) call resize(sym)
               iatom = iatom + 1
               read(line, *, iostat=error) coord(:, iatom), sym(iatom)
            enddo coord_group

         elseif (.not.has_periodic .and. index(line, 'periodic') == 2) then
            has_periodic = .true.
            ! $periodic 0/1/2/3
            read(line(10:), *, iostat=error) periodic

         elseif (.not.has_lattice .and. index(line, 'lattice') == 2) then
            has_lattice = .true.
            ! $lattice bohr / $lattice angs
            call select_unit(line, lattice_in_bohr)
            cell_vectors = 0
            lattice_string = ''
            lattice_group: do while(error == 0)
               call getline(unit, line, error)
               if (index(line, flag) == 1) then
                  backspace(unit)
                  exit lattice_group
               endif
               cell_vectors = cell_vectors + 1
               lattice_string = lattice_string // ' ' // line
            enddo lattice_group

         elseif (.not.has_cell .and. index(line, 'cell') == 2) then
            has_cell = .true.
            ! $cell bohr / $cell angs
            call select_unit(line, lattice_in_bohr)
            call getline(unit, cell_string, error)
            if (debug) print*, cell_string

         endif
      endif
   enddo

   if (.not.has_coord .or. iatom == 0) then
      iomsg = "coordinates not present, cannot work without coordinates"
      return
   endif

   if (has_cell .and. has_lattice) then
      iomsg = "both lattice and cell group are present"
      return
   endif

   if (.not.has_periodic .and. (has_cell .or. has_lattice)) then
      iomsg = "cell and lattice definition is present, but periodicity is not given"
      return
   endif

   if (periodic > 0 .and. .not.(has_cell .or. has_lattice)) then
      iomsg = "system is periodic but definition of lattice/cell is missing"
      return
   endif

   if (.not.cartesian .and. periodic == 0) then
      iomsg = "fractional coordinates do not work for molecular systems"
      return
   endif

   call mol%allocate(iatom)

   mol%sym = sym(:len(mol))
   mol%at = toNumber(sym(:len(mol)))
   if (any(mol%at == 0)) then
      iomsg = "unknown element '"//sym(minval(mol%at))//"' present"
      return
   endif

   mol%npbc = periodic
   if (periodic > 0) then
      mol%pbc(:periodic) = .true.
      mol%boundaryCondition = boundaryCondition%pbc3d
   else
      mol%boundaryCondition = boundaryCondition%cluster
   end if

   if (has_cell) then
      read(cell_string, *, iostat=error) latvec(:p_ncp(periodic))
      if (debug) print*, latvec(:p_ncp(periodic))
      if (lattice_in_bohr) then
         conv = 1.0_wp
      else
         conv = aatoau
      endif
      select case(periodic)
      case(1)
         mol%cellpar = [latvec(1)*conv, 1.0_wp, 1.0_wp, &
            &           pi/2, pi/2, pi/2]
      case(2)
         mol%cellpar = [latvec(1)*conv, latvec(2)*conv, 1.0_wp, &
            &           pi/2, pi/2, latvec(3)*pi/180.0_wp]
      case(3)
         mol%cellpar = [latvec(1:3)*conv, latvec(4:6)*pi/180.0_wp]
      end select
      call cell_to_dlat(mol%cellpar,mol%lattice)
   endif

   if (has_lattice) then
      if (cell_vectors /= periodic) then
         iomsg = "number of cell vectors does not match periodicity"
         return
      endif
      read(lattice_string, *, iostat=error) latvec(:p_nlv(periodic))
      if (lattice_in_bohr) then
         conv = 1.0_wp
      else
         conv = aatoau
      endif
      j = 0
      do i = 1, periodic
         mol%lattice(:periodic,i) = latvec(j+1:j+periodic) * conv
         j = j + periodic
      enddo
   endif

   if (cartesian) then
      if (coord_in_bohr) then
         conv = 1.0_wp
      else
         conv = aatoau
      endif
      mol%xyz = coord(:, :len(mol)) * conv
      if (periodic > 0) &
         &call xyz_to_abc(len(mol),mol%lattice,mol%xyz,mol%abc,mol%pbc)
   else
      mol%abc = coord
      call abc_to_xyz(len(mol),mol%lattice,coord,mol%xyz)
   endif

   ! save data on input format
   mol%turbo = turbo_info(cartesian=cartesian, lattice=has_lattice, &
      &                   angs_lattice=.not.lattice_in_bohr, &
      &                   angs_coord=.not.coord_in_bohr)

   status = .true.

contains

   subroutine select_unit(line, in_bohr, cartesian)
      character(len=*), intent(in) :: line
      logical, intent(out) :: in_bohr
      logical, intent(out), optional :: cartesian
      in_bohr = index(line, ' angs') == 0
      if (present(cartesian)) cartesian = index(line, ' frac') == 0
   end subroutine select_unit

end subroutine readMoleculeCoord


end module xtb_io_reader_turbomole
