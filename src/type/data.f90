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

module xtb_type_data
   use xtb_mctc_accuracy, only : wp
   use xtb_type_pcem
   use xtb_iff_data, only : TIFFData

   implicit none

   public :: scc_results
   public :: freq_results

   private

   type :: TIFFResults
      !> Number of atoms
      integer :: n       
      !> Ordinary numbers
      integer, allocatable :: at(:)
      !> Coordinates 
      real(wp), allocatable :: xyz(:, :)
      !> Charges
      real(wp), allocatable :: q(:)
      !> Number of LMOs 
      integer :: nlmo
      !> LMO positions 
      real(wp), allocatable :: rlmo(:, :)
      !> LMO values
      integer, allocatable :: lmo(:)
      !> Charge related stuff
      real(wp), allocatable :: qct(:, :)
      !> HOMO, LUMO, Dipol
      real(wp) :: elumo
      real(wp) :: ehomo
      real(wp) :: dipol

   contains

      procedure :: delete
      procedure :: allocateIFFResults

   end type TIFFResults

   type :: scc_results
      real(wp) :: e_atom = 0.0_wp
      real(wp) :: e_elec = 0.0_wp
      real(wp) :: e_total = 0.0_wp
      real(wp) :: e_rep = 0.0_wp
      real(wp) :: e_es = 0.0_wp
      real(wp) :: e_aes = 0.0_wp
      real(wp) :: e_axc = 0.0_wp
      real(wp) :: e_disp = 0.0_wp
      real(wp) :: e_xb = 0.0_wp
      real(wp) :: g_born = 0.0_wp
      real(wp) :: g_sasa = 0.0_wp
      real(wp) :: g_hb = 0.0_wp
      real(wp) :: g_shift = 0.0_wp
      real(wp) :: hl_gap = 0.0_wp
      real(wp) :: dipole(3) = (/0.0_wp,0.0_wp,0.0_wp/)
      real(wp) :: alpha(3,3) = reshape((/0.0_wp, 0.0_wp, 0.0_wp, &
         & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp/), (/3, 3/))
      real(wp) :: quadrupole(6) = (/0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp/)
      real(wp) :: molpol = 0.0_wp
      real(wp) :: g_solv = 0.0_wp
      real(wp) :: g_total = 0.0_wp
      real(wp) :: gnorm = 0.0_wp
      logical  :: converged = .true.
      type(tb_pcem) :: pcem = tb_pcem()
      real(wp) :: e_bond = 0.0_wp
      real(wp) :: e_angl = 0.0_wp
      real(wp) :: e_tors = 0.0_wp
      real(wp) :: e_hb = 0.0_wp
      real(wp) :: e_batm = 0.0_wp
      real(wp) :: e_ext = 0.0_wp
      type(TIFFResults), allocatable :: iff_results
   contains 
      procedure :: print => print_scc_results
   end type scc_results

   type freq_results
      integer  :: n = 0
      integer  :: n3 = 0
      integer  :: n3true = 0
      integer  :: lowmode = 0
      integer  :: nimag = 0
      logical  :: linear = .false.
      logical  :: mweighted = .false.
      real(wp) :: zp = 0.0_wp
      real(wp) :: temp = 0.0_wp
      real(wp) :: etot = 0.0_wp
      real(wp) :: htot = 0.0_wp
      real(wp) :: gtot = 0.0_wp
      real(wp) :: gnorm = 0.0_wp
      real(wp),allocatable :: grad(:,:)
      real(wp),allocatable :: freq(:)
      real(wp),allocatable :: hess(:,:)
      real(wp),allocatable :: rmass(:)
      real(wp),allocatable :: dipt(:)
      real(wp),allocatable :: polt(:)
      character(len=:),allocatable :: pg
   contains
      procedure :: allocate => allocate_freq_results
      procedure :: deallocate => deallocate_freq_results
   end type freq_results

contains
   
subroutine delete(self)
   class(TIFFResults), intent(out) :: self
   if (allocated(self%at)) deallocate (self%at)
   if (allocated(self%xyz)) deallocate (self%xyz)
   if (allocated(self%q)) deallocate (self%q)
   if (allocated(self%lmo)) deallocate (self%lmo)
   if (allocated(self%rlmo)) deallocate (self%rlmo)
   if (allocated(self%qct)) deallocate (self%qct)
end subroutine delete

subroutine allocateIFFResults(self, n)
   use xtb_type_molecule, only: TMolecule
   class(TIFFResults), intent(out) :: self
   integer, intent(in) :: n
   self%n = n
   allocate(self%at(n), self%lmo(10*n), source = 0)

   allocate(self%xyz(3,n), self%q(n), self%rlmo(4,10*n),&
           & self%qct(n,2), source= 0.0_wp)
end subroutine allocateIFFResults

subroutine allocate_freq_results(self,n)
   implicit none
   class(freq_results) :: self
   integer,intent(in)  :: n
   call self%deallocate
   self%n      = n
   self%n3     = 3*n
   self%n3true = 3*n
   self%pg = 'c1'
   allocate( self%grad (3,n),     source = 0.0_wp )
   allocate( self%hess (3*n,3*n), source = 0.0_wp )
   allocate( self%freq (3*n),     source = 0.0_wp )
   allocate( self%rmass(3*n),     source = 0.0_wp )
   allocate( self%dipt (3*n),     source = 0.0_wp )
   allocate( self%polt (3*n),     source = 0.0_wp )
end subroutine allocate_freq_results

subroutine deallocate_freq_results(self)
   implicit none
   class(freq_results) :: self
   self%n = 0
   self%n3 = 0
   self%n3true = 0
   self%nimag = 0
   if (allocated( self%grad )) deallocate( self%grad )
   if (allocated( self%hess )) deallocate( self%hess )
   if (allocated( self%freq )) deallocate( self%freq )
   if (allocated( self%rmass)) deallocate( self%rmass)
   if (allocated( self%dipt )) deallocate( self%dipt )
   if (allocated( self%polt )) deallocate( self%polt )
   if (allocated( self%pg   )) deallocate( self%pg   )
end subroutine deallocate_freq_results

!> print SCC results (for debug)
subroutine print_scc_results(self, unit)
   use iso_fortran_env, only : output_unit
   class(scc_results) :: self
   integer, optional :: unit

   integer :: out ! output unit holder
   character(len=*), parameter :: dfmt = '(3x,a,2x,f12.6)' ! format for double precision
   integer :: i
   
   if (present(unit)) then
      out = unit
   else
      out = output_unit
   endif

   write(out, '(3x,a,/,2x,11("="))') 'SCC Results'

   write(out, dfmt) 'e_total   :', self%e_total
   write(out, dfmt) 'hl_gap    :', self%hl_gap
   do i=1,3
      write(out, dfmt) 'dipole    :', self%dipole(i)
   enddo

end subroutine print_scc_results

end module xtb_type_data
