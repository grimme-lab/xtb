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

module tbdef_data
   use iso_fortran_env, wp => real64

   implicit none

   public :: scc_results
   public :: freq_results

   private

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
      real(wp) :: molpol = 0.0_wp
      real(wp) :: g_solv = 0.0_wp
      real(wp) :: g_total = 0.0_wp
      real(wp) :: gnorm = 0.0_wp
      logical  :: converged = .true.
   end type scc_results

   type freq_results
      integer  :: n = 0
      integer  :: n3 = 0
      integer  :: n3true = 0
      integer  :: lowmode = 0
      integer  :: nimag = 0
      logical  :: linear = .false.
      logical  :: mweighted = .false.
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

end module tbdef_data
