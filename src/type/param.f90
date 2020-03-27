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

module xtb_type_param
   use xtb_mctc_accuracy, only : wp
   implicit none

   public :: dftd_parameter
   public :: TxTBParameter
   public ::  scc_parameter

   public :: chrg_parameter
   private

   type :: TxTBParameter
      real(wp) :: kshell(0:3) = 0.0_wp
      real(wp) :: ksp = 0.0_wp
      real(wp) :: ksd = 0.0_wp
      real(wp) :: kpd = 0.0_wp
      real(wp) :: kdiff = 0.0_wp
      real(wp) :: kdiffa = 0.0_wp
      real(wp) :: kdiffb = 0.0_wp
      real(wp) :: enshell(0:3) = 0.0_wp
      real(wp) :: enscale4 = 0.0_wp
      real(wp) :: cnshell(2, 0:3) = 0.0_wp
      real(wp) :: gam3shell(2, 0:3) = 0.0_wp
      real(wp) :: srbshift = 0.0_wp
      real(wp) :: srbpre = 0.0_wp
      real(wp) :: srbexp = 0.0_wp
      real(wp) :: srbken = 0.0_wp
      real(wp) :: wllscal = 0.0_wp
      real(wp) :: gscal = 0.0_wp
      real(wp) :: zcnf = 0.0_wp
      real(wp) :: tscal = 0.0_wp
      real(wp) :: kcn = 0.0_wp
      real(wp) :: fpol = 0.0_wp
      real(wp) :: ken = 0.0_wp
      real(wp) :: lshift = 0.0_wp
      real(wp) :: lshifta = 0.0_wp
      real(wp) :: split = 0.0_wp
      real(wp) :: zqf = 0.0_wp
      real(wp) :: alphaj = 0.0_wp
      real(wp) :: kexpo = 0.0_wp
      real(wp) :: dispa = 0.0_wp
      real(wp) :: dispb = 0.0_wp
      real(wp) :: dispc = 0.0_wp
      real(wp) :: dispatm = 0.0_wp
      real(wp) :: xbdamp = 0.0_wp
      real(wp) :: xbrad = 0.0_wp
      real(wp) :: aesshift = 0.0_wp
      real(wp) :: aesexp = 0.0_wp
      real(wp) :: aesrmax = 0.0_wp
      real(wp) :: aesdmp3 = 0.0_wp
      real(wp) :: aesdmp5 = 0.0_wp
      real(wp) :: ipeashift = 0.0_wp
      real(wp) :: renscale = 0.0_wp
   end type TxTBParameter

   type :: dftd_parameter
      real(wp) :: s6  = -1.0_wp
      real(wp) :: s8  = -1.0_wp
      real(wp) :: s10 =  0.0_wp
      real(wp) :: a1  = -1.0_wp
      real(wp) :: a2  = -1.0_wp
      real(wp) :: s9  =  1.0_wp
      integer  :: alp = 16
      ! for MBD@rsSCS
      real(wp) :: beta = 1.0_wp
   end type dftd_parameter

   type :: scc_parameter
      real(wp) :: gam3l(0:3)
      real(wp) :: ipshift
      real(wp) :: eashift
      real(wp) :: g_a
      real(wp) :: g_c
      real(wp) :: wf
   end type scc_parameter

   type chrg_parameter
      integer  :: n
      real(wp),allocatable :: en(:)
      real(wp),allocatable :: gam(:)
      real(wp),allocatable :: kappa(:)
      real(wp),allocatable :: alpha(:)
      real(wp),allocatable :: dpol(:)
      real(wp),allocatable :: beta(:)
   contains
      procedure :: allocate => allocate_chrgeq
      procedure :: deallocate => deallocate_chrgeq
   end type chrg_parameter

contains

subroutine allocate_chrgeq(self,n,extended)
   implicit none
   class(chrg_parameter) :: self
   integer,intent(in) :: n
   logical,intent(in),optional :: extended
   logical :: multipoles
   if (present(extended)) then
      multipoles = extended
   else
      multipoles = .false.
   endif
   call self%deallocate
   allocate( self%en(n),    source = 0.0_wp )
   allocate( self%gam(n),   source = 0.0_wp )
   allocate( self%kappa(n), source = 0.0_wp )
   allocate( self%alpha(n), source = 0.0_wp )
   if (multipoles) then
   allocate( self%dpol(n),  source = 0.0_wp )
   allocate( self%beta(n),  source = 0.0_wp )
   endif
end subroutine allocate_chrgeq

subroutine deallocate_chrgeq(self)
   implicit none
   class(chrg_parameter) :: self
   if (allocated(self%en))    deallocate(self%en)
   if (allocated(self%gam))   deallocate(self%gam)
   if (allocated(self%kappa)) deallocate(self%kappa)
   if (allocated(self%alpha)) deallocate(self%alpha)
   if (allocated(self%dpol))  deallocate(self%dpol)
   if (allocated(self%beta))  deallocate(self%beta)
end subroutine deallocate_chrgeq

end module xtb_type_param
