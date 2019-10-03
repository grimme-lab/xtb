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

!> provides definition of calculation options type
module tbdef_options
   use iso_fortran_env, wp => real64
   use iso_c_binding
   implicit none

   public :: tb_environment

   public :: tb_options

   public :: peeq_options
   public :: c_peeq_options
   public :: scc_options
   public :: c_scc_options
   public :: eeq_options
   public :: c_eeq_options
   public :: dftd_options
   public :: c_dftd_options
   public :: assignment(=)

   private

!> bundles all calculation options in one type
   type :: tb_options
      character(len=:),allocatable :: fname    !< geometry file
      character(len=:),allocatable :: xcontrol !< instruction file
      character(len=:),allocatable :: fnv      !< parameter file
      real(wp) :: acc = 0.0_wp                 !< calculation accuracy
      real(wp) :: etemp = 0.0_wp
      integer  :: prlevel = 0
      integer  :: maxiter = 0
      logical  :: copycontrol = .false.        !< copy control file
      logical  :: restart = .false.            !< print restart file
      logical  :: grad = .false.               !< calculate gradient
      logical  :: define = .false.             !< check the current run
      logical  :: verbose = .false.            !< print more information
      logical  :: veryverbose = .false.        !< clutter the screen more
      logical  :: silent = .false.             !< clutter the screen less
   contains
   procedure :: allocate => allocate_options
   procedure :: deallocate => deallocate_options
   procedure :: write => write_options
   procedure :: default_options
   procedure :: export_peeq => export_peeq_options
   procedure :: export_scc  => export_scc_options
   end type tb_options

   type :: tb_environment
      character(len=:),allocatable :: whoami
      character(len=:),allocatable :: hostname
      character(len=:),allocatable :: home
      character(len=:),allocatable :: path
      character(len=:),allocatable :: xtbpath
      character(len=:),allocatable :: xtbhome
   contains
      procedure :: setup => setup_environment
   end type tb_environment

   type :: scc_options
      sequence
      integer  :: prlevel = 0
      integer  :: parallel = 0
      real(wp) :: acc = 0.0_wp
      real(wp) :: etemp = 0.0_wp
      logical  :: grad = .false.
      logical  :: restart = .false.
      logical  :: ccm = .false.
      integer  :: maxiter = 0
      character(len=20) :: solvent = "none"
   end type scc_options

   type,bind(C) :: c_scc_options
      integer(c_int)  :: prlevel = 0_c_int
      integer(c_int)  :: parallel = 0_c_int
      real(c_double)  :: acc = 0.0_c_double
      real(c_double)  :: etemp = 0.0_c_double
      logical(c_bool) :: grad = .false._c_bool
      logical(c_bool) :: restart = .false._c_bool
      logical(c_bool) :: ccm = .false._c_bool
      integer(c_int)  :: maxiter = 0_c_int
      character(c_char) :: solvent(20) = c_null_char
   end type c_scc_options

   type :: peeq_options
      sequence
      integer  :: prlevel = 0
      integer  :: parallel = 0
      real(wp) :: acc = 0.0_wp
      real(wp) :: etemp = 0.0_wp
      logical  :: grad = .false.
      logical  :: ccm = .false.
      character(len=20) :: solvent = "none"
   end type peeq_options

   type,bind(C) :: c_peeq_options
      integer(c_int)  :: prlevel = 0_c_int
      integer(c_int)  :: parallel = 0_c_int
      real(c_double)  :: acc = 0.0_c_double
      real(c_double)  :: etemp = 0.0_c_double
      logical(c_bool) :: grad = .false._c_bool
      logical(c_bool) :: ccm = .false._c_bool
      character(c_char) :: solvent(20) = c_null_char
   end type c_peeq_options

   type :: eeq_options
      sequence
      integer  :: prlevel = 0
      integer  :: param = 0
      logical  :: grad = .false.
   end type eeq_options

   type,bind(C) :: c_eeq_options
      integer(c_int)  :: prlevel = 0_c_int
      integer(c_int)  :: param = 0_c_int
      logical(c_bool) :: grad = .false._c_bool
   end type c_eeq_options

!> calculation setup
   type :: dftd_options
      sequence
      integer  :: lmbd = -1                 !< kind of non-additivity correction
      integer  :: refq = -1                 !< kind of charge model
      real(wp) :: wf = 0.0_wp               !< weighting factor
      real(wp) :: g_a = 0.0_wp              !< charge scale height
      real(wp) :: g_c = 0.0_wp              !< charge scale steepness
      logical  :: lmolpol = .false.         !< calculate molecular properties?
      logical  :: lenergy = .false.         !< calculate dispersion energy?
      logical  :: lgradient = .false.       !< calculate dispersion gradient?
      logical  :: verbose = .false.         !< print more information
      logical  :: veryverbose = .false.     !< clutter the screen more
      logical  :: silent = .false.          !< clutter the screen less
   end type dftd_options

!> calculation setup
   type,bind(C) :: c_dftd_options
      integer(c_int)  :: lmbd = -1             !< kind of non-additivity correction
      integer(c_int)  :: refq = -1             !< kind of charge model
      real(c_double)  :: wf = 0.0_wp           !< weighting factor
      real(c_double)  :: g_a = 0.0_wp          !< charge scale height
      real(c_double)  :: g_c = 0.0_wp          !< charge scale steepness
      logical(c_bool) :: lmolpol = .false.     !< calculate molecular properties?
      logical(c_bool) :: lenergy = .false.     !< calculate dispersion energy?
      logical(c_bool) :: lgradient = .false.   !< calculate dispersion gradient?
      logical(c_bool) :: verbose = .false.     !< print more information
      logical(c_bool) :: veryverbose = .false. !< clutter the screen more
      logical(c_bool) :: silent = .false.      !< clutter the screen less
   end type c_dftd_options

   interface assignment(=)
      module procedure :: convert_peeq_options_c_to_f
      module procedure :: convert_peeq_options_f_to_c
      module procedure :: convert_scc_options_c_to_f
      module procedure :: convert_scc_options_f_to_c
      module procedure :: convert_eeq_options_c_to_f
      module procedure :: convert_eeq_options_f_to_c
      module procedure :: convert_dftd_options_c_to_f
      module procedure :: convert_dftd_options_f_to_c
   end interface

contains

subroutine setup_environment(self)
   use mctc_systools
   implicit none
   class(tb_environment), intent(inout) :: self
   integer :: err
   call rdarg(0,self%whoami,err)
   call rdvar('HOSTNAME',self%hostname,err)
   call rdvar('HOME',self%home,err)
   call rdvar('PATH',self%path,err)
   call rdvar('XTBHOME',self%xtbhome,err)
   if (err.ne.0) self%xtbhome = self%home
   call rdvar('XTBPATH',self%xtbpath,err)
   if (err.ne.0) self%xtbpath = self%xtbhome
end subroutine setup_environment

subroutine default_options(self)
   implicit none
   class(tb_options),intent(out) :: self

   self%acc = 1.0_wp
   self%restart = .true.
   self%copycontrol = .false.
   self%grad = .false.
   self%etemp = 300.0_wp
   self%maxiter = 250

end subroutine default_options

!> constructor for calculation options
subroutine allocate_options(self)
   implicit none
   class(tb_options),intent(inout) :: self
   call self%deallocate
end subroutine allocate_options

!> deconstructor for calculation options
subroutine deallocate_options(self)
   implicit none
   class(tb_options),intent(inout) :: self
   if (allocated(self%fname))    deallocate(self%fname)
   if (allocated(self%xcontrol)) deallocate(self%xcontrol)
   if (allocated(self%fnv))      deallocate(self%fnv)
end subroutine deallocate_options

!> print information about current calculation options to unit
subroutine write_options(self,iunit,comment)
   implicit none
   class(tb_options),intent(in) :: self    !< calculation options
   integer,          intent(in) :: iunit   !< file handle
   character(len=*), intent(in) :: comment !< name of the variable
   character(len=*),parameter :: dfmt = '(1x,a,1x,"=",1x,g0)'

   write(iunit,'(72(">"))')
   write(iunit,'(1x,"*",1x,a)') "Writing 'tb_options' class"
   write(iunit,'(  "->",1x,a)') comment
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "status of the fields"
   write(iunit,dfmt) "real    :: acc         ",self%acc
   write(iunit,dfmt) "logical :: copycontrol ",self%copycontrol
   write(iunit,dfmt) "logical :: restart     ",self%restart
   write(iunit,dfmt) "logical :: grad        ",self%grad
   write(iunit,dfmt) "logical :: define      ",self%define
   write(iunit,dfmt) "logical :: verbose     ",self%verbose
   write(iunit,dfmt) "logical :: veryverbose ",self%veryverbose
   write(iunit,dfmt) "logical :: silent      ",self%silent
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "allocation status"
   write(iunit,dfmt) "allocated? fname       ",allocated(self%fname)
   write(iunit,dfmt) "allocated? xcontrol    ",allocated(self%xcontrol)
   write(iunit,dfmt) "allocated? fnv         ",allocated(self%fnv)
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "size of memory allocation"
   if (allocated(self%fname)) then
   write(iunit,dfmt) "length  :: fname       ",len(self%fname)
   endif
   if (allocated(self%xcontrol)) then
   write(iunit,dfmt) "length  :: xcontrol    ",len(self%xcontrol)
   endif
   if (allocated(self%fnv)) then
   write(iunit,dfmt) "length  :: fnv         ",len(self%fnv)
   endif
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "string content"
   if (allocated(self%fname)) then
   write(iunit,dfmt) "string  :: fname       ",self%fname
   endif
   if (allocated(self%xcontrol)) then
   write(iunit,dfmt) "string  :: xcontrol    ",self%xcontrol
   endif
   if (allocated(self%fnv)) then
   write(iunit,dfmt) "string  :: fnv         ",self%fnv
   endif
   write(iunit,'(72("<"))')
end subroutine write_options

pure function export_scc_options(self) result(opt)
   implicit none
   class(tb_options), intent(in) :: self
   type(scc_options) :: opt
   opt%prlevel = self%prlevel
   opt%maxiter = self%maxiter
   opt%acc     = self%acc
   opt%etemp   = self%etemp
   opt%grad    = self%grad
end function export_scc_options

pure elemental subroutine convert_scc_options_c_to_f &
      (f_opt,c_opt)
   use iso_c_binding, only : c_null_char
   implicit none
   type(c_scc_options), intent(in)  :: c_opt
   type(scc_options),   intent(out) :: f_opt
   integer :: i
   f_opt%prlevel = c_opt%prlevel
   f_opt%maxiter = c_opt%maxiter
   f_opt%acc     = c_opt%acc
   f_opt%etemp   = c_opt%etemp
   f_opt%grad    = c_opt%grad
   f_opt%restart = c_opt%restart
   f_opt%parallel= c_opt%parallel
   f_opt%ccm     = c_opt%ccm

   f_opt%solvent = ''
   do i = 1, 20
      if (c_opt%solvent(i).eq.c_null_char) exit
      f_opt%solvent(i:i) = c_opt%solvent(i)
   enddo

end subroutine convert_scc_options_c_to_f

pure elemental subroutine convert_scc_options_f_to_c &
      (c_opt,f_opt)
   use iso_c_binding, only : c_null_char
   implicit none
   type(scc_options),  intent(in)  :: f_opt
   type(c_scc_options),intent(out) :: c_opt
   integer :: i
   c_opt%prlevel = f_opt%prlevel
   c_opt%maxiter = f_opt%maxiter
   c_opt%acc     = f_opt%acc
   c_opt%etemp   = f_opt%etemp
   c_opt%grad    = f_opt%grad
   c_opt%restart = f_opt%restart
   c_opt%parallel= f_opt%parallel
   c_opt%ccm     = f_opt%ccm

   do i = 1, len_trim(f_opt%solvent)
      c_opt%solvent(i) = f_opt%solvent(i:i)
   enddo
   c_opt%solvent(len_trim(f_opt%solvent)+1) = c_null_char

end subroutine convert_scc_options_f_to_c

pure function export_peeq_options(self) result(opt)
   implicit none
   class(tb_options), intent(in) :: self
   type(peeq_options) :: opt
   opt%prlevel = self%prlevel
   opt%acc     = self%acc
   opt%etemp   = self%etemp
   opt%grad    = self%grad
end function export_peeq_options

pure elemental subroutine convert_peeq_options_c_to_f &
      (f_opt,c_opt)
   implicit none
   type(c_peeq_options), intent(in)  :: c_opt
   type(peeq_options),   intent(out) :: f_opt
   integer :: i
   f_opt%prlevel = c_opt%prlevel
   f_opt%ccm     = c_opt%ccm
   f_opt%acc     = c_opt%acc
   f_opt%etemp   = c_opt%etemp
   f_opt%grad    = c_opt%grad
   f_opt%parallel= c_opt%parallel

   f_opt%solvent = ''
   do i = 1, 20
      if (c_opt%solvent(i).eq.c_null_char) exit
      f_opt%solvent(i:i) = c_opt%solvent(i)
   enddo
end subroutine convert_peeq_options_c_to_f

pure elemental subroutine convert_peeq_options_f_to_c &
      (c_opt,f_opt)
   implicit none
   type(peeq_options),  intent(in)  :: f_opt
   type(c_peeq_options),intent(out) :: c_opt
   integer :: i
   c_opt%prlevel = f_opt%prlevel
   c_opt%ccm     = f_opt%ccm
   c_opt%acc     = f_opt%acc
   c_opt%etemp   = f_opt%etemp
   c_opt%grad    = f_opt%grad
   c_opt%parallel= f_opt%parallel

   do i = 1, len_trim(f_opt%solvent)
      c_opt%solvent(i) = f_opt%solvent(i:i)
   enddo
   c_opt%solvent(len_trim(f_opt%solvent)+1) = c_null_char
end subroutine convert_peeq_options_f_to_c

pure elemental subroutine convert_eeq_options_c_to_f &
      (f_opt,c_opt)
   implicit none
   type(c_eeq_options), intent(in)  :: c_opt
   type(eeq_options),   intent(out) :: f_opt
   f_opt%prlevel = c_opt%prlevel
   f_opt%param   = c_opt%param
   f_opt%grad    = c_opt%grad
end subroutine convert_eeq_options_c_to_f

pure elemental subroutine convert_eeq_options_f_to_c &
      (c_opt,f_opt)
   implicit none
   type(eeq_options),  intent(in)  :: f_opt
   type(c_eeq_options),intent(out) :: c_opt
   c_opt%prlevel = f_opt%prlevel
   c_opt%param   = f_opt%param
   c_opt%grad    = f_opt%grad
end subroutine convert_eeq_options_f_to_c

pure elemental subroutine convert_dftd_options_c_to_f &
      (f_dopt,c_dopt)
   implicit none
   type(c_dftd_options),intent(in)  :: c_dopt
   type(dftd_options),  intent(out) :: f_dopt
   f_dopt%lmbd        = c_dopt%lmbd
   f_dopt%refq        = c_dopt%refq
   f_dopt%wf          = c_dopt%wf
   f_dopt%g_a         = c_dopt%g_a
   f_dopt%g_c         = c_dopt%g_c
   f_dopt%lmolpol     = c_dopt%lmolpol
   f_dopt%lenergy     = c_dopt%lenergy
   f_dopt%lgradient   = c_dopt%lgradient
   f_dopt%verbose     = c_dopt%verbose
   f_dopt%veryverbose = c_dopt%veryverbose
   f_dopt%silent      = c_dopt%silent
end subroutine convert_dftd_options_c_to_f

pure elemental subroutine convert_dftd_options_f_to_c &
      (c_dopt,f_dopt)
   implicit none
   type(c_dftd_options),intent(out) :: c_dopt
   type(dftd_options),  intent(in)  :: f_dopt
   c_dopt%lmbd        = f_dopt%lmbd
   c_dopt%refq        = f_dopt%refq
   c_dopt%wf          = f_dopt%wf
   c_dopt%g_a         = f_dopt%g_a
   c_dopt%g_c         = f_dopt%g_c
   c_dopt%lmolpol     = f_dopt%lmolpol
   c_dopt%lenergy     = f_dopt%lenergy
   c_dopt%lgradient   = f_dopt%lgradient
   c_dopt%verbose     = f_dopt%verbose
   c_dopt%veryverbose = f_dopt%veryverbose
   c_dopt%silent      = f_dopt%silent
end subroutine convert_dftd_options_f_to_c

end module tbdef_options
