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

module xtb_type_anc
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas
   use xtb_mctc_lapack, only : lapack_spevd, lapack_syevd
   implicit none

   public :: tb_anc
   private

   type tb_anc
      integer  :: n    !< number of atoms
      integer  :: n3   !< dimension of hessian
      integer  :: nvar !< actual dimension
      real(wp) :: hlow
      real(wp) :: hmax
      real(wp),allocatable :: hess(:)
      real(wp),allocatable :: B(:,:)
      real(wp),allocatable :: eigv(:)
      real(wp),allocatable :: coord(:)
      real(wp),allocatable :: xyz(:,:)
   contains
   procedure :: allocate => allocate_anc
   procedure :: deallocate => deallocate_anc
   procedure :: write => write_anc
   procedure :: new => generate_anc_blowup
   procedure :: get_cartesian
   end type tb_anc

contains

subroutine allocate_anc(self,n,nvar,hlow,hmax)
   implicit none
   class(tb_anc),intent(inout)  :: self
   integer,      intent(in)     :: n
   integer,      intent(in)     :: nvar
   integer                      :: n3
   real(wp),intent(in),optional :: hlow
   real(wp),intent(in),optional :: hmax
   n3 = 3*n
   call self%deallocate
   self%n    = n
   self%n3   = 3*n
   self%nvar = nvar
   if (present(hlow)) self%hlow = hlow
   if (present(hmax)) self%hmax = hmax
   allocate( self%hess(nvar*(nvar+1)/2), source = 0.0_wp )
   allocate( self%B(n3,n3),              source = 0.0_wp )
   allocate( self%eigv(n3),              source = 0.0_wp )
   allocate( self%coord(nvar),           source = 0.0_wp )
   allocate( self%xyz(3,n),              source = 0.0_wp )
end subroutine allocate_anc

subroutine deallocate_anc(self)
   implicit none
   class(tb_anc),intent(inout) :: self
   self%n    = 0
   self%n3   = 0
   self%nvar = 0
   if (allocated(self%hess )) deallocate( self%hess  )
   if (allocated(self%B    )) deallocate( self%B     )
   if (allocated(self%eigv )) deallocate( self%eigv  )
   if (allocated(self%coord)) deallocate( self%coord )
   if (allocated(self%xyz  )) deallocate( self%xyz   )
end subroutine deallocate_anc

!> @brief print information about current approximate normal coordinates to unit
subroutine write_anc(self,iunit,comment)
   implicit none
   class(tb_anc),   intent(in) :: self    !< approximate normal coordinates
   integer,         intent(in) :: iunit   !< file handle
   character(len=*),intent(in) :: comment !< name of the variable
   character(len=*),parameter  :: dfmt = '(1x,a,1x,"=",1x,g0)'

   write(iunit,'(72(">"))')
   write(iunit,'(1x,"*",1x,a)') "Writing 'tb_anc' class"
   write(iunit,'(  "->",1x,a)') comment
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "status of the fields"
   write(iunit,dfmt) "integer :: n           ",self%n
   write(iunit,dfmt) "integer :: n3          ",self%n3
   write(iunit,dfmt) "integer :: nvar        ",self%nvar
   write(iunit,dfmt) "real    :: hlow        ",self%hlow
   write(iunit,dfmt) "real    :: hmax        ",self%hmax
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "allocation status"
   write(iunit,dfmt) "allocated? hess(:)     ",allocated(self%hess)
   write(iunit,dfmt) "allocated? B(:)        ",allocated(self%B)
   write(iunit,dfmt) "allocated? eigv(:)     ",allocated(self%eigv)
   write(iunit,dfmt) "allocated? coord(:)    ",allocated(self%coord)
   write(iunit,dfmt) "allocated? xyz(:,:)    ",allocated(self%xyz)
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "size of memory allocation"
   if (allocated(self%hess)) then
   write(iunit,dfmt) "size(1) :: hess(*)     ",size(self%hess,1)
   endif
   if (allocated(self%B)) then
   write(iunit,dfmt) "size(1) :: B(*,:)      ",size(self%B,1)
   write(iunit,dfmt) "size(2) :: B(:,*)      ",size(self%B,2)
   endif
   if (allocated(self%eigv)) then
   write(iunit,dfmt) "size(1) :: eigv(*)     ",size(self%eigv,1)
   endif
   if (allocated(self%coord)) then
   write(iunit,dfmt) "size(1) :: coord(*)    ",size(self%coord,1)
   endif
   if (allocated(self%xyz)) then
   write(iunit,dfmt) "size(1) :: xyz(*,:)    ",size(self%xyz,1)
   write(iunit,dfmt) "size(2) :: xyz(:,*)    ",size(self%xyz,2)
   endif
   write(iunit,'(72("<"))')
end subroutine write_anc

subroutine generate_anc_blowup(self,iunit,xyz,hess,pr,linear)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_la
   use xtb_detrotra, only : detrotra8
   implicit none
   class(tb_anc),intent(inout) :: self
   integer,      intent(in)    :: iunit
   real(wp),     intent(in)    :: xyz(3,self%n)
   real(wp),     intent(inout) :: hess(self%n3,self%n3)
   logical,      intent(in)    :: pr
   logical,      intent(in)    :: linear

   real(wp),parameter   :: thr1 = 1.0e-10_wp
   real(wp),parameter   :: thr2 = 1.0e-11_wp
   integer, parameter   :: maxtry = 4
   integer  :: i,itry
   integer  :: nvar
   integer  :: info
   integer  :: lwork
   integer  :: liwork
   logical  :: fail
   integer, allocatable :: iwork(:)
   real(wp) :: elow,damp,thr
   real(wp),allocatable :: aux(:)

   self%xyz = xyz

   thr = thr2
   lwork  = 1 + 6*self%n3 + 2*self%n3**2
   liwork = 8 * self%n3

   allocate(iwork(liwork), source = 0 )
   allocate(aux(lwork), source = 0.0_wp )

   call lapack_syevd('V','U',self%n3,hess,self%n3,self%eigv, &
      &        aux,lwork,iwork,liwork,info)

   call detrotra8(linear,self%n,self%xyz,hess,self%eigv) 

   !elow = 1.0e+99_wp
   elow = minval(self%eigv,mask=(abs(self%eigv) > thr1))
!   do i = 1, self%n3
!      if (abs(self%eigv(i)) > thr1 ) elow = min(elow,self%eigv(i))
!   enddo

   damp = max(self%hlow - elow,0.0_wp) 
   where(abs(self%eigv) > thr2) self%eigv = self%eigv + damp
!   do i = 1, self%n3
!      if (abs(self%eigv(i)) > thr2 ) self%eigv(i) = self%eigv(i) + damp
!   enddo

   if(pr)then
      write(iunit,*) 'Shifting diagonal of input Hessian by ', damp
      write(iunit,*) 'Lowest  eigenvalues of input Hessian'
      write(iunit,'(6F12.6)')(self%eigv(i),i=1,min(18,self%n3))
      write(iunit,*) 'Highest eigenvalues'
      write(iunit,'(6F12.6)')(self%eigv(i),i=self%n3-5,self%n3)
      write(iunit,*)
   endif

   fail = .true.
   get_anc: do itry = 1, maxtry
      self%B = 0.0_wp
      self%hess = 0.0_wp
      nvar = 0
      ! take largest (positive) first
      do i = self%n3, 1, -1
         if (abs(self%eigv(i)) > thr .and. nvar < self%nvar) then
            nvar = nvar+1
            self%B(:,nvar) = hess(:,i)
            self%hess(nvar+nvar*(nvar-1)/2) = &
               min(max(self%eigv(i),self%hlow),self%hmax)
         endif
      enddo

      if (nvar.ne.self%nvar) then
         thr = thr * 0.1_wp
         cycle get_anc
      endif

      fail = .false.
      exit get_anc

   enddo get_anc

   if (fail) then
      write(*,*) 'nvar, selv%nvar',nvar,self%nvar
      call raise('E',"ANC generation failed!")
   end if

   call sort(self%n3,self%nvar,self%hess,self%B)

   self%coord = 0.0_wp
end subroutine generate_anc_blowup

subroutine generate_anc_packed(self,xyz,hess,pr)
   use xtb_mctc_io, only : stdout
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_la
   implicit none
   class(tb_anc),intent(inout) :: self
   real(wp),     intent(in)    :: xyz(3,self%n)
   real(wp),     intent(inout) :: hess(self%n3*(self%n3+1)/2)
   logical,      intent(in)    :: pr

   real(wp),parameter   :: thr1 = 1.0e-10_wp
   real(wp),parameter   :: thr2 = 1.0e-11_wp
   integer, parameter   :: maxtry = 4
   integer  :: i,itry
   integer  :: nvar
   integer  :: info
   integer  :: lwork
   integer  :: liwork
   logical  :: fail
   integer, allocatable :: iwork(:)
   real(wp) :: elow,damp,thr
   real(wp),allocatable :: aux(:)
   real(wp),allocatable :: u(:,:)

   self%xyz = xyz

   thr = thr2
   lwork  = 1 + 6*self%n3 + 2*self%n3**2
   liwork = 8 * self%n3

   allocate(iwork(liwork), source = 0 )
   allocate(aux(lwork), source = 0.0_wp )
   allocate(u(self%n3,self%n3), source = 0.0_wp )

   call lapack_spevd('V','U',self%n3,hess,self%eigv,u,self%n3, &
      &        aux,lwork,iwork,liwork,info)

   !elow = 1.0e+99_wp
   elow = minval(self%eigv,mask=(abs(self%eigv) > thr1))
   !do i = 1, self%n3
   !   if (abs(self%eigv(i)) > thr1 ) elow = min(elow,self%eigv(i))
   !enddo

   damp = max(self%hlow - elow,0.0_wp)
   where(abs(self%eigv) > thr2) self%eigv = self%eigv + damp
!   do i = 1, self%n3
!      if (abs(self%eigv(i)) > thr2 ) self%eigv(i) = self%eigv(i) + damp
!   enddo

   if(pr)then
      write(stdout,*) 'Shifting diagonal of input Hessian by ', damp
      write(stdout,*) 'Lowest  eigenvalues of input Hessian'
      write(stdout,'(6F12.6)')(self%eigv(i),i=1,min(18,self%n3))
      write(stdout,*) 'Highest eigenvalues'
      write(stdout,'(6F12.6)')(self%eigv(i),i=self%n3-5,self%n3)
      write(stdout,*)
   endif

   fail = .true.
   get_anc: do itry = 1, maxtry
      self%B = 0.0_wp
      self%hess = 0.0_wp
      nvar = 0
      ! take largest (positive) first
      do i = self%n3, 1, -1
         if (abs(self%eigv(i)) > thr .and. nvar < self%nvar) then
            nvar = nvar+1
            self%B(:,nvar) = u(:,i)
            self%hess(nvar+nvar*(nvar-1)/2) = &
               min(max(self%eigv(i),self%hlow),self%hmax)
         endif
      enddo

      if (nvar.ne.self%nvar) then
         thr = thr * 0.1_wp
         cycle get_anc
      endif

      fail = .false.
      exit get_anc

   enddo get_anc

   if (fail) &
      call raise('E',"ANC generation failed!")

   call sort(self%n3,self%nvar,self%hess,self%B)

   self%coord = 0.0_wp
end subroutine generate_anc_packed

pure subroutine sort(nat3,nvar,hess,b)
   implicit none
   integer :: ii,k,j,m,i
   integer, intent(in)    :: nat3,nvar
   real(wp),intent(inout) :: hess(nvar*(nvar+1)/2)
   real(wp),intent(inout) :: b(nat3,nat3)
   real(wp) :: pp,sc1
   real(wp),allocatable   :: edum(:)
   allocate( edum(nvar), source = 0.0_wp )

   do k=1,nvar
      edum(k)=hess(k+k*(k-1)/2)
   enddo
!  sort
   do ii = 2, nvar
      i = ii - 1
      k = i
      pp= edum(i)
      do j = ii, nvar
         if (edum(j) .gt. pp) cycle
         k = j
         pp= edum(j)
      enddo
      if (k .eq. i) cycle
      edum(k) = edum(i)
      edum(i) = pp
      do m=1,nat3
         sc1=b(m,i)
         b(m,i)=b(m,k)
         b(m,k)=sc1
      enddo
   enddo

   do k=1,nvar
      hess(k+k*(k-1)/2)=edum(k)
   enddo

end subroutine sort

subroutine get_cartesian(self,xyz)
   use xtb_mctc_accuracy, only : wp
   implicit none
   class(tb_anc),intent(in) :: self
   integer :: m,i,j,k
   real(wp),intent(out) :: xyz (3,self%n)
   real(wp) :: dum

!  generate cartesian displacement vector
   xyz = self%xyz
   call dgemv('n',self%n3,self%nvar,1.0_wp,self%B,self%n3,self%coord,1,1.0_wp,xyz,1)

end subroutine get_cartesian

end module xtb_type_anc
