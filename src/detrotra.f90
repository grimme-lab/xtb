! try to determine rot/tra vibrational modes
! even if contaminated from real*4 diag
! SG 5/20

subroutine detrotra4(linear,n,at,xyz,h,eig)
      use xtb_mctc_accuracy, only : sp, wp
      implicit none
      integer, intent(in)     :: n
      integer, intent(in)     :: at(n)
      real(wp), intent(in)    :: xyz(3,n)
      real(sp), intent(in)    :: h(3*n,3*n)    ! values from projected Lindh diag
      real(sp), intent(inout) :: eig(3*n)      ! eigenvectors from projected Lindh diag
      logical, intent(in)     :: linear

      integer               :: i,j,k,kk,ii,nn,n3,nend
      integer,allocatable   :: ind(:)
      real(wp), allocatable :: tmp(:,:)
      real(wp), allocatable :: e(:)
      real(wp)              :: a0,b0,c0

      allocate(tmp(3,n),e(3*n),ind(3*n))

      n3 = 3*n
      nn = 0
      do ii=1, n3
         if(eig(ii).gt.0.05) cycle  ! only lowest checked

         kk=0
         do j=1,n
            do k=1,3
               kk=kk+1
               tmp(k,j)=xyz(k,j)+h(kk,ii) ! distort along mode ii
            enddo
         enddo

         c0=0            ! compared all interatomic distances of original and distortet geom.
         do i=2,n
            do j=1,i-1
               a0 = sqrt((xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2)
               b0 = sqrt((tmp(1,i)-tmp(1,j))**2+(tmp(2,i)-tmp(2,j))**2+(tmp(3,i)-tmp(3,j))**2)
               c0 = c0+(a0-b0)**2   
            enddo
         enddo
         nn = nn + 1
         e(nn)=sqrt(c0/n)*abs(eig(ii))  ! weight by Lindh eigenvalue
         ind(nn)=nn

      enddo

      call qsort(e, 1, nn, ind)   ! sort

      nend = 6
      if (linear) nend = 5

      do i=1,nend
!        write(*,'(2i3,e16.6,6f12.4)') i,ind(i),eig(ind(i)),e(i)   
         eig(ind(i)) = 0.0  ! identifier for rot/tra
      enddo

!     eig = eig * 100000.
!     call g98fakex('lindh.out',n,at,xyz,eig,h)

      end subroutine detrotra4
