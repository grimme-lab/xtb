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

module xtb_detrotra
  implicit none
  private
  public :: detrotra4, detrotra8

  contains

  subroutine detrotra4(linear,mol,h,eig)
    use xtb_mctc_accuracy, only : sp, wp
    use xtb_type_molecule
    implicit none
    type(TMolecule), intent(in) :: mol
    real(sp), intent(in)    :: h(3*mol%n,3*mol%n)    ! values from projected Lindh diag
    real(sp), intent(inout) :: eig(3*mol%n)          ! eigenvectors from projected Lindh diag
    logical, intent(in)     :: linear
  
    integer               :: i,j,k,kk,ii,nn,n3,nend
    integer,allocatable   :: ind(:)
    real(wp), allocatable :: tmp(:,:)
    real(wp), allocatable :: e(:)
    real(wp)              :: a0,b0,c0
  
    n3 = 3*mol%n
  
    allocate(tmp(3,mol%n),e(n3),ind(n3))
  
    nn = 0
    do ii=1, n3
       if(eig(ii).gt.0.05) cycle  ! only lowest checked
  
       kk=0
       do j=1,mol%n
          do k=1,3
             kk=kk+1
             tmp(k,j)=mol%xyz(k,j)+h(kk,ii) ! distort along mode ii
          enddo
       enddo
  
       c0=0            ! compared all interatomic distances of original and distortet geom.
       do i=2,mol%n
          do j=1,i-1
             a0 = sqrt((mol%xyz(1,i)-mol%xyz(1,j))**2+(mol%xyz(2,i)-mol%xyz(2,j))**2+(mol%xyz(3,i)-mol%xyz(3,j))**2)
             b0 = sqrt((tmp(1,i)-tmp(1,j))**2+(tmp(2,i)-tmp(2,j))**2+(tmp(3,i)-tmp(3,j))**2)
             c0 = c0+(a0-b0)**2   
          enddo
       enddo
       nn = nn + 1
       e(nn)=sqrt(c0/mol%n)*abs(eig(ii))  ! weight by Lindh eigenvalue
       ind(nn)=nn
  
    enddo
  
    call qsort(e, 1, nn, ind)   ! sort
  
    nend = 6
    if (linear) nend = 5
  
    do i=1,nend
       eig(ind(i)) = 0.0  ! identifier for rot/tra
    enddo
  
  end subroutine detrotra4

!> determine rotational and translational modes 
subroutine detrotra8(linear,n,xyz,h,eig)
   
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule
   implicit none
    
   integer, intent(in)     :: n
   real(wp), intent(in)    :: xyz(3,n)    ! values from projected Lindh diag
   real(wp), intent(in)    :: h(3*n,3*n)    ! values from projected Lindh diag
   real(wp), intent(inout) :: eig(3*n)          ! eigenvalues from projected Lindh diag
   logical, intent(in)     :: linear
  
   integer               :: i,j,k,kk,ii,nn,n3,nend
   integer,allocatable   :: ind(:)
   real(wp), allocatable :: tmpxyz(:,:)
   real(wp), allocatable :: e(:)
   real(wp)              :: a0,b0,c0

   n3 = 3*n

   allocate(tmpxyz(3,n),e(n3),ind(n3))

   nn = 0
   do ii=1, n3

      if(eig(ii).gt.0.05) cycle  ! check only low-lying modes 
      
      ! distort initial geometry along ii-th mode !
      kk=0
      do j=1,n
         do k=1,3
            kk=kk+1
            tmpxyz(k,j)=xyz(k,j)+h(kk,ii)
         enddo
      enddo

      ! compare all interatomic distances of original and distortet geom. !
      c0=0           
      do i=2,n
         do j=1,i-1
            a0 = sqrt((xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2)
            b0 = sqrt((tmpxyz(1,i)-tmpxyz(1,j))**2+(tmpxyz(2,i)-tmpxyz(2,j))**2+(tmpxyz(3,i)-tmpxyz(3,j))**2)
            c0 = c0+(a0-b0)**2 ! sum of squared differences  
         enddo
      enddo

      nn = nn + 1 
      e(nn)=sqrt(c0/n)*abs(eig(ii))  ! weight by Lindh eigenvalue
      ind(nn)=nn

   enddo

   ! sort in ascending order !
   call qsort(e,1,nn,ind)   

   nend = 6
   if (linear) nend = 5

   ! set lowest eigenvalues to zero !
   do i=1,nend
      eig(ind(i)) = 0.0  
   enddo

end subroutine detrotra8
  
end module xtb_detrotra
