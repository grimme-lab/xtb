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

!----------------------------------------------------------------------------------------------
! module lidep
!
! Contains some routines related to problems that can arise due to
! linear dependencies in the basis set, i.e., if any of the eigenvalues
! of the diagonalized overlap matrix approach zero.
! With these routines we can get cut-off eigenvalues of the overlap that are too small
! and generate a transformation matrix X. With this matrix X we
! can then transform our Hamiltionian to H'=X**(T)*H*X, which can then be diagonalized
! without need for the Overlap.
!
! For more information see Chapter 3.4.5 "Orthogonalization of the Basis" in Szabo/Ostlund
!
! P.Pracht, May 2019
!----------------------------------------------------------------------------------------------
module lidep
   use iso_fortran_env, only : wp => real64
   use mctc_la, only : sygvd,gemm,symm
   use setparam , only: lidethr
   implicit none

   !real(wp) :: lidethr    ! cut-off threshold for small overlap eigenvalues
   logical  :: orthog = .false.       ! global logical if the orthogonal basis should be used

   contains
!----------------------------------------------------------------------------------------------
!  subroutine cholesky
!
!  Try to perform a Cholesky-factorization A = U**T*U, where
!  A has to be a symmetric positive definite matrix and U is  
!  the normalized upper triangular matrix.
!  The Cholesky decomposition is used here ONLY to test if
!  A is positive definite. 
!  This check is usefull to pre-test if the overlap will
!  have near linear dependencies, and is faster than checking the
!  eigenvalues of the overlap.
!
!  INPUT  - ndim  dimension of the matrix A
!         - A     symmetric matrix with the dimension A(ndim,ndim)
!                 A is not modified and will be the same on In- and Output
!
!  OUTPUT - fail  logical to indicate whether the Cholesky-factorization
!                 was successful
!
!----------------------------------------------------------------------------------------------
subroutine cholesky(iunit,printl,ndim,A,fail)
    implicit none
  !---
    integer,intent(in)  :: iunit
    integer,intent(in)  :: ndim
    real(wp),intent(in) :: A(ndim,ndim)
    logical,intent(out) :: fail
    logical,intent(in)  :: printl
  !---
    integer :: ierr
    real(wp),allocatable  :: U(:,:)

    if(printl) &
    & write(iunit,'(2x,a)',advance='no')'Checking positiv definite overlap ...'

  !---
    allocate(U(ndim,ndim))
    fail=.false.
    U=A
    call dpotrf('U',ndim,U,ndim,ierr) !LAPACK routine  A = U**T*U
                                      !on input U=A, retruns U

    if(printl) then
      write(iunit,'(2x,a)')'done.'
      write(iunit,*)
    endif

    if(ierr /= 0) then
      fail=.true.
      if(printl)then
      write(iunit,'(a)') '**** WARNING ****'
      write(iunit,'(a)') ' Cholesky factorization of the overlap failed!'
      write(iunit,'(a)') ' It is possible that there are near linear dependecies in the basis.'
      write(iunit,'(a)') ' Therefore PEEQ will work with an orthogonal basis by using'
      write(iunit,'(a)') ' canonical orthogonalization.'
      write(iunit,'(a)') '**** WARNING ****'
      write(iunit,'(a)')
      endif
    endif
    deallocate(U)
    return
end subroutine cholesky

!----------------------------------------------------------------------------------------------
!  subroutine renorm
!
!  Re-normalize the matirx A according to the diagonal elements.
!
!  INPUT  - ndim  dimension of the matrix A
!         - A     symmetric matrix with the dimension A(ndim,ndim)
!
!  OUTPUT - 
!
!----------------------------------------------------------------------------------------------
subroutine renorm(ndim,A)
    implicit none
  !---
    integer,intent(in)  :: ndim
    real(wp),intent(inout) :: A(ndim,ndim)
  !---
    integer :: ierr,i,j,k,l

    !do i=1,ndim
    !   A(:,i)=A(:,i)/A(i,i)
    !enddo
    do i=1,ndim
      do j=1,i-1
        A(i,j)=A(i,j)/A(i,i)
        A(j,i)=A(i,j)
      enddo
    enddo
    return
end subroutine renorm


!----------------------------------------------------------------------------------------------
!  subroutine canorthog
!
!  In order to get the orthogonal basis we need a transformation
!  matrix X. In practice, this X can be modified to circumvent 
!  near linear dependencies. This can hoever only be done if the 
!  so-called 'canonical orthogonalization' is used, i.e.,
!  X is initially determined as X = U*s**(-1/2), where
!  U is a unitary matix and s**(-1/2) is the inverse square root
!  of the diagonalized overlap matrix S.
!  As a sanity check, X**T*S*X = 1 has to be valid.
!  See Szabo/Ostlund, equation 3.169 and following.
!  This routine determines X = U*s**(-1/2) from the overlap S.
!
!  INPUT  - ndim    dimension of the matrices S and X
!         - S       overlap matrix with the dimension S(ndim,ndim)
!
!  OUTPUT - X       transformation matrix X, with the dimension X(ndim,newdim) 
!         - newdim  new dimension of trafo X
!         - fail    logical to indicate whether the setup of matrix X
!                   was successful
!
!----------------------------------------------------------------------------------------------
subroutine canorthog(iunit,ndim,S,X,newdim,printl,fail)
    implicit none
   !---   
    integer,intent(in)   :: iunit
    integer,intent(in)   :: ndim
    real(wp),intent(in)  :: S(ndim,ndim)
    real(wp),intent(out) :: X(ndim,ndim)
    logical,intent(out)  :: fail
    integer,intent(out)  :: newdim 
    logical,intent(in)   :: printl
   !---
    integer :: i
    real(wp),allocatable :: U(:,:)
    real(wp),allocatable :: ssq(:)
    real(wp),allocatable :: sisq(:)
    real(wp),allocatable :: P(:,:),P2(:,:)
   !---
    allocate(U(ndim,ndim),ssq(ndim),sisq(ndim))
    fail=.false.
    X=0.0d0
   !--- diagonalize overlap and get eigenvectors U and eigenvalues s
    if(printl) &
    & write(iunit,'(2x,a)',advance='no')'Diagonalization of the Overlap   ...'
    call canorthog2(ndim,S,U,ssq,fail)
    if(fail)then
       if(printl) write(iunit,'(2x,a)')'failed.'
       return
    else
      if(printl) write(iunit,'(2x,a)')'done.'
    endif

   !--- restructure U and ssq
    call sorteigen(ndim,ssq,U)

   !--- cut-off small eigenvalues and build the inverse square root s**(-1/2)
    newdim=ndim
    call lidepcut(iunit,ndim,ssq,U,sisq,newdim,printl)
    !call prmat(6,e,nbf,1,"eigenvalues")
    !call prmat(6,U,ndim,ndim,"U")
   !--- sanity check stuff
    !newdim=ndim
    !do i=1,ndim
    !   sisq(i)=1.0d0 / sqrt(ssq(i))
    !enddo
 
   !--- get X and its new dimensions
    if(printl) &
    & write(iunit,'(2x,a)',advance='no')'Building transformation matrix X ...'
    call buildtrafoX(ndim,ssq,U,sisq,newdim,X)
    !call prmat(6,X,ndim,ndim,"X")
    if(printl) then
     write(iunit,'(2x,a)')'done.'
     write(iunit,'(a)')
    endif
   !--- more sanity check stuff
    !call prmat(6,X,ndim,ndim,"Trafo X")
    !allocate(P(ndim,ndim),P2(ndim,ndim))
    !call gemm('N','N',ndim,newdim,ndim,1.0d0,S,ndim,X,ndim,0.0d0,P,ndim)
    !call gemm('T','N',ndim,ndim,newdim,1.0d0,X,ndim,P,ndim,0.0d0,P2,ndim)
    !call prmat(6,P2,ndim,ndim,"X**T*S*X")
    !deallocate(P2,P)
    !stop

    deallocate(sisq,ssq,U)
    return
end subroutine canorthog

!----------------------------------------------------------------------------------------------
!  subroutine canorthog2
!
!  Diagonalizes the overlap matrix S and returns the unitary transformation
!  matrix U and the eigenvalues s, i.e., U**T*S*U = s is performed.
!  
!
!  INPUT  - ndim  dimension of the matrices S and X
!         - S     overlap matrix with the dimension S(ndim,ndim)
!
!  OUTPUT - X     transformation matrix U, with the dimension U(ndim,ndim) 
!         - ssq   eigenvalues s of matrix S
!         - fail  logical to indicate whether the setup of U and 
!                 was successful
!
!----------------------------------------------------------------------------------------------
subroutine canorthog2(ndim,S,U,ssq,fail)
    implicit none
   !---   
    integer,intent(in)   :: ndim
    real(wp),intent(in)  :: S(ndim,ndim)
    real(wp),intent(out) :: U(ndim,ndim)
    real(wp),intent(out) :: ssq(ndim)
    logical,intent(out)  :: fail
   !---
    integer :: ierr
    integer :: lwork
    real(wp),dimension(1) :: wdim
    real(wp),allocatable  :: work(:)

   !---
    fail=.false.
    U=S
    ssq=0.0d0
   !--- determine optimal work space
   !--- only calculates the optimal size of the 'work' array
    call dsyev('v','l',ndim,U,ndim,ssq,wdim,-1,ierr)
    if(ierr /= 0)then
       fail=.true.
       return
    endif
    lwork=int(wdim(1))
    allocate(work(lwork))

   !--- diagonalize overlap s = U**T*S*U
    call dsyev('v','l',ndim,U,ndim,ssq,work,lwork,ierr)
    if(ierr /= 0)then
       fail=.true.
       return
    endif
    !!write(iunit,*) ssq
    return
end subroutine canorthog2


!----------------------------------------------------------------------------------------------
!  subroutine sorteigen
!
!  We are allowed to order the eigenvalues in any way in the diagonal matrix s,
!  provided we order the columns of U in the same way.
!  Sort the overlap eigenvalues form highest to lowest and also
!  sort the eigenvectors accordingly. 
!  s is typically structured in ascending order from DSYEV, 
!  hence, this order only has to be reversed.
!  
!
!  INPUT     - ndim  dimension of eigenvalue array s and matrix U
!
!  IN/OUTPUT - s     eigenvalues of the overlap with the dimension s(ndim)
!            - U     eigenvectors of the overlap with the dimension U(ndim,ndim)
!
!----------------------------------------------------------------------------------------------
subroutine sorteigen(ndim,s,U)
    implicit none
   !---
    integer,intent(in)     :: ndim
    real(wp),intent(inout) :: s(ndim)
    real(wp),intent(inout) :: U(ndim,ndim)
   !---
    integer  :: i,j,k,l,imax
    real(wp) :: smax
    real(wp) :: Umax
   !---
    imax=ndim-1
    do i=1,imax
      smax=s(i)
      k=i
      l=i+1
      do j=l,ndim
         if(s(j).le.smax) cycle
         k=j
         smax=s(j)
      enddo
      Umax=s(i)
      s(i)=s(k)
      s(k)=Umax
      do j=1,ndim
         Umax=U(j,i)
         U(j,i)=U(j,k)
         U(j,k)=Umax
      enddo
    enddo
  !---
    return
end subroutine sorteigen



!----------------------------------------------------------------------------------------------
!  subroutine lidepcut
!
!  Cut-off small eigenvalues, which can be necessary to prevent linear dependencies
!  
!
!  INPUT     - ndim    dimension of eigenvalue array s and matrix U
!            - printl  print some output?
!
!  IN/OUTPUT - s     eigenvalues of the overlap with the dimension s(ndim)
!                    on OUTPUT this is s**(1/2)
!            - U     eigenvectors of the overlap with the dimension U(ndim,ndim)
!
!  OUTPUT    - sisq  array that contains inverse square root of the eigenvalues s**(-1/2)
!
!----------------------------------------------------------------------------------------------
subroutine lidepcut(iunit,ndim,s,U,sisq,newdim,printl)
    implicit none
   !---
    integer,intent(in)     :: iunit
    integer,intent(in)     :: ndim
    real(wp),intent(inout) :: s(ndim)
    real(wp),intent(inout) :: U(ndim,ndim)
    real(wp),intent(out)   :: sisq(ndim)
    integer,intent(out)    :: newdim
    logical,intent(in)     :: printl
   !--- 
    integer  :: i,j,imax
    real(wp) :: eigen
    real(wp) :: ehig,elow
   !---
    !lidethr=1.0d-4

    sisq=0.0d0
    newdim=ndim
   !--- min and max eigenvalues of the overlap
    ehig = s(1)
    elow = s(ndim)
   !--- cut eigenvalues and eigenvectors that are below the predefined threshold
    if(printl) &
    & write(iunit,'(2x,a)',advance='no')'Cutting off small eigenvalues    ...'
    do i=1,ndim
      eigen = s(i)
      if(eigen > lidethr) then
        eigen = sqrt(eigen)
        s(i) = eigen
        sisq(i) = 1.0d0 / eigen
      else
        s(i) = 0.0d0
        sisq(i) = 0.0d0
        U(1:ndim,i) = 0.0d0
        newdim=newdim-1
      endif
    enddo
    if(printl)then
    write(iunit,'(2x,a)')'done.'
   !--- some printout
    write(iunit,'(2x,a)')'Maximum eigenvalues of the overlap:'
    write(iunit,'(4x,a,f10.4)')'Largest eigenvalue              : ',ehig
    write(iunit,'(4x,a,f10.4)')'Smallest eigenvalue             : ',elow
    write(iunit,'(2x,a,e14.4)')'Eigenvalue cut-off threshold      : ',lidethr
    write(iunit,'(2x,a,i6)')   'Initial number of eigenvectors    : ',ndim
    write(iunit,'(2x,a,i6)')   'Removed eigenvectors              : ',ndim-newdim
    write(iunit,'(2x,a,i6)')   'Number of remaining eigenvectors  : ',newdim
    write(iunit,'(4x,a,f10.4)')'New smallest eigenvalue         : ',s(newdim)**2
    endif
   !---
    return
end subroutine lidepcut


!----------------------------------------------------------------------------------------------
!  subroutine buildtrafoX
!
!  Generate the transofrmation matrix X from U and s**(-1/2).
!  See Szabo/Ostlund Eq. 3.171
!  
!
!  INPUT     - ndim  dimension of eigenvalue array s and matrix U
!            - ssq   square-root of the overlap eigenvalues with the dimension s(ndim)
!            - U     eigenvectors of the overlap with the dimension U(ndim,ndim)
!            - sisq  inverse square root of the overlap eigenvalues s**(-1/2)
!
!  OUTPUT    - X     transformation matrix X=U*s**(-1/2)
!
!----------------------------------------------------------------------------------------------
subroutine buildtrafoX(ndim,ssq,U,sisq,newdim,X)
    implicit none
   !---
    integer,intent(in)     :: ndim
    integer,intent(in)     :: newdim
    real(wp),intent(in)    :: ssq(ndim)
    real(wp),intent(in)    :: U(ndim,ndim)
    real(wp),intent(in)    :: sisq(ndim)
    real(wp),intent(out)   :: X(ndim,ndim)
   !---
    integer :: i,j,k,l
   !---
    X=0.0d0
    do i=1,ndim
     do j=1,newdim
      X(i,j)=U(i,j) * sisq(j)
     enddo
    enddo
   !---
    return
end subroutine buildtrafoX



!----------------------------------------------------------------------------------------------
!  subroutine orthgsolve
!
!  Modification of the solv routine.
!  Diagonalization of H is done with the trafo matrix X = U*s**(-1/2), i.e.,
!  in the orthogonal basis.
!  See Szabo/Ostlund equation 3.176:
!  (X**T*F*X)*C' = (X**T*S*X)*C'*e
!  which is possible since X**T*S*X=1
!
!  INPUT  - full    some logical to indicate was should be done in this routine
!         - ndim    original dimension of the matrices H,S,X,P and vector e
!         - cutdim  new second dimension of the rectangular trafo matrix X
!         - ihomo 
!         - acc     accuracy setting
!         - H       Hamiltonian matrix
!         - X       trafo matrix with dimensions X(ndim,cutdim)
!         - S       overlap matrix with the dimension S(ndim,ndim)
!         - e       eigenvalues
!         
!
!  OUTPUT - fail    logical to indicate whether the diagonalization
!                   was successful
!         - 
!
!----------------------------------------------------------------------------------------------
subroutine orthgsolve(full,ndim,cutdim,ihomo,acc,H,S,X,P,e,fail)
   implicit none
   integer, intent(in)   :: ndim
   integer, intent(in)   :: cutdim
   logical, intent(in)   :: full
   integer, intent(in)   :: ihomo
   real(wp),intent(inout):: H(ndim,ndim)
   real(wp),intent(in)   :: S(ndim,ndim)
   real(wp),intent(inout):: X(ndim,ndim)
   real(wp),intent(out)  :: P(ndim,ndim)
   real(wp),intent(out)  :: e(ndim)
   real(wp),intent(in)   :: acc
   logical, intent(out)  :: fail

   integer :: i,j,info,lwork,liwork,nfound,iu
   integer :: nbf,xbf
   integer, allocatable :: iwork(:),ifail(:)
   real(wp),allocatable :: aux  (:)
   real(wp),allocatable :: Haux(:,:),Paux(:,:)
   fail =.false.


  !--- DIAG IN ORTHOGONAL BASIS WITH X=U*s^-1/2 TRAFO
    nbf = ndim
    xbf = cutdim
    lwork  = 1 + 6*nbf + 2*nbf**2
    allocate (aux(lwork))
  !--- calculate F'=X**T*F*X  

    call dgemm('n','n',nbf,xbf,nbf,1.0d0,H,nbf,X,nbf,0.0d0,P,nbf)
    call dgemm('t','n',xbf,xbf,nbf,1.0d0,X,nbf,P,nbf,0.0d0,H,xbf)
    
  !--- Caluclate C' and ε from F'*C' = C'*ε
    call dsyev('v','u',xbf,H,xbf,e,aux,lwork,info)
    if(info.ne.0)then
       fail=.true.
       return
    endif
  !--- calculate C = X*C'
    P=0.0d0
    call dgemm('N','N',nbf,xbf,xbf,1.0d0,X,nbf,H,xbf,0.0d0,P,nbf)
    H=P
    deallocate(aux)
    return
end subroutine orthgsolve


subroutine orthgsolve2(full,ndim,cutdim,ihomo,acc,H,S,X,P,e,fail)
   implicit none
   integer, intent(in)   :: ndim
   integer, intent(in)   :: cutdim
   logical, intent(in)   :: full
   integer, intent(in)   :: ihomo
   real(wp),intent(inout):: H(ndim,ndim)
   real(wp),intent(in)   :: S(ndim,ndim)
   real(wp),intent(inout):: X(ndim,ndim)
   real(wp),intent(out)  :: P(ndim,ndim)
   real(wp),intent(out)  :: e(ndim)
   real(wp),intent(in)   :: acc
   logical, intent(out)  :: fail

   integer :: i,j,info,lwork,liwork,nfound,iu
   integer :: nbf,xbf
   integer, allocatable :: iwork(:),ifail(:)
   real(wp),allocatable :: aux  (:)
   real(wp),allocatable :: Haux(:,:),Paux(:,:),Xaux(:,:)
   fail =.false.


  !--- DIAG IN ORTHOGONAL BASIS WITH X=U*s^-1/2 TRAFO
    nbf = ndim
    xbf = cutdim
    lwork  = 1 + 6*nbf + 2*nbf**2
    allocate (aux(lwork))
  !--- calculate F'=X**T*F*X  
    allocate(Xaux(nbf,xbf),Paux(nbf,xbf))
    Xaux=X(:,1:xbf)
    call dgemm('n','n',nbf,xbf,nbf,1.0d0,H,nbf,Xaux,nbf,0.0d0,Paux,nbf)
    allocate(Haux(xbf,xbf))
    call dgemm('t','n',xbf,xbf,nbf,1.0d0,Xaux,nbf,Paux,nbf,0.0d0,Haux,xbf)

  !--- Caluclate C' and ε from F'*C' = C'*ε
    call dsyev('v','u',xbf,Haux,xbf,e,aux,lwork,info)
    !call prmat(6,e,nbf,1,"eigenvalues")

    if(info.ne.0)then
       fail=.true.
       return
    endif
  !--- calculate C = X*C'
    P=0.0d0
    call dgemm('N','N',nbf,xbf,xbf,1.0d0,Xaux,nbf,Haux,xbf,0.0d0,P,nbf)
    H=P

    deallocate(Haux,Paux,Xaux)
    deallocate(aux)
    return
end subroutine orthgsolve2


!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
end module lidep
