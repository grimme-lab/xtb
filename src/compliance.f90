module xtb_compliance
   use xtb_mctc_accuracy, only : wp
   !---------------------------------------------------------------------------
   ! C = F^{-1} in internal coordinates.
   ! Works for any nint (diatomic, linear, mixed, general).
   ! nint passed explicitly — no hardcoded 3N-6.
   !---------------------------------------------------------------------------
   implicit none
   private
   public :: compute_compliance

contains

   subroutine compute_compliance(unit, H, B, natoms, nint, C, stat)
      !----------------------------------------------------------------------
      ! F = G^{-1} B H B^T G^{-1},   C = F^{-1}
      !
      ! G = B B^T may be singular when nint > rank(B), e.g. for molecules
      ! with symmetric tops (CH3, NH2) or linear functional groups where
      ! the Z-matrix coordinate set is slightly redundant.  In that case
      ! the Moore-Penrose pseudoinverse G^+ is used instead of G^{-1},
      ! computed via SVD with threshold tol = eps_svd * sigma_max.
      ! The compliance matrix C is then the pseudoinverse of F.
      ! Rows/columns of C corresponding to zero singular values of F are
      ! set to zero (they describe redundant, unphysical coordinates).
      !----------------------------------------------------------------------
      integer, intent(in)  :: unit
      integer,  intent(in)  :: natoms, nint
      real(wp),  intent(in)  :: H(3*natoms, 3*natoms)
      real(wp),  intent(in)  :: B(nint, 3*natoms)
      real(wp),  intent(out) :: C(nint, nint)
      integer,  intent(out) :: stat

      integer  :: ndim, lwork, info, i, rank_g, rank_f
      real(wp)  :: tol_g, tol_f
      real(wp), parameter :: eps_svd = 1.0d-10
      real(wp), allocatable :: G(:,:), Ginv(:,:), BH(:,:), BHBt(:,:), F(:,:)
      real(wp), allocatable :: work(:), U(:,:), VT(:,:), S(:), Utmp(:,:)

      stat=0; ndim=3*natoms
      allocate(G(nint,nint), Ginv(nint,nint), BH(nint,ndim), &
               BHBt(nint,nint), F(nint,nint))
      allocate(U(nint,nint), VT(nint,nint), S(nint))

      !--- G = B B^T ---
      call dgemm('N','T',nint,nint,ndim,1d0,B,nint,B,nint,0d0,G,nint)

      !--- Ginv = G^+ via SVD (handles redundant coordinate sets) ---
      Ginv = G
      lwork = -1; allocate(work(1))
      call dgesvd('A','A',nint,nint,Ginv,nint,S,U,nint,VT,nint,work,lwork,info)
      lwork = int(work(1)); deallocate(work); allocate(work(lwork))
      call dgesvd('A','A',nint,nint,Ginv,nint,S,U,nint,VT,nint,work,lwork,info)
      if (info/=0) then
         write(unit, '(A,I0)') 'compute_compliance: DGESVD(G) info=',info
         stat=info; return
      end if
      tol_g = eps_svd * S(1)
      rank_g = count(S > tol_g)
      if (rank_g < nint) &
         write(unit, '(A,I0,A,I0)') &
            '  Note: G matrix rank=', rank_g, ' < nint=', nint
      ! Ginv = V * S^{-1} * U^T  (only for singular values > tol_g)
      Ginv = 0d0
      allocate(Utmp(nint,nint), source=0.0_wp)
      do i = 1, nint
         if (S(i) <= tol_g) cycle
         Utmp(:,i) = U(:,i) / S(i)
      end do
      call dgemm('T','T',nint,nint,nint,1d0,VT,nint,Utmp,nint,0d0,Ginv,nint)
      deallocate(Utmp)
      call symmetrise(Ginv, nint)

      !--- BH = B H,  BHBt = B H B^T ---
      call dgemm('N','N',nint,ndim,ndim,1d0,B,nint,H,ndim,0d0,BH,nint)
      call dgemm('N','T',nint,nint,ndim,1d0,BH,nint,B,nint,0d0,BHBt,nint)

      !--- F = G^+ B H B^T G^+ ---
      call dgemm('N','N',nint,nint,nint,1d0,Ginv,nint,BHBt,nint,0d0,BH(1:nint,1:nint),nint)
      call dgemm('N','N',nint,nint,nint,1d0,BH(1:nint,1:nint),nint,Ginv,nint,0d0,F,nint)
      call symmetrise(F,nint)

      !--- C = F^+ via SVD ---
      S = 0d0
      C = F
      call dgesvd('A','A',nint,nint,C,nint,S,U,nint,VT,nint,work,lwork,info)
      if (info/=0) then
         write(unit, '(A,I0)') 'compute_compliance: DGESVD(F) info=',info
         stat=info; return
      end if
      tol_f = eps_svd * S(1)
      rank_f = count(S > tol_f)
      if (rank_f < nint) &
         write(unit, '(A,I0,A,I0)') &
            '  Note: F matrix rank=', rank_f, ' < nint=', nint
      C = 0d0
      allocate(Utmp(nint,nint), source=0.0_wp)
      do i = 1, nint
         if (S(i) <= tol_f) cycle
         Utmp(:,i) = U(:,i) / S(i)
      end do
      call dgemm('T','T',nint,nint,nint,1d0,VT,nint,Utmp,nint,0d0,C,nint)
      deallocate(Utmp)
      call symmetrise(C, nint)

      deallocate(G,Ginv,BH,BHBt,F,U,VT,S,work)

   end subroutine compute_compliance

   subroutine symmetrise(A,n)
      integer, intent(in)    :: n
      real(wp), intent(inout) :: A(n,n)
      integer :: i,j
      do i=1,n; do j=i+1,n; A(j,i)=A(i,j); end do; end do
   end subroutine symmetrise

end module xtb_compliance
