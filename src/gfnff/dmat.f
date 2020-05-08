cccccccccccccccccccccccccccccccccccccccccccccc
c density matrix
c C: MO coefficient
c X: scratch
c P  dmat
cccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dmat(ndim,focc,C,P)
      implicit none
      integer ndim
      real*8 focc(*)
      real*8 C(ndim,ndim)
      real*8 P(ndim,ndim)
      integer i,m
      real*8,allocatable ::Ptmp(:,:)
              
      allocate(Ptmp(ndim,ndim))                  
      do m=1,ndim  
         do i=1,ndim
            Ptmp(i,m)=C(i,m)*focc(m)
         enddo
      enddo
      call DGEMM('N','T',ndim,ndim,ndim,1.0d0,C,
     .                   ndim,Ptmp,ndim,0.0d0,P,ndim)
      deallocate(Ptmp)

      end

