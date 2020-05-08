!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
      implicit none
      integer nel,nopen,ndim,ihomoa,ihomob
      real*8 focca(ndim), foccb(ndim)
      integer focc(ndim)
      integer i,na,nb,ihomo

      focc=0
      focca=0
      foccb=0
! even nel      
      if(mod(nel,2).eq.0)then
      ihomo=nel/2
      do i=1,ihomo 
         focc(i)=2
      enddo
      if(2*ihomo.ne.nel) then
         ihomo=ihomo+1
         focc(ihomo)=1
         if(nopen.eq.0)nopen=1
      endif
      if(nopen.gt.1)then
         do i=1,nopen/2
            focc(ihomo-i+1)=focc(ihomo-i+1)-1
            focc(ihomo+i)=focc(ihomo+i)+1
         enddo
      endif
! odd nel      
      else
      na=nel/2+(nopen-1)/2+1
      nb=nel/2-(nopen-1)/2
      do i=1,na             
         focc(i)=focc(i)+1
      enddo
      do i=1,nb             
         focc(i)=focc(i)+1
      enddo
      endif

      do i=1,ndim
         if(focc(i).eq.2)then
            focca(i)=1.0d0
            foccb(i)=1.0d0
         endif
         if(focc(i).eq.1)focca(i)=1.0d0
      enddo

      ihomoa=0
      ihomob=0
      do i=1,ndim
         if(focca(i).gt.0.99) ihomoa=i
         if(foccb(i).gt.0.99) ihomob=i
      enddo

      if(ihomoa.lt.1) stop 'internal error in occu'
      end

