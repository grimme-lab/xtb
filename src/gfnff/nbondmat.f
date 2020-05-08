      subroutine ff_bond(n,nb,pair)
      implicit none
      integer n,nb(20,n),pair(n*(n+1)/2)

      integer nn(n),list(20,n),jj,nlist(20,n),i2,k,iadr

      nn(1:n)=nb(20,1:n)

      pair=0
      list=0
      do i=1,n
         ni=nn(i)
         list(1:ni,i)=nb(1:ni,i)
      enddo

      nlist=list

c one bond, tag=1      
      tag=1
      call pairs(n,nn,list,pair,tag)

c determine up to 4 bonds in between      
      do d=1,3

c loop over atoms      
      do i=1,n
         ni=nn(i)
         newi=ni
c all neighbors of i         
         do ii=1,ni
            i1=list(ii,i)
            ni1=nb(20,i1)
c all neighbors of neighbors of i         
            do iii=1,ni1
               newatom=nb(iii,i1)
               da=.false.
               do j=1,newi
                  if(newatom.eq.list(j,i))da=.true.
               enddo
               if(.not.da)then
                 newi=newi+1
                 nlist(newi,i)=newatom           
               endif
            enddo
         enddo
         nnn(i)=newi
      enddo

      list=nlist
      nn  =nnn

c one bond more
      tag=tag+1
      call pairs(n,nn,list,pair,tag)

      enddo

      end


      subroutine pairs(n,nn,list,pair,tag)
      implicit none
      integer n,nn(n),list(20,n),tag
      integer i,j,k,ni,nj,ii,jj
      integer pair(n*(n+1)/2)
      logical dai,daj

      k=0
      do i=1,n
         ni=nn(i)
         do j=1,i     
            k=k+1 
            if(i.eq.j)cycle
            nj=nn(j)
            dai=.false.
            daj=.false.
            do ii=1,ni
               if(list(ii,i).eq.j)daj=.true.
            enddo
            do jj=1,nj
               if(list(jj,j).eq.i)dai=.true.
            enddo
            if(dai.and.daj.and.pair(k).eq.0) then
               pair(k)=tag
            endif
         enddo
      enddo

      end

