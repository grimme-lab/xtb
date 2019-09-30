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

c--------1---------2---------3---------4---------5---------6---------7--
c reference : g. cisneros, m .berrondo, c. f. bunge
c             computers and chemistry vol. 4, p 281 (1986)
c for ibm rs6000:
c                 s. grimme, uni bonn, 9.92
c erweitert auf gleichzeitige behandlung mehrerer wurzeln
c nach einem verfahren von
c       b. liu
c       numerical algorithms in chemistry:algebraic methods,
c       lbl-8158 lawrence berkeley laboratory
c       eds.: c. moler and i. shavitt
c                 b. engels, uni bonn, 12.93
c direct version and several modifications s.grimme, 6.2.99
c new general (BLAS) version SG, 12/18
c note: the current version only works for the LOWEST root i.e. nr=1 !!!
c-----------------------------------------------------------------------
 
      subroutine ddavid(pr,ini,n,nr,crite,H,C,e)
      implicit none                        
      logical pr      ! print logical
      logical ini     ! initialize start vector if .true., if false starts from previous vector
      integer n       ! dimension
      integer nr      ! # roots
      real*8  crite   ! eigenvalue convergence threshold
      real*8  H(n,n)  ! matrix to be diagonalized
      real*8  C(n,nr) ! eigenvectors
      real*8  e(nr)   ! eigenvalues

c local
      integer maxiter        ! maximum # of iterations
      parameter (maxiter=100)      
      integer iter,iconf(nr),ineue(nr),janf,jneu1,lun1,lun2
      integer iideks(nr*maxiter),mx2,idum,j,jalt,ilauf,jneu,nneue
      integer lauf,l1,l2,k,LWORK,LIWORK,INFO,ianf,i,ien,ico,ialt,memlun2
      integer,allocatable :: iwork(:)

      real*8 valn(nr),uim,s,ddot,one,zero,denerg
      real*8, allocatable :: adiag(:),vecf1(:),vecf2(:),w(:)
      real*8, allocatable :: Uaug(:,:),d(:),aux(:)
      real*8, allocatable :: AB(:,:),av(:),tmpav(:,:)
      real*8, allocatable :: HP(:)
      parameter (one =1)      
      parameter (zero=0)      

      mx2=nr*maxiter

      if(pr)then
      write(*,'(/,10x,''******************************************'')')
      write(*,'(10x,''*            multi-root davidson (R8)    *'')')
      write(*,'(10x,''******************************************'',/)')
      write(*,*) 'dim ',n,' # roots ',nr
      endif

      allocate(adiag(n),vecf1(n*nr),vecf2(n*nr),w(n),av(mx2*(mx2+1)/2),
     .         HP(n*(n+1)/2))

c IO
      lun1=86
      lun2=87
      open(unit=lun1,recl=n*8,access='direct',
     .        file='david1.tmp',
     .        form='unformatted')    
      open(unit=lun2,recl=n*8,access='direct',
     .        file='david2.tmp',
     .        form='unformatted')    

      k=0
      do i=1,n
         do j=1,i
            k=k+1
            HP(k)=H(j,i)
         enddo
      enddo

C H * C for initialization
      ianf = 1
      do i = 1,nr
         if(ini) then
          do k=1,n
            call random_number(s)
            C(k,i)=s
          enddo
          s=ddot(n,C(1,i),1,C(1,i),1)
          C(1:n,i)=C(1:n,i)/sqrt(s)
         endif
         call dmwrite(n,lun1,C(1,i),i)
c        call dgemv('N',n,n,ONE,H,n,C(1,i),1,ZERO,vecf2(ianf),1)
         call dspmv('U',n, ONE,HP,  C(1,i),1,ZERO,vecf2(ianf),1)
         call dmwrite(n,lun2,vecf2(ianf),i)   
         ianf = ianf + n
      enddo

c aufbau des iideks feldes
      iideks(1) = 1
      do idum = 2,mx2
         iideks(idum) = iideks(idum - 1) + idum
      enddo
      valn = 0
      iconf= 0
      e    = 0

      do i=1,n
         adiag(i)=H(i,i)
      enddo
    
      if(nr.eq.1)then
         av(1)=ddot(n,C(1,1),1,vecf2(1),1)
      else
         stop '# roots > 1 not implemented'
c aufbau der startmatrix av = bi*a*bj mit bi, bj startvektoren
c     allocate(tmpav(nr,nr),AB(n,nr))
c     k=0
c     do l1=1,nr
c        do l2=1,l1
c           k=k+1
c           av(k)=tmpav(l2,l1)
c        enddo
c     enddo
c     deallocate(tmpav,AB)
      endif
c done

      j = nr 
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c eigentliche schleife im davidson
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do iter = 1, maxiter

      lwork  = 1 + 6*j + 2*j**2
      liwork = 8*j
      allocate(Uaug(j,j),d(j),iwork(liwork),aux(lwork))

      k=0
      do l1=1,j
         do l2=1,l1
            k=k+1
            Uaug(l2,l1) = av(k)
            Uaug(l1,l2) = av(k)
         enddo
      enddo
      call dsyevd('V','U',j,Uaug,j,d,aux,LWORK,IWORK,LIWORK,INFO)
      valn(1:nr) = d(1:nr)

c aufbau der eigentlichen vektoren, die stehen dann auf vecf1
c multiplikation der vorherigen entwicklungsvektoren mit den entwicklungskoeffizienten
      vecf1=0     
      do i=1,j
         call dmread(n,lun1,w,i)
         ianf = 1
         do lauf = 1,nr    
            uim = Uaug(i,lauf)
            call daxpy(n,uim,w,1,vecf1(ianf),1)
            ianf = ianf + n
         enddo
      enddo

c aufbau -E*bi: vecf1 ist bi; vecf2 dann E*bi
      ianf = 1
      do lauf = 1,nr    
         ien = ianf + n -1
         vecf2(ianf:ien) = -valn(lauf) * vecf1(ianf:ien)
         ianf = ianf + n
      enddo
c aufbau des residuen vektors (h*bi-e*bi)
c steht dann auf vecf2 ueberschreibt also -e*bi
c h*bi-e*bi steht jetzt auf vecf2
      do i=1,j
         call dmread(n,lun2,w,i)
         memlun2 = i
         ianf = 1
         do lauf = 1,nr   
            uim = Uaug(i,lauf)
            call daxpy(n,uim,w,1,vecf2(ianf),1)
            ianf = ianf + n
         enddo
      enddo
      deallocate(Uaug,d,iwork,aux)

      ianf = 1
      do lauf = 1,nr            
         C(1:n,lauf)=vecf1(ianf:ianf+n-1)
         ianf = ianf + n
      enddo

c aufbau (h*bi - e*bi)/(e - haa); steht danach auf vecf2
      ianf = 1
      do lauf = 1,nr   
         call dvdssvd(n,valn(lauf),vecf2(ianf),adiag,vecf1(ianf))
         ianf = ianf + n
      enddo

      denerg = 0.0d0
      do lauf = 1,nr
         denerg = denerg + abs(valn(lauf) - e(lauf))
         if ( abs(valn(lauf) - e(lauf)) .lt. crite) iconf(lauf) = 1
      enddo
      denerg = denerg / nr   
      ico = sum(iconf)
      if(pr) write(*,*) iter,ico,denerg,valn(1:nr)

      if (ico .eq. nr) then
          if(pr) write(*,*) 'all roots converged'
          go to 9999
      endif

      if(j.gt.0) then
c-- mit allen alten
         nneue = 0
         ialt = j
         ianf = 1
         do lauf = 1,nr    
c orthogonaliesung des betrachteten auf die alten
           do jalt = 1,ialt
              call dmread(n,lun1,w,jalt)
              s=-ddot(n,w,1,vecf1(ianf),1)
              call daxpy(n,s,w,1,vecf1(ianf),1)
           enddo
c normierung dessen was vom betrachteten uebrig bleibt
           s=ddot(n,vecf1(ianf),1,vecf1(ianf),1)
           if (s.gt.0.00000001)then
c neuer wird mitgenommen
              s = ONE /sqrt(s)
              vecf1(ianf:ianf+n-1)= vecf1(ianf:ianf+n-1) * s
              ialt = ialt + 1
              nneue = nneue + 1
              ineue(nneue) = ianf
c wegschreiben des neuen zu den alten
              call dmwrite(n,lun1,vecf1(ianf),jalt)
           else
              goto 9999
           endif
           ianf = ianf + n
         enddo
      endif

c umspeichern der ueberlebenden vektoren auf vecf1
      ianf = 1
      do lauf = 1,nneue
         call dcopy(n,vecf1(ineue(lauf)),1,vecf1(ianf),1)
         ianf = ianf + n
      enddo

C H * C
      ianf = 1
      do i = 1,nneue
c        call dgemv('N',n,n,ONE,H,n,vecf1(ianf),1,ZERO,vecf2(ianf),1)
         call dspmv('U',n, ONE,HP,  vecf1(ianf),1,ZERO,vecf2(ianf),1)
         ianf = ianf + n
      enddo

      ianf = 1
      do i = 1,nneue
         call dmwrite(n,lun2,vecf2(ianf),memlun2+i)   
         ianf = ianf + n
      enddo

c berechnung der neuen matrixelemente der davidson-matrix
c zunaechst mit den alten
      do jalt = 1,j 
        call dmread(n,lun1,w,jalt)
        ianf = 1 
        ilauf = iideks(j) + jalt
        do jneu = 1,nneue 
           av(ilauf) = ddot(n,w,1,vecf2(ianf),1)
           ilauf = ilauf + jneu + j 
           ianf = ianf + n 
        enddo  
      enddo
c dann mit den neuen
      ianf = 1 
      do jneu = 1,nneue
         janf = 1 
         ilauf = iideks(j+jneu) - jneu + 1 
         do jneu1 = 1,jneu
            av(ilauf) = ddot(n,vecf2(ianf),1,vecf1(janf),1)
            janf = janf + n 
            ilauf = ilauf + 1 
        enddo
        ianf = ianf + n 
      enddo

c increase expansion space and iterate further
      e = valn 
      j = j + nneue

      enddo   

      if(pr) write(*,*) 'Warning: davidson not properly converged'
c exit

9999  continue

      deallocate(adiag,vecf1,vecf2,w,av,HP)
      close (lun1,status='delete')
      close (lun2,status='delete')

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dvdssvd(n,sl,v1,v2,v)
      implicit none                          
      real*8 sl, v1(*), v(*), v2(*)
      integer n

      integer i,m,mp1

      do i = 1,n
        v(i)=v1(i)/(sl-v2(i))
      enddo

      return
      end 

***********************************************************************

      subroutine dmwrite(n,iwo,v,irec)
      implicit none
      real*8 v(n)
      integer n,iwo,irec
      write(iwo,rec=irec) v
      return
      end 


      subroutine dmread(n,iwo,v,irec)
      implicit none
      real*8 v(n)
      integer n,iwo,irec
      read(iwo,rec=irec) v
      return
      end 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c real*4 version
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sdavid(pr,ini,n,nr,crite,H,C,e)
      implicit none                        
      logical pr      ! print logical
      logical ini     ! initialize start vector if .true., if false starts from previous vector
      integer n       ! dimension
      integer nr      ! # roots
      real*4  crite   ! eigenvalue convergence threshold
      real*4  H(n,n)  ! matrix to be diagonalized
      real*4  C(n,nr) ! eigenvectors
      real*4  e(nr)   ! eigenvalues

c local
      integer maxiter        ! maximum # of iterations
      parameter (maxiter=100)      
      integer iter,iconf(nr),ineue(nr),janf,jneu1,lun1,lun2
      integer iideks(nr*maxiter),mx2,idum,j,jalt,ilauf,jneu,nneue
      integer lauf,l1,l2,k,LWORK,LIWORK,INFO,ianf,i,ien,ico,ialt,memlun2
      integer,allocatable :: iwork(:)

      real*4 valn(nr),uim,s,sdot,one,zero,denerg
      real*4, allocatable :: adiag(:),vecf1(:),vecf2(:),w(:)
      real*4, allocatable :: Uaug(:,:),d(:),aux(:)
      real*4, allocatable :: AB(:,:),av(:),tmpav(:,:)
      real*4, allocatable :: HP(:)
      parameter (one =1)      
      parameter (zero=0)      

      mx2=nr*maxiter

      if(pr)then
      write(*,'(/,10x,''******************************************'')')
      write(*,'(10x,''*            multi-root davidson (R4)    *'')')
      write(*,'(10x,''******************************************'',/)')
      write(*,*) 'dim ',n,' # roots ',nr
      endif

      allocate(adiag(n),vecf1(n*nr),vecf2(n*nr),w(n),av(mx2*(mx2+1)/2),
     .         HP(n*(n+1)/2))

c IO
      lun1=86
      lun2=87
      open(unit=lun1,recl=n*4,access='direct',
     .        file='david1.tmp',
     .        form='unformatted')    
      open(unit=lun2,recl=n*4,access='direct',
     .        file='david2.tmp',
     .        form='unformatted')    

      k=0
      do i=1,n
         do j=1,i
            k=k+1
            HP(k)=H(j,i)
         enddo
      enddo

C H * C for initialization
      ianf = 1
      do i = 1,nr
         if(ini) then
          do k=1,n
            call random_number(s)
            C(k,i)=s
          enddo
          s=sdot(n,C(1,i),1,C(1,i),1)
          C(1:n,i)=C(1:n,i)/sqrt(s)
         endif
         call smwrite(n,lun1,C(1,i),i)
c        call sgemv('N',n,n,ONE,H,n,C(1,i),1,ZERO,vecf2(ianf),1)
         call sspmv('U',n,  ONE,HP, C(1,i),1,ZERO,vecf2(ianf),1)
         call smwrite(n,lun2,vecf2(ianf),i)   
         ianf = ianf + n
      enddo

c aufbau des iideks feldes
      iideks(1) = 1
      do idum = 2,mx2
         iideks(idum) = iideks(idum - 1) + idum
      enddo
      valn = 0
      iconf= 0
      e    = 0

      do i=1,n
         adiag(i)=H(i,i)
      enddo
    
      if(nr.eq.1)then
         av(1)=sdot(n,C(1,1),1,vecf2(1),1)
      else
         stop '# roots > 1 not implemented'
c aufbau der startmatrix av = bi*a*bj mit bi, bj startvektoren
c     allocate(tmpav(nr,nr),AB(n,nr))
c     k=0
c     do l1=1,nr
c        do l2=1,l1
c           k=k+1
c           av(k)=tmpav(l2,l1)
c        enddo
c     enddo
c     deallocate(tmpav,AB)
      endif
c done

      j = nr 
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c eigentliche schleife im davidson
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do iter = 1, maxiter

      lwork  = 1 + 6*j + 2*j**2
      liwork = 8*j
      allocate(Uaug(j,j),d(j),iwork(liwork),aux(lwork))

      k=0
      do l1=1,j
         do l2=1,l1
            k=k+1
            Uaug(l2,l1) = av(k)
            Uaug(l1,l2) = av(k)
         enddo
      enddo
      call ssyevd('V','U',j,Uaug,j,d,aux,LWORK,IWORK,LIWORK,INFO)
      valn(1:nr) = d(1:nr)

c aufbau der eigentlichen vektoren, die stehen dann auf vecf1
c multiplikation der vorherigen entwicklungsvektoren mit den entwicklungskoeffizienten
      vecf1=0     
      do i=1,j
         call smread(n,lun1,w,i)
         ianf = 1
         do lauf = 1,nr    
            uim = Uaug(i,lauf)
            call saxpy(n,uim,w,1,vecf1(ianf),1)
            ianf = ianf + n
         enddo
      enddo

c aufbau -E*bi: vecf1 ist bi; vecf2 dann E*bi
      ianf = 1
      do lauf = 1,nr    
         ien = ianf + n -1
         vecf2(ianf:ien) = -valn(lauf) * vecf1(ianf:ien)
         ianf = ianf + n
      enddo
c aufbau des residuen vektors (h*bi-e*bi)
c steht dann auf vecf2 ueberschreibt also -e*bi
c h*bi-e*bi steht jetzt auf vecf2
      do i=1,j
         call smread(n,lun2,w,i)
         memlun2 = i
         ianf = 1
         do lauf = 1,nr   
            uim = Uaug(i,lauf)
            call saxpy(n,uim,w,1,vecf2(ianf),1)
            ianf = ianf + n
         enddo
      enddo
      deallocate(Uaug,d,iwork,aux)

      ianf = 1
      do lauf = 1,nr            
         C(1:n,lauf)=vecf1(ianf:ianf+n-1)
         ianf = ianf + n
      enddo

c aufbau (h*bi - e*bi)/(e - haa); steht danach auf vecf2
      ianf = 1
      do lauf = 1,nr   
         call svdssvd(n,valn(lauf),vecf2(ianf),adiag,vecf1(ianf))
         ianf = ianf + n
      enddo

      denerg = 0.0d0
      do lauf = 1,nr
         denerg = denerg + abs(valn(lauf) - e(lauf))
         if ( abs(valn(lauf) - e(lauf)) .lt. crite) iconf(lauf) = 1
      enddo
      denerg = denerg / nr   
      ico = sum(iconf)
      if(pr) write(*,*) iter,ico,denerg,valn(1:nr)

      if (ico .eq. nr) then
          if(pr) write(*,*) 'all roots converged'
          go to 9999
      endif

      if(j.gt.0) then
c-- mit allen alten
         nneue = 0
         ialt = j
         ianf = 1
         do lauf = 1,nr    
c orthogonaliesung des betrachteten auf die alten
           do jalt = 1,ialt
              call smread(n,lun1,w,jalt)
              s=-sdot(n,w,1,vecf1(ianf),1)
              call saxpy(n,s,w,1,vecf1(ianf),1)
           enddo
c normierung dessen was vom betrachteten uebrig bleibt
           s=sdot(n,vecf1(ianf),1,vecf1(ianf),1)
           if (s.gt.0.00000001)then
c neuer wird mitgenommen
              s = ONE /sqrt(s)
              vecf1(ianf:ianf+n-1)= vecf1(ianf:ianf+n-1) * s
              ialt = ialt + 1
              nneue = nneue + 1
              ineue(nneue) = ianf
c wegschreiben des neuen zu den alten
              call smwrite(n,lun1,vecf1(ianf),jalt)
           else
              goto 9999
           endif
           ianf = ianf + n
         enddo
      endif

c umspeichern der ueberlebenden vektoren auf vecf1
      ianf = 1
      do lauf = 1,nneue
         call scopy(n,vecf1(ineue(lauf)),1,vecf1(ianf),1)
         ianf = ianf + n
      enddo

C H * C
      ianf = 1
      do i = 1,nneue
c        call sgemv('N',n,n,ONE,H,n,vecf1(ianf),1,ZERO,vecf2(ianf),1)
         call sspmv('U',n,  ONE,HP, vecf1(ianf),1,ZERO,vecf2(ianf),1)
         ianf = ianf + n
      enddo

      ianf = 1
      do i = 1,nneue
         call smwrite(n,lun2,vecf2(ianf),memlun2+i)   
         ianf = ianf + n
      enddo

c berechnung der neuen matrixelemente der davidson-matrix
c zunaechst mit den alten
      do jalt = 1,j 
        call smread(n,lun1,w,jalt)
        ianf = 1 
        ilauf = iideks(j) + jalt
        do jneu = 1,nneue 
           av(ilauf) = sdot(n,w,1,vecf2(ianf),1)
           ilauf = ilauf + jneu + j 
           ianf = ianf + n 
        enddo  
      enddo
c dann mit den neuen
      ianf = 1 
      do jneu = 1,nneue
         janf = 1 
         ilauf = iideks(j+jneu) - jneu + 1 
         do jneu1 = 1,jneu
            av(ilauf) = sdot(n,vecf2(ianf),1,vecf1(janf),1)
            janf = janf + n 
            ilauf = ilauf + 1 
        enddo
        ianf = ianf + n 
      enddo

c increase expansion space and iterate further
      e = valn 
      j = j + nneue

      enddo   

      if(pr) write(*,*) 'Warning: davidson not properly converged'
c exit

9999  continue

      deallocate(adiag,vecf1,vecf2,w,av,HP)
      close (lun1,status='delete')
      close (lun2,status='delete')

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine svdssvd(n,sl,v1,v2,v)
      implicit none                          
      real*4 sl, v1(*), v(*), v2(*)
      integer n

      integer i,m,mp1

      do i = 1,n
        v(i)=v1(i)/(sl-v2(i))
      enddo

      return
      end 

***********************************************************************

      subroutine smwrite(n,iwo,v,irec)
      implicit none
      real*4 v(n)
      integer n,iwo,irec
      write(iwo,rec=irec) v
      return
      end 


      subroutine smread(n,iwo,v,irec)
      implicit none
      real*4 v(n)
      integer n,iwo,irec
      read(iwo,rec=irec) v
      return
      end 





