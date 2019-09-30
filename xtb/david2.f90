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

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! real*4 version
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine sdavid2(n,crite,H,C,e,fail)
   use iso_fortran_env, wp => real32
   implicit none
   logical,parameter :: pr = .false.
   logical,parameter :: ini = .false.
   integer n       ! dimension
   integer,parameter :: nr = 1
   real(wp)  crite   ! eigenvalue convergence threshold
   real(wp)  H(n,n)  ! matrix to be diagonalized
   real(wp)  C(n,nr) ! eigenvectors
   real(wp)  e(nr)   ! eigenvalues
   logical,intent(out) :: fail

   ! local
   integer maxiter        ! maximum # of iterations
   parameter (maxiter=100)
   integer iter,ineue(1),janf,lun1,lun2
   integer iideks(maxiter),idum,j,jalt,ilauf,jneu
   integer l1,l2,k,LWORK,LIWORK,INFO,i,ien,ialt,memlun2
   integer,allocatable :: iwork(:)
   logical lconf

   real(wp) valn(1),uim,s,sdot,one,zero,denerg
   real(wp), allocatable :: adiag(:),vecf1(:),vecf2(:),w(:)
   real(wp), allocatable :: Uaug(:,:),d(:),aux(:)
   real(wp), allocatable :: AB(:,:),av(:),tmpav(:,:)
   real(wp), allocatable :: HP(:)
   parameter (one =1)
   parameter (zero=0)

   fail = .true.

   if(pr)then
      write(*,'(/,10x,''******************************************'')')
      write(*,'(10x,''*            multi-root davidson (R4)    *'')')
      write(*,'(10x,''******************************************'',/)')
      write(*,*) 'dim ',n,' # roots ',1
   endif

   allocate(adiag(n),vecf1(n),vecf2(n),w(n),av(maxiter*(maxiter+1)/2), &
      &     HP(n*(n+1)/2))

   ! IO
   lun1=86
   lun2=87
   open(unit=lun1,recl=n*4,access='direct',file='david1.tmp',form='unformatted')
   open(unit=lun2,recl=n*4,access='direct',file='david2.tmp',form='unformatted')

   k=0
   do i=1,n
      do j=1,i
         k=k+1
         HP(k)=H(j,i)
      enddo
   enddo

   ! H * C for initialization
   call smwrite(n,lun1,C(1,1),1)
   call sspmv('U',n,  ONE,HP, C(1,1),1,ZERO,vecf2,1)
   call smwrite(n,lun2,vecf2,1)

   ! aufbau des iideks feldes
   iideks(1) = 1
   do idum = 2,maxiter
      iideks(idum) = iideks(idum - 1) + idum
   enddo
   valn = 0
   lconf= .false.
   e    = 0

   do i=1,n
      adiag(i)=H(i,i)
   enddo

   av(1)=sdot(n,C(1,1),1,vecf2,1)
   ! done

   j = 1

   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   ! eigentliche schleife im davidson
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   do iter = 1, maxiter-1

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
      valn(1:1) = d(1:1)

      ! aufbau der eigentlichen vektoren, die stehen dann auf vecf1
      ! multiplikation der vorherigen entwicklungsvektoren mit den entwicklungskoeffizienten
      vecf1=0.0_wp
      do i=1,j
         call smread(n,lun1,w,i)
         uim = Uaug(i,1)
         call saxpy(n,uim,w,1,vecf1,1)
      enddo

      ! aufbau -E*bi: vecf1 ist bi; vecf2 dann E*bi
      vecf2 = -valn(1) * vecf1
      ! aufbau des residuen vektors (h*bi-e*bi)
      ! steht dann auf vecf2 ueberschreibt also -e*bi
      ! h*bi-e*bi steht jetzt auf vecf2
      do i=1,j
         call smread(n,lun2,w,i)
         memlun2 = i
         uim = Uaug(i,1)
         call saxpy(n,uim,w,1,vecf2,1)
      enddo
      deallocate(Uaug,d,iwork,aux)

      C(1:n,1)=vecf1

      ! aufbau (h*bi - e*bi)/(e - haa); steht danach auf vecf2
      vecf1=vecf2/(valn(1)-adiag)

      denerg = abs(valn(1) - e(1))
      lconf = abs(valn(1) - e(1)) .lt. crite
      if (pr) write(*,*) iter,lconf,denerg,valn(1:1)

      if (lconf) then
         if (pr) write(*,*) 'all roots converged'
         fail = .false.
         go to 99
      endif

      if (j.gt.0) then
         !-- mit allen alten
         ialt = j
         ! orthogonaliesung des betrachteten auf die alten
         do jalt = 1,ialt
            call smread(n,lun1,w,jalt)
            s=-sdot(n,w,1,vecf1,1)
            call saxpy(n,s,w,1,vecf1,1)
         enddo
         ! normierung dessen was vom betrachteten uebrig bleibt
         s=sdot(n,vecf1,1,vecf1,1)
         if (s.gt.0.00000001) then
            ! neuer wird mitgenommen
            s = ONE /sqrt(s)
            vecf1= vecf1 * s
            ialt = ialt + 1
            ! wegschreiben des neuen zu den alten
            call smwrite(n,lun1,vecf1,jalt)
         else
            fail = .false.
            goto 99
         endif
      endif

      ! H * C
      call sspmv('U',n,  ONE,HP, vecf1,1,ZERO,vecf2,1)

      call smwrite(n,lun2,vecf2,memlun2+1)

      ! berechnung der neuen matrixelemente der davidson-matrix
      ! zunaechst mit den alten
      do jalt = 1,j
         call smread(n,lun1,w,jalt)
         ilauf = iideks(j) + jalt
         av(ilauf) = sdot(n,w,1,vecf2,1)
         ilauf = ilauf + 1 + j
      enddo
      ! dann mit den neuen
      av(iideks(j+1)) = sdot(n,vecf2,1,vecf1,1)

      ! increase expansion space and iterate further
      e = valn
      j = j + 1

   enddo

   if(pr) write(*,*) 'Warning: davidson not properly converged'
   ! exit

99 continue

   deallocate(adiag,vecf1,vecf2,w,av,HP)
   close (lun1,status='delete')
   close (lun2,status='delete')

end subroutine sdavid2

subroutine solver_sdavidson(n,crite,Hp,C,e,fail,pr)
   use iso_fortran_env, wp => real32
   implicit none
   logical, intent(in) :: pr
   logical,parameter :: ini = .false.
   integer n       ! dimension
   integer,parameter :: nr = 1
   real(wp)  crite   ! eigenvalue convergence threshold
   real(wp)  Hp(n*(n+1)/2)  ! matrix to be diagonalized
   real(wp)  C(n,nr) ! eigenvectors
   real(wp)  e(nr)   ! eigenvalues
   logical,intent(out) :: fail

   ! local
   integer maxiter        ! maximum # of iterations
   parameter (maxiter=100)
   integer iter,ineue(1),janf!,lun1,lun2
   integer iideks(maxiter),idum,j,jalt,ilauf,jneu
   integer l1,l2,k,LWORK,LIWORK,INFO,i,ien,ialt,memlun2
   integer,allocatable :: iwork(:)
   logical lconf
   real(wp), allocatable :: lun1(:,:),lun2(:,:)
   integer, parameter :: initial_dyn_array_size = 10

   real(wp) valn(1),uim,s,sdot,denerg
   real(wp), allocatable :: adiag(:),vecf1(:),vecf2(:),w(:)
   real(wp), allocatable :: Uaug(:,:),d(:),aux(:)
   real(wp), allocatable :: AB(:,:),av(:),tmpav(:,:)

   fail = .true.

   if (pr) then
      write(*,'(/,10x,''******************************************'')')
      write(*,'(10x,''*            multi-root davidson (R4)    *'')')
      write(*,'(10x,''******************************************'',/)')
      write(*,*) 'dim ',n,' # roots ',1
   endif

   allocate(adiag(n),vecf1(n),vecf2(n),w(n),av(maxiter*(maxiter+1)/2))

   allocate(lun1(n,initial_dyn_array_size), lun2(n,initial_dyn_array_size), &
      &     source = 0.0_wp)
   ! IO
   !lun1=86
   !lun2=87
   !open(unit=lun1,recl=n*4,access='direct',status='scratch',form='unformatted')
   !open(unit=lun2,recl=n*4,access='direct',status='scratch',form='unformatted')

   ! H * C for initialization
   call smwrite(n,lun1,C(1,1),1)
   call sspmv('U',n,  1.0_wp,HP, C(1,1),1,0.0_wp,vecf2,1)
   call smwrite(n,lun2,vecf2,1)

   ! aufbau des iideks feldes
   iideks(1) = 1
   do idum = 2,maxiter
      iideks(idum) = iideks(idum - 1) + idum
   enddo
   valn = 0
   lconf= .false.
   e    = 0

   do i=1,n
      adiag(i)=HP(i*(i-1)/2)
   enddo

   av(1)=sdot(n,C(1,1),1,vecf2,1)
   ! done

   j = 1

   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   ! eigentliche schleife im davidson
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   do iter = 1, maxiter-1

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
      valn(1:1) = d(1:1)

      ! aufbau der eigentlichen vektoren, die stehen dann auf vecf1
      ! multiplikation der vorherigen entwicklungsvektoren mit den entwicklungskoeffizienten
      vecf1=0.0_wp
      do i=1,j
         call smread(n,lun1,w,i)
         uim = Uaug(i,1)
         call saxpy(n,uim,w,1,vecf1,1)
      enddo

      ! aufbau -E*bi: vecf1 ist bi; vecf2 dann E*bi
      vecf2 = -valn(1) * vecf1
      ! aufbau des residuen vektors (h*bi-e*bi)
      ! steht dann auf vecf2 ueberschreibt also -e*bi
      ! h*bi-e*bi steht jetzt auf vecf2
      do i=1,j
         call smread(n,lun2,w,i)
         memlun2 = i
         uim = Uaug(i,1)
         call saxpy(n,uim,w,1,vecf2,1)
      enddo
      deallocate(Uaug,d,iwork,aux)

      C(1:n,1)=vecf1

      ! aufbau (h*bi - e*bi)/(e - haa); steht danach auf vecf2
      vecf1=vecf2/(valn(1)-adiag)

      denerg = abs(valn(1) - e(1))
      lconf = abs(valn(1) - e(1)) .lt. crite
      if (pr) write(*,*) iter,lconf,denerg,valn(1:1)

      if (lconf) then
         if (pr) write(*,*) 'all roots converged'
         fail = .false.
         go to 99
      endif

      if (j.gt.0) then
         !-- mit allen alten
         ialt = j
         ! orthogonaliesung des betrachteten auf die alten
         do jalt = 1,ialt
            call smread(n,lun1,w,jalt)
            s=-sdot(n,w,1,vecf1,1)
            call saxpy(n,s,w,1,vecf1,1)
         enddo
         ! normierung dessen was vom betrachteten uebrig bleibt
         s=sdot(n,vecf1,1,vecf1,1)
         if (s.gt.0.00000001) then
            ! neuer wird mitgenommen
            s = 1.0_wp /sqrt(s)
            vecf1= vecf1 * s
            ialt = ialt + 1
            ! wegschreiben des neuen zu den alten
            call smwrite(n,lun1,vecf1,jalt)
         else
            fail = .false.
            goto 99
         endif
      endif

      ! H * C
      call sspmv('U',n,  1.0_wp,HP, vecf1,1,0.0_wp,vecf2,1)

      call smwrite(n,lun2,vecf2,memlun2+1)

      ! berechnung der neuen matrixelemente der davidson-matrix
      ! zunaechst mit den alten
      do jalt = 1,j
         call smread(n,lun1,w,jalt)
         ilauf = iideks(j) + jalt
         av(ilauf) = sdot(n,w,1,vecf2,1)
         ilauf = ilauf + 1 + j
      enddo
      ! dann mit den neuen
      av(iideks(j+1)) = sdot(n,vecf2,1,vecf1,1)

      ! increase expansion space and iterate further
      e = valn
      j = j + 1

   enddo

   if(pr) write(*,*) 'Warning: davidson not properly converged'
   ! exit

99 continue

   deallocate(adiag,vecf1,vecf2,w,av,lun1,lun2)
   !close (lun1,status='delete')
   !close (lun2,status='delete')

contains

subroutine smwrite(n,iwo,v,irec)
   implicit none
   real(wp), intent(inout), allocatable :: iwo(:,:)
   real(wp), intent(in)  :: v(n)
   integer,  intent(in)  :: n,irec
   real(wp), allocatable :: tmp(:,:)
   integer :: d2,dn
   d2 = size(iwo,2)
   if (irec > d2) then
      dn = d2 + d2/2 + 1
      allocate(tmp(n,dn))
      tmp(:,:d2) = iwo
      deallocate(iwo)
      call move_alloc(tmp,iwo)
   endif
   iwo(:,irec) = v
   !write(iwo,rec=irec) v
   !return
end subroutine smwrite

subroutine smread(n,iwo,v,irec)
   implicit none
   real(wp), intent(out) :: v(n)
   real(wp), intent(in)  :: iwo(:,:)
   integer,  intent(in)  :: n,irec
   v = iwo(:,irec)
   !read(iwo,rec=irec) v
   !return
end subroutine smread

end subroutine solver_sdavidson
