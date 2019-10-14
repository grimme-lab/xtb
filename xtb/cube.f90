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

!cccccccccccccccccccccccccccccccc
!    write out cube file        c
!cccccccccccccccccccccccccccccccc
! ncent  : # atoms
! nmo    : # MOs
! nbf    : # AOs
! nprims : # primitives (in total)
! xyz(4,ncent) : Cartesian coordinates & nuclear charge
! cxip(nprims) : contraction coefficients of primitives
! exip(nprims) : exponents of primitives
! cmo(nbf,nmo) : LCAO-MO coefficients
! eval(nmo)    : orbital eigenvalues
! occ(nmo)     : occupation # of MO
! ipty(nprims) : angular momentum of primitive function
! ipao(nbf)    : # primitives in contracted AO
! ibf(ncent)   : # of contracted AOs on atom

subroutine cube(n,nmo,nbf,xyz,at,cmo,eval,occ,fname,basis)
   use tbdef_basisset
   use setparam
   implicit none
   type(tb_basisset), intent(in) :: basis

   real*8, intent ( in ) :: xyz(3,n)
   real*8, intent ( in ) :: eval(nmo)
   real*8, intent ( in ) :: occ (nmo)
   real*8, intent ( in ) :: cmo(nbf,nmo)
   integer, intent( in ) :: at(n)
   integer, intent( in ) :: n,nmo,nbf
   character*(*) fname

   real*8 ,allocatable  ::Ptmp (:,:)
   real*8 ,allocatable  ::P    (:,:)
   real*8 ,allocatable  ::C    (:,:)
   real*4 ,allocatable  ::cb   (:,:,:)
   real*8 ,allocatable  ::array(:)
   integer,allocatable  ::matlist (:,:)

   integer i,j,k,m,nm,ii,jj,iii,jjj,npri,nprj,primstart(nbf),nexp
   integer iat,jat,xst,yst,zst,cou,iprimcount,jprimcount,iiii,jjjj
   real*8 thr,thr3,intcut,intcut2
   real*8 w0,w1,t0,t1,r,r2,val,rab,ccc,gridp(3),xyza(3),xyzb(3),dr3
   real*8 step,px,py,pz,xinc,yinc,zinc,nx,ny,nz,vv(1),v,est,nfod
   real*8 f1,f2,dx1,dy1,dz1,dx2,dy2,dz2,r1,dum,r1xy,r2xy
   real*8 dxx1,dyy1,dzz1,dxx2,dyy2,dzz2,ar1,ar2,cc
   real*8 fastexp
   integer ifile

   call timing(t0,w0)
   write(*,*)
   write(*,*)'cube file module (SG, 7/16)'
   thr = cube_pthr ! Dmat pre-screen
   step= cube_step ! grid step (Bohr)
   intcut=8.00d0   ! primitive cut
   intcut2=2.0d0*intcut

   nexp=100*int(intcut2) ! size of exp array ie arguments 0:nexp
   allocate(array(0:nexp))
   do i=0,nexp
      dum=float(i)/100.0
      array(i)=exp(-dum)
   enddo

   !     do i=0,160
   !        dum=float(i)/10.0
   !        write(*,*) int(10*dum)
   !        write(*,*) exp(-dum),array(int(10*dum)),fastexp(nexp,array,dum)
   !     enddo
   !     stop

   allocate(Ptmp(nbf,nbf),P(nbf,nbf),C(nbf,nbf),matlist(2,nbf*(nbf+1)/2))
   C(1:nbf,1:nmo) = cmo(1:nbf,1:nmo)
   do m=1,nmo
      do i=1,nbf
         Ptmp(i,m)=cmo(i,m)*occ(m)
      enddo
   enddo
   if(nmo.ne.nbf)then
      Ptmp(1:nbf,nmo+1:nbf)=0
      C   (1:nbf,nmo+1:nbf)=0
   endif
   call DGEMM('N','T',nbf,nbf,nbf, 1.0d0,C,nbf,Ptmp,nbf,0.0d0,P,nbf)
   deallocate(Ptmp,C)

   nm=0
   npri=0
   do i=1,nbf
      primstart(i)=npri  ! start index for prims of AO i
      npri=npri+basis%nprim(i)
      do j=1,i
         if(abs(P(j,i)).gt.thr)then
            nm=nm+1
            matlist(1,nm)=i
            matlist(2,nm)=j
            if(i.ne.j) P(j,i)=2.0d0*P(j,i)
         endif
      enddo
   enddo
   write(*,'('' cube_pthr     : '',f7.3)')cube_pthr
   write(*,'('' cube_step     : '',f7.3)')cube_step
   write(*,'('' non-zero P (%): '',f7.3,''   nmat:'',i8)') &
   & 100.*float(nm)/float(nbf*(nbf+1)/2),nm

   px=maxval(xyz(1,1:n))+3.0
   py=maxval(xyz(2,1:n))+3.0
   pz=maxval(xyz(3,1:n))+3.0
   nx=minval(xyz(1,1:n))-3.0
   ny=minval(xyz(2,1:n))-3.0
   nz=minval(xyz(3,1:n))-3.0
   write(*,*)'Grid Boundaries (x y z) :'
   write(*,*)px,py,pz
   write(*,*)nx,ny,nz

   ! calculate step size and number of steps (step size approx 0.2-0.5)
   xst=floor((abs(px)+abs(nx))/step)
   xinc=(abs(px)+abs(nx))/xst
   yst=floor((abs(py)+abs(ny))/step)
   yinc=(abs(py)+abs(ny))/yst
   zst=floor((abs(pz)+abs(nz))/step)
   zinc=(abs(pz)+abs(nz))/zst
   write(*,*)'Total # of points', (xst+1)*(yst+1)*(zst+1)
   dr3=xinc*yinc*zinc

   write(*,*)'writing ',trim(fname)
   call open_file(ifile,fname,'w')
   write(ifile,*)'xtb spin/fo density'
   write(ifile,*)'By S.G.'
   write(ifile,101)n,nx,ny,nz
   write(ifile,101)xst+1,xinc,0.0,0.0
   write(ifile,101)yst+1,0.0,yinc,0.0
   write(ifile,101)zst+1,0.0,0.0,zinc
   do i=1,n
      write(ifile,102)at(i),0.0,xyz(1,i),xyz(2,i),xyz(3,i)
   enddo

   allocate(cb(0:zst,0:yst,0:xst))

   nfod=0
   cb  =0

   !     Dmat loop         -----------------------------------
   do m=1,nm
      ii=matlist(1,m)
      jj=matlist(2,m)
      iat=basis%aoat(ii)
      jat=basis%aoat(jj)
      xyza(1:3)=xyz(1:3,iat)
      xyzb(1:3)=xyz(1:3,jat)
      npri=basis%nprim(ii)
      nprj=basis%nprim(jj)
      iprimcount=primstart(ii)
      jprimcount=primstart(jj)
      rab=(xyza(1)-xyzb(1))**2  &
      &      +(xyza(2)-xyzb(2))**2  &
      &      +(xyza(3)-xyzb(3))**2
      !        prim loop
      do iiii=1,npri
         iii=iprimcount+iiii
         ccc=P(jj,ii)*basis%cont(iii)
         dum=rab*basis%alp(iii)
         do jjjj=1,nprj
            jjj=jprimcount+jjjj
            est=dum*basis%alp(jjj)/(basis%alp(iii)+basis%alp(jjj))
            if(est.lt.intcut)then
               cc=ccc*basis%cont(jjj)
               !                 grid loops
               gridp(1)=nx
               do i=0,xst
                  gridp(1)=nx+(xinc*i)
                  dx1=xyza(1)-gridp(1)
                  dx2=xyzb(1)-gridp(1)
                  dxx1=dx1*dx1
                  dxx2=dx2*dx2
                  gridp(2)=ny
                  do j=0,yst
                     gridp(2)=ny+(yinc*j)
                     dy1=xyza(2)-gridp(2)
                     dy2=xyzb(2)-gridp(2)
                     dyy1=dy1*dy1
                     dyy2=dy2*dy2
                     gridp(3)=nz
                     r1xy=dxx1+dyy1
                     r2xy=dxx2+dyy2
                     do k=0,zst
                        gridp(3)=nz+(zinc*k)
                        dz1=xyza(3)-gridp(3)
                        dz2=xyzb(3)-gridp(3)
                        dzz1=dz1*dz1
                        dzz2=dz2*dz2
                        r1=r1xy+dzz1
                        !                    primitive function value
                        ar1=basis%alp(iii)*r1
                        if(ar1.lt.intcut2)then  ! exp(-16) < 1.d-7 i.e. zero
                           call primvalf(dx1,dy1,dz1,dxx1,dyy1,dzz1,ar1, &
                              &             basis%lao(ii),nexp,array,f1)
                           r2=r2xy+dzz2
                           ar2=basis%alp(jjj)*r2
                           if(ar2.lt.intcut2)then
                              call primvalf(dx2,dy2,dz2,dxx2,dyy2,dzz2,ar2, &
                                 &             basis%lao(jj),nexp,array,f2)
                              cb(k,j,i)=cb(k,j,i)+cc*f1*f2
                           endif
                        endif
                     enddo
                  enddo
               enddo
            endif
         enddo
      enddo
   enddo
   !     Dmat loop end     -----------------------------------

   ! write
   cou=1
   do i=0,xst
      do j=0,yst
         do k=0,zst
            if (cou.lt.6) then
               write(ifile,'(E14.8,1X)',advance='no')cb(k,j,i)
               cou=cou+1
            else
               write(ifile,'(E14.8)')cb(k,j,i)
               cou=1
            endif
            nfod=nfod+cb(k,j,i)*dr3
         enddo
      enddo
   enddo
   call close_file(ifile)

   101   format(I5,3F16.6)
   102   format(I5,4F16.6)

   write(*,'('' N    (numint) : '',f6.3)')nfod
   call timing(t1,w1)
   call prtime(6,t1-t0,w1-w0,'cube')
end subroutine cube


subroutine primval(dx,dy,dz,dx2,dy2,dz2,alpr2,lao,f)
   implicit none
   integer lao
   real*8 f,dx2,dy2,dz2
   real*8 dx,dy,dz,alpr2

   goto (100,201,202,203,301,302,303,304,305,306) lao

   100  f=exp(-alpr2)
   return
   201  f=exp(-alpr2)*dx
   return
   202  f=exp(-alpr2)*dy
   return
   203  f=exp(-alpr2)*dz
   return
   301  f=exp(-alpr2)*dx2
   return
   302  f=exp(-alpr2)*dy2
   return
   303  f=exp(-alpr2)*dz2
   return
   304  f=exp(-alpr2)*dx*dy
   return
   305  f=exp(-alpr2)*dx*dz
   return
   306  f=exp(-alpr2)*dy*dz
   return

end subroutine primval

! routine with exp table lookup
subroutine primvalf(dx,dy,dz,dx2,dy2,dz2,alpr2,lao,nexp,a,f)
   implicit none
   integer lao,nexp
   real*8 a(0:nexp)
   real*8 f,dx2,dy2,dz2
   real*8 dx,dy,dz,alpr2
   real*8 fastexp

   goto (100,201,202,203,301,302,303,304,305,306) lao

   100  f=fastexp(nexp,a,alpr2)
   return
   201  f=fastexp(nexp,a,alpr2)*dx
   return
   202  f=fastexp(nexp,a,alpr2)*dy
   return
   203  f=fastexp(nexp,a,alpr2)*dz
   return
   301  f=fastexp(nexp,a,alpr2)*dx2
   return
   302  f=fastexp(nexp,a,alpr2)*dy2
   return
   303  f=fastexp(nexp,a,alpr2)*dz2
   return
   304  f=fastexp(nexp,a,alpr2)*dx*dy
   return
   305  f=fastexp(nexp,a,alpr2)*dx*dz
   return
   306  f=fastexp(nexp,a,alpr2)*dy*dz
   return

end subroutine primvalf

real*8 function fastexp(nexp,array,arg)
   integer nexp
   real*8 arg
   real*8 array(0:nexp)

   fastexp=array(int(100.0d0*arg))

end function fastexp
