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


c  .....................................................................
c
c   Calculates the pointgroup, the nuclear exchange table and the
c   transfomation matrices
c
c  .....................................................................

      Subroutine symtrans(sflies,thrsym,xyz,natoms,ict,trans,ntrans)
      use iso_fortran_env, only : wp => real64
      
c  .....................................................................
c
c   Input
c   -----
c
c    xyz     : cartesian coordinates of the atoms
c    
c    natoms  : number of the atoms
c
c    thrsym  : threshhold
c
c
c   Output
c   ------
c
c    sflies  : Schoenflies symbol
c
c    ict     : nuclear exchange table 
c 
c    ntrans  : number of transformations
c
c    trans   : transformation matrices
c
c  .....................................................................

      implicit none

      integer ndi14
      parameter (ndi14=120)

      logical
     ,            show 

      character
     ,            sflies*(*),csf*1,cla*1

      integer
     ,            natoms,nn,nprisy,iback,ngen,ntrans,mi(ndi14),
     ,            mj(ndi14),natomsred,nhessred,
     ,            invt(ndi14),ict(natoms,ndi14),ierr,i,j,k,l,
     ,            jkseq(2,natoms),ictnr(natoms,ndi14),ictt,
     ,            lieff(natoms)

      real(wp)
     ,            xyz(3,natoms),gen(9,3),trans(9,ndi14),thrsym


      nprisy=-1
      show=.false.
csg
      ict = 0

c ... Decompose input Schoenflies symbol
      call grpsmb(sflies,csf,cla,nn,nprisy)

c ... Check validity of Schoenflies symbol
      if (.not.(csf.eq.'s'.and.cla.eq.' '
     ,     .or. csf.eq.'c'.and.(cla.eq.'h'.or.cla.eq.'v'.or.cla.eq.' ')
     ,     .or. csf.eq.'d'.and.(cla.eq.'h'.or.cla.eq.'d'.or.cla.eq.' ')
     ,     .or. csf.eq.'t'.and.(cla.eq.'h'.or.cla.eq.'d'.or.cla.eq.' ')
     ,     .or. csf.eq.'o'.and.(cla.eq.'h'.or.cla.eq.' ')
     ,     .or. csf.eq.'i'.and.(cla.eq.'h'.or.cla.eq.' '))) then
c        write (*,'(a)') '<symmetry> INVALID SCHOENFLIES SYMBOL ! '
c        write (*,*) '*',sflies,'*',csf,cla,nn
         call flush(6)
c        error STOP 
       end if

c ... Form generators of group corresponding to Schoenflies symbol
      call getgen(csf,nn,cla,gen,ngen,iback,nprisy)
      if(iback.ne.0) then
        write (*,'(a)') '<symmetry> : Error in getgen' 
        call flush(6)
        error STOP 
      end if

c ... Generate transformation matrices 'trans' of the symmetry group
c ... Find inverse group operations 'invt' and memorize group
c ... generation path 'mi','mj' for subsequent construction
c ... of representations
      call group(gen,ngen,trans,ntrans,ndi14,thrsym,invt,mi,mj,iback,
     ,           nprisy)
 
c ...
c ... determines the nuclear exchange group
c ... ict(i,ig) is the nucleus obtained on nucleus i
c ... by means of the ig's symmetry operation
c ... symmetries specified as 3*3 matrices on array trans
c ... ntrans = group order
c ... natoms = number of atoms
c ... xyz = nuclear coordinates
c ...
      call nucex(ict,natoms,trans,ntrans,thrsym,xyz,natoms,ierr,show)

c ... delete symmetry redundant atoms in ict
      ictnr=ict
      do i=1,natoms
        do k=2,ntrans
          ictt=ict(i,k)
          do l=1,(k-1)
            if (ictt.eq.ict(i,l)) ictnr(i,k)=0 
          end do
        end do
      end do  
    
c ... get list of sym.equivalent atoms jkseq(.,i), i=1,natoms
      do i=1,natoms
        jkseq(1,i)=0
        jkseq(2,i)=0
        do j=1,(i-1)
          do k=2,ntrans
            if (i.eq.ict(j,k).and.jkseq(1,i).eq.0) then
              jkseq(1,i)=j
              jkseq(2,i)=k
              exit
            end if
          end do
        end do 
      end do
      natomsred=0
      do i=1,natoms
        if (jkseq(1,i).eq.0) then
          natomsred=natomsred+1
          lieff(natomsred)=i
        end if
      end do
      nhessred=3*natomsred
cc    write (*,*) 'natomsred=',natomsred
cc    write (*,*) 'lieff(.)=',(lieff(i),i=1,natomsred)

cc    write(*,*) 'sflies',sflies,'/'
cc    write(*,*) 'ict ntrans=',ntrans
cc    do i=1,natoms
cc      write(*,*)(ict(i,j),j=1,ntrans)
cc    end do
cc    write(*,*) 'trans'
cc    do j=1,ntrans
cc      write(*,'(9f7.4)')(trans(i,j),i=1,9)
cc    end do
cc    write(*,*) 'ictnr',ntrans
cc    do i=1,natoms
cc      write(*,*)(ictnr(i,j),j=1,ntrans)
cc    end do
cc    write(*,*) 'jkseq'
cc    do i=1,natoms
cc      write(*,*)(jkseq(j,i),j=1,2)
cc    end do

      End


c RCS ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Id: grpsmb.f,v 1.2 1998/11/03 15:21:49 klaus Exp $
c $Log: grpsmb.f,v $
c Revision 1.2  1998/11/03 15:21:49  klaus
c change getcor into allocate (F90)
c
c Revision 1.1  1992/09/10 13:41:20  tomjones
c Initial revision
c
c RCS ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c--------------------------------------------------------------------
      subroutine grpsmb(sflies,csf,cla,nn,nprisy)
      use iso_fortran_env, only : wp => real64
      implicit real(wp)(a-h,o-z)
      character csf,cla,zn(0:9),sflies*4
      data zn /'0','1','2','3','4','5','6','7','8','9'/
c
      nn=0
c     nn is the symmetry number of the z-axis
      csf=sflies(1:1)
      do 5 i=1,9
         if (sflies(2:2).eq.zn(i)) nn=i
    5    continue
      if (nn.eq.0) then
         cla=sflies(2:2)
c        if (nprisy.ge.-1) write (6,'(/,
c    1      '' symmetry group of the molecule :   '',2a1)') csf,cla
         if (csf.eq.'c'.and.(cla.eq.'s'.or.(nn.eq.0.and.cla.eq.'h')))
     1   then
            nn=1
            cla='h'
         elseif (csf.eq.'c' .and. cla.eq.'i') then
            nn=2
            cla=' '
            csf='s'
         endif
      else
         np=-1
         do 10 i=0,9
            if (sflies(3:3).eq.zn(i)) np=i
   10       continue
         if (np.ge.0) then
            nn=10*nn+np
            cla=sflies(4:4)
c           if (nprisy.ge.-1)
c    1        write (6,'(/,'' symmetry group of the molecule :   '',
c    2           a1,i2,a1)') csf,nn,cla
         else
            cla=sflies(3:3)
c           if (nprisy.ge.-1)
c    1        write (6,'(/,'' symmetry group of the molecule :   '',
c    2           a1,i1,a1)') csf,nn,cla
         endif
      endif
      return
      end
c RCS ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Id: group.f,v 1.3 1998/11/03 15:21:36 klaus Exp $
c $Log: group.f,v $
c Revision 1.3  1998/11/03 15:21:36  klaus
c change getcor into allocate (F90)
c
c Revision 1.2  1997/10/01 13:07:02  marcok
c The hp780 ld complained about not initialized variables. Done now.
c
c Revision 1.1  1992/09/10 13:41:12  tomjones
c Initial revision
c
c RCS ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c--------------------------------------------------------------------
      subroutine group(gen,ngen,u,nt,max14,thrsym,inv,ipath,jpath,ierr,
     1                 nprisy)
      use iso_fortran_env, only : wp => real64
      implicit real(wp)(a-h,o-z)
      dimension u(9,max14),gen(9,ngen),inv(max14),scr(9)
      dimension ipath(max14),jpath(max14)
c      generates 3x3 matrices forming group
c      from ngen generators on gen
c      output:
c      u = 3x3 matrices
c      nt = number of symmetry operations
c      max14 = maximum order of group
c      inv  containes index of invers operation
      it=0
      a0=0
      a1=1
      ierr=0
c      set neutral element
c      (every point group has the neutral element as first operation)
      nt = 1
      do 200 i=1,9
 200  u(i,nt) = a0
      u(1,nt) = a1
      u(5,nt) = a1
      u(9,nt) = a1
      ipath(1)=0
      jpath(1)=0
      if (ngen.eq.0) goto 700
      do 500 igen=1,ngen
c         loop over generators
         ncoset = 1
         nto = nt
         do 220 i=1,9
 220        scr(i) = gen(i,igen)
c        check if generator is already contained in (sub)group
c        of order nto :
         do 240 j=1,nt
            if(dif(u(1,j),scr(1),9).lt.thrsym) go to 500
 240        continue
         ipath(nt+1)=-igen
         jpath(nt+1)=0
c        for new group element form coset  :
 245     do 250 it=1,nto
            nt = nt + 1
            if (nt.gt.max14) goto 900
            if (it.eq.1) goto 250
            ipath(nt)=nto*ncoset+1
            jpath(nt)=it
 250        call mult3 (u(1,nto*ncoset+it),scr,u(1,it),3)
         ncoset = ncoset + 1
c        check : do subgroup plus cosets form a group ?
         do 330 it = 1,nt
            do 320 jcoset = 1,ncoset-1
               call mult3 (scr,u(1,it),u(1,1+nto*jcoset),3)
               do 310 kt = 1,nt
                  if (dif(u(1,kt),scr(1),9).lt.thrsym) goto 320
 310              continue
c              new symmetry operation found
               ipath(nt+1)=it
               jpath(nt+1)=1+nto*jcoset
               goto 245
 320           continue
 330        continue
c        cosets now form group - take next generator
 500     continue
      if (nprisy.ge.2) then
       write(6,'(/,i5,a)') nt,' symmetry operations found :'
      elseif (nprisy.ge.0) then
       write(6,'(/,i5,a)') nt,' symmetry operations found'
      endif
c
c     find inverse operators
      inv(1) = 1
      if (nt.eq.1) goto 700
      do 660 i=2,nt
      inv(i) = 0
      do 640 j=2,nt
      call mult3 (scr,u(1,i),u(1,j),3)
      s = dif (scr,u(1,1),9)
      if (s.gt.thrsym) go to 640
      inv(i) = j
      go to 660
  640 continue
      go to 950
  660 continue
  700 continue
c
      if (nprisy.ge.2) then
         write(6,750)
  750    format(/,10x,'xx',6x,'yx',6x,'zx',6x,'xy',6x,'yy',6x,
     1          'zy',6x,'xz',6x,'yz',6x,'zz')
         do 5 it=1,nt
    5       write(6,'(i4,1x,9f8.4)') it,(u(i,it),i=1,9)
         write(6,*)
         write(6,'(/,'' group element generation record'',
     1          '' (negative numbers correspond to generators) :'')')
         write (6,'(6(2x,i3,''='',i3,''*'',i3))')
     1            (it,ipath(it),jpath(it),it=1,nt)
      end if
c
      return
c
  900 write(6,'(a,i4)') ' size of group larger than parameter ndi14 =',
     1                  max14
      ierr=1
      return
c
  950 write(6,*) ' inversion problem '
      ierr=1
      return
      end
c RCS cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Id: nucex.f,v 1.4 2000/10/20 13:51:18 haettig Exp $
c $Log: nucex.f,v $
c Revision 1.4  2000/10/20 13:51:18  haettig
c replaced fortran STOPs by call to quit routine.
c
c Revision 1.3  1998/11/03 15:22:26  klaus
c change getcor into allocate (F90)
c
c Revision 1.2  1995/04/05 17:14:37  holger
c Source formatting.
c
c Revision 1.1  1992/09/10  13:42:02  tomjones
c Initial revision
c RCS cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nucex(ict,max10,trans,ntrans,thrsym,xyz,natoms,
     1                 ierr,show)
      use iso_fortran_env, only : wp => real64
c
c     determines the nuclear exchange group
c     ict(i,ig) is the nucleus obtained on nucleus i
c     by means of the ig's symmetry operation
c     symmetries specified as 3*3 matrices on array trans
c     ntrans = group order
c     natoms = number of atoms
c     xyz = nuclear coordinates
c
      implicit real(wp)(a-h,o-z)
      dimension ict(max10,ntrans),trans(9,ntrans),scr(3),xyz(3,natoms)
      logical show
c
c     set neutral operation
c
      ierr=0 
      imotz=ierr
csg   ierr=0
      do 120 ic=1,natoms
         ict(ic,1) = ic
 120  continue
      if (ntrans.eq.1) goto 400
      do 300 ic=1,natoms
         do 250 it=2,ntrans
            ict(ic,it) = 0
            call mult3 (scr,trans(1,it),xyz(1,ic),1)
            do 200 jc=1,natoms
               if (dif(scr,xyz(1,jc),3).gt.thrsym) goto 200
               ict(ic,it) = jc
               goto 250
 200        continue
            goto 500
 250     continue
 300  continue
 400  continue
csg
c     open(unit=66,file='trans',form='unformatted')
c     write(66)ntrans
c     do it=1,ntrans
c        write(66) (trans(k,it),k=1,9)
c     enddo    
c     do iatom=1,natoms
c        write(66) (ict(iatom,it),it=1,ntrans)
c     enddo       
c     close(66)
csg
      if (show) then
         write(6,'(/,10x,a,/)') 'nuclear exchange table :'
         write(6,'(1x,a2,2x,24i3)') 'it',(it,it=1,ntrans)
         write(6,*)
         do 30 iatom=1,natoms
            write(6,'(5x,24i3)') (ict(iatom,it),it=1,ntrans)
 30      continue
      endif
      return
 500  ierr=1
      if (imotz.eq.-1) return
      write(6,*) ' error in nuclear exchange group determination'
      write(6,*) ' coordinates '
      do 600 ic=1,natoms
         write(6,'(1x,i5,5x,3f14.8)') ic,(xyz(k,ic),k=1,3)
 600  continue
      write(6,*) ' transformation matrices '
      do 700 it=1,ntrans
         write(6,'(1x,i5)') it
         write(6,'(1x,9f8.4)') (trans(k,it),k=1,9)
 700  continue
      error STOP 'fatal error in nucex.'
      end
c RCS ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Id: getgen.f,v 1.1 1992/09/10 13:40:37 tomjones Exp $
c $Log: getgen.f,v $
c Revision 1.1  1992/09/10 13:40:37  tomjones
c Initial revision
c
c RCS ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c--------------------------------------------------------------------
      subroutine getgen (csf,nn,cla,gen,ngen,ierr,nprisy)
      use iso_fortran_env, only : wp => real64
      implicit real(wp)(a-h,o-y)
      implicit character (z)
      dimension gen(9,3),unit(9)
      character csf,cla
      data zc,zd,zs,zt,zo,zi,zv,zh/'c','d','s','t','o','i','v','h'/
      a0=0
      a1=1
      a2=2
      a3=3
      a4=4
      a5=5
      ierr=0
      pi=a4*atan(a1)
      twopi=pi+pi
c       generate unit matrix
      do 10 i=2,8
  10  unit(i) = a0
      unit(1) = a1
      unit(5) = a1
      unit(9) = a1
      ngen=0
c     write(6,*) ' this program sets generators for point groups'
c     write(6,*) ' specified by schoenflies symbols, e.g. d2d'
c     write(6,*) ' please use s1 for cs  and s2 for ci'
c     write(6,*) ' main axis is always the  z - axis'
c     write(6,*) ' secondary axis is always x, e.g. for d3'
c     write(6,*) ' sigma(v) plane is always xz, e.g. for c3v'
c     write(6,*) ' input  schoenflies symbol: c2v  or td etc'
      if (nprisy.ge.0) write(6,20)
   20 format(/,' the group has the following generators :')
c     note that in ih the mirror plane is perpendicular to the c2-axis
c     but not to the c5- or c3-axis
      if (csf.eq.zi.and.cla.eq.zh) cla=zv
c
c       set generators for input group
c
c       cn(z)
c
      n=0
      if (nn.gt.0.and.(csf.eq.zc.or.csf.eq.zd).and.cla.ne.zd) n=nn
      if (csf.eq.zt.and.cla.ne.zd) n=2
      if (csf.eq.zo) n=4
      if(csf.eq.zi) n=5
      if (n.le.0) go to 250
      if (n.le.9 .and. nprisy.ge.0) write(6,30) n
      if (n.gt.9 .and. nprisy.ge.0) write(6,31) n
   30 format(3x,'c',i1,'(z)')
   31 format(3x,'c',i2,'(z)')
      ngen = ngen +1
      if (ngen.gt.3) go to 1000
      do 200 i=1,9
  200 gen(i,ngen) = unit(i)
      an = n
      alpha = twopi/an
      x = sin(alpha)
      y = cos(alpha)
      gen(1,ngen) = y
      gen(5,ngen) = y
      gen(4,ngen) = -x
      gen(2,ngen) = x
  250 continue
c
c       sn(z)
c
      n=0
      if (csf.eq.zs) n=nn
      if (csf.eq.cla) n=2*nn
      if (csf.eq.zt.and.cla.eq.zd) n=4
      if (n.le.0) go to 400
      if (n.le.9 .and. nprisy.ge.0) write(6,300) n
      if (n.gt.9 .and. nprisy.ge.0) write(6,301) n
  300 format(3x,'s',i1,'(z)')
  301 format(3x,'s',i2,'(z)')
      ngen = ngen +1
      if (ngen.gt.3) go to 1000
      do 350 i=1,9
  350 gen(i,ngen) =-unit(i)
      an = n
      alpha = twopi/an
      x = sin(alpha)
      y = cos(alpha)
      gen(1,ngen) = y
      gen(5,ngen) = y
      gen(4,ngen) = -x
      gen(2,ngen) = x
  400 continue
c
c        c3(1,1,1)  for t,td,th,o,oh
c
      if (csf.ne.zt.and.csf.ne.zo) go to 450
      if (nprisy.ge.0) write(6,410)
  410 format(3x,'c3(1,1,1)')
      ngen = ngen +1
      if (ngen.gt.3) go to 1000
      do 420 i=1,9
  420 gen(i,ngen) = a0
      gen(2,ngen) = a1
      gen(6,ngen) = a1
      gen(7,ngen) = a1
  450 continue
c
c       c3 for i or ih
c
      if (csf.ne.zi) go to 500
      if (nprisy.ge.0) write(6,480)
  480 format(3x,'c3 for i or ih')
      ngen = ngen + 1
      if (ngen.gt.3) go to 1000
      sqr3 = sqrt(a3)
      cosd = sin(a2*pi/a5)/((a1-cos(a2*pi/a5))*sqr3)
      sind = sqrt(a1-cosd**2)
      gen(1,ngen) = a1 - a3*cosd**2/a2
      gen(2,ngen) = sqr3*cosd/a2
      gen(3,ngen) = a3*sind*cosd/a2
      gen(4,ngen) = -gen(2,ngen)
      gen(5,ngen) = -a1/a2
      gen(6,ngen) = sqr3*sind/a2
      gen(7,ngen) = gen(3,ngen)
      gen(8,ngen) = -gen(6,ngen)
      gen(9,ngen) = a1-a3*sind**2/a2
  500 continue
c
c       c2(x)
c
      if(csf.ne.zd) go to 600
      if (nprisy.ge.0) write(6,520)
  520 format(3x,'c2(x)')
      i=1
      ngen = ngen +1
      if (ngen.gt.3) go to 1000
      do 550 j=1,9
  550 gen(j,ngen) = -unit(j)
      gen(4*i-3,ngen) = a1
  600 continue
c
c       mirror plane sigma(xz) or sigma(xy)
c
      if (cla.ne.zv.and.cla.ne.zh) go to 750
      if (cla.eq.zv) then
         if (nprisy.ge.0) write(6,610)
  610    format(3x,'mirror plane sigma(xz)')
         i=2
      else
         if (nprisy.ge.0) write(6,620)
  620    format(3x,'mirror plane sigma(xy)')
         i=3
      endif
      ngen = ngen + 1
      if (ngen.gt.3) go to 1000
      do 700 j=1,9
  700 gen(j,ngen) = unit(j)
      gen(4*i-3,ngen) = -a1
  750 return
 1000 write(6,*) ' more than ',3,' generators'
      ierr=1
      end


      subroutine symmetry(natoms,nat3,xyz0,h,eig,totsym,nvar)
      use iso_fortran_env, only : wp => real64
      implicit none
      integer natoms,nat3,totsym(nat3),nvar
      real(wp) h(nat3,nat3),xyz0(3,natoms),eig(nat3)
      real(wp) :: xyz(3,natoms)

      real(wp),allocatable :: trans(:,:)             
      real(wp) thrsym
      integer ntrans,i,j,nat,it,k,m
      logical lok
      integer :: ifile

      call open_binary(ifile,'trans','r')
      read(ifile)ntrans
      allocate(trans(9,ntrans))
      do it=1,ntrans
         read(ifile) (trans(i,it),i=1,9)
      enddo   
      call close_file(ifile)

      thrsym=1.d-6*natoms
      totsym = 0
      nvar = 0
      do i=1,nat3
         xyz = xyz0
         m = 0
         do j=1,natoms
            do k=1,3
               m = m + 1
               xyz(k,j)=xyz(k,j)+h(m,i)*0.005
            enddo
         enddo
         call checksym(thrsym,natoms,ntrans,trans,xyz,lok)
         if(lok.and.abs(eig(i)).gt.1.d-8)then
            totsym(i)=1
            nvar = nvar + 1
         endif
      enddo
      write(*,*)'Number of symmetry operations ',ntrans
      write(*,*)'Nvar, symmetry restricted     ',nvar

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine grdsym(grd,natoms)
      use iso_fortran_env, only : wp => real64
      use symparam
c,ntrans,ict,trans)

      implicit real(wp) (a-h,o-z)
      dimension grd(3,natoms)
c     dimension ict(natoms,ntrans),trans(9,ntrans)
      data garnix/1.d-14/, zero/0.d0/, one/1.d0/
      real(wp) :: det(3,natoms),scr(3)

c     ----- symmetrize cartesian gradient/displ vector -----
c     using transformation matrices and nuclear exchange table ict
c     ict = transformation table "atoms versus symmetry operations"

      if (ntrans.gt.1) then
         do 200 iat = 1,natoms
            det(1,iat) = zero
            det(2,iat) = zero
            det(3,iat) = zero
  200    continue
         do 300 iat=1,natoms
            do 400 itrans=1,ntrans
               call mult3(scr,trans(1,itrans),grd(1,iat),1)
               isymat=ict(iat,itrans)
               det(1,isymat)=det(1,isymat) + scr(1)
               det(2,isymat)=det(2,isymat) + scr(2)
               det(3,isymat)=det(3,isymat) + scr(3)
  400       continue
  300    continue
         ogroup=one/dble(ntrans)
         do 500 iat=1,natoms
            grd(1,iat)=det(1,iat)*ogroup
            grd(2,iat)=det(2,iat)*ogroup
            grd(3,iat)=det(3,iat)*ogroup
  500    continue
      end if

c     brush away numerical zeros
      do 575 iat=1,natoms
      do 575 j=1,3
      if (dabs(grd(j,iat)).lt.garnix) grd(j,iat) = zero
  575 continue

      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine checksym(thrsym,natoms,ntrans,trans,xyz,same)
      use iso_fortran_env, only : wp => real64
      implicit none
      real(wp) xyz(3,natoms),trans(9,ntrans)
      integer natoms,ntrans
      integer ic,it,jc
      logical same
      real(wp) scr(3),thrsym
      real(wp) dif

      same=.true.
      do 300 ic=1,natoms
         do 250 it=2,ntrans
            call mult3 (scr,trans(1,it),xyz(1,ic),1)
            do 200 jc=1,natoms
               if (dif(scr,xyz(1,jc),3).gt.thrsym) goto 200
               goto 250
 200        continue
            same = .false.
            return
 250     continue
 300  continue

      end

c ---------------------------------------------------------------------

      function dif(a,b,n)
      use iso_fortran_env, only : wp => real64
      implicit real(wp) (a-h,o-z)
      dimension a(n),b(n)
      s = 0.d0
      do 10 k=1,n
  10  s = s + (b(k)-a(k))*(b(k)-a(k))
      dif = dsqrt(s)
      return
      end

c ---------------------------------------------------------------------

      subroutine mult3 (a,b,c,n)
      use iso_fortran_env, only : wp => real64
      implicit real(wp) (a-h,o-z)
c ---------------------------------------------------------------------
c     multiplication of two matrices b,c
c      a = b*c     dim: a(3,n),b(3,3),c(3,n)
c ---------------------------------------------------------------------
      dimension a(3,n),b(3,3),c(3,n)
      do 30 j=1,n
      do 20 i=1,3
      s = 0.d0
      do 10 k=1,3
  10  s = s + b(i,k)*c(k,j)
  20  a(i,j) = s
  30  continue

      return
      end

