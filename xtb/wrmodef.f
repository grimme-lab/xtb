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

      subroutine wrmodef(typ,n,at,xyzin,wbo,rmass,freq,u,udum,vthr,
     &                   linear)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in)    :: typ,n,at(n)
      real(wp),intent(in)    :: u(3*n,3*n),freq(3*n)
      real(wp),intent(in)    :: xyzin(3,n),wbo(n,n),rmass(3*n)
      real(wp),intent(inout) :: udum(3*n,3*n)
      real(wp),intent(in)    :: vthr
      logical, intent(in)    :: linear

      real(wp) :: one
      real(wp) :: norm,w
      real(wp),allocatable :: bmat(:,:),xyz(:,:),ttyp(:),rdum(:)
      real(wp),allocatable :: tmp(:,:),pop(:),deloc(:),xtyp(:,:)
      real(wp),allocatable :: dit (:),geo(:,:),vtyp(:),fdum(:),tors(:)
      parameter (one=1.0d0)
      integer, allocatable :: ind(:)
      integer, allocatable :: na(:),nb(:),nc(:)
      integer i,k,j,ik,kl,l,root,lend,b1,b2,n3
      integer ich ! file handle
      integer :: nfree,nnull

      if (linear) then
         nnull = 5
      else
         nnull = 6
      endif
      nfree = 3*n-nnull

      allocate( bmat(nfree,3*n), xyz(3,n), ttyp(3*n), rdum(3*n), 
     &          tmp(3*n,3*n), pop(n), deloc(3*n), xtyp(3*n,3),
     &          dit(nfree), geo(3,n), vtyp(3), fdum(3*n), tors(3*n))
      allocate( ind(3*n), na(n), nb(n), nc(n) )

      geo=0
      ttyp(1:nnull)=0
      fdum(1:nnull)=0
      rdum(1:nnull)=0
      xyz=xyzin*0.52916790d0
      n3=3*n

      call XYZINT(XYZ,n,NA,NB,NC,one,GEO)

      call bzmat(n,at,xyz,bmat)

      do root=nnull+1,n3

         dit = 0
         do k=1,n3
            do j=1,3*n-nnull
               dit(j)=dit(j)+bmat(j,k)*u(k,root)
            enddo
         enddo

         do k=1,n
            norm=0
            do ik=1,3
               norm=norm+u((k-1)*3+ik,root)**2
            enddo
            pop(k)=norm
         enddo

         norm=0
         do k=1,n
            norm=norm+pop(k)**2
         enddo
         deloc(root)=1.d0/(norm+1.d-8)

         vtyp = 0
         kl=0
         do k=2,n
            lend=3
            if(k.eq.2) lend=1
            if(k.eq.3) lend=2
            do l=1,lend
               kl=kl+1
c              b1=nc(k)
c              b2=nb(k)
c              w=1.
c              if(l.eq.3.and.dit(kl).gt.0.1)then
c                write(*,*) b1,b2,wbo(b1,b2)
c                if(b1.gt.0.and.b2.gt.0.and.wbo(b1,b2).gt.0.5)
c    .           w=1.0d0/wbo(b1,b2)**8
c              endif
               vtyp(l)=vtyp(l)+dit(kl)**2
            enddo
         enddo

         norm=sum(vtyp)
         vtyp=vtyp/norm

c        write(*,'(i5,2f10.2,5x,4f7.2)')
c    .   root,freq(root),rmass(root),vtyp(1:3)

         ttyp(root)=vtyp(3)**2
         tors(root)=vtyp(3)
      enddo

      write(*,*)
      write(*,*) 'normal mode %tage of torsion character'
      write(*,'(8(i4,'':'',f6.2))') (i,ttyp(i),i=1,n3)

      do root=nnull+1,n3
         if(freq(root).gt.0) then
            ttyp(root)=freq(root)+(1.-ttyp(root))*freq(root)*2.
         else
            ttyp(root)=freq(root)
         endif
      enddo
      if(typ.eq.1)then
      do root=nnull+1,n3
         if(freq(root).gt.0) then
            ttyp(root)=ttyp(root)+10.*deloc(root)**2
         else
            ttyp(root)=freq(root)
         endif
      enddo
      endif
      do i=1,nnull
         ttyp(i)=-1.d+42
      enddo

c sort modes according to typ (and freq)
      do i=1,n3
         ind(i)=i
      enddo
      call Qsort(ttyp, 1, n3, ind)

      j=0
      do i=nnull+1,n3
         udum(1:n3,i)=u(1:n3,ind(i))
         fdum(i)=freq(ind(i))
         rdum(i)=rmass(ind(i))
c        write(*,'(i5,2f10.2,5x,4f7.2)')
c    .         i,fdum(i),tors(ind(i)),rdum(i)
         if(fdum(i).lt.vthr.and.tors(ind(i)).gt.0.5.and.i.gt.j) j=i
c        write(*,'(i4,5f10.2)') i,fdum(i),deloc(ind(i))!,xtyp(ind(i),1:3)
      enddo
      write(*,*)
     &   'recommended # of modes for mode following',max(j-nnull,4)


      if(typ.eq.0)then
         call open_file(ich,'tmpxx','w')
         write(ich,*) max(j-6,4)
         call close_file(ich)
      endif


      if(typ.eq.0)then
         call open_binary(ich,'xtb_normalmodes','w')
         call g98fake('g98_canmode.out',n,at,xyzin,fdum,udum,tmp)
      else
         call open_file(ich,'xtb_localmodes','w')
         call g98fake('g98_locmode.out',n,at,xyzin,fdum,udum,tmp)
      endif
      write(ich)n3
      write(ich)fdum
      write(ich)rdum
      write(ich)udum
      call close_file(ich)

      end subroutine wrmodef

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine g98fake(fname,n,at,xyz,freq,u2,u)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in) :: n,at(n)
      real(wp),intent(in) :: freq(3*n),xyz(3,n),u2(3*n,3*n)
      real(wp),intent(inout) :: u(3*n,3*n)
      character(len=*),intent(in) :: fname

      integer gu,i,j,ka,kb,kc,la,lb,k
      character(len=2) irrep
      real(wp),allocatable :: red_mass(:)
      real(wp),allocatable :: force   (:)
      real(wp),allocatable :: ir_int  (:)
      real(wp),allocatable :: f2      (:)
      real(wp),allocatable :: zero
      allocate( red_mass(3*n),
     &          force   (3*n),
     &          ir_int  (3*n),
     &          f2      (3*n) )

      irrep='a'
      red_mass=99.0
      force   =99.0
      ir_int  =99.0
      zero    =0.0

      k=0
      do i=1,3*n
         if(abs(freq(i)).gt.1.d-1)then
            k=k+1
            u(1:3*n,k)=u2(1:3*n,i)
            f2(k)=freq(i)
         endif
      enddo

      gu=55
      call open_file(gu,fname,'w')
      write (gu,'('' Entering Gaussian System'')')
      write (gu,'('' *********************************************'')')
      write (gu,'('' Gaussian 98:'')')
      write (gu,'('' frequency output generated by the xtb code'')')
      write (gu,'('' *********************************************'')')

      write (gu,*) '                        Standard orientation:'
      write (gu,*) '---------------------------------------------',
     .             '-----------------------'
      write (gu,*) ' Center     Atomic     Atomic',
     .             '              Coordinates (Angstroms)'
      write (gu,*) ' Number     Number      Type ',
     .             '             X           Y           Z'
      write (gu,*) '-----------------------',
     .             '---------------------------------------------'
      j=0
      do i=1,n
       write(gu,111) i,at(i),j,xyz(1:3,i)*0.52917726
      enddo
      write (gu,*) '----------------------',
     .              '----------------------------------------------'
      write (gu,*) '    1 basis functions        1 primitive gaussians'
      write (gu,*) '    1 alpha electrons        1 beta electrons'
      write (gu,*)
111   format(i5,i11,i14,4x,3f12.6)

      write (gu,*) 'Harmonic frequencies (cm**-1), IR intensities',
     .              ' (KM/Mole),'
      write (gu,*) 'Raman scattering activities (A**4/amu),',
     .              ' Raman depolarization ratios,'
      write (gu,*) 'reduced masses (AMU), force constants (mDyne/A)',
     .             ' and normal coordinates:'

      ka =1
      kc=3
60    kb=min0(kc,k)
      write (gu,100) (j,j=ka,kb)
      write (gu,105) (irrep,j=ka,kb)
      write (gu,110) ' Frequencies --',(f2(j),j=ka,kb)
      write (gu,110) ' Red. masses --',(red_mass(j),j=ka,kb)
      write (gu,110) ' Frc consts  --',(force(j),j=ka,kb)
      write (gu,110) ' IR Inten    --',(ir_int(j),j=ka,kb)
      write (gu,110) ' Raman Activ --',(zero,j=ka,kb)
      write (gu,110) ' Depolar     --',(zero,j=ka,kb)
      write (gu,*)'Atom AN      X      Y      Z        X      Y',
     .            '      Z        X      Y      Z'
      la=1
70    lb=n
      do  i=la,lb
        write (gu,130) i,at(i),
     .              (u(i*3-2,j),
     .               u(i*3-1,j),
     .               u(i*3  ,j),j=ka,kb)
      enddo
      if (lb.eq.n) go to 90
      go to 70
90    if (kb.eq.k) then
       return
      endif
      ka=kc+1
      kc=kc+3
      go to 60

100   format (3(20x,i3))
105   format (3x,3(18x,a5))
110   format (a15,f11.4,12x,f11.4,12x,f11.4)
130   format (2i4,3(2x,3f7.2))

      call close_file(gu)
      return

      end subroutine g98fake

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
