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

!******************************************
! rotate rigid fragment in binding "pocket"
! of other fragment
!****************************************** 
      subroutine pocketscan(n,at,xyz,nout) 
      use splitparam
      use setparam
      use ncoord, only : ncoord_d3
      implicit none
      integer n,at(n),i,j,k,ios,ierr,npath,nout
      real*8 xyz(3,n),cn0(n),cn(n),coord(3,n),f,dum,cthr
      real*8 cma1(3),cma2(3),trfm1(3,3),trfm2(3,3)
      real*8, allocatable :: xyzpath(:,:,:)
      character*2 dummy,asym
      integer :: ich

      npath=mode_nscan
      if(mod(npath,2).eq.0) npath=npath+1 ! make it odd
      ierr=0
      nout=0
      f=0.5291770d0
      cn0=0.0d0
      cthr=1.0d0
!     call ncoord2(.false.,n,at,xyz,cn0,500.0d0)
      call ncoord_d3(n,at,xyz,cn0,500.0d0)
      ! obtain cma's and principal axes of 1 and 2
      call ifaceaxis(n,at,xyz,cma1,trfm1,cma2,trfm2)         

      allocate(xyzpath(3,n,3*npath), stat=ios)
      if(ios.ne.0) stop 'error in mode 0 alloc'
      xyzpath=0.0d0
      !***************************************
      ! now rotate fragment 1 around axes of 1 
      !***************************************
      call pocketrotation(n,at,xyz,cma1,trfm1,1,npath,xyzpath)

      call open_file(ich,'xtb_modescan_0.xyz','w')
      ! check for clashes and do reoptimization
      do k=1,3*npath
        coord(1:3,1:n)=xyzpath(1:3,1:n,k)
        cn=0.0d0
!       call ncoord2(.false.,n,at,coord,cn,500.0d0)
        call ncoord_d3(n,at,coord,cn,500.0d0)
        dum=-1.0d0
        do i=1,n
           if(abs(cn(i)-cn0(i)).gt.dum) dum=abs(cn(i)-cn0(i)) ! clash check
        enddo
        if(dum.gt.cthr) then
           ierr=ierr+1
           cycle
        endif
        write(ich,*)n
        write(ich,*)
        do i=1,n
           write(ich,'(1x,a2,1x,3f18.8)')
     .     asym(at(i)),xyzpath(1:3,i,k)*f
        enddo
      enddo

      deallocate(xyzpath)

      ! ***
      allocate(xyzpath(3,n,3*npath), stat=ios)
      if(ios.ne.0) stop 'error in mode 0 alloc'
      xyzpath=0.0d0

      !***************************************
      ! now rotate fragment 2 around axes of 2
      !***************************************
      call pocketrotation(n,at,xyz,cma2,trfm2,2,npath,xyzpath)
      ! check for clashes and do reoptimization
      do k=1,3*npath
        coord(1:3,1:n)=xyzpath(1:3,1:n,k)
        cn=0.0d0
!       call ncoord2(.false.,n,at,coord,cn,500.0d0)
        call ncoord_d3(n,at,coord,cn,500.0d0)
        dum=-1.0d0
        do i=1,n
           if(abs(cn(i)-cn0(i)).gt.dum) dum=abs(cn(i)-cn0(i)) ! clash check
        enddo
        if(dum.gt.cthr) then
           ierr=ierr+1
           cycle
        endif
        write(ich,*)n
        write(ich,*)
        do i=1,n
           write(ich,'(1x,a2,1x,3f18.8)')
     .     asym(at(i)),xyzpath(1:3,i,k)*f
        enddo
      enddo
      call close_file(ich)

      deallocate(xyzpath)


      nout=6*npath - ierr

      return 
      end

      subroutine pocketrotation(n,at,xyz,cma,trfm,ifrag,npath,xyzpath)  
      use splitparam
      implicit none
      real*8 xyz(3,n),coord(3,n),turncrd(3,n),addcrd(3,n),turncrd0(3,n)
      integer n,at(n),ifrag,npath
      real*8 cma(3),trfm(3,3),rotm(3,3),dblepi,step,f
      real*8 xyzpath(3,n,3*npath)
      integer i,j,k
      character*2 asym

      f=0.5291770d0
      dblepi =2.0d0*3.14159265358979D0
      if(npath.lt.1) then
        step=dblepi/12.0d0
        npath=11
      else 
        step=dblepi/dble(npath+1)
      endif
 
      coord=0.0d0
      turncrd=0.0d0
      turncrd0=0.0d0
      addcrd=0.0d0
      do i=1,n
         if(splitlist(i).eq.ifrag) then 
           turncrd0(1:3,i)=xyz(1:3,i)-cma(1:3) ! part that is rotated
         else
           addcrd(1:3,i)=xyz(1:3,i)-cma(1:3) ! part that is fixed and simply added to the zeros of turncrd to yield the new coord
           ! cma is added here, simply to have a simpler addition later on
         endif
      enddo
      turncrd0=matmul(trfm,turncrd0)
      !*******************
      ! rotation around x
      !*******************
      turncrd=turncrd0
      rotm(1,1)=1.0d0
      rotm(2,1)=0.0d0
      rotm(3,1)=0.0d0
      rotm(1,2)=0.0d0
      rotm(2,2)=cos(step)
      rotm(3,2)=sin(step)
      rotm(1,3)=0.0d0
      rotm(2,3)=-sin(step)
      rotm(3,3)=cos(step)
      j=0
      ! now rotate stepwise 
      do k=1,npath
         j=j+1
         coord=0.0d0
         turncrd=matmul(rotm,turncrd)
         coord=matmul(transpose(trfm),turncrd)
         do i=1,n
            xyzpath(1:3,i,j)=coord(1:3,i)+cma(1:3)+addcrd(1:3,i)
         enddo
      enddo

      !*******************
      ! rotation around y
      !*******************
      turncrd=turncrd0
      rotm(1,1)=cos(step)
      rotm(2,1)=0.0d0
      rotm(3,1)=sin(step)
      rotm(1,2)=0.0d0
      rotm(2,2)=1.0d0
      rotm(3,2)=0.0d0
      rotm(1,3)=-sin(step)
      rotm(2,3)=0.0d0
      rotm(3,3)=cos(step)

      ! now rotate stepwise 
      do k=1,npath
         j=j+1
         coord=0.0d0
         turncrd=matmul(rotm,turncrd)
         coord=matmul(transpose(trfm),turncrd)
         do i=1,n
            xyzpath(1:3,i,j)=coord(1:3,i)+cma(1:3)+addcrd(1:3,i)
         enddo
      enddo

      !*******************
      ! rotation around z
      !*******************
      turncrd=turncrd0
      rotm(1,1)=cos(step)
      rotm(2,1)=sin(step)
      rotm(3,1)=0.0d0
      rotm(1,2)=-sin(step)
      rotm(2,2)=cos(step)
      rotm(3,2)=0.0d0
      rotm(1,3)=0.0d0
      rotm(2,3)=0.0d0
      rotm(3,3)=1.0d0

      ! now rotate stepwise 
      do k=1,npath
         j=j+1
         coord=0.0d0
         turncrd=matmul(rotm,turncrd)
         coord=matmul(transpose(trfm),turncrd)
         do i=1,n
            xyzpath(1:3,i,j)=coord(1:3,i)+cma(1:3)+addcrd(1:3,i)
         enddo
      enddo


      return
      end     


      !*********************************************
      ! principal rotation axes for pocket rotation
      !*********************************************
      subroutine ifaceaxis(numat,nat,coord,cma1,trfm1,cma2,trfm2)    
      use splitparam
      implicit none
      real*8 coord(3,numat)
      integer numat,nat(numat)
      real*8 cma1(3),cma2(3),trfm1(3,3),trfm2(3,3)
      real*8 t1(6),t2(6),eig(3),x(numat),y(numat),z(numat)
      real*8 sum1w,sum1wx,sum1wy,sum1wz,atmss,eps
      real*8 sum2w,sum2wx,sum2wy,sum2wz,xsum
      integer i,j,k

      sum1w=1.d-20                                                  
      sum1wx=0.d0                                                   
      sum1wy=0.d0                                                   
      sum1wz=0.d0                                                   
      sum2w=1.d-20                                                  
      sum2wx=0.d0                                                   
      sum2wy=0.d0                                                   
      sum2wz=0.d0                                                   

      do i=1,numat         
         atmss=atmass(i)                                  
         if(splitlist(i).eq.1)then
           sum1w =sum1w+atmss                                       
           sum1wx=sum1wx+atmss*coord(1,i)                          
           sum1wy=sum1wy+atmss*coord(2,i)                          
           sum1wz=sum1wz+atmss*coord(3,i)     
         else 
           sum2w =sum2w+atmss             
           sum2wx=sum2wx+atmss*coord(1,i) 
           sum2wy=sum2wy+atmss*coord(2,i) 
           sum2wz=sum2wz+atmss*coord(3,i)  
         endif             
      enddo                                                  
                                                                   
      eps=1.d-3
      sum1wx=sum1wx/sum1w                                             
      sum1wy=sum1wy/sum1w                                             
      sum1wz=sum1wz/sum1w                                             
      sum2wx=sum2wx/sum2w                                             
      sum2wy=sum2wy/sum2w                                             
      sum2wz=sum2wz/sum2w                                             
      cma1=0.0d0
      cma2=0.0d0
      cma1(1)=sum1wx
      cma1(2)=sum1wy
      cma1(3)=sum1wz
      cma2(1)=sum2wx
      cma2(2)=sum2wy
      cma2(3)=sum2wz

      ! shift to respective CMA
      do i=1,numat 
         if(splitlist(i).eq.1)then                              
             x(i)=coord(1,i)-sum1wx 
             y(i)=coord(2,i)-sum1wy 
             z(i)=coord(3,i)-sum1wz 
         else
             x(i)=coord(1,i)-sum2wx
             y(i)=coord(2,i)-sum2wy
             z(i)=coord(3,i)-sum2wz
         endif
      enddo

      ! initialize moment of inertia tensors      
      do i=1,6                                                  
         t1(i)=dble(i)*1.0d-10                                         
         t2(i)=dble(i)*1.0d-10
      enddo
                                                                   
      do i=1,numat                                           
         atmss=atmass(i)         
         if(splitlist(i).eq.1)then                            
            t1(1)=t1(1)+atmss*(y(i)**2+z(i)**2)+eps                 
            t1(2)=t1(2)-atmss*x(i)*y(i)
            t1(3)=t1(3)+atmss*(z(i)**2+x(i)**2)+eps                 
            t1(4)=t1(4)-atmss*z(i)*x(i)
            t1(5)=t1(5)-atmss*y(i)*z(i)
            t1(6)=t1(6)+atmss*(x(i)**2+y(i)**2)+eps                
         else
            t2(1)=t2(1)+atmss*(y(i)**2+z(i)**2)+eps
            t2(2)=t2(2)-atmss*x(i)*y(i)
            t2(3)=t2(3)+atmss*(z(i)**2+x(i)**2)+eps
            t2(4)=t2(4)-atmss*z(i)*x(i)
            t2(5)=t2(5)-atmss*y(i)*z(i)
            t2(6)=t2(6)+atmss*(x(i)**2+y(i)**2)+eps
         endif 
      enddo                                                


      trfm1=0.0d0                                   
      trfm2=0.0d0
      eig=0.0d0
      ! solve eigenvalue problem for iface 1
      call rsp(t1,3,3,eig,trfm1)
      eig=0.0d0
      ! solve eigenvalue problem for iface 2
      call rsp(t2,3,3,eig,trfm2)

                                     

c   check for right-handedness and flip otherwise
      xsum=0.0d0
      xsum=trfm1(1,1)*(trfm1(2,2)*trfm1(3,3)-trfm1(3,2)*trfm1(2,3)) +   
     1     trfm1(1,2)*(trfm1(2,3)*trfm1(3,1)-trfm1(2,1)*trfm1(3,3)) +   
     2     trfm1(1,3)*(trfm1(2,1)*trfm1(3,2)-trfm1(2,2)*trfm1(3,1))     
      if( xsum .lt. 0.0d0) then                                         
         do j=1,3                                               
           trfm1(j,1)=-trfm1(j,1)                                      
         enddo
      endif                                                        

      xsum=0.0d0
      xsum=trfm2(1,1)*(trfm2(2,2)*trfm2(3,3)-trfm2(3,2)*trfm2(2,3)) +   
     1     trfm2(1,2)*(trfm2(2,3)*trfm2(3,1)-trfm2(2,1)*trfm2(3,3)) +   
     2     trfm2(1,3)*(trfm2(2,1)*trfm2(3,2)-trfm2(2,2)*trfm2(3,1))     
      if( xsum .lt. 0.0d0) then                                         
         do j=1,3                                               
           trfm2(j,1)=-trfm2(j,1)                                      
         enddo
      endif                                                        

                                                                   
      return               
      end                 
