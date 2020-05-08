
      subroutine mmass(n,at,mass)                
      implicit none 
      real*8 mass
      integer n,at(n)
      integer i
      real*8 ams(107),amutokg
      data  ams /  1.00790d0,  4.00260d0,  6.94000d0,  9.01218d0,
     .10.81000d0, 12.01100d0, 14.00670d0, 15.99940d0, 18.99840d0,
     .20.17900d0, 22.98977d0, 24.30500d0, 26.98154d0, 28.08550d0,
     .30.97376d0, 32.06000d0, 35.45300d0, 39.94800d0, 39.09830d0,
     .40.08000d0, 44.95590d0, 47.90000d0, 50.94150d0, 51.99600d0,
     .54.93800d0, 55.84700d0, 58.93320d0, 58.71000d0, 63.54600d0,
     .65.38000d0, 69.73500d0, 72.59000d0, 74.92160d0, 78.96000d0,
     .79.90400d0, 83.80000d0, 85.46780d0, 87.62000d0, 88.90590d0,
     .91.22000d0, 92.90640d0, 95.94000d0, 98.90620d0, 101.0700d0,
     .102.9055d0, 106.4000d0, 107.8680d0, 112.4100d0, 114.8200d0,
     .118.6900d0, 121.7500d0, 127.6000d0, 126.9045d0, 131.3000d0,
     .132.9054d0, 137.3300d0, 
     .138.91,140.12,140.91,144.24,147.00,150.36,151.97,157.25,
     .158.93,162.50,164.93,167.26,168.93,173.04,174.97,
     .178.4900d0, 180.9479d0,
     .183.8500d0, 186.2070d0, 190.2000d0, 192.2200d0, 195.0900d0,
     .196.9665d0, 200.5900d0, 204.3700d0, 207.2000d0, 208.9804d0,
     .209.,210.,222.,21*0.000d0/
      parameter (amutokg=1.66053886E-27)
      mass=0
      do i=1,n
         mass=mass+ams(at(i))
      enddo
      mass=mass*amutokg
      end

      subroutine getvolume(n,at,xyz,v)
      implicit none
      real*8 xyz(3,n),v
      integer n,at(n)

      real*8 dum(8),rx,ry,rz 
      real x,y,z,rr,r
      character*2 asym
      integer i,j

!     dum=0
!     do i=1,n
!        rx=xyz(1,i)
!        ry=xyz(2,i)
!        rz=xyz(3,i)
!        rr=sqrt(rx**2+ry**2+rz**2)         
!        if(rx.gt.0.and.ry.gt.0.and.rz.gt.0.and.rr.gt.dum(1))dum(1)=rr  ! + + +
!        if(rx.le.0.and.ry.gt.0.and.rz.gt.0.and.rr.gt.dum(2))dum(2)=rr  ! - + +
!        if(rx.gt.0.and.ry.le.0.and.rz.gt.0.and.rr.gt.dum(3))dum(3)=rr  ! + - +
!        if(rx.gt.0.and.ry.gt.0.and.rz.le.0.and.rr.gt.dum(4))dum(4)=rr  ! + + -
!        if(rx.le.0.and.ry.le.0.and.rz.gt.0.and.rr.gt.dum(5))dum(5)=rr  ! - - +
!        if(rx.le.0.and.ry.gt.0.and.rz.le.0.and.rr.gt.dum(6))dum(6)=rr  ! - + -
!        if(rx.gt.0.and.ry.le.0.and.rz.le.0.and.rr.gt.dum(7))dum(7)=rr  ! + - -
!        if(rx.le.0.and.ry.le.0.and.rz.le.0.and.rr.gt.dum(8))dum(8)=rr  ! - - -
!     enddo
      rx=maxval(xyz(1,:))-minval(xyz(1,:))
      ry=maxval(xyz(2,:))-minval(xyz(2,:))
      rz=maxval(xyz(3,:))-minval(xyz(3,:))

!     r=sum(dum)/8.0d0+1.0d0
      v=(rx+0.)*(ry+0.)*(rz+0.)*(0.52917726D-10**3) !*4.*3.14159265/3.)

      end

      subroutine axis(numat,nat,xyz,aa,bb,cc)        
      implicit double precision (a-h,o-z)                        
      dimension xyz(3,*), ams(107)                                 
      integer nat(*)

      data  ams /  1.00790d0,  4.00260d0,  6.94000d0,  9.01218d0,
     .10.81000d0, 12.01100d0, 14.00670d0, 15.99940d0, 18.99840d0,
     .20.17900d0, 22.98977d0, 24.30500d0, 26.98154d0, 28.08550d0,
     .30.97376d0, 32.06000d0, 35.45300d0, 39.94800d0, 39.09830d0,
     .40.08000d0, 44.95590d0, 47.90000d0, 50.94150d0, 51.99600d0,
     .54.93800d0, 55.84700d0, 58.93320d0, 58.71000d0, 63.54600d0,
     .65.38000d0, 69.73500d0, 72.59000d0, 74.92160d0, 78.96000d0,
     .79.90400d0, 83.80000d0, 85.46780d0, 87.62000d0, 88.90590d0,
     .91.22000d0, 92.90640d0, 95.94000d0, 98.90620d0, 101.0700d0,
     .102.9055d0, 106.4000d0, 107.8680d0, 112.4100d0, 114.8200d0,
     .118.6900d0, 121.7500d0, 127.6000d0, 126.9045d0, 131.3000d0,
     .132.9054d0, 137.3300d0, 
     .138.91,140.12,140.91,144.24,147.00,150.36,151.97,157.25,
     .158.93,162.50,164.93,167.26,168.93,173.04,174.97,
     .178.4900d0, 180.9479d0,
     .183.8500d0, 186.2070d0, 190.2000d0, 192.2200d0, 195.0900d0,
     .196.9665d0, 200.5900d0, 204.3700d0, 207.2000d0, 208.9804d0,
     .209.,210.,222.,21*0.000d0/
                                                                                
      PARAMETER (BOHR=0.52917726)                                               
      dimension t(6), rot(3), xyzmom(3), eig(3), evec(3,3)                      
      dimension x(numat),y(numat),z(numat),coord(3,numat)
      data t /6*0.d0/                                                           
************************************************************************        
*     const1 =  10**40/(n*a*a)                                                  
*               n = avergadro's number                                          
*               a = cm in an angstrom                                           
*               10**40 is to allow units to be 10**(-40)gram-cm**2              
*                                                                               
************************************************************************        
      const1 = 1.66053d0                                                        
************************************************************************        
*                                                                               
*     const2 = conversion factor from angstrom-amu to cm**(-1)                  
*                                                                               
*            = (planck's constant*n*10**16)/(8*pi*pi*c)                         
*            = 6.62618*10**(-27)[erg-sec]*6.02205*10**23*10**16/                
*              (8*(3.1415926535)**2*2.997925*10**10[cm/sec])                    
*                                                                               
************************************************************************        
      const2=16.8576522d0                                                       
c                                                                               
      sumw=1.d-20                                                               
      sumwx=0.d0                                                                
      sumwy=0.d0                                                                
      sumwz=0.d0                                                                

      coord(1:3,1:numat)=xyz(1:3,1:numat)*bohr
                                                                                
         do 10 i=1,numat                                                        
            atmass=ams(nat(i))                                                  
            sumw=sumw+atmass                                                    
            sumwx=sumwx+atmass*coord(1,i)                                       
            sumwy=sumwy+atmass*coord(2,i)
            sumwz=sumwz+atmass*coord(3,i)
   10    continue                                                               
c                                                                               
      sumwx=sumwx/sumw                                                          
      sumwy=sumwy/sumw                                                          
      sumwz=sumwz/sumw                                                          
      f=1.0d0/bohr        
      do 30 i=1,numat                                                           
         x(i)=coord(1,i)-sumwx                                                  
         y(i)=coord(2,i)-sumwy                                                  
   30 z(i)=coord(3,i)-sumwz                                                     

************************************************************************        
*    matrix for moments of inertia is of form                                   
*                                                                               
*           |   y**2+z**2                         |                             
*           |    -y*x       z**2+x**2             | -i =0                       
*           |    -z*x        -z*y       x**2+y**2 |                             
*                                                                               
************************************************************************        
      do 40 i=1,6                                                               
   40 t(i)=dble(i)*1.0d-10                                                      
         do 50 i=1,numat                                                        
            atmass=ams(nat(i))                                                  
            t(1)=t(1)+atmass*(y(i)**2+z(i)**2)                                  
            t(2)=t(2)-atmass*x(i)*y(i)                                          
            t(3)=t(3)+atmass*(z(i)**2+x(i)**2)                                  
            t(4)=t(4)-atmass*z(i)*x(i)                                          
            t(5)=t(5)-atmass*y(i)*z(i)                                          
            t(6)=t(6)+atmass*(x(i)**2+y(i)**2)                                  
   50    continue                                                               
      call rsp(t,3,3,eig,evec)                                                  
         do 70 i=1,3                                                            
            if(eig(i).lt.3.d-4) then                                            
               rot(i)=0.d0                                                      
            else                                                                
               rot(i)=2.9979245d+4*const2/eig(i)  
            endif                                                               
   70    continue                                                               

      aa=rot(3)
      bb=rot(2)
      cc=rot(1)
                                                                                
      return                                                                    
      end                                                                       

