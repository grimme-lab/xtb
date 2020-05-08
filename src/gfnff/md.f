      subroutine md(n,ichrg,at,xyz,q,grd,pdb) 
      use gff_mdcom               
      use gff_shake_module
      use gff_pdbio
      implicit none
      integer n,at(n),ichrg
      real*8 xyz(3,n)
      real*8 grd(3,n),q(n)

      real*8  step,amu2au,eel1,kboltz,taut,xlam2,accu,driftthr
      real*8  Ekin,fstoau,amutoau,angtoau,tmass,f,autokcal,mintime
      real*8  Tinit,Tav,T,ekav,dav,cma(3),dum2,maxtime,epav 
      real*8  dum,edum,eerror,xx(10),molmass,amutokg,metadynk
      real*8  nfreedom,t0,w0,t1,w1,eav
      real*8  pi,be(3),b0(3),tor,dip(3),epot,tstep
      logical ex,ldum
      parameter (pi=3.14159265358979D0)
      parameter (autokcal=627.509541d0)
      parameter (fstoau =0.413413733365614D+02)
      parameter (angtoau=1.0d0/0.529177260d0)
      parameter (amutoau=1.66053886E-27/9.10938215E-31)
      parameter (amutokg=1.66053886E-27)
      parameter (kboltz =0.316681534524639E-05)

      integer i,j,k,ic,jc,ia,ja,ii,jj,ndum,nmax,iav
      integer nstep,ndump,mdump,dumpstep,screendump
      integer nref,metadynnew
      real*8 ,allocatable ::velo(:,:)
      real*8 ,allocatable ::vel (:,:)
      real*8 ,allocatable ::veln(:,:)
      real*8 ,allocatable ::xyzo(:,:)
      real*8 ,allocatable ::acc (:,:)
      real*8 ,allocatable ::tmpx(:,:)
      real*8 ,allocatable ::mass(:)
      character*2 asym         
      character*80 atmp
      type(pdb_atom) pdb        
! mass
      real*8   :: ams(107)
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

      call execute_command_line('rm xtbmdok 2>/dev/null')
      write(*,'(13x,''======================================='')')
      write(*,'(13x,''|                                     |'')')
      write(*,'(13x,''|        Molecular Dynamics           |'')')
      write(*,'(13x,''|                                     |'')')
      write(*,'(13x,''======================================='')')
      write(*,*)'trajectory on xtb.trj'


c just screen      
      screendump=200
      dumpstep=dumpmd/tstep0

      allocate(velo(3,n),vel(3,n),veln(3,n),
     .         xyzo(3,n),acc(3,n),mass(n))

      molmass=0
      do i=1,n
         molmass=molmass+ams(at(i))
         mass(i)=ams(at(i))*amutoau
      enddo
      tmass=molmass*amutoau
      do i=1,n
         if(at(i).eq.1.and.hmass.gt.0) then
            mass(i)=hmass*amutoau
         endif
      enddo
      molmass=molmass*amutokg                      

      nmax=tmax*1000/tstep0

      write(*,*) 
      write(*,'('' MD time /ps        :'',f8.2)')tmax  
      write(*,'('' dt /fs             :'',f8.2)')tstep0
      write(*,'('' temperature /K     :'',f8.2)')Tsoll       
      write(*,'('' max steps          :'',i6  )')nmax 
      write(*,'('' dumpstep(trj) /fs  :'',f8.2,i6)')dumpmd,dumpstep
      write(*,'('' H atoms mass (amu) :'',f8.2)')hmass

      nfreedom = 3.*dble(n)
c     if(mdrtrconstr ) nfreedom = nfreedom - 6.0d0
      write(*,'('' # deg. of freedom  :'',i6  )')idint(nfreedom)

      taut=500. ! damping of thermostat, 1000 is slow heating

      Tinit=Tsoll

      f=1 
      if(.not.thermostat) f = 1
      edum=f*Tinit*0.5*kboltz*nfreedom             
      call mdinitu(n,at,velo,mass,edum)

      if(thermostat)   write(*,*) 'Berendsen THERMOSTAT on'

c     if(metadyn_k(1).ne.0)then ! make first = start reference geometry
c        write(*,'('' META-DYNAMICS on, k,alpha,update:'',2f9.3,i4)')
c    .   metadyn_k(1),metadyn_a,metadyn_n
c        metadynk = metadyn_k(1) ! save it for time dep. scaling
c        nmd=metadyn_n       ! save how many refs are max taken
c        do i=1,n
c           do j=1,3
c              call random_number(dum2)
c              metadyn_xyz(j,i,1)=xyz(j,i)+1.d-6*dum2
c           enddo
c        enddo
c        metadyn_n=1 ! initialize set counter in common
c     endif
c     metadynnew=0

      atmp='xtb.trj'
      open(unit=1,file=atmp)
      if(pdb%mode) call pdb%export('xtb.trj.pdb',xyz,n,3)
      if(.not.thermostat)then
      write(*,'(9x,''time (ps)'',4x,''<Epot>'',6x,''Ekin   <T>   T'',5x,
     .         ''Etot'',7x,''error'','' '')')
      else
      write(*,'(9x,''time (ps)'',4x,''<Epot>'',6x,''Ekin   <T>   T'',5x,
     .              ''Etot'')')
      endif

      tstep=tstep0*fstoau
       
      mdump=screendump
      ndump=dumpstep

      call ekinet(n,velo,mass,Ekin)

      call gfnff_eg(.false.,n,ichrg,at,xyz,.true.,grd,epot)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c MD loop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Edum   =Epot+Ekin
      Eerror =0
      Tav    =0
      Ekav   =0
      Epav   =0
      Eav    =0
      T      =0
      nstep  =0  
      iav    =0
 
                                                      call timing(t0,w0)
      do k=1,nmax

cccccccccccccccccccc      
c calc E and forces                                       STEP 1
cccccccccccccccccccc
      call gfnff_eg(.false.,n,ichrg,at,xyz,.true.,grd,epot)

c new metadynamics stuff, i.e., take every mdump fs a new ref structure
c     metadynnew=metadynnew+1
c     if(metadyn_k(1).ne.0) then  ! damp every new potential in the first steps
c        metadyn_k(1:metadyn_n) = metadynk
c         metadyn_k(metadyn_n)  = metadyn_k(metadyn_n) 
c    .    * (2.0/(1.0+exp(-0.03*dble(metadynnew)))-1.) ! damp the Vbias, ie 90 % in the first 50-100 steps
c     endif
c     if(cdump.gt.cdump0.and.metadyn_k(1).ne.0) then  ! every mdump fs
c        if(metadyn_n.lt.nmd)then ! fill it
c           metadyn_n = metadyn_n +1
c           metadynnew= 0         ! new one
c           metadyn_xyz(1:3,1:n,metadyn_n)=xyz(1:3,1:n)
c        else                     ! kick out oldest
c           do i=2,nmd
c              metadyn_xyz(1:3,1:n,i-1)=metadyn_xyz(1:3,1:n,i)
c           enddo
c           metadyn_xyz(1:3,1:n,nmd)=xyz(1:3,1:n)
c        endif
c     endif

      if(epot.ne.epot.or.
     .   epot.gt.1.d+5.or.
     .   epot.lt.-1.d+5.or.
     .   T.gt.10000) then 
         write(*,*) epot,T                              
         write(*,*) 'MD is unstable, emergency exit '
         write(*,*) 'but still taking it as converged!'
         goto 1000 
      endif
                                                      
      if(k.eq.100)then
                                                      call timing(t1,w1)
         write(*,'(''est. speed in wall clock h for 100 ps :'',f6.2)') 
     .   100.0*(1000./tstep0)*(w1-w0)/100.0d0/3600.
      endif

c     dump coords to xyz and scoord
c screen
      if(mdump.gt.screendump-1)then
         if(.not.thermostat)then
         write(*,'(i7,f8.2,F13.5,F9.4,2F6.0,F12.5,4F10.4)')
     .   nstep,0.001*nstep*tstep/fstoau,Epav/k,Ekin,Tav/k,T,Epot+Ekin,
     .   Edum/k-Epot-Ekin
         else
         write(*,'(i7,f8.2,F13.5,F9.4,2F6.0,F12.5,E14.6)')
     .   nstep,0.001*nstep*tstep/fstoau,Epav/k,Ekin,Tav/k,T,Epot+Ekin
         endif
         mdump=0
      endif
c dump xyz (trj) 
      if(ndump.gt.dumpstep-1)then
         ndump=0
      if(pdb%mode) then
      call pdb%export('xtb.trj.pdb',xyz,n,2)
      else
         write(1,*)n
         write(1,'(F20.12,f20.5,''  ##'')')Epot,0.001*nstep*tstep/fstoau
         do i=1,n
            write(1,'(1x,a2,1x,3f18.9)')
     .      asym(at(i)),xyz(1:3,i)/angtoau
         enddo
      endif
      endif

c     compute the accelaration at t
      do i=1,n
         do j=1,3
            acc(j,i)=-grd(j,i)/mass(i)
         enddo
      enddo

c     store positions (at t); velocities are at t-1/2dt
      xyzo = xyz

cccccccccccccccccccc      
c temperature and pressure/density control               STEP 2
cccccccccccccccccccc      

c     estimate velocities at t
      veln = velo + 0.5d0*tstep*acc

c     compute kinetic energy
      call ekinet(n, veln, mass, Ekin)

c     compute temperature
      T = 2.d0*Ekin/nfreedom/kboltz

c     compute temperature scaling factors
      xlam2 = dsqrt(1.0d0 + (tstep0/taut)*(Tinit/T-1.0d0))
      
cccccccccccccccccccc      
c velocity and position update                           STEP 3
cccccccccccccccccccc 

       if(thermostat) then
         do i=1,n
          do j=1,3
            vel(j,i)=xlam2*(velo(j,i)+acc(j,i)*tstep)
          enddo
         enddo
       else
        do i=1,n
          do j=1,3
           vel(j,i)=velo(j,i)+acc(j,i)*tstep
          enddo
        enddo
       endif

c     update positions to t+dt
        do i=1,n
          do j=1,3
             xyz(j,i) = xyzo(j,i) + vel(j,i)*tstep
          enddo
        enddo

ccccccccccccccccccc
c new temperature and pressure scaling factors
ccccccccccccccccccc

c     estimate velocities at t
      veln = 0.5d0*(velo + vel)

c     compute kinetic energy
      call ekinet(n, veln, mass, Ekin)

c     compute temperature
      T = 2.d0*Ekin/nfreedom/kboltz

c     compute temperature scaling factors
      if(thermostat) then
        xlam2 = dsqrt(1.0d0 + (tstep0/taut)*(Tinit/T-1.0d0))
      endif

cccccccccccccccccccc
c SHAKE (apply constraint at t+dt)
cccccccccccccccccccc

c     call do_shake(n,xyzo,xyz,vel,acc,mass,tstep)

c     update velocities
      velo = vel

c remove trans/rot velocities 
      call rmrottr(n,mass,velo,xyz)

c compute averages              
      nstep=nstep+1
      ndump=ndump+1
      mdump=mdump+1
      Edum=Edum+Epot+Ekin
      Eerror=Edum/k-Epot-Ekin
      Tav =Tav+T
      Ekav=Ekav+ekin
      Epav=Epav+epot
      if(k.gt.nmax/3) then
         iav=iav+1
         Eav=Eav+epot
      endif

c end MD loop
      enddo

c     exit
1000  close(1)

      write(*,*) 'average properties '
      write(*,*) 'Epot               :',Epav/k
      write(*,*) 'Epot (last 2/3)    :',Eav/iav
      write(*,*) 'Ekin               :',Ekav/k
      write(*,*) 'Etot               :',(Ekav+Epav)/k
      write(*,*) 'T                  :',Tav/k

      call execute_command_line('touch xtbmdok')
      write(*,*) 'normal exit of md()'

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ekinet(n,velo,mass,E)
      implicit none
      real*8 velo(3,n),mass(n),E
      integer n,i
      e=0
      do i=1,n
         e=e+mass(i)*(velo(1,i)**2
     .               +velo(2,i)**2
     .               +velo(3,i)**2)
      enddo
      e=e*0.5
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c rescale velos to get temperature tsoll
c H's get half of the velo 

      subroutine veloscal(n,iat,velo,mass,tsoll,t)
      implicit none
      integer n
      real*8 velo (3,n)
      real*8 mass(n)
      integer iat(n),nh,i,j
      real*8 velo2(3,n)
      real*8 tsoll,t,edum,f

      velo = velo * sqrt(tsoll/t)

      velo2 = velo
      j=0
      f=0.5d0
 10   j=j+1  
      do i=1,n
         if(iat(i).eq.1)then
            velo(1:3,i)=velo2(1:3,i)*f
         else
            velo(1:3,i)=velo2(1:3,i)
         endif
      enddo

 20   call ekinet(n,velo,mass,edum)
      t=edum/(0.5*3*n*0.316681534524639E-05)
      f=f+0.0001
      if(abs(t-tsoll).gt.1.and.j.lt.100000) goto 10
      
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine wrmdrestart(n,xyz,velo)
      implicit none
      integer n,i
      real*8 xyz(3,n),velo(3,n)

      open(unit=112,file='mdrestart')
      write(112,*) '-1.0'
      do i=1,n
         write(112,'(6D22.14)')xyz(1:3,i),velo(1:3,i)
      enddo
      close(112)

      end

      subroutine rdmdrestart(n,xyz,velo)
      implicit none
      integer n,i
      real*8 xyz(3,n),velo(3,n),dum

      open(unit=112,file='mdrestart')
      read(112,*) dum
      do i=1,n
         read (112,'(6D22.14)')xyz(1:3,i),velo(1:3,i)
      enddo
      close(112)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c initialize velocities uniformly
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine mdinitu(n,iat,velo,mass,Ekin)
      implicit none
      integer n
      real*8 velo(3,n)
      integer iat(n)
      real*8 mass(n)
      real*8 Ekin   
      real x,ranf

      real*8 eperat,v,f,t,edum,f2
      integer i,j
      integer, allocatable :: iseed(:)

c     if(samerand)then
c       call random_seed(size=j)
c       allocate(iseed(j))
c       iseed(1:j)=41                    ! start random number generator for same sequenz
c       do i=1,j
c          iseed(i)=iseed(i)+j
c       enddo
c       call random_seed(put=iseed)
c       deallocate(iseed)
c     else
        call random_seed()
c     endif

      eperat=Ekin/(3.0*n)

      do i=1,n
         f2=1
         if(iat(i).eq.1) f2=2
         v=sqrt(2*eperat/mass(i))
         f=1.0d0
         call random_number(x)
         if(x.gt.0.5)f=-1.0d0
         velo(1,i)=v*f*f2
         f=1.0d0
         call random_number(x)
         if(x.gt.0.5)f=-1.0d0
         velo(2,i)=v*f*f2
         f=1.0d0
         call random_number(x)
         if(x.gt.0.5)f=-1.0d0
         velo(3,i)=v*f*f2
      enddo 

      call ekinet(n,velo,mass,edum)

      t=edum/(0.5*3*n*0.316681534524639E-05)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c integer random number n=<irand<=1

      integer function irand(n)
      integer n
      real x,ranf
      call random_number(x)
      irand=n*x
      if(irand.lt.1) irand=1
      if(irand.gt.n) irand=n
      end
