! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

module xtb_dynamic
   use xtb_io_writer, only : writeMolecule
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_filetypes, only : fileType
   use xtb_single, only : singlepoint
   use xtb_intmodes, only : xyzgeo
   use xtb_metadynamic

contains
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

pure subroutine xyzsort(n,lab,ew,xyz,list)
   implicit none
   integer, intent(in) :: n
   integer, intent(in) :: lab
   integer, intent(inout) :: list(*)
   real(wp),intent(inout) :: ew(*)
   real(wp),intent(inout) :: xyz(3,n,*)
   integer  :: i,ii,j,k,m
   integer  :: ihilf
   real(wp) :: pp,hilf

   do ii = 2, lab
      i = ii - 1
      k = i
      pp= ew(i)
      do j = ii, lab
         if (ew(j) .gt. pp) cycle
         k = j
         pp= ew(j)
      enddo
      if (k .eq. i) cycle
      ew(k) = ew(i)
      ew(i) = pp

      ihilf=list(i)
      list(i)=list(k)
      list(k)=ihilf
      do m=1,n
         do j=1,3
            hilf=xyz(j,m,i)
            xyz(j,m,i)=xyz(j,m,k)
            xyz(j,m,k)=hilf
         enddo
      enddo
   enddo

end subroutine xyzsort

pure subroutine xyzsort2(ndim,n,lab,ew,xyz)
   implicit none
   integer, intent(in) :: ndim
   integer, intent(in) :: n
   integer, intent(in) :: lab
   real(wp),intent(inout) :: ew(0:ndim)
   real(wp),intent(inout) :: xyz(3,n,*)
   integer  :: i,ii,j,k,m
   real(wp) :: pp, hilf

   do ii = 2, lab
      i = ii - 1
      k = i
      pp= ew(i)
      do j = ii, lab
         if (ew(j) .gt. pp) cycle
         k = j
         pp= ew(j)
      enddo
      if (k .eq. i) cycle
      ew(k) = ew(i)
      ew(i) = pp

      do m=1,n
         do j=1,3
            hilf=xyz(j,m,i)
            xyz(j,m,i)=xyz(j,m,k)
            xyz(j,m,k)=hilf
         enddo
      enddo
   enddo

end subroutine xyzsort2

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine wrc(fname,n,xyz,at)
   use xtb_setparam
   use xtb_mctc_symbols, only : toSymbol
   implicit none
   character*(*) fname
   integer n,at(n),j
   real(wp) xyz(3,n)
   integer :: ich ! file handle

   call open_file(ich,fname,'w')
   write(ich,'(''$coord'')')
   do j=1,n
      write(ich,'(3F24.10,5x,a2)') xyz(1:3,j),toSymbol(at(j))
   enddo
   write(ich,'(''$end'')')
   write(ich,'(''$set'')')
   write(ich,'('' chrg     '',i2)') set%ichrg
   write(ich,'('' uhf      '',i2)') set%nalphabeta
   write(ich,'(''$end'')')
   call close_file(ich)

end subroutine wrc

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine boltz(n,t,e,p)
   implicit none
   integer,  intent(in)  :: n
   real(wp), intent(in)  :: t
   real(wp), intent(in)  :: e(n)
   real(wp), intent(out) :: p(n)
   real(wp) f
   f = 8.31451_wp * t / 4.184e+3_wp
   p = exp(-e/f)
   p = p / sum(p)
end subroutine boltz

subroutine md(env,mol,chk,calc, &
      &       egap,et,maxiter,epot,grd,sigma,icall,Tsoll,cdump2)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autokcal, aatoau, amutokg, amutoau, fstoau
   use xtb_mctc_constants, only : pi, kB
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_calculator
   use xtb_type_restart
   use xtb_type_data
   use xtb_shake, only: do_shake,ncons
   use xtb_setparam
   use xtb_fixparam
   use xtb_scanparam
   use xtb_splitparam
   use xtb_type_setvar, only: metadyn_setvar
   implicit none

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   type(TMolecule),intent(inout) :: mol
   type(TRestart),intent(inout) :: chk
   class(TCalculator), intent(inout) :: calc
   integer  :: icall
   integer, intent(in) :: maxiter
   integer, intent(inout) :: cdump2
   real(wp), intent(inout) :: epot
   real(wp), intent(in) :: et
   real(wp), intent(inout) :: egap
   real(wp), intent(inout) :: grd(3,mol%n)
   real(wp), intent(inout) :: sigma(3,3)
   real(wp), intent(in) :: Tsoll
   type(scc_results) :: res

   real(wp) :: step,eel1,tstep,taut,xlam2,accu,driftthr
   real(wp) :: Ekin,tmass,f,mintime
   real(wp) :: Tinit,Tav,T,epav,ekav,dav,cma(3),bave,bavt
   real(wp) :: dum,edum,eerror,xx(10),molmass,slope,maxtime
   real(wp) :: tstep0,tmax,nfreedom,t0,w0,t1,w1,ep_prec,rege(4)
   real(wp) :: tors(mol%n),be(3),b0(3),tor,dip(3)
   real(wp) :: rcoord(3)
   logical  :: ex,thermostat,restart,confdump,equi,gmd,ldum

   integer :: i,j,k,ic,jc,ia,ja,ii,jj,ndum,cdump,nmax,ibin
   integer :: nstep,ndump,mdump,dumpstep,screendump,acount
   integer :: cdump0,nbo(0:20,mol%n),k3,nreg,ngeoav
   integer :: blockl,iblock,nblock
   real(wp),allocatable :: velo(:,:)
   real(wp),allocatable :: vel (:,:)
   real(wp),allocatable :: veln(:,:)
   real(wp),allocatable :: xyzo(:,:)
   real(wp),allocatable :: acc (:,:)
   real(wp),allocatable :: tmpx(:,:)
   real(wp),allocatable :: mass(:)
   real(wp),allocatable :: blocke(:),blockt(:)
   real(wp),allocatable :: intgeo(:,:),intgeo0(:,:)
   real(wp),allocatable :: intgeoav(:,:),tmpg(:,:)
   integer, allocatable :: na(:),nb(:),nc(:)
   character(len=80) :: atmp
   character(len=:),allocatable :: fname
   integer :: ich,trj,pdb,imdl
   logical :: exist

   ! Displace reference geometry by 1e-6
   real(wp), parameter :: atom_displacement = 1.0e-6_wp

   type(metadyn_setvar) :: metasetlocal
   real(wp) :: emtd
   real(wp) :: metatime
   metatime = 0.0_wp

   call delete_file('xtbmdok')

   if(icall.eq.0)then
      write(*,*)'trajectories on xtb.trj or xtb.trj.<n>'
   endif

   ! special equi and GMD settings
   tmax=set%time_md
   equi=.false.
   gmd=.false.
   if(icall.eq.-1) equi = .true.   ! equilibration
   if(icall.eq.-2)  gmd = .true.   ! GMD production run

   if(equi) then
      mintime = 5.0d0              ! minimum time in ps
      maxtime = tmax               ! maximum time in ps
      driftthr=2.d-3               ! stop MD if drift = slope Epot (block average) regression (4 points) below this value
   endif
   if(gmd) then
      mintime = tmax               ! minimum time in ps
      maxtime = tmax * 4           ! maximum time in ps
      driftthr=1.d-3
   endif

   ! real coord dump to e.g. scoord.n in siman
   confdump=.false.
   if(nscan.eq.0.and.cdump2.ge.0) confdump=.true.
   if(set%ceasefiles) confdump=.false.
   ! just screen
   screendump=200

   ! blocklength for SD of energy and T
   blockl=min(5000,idint(5000/set%tstep_md))
   allocate(blocke(blockl),blockt(blockl))

   ! take paramters from common
   tstep0=set%tstep_md
   cdump0=set%dump_md/set%tstep_md    ! scoord
   dumpstep=set%dump_md2/set%tstep_md ! xyz
   thermostat=set%nvt_md
   restart=set%restart_md
   accu=set%accu_md

   allocate(velo(3,mol%n),vel(3,mol%n),veln(3,mol%n),xyzo(3,mol%n),acc(3,mol%n),mass(mol%n))

   call neighbor(mol%n,mol%xyz,mol%at,nbo) ! neighbor list
   molmass=0
   do i=1,mol%n
      molmass=molmass+atmass(i)
      mass(i)=atmass(i)*amutoau
   enddo
   tmass=molmass*amutoau
   do i=1,mol%n
      if(mol%at(i).eq.1.and.set%md_hmass.gt.0) then
         mass(i)=dble(set%md_hmass)*amutoau
         ! k=nbo(1,i) ! atom to which H is bonded
         ! mass(k)=mass(k)-mass(i)+ams(1)*amutoau ! reduce by the increase H mass
      endif
   enddo
   molmass=molmass*amutokg

   nmax=tmax*1000/tstep0
   if(gmd.or.equi) nmax=maxtime*1000/tstep0

   write(*,*)
   write(*,'('' MD time /ps        :'',f8.2)')tmax
   write(*,'('' dt /fs             :'',f8.2)')tstep0
   write(*,'('' SCC accuracy       :'',f8.2)')accu
   write(*,'('' temperature /K     :'',f8.2)')Tsoll
   write(*,'('' max steps          :'',i6  )')nmax
   write(*,'('' block length (av.) :'',i6  )')blockl
   write(*,'('' dumpstep(trj) /fs  :'',f8.2,i6)')set%dump_md2,dumpstep
   write(*,'('' dumpstep(coords)/fs:'',f8.2,i6)')set%dump_md,cdump0
   if(equi.or.gmd)then
      write(*,'('' minimum runtime(ps):'',f6.1)')mintime
      write(*,'('' maximum runtime(ps):'',f6.1)')maxtime
   endif
   if(set%md_hmass.gt.0) &
      &write(*,'('' H atoms mass (amu) :'',i6  )')set%md_hmass

   nfreedom = 3._wp*real(mol%n,wp)
   if(set%mdrtrconstr ) nfreedom = nfreedom - 6.0d0
   if(set%shake_md    ) nfreedom = nfreedom - dble(ncons)
   if(zconstr.eq.1) nfreedom = nfreedom - dble(iatf1)  ! fragment 1 in Z plane
   write(*,'('' # deg. of freedom  :'',i6  )')idint(nfreedom)

   taut=500. ! damping of thermostat, 1000 is slow heating
   if(equi) taut=100.

   Tinit=Tsoll
   if (equi) Tinit=0.3*Tsoll   ! slow equi

   if(.not.restart)then
      f=1
      if(.not.thermostat) f = 2
      edum=f*Tinit*0.5*kB*nfreedom
      call mdinitu(mol%n,mol%at,velo,mass,edum)
   else
      call rdmdrestart(mol%n,mol%xyz,velo)
   endif

   if(set%shake_md) then
      write(*,*) 'SHAKE on. # bonds  :',ncons,' all:',.not.set%xhonly
   else
      write(*,*) 'SHAKE off'
   endif
   if(thermostat)   write(*,*) 'Berendsen THERMOSTAT on'
   if(equi)       write(*,*) 'EQUILIBRATION mode'
   if( gmd)       write(*,*) 'GMD mode'
   if(restart)    write(*,*) 'RESTART'


   !--- For "true" metadynamics an independent RMSD potential
   !    is constructed for this MD. "metaset" contains alpha and k
   !    "metasetlocal" is the dynamic potential for the MTD.
   !    If an ensemble was read into "metaset", only a static RMSD bias is applied.

   metasetlocal = metaset
   if((metaset%nstruc > 0).and.(metaset%static))then !if >0, an ensemble was read --> static potential
       metasetlocal%maxsave=0 !save nothing new --> deactivate dynamic potential
   else
       metaset%nstruc=0  !avoid calculation of RMSD potential within singlepoint routine
                         !--> would lead to double counting with "metasetlocal"
   endif

   if ((metasetlocal%maxsave.gt.0) .or. (metaset%nstruc > 0)) then
      write(env%unit,'(" --- metadynamics parameter ---")')
      write(env%unit,'(" kpush  :",f9.3)') metasetlocal%global_factor
      write(env%unit,'(" alpha  :",f9.3)') metasetlocal%global_width
      write(env%unit,'(" update :",i8)')   metasetlocal%maxsave
      if((metasetlocal%ramp - 0.03_wp) .ne. 0.0_wp)then
      write(env%unit,'(" ramp   :",f9.3)') metasetlocal%ramp
      endif
      if (metasetlocal%nstruc.eq.0) then
         do i = 1, mol%n
            ! Generate randomly displaced geometry
            do
               ! Generate numbers in [0,1]
               call random_number(rcoord)
               ! Convert numbers to [-1, 1]
               rcoord = 2.0_wp*rcoord-1.0_wp
               ! Ensure that displacement is large enough so that it
               ! can be normalized to the unit sphere
               if(norm2(rcoord) >= 1e-8) then
                  exit
               endif
            enddo
            ! Normalize displacement to the unit sphere
            rcoord = rcoord/norm2(rcoord)
            ! Assign displaced geometry
            metasetlocal%xyz(:,i,1) = mol%xyz(:,i) + atom_displacement*rcoord
         enddo
         metasetlocal%nstruc = 1
      else
      write(env%unit,'(" number of input RMSDs :",i4)') metasetlocal%nstruc
      endif
   endif

   atmp='xtb.trj'
   if(icall.gt.0)then
      write(atmp,'(''xtb.trj.'',i0)')icall
   endif
   call open_file(trj,trim(atmp),'w')
   pdb = -1
   if (allocated(mol%pdb) .and. icall == 0) then
      call open_file(pdb, 'xtb-trj.pdb', 'w')
      imdl = 0
   end if

   if(.not.thermostat)then
      write(*,'(9x,''time (ps)'',4x,''<Epot>'',6x,''Ekin   <T>   T'',5x, &
         &         ''Etot'',7x,''error'','' '')')
   else
      write(*,'(9x,''time (ps)'',4x,''<Epot>'',6x,''Ekin   <T>   T'',5x, &
         &              ''Etot'')')
   endif

   tstep=tstep0*fstoau

   mdump=screendump
   ndump=dumpstep

   call ekinet(mol%n,velo,mass,Ekin)

   grd=0.0_wp
   epot=0.0_wp
   call singlepoint &
      &     (env,mol,chk,calc, &
      &      egap,et,maxiter,0,.true.,.false.,1.0_wp,epot,grd,sigma,res)

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! MD loop
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   Edum   =Epot+Ekin
   Eerror =0
   Tav    =0
   Epav   =0
   Ekav   =0
   iblock =0
   nblock =0
   acount =0
   k3     =0
   ep_prec=0
   nreg   =0
   T      =0
   cdump  =0
   nstep  =0
   ngeoav =0

   call zeroz(mol%n,velo) ! optional z-plane fix
   !     call zero3n6(mol%n,velo,acc) ! fix 6 deg. of freedom
   call timing(t0,w0)
   do k=1,nmax

      !ccccccccccccccccccc
      ! calc E and forces                                       STEP 1
      !ccccccccccccccccccc
      accu=set%accu_md
      if(acount.eq.10)then  ! accurate SCC
         accu=1.0d0
         acount=0
      else
         acount=acount+1
      endif
      epot=0.0_wp
      grd = 0.0_wp
      call singlepoint &
         &     (env,mol,chk,calc, &
         &      egap,et,maxiter,0,.true.,.true.,accu,epot,grd,sigma,res)

      if (metasetlocal%maxsave.ne.0) then
         metatime = metatime + 1.0_wp
         metasetlocal%factor(1:metasetlocal%nstruc) = metasetlocal%global_factor
         metasetlocal%factor(metasetlocal%nstruc) = metasetlocal%factor(metasetlocal%nstruc) &
            &    * (2.0_wp/(1.0_wp+exp(-metasetlocal%ramp*metatime))-1.0_wp)
         if (cdump.gt.cdump0) then
            if (metasetlocal%nstruc.lt.metasetlocal%maxsave) then
               metatime = 0.0
               metasetlocal%nstruc = metasetlocal%nstruc + 1
               metasetlocal%xyz(:,:,metasetlocal%nstruc) = mol%xyz
            else
               do i = 2, metasetlocal%maxsave
                  metasetlocal%xyz(:,:,i-1) = metasetlocal%xyz(:,:,i)
               enddo
               metasetlocal%xyz(:,:,metasetlocal%maxsave) = mol%xyz
            endif
            write(*,'(2x,"adding snapshot to metadynamics bias")')
         endif
      !-------------------------------------------------------------------------
         emtd = 0.0d0
         call metadynamic (metasetlocal,mol%n,mol%at,mol%xyz,emtd,grd)
         epot = epot + emtd
      endif

      if(acount.eq.0)then  ! take only accurate Epot for average
         k3=k3+1
         ep_prec=ep_prec+epot
      endif

      ! okay, but you can also check for res%converged, can't you?
      if(epot.ne.epot.or. &
         &   epot.gt.1.d+5.or. &
         &   epot.lt.-1.d+5.or. & ! WHY?! this is completely resonable for >10000 atoms
         &   T.gt.10000) then
         write(*,*) epot,T
         write(*,*) 'MD is unstable, emergency exit '
         write(*,*) 'but still taking it as converged!'
         goto 1000
      endif

      if(k.eq.100)then
         call timing(t1,w1)
         write(*,'(''est. speed in wall clock h for 100 ps :'',f6.2)')  &
            &   100.0*(1000./tstep0)*(w1-w0)/100.0d0/3600.
      endif

      if(iblock.eq.blockl)then
         if(equi)then
            Tinit=Tinit*1.5
            Tinit=min(Tsoll,Tinit)
         endif
         nblock=nblock+1
         iblock=0
         call blocksd(mol%n,blockl,blocke,blockt,bave,bavt)
         nreg = nreg + 1
         rege(:) = [rege(2:), bave]
         if(nreg.ge.4)then
            call regress(1,4,rege,slope)
         else
            slope=99.
         endif
         write(*,'(''block <Epot> / <T> :'',f14.5,f6.0,5x, &
            &             ''drift:'',d10.2,3x,''Tbath :'',f5.0)')  &
            &             bave,bavt,slope,Tinit
      else
         iblock=iblock+1
         blocke(iblock)=epot
         blockt(iblock)=T
      endif

      !     dump coords to xyz and scoord
      !! ========================================================================
      ! dump for conformational sampling (scoor)
      if(cdump.gt.cdump0.and.confdump)then
         cdump2=cdump2+1
         call getname1(cdump2,atmp)
         call open_file(ich,trim(atmp),'w')
         call writeMolecule(mol, ich, fileType%tmol)
         call close_file(ich)
         cdump=0
      endif
      ! screen
      if(mdump.gt.screendump-1)then
         if(.not.thermostat)then
            write(*,'(i7,f8.2,F13.5,F9.4,2F6.0,F12.5,4F10.4)') &
               &   nstep,0.001*nstep*tstep/fstoau,Epav/k,Ekin,Tav/k,T,Epot+Ekin, &
               &   Edum/k-Epot-Ekin
         else
            write(*,'(i7,f8.2,F13.5,F9.4,2F6.0,F12.5,E14.6)') &
               &   nstep,0.001*nstep*tstep/fstoau,Epav/k,Ekin,Tav/k,T,Epot+Ekin
         endif
         mdump=0
      endif
      ! dump xyz (trj)
      if(ndump.gt.dumpstep-1)then
         ndump=0
         call writeMolecule(mol, trj, fileType%xyz, energy=epot, gnorm=res%gnorm)
         if(set%velodump)then
            do i=1,mol%n
               write(trj,'(3f20.14)')velo(1:3,i)
            enddo
         endif
         if (pdb /= -1) then
            imdl = imdl+1
            call writeMolecule(mol, pdb, fileType%pdb, number=imdl, &
               & energy=epot, gnorm=res%gnorm)
         end if
         if (set%forcewrrestart) then
            call wrmdrestart(mol%n,mol%xyz,velo)
         endif
      endif
      !! ========================================================================

      !     compute the accelaration at t
      do i=1,mol%n
         acc(:,i)=-grd(:,i)/mass(i)
      enddo

      call zeroz(mol%n,acc) ! z-plane fix
      ! call zero3n6(mol%n,velo,acc) ! fix 6 deg. of freedom

      ! store positions (at t); velocities are at t-1/2dt
      xyzo = mol%xyz

      !ccccccccccccccccccc
      ! temperature and pressure/density control               STEP 2
      !ccccccccccccccccccc

      ! estimate velocities at t
      veln = velo + 0.5d0*tstep*acc

      ! compute kinetic energy
      call ekinet(mol%n, veln, mass, Ekin)

      ! compute temperature
      T = 2.d0*Ekin/nfreedom/kB

      ! compute temperature scaling factors
      xlam2 = dsqrt(1.0d0 + (tstep0/taut)*(Tinit/T-1.0d0))

      !ccccccccccccccccccc
      ! velocity and position update                           STEP 3
      !ccccccccccccccccccc

      if(thermostat) then
         vel = xlam2 * (velo + acc*tstep)
      else
         vel = velo + acc*tstep
      endif

      ! update positions to t+dt
      mol%xyz = xyzo + vel*tstep

      !cccccccccccccccccc
      ! new temperature and pressure scaling factors
      !cccccccccccccccccc

      ! estimate velocities at t
      veln = 0.5d0*(velo + vel)

      ! compute kinetic energy
      call ekinet(mol%n, veln, mass, Ekin)

      ! compute temperature
      T = 2.d0*Ekin/nfreedom/kB

      ! compute temperature scaling factors
      if(thermostat) then
         xlam2 = dsqrt(1.0d0 + (tstep0/taut)*(Tinit/T-1.0d0))
      endif

      !ccccccccccccccccccc
      ! SHAKE (apply constraint at t+dt)
      !ccccccccccccccccccc

      if(set%shake_md) call do_shake(mol%n,xyzo,mol%xyz,vel,acc,mass,tstep)

      ! update velocities
      velo = vel

      ! remove trans/rot velocities
      call rmrottr(mol%n,mass,velo,mol%xyz)

      ! average internal coords
      if(gmd)then
         ngeoav=ngeoav+1
         tmpg=mol%xyz/aatoau
         call xyzgeo(tmpg,mol%n,na,nb,nc,1.d0,intgeo)
         do j=4,mol%n
            !          if(j.eq.20.or.j.eq.8)  write(*,*) j,tors(j),intgeo(3,j)
            if(tors(j)-intgeo(3,j).gt.pi)intgeo(3,j)=intgeo(3,j)+2.d0*pi ! check 0/180 torsion changes
            if(abs(tors(j)-intgeo(3,j)).gt.pi)then   ! check
               intgeo(3,j)=-2.0*pi-intgeo(3,j)
            endif
            if(abs(tors(j)-intgeo(3,j)).gt.pi) then  ! check
               intgeo(3,j)=-intgeo(3,j)
            endif
            if(abs(tors(j)-intgeo(3,j)).gt.pi) then  ! check
               intgeo(3,j)=2.*intgeo(3,j)
            endif
            if(abs(tors(j)-intgeo(3,j)).gt.pi) then  ! still in error
               write(*,*) j,tors(j),intgeo(3,j)
               stop 'error in md/zmat'
            endif
         enddo
         tors(1:mol%n)=intgeo(3,1:mol%n)
         intgeoav = intgeoav + intgeo-intgeo0
      endif

      ! compute averages
      nstep=nstep+1
      ndump=ndump+1
      mdump=mdump+1
      cdump=cdump+1
      Edum=Edum+Epot+Ekin
      Eerror=Edum/k-Epot-Ekin
      Tav =Tav+T
      Epav=Epav+epot
      Ekav=Ekav+ekin
      ! pmfav(k)=rcma

      if((equi.or.gmd).and.0.001*nstep*tstep/fstoau.gt.mintime)then ! check for equilibration end exit
         if( nblock.gt.1 .and. abs(bavt-Tinit)/Tinit.lt.0.02.and. &
            &       abs(slope).lt.driftthr) then
            write(*,*) 'GOOD AVERAGES REACHED'
            goto 1000
         endif
      endif
      if((equi.or.gmd).and.0.001*nstep*tstep/fstoau.gt.maxtime)then ! check for equilibration end exit
         write(*,*) 'MAXIMUM RUN TIME EXCEEDED'
         goto 1000
      endif

      ! end MD loop
   enddo

   ! exit
1000 call close_file(trj)
   if (pdb /= -1) then
      call close_file(pdb)
   end if

   write(*,*) 'average properties '
   write(*,*) 'Epot               :',Epav/k
   write(*,*) 'Epot (accurate SCC):',ep_prec/k3
   write(*,*) 'Ekin               :',Ekav/k
   write(*,*) 'Etot               :',(Ekav+Epav)/k
   write(*,*) 'T                  :',Tav/k

   if (abs(Tav/k-Tinit).gt.0.02*Tinit .and. &
      &   k.gt.500 .and. thermostat .and. (.not.equi)) &
      &   write(*,*)'thermostating problem'

   call wrmdrestart(mol%n,mol%xyz,velo)

   call touch_file('xtbmdok')

   write(*,*) 'normal exit of md()'

end subroutine md

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      pure subroutine regress(n1,n2,rege,slope)
      implicit none
      real(wp), intent(in) :: rege(:)
      real(wp), intent(out) :: slope
      integer, intent(in) :: n1,n2
      integer :: n,i
      real(wp) :: sx,sy,sxx,sxy,x

      n=n2-n1+1
      sx=0
      sy=0
      sxy=0
      sxx=0
      x=0
      do i=n1,n2
         x=x+1.
         sx=sx+x
         sy=sy+rege(i)
         sxx=sxx+x**2
         sxy=sxy+x*rege(i)
      enddo

      slope=(real(n, wp)*sxy-sx*sy)/(real(n, wp)*sxx-sx*sx)

      end subroutine regress

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine blocksd(n,nbl,ebl,tbl,esd,tsd)
      implicit none
      integer n,nbl
      real(wp) ebl(nbl),tbl(nbl),dum,av,esd,tsd
      integer i

      dum=0
      do i=1,nbl
         dum=dum+ebl(i)
      enddo
      av=dum/dble(nbl)
      esd=av

!     dum=0
!     do i=1,nbl
!        dum=dum+(ebl(i)-av)**2
!     enddo
!     esd=sqrt(dum/dble(nbl-1))/dble(n)

      dum=0
      do i=1,nbl
         dum=dum+tbl(i)
      enddo
      av=dum/dble(nbl)
      tsd=av
!     dum=0
!     do i=1,nbl
!        dum=dum+(tbl(i)-av)**2
!     enddo
!     tsd=sqrt(dum/dble(nbl-1))

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine ekinet(n,velo,mass,E)
   implicit none
   real(wp) velo(3,n),mass(n),E
   integer n,i
   e=0
   do i=1,n
      e=e+mass(i)*(velo(1,i)**2+velo(2,i)**2+velo(3,i)**2)
   enddo
   e=e*0.5
end subroutine ekinet

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! rescale velos to get temperature tsoll
! H's get half of the velo

      subroutine veloscal(n,iat,velo,mass,tsoll,t)
      implicit none
      integer n
      real(wp) velo (3,n)
      real(wp) mass(n)
      integer iat(n),nh,i,j
      real(wp) velo2(3,n)
      real(wp) tsoll,t,edum,f

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

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine wrmdrestart(n,xyz,velo)
   implicit none
   integer n,i
   real(wp) xyz(3,n),velo(3,n)
   integer :: ich

   call open_file(ich,'mdrestart','w')
   write(ich,*) '-1.0'
   do i=1,n
      write(ich,'(6D22.14)')xyz(1:3,i),velo(1:3,i)
   enddo
   call close_file(ich)

end subroutine wrmdrestart

subroutine rdmdrestart(n,xyz,velo)
   use xtb_setparam, only : get_namespace
   implicit none
   integer n,i
   real(wp) xyz(3,n),velo(3,n),dum
   integer :: ich
   character(len=:),allocatable :: fname

   call open_file(ich,'mdrestart','r')
   read(ich,*) dum
   do i=1,n
      read (ich,'(6D22.14)')xyz(1:3,i),velo(1:3,i)
   enddo
   call close_file(ich)

end subroutine rdmdrestart

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! initialize velocities uniformly
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine mdinitu(n,iat,velo,mass,Ekin)
   use xtb_setparam
   implicit none
   integer n
   real(wp) velo(3,n)
   integer iat(n)
   real(wp) mass(n)
   real(wp) Ekin
   real x,ranf

   real(wp) eperat,v,f,t,edum,f2
   integer i,j
   integer, allocatable :: iseed(:)

   call initrand

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

end subroutine mdinitu

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! integer random number n=<irand<=1

integer function irand(n)
   integer n
   real x,ranf
   call random_number(x)
   irand=n*x
   if(irand.lt.1) irand=1
   if(irand.gt.n) irand=n
end function irand

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine zeroz(n,array)
   use xtb_splitparam
   use xtb_scanparam
   implicit none
   integer n
   real(wp) array(3,*)
   integer i

   if(zconstr.eq.0) return

   do i=1,n
      if(splitlist(i).eq.1) array(3,i) = 0
   enddo

   return
end subroutine zeroz

subroutine zero3n6(n,v,a)
   use xtb_setparam
   implicit none
   integer n
   real(wp) a(3,n),v(3,n)

   if(.not.set%mdrtrconstr) return

   a(1:3,n)=0.d0
   a(1:2,n-1)=0.d0
   a(1,n-2)=0.d0
   v(1:3,n)=0.d0
   v(1:2,n-1)=0.d0
   v(1,n-2)=0.d0

   return
end subroutine zero3n6

end module xtb_dynamic
