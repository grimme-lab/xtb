! This file is part of xtb.
!
! Copyright (C) 2019-2020 Stefan Grimme
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
module xtb_gfnff_ini

contains

subroutine gfnff_ini(env,pr,makeneighbor,mol,gen,param,topo,accuracy)
      use xtb_mctc_accuracy, only : wp, sp
      use xtb_type_molecule
      use xtb_type_environment, only : TEnvironment
      use xtb_gfnff_param, only : efield, gfnff_thresholds
      use xtb_gfnff_data, only : TGFFData
      use xtb_gfnff_topology, only : TGFFTopology
      use xtb_gfnff_generator, only : TGFFGenerator
      use xtb_gfnff_ini2
      use xtb_gfnff_eg, only : gfnff_dlogcoord
      use xtb_gfnff_mrec, only : mrecgff
      use xtb_gfnff_fraghess
      use xtb_restart
      use xtb_mctc_constants

      implicit none
      character(len=*), parameter :: source = 'gfnff_ini'
!--------------------------------------------------------------------------------------------------
      type(TEnvironment), intent(inout) :: env
      type(TMolecule), intent(in) :: mol   ! # molecule type
      type(TGFFTopology), intent(inout) :: topo
      type(TGFFGenerator), intent(in) :: gen
      type(TGFFData), intent(in) :: param
      real(wp), intent(in) :: accuracy

      logical, intent(in) :: pr            ! print flag
      logical, intent(in) :: makeneighbor  ! make a neigbor list or use existing one?
!--------------------------------------------------------------------------------------------------

      integer ati,atj,atk,i,j,k,l,lin,nn,ii,jj,kk,ll,m,rings,ia,ja,ij,ix,nnn,idum,ip,ji,no
      integer ineig,jneig,nrot,bbtyp,ringtyp,nn1,nn2,hybi,hybj,pis,ka,nh,jdum,hcalc,nc
      integer ringsi,ringsj,ringsk,ringl,npi,nelpi,picount,npiall,maxtors,rings4,nheav
      integer nm,maxhb,ki,n13,current,ncarbo,mtyp1,mtyp2
      integer ind3(3),sr(20),cr(10,20),niel(86)
      integer qloop_count,nf,nsi,nmet,nhi,nhj,ifrag
      integer hbA,hbH,Bat,atB,Aat,Hat
      integer AHB_nr
      integer bond_hbn
      interface
         integer function itabrow6(i)
            integer i
         end function
      end interface

      real(wp) r0,ff,omega,f1,f2,phi,valijklff,ringf,fcn
      real(wp) shift,dum,dum1,dum2,dum4,qafac,fqq,feta
      real(wp) sumppi,fpi,fxh,fijk,fsrb2,ees
      real(wp) fheavy,fn,eold,fctot,fij
      real(wp) hbpi(2),hbpj(2),sdum3(3)
      real(wp) bstrength
      real(wp) xx(20)
      real(wp) fkl,qreps,fbsmall,bohr

      parameter (bohr=1.0_wp/0.52917726_wp)

      real(wp), parameter :: rabd_cutoff = 13.0_wp

      logical lring,picon,notpicon,bridge,sp3ij,ccij,success
      logical heavy,triple,piat,sp3kl,ex,cnij,frag_charges_known

      integer,allocatable :: btyp(:),imetal(:),nbm(:,:),nbf(:,:)
      integer,allocatable :: itag(:)
      integer,allocatable :: piadr(:),piadr2(:),piadr3(:),piadr4(:)
      integer,allocatable :: itmp(:),sring(:,:),cring(:,:,:)
      integer,allocatable :: ipis(:),pimvec(:),nbpi(:,:),piel(:)
      integer,allocatable :: lin_AHB(:)
      integer,allocatable :: bond_hbl(:,:)
      real(wp),allocatable:: rab  (:)
      real(wp),allocatable:: sqrab(:)
      real(wp),allocatable:: cn   (:)
      real(wp),allocatable:: dcn(:,:,:)
      real(wp),allocatable:: dgam(:), dxi(:)
      real(wp),allocatable:: mchar(:)
      real(wp),allocatable:: rtmp (:)
      real(wp),allocatable:: pbo  (:)
      real(wp),allocatable:: qtmp (:), dqa(:), qah(:)
      real(wp),allocatable:: Api(:,:),S(:,:),Pold(:,:),pibo(:),occ(:),eps(:)
      real(wp),allocatable:: pispop(:),pisea(:),pisip(:),apisave(:,:)
      real(sp),allocatable:: rabd(:,:)

      character(len=255) atmp
      integer  :: ich, err
      real(wp) :: dispthr, cnthr, repthr, hbthr1, hbthr2
      logical :: exitRun

      call gfnff_thresholds(accuracy, dispthr, cnthr, repthr, hbthr1, hbthr2)

      if (pr) then
         write(env%unit,*)
         write(env%unit,'(10x,"entering GFN-FF setup routine... ",i0)') mol%n
      endif

      write(env%unit,*)
      write(env%unit,'(10x,"==================== Thresholds ====================")')
      write(env%unit,'(10x,"CN  :",f12.5)')   cnthr
      write(env%unit,'(10x,"rep :",f12.5)')   repthr
      write(env%unit,'(10x,"disp:",f12.5)')   dispthr
      write(env%unit,'(10x,"HB1 :",f12.5)')   hbthr1
      write(env%unit,'(10x,"HB2 :",f12.5)')   hbthr2
      write(env%unit,*)

      allocate( rab(mol%n*(mol%n+1)/2), source = 0.0d0 )
      allocate( cn(mol%n), source = 0.0d0 )
      allocate( sqrab(mol%n*(mol%n+1)/2), source = 0.0d0 )
      allocate( topo%hyb(mol%n), source = 0 )
      allocate( topo%alphanb(mol%n*(mol%n+1)/2), source = 0.0d0 )
      allocate( rtmp(mol%n*(mol%n+1)/2), source = 0.0d0 )
      allocate( pbo(mol%n*(mol%n+1)/2), source = 0.0d0 )
      allocate( piadr(mol%n), source = 0 )
      allocate( piadr2(mol%n), source = 0 )
      allocate( topo%bpair(mol%n*(mol%n+1)/2), source = 0 )
      allocate( itmp(mol%n), source = 0 )
      allocate( itag(mol%n), source = 0 )
      allocate( sring(20,mol%n), source = 0 )
      allocate( cring(10,20,mol%n), source = 0 )
      allocate( piadr3(mol%n), source = 0 )
      allocate( piadr4(mol%n), source = 0 )
      allocate( qtmp(mol%n), source = 0.0d0 )
      allocate( dxi(mol%n), source = 0.0d0 )
      allocate( dgam(mol%n), source = 0.0d0 )
      allocate( topo%chieeq(mol%n), source = 0.0d0 )
      allocate( topo%gameeq(mol%n), source = 0.0d0 )
      allocate( topo%alpeeq(mol%n), source = 0.0d0 )
      allocate( topo%qa(mol%n), source = 0.0d0 )
      allocate( dqa(mol%n), source = 0.0d0)
      allocate( qah(mol%n), source = 0.0d0 )
      allocate( nbm(20,mol%n), source = 0 )
      allocate( mchar(mol%n), source = 0.0d0 )
      allocate( imetal(mol%n), source = 0 )
      allocate( topo%zetac6(mol%n*(mol%n+1)/2), source = 0.0d0 )
      allocate( topo%xyze0(3,mol%n), source = 0.0d0 )
      allocate( nbf(20,mol%n), source = 0 )

      niel=0
      do i=1,mol%n
         niel(mol%at(i))=niel(mol%at(i))+1
      enddo

      write(env%unit,'(10x,"Pauling EN used:")')
      do i=1,86
         if(niel(i).gt.0) write(env%unit,'(10x,"Z :",i2,"  EN :",f6.2)') i,param%en(i)
      enddo

      dum = sqrt(sum(efield**2))
      write(env%unit,'(10x,"electric field strengths (au):",f6.3)') dum
!     alp = alp *(1.+0.0*dum)

      write(env%unit,*)
      write(env%unit,'(10x," ------------------------------------------------- ")')
      write(env%unit,'(10x,"|           Force Field Initialization            |")')
      write(env%unit,'(10x," ------------------------------------------------- ")')
      write(env%unit,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! distances and bonds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      topo%xyze0 = mol%xyz ! initial geom

      write(env%unit,'(10x,"distances ...")')
      pbo   = 0
      rab   = 0
      sqrab = 0
      do i=1,mol%n
         ati=mol%at(i)
         kk=i*(i-1)/2
         do j=1,i-1
            atj=mol%at(j)
            k=kk+j
            sqrab(k)=(mol%xyz(1,i)-mol%xyz(1,j))**2+(mol%xyz(2,i)-mol%xyz(2,j))**2+(mol%xyz(3,i)-mol%xyz(3,j))**2
            rab(k)  =sqrt(sqrab(k))
            if(rab(k).lt.1.d-3) then
               write(env%unit,*) i,j,ati,atj,rab(k)
               call env%error("Particular close distance present", source)
               exit
            endif
         enddo
      enddo

      call env%check(exitRun)
      if (exitRun) then
         return
      end if

!     Calculate CN and derivative
      allocate(dcn(3,mol%n,mol%n), source = 0.0d0 )
      call gfnff_dlogcoord(mol%n,mol%at,mol%xyz,rab,cn,dcn,cnthr,param) ! dcn needed
      do i=1,mol%n
         dum2=0
         do j=1,mol%n
            dum2=dum2+sqrt(dcn(1,j,i)**2+dcn(2,j,i)**2+dcn(3,j,i)**2)
         enddo
         mchar(i) = exp(-0.005d0*param%en(mol%at(i))**8)*dum2/(cn(i)+1.0d0)     ! estimated metallic character as ratio of av. dCN and CN
                                                                      ! and an EN cut-off function, used in neigbor routinen and for BS estimate
      enddo
      deallocate(dcn)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! neighbor list, hyb and ring info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      topo%qa = 0
      qloop_count = 0

!111   continue
!  do the loop only if factor is significant
   do while (qloop_count.lt.2.and.gen%rqshrink.gt.1.d-3)

      write(env%unit,'(10x,"----------------------------------------")')
      write(env%unit,'(10x,"generating topology and atomic info file ...")')
      call gfnff_neigh(env,makeneighbor,mol%n,mol%at,mol%xyz,rab,gen%rqshrink, &
         & gen%rthr,gen%rthr2,gen%linthr,mchar,topo%hyb,itag,nbm,nbf,param,topo)

      do i=1,mol%n
         imetal(i)=param%metal(mol%at(i))
         if(topo%nb(20,i).le.4.and.param%group(mol%at(i)).gt.3) imetal(i)=0 ! Sn,Pb,Bi, with small CN are better described as non-metals
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! bonds (non bonded directly in EG)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      topo%bpair=0
      do i=1,mol%n
         do j=1,topo%nb(20,i)
            k=topo%nb(j,i)
            topo%bpair(lin(k,i))=1
         enddo
      enddo
      topo%nbond = sum(topo%bpair)
      topo%nbond_blist = topo%nbond
      allocate( topo%blist(2,topo%nbond), source = 0 )
      allocate( btyp(topo%nbond), source = 0 )
      allocate( pibo(topo%nbond), source = 0.0d0 )

      pibo  = -99.
      topo%nbond = 0
      do i=1,mol%n
         kk=i*(i-1)/2
         do j=1,i-1
            k=kk+j
            if ( topo%bpair(k) .eq. 1 ) then  ! bonds
                topo%nbond = topo%nbond +1
                topo%blist(1,topo%nbond)=i
                topo%blist(2,topo%nbond)=j
            endif
         enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hueckel setup for all first-row sp2 and sp atoms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! setup list of all possible pi atoms
      k=0
      piadr =0
      piadr2=0
      do i=1,mol%n ! setup loop
         piat =(topo%hyb(i).eq.1.or.topo%hyb(i).eq.2).and.pilist(mol%at(i)) ! sp or sp2 and CNOFS
         kk=0
         do j=1,topo%nb(20,i)
            jj=topo%nb(j,i)
            if(mol%at(i).eq.8.and.mol%at(jj).eq.16.and.topo%hyb(jj).eq.5) then
                                                             piat=.false.
                                                             cycle ! SO3   is not a pi
                                                             endif
            if(topo%hyb(jj).eq.1.or.topo%hyb(jj).eq.2)  kk=kk+1         ! attached to sp2 or sp
         enddo
         picon=kk.gt.0.and.nofs(mol%at(i))                     ! an N,O,F (sp3) on sp2
         if(mol%at(i).eq. 7.and.topo%nb(20,i).gt.3) cycle           ! NR3-X is not a pi
         if(mol%at(i).eq.16.and.topo%hyb(i).eq.5  ) cycle           ! SO3   is not a pi
         if(picon.or.piat) then
            k=k+1
            piadr (k)=i
            piadr2(i)=k
         endif
      enddo
      npiall=k
! make pi neighbor list
      allocate( nbpi(20,npiall),pimvec(npiall), source = 0 )
      nbpi=0
      do i=1,mol%n
         if(piadr2(i).eq.0) cycle
         ii=piadr2(i)
         nbpi(20,ii)=0
         do j=1,topo%nb(20,i)
            k=topo%nb(j,i)
            if(piadr2(k).gt.0)then
               nbpi(20,ii)=nbpi(20,ii)+1
               nbpi(nbpi(20,ii),ii)=piadr2(k)
            endif
         enddo
      enddo

! assign pi atoms to fragments
      call mrecgff(npiall,nbpi,picount,pimvec)
      deallocate(nbpi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! setup xi correction for EEQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      dxi = 0 ! default none
      do i=1,mol%n
         ati=mol%at(i)
         nn =topo%nb(20,i)
         if(nn.eq.0) cycle
         ip =piadr2(i)
         ji =topo%nb(1,i) ! first neighbor
         nh =0
         nm =0
         do j=1,nn
            if(mol%at(topo%nb(j,i)).eq.1)        nh=nh+1
            if(imetal(  topo%nb(j,i) ).ne.0) nm=nm+1
         enddo
!     hydrogen
!        if(ati.eq.1.and.nn.gt.1)                                                 dxi(i)=dxi(i)-nn*0.01
!     boron
         if(ati.eq.5)                                                             dxi(i)=dxi(i)+nh*0.015
!     carbon
         if(ati.eq.6.and.nn.eq.2.and.itag(i).eq.1)                                dxi(i)=-0.15 ! make carbene more negative
!        if(ati.eq.6.and.nn.eq.2)then
!           ki=topo%nb(2,i)
!           if(mol%at(ki).eq.8.and.mol%at(ji).eq.8.and.topo%nb(20,ji).eq.1.and.topo%nb(20,ki).eq.1)then          ! free CO2
!                                                                                 dxi(ki)=0.19 ! lower EN for O
!                                                                                 dxi(ji)=0.19 !  "    "   "  "
!           endif
!        endif
         if(ati.eq.6.and.nn.eq.1.and.mol%at(ji).eq.8.and.topo%nb(20,ji).eq.1)              dxi(ji)=0.15! free CO
!     nitrogen
!        if(ati.eq.7.and.nn.eq.1.and.mol%at(ji).eq.6)                                 dxi(i)=0.00  !CN
!     oxygen / group 6
         if(ati.eq.8.and.nn.eq.1.and.ip.ne.0.and.mol%at(ji).eq.7.and.piadr2(ji).ne.0) dxi(i)= 0.05    ! nitro oxygen, otherwise NO2 HBs are too strong
         if(ati.eq.8.and.nn.eq.2.and.nh.eq.2)                                     dxi(i)=-0.02    ! H2O
         if(param%group(ati).eq.6.and.nn.gt.2)                                          dxi(i)=dxi(i)+nn*0.005! good effect
         if(ati.eq.8.or.ati.eq.16)                                                dxi(i)=dxi(i)-nh*0.005
!    fluorine / group 7
         if(param%group(ati).eq.7.and.ati.gt.9.and.nn.gt.1) then ! polyvalent Cl,Br ...
                                         if(nm.eq.0)then
                                                                                  dxi(i)=dxi(i)-nn*0.021! good effect
                                         else
                                                                                  dxi(i)=dxi(i)+nn*0.05 ! good effect for TMs
                                         endif
         endif
      enddo

!     prepare EEQ xi ATOMIC parameter
!     at this point for the non-geom. dep. charges qa with CN = nb
      do i=1,mol%n
         ati=mol%at(i)
         dum =min(dble(topo%nb(20,i)),gen%cnmax)  ! limits it
!                   base val  spec. corr.    CN dep.
         topo%chieeq(i)=-param%chi(ati) + dxi(i) + param%cnf(ati)*sqrt(dum)
         topo%gameeq(i)= param%gam(ati)
         if(imetal(i).eq.2)then           ! the "true" charges for the TM metals are small (for various reasons)
            topo%chieeq(i)=topo%chieeq(i)-gen%mchishift ! so take for the non-geom. dep. ones less electronegative metals yield more q+
         endif                            ! which reflect better the true polarity used for guessing various
                                          ! potential terms. The positive effect of this is big.
         topo%alpeeq(i)= param%alp(ati)**2
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! topology based charges
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(env%unit,'(10x,"pair mat ...")')
      call nbondmat(mol%n,topo%nb,topo%bpair)  ! get number of cov. bonds between atoms up to 4 bonds
      write(env%unit,'(10x,"computing topology distances matrix with Floyd-Warshall algo ...")')
      allocate( rabd(mol%n,mol%n), source = 0.0e0_sp)
      rabd = rabd_cutoff
!     determine topology distances by Floyd-Warshall algo
!     they are used in the EEQ to determine qa (approximate topology charges)
      do i = 1, mol%n
        rabd(i, i) = 0.0
        do j = 1, topo%nb(20, i)
           k=topo%nb(j,i)
           rabd(k, i) = param%rad(mol%at(i)) + param%rad(mol%at(k))
           rabd(i, k) = rabd(k, i)
        end do
      end do
      do k = 1, mol%n
      do i = 1, mol%n
         if (rabd(i, k) > gen%tdist_thr) cycle
         do j = 1, mol%n
            if (rabd(k, j) >  gen%tdist_thr) cycle
            if (rabd(i, j) > (rabd(i, k) + rabd(k, j))) then
                rabd(i, j) =  rabd(i, k) + rabd(k, j)
            end if
         end do
      end do
      end do

      do i=1,mol%n
         do j=1,i-1
            ij=lin(j,i)
            if(rabd(j,i).gt.gen%tdist_thr) rabd(j,i)=rabd_cutoff ! values not properly considered
            rtmp(ij) = gen%rfgoed1* rabd(j,i) / 0.52917726d0
         enddo
      enddo
      deallocate(rabd)

      frag_charges_known=.false.
      write(env%unit,'(10x,"making topology EEQ charges ...")')
      if(topo%nfrag.le.1) then                           ! nothing is known
!     first check for fragments
      call mrecgff(mol%n,nbf,topo%nfrag,topo%fraglist)
      write(env%unit,'(10x,"#fragments for EEQ constrain: ",i0)') topo%nfrag
!     read QM info if it exists
      call open_file(ich, 'charges', 'r')
      if (ich /= -1) then
         qtmp = 0
         err = 0
         i = 0
         do while(err == 0)
            read(ich,*,iostat=err) dum
            if (err /= 0) exit
            if (i < mol%n) then
               i = i+1
               qtmp(topo%fraglist(i))=qtmp(topo%fraglist(i))+dum
            else
               call env%warning("More charges than atoms present, assuming missmatch", source)
               err = 1
            end if
         enddo
         if (is_iostat_end(err) .and. i == mol%n) err = 0
         call close_file(ich)
         if (err == 0) then
            if (i < mol%n .or. abs(sum(qtmp) - mol%chrg) > 1.0e-3_wp) then
               call env%warning("Rejecting external charges input due to missmatch", source)
            else
               topo%qfrag=dnint(qtmp)
               write(env%unit,'(10x,"fragment charges from <charges> :",10F7.3)') topo%qfrag(1:topo%nfrag)
            end if
         else
            call env%warning("Could not initialize fragment charges from file", source)
         endif

         call env%check(exitRun)
         if (exitRun) then
            return
         end if
      endif
      if(mol%n.lt.100.and.topo%nfrag.gt.2.and.mol%chrg.ne.0.and.sum(topo%qfrag(2:topo%nfrag)).gt.999) then
         itmp=0
         do i=1,mol%n
            itmp(topo%fraglist(i))=itmp(topo%fraglist(i))+1
         enddo
         do i=1,topo%nfrag
            write(env%unit,*)i,itmp(i)
         enddo
         call env%error('fragment charge input required', source)
         return
      endif
      if(mol%n.ge.100.and.topo%nfrag.gt.2.and.mol%chrg.ne.0.and.sum(topo%qfrag(2:topo%nfrag)).gt.999) then
         topo%qfrag(1)=mol%chrg
         topo%qfrag(2:topo%nfrag)=0
      endif
      if(topo%nfrag.eq.2.and.mol%chrg.ne.0.and.sum(topo%qfrag(2:topo%nfrag)).gt.999) then
         write(env%unit,*) 'trying auto detection of charge on 2 fragments:'
         topo%qfrag(1)=0
         topo%qfrag(2)=mol%chrg
         call goedeckera(env,mol%n,mol%at,topo%nb,rtmp,topo%qa,dum1,topo)
         call env%check(exitRun)
         if (exitRun) then
            call env%error("Failed to generate charges", source)
            return
         end if
         topo%qfrag(2)=0
         topo%qfrag(1)=mol%chrg
         call goedeckera(env,mol%n,mol%at,topo%nb,rtmp,topo%qa,dum2,topo)
         call env%check(exitRun)
         if (exitRun) then
            call env%error("Failed to generate charges", source)
            return
         end if
         if(dum1.lt.dum2) then
            topo%qfrag(1)=0
            topo%qfrag(2)=mol%chrg
         endif
         write(env%unit,*) 'dEes      :',dum1-dum2
         write(env%unit,*) 'charge 1/2:',topo%qfrag(1:2)
      endif
      else if (allocated(mol%pdb)) then ! frag_charges_known
         write(env%unit,'(10x,"#fragments for EEQ constrain from pdb file: ",i0)') topo%nfrag
         frag_charges_known=.true.
      endif

!     make estimated, topology only EEQ charges from rabd values, including "right" fragment charge
      call goedeckera(env,mol%n,mol%at,topo%nb,rtmp,topo%qa,ees,topo)
      call env%check(exitRun)
      if (exitRun) then
         call env%error("Failed to generate charges", source)
         return
      end if

!     estimate how much of the frag charge is on the pi sub systems
      if(picount.gt.0.and.qloop_count.gt.0) then
         allocate( ipis(picount), source = 0 )
         if(frag_charges_known) then                        ! PDB case
         ipis = 0
         do pis=1,picount ! loop over pi systems
            do k=1,npiall
               if(pimvec(k).eq.pis) then
                  ipis(pis)=ipis(pis)+topo%qpdb(piadr(k))
                  topo%qpdb(piadr(k))=0
               endif
            enddo
         enddo
         else                                               ! general case
         qtmp = topo%qa ! save the "right" ones
         qah  = topo%qa
         call qheavy(mol%n,mol%at,topo%nb,qah) ! heavy atoms only ie H condensed to neighbor
         do pis=1,picount ! loop over pi systems
            do k=1,npiall
               if(pimvec(k).eq.pis) then
                  kk=piadr(k)
                  ifrag=topo%fraglist(kk) !the pi atom of this pi fragment is in EEQ fragment ifrag
                  exit
               endif
            enddo
            dum2=topo%qfrag(ifrag) ! save
            topo%qfrag(ifrag) = 0 ! make only this EEQ fragment neutral
            call goedeckera(env,mol%n,mol%at,topo%nb,rtmp,topo%qa,ees,topo) ! for neutral
            call env%check(exitRun)
            if (exitRun) then
               call env%error("Failed to generate charges", source)
               return
            end if
            topo%qfrag(ifrag) = dum2 ! back
            call qheavy(mol%n,mol%at,topo%nb,topo%qa)
            dqa =qah-topo%qa ! difference charges upon ionization
            dum1=0
            dum=0
            do k=1,npiall
               if(pimvec(k).eq.pis) dum=dum+dqa(piadr(k)) ! only pi atoms
            enddo
            dum = dum * 1.1 !charges tend to be slightly too small 1.1-1.2
            ipis(pis)=idnint(dum)
            dum1=dum1+dum
         enddo
         topo%qa = qtmp ! put "right" charges used in FF construction and for HB/XB in place
         endif
      endif

      if(qloop_count.eq.0) itmp(1:mol%n)=topo%nb(20,1:mol%n)
      qloop_count=qloop_count+1
      if(qloop_count.lt.2.and.gen%rqshrink.gt.1.d-3) then  ! do the loop only if factor is significant

         deallocate(topo%blist,btyp,pibo,pimvec)
!         goto 111
      endif
   end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! change EEQ J with estimated q
! which is a kind of third-order term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=1,mol%n
         ff=0                           ! do nothing
         if(mol%at(i).eq. 1)                 ff=-0.08 ! H
         if(mol%at(i).eq. 5)                 ff=-0.05 ! B
         if(mol%at(i).eq. 6)                 then
                                         ff=-0.27 ! C
                         if(topo%hyb(i).lt.3) ff=-0.45 ! unsat
                         if(topo%hyb(i).lt.2) ff=-0.34 ! unsat
         endif
         if(mol%at(i).eq. 7)                 then
                                         ff=-0.13 ! N
           if(piadr(i).ne.0)             ff=-0.14
           if(amide(mol%n,mol%at,topo%hyb,topo%nb,piadr,i))ff=-0.16
         endif
         if(mol%at(i).eq. 8)                 then
                                         ff=-0.15 ! O
                         if(topo%hyb(i).lt.3) ff=-0.08 ! unsat
         endif
         if(mol%at(i).eq. 9)                 ff= 0.10 ! F
         if(mol%at(i).gt.10)                 ff=-0.02 ! heavy
         if(mol%at(i).eq.17)                 ff=-0.02 ! Cl
         if(mol%at(i).eq.35)                 ff=-0.11 ! Br
         if(mol%at(i).eq.53)                 ff=-0.07 ! I
         if(imetal(i).eq.1)              ff=-0.08 ! M maing
         if(imetal(i).eq.2)              ff=-0.9  ! M TM    ??? too large
         if(param%group(mol%at(i)).eq.8)           ff= 0.0  ! RG
         dgam(i)=topo%qa(i)*ff
      enddo

!     prepare true EEQ parameter, they are ATOMIC not element specific!
      do i=1,mol%n
!                   base val   spec. corr.
         topo%chieeq(i)=-param%chi(mol%at(i)) + dxi(i)
         topo%gameeq(i)= param%gam(mol%at(i)) +dgam(i)
         if (amideH(mol%n,mol%at,topo%hyb,topo%nb,piadr2,i)) topo%chieeq(i) = topo%chieeq(i) - 0.02
         ff = 0
         if(mol%at(i).eq.6)       ff= 0.09
         if(mol%at(i).eq.7)       ff=-0.21
         if(param%group(mol%at(i)).eq.6)ff=-0.03
         if(param%group(mol%at(i)).eq.7)ff= 0.50
         if(imetal(i).eq.1)   ff= 0.3
         if(imetal(i).eq.2)   ff=-0.1
         topo%alpeeq(i) = (param%alp(mol%at(i))+ff*topo%qa(i))**2
      enddo
      deallocate(dgam,dxi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get ring info (smallest ring size)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(env%unit,'(10x,"rings ...")')
!$omp parallel default(none) private(i,cr,sr) shared(mol,nbm,cring,sring)
!$omp do
      do i=1,mol%n
         call getring36(mol%n,mol%at,nbm,i,cr,sr)
         cring(1:10,1:20,i)=cr(1:10,1:20)
         sring(     1:20,i)=sr(1:20)
      enddo
!$omp end do
!$omp end parallel
      deallocate(nbm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! bonded atom triples not included in
! bend and tors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      idum=1000*mol%n
      allocate( topo%b3list(3,idum), source = 0 )
      topo%nbatm=0
      do i=1,mol%n
         do j=1,i-1
            ij=lin(j,i)
            if(topo%bpair(ij).eq.3) then  ! 1,4 exclusion of back-pair makes it worse, 1,3 makes little effect
            do m=1,topo%nb(20,j)
               k=topo%nb(m,j)
               topo%nbatm=topo%nbatm+1
               topo%b3list(1,topo%nbatm)=i
               topo%b3list(2,topo%nbatm)=j
               topo%b3list(3,topo%nbatm)=k
            enddo
            do m=1,topo%nb(20,i)
               k=topo%nb(m,i)
               topo%nbatm=topo%nbatm+1
               topo%b3list(1,topo%nbatm)=i
               topo%b3list(2,topo%nbatm)=j
               topo%b3list(3,topo%nbatm)=k
            enddo
            endif
         enddo
      enddo
      if(topo%nbatm.gt.idum) then
         write(env%unit,*) idum,topo%nbatm
         call env%error('overflow in ini', source)
         return
      endif
      write(env%unit,'(10x,"# BATM",3x,i0)') topo%nbatm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! non bonded pair exponents
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=1,mol%n
         ati=mol%at(i)
         fn=1.0d0 + gen%nrepscal/(1.0d0+dble(topo%nb(20,i))**2)
         dum1=param%repan(ati)*(1.d0 + topo%qa(i)*gen%qrepscal)*fn ! a small but physically correct decrease of repulsion with q
         f1=zeta(ati,topo%qa(i))
         do j=1,i-1
            atj=mol%at(j)
            fn=1.0d0 + gen%nrepscal/(1.0d0+dble(topo%nb(20,j))**2)
            dum2=param%repan(atj)*(1.d0 + topo%qa(j)*gen%qrepscal)*fn
            f2=zeta(atj,topo%qa(j))
            ij=lin(j,i)
            ff = 1.0d0
            if(ati.eq.1.and.atj.eq.1) then
               ff = 1.0d0*gen%hhfac                     ! special H ... H case (for other pairs there is no good effect of this)
               if(topo%bpair(ij).eq.3) ff=ff*gen%hh14rep     ! 1,4 case important for right torsion pot.
               if(topo%bpair(ij).eq.2) ff=ff*gen%hh13rep     ! 1,3 case
            endif
            if((ati.eq.1.and.param%metal(atj).gt.0).or.(atj.eq.1.and.param%metal(ati).gt.0)) ff=0.85 ! M...H
            if((ati.eq.1.and.atj.eq.6).or.(atj.eq.1.and.ati.eq.6))               ff=0.91 ! C...H, good effect
            if((ati.eq.1.and.atj.eq.8).or.(atj.eq.1.and.ati.eq.8))               ff=1.04 ! O...H, good effect
            topo%alphanb(ij)=sqrt(dum1*dum2)*ff
            topo%zetac6(ij)=f1*f2  ! D4 zeta scaling using qref=0
         enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make list of HB donor bascity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !atom specific (not element) basicity parameters
      allocate( topo%hbbas(mol%n), source =1.0d0 )
      do i = 1,mol%n
         nn=topo%nb(20,i)
         ati=mol%at(i)
         topo%hbbas(i)=param%xhbas(mol%at(i))
         ! Carbene:
         if(ati.eq.6.and.nn.eq.2.and.itag(i).eq.1) topo%hbbas(i) = 1.46
         ! Carbonyl R-C=O
         if(ati.eq.8.and.nn.eq.1.and.mol%at(topo%nb(nn,i)).eq.6) topo%hbbas(i) = 0.68
         ! Nitro R-N=O
         if(ati.eq.8.and.nn.eq.1.and.mol%at(topo%nb(nn,i)).eq.7) topo%hbbas(i) = 0.47
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make list of HB donor acidity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !atom specific (not element) basicity parameters
      allocate( topo%hbaci(mol%n), source =1.0d0 )
      do i = 1,mol%n
         topo%hbaci(i)=param%xhaci(mol%at(i))
      end do
      do i = 1,mol%n
         nn=topo%nb(1,i)
         topo%hbaci(i)=param%xhaci(mol%at(i))
         ! AmideH:
         if (amideH(mol%n,mol%at,topo%hyb,topo%nb,piadr2,i)) topo%hbaci(nn) = topo%hbaci(nn) * 0.80
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make list of ABs for HAB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      allocate( topo%hbatHl(mol%n),topo%hbatABl(2,mol%n*(mol%n+1)/2), source = 0 )

      topo%nathbH=0
      do i=1,mol%n
         if(mol%at(i).ne.1)  cycle
         if(topo%hyb(i).eq.1) cycle      ! exclude bridging hydrogens from HB correction
         ff=gen%hqabthr
         j=topo%nb(1,i)
         if(j.le.0) cycle
         if(mol%at(j).gt.10) ff=ff-0.20                ! H on heavy atoms may be negatively charged
         if(mol%at(j).eq.6.and.topo%hyb(j).eq.3) ff=ff+0.05 ! H on sp3 C must be really positive 0.05
         if(topo%qa(i).gt.ff)then                       ! make list of HB H atoms but only if they have a positive charge
            topo%nathbH=topo%nathbH+1
            topo%hbatHl(topo%nathbH)=i
         endif
      enddo
      write(env%unit,'(10x,"# H in HB",3x,i0)') topo%nathbH

      topo%nathbAB=0
      do i=1,mol%n
         if(mol%at(i).eq. 6.and.piadr2(i).eq.0) cycle ! C sp or sp2 pi
         ff=gen%qabthr
         if(mol%at(i).gt.10) ff=ff+0.2   ! heavy atoms may be positively charged
         if(topo%qa(i).gt.ff) cycle
         do j=1,i-1
            ff=gen%qabthr
            if(mol%at(j).gt.10) ff=ff+0.2  ! heavy atoms may be positively charged
            if(topo%qa(j).gt.ff) cycle
            call hbonds(i,j,hbpi,hbpj,param,topo)
            if(hbpi(1)*hbpj(2).lt.1.d-6.and.hbpi(2)*hbpj(1).lt.1.d-6)cycle
            if(mol%at(j).eq. 6.and.piadr2(j).eq.0) cycle ! C sp or sp2 pi
            topo%nathbAB = topo%nathbAB + 1
            topo%hbatABl(1,topo%nathbAB)=i
            topo%hbatABl(2,topo%nathbAB)=j
         enddo
      enddo

! make ABX list
      m=0
      do i=1,mol%n
         do ia=1,topo%nb(20,i)
            ix=topo%nb(ia,i)
            if(xatom(mol%at(ix))) then
            if(mol%at(ix).eq.16.and.topo%nb(20,ix).gt.2) cycle ! no sulphoxide etc S
            do j=1,mol%n
               if(i.eq.j.or.j.eq.ix) cycle
               if(topo%bpair(lin(j,ix)).le.3) cycle   ! must be A...B and not X-B i.e. A-X...B
               if(param%xhbas(mol%at(j)).lt.1.d-6) cycle   ! B must be O,N,...
               if(param%group(mol%at(j)).eq.4    ) then
                  if(piadr2(j).eq.0.or.topo%qa(j).gt.0.05) cycle   ! must be a (pi)base
               endif
               m=m+1
            enddo
            endif
         enddo
      enddo
      topo%natxbAB=m
      allocate(topo%xbatABl(3,topo%natxbAB), source = 0 )
      m=0
      do i=1,mol%n
         do ia=1,topo%nb(20,i)
            ix=topo%nb(ia,i)
            if(xatom(mol%at(ix))) then
            if(mol%at(ix).eq.16.and.topo%nb(20,ix).gt.2) cycle ! no sulphoxide etc S
            do j=1,mol%n
               if(i.eq.j.or.j.eq.ix) cycle
               if(topo%bpair(lin(j,ix)).le.3) cycle  ! must be A...B and not X-B i.e. A-X...B
               if(param%xhbas(mol%at(j)).lt.1.d-6) cycle  ! B must be O,N,...
               if(param%group(mol%at(j)).eq.4    ) then
                  if(piadr2(j).eq.0.or.topo%qa(j).gt.0.05) cycle   ! must be a (pi)base
               endif
               m=m+1
               topo%xbatABl(1,m)=i
               topo%xbatABl(2,m)=j
               topo%xbatABl(3,m)=ix
            enddo
            endif
         enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! do Hueckel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(picount.gt.0) then
      write(env%unit,'(10x,"doing iterative Hueckel for ",i0," subsystem(s) ...")') picount
      allocate( pispop(picount),pisip(picount),pisea(picount), source = 0.0d0 )
      allocate( piel(mol%n), source = 0 )
      itmp = 0 ! save pi atom info
      hcalc= 0
      pisip= 0
      pisea= 0

      if(pr) then
         write(env%unit,'(10x,"iterative Hueckel run to get P ...")')
      endif
      do pis=1,picount ! loop over pi systems
      npi   =0
      nelpi =0
      piadr3=0
      piadr4=0
      piel  =0

      do k=1,npiall
         if(pimvec(k).eq.pis) then
            npi=npi+1
            ati =mol%at (piadr(k))
            hybi=topo%hyb(piadr(k))
            ii  =nelpi
            if(ati.eq.5.and.hybi.eq.1)           nelpi=nelpi+1  ! B in borine
            if(ati.eq.6.and.itag(piadr(k)).ne.1) nelpi=nelpi+1  ! skip if its a carbene (tag itag=1)
            if(ati.eq.7.and.hybi.eq.2.and.itag(piadr(k)).eq.1) &
     &                                           nelpi=nelpi+1  ! the itag=1 avoids an odd el number for the nitro group (its 4)
            if(ati.eq.7.and.hybi.le.2)           nelpi=nelpi+1
            if(ati.eq.7.and.hybi.eq.3)           nelpi=nelpi+2
            if(ati.eq.8.and.hybi.eq.1)           nelpi=nelpi+1
            if(ati.eq.8.and.hybi.eq.2)           nelpi=nelpi+1
            if(ati.eq.8.and.hybi.eq.3)           nelpi=nelpi+2
            if(ati.eq.9.and.hybi.ne.1)           nelpi=nelpi+2
            if(ati.eq.9.and.hybi.eq.1)           nelpi=nelpi+3 !??? otherwise fluor-furan+ is wrong
            if(ati.eq.16.and.hybi.eq.1)          nelpi=nelpi+1
            if(ati.eq.16.and.hybi.eq.2)          nelpi=nelpi+1
            if(ati.eq.16.and.hybi.eq.3)          nelpi=nelpi+2
            if(ati.eq.17.and.hybi.eq.0)          nelpi=nelpi+2
            if(ati.eq.17.and.hybi.eq.1)          nelpi=nelpi+3
            piadr3(npi)=piadr(k) ! map to original, full atom set
            piadr4(piadr(k))=npi
            piel(piadr(k))=nelpi-ii
            if(piel(piadr(k)).gt.2)piel(piadr(k))=2
         endif
      enddo
      nelpi=nelpi-ipis(pis)
      if(npi.lt.2.or.nelpi.lt.1) cycle
      allocate(Api(npi,npi),apisave(npi,npi),Pold(npi,npi),S(npi,npi),occ(npi),eps(npi)) ! S is just scratch here

      eold= 0
      Pold= 2.d0/3.d0
! iterative Hueckel loop, off-diag terms are reduced depending on P to avoid overdelocalization
      do nn=1,nint(gen%maxhiter)      ! just some iterations
      Api = 0
      do i=1,npi
         ii=piadr3(i)
         Api(i,i)=gen%hdiag(mol%at(ii))+topo%qa(ii)*gen%hueckelp3-dble(piel(ii)-1)*gen%pilpf
      enddo
!     loop over bonds for pair interactions
      do i=1,topo%nbond
         ii=topo%blist(1,i)
         jj=topo%blist(2,i)
         ia=piadr4(ii)
         ja=piadr4(jj)
         if(ia.gt.0.and.ja.gt.0) then
            dum=1.d-9*rab(lin(ii,jj))                                 ! distort so that Huckel for e.g. COT localizes to right bonds
            dum=sqrt(gen%hoffdiag(mol%at(ii))*gen%hoffdiag(mol%at(jj)))-dum           ! better than arithmetic
            dum2=gen%hiter
            if(topo%hyb(ii).eq.1)                 dum2=dum2*gen%htriple        ! triple bond is different
            if(topo%hyb(jj).eq.1)                 dum2=dum2*gen%htriple        ! triple bond is different
            Api(ja,ia)=-dum  * (1.0d0-dum2*(2.0d0/3.0d0-Pold(ja,ia))) ! Pmat scaling with benzene as reference
            Api(ia,ja)=Api(ja,ia)
         endif
      enddo

      apisave = Api
      call gfnffqmsolve(.false.,Api,S,.false.,4000.0d0,npi,0,nelpi,dum,occ,eps)  !diagonalize, 4000 better than 300

      do i=1,npi  ! save IP/EA
         if(occ(i).gt.0.5) then
                       pisip(pis)=eps(i)   ! IP
         if(i+1.lt.npi)pisea(pis)=eps(i+1) ! EA
         endif
      enddo
      if(abs(dum-eold).lt.1.d-4) exit  ! end of iterations
      Pold = Api
      eold = dum
      enddo
! end of iterative loop
      if(pr)then
         write(env%unit,'(''Hueckel system :'',i3,'' charge : '',i3,'' ndim/Nel :'',2i5, &
     &         3x, ''eps(HOMO/LUMO)'',2f12.6)')pis,ipis(pis),npi,nelpi,pisip(pis),pisea(pis)
         if(pisip(pis).gt.0.40) then
            write(env%unit,*)'WARNING: probably wrong pi occupation. Second attempt with Nel=Nel-1!'
            do i=1,mol%n
               if(piadr4(i).ne.0) write(env%unit,*) 'at,nb,topo%hyb,Npiel:', i,mol%sym(i),topo%nb(20,i),topo%hyb(i),piel(i)
            enddo
            nelpi=nelpi-1
            Api = Apisave
            call gfnffqmsolve(.false.,Api,S,.false.,300.0d0,npi,0,nelpi,dum,occ,eps)  !diagonalize
            call PREIG(6,occ,1.0d0,eps,1,npi)
         endif
      endif
! save BO
      do i=1,topo%nbond
         ii=topo%blist(1,i)
         jj=topo%blist(2,i)
         ia=piadr4(ii)
         ja=piadr4(jj)
         if(ia.gt.0.and.ja.gt.0) then
            pibo(i)=Api(ja,ia)
            pbo(lin(ii,jj))=Api(ja,ia)
            itmp(ii)=1
            itmp(jj)=1
         endif
      enddo
      deallocate(Api,apisave,Pold,S,occ,eps)
      enddo
! end of pi system loop
      piadr = itmp  ! array used for identifying pi atoms in following codes
      deallocate(pispop,pisip,pisea,ipis,pimvec,piel)
      endif
!----------- end Hueckel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! modify hyb due to pi assignment
! and output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=1,mol%n
         if(topo%hyb(i).eq.2 .and. piadr(i).eq.0 .and. topo%nb(20,i).eq.3 .and. param%group(mol%at(i)).eq.4)then ! C,Si,Ge... CN=3, no pi
            jj=topo%nb(1,i)
            kk=topo%nb(2,i)
            ll=topo%nb(3,i)
            phi=omega(mol%n,mol%xyz,i,jj,kk,ll)  ! the shitty second geom. dep. term GEODEP
            if(abs(phi)*180./pi.gt.40.d0) topo%hyb(i) = 3  ! change to sp^3
         endif
      enddo

!     if(pr)then
      write(env%unit,*)
      write(env%unit,'(2x,"atom   neighbors  erfCN metchar sp-hybrid imet pi  qest     coordinates")')
      do i=1,mol%n
         j = topo%hyb(i)
         if(amide(mol%n,mol%at,topo%hyb,topo%nb,piadr,i))  j=-topo%hyb(i)
         if(mol%at(i).eq.6.and.itag(i).eq.1) j=-topo%hyb(i)
         write(env%unit,'(i5,2x,a2,3x,i4,3x,f5.2,2x,f5.2,8x,i2,3x,i2,3x,i2,2x,f6.3,3f12.6)') &
     &             i,mol%sym(i),topo%nb(20,i),cn(i),mchar(i),j,imetal(i),piadr(i),topo%qa(i),mol%xyz(1:3,i)
      enddo

!     compute fragments and charges for output (check for CT)
!     call mrecgff(mol%n,topo%nb,nmol,piadr3)
!     write(env%unit,*) 'Nmol',nmol
      if(pr)then
      write(env%unit,'(/,''molecular fragment  # atoms  topo charge'')')
      do i=1,topo%nfrag
         dum=0
           m=0
         do k=1,mol%n
            if(topo%fraglist(k).eq.i) then
               m=m+1
               dum=dum+topo%qa(k)
            endif
         enddo
         write(env%unit,'(5x,i3,10x,i4,10x,f8.3)')i,m,dum
      enddo
      write(env%unit,*)
      endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!            bonds
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call gfnffrab(mol%n,mol%at,cn,rtmp)           ! guess RAB for output

      topo%nbond_vbond = topo%nbond
      allocate( topo%vbond(3,topo%nbond), source = 0.0d0 )

      write(env%unit,*)
      write(env%unit,'(10x,"#atoms :",3x,i0)') mol%n
      write(env%unit,'(10x,"#bonds :",3x,i0)') topo%nbond
      if(pr)then
      write(env%unit,*)
      write(env%unit,*) 'bond atoms        type  in ring    R      R0    piBO    fqq  kbond(tot)  alp'
      endif

      do i=1,topo%nbond
         ii=topo%blist(1,i)
         jj=topo%blist(2,i)
         ij=lin(ii,jj)
         ia=mol%at(ii)
         ja=mol%at(jj)
         call ringsbond(mol%n,ii,jj,cring,sring,rings)
         shift =0.d0
         fxh   =1.d0
         ringf =1.d0
         fqq   =1.d0
         fpi   =1.d0
         fheavy=1.d0
         fheavy=1.d0
         fcn   =1.d0
         fsrb2 =gen%srb2
         bridge=.false.
         shift =0.d0
! assign bond type
                                                                  btyp(i)=1 ! single
         if(topo%hyb(ii).eq.2.and.topo%hyb(jj).eq.2)                        btyp(i)=2 ! sp2-sp2 = pi
         if(topo%hyb(ii).eq.3.and.topo%hyb(jj).eq.2.and.ia.eq.7)            btyp(i)=2 ! N-sp2
         if(topo%hyb(jj).eq.3.and.topo%hyb(ii).eq.2.and.ja.eq.7)            btyp(i)=2 ! N-sp2
         if(topo%hyb(ii).eq.1.or. topo%hyb(jj).eq.1)                        btyp(i)=3 ! sp-X i.e. no torsion
         if((param%group(ia).eq.7.or.ia.eq.1).and.topo%hyb(ii).eq.1)then
                                                                  btyp(i)=3 ! linear halogen i.e. no torsion
                                                                  bridge=.true.
         endif
         if((param%group(ja).eq.7.or.ja.eq.1).and.topo%hyb(jj).eq.1)then
                                                                  btyp(i)=3 ! linear halogen i.e. no torsion
                                                                  bridge=.true.
         endif
         if(topo%hyb(ii).eq.5.or. topo%hyb(jj).eq.5)                        btyp(i)=4 ! hypervalent
         if(imetal(ii).gt.0.or.imetal(jj).gt.0)                   btyp(i)=5 ! metal
         if(imetal(ii).eq.2.and.imetal(jj).eq.2)                  btyp(i)=7 ! TM metal-metal
         if(imetal(jj).eq.2.and.itag(ii).eq.-1.and.piadr(ii).gt.0)btyp(i)=6 ! eta
         if(imetal(ii).eq.2.and.itag(jj).eq.-1.and.piadr(jj).gt.0)btyp(i)=6 ! eta
         bbtyp=btyp(i)
! normal bond
         if(bbtyp.lt.5)then
            hybi=max(topo%hyb(ii),topo%hyb(jj))
            hybj=min(topo%hyb(ii),topo%hyb(jj))
            if(hybi.eq.5.or.hybj.eq.5) then
            bstrength=gen%bstren(4)                                       ! base value hypervalent
            else
            bstrength=gen%bsmat(hybi,hybj)                                ! base value normal hyb
            endif
            if(hybi.eq.3.and.hybj.eq.2.and.(ia.eq.7.or.ja.eq.7)) &
     &                                      bstrength=gen%bstren(2)*1.04   ! N-sp2

            if(bridge)then
               if(param%group(ia).eq.7)           bstrength=gen%bstren(1)*0.50d0 ! bridging X
               if(param%group(ja).eq.7)           bstrength=gen%bstren(1)*0.50d0 ! bridging X
               if(ia.eq.1.or.ia.eq.9)       bstrength=gen%bstren(1)*0.30d0 ! bridging H/F
               if(ja.eq.1.or.ja.eq.9)       bstrength=gen%bstren(1)*0.30d0 ! bridging H/F
            endif
            if(bbtyp.eq.4)                  shift=gen%hyper_shift          ! hypervalent
            if(ia.eq.1.or.ja.eq.1)          shift=gen%rabshifth            ! XH
            if(ia.eq.9.and.ja.eq.9)         shift=0.22                 ! f2
            if(topo%hyb(ii).eq.3.and.topo%hyb(jj).eq.0)shift=shift-0.022         ! X-sp3
            if(topo%hyb(ii).eq.0.and.topo%hyb(jj).eq.3)shift=shift-0.022         ! X-sp3
            if(topo%hyb(ii).eq.1.and.topo%hyb(jj).eq.0)shift=shift+0.14          ! X-sp
            if(topo%hyb(ii).eq.0.and.topo%hyb(jj).eq.1)shift=shift+0.14          ! X-sp
            if( (ia.eq.1.and.ja.eq.6) )then
                           call ringsatom(mol%n,jj,cring,sring,ringsj)
                           if(ringsj.eq.3)                 fxh=1.05    ! 3-ring CH
                           if(ctype(mol%n,mol%at,topo%nb,piadr,jj).eq.1)fxh=0.95    ! aldehyd CH
            endif
            if( (ia.eq.6.and.ja.eq.1) )then
                           call ringsatom(mol%n,ii,cring,sring,ringsi)
                           if(ringsi.eq.3)                 fxh=1.05    ! 3-ring CH
                           if(ctype(mol%n,mol%at,topo%nb,piadr,ii).eq.1)fxh=0.95    ! aldehyd CH
            endif
            if( (ia.eq.1.and.ja.eq.5) )                    fxh=1.10    ! BH
            if( (ja.eq.1.and.ia.eq.5) )                    fxh=1.10    !
            if( (ia.eq.1.and.ja.eq.7) )                    fxh=1.06    ! NH
            if( (ja.eq.1.and.ia.eq.7) )                    fxh=1.06    !
            if( (ia.eq.1.and.ja.eq.8) )                    fxh=0.93    ! OH
            if( (ja.eq.1.and.ia.eq.8) )                    fxh=0.93    !
            if(bbtyp.eq.3.and.ia.eq.6.and.ja.eq.8) bstrength=gen%bstren(3)*0.90d0 ! makes CO right and M-CO reasonable
            if(bbtyp.eq.3.and.ia.eq.8.and.ja.eq.6) bstrength=gen%bstren(3)*0.90d0 !
!           modify locally for triple bonds
            if( bbtyp.eq.3.and. (topo%hyb(ii).eq.0.or.topo%hyb(jj).eq.0) ) bbtyp=1 ! sp-sp3
            if( bbtyp.eq.3.and. (topo%hyb(ii).eq.3.or.topo%hyb(jj).eq.3) ) bbtyp=1 ! sp-sp3
            if( bbtyp.eq.3.and. (topo%hyb(ii).eq.2.or.topo%hyb(jj).eq.2) ) bbtyp=2 ! sp-sp2
!           Pi stuff
            if(pibo(i).gt.0) then
                          shift=gen%hueckelp*(gen%bzref - pibo(i)) ! ref value = no correction is benzene, P=2/3
                          if(bbtyp.ne.3.and.pibo(i).gt.0.1) then
                                                            btyp(i)=2
                                                            bbtyp=2
                          endif
                          fpi=1.0d0-gen%hueckelp2*(gen%bzref2 - pibo(i)) ! deepness
            endif
            if(ia.gt.10.and.ja.gt.10)then
              fcn=fcn/(1.0d0+0.007*dble(topo%nb(20,ii))**2)
              fcn=fcn/(1.0d0+0.007*dble(topo%nb(20,jj))**2)
            endif
            qafac=topo%qa(ii)*topo%qa(jj)*70.0d0
            fqq=1.0d0+gen%qfacbm0*exp(-15.d0*qafac)/(1.0d0+exp(-15.d0*qafac))
! metal involed
         else
            shift=0
            bstrength=gen%bstren(bbtyp)
            if(bbtyp.eq.7)then ! TM-TM
               if(itabrow6(ia).gt.4.and.itabrow6(ja).gt.4) bstrength=gen%bstren(8) ! 4/5d-4/5d
               if(itabrow6(ia).eq.4.and.itabrow6(ja).gt.4) bstrength=gen%bstren(9) ! 3d-4/5d
               if(itabrow6(ja).eq.4.and.itabrow6(ia).gt.4) bstrength=gen%bstren(9) ! 3d-4/5d
               dum=2.0d0*mchar(ii)+2.0d0*mchar(jj)
               dum=min(dum,0.5d0)  ! limit the "metallic" correction
               bstrength=bstrength*(1.0d0-dum)
            endif
                                                  mtyp1=0  ! no metal
                                                  mtyp2=0
            if(param%group(ia).eq.1)                    mtyp1=1  ! Li...
            if(param%group(ia).eq.2)                    mtyp1=2  ! Be...
            if(param%group(ia).gt.2.and.imetal(ii).eq.1)mtyp1=3  ! main group
            if(imetal(ii).eq.2)                   mtyp1=4  ! TM
            if(param%group(ja).eq.1)                    mtyp2=1  ! Li...
            if(param%group(ja).eq.2)                    mtyp2=2  ! Be...
            if(param%group(ja).gt.2.and.imetal(jj).eq.1)mtyp2=3  ! main group
            if(imetal(jj).eq.2)                   mtyp2=4  ! TM
            qafac=topo%qa(ii)*topo%qa(jj)*25.0d0
            dum=exp(-15.d0*qafac)/(1.0d0+exp(-15.d0*qafac))
            fqq=1.0d0+dum * (gen%qfacbm(mtyp1)+gen%qfacbm(mtyp2))*0.5   ! metal charge corr.
            if(imetal(ii).eq.2.and.ja.gt.10)      fheavy=0.65d0 ! heavy gen. ligand
            if(imetal(jj).eq.2.and.ia.gt.10)      fheavy=0.65d0
            if(imetal(ii).eq.2.and.ja.eq.15)      fheavy=1.60d0 ! P ligand
            if(imetal(jj).eq.2.and.ia.eq.15)      fheavy=1.60d0
            if(imetal(ii).eq.2.and.param%group(ja).eq.6)fheavy=0.85d0 ! chalcogen ligand
            if(imetal(jj).eq.2.and.param%group(ia).eq.6)fheavy=0.85d0
            if(imetal(ii).eq.2.and.param%group(ja).eq.7)fheavy=1.30d0 ! halogen ligand
            if(imetal(jj).eq.2.and.param%group(ia).eq.7)fheavy=1.30d0
            if(imetal(ii).eq.2.and.ja.eq.1.and.itabrow6(ia).le.5) fxh =0.80d0 ! hydrogen 3d/4d
            if(imetal(jj).eq.2.and.ia.eq.1.and.itabrow6(ja).le.5) fxh =0.80d0 ! hydrogen 3d/4d
            if(imetal(ii).eq.2.and.ja.eq.1.and.itabrow6(ia).gt.5) fxh =1.00d0 ! hydrogen 5d
            if(imetal(jj).eq.2.and.ia.eq.1.and.itabrow6(ja).gt.5) fxh =1.00d0 ! hydrogen 5d
            if(imetal(ii).eq.1.and.ja.eq. 1)      fxh   =1.20d0
            if(imetal(jj).eq.1.and.ia.eq. 1)      fxh   =1.20d0
            if(imetal(jj).eq.2.and.topo%hyb(ii).eq.1) then !CO/CN/NC...
                                       if(ia.eq.6)then
                                                  fpi   =1.5d0
                                                  shift=-0.45d0
                                       endif
                                       if(ia.eq.7.and.topo%nb(20,ii).ne.1)then
                                                  fpi   =0.4d0
                                                  shift= 0.47d0
                                       endif
            endif
            if(imetal(ii).eq.2.and.topo%hyb(jj).eq.1) then !CO/CN/NC...
                                       if(ja.eq.6)then
                                                  fpi   =1.5d0
                                                  shift=-0.45d0
                                       endif
                                       if(ja.eq.7.and.topo%nb(20,jj).ne.1)then
                                                  fpi   =0.4d0
                                                  shift= 0.47d0
                                       endif
            endif
            if(imetal(ii).eq.2)                   shift=shift+gen%metal2_shift   ! metal shift TM
            if(imetal(jj).eq.2)                   shift=shift+gen%metal2_shift   !
            if(imetal(ii).eq.1.and.param%group(ia).le.2)shift=shift+gen%metal1_shift   ! metal shift group 1+2
            if(imetal(jj).eq.1.and.param%group(ja).le.2)shift=shift+gen%metal1_shift   !
            if(mtyp1     .eq.3)                   shift=shift+gen%metal3_shift   ! metal shift MG
            if(mtyp2     .eq.3)                   shift=shift+gen%metal3_shift   !
            if(bbtyp.eq.6.and.param%metal(ia).eq.2)     shift=shift+gen%eta_shift*topo%nb(20,ii)! eta coordinated
            if(bbtyp.eq.6.and.param%metal(ja).eq.2)     shift=shift+gen%eta_shift*topo%nb(20,jj)! eta coordinated
            if(mtyp1.gt.0.and.mtyp1.lt.3) fcn=fcn/(1.0d0+0.100*dble(topo%nb(20,ii))**2)
            if(mtyp2.gt.0.and.mtyp2.lt.3) fcn=fcn/(1.0d0+0.100*dble(topo%nb(20,jj))**2)
            if(mtyp1.eq.3)                fcn=fcn/(1.0d0+0.030*dble(topo%nb(20,ii))**2)
            if(mtyp2.eq.3)                fcn=fcn/(1.0d0+0.030*dble(topo%nb(20,jj))**2)
            if(mtyp1.eq.4)                fcn=fcn/(1.0d0+0.036*dble(topo%nb(20,ii))**2)
            if(mtyp2.eq.4)                fcn=fcn/(1.0d0+0.036*dble(topo%nb(20,jj))**2)
            if(mtyp1.eq.4.or.mtyp2.eq.4)then
              fsrb2=-gen%srb2*0.22! weaker, inverse EN dep. for TM metals
            else
              fsrb2= gen%srb2*0.28! "normal" for other metals
            endif
         endif

         if(ia.gt.10.and.ja.gt.10) then  ! both atoms are heavy
             shift = shift + gen%hshift3
             if(ia.gt.18) shift = shift + gen%hshift4
             if(ja.gt.18) shift = shift + gen%hshift4
             if(ia.gt.36) shift = shift + gen%hshift5
             if(ja.gt.36) shift = shift + gen%hshift5
         endif

! shift
         topo%vbond(1,i) = gen%rabshift + shift   ! value for all bonds + special part

! RINGS prefactor
         if(rings.gt.0) ringf = 1.0d0 + gen%fringbo*(6.0d0-dble(rings))**2  ! max ring size is 6

! steepness
         topo%vbond(2,i) =  gen%srb1*( 1.0d0 + fsrb2*(param%en(ia)-param%en(ja))**2 + gen%srb3*bstrength )

! tot prefactor        atoms              spec     typ       qterm    heavy-M  pi   XH(3ring,OH...) CN for M
         topo%vbond(3,i) = -param%bond(ia)*param%bond(ja) * ringf * bstrength * fqq * fheavy * fpi * fxh * fcn
!        write(env%unit,*) bond(ia),bond(ja),ringf,bstrength,fqq,fheavy,fpi,fxh

! output
         r0 = (rtmp(ij)+topo%vbond(1,i))*0.529167
         if(pr) write(env%unit,'(2a3,2i5,2x,2i5,2x,6f8.3)') &
     &   mol%sym(ii),mol%sym(jj),ii,jj,bbtyp,rings,0.529167*rab(ij),r0,pibo(i),fqq,topo%vbond(3,i),topo%vbond(2,i)
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     scale FC if bond is part of hydrogen bridge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      do i=1,topo%nbond
!         ii=topo%blist(1,i)
!         jj=topo%blist(2,i)
!         ia=mol%at(ii)
!         ja=mol%at(jj)
!         do j=1,topo%nathbAB
!            hbA=topo%hbatABl(2,j)
!            !O-H case:
!            if (ia.eq.8.and.ja.eq.1.or.ia.eq.1.and.ja.eq.8) then
!               if (ii.eq.hbA.or.jj.eq.hbA) then
!                  topo%vbond(2,i) = topo%vbond(2,i) * 1.00
!                  !topo%vbond(1,i) = topo%vbond(1,i) * 1.00
!               end if
!            !N-H case
!          else if (ia.eq.7.and.ja.eq.1.or.ia.eq.1.and.ja.eq.7) then
!               if (ii.eq.hbA.or.jj.eq.hbA) then
!                  topo%vbond(2,i) = topo%vbond(2,i) * 1.00
!                  !topo%vbond(1,i) = topo%vbond(1,i) * 1.00
!               end if
!            end if
!         end do
!      end do

      deallocate(rtmp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     scale FC if bond is part of hydrogen bridge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Set up fix hblist just like for the HB term
      call bond_hbset0(mol%n,mol%at,mol%xyz,sqrab,bond_hbn,topo,hbthr1,hbthr2)
      allocate(bond_hbl(3,bond_hbn))
      allocate(topo%nr_hb(topo%nbond), source=0)
      call bond_hbset(mol%n,mol%at,mol%xyz,sqrab,bond_hbn,bond_hbl,topo,hbthr1,hbthr2)

      !Set up AH, B and nr. of B list
      call bond_hb_AHB_set0(mol%n,mol%at,topo%nbond,bond_hbn,bond_hbl,AHB_nr,topo)
      allocate( lin_AHB(0:AHB_nr), source=0  )
      call bond_hb_AHB_set1(mol%n,mol%at,topo%nbond,bond_hbn,bond_hbl,AHB_nr,lin_AHB,topo%bond_hb_nr,topo%b_max,topo)
      allocate( topo%bond_hb_AH(2,topo%bond_hb_nr), source = 0 )
      allocate( topo%bond_hb_B(topo%b_max,topo%bond_hb_nr), source = 0 )
      allocate( topo%bond_hb_Bn(topo%bond_hb_nr), source = 0 )
      call bond_hb_AHB_set(mol%n,mol%at,topo%nbond,bond_hbn,bond_hbl,AHB_nr,lin_AHB,topo)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!               bend
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      topo%nangl=0
      do i=1,mol%n
         nn=topo%nb(20,i)                 ! take full set to include M-X-Y
         if(nn.le.1) cycle
         if(topo%nb(20,i).gt.6) cycle     ! no highly coordinated atom
         ati=mol%at(i)
         do j=1,nn
            do k=1,j-1
               jj=topo%nb(j,i)
               kk=topo%nb(k,i)
               atj=mol%at(jj)
               atk=mol%at(kk)
               fijk=param%angl(ati)*param%angl2(atj)*param%angl2(atk)
               if(fijk.lt.gen%fcthr) cycle     ! too small
               topo%nangl=topo%nangl+1
            enddo
         enddo
      enddo

      write(env%unit,'(10x,"#angl  :",3x,i0)') topo%nangl
      if(pr)then
      write(env%unit,*)
      write(env%unit,*) 'angle atoms        phi0    phi      FC  pi rings'
      endif

      topo%nangl_alloc = topo%nangl
      allocate( topo%alist(3,topo%nangl), source = 0 )
      allocate( topo%vangl(2,topo%nangl), source = 0.0d0 )
      topo%nangl=0
      do i=1,mol%n
         nn=topo%nb(20,i)
         if(nn.le.1) cycle
         if(topo%nb(20,i).gt.6) cycle
         ii=i
         ati=mol%at(i)
         do j=1,nn
            do k=1,j-1
               jj=topo%nb(j,i)
               kk=topo%nb(k,i)
               atj=mol%at(jj)
               atk=mol%at(kk)
               fijk=param%angl(ati)*param%angl2(atj)*param%angl2(atk)
               if(fijk.lt.gen%fcthr) cycle     ! too small
               call bangl(mol%xyz,jj,i,kk,phi)
               if(param%metal(ati).gt.0.and.phi*180./pi.lt.60.) cycle ! skip eta cases even if CN < 6 (e.g. CaCp+)
               feta=1.0d0
               if(imetal(ii).eq.2.and.itag(jj).eq.-1.and.piadr(jj).gt.0) feta=0.3d0       ! eta coord.
               if(imetal(ii).eq.2.and.itag(kk).eq.-1.and.piadr(kk).gt.0) feta=feta*0.3d0  !
               nh =0
               if(atj.eq.1)nh=nh+1
               if(atk.eq.1)nh=nh+1
               nnn=0
               if(atj.eq.7)nnn=nnn+1
               if(atk.eq.7)nnn=nnn+1
               no=0
               if(atj.eq.8)no=no+1
               if(atk.eq.8)no=no+1
               nheav=0
               if(atj.gt.14)nheav=nheav+1
               if(atk.gt.14)nheav=nheav+1
               nsi=0
               if(atj.eq.14)nsi=nsi+1
               if(atk.eq.14)nsi=nsi+1
               nc=0
               if(atj.eq.6)nc=nc+1
               if(atk.eq.6)nc=nc+1
               nmet=0
               if(param%metal(atj).ne.0)nmet=nmet+1
               if(param%metal(atk).ne.0)nmet=nmet+1
               npi=0
               if(piadr(jj).ne.0)npi=npi+1
               if(piadr(kk).ne.0)npi=npi+1
               topo%nangl=topo%nangl+1
               topo%alist(1,topo%nangl)=ii
               topo%alist(2,topo%nangl)=jj
               topo%alist(3,topo%nangl)=kk
               call ringsbend(mol%n,ii,jj,kk,cring,sring,rings)
               triple=(topo%hyb(ii).eq.1 .or. topo%hyb(jj).eq.1) .or. &
     &                (topo%hyb(ii).eq.1 .or. topo%hyb(kk).eq.1)
               if(imetal(ii).eq.0.and.imetal(jj).eq.0.and.imetal(kk).eq.0) then
               fqq=1.0d0-(topo%qa(ii)*topo%qa(jj)+topo%qa(ii)*topo%qa(kk))*gen%qfacBEN      ! weaken it
               else
               fqq=1.0d0-(topo%qa(ii)*topo%qa(jj)+topo%qa(ii)*topo%qa(kk))*gen%qfacBEN*2.5
               endif
               f2 =1.0d0
               fn =1.0d0

!-------------------------
! definitions come here
!-------------------------

!!!!!!!!!!
! DEFAULT
!!!!!!!!!!
               r0=100.0

               if(topo%hyb(i).eq.1)                                   r0=180.
               if(topo%hyb(i).eq.2)                                   r0=120.
               if(topo%hyb(i).eq.3)                                   r0=109.5
               if(topo%hyb(i).eq.3.and.mol%at(i).gt.10) then
                                               if(nn.le.3)       r0=gen%aheavy3    ! heavy maingroup three coordinated
                                               if(nn.ge.4)       r0=gen%aheavy4    ! heavy maingroup four  coordinated
                                if(nn.eq.4.and.param%group(ati).eq.5)  r0=109.5      ! four coordinated group 5
                  if(nn.eq.4.and.param%group(ati).eq.4.and.ati.gt.49)  r0=109.5      ! four coordinated Sn, Pb
                                           if(param%group(ati).eq.4)   r0=r0-nh*5.   ! smaller angles for XHn Si...
                                           if(param%group(ati).eq.5)   r0=r0-nh*5.   ! smaller angles for XHn P..
                                           if(param%group(ati).eq.6)   r0=r0-nh*5.   ! smaller angles for XHn S..
               endif
               if(topo%hyb(i).eq.5)                                   then
                                                                 r0=90.
                                                                 f2=0.11       ! not very important
                                       if(phi*180./pi.gt.gen%linthr) r0=180.       ! hypervalent coordination can be linear GEODEP
               endif
!!!!!!!!!!
! B
!!!!!!!!!!
               if(ati.eq.5)then
                  if(topo%hyb(i).eq.3)                                r0=115.
                  if(topo%hyb(i).eq.2)                                r0=115.
               endif
!!!!!!!!!!
! C cases
!!!!!!!!!!
               if(ati.eq.6)then
                  if(topo%hyb(i).eq.3.and.nh.eq.2)                    r0=108.6  ! CHH
                  if(topo%hyb(i).eq.3.and.no.eq.1)                    r0=108.5  ! COR
                  if(topo%hyb(i).eq.2.and.no.eq.2)                    r0=122.   ! COO
                  if(topo%hyb(i).eq.2.and.no.eq.1)                    f2=0.7    ! C=O
                  if(topo%hyb(i).eq.1.and.no.eq.2)                    then
                                                                 triple=.false.   ! CO2
                                                                 f2=2.0
                  endif
                  if(topo%hyb(i).eq.3.and.nn.gt.4)                    then
                                       if(phi*180./pi.gt.gen%linthr) r0=180.       ! hypervalent coordination can be linear GEODEP
                  endif
               endif
!!!!!!!!!!
! O cases
!!!!!!!!!!
               if(ati.eq.8.and.nn.eq.2) then
                                                                 r0=104.5
!                   H2O
                    if(nh.eq.2)                                  then
                                                                 r0=100. ! compensate ES of the Hs
                                                                 f2=1.20 ! H2O is better with 1.2-1.3 but H2O in fit behaves differently
                                                                 endif
                                                                 r0 = r0 + 7.*nsi   ! O angles widen with Si attached
                                                                 r0 = r0 +14.*nmet  ! O angles widen with M attached
                    if(npi.eq.2)                                 then
                                                                 r0=109. ! e.g. Ph-O-Ph
                                                                 endif
                    if(nmet.gt.0.and.phi*180./pi.gt.gen%linthr)      then
                                                                 r0=180. ! metal coordination can be linear GEODEP
                                                                 f2 = 0.3
                                                                 endif
               endif
!!!!!!!!!!
! N cases
!!!!!!!!!!
               if(ati.eq.7.and.nn.eq.2) then
                                                                  f2=1.4
                                                                  r0=115.
                                    if(rings.ne.0)                r0=105.
                                    if(mol%at(kk).eq.8.or.mol%at(jj).eq.8)r0=103.
                                    if(mol%at(kk).eq.9.or.mol%at(jj).eq.9)r0=102.
                                    if(topo%hyb(i).eq.1)               r0=180.   ! NC or NNN
               if(imetal(jj).eq.2.and.topo%hyb(i).eq.1.and.mol%at(kk).eq.7)r0=135.   ! NN on M
               if(imetal(kk).eq.2.and.topo%hyb(i).eq.1.and.mol%at(jj).eq.7)r0=135.   ! NN on M
               endif
! NR3
               if(ati.eq.7.and.topo%hyb(i).eq.3)then
!                 in pi system
                  if(npi.gt.0)then
                                                                 if(amide(mol%n,mol%at,topo%hyb,topo%nb,piadr,i))then
                                                                   r0=115.
                                                                   f2=1.2d0
                                                                 else
                                                                   sumppi=pbo(lin(ii,jj))+pbo(lin(ii,kk))
                                                                   r0=113.
                                                                   f2=1.d0-sumppi*0.7d0 ! must be -!
                                                                 endif
                  else
                                                                 r0=104. ! sat. pyr. N, steep around 106
                                                                 f2=0.40 ! 1.0 is better for NH3
                                                                 f2=f2+nh*0.19
                                                                 f2=f2+no*0.25
                                                                 f2=f2+nc*0.01
                  endif
               endif
!!!!!!!!!!
! RING < 5
!!!!!!!!!!
               if(rings.eq.3)                                    r0=82. ! 60 gives too little strain
               if(rings.eq.4)                                    r0=96.
               if(rings.eq.5.and.ati.eq.6)                       r0=109.
!!!!!!!!!!
! specials
!!!!!!!!!!
! R-X in 3-rings e.g. cyclopropene
               if(rings.eq.0)then
                  call ringsatom(mol%n,i,cring,sring,idum)
                  if(idum.eq.3) then
                     call ringsatom(mol%n,jj,cring,sring,ringsj)
                     call ringsatom(mol%n,kk,cring,sring,ringsk)
                     if(ringsj+ringsk.eq.102)                    r0=r0+4.d0
                  endif
               endif

! triple bonds
               if (triple) then
                                                                 f2 = 0.60d0  ! complex 7 in S30L makes artificial torsions if this is 0.4 which is
                                                                              ! slightly better for the phenylmethylethyne bending pot.
                 if(atj.eq.7.or.atk.eq.7)                        f2 = 1.00d0
                 if((imetal(jj).eq.2.or.imetal(kk).eq.2).and.phi*180./pi.gt.gen%linthr) then
                                 if(ati.eq.6.and.atj.eq.6)       f2=3.   ! M-CC
                                 if(ati.eq.6.and.atk.eq.6)       f2=3.   ! M-CC
                                 if(ati.eq.6.and.atj.eq.7)       f2=3.   ! M-CN
                                 if(ati.eq.6.and.atk.eq.7)       f2=3.   ! M-CN
                                 if(ati.eq.6.and.param%group(atj).eq.6)f2=14.  ! M-CO or CS
                                 if(ati.eq.6.and.param%group(atk).eq.6)f2=14.  ! M-CO or CS
                                 if(ati.eq.7.and.atj.eq.7)       f2=10.  ! M-NN
                                 if(ati.eq.7.and.atj.eq.6)       f2=10.  ! M-NC
                                 if(ati.eq.7.and.atk.eq.6)       f2=10.  ! M-NC
                                 if(ati.eq.7.and.atj.eq.8)then; r0=180.;f2=12.; endif  ! M-NO
                                 if(ati.eq.7.and.atk.eq.8)then; r0=180.;f2=12.; endif  ! M-NO
                                 endif
               endif
! carbene analogous
               if(param%group(ati).eq.4.and.nn.eq.2.and.itag(i).eq.1)  then
                                                                 if(ati.eq.6) r0=145.
                                                                 if(ati.gt.6) r0= 90.
               endif
! SO3X
               if(param%group(ati).eq.6.and.nn.eq.4.and.no.ge.1)       r0=115.
! halogens CN=2
               if(param%group(ati).eq.7.and.topo%hyb(i).eq.1)               then
                                                                 if(ati.eq.9)  r0= 90.
                                                                 if(ati.eq.17) r0= 90.
                                                                 if(ati.eq.35) r0= 90.
                                                                 if(ati.eq.53) r0= 90.
                                        if(ati.gt.9.and.phi*180./pi.gt.gen%linthr) r0=180. ! change to linear if linear coordinated, GEODEP
                                                                 f2=0.6/dble(ati)**0.15
               endif
! PB or Sn can be pyramidal
               if(topo%hyb(i).eq.3.and.param%group(ati).eq.4.and.ati.gt.32.and.topo%qa(i).gt.0.4)then
                  if(phi*180./pi.gt.140.) then
                     r0=180. ! change to linear
                  endif
                  if(phi*180./pi.lt.100.) then
                     r0=90.
                  endif
                  f2=1.0
               endif
! METAL
               if(imetal(ii).gt.0) then
                                   if(topo%hyb(i).eq.0)               then
                                                                 r0=90.
                                                                 f2=1.35  ! important difference to other bends, big effect 1.15,1.25,1.35
                                   endif
                                   if(topo%hyb(i).eq.1)               r0=180.
                                   if(topo%hyb(i).eq.2)               r0=120.
                                   if(topo%hyb(i).eq.3)               r0=109.5
                                   if(phi*180./pi.gt.gen%linthr)     r0=180. ! change to linear
               endif

               fn=1.0d0 - 2.36d0/dble(nn)**2

!----------------------
! end of definitions
!----------------------
               topo%vangl(1,topo%nangl)=r0*pi/180.
               fbsmall=(1.0d0-gen%fbs1*exp(-0.64*(topo%vangl(1,topo%nangl)-pi)**2))

!                          central*neigbor charge spec. met.  small angle corr.
               topo%vangl(2,topo%nangl)= fijk * fqq * f2 * fn * fbsmall * feta
!              write(env%unit,*) param%angl(ati),param%angl2(atj),param%angl2(atk), param%angl(ati)*param%angl2(atj)*param%angl2(atk), fqq,f2,fn,fbsmall
               if(pr)write(env%unit,'(3i5,2x,3f8.3,l2,i4)') ii,jj,kk,r0,phi*180./pi,topo%vangl(2,topo%nangl),picon,rings
            enddo
         enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!              torsion
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      topo%ntors=sum(piadr)+mol%n
      do m=1,topo%nbond
         ii=topo%blist(1,m)
         jj=topo%blist(2,m)
         if(btyp(m).eq.3.or.btyp(m).eq.6)           cycle ! no sp-sp or metal eta
         if(param%tors(mol%at(ii)).lt.0.or.param%tors(mol%at(jj)).lt.0) cycle ! no negative values
         if(param%tors(mol%at(ii))*param%tors(mol%at(jj)).lt.1.d-3)     cycle ! no small values
         if(param%metal(mol%at(ii)).gt.1.and.topo%nb(20,ii).gt.4)  cycle ! no HC metals
         if(param%metal(mol%at(jj)).gt.1.and.topo%nb(20,jj).gt.4)  cycle !
         topo%ntors=topo%ntors+topo%nb(20,ii)*topo%nb(20,jj)*2 ! upper limit
      enddo
      maxtors=topo%ntors
      if(pr) write(env%unit,*) 'torsion atoms        nrot   rings    phi0    phi      FC'

      topo%ntors_alloc = topo%ntors
      allocate( topo%tlist(5,topo%ntors), source = 0 )
      allocate( topo%vtors(2,topo%ntors), source = 0.0d0 )
      topo%ntors=0
      do m=1,topo%nbond
         ii=topo%blist(1,m)
         jj=topo%blist(2,m)
         if(btyp(m).eq.3.or.btyp(m).eq.6) cycle    ! metal eta or triple
         fij=param%tors(mol%at(ii))*param%tors(mol%at(jj))             ! atom contribution, central bond
         if(fij.lt.gen%fcthr)                 cycle
         if(param%tors(mol%at(ii)).lt.0.or.param%tors(mol%at(jj)).lt.0) cycle ! no negative values
         if(param%metal(mol%at(ii)).gt.1.and.topo%nb(20,ii).gt.4)  cycle ! no HC metals
         if(param%metal(mol%at(jj)).gt.1.and.topo%nb(20,jj).gt.4)  cycle !
         fqq=1.0d0+abs(topo%qa(ii)*topo%qa(jj))*gen%qfacTOR      ! weaken it for e.g. CF-CF and similar
         call ringsbond(mol%n,ii,jj,cring,sring,rings) ! i and j in same ring
         lring=.false.
         ccij =.false.
         if(rings.gt.0) lring=.true.
         sp3ij=topo%hyb(ii).eq.3.and.topo%hyb(jj).eq.3
         if(mol%at(ii).eq.6.and.mol%at(jj).eq.6) ccij=.true.
         nhi=1
         nhj=1
         do ineig=1,topo%nb(20,ii)
            if(mol%at(topo%nb(ineig,ii)).eq.1) nhi=nhi+1
         enddo
         do jneig=1,topo%nb(20,jj)
            if(mol%at(topo%nb(jneig,jj)).eq.1) nhj=nhj+1
         enddo
         fij=fij*(dble(nhi)*dble(nhj))**0.07 ! n H term
         ! amides and alpha carbons in peptides/proteins
         if (alphaCO(mol%n,mol%at,topo%hyb,topo%nb,piadr,ii,jj)) fij = fij * 1.3d0
         if (amide(mol%n,mol%at,topo%hyb,topo%nb,piadr,ii).and.topo%hyb(jj).eq.3.and.mol%at(jj).eq.6) fij = fij * 1.3d0
         if (amide(mol%n,mol%at,topo%hyb,topo%nb,piadr,jj).and.topo%hyb(ii).eq.3.and.mol%at(ii).eq.6) fij = fij * 1.3d0
         ! hypervalent
         if(btyp(m).eq.4)                                                                   fij = fij * 0.2d0
!        loop over neighbors of ij
         do ineig=1,topo%nb(20,ii)
            kk=topo%nb(ineig,ii)
            if(kk.eq.jj) cycle
            do jneig=1,topo%nb(20,jj)
               ll=topo%nb(jneig,jj)
               if(ll.eq.ii) cycle
               if(ll.eq.kk) cycle
               if(chktors(mol%n,mol%xyz,ii,jj,kk,ll)) cycle  ! near 180
               fkl=param%tors2(mol%at(kk))*param%tors2(mol%at(ll))       ! outer kl term
               if(mol%at(kk).eq.7.and.piadr(kk).eq.0) fkl=param%tors2(mol%at(kk))*param%tors2(mol%at(ll))*0.5
               if(mol%at(ll).eq.7.and.piadr(ll).eq.0) fkl=param%tors2(mol%at(kk))*param%tors2(mol%at(ll))*0.5
               if(fkl.lt.gen%fcthr)               cycle
               if(param%tors(mol%at(kk)).lt.0.or.param%tors(mol%at(ll)).lt.0) cycle ! no negative values
               f1 = gen%torsf(1)
               f2 = 0.0d0
               fkl=fkl*(dble(topo%nb(20,kk))*dble(topo%nb(20,ll)))**(-0.14)  ! CN term

!-----------------------
! definitions come here
!-----------------------
               if ( lring ) then
                    if(rings.gt.3)then
                    call ringstors(mol%n,ii,jj,kk,ll,cring,sring,rings4) ! smallest ring in which i,j,k,l are
                    else
                    rings4=3 ! the 3-ring is special
                    endif
! RING CASE
                    nrot = 1
                    if(btyp(m).eq.2) nrot = 2 ! max at 90 for pi and symmetric at 0,-180,180
                    phi  = 0  ! cis
                    if ( btyp(m).eq.1 .and. rings4 .gt. 0 )then
                      call ringstorl(mol%n,ii,jj,kk,ll,cring,sring,ringl)  ! largest ring in which i,j,k,l are
                      notpicon=piadr(kk).eq.0.and.piadr(ll).eq.0                      ! do it only for sat. rings
                      if ( rings4 .eq. 3 .and.                      notpicon) then; nrot=1; phi = 0.d0; f1=gen%fr3; endif
                      if ( rings4 .eq. 4 .and. ringl.eq.rings4 .and.notpicon) then; nrot=6; phi =30.d0; f1=gen%fr4; endif
                      if ( rings4 .eq. 5 .and. ringl.eq.rings4 .and.notpicon) then; nrot=6; phi =30.d0; f1=gen%fr5; endif
                      if ( rings4 .eq. 6 .and. ringl.eq.rings4 .and.notpicon) then; nrot=3; phi =60.d0; f1=gen%fr6; endif
                    endif
                    if ( rings4.eq.0 .and. btyp(m).eq.1 .and. &
                                 topo%nb(20,kk).eq.1.and.topo%nb(20,ll).eq.1) then; nrot=6; phi =30.d0; f1=0.30; endif
                    if(btyp(m).eq.2 .and. rings.eq.5 .and. mol%at(ii)*mol%at(jj).eq.42) then
                       if(amide(mol%n,mol%at,topo%hyb,topo%nb,piadr,ii).or.amide(mol%n,mol%at,topo%hyb,topo%nb,piadr,jj)) f1=5.  ! improving CB7
                    endif
               else
! ACYCLIC
                                                      phi  = 180.d0 ! trans
                                                      nrot = 1
                    if(topo%hyb(ii).eq.3.and.topo%hyb(jj).eq.3) nrot = 3 ! Me case
                    if(btyp(m).eq.2)                  nrot = 2 ! max at 90 for pi and symmetric at 0,-180,180
                    if(piadr(ii).gt.0.and.(piadr(jj).eq.0.and.topo%hyb(jj).eq.3))then  ! pi-sp3
                                                         f1=0.5d0
                                         if(mol%at(ii).eq.7) f1=0.2d0 ! important for CB7 conf.
                                                         phi  =180.d0
                                                         nrot = 3
                    endif
                    if(piadr(jj).gt.0.and.(piadr(ii).eq.0.and.topo%hyb(ii).eq.3))then
                                                         f1=0.5d0
                                         if(mol%at(jj).eq.7) f1=0.2d0 ! important for CB7 conf.
                                                         phi  =180.d0
                                                         nrot = 3
                    endif
               endif
! SP3 specials
               if(topo%hyb(ii).eq.3.and.topo%hyb(jj).eq.3) then
! N-N, P-P ...
               if(param%group(mol%at(ii)).eq.5.and.param%group(mol%at(jj)).eq.5) then
                                                             nrot=3
                                                             phi=60.d0
                                                             f1= 3.0d0
               endif
! 5-6
               if((param%group(mol%at(ii)).eq.5.and.param%group(mol%at(jj)).eq.6).or. &
     &            (param%group(mol%at(ii)).eq.6.and.param%group(mol%at(jj)).eq.5))then
                                                             nrot=2
                                                             phi=90.d0
                                                             f1=1.0d0
                            if(mol%at(ii).ge.15.and.mol%at(jj).ge.15)f1=20.0d0
               endif
! O-O, S-S ...
               if(param%group(mol%at(ii)).eq.6.and.param%group(mol%at(jj)).eq.6) then
                                                             nrot=2
                                                             phi=90.d0
                                                             f1=5.0d0
                            if(mol%at(ii).ge.16.and.mol%at(jj).ge.16)f1=25.0d0 ! better for h2s2
               endif
               endif
! pi system
               if(pibo(m).gt.0) then
                                f2=pibo(m)*exp(-2.5d0*(1.24d0-pibo(m))**14)  ! decrease to very small values for P < 0.3
                                                                             ! values of 2.5 instead of 2.4 give larger tangles
                                                                             ! the parameter 1.24 is very sensitive ie 1.25 yield 5 deg more in 1,3cB
                                if(piadr(kk).eq.0.and.mol%at(kk).gt.10) f2=f2*1.3! the pi BO becomes more significant if heavies are attached
                                if(piadr(ll).eq.0.and.mol%at(ll).gt.10) f2=f2*1.3
                                f1 = f1 * 0.55
               endif

               if(topo%hyb(kk).eq.5.or.topo%hyb(ll).eq.5) fkl = fkl * 1.5 ! hypervalent corr.
!--------------------
! end of definitions
!-------------------

! total FC            sigma       pi             charge central outer kl
               fctot = (f1 + 10.d0*gen%torsf(2)*f2) * fqq * fij * fkl

               if(fctot.gt.gen%fcthr) then ! avoid tiny potentials
                  topo%ntors=topo%ntors+1
                  if (topo%ntors.gt.maxtors) then
                     call env%error('internal (torsion setup) error', source)
                     return
                  end if
                  topo%tlist(1,topo%ntors)=ll
                  topo%tlist(2,topo%ntors)=ii
                  topo%tlist(3,topo%ntors)=jj
                  topo%tlist(4,topo%ntors)=kk
                  topo%tlist(5,topo%ntors)=nrot
                  topo%vtors(1,topo%ntors)=phi*pi/180.0d0
                  topo%vtors(2,topo%ntors)=fctot
!                 printout
                  phi=valijklff(mol%n,mol%xyz,ll,ii,jj,kk)
                  if(pr)write(env%unit,'(4i5,2x,i2,5x,i2,4x,3f8.3)') &
     &            ii,jj,kk,ll,topo%tlist(5,topo%ntors),rings,topo%vtors(1,topo%ntors)*180./pi,phi*180./pi,topo%vtors(2,topo%ntors)
               endif

! extra rot=1 torsion potential for sp3-sp3 to get gauche conf energies well
               sp3kl=topo%hyb(kk).eq.3.and.topo%hyb(ll).eq.3
               if(sp3kl.and.sp3ij.and.(.not.lring).and.btyp(m).lt.5) then
                  topo%ntors=topo%ntors+1
                  if (topo%ntors.gt.maxtors) then
                     call env%error('internal (torsion setup) error', source)
                     return
                  end if
                  ff = gen%torsf(6)
                  if(mol%at(ii).eq.7.or.mol%at(jj).eq.7) ff = gen%torsf(7)
                  if(mol%at(ii).eq.8.or.mol%at(jj).eq.8) ff = gen%torsf(8)
                  topo%tlist(1,topo%ntors)=ll
                  topo%tlist(2,topo%ntors)=ii
                  topo%tlist(3,topo%ntors)=jj
                  topo%tlist(4,topo%ntors)=kk
                  topo%tlist(5,topo%ntors)=1
                  topo%vtors(1,topo%ntors)=pi
                  topo%vtors(2,topo%ntors)= ff * fij * fkl *fqq
                  if(pr)write(env%unit,'(4i5,2x,i2,5x,i2,4x,3f8.3)') &
     &            ii,jj,kk,ll,topo%tlist(5,topo%ntors),rings,topo%vtors(1,topo%ntors)*180./pi,phi*180./pi,topo%vtors(2,topo%ntors)
               endif

            enddo ! neighbors ij
         enddo
      enddo ! bond loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! out-of-plane, improper (three-fold coordinated central pi atom i or an N)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(pr) write(env%unit,*) 'out-of-plane atoms          phi0    phi      FC'
      do i=1,mol%n
         if(topo%nb(20,i).ne.3) cycle
         if(piadr(i).eq.0) then
            if(      mol%at(i) .ne.7) cycle
         endif
         topo%ntors=topo%ntors+1
         jj=topo%nb(1,i)
         kk=topo%nb(2,i)
         ll=topo%nb(3,i)
!        sort atoms according to distance to central atom such that the same inversion angle def. always results
         sdum3(1)=rab(lin(i,jj))
         sdum3(2)=rab(lin(i,kk))
         sdum3(3)=rab(lin(i,ll))
         ind3(1)=jj
         ind3(2)=kk
         ind3(3)=ll
         call ssort(3,sdum3,ind3)
         jj=ind3(1)
         kk=ind3(2)
         ll=ind3(3)
         topo%tlist(1,topo%ntors)=i
         topo%tlist(2,topo%ntors)=jj
         topo%tlist(3,topo%ntors)=kk
         topo%tlist(4,topo%ntors)=ll
         if(piadr(i).eq.0.and.mol%at(i).eq.7) then  ! sat N case
           r0=80.0d0
           ff=0.60d0
           topo%tlist(5,topo%ntors)=-1
           topo%vtors(1,topo%ntors)=r0*pi/180. ! double min at +/- phi0
           topo%vtors(2,topo%ntors)=0.0d0
           do m=1,topo%nb(20,i)
              idum=topo%nb(m,i)
              topo%vtors(2,topo%ntors)=topo%vtors(2,topo%ntors)+ff*sqrt(param%repz(mol%at(idum)))  ! NX3 has higher inv barr. than NH3
           enddo
         else
           ncarbo=0
           nf    =0
           do m=1,topo%nb(20,i)
              idum=topo%nb(m,i)
              if(mol%at(idum).eq.8.or.mol%at(idum).eq.16) ncarbo=ncarbo+1
              if(param%group(mol%at(idum)).eq.7           ) nf    =nf    +1
           enddo
           fqq=1.0d0+topo%qa(i)*5.0d0
           topo%tlist(5,topo%ntors)=0         ! phi0=0 case (pi)
           topo%vtors(1,topo%ntors)=0.0d0     !  "      "
           sumppi=pbo(lin(i,jj))+pbo(lin(i,kk))+pbo(lin(i,ll))
           f2=1.0d0-sumppi*gen%torsf(5)
!                         base val  piBO  charge term
           topo%vtors(2,topo%ntors)=gen%torsf(3) * f2 * fqq
!          carbonyl corr.
           if(mol%at(i).eq.5.and.ncarbo.gt.0)             topo%vtors(2,topo%ntors)=topo%vtors(2,topo%ntors)*38.
           if(mol%at(i).eq.6.and.ncarbo.gt.0)             topo%vtors(2,topo%ntors)=topo%vtors(2,topo%ntors)*38.
           if(mol%at(i).eq.6.and.nf.gt.0.and.ncarbo.eq.0) topo%vtors(2,topo%ntors)=topo%vtors(2,topo%ntors)*10.
           if(mol%at(i).eq.7.and.ncarbo.gt.0)             topo%vtors(2,topo%ntors)=topo%vtors(2,topo%ntors)*10./f2 ! no pi dep
         endif
!        printout
         phi=omega(mol%n,mol%xyz,i,jj,kk,ll)
         if(pr)write(env%unit,'(4i5,7x,3f8.3)') i,jj,kk,ll,topo%vtors(1,topo%ntors)*180./pi,phi*180./pi,topo%vtors(2,topo%ntors)
      enddo

      write(env%unit,'(10x,"#tors  :",3x,i0)') topo%ntors
      write(env%unit,'(10x,"#nmol  :",3x,i0)') topo%nfrag

! all done

      topo%maxsystem = 5000
      !if(mol%n.gt.800)then
      block
         use xtb_setparam, only : set, p_modh_old

         if (set%mhset%model == p_modh_old) then
            call fragmentize(mol%n,mol%at,mol%xyz,topo%maxsystem,500,rab,topo%nb, &
               & topo%ispinsyst,topo%nspinsyst,topo%nsystem,env)
         else
            topo%nsystem=1
         endif
      end block


      write(env%unit,'(10x,"#optfrag :",3x,i0)') topo%nfrag
      !write(env%unit,*) '#optfrag :',nsystem

      if(pr)then
      write(env%unit,*)
      write(env%unit,*) 'GFN-FF setup done.'
      write(env%unit,*)
      endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz)

!> @brief charge scaling function
pure elemental function zeta(at,q)
use xtb_mctc_accuracy, only : wp
   implicit none
   integer ,intent(in) :: at
   real(wp),intent(in) :: q

   real(wp)           :: zeta,qmod
   real(wp),parameter :: zeff(86) = (/ &
   &   1,                                                 2,  & ! H-He
   &   3, 4,                               5, 6, 7, 8, 9,10,  & ! Li-Ne
   &  11,12,                              13,14,15,16,17,18,  & ! Na-Ar
   &  19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,  & ! K-Kr
   &   9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,  & ! Rb-Xe
   &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Cs-Lu
   &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26/)  ! Hf-Rn
!! Semiempirical Evaluation of the GlobalHardness of the Atoms of 103
!! Elements of the Periodic Table Using the Most Probable Radii as
!! their Size Descriptors DULAL C. GHOSH, NAZMUL ISLAM 2009 in
!! Wiley InterScience (www.inte"rscience.wiley.com).
!! DOI 10.1002/qua.22202
!! values in the paper multiplied by two because
!! (ii:ii)=(IP-EA)=d^2 E/dN^2 but the hardness
!! definition they use is 1/2d^2 E/dN^2 (in Eh)
   real(wp),parameter :: c(1:86) = (/ &
  &0.47259288_wp,0.92203391_wp,0.17452888_wp,0.25700733_wp,0.33949086_wp,0.42195412_wp, & ! H-C
  &0.50438193_wp,0.58691863_wp,0.66931351_wp,0.75191607_wp,0.17964105_wp,0.22157276_wp, & ! N-Mg
  &0.26348578_wp,0.30539645_wp,0.34734014_wp,0.38924725_wp,0.43115670_wp,0.47308269_wp, & ! Al-Ar
  &0.17105469_wp,0.20276244_wp,0.21007322_wp,0.21739647_wp,0.22471039_wp,0.23201501_wp, & ! Ca-Cr
  &0.23933969_wp,0.24665638_wp,0.25398255_wp,0.26128863_wp,0.26859476_wp,0.27592565_wp, & ! Mn-Zn
  &0.30762999_wp,0.33931580_wp,0.37235985_wp,0.40273549_wp,0.43445776_wp,0.46611708_wp, & ! Ga-Kr
  &0.15585079_wp,0.18649324_wp,0.19356210_wp,0.20063311_wp,0.20770522_wp,0.21477254_wp, & ! Rb-Mo
  &0.22184614_wp,0.22891872_wp,0.23598621_wp,0.24305612_wp,0.25013018_wp,0.25719937_wp, & ! Tc-Cd
  &0.28784780_wp,0.31848673_wp,0.34912431_wp,0.37976593_wp,0.41040808_wp,0.44105777_wp, & ! In-Xe
  &0.05019332_wp,0.06762570_wp,0.08504445_wp,0.10247736_wp,0.11991105_wp,0.13732772_wp, & ! Cs-Nd
  &0.15476297_wp,0.17218265_wp,0.18961288_wp,0.20704760_wp,0.22446752_wp,0.24189645_wp, & ! Pm-Dy
  &0.25932503_wp,0.27676094_wp,0.29418231_wp,0.31159587_wp,0.32902274_wp,0.34592298_wp, & ! Ho-Hf
  &0.36388048_wp,0.38130586_wp,0.39877476_wp,0.41614298_wp,0.43364510_wp,0.45104014_wp, & ! Ta-Pt
  &0.46848986_wp,0.48584550_wp,0.12526730_wp,0.14268677_wp,0.16011615_wp,0.17755889_wp, & ! Au-Po
  &0.19497557_wp,0.21240778_wp/)

   intrinsic :: exp

   qmod = zeff(at) + q
   if (qmod.lt.0._wp) then
      zeta = exp( 3.0d0 )
   else
      zeta = exp( 3.0d0  * ( 1._wp - exp( c(at) * ( 1._wp - zeff(at) / qmod ) ) ) )
   endif

end function zeta

end subroutine gfnff_ini

end module xtb_gfnff_ini
