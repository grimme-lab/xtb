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

!! ========================================================================
!  RMSD biased push/pull path finder (RMSD-BPP)
!  SG 12/18
!! ========================================================================
subroutine bias_path(mol,wfx,calc,egap,et,maxiter,epot,grd,sigma)
   use iso_fortran_env, wp => real64, istdout => output_unit

   use mctc_econv

   use tbdef_molecule
   use tbdef_calculator
   use tbdef_wavefunction

   use setparam
   use fixparam

   use optimizer
   use ls_rmsd

   implicit none

   type(tb_molecule),    intent(inout) :: mol
   type(tb_wavefunction),intent(inout) :: wfx
   type(tb_calculator),  intent(in) :: calc
   integer, intent(in)    :: maxiter
   real(wp),intent(in)    :: epot
   real(wp),intent(in)    :: et
   real(wp),intent(inout) :: egap
   real(wp),intent(inout) :: grd(3,mol%n)
   real(wp),intent(inout) :: sigma(3,3)

   real(wp),allocatable :: xyz0(:,:)
   real(wp),allocatable :: xyz1(:,:)
   real(wp),allocatable :: g0(:,:),g1(:,:) 
   real(wp),allocatable :: xyznew(:,:,:)
   real(wp),allocatable :: enew(:)
   real(wp),allocatable :: dat(:,:,:)         
   integer, allocatable :: at1(:)

   integer maxpath
   parameter (maxpath=2000) ! max # points on one path
   integer i,j,k,nn,npath(100),olev,skip,irun,maxoptiter
   real(wp) e,U(3,3),x(3),y(3),rms,barr,drms,gnew,e0,e1,dum0
   real(wp) xx(10),dum,factor,mem(5,100),distort,gnorm,dum1,ang
   real(wp) ma,mfac,gthr,fnat
   parameter (distort=0.005)  ! set such that symmetric HCN and trans-dichlorethen react
   character(len=80) :: atmp,btmp
   logical fail,ex
   integer :: ilog,idat ! file handle

   if(pathset%nrun.gt.100) call raise('E','too big path_nrun (max 100)',1)
   call metaset%allocate(mol%n,2)
   allocate(enew(maxpath),xyznew(3,mol%n,maxpath),dat(2,maxpath,100))
   allocate(xyz0(3,mol%n),xyz1(3,mol%n),g0(3,mol%n),g1(3,mol%n),at1(mol%n))

   if (.not.allocated(pathset%fname)) then
      call raise('E','No product structure given!',1)
   endif
   write(istdout,*) 'reading reference structures from '//pathset%fname//' ...'
   call rdcoord(pathset%fname,mol%n,xyz1,at1)
   if (any(at1.ne.mol%at)) then
      call raise('E','Atom type missmatch between reactant and product!',1)
   endif

   if (pathset%nat.gt.0) then
      metaset%nat = pathset%nat
      metaset%atoms = pathset%atoms
      fnat=real(metaset%nat,wp)   ! # atoms used in RMSD
   else
      metaset%nat = 0
      fnat= mol%n  ! all atoms used in RMSD
   endif
   metaset%width=pathset%alp   ! copy width to metadyn common block
   metaset%nstruc = 2

   write(istdout,'("initial k push/pull              :",2f9.3)') &
      &        pathset%kpush,pathset%kpull
   write(istdout,'("Gaussian width (1/Bohr)          :",f9.3)')metaset%width
   write(istdout,'("# random initial distor incr +/- :",f9.3)')distort  
   write(istdout,'("# optlevel                       :",i4)')optset%optlev   
   write(istdout,'("# of runs                        :",i4)')pathset%nrun
   write(istdout,'("approx # structures opt on path  :",i4)')pathset%nopt
   write(istdout,'("# of ''an''-optimization steps     :",i4)')pathset%anopt
   factor = 1.0
   xyz0 = mol%xyz
   metaset%xyz(:,:,1) = xyz0
   metaset%xyz(:,:,2) = xyz1

   metaset%width =1.0 ! start value
   write(istdout,*) 
   write(istdout,*) ' k_push/k_pull internally multiplied by Natoms!'
   write(istdout,*) ' user info on effect of push/pull parameters'
   write(istdout,*) ' ang. ed-prod G should be small, Govlp large'
   write(istdout,*)  &
      &         '  width   k_push  k_pull  |Gbias| ang. ed-prod G  Govlp'
   do i=1,9
      ! forward
      metaset%factor(1) = pathset%kpush * fnat * factor
      metaset%factor(2) = pathset%kpull * fnat * factor
      g0 = 0
      call metadynamic(metaset,mol%n,mol%at,xyz0,e0,g0)
      ! backward
      metaset%factor(2) = pathset%kpush * fnat * factor
      metaset%factor(1) = pathset%kpull * fnat * factor
      g1 = 0
      call metadynamic(metaset,mol%n,mol%at,xyz1,e1,g1)

      dum=1.d-9
      dum0=1.d-9
      dum1=1.d-9
      do j=1,mol%n
         do k=1,3
            dum0=dum0+g0(k,j)*g0(k,j)
            dum1=dum1+g1(k,j)*g1(k,j)
            dum=dum  +g0(k,j)*g1(k,j)
         enddo
      enddo
      ang =acos(-dum/(sqrt(dum0*dum1)))*360./3.14159
      dum1=0.5*(sum(abs(g0))+sum(abs(g1)))
      write(istdout,'(3f8.3,f12.6,f8.1,f12.5)')metaset%width, &
         &   metaset%factor(2)/fnat, &
         &   metaset%factor(1)/fnat, &
         &   dum1,ang,-dum
      metaset%width = metaset%width - 0.1
      factor = factor * 1.5
   enddo

   call rmsd(mol%n,xyz0,xyz1,0,U,x,y,rms,.false.,g0)
   write(istdout,'("input reactant/product RMSD (Bohr)  :",f9.3,/)')rms
   if(rms.gt.4)  &
      & write(istdout,*) 'WARNING: for large reactant/product RMSD,', &
   &            ' small $path/alp values (width) should be tested.'
   if(optset%optlev.lt.1) &
      &  write(istdout,*) 'WARNING: $opt/level=tight or vtight recommended!'

!! ------------------------------------------------------------------------
!  loop over runs with increasing Vbias
!! ------------------------------------------------------------------------
   metaset%width = pathset%alp 
   factor = 1
   do irun=1,pathset%nrun
      ! distort slightly to avoid trapping
      xyz0 = metaset%xyz(:,:,1)
      do i=1,mol%n
         do j=1,3
            call random_number(dum)
            if(mod(i,2).eq.0)then
               xyz0(j,i)=xyz0(j,i)+distort*dum 
            else
               xyz0(j,i)=xyz0(j,i)-distort*dum 
            endif
         enddo
      enddo
      ! set for metadyn routine (common)
      metaset%factor(1) = pathset%kpush * factor * fnat 
      metaset%factor(2) = pathset%kpull * factor * fnat 

!! ------------------------------------------------------------------------
      ! make path
!! ------------------------------------------------------------------------
      write(istdout,*)
      write(istdout,*)'optimizing inital path with Vbias ...'
      write(istdout,'("actual k push/pull :",2f8.3)')  &
         &          metaset%factor(1:2)/fnat           
      maxoptiter=0         

      mol%xyz = xyz0
      call geometry_optimization &
         &          (mol,wfx,calc, &
         &           egap,et,maxiter,maxoptiter,e,grd,sigma,optset%optlev, &
         &           .true.,.false.,fail)
      write(atmp,'("mv xtbopt.log xtbpath_biasopt_",i0,".xyz")') &
         &      irun 
      call execute_command_line(atmp)                    

!! ------------------------------------------------------------------------
      ! read opt history file (kept)
!! ------------------------------------------------------------------------
      write(atmp,'("xtbpath_biasopt_",i0,".xyz")') irun 
      write(btmp,'("xtbpath_",i0,".xyz")') irun 
      open(newunit=ilog,file=atmp)
      j=1 
 10   continue
      read(ilog,'(a)',end=100) atmp
      read(ilog,'(a)')         atmp
      call readl(atmp,xx,nn)
      gnew=xx(2) ! gradient 
      do i=1,mol%n
         read(ilog,'(a)') atmp
         call readl(atmp,xx,nn)
         xyznew(1:3,i,j)=xx(1:3)*aatoau
      enddo
      if(gnew.gt.1.d-8) j=j+1 ! discard ancopt restart points
      goto 10
100   continue
      if(j-1.gt.maxpath) call raise('E','internal error in path.f',1)
      close(ilog,status='delete')
      npath(irun) = j-1 

      skip = max(1,npath(irun) / pathset%nopt) ! this leads approx to the requested # points

!! ------------------------------------------------------------------------
      ! opt all without Vbias
!! ------------------------------------------------------------------------
      metaset%nstruc = 0   ! no constraint
      olev=p_olev_crude    ! just steps
      maxoptiter=pathset%anopt! max steps for opt
      write(istdout,*) ' # of points on xtbopt.log: ',npath(irun)
      write(istdout,*) ' doing short optimizations without Vbias ... '
      k=0
      do i=1,npath(irun)
         if(mod(i,skip).ne.0) cycle
         k=k+1
         mol%xyz=xyznew(:,:,i)
         call geometry_optimization &
            &       (mol,wfx,calc, &
            &        egap,et,maxiter,maxoptiter,enew(k),grd,sigma,olev, &
            &        .false.,.true.,fail)
         if(k.gt.1.and.mod(k,10).eq.0)  &
            &      write(istdout,'(i5,"  unconstrained dE :",f14.6)')k, &
            &      enew(k)-enew(1)
         xyznew(:,:,k)=mol%xyz
      enddo
      npath(irun)=k
      metaset%nstruc = 2 ! reactivate

!! ------------------------------------------------------------------------
      ! check RMSD if something happened
!! ------------------------------------------------------------------------
      call rmsd(mol%n,xyznew(:,:,1),xyznew(:,:,k),0,U,x,y,rms,.false.,g0)
      mem(3,irun)=rms

!! ------------------------------------------------------------------------
      ! save results
!! ------------------------------------------------------------------------
      open(newunit=ilog,file=btmp)
      barr=-1.d+42
      rms =0
      do i=1,npath(irun) 
         if(i.gt.1)then
            call rmsd(mol%n,xyznew(:,:,i-1),xyznew(:,:,i),0,U,x,y,drms, &
                      .false.,g0)
         endif
         dum=autokcal*(enew(i)-enew(1))
         if(dum.gt.barr) barr=dum
         call wrlog2(ilog,mol%n,xyznew(1,1,i),mol%at,dum)
         dat(1,i,irun) = rms 
         dat(2,i,irun) = dum 
         rms=rms+drms ! the plot coordinate in dat file is the sum of dRMSD with previous structure
      enddo
      close(ilog)
      write(istdout,*) 'final path written to ',trim(btmp)  
      mem(1,irun)=barr
      mem(2,irun)=autokcal*(enew(npath(irun))-enew(1))
      call rmsd(mol%n,metaset%xyz(:,:,1),xyznew(:,:,1),0,U,x,y,rms,.false.,g0)
      mem(4,irun)=rms
      call rmsd(mol%n,metaset%xyz(:,:,2),xyznew(:,:,npath(irun)),0,U,x,y,rms, &
                .false.,g0)
      mem(5,irun)=rms
      ! increase power (but keep width)
      factor = factor * 1.5

   enddo
!! ------------------------------------------------------------------------
   ! done
!! ------------------------------------------------------------------------

!! ------------------------------------------------------------------------
   ! final output
!! ------------------------------------------------------------------------
   open(unit=idat,file='xtbpath.dat') ! for comparison all pot curves
   write(istdout,*) 
   write(istdout,*) 'energies in kcal/mol, RMSD in Bohr'
   do irun=1,pathset%nrun
      if(mem(3,irun).lt.0.5.or.abs(mem(1,irun)).lt.1) &
         & write(istdout,*) 'WARNING: may be no reaction so increase |kpush/pull|' &
         &           ,' or decrease kpath_alp'
      write(istdout,'("run",i2,"  barrier:", &
         &       f7.2,"  dE:",f7.2, &
         &   "  ed-prod RMSD:",f6.2, &
         &   "  ed-startpath RMSD:",f6.2, &
         &   "  prod-endpath RMSD:",f6.2 &
         &       )') irun,mem(1:5,irun)
      drms=dat(1,npath(irun),irun)/mem(3,irun)    ! normalize path length
      do j=1,npath(irun)
         write(idat,*) dat(1,j,irun)/drms,dat(2,j,irun)
      enddo
      write(idat,*)
   enddo
   close(idat)
   write(istdout,*)  &
      'check potential xtbpath.dat (cummu. dRMSD in Bohr vs. E in kcal)'

   call execute_command_line('rm -rf xtbopt.log xtbopt.coord') ! just crab from last short opt 

end subroutine bias_path


