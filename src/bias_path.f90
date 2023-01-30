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
module xtb_biaspath
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: bias_path

contains

!! ========================================================================
!  RMSD biased push/pull path finder (RMSD-BPP)
!  SG 12/18, modified 06/20
!! ========================================================================
subroutine bias_path(env, mol, chk, calc, egap, et, maxiter, epot, grd, sigma)
   use xtb_mctc_filetypes, only : fileType
   use xtb_mctc_convert
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_calculator
   use xtb_type_restart
   use xtb_type_data, only : scc_results
   use xtb_single, only : singlepoint
   use xtb_io_writer, only : writeMolecule
   use xtb_setparam
   use xtb_fixparam
   use xtb_cqpath
   use xtb_geoopt
   use xtb_optimizer
   use xtb_lsrmsd
   use xtb_metadynamic
   use xtb_readin
   use xtb_axis, only : axis
   use xtb_param_atomicRad, only : atomicRad

   implicit none

   character(len=*), parameter :: source = 'xtb_path'

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   type(TMolecule),    intent(inout) :: mol
   type(TRestart),intent(inout) :: chk
   class(TCalculator), intent(inout) :: calc
   integer, intent(in)    :: maxiter
   real(wp),intent(in)    :: epot
   real(wp),intent(in)    :: et
   real(wp),intent(inout) :: egap
   real(wp),intent(inout) :: grd(3,mol%n)
   real(wp),intent(inout) :: sigma(3,3)

   real(wp),allocatable :: xyze(:,:)
   real(wp),allocatable :: xyzp(:,:)
   real(wp),allocatable :: gtmp(:,:)
   real(wp),allocatable :: xyzpath(:,:,:)
   real(wp),allocatable :: xyzdum (:,:,:)
   real(wp),allocatable :: bo(:,:), cn(:)
   real(wp),allocatable :: pvec1(:), pvec2(:)
   real(wp),allocatable :: epath(:), mempath(:,:)
   integer, allocatable :: atp(:)
   integer, allocatable :: splitl(:)
   integer, allocatable :: npath(:)

   integer i,j,k,ii,nn,olev,skip,irun,maxoptiter,nk,maxcycle,idum,maxpath
   integer np,iter,its,modef,run,npathk,tsmodeav,ialp,nrefine,nalp,imem(100)
   integer ilog,ilog2,idat ! file handle

   real(wp) U(3,3),x(3),y(3),rms,barr,drms,gnew,e1,dum0,r,rco,rms2,factor2
   real(wp) x3(3),t3(3),xx(10),dum,factor,distort,gnorm,dum1,thr,scc,w1,w2
   real(wp) fnat,rmsd_prod_thr,alp_change,k_change,edrmsd,ppull
   real(wp) gdum,gts,emax,eback,eforw,eforw_prev,gsum_prev,pthr
   real(wp) de_min,grad_min,e_ed,e_pr,de,ddot,gthr

   character(len=80) :: atmp,btmp
   logical fail,ex,done,degen
   type(scc_results) :: res

   !! ------------------------------------------------------------------------
   !  some parameters
   !! ------------------------------------------------------------------------

   nrefine =set%pathset%nrun! # of refinement (full path opt) steps (usually 1)
   nalp    =4           ! # of test runs with different alpha bias
   nk      =6           ! # of test runs with different k bias
   maxpath =nalp*nk

   ppull= set%pathset%ppull ! default pull strength on path point

   pthr = 0.00_wp       ! include points on path depending on overlap of approx TS modes
   gthr = 0.01_wp       ! Gnorm at TS exit thr
   rmsd_prod_thr=0.2_wp ! RMSD between product and last point on path found
   ! if less, loop terminates and path is considered as ok
   k_change  = 1.5_wp   ! increase k push/pull in each try by this
   alp_change= 0.2_wp   ! decrease alpha in each alp try by this
   scc       = 2.0_wp   ! high accuracy not important in SP

   set%optset%maxdispl_opt=0.50_wp ! don't make it too small (more points, erratic behavior)
   set%optset%hlow_opt    =0.03_wp ! too small values can cause spikes

   ilog =14
   ilog2=142

   call metaset%allocate(mol%n,2)
   allocate(splitl(mol%n),gtmp(3,mol%n),pvec1(3*mol%n),pvec2(3*mol%n))
   allocate(xyze(3,mol%n),xyzp(3,mol%n),bo(mol%n,mol%n),atp(mol%n),cn(mol%n))
   allocate(npath(maxpath),mempath(4,maxpath))

   if (.not.allocated(set%pathset%fname)) then
      call env%error('No product structure given!', source)
      return
   endif
   write(env%unit,*) 'reading reference structures from '//set%pathset%fname//' ...'
   call rdcoord(set%pathset%fname,mol%n,xyzp,atp) ! product
   if (any(atp.ne.mol%at)) then
      call env%error('Atom type missmatch between reactant and product!', source)
   endif
   xyze=mol%xyz ! reactant

   call rmsd(mol%n,xyze,xyzp,0,U,x,y,rms,.false.,gtmp)
   write(env%unit,'("reactant product RMSD :",f9.3)')rms

   if (set%pathset%nat.gt.0) then
      metaset%nat = set%pathset%nat
      metaset%atoms = set%pathset%atoms
      fnat=real(metaset%nat,wp)   ! # atoms used in RMSD
   else
      metaset%nat = 0
      fnat= mol%n  ! all atoms used in RMSD
   endif

   write(env%unit,'("initial k push/pull (in code xNat) :",2f9.3)')set%pathset%kpush,set%pathset%kpull
   write(env%unit,'("initial Gaussian width (1/Bohr)    :",f9.3)')set%pathset%alp
   write(env%unit,'("# refinement runs                  :",i4)')nrefine
   write(env%unit,'("# of ''an''-optimization steps       :",i4)')set%pathset%anopt
   write(env%unit,'("# optlevel                         :",i4)')set%optset%optlev
   write(env%unit,*)

   !! ------------------------------------------------------------------------
   !  check if its poly-molecular reaction (in case increas rmsd_prod_thr)
   !! ------------------------------------------------------------------------

   bo = 0
   do i=1,mol%n
      cn(i)=0
      do j=1,mol%n
         r=sqrt(sum(xyzp(:,i)-xyzp(:,j)))
         rco=(atomicRad(mol%at(i))+atomicRad(mol%at(j)))*autoaa
         if(r.lt.2.5*rco) then
            bo(j,i)=1
            cn(i)=cn(i)+1
         endif
      enddo
      bo(i,i)=0
   enddo
   call mrec(i,xyzp,cn,bo,mol%n,mol%at,splitl)
   rmsd_prod_thr = rmsd_prod_thr * float(i) ! RMSD can be large if fragments are involved
   deallocate(bo,cn,splitl)

   call axis(mol%n,mol%at,xyze,x(1),x(2),x(3))
   call axis(mol%n,mol%at,xyzp,y(1),y(2),y(3))

   degen = .false.
   dum = 2.0_wp*abs(sum(x)-sum(y))/(sum(x)+sum(y))
   if(dum.lt.0.01) degen = .true.  ! possibly degenrate E/P
   write(env%unit,'("degenerate system :",l,2f9.6)')degen,sum(x),sum(y)

   !! ------------------------------------------------------------------------
   !  loop over runs with increasing Vbias
   !! ------------------------------------------------------------------------

   metaset%xyz(:,:,1) = xyze
   metaset%xyz(:,:,2) = xyzp
   done  = .false.
   irun  = 0
   factor2= 0.0_wp
   call delete_file('.xtbtmpmode')
   call delete_file('hessian')

   bias_loop: do ialp=1,nalp ! alp
      factor = 1.0_wp
      do run=1,nk  ! k_push/pull

         irun = irun + 1
         ! set for metadyn routine (common)
         metaset%factor(1) = set%pathset%kpush * factor * fnat
         metaset%factor(2) = set%pathset%kpull * factor * fnat
         metaset%width(:) = set%pathset%alp - factor2
         metaset%width(:) = max(metaset%width(:),0.2_wp)

         !! ------------------------------------------------------------------------
         ! make path
         !! ------------------------------------------------------------------------

         if(degen)then ! avoid symmetry trapping
            do i=1,mol%n
               do j=1,3
                  call random_number(dum)
                  if(mod(i,2).eq.0)then
                     mol%xyz(j,i)=mol%xyz(j,i)+0.05*dum
                  else
                     mol%xyz(j,i)=mol%xyz(j,i)-0.05*dum
                  endif
               enddo
            enddo
         endif

         metaset%nstruc = 2
         maxoptiter     = 0
         olev           = set%optset%optlev       ! opt level
         mol%xyz        = xyze

         atmp='xtbopt.log'
         call open_file(ilog,atmp,'w')
         call ancopt (env,ilog,mol,chk,calc,egap,et,maxiter,maxoptiter,dum,grd,sigma,olev,.false.,fail)
         call close_file(ilog)

         !! ------------------------------------------------------------------------
         ! read opt history file (kept)
         !! ------------------------------------------------------------------------

         call cqpath_read_pathfile_parameter(atmp,j,k,npath(irun))
         if(k.ne.mol%n) call env%terminate('read error in bias_path, 1')
         allocate(xyzpath(3,mol%n,npath(irun)),epath(npath(irun)))
         call cqpath_read_pathfile(atmp,j,k,npath(irun),xyzpath,mol%at,epath)
         xyzpath=xyzpath/0.52917726_wp

         ! check RMSD if something happened
         mol%xyz(:, :) = xyzpath(:, :, npath(irun))
         call rmsd(mol%n,mol%xyz,xyzp,0,U,x,y,rms,.false.,gtmp)

         write(env%unit,'(i3," # points, run ",i3," for k push/pull/alpha :",3f8.3,5x," prod-ed RMSD:",f8.3)')  &
            &  npath(irun),run,metaset%factor(1:2)/fnat, metaset%width(1), rms

         !! ------------------------------------------------------------------------
         ! save results
         !! ------------------------------------------------------------------------

         write(btmp,'("xtbpath_",i0,".xyz")') irun
         call open_file(ilog2,btmp,'w')
         metaset%nstruc = 0 ! no bias
         barr=-1.d+42
         do i=1,npath(irun)
            mol%xyz(:, :) = xyzpath(:, :, i)
            call singlepoint(env,mol,chk,calc,egap,et,maxiter,0,.true.,.false.,scc,epath(i),grd,sigma,res)
            dum=autokcal*(epath(i)-epath(1))
            if(dum.gt.barr) barr=dum
            call writeMolecule(mol, ilog2, fileType%xyz, energy=dum)
         enddo
         call close_file(ilog2)
         mempath(1,irun)=barr
         mempath(2,irun)=autokcal*(epath(npath(irun))-epath(1))
         mempath(3,irun)=rms
         mempath(4,irun)=metaset%width(1)
         deallocate(xyzpath,epath)

         ! increase power
         factor = factor * k_change

         if(done) exit bias_loop

         if(rms.lt.rmsd_prod_thr) then ! found product
            factor = factor / k_change / 1.2_wp  ! test a bit smaller push/pull
            done = .true. ! close enough to product, add one more (with less bias) and then exit
         endif

      enddo
      factor2= factor2 + alp_change
   enddo bias_loop
   
   !! ------------------------------------------------------------------------
   ! output and find path yielding product
   !! ------------------------------------------------------------------------

   write(env%unit,*)
   write(env%unit,*) 'path trials (see xtbpath_*.xyz), energies in kcal/mol'
   j = 0
   do i=1,irun
      write(env%unit,'("run",i2,"  barrier:",f7.2,"  dE:",f7.2, "  product-end path RMSD:",f8.3 &
         &       )') i,mempath(1:3,i)
      if(mempath(3,i).lt.rmsd_prod_thr) then
         j = j + 1
         imem(j) = i
      endif
   enddo
   if(j.eq.0) then
      call env%error('No product! Decrease alp in $path and/or optlevel.', source)
      return
   end if
   dum = 1.d+42
   k = 1
   do i=1,j   ! take from those the one with lowest barrier
      ii=imem(i)
      if(mempath(1,ii).lt.dum) then
         dum = mempath(1,ii)
         k = ii
      endif
   enddo
   npathk = k

   ! important modification of the pull strength on the path
   ppull = ppull * (mempath(4,npathk)/1.2_wp)**2  ! if the opt. path alpha is small, the
   ! kpull for the constrained opt. must be small
   ! e.g. for helicene or BCF racem.

   !! ------------------------------------------------------------------------
   !  take best path for refinement
   !! ------------------------------------------------------------------------

   write(btmp,'("xtbpath_",i0,".xyz")') npathk
   call cqpath_read_pathfile_parameter(btmp,j,k,np)
   allocate(xyzpath(3,mol%n,np),epath(np))
   call cqpath_read_pathfile(btmp,j,k,np,xyzpath,mol%at,epath)
   xyzpath=xyzpath/0.52917726_wp
   write(env%unit,'("path ",i2," taken with ",i4," points.")') npathk, np

   deallocate(npath,mempath)

   !! ------------------------------------------------------------------------
   !  make the path nicer for optimizations
   !! ------------------------------------------------------------------------

   write(env%unit,'("screening points ...")')
   if(np.gt.set%pathset%nopt)then
      allocate(xyzdum(3,mol%n,np))
      call screenpath(np,set%pathset%nopt,idum,mol%n,xyzpath,epath,xyzdum)
      xyzdum(:,:,1:idum)=xyzpath(1:,:,1:idum)
      deallocate(xyzpath,epath)
      np=idum
      allocate(xyzpath(3,mol%n,np),epath(np))
      xyzpath(:,:,1:np)=xyzdum(1:,:,1:np)
      deallocate(xyzdum)
      write(env%unit,'("new # points :",i3)') np
   endif

   !  call cleanpath(np,mol%n,xyzpath,xyze,xyzp,pthr) ! not good

   epath = 0
   btmp='xtbpath_0.xyz' ! write the path before opt.
   write(env%unit,'("start path on file ",a)') btmp
   call open_file(ilog2,btmp,'w')
   metaset%nstruc = 0 ! unbiased energy and gradient
   mol%xyz=xyzpath(:,:,1)
   call singlepoint(env,mol,chk,calc,egap,et,maxiter,0,.true.,.false.,scc,epath(1),grd,sigma,res)
   call writeMolecule(mol, ilog2, fileType%xyz, energy=0.0_wp)
   do i=2,np-1
      mol%xyz=xyzpath(:,:,i)
      call singlepoint(env,mol,chk,calc,egap,et,maxiter,0,.true.,.false.,scc,epath(i),grd,sigma,res)
      call writeMolecule(mol, ilog2, fileType%xyz, energy=autokcal*(epath(i)-epath(1)))
   end do
   mol%xyz=xyzpath(:,:,np)
   call singlepoint(env,mol,chk,calc,egap,et,maxiter,0,.true.,.false.,scc,epath(np),grd,sigma,res)
   call writeMolecule(mol, ilog2, fileType%xyz, energy=autokcal*(epath(np)-epath(1)))
   call close_file(ilog2)

   !! ------------------------------------------------------------------------
   !  new section for path refinment by restricted (path biased) opt.
   !! ------------------------------------------------------------------------

   maxoptiter=set%pathset%anopt

   allocate(xyzdum(3,mol%n,np))
   call open_file(ilog,'xtbopt.log','w')
   eforw= 1.0d+42
   gts  = 1.0d+42
   do iter=1,nrefine
      write(env%unit,'(''refinement cycle'',i4)')iter
      call alignpath(mol%n,np,xyzpath)

      xyzdum(:,:,1:np) = xyzpath(:,:,1:np)      ! save coords for path mode projection
      do i=2,np-1                               ! one opt for all points
         if(mod(i,10).eq.0.or.i.eq.2) write(*,*) 'optimizing points ',i,' ...'
         mol%xyz=xyzpath(:,:,i)
         call metadyn_tsmode(mol%n,np,i,xyzdum,mol%xyz,ppull) ! add bias on path (and constrain optionally)
         call ancopt (env,ilog,mol,chk,calc,egap,et,maxiter,maxoptiter,epath(i),grd,sigma,olev,.false.,fail)
         xyzpath(:,:,i)=mol%xyz
      end do

      call open_file(ilog2,'xtbpath.xyz','w')
      mol%xyz=xyzpath(:,:,1) ! first point
      call writeMolecule(mol, ilog2, fileType%xyz, energy=0.0_wp)
      write(env%unit,*)
      dum =-1.d+42
      metaset%nstruc = 0 ! unbiased energy and gradient
      do i=2,np-1
         mol%xyz=xyzpath(:,:,i)
         call singlepoint(env,mol,chk,calc,egap,et,maxiter,0,.true.,.true.,scc,epath(i),grd,sigma,res)
         call writeMolecule(mol, ilog2, fileType%xyz, energy=autokcal*(epath(i)-epath(1)))
         if(epath(i).gt.dum)then
            dum=epath(i)
            gts=sqrt(sum(grd**2))
            its=i
         endif
      enddo
      mol%xyz=xyzpath(:,:,np)  ! last point
      call writeMolecule(mol, ilog2, fileType%xyz, energy=autokcal*(epath(np)-epath(1)))
      call close_file(ilog2)

      gsum_prev=gts
      eforw_prev=eforw
      emax=maxval(epath(1:np))
      eforw=(emax-epath(1) )*autokcal
      eback=(emax-epath(np))*autokcal
      de   =(epath(np)-epath(1))*autokcal
      write(env%unit,'(''forward  barrier (kcal)  :'',f10.3)') eforw
      write(env%unit,'(''backward barrier (kcal)  :'',f10.3)') eback
      write(env%unit,'(''reaction energy  (kcal)  :'',f10.3)') de
      write(env%unit,'(''opt. pull strength       :'',f10.3)') ppull
      write(env%unit,'(''norm(g) at est. TS, point:'', f8.5,i4)') gts, its

      if(gts.lt.gthr) then
         write(env%unit,'("terminated because grad at TS <threshold")')
         exit
      endif

   end do
   write(env%unit,*)
   if(iter.ge.nrefine)write(env%unit,'("terminated because max. # cycles reached")')

   !! ------------------------------------------------------------------------
   !  output
   !! ------------------------------------------------------------------------

   btmp='xtbpath_ts.xyz'
   write(env%unit,'("estimated TS on file ",a)') btmp
   call open_file(ilog,btmp,'w')
   mol%xyz=xyzpath(:,:,its)
   call writeMolecule(mol, ilog, fileType%xyz, energy=epath(its))
   close(ilog)

   dum=0
   write(env%unit,'("path data (pmode=approx. path mode):")')
   write(env%unit,'("point",5x,"drms",5x,"energy",1x,"pmode ovlp",1x,"pmode grad")')
   do i=2,np
      call guess_tsmode(mol%n,np,i,xyzpath,pvec1,pvec2)
      gdum=ddot(3*mol%n,pvec1,1,pvec2,1)
      mol%xyz=xyzpath(:,:,i)
      call rmsd(mol%n,xyzpath(:,:,i),xyzpath(:,:,i-1),0,U,x,y,rms,.false.,gtmp)
      dum1=(epath(i)-epath(i-1))/(rms*real(mol%n))  ! gradient along path
      write(env%unit,'(i4,3f10.3,f10.5)') i,dum,(epath(i)-epath(1))*autokcal,gdum,dum1
      dum=dum+rms
   enddo

end subroutine bias_path

!! ========================================================================
!  remove unecessary points
!! ========================================================================
subroutine screenpath(np,npwanted,npnew,n,xyz,e,xyzdum)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment
   implicit none
   integer   n,np,npwanted,npnew
   real(wp)  xyz   (3,n,np)       ! path coords
   real(wp)  e     (np)           ! energies
   real(wp)  xyzdum(3,n,np)       ! path coords

   integer i, j, k
   real(wp) :: dum1,dum2,thr
   real(wp), allocatable :: grad(:)

   allocate(grad(np))

   xyzdum = xyz

   do i=2,np-1
      dum2 = sum((xyz(:,:,i+1)-xyz(:,:,i-1))**2)+1.0_wp
      dum1 = e(i+1)-e(i-1)
      grad(i) = dum1/dum2
   enddo

   thr=0.01_wp
   k  = 0
10 continue
   k = k + 1
   j = 1
   do i=2,np-1
      if(abs(grad(i)).gt.thr)then
         j = j + 1
         xyz(:,:,j) = xyzdum(:,:,i)
      endif
   enddo
   if(j.ge.npwanted.and.k.lt.100)then
      thr = thr * 1.3
      goto 10
   endif

   npnew = j + 1
   xyz(:,:,1) = xyzdum(:,:,1)
   xyz(:,:,npnew) = xyzdum(:,:,np)

end subroutine


!! ========================================================================
!  put a bias on structure of actual point on path (and project out mode)
!! ========================================================================

subroutine metadyn_tsmode(n,sn,its,xyzpath,xyzact,kpull)
   use xtb_mctc_accuracy, only : wp
   use xtb_metadynamic
   use xtb_fixparam
   implicit none
   integer   n,its,sn
   real(wp)  xyzpath(3,n,sn)    ! path coords
   real(wp)  xyzact (3,n)       ! actual "
   real(wp)  kpull

   real(wp), allocatable :: pvec1(:), pvec2(:)
   integer :: i,j,k
   integer :: modef, idat


   metaset%nstruc = 1

   metaset%factor(1) = -kpull * real(n)  ! pull point on path, 0.05 good choice
   metaset%width(1) = 1.3_wp

   metaset%xyz(:,:,1) = xyzact(:,:)

   ! do not project, is unstable

   !   allocate(pvec1(3*n),pvec2(3*n))
   !   call guess_tsmode(n,sn,its,xyzpath,pvec1,pvec2)
   !   modef = 1
   !   idat  =15
   !   call open_binary(idat,'.xtbtmpmode','w')
   !   write (idat) modef
   !   write (idat) 0.5*(pvec1+pvec2)
   !   call close_file(idat)
   !   deallocate(pvec1,pvec2)

end subroutine

!! ========================================================================
!  get TS mode(s) (two directions) at path point its
!! ========================================================================

subroutine guess_tsmode(n,sn,its,xyz,pvec1,pvec2)
   use xtb_mctc_accuracy, only : wp
   implicit none
   integer   n,its,sn
   real(wp)  xyz(3,n,sn)
   real(wp)  pvec1(3*n)
   real(wp)  pvec2(3*n)

   integer :: i,j,k
   real(wp)  rms1,rms2,dum,dum2,third1,third2
   real(wp), allocatable, dimension(:) :: v
   parameter (third1 = 1.0_wp/3.0_wp)
   parameter (third2 = 2.0_wp/3.0_wp)

   i   = 0
   rms1= 0
   rms2= 0
   do k=1,n
      do j=1,3
         i = i + 1
         dum = (xyz(j,k,its+1) - xyz(j,k,its))
         if(its+2.le.sn)then
            dum2 = (xyz(j,k,its+2) - xyz(j,k,its+1))
            dum  = third2*dum+third1*dum2
         endif
         pvec1(i) = dum
         rms1= rms1 + dum**2
         dum = (xyz(j,k,its  ) - xyz(j,k,its-1))
         if(its.gt.2)then
            dum2 = (xyz(j,k,its-1) - xyz(j,k,its-2))
            dum  = third2*dum+third1*dum2
         endif
         pvec2(i) = dum
         rms2= rms2 + dum**2
      enddo
   enddo
   pvec1(:) = pvec1(:) / dsqrt(rms1) ! the vectors point in the same direction
   pvec2(:) = pvec2(:) / dsqrt(rms2)

end subroutine

subroutine alignpath(nat,sn,xyz)
   implicit none
   integer nat,sn
   real(wp) xyz(3,nat,sn)

   real(wp), allocatable, dimension(:,:) :: sp_atom
   real(wp), allocatable, dimension(:) :: x,y,x2,y2,strecke_X
   real(wp) :: strecke_X_inc,temp_strecke_X,xyz_0(3),xyz_1(3),xyz_diff(3)
   integer :: i,j,k

   allocate (strecke_X(sn),sp_atom(3,sn),x(sn),y(sn),x2(sn),y2(sn))

   ! loop over all atoms to get an array of all points of each atom along the path
   do i=1,nat
      sp_atom(:,:) = xyz(:,i,:)
      ! loop over path
      xyz_0(:) = sp_atom(:,1)
      strecke_X(1) = 0.0
      do j=2,sn
         xyz_1(:) = sp_atom(:,j)
         xyz_diff = xyz_1 - xyz_0
         xyz_0 = xyz_1
         strecke_X(j) = strecke_X(j-1) + sqrt(sum(xyz_diff**2))
      end do
      ! make a cubic spline for x,y,z dependent on strecke_X
      !       X = strecke_X , Y = [x|y|z]
      do k=1,3
         x(:) = strecke_X(:)
         y(:) = sp_atom(k,:)
         ! reposition on x
         call cqpath_repos_sym(sn,x,x2)               ! symmetric
         call cqpath_damped_cube_spline(sn,x,y,x2,y2) ! much better than non-damped
         sp_atom(k,:) = y2(:)
      end do
      ! save point
      xyz(:,i,:) = sp_atom(:,:)
   end do

   deallocate (strecke_X,sp_atom,x,y,x2,y2)
end subroutine

subroutine cqpath_repos_sym(n,x,x2)
   ! in    n  - dimension of array x,y
   !       x  - array x
   ! out   x2 - symmetric points from array x
   !
   implicit none
   integer, intent(in)  :: n
   real(wp), intent(in)  :: x(n)
   real(wp), intent(out) :: x2(n)

   real(wp) x_wert, x_anfang, x_ende, x_schritt
   integer :: i

   ! symmetric distribution
   ! dx(i) = distance/points
   x_anfang = x(1)
   x_ende = x(n)
   x_schritt = (x_ende-x_anfang)/(REAL(n,wp)-1.0)

   x_wert = x_anfang
   x2(1) = x_anfang
   do i=2,n-1
      x_wert = x_wert + x_schritt
      x2(i) = x_wert
   end do
   x2(n) = x_ende

end subroutine

subroutine cqpath_repos_asym(n,x1,x2)
   implicit none
   integer, intent(in)  :: n
   real(wp), intent(in)  :: x1(n)
   real(wp), intent(out) :: x2(n)

   real(wp) x_schritt, x_anfang, x_ende, x, ynorm
   real(wp), dimension(:), allocatable :: y
   integer :: i

   allocate(y(n))

   ynorm = 0.
   x = -1.
   do i=1,n
      y(i)=1./exp(-2.0*x**4)
      ynorm = ynorm + y(i)
      x = x + 2./float(n-1)
   enddo
   y = real(n,8) * y / ynorm

   x_anfang = x1(1)
   x_ende   = x1(n)
   x_schritt= (x_ende-x_anfang)/(REAL(n,wp)-2.0)

   x = x_anfang
   do i=1,n
      x2(i) = x
      x = x + x_schritt*y(i)
   end do

end subroutine

subroutine cqpath_cube_spline(n,x,y,x2,y2)
   ! in    n  - dimension of array x,y
   !       x  - array x
   !       y  - array y
   ! out   x2 - symmetric points from x
   !       y2 - symmetric points from y
   implicit none
   integer, intent(in) :: n
   real(wp), intent(in) :: x(n),y(n),x2(n)
   real(wp), intent(out) :: y2(n)

   real(wp), dimension(:), allocatable :: ypp
   real(wp) y_wert, y_wert_ypval, y_wert_yppval
   integer :: i

   allocate (ypp(n))
   call spline_cubic_set ( n, x, y, 0, 0.0d0, 0, 0.0d0, ypp )
   y2(1) = y(1)
   do i=2,n-1
      call spline_cubic_val ( n, x, y, ypp, x2(i), y_wert, y_wert_ypval, y_wert_yppval )
      y2(i) = y_wert
   end do
   y2(n) = y(n)
   deallocate (ypp)
   return
end subroutine

subroutine cqpath_damped_cube_spline(n,x,y,x2,y2)
   ! in    n  - dimension of array x,y
   !       x  - array x
   !       y  - array y
   ! out   x2 - symmetric points from x
   !       y2 - symmetric points from y
   implicit none
   integer, intent(in) :: n
   real(wp), intent(in) :: x(n),y(n),x2(n)
   real(wp), intent(out) :: y2(n)

   real(wp), dimension(:,:), allocatable :: abcd
   real(wp) y_wert
   integer :: i

   allocate (abcd(n,4))
   call cqpath_damped_cube_spline_set(n,x,y,abcd)

   y2(1) = y(1)
   do i=2,n-1
      call cqpath_damped_cube_spline_val(n,abcd,x,x2(i),y_wert)
      y2(i) = y_wert
   end do
   y2(n) = y(n)

   deallocate (abcd)
   return
end subroutine

subroutine cqpath_damped_cube_spline_set(n,x,y,abcd)
   ! in    n  - dimension von array x,y
   !       x  - array x
   !       y  - array y
   ! out   abcd - array of polynom coefficients
   !
   ! Based on Constrained Cubic Spline Interpolation for Chemical Engineering Application
   ! by CJC Kruger
   !
   ! Rewritten and optimized in fortran
   implicit none
   integer, intent(in) :: n
   real(wp), intent(in) :: x(n)
   real(wp), intent(in) :: y(n)
   real(wp), intent(out) :: abcd(n,4)

   real(wp) fs_x0,fs_x1,fss_x0,fss_x1
   real(wp) x0,x1,x2,y0,y1,y2
   real(wp) x1_x0,y1_y0,x0_2,x0_3,x1_x0_2,x1_2
   real(wp) steigung
   real(wp), dimension(:), allocatable :: dx,dy,f1

   integer :: i,j,nf

   nf = n-1
   allocate (dx(nf),dy(nf),f1(n))

   do i=1,nf
      dx(i) = x(i+1)-x(i)
      dy(i) = y(i+1)-y(i)
   end do

   do i=2,nf
      steigung = dy(i-1)*dy(i)
      if (steigung > 0.0) then
         f1(i) = 2.0/(dx(i)/dy(i)+dx(i-1)/dy(i-1))
      else if (steigung <= 0.0) then
         f1(i) = 0.0
      end if
   end do

   f1(1) = 3.0*dy(1)/(2.0*dx(1)) - f1(2)/2.0
   f1(n) = 3.0*dy(n-1)/(2.0*dx(n-1)) - f1(n-1)/2.0

   do i=2,n
      x1_x0 = dx(i-1)
      y1_y0 = dy(i-1)
      x0_2 = x(i-1)*x(i-1)
      x0_3 = x0_2*x(i-1)
      x1_x0_2 = x1_x0*x1_x0
      x1_2 = x(i)*x(i)

      fss_x0 = -2.0*(f1(i) + 2.0*f1(i-1))/(x1_x0)+6.0*(y1_y0)/x1_x0_2
      fss_x1 = 2.0*(2.0*f1(i) + f1(i-1))/(x1_x0)-6.0*(y1_y0)/x1_x0_2

      abcd(i,4) = (fss_x1 - fss_x0)/(6.0*(x1_x0))
      abcd(i,3) = (x(i)*fss_x0-x(i-1)*fss_x1)/(2.0*x1_x0)
      abcd(i,2) = ((y1_y0)-abcd(i,3)*(x1_2-x0_2)-abcd(i,4)*(x1_2*x(i)-x0_3))/(x1_x0)
      abcd(i,1) = y(i-1) - abcd(i,2)*x(i-1) - abcd(i,3) * x0_2 - abcd(i,4) * x0_3
   end do

   deallocate(dx,dy,f1)
end subroutine

subroutine cqpath_damped_cube_spline_val(n,abcd,x,x1,y1)
   ! in    n    - dimension von array x,y
   !       abcd - array of polynom coefficients
   !       x    - x
   ! out   y    - f(x,abcd)
   !
   ! Based on Constrained Cubic Spline Interpolation for Chemical Engineering Application
   ! by CJC Kruger
   !
   ! Rewritten and optimized in fortran
   integer, intent(in) :: n
   real(wp), intent(in) :: abcd(n,4)
   real(wp), intent(in) :: x(n)
   real(wp), intent(in) :: x1
   real(wp), intent(out) :: y1

   integer :: j
   real(wp) :: xx

   do j=2,n
      if (x1 < x(j)) then
         if (x1 >= x(j-1)) then
            xx = x1*x1
            y1 = abcd(j,1)+abcd(j,2)*x1+abcd(j,3)*xx+abcd(j,4)*xx*x1
         end if
      end if
   end do
   return
end subroutine

!! ========================================================================
!  remove bad points, not used
!! ========================================================================

subroutine cleanpath(np,n,xyz,xyz0,xyz1,pthr)
   use xtb_type_environment
   implicit none
   integer   n,np
   real(wp)  xyz(3,n,np)       ! path coords
   real(wp)  xyz0(3,n)         ! path coords
   real(wp)  xyz1(3,n)         ! path coords
   real(wp)  pthr              ! threshold

   real(wp), allocatable :: pvec1(:), pvec2(:), ovlp(:), delta(:,:), xyzdum(:,:,:)
   real(wp)              :: dum, ddot
   integer               :: i,j,k,l
   integer               :: ib, ie

   write(*,'(''input path length :'',i3)') np

   allocate(pvec1(3*n),pvec2(3*n),ovlp(np),delta(3,n),xyzdum(3,n,np))

   ovlp = 1
   do i=2,np-1
      call guess_tsmode(n,np,i,xyz,pvec1,pvec2)
      ovlp(i)=ddot(3*n,pvec1,1,pvec2,1)
   enddo
   write(*,*) 'path overlaps before alignment'
   write(*,'(20f6.2)') ovlp(1:np)
   call alignpath(n,np,xyz)
   do i=2,np-1
      call guess_tsmode(n,np,i,xyz,pvec1,pvec2)
      ovlp(i)=ddot(3*n,pvec1,1,pvec2,1)
   enddo
   write(*,*) 'path overlaps after alignment'
   write(*,'(20f6.2)') ovlp(1:np)

   ! get bad points and replace them by linear interpolation
   do i=2,np-1
      if(ovlp(i).lt.pthr)then
         write(*,*) 'bad point ',i
         do k=i-1,1,-1
            ib=k
            if(ovlp(k).gt.pthr) then
               ib = k
               exit
            endif
         enddo
         do k=i+1,np
            ie=k
            if(ovlp(k).gt.pthr) then
               ie = k
               exit
            endif
         enddo
         !         write(*,*) i,ib,ie,ovlp(i)
         ! structures ib+1...i-1 are unreliable
         delta (:,:) = (xyz(:,:,ie) - xyz(:,:,ib))/float(ib-ie)
         l = 0
         do k=ib+1,ie-1
            l = l + 1
            xyzdum(:,:,k) = xyz(:,:,ib) - delta(:,:)*float(l)
         enddo
      endif
   enddo

   do i=2,np-1
      if(ovlp(i).lt.pthr) xyz(:,:,i) = xyzdum(:,:,i)
   enddo

   call alignpath(n,np,xyz)
   do i=2,np-1
      call guess_tsmode(n,np,i,xyz,pvec1,pvec2)
      ovlp(i)=ddot(3*n,pvec1,1,pvec2,1)
   enddo
   write(*,*) 'path overlaps after cleaning and alignment'
   write(*,'(20f6.2)') ovlp(1:np)

end subroutine

end module xtb_biaspath
