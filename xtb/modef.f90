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

module modef
   use iso_fortran_env, wp => real64
   use single, only : singlepoint
contains
   ! mode following routine for conformational searches or 1D-anharmonic calc
   ! thres types of operation:
   ! mode_local = 0 : conformational PES scan based on normal coords
   ! mode_local = 1 : conformational PES scan based on localized normal coords
   ! mode_local =-1 : PES scan for anharmonic corrections (ie no minima opt.)

subroutine modefollow(mol,wfn,calc,egap,et,maxiter,epot,grd,sigma)
   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding, only : c_null_char

   use mctc_econv, only : autokcal, aatoau, autorcm, amutoau

   use tbdef_molecule
   use tbdef_calculator
   use tbdef_wavefunction
   use tbdef_data

   use setparam
   use splitparam

   implicit none
   intrinsic date_and_time
   type(tb_molecule), intent(inout) :: mol
   type(tb_wavefunction),intent(inout) :: wfn
   type(tb_calculator),intent(in) :: calc
   integer :: icall,maxiter
   real(wp) :: epot,et,egap
   real(wp), intent(inout) :: grd(3,mol%n)
   real(wp), intent(inout) :: sigma(3,3)
   type(scc_results) :: res
   interface
      function rename(oldpath, newpath) bind(c) result(stat)
         use iso_c_binding
         character(len=*,kind=c_char), intent(in) :: oldpath, newpath
         integer(c_int) :: stat
      end function rename
   end interface

   real(wp),allocatable:: freq(:), u(:,:), uu(:,:)
   real(wp),allocatable:: xyza (:,:,:), rmass(:), e(:), e2(:), geo(:,:)
   real(wp),allocatable:: coord0(:,:), displ(:), xyzmin(:,:), bmat(:,:)
   real(wp),allocatable:: x(:),b(:),c(:),d(:),xx(:),yy(:),yy1(:),yy2(:)

   real(wp) :: step,step2,x1,x2,dx,tmpx,yyy,yyy1,yyy2,emin,efreq,eps,dum
   real(wp) :: norm,fc,updscal,aa1,bb1,cc1,aa2,bb2,cc2,ee,thr,dip(3)
   integer :: ii,i,j,k,n3,maxoptiter,ia,ic,ij,ss,nstep,mode,np,kk
   integer :: nsplinep,nstp,istart,iend,nporg,nprj
   integer :: m,ifind,isave,imod,kkk,minpos,tag,n36
   !integer :: map(1000),stplist(100),stptyp(100),na(n),nb(n),nc(n)
   integer,allocatable :: map(:),stplist(:),stptyp(:)
   integer,allocatable :: na(:),nb(:),nc(:)
   character(len=2) :: asym,adum
   character(len=30) :: ctmp,fname
   logical :: ex,fail,vecupdate,scfonly,ldum,intermol
   integer :: ich,val(8)

   write(*,*)
   write(*,'(7x,''======================================='')')
   write(*,'(7x,''|     M O D E   F O L L O W I N G     |'')')
   write(*,'(7x,''======================================='')')
   write(*,*)

   !     generate a tag to identify in confg the scan
   !call system("date | gawk '{print $4}' | sed s/':'/''/g > .tmpxxx")
   !open(newunit=ich,file='.tmpxxx')
   !read(ich,*)tag
   !close(ich,status='delete')
   !ppracht 09/2018 begin
   allocate(map(1000),stplist(100),stptyp(100),na(mol%n),nb(mol%n),nc(mol%n))

   call date_and_time(values=val)
   tag=365*val(1) + 30.5*val(2) + val(3) ! actual date
   !ppracht 09/2018 end

   call delete_file('xtbmodefok')

   mode=mode_follow
   tag=tag+100000*mode
   intermol=iatf2.gt.0.and.mode.lt.7

   imod=mode
   fname='xtb_normalmodes'
   if(mode_local.gt.0) then    ! file changed to local modes
      fname='xtb_localmodes'
      imod=mode+1000
   endif

   nprj=1
   if(mode_prj.gt.0) nprj=2

   write(*,*) 'intermolecular mode : ',intermol
   if(mode_follow.ne.0)then
      write(*,*) 'normal mode file    : ',trim(fname)
      write(*,*) '# of projected modes: ',nprj
   endif

   ! error section
   inquire(file=fname,exist=ex)
   if(.not.ex.and.mode_follow.ne.0) &
      & call raise('E','Hessian not found! run -hess first!',1)

   if(mode.lt.6.and.(.not.intermol)) call raise('E','mode <7(6) makes no sense!',1)
   if(mode_follow.gt.3*mol%n) call raise('E','mode >3N makes no sense!',1)

   !call system('mv hessian hessian.xtbtmp 2>/dev/null') ! don't use hessian in ancopt
   i = rename('hessian'//c_null_char,'hessian.xtbtmp'//c_null_char)              ! don't use hessian in ancopt

   n3=3*mol%n
   allocate(freq(n3),u(n3,n3),uu(n3,nprj),xyzmin(3,mol%n), &
   &         geo(3,mol%n),coord0(3,mol%n),displ(n3),xyza(3,mol%n,1000),rmass(n3))


   ! mode_follow=0: do rotations within pockets
   if(mode_follow.eq.0) then
      call pocketscan(mol%n,mol%at,mol%xyz,nstp)
      xyza=0.0d0
      ! read and re-opt
      call open_file(ich,'xtb_modescan_0.xyz','r')
      rewind(ich)
      do k=1,nstp
         read(ich,*)
         read(ich,*)
         do i=1,mol%n
            read(ich,*)adum,xyza(1:3,i,k)
         enddo
      enddo
      call close_file(ich)
      imod=0
      k=0
      do ii=1,nstp
         k=k+1
         mol%xyz(1:3,1:mol%n)=xyza(1:3,1:mol%n,ii)/aatoau
         call geometry_optimization &
         &          (mol,wfn,calc, &
         &           egap,et,maxiter,maxoptiter,ee,grd,sigma,optset%optlev, &
         &           .false.,.true.,fail)
         kk=k
         call modeminname(imod,kk,ctmp)
         inquire(file=ctmp,exist=ex)
         if(ex)then
            444            kk=kk+1
            call modeminname(imod,kk,ctmp)
            inquire(file=ctmp,exist=ex)
            if(ex) goto 444
         endif
         write(*,'(''Erel (w.r.t input min / kcal :'' f8.3)') &
         &       (ee-epot)*autokcal
         call open_file(ich,ctmp,'w')
         write(ich,*)mol%n
         write(ich,'(F16.8,1x,i9)')ee,tag
         do j=1,mol%n
            write(ich,'(1x,a2,1x,3f18.8)') &
            &          asym(mol%at(j)),mol%xyz(1:3,j)/aatoau
         enddo
         call close_file(ich)
      enddo
      goto 999
   endif


   ! read normal modes
   ! normal case
   call open_binary(ich,fname,'r')
   read (ich)k
   if(k.ne.n3) call raise('E','severe read error on modes',1)
   read (ich)freq
   read (ich)rmass
   read (ich)u
   call close_file(ich)
   ! remove masses from normal modes and re-normalize
   do m=1,n3
      norm=0
      do ia=1,mol%n
         do ic=1,3
            ii = (ia-1)*3+ic
            u(ii,m)=u(ii,m)*sqrt(atmass(ia))  ! this is correct... now it is
            norm=norm+u(ii,m)**2
         enddo
      enddo
      u(1:n3,m)=u(1:n3,m)/(sqrt(norm)+1.d-12)
   enddo

   if(intermol.and.mode_follow.lt.7) then
      write(*,*) 'normal mode file    : xtb_inormalmodes'
      call open_binary(ich,'xtb_inormalmodes','r')
      !        read intermolecular modes and put them on 1-6
      read(ich)freq(1:6)
      read(ich)rmass(1:6)
      read(ich)u(1:3*mol%n,1:6)
      call close_file(ich)
   endif

   if(nprj.eq.2)then
      uu(1:3*mol%n,2)=u(1:3*mol%n,mode_prj) ! assume that the additional mode which should be projected out is nr. mode_prj
      if(mode.eq.mode_prj) call raise('E','bad choice of options',1) ! ie, normally the TS one
   endif

   maxoptiter=optset%maxoptcycle
   isave=optset%micro_opt
   optset%micro_opt=1000 ! no coordinate update in ancopt

   if(mod(mode_nscan,2).eq.0) mode_nscan=mode_nscan+1 ! make it odd
   nstep=(mode_nscan-1)/2

   fc=0.5*(rmass(mode)*amutoau)*(freq(mode)/autorcm)**2  ! force constant
   step=mode_step*sqrt(2*0.04/fc)/nstep               ! to get 0.04 Eh energy i.e. consistent paths
   step=min(step,3.0d0)                               ! for stiff and floppy systems
   step=max(step,0.08d0)
   if(mode_local.gt.0) step=mode_step                 ! freq. has no meaning for loc modes

   updscal=mode_updat

   scfonly=.false.
   if(mol%n.eq.2.or.maxoptiter.lt.1) scfonly=.true.
   write(*,'(''# steps              :'',i5   )')mode_nscan
   write(*,'(''along normal mode    :'',i5   )')mode
   write(*,'(''frequency /cm-1      :'',F12.3)')freq(mode)
   write(*,'(''force constant       :'',F10.5)')fc
   write(*,'(''mode vector mixing   :'',F10.5)')updscal
   write(*,'(''step length          :'',F10.5)')step
   write(*,'(''maxiter (opt)        :'',i5   )')maxoptiter
   write(*,'(''optlevel for minima  :'',i5   )')optset%optlev


   np=2*nstep+1    ! total # of points on path
   allocate(e(np))

   !     check if input is ok
   call singlepoint(istdout,mol,wfn,calc, &
   &         egap,et,maxiter,0,.true.,.true.,1.0d0,e(nstep+1),grd,sigma,res)
   dum=sqrt(sum(grd**2))
   write(*,'(''RMS gradient         :'',F10.5)')dum
   write(*,*)
   !if(dum.gt.0.005.and.nprj.eq.1.and.freq(mode).gt.0) then
   !  write(*,*)'WARNING WARNING WARNING WARNING WARNING WARNING'
   !  write(*,*)'mode following on incompletely optimized geometry!'
   !  write(*,*)'makes no sense and hence exiting.'
   !  goto 999
   !endif
   e(1:np)=1000. ! if points are skipped
   efreq=max(10.0d0,50*freq(mode)*2.8591D-3) ! in kcal
   efreq=min(efreq,50.0d0)                   ! in kcal
   !     the minimum (which should be the input struc)
   e(nstep+1)=epot
   xyza(1:3,1:mol%n,nstep+1)=mol%xyz(1:3,1:mol%n)

   !     n36=3*n-6
   !     allocate(bmat(n36,n3))
   !     call bzmat(n,at,xyz,bmat)
   !     call xyzint(xyz,n,na,nb,nc,1.0d0,geo)
   !     call gmetry(n,geo,mol%xyz,na,nb,nc)
   !     xyza(1:3,1:n,nstep+1)=mol%xyz(1:3,1:n)

   ! -direction
   uu(1:n3,1)=u(1:n3,mode) ! start org mode
   call open_binary(ich,'.xtbtmpmode','w')
   write(ich) nprj
   write(ich) uu
   call close_file(ich)
   xyzmin = mol%xyz
   do ss=1,nstep
      do ia=1,mol%n
         do ic=1,3
            ij = (ia-1)*3+ic
            mol%xyz(ic,ia)=mol%xyz(ic,ia)-step*uu(ij,1)
         enddo
      enddo
      !        call intmodestep(n,bmat,uu,-step,geo,na,nb,nc,mol%xyz)
      coord0 = mol%xyz
      if(scfonly)then
         call singlepoint(istdout,mol,wfn,calc, &
         &         egap,et,maxiter,0,.true.,.false.,1.0d0,e(nstep-ss+1),grd,sigma,res)
      else
         call geometry_optimization &
         &       (mol,wfn,calc, &
         &        egap,et,maxiter,maxoptiter,e(nstep-ss+1),grd,sigma,optset%optlev, &
         &        .false.,.true.,fail)
      endif
      xyza(1:3,1:mol%n,nstep-ss+1)=mol%xyz(1:3,1:mol%n)
      write(*,'(''step length :'',f6.2,'' dE /kcal :'',F8.3)') &
      &   -step*ss,autokcal*(e(nstep-ss+1)-e(nstep+1))
      !        if(ss.eq.1.and.e(nstep-ss+1)-e(nstep+1).gt.0.1)then
      if(e(nstep-ss+1)-e(nstep+1).gt.0.1)then
         write(*,*) 'exit because energy is too high',nstep-ss+1
         goto 888
      endif
      do ia=1,mol%n
         do ic=1,3
            displ((ia-1)*3+ic)=(mol%xyz(ic,ia) - coord0(ic,ia))*updscal
         enddo
      enddo
      if(updscal.gt.1.d-6)call writeuu(.true.,n3,displ,uu,nprj)  ! new search vector
   enddo

   ! +direction
   888  write(*,*)
   !     call xyzint(xyz,n,na,nb,nc,1.0d0,geo)
   mol%xyz = xyzmin
   uu(1:n3,1)=u(1:n3,mode) ! start org mode
   call open_binary(ich,'.xtbtmpmode','w')  ! important: switch on projection mode in opt
   write(ich) nprj
   write(ich) uu
   call close_file(ich)
   do ss=1,nstep
      do ia=1,mol%n
         do ic=1,3
            ij = (ia-1)*3+ic
            mol%xyz(ic,ia)=mol%xyz(ic,ia)+step*uu(ij,1)
         enddo
      enddo
      !        call intmodestep(n,bmat,uu,step,geo,na,nb,nc,mol%xyz)
      coord0 = mol%xyz
      if(scfonly)then
         call singlepoint(istdout,mol,wfn,calc, &
         &         egap,et,maxiter,0,.true.,.false.,1.0d0,e(nstep+ss+1),grd,sigma,res)
      else
         call geometry_optimization &
         &       (mol,wfn,calc, &
         &        egap,et,maxiter,maxoptiter,e(nstep+ss+1),grd,sigma,optset%optlev, &
         &        .false.,.true.,fail)
      endif
      xyza(:,:,nstep+ss+1)=mol%xyz
      write(*,'(''step length :'',f6.2,'' dE /kcal :'',F8.3)') &
      &   step*ss,autokcal*(e(nstep+ss+1)-e(nstep+1))
      if(e(nstep+ss+1)-e(nstep+1).gt.0.1)then
         write(*,*) 'exit because energy is too high',nstep+ss+1
         goto 889
      endif
      do ia=1,mol%n
         do ic=1,3
            displ((ia-1)*3+ic)=(mol%xyz(ic,ia) - coord0(ic,ia))*updscal
         enddo
      enddo
      displ=-displ
      if(updscal.gt.1.d-6)call writeuu(.true.,n3,displ,uu,nprj)  ! new search vector
   enddo

   889  write(*,*) 'writing path to file xtb_modescan_*'
   write(ctmp,'(''xtb_modescan_'',i0,''.xyz'')') imod
   call open_file(ich,ctmp,'w')

   emin=minval(e)
   e=e-emin    ! set min to zero
   !     if(mode_local.gt.0) efreq=50.
   !     remove very high-lying points
   istart=1
   iend  =np
   do i=1,np
      if(autokcal*e(i).lt.efreq)then
         istart=i
         exit
      endif
   enddo
   do i=np,1,-1
      if(autokcal*e(i).lt.efreq)then
         iend=i
         exit
      endif
   enddo
   nporg=np
   np=iend-istart+1
   if(np.le.1) goto 999 ! curve is bogous
   allocate(e2(np),b(np),c(np),d(np),x(np))
   e2(1:np)=e(istart:iend)
   write(*,*) 'removing high energy points ',nporg-np
   !     write the xyz file
   step2=-nstep*step
   kk=0
   do k=1,nporg
      if(k.ge.istart.and.k.le.iend) then
         kk=kk+1
         map(kk)=k
         x(kk)=step2
         write(ich,*)mol%n
         write(ich,'(F16.8,1x,3F9.3)')e(k),freq(mode),rmass(mode),step2
         write(143,*)step2,627.51*e(k)
         do i=1,mol%n
            write(ich,'(1x,a2,1x,3f18.8)') &
            &      asym(mol%at(i)),xyza(1:3,i,k)/aatoau
         enddo
      endif
      step2=step2+step
   enddo
   call close_file(ich)

   ! spline the potential
   !                  x= x input, e2= y input, b,c,d are the spline coeffcients
   !                  dimensions are e2(np),b(np),c(np),d(np),x(np)
   !                  np= # of input points
   call spline2(x, e2, b, c, d, np)

   nsplinep=50000 ! # of points to be generated (interpolated)
   allocate(xx(nsplinep), yy(nsplinep), yy1(nsplinep), yy2(nsplinep))

   x1=x(1) +0.05*x(1)
   x2=x(np)+0.05*x(np)
   dx=(x2-x1)/(nsplinep-1)
   tmpx=x1
   do i=1,nsplinep
      xx(i)=tmpx
      !                      x, original points, spline coefficients, # input points,
      call ispline(tmpx, x, e2, b, c, d, np, yyy,yyy1,yyy2)
      yy (i)=yyy  ! splined value at tmpx
      yy1(i)=yyy1 ! first derivative at tmpx
      yy2(i)=yyy2 ! second derivative at tmpx
      tmpx=tmpx+dx
   enddo
   !     find minima and saddle points
   call curveanal(nsplinep,xx,yy,yy1,yy2,stplist,stptyp,nstp,eps)

   !call system('rm .xtbtmpmode') ! switch off projection mode in opt
   call delete_file('.xtbtmpmode')
   optset%micro_opt =isave   ! back to default
   maxoptiter=0       ! because we need accurate structures which allow comparison

   if(nstp.gt.1.and.mode_local.ge.0)then  ! good, there are new minima and its not a anharm run
      !        at which point on the curve is the min?
      k=0
      thr=0.1 ! different from input ie RK off the min at 0 ?
      if(mode.eq.7) thr=-99 ! include in a single case the input min.
      ! so that the ensemble contains it
      do i=1,nstp
         if(stptyp(i).ge.2)then  ! min?
            dx=xx(stplist(i))
            if(abs(dx).gt.thr)then
               ii=ifind(dx,np,x)
               k=k+1
               write(*,'(''minimum '',i2,'' at point '',i2,  &
               &                 '' re-opt:'')') k,ii
               mol%xyz = xyza(:,:,map(ii))
               call geometry_optimization &
               &          (mol,wfn,calc, &
               &           egap,et,maxiter,maxoptiter,ee,grd,sigma,optset%optlev, &
               &           .false.,.true.,fail)
               kk=k
               call modeminname(imod,kk,ctmp)
               inquire(file=ctmp,exist=ex)
               if(ex)then
   44             kk=kk+1
                  call modeminname(imod,kk,ctmp)
                  inquire(file=ctmp,exist=ex)
                  if(ex) goto 44
               endif
               write(*,'(''Erel (w.r.t input min / kcal :'' f8.3)')  &
               &       (ee-epot)*autokcal
               call open_file(ich,ctmp,'w')
               write(ich,*)mol%n
               write(ich,'(F16.8,1x,i9)')ee,tag
               do j=1,mol%n
                  write(ich,'(1x,a2,1x,3f18.8)') &
                  &          asym(mol%at(j)),mol%xyz(1:3,j)/aatoau
               enddo
               call close_file(ich)
            endif
         endif
      enddo
   endif
   ! local searches often generate no new minima but lets try anyway a few opts
   ! if no minimum is found, the endpoints sometimes lead to one in a full opt
   ldum=nstp.le.2.and.mode_local.ge.0.and.mode.gt.6 ! exclude anharm case and intermol one
   if(mode_local.gt.0.or.ldum) then
      minpos=nstep+1
      do k=1,mode_nscan,5
         if(abs(k-minpos).lt.5) cycle
         write(*,'(''re-opt at point '',i2)') k
         mol%xyz=xyza(:,:,k)
         call geometry_optimization &
         &       (mol,wfn,calc, &
         &        egap,et,maxiter,maxoptiter,ee,grd,sigma,optset%optlev, &
         &        .false.,.true.,fail)
         kk=k
         call modeminname(imod,kk,ctmp)
         inquire(file=ctmp,exist=ex)
         if(ex)then
            45           kk=kk+1
            call modeminname(imod,kk,ctmp)
            inquire(file=ctmp,exist=ex)
            if(ex) goto 45
         endif
         write(*,'(''Erel (w.r.t input min / kcal :'' f8.3)')  &
         &           (ee-epot)*autokcal
         call open_file(ich,ctmp,'w')
         write(ich,*)mol%n
         write(ich,'(F16.8,1x,i9)')ee,tag
         do j=1,mol%n
            write(ich,'(1x,a2,1x,3f18.8)') &
            &         asym(mol%at(j)),mol%xyz(1:3,j)/aatoau
         enddo
         call close_file(ich)
      enddo
   endif

   999  continue
   !call system('mv hessian.xtbtmp hessian 2>/dev/null') ! restore
   i = rename('hessian.xtbtmp'//c_null_char,'hessian'//c_null_char)              ! restore
   !call system('touch xtbmodefok')                     ! for parallel runs
   call touch_file('xtbmodefok')                    ! for parallel runs

   !ppracht 09/2018 begin ! saw: some are not always allocated
   if(allocated(yy2))    deallocate(yy2)
   if(allocated(yy1))    deallocate(yy1)
   if(allocated(xx))     deallocate(xx)
   if(allocated(yy))     deallocate(yy)
   if(allocated(x))      deallocate(x)
   if(allocated(d))      deallocate(d)
   if(allocated(c))      deallocate(c)
   if(allocated(b))      deallocate(b)
   if(allocated(e))      deallocate(e)
   if(allocated(e2))     deallocate(e2)
   if(allocated(rmass))  deallocate(rmass)
   if(allocated(xyza))   deallocate(xyza)
   if(allocated(displ))  deallocate(displ)
   if(allocated(coord0)) deallocate(coord0)
   if(allocated(geo))    deallocate(geo)
   if(allocated(xyzmin)) deallocate(xyzmin)
   if(allocated(uu))     deallocate(uu)
   if(allocated(u))      deallocate(u)
   if(allocated(freq))   deallocate(freq)
   if(allocated(nc))     deallocate(nc)
   if(allocated(nb))     deallocate(nb)
   if(allocated(na))     deallocate(na)
   if(allocated(stptyp)) deallocate(stptyp)
   if(allocated(stplist))deallocate(stplist)
   if(allocated(map))    deallocate(map)
   !ppracht 09/2018 end

end subroutine modefollow

subroutine modeminname(mode,mini,nam)
   implicit none
   integer mode,mini
   character(len=*) :: nam
   write(nam,'(''xtb_modemin_'',i0,''_'',i0,''.xyz'')')mini,mode
end subroutine modeminname

subroutine writeuu(wr,n3,d,uu,nprj)
   implicit none
   logical wr
   integer n3,i,nprj
   real(wp) d(n3),uu(n3,nprj),norm
   integer ich

   uu(1:n3,1) = uu(1:n3,1) - d(1:n3)
   norm=sum(uu(1:n3,1)**2)
   uu(1:n3,1) = uu(1:n3,1) / sqrt(norm)
   if(wr)then
      call open_binary(ich,'.xtbtmpmode','w')
      write(ich) nprj
      write(ich) uu
      call close_file(ich)
   endif

end subroutine writeuu

end module modef

