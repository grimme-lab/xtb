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

module xtb_hessian
   use xtb_mctc_accuracy, only : wp

contains

subroutine numhess( &
      & env,mol,chk0,calc, &
      & egap,et,maxiter,etot,gr,sr,res)
   use xtb_mctc_accuracy, only : wp
!$ use omp_lib

   use xtb_mctc_convert

!! ========================================================================
!  type definitions
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_calculator
   use xtb_type_data

   use xtb_setparam
   use xtb_splitparam
   use xtb_fixparam

   use xtb_single, only : singlepoint
   use xtb_axis, only : axis

   use xtb_gfnff_calculator, only : TGFFCalculator

   implicit none

   !> Source of errors in the main program unit
   character(len=*), parameter :: source = "hessian_numhess"

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(inout) :: mol
   integer, intent(in)    :: maxiter
   type(TRestart),intent(inout) :: chk0
   class(TCalculator), intent(inout) :: calc
   real(wp) :: eel
   real(wp),intent(inout) :: etot
   real(wp),intent(in)    :: et
   real(wp),intent(inout) :: egap
   real(wp),intent(inout) :: gr(3,mol%n)
   real(wp),intent(inout) :: sr(3,3)
   type(freq_results),intent(out) :: res

   type(TRestart) :: chk
   type(scc_results) :: sccr,sccl
   real(wp) :: rij(3),step,zpve,t1,t0,dumi,dum,xsum
   real(wp) :: dumj,acc,w0,w1,step2,aa,bb,cc,scalh,hof,h298
   real(wp) :: sum1,sum2,trdip(3),dipole(3)
   real(wp) :: trpol(3),sl(3,3)
   integer  :: n3,i,j,k,ic,jc,ia,ja,ii,jj,info,lwork,a,b,ri,rj
   integer  :: nread,kend,lowmode
   integer  :: nonfrozh,izero(6)
   integer, allocatable :: nb(:,:)
   integer, allocatable :: indx(:),molvec(:)
   real(wp),allocatable :: bond(:,:)

!$ integer  :: nproc
   logical :: parallize

   real(wp),allocatable :: h (:,:)
   real(wp),allocatable :: hss(:)
   real(wp),allocatable :: aux (:)
   real(wp),allocatable :: isqm(:)
   real(wp),allocatable :: gl  (:,:)
   real(wp),allocatable :: xyzsave(:,:)
   real(wp),allocatable :: pold(:)
   real(wp),allocatable :: dipd(:,:)

   type(TMolecule) :: tmol

   logical :: ex,rd
   integer :: ich ! file handle
   integer :: err
   character(len=:),allocatable :: hname,fname
   logical :: exitRun
   character(len=128) :: errStr

   parameter (scalh =1.00d0)

   n3=3*mol%n
   call res%allocate(mol%n)
   res%n3true = n3-3*freezeset%n

   allocate(hss(n3*(n3+1)/2),h(n3,n3), &
      & gl(3,mol%n),isqm(n3),xyzsave(3,mol%n),dipd(3,n3), &
      & pold(n3),nb(20,mol%n),indx(mol%n),molvec(mol%n),bond(mol%n,mol%n))

   rd=.false.
   xyzsave = mol%xyz

   step=0.0001_wp
   call rotmol(mol%n,mol%xyz,step,2.*step,3.*step)

   ! step length
   step=step_hess
   if(extcode.eq.5) step=step*2.0_wp ! MOPAC is not very accurate
   ! SCC accuraccy
   acc=accu_hess

   call singlepoint &
      & (env,mol,chk0,calc, &
      &  egap,et,maxiter,0,.true.,.true.,acc,res%etot,res%grad,sr,sccr)

   write(env%unit,'(''step length          :'',F10.5)') step
   write(env%unit,'(''SCC accuracy         :'',F10.5)') acc
   write(env%unit,'(''Hessian scale factor :'',F10.5)') scalh
   write(env%unit,'(''frozen atoms in %    :'',F10.5,i5)') &
      & real(freezeset%n,wp)/real(mol%n,wp)*100,freezeset%n

   res%gnorm = norm2(res%grad)
   write(env%unit,'(''RMS gradient         :'',F10.5)', advance='no') res%gnorm
   if(res%gnorm.gt.0.002_wp) then
      write(env%unit, '(1x,"!! INCOMPLETELY OPTIMIZED GEOMETRY !!")')
      call env%warning('Hessian on incompletely optimized geometry!', source)
      call env%check(exitRun)
      if (exitRun) return
   else
      write(env%unit, '(a)')
   end if

   res%linear=.false.
   call axis(mol%n,mol%at,mol%xyz,aa,bb,cc)
   if(cc.lt.1.d-10) res%linear=.true.
   call timing(t0,w0)
   step2=0.5d0/step

   h = 0.0_wp

!! ========================================================================
!  Hessian part -----------------------------------------------------------

   parallize = .true.
   select type(calc)
   type is (TGFFCalculator)
      calc%update = .false.
      parallize = .false.
   end select

   omp_parallize: if (parallize) then

   if(freezeset%n.gt.0) then
      ! for frozfc of about 10 the frozen modes
      ! approach 5000 cm-1, i.e., come too close to
      ! the real ones
      nonfrozh=mol%n-freezeset%n
      do a = 1,mol%n
         res%freq(a)=float(a)
      enddo
      do a=1,freezeset%n
         res%freq(freezeset%atoms(a))=freezeset%atoms(a)*100000d0
      enddo
      call sortind(mol%n,res%freq)
      do a=1,nonfrozh
         indx(a)=idint(res%freq(a))
      enddo
      do a=nonfrozh+1,mol%n
         indx(a)=idint(res%freq(a)/100000d0)
      enddo
      write(*,'(''atoms frozen in Hessian calc.:'',10i4)') &
         & indx(nonfrozh+1:mol%n)
      ! the index array indx(1:n) contains from 1:nonfrozh in ascending order
      ! the non-frozen atoms and from nonfrozh+1:n the frozen ones also in
      ! ascending order
      ! now compute a subblock of the Hessian
      !$ nproc = omp_get_num_threads()
      !$omp parallel default(shared) &
      !$omp firstprivate(et,maxiter,acc) &
      !$omp private(ia,ic,ii,ja,jc,jj,eel,gr,gl,egap,sccr,sccl,sr,sl,chk,tmol) &
      !$omp shared (mol,h,dipd,pold,step,step2,t1,t0,w1,w0,indx,nonfrozh,env,calc,chk0)
      !$ call omp_set_num_threads(1)
#ifdef WITH_MKL
      !$ call mkl_set_num_threads(1)
#endif
      !$omp do schedule(dynamic)
      do a = 1, nonfrozh
         ia=indx(a)
         do ic = 1, 3
            ii = (ia-1)*3+ic
            tmol = mol
            chk = chk0
            tmol%xyz(ic,ia)=tmol%xyz(ic,ia)+step
            call singlepoint &
               & (env,tmol,chk,calc, &
               &  egap,et,maxiter,0,.true.,.true.,acc,eel,gr,sr,sccr)
            tmol = mol
            chk = chk0
            dipd(1:3,ii)=sccr%dipole(1:3)
            pold(ii)=sccr%molpol
            mol%xyz(ic,ia)=mol%xyz(ic,ia)-2.*step
            call singlepoint &
               & (env,tmol,chk,calc, &
               &  egap,et,maxiter,0,.true.,.true.,acc,eel,gl,sl,sccl)
            tmol%xyz(ic,ia)=tmol%xyz(ic,ia)+step
            dipd(1:3,ii)=(dipd(1:3,ii)-sccl%dipole(1:3))*step2
            pold(ii)=(pold(ii)-sccl%molpol)*step2
            do b = 1, mol%n
               ja = indx(b)
               do jc = 1, 3
                  jj = (ja-1)*3 + jc
                  h(ii,jj) =(gr(jc,ja) - gl(jc,ja)) * step2
                  h(jj,ii) = h(ii,jj)  ! symmetric coupling part
               enddo
            enddo
         enddo
         if(a.eq.3)then
            call timing(t1,w1)
            write(*,'(''estimated CPU  time'',F10.2,'' min'')') &
               & 0.3333333d0*nonfrozh*(t1-t0)/60.
            write(*,'(''estimated wall time'',F10.2,'' min'')') &
               & 0.3333333d0*nonfrozh*(w1-w0)/60.
         endif
      enddo
      !$omp end do
      !$omp end parallel
      !$ call omp_set_num_threads(nproc)
#ifdef WITH_MKL
      !$ call mkl_set_num_threads(nproc)
#endif

   else
!! ------------------------------------------------------------------------
!  normal case
!! ------------------------------------------------------------------------
      !$ nproc = omp_get_num_threads()
      !$omp parallel default(shared) &
      !$omp firstprivate(et,maxiter,acc) &
      !$omp private(ia,ic,ii,ja,jc,jj,eel,gr,gl,egap,sccr,sccl,sr,sl,chk,tmol) &
      !$omp shared (mol,h,dipd,pold,step,step2,t1,t0,w1,w0,xyzsave,env,calc,chk0)
      !$ call omp_set_num_threads(1)
#ifdef WITH_MKL
      !$ call mkl_set_num_threads(1)
#endif
      !$omp do schedule(dynamic)
      do ia = 1, mol%n
         do ic = 1, 3
            ii = (ia-1)*3+ic

            tmol = mol
            chk = chk0 ! initialize wavefunction
            tmol%xyz(ic,ia)=xyzsave(ic,ia)+step

            gr = 0.0_wp
            eel = 0.0_wp
            call singlepoint &
               & (env,tmol,chk,calc, &
               &  egap,et,maxiter,-1,.true.,.true.,acc,eel,gr,sr,sccr)
            dipd(1:3,ii)=sccr%dipole(1:3)
            pold(ii)=sccr%molpol

            tmol = mol
            chk = chk0 ! reset wavefunction
            tmol%xyz(ic,ia)=xyzsave(ic,ia)-step

            gl = 0.0_wp
            eel = 0.0_wp
            call singlepoint &
               & (env,tmol,chk,calc, &
               &  egap,et,maxiter,-1,.true.,.true.,acc,eel,gl,sl,sccl)
            tmol%xyz(ic,ia)=xyzsave(ic,ia)
            dipd(1:3,ii)=(dipd(1:3,ii)-sccl%dipole(1:3))*step2
            pold(ii)=(pold(ii)-sccl%molpol)*step2

            do ja= 1, mol%n
               do jc = 1, 3
                  jj = (ja-1)*3 + jc
                  h(ii,jj) =(gr(jc,ja) - gl(jc,ja)) * step2
               enddo
            enddo

            call tmol%deallocate
         enddo

         if(ia.eq.3)then
            call timing(t1,w1)
            write(*,'(''estimated CPU  time'',F10.2,'' min'')') &
               & 0.3333333d0*mol%n*(t1-t0)/60.
            write(*,'(''estimated wall time'',F10.2,'' min'')') &
               & 0.3333333d0*mol%n*(w1-w0)/60.
         endif
      enddo
      !$omp end do
      !$omp end parallel
      !$ call omp_set_num_threads(nproc)
#ifdef WITH_MKL
      !$ call mkl_set_num_threads(nproc)
#endif

   endif

   else omp_parallize

   if(freezeset%n.gt.0) then
      ! for frozfc of about 10 the frozen modes
      ! approach 5000 cm-1, i.e., come too close to
      ! the real ones
      nonfrozh=mol%n-freezeset%n
      do a = 1,mol%n
         res%freq(a)=float(a)
      enddo
      do a=1,freezeset%n
         res%freq(freezeset%atoms(a))=freezeset%atoms(a)*100000d0
      enddo
      call sortind(mol%n,res%freq)
      do a=1,nonfrozh
         indx(a)=idint(res%freq(a))
      enddo
      do a=nonfrozh+1,mol%n
         indx(a)=idint(res%freq(a)/100000d0)
      enddo
      write(*,'(''atoms frozen in Hessian calc.:'',10i4)') &
         & indx(nonfrozh+1:mol%n)
      ! the index array indx(1:n) contains from 1:nonfrozh in ascending order
      ! the non-frozen atoms and from nonfrozh+1:n the frozen ones also in
      ! ascending order
      ! now compute a subblock of the Hessian
      do a = 1, nonfrozh
         ia=indx(a)
         do ic = 1, 3
            ii = (ia-1)*3+ic
            tmol = mol
            chk = chk0
            tmol%xyz(ic,ia)=tmol%xyz(ic,ia)+step
            call singlepoint &
               & (env,tmol,chk,calc, &
               &  egap,et,maxiter,0,.true.,.true.,acc,eel,gr,sr,sccr)
            tmol = mol
            chk = chk0
            dipd(1:3,ii)=sccr%dipole(1:3)
            pold(ii)=sccr%molpol
            mol%xyz(ic,ia)=mol%xyz(ic,ia)-2.*step
            call singlepoint &
               & (env,tmol,chk,calc, &
               &  egap,et,maxiter,0,.true.,.true.,acc,eel,gl,sl,sccl)
            tmol%xyz(ic,ia)=tmol%xyz(ic,ia)+step
            dipd(1:3,ii)=(dipd(1:3,ii)-sccl%dipole(1:3))*step2
            pold(ii)=(pold(ii)-sccl%molpol)*step2
            do b = 1, mol%n
               ja = indx(b)
               do jc = 1, 3
                  jj = (ja-1)*3 + jc
                  h(ii,jj) =(gr(jc,ja) - gl(jc,ja)) * step2
                  h(jj,ii) = h(ii,jj)  ! symmetric coupling part
               enddo
            enddo
         enddo
         if(a.eq.3)then
            call timing(t1,w1)
            write(*,'(''estimated CPU  time'',F10.2,'' min'')') &
               & 0.3333333d0*nonfrozh*(t1-t0)/60.
            write(*,'(''estimated wall time'',F10.2,'' min'')') &
               & 0.3333333d0*nonfrozh*(w1-w0)/60.
         endif
      enddo

   else
!! ------------------------------------------------------------------------
!  normal case
!! ------------------------------------------------------------------------
      do ia = 1, mol%n
         do ic = 1, 3
            ii = (ia-1)*3+ic

            tmol = mol
            chk = chk0 ! initialize wavefunction
            tmol%xyz(ic,ia)=xyzsave(ic,ia)+step

            gr = 0.0_wp
            eel = 0.0_wp
            call singlepoint &
               & (env,tmol,chk,calc, &
               &  egap,et,maxiter,-1,.true.,.true.,acc,eel,gr,sr,sccr)
            dipd(1:3,ii)=sccr%dipole(1:3)
            pold(ii)=sccr%molpol

            tmol = mol
            chk = chk0 ! reset wavefunction
            tmol%xyz(ic,ia)=xyzsave(ic,ia)-step

            gl = 0.0_wp
            eel = 0.0_wp
            call singlepoint &
               & (env,tmol,chk,calc, &
               &  egap,et,maxiter,-1,.true.,.true.,acc,eel,gl,sl,sccl)
            tmol%xyz(ic,ia)=xyzsave(ic,ia)
            dipd(1:3,ii)=(dipd(1:3,ii)-sccl%dipole(1:3))*step2
            pold(ii)=(pold(ii)-sccl%molpol)*step2

            do ja= 1, mol%n
               do jc = 1, 3
                  jj = (ja-1)*3 + jc
                  h(ii,jj) =(gr(jc,ja) - gl(jc,ja)) * step2
               enddo
            enddo

            call tmol%deallocate
         enddo

         if(ia.eq.3)then
            call timing(t1,w1)
            write(*,'(''estimated CPU  time'',F10.2,'' min'')') &
               & 0.3333333d0*mol%n*(t1-t0)/60.
            write(*,'(''estimated wall time'',F10.2,'' min'')') &
               & 0.3333333d0*mol%n*(w1-w0)/60.
         endif
      enddo

   endif

   end if omp_parallize

!  Hessian done -----------------------------------------------------------
!! ========================================================================

   if(freezeset%n.gt.0)then
      ! inverse mass array
      do a = 1, mol%n
         ia = indx(a)
         do ic = 1, 3
            ii = (ia-1)*3+ic
            isqm(ii)=1.0_wp/sqrt(atmass(ia))
         enddo
      enddo
      ! at this point the H matrix has zeros in the frozen atom block
      do a = nonfrozh+1,mol%n
         ia = indx(a)
         do ic = 1, 3
            ii = (ia-1)*3+ic
            h(ii,ii)=freezeset%fc   ! fill frozen diagonal block only
         enddo
      enddo
      res%hess = h ! copy
   else
      ! symmetrize
      do i=1,n3
         do j=1,n3
            res%hess(j,i)=(h(i,j)+h(j,i))*0.5
            if(abs(h(i,j)-h(j,i)).gt.1.d-2) then
               write(errStr,'(a,1x,i0,1x,i0,1x,a,1x,es14.6,1x,es14.6)') &
                  & 'Hessian element ',i,j,' is not symmetric:',h(i,j),h(j,i)
               call env%warning(trim(errStr), source)
            endif
         enddo
      enddo
      ! inverse mass array
      do ia = 1, mol%n
         do ic = 1, 3
            ii = (ia-1)*3+ic
            isqm(ii)=1.0_wp/sqrt(atmass(ia))
         enddo
      enddo
   endif

   if (pr_dftbp_hessian_out) then
      call writeHessianOut('hessian.out', res%hess)
      write(env%unit, '(A)') "DFTB+ style hessian.out written"
   end if

   ! prepare all for diag
   ! copy
   k=0
   do i=1,n3
      do j=1,i
         k=k+1
         hss(k)=res%hess(j,i)
      enddo
   enddo
   ! project
   if(.not.res%linear)then ! projection does not work for linear mol.
      call trproj(mol%n,n3,mol%xyz,hss,.false.,0,res%freq,1) ! freq is dummy
   endif
   ! non mass weigthed Hessian in hss
   hname = 'hessian'
   write(env%unit,'(a)')
   write(env%unit,'("writing file <",a,">.")') hname
   call wrhess(n3,hss,hname)

   ! include masses
   k=0
   do i=1,n3
      do j=1,i
         k=k+1
         res%hess(j,i)=hss(k)*isqm(i)*isqm(j)*scalh
         res%hess(i,j)=res%hess(j,i)
      enddo
   enddo
   ! diag
   lwork  = 1 + 6*n3 + 2*n3**2
   allocate(aux(lwork))
   call dsyev ('V','U',n3,res%hess,n3,res%freq,aux,lwork,info)
   if(info.ne.0) then
      call env%error('Diagonalization of hessian failed', source)
      return
   end if

   write(env%unit,'(a)')
   if(res%linear)then
      write(env%unit,'(1x,a)') 'vibrational frequencies (cm-1)'
   else
      write(env%unit,'(1x,a)') 'projected vibrational frequencies (cm-1)'
   endif
   k=0
   do i=1,n3
      res%freq(i)=autorcm*sign(sqrt(abs(res%freq(i))),res%freq(i))/sqrt(amutoau)
      if(abs(res%freq(i)).lt.0.01_wp) then
         k=k+1
         izero(k)=i
      endif
   enddo

   ! sort such that rot/trans are modes 1:6, H/isqm are scratch
   kend=6
   if(res%linear)then
      kend=5
      do i=1,kend
         izero(i)=i
      enddo
      res%freq(1:5)=0
   endif
   do k=1,kend
      h(1:n3,k)=res%hess(1:n3,izero(k))
      isqm(  k)=res%freq(izero(k))
   enddo
   j=kend
   do k=1,n3
      if(abs(res%freq(k)).gt.0.01_wp)then
         j=j+1
         if(j.gt.n3) then
            call env%error('internal error while sorting hessian', source)
            return
         end if
         h(1:n3,j)=res%hess(1:n3,k)
         isqm(  j)=res%freq(   k)
      endif
   enddo
   res%hess = h
   res%freq = isqm
   call PREIGF(env%unit,res%freq,res%n3true)

   ! reduced mass
   res%lowmode=1
   k=0
   do i=1,n3
      if(res%freq(i).lt.mode_vthr) res%lowmode=i
      xsum=0
      k=k+1
      do ia=1,mol%n
         do ic=1,3
            ii = (ia-1)*3+ic
            xsum=xsum+atmass(ia)*res%hess(ii,i)**2
         enddo
      enddo
      res%rmass(i)=xsum
   enddo
   ! IR intensity
   do i = 1, n3
      do k = 1, 3
         sum2 = 0.0_wp
         do j = 1, n3
            sum2 = sum2 + res%hess(j,i)*dipd(k,j)
         end do
         trdip(k) = sum2
      end do
      res%dipt(i) = sqrt (trdip(1)**2+trdip(2)**2+trdip(3)**2)
   end do
   ! Raman intensity
   do i = 1, n3
      sum2 = 0.0_wp
      do j = 1, n3
         sum2 = sum2 + res%hess(j,i)*pold(j)
      end do
      res%polt(i) = abs(sum2)
   end do

end subroutine numhess

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine rotmol(n,xyz,xrot,yrot,zrot)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only: pi
   implicit none
   integer :: n,i
   real(wp) :: xrot,yrot,zrot,xyz(3,n)
   real(wp) :: ang,xo,yo

   ang=xrot*pi/180.0_wp
   do i=1,n
      xo=xyz(2,i)
      yo=xyz(3,i)
      xyz(2,i)= xo*cos(ang)+yo*sin(ang)
      xyz(3,i)=-xo*sin(ang)+yo*cos(ang)
   enddo
   ang=yrot*pi/180.0_wp
   do i=1,n
      xo=xyz(1,i)
      yo=xyz(3,i)
      xyz(1,i)= xo*cos(ang)+yo*sin(ang)
      xyz(3,i)=-xo*sin(ang)+yo*cos(ang)
   enddo
   ang=zrot*pi/180.0_wp
   do i=1,n
      xo=xyz(1,i)
      yo=xyz(2,i)
      xyz(1,i)= xo*cos(ang)+yo*sin(ang)
      xyz(2,i)=-xo*sin(ang)+yo*cos(ang)
   enddo

end subroutine rotmol

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine distort(mol,freq,u)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_filetypes, only : generateFileName
   use xtb_type_molecule
   use xtb_io_writer, only : writeMolecule
   implicit none
   type(TMolecule), intent(inout) :: mol
   real(wp) u(:,:),freq(:)

   integer n3,i,imag,jj,ja,jc,ich
   real(wp) f,thr
   real(wp), allocatable :: xyz2(:,:)
   character(len=:), allocatable :: fname

   xyz2 = mol%xyz

   n3=3*len(mol)
   ! cut-off for what is considered to be imag
   thr=5.0

   imag=0
   do i=1,n3
      if(freq(i).lt.0.and.abs(freq(i)).gt.thr)then
         imag=imag+1
      endif
   enddo
   if(imag.eq.0) return
   write(*,'(''imag cut-off (cm-1) :'',f8.2)') thr
   if(imag.eq.1) &
      & write(*,*) 'found ',imag,' significant imaginary frequency'
   if(imag.gt.1) &
      & write(*,*) 'found ',imag,' significant imaginary frequencies'

   ! magnitude of distortion
   f=0.5/float(imag)

   do i=1,n3
      if(freq(i).lt.0.and.abs(freq(i)).gt.thr)then
         do ja=1,len(mol)
            do jc=1,3
               jj = (ja-1)*3 + jc
               mol%xyz(jc,ja)=mol%xyz(jc,ja)+f*u(jj,i)
            enddo
         enddo
      endif
   enddo

   call generateFileName(fname, 'xtbhess', '', mol%ftype)
   write(*,*) 'writing imag mode distorted coords to '//fname
   write(*,*) 'for further optimization.'
   call open_file(ich, fname, 'w')
   call writeMolecule(mol, ich)
   call close_file(ich)

   mol%xyz = xyz2

end subroutine distort

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine write_tm_vibspectrum(ich,n3,freq,ir_int)
   integer, intent(in)  :: ich ! file handle
   integer, intent(in)  :: n3
   real(wp),intent(in)  :: freq(n3)
   real(wp),intent(in)  :: ir_int(n3)
   integer  :: i
   real(wp) :: thr=0.01_wp
   write(ich,'("$vibrational spectrum")')
   write(ich,'("#  mode     symmetry     wave number   IR intensity    selection rules")')
   write(ich,'("#                         cm**(-1)        (amu)          IR     RAMAN")')
   do i = 1, n3
      if (abs(freq(i)).lt.thr) then
      write(ich,'(i6,9x,    f18.2,f16.5,7x," - ",5x," - ")') &
         i,freq(i),0.0_wp
      else
      write(ich,'(i6,8x,"a",f18.2,f16.5,7x,"YES",5x,"YES")') &
         i,freq(i),ir_int(i)
      endif
   enddo
   write(ich,'("$end")')
end subroutine

subroutine g98fake2(fname,n,at,xyz,freq,red_mass,ir_int,u2)
   use xtb_mctc_accuracy, only : wp
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: freq(3*n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: u2(3*n,3*n)
   character(len=*),intent(in) :: fname
   real(wp),intent(in)  :: red_mass(3*n)
   real(wp),intent(in)  :: ir_int  (3*n)

   integer  :: gu,i,j,ka,kb,kc,la,lb,k
   character(len=2) :: irrep
   real(wp),allocatable :: u(:,:)
   real(wp),allocatable :: red(:)
   real(wp),allocatable :: f2 (:)
   real(wp),allocatable :: ir (:)
   real(wp) :: zero

   allocate( u(3*n,3*n), red(3*n), f2(3*n), ir(3*n), source = 0.0_wp )

   irrep='a'
   zero    =0.0

   k=0
   do i=1,3*n
      if(abs(freq(i)).gt.1.d-1)then
         k=k+1
         u(1:3*n,k)=u2(1:3*n,i)
         f2(k)=freq(i)
         ir(k)=ir_int(i)
         red(k)=red_mass(i)
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
   write (gu,*) '---------------------------------------------', &
      & '-----------------------'
   write (gu,*) ' Center     Atomic     Atomic', &
      & '              Coordinates (Angstroms)'
   write (gu,*) ' Number     Number      Type ', &
      & '             X           Y           Z'
   write (gu,*) '-----------------------', &
      & '---------------------------------------------'
   j=0
   do i=1,n
      write(gu,111) i,at(i),j,xyz(1:3,i)*0.52917726
   enddo
   write (gu,*) '----------------------', &
      & '----------------------------------------------'
   write (gu,*) '    1 basis functions        1 primitive gaussians'
   write (gu,*) '    1 alpha electrons        1 beta electrons'
   write (gu,*)
   111 format(i5,i11,i14,4x,3f12.6)

   write (gu,*) 'Harmonic frequencies (cm**-1), IR intensities', &
      & ' (KM/Mole),'
   write (gu,*) 'Raman scattering activities (A**4/amu),', &
      & ' Raman depolarization ratios,'
   write (gu,*) 'reduced masses (AMU), force constants (mDyne/A)', &
      & ' and normal coordinates:'

   ka=1
   kc=3
   60  kb=min0(kc,k)
   write (gu,100) (j,j=ka,kb)
   write (gu,105) (irrep,j=ka,kb)
   write (gu,110) ' Frequencies --',(f2(j),j=ka,kb)
   write (gu,110) ' Red. masses --',(red(j),j=ka,kb)
   write (gu,110) ' Frc consts  --',(zero,j=ka,kb)
   write (gu,110) ' IR Inten    --',(ir(j),j=ka,kb)
   write (gu,110) ' Raman Activ --',(zero,j=ka,kb)
   write (gu,110) ' Depolar     --',(zero,j=ka,kb)
   write (gu,*)'Atom AN      X      Y      Z        X      Y', &
      & '      Z        X      Y      Z'
   la=1
   70  lb=n
   do  i=la,lb
      write (gu,130) i,at(i), (u(i*3-2,j),  u(i*3-1,j),  u(i*3  ,j),j=ka,kb)
   enddo
   if (lb.eq.n) go to 90
   go to 70
   90  if (kb.eq.k) then
      return
   endif
   ka=kc+1
   kc=kc+3
   go to 60

   100 format (3(20x,i3))
   105 format (3x,3(18x,a5))
   110 format (a15,f11.4,12x,f11.4,12x,f11.4)
   130 format (2i4,3(2x,3f7.2))

   write(gu,'(''end of file'')')
   call close_file(gu)
   return

end subroutine g98fake2

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine g98fake(fname,n,at,xyz,freq,u2)
   use xtb_mctc_accuracy, only : wp
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: freq(3*n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp),intent(in)  :: u2(3*n,3*n)
   character(len=*),intent(in) :: fname

   integer  :: gu,i,j,ka,kb,kc,la,lb,k
   character(len=2) :: irrep
   real(wp),allocatable :: u(:,:)
   real(wp),allocatable :: red_mass(:)
   real(wp),allocatable :: force   (:)
   real(wp),allocatable :: ir_int  (:)
   real(wp),allocatable :: f2      (:)
   real(wp) :: zero

   allocate( u(3*n,3*n), red_mass(3*n), force(3*n), ir_int(3*n), f2(3*n), &
             source = 0.0_wp )

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
   write (gu,*) '---------------------------------------------', &
      & '-----------------------'
   write (gu,*) ' Center     Atomic     Atomic', &
      & '              Coordinates (Angstroms)'
   write (gu,*) ' Number     Number      Type ', &
      & '             X           Y           Z'
   write (gu,*) '-----------------------', &
      & '---------------------------------------------'
   j=0
   do i=1,n
      write(gu,111) i,at(i),j,xyz(1:3,i)*0.52917726
   enddo
   write (gu,*) '----------------------', &
      & '----------------------------------------------'
   write (gu,*) '    1 basis functions        1 primitive gaussians'
   write (gu,*) '    1 alpha electrons        1 beta electrons'
   write (gu,*)
   111 format(i5,i11,i14,4x,3f12.6)

   write (gu,*) 'Harmonic frequencies (cm**-1), IR intensities',' (KM/Mole),'
   write (gu,*) 'Raman scattering activities (A**4/amu),', &
      & ' Raman depolarization ratios,'
   write (gu,*) 'reduced masses (AMU), force constants (mDyne/A)', &
      & ' and normal coordinates:'

   ka=1
   kc=3
   60  kb=min0(kc,k)
   write (gu,100) (j,j=ka,kb)
   write (gu,105) (irrep,j=ka,kb)
   write (gu,110) ' Frequencies --',(f2(j),j=ka,kb)
   write (gu,110) ' Red. masses --',(red_mass(j),j=ka,kb)
   write (gu,110) ' Frc consts  --',(force(j),j=ka,kb)
   write (gu,110) ' IR Inten    --',(ir_int(j),j=ka,kb)
   write (gu,110) ' Raman Activ --',(zero,j=ka,kb)
   write (gu,110) ' Depolar     --',(zero,j=ka,kb)
   write (gu,*)'Atom AN      X      Y      Z        X      Y', &
      & '      Z        X      Y      Z'
   la=1
   70  lb=n
   do  i=la,lb
      write (gu,130) i,at(i), (u(i*3-2,j),  u(i*3-1,j),  u(i*3  ,j),j=ka,kb)
   enddo
   if (lb.eq.n) go to 90
   go to 70
   90  if (kb.eq.k) then
      return
   endif
   ka=kc+1
   kc=kc+3
   go to 60

   100 format (3(20x,i3))
   105 format (3x,3(18x,a5))
   110 format (a15,f11.4,12x,f11.4,12x,f11.4)
   130 format (2i4,3(2x,3f7.2))

   write(gu,'(''end of file'')')
   call close_file(gu)
   return

end subroutine g98fake

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine trproj(natoms,nat3,xyz,hess,ldebug,nmode,mode,ndim)

   !----------------------------------------------------------------------
   !  subroutine trproj drives projection of hessian out of the
   !  space of translational and rotational motions:
   !  first get xyz c.m.; second get transl. and rot. projection matrix
   !
   !  get center of mass coordinates with unit mass
   !
   ! Input
   !   natoms  = number of atoms
   !   nat3    = 3*natoms
   !   xyz     = cartesian coordinates
   !       ldebug  = debug flag = .true. for debugging
   !
   ! Ouput fo gtrprojm.f
   !   xyzucm  = temporary c.m. coordinates
   !
   !       hess        = projected hessian out of space of transl. and rot.
   !                 motion
   !
   !----------------------------------------------------------------------

   implicit none

   ! Input
   logical, intent(in) :: ldebug
   integer, intent(in) :: natoms,nat3,nmode,ndim
   real(8), dimension(3,natoms) :: xyz
   real(8), dimension(nat3,ndim):: mode
   ! Ouput
   real(8), dimension(nat3*(nat3+1)/2) :: hess
   ! Local
   integer :: i
   real(8) :: xm,ym,zm
   real(8), dimension(3,natoms) ::xyzucm

   xyzucm(:,:) = xyz(:,:)

   xm = 0.0d0
   ym = 0.0d0
   zm = 0.0d0

   do i=1,natoms
      xm = xm + xyzucm(1,i)
      ym = ym + xyzucm(2,i)
      zm = zm + xyzucm(3,i)
   end do

   xm = xm/natoms
   ym = ym/natoms
   zm = zm/natoms


   do i=1,natoms
      xyzucm(1,i) = xyzucm(1,i) - xm
      xyzucm(2,i) = xyzucm(2,i) - ym
      xyzucm(3,i) = xyzucm(3,i) - zm
   end do

   ! get translational and rotational projection matrix

   call gtrprojm(natoms,nat3,xyzucm,hess,ldebug,nmode,mode,ndim)

end subroutine trproj

!----------------------------------------------------------------------
subroutine gtrprojm(natoms,nat3,xyzucm,hess,ldebug,nmode,mode,ndim)
   !----------------------------------------------------------------------
   ! calculating the translational-rotational projection matrix
   !
   ! Input
   !   natoms  = number of atoms
   !   nat3    = 3*natoms
   !   xyzucm  = coords c.m. from gxyzucm.f
   !           hess    = hessian
   !       ldebug      = debug flag = .true. for debugging
   !
   ! Ouput
   !   fmat    = F-matrix with translational and rotational vectors
   !       pmat        = projection matrix P = (1-FFt)
   !   hess    = projected hessian
   !----------------------------------------------------------------------
   use xtb_fixparam
   use xtb_mctc_la, only : blckmgs,syprj

   implicit none

   ! Input
   logical, intent(in) :: ldebug
   integer, intent(in) :: natoms,nat3,nmode,ndim
   real(8), dimension(3,natoms) :: xyzucm
   real(8), dimension(nat3,ndim):: mode
   ! Ouput
   real(8), dimension(nat3*(nat3+1)/2) :: hess

   ! Local
   integer :: i,ii,iii
   real(8), allocatable :: fmat(:,:)
   integer :: nprj

   nprj=6
   if(nmode.gt.0) nprj=nprj+nmode
   if(nmode.lt.0) nprj=nprj+fixset%n*3
   allocate(fmat(nat3,nprj))
   fmat(:,:) = 0.0d0

   if(nmode.ge.0) then
      do i=1,natoms
         do ii=1,3
            !        translation vectors
            fmat(3*(i-1)+ii,ii) = 1.0d0
         end do
         !        rotational vectors
         fmat(3*(i-1)+1,4) =  0.0d0
         fmat(3*(i-1)+2,4) = -xyzucm(3,i)
         fmat(3*(i-1)+3,4) =  xyzucm(2,i)
         fmat(3*(i-1)+1,5) =  xyzucm(3,i)
         fmat(3*(i-1)+2,5) =  0.0d0
         fmat(3*(i-1)+3,5) = -xyzucm(1,i)
         fmat(3*(i-1)+1,6) = -xyzucm(2,i)
         fmat(3*(i-1)+2,6) =  xyzucm(1,i)
         fmat(3*(i-1)+3,6) =  0.0d0
      end do
   endif

   if(nmode.gt.0) then  ! NMF
      do i=1,nmode
         fmat(1:nat3,6+i)=mode(1:nat3,i)
      enddo
   endif

   if(nmode.lt.0) then ! exact fixing
      do i=1,natoms
         !        rotational vectors
         fmat(3*(i-1)+1,1) =  0.0d0
         fmat(3*(i-1)+2,1) = -xyzucm(3,i)
         fmat(3*(i-1)+3,1) =  xyzucm(2,i)
         fmat(3*(i-1)+1,2) =  xyzucm(3,i)
         fmat(3*(i-1)+2,2) =  0.0d0
         fmat(3*(i-1)+3,2) = -xyzucm(1,i)
         fmat(3*(i-1)+1,3) = -xyzucm(2,i)
         fmat(3*(i-1)+2,3) =  xyzucm(1,i)
         fmat(3*(i-1)+3,3) =  0.0d0
      enddo
      do i=1,fixset%n
         iii=fixset%atoms(i)
         do ii=1,3
            fmat(3*(iii-1)+ii,3+(i-1)*3+ii) = 1.0d0
         end do
      enddo
   endif

   if(ldebug) then
      write(*,'(a)')
      write(*,'(a)') ' Basis vectors before orthonormalization'
      write(*,'(3e22.14)') fmat
   end if

   ! do orthogonalization

   call  blckmgs(nat3,nprj,nat3,fmat)

   ! do projection

   call syprj(nat3,nprj,fmat,nat3,hess)

   deallocate(fmat)

end subroutine gtrprojm


subroutine wrhess(nat3,h,fname)
   use xtb_lin, only : lin
   implicit none
   integer, intent(in) :: nat3
   real(wp),intent(in) :: h(nat3*(nat3+1)/2)
   character(len=*),intent(in) :: fname
   integer iunit,i,j,mincol,maxcol,k
   character(len=5)  :: adum
   character(len=80) :: a80

   adum='   '
   call open_file(iunit,fname,'w')
   a80='$hessian'
   write(iunit,'(a)')a80
   do i=1,nat3
      maxcol = 0
      k=0
      200    mincol = maxcol + 1
      k=k+1
      maxcol = min(maxcol+5,nat3)
      write(iunit,'(a5,5f15.10)')adum,(h(lin(i,j)),j=mincol,maxcol)
      if (maxcol.lt.nat3) goto 200
   enddo
   call close_file(iunit)

end subroutine wrhess

subroutine rdhess(nat3,h,fname)
   implicit none
   integer, intent(in)  :: nat3
   real(wp),intent(out) :: h(nat3,nat3)
   character(len=*),intent(in) :: fname
   integer  :: iunit,i,j,mincol,maxcol
   character(len=5)  :: adum
   character(len=80) :: a80

   !     write(*,*) 'Reading Hessian <',trim(fname),'>'
   call open_file(iunit,fname,'r')
   50  read(iunit,'(a)')a80
   if(index(a80,'$hessian').ne.0)then
      do i=1,nat3
         maxcol = 0
         200       mincol = maxcol + 1
         maxcol = min(maxcol+5,nat3)
         read(iunit,*)(h(j,i),j=mincol,maxcol)
         if (maxcol.lt.nat3) goto 200
      enddo
      call close_file(iunit)
      goto 300
   endif
   goto 50

   300 return
end subroutine rdhess

pure subroutine sortind(nvar,edum)
   implicit none
   integer, intent(in)    :: nvar
   real(wp),intent(inout) :: edum(nvar)
   integer  :: ii,k,i,j
   real(wp) :: pp

   do ii = 2, nvar
      i = ii - 1
      k = i
      pp= edum(i)
      do j = ii, nvar
         if (edum(j) .gt. pp) cycle
         k = j
         pp= edum(j)
      enddo
      if (k .eq. i) cycle
      edum(k) = edum(i)
      edum(i) = pp
   enddo

end subroutine sortind


!> Write the second derivative matrix
subroutine writeHessianOut(fileName, pDynMatrix)

   !> File name
   character(*), intent(in) :: fileName

   !> Dynamical (Hessian) matrix
   real(wp), intent(in) :: pDynMatrix(:,:)

   !> Format string for energy second derivative matrix
   character(len=*), parameter :: formatHessian = '(4f16.10)'

   integer :: ii, fd

   call open_file(fd, fileName, 'w')
   do ii = 1, size(pDynMatrix, dim=2)
      write(fd, formatHessian) pDynMatrix(:, ii)
   end do
   call close_file(fd)

end subroutine writeHessianOut

end module xtb_hessian
