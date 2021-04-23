! This file is part of xtb.
!
! Copyright (C) 2017-2021 Stefan Grimme
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
   use xtb_freq_io, only : rdhess, wrhess, writeHessianOut, &
      & write_tm_vibspectrum, g98fake, g98fake2
   use xtb_freq_numdiff, only : numdiff2
   use xtb_freq_project, only : trproj
   use xtb_freq_turbomole, only : aoforce_hessian

contains

subroutine numhess( &
      & env,mol,chk0,calc, &
      & egap,et,maxiter,etot,gr,sr,res)
   use xtb_mctc_accuracy, only : wp
!$ use omp_lib

   use xtb_mctc_convert
   use xtb_mctc_blas

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
   use xtb_metadynamic

   use xtb_single, only : singlepoint
   use xtb_axis, only : axis

   use xtb_gfnff_calculator, only : TGFFCalculator
   use xtb_type_dummycalc, only : TDummyCalculator

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
   real(wp) :: ebias
   real(wp) :: alp1,alp2
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
   real(wp),allocatable :: htb (:,:)
   real(wp),allocatable :: hbias (:,:)
   real(wp),allocatable :: hss(:)
   real(wp),allocatable :: hsb(:)
   real(wp),allocatable :: fc_tb(:)
   real(wp),allocatable :: fc_bias(:)
   real(wp),allocatable :: v(:)
   real(wp),allocatable :: fc_tmp(:)
   real(wp),allocatable :: freq_scal(:)
   real(wp),allocatable :: aux (:)
   real(wp),allocatable :: isqm(:)
   real(wp),allocatable :: gl  (:,:)
   real(wp),allocatable :: xyzsave(:,:)
   real(wp),allocatable :: pold(:)
   real(wp),allocatable :: dipd(:,:)
   real(wp),allocatable :: amass(:)

   type(TMolecule) :: tmol

   logical :: ex,rd
   integer :: ich ! file handle
   integer :: err
   character(len=:),allocatable :: hname,fname
   logical :: exitRun
   character(len=128) :: errStr

   n3=3*mol%n
   call res%allocate(mol%n)
   res%n3true = n3-3*freezeset%n

   allocate(hss(n3*(n3+1)/2),hsb(n3*(n3+1)/2),h(n3,n3),htb(n3,n3),hbias(n3,n3), &
      & gl(3,mol%n),isqm(n3),xyzsave(3,mol%n),dipd(3,n3), &
      & pold(n3),nb(20,mol%n),indx(mol%n),molvec(mol%n),bond(mol%n,mol%n), &
      & v(n3),fc_tmp(n3),freq_scal(n3),fc_tb(n3),fc_bias(n3),amass(n3))

   rd=.false.
   xyzsave = mol%xyz

   step=0.0001_wp
   call rotmol(mol%n,mol%xyz,step,2.*step,3.*step)

   ! step length
   step=step_hess
   if(extcode.eq.5) step=step*2.0_wp ! MOPAC is not very accurate
   ! SCC accuraccy
   acc=accu_hess
   scalh=scale_hess

   call singlepoint &
      & (env,mol,chk0,calc, &
      &  egap,et,maxiter,0,.true.,.true.,acc,res%etot,res%grad,sr,sccr)

   if (runtyp.eq.p_run_bhess) then
   write(env%unit,'(''kpush                :'',F10.5)') metaset%factor(metaset%nstruc)
   write(env%unit,'(''alpha                :'',F10.5)') metaset%global_width
   end if
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
   step2=0.5_wp/step

   h = 0.0_wp
   htb = 0.0_wp
   hbias = 0.0_wp

!! ========================================================================
!  Hessian part -----------------------------------------------------------

   !analytical hessian calculation
   if(mode_extrun .eq. p_ext_turbomole) then    
           dipd=0.0_wp
           call aoforce_hessian(env,mol,h,dipd) 
   else !numerical hessian calculation
   parallize = .true.
   select type(calc)
   type is (TDummyCalculator)
      parallize = .false.
   end select

   if(freezeset%n.gt.0) then
      ! for frozfc of about 10 the frozen modes
      ! approach 5000 cm-1, i.e., come too close to
      ! the real ones
      nonfrozh=mol%n-freezeset%n
      do a = 1,mol%n
         res%freq(a)=float(a)
      enddo
      do a=1,freezeset%n
         res%freq(freezeset%atoms(a))=freezeset%atoms(a)*100000_wp
      enddo
      call sortind(mol%n,res%freq)
      do a=1,nonfrozh
         indx(a)=idint(res%freq(a))
      enddo
      do a=nonfrozh+1,mol%n
         indx(a)=idint(res%freq(a)/100000_wp)
      enddo
      write(*,'(''atoms frozen in Hessian calc.:'',10i4)') &
         & indx(nonfrozh+1:mol%n)

      h = 0.0_wp
      dipd = 0.0_wp
      pold = 0.0_wp
      call numdiff2(env, mol, chk0, calc, indx(:nonfrozh), step, h, dipd, parallize)

   else
!! ------------------------------------------------------------------------
!  normal case
!! ------------------------------------------------------------------------
      h = 0.0_wp
      dipd = 0.0_wp
      pold = 0.0_wp
      call numdiff2(env, mol, chk0, calc, step, h, dipd, parallize)
   endif

   endif !From turbomole_exception


!  Hessian done -----------------------------------------------------------
!! ========================================================================

   if (runtyp.eq.p_run_bhess) call numhess_rmsd(env,mol,hbias)

   if(mode_extrun .eq. p_ext_turbomole .AND. runtyp.eq.p_run_bhess) then 
        h = h + hbias !h is biased
   end if

   if(freezeset%n.gt.0)then
      ! inverse mass array
      do a = 1, mol%n
         ia = indx(a)
         do ic = 1, 3
            ii = (ia-1)*3+ic
            isqm(ii)=1.0_wp/sqrt(atmass(ia))
            amass(ii)=isqm(ii)/sqrt(amutoau)
         enddo
      enddo
      do a = 1, nonfrozh
         ia = indx(a)
         do ic = 1, 3
            ii = (ia-1)*3+ic
            do b = 1, nonfrozh
               ja = indx(b)
               do jc = 1, 3
                  jj = (ja-1)*3+jc
                  if(abs(h(ii,jj)-h(jj,ii)).gt.1.d-2) then
                     write(errStr,'(a,1x,i0,1x,i0,1x,a,1x,es14.6,1x,es14.6)') &
                        & 'Hessian element ',i,j,' is not symmetric:',h(i,j),h(j,i)
                     call env%warning(trim(errStr), source)
                  endif
                  h(jj,ii) = 0.5_wp*(h(ii,jj)+h(jj,ii))
               enddo
            enddo
         enddo
      enddo
      ! at this point the H matrix has zeros in the frozen atom block
      do a = nonfrozh+1,mol%n
         ia = indx(a)
         do ic = 1, 3
            ii = (ia-1)*3+ic
            h(:, ii) = 0.0_wp
            h(ii, :) = 0.0_wp
            h(ii,ii)=freezeset%fc   ! fill frozen diagonal block only
         enddo
      enddo
      res%hess = h ! copy
   else
      ! symmetrize
      do i=1,n3
         do j=1,n3
            res%hess(j,i)=(h(i,j)+h(j,i))*0.5_wp
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
            amass(ii)=isqm(ii)/sqrt(amutoau)
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
   ! same for bhess run
   if (runtyp.eq.p_run_bhess) then
      k=0
      do i=1,n3
         do j=1,i
            k=k+1
            hsb(k)=hbias(j,i)
         enddo
      enddo
   end if
   ! project
   if(.not.res%linear)then ! projection does not work for linear mol.
      if (runtyp.eq.p_run_bhess) then
         call trproj(mol%n,n3,mol%xyz,hsb,.false.,0,res%freq,1) ! freq is dummy
      end if
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
   ! same for bhess run
   if (runtyp.eq.p_run_bhess) then
      k=0
      do i=1,n3
         do j=1,i
            k=k+1
            hbias(j,i)=hsb(k)*isqm(i)*isqm(j)*scalh
            hbias(i,j)=hbias(j,i)
         enddo
      enddo
   end if
   ! calcualte htb without RMSD bias
   if (runtyp.eq.p_run_bhess) htb=res%hess-hbias
   ! diag
   lwork  = 1 + 6*n3 + 2*n3**2
   allocate(aux(lwork))
   call dsyev ('V','U',n3,res%hess,n3,res%freq,aux,lwork,info)
   if(info.ne.0) then
      call env%error('Diagonalization of hessian failed', source)
      return
   end if

   ! calculate fc_tb and fc_bias
   alp1=1.27_wp
   alp2=1.5d-4
   if (runtyp.eq.p_run_bhess) then
      do j=1,n3
         v(1:n3) = res%hess(1:n3,j) ! modes
         call mctc_gemv(htb,v,fc_tmp)
         fc_tb(j) = mctc_dot(v,fc_tmp)
         call mctc_gemv(hbias,v,fc_tmp)
         fc_bias(j) = mctc_dot(v,fc_tmp)
         if (abs(res%freq(j)).gt.1.0d-6) then
            freq_scal(j) = sqrt( (fc_tb(j)+alp2) / ( (fc_tb(j)+alp2) +  alp1*fc_bias(j) ) )
            if (fc_tb(j).lt.0.and.fc_bias(j).ne.0) then
               freq_scal(j) = -sqrt( (abs(fc_tb(j))+alp2) / ( (abs(fc_tb(j))+alp2) + alp1*fc_bias(j) ) )
            end if
         else
            freq_scal(j) = 1.0_wp
         end if
      end do
   end if

   write(env%unit,'(a)')
   if(res%linear)then
      write(env%unit,'(1x,a)') 'vibrational frequencies (cm⁻¹)'
   else
      write(env%unit,'(1x,a)') 'projected vibrational frequencies (cm⁻¹)'
   endif
   k=0
   do i=1,n3
      ! Eigenvalues in atomic units, convert to wavenumbers
      res%freq(i)=autorcm*sign(sqrt(abs(res%freq(i))),res%freq(i))/sqrt(amutoau)
      if(abs(res%freq(i)).lt.0.01_wp) then
         k=k+1
         izero(k)=i
      endif
   enddo

   ! scale frequencies
   if (runtyp.eq.p_run_bhess) then
      do j=1,n3
         res%freq(j)=freq_scal(j)*res%freq(j)
      end do
   end if

   if (verbose.and.runtyp.eq.p_run_bhess) then
      write(env%unit,'(4x,"freq   fc_tb      fc_bias    scal")')
      do i=1,n3
         write(env%unit,'(f8.2,2x,f9.6,2x,f9.6,2x,f7.4)') &
         res%freq(i),fc_tb(i),fc_bias(i),freq_scal(i)
      end do
      write(env%unit,*)
   end if

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

   !--- IR intensity ---!
   !  1. res%hess corresponds to the orthonormal eigenvectors of the hessian
   !     matrix (-> normal modes of vibration). Mass-weighting is introduced
   !     back again via multiplying with amass(j).
   !
   !  2. res%hess(j,i) is the matrix which transforms a derivative with
   !     respect to the j-th cartesian coordinate ("dipd") into a derivative with
   !     respect to the i-th (mass-weighted) normal coordinate.
   !
   !  3. amass(j) = 1/sqrt(m(j)); m(j) is given in atomic units (a.u.).
   !
   !  4. matmul(D x H) = U
   !
   !  5. D = dipd(3,n3); H = res%hess(n3:n3); U = Matrix with dipol derivatives
   !                                              in x, y and z direction per mode

   do i = 1, n3
      do k = 1, 3
         sum2 = 0.0_wp
         do j = 1, n3
            sum2 = sum2 + dipd(k,j)*(res%hess(j,i)*amass(j))
         end do
         trdip(k) = sum2
      end do
      res%dipt(i) = autokmmol*(trdip(1)**2+trdip(2)**2+trdip(3)**2)
   end do
   ! Raman intensity
   do i = 1, n3
      sum2 = 0.0_wp
      do j = 1, n3
         sum2 = sum2 + pold(j)*(res%hess(j,i)*amass(j))
      end do
      res%polt(i) = abs(sum2)
   end do

end subroutine numhess

subroutine numhess_rmsd( &
      & env,mol,hbias)
   use xtb_mctc_accuracy, only : wp
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
   use xtb_metadynamic

   implicit none
   !> Dummy
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(inout) :: mol
   real(wp),intent(inout)         :: hbias(mol%n*3,mol%n*3)
   !> Stack
   type(TMolecule) :: tmol
   real(wp) :: ebias
   real(wp) :: step,step2
   integer  :: n3,i,j,k,ic,jc,ia,ja,ii,jj,a,b
   real(wp),allocatable :: gr(:,:)
   real(wp),allocatable :: gl(:,:)
   real(wp),allocatable :: xyzsave(:,:)

   n3=3*mol%n

   allocate(gr(3,mol%n),gl(3,mol%n),xyzsave(3,mol%n))

   xyzsave = mol%xyz

   ! step length
   step=0.0001_wp
   step=step_hess
   step2=0.5_wp/step

!! ========================================================================
!  RMSD part -----------------------------------------------------------

   do ia = 1, mol%n
      do ic = 1, 3
         ii = (ia-1)*3+ic

         tmol=mol
         tmol%xyz(ic,ia)=xyzsave(ic,ia)+step

         gr = 0.0_wp
         ebias = 0.0_wp
         call metadynamic(metaset,tmol%n,tmol%at,tmol%xyz,ebias,gr)

         tmol=mol
         tmol%xyz(ic,ia)=xyzsave(ic,ia)-step

         gl = 0.0_wp
         ebias = 0.0_wp
         call metadynamic(metaset,tmol%n,tmol%at,tmol%xyz,ebias,gl)

         tmol%xyz(ic,ia)=xyzsave(ic,ia)

         do ja= 1, mol%n
            do jc = 1, 3
               jj = (ja-1)*3 + jc
               hbias(ii,jj) =(gr(jc,ja) - gl(jc,ja)) * step2
            enddo
         enddo

         call tmol%deallocate
      enddo

   end do

!  RMSD done -----------------------------------------------------------
!! ========================================================================

end subroutine numhess_rmsd

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

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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

end module xtb_hessian
