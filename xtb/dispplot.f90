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

subroutine dispersion_surfaceplot(nat,at,xyz,q,wf,g_a,g_c,lmbd,dfparam,surface)
   use iso_fortran_env, wp => real64
   use setparam
   use tbdef_param
   use grid_module
   use dftd4
   use printout
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: wf,g_a,g_c
   real(wp),intent(in)  :: q(nat)
   type(dftd_parameter),intent(in) :: dfparam
   integer, intent(in)  :: lmbd
   type(tb_grid),intent(inout) :: surface

   integer  :: i,j,k
   ! local copies
   integer  :: ndim
   real(wp) :: etmp,embd,save_tmp,save_mbd
   integer, allocatable :: at_probe(:)
   real(wp),allocatable :: xyz_probe(:,:)
   real(wp),allocatable :: covcn(:)
   real(wp),allocatable :: q_probe(:)
   real(wp),allocatable :: gw(:)
   real(wp),allocatable :: c6abns(:,:)
   real(wp),allocatable :: dummy(:,:)
!  real(wp),allocatable :: disp(:)

   ! parameters
   integer, parameter   :: probe_atom(1) = (/1/) ! Ar

   if (allocated(surface%rho)) deallocate(surface%rho)
   allocate(surface%rho(surface%n,2), source = 0.0_wp)

   ! make a local copy with extended size
   allocate( at_probe(nat+1), source = 0 )
   allocate( xyz_probe(3,nat+1),q_probe(nat+1),covcn(nat+1), source = 0.0_wp )
   at_probe(1:nat) = at
   at_probe(nat+1) = probe_atom(1)
   xyz_probe(:,1:nat) = xyz
   xyz_probe(:,nat+1) = (/(1.0e10_wp,k=1,3)/)
   q_probe(1:nat) = q
   q_probe(nat+1) = 0.0_wp

   call d4init(nat+1,at_probe,g_a,g_c,p_refq_gfn2xtb,ndim)
   allocate( c6abns(ndim,ndim),gw(ndim),dummy(ndim,ndim), source = 0.0_wp )
   call covncoord(nat,at,xyz,covcn,1600.0_wp)
   covcn(nat+1) = 0.0_wp
   call d4(nat+1,ndim,at_probe,wf,g_a,g_c,covcn,gw,c6abns)
   deallocate(dummy)
   xyz_probe(:,nat+1) = (/(1.0e10_wp,k=1,3)/)
   call edisp(nat+1,ndim,at_probe,q_probe,xyz_probe,dfparam,g_a,g_c, &
   &          gw,c6abns,lmbd,save_tmp,emany=save_mbd)

   !$omp parallel default(none) &
   !$omp shared(at_probe,nat,q_probe,dfparam,wf,g_a,g_c,c6abns,surface,ndim, &
   !$omp        gw,lmbd,save_tmp,save_mbd,xyz) &
   !$omp firstprivate(xyz_probe) &
   !$omp private(k,etmp,embd)
   !$omp do
   do k = 1, surface%n
      xyz_probe(:,:nat) = xyz
      xyz_probe(:,nat+1) = surface%x(:,k)
      call edisp(nat+1,ndim,at_probe,q_probe,xyz_probe,dfparam,g_a,g_c, &
      &          gw,c6abns,lmbd,Etmp,emany=embd)
      surface%rho(k,1) = etmp-save_tmp
      surface%rho(k,2) = embd-save_mbd
!      call dispgrad(nat+1,ndim,at_probe,q,xyz_probe,dfparam,wf,g_a,g_c, &
!      &             c6abns,3,g,etmp)
      
   enddo
   !$omp enddo
   !$omp end parallel

end subroutine dispersion_surfaceplot
!subroutine dispersion_surfaceplot
!   use iso_fortran_env, wp => real64
!   use grid_module
!   use printout, only : writecosmofile
!   implicit none
!   type(tb_grid) :: surface
!   integer  :: nat
!   integer  :: at(nat)
!   real(wp) :: xyz(3,nat)
!
!   ! loop indices
!   integer  :: i,ii,ij,ia,j,jj,ja
!   integer  :: iref,jref
!   integer  :: k,kk,l,ll
!
!   real(wp),allocatable :: c(:),r(:),c6(:),c6_0(:)
!
!   allocate( c(nat), r(nat), c6(nat), c6(nat*(1+nat)/2), source = 0.0_wp )
!   c = 0.0_wp
!   r = 0.0_wp
!   c6 = 0.0_wp
!   c6_0 = 0.0_wp
!
!   ! precalculate all information which not explicitly depend on the grid point
!   k = 0
!   do i = 1, nat
!      ia = at(i)
!      ii = i*(i-1)/2
!      c(ii) = 3._wp*r4r2(ia)*r4r2(iprobe)
!      c6tmp = 0.0_wp
!      c6tmp0 = 0.0_wp
!      do iref = 1, refn(ia)
!         k = k+1
!         itbl(iref,i) = k
!         alpha = refal(:,iref,ia)*alpha_probe
!         c6tmp = c6tmp + zetvec(k)*thopi*trapzd(alpha)
!         c6tmp0 = c6tmp0 + zerovec(k)*thopi*trapzd(alpha)
!      enddo
!      c6(i) = c6tmp
!      c6_0(ii) = c6tmp0
!      ! get all pairs for the ATM calculation
!      do j = 1, i-1
!         ja = at(j)
!         ij = ii + j
!         r4r2ij = 3._wp*r4r2(ia)*r4r2(ja)
!         r(ij) = norm2(xyz(:,i)-xyz(:,j))
!         c(ij) = par%a1*sqrt(r4r2ij)+par%a2
!         c6tmp = 0.0_wp
!         do iref = 1, refn(ia)
!            kk = itbl(iref,i)
!            do jref = 1, refn(ja)
!               ll = itbl(jref,j)
!               alpha = refal(:,iref,ia)*refal(:,jref,ja)
!               c6tmp = c6tmp + zerovec(kk)*zerovec(ll)*thopi*trapzd(alpha)
!            enddo
!         enddo
!         c6_0(ij) = c6tmp
!      enddo
!   enddo
!
!   do k = 1, surface%n
!      do i = 1, nat
!         ii = i*(i-1)/2
!         ia = at(i)
!         rik = r(ii)
!         r4r2ij = 3.0_wp * r4r2(ia)*r4r2(iprobe)
!         cik = par%a1 * sqrt(r4r2ij) + par%a2
!         r6  = fdmpr_bj( 6,rik,cik)
!         r8  = fdmpr_bj( 8,rik,cik)
!         r10 = fdmpr_bj(10,rik,cik)
!         ! pairwise terms
!         Epoint = c6(i) * (par%s6*r6 + par%s8*r4r2ij*r8 &
!         &                           + par%s10*r4r2ij**2*49.0_wp/40.0_wp*r10 )
!         ! nonadditive terms
!         do j = 1, i-1
!            ! get indices
!            ij = ii + j
!            jj = j*(j-1)/2
!            ! store in temps
!            r2ij = r(ij)**2
!            r2jk = r(jj)
!            cij = c(ij)
!            cjk = c(jj)
!            atm = ( 0.375_wp * (r2ij+r2jk-r2ik) &
!            &                * (r2ij+r2ik-r2jk) &
!            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp
!            fdmp = one/(one+six*((cijk/rijk)**oth)**par%alp)
!            oor9ijk = atm/rijk**3*fdmp
!            E = E + c9ijk * oor9ijk
!         enddo
!      enddo
!   enddo
!
!end subroutine dispersion_surfaceplot
