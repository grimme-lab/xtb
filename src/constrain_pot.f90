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
module xtb_constrainpot
contains

subroutine constrain_zaxis(fix,n,at,xyz,g,e)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_setvar
   implicit none
   type(fix_setvar),intent(in) :: fix
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: e
   real(wp),intent(inout) :: g(3,n)

   integer  :: i,ii
   real(wp) :: z0,dz

   if (fix%n.eq.0) return

   do i = 1, fix%n
      ii = fix%atoms(i)
      z0 = fix%val(i)
      dz = xyz(3,ii) - z0
      e  = e + fix%fc*dz**2
      g(3,ii) = g(3,ii) + 2.0d0*fix%fc*dz
   enddo

end subroutine constrain_zaxis

subroutine constrain_pos(fix,n,at,xyz,g,e)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_setvar
   implicit none
   type(fix_setvar),intent(in) :: fix
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: e
   real(wp),intent(inout) :: g(3,n)

   integer  :: i,j,ii,jj,k
   real(wp) :: r,rij(3),tmp,d

   if (fix%n.eq.0) return

   k = 0
   do i = 1, fix%n
      ii = fix%atoms(i)
      do j = 1, i-1
         jj = fix%atoms(j)
         k = k+1
         rij = xyz(:,jj) - xyz(:,ii)
         r = norm2(rij)
         d = r - fix%val(k)
         e = e + fix%fc * d**2
         tmp = fix%fc * 2.0_wp * d
         g(:,jj) = g(:,jj) + tmp*rij/r
         g(:,ii) = g(:,ii) - tmp*rij/r
      enddo
   enddo
contains
pure elemental integer function lin(i1,i2)
   integer,intent(in) :: i1,i2
   integer :: idum1,idum2
   idum1=max(i1,i2)
   idum2=min(i1,i2)
   lin=idum2+idum1*(idum1-1)/2
end function lin

end subroutine constrain_pos

subroutine qpothess2(fix,n,at,xyz,h)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_setvar
   implicit none
   type(fix_setvar),intent(in) :: fix
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: h(n*3*(n*3+1)/2)

   real(wp) :: rij(3),dr,r0
   real(wp) :: t1,t0,w0,w1,dx,dx2,dy,dy2,r,r2,r3
   integer  :: i,j,k,l,m,ii,jj,iik,iikl,ikjm,ijk,ij

   if(fix%n.eq.0) return

   ! calculate hessian
   do i = 1, fix%n !loop over all atoms
      ii = (fix%atoms(i)-1)*3
      do j = 1, fix%n !inner loop for sum over all atoms
         if (i.ne.j) then
            rij = xyz(:,fix%atoms(j))-xyz(:,fix%atoms(i))
            r  = norm2(rij)
            r2 = r*r
            r3 = r2*r
            ij = lin(i-1,j-1)+1
            r0 = fix%val(ij)
            dr = r-r0
            do k = 1, 3  !loop diagonal elements (xyz components)
               dx = xyz(k,fix%atoms(i))-xyz(k,fix%atoms(j))
               dx2 = dx*dx
               iik = lin(ii+k,ii+k)
               h(iik) = h(iik) + 2.0_wp*fix%fc*(1.0_wp+(dx2/r2)-(dx2*dr/r3)-(r0/r))
               do l = k+1, 3 !loop same-atom block-diagonal elements
                  dy = xyz(l,fix%atoms(i))-xyz(l,fix%atoms(j))
                  iikl = lin(ii+k,ii+l)
                  h(iikl) = h(iikl) + 2.0_wp*fix%fc*r0*dx*dy/r3
               enddo !end loop same-atom block-diagonal elements
            enddo !end loop diagonal elements (xyz components)
         endif
      enddo !end inner loop for sum over all atoms
      do j = i+1, fix%n !loop over the rest (mixed atoms)
         rij = xyz(:,fix%atoms(j))-xyz(:,fix%atoms(i))
         r  = norm2(rij)
         r2 = r*r
         r3 = r2*r
         ij = lin(i-1,j-1)+1
         r0 = fix%val(ij)
         dr = r-r0
         jj = (fix%atoms(j)-1)*3
         do k = 1, 3  !loop diagonal elements (xyz components)
            dx = xyz(k,fix%atoms(i))-xyz(k,fix%atoms(j))
            do m = 1, 3
               if (k.eq.m) then !same component case
                  dx2 = dx*dx
                  ijk = lin(ii+k,jj+k)
                  h(ijk) = h(ijk) - 2.0_wp*fix%fc*(1.0_wp+dx2/r2-dx2*dr/r3-r0/r)
               else !different component case
                  dy = xyz(m,fix%atoms(i))-xyz(m,fix%atoms(j))
                  ikjm = lin(ii+k,jj+m)
                  h(ikjm) = h(ikjm) - 2.0_wp*fix%fc*r0*dx*dy/r3
               endif
            enddo
         enddo !end loop diagonal elements (xyz components)
       enddo
   enddo !end loop atoms

contains
pure elemental integer function lin(i1,i2)
   integer,intent(in) :: i1,i2
   integer :: idum1,idum2
   idum1=max(i1,i2)
   idum2=min(i1,i2)
   lin=idum2+idum1*(idum1-1)/2
end function lin

end subroutine qpothess2

subroutine constrain_dist(fix,n,at,xyz,g,e)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_setvar
   implicit none
   type(fix_setvar),intent(in) :: fix
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: e
   real(wp),intent(inout) :: g(3,n)

   integer i,j,k,l,m,mm
   real(wp)rij(3),dum,r0,r,vp(3),ra(3),rb(3)
   real(wp)va(3),vb(3),vc(3),vab(3),vcb(3),deda(3),dedc(3),dedb(3)
   real(wp)dda(3),ddb(3),ddc(3),ddd(3)
   real(wp)c0,rab2,rcb2,theta,dt,rmul1,rmul2,deddt,rp,cosa,ea,vlen,d
   real(wp)c1,phi0,phi,dphi1,x1cos,x1sin,dij,valijkl,ff,dum1,dum2
   real(wp)dx,dy,dz,termx,termy,termz,abcx,abcy,abcz,a,z0,expo,expom1

   if(fix%n.eq.0) return

   do m = 1, fix%n
      mm = 2*m-1
      i = fix%atoms(mm)
      j = fix%atoms(mm+1)
      if (i == j) cycle ! workaround for ifx 2024.1.0
      r0= fix%val(m)
      rij=xyz(:,j)-xyz(:,i)
      r = norm2(rij)
      d=r-r0
!     e=e+fix%fc*d*d
!     ff=fix%fc*2.0d0*d
      dum= d**fix%expo(m)
      dum2=d**(fix%expo(m)-1.0_wp)
      e=e+fix%fc*dum
      ff=fix%fc*fix%expo(m)*dum2
      dum=ff/r
      g(:,j)=g(:,j)+dum*rij
      g(:,i)=g(:,i)-dum*rij

   enddo

end subroutine constrain_dist

subroutine constrain_angle(fix,n,at,xyz,g,e)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_setvar
   use xtb_basic_geo
   implicit none
   type(fix_setvar),intent(in) :: fix
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: e
   real(wp),intent(inout) :: g(3,n)

   integer i,j,k,l,m,mm
   real(wp)rij(3),dum,r0,ffc,vp(3),ra(3),rb(3)
   real(wp)va(3),vb(3),vc(3),vab(3),vcb(3),deda(3),dedc(3),dedb(3)
   real(wp)dda(3),ddb(3),ddc(3),ddd(3)
   real(wp)c0,rab2,rcb2,theta,dt,rmul1,rmul2,deddt,rp,cosa,ea,vlen,d
   real(wp)dx,dy,dz,termx,termy,termz,abcx,abcy,abcz,a,z0,expo,expom1

   if(fix%n.eq.0) return

   do m = 1, fix%n
      mm = 3*m-2
      i = fix%atoms(mm)
      j = fix%atoms(mm+1)
      k = fix%atoms(mm+2)
      c0 = fix%val(m)
      va = xyz(1:3,i)
      vb = xyz(1:3,j)
      vc = xyz(1:3,k)
      vab = va-vb
      vcb = vc-vb
      rab2 = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
      rcb2 = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
      call crprod(vcb,vab,vp)
      rp = norm2(vp)+1.d-14
      call impsc(vab,vcb,cosa)
      cosa = dble(min(1.0d0,max(-1.0d0,cosa)))
      theta= dacos(cosa)
      dt  = theta - c0
      ea  = fix%fc* dt**2
      deddt = 2.d0 * fix%fc * dt
      e = e + ea
      call crprod(vab,vp,deda)
      rmul1 = -deddt / (rab2*rp)
      deda = deda*rmul1
      call crprod(vcb,vp,dedc)
      rmul2 =  deddt / (rcb2*rp)
      dedc = dedc*rmul2
      dedb = deda+dedc
      g(:,i) = g(:,i) + deda
      g(:,j) = g(:,j) - dedb
      g(:,k) = g(:,k) + dedc
   enddo

end subroutine constrain_angle

subroutine constrain_dihedral(fix,n,at,xyz,g,e)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_type_setvar
   use xtb_basic_geo
   implicit none
   type(fix_setvar),intent(in) :: fix
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: e
   real(wp),intent(inout) :: g(3,n)

   integer i,j,k,l,m,mm
   real(wp)rij(3),dum,r0,ffc,vp(3),ra(3),rb(3)
   real(wp)va(3),vb(3),vc(3),vab(3),vcb(3),deda(3),dedc(3),dedb(3)
   real(wp)dda(3),ddb(3),ddc(3),ddd(3)
   real(wp)c0,rab2,rcb2,theta,dt,rmul1,rmul2,deddt,rp,cosa,ea,vlen,d
   real(wp)c1,phi0,phi,dphi1,x1cos,x1sin,dij,ff,dum1,dum2
   real(wp)dx,dy,dz,termx,termy,termz,abcx,abcy,abcz,a,z0,expo,expom1

   if(fix%n.eq.0) return

   do m = 1, fix%n
      mm = 4*m-3
      i= fix%atoms(mm)
      j= fix%atoms(mm+1)
      k= fix%atoms(mm+2)
      l= fix%atoms(mm+3)
      phi0 =fix%val(m)
      phi=valijkl(n,xyz,i,j,k,l)
      if(abs(phi-pi).lt.1.d-8.or.abs(phi).lt.1.d-8) phi=phi+1.d-8
      call dphidr(n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)
      dphi1=phi0-phi+pi
      x1cos=cos(dphi1)
      x1sin=sin(dphi1)
      dij  =fix%fc*x1sin
      g(:,i)=g(:,i)+dij*dda
      g(:,j)=g(:,j)+dij*ddb
      g(:,k)=g(:,k)+dij*ddc
      g(:,l)=g(:,l)+dij*ddd
      e=e+fix%fc*(1.0_wp+x1cos)
   enddo

end subroutine constrain_dihedral

subroutine constrain_pot(fix,n,at,xyz,g,e)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_setvar
   implicit none
   type(constr_setvar),intent(in) :: fix
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz(3,n)
   real(wp),intent(inout) :: e
   real(wp),intent(inout) :: g(3,n)

   call constrain_pos     (fix%pos,     n,at,xyz,g,e)
   call constrain_dist    (fix%dist,    n,at,xyz,g,e)
   call constrain_angle   (fix%angle,   n,at,xyz,g,e)
   call constrain_dihedral(fix%dihedral,n,at,xyz,g,e)

end subroutine constrain_pot

subroutine constrain_hess(fix,n,at,xyz0,Hess)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_setvar
   implicit none
   type(constr_setvar),intent(in) :: fix
   integer, intent(in)    :: n
   integer, intent(in)    :: at(n)
   real(wp),intent(in)    :: xyz0(3,n)
   real(wp),intent(inout) :: Hess((n*3)*((n*3)+1)/2)

   integer ia,ic,ja,jc,i,j,k,n3,ii,jj
   real(wp) step,step2,e
   real(wp),allocatable :: xyz(:,:)
   real(wp),allocatable :: h(:,:),gr(:,:),gl(:,:)

   if(fix%n.eq.0) return

   allocate( h(n*3,n*3),gr(3,n),gl(3,n),xyz(3,n), source = 0.0_wp )

   xyz = xyz0

   n3=3*n
   step =0.0001_wp
   step2=0.5_wp/step
   e=0

   do ia = 1, n
      do ic = 1, 3
         ii = (ia-1)*3+ic
         xyz(ic,ia)=xyz0(ic,ia)+step
         gr=0.0_wp
         call constrpot(n,at,xyz,gr,e)
         xyz(ic,ia)=xyz0(ic,ia)-step
         gl=0.0_wp
         call constrpot(n,at,xyz,gl,e)
         xyz(ic,ia)=xyz0(ic,ia)
         do ja = 1, n
            do jc = 1, 3
               jj = (ja-1)*3 + jc
               h(jj,ii) =(gr(jc,ja) - gl(jc,ja)) * step2
            enddo
         enddo
      enddo
   enddo

!  call prmat(6,h,3*n,3*n,'H')

   k=0
   do i=1,n3
      do j=1,i
         k=k+1
         hess(k)=hess(k)+0.5_wp*(h(j,i)+h(i,j))
      enddo
   enddo

end subroutine constrain_hess

end module xtb_constrainpot
