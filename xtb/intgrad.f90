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

! --------------------------------------------------------------[SAW1907]-
!> Backend implementation for all overlap distribution based integrals,
!  we use a hardcoded horizontal Obara--Saika recursion relation to get
!  the job done, this is working code, so think twice before modifying it!
module intgrad
   use iso_fortran_env, only : wp => real64
   implicit none

   integer,private,parameter :: lx(84) = (/ &
      & 0, &
      & 1,0,0, &
      & 2,0,0,1,1,0, &
      & 3,0,0,2,2,1,0,1,0,1, &
      & 4,0,0,3,3,1,0,1,0,2,2,0,2,1,1, &
      & 5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1, &
      & 6,0,0,3,3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/)
   integer,private,parameter :: ly(84) = (/ &
      & 0, &
      & 0,1,0, &
      & 0,2,0,1,0,1, &
      & 0,3,0,1,0,2,2,0,1,1, &
      & 0,4,0,1,0,3,3,0,1,2,0,2,1,2,1, &
      & 0,5,0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2, &
      & 0,6,0,3,0,3,1,0,0,1,5,5,2,0,0,2,4,4,2,1,3,1,3,2,1,4,1,2/)
   integer,private,parameter :: lz(84) = (/ &
      & 0, &
      & 0,0,1, &
      & 0,0,2,0,1,1, &
      & 0,0,3,0,1,0,1,2,2,1, &
      & 0,0,4,0,1,0,1,3,3,0,2,2,1,1,2, &
      & 0,0,5,0,2,0,3,2,3,0,1,0,1,4,4,3,1,1,1,2,2, &
      & 0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,1,1,4,2/)
   !     integer,parameter :: iall(4,4) = reshape( (/
   !     1,4,10,20,4,10,20,35,10,20,35,56,20,35,56,84 /), shape(iall) )
   !     s    px   py   pz    dx²   dy²   dz²   dxy   dxz   dyz
   !     1    2    3    4     5     6     7     8     9     10
   !     fx³  fy³  fz³  fx²y  fx²z  fy²x  fy²z  fxz²  fyz²  fxyz
   !     11   12   13   14    15    16    17    18    19    20
   !     gx⁴  gy⁴  gz⁴  gx³y  gx³z  gy³x  gy³z  gz³x  gz³y  gx²y²
   !     21   22   23   24    25    26    27    28    29    30
   !     gx²z²     gy²z²      gx²yz       gy²xz       gz²xy
   !     31        32         33          34          35


contains

! --------------------------------------------------------------[SAW1907]-
!> calculates a partial overlap in one cartesian direction
pure elemental function olapp(l,gama) result(s)
   implicit none
   integer,intent(in) :: l
   real(wp),intent(in) :: gama
   real(wp) :: s
   real(wp) :: gm
   integer  :: lh
   real(wp),parameter :: dftr(0:7) = & ! see OEIS A001147
      & [1._wp,1._wp,3._wp,15._wp,105._wp,945._wp,10395._wp,135135._wp]
   if (mod(l,2).ne.0) then
      s=0._wp
   else
      lh=l/2
      gm=0.5_wp/gama
      s=gm**lh*dftr(lh)
   endif
end function olapp

! --------------------------------------------------------------[SAW1907]-
!> returns center of product Gaussian from two Gaussians by GPT
pure function gpcenter(alp,ra,bet,rb) result(rp)
   implicit none
   real(wp),intent(in) :: alp,bet
   real(wp),intent(in) :: ra(3),rb(3)
   real(wp) :: rp(3)

   rp = (alp*ra + bet*rb)/(alp+bet)

end function gpcenter

! --------------------------------------------------------------[SAW1907]-
pure subroutine build_kab(ra,alp,rb,bet,gama,kab)
   use mctc_constants
   implicit none
   !     this computes the center, exponent, and multiplying factor of
   !     a single gaussian which can replace the product of two gaussian
   !     centers a and b, and exponents alpha and beta.
   real(wp),intent(in)  :: ra(3),alp,rb(3),bet
   real(wp),intent(out) :: gama,kab
   real(wp) :: rab2,rab(3),est,gm
   gama  = alp+bet
   gm    = 1.0_wp/gama
   rab   = ra - rb
   rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
   est   = rab2*alp*bet*gm
   kab   = exp(-est)*(sqrtpi*sqrt(gm))**3
end subroutine build_kab

! --------------------------------------------------------------[SAW1712]-
!> not needed anymore, but the math was lengthy, so I rather keep it
pure subroutine rhftce2(cfs,a,e,iff)
   implicit none
   integer, intent(in)  :: iff
   real(wp),intent(in)  :: a(*),e(*)
   real(wp),intent(inout) :: cfs(*)
   real(wp),parameter   :: c2 = 2.0_wp
   real(wp),parameter   :: c3 = 3.0_wp
   real(wp)  :: aex,aey,aez
   ! ---- e = center of product function, a = center of single gaussian
   aex = e(1)-a(1)
   aey = e(2)-a(2)
   aez = e(3)-a(3)
   select case(iff)
   case(1) ! s
      continue
   case(2) ! x·s + px
      cfs( 1)=aex*cfs(2)
   case(3) ! y·s + py
      cfs( 1)=aey*cfs(3)
   case(4) ! z·s + pz
      cfs( 1)=aez*cfs(4)
   case(5) ! x²·s + 2x·px + dx²
      cfs( 1)=aex*aex*cfs(5)
      cfs( 2)=c2*aex*cfs(5)
   case(6) ! y²·s + 2y·py + dy²
      cfs( 1)=aey*aey*cfs(6)
      cfs( 3)=c2*aey*cfs(6)
   case(7) ! z²·s + 2z·pz + dz²
      cfs( 1)=aez*aez*cfs(7)
      cfs( 4)=c2*aez*cfs(7)
   case(8) ! xy·s + y·px + x·py + dxy
      cfs( 1)=aex*aey*cfs(8)
      cfs( 2)=aey*cfs(8)
      cfs( 3)=aex*cfs(8)
   case(9) ! xz·s + z·px + x·pz + dxz
      cfs( 1)=aex*aez*cfs(9)
      cfs( 2)=aez*cfs(9)
      cfs( 4)=aex*cfs(9)
   case(10) ! yz·s + z·py + y·pz + dyz
      cfs( 1)=aey*aez*cfs(10)
      cfs( 3)=aez*cfs(10)
      cfs( 4)=aey*cfs(10)
   case(11) ! x³·s + 3x²·px + 3x·dx² + fx³
      cfs( 1)=aex*aex*aex*cfs(11)
      cfs( 2)=c3*aex*aex*cfs(11)
      cfs( 5)=c3*aex*cfs(11)
   case(12) ! y³·s + 3y²·py + 3y·dy² + fy³
      cfs( 1)=aey*aey*aey*cfs(12)
      cfs( 3)=c3*aey*aey*cfs(12)
      cfs( 6)=c3*aey*cfs(12)
   case(13) ! z³·s + 3z²·pz + 3z·dz² + fz³
      cfs( 1)=aez*aez*aez*cfs(13)
      cfs( 4)=c3*aez*aez*cfs(13)
      cfs( 7)=c3*aez*cfs(13)
   case(14) ! x²y·s + 2xy·px + x²·py + x·dx² + 2x·dxy + fx²y
      cfs( 1)=aex*aex*aey*cfs(14)
      cfs( 2)=c2*aex*aey*cfs(14)
      cfs( 3)=aex*aex*cfs(14)
      cfs( 5)=aey*cfs(14)
      cfs( 8)=c2*aex*cfs(14)
   case(15) ! x²z·s + 2xz·px + x²·pz + z·dx² + 2x·dxz + fx²z
      cfs( 1)=aex*aex*aez*cfs(15)
      cfs( 2)=c2*aex*aez*cfs(15)
      cfs( 4)=aex*aex*cfs(15)
      cfs( 5)=aez*cfs(15)
      cfs( 9)=c2*aex*cfs(15)
   case(16) ! y²x·s + y²·px + 2xy·py + x·dy² + 2y·dxy + fy²x
      cfs( 1)=aey*aey*aex*cfs(16)
      cfs( 2)=aey*aey*cfs(16)
      cfs( 3)=c2*aey*aex*cfs(16)
      cfs( 6)=aex*cfs(16)
      cfs( 8)=c2*aey*cfs(16)
   case(17) ! y²z·s + 2yz·py + y²·pz + z·dy² + 2y·dyz + fy²z
      cfs( 1)=aey*aey*aez*cfs(17)
      cfs( 3)=c2*aey*aez*cfs(17)
      cfs( 4)=aey*aey*cfs(17)
      cfs( 6)=aez*cfs(17)
      cfs(10)=c2*aey*cfs(17)
   case(18) ! xz²·s + z²·px + 2xz·pz + x·dz² + 2z·dxz + fxz²
      cfs( 1)=aez*aez*aex*cfs(18)
      cfs( 2)=aez*aez*cfs(18)
      cfs( 4)=c2*aez*aex*cfs(18)
      cfs( 7)=aex*cfs(18)
      cfs( 9)=c2*aez*cfs(18)
   case(19) ! yz²·s + z²·py + 2yz·pz + y·dz² + 2z·dyz + fyz²
      cfs( 1)=aez*aez*aey*cfs(19)
      cfs( 3)=aez*aez*cfs(19)
      cfs( 4)=c2*aez*aey*cfs(19)
      cfs( 7)=aey*cfs(19)
      cfs(10)=c2*aez*cfs(19)
   case(20) ! xyz·s + yz·px + xz·py + xy·pz + z·dxy + y·dxz + x·dyz + fxyz
      cfs( 1)=aex*aey*aez*cfs(20)
      cfs( 2)=aez*aey*cfs(20)
      cfs( 3)=aex*aez*cfs(20)
      cfs( 4)=aex*aey*cfs(20)
      cfs( 8)=aez*cfs(20)
      cfs( 9)=aey*cfs(20)
      cfs(10)=aex*cfs(20)
   end select
   return
end subroutine rhftce2

! --------------------------------------------------------------[SAW1801]-
pure subroutine dtrf2(s,li,lj)
   use mctc_la, only : gemm
   implicit none
   real(wp),intent(inout) :: s(6,6)
   integer, intent(in)    :: li,lj
   ! CAO-AO transformation
   real(wp) :: trafo(6,6)
   parameter (trafo = reshape((/ & ! copied from scf.f, simplyfied
      ! --- dS
      & sqrt(1.0_wp/5.0_wp), &
      & sqrt(1.0_wp/5.0_wp), &
      & sqrt(1.0_wp/5.0_wp), &
      & 0.0_wp,0.0_wp,0.0_wp, &
      ! --- dx²-y²
      & 0.5_wp*sqrt(3.0_wp), &
      &-0.5_wp*sqrt(3.0_wp), &
      & 0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
      ! --- dz²
      & 0.5_wp,0.5_wp,-1.0_wp, &
      & 0.0_wp,0.0_wp, 0.0_wp, &
      ! --- rest
      & 0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp,0.0_wp, &
      & 0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp, &
      & 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp /), shape(trafo)))

   real(wp) s2(6,6),sspher,dum(6,6)
   integer ii,jj,m,n

   !     transformation not needed for pure s/p overlap -> do nothing
   if (li.lt.2.and.lj.lt.2) return

   ! --- means li.ge.2.or.lj.ge.2, so one of them is a d-shell
   !     assuming its on jat ... a wild guess
   select case(li)
   case(0) ! s-d
      do jj=1,6
         sspher=0
         do m=1,6
            sspher=sspher+trafo(m,jj)*s(m,1)
         enddo
         s2(jj,1)=sspher
      enddo
      s(1:5,1) = s2(2:6,1)
      return
   case(1) ! p-d
      do ii=1,3
         do jj=1,6
            sspher=0
            do m=1,6
               sspher=sspher+trafo(m,jj)*s(m,ii)
            enddo
            s2(jj,ii)=sspher
         enddo
         s(1:5,ii) = s2(2:6,ii)
      enddo
      return
   end select
   !     wasn't there, then try iat ...
   select case(lj)
   case(0) ! d-s
      do jj=1,6
         sspher=0
         do m=1,6
            sspher=sspher+trafo(m,jj)*s(1,m)
         enddo
         s2(1,jj)=sspher
      enddo
      s(1,1:5) = s2(1,2:6)
      return
   case(1) ! d-p
      do ii=1,3
         do jj=1,6
            sspher=0
            do m=1,6
               sspher=sspher+trafo(m,jj)*s(ii,m)
            enddo
            s2(ii,jj)=sspher
         enddo
         s(ii,1:5) = s2(ii,2:6)
      enddo
      return
   end select
   !     if not returned up to here -> d-d
   ! CB: transposing s in first dgemm is important for integrals other than S
   CALL GEMM('T','N',6,6,6,1._wp,s,6,trafo,6,0._wp,dum,6)
   CALL GEMM('T','N',6,6,6,1._wp,dum,6,trafo,6,0._wp,s2,6)
   s(1:5,1:5) = s2(2:6,2:6)
   return

end subroutine dtrf2

! --------------------------------------------------------------[SAW1801]-
pure subroutine build_hshift(cfs,a,e,l)
   implicit none
   integer,intent(in)  :: l(3)
   real(wp), intent(in)  :: a(3),e(3)
   real(wp), intent(inout) :: cfs(3,*)
   real(wp), parameter   :: c2 = 2.0_wp
   real(wp), parameter   :: c3 = 3.0_wp
   real(wp), parameter   :: c4 = 4.0_wp
   real(wp), parameter   :: c6 = 6.0_wp
   integer :: i
   real(wp)  :: ae
   ! --- e = center of product function, a = center of single gaussian
   do i = 1, 3
      ae = e(i)-a(i)
      select case(l(i))
      case(0) ! s
         continue
      case(1) ! p
         cfs(i,1)=ae*cfs(i,2)
      case(2) ! d
         cfs(i,1)=ae*ae*cfs(i,3)
         cfs(i,2)=c2*ae*cfs(i,3)
      case(3) ! f
         cfs(i,1)=ae*ae*ae*cfs(i,4)
         cfs(i,2)=c3*ae*ae*cfs(i,4)
         cfs(i,3)=c3*ae*cfs(i,4)
      case(4) ! g
         cfs(i,1)=ae*ae*ae*ae*cfs(i,5)
         cfs(i,2)=c4*ae*ae*ae*cfs(i,5)
         cfs(i,3)=c6*ae*ae*cfs(i,5)
         cfs(i,4)=c4*ae*cfs(i,5)
      end select
   enddo
end subroutine build_hshift

! --------------------------------------------------------------[SAW1801]-
pure subroutine build_hshift2(cfs,a,e,l)
   implicit none
   integer,intent(in)  :: l
   real(wp), intent(in)  :: a,e
   real(wp), intent(inout) :: cfs(*)
   integer :: i
   real(wp)  :: ae
   ! --- e = center of product function, a = center of single gaussian
   ae = e-a
   select case(l)
   case(0) ! s
      continue
   case(1) ! p
      cfs(1)=ae*cfs(2)
   case(2) ! d
      cfs(1)=ae*ae*cfs(3)
      cfs(2)= 2*ae*cfs(3)
   case(3) ! f
      cfs(1)=ae*ae*ae*cfs(4)
      cfs(2)= 3*ae*ae*cfs(4)
      cfs(3)= 3*ae*cfs(4)
   case(4) ! g
      cfs(1)=ae*ae*ae*ae*cfs(5)
      cfs(2)= 4*ae*ae*ae*cfs(5)
      cfs(3)= 6*ae*ae*cfs(5)
      cfs(4)= 4*ae*cfs(5)
   end select
end subroutine build_hshift2

! --------------------------------------------------------------[SAW1801]-
pure subroutine prod3(a,b,d,la,lb)
   implicit none
   integer,intent(in)    :: la,lb
   real(wp), intent(in)    :: a(*),b(*)
   real(wp), intent(inout) :: d(*)
   integer :: i
   if(la.ge.4.or.lb.ge.4) goto 40
   if(la.ge.3.or.lb.ge.3) goto 30
   if(la.ge.2.or.lb.ge.2) goto 20
   ! <s|s> = <s>
   d(1)=a(1)*b(1)
   if(la.eq.0.and.lb.eq.0) return
   ! <s|p> = <s|*(|s>+|p>)
   !       = <s> + <p>
   d(2)=a(1)*b(2)+a(2)*b(1)
   if(la.eq.0.or.lb.eq.0) return
   ! <p|p> = (<s|+<p|)*(|s>+|p>)
   !       = <s> + <p> + <d>
   d(3)=a(2)*b(2)
   return
20 continue
   ! <s|d> = <s|*(|s>+|p>+|d>)
   !       = <s> + <p> + <d>
   d(1)=a(1)*b(1)
   d(2)=a(1)*b(2)+a(2)*b(1)
   d(3)=a(1)*b(3)+a(3)*b(1)
   if(la.eq.0.or.lb.eq.0) return
   ! <p|d> = (<s|+<p|)*(|s>+|p>+|d>)
   !       = <s> + <p> + <d> + <f>
   d(3)=d(3)+a(2)*b(2)
   d(4)=a(2)*b(3)+a(3)*b(2)
   if(la.le.1.or.lb.le.1) return
   ! <d|d> = (<s|+<p|+<d|)*(|s>+|p>+|d>)
   !       = <s> + <p> + <d> + <f> + <g>
   d(5)=a(3)*b(3)
   return
30 continue
   ! <s|f> = <s|*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f>
   d(1)=a(1)*b(1)
   d(2)=a(1)*b(2)+a(2)*b(1)
   d(3)=a(1)*b(3)+a(3)*b(1)
   d(4)=a(1)*b(4)+a(4)*b(1)
   if(la.eq.0.or.lb.eq.0) return
   ! <p|f> = (<s|+<p|)*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f> + <g>
   d(3)=d(3)+a(2)*b(2)
   d(4)=d(4)+a(2)*b(3)+a(3)*b(2)
   d(5)=a(2)*b(4)+a(4)*b(2)
   if(la.le.1.or.lb.le.1) return
   ! <d|f> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f> + <g> + <h>
   d(5)=d(5)+a(3)*b(3)
   d(6)=a(3)*b(4)+a(4)*b(3)
   if(la.le.2.or.lb.le.2) return
   ! <f|f> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
   d(7)=a(4)*b(4)
   return
40 continue
   ! <s|g> = <s|*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g>
   d(1)=a(1)*b(1)
   d(2)=a(1)*b(2)+a(2)*b(1)
   d(3)=a(1)*b(3)+a(3)*b(1)
   d(4)=a(1)*b(4)+a(4)*b(1)
   d(5)=a(1)*b(5)+a(5)*b(1)
   if(la.eq.0.or.lb.eq.0) return
   ! <p|g> = (<s|+<p|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h>
   d(3)=d(3)+a(2)*b(2)
   d(4)=d(4)+a(2)*b(3)+a(3)*b(2)
   d(5)=d(5)+a(2)*b(4)+a(4)*b(2)
   d(6)=a(2)*b(5)+a(5)*b(2)
   if(la.le.1.or.lb.le.1) return
   ! <d|g> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
   d(5)=d(5)+a(3)*b(3)
   d(6)=d(5)+a(3)*b(4)+a(4)*b(3)
   d(7)=a(3)*b(5)+a(5)*b(3)
   if(la.le.2.or.lb.le.2) return
   ! <f|g> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k>
   d(7)=d(7)+a(4)*b(4)
   d(8)=a(4)*b(5)+a(5)*b(4)
   if(la.le.3.or.lb.le.3) return
   ! <g|g> = (<s|+<p|+<d|+<f|+<g|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k> + <l>
   d(9)=a(5)*b(5)

end subroutine prod3

! --------------------------------------------------------------[SAW1801]-
pure subroutine prod2(a,b,d,la,lb)
   implicit none
   integer,intent(in)    :: la(3),lb(3)
   real(wp), intent(in)    :: a(3,*),b(3,*)
   real(wp), intent(inout) :: d(3,*)
   integer :: i
   do i = 1, 3
      if(la(i).ge.4.or.lb(i).ge.4) goto 40
      if(la(i).ge.3.or.lb(i).ge.3) goto 30
      if(la(i).ge.2.or.lb(i).ge.2) goto 20
      ! <s|s> = <s>
      d(i,1)=a(i,1)*b(i,1)
      if(la(i).eq.0.and.lb(i).eq.0) cycle
      ! <s|p> = <s|*(|s>+|p>)
      !       = <s> + <p>
      d(i,2)=a(i,1)*b(i,2)+a(i,2)*b(i,1)
      if(la(i).eq.0.or.lb(i).eq.0) cycle
      ! <p|p> = (<s|+<p|)*(|s>+|p>)
      !       = <s> + <p> + <d>
      d(i,3)=a(i,2)*b(i,2)
      cycle
   20 continue
      ! <s|d> = <s|*(|s>+|p>+|d>)
      !       = <s> + <p> + <d>
      d(i,1)=a(i,1)*b(i,1)
      d(i,2)=a(i,1)*b(i,2)+a(i,2)*b(i,1)
      d(i,3)=a(i,1)*b(i,3)+a(i,3)*b(i,1)
      if(la(i).eq.0.or.lb(i).eq.0) cycle
      ! <p|d> = (<s|+<p|)*(|s>+|p>+|d>)
      !       = <s> + <p> + <d> + <f>
      d(i,3)=d(i,3)+a(i,2)*b(i,2)
      d(i,4)=a(i,2)*b(i,3)+a(i,3)*b(i,2)
      if(la(i).le.1.or.lb(i).le.1) cycle
      ! <d|d> = (<s|+<p|+<d|)*(|s>+|p>+|d>)
      !       = <s> + <p> + <d> + <f> + <g>
      d(i,5)=a(i,3)*b(i,3)
      cycle
   30 continue
      ! <s|f> = <s|*(|s>+|p>+|d>+|f>)
      !       = <s> + <p> + <d> + <f>
      d(i,1)=a(i,1)*b(i,1)
      d(i,2)=a(i,1)*b(i,2)+a(i,2)*b(i,1)
      d(i,3)=a(i,1)*b(i,3)+a(i,3)*b(i,1)
      d(i,4)=a(i,1)*b(i,4)+a(i,4)*b(i,1)
      if(la(i).eq.0.or.lb(i).eq.0) cycle
      ! <p|f> = (<s|+<p|)*(|s>+|p>+|d>+|f>)
      !       = <s> + <p> + <d> + <f> + <g>
      d(i,3)=d(i,3)+a(i,2)*b(i,2)
      d(i,4)=d(i,4)+a(i,2)*b(i,3)+a(i,3)*b(i,2)
      d(i,5)=a(i,2)*b(i,4)+a(i,4)*b(i,2)
      if(la(i).le.1.or.lb(i).le.1) cycle
      ! <d|f> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>)
      !       = <s> + <p> + <d> + <f> + <g> + <h>
      d(i,5)=d(i,5)+a(i,3)*b(i,3)
      d(i,6)=a(i,3)*b(i,4)+a(i,4)*b(i,3)
      if(la(i).le.2.or.lb(i).le.2) cycle
      ! <f|f> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>)
      !       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
      d(i,7)=a(i,4)*b(i,4)
      cycle
   40 continue
      ! <s|g> = <s|*(|s>+|p>+|d>+|f>+|g>)
      !       = <s> + <p> + <d> + <f> + <g>
      d(i,1)=a(i,1)*b(i,1)
      d(i,2)=a(i,1)*b(i,2)+a(i,2)*b(i,1)
      d(i,3)=a(i,1)*b(i,3)+a(i,3)*b(i,1)
      d(i,4)=a(i,1)*b(i,4)+a(i,4)*b(i,1)
      d(i,5)=a(i,1)*b(i,5)+a(i,5)*b(i,1)
      if(la(i).eq.0.or.lb(i).eq.0) cycle
      ! <p|g> = (<s|+<p|)*(|s>+|p>+|d>+|f>+|g>)
      !       = <s> + <p> + <d> + <f> + <g> + <h>
      d(i,3)=d(i,3)+a(i,2)*b(i,2)
      d(i,4)=d(i,4)+a(i,2)*b(i,3)+a(i,3)*b(i,2)
      d(i,5)=d(i,5)+a(i,2)*b(i,4)+a(i,4)*b(i,2)
      d(i,6)=a(i,2)*b(i,5)+a(i,5)*b(i,2)
      if(la(i).le.1.or.lb(i).le.1) cycle
      ! <d|g> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>+|g>)
      !       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
      d(i,5)=d(i,5)+a(i,3)*b(i,3)
      d(i,6)=d(i,5)+a(i,3)*b(i,4)+a(i,4)*b(i,3)
      d(i,7)=a(i,3)*b(i,5)+a(i,5)*b(i,3)
      if(la(i).le.2.or.lb(i).le.2) cycle
      ! <f|g> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>+|g>)
      !       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k>
      d(i,7)=d(i,7)+a(i,4)*b(i,4)
      d(i,8)=a(i,4)*b(i,5)+a(i,5)*b(i,4)
      if(la(i).le.3.or.lb(i).le.3) cycle
      ! <g|g> = (<s|+<p|+<d|+<f|+<g|)*(|s>+|p>+|d>+|f>+|g>)
      !       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k> + <l>
      d(i,9)=a(i,5)*b(i,5)
   enddo
end subroutine prod2

! --------------------------------------------------------------[SAW1801]-
pure subroutine dsawab(l,ga,v,d)
   implicit none
   integer,intent(in)  :: l
   real(wp), intent(in)  :: ga,d
   real(wp), intent(out) :: v(*)
   real(wp)  :: t(3)
   ! --- for overlap, first moment and second moment
   t(1:3)=olapp((/l,l+1,l+2/),ga)
   ! --- build everything together
   v(1) = t(1)
   v(2) = t(2)+d*t(1)
   v(3) = t(3)+2*d*t(2)+d**2*t(1)
end subroutine dsawab

! --------------------------------------------------------------[SAW1907]-
!     a: center of first gaussian
!     b: center of second gaussian
!     c: aufpunkt of moment operator
!     alpi: alpha/exponent of a
!     alpj: beta/exponent of b
!     la/lb: defines l
!     nt: dimension of g
pure subroutine build_sdq_ints(a,b,c,alpi,alpj,la,lb,v)
   implicit none
   !     aufpunkte,ref point,intarray
   integer,intent(in)  :: la,lb
   real(wp), intent(in)  :: alpi,alpj
   real(wp), intent(in)  :: a(3),b(3),c(3)
   real(wp), intent(out) :: v(10)
   !     local variables
   real(wp)  :: d(3),dd(0:8,3),va(3),val(3,3)
   real(wp)  :: e(3),aa(0:3,3),bb(0:3,3),t(0:8)
   real(wp)  :: gama,kab
   integer :: i,j,ij(3),ii(3),jj(3),lmax

   val = 0

   aa = 0
   bb = 0
   dd = 0
   ii = [lx(la),ly(la),lz(la)]
   jj = [lx(lb),ly(lb),lz(lb)]
   ij = ii+jj
   aa(ii(1),1)=1.0_wp
   aa(ii(2),2)=1.0_wp
   aa(ii(3),3)=1.0_wp
   bb(jj(1),1)=1.0_wp
   bb(jj(2),2)=1.0_wp
   bb(jj(3),3)=1.0_wp
   ! apply product theorem
   ! e is center of product gaussian with exponent gama
   e = gpcenter(alpi,a,alpj,b)
   ! c is reference point
   d = e - c
   call build_kab(a,alpi,b,alpj,gama,kab)
   do i = 1, 3
      ! calculate cartesian prefactor for first gaussian
      call build_hshift2(aa(:,i),a(i),e(i),ii(i))    ! <a|
      ! calculate cartesian prefactor for second gaussian
      call build_hshift2(bb(:,i),b(i),e(i),jj(i))    ! |b>
      ! form their product
      call prod3(aa(:,i),bb(:,i),dd(:,i),ii(i),jj(i))
      lmax = ij(i)
      t(:lmax+2) = olapp([(j,j=0,lmax+2)],gama)
      do j = 0, lmax
         !    <a|b> <a|x|b>             <a|x²|b>
         va = [t(j), t(j+1) + d(i)*t(j), t(j+2) + 2*d(i)*t(j+1) + d(i)*d(i)*t(j)]
         val(i,1:3) = val(i,1:3) + dd(j,i)*va(1:3)
      enddo
   enddo
   v( 1)=kab*(val(1,1)*val(2,1)*val(3,1))
   v( 2)=kab*(val(1,2)*val(2,1)*val(3,1))
   v( 3)=kab*(val(1,1)*val(2,2)*val(3,1))
   v( 4)=kab*(val(1,1)*val(2,1)*val(3,2))
   v( 5)=kab*(val(1,3)*val(2,1)*val(3,1))
   v( 6)=kab*(val(1,1)*val(2,3)*val(3,1))
   v( 7)=kab*(val(1,1)*val(2,1)*val(3,3))
   v( 8)=kab*(val(1,2)*val(2,2)*val(3,1))
   v( 9)=kab*(val(1,2)*val(2,1)*val(3,2))
   v(10)=kab*(val(1,1)*val(2,2)*val(3,2))

end subroutine build_sdq_ints

! --------------------------------------------------------------[SAW1907]-
!     a: center of first gaussian
!     b: center of second gaussian
!     c: aufpunkt of moment operator
!     alpi: alpha/exponent of a
!     alpj: beta/exponent of b
!     la/lb: defines l
!     g: gradient
!     nt: dimension of g
pure subroutine build_dsdq_ints(a,b,c,alpi,alpj,la,lb,v,g)
   implicit none
   !     aufpunkte,ref point,intarray
   integer,intent(in)  :: la,lb
   real(wp), intent(in)  :: alpi,alpj
   real(wp), intent(in)  :: a(3),b(3),c(3)
   real(wp), intent(out) :: v(10),g(3,10)
   !     local variables
   real(wp)  :: d(3),dd(0:8,3),gg(0:8,3),va(3),val(3,3),gra(3,3)
   real(wp)  :: e(3),aa(0:3,3),aap(0:4,3),aam(0:4,3),bb(0:3,3),t(0:8)
   real(wp)  :: gama,kab
   integer :: i,j,ij(3),ii(3),jj(3),lmax

   val = 0
   gra = 0

   aa = 0
   bb = 0
   dd = 0
   gg = 0
   ii = [lx(la),ly(la),lz(la)]
   jj = [lx(lb),ly(lb),lz(lb)]
   ij = ii+jj
   aa(ii(1),1)=1.0_wp
   aa(ii(2),2)=1.0_wp
   aa(ii(3),3)=1.0_wp
   bb(jj(1),1)=1.0_wp
   bb(jj(2),2)=1.0_wp
   bb(jj(3),3)=1.0_wp
   !     d/dX<a|b> = alpi<a+1|b> - i<a-1|b> = -alpj<a|b+1> + j<a|b-1>
   aap=0
   aam=0
   aap(ii(1)+1,1)=2*alpi
   aap(ii(2)+1,2)=2*alpi
   aap(ii(3)+1,3)=2*alpi
   if(ii(1).gt.0) aam(ii(1)-1,1)=-ii(1)
   if(ii(2).gt.0) aam(ii(2)-1,2)=-ii(2)
   if(ii(3).gt.0) aam(ii(3)-1,3)=-ii(3)
   ! apply product theorem
   ! e is center of product gaussian with exponent gama
   e = gpcenter(alpi,a,alpj,b)
   ! c is reference point
   d = e - c
   call build_kab(a,alpi,b,alpj,gama,kab)
   do i = 1, 3
      ! calculate cartesian prefactor for first gaussian
      call build_hshift2(aa(:,i),a(i),e(i),ii(i))    ! <a|
      call build_hshift2(aap(:,i),a(i),e(i),ii(i)+1) ! <a+1|
      call build_hshift2(aam(:,i),a(i),e(i),ii(i)-1) ! <a-1|
      ! calculate cartesian prefactor for second gaussian
      call build_hshift2(bb(:,i),b(i),e(i),jj(i))    ! |b>
      ! form their product
      call prod3(aa(:,i),bb(:,i),dd(:,i),ii(i),jj(i))
      aap(:,i) = aap(:,i) + aam(:,i)
      call prod3(aap(:,i),bb(:,i),gg(:,i),ii(i)+1,jj(i))
      lmax = ij(i)+1
      t(:lmax+2) = olapp([(j,j=0,lmax+2)],gama)
      do j = 0, lmax
         !    <a|b> <a|x|b>             <a|x²|b>
         va = [t(j), t(j+1) + d(i)*t(j), t(j+2) + 2*d(i)*t(j+1) + d(i)*d(i)*t(j)]
         val(i,1:3) = val(i,1:3) + dd(j,i)*va(1:3)
         gra(i,1:3) = gra(i,1:3) + gg(j,i)*va(1:3)
      enddo
   enddo
   v( 1)=kab*(val(1,1)*val(2,1)*val(3,1))
   v( 2)=kab*(val(1,2)*val(2,1)*val(3,1))
   v( 3)=kab*(val(1,1)*val(2,2)*val(3,1))
   v( 4)=kab*(val(1,1)*val(2,1)*val(3,2))
   v( 5)=kab*(val(1,3)*val(2,1)*val(3,1))
   v( 6)=kab*(val(1,1)*val(2,3)*val(3,1))
   v( 7)=kab*(val(1,1)*val(2,1)*val(3,3))
   v( 8)=kab*(val(1,2)*val(2,2)*val(3,1))
   v( 9)=kab*(val(1,2)*val(2,1)*val(3,2))
   v(10)=kab*(val(1,1)*val(2,2)*val(3,2))
   g(1, 1)=kab*(gra(1,1)*val(2,1)*val(3,1))
   g(2, 1)=kab*(val(1,1)*gra(2,1)*val(3,1))
   g(3, 1)=kab*(val(1,1)*val(2,1)*gra(3,1))
   g(1, 2)=kab*(gra(1,2)*val(2,1)*val(3,1))
   g(2, 2)=kab*(val(1,2)*gra(2,1)*val(3,1))
   g(3, 2)=kab*(val(1,2)*val(2,1)*gra(3,1))
   g(1, 3)=kab*(gra(1,1)*val(2,2)*val(3,1))
   g(2, 3)=kab*(val(1,1)*gra(2,2)*val(3,1))
   g(3, 3)=kab*(val(1,1)*val(2,2)*gra(3,1))
   g(1, 4)=kab*(gra(1,1)*val(2,1)*val(3,2))
   g(2, 4)=kab*(val(1,1)*gra(2,1)*val(3,2))
   g(3, 4)=kab*(val(1,1)*val(2,1)*gra(3,2))
   g(1, 5)=kab*(gra(1,3)*val(2,1)*val(3,1))
   g(2, 5)=kab*(val(1,3)*gra(2,1)*val(3,1))
   g(3, 5)=kab*(val(1,3)*val(2,1)*gra(3,1))
   g(1, 6)=kab*(gra(1,1)*val(2,3)*val(3,1))
   g(2, 6)=kab*(val(1,1)*gra(2,3)*val(3,1))
   g(3, 6)=kab*(val(1,1)*val(2,3)*gra(3,1))
   g(1, 7)=kab*(gra(1,1)*val(2,1)*val(3,3))
   g(2, 7)=kab*(val(1,1)*gra(2,1)*val(3,3))
   g(3, 7)=kab*(val(1,1)*val(2,1)*gra(3,3))
   g(1, 8)=kab*(gra(1,2)*val(2,2)*val(3,1))
   g(2, 8)=kab*(val(1,2)*gra(2,2)*val(3,1))
   g(3, 8)=kab*(val(1,2)*val(2,2)*gra(3,1))
   g(1, 9)=kab*(gra(1,2)*val(2,1)*val(3,2))
   g(2, 9)=kab*(val(1,2)*gra(2,1)*val(3,2))
   g(3, 9)=kab*(val(1,2)*val(2,1)*gra(3,2))
   g(1,10)=kab*(gra(1,1)*val(2,2)*val(3,2))
   g(2,10)=kab*(val(1,1)*gra(2,2)*val(3,2))
   g(3,10)=kab*(val(1,1)*val(2,2)*gra(3,2))

end subroutine build_dsdq_ints

! --------------------------------------------------------------[SAW1801]-
!> move gradient operator from center a to center b
!  might look complicated, but take it from me: integrals are usually complicated.
pure subroutine shiftintg(g,s,r)
   implicit none
   real(wp),intent(inout) :: g(3,19)
   real(wp),intent(in)    :: s(10),r(3)
   g(:,11)=g(:,2)-r(1)*g(:,1)
   g(:,12)=g(:,3)-r(2)*g(:,1)
   g(:,13)=g(:,4)-r(3)*g(:,1)
   g(:,14)=g(:,5)-2*r(1)*g(:,2)+r(1)**2*g(:,1)
   g(:,15)=g(:,6)-2*r(2)*g(:,3)+r(2)**2*g(:,1)
   g(:,16)=g(:,7)-2*r(3)*g(:,4)+r(3)**2*g(:,1)
   g(:,17)=g(:,8)-r(1)*g(:,3)-r(2)*g(:,2)+r(1)*r(2)*g(:,1)
   g(:,18)=g(:,9)-r(1)*g(:,4)-r(3)*g(:,2)+r(1)*r(3)*g(:,1)
   g(:,19)=g(:,10)-r(2)*g(:,4)-r(3)*g(:,3)+r(2)*r(3)*g(:,1)
   g(1,11)=g(1,11)-s(1)
   g(2,12)=g(2,12)-s(1)
   g(3,13)=g(3,13)-s(1)
   g(1,14)=g(1,14)-2*s(2)+2*r(1)*s(1)
   g(2,15)=g(2,15)-2*s(3)+2*r(2)*s(1)
   g(3,16)=g(3,16)-2*s(4)+2*r(3)*s(1)
   g(1,17)=g(1,17)-s(3)+r(2)*s(1)
   g(2,17)=g(2,17)-s(2)+r(1)*s(1)
   g(1,18)=g(1,18)-s(4)+r(3)*s(1)
   g(3,18)=g(3,18)-s(2)+r(1)*s(1)
   g(2,19)=g(2,19)-s(4)+r(3)*s(1)
   g(3,19)=g(3,19)-s(3)+r(2)*s(1)
end subroutine shiftintg

! --------------------------------------------------------------[SAW1805]-
!     a: center of first gaussian
!     b: center of second gaussian
!     c: aufpunkt of moment operator
!     alpi: alpha/exponent of a
!     alpj: beta/exponent of b
!     la/lb: defines l
!     g: gradient
!     nt: dimension of g
pure subroutine build_ds_ints(a,b,alpi,alpj,la,lb,v,g)
   implicit none
   !     aufpunkte,ref point,intarray
   integer,intent(in)  :: la,lb
   real(wp), intent(in)  :: alpi,alpj
   real(wp), intent(in)  :: a(3),b(3)
   real(wp), intent(out) :: v,g(3)
   !     local variables
   real(wp)  :: d(3),dd(0:8,3),gg(0:8,3),va(3),val(3),gra(3),t
   real(wp)  :: e(3),aa(0:3,3),aap(0:4,3),aam(0:4,3),bb(0:3,3)
   real(wp)  :: gama,kab
   integer :: i,j,ij(3),ii(3),jj(3)

   val = 0
   gra = 0

   aa = 0
   bb = 0
   dd = 0
   gg = 0
   ii = [lx(la),ly(la),lz(la)]
   jj = [lx(lb),ly(lb),lz(lb)]
   ij = ii+jj
   aa(ii(1),1)=1.0_wp
   aa(ii(2),2)=1.0_wp
   aa(ii(3),3)=1.0_wp
   bb(jj(1),1)=1.0_wp
   bb(jj(2),2)=1.0_wp
   bb(jj(3),3)=1.0_wp
   !     d/dX<a|b> = alpi<a+1|b> - i<a-1|b> = -alpj<a|b+1> + j<a|b-1>
   aap=0
   aam=0
   aap(ii(1)+1,1)=2*alpi
   aap(ii(2)+1,2)=2*alpi
   aap(ii(3)+1,3)=2*alpi
   if(ii(1).gt.0) aam(ii(1)-1,1)=-ii(1)
   if(ii(2).gt.0) aam(ii(2)-1,2)=-ii(2)
   if(ii(3).gt.0) aam(ii(3)-1,3)=-ii(3)
   ! --- apply product theorem, not affected by gradient?
   e = gpcenter(alpi,a,alpj,b)
   call build_kab(a,alpi,b,alpj,gama,kab)
   do i = 1, 3
      ! --- calculate cartesian prefactor for first gaussian
      call build_hshift2(aa(:,i),a(i),e(i),ii(i))    ! <a|
      call build_hshift2(aap(:,i),a(i),e(i),ii(i)+1) ! <a+1|
      call build_hshift2(aam(:,i),a(i),e(i),ii(i)-1) ! <a-1|
      ! --- calculate cartesian prefactor for second gaussian
      call build_hshift2(bb(:,i),b(i),e(i),jj(i))    ! |b>
      ! --- form their product
      call prod3(aa(:,i),bb(:,i),dd(:,i),ii(i),jj(i))
      aap(:,i) = aap(:,i) + aam(:,i)
      call prod3(aap(:,i),bb(:,i),gg(:,i),ii(i)+1,jj(i))
      ! --- e is center of product gaussian with exponent gama
      do j=0,ij(i)+1
         t = olapp(j,gama)
         val(i) = val(i) + dd(j,i)*t
         gra(i) = gra(i) + gg(j,i)*t
      enddo
   enddo
   v=kab*(val(1)*val(2)*val(3))
   g(1)=kab*(gra(1)*val(2)*val(3))
   g(2)=kab*(val(1)*gra(2)*val(3))
   g(3)=kab*(val(1)*val(2)*gra(3))

end subroutine build_ds_ints

end module intgrad
