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

! --------------------------------------------------------------[SAW1907]-
!> Backend implementation for all overlap distribution based integrals,
!  we use a hardcoded horizontal Obara--Saika recursion relation to get
!  the job done, this is working code, so think twice before modifying it!
module xtb_intgrad
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   implicit none

   integer, parameter :: maxl = 6
   integer, parameter :: maxl2 = maxl*2

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
   integer, private, parameter :: lxyz(3, 84) = reshape(&
      & [lx, ly, lz], shape(lxyz), order=[2, 1])
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

   integer,private, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer,private, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

   real(wp),private, parameter :: trafo(6,6) = reshape((/ & ! copied from scf.f, simplyfied
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
      & 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp /), shape(trafo))

   integer, parameter, private :: itt(0:3) = [0,1,4,10]

contains


! --------------------------------------------------------------[SAW1907]-
!> calculates a partial overlap in one cartesian direction
pure elemental function olapp(l,gama) result(s)
   !$acc routine seq
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
  !$acc routine seq
   real(wp),intent(in) :: alp,bet
   real(wp),intent(in) :: ra(3),rb(3)
   real(wp) :: rp(3)

   rp = (alp*ra + bet*rb)/(alp+bet)

end function gpcenter

! --------------------------------------------------------------[SAW1907]-
pure subroutine build_kab(ra,alp,rb,bet,gama,kab)
  !$acc routine seq
   use xtb_mctc_constants
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

! --------------------------------------------------------------[SAW1801]-
pure subroutine dtrf2(s,li,lj)
   use xtb_mctc_blas, only : mctc_gemm
   real(wp),intent(inout) :: s(6,6)
   integer, intent(in)    :: li,lj
   ! CAO-AO transformation
   real(wp), parameter :: trafo(6,6) = reshape((/ & ! copied from scf.f, simplyfied
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
      & 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp /), shape(trafo))

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
   CALL mctc_gemm(trafo, s, dum, transa='T')
   CALL mctc_gemm(dum, trafo, s2, transb='N')
   s(1:5,1:5) = s2(2:6,2:6)
   return

end subroutine dtrf2

! --------------------------------------------------------------[SAW1801]-
pure subroutine build_hshift2(cfs,a,e,l)
   !$acc routine seq
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
   !$acc routine seq
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

pure subroutine horizontal_shift(ae, l, cfs)
   integer, intent(in) :: l
   real(wp), intent(in) :: ae
   real(wp), intent(inout) :: cfs(*)
   integer :: i
   select case(l)
   case(0) ! s
      continue
   case(1) ! p
      cfs(1)=cfs(1)+ae*cfs(2)
   case(2) ! d
      cfs(1)=cfs(1)+ae*ae*cfs(3)
      cfs(2)=cfs(2)+ 2*ae*cfs(3)
   case(3) ! f
      cfs(1)=cfs(1)+ae*ae*ae*cfs(4)
      cfs(2)=cfs(2)+ 3*ae*ae*cfs(4)
      cfs(3)=cfs(3)+ 3*ae*cfs(4)
   case(4) ! g
      cfs(1)=cfs(1)+ae*ae*ae*ae*cfs(5)
      cfs(2)=cfs(2)+ 4*ae*ae*ae*cfs(5)
      cfs(3)=cfs(3)+ 6*ae*ae*cfs(5)
      cfs(4)=cfs(4)+ 4*ae*cfs(5)
   end select
end subroutine horizontal_shift

pure subroutine form_product(a, b, la, lb, d)
   integer, intent(in) :: la, lb
   real(wp), intent(in) :: a(*), b(*)
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

end subroutine form_product


pure subroutine overlap_3d(ri, rj, rp, ai, aj, li, lj, s1d, s3d)
   real(wp), intent(in) :: ri(3)
   real(wp), intent(in) :: rj(3)
   real(wp), intent(in) :: rp(3)
   real(wp), intent(in) :: ai
   real(wp), intent(in) :: aj
   integer, intent(in) :: li(3)
   integer, intent(in) :: lj(3)
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: s3d

   integer :: i, j, k, l
   real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), val(3)

   val(:) = 0.0_wp

   do k = 1, 3
      vv(:) = 0.0_wp
      vi(:) = 0.0_wp
      vj(:) = 0.0_wp
      vi(li(k)) = 1.0_wp
      vj(lj(k)) = 1.0_wp

      call horizontal_shift(rp(k) - ri(k), li(k), vi)
      call horizontal_shift(rp(k) - rj(k), lj(k), vj)
      call form_product(vi, vj, li(k), lj(k), vv)
      do l = 0, li(k) + lj(k)
         val(k) = val(k) + s1d(l) * vv(l)
      end do
   end do

   s3d = val(1) * val(2) * val(3)

end subroutine overlap_3d

pure subroutine multipole_3d(ri, rj, rc, rp, ai, aj, li, lj, s1d, s3d)
   real(wp), intent(in) :: ri(3)
   real(wp), intent(in) :: rj(3)
   real(wp), intent(in) :: rp(3)
   real(wp), intent(in) :: rc(3)
   real(wp), intent(in) :: ai
   real(wp), intent(in) :: aj
   integer, intent(in) :: li(3)
   integer, intent(in) :: lj(3)
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: s3d(10)

   integer :: k, l
   real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), val(3, 3), rpc

   val(:, :) = 0.0_wp

   do k = 1, 3
      vv(:) = 0.0_wp
      vi(:) = 0.0_wp
      vj(:) = 0.0_wp
      vi(li(k)) = 1.0_wp
      vj(lj(k)) = 1.0_wp
      rpc = rp(k) - rc(k)

      call horizontal_shift(rp(k) - ri(k), li(k), vi)
      call horizontal_shift(rp(k) - rj(k), lj(k), vj)
      call form_product(vi, vj, li(k), lj(k), vv)
      do l = 0, li(k) + lj(k)
         val(k, 1) = val(k, 1) + s1d(l) * vv(l)
         val(k, 2) = val(k, 2) + (s1d(l+1) + rpc*s1d(l)) * vv(l)
         val(k, 3) = val(k, 3) + (s1d(l+2) + 2*rpc*s1d(l+1) + rpc*rpc*s1d(l)) * vv(l)
      end do
   end do

   s3d( 1) = val(1, 1) * val(2, 1) * val(3, 1)
   s3d( 2) = val(1, 2) * val(2, 1) * val(3, 1)
   s3d( 3) = val(1, 1) * val(2, 2) * val(3, 1)
   s3d( 4) = val(1, 1) * val(2, 1) * val(3, 2)
   s3d( 5) = val(1, 3) * val(2, 1) * val(3, 1)
   s3d( 6) = val(1, 1) * val(2, 3) * val(3, 1)
   s3d( 7) = val(1, 1) * val(2, 1) * val(3, 3)
   s3d( 8) = val(1, 2) * val(2, 2) * val(3, 1)
   s3d( 9) = val(1, 2) * val(2, 1) * val(3, 2)
   s3d(10) = val(1, 1) * val(2, 2) * val(3, 2)

end subroutine multipole_3d

pure subroutine overlap_grad_3d(ri, rj, rp, ai, aj, li, lj, s1d, s3d, ds3d)
   real(wp), intent(in) :: ri(3)
   real(wp), intent(in) :: rj(3)
   real(wp), intent(in) :: rp(3)
   real(wp), intent(in) :: ai
   real(wp), intent(in) :: aj
   integer, intent(in) :: li(3)
   integer, intent(in) :: lj(3)
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: s3d
   real(wp), intent(out) :: ds3d(3)

   integer :: k, l
   real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), val(3)
   real(wp) :: gi(0:maxl), gg(0:maxl2), gra(3)

   val(:) = 0.0_wp
   gra(:) = 0.0_wp

   do k = 1, 3
      vv(:) = 0.0_wp
      gg(:) = 0.0_wp
      vi(:) = 0.0_wp
      vj(:) = 0.0_wp
      gi(:) = 0.0_wp

      vi(li(k)) = 1.0_wp
      vj(lj(k)) = 1.0_wp
      gi(li(k)+1) = 2*ai
      if (li(k) > 0) gi(li(k)-1) = -li(k)

      call horizontal_shift(rp(k) - ri(k), li(k)-1, gi)
      call horizontal_shift(rp(k) - ri(k), li(k)+1, gi)
      call horizontal_shift(rp(k) - ri(k), li(k), vi)
      call horizontal_shift(rp(k) - rj(k), lj(k), vj)
      call form_product(vi, vj, li(k), lj(k), vv)
      call form_product(gi, vj, li(k)+1, lj(k), gg)
      do l = 0, li(k) + lj(k) + 1
         val(k) = val(k) + s1d(l) * vv(l)
         gra(k) = gra(k) + s1d(l) * gg(l)
      end do
   end do

   s3d = val(1) * val(2) * val(3)
   ds3d(1) = gra(1) * val(2) * val(3)
   ds3d(2) = val(1) * gra(2) * val(3)
   ds3d(3) = val(1) * val(2) * gra(3)

end subroutine overlap_grad_3d

pure subroutine multipole_grad_3d(ri, rj, rc, rp, ai, aj, li, lj, s1d, s3d, ds3d)
   real(wp), intent(in) :: ri(3)
   real(wp), intent(in) :: rj(3)
   real(wp), intent(in) :: rc(3)
   real(wp), intent(in) :: rp(3)
   real(wp), intent(in) :: ai
   real(wp), intent(in) :: aj
   integer, intent(in) :: li(3)
   integer, intent(in) :: lj(3)
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: s3d(10)
   real(wp), intent(out) :: ds3d(3, 10)

   integer :: k, l
   real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), val(3, 3), rpc
   real(wp) :: gi(0:maxl), gg(0:maxl2), gra(3, 3)

   val(:, :) = 0.0_wp
   gra(:, :) = 0.0_wp

   do k = 1, 3
      vv(:) = 0.0_wp
      gg(:) = 0.0_wp
      vi(:) = 0.0_wp
      vj(:) = 0.0_wp
      gi(:) = 0.0_wp
      rpc = rp(k) - rc(k)

      vi(li(k)) = 1.0_wp
      vj(lj(k)) = 1.0_wp
      gi(li(k)+1) = 2*ai
      if (li(k) > 0) gi(li(k)-1) = -li(k)

      call horizontal_shift(rp(k) - ri(k), li(k)-1, gi)
      call horizontal_shift(rp(k) - ri(k), li(k)+1, gi)
      call horizontal_shift(rp(k) - ri(k), li(k), vi)
      call horizontal_shift(rp(k) - rj(k), lj(k), vj)
      call form_product(vi, vj, li(k), lj(k), vv)
      call form_product(gi, vj, li(k)+1, lj(k), gg)
      do l = 0, li(k) + lj(k) + 1
         val(k, 1) = val(k, 1) + s1d(l) * vv(l)
         val(k, 2) = val(k, 2) + (s1d(l+1) + rpc*s1d(l)) * vv(l)
         val(k, 3) = val(k, 3) + (s1d(l+2) + 2*rpc*s1d(l+1) + rpc*rpc*s1d(l)) * vv(l)
         gra(k, 1) = gra(k, 1) + s1d(l) * gg(l)
         gra(k, 2) = gra(k, 2) + (s1d(l+1) + rpc*s1d(l)) * gg(l)
         gra(k, 3) = gra(k, 3) + (s1d(l+2) + 2*rpc*s1d(l+1) + rpc*rpc*s1d(l)) * gg(l)
      end do
   end do

   s3d( 1) = val(1, 1) * val(2, 1) * val(3, 1)
   s3d( 2) = val(1, 2) * val(2, 1) * val(3, 1)
   s3d( 3) = val(1, 1) * val(2, 2) * val(3, 1)
   s3d( 4) = val(1, 1) * val(2, 1) * val(3, 2)
   s3d( 5) = val(1, 3) * val(2, 1) * val(3, 1)
   s3d( 6) = val(1, 1) * val(2, 3) * val(3, 1)
   s3d( 7) = val(1, 1) * val(2, 1) * val(3, 3)
   s3d( 8) = val(1, 2) * val(2, 2) * val(3, 1)
   s3d( 9) = val(1, 2) * val(2, 1) * val(3, 2)
   s3d(10) = val(1, 1) * val(2, 2) * val(3, 2)
   ds3d(1, 1) = gra(1, 1) * val(2, 1) * val(3, 1)
   ds3d(2, 1) = val(1, 1) * gra(2, 1) * val(3, 1)
   ds3d(3, 1) = val(1, 1) * val(2, 1) * gra(3, 1)
   ds3d(1, 2) = gra(1, 2) * val(2, 1) * val(3, 1)
   ds3d(2, 2) = val(1, 2) * gra(2, 1) * val(3, 1)
   ds3d(3, 2) = val(1, 2) * val(2, 1) * gra(3, 1)
   ds3d(1, 3) = gra(1, 1) * val(2, 2) * val(3, 1)
   ds3d(2, 3) = val(1, 1) * gra(2, 2) * val(3, 1)
   ds3d(3, 3) = val(1, 1) * val(2, 2) * gra(3, 1)
   ds3d(1, 4) = gra(1, 1) * val(2, 1) * val(3, 2)
   ds3d(2, 4) = val(1, 1) * gra(2, 1) * val(3, 2)
   ds3d(3, 4) = val(1, 1) * val(2, 1) * gra(3, 2)
   ds3d(1, 5) = gra(1, 3) * val(2, 1) * val(3, 1)
   ds3d(2, 5) = val(1, 3) * gra(2, 1) * val(3, 1)
   ds3d(3, 5) = val(1, 3) * val(2, 1) * gra(3, 1)
   ds3d(1, 6) = gra(1, 1) * val(2, 3) * val(3, 1)
   ds3d(2, 6) = val(1, 1) * gra(2, 3) * val(3, 1)
   ds3d(3, 6) = val(1, 1) * val(2, 3) * gra(3, 1)
   ds3d(1, 7) = gra(1, 1) * val(2, 1) * val(3, 3)
   ds3d(2, 7) = val(1, 1) * gra(2, 1) * val(3, 3)
   ds3d(3, 7) = val(1, 1) * val(2, 1) * gra(3, 3)
   ds3d(1, 8) = gra(1, 2) * val(2, 2) * val(3, 1)
   ds3d(2, 8) = val(1, 2) * gra(2, 2) * val(3, 1)
   ds3d(3, 8) = val(1, 2) * val(2, 2) * gra(3, 1)
   ds3d(1, 9) = gra(1, 2) * val(2, 1) * val(3, 2)
   ds3d(2, 9) = val(1, 2) * gra(2, 1) * val(3, 2)
   ds3d(3, 9) = val(1, 2) * val(2, 1) * gra(3, 2)
   ds3d(1,10) = gra(1, 1) * val(2, 2) * val(3, 2)
   ds3d(2,10) = val(1, 1) * gra(2, 2) * val(3, 2)
   ds3d(3,10) = val(1, 1) * val(2, 2) * gra(3, 2)

end subroutine multipole_grad_3d


! --------------------------------------------------------------[SAW1801]-
!> move gradient operator from center a to center b
!  might look complicated, but take it from me: integrals are usually complicated.
pure subroutine shiftintg(g,s,r)
   !$acc routine seq
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


pure subroutine get_overlap(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj,point,intcut, &
      &                nprim,primcount,alp,cont,sint)
   integer, intent(in)  :: icao
   integer, intent(in)  :: jcao
   integer, intent(in)  :: naoi
   integer, intent(in)  :: naoj
   integer, intent(in)  :: ishtyp
   integer, intent(in)  :: jshtyp
   real(wp),intent(in)  :: ri(3)
   real(wp),intent(in)  :: rj(3)
   real(wp),intent(in)  :: point(3)
   real(wp),intent(in)  :: intcut
   real(wp),intent(out) :: sint(:,:)

   integer, intent(in)  :: nprim(:)
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)

   integer  :: ip,iprim,mli,jp,jprim,mlj,k,iptyp,jptyp
   real(wp) :: rij(3),rij2,alpi,alpj,ci,cj,cc,kab,rp(3),t(0:8)
   real(wp) :: ab,est,saw(10)

   real(wp),parameter :: max_r2 = 2000.0_wp
   real(wp),parameter :: sqrtpi = sqrt(pi)

   sint = 0.0_wp
   iptyp = itt(ishtyp)
   jptyp = itt(jshtyp)

   rij = ri - rj
   rij2 = rij(1)**2 + rij(2)**2 + rij(3)**2

   if(rij2.gt.max_r2) return

   do ip = 1, nprim(icao+1)
      iprim = ip + primcount(icao+1)
      ! exponent the same for each l component
      alpi = alp(iprim)
      do jp = 1, nprim(jcao+1)
         jprim=jp+primcount(jcao+1)
         ! exponent the same for each l component
         alpj=alp(jprim)
         ab=1.0_wp/(alpi+alpj)
         est=rij2*alpi*alpj*ab
         if(est.gt.intcut) cycle
         kab = exp(-est)*(sqrtpi*sqrt(ab))**3
         rp = (alpi*ri + alpj*rj)*ab
         do k = 0, ishtyp + jshtyp
            t(k) = olapp(k, alpi+alpj)
         end do
         ! now compute integrals
         do mli=1,naoi
            iprim = ip + primcount(icao+mli)
            ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
            ci = cont(iprim)
            do mlj=1,naoj
               jprim = jp + primcount(jcao+mlj)
               saw = 0.0_wp
               ! prim-prim  integrals
               call multipole_3d(ri,rj,point,rp,alpi,alpj, &
                  & lxyz(:,iptyp+mli),lxyz(:,jptyp+mlj),t,saw)
               cc = kab*cont(jprim)*ci
               sint(mlj,mli) = sint(mlj,mli)+saw(1)*cc! pbc_w(jat,iat)
            enddo ! mlj
         enddo ! mli
      enddo ! jp
   enddo ! ip

end subroutine get_overlap

pure subroutine get_grad_overlap(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj,point,intcut, &
      &                     nprim,primcount,alp,cont,sdq,sdqg)
   integer, intent(in)  :: icao
   integer, intent(in)  :: jcao
   integer, intent(in)  :: naoi
   integer, intent(in)  :: naoj
   integer, intent(in)  :: ishtyp
   integer, intent(in)  :: jshtyp
   real(wp),intent(in)  :: ri(3)
   real(wp),intent(in)  :: rj(3)
   real(wp),intent(in)  :: point(3)
   real(wp),intent(in)  :: intcut
   real(wp),intent(out) :: sdq(:,:)
   real(wp),intent(out) :: sdqg(:,:,:)

   integer, intent(in)  :: nprim(:)
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)

   integer  :: ip,iprim,mli,jp,jprim,mlj,k,iptyp,jptyp
   real(wp) :: rij(3),rij2,alpi,alpj,ci,cj,cc,kab,rp(3),t(0:8)
   real(wp) :: ab,est,saw,sawg(3)

   real(wp),parameter :: max_r2 = 2000.0_wp
   real(wp),parameter :: sqrtpi = sqrt(pi)

   sdqg = 0.0_wp
   sdq  = 0.0_wp
   iptyp = itt(ishtyp)
   jptyp = itt(jshtyp)

   rij = ri - rj
   rij2 = rij(1)**2 + rij(2)**2 + rij(3)**2

   if(rij2.gt.max_r2) return

   do ip=1,nprim(icao+1)
      iprim=ip+primcount(icao+1)
      ! exponent the same for each l component
      alpi=alp(iprim)
      do jp=1,nprim(jcao+1)
         jprim=jp+primcount(jcao+1)
         ! exponent the same for each l component
         alpj=alp(jprim)
         ab = 1.0_wp/(alpi+alpj)
         est=alpi*alpj*rij2*ab
         if(est.gt.intcut) cycle
         kab = exp(-est)*(sqrtpi*sqrt(ab))**3
         rp = (alpi*ri + alpj*rj)*ab
         do k = 0, ishtyp + jshtyp + 1
            t(k) = olapp(k, alpi+alpj)
         end do
         !--------------- compute gradient ----------
         ! now compute integrals  for different components of i(e.g., px,py,pz)
         do mli=1,naoi
            iprim=ip+primcount(icao+mli)
            ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
            ci=cont(iprim)
            do mlj=1,naoj
               jprim=jp+primcount(jcao+mlj)
               cc=kab*cont(jprim)*ci
               saw=0;sawg=0
               call overlap_grad_3d(ri,rj,rp,alpi,alpj,lxyz(:,iptyp+mli),lxyz(:,jptyp+mlj),t,saw,sawg)
               sdq(mlj,mli) = sdq(mlj,mli)+saw*cc
               sdqg(:,mlj,mli) = sdqg(:,mlj,mli) &
                  & + sawg(:)*cc
            enddo ! mlj : Cartesian component of j prims
         enddo  ! mli : Cartesian component of i prims
      enddo ! jp : loop over j prims
   enddo  ! ip : loop over i prims
end subroutine get_grad_overlap


pure subroutine get_multiints(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj,point, &
      &                           intcut,nprim,primcount,alp,cont,ss,dd,qq)
   integer, intent(in)  :: icao
   integer, intent(in)  :: jcao
   integer, intent(in)  :: naoi
   integer, intent(in)  :: naoj
   integer, intent(in)  :: ishtyp
   integer, intent(in)  :: jshtyp
   real(wp),intent(in)  :: ri(3)
   real(wp),intent(in)  :: rj(3)
   real(wp),intent(in)  :: point(3)
   real(wp),intent(in)  :: intcut
   real(wp),intent(out) :: ss(:,:)
   real(wp),intent(out) :: dd(:,:,:)
   real(wp),intent(out) :: qq(:,:,:)

   integer, intent(in)  :: nprim(:)
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)

   integer  :: ip,iprim,mli,jp,jprim,mlj,k,iptyp,jptyp
   real(wp) :: rij(3),rp(3),rij2,alpi,alpj,ci,cj,cc,kab,t(0:8)
   real(wp) :: ab,est,saw(10),sawg(3,10)

   real(wp),parameter :: max_r2 = 2000.0_wp
   real(wp),parameter :: sqrtpi = sqrt(pi)

   ss = 0.0_wp
   dd = 0.0_wp
   qq = 0.0_wp
   iptyp = itt(ishtyp)
   jptyp = itt(jshtyp)

   rij = ri - rj
   rij2 = rij(1)**2 + rij(2)**2 + rij(3)**2

   if(rij2.gt.max_r2) return

   ! we go through the primitives (because the screening is the same for all of them)
   do ip = 1,nprim(icao+1)
      iprim = ip+primcount(icao+1)
      ! exponent the same for each l component
      alpi = alp(iprim)
      do jp = 1,nprim(jcao+1)
         jprim = jp+primcount(jcao+1)
         ! exponent the same for each l component
         alpj = alp(jprim)
         ab = 1.0_wp/(alpi+alpj)
         est = alpi*alpj*rij2*ab
         if(est.gt.intcut) cycle
         kab = exp(-est)*(sqrtpi*sqrt(ab))**3
         rp = (alpi*ri + alpj*rj)*ab
         do k = 0, ishtyp + jshtyp + 2
            t(k) = olapp(k, alpi+alpj)
         end do
         !--------------- compute gradient ----------
         ! now compute integrals  for different components of i(e.g., px,py,pz)
         do mli = 1,naoi
            iprim = ip+primcount(icao+mli)
            ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
            ci = cont(iprim)
            do mlj = 1,naoj
               jprim = jp+primcount(jcao+mlj)
               cc = kab*cont(jprim)*ci
               saw = 0;sawg = 0
               call multipole_3d(ri,rj,point,rp,alpi,alpj, &
                  & lxyz(:,iptyp+mli),lxyz(:,jptyp+mlj),t,saw)
               ss(mlj,mli) = ss(mlj,mli) + saw(1)*cc
               dd(:,mlj,mli) = dd(:,mlj,mli) + saw(2:4)*cc
               qq(:,mlj,mli) = qq(:,mlj,mli) + saw(5:10)*cc
            enddo ! mlj : Cartesian component of j prims
         enddo  ! mli : Cartesian component of i prims
      enddo ! jp : loop over j prims
   enddo  ! ip : loop over i prims
end subroutine get_multiints


pure subroutine get_grad_multiint(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj, &
      &                           intcut,nprim,primcount,alp,cont,sdq,sdqg)
   integer, intent(in)  :: icao
   integer, intent(in)  :: jcao
   integer, intent(in)  :: naoi
   integer, intent(in)  :: naoj
   integer, intent(in)  :: ishtyp
   integer, intent(in)  :: jshtyp
   real(wp),intent(in)  :: ri(3)
   real(wp),intent(in)  :: rj(3)
   real(wp),intent(in)  :: intcut
   real(wp),intent(out) :: sdq(:,:,:)
   real(wp),intent(out) :: sdqg(:,:,:,:)

   integer, intent(in)  :: nprim(:)
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)

   integer  :: ip,iprim,mli,jp,jprim,mlj,k,iptyp,jptyp
   real(wp) :: rij(3),rp(3),rij2,alpi,alpj,ci,cj,cc,kab,t(0:8)
   real(wp) :: ab,est,saw(10),sawg(3,10)

   real(wp),parameter :: max_r2 = 2000.0_wp
   real(wp),parameter :: sqrtpi = sqrt(pi)

   sdqg = 0.0_wp
   sdq  = 0.0_wp
   iptyp = itt(ishtyp)
   jptyp = itt(jshtyp)

   rij = ri - rj
   rij2 = rij(1)**2 + rij(2)**2 + rij(3)**2

   if(rij2.gt.max_r2) return

   ! we go through the primitives (because the screening is the same for all of them)
   do ip = 1,nprim(icao+1)
      iprim = ip+primcount(icao+1)
      ! exponent the same for each l component
      alpi = alp(iprim)
      do jp = 1,nprim(jcao+1)
         jprim = jp+primcount(jcao+1)
         ! exponent the same for each l component
         alpj = alp(jprim)
         ab = 1.0_wp/(alpi+alpj)
         est = alpi*alpj*rij2*ab
         if(est.gt.intcut) cycle
         kab = exp(-est)*(sqrtpi*sqrt(ab))**3
         rp = (alpi*ri + alpj*rj)*ab
         do k = 0, ishtyp + jshtyp + 3
            t(k) = olapp(k, alpi+alpj)
         end do
         !--------------- compute gradient ----------
         ! now compute integrals  for different components of i(e.g., px,py,pz)
         do mli = 1,naoi
            iprim = ip+primcount(icao+mli)
            ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
            ci = cont(iprim)
            do mlj = 1,naoj
               jprim = jp+primcount(jcao+mlj)
               cc = kab*cont(jprim)*ci
               saw = 0;sawg = 0
               call multipole_grad_3d(ri,rj,rj,rp,alpi,alpj, &
                  & lxyz(:,iptyp+mli),lxyz(:,jptyp+mlj),t,saw,sawg)
               sdq(:,mlj,mli) = sdq(:,mlj,mli) + saw*cc
               sdqg(:,:10,mlj,mli) = sdqg(:,:10,mlj,mli) + sawg*cc
            enddo ! mlj : Cartesian component of j prims
         enddo  ! mli : Cartesian component of i prims
      enddo ! jp : loop over j prims
   enddo  ! ip : loop over i prims
   do mli = 1,naoi
      do mlj = 1,naoj
         call shiftintg(sdqg(:,:,mlj,mli),sdq(:,mlj,mli),rij)
      enddo
   enddo
end subroutine get_grad_multiint


!> computes the dipole and quadrupole integrals and performs screening to
!  determine, which contribute to potential
subroutine sdqint(nShell,angShell,nat,at,nbf,nao,xyz,intcut,caoshell,saoshell, &
      &           nprim,primcount,alp,cont,sint,dpint,qpint)
   integer, intent(in) :: nShell(:)
   integer, intent(in) :: angShell(:,:)
   !> # of atoms
   integer, intent(in)  :: nat
   !> # of spherical AOs (SAOs)
   integer, intent(in)  :: nao
   !> # of Cartesian AOs (CAOs)
   integer, intent(in)  :: nbf
   integer, intent(in)  :: at(nat)
   !> cartesian coordinates
   real(wp),intent(in)  :: xyz(3,nat)
   !> integral cutoff according to prefactor from Gaussian product theorem
   real(wp),intent(in)  :: intcut
   !> map shell of atom to index in CAO space (lowest Cart. component is taken)
   integer, intent(in)  :: caoshell(:,:)
   !> map shell of atom to index in SAO space (lowest m_l component is taken)
   integer, intent(in)  :: saoshell(:,:)
   integer, intent(in)  :: nprim(:)
   !> index of first primitive (over entire system) of given CAO
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)
   !> overlap integral matrix
   real(wp),intent(out) :: sint(nao,nao)
   !> dipole integral matrix
   real(wp),intent(out) :: dpint(3,nao,nao)
   !> quadrupole integral matrix
   real(wp),intent(out) :: qpint(6,nao,nao)


   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,km,mi,mj,ij
   real(wp) tmp1,tmp2,tmp3,tmp4,step,step2
   real(wp) dx,dy,dz,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cc,cj,alpi,rab2,ab,est
   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

   real(wp)  ra(3),rb(3),f1,f2,point(3)
   real(wp) dtmp(3),qtmp(6),ss(6,6),dd(3,6,6),qq(6,6,6),tmp(6,6)
   integer ip,jp,iat,jat,izp,jzp,ish,jsh,icao,jcao,iao,jao,jshmax
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer itt(0:3)
   parameter(itt  =(/0,1,4,10/))
   real(wp) :: saw(10)


   ! integrals
   sint = 0.0_wp
   dpint = 0.0_wp
   qpint = 0.0_wp
   ! --- Aufpunkt for moment operator
   point = 0.0_wp

   !$OMP PARALLEL DO schedule(runtime) &
   !$omp PRIVATE (iat,jat,izp,cc,ci,ra,rb,saw, &
   !$omp& rab2,jzp,ish,ishtyp,icao,naoi,iptyp, &
   !$omp& jsh,jshmax,jshtyp,jcao,naoj,jptyp,ss,dd,qq, &
   !$omp& est,alpi,alpj,ab,iprim,jprim,ip,jp, &
   !$omp& mli,mlj,tmp,tmp1,tmp2,iao,jao,ii,jj,k,ij) &
   !$omp shared (sint,dpint,qpint)
   do iat = 1,nat
      ra(1:3) = xyz(1:3,iat)
      izp = at(iat)
      do jat = 1,iat-1
         rb(1:3) = xyz(1:3,jat)
         jzp = at(jat)
         rab2 = sum( (rb-ra)**2 )
         !           ints < 1.d-9 for RAB > 40 Bohr
         if(rab2.gt.2000) cycle
         do ish = 1,nShell(izp)
            ishtyp = angShell(ish,izp)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            do jsh = 1,nShell(jzp)
               jshtyp = angShell(jsh,jzp)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)
               ss = 0.0_wp
               dd = 0.0_wp
               qq = 0.0_wp
               call get_multiints(icao,jcao,naoi,naoj,ishtyp,jshtyp,ra,rb,point, &
                  &               intcut,nprim,primcount,alp,cont,ss,dd,qq)
               !transform from CAO to SAO
               call dtrf2(ss,ishtyp,jshtyp)
               do k = 1,3
                  tmp(1:6,1:6) = dd(k,1:6,1:6)
                  call dtrf2(tmp,ishtyp,jshtyp)
                  dd(k,1:6,1:6) = tmp(1:6,1:6)
               enddo
               do k = 1,6
                  tmp(1:6,1:6) = qq(k,1:6,1:6)
                  call dtrf2(tmp,ishtyp,jshtyp)
                  qq(k,1:6,1:6) = tmp(1:6,1:6)
               enddo
               do ii = 1,llao2(ishtyp)
                  iao = ii+saoshell(ish,iat)
                  do jj = 1,llao2(jshtyp)
                     jao = jj+saoshell(jsh,jat)
                     sint(iao, jao) = sint(iao, jao) + ss(jj, ii)
                     sint(jao, iao) = sint(jao, iao) + ss(jj, ii)
                     dpint(1:3, iao, jao) = dpint(1:3, iao, jao) + dd(1:3, jj, ii)
                     dpint(1:3, jao, iao) = dpint(1:3, jao, iao) + dd(1:3, jj, ii)
                     qpint(1:6, iao, jao) = qpint(1:6, iao, jao) + qq(1:6, jj, ii)
                     qpint(1:6, jao, iao) = qpint(1:6, jao, iao) + qq(1:6, jj, ii)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

   ! diagonal elements
   do iat = 1, nat
      ra = xyz(:, iat)
      izp = at(iat)
      do ish = 1, nShell(izp)
         ishtyp = angShell(ish,izp)
         do iao = 1, llao2(ishtyp)
            i = iao+saoshell(ish,iat)
            sint(i,i) = 1.0_wp + sint(i,i)
         end do

         icao = caoshell(ish,iat)
         naoi = llao(ishtyp)
         iptyp = itt(ishtyp)
         do jsh = 1, ish
            jshtyp = angShell(jsh,izp)
            jcao = caoshell(jsh,iat)
            naoj = llao(jshtyp)
            jptyp = itt(jshtyp)
            ss = 0.0_wp
            dd = 0.0_wp
            qq = 0.0_wp
            call get_multiints(icao,jcao,naoi,naoj,ishtyp,jshtyp,ra,ra,point, &
               &               intcut,nprim,primcount,alp,cont,ss,dd,qq)
            !transform from CAO to SAO
            !call dtrf2(ss,ishtyp,jshtyp)
            do k = 1,3
               tmp(1:6, 1:6) = dd(k,1:6, 1:6)
               call dtrf2(tmp, ishtyp, jshtyp)
               dd(k, 1:6, 1:6) = tmp(1:6, 1:6)
            enddo
            do k = 1,6
               tmp(1:6, 1:6) = qq(k, 1:6, 1:6)
               call dtrf2(tmp, ishtyp, jshtyp)
               qq(k, 1:6, 1:6) = tmp(1:6, 1:6)
            enddo
            do ii = 1, llao2(ishtyp)
               iao = ii + saoshell(ish,iat)
               do jj = 1, llao2(jshtyp)
                  jao = jj + saoshell(jsh,iat)
                  if (jao > iao .and. ish ==  jsh) cycle
                  dpint(1:3, iao, jao) = dpint(1:3, iao, jao) + dd(1:3, jj, ii)
                  if (iao /= jao) then
                     dpint(1:3, jao, iao) = dpint(1:3, jao, iao) + dd(1:3, jj, ii)
                  end if
                  qpint(1:6, iao, jao) = qq(1:6, jj, ii)
                  if (jao /= iao) then
                     qpint(1:6, jao, iao) = qq(1:6, jj, ii)
                  end if
               end do
            end do
         end do
      end do
   end do

end subroutine sdqint


pure subroutine build_sdq_ints_gpu(a,b,c,alpi,alpj,la,lb,kab,t,e,lx,ly,lz,v)
   !$acc routine seq
   !     aufpunkte,ref point,intarray
   integer,intent(in)  :: la,lb
   real(wp), intent(in)  :: alpi,alpj
   real(wp), intent(in)  :: a(3),b(3),c(3)
   real(wp), intent(in)  :: kab,t(0:8),e(3)
   integer, intent(in)  :: lx(84),ly(84),lz(84)
   real(wp), intent(out) :: v(10)
   !     local variables
   real(wp)  :: dd(0:8),va(3),val(3,3)
   real(wp)  :: aa(0:3,3),bb(0:3,3)
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

   do i = 1, 3
      ! calculate cartesian prefactor for first gaussian
      call build_hshift2(aa(:,i),a(i),e(i),ii(i))    ! <a|
      ! calculate cartesian prefactor for second gaussian
      call build_hshift2(bb(:,i),b(i),e(i),jj(i))    ! |b>
      ! form their product
      call prod3(aa(:,i),bb(:,i),dd(:),ii(i),jj(i))
      lmax = ij(i)
      do j = 0, lmax
         !    <a|b> <a|x|b>             <a|x²|b>
         va = [t(j), t(j+1) + e(i)*t(j), t(j+2) + 2*e(i)*t(j+1) + e(i)*e(i)*t(j)]
         val(i,1:3) = val(i,1:3) + dd(j)*va(1:3)
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
end subroutine build_sdq_ints_gpu


!> Computes the overlap dipole and quadrupole integrals
subroutine sdqint_gpu(nShell, angShell, nat, at, nbf, nao, xyz, trans, &
      & intcut, caoshell, saoshell, nprim, primcount, alp, cont, &
      & sint, dpint, qpint)
   integer, intent(in) :: nShell(:)
   integer, intent(in) :: angShell(:,:)
   !> # of atoms
   integer, intent(in)  :: nat
   !> # of spherical AOs (SAOs)
   integer, intent(in)  :: nao
   !> # of Cartesian AOs (CAOs)
   integer, intent(in)  :: nbf
   integer, intent(in)  :: at(:)
   !> Cartesian coordinates
   real(wp),intent(in)  :: xyz(:, :)
   !> Integral cutoff according to prefactor from Gaussian product theorem
   real(wp),intent(in)  :: intcut
   real(wp),intent(in)  :: trans(:, :)
   !> Map shell of atom to index in CAO space (lowest Cart. component is taken)
   integer, intent(in)  :: caoshell(:,:)
   !> Map shell of atom to index in SAO space (lowest m_l component is taken)
   integer, intent(in)  :: saoshell(:,:)
   integer, intent(in)  :: nprim(:)
   !> Index of first primitive (over entire system) of given CAO
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)
   !> Overlap integral matrix
   real(wp),intent(out) :: sint(:, :)
   !> Dipole integral matrix
   real(wp),intent(out) :: dpint(:, :, :)
   !> Quadrupole integral matrix
   real(wp),intent(out) :: qpint(:, :, :)


   integer i,j,k,l,m,ii,jj,ij,ioff,joff,kk
   real(wp) thr2,ci,cc,cj,alpi,rab2,ab,est,alpj

   real(wp) ra(3),rb(3),point(3)
   real(wp) dtmp(3),qtmp(6),ss(6,6),dd(6,6,3),qq(6,6,6),tmp(6,6),dd2(3,6,6)
   integer ip,jp,iat,jat,izp,jzp,ish,jsh,icao,jcao,iao,jao,jshmax
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer :: il, jl, itr
   real(wp) :: zi, zj, zetaij, km, hii, hjj, hav, shpoly
   integer itt(0:3)
   parameter(itt  =(/0,1,4,10/))
   real(wp) :: saw(10)

   real(wp) :: gama,kab,t(0:8),e(3)
   real(wp),parameter :: sqrtpi = sqrt(pi)

   real(wp) s2(6,6),dum(6,6),sspher

   !$acc enter data create(sint(:, :), dpint(:, :, :), qpint(:, :, :)) &
   !$acc copyin(xyz, at, nShell, trans, caoshell, saoshell, &
   !$acc& nprim, primcount, alp, cont, intcut, point, angShell(:, :))

   !$acc kernels default(present)
   ! integrals
   sint = 0.0_wp
   dpint = 0.0_wp
   qpint = 0.0_wp
   !$acc end kernels
   ! --- Aufpunkt for moment operator
   point = 0.0_wp

   !$acc parallel vector_length(32) private(ss,dd,qq,rb,iat,jat,izp,ci,ra,& 
   !$acc& rab2,jzp,ish,ishtyp,icao,naoi,iptyp, &
   !$acc& jsh,jshmax,jshtyp,jcao,naoj,jptyp,shpoly, &
   !$acc& est,alpi,alpj,ab,iprim,jprim,ip,jp,il,jl,hii,hjj,km,zi,zj,zetaij,hav, &
   !$acc& mli,mlj,iao,jao,ii,jj,k,ij,itr,ioff,joff,t,e,saw,s2,dum)

   !$acc loop gang collapse(2)
   do iat = 1, nat
      do jat = 1, nat
         if (jat.ge.iat) cycle
         ra(1:3) = xyz(1:3,iat)
         izp = at(iat)
         jzp = at(jat)
         do ish = 1, nShell(izp)
            ishtyp = angShell(ish,izp)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            do jsh = 1, nShell(jzp)
               jshtyp = angShell(jsh,jzp)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               il = ishtyp+1
               jl = jshtyp+1

               !$acc loop private(ss,dd,qq,s2,dum)
               do itr = 1, size(trans, dim=2)
                  !$acc cache(ss,dd,qq,s2,dum)
                  rb(1:3) = xyz(1:3,jat) + trans(:, itr)
                  rab2 = sum( (rb-ra)**2 )

                  ss = 0.0_wp
                  dd = 0.0_wp
                  qq = 0.0_wp

                  !call get_multiints(icao,jcao,naoi,naoj,ishtyp,jshtyp,ra,rb,point, &
                  !   &               intcut,nprim,primcount,alp,cont,ss,dd,qq)
                  !$acc loop vector private(saw,t,e) independent collapse(2)
                  do ip = 1,nprim(icao+1)
                    do jp = 1,nprim(jcao+1)

                      iprim = ip+primcount(icao+1)
                      jprim = jp+primcount(jcao+1)
                      alpi = alp(iprim)
                      alpj = alp(jprim)

                      gama = alpi+alpj
                      ab = 1.0_wp/(gama)
                      e = (alpi * ra + alpj * rb)/(alpi + alpj)
                      t(0:8) = olapp([(j,j=0,8)],gama)
                      est = rab2*alpi*alpj*ab
                      kab = exp(-est)*(sqrtpi*sqrt(ab))**3

                      do mli = 1,naoi
                        do mlj = 1,naoj
                          saw = 0.0_wp

                          ! prim-prim quadrupole and dipole integrals
                          call build_sdq_ints_gpu(ra,rb,point,alpi,alpj, &
                            & iptyp+mli,jptyp+mlj,kab,t,e,lx,ly,lz,saw)

                          iprim = ip+primcount(icao+mli)
                          jprim = jp+primcount(jcao+mlj)
                          ci = cont(iprim) ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
                          cc = cont(jprim)*ci
                          ! from primitive integrals fill CAO-CAO matrix for ish-jsh block
                          !                             ! overlap
                          !$acc atomic
                          ss(mlj,mli) = ss(mlj,mli)+saw(1)*cc
                          !$acc end atomic
                          ! dipole
                          do k=1,3
                            !$acc atomic
                            dd(mlj,mli,k) = dd(mlj,mli,k)+saw(k+1)*cc
                            !$acc end atomic
                          enddo
                          do k = 1,6
                            ! quadrupole
                            !$acc atomic
                            qq(mlj,mli,k) = qq(mlj,mli,k)+saw(k+4)*cc
                            !$acc end atomic
                          enddo
                        enddo ! mlj
                      enddo ! mli
                    enddo ! jp
                  enddo ! ip

                  !transform from CAO to SAO
                  !call dtrf2(ss,ishtyp,jshtyp)
                  ! DTRF2(S)
                  !select case(ishtyp)
                  !case(0) ! s-d
                  if (ishtyp.ge.2.or.jshtyp.ge.2) then
                    if (ishtyp.eq.0) then
                      !ss
                      do jj=1,6
                        sspher=0
                        do m=1,6
                          sspher=sspher+trafo(m,jj)*ss(m,1)
                        enddo
                        s2(jj,1)=sspher
                      enddo
                      ss(1:5,1) = s2(2:6,1)

                      !dd
                      do k = 1,3
                        do jj=1,6
                          sspher=0
                          do m=1,6
                            sspher=sspher+trafo(m,jj)*dd(m,1,k)
                          enddo
                          s2(jj,1)=sspher
                        enddo
                        dd(1:5,1,k) = s2(2:6,1)
                      enddo

                      !qq
                      do k = 1,6
                        do jj=1,6
                          sspher=0
                          do m=1,6
                            sspher=sspher+trafo(m,jj)*qq(m,1,k)
                          enddo
                          s2(jj,1)=sspher
                        enddo
                        qq(1:5,1,k) = s2(2:6,1)
                      enddo

                    !case(1) ! p-d
                    else if (ishtyp.eq.1) then
                      !ss
                      do ii=1,3
                        do jj=1,6
                          sspher=0
                          do m=1,6
                            sspher=sspher+trafo(m,jj)*ss(m,ii)
                          enddo
                          s2(jj,ii)=sspher
                        enddo
                        ss(1:5,ii) = s2(2:6,ii)
                      enddo

                      ! dd
                      do k = 1,3
                        do ii=1,3
                          do jj=1,6
                            sspher=0
                            do m=1,6
                              sspher=sspher+trafo(m,jj)*dd(m,ii,k)
                            enddo
                            s2(jj,ii)=sspher
                          enddo
                          dd(1:5,ii,k) = s2(2:6,ii)
                        enddo
                      enddo

                      !qq
                      do k = 1,6
                        do ii=1,3
                          do jj=1,6
                            sspher=0
                            do m=1,6
                              sspher=sspher+trafo(m,jj)*qq(m,ii,k)
                            enddo
                            s2(jj,ii)=sspher
                          enddo
                          qq(1:5,ii,k) = s2(2:6,ii)
                        enddo
                      enddo
                      !return
                      !end select
                      !     wasn't there, then try iat ...
                      !select case(jshtyp)
                      !case(0) ! d-s
                    else if (jshtyp.eq.0) then
                      !ss
                      do jj=1,6
                        sspher=0
                        do m=1,6
                          sspher=sspher+trafo(m,jj)*ss(1,m)
                        enddo
                        s2(1,jj)=sspher
                      enddo
                      ss(1,1:5) = s2(1,2:6)

                      !dd
                      do k = 1,3
                        do jj=1,6
                          sspher=0
                          do m=1,6
                            sspher=sspher+trafo(m,jj)*dd(1,m,k)
                          enddo
                          s2(1,jj)=sspher
                        enddo
                        dd(1,1:5,k) = s2(1,2:6)
                      enddo

                      !qq
                      do k = 1,6
                        do jj=1,6
                          sspher=0
                          do m=1,6
                            sspher=sspher+trafo(m,jj)*qq(1,m,k)
                          enddo
                          s2(1,jj)=sspher
                        enddo
                        qq(1,1:5,k) = s2(1,2:6)
                      enddo
                      !return
                      !case(1) ! d-p
                    else if (jshtyp.eq.1) then

                      !ss
                      do ii=1,3
                        do jj=1,6
                          sspher=0
                          do m=1,6
                            sspher=sspher+trafo(m,jj)*ss(ii,m)
                          enddo
                          s2(ii,jj)=sspher
                        enddo
                        ss(ii,1:5) = s2(ii,2:6)
                      enddo

                      !dd
                      do k = 1,3
                        do ii=1,3
                          do jj=1,6
                            sspher=0
                            do m=1,6
                              sspher=sspher+trafo(m,jj)*dd(ii,m,k)
                            enddo
                            s2(ii,jj)=sspher
                          enddo
                          dd(ii,1:5,k) = s2(ii,2:6)
                        enddo
                      enddo

                      !qq
                      do k = 1,6
                        do ii=1,3
                          do jj=1,6
                            sspher=0
                            do m=1,6
                              sspher=sspher+trafo(m,jj)*qq(ii,m,k)
                            enddo
                            s2(ii,jj)=sspher
                          enddo
                          qq(ii,1:5,k) = s2(ii,2:6)
                        enddo
                      enddo

                    else
                      !  return
                      !end select
                      !     if not returned up to here -> d-d
                      ! CB: transposing s in first dgemm is important for integrals other than S

                      !CALL mctc_gemm(trafo, s, dum, transa='T')

                      !ss -----------
                      do i=1,6
                        do j=1,6
                          sspher = 0.0_wp
                          do k=1,6
                            sspher = sspher + trafo(k,i) * ss(k,j)
                          enddo
                          dum(i,j) = sspher
                        enddo
                      enddo

                      !CALL mctc_gemm(dum, trafo, s2, transb='N')

                      do i=1,6
                        do j=1,6
                          sspher = 0.0_wp
                          do k=1,6
                            sspher = sspher + dum(i,k) * trafo(k,j)
                          enddo
                          s2(i,j) = sspher
                        enddo
                      enddo

                      do i=1,5
                        do j=1,5
                          ss(j,i) = s2(j+1,i+1)
                        enddo
                      enddo
                      ! end ss -----------

                      !dd --------
                      do kk = 1,3
                        do i=1,6
                          do j=1,6
                            sspher = 0.0_wp
                            do k=1,6
                              sspher = sspher + trafo(k,i) * dd(k,j,kk)
                            enddo
                            dum(i,j) = sspher
                          enddo
                        enddo

                        !CALL mctc_gemm(dum, trafo, s2, transb='N')

                        do i=1,6
                          do j=1,6
                            sspher = 0.0_wp
                            do k=1,6
                              sspher = sspher + dum(i,k) * trafo(k,j)
                            enddo
                            s2(i,j) = sspher
                          enddo
                        enddo

                        do i=1,5
                          do j=1,5
                            dd(j,i,kk) = s2(j+1,i+1)
                          enddo
                        enddo
                      enddo
                      !end dd --------------

                      !qq ---------------
                      do kk = 1,6
                        do i=1,6
                          do j=1,6
                            sspher = 0.0_wp
                            do k=1,6
                              sspher = sspher + trafo(k,i) * qq(k,j,kk)
                            enddo
                            dum(i,j) = sspher
                          enddo
                        enddo

                        !CALL mctc_gemm(dum, trafo, s2, transb='N')

                        do i=1,6
                          do j=1,6
                            sspher = 0.0_wp
                            do k=1,6
                              sspher = sspher + dum(i,k) * trafo(k,j)
                            enddo
                            s2(i,j) = sspher
                          enddo
                        enddo

                        do i=1,5
                          do j=1,5
                            qq(j,i,kk) = s2(j+1,i+1)
                          enddo
                        enddo
                      enddo
                      !end qq -----------------

                    endif
                  endif

                  ioff = saoshell(ish,iat)
                  joff = saoshell(jsh,jat)

                  ! would be good to change data order of sint,dpint,qpint for 
                  ! memory coalescence
                  !$acc loop vector collapse(2)
                  do ii = 1,llao2(ishtyp)
                     do jj = 1,llao2(jshtyp)
                       iao = ii + ioff
                       jao = jj + joff
                       !$acc atomic
                       sint(iao, jao) = sint(iao, jao) + ss(jj, ii)
                       !$acc atomic
                       sint(jao, iao) = sint(jao, iao) + ss(jj, ii)
                       do k = 1,3
                         !$acc atomic
                         dpint(k, iao, jao) = dpint(k, iao, jao) + dd(jj, ii, k)
                         !$acc atomic
                         dpint(k, jao, iao) = dpint(k, jao, iao) + dd(jj, ii, k)
                       enddo
                       do k = 1,6
                         !$acc atomic
                         qpint(k, iao, jao) = qpint(k, iao, jao) + qq(jj, ii, k)
                         !$acc atomic
                         qpint(k, jao, iao) = qpint(k, jao, iao) + qq(jj, ii, k)
                       enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$acc end parallel

   !$acc exit data copyout(sint(:, :), dpint(:, :, :), qpint(:, :, :))

   ! dd array is transposed in ACC region, so dd2 is created for the diagonal
   ! elements, which is still running on CPU

   ! diagonal elements
   do iat = 1, nat
      ra = xyz(:, iat)
      izp = at(iat)
      do ish = 1, nShell(izp)
         ishtyp = angShell(ish,izp)
         do iao = 1, llao2(ishtyp)
            i = iao+saoshell(ish,iat)
            ii = i*(1+i)/2
            sint(i,i) = 1.0_wp + sint(i,i)
         end do

         icao = caoshell(ish,iat)
         naoi = llao(ishtyp)
         iptyp = itt(ishtyp)
         do jsh = 1, ish
            jshtyp = angShell(jsh,izp)
            jcao = caoshell(jsh,iat)
            naoj = llao(jshtyp)
            jptyp = itt(jshtyp)
            ss = 0.0_wp
            dd2 = 0.0_wp
            qq = 0.0_wp
            call get_multiints(icao,jcao,naoi,naoj,ishtyp,jshtyp,ra,ra,point, &
               &               intcut,nprim,primcount,alp,cont,ss,dd2,qq)
            !transform from CAO to SAO
            !call dtrf2(ss,ishtyp,jshtyp)
            do k = 1,3
               tmp(1:6, 1:6) = dd2(k,1:6, 1:6)
               call dtrf2(tmp, ishtyp, jshtyp)
               dd2(k, 1:6, 1:6) = tmp(1:6, 1:6)
            enddo
            do k = 1,6
               tmp(1:6, 1:6) = qq(k, 1:6, 1:6)
               call dtrf2(tmp, ishtyp, jshtyp)
               qq(k, 1:6, 1:6) = tmp(1:6, 1:6)
            enddo
            do ii = 1, llao2(ishtyp)
               iao = ii + saoshell(ish,iat)
               do jj = 1, llao2(jshtyp)
                  jao = jj + saoshell(jsh,iat)
                  if (jao > iao .and. ish ==  jsh) cycle
                  dpint(1:3, iao, jao) = dpint(1:3, iao, jao) + dd2(1:3, jj, ii)
                  if (iao /= jao) then
                     dpint(1:3, jao, iao) = dpint(1:3, jao, iao) + dd2(1:3, jj, ii)
                  end if
                  qpint(1:6, iao, jao) = qq(1:6, jj, ii)
                  if (jao /= iao) then
                     qpint(1:6, jao, iao) = qq(1:6, jj, ii)
                  end if
               end do
            end do
         end do
      end do
   end do

   !$acc exit data delete(xyz, at, nShell, caoshell, saoshell, &
   !$acc& nprim, primcount, alp, cont, intcut, trans, point, angShell(:, :))

end subroutine sdqint_gpu


end module xtb_intgrad
