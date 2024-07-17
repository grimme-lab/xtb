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

module xtb_approxrab
   use xtb_mctc_accuracy, only : wp
  implicit none
  private
  public :: pbc_approx_rab, approx_rab, approx_bonds

! parameter blocks
  real(wp),private, dimension(103) :: cnfak
  real(wp),private, dimension(103) :: r0
  real(wp),private, dimension(103) :: en
  real(wp),private, dimension(4,2):: p
! START PARAMETER--------------------------------------------------
      data en /&
      2.30085633, 2.78445145, 1.52956084, 1.51714704, 2.20568300,&
      2.49640820, 2.81007174, 4.51078438, 4.67476223, 3.29383610,&
      2.84505365, 2.20047950, 2.31739628, 2.03636974, 1.97558064,&
      2.13446570, 2.91638164, 1.54098156, 2.91656301, 2.26312147,&
      2.25621439, 1.32628677, 2.27050569, 1.86790977, 2.44759456,&
      2.49480042, 2.91545568, 3.25897750, 2.68723778, 1.86132251,&
      2.01200832, 1.97030722, 1.95495427, 2.68920990, 2.84503857,&
      2.61591858, 2.64188286, 2.28442252, 1.33011187, 1.19809388,&
      1.89181390, 2.40186898, 1.89282464, 3.09963488, 2.50677823,&
      2.61196704, 2.09943450, 2.66930105, 1.78349472, 2.09634533,&
      2.00028974, 1.99869908, 2.59072029, 2.54497829, 2.52387890,&
      2.30204667, 1.60119300, 2.00000000, 2.00000000, 2.00000000,& ! 60
      2.00000000, 2.00000000, 2.00000000, 2.00000000, 2.00000000,&
      2.00000000, 2.00000000, 2.00000000, 2.00000000, 2.00000000,& ! 70
      2.00000000, 2.30089349, 1.75039077, 1.51785130, 2.62972945,&
      2.75372921, 2.62540906, 2.55860939, 3.32492356, 2.65140898,& ! 80
      1.52014458, 2.54984804, 1.72021963, 2.69303422, 1.81031095,&
      2.34224386,&
      2.52387890,& !@thomas TODO
      2.30204667, 1.60119300, 2.00000000, 2.00000000, 2.00000000,& ! 92
      2.00000000, 2.00000000, 2.00000000, 2.00000000, 2.00000000,&
      2.00000000, 2.00000000, 2.00000000, 2.00000000, 2.00000000,& ! 102
      2.00000000&
      /
      data r0 /&
      0.55682207, 0.80966997, 2.49092101, 1.91705642, 1.35974851,&
      0.98310699, 0.98423007, 0.76716063, 1.06139799, 1.17736822,& ! 10
      2.85570926, 2.56149012, 2.31673425, 2.03181740, 1.82568535,&
      1.73685958, 1.97498207, 2.00136196, 3.58772537, 2.68096221,& ! 20
      2.23355957, 2.33135502, 2.15870365, 2.10522128, 2.16376162,&
      2.10804037, 1.96460045, 2.00476257, 2.22628712, 2.43846700,& ! 30
      2.39408483, 2.24245792, 2.05751204, 2.15427677, 2.27191920,&
      2.19722638, 3.80910350, 3.26020971, 2.99716916, 2.71707818,& ! 40
      2.34950167, 2.11644818, 2.47180659, 2.32198800, 2.32809515,&
      2.15244869, 2.55958313, 2.59141300, 2.62030465, 2.39935278,& ! 50
      2.56912355, 2.54374096, 2.56914830, 2.53680807, 4.24537037,&
      3.66542289, 3.22480000, 3.21280000, 3.10550000, 3.10200000,& ! 60
      3.10840000, 3.14030000, 3.06390000, 3.10730000, 3.10000000,&
      3.11910000, 3.10760000, 3.13740000, 3.09740000, 2.92860000,&
      3.05880000, 2.34880037, 2.37597108, 2.49067697, 2.14100577,&
      2.33473532, 2.19498900, 2.12678348, 2.34895048, 2.33422774,&
      2.86560827, 2.62488837, 2.88376127, 2.75174124, 2.83054552,&
      2.63264944,&
      4.24537037,& !@thomas TODO
      3.66542289, 4.20000000, 4.20000000, 4.20000000, 4.20000000,& ! 92
      4.20000000, 4.20000000, 4.20000000, 4.20000000, 4.20000000,&
      4.20000000, 4.20000000, 4.20000000, 4.20000000, 4.20000000,& ! 102
      4.20000000&
      /
      data cnfak /&
      0.17957827, 0.25584045,-0.02485871, 0.00374217, 0.05646607,&
      0.10514203, 0.09753494, 0.30470380, 0.23261783, 0.36752208,&
      0.00131819,-0.00368122,-0.01364510, 0.04265789, 0.07583916,&
      0.08973207,-0.00589677, 0.13689929,-0.01861307, 0.11061699,&
      0.10201137, 0.05426229, 0.06014681, 0.05667719, 0.02992924,&
      0.03764312, 0.06140790, 0.08563465, 0.03707679, 0.03053526,&
     -0.00843454, 0.01887497, 0.06876354, 0.01370795,-0.01129196,&
      0.07226529, 0.01005367, 0.01541506, 0.05301365, 0.07066571,&
      0.07637611, 0.07873977, 0.02997732, 0.04745400, 0.04582912,&
      0.10557321, 0.02167468, 0.05463616, 0.05370913, 0.05985441,& ! 60
      0.02793994, 0.02922983, 0.02220438, 0.03340460,-0.04110969,&
     -0.01987240, 0.07260201, 0.07700000, 0.07700000, 0.07700000,& ! 70
      0.07700000, 0.07700000, 0.07700000, 0.07700000, 0.07700000,&
      0.07700000, 0.07700000, 0.07700000, 0.07700000, 0.07700000,& ! 80
      0.07700000, 0.08379100, 0.07314553, 0.05318438, 0.06799334,&
      0.04671159, 0.06758819, 0.09488437, 0.07556405, 0.13384502,& ! 90
      0.03203572, 0.04235009, 0.03153769,-0.00152488, 0.02714675,&
      0.04800662,&
      0.04582912,& !@thomas TODO
      0.10557321, 0.02167468, 0.05463616, 0.05370913, 0.05985441,& ! 92
      0.02793994, 0.02922983, 0.02220438, 0.03340460,-0.04110969,&
     -0.01987240, 0.07260201, 0.07700000, 0.07700000, 0.07700000,& ! 102
      0.07700000&
     /
! END PARAMETER-----------------------------------------------------

! global EN polynomial parameter (NOTE: x 10^3)
   data p /29.84522887_wp, -1.70549806_wp,  6.54013762_wp,  6.39169003_wp, &
           -8.87843763_wp,  2.10878369_wp,  0.08009374_wp, -0.85808076_wp/

contains

pure subroutine pbc_approx_rab(n,at,xyz,cn,dcndr,dcndL,nsrb,srblist,shift, &
      &                        rab,drabdr,drabdL)
   implicit none
   ! intent in
   integer,                    intent(in) :: n
   integer,                    intent(in) :: nsrb
   integer, dimension(n),      intent(in) :: at
   integer, dimension(2,nsrb), intent(in) :: srblist
   real(wp), dimension(3,n),   intent(in) :: xyz
   real(wp), dimension(n),     intent(in) :: cn
   real(wp), dimension(3,n,n), intent(in) :: dcndr
   real(wp), dimension(3,3,n), intent(in) :: dcndL
   real(wp), intent(in) :: shift
   ! intent inout
   real(wp), dimension(nsrb),     intent(out) :: rab
   real(wp), dimension(3,n,nsrb), intent(out) :: drabdr
   real(wp), dimension(3,3,nsrb), intent(out) :: drabdL
   ! local variable
   integer :: i,j,k,m
   integer :: iat,jat,ati,atj
   integer :: ir,jr
   real(wp) :: ra,rb
   real(wp) :: den
   real(wp) :: k1,k2
   real(wp) :: ff

   rab  = 0.0_wp
   drabdr = 0.0_wp
   drabdL = 0.0_wp

   do concurrent(k = 1:nsrb)
      ! enroll srblist
      i = srblist(1,k)
      j = srblist(2,k)
      ati = at(i)
      atj = at(j)
      ir = itr(ati)
      jr = itr(atj)
      ra=r0(ati)+cnfak(ati)*cn(i)+shift
      rb=r0(atj)+cnfak(atj)*cn(j)+shift
      den=abs(en(ati) - en(atj))
      k1=0.005_wp*(p(ir,1) + p(jr,1))
      k2=0.005_wp*(p(ir,2) + p(jr,2))
      ff=1.0_wp - k1*den - k2*den**2
      ! save distances => rab
      rab(k) = (ra + rb)*ff
      ! save gradient => drabdr
      drabdr(:,:,k) = ff*(cnfak(ati)*dcndr(:,:,i) + cnfak(atj)*dcndr(:,:,j))
      drabdL(:,:,k) = ff*(cnfak(ati)*dcndL(:,:,i) + cnfak(atj)*dcndL(:,:,j))
   enddo ! k

end subroutine pbc_approx_rab

pure subroutine approx_bonds(n,at,xyz,cn,bonds,shift)
   implicit none
   ! intent in
   integer,                  intent(in) :: n
   integer,  dimension(n),   intent(in) :: at
   real(wp), dimension(3,n), intent(in) :: xyz
   real(wp), dimension(n),   intent(in) :: cn
   real(wp),                 intent(in) :: shift
   ! intent inout
   integer,  dimension(n,n), intent(out) :: bonds
   ! local variable
   integer :: i,j,k,m
   integer :: iat,jat,ati,atj
   integer :: ir,jr
   real(wp) :: ra,rb
   real(wp) :: den
   real(wp) :: k1,k2
   real(wp) :: ff,tmp
   real(wp) :: r2,rab,rij(3)

   bonds = 0

   do i = 1, n
      bonds(i,i) = ceiling(cn(i))
      do j = 1, i-1
         rij = xyz(:,i) - xyz(:,j)
         r2 = sum(rij**2)
         if (r2.gt.200.0_wp) cycle
         ati = at(i)
         atj = at(j)
         ir = itr(ati)
         jr = itr(atj)
         ra=r0(ati)+cnfak(ati)*cn(i)+shift
         rb=r0(atj)+cnfak(atj)*cn(j)+shift
         den=abs(en(ati) - en(atj))
         k1=0.005_wp*(p(ir,1) + p(jr,1))
         k2=0.005_wp*(p(ir,2) + p(jr,2))
         ff=1.0_wp - k1*den - k2*den**2
         ! save distances => rab
         rab = ((ra + rb)*ff)**2
         if (r2 < rab) then
            bonds(i,j) = 1
            bonds(j,i) = 1
         endif
      enddo ! k
   enddo ! k

end subroutine approx_bonds

pure subroutine approx_rab(n,at,xyz,cn,dcndr,nsrb,srblist,shift,rab,grab)
   implicit none
   ! intent in
   integer,                    intent(in) :: n
   integer,                    intent(in) :: nsrb
   integer, dimension(n),      intent(in) :: at
   integer, dimension(2,nsrb), intent(in) :: srblist
   real(wp), dimension(3,n),   intent(in) :: xyz
   real(wp), dimension(n),     intent(in) :: cn
   real(wp), dimension(3,n,n), intent(in) :: dcndr
   real(wp), intent(in) :: shift
   ! intent inout
   real(wp), dimension(nsrb),     intent(out) :: rab
   real(wp), dimension(3,n,nsrb), intent(out) :: grab
   ! local variable
   integer :: i,j,k,m
   integer :: iat,jat,ati,atj
   integer :: ir,jr
   real(wp) :: ra,rb
   real(wp) :: den
   real(wp) :: k1,k2
   real(wp) :: ff

   rab  = 0.0_wp
   grab = 0.0_wp

   do concurrent(k = 1:nsrb)
      ! enroll srblist
      i = srblist(1,k)
      j = srblist(2,k)
      ati = at(i)
      atj = at(j)
      ir = itr(ati)
      jr = itr(atj)
      ra=r0(ati)+cnfak(ati)*cn(i)+shift
      rb=r0(atj)+cnfak(atj)*cn(j)+shift
      den=abs(en(ati) - en(atj))
      k1=0.005_wp*(p(ir,1) + p(jr,1))
      k2=0.005_wp*(p(ir,2) + p(jr,2))
      ff=1.0_wp - k1*den - k2*den**2
      ! save distances => rab
      rab(k) = (ra + rb)*ff
      ! save gradient => grab
      do m = 1, n
         grab(:,m,k)=ff*(cnfak(ati)*dcndr(:,m,i)&
            &          + cnfak(atj)*dcndr(:,m,j))
      enddo ! m
   enddo ! k

end subroutine approx_rab

! row in PSE for given ordinal number; (values>4) => 4
pure elemental integer function itr(i)
   implicit none
   integer, intent(in) :: i
   itr = 0
   if(i.gt.0.and.i.le.2) then
      itr = 1
   elseif(i.gt.2.and.i.le.10) then
      itr = 2
   elseif(i.gt.10.and.i.le.18) then
      itr = 3
   elseif(i.gt.18) then
      itr = 4
   endif

   return
end function itr

end module xtb_approxrab

