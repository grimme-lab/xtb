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

module approxrab
  use iso_fortran_env, wp => real64
  implicit none

! parameter blocks
  real(wp),private, dimension(86) :: cnfak
  real(wp),private, dimension(86) :: r0
  real(wp),private, dimension(86) :: en
  real(wp),private, dimension(4,2):: p
! START PARAMETER--------------------------------------------------
  data en /&
       2.30085633_wp, 2.78445145_wp, 1.52956084_wp, 1.51714704_wp, 2.20568300_wp,&
       2.49640820_wp, 2.81007174_wp, 4.51078438_wp, 4.67476223_wp, 3.29383610_wp,&
       2.84505365_wp, 2.20047950_wp, 2.31739628_wp, 2.03636974_wp, 1.97558064_wp,&
       2.13446570_wp, 2.91638164_wp, 1.54098156_wp, 2.91656301_wp, 2.26312147_wp,&
       2.25621439_wp, 1.32628677_wp, 2.27050569_wp, 1.86790977_wp, 2.44759456_wp,&
       2.49480042_wp, 2.91545568_wp, 3.25897750_wp, 2.68723778_wp, 1.86132251_wp,&
       2.01200832_wp, 1.97030722_wp, 1.95495427_wp, 2.68920990_wp, 2.84503857_wp,&
       2.61591858_wp, 2.64188286_wp, 2.28442252_wp, 1.33011187_wp, 1.19809388_wp,&
       1.89181390_wp, 2.40186898_wp, 1.89282464_wp, 3.09963488_wp, 2.50677823_wp,&
       2.61196704_wp, 2.09943450_wp, 2.66930105_wp, 1.78349472_wp, 2.09634533_wp,&
       2.00028974_wp, 1.99869908_wp, 2.59072029_wp, 2.54497829_wp, 2.52387890_wp,&
       2.30204667_wp, 1.60119300_wp, 2.00000000_wp, 2.00000000_wp, 2.00000000_wp,&
       2.00000000_wp, 2.00000000_wp, 2.00000000_wp, 2.00000000_wp, 2.00000000_wp,&
       2.00000000_wp, 2.00000000_wp, 2.00000000_wp, 2.00000000_wp, 2.00000000_wp,&
       2.00000000_wp, 2.30089349_wp, 1.75039077_wp, 1.51785130_wp, 2.62972945_wp,&
       2.75372921_wp, 2.62540906_wp, 2.55860939_wp, 3.32492356_wp, 2.65140898_wp,&
       1.52014458_wp, 2.54984804_wp, 1.72021963_wp, 2.69303422_wp, 1.81031095_wp,&
       2.34224386_wp&
       /
  data r0 /&
       0.55682207_wp, 0.80966997_wp, 2.49092101_wp, 1.91705642_wp, 1.35974851_wp,&
       0.98310699_wp, 0.98423007_wp, 0.76716063_wp, 1.06139799_wp, 1.17736822_wp,&
       2.85570926_wp, 2.56149012_wp, 2.31673425_wp, 2.03181740_wp, 1.82568535_wp,&
       1.73685958_wp, 1.97498207_wp, 2.00136196_wp, 3.58772537_wp, 2.68096221_wp,&
       2.23355957_wp, 2.33135502_wp, 2.15870365_wp, 2.10522128_wp, 2.16376162_wp,&
       2.10804037_wp, 1.96460045_wp, 2.00476257_wp, 2.22628712_wp, 2.43846700_wp,&
       2.39408483_wp, 2.24245792_wp, 2.05751204_wp, 2.15427677_wp, 2.27191920_wp,&
       2.19722638_wp, 3.80910350_wp, 3.26020971_wp, 2.99716916_wp, 2.71707818_wp,&
       2.34950167_wp, 2.11644818_wp, 2.47180659_wp, 2.32198800_wp, 2.32809515_wp,&
       2.15244869_wp, 2.55958313_wp, 2.59141300_wp, 2.62030465_wp, 2.39935278_wp,&
       2.56912355_wp, 2.54374096_wp, 2.56914830_wp, 2.53680807_wp, 4.24537037_wp,&
       3.66542289_wp, 3.19903011_wp, 2.80000000_wp, 2.80000000_wp, 2.80000000_wp,&
       2.80000000_wp, 2.80000000_wp, 2.80000000_wp, 2.80000000_wp, 2.80000000_wp,&
       2.80000000_wp, 2.80000000_wp, 2.80000000_wp, 2.80000000_wp, 2.80000000_wp,&
       2.80000000_wp, 2.34880037_wp, 2.37597108_wp, 2.49067697_wp, 2.14100577_wp,&
       2.33473532_wp, 2.19498900_wp, 2.12678348_wp, 2.34895048_wp, 2.33422774_wp,&
       2.86560827_wp, 2.62488837_wp, 2.88376127_wp, 2.75174124_wp, 2.83054552_wp,&
       2.63264944_wp&
       /
  data cnfak /&
       0.17957827_wp, 0.25584045_wp,-0.02485871_wp, 0.00374217_wp, 0.05646607_wp,&
       0.10514203_wp, 0.09753494_wp, 0.30470380_wp, 0.23261783_wp, 0.36752208_wp,&
       0.00131819_wp,-0.00368122_wp,-0.01364510_wp, 0.04265789_wp, 0.07583916_wp,&
       0.08973207_wp,-0.00589677_wp, 0.13689929_wp,-0.01861307_wp, 0.11061699_wp,&
       0.10201137_wp, 0.05426229_wp, 0.06014681_wp, 0.05667719_wp, 0.02992924_wp,&
       0.03764312_wp, 0.06140790_wp, 0.08563465_wp, 0.03707679_wp, 0.03053526_wp,&
      -0.00843454_wp, 0.01887497_wp, 0.06876354_wp, 0.01370795_wp,-0.01129196_wp,&
       0.07226529_wp, 0.01005367_wp, 0.01541506_wp, 0.05301365_wp, 0.07066571_wp,&
       0.07637611_wp, 0.07873977_wp, 0.02997732_wp, 0.04745400_wp, 0.04582912_wp,&
       0.10557321_wp, 0.02167468_wp, 0.05463616_wp, 0.05370913_wp, 0.05985441_wp,&
       0.02793994_wp, 0.02922983_wp, 0.02220438_wp, 0.03340460_wp,-0.04110969_wp,&
      -0.01987240_wp, 0.07260201_wp, 0.07700000_wp, 0.07700000_wp, 0.07700000_wp,&
       0.07700000_wp, 0.07700000_wp, 0.07700000_wp, 0.07700000_wp, 0.07700000_wp,&
       0.07700000_wp, 0.07700000_wp, 0.07700000_wp, 0.07700000_wp, 0.07700000_wp,&
       0.07700000_wp, 0.08379100_wp, 0.07314553_wp, 0.05318438_wp, 0.06799334_wp,&
       0.04671159_wp, 0.06758819_wp, 0.09488437_wp, 0.07556405_wp, 0.13384502_wp,&
       0.03203572_wp, 0.04235009_wp, 0.03153769_wp,-0.00152488_wp, 0.02714675_wp,&
       0.04800662_wp&
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
      do m = 1, n
         drabdr(:,m,k)=ff*(cnfak(ati)*dcndr(:,m,i)&
            &            + cnfak(atj)*dcndr(:,m,j))
      enddo ! m
      drabdL(:,:,k) = ff*(cnfak(ati)*dcndL(:,:,i)&
         &              + cnfak(atj)*dcndL(:,:,j))
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

end module approxrab

