!! ------------------------------------------------------------------------
!  reimplementation of the CM5 charges
!! ------------------------------------------------------------------------
subroutine calc_cm5(nat,at,xyz,q,cm5,cm5a,dcm5)
use iso_fortran_env, wp => real64
use mctc_econv
implicit none
integer, intent(in)  :: nat
integer, parameter   :: mz=118
real(wp),intent(in)  :: xyz(3,nat)
integer, intent(in)  :: at(nat)
real(wp),intent(in)  :: q(nat)
real(wp),intent(out) :: cm5(nat)
real(wp),intent(out) :: cm5a(nat)
real(wp),intent(out) :: dcm5(3,nat,nat)
real(wp) :: d(mz,mz)
real(wp) :: dist,bkk,bkkd,bkkda
integer  :: i,j,k
real(wp) :: alp,rab(3)
!
! COVALENT RADII
!
! based on "Atomic Radii of the Elements," M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
! in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
! edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
! corrected Nov. 17, 2010 for the 92nd edition.
!
real(wp),parameter :: rad(mz) = (/ &
   & 0.32_wp,0.37_wp,1.30_wp,0.99_wp,0.84_wp,0.75_wp,0.71_wp,0.64_wp, &
   & 0.60_wp,0.62_wp,1.60_wp,1.40_wp,1.24_wp,1.14_wp,1.09_wp,1.04_wp, &
   & 1.00_wp,1.01_wp,2.00_wp,1.74_wp,1.59_wp,1.48_wp,1.44_wp,1.30_wp, &
   & 1.29_wp,1.24_wp,1.18_wp,1.17_wp,1.22_wp,1.20_wp,1.23_wp,1.20_wp, &
   & 1.20_wp,1.18_wp,1.17_wp,1.16_wp,2.15_wp,1.90_wp,1.76_wp,1.64_wp, &
   & 1.56_wp,1.46_wp,1.38_wp,1.36_wp,1.34_wp,1.30_wp,1.36_wp,1.40_wp, &
   & 1.42_wp,1.40_wp,1.40_wp,1.37_wp,1.36_wp,1.36_wp,2.38_wp,2.06_wp, &
   & 1.94_wp,1.84_wp,1.90_wp,1.88_wp,1.86_wp,1.85_wp,1.83_wp,1.82_wp, &
   & 1.81_wp,1.80_wp,1.79_wp,1.77_wp,1.77_wp,1.78_wp,1.74_wp,1.64_wp, &
   & 1.58_wp,1.50_wp,1.41_wp,1.36_wp,1.32_wp,1.30_wp,1.30_wp,1.32_wp, &
   & 1.44_wp,1.45_wp,1.50_wp,1.42_wp,1.48_wp,1.46_wp,2.42_wp,2.11_wp, &
   & 2.01_wp,1.90_wp,1.84_wp,1.83_wp,1.80_wp,1.80_wp,1.73_wp,1.68_wp, &
   & 1.68_wp,1.68_wp,1.65_wp,1.67_wp,1.73_wp,1.76_wp,1.61_wp,1.57_wp, &
   & 1.49_wp,1.43_wp,1.41_wp,1.34_wp,1.29_wp,1.28_wp,1.21_wp,1.22_wp, &
   & 1.36_wp,1.43_wp,1.62_wp,1.75_wp,1.65_wp,1.57_wp /)*aatoau
!
! cm5 model parameters
!
! atomwise parameters
real(wp),parameter :: a0(mz) = (/ &
   & 0.0056_wp,-0.1543_wp, 0.0000_wp, 0.0333_wp,-0.1030_wp,-0.0446_wp, &
   &-0.1072_wp,-0.0802_wp,-0.0629_wp,-0.1088_wp, 0.0184_wp, 0.0000_wp, &
   &-0.0726_wp,-0.0790_wp,-0.0756_wp,-0.0565_wp,-0.0444_wp,-0.0767_wp, &
   & 0.0130_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &-0.0512_wp,-0.0557_wp,-0.0533_wp,-0.0399_wp,-0.0313_wp,-0.0541_wp, &
   & 0.0092_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &-0.0361_wp,-0.0393_wp,-0.0376_wp,-0.0281_wp,-0.0220_wp,-0.0381_wp, &
   & 0.0065_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   & 0.0000_wp, 0.0000_wp,-0.0255_wp,-0.0277_wp,-0.0265_wp,-0.0198_wp, &
   &-0.0155_wp,-0.0269_wp, 0.0046_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp,-0.0179_wp,-0.0195_wp, &
   &-0.0187_wp,-0.0140_wp,-0.0110_wp,-0.0189_wp /)
!
d = 0.0_wp
do i=1,mz
   do j=i+1,mz
      d(i,j)=a0(i)-a0(j)
   enddo
enddo
! pairwise parameters
d( 1, 6)= 0.0502_wp
d( 1, 7)= 0.1747_wp
d( 1, 8)= 0.1671_wp
d( 6, 7)= 0.0556_wp
d( 6, 8)= 0.0234_wp
d( 7, 8)=-0.0346_wp
!
do i=1,mz
   do j=i+1,mz
      d(j,i)=-d(i,j)
   enddo
enddo
! alpha
alp=2.4740_wp/aatoau
! c-coefficient: 0.7050   ! already included in a0
dcm5 = 0.0_wp
cm5a = 0.0_wp

do k=1,nat
   do i=1,nat
      if (at(k).ne.at(i)) then
         rab=(xyz(:,k)-xyz(:,i))
         dist=norm2(rab)
         bkk=exp(-alp*(dist-rad(at(k))-rad(at(i))))
         bkkd=bkk*d(at(k),at(i))
         cm5a(k)=cm5a(k)+bkkd
         bkkda=bkkd*alp/dist
         dcm5(:,i,k)=bkkda*rab
      endif
   enddo
enddo
do k = 1, nat
   cm5(k)=q(k)+cm5a(k)
enddo

end subroutine calc_cm5
