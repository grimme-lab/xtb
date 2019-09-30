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

! SAW ------------------------------------------------------------- 1803
!     STO-NG by R. F. Stewart, J. Chem. Phys., 52 431-438, 1970.
! ng: number of primitive gaussians per slater function
! n:  principal quantum number (max. 5)
! l:  angular quantum number (s=0, p=1, d=2, f=3, g=4)
! alpha/coeff: exponents and coeffients of primitive gaussians

!* this is just a driver for the STO-NG routines
pure subroutine slater(ng,n,l,zeta,alpha,coeff,info)
   implicit none
   interface
      pure subroutine sto1g(ityp,zeta,alpha,coeff)
      integer, parameter   :: dp=kind(1.0d0)
      integer, intent(in)  :: ityp
      real(dp),intent(in)  :: zeta
      real(dp),intent(out) :: alpha(*)
      real(dp),intent(out) :: coeff(*)
      end subroutine sto1g
      pure subroutine sto2g(ityp,zeta,alpha,coeff)
      integer, parameter   :: dp=kind(1.0d0)
      integer, intent(in)  :: ityp
      real(dp),intent(in)  :: zeta
      real(dp),intent(out) :: alpha(*)
      real(dp),intent(out) :: coeff(*)
      end subroutine sto2g
      pure subroutine sto3g(ityp,zeta,alpha,coeff)
      integer, parameter   :: dp=kind(1.0d0)
      integer, intent(in)  :: ityp
      real(dp),intent(in)  :: zeta
      real(dp),intent(out) :: alpha(*)
      real(dp),intent(out) :: coeff(*)
      end subroutine sto3g
      pure subroutine sto4g(ityp,zeta,alpha,coeff)
      integer, parameter   :: dp=kind(1.0d0)
      integer, intent(in)  :: ityp
      real(dp),intent(in)  :: zeta
      real(dp),intent(out) :: alpha(*)
      real(dp),intent(out) :: coeff(*)
      end subroutine sto4g
      pure subroutine sto5g(ityp,zeta,alpha,coeff)
      integer, parameter   :: dp=kind(1.0d0)
      integer, intent(in)  :: ityp
      real(dp),intent(in)  :: zeta
      real(dp),intent(out) :: alpha(*)
      real(dp),intent(out) :: coeff(*)
      end subroutine sto5g
      pure subroutine sto6g(ityp,zeta,alpha,coeff)
      integer, parameter   :: dp=kind(1.0d0)
      integer, intent(in)  :: ityp
      real(dp),intent(in)  :: zeta
      real(dp),intent(out) :: alpha(*)
      real(dp),intent(out) :: coeff(*)
      end subroutine sto6g
   end interface
   integer, parameter   :: dp=kind(1.0d0)
   integer, intent(in)  :: ng
   integer, intent(in)  :: n
   integer, intent(in)  :: l
   real(dp),intent(in)  :: zeta
   real(dp),intent(out) :: alpha(*)
   real(dp),intent(out) :: coeff(*)
   integer, intent(out) :: info
   integer :: ityp

   ! basic checks first, we cannot support principal QN of six or higher
   ! also you should not violate l ∊ [n-1, n-2, ..., 1, 0]
   if ((n.gt.5).or.(n.le.l)) then
      info = 2
      return
   endif
   ! Who would actually try this? Better save, than sorry.
   if (zeta.le.0.0_dp) then
      info = 4
      return
   endif

   ! we have to use a little hack here, 
   ! if you pass n and l correctly, everything is fine
   ! ityp: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 
   !    n: 1 2 3 4 5 2 3 4 5  3  4  5  4  5  5
   !    l: 0 0 0 0 0 1 1 1 1  2  2  2  3  3  4
   select case(l) ! integer hack:
   case(0); ityp = n    ! s
   case(1); ityp = 4+n  ! p
   case(2); ityp = 7+n  ! d
   case(3); ityp = 9+n  ! f
   case(4); ityp = 10+n ! g
   case default ! I'm sorry, no h-functions for you today
      info = 3
      return
   end select

   select case(ng) ! select right STO-NG representation
   case(1); call sto1g(ityp,zeta,alpha,coeff)
   case(2); call sto2g(ityp,zeta,alpha,coeff)
   case(3); call sto3g(ityp,zeta,alpha,coeff)
   case(4); call sto4g(ityp,zeta,alpha,coeff)
   case(5); call sto5g(ityp,zeta,alpha,coeff)
   case(6); call sto6g(ityp,zeta,alpha,coeff)
   case default ! currently we cannot go beyond 6 primitives
      info = 1
      return
   end select

   info = 0 ! success, hurray

end subroutine slater

pure subroutine sto1g(ityp,zeta,alpha,coeff)
   implicit none
   integer, parameter   :: dp=kind(1.0d0)
   integer, intent(in)  :: ityp
   real(dp),intent(in)  :: zeta
   real(dp),intent(out) :: alpha(*)
   real(dp),intent(out) :: coeff(*)

   integer,parameter :: ng=1,nf=15,l(*)=(/0,0,0,0,0,1,1,1,1,2,2,2,3,3,4/)
   real(dp),parameter :: top = 0.5_dp/atan(1.0_dp) ! 0.636619772367582_dp
   real(dp),parameter :: dfactorial(*) = & ! see OEIS A001147
   &  (/1._dp,1._dp,3._dp,15._dp,105._dp,945._dp,10395._dp,135135._dp/)
   real(dp),parameter :: palpha(ng,nf) = reshape((/ &
!* 1s
   2.709498091e-1_dp, &
!* 2s
   1.012151084e-1_dp, &
!* 3s
   5.296881757e-2_dp, &
!* 4s
   3.264600274e-2_dp, &
!* 5s
   2.216912938e-2_dp, &
!* 2p
   1.759666885e-1_dp, &
!* 3p
   9.113614253e-2_dp, &
!* 4p
   5.578350235e-2_dp, &
!* 5p
   3.769845216e-2_dp, &
!* 3d
   1.302270363e-1_dp, &
!* 4d
   7.941656339e-2_dp, &
!* 5d
   5.352200793e-2_dp, &
!* 4f
   1.033434062e-1_dp, &
!* 5f
   6.952785407e-2_dp, &
!* 5g
   8.565417784e-2_dp /),shape(palpha))
   
!* <φ|φ> = (2i-1)!!(2j-1)!!(2k-1)!!/(4α)^(i+j+k) · sqrt(π/2α)³
!  N² = (4α)^(i+j+k)/((2i-1)!!(2j-1)!!(2k-1)!!)  · sqrt(2α/π)³
!  N = (4α)^((i+j+k)/2) / sqrt((2i-1)!!(2j-1)!!(2k-1)!!) · (2α/π)^(3/4)
   alpha(:ng) = palpha(:,ityp) * zeta**2
   coeff(:ng) = (top*alpha(:ng))**0.75_dp &
   &  * (4*alpha(:ng))**(l(ityp)) / sqrt(dfactorial(l(ityp)+1))

end subroutine sto1g

pure subroutine sto2g(ityp,zeta,alpha,coeff)
   implicit none
   integer, parameter   :: dp=kind(1.0d0)
   integer, intent(in)  :: ityp
   real(dp),intent(in)  :: zeta
   real(dp),intent(out) :: alpha(*)
   real(dp),intent(out) :: coeff(*)

   integer,parameter :: ng=2,nf=15,l(*)=(/0,0,0,0,0,1,1,1,1,2,2,2,3,3,4/)
   real(dp),parameter :: top = 0.5_dp/atan(1.0_dp) ! 0.636619772367582_dp
   real(dp),parameter :: dfactorial(*) = & ! see OEIS A001147
   &  (/1._dp,1._dp,3._dp,15._dp,105._dp,945._dp,10395._dp,135135._dp/)
   real(dp),parameter :: palpha(ng,nf) = reshape((/ &
!* 1s
   8.518186635e-1_dp,1.516232927e-1_dp, &
!* 2s
   1.292278611e-1_dp,4.908584205e-2_dp, &
!* 3s
   6.694095822e-1_dp,5.837135094e-2_dp, &
!* 4s
   2.441785453e-1_dp,4.051097664e-2_dp, &
!* 5s
   1.213425654e-1_dp,3.133152144e-2_dp, &
!* 2p
   4.323908358e-1_dp,1.069139065e-1_dp, &
!* 3p
   1.458620964e-1_dp,5.664210742e-2_dp, &
!* 4p
   6.190052680e-2_dp,2.648418407e-2_dp, &
!* 5p
   2.691294191e-1_dp,3.980805011e-2_dp, &
!* 3d
   2.777427345e-1_dp,8.336507714e-2_dp, &
!* 4d
   1.330958892e-1_dp,5.272119659e-2_dp, &
!* 5d
   6.906014388e-2_dp,3.399457777e-2_dp, &
!* 4f
   2.006693538e-1_dp,6.865384900e-2_dp, &
!* 5f
   1.156094555e-1_dp,4.778940916e-2_dp, &
!* 5g
   1.554531559e-1_dp,5.854079811e-2_dp /),shape(palpha))
   real(dp),parameter :: pcoeff(ng,nf) = reshape((/ &
!* 1s
   4.301284983e-1_dp,6.789135305e-1_dp, &
!* 2s
   7.470867124e-1_dp,2.855980556e-1_dp, &
!* 3s
  -1.529645716e-1_dp,1.051370110e+0_dp, &
!* 4s
  -3.046656896e-1_dp,1.146877294e+0_dp, &
!* 5s
  -5.114756049e-1_dp,1.307377277e+0_dp, &
!* 2p
   4.522627513e-1_dp,6.713122642e-1_dp, &
!* 3p
   5.349653114e-1_dp,5.299607212e-1_dp, &
!* 4p
   8.743116767e-1_dp,1.513640107e-1_dp, &
!* 5p
  -1.034227010e-1_dp,1.033376378e+0_dp, &
!* 3d
   4.666137923e-1_dp,6.644706516e-1_dp, &
!* 4d
   4.932764167e-1_dp,5.918727866e-1_dp, &
!* 5d
   6.539405185e-1_dp,3.948945302e-1_dp, &
!* 4f
   4.769346276e-1_dp,6.587383976e-1_dp, &
!* 5f
   4.856637346e-1_dp,6.125980914e-1_dp, &
!* 5g
   4.848298074e-1_dp,6.539381621e-1_dp /),shape(pcoeff))
   
!* <φ|φ> = (2i-1)!!(2j-1)!!(2k-1)!!/(4α)^(i+j+k) · sqrt(π/2α)³
!  N² = (4α)^(i+j+k)/((2i-1)!!(2j-1)!!(2k-1)!!)  · sqrt(2α/π)³
!  N = (4α)^((i+j+k)/2) / sqrt((2i-1)!!(2j-1)!!(2k-1)!!) · (2α/π)^(3/4)
   alpha(:ng) = palpha(:,ityp) * zeta**2
   coeff(:ng) = pcoeff(:,ityp) * (top*alpha(:ng))**0.75_dp &
   &  * sqrt(4*alpha(:ng))**(l(ityp)) / sqrt(dfactorial(l(ityp)+1))

end subroutine sto2g

pure subroutine sto3g(ityp,zeta,alpha,coeff)
   implicit none
   integer, parameter   :: dp=kind(1.0d0)
   integer, intent(in)  :: ityp
   real(dp),intent(in)  :: zeta
   real(dp),intent(out) :: alpha(*)
   real(dp),intent(out) :: coeff(*)

   integer,parameter :: ng=3,nf=15,l(*)=(/0,0,0,0,0,1,1,1,1,2,2,2,3,3,4/)
   real(dp),parameter :: top = 0.5_dp/atan(1.0_dp) ! 0.636619772367582_dp
   real(dp),parameter :: dfactorial(*) = & ! see OEIS A001147
   &  (/1._dp,1._dp,3._dp,15._dp,105._dp,945._dp,10395._dp,135135._dp/)
   real(dp),parameter :: palpha(ng,nf) = reshape((/ &
!* 1s
   2.227660584e+0_dp, 4.057711562e-1_dp, 1.098175104e-1_dp, &
!* 2s
   2.581578398e+0_dp, 1.567622104e-1_dp, 6.018332272e-2_dp, &
!* 3s
   5.641487709e-1_dp, 6.924421391e-2_dp, 3.269529097e-2_dp, &
!* 4s
   2.267938753e-1_dp, 4.448178019e-2_dp, 2.195294664e-2_dp, &
!* 5s
   1.080198458e-1_dp, 4.408119382e-2_dp, 2.610811810e-2_dp, &
!* 2p
   9.192379002e-1_dp, 2.359194503e-1_dp, 8.009805746e-2_dp, &
!* 3p
   2.692880368e+0_dp, 1.489359592e-1_dp, 5.739585040e-2_dp, &
!* 4p
   4.859692220e-1_dp, 7.430216918e-2_dp, 3.653340923e-2_dp, &
!* 5p
   2.127482317e-1_dp, 4.729648620e-2_dp, 2.604865324e-2_dp, &
!* 3d
   5.229112225e-1_dp, 1.639595876e-1_dp, 6.386630021e-2_dp, &
!* 4d
   1.777717219e-1_dp, 8.040647350e-2_dp, 3.949855551e-2_dp, &
!* 5d
   4.913352950e-1_dp, 7.329090601e-2_dp, 3.594209290e-2_dp, &
!* 4f
   3.483826963e-1_dp, 1.249380537e-1_dp, 5.349995725e-2_dp, &
!* 5f
   1.649233885e-1_dp, 7.487066646e-2_dp, 3.135787219e-2_dp, &
!* 5g
   2.545432122e-1_dp, 1.006544376e-1_dp, 4.624463922e-2_dp /),shape(palpha))
   real(dp),parameter :: pcoeff(ng,nf) = reshape((/ &
!* 1s
   1.543289673e-1_dp, 5.353281423e-1_dp, 4.446345422e-1_dp, &
!* 2s
  -5.994474934e-2_dp, 5.960385398e-1_dp, 4.581786291e-1_dp, &
!* 3s
  -1.782577972e-1_dp, 8.612761663e-1_dp, 2.261841969e-1_dp, &
!* 4s
  -3.349048323e-1_dp, 1.056744667e+0_dp, 1.256661680e-1_dp, &
!* 5s
  -6.617401158e-1_dp, 7.467595004e-1_dp, 7.146490945e-1_dp, &
!* 2p
   1.623948553e-1_dp, 5.661708862e-1_dp, 4.223071752e-1_dp, &
!* 3p
  -1.061945788e-2_dp, 5.218564264e-1_dp, 5.450015143e-1_dp, &
!* 4p
  -6.147823411e-2_dp, 6.604172234e-1_dp, 3.932639495e-1_dp, &
!* 5p
  -1.389529695e-1_dp, 8.076691064e-1_dp, 2.726029342e-1_dp, &
!* 3d
   1.686596060e-1_dp, 5.847984817e-1_dp, 4.056779523e-1_dp, &
!* 4d
   2.308552718e-1_dp, 6.042409177e-1_dp, 2.595768926e-1_dp, &
!* 5d
  -2.010175008e-2_dp, 5.899310608e-1_dp, 4.658445960e-1_dp, &
!* 4f
   1.737856685e-1_dp, 5.973380628e-1_dp, 3.929395614e-1_dp, &
!* 5f
   1.909729355e-1_dp, 6.146060459e-1_dp, 3.059611271e-1_dp, &
!* 5g
   1.780980905e-1_dp, 6.063757846e-1_dp, 3.828552923e-1_dp /),shape(pcoeff))
   
!* <φ|φ> = (2i-1)!!(2j-1)!!(2k-1)!!/(4α)^(i+j+k) · sqrt(π/2α)³
!  N² = (4α)^(i+j+k)/((2i-1)!!(2j-1)!!(2k-1)!!)  · sqrt(2α/π)³
!  N = (4α)^((i+j+k)/2) / sqrt((2i-1)!!(2j-1)!!(2k-1)!!) · (2α/π)^(3/4)
   alpha(:ng) = palpha(:,ityp) * zeta**2
   coeff(:ng) = pcoeff(:,ityp) * (top*alpha(:ng))**0.75_dp &
   &  * sqrt(4*alpha(:ng))**(l(ityp)) / sqrt(dfactorial(l(ityp)+1))

end subroutine sto3g

pure subroutine sto4g(ityp,zeta,alpha,coeff)
   implicit none
   integer, parameter   :: dp=kind(1.0d0)
   integer, intent(in)  :: ityp
   real(dp),intent(in)  :: zeta
   real(dp),intent(out) :: alpha(*)
   real(dp),intent(out) :: coeff(*)

   integer,parameter :: ng=4,nf=15,l(*)=(/0,0,0,0,0,1,1,1,1,2,2,2,3,3,4/)
   real(dp),parameter :: top = 0.5_dp/atan(1.0_dp) ! 0.636619772367582_dp
   real(dp),parameter :: dfactorial(*) = & ! see OEIS A001147
   &  (/1._dp,1._dp,3._dp,15._dp,105._dp,945._dp,10395._dp,135135._dp/)
   real(dp),parameter :: palpha(ng,nf) = reshape((/ &
!* 1s
   5.216844534e+0_dp, 9.546182760e-1_dp, &
   2.652034102e-1_dp, 8.801862774e-2_dp, &
!* 2s
   1.161525551e+1_dp, 2.000243111e+0_dp, &
   1.607280687e-1_dp, 6.125744532e-2_dp, &
!* 3s
   1.513265591e+0_dp, 4.262497508e-1_dp, &
   7.643320863e-2_dp, 3.760545063e-2_dp, &
!* 4s
   3.242212833e-1_dp, 1.663217177e-1_dp, &
   5.081097451e-2_dp, 2.829066600e-2_dp, &
!* 5s
   8.602284252e-1_dp, 1.189050200e-1_dp, &
   3.446076176e-2_dp, 1.974798796e-2_dp, &
!* 2p
   1.798260992e+0_dp, 4.662622228e-1_dp, &
   1.643718620e-1_dp, 6.543927065e-2_dp, &
!* 3p
   1.853180239e+0_dp, 1.915075719e-1_dp, &
   8.655487938e-2_dp, 4.184253862e-2_dp, &
!* 4p
   1.492607880e+0_dp, 4.327619272e-1_dp, &
   7.553156064e-2_dp, 3.706272183e-2_dp, &
!* 5p
   3.962838833e-1_dp, 1.838858552e-1_dp, &
   4.943555157e-2_dp, 2.750222273e-2_dp, &
!* 3d
   9.185846715e-1_dp, 2.920461109e-1_dp, &
   1.187568890e-1_dp, 5.286755896e-2_dp, &
!* 4d
   1.995825422e+0_dp, 1.823461280e-1_dp, &
   8.197240896e-2_dp, 4.000634951e-2_dp, &
!* 5d
   4.230617826e-1_dp, 8.293863702e-2_dp, &
   4.590326388e-2_dp, 2.628744797e-2_dp, &
!* 4f
   5.691670217e-1_dp, 2.074585819e-1_dp, &
   9.298346885e-2_dp, 4.473508853e-2_dp, &
!* 5f
   2.017831152e-1_dp, 1.001952178e-1_dp, &
   5.441006630e-2_dp, 3.037569283e-2_dp, &
!* 5g
   3.945205573e-1_dp, 1.588100623e-1_dp, &
   7.646521729e-1_dp, 3.898703611e-2_dp /),shape(palpha))
   real(dp),parameter :: pcoeff(ng,nf) = reshape((/ &
!* 1s
   5.675242080e-2_dp, 2.601413550e-1_dp, &
   5.328461143e-1_dp, 2.916254405e-1_dp, &
!* 2s
  -1.198411747e-2_dp,-5.472052539e-2_dp, &
   5.805587176e-1_dp, 4.770079976e-1_dp, &
!* 3s
  -3.295496352e-2_dp,-1.724516959e-1_dp, &
   7.518511194e-1_dp, 3.589627317e-1_dp, &
!* 4s
  -1.120682822e-1_dp,-2.845426863e-1_dp, &
   8.909873788e-1_dp, 3.517811205e-1_dp, &
!* 5s
   1.103657561e-2_dp,-5.606519023e-1_dp, &
   1.179429987e+0_dp, 1.734974376e-1_dp, &
!* 2p
   5.713170255e-2_dp, 2.857455515e-1_dp, &
   5.517873105e-1_dp, 2.632314924e-1_dp, &
!* 3p
  -1.434249391e-2_dp, 2.755177589e-1_dp, &
   5.846750879e-1_dp, 2.144986514e-1_dp, &
!* 4p
  -6.035216774e-3_dp,-6.013310874e-2_dp, &
   6.451518200e-1_dp, 4.117923820e-1_dp, &
!* 5p
  -1.801459207e-2_dp,-1.360777372e-1_dp, &
   7.533973719e-1_dp, 3.409304859e-1_dp, &
!* 3d
   5.799057705e-2_dp, 3.045581349e-1_dp, &
   5.601358038e-1_dp, 2.432423313e-1_dp, &
!* 4d
  -2.816702620e-3_dp, 2.177095871e-1_dp, &
   6.058047348e-1_dp, 2.717811257e-1_dp, &
!* 5d
  -2.421626009e-2_dp, 3.937644956e-1_dp, &
   5.489520286e-1_dp, 1.190436963e-1_dp, &
!* 4f
   5.902730589e-2_dp, 3.191828952e-1_dp, &
   5.639423893e-1_dp, 2.284796537e-1_dp, &
!* 5f
   9.174268830e-2_dp, 4.023496947e-1_dp, &
   4.937432100e-1_dp, 1.254001522e-1_dp, &
!* 5g
   6.010484250e-2_dp, 3.309738329e-1_dp, &
   5.655207585e-1_dp, 2.171122608e-1_dp /),shape(pcoeff))
   
!* <φ|φ> = (2i-1)!!(2j-1)!!(2k-1)!!/(4α)^(i+j+k) · sqrt(π/2α)³
!  N² = (4α)^(i+j+k)/((2i-1)!!(2j-1)!!(2k-1)!!)  · sqrt(2α/π)³
!  N = (4α)^((i+j+k)/2) / sqrt((2i-1)!!(2j-1)!!(2k-1)!!) · (2α/π)^(3/4)
   alpha(:ng) = palpha(:,ityp) * zeta**2
   coeff(:ng) = pcoeff(:,ityp) * (top*alpha(:ng))**0.75_dp &
   &  * sqrt(4*alpha(:ng))**(l(ityp)) / sqrt(dfactorial(l(ityp)+1))

end subroutine sto4g

pure subroutine sto5g(ityp,zeta,alpha,coeff)
   implicit none
   integer, parameter   :: dp=kind(1.0d0)
   integer, intent(in)  :: ityp
   real(dp),intent(in)  :: zeta
   real(dp),intent(out) :: alpha(*)
   real(dp),intent(out) :: coeff(*)

   integer,parameter :: ng=5,nf=15,l(*)=(/0,0,0,0,0,1,1,1,1,2,2,2,3,3,4/)
   real(dp),parameter :: top = 0.5_dp/atan(1.0_dp) !0.636619772367582_dp
   real(dp),parameter :: dfactorial(*) = & ! see OEIS A001147
   &  (/1._dp,1._dp,3._dp,15._dp,105._dp,945._dp,10395._dp,135135._dp/)
   real(dp),parameter :: palpha(ng,nf) = reshape((/ &
!* 1s
   1.130563696e+1_dp, 2.071728178e+0_dp, 5.786484833e-1_dp, &
   1.975724573e-1_dp, 7.445271746e-2_dp, &
!* 2s
   8.984956862e+0_dp, 1.673710636e+0_dp, 1.944726668e-1_dp, &
   8.806345634e-2_dp, 4.249068522e-2_dp, &
!* 3s
   4.275877914e+0_dp, 1.132409433e+0_dp, 4.016256968e-1_dp, &
   7.732370620e-2_dp, 3.800708627e-2_dp, &
!* 4s
   2.980263783e+0_dp, 3.792228833e-1_dp, 1.789717224e-1_dp, &
   5.002110360e-2_dp, 2.789361681e-2_dp, &
!* 5s
   7.403763257e-1_dp, 1.367990863e-1_dp, 9.135301779e-2_dp, &
   3.726907315e-2_dp, 2.241490836e-2_dp, &
!* 2p
   3.320386533e+0_dp, 8.643257633e-1_dp, 3.079819284e-1_dp, &
   1.273309895e-1_dp, 5.606243164e-2_dp, &
!* 3p
   6.466803859e+0_dp, 1.555914802e+0_dp, 1.955925255e-1_dp, &
   8.809647701e-2_dp, 4.234835707e-2_dp, &
!* 4p
   1.091977298e+0_dp, 3.719985051e-1_dp, 8.590019352e-2_dp, &
   4.786503860e-2_dp, 2.730479990e-2_dp, &
!* 5p
   3.422168934e-1_dp, 1.665099900e-1_dp, 5.443732013e-2_dp, &
   3.367775277e-2_dp, 2.091949042e-2_dp, &
!* 3d
   1.539033958e+0_dp, 4.922090297e-1_dp, 2.029756928e-1_dp, &
   9.424112917e-2_dp, 4.569058269e-2_dp, &
!* 4d
   1.522122079e+0_dp, 2.173041823e-1_dp, 1.084876577e-1_dp, &
   5.836797641e-2_dp, 3.206682246e-2_dp, &
!* 5d
   9.702946470e-1_dp, 3.603270196e-1_dp, 8.668717752e-2_dp, &
   4.833708379e-2_dp, 2.751899341e-2_dp, &
!* 4f
   8.925960415e-1_dp, 3.277589120e-1_dp, 1.492869962e-1_dp, &
   7.506099109e-2_dp, 3.892475795e-2_dp, &
!* 5f
   1.670735676e+0_dp, 2.072477219e-1_dp, 1.024709357e-1_dp, &
   5.531913898e-2_dp, 3.072866652e-2_dp, &
!* 5g
   5.895429375e-1_dp, 2.393343780e-1_dp, 1.172646904e-1_dp, &
   6.254074479e-2_dp, 3.411243214e-2_dp /),shape(palpha))
   real(dp),parameter :: pcoeff(ng,nf) = reshape((/ &
!* 1s
   2.214055312e-2_dp, 1.135411520e-1_dp, 3.318161484e-1_dp, &
   4.825700713e-1_dp, 1.935721966e-1_dp, &
!* 2s
  -1.596349096e-2_dp,-5.685884883e-2_dp, 3.698265599e-1_dp, &
   5.480512593e-1_dp, 1.472634893e-1_dp, &
!* 3s
  -3.920358850e-3_dp,-4.168430506e-2_dp,-1.637440990e-1_dp, &
   7.419373723e-1_dp, 3.724364929e-1_dp, &
!* 4s
   1.513948997e-3_dp,-7.316801518e-2_dp,-3.143703799e-1_dp, &
   9.032615169e-1_dp, 3.294210848e-1_dp, &
!* 5s
   1.375523371e-2_dp,-3.097344179e-1_dp,-3.199192259e-1_dp, &
   1.084547038e+0_dp, 3.345288361e-1_dp, &
!* 2p
   2.079051117e-2_dp, 1.235472099e-1_dp, 3.667738986e-1_dp, &
   4.834930290e-1_dp, 1.653444074e-1_dp, &
!* 3p
  -2.329023747e-3_dp,-1.357395221e-2_dp, 2.632185383e-1_dp, &
   5.880427024e-1_dp, 2.242794445e-1_dp, &
!* 4p
  -1.143929558e-2_dp,-6.322651538e-2_dp, 4.398907721e-1_dp, &
   5.245859166e-1_dp, 1.017072253e-1_dp, &
!* 5p
  -3.113958289e-2_dp,-1.374007017e-1_dp, 5.573881018e-1_dp, &
   4.855428100e-1_dp, 6.605423564e-2_dp, &
!* 3d
   2.020869128e-2_dp, 1.321157923e-1_dp, 3.911240346e-1_dp, &
   4.779609701e-1_dp, 1.463662294e-1_dp, &
!* 4d
  -3.673711876e-3_dp, 1.167122499e-1_dp, 4.216476416e-1_dp, &
   4.547673415e-1_dp, 1.037803318e-1_dp, &
!* 5d
  -3.231527611e-3_dp,-2.434931372e-2_dp, 3.440817054e-1_dp, &
   5.693674316e-1_dp, 1.511340183e-1_dp, &
!* 4f
   1.999839052e-2_dp, 1.395427440e-1_dp, 4.091508237e-1_dp, &
   4.708252119e-1_dp, 1.328082566e-1_dp, &
!* 5f
  -7.301193568e-4_dp, 8.414991343e-2_dp, 3.923683153e-1_dp, &
   5.040033146e-1_dp, 1.328979300e-1_dp, &
!* 5g
   1.998085812e-2_dp, 1.460384050e-1_dp, 4.230565459e-1_dp, &
   4.635699665e-1_dp, 1.226411691e-1_dp /),shape(pcoeff))
   
!* <φ|φ> = (2i-1)!!(2j-1)!!(2k-1)!!/(4α)^(i+j+k) · sqrt(π/2α)³
!  N² = (4α)^(i+j+k)/((2i-1)!!(2j-1)!!(2k-1)!!)  · sqrt(2α/π)³
!  N = (4α)^((i+j+k)/2) / sqrt((2i-1)!!(2j-1)!!(2k-1)!!) · (2α/π)^(3/4)
   alpha(:ng) = palpha(:,ityp) * zeta**2
   coeff(:ng) = pcoeff(:,ityp) * (top*alpha(:ng))**0.75_dp &
   &  * sqrt(4*alpha(:ng))**(l(ityp)) / sqrt(dfactorial(l(ityp)+1))

end subroutine sto5g

pure subroutine sto6g(ityp,zeta,alpha,coeff)
   implicit none
   integer, parameter   :: dp=kind(1.0d0)
   integer, intent(in)  :: ityp
   real(dp),intent(in)  :: zeta
   real(dp),intent(out) :: alpha(*)
   real(dp),intent(out) :: coeff(*)

   integer,parameter :: ng=6,nf=15,l(*)=(/0,0,0,0,0,1,1,1,1,2,2,2,3,3,4/)
   real(dp),parameter :: top = 0.5_dp/atan(1.0_dp) !0.636619772367582_dp
   real(dp),parameter :: dfactorial(*) = & ! see OEIS A001147
   &  (/1._dp,1._dp,3._dp,15._dp,105._dp,945._dp,10395._dp,135135._dp/)
   real(dp),parameter :: palpha(ng,nf) = reshape((/ &
!* 1s
   2.310303149e+1_dp, 4.235915534e+0_dp, 1.185056519e+0_dp, &
   4.070988982e-1_dp, 1.580884151e-1_dp, 6.510953954e-2_dp, &
!* 2s
   2.768496241e+1_dp, 5.077140627e+0_dp, 1.426786050e+0_dp, &
   2.040335729e-1_dp, 9.260298399e-2_dp, 4.416183978e-2_dp, &
!* 3s
   3.273031938e+0_dp, 9.200611311e-1_dp, 3.593349765e-1_dp, &
   6.636686991e-2_dp, 4.797373812e-2_dp, 2.724741144e-2_dp, &
!* 4s
   3.232838646e+0_dp, 3.605788802e-1_dp, 1.717905487e-1_dp, &
   5.277666487e-2_dp, 3.163400284e-2_dp, 1.874093091e-2_dp, &
!* 5s
   1.410128298e+0_dp, 5.077878915e-1_dp, 1.847926858e-1_dp, &
   1.061070594e-1_dp, 3.669584901e-2_dp, 2.213558430e-2_dp, &
!* 2p
   5.868285913e+0_dp, 1.530329631e+0_dp, 5.475665231e-1_dp, &
   2.288932733e-1_dp, 1.046655969e-1_dp, 4.948220127e-2_dp, &
!* 3p
   5.077973607e+0_dp, 1.340786940e+0_dp, 2.248434849e-1_dp, &
   1.131741848e-1_dp, 6.076408893e-2_dp, 3.315424265e-2_dp, &
!* 4p
   2.389722618e+0_dp, 7.960947826e-1_dp, 3.415541380e-1_dp, &
   8.847434525e-2_dp, 4.958248334e-2_dp, 2.816929784e-2_dp, &
!* 5p
   3.778623374e+0_dp, 3.499121109e-1_dp, 1.683175469e-1_dp, &
   5.404070736e-2_dp, 3.328911801e-2_dp, 2.063815019e-2_dp, &
!* 3d
   2.488296923e+0_dp, 7.981487853e-1_dp, 3.311327490e-1_dp, &
   1.559114463e-1_dp, 7.817734732e-2_dp, 4.058484363e-2_dp, &
!* 4d
   4.634239420e+0_dp, 1.341648295e+0_dp, 2.209593028e-1_dp, &
   1.101467943e-1_dp, 5.904190370e-2_dp, 3.232628887e-2_dp, &
!* 5d
   8.820520428e-1_dp, 3.410838409e-1_dp, 9.204308840e-2_dp, &
   5.472831774e-2_dp, 3.391202830e-2_dp, 2.108227374e-2_dp, &
!* 4f
   1.357718039e+0_dp, 5.004907278e-1_dp, 2.296565064e-1_dp, &
   1.173146814e-1_dp, 6.350097171e-2_dp, 3.474556673e-2_dp, &
!* 5f
   1.334096840e+0_dp, 2.372312347e-1_dp, 1.269485144e-1_dp, &
   7.290318381e-2_dp, 4.351355997e-2_dp, 2.598071843e-2_dp, &
!* 5g
   8.574668996e-1_dp, 3.497184772e-1_dp, 1.727917060e-1_dp, &
   9.373643151e-2_dp, 5.340032759e-2_dp, 3.051364464e-2_dp /),shape(palpha))
   real(dp),parameter :: pcoeff(ng,nf) = reshape((/ &
!* 1s
   9.163596280e-3_dp, 4.936149294e-2_dp, 1.685383049e-1_dp, &
   3.705627997e-1_dp, 4.164915298e-1_dp, 1.303340841e-1_dp, &
!* 2s
  -4.151277819e-3_dp,-2.067024148e-2_dp,-5.150303337e-2_dp, &
   3.346271174e-1_dp, 5.621061301e-1_dp, 1.712994697e-1_dp, &
!* 3s
  -6.775596947e-3_dp,-5.639325779e-2_dp,-1.587656086e-1_dp, &
   5.534527651e-1_dp, 5.015351020e-1_dp, 7.223633674e-2_dp, &
!* 4s
   1.374817488e-3_dp,-8.666390043e-2_dp,-3.130627309e-1_dp, &
   7.812787397e-1_dp, 4.389247988e-1_dp, 2.487178756e-2_dp, &
!* 5s
   2.695439582e-3_dp, 1.850157487e-2_dp,-9.588628125e-2_dp, &
  -5.200673560e-1_dp, 1.087619490e+0_dp, 3.103964343e-1_dp, &
!* 2p
   7.924233646e-3_dp, 5.144104825e-2_dp, 1.898400060e-1_dp, &
   4.049863191e-1_dp, 4.012362861e-1_dp, 1.051855189e-1_dp, &
!* 3p
  -3.329929840e-3_dp,-1.419488340e-2_dp, 1.639395770e-1_dp, &
   4.485358256e-1_dp, 3.908813050e-1_dp, 7.411456232e-2_dp, &
!* 4p
  -1.665913575e-3_dp,-1.657464971e-2_dp,-5.958513378e-2_dp, &
   4.053115554e-1_dp, 5.433958189e-1_dp, 1.204970491e-1_dp, &
!* 5p
   1.163246387e-4_dp,-2.920771322e-2_dp,-1.381051233e-1_dp, &
   5.706134877e-1_dp, 4.768808140e-1_dp, 6.021665516e-2_dp, &
!* 3d
   7.283828112e-3_dp, 5.386799363e-2_dp, 2.072139149e-1_dp, &
   4.266269092e-1_dp, 3.843100204e-1_dp, 8.902827546e-2_dp, &
!* 4d
  -4.749842876e-4_dp,-3.566777891e-3_dp, 1.108670481e-1_dp, &
   4.159646930e-1_dp, 4.621672517e-1_dp, 1.081250196e-1_dp, &
!* 5d
  -4.097311019e-3_dp,-2.508271857e-2_dp, 2.648458555e-1_dp, &
   5.097437054e-1_dp, 2.654483467e-1_dp, 2.623132212e-2_dp, &
!* 4f
   6.930234381e-3_dp, 5.634263745e-2_dp, 2.217065797e-1_dp, &
   4.411388883e-1_dp, 3.688112625e-1_dp, 7.787514504e-2_dp, &
!* 5f
  -9.486751531e-4_dp, 4.624275998e-2_dp, 2.373699784e-1_dp, &
   4.589112231e-1_dp, 3.205010548e-1_dp, 5.077063693e-2_dp, &
!* 5g
   6.729778096e-3_dp, 5.874145170e-2_dp, 2.339955227e-1_dp, &
   4.512983737e-1_dp, 3.552053926e-1_dp, 6.974153145e-2_dp /),shape(pcoeff))
   
!* <φ|φ> = (2i-1)!!(2j-1)!!(2k-1)!!/(4α)^(i+j+k) · sqrt(π/2α)³
!  N² = (4α)^(i+j+k)/((2i-1)!!(2j-1)!!(2k-1)!!)  · sqrt(2α/π)³
!  N = (4α)^((i+j+k)/2) / sqrt((2i-1)!!(2j-1)!!(2k-1)!!) · (2α/π)^(3/4)
   alpha(:ng) = palpha(:,ityp) * zeta**2
   coeff(:ng) = pcoeff(:,ityp) * (top*alpha(:ng))**0.75_dp &
   &  * sqrt(4*alpha(:ng))**(l(ityp)) / sqrt(dfactorial(l(ityp)+1))

end subroutine sto6g

