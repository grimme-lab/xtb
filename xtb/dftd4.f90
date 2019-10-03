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

   module dftd4
   use iso_fortran_env, only : wp => real64
   use mctc_constants, only : pi
   !! ========================================================================
   !  mix in the covalent coordination number from the ncoord module
   !  also get the CN-Parameters to inline the CN-derivative in the gradient
   use ncoord, only : covncoord => ncoord_d4, kn,k1,k4,k5,k6
   use tbdef_param, only : dftd_parameter
   implicit none

      real(wp) :: thopi,ootpi
      parameter ( thopi = 3._wp/pi )
      parameter ( ootpi = 0.5_wp/pi )

      integer,parameter :: p_refq_gfn2xtb          = 0
      integer,parameter :: p_refq_gasteiger        = 1
      integer,parameter :: p_refq_hirshfeld        = 2
      integer,parameter :: p_refq_periodic         = 3
      integer,parameter :: p_refq_gfn2xtb_gbsa_h2o = 4
      integer,parameter :: p_refq_goedecker        = 5

      integer,parameter :: p_mbd_none       = 0
      integer,parameter :: p_mbd_rpalike    = 1
      integer,parameter :: p_mbd_exact_atm  = 2
      integer,parameter :: p_mbd_approx_atm = 3

      integer,private,parameter :: max_elem = 118
      real(wp),parameter :: zeff(max_elem) = (/ &
      &   1,                                                 2,  & ! H-He
      &   3, 4,                               5, 6, 7, 8, 9,10,  & ! Li-Ne
      &  11,12,                              13,14,15,16,17,18,  & ! Na-Ar
      &  19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,  & ! K-Kr
      &   9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,  & ! Rb-Xe
      &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Cs-Lu
      &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26, & ! Hf-Rn
      !  just copy & paste from above
      &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Fr-Lr
      &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 /) ! Rf-Og

   !  Semiempirical Evaluation of the GlobalHardness of the Atoms of 103
   !  Elements of the Periodic Table Using the Most Probable Radii as
   !  their Size Descriptors DULAL C. GHOSH, NAZMUL ISLAM 2009 in 
   !  Wiley InterScience (www.interscience.wiley.com).
   !  DOI 10.1002/qua.22202
   !  values in the paper multiplied by two because 
   !  (ii:ii)=(IP-EA)=d^2 E/dN^2 but the hardness
   !  definition they use is 1/2d^2 E/dN^2 (in Eh)
      real(wp),parameter :: gam(1:max_elem) = (/ &
     &0.47259288,0.92203391,0.17452888,0.25700733,0.33949086,0.42195412, & ! H-C
     &0.50438193,0.58691863,0.66931351,0.75191607,0.17964105,0.22157276, & ! N-Mg
     &0.26348578,0.30539645,0.34734014,0.38924725,0.43115670,0.47308269, & ! Al-Ar
     &0.17105469,0.20276244,0.21007322,0.21739647,0.22471039,0.23201501, & ! Ca-Cr
     &0.23933969,0.24665638,0.25398255,0.26128863,0.26859476,0.27592565, & ! Mn-Zn
     &0.30762999,0.33931580,0.37235985,0.40273549,0.43445776,0.46611708, & ! Ga-Kr
     &0.15585079,0.18649324,0.19356210,0.20063311,0.20770522,0.21477254, & ! Rb-Mo
     &0.22184614,0.22891872,0.23598621,0.24305612,0.25013018,0.25719937, & ! Tc-Cd
     &0.28784780,0.31848673,0.34912431,0.37976593,0.41040808,0.44105777, & ! In-Xe
     &0.05019332,0.06762570,0.08504445,0.10247736,0.11991105,0.13732772, & ! Cs-Nd
     &0.15476297,0.17218265,0.18961288,0.20704760,0.22446752,0.24189645, & ! Pm-Dy
     &0.25932503,0.27676094,0.29418231,0.31159587,0.32902274,0.34592298, & ! Ho-Hf
     &0.36388048,0.38130586,0.39877476,0.41614298,0.43364510,0.45104014, & ! Ta-Pt
     &0.46848986,0.48584550,0.12526730,0.14268677,0.16011615,0.17755889, & ! Au-Po
     &0.19497557,0.21240778,0.07263525,0.09422158,0.09920295,0.10418621, & ! At-Th
     &0.14235633,0.16394294,0.18551941,0.22370139,0.00000000,0.00000000, & ! Pa-Cm
     &0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000, & ! Bk-No
     &0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000, & ! Rf-Mt
     &0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000, & ! Ds-Mc
     &0.00000000,0.00000000,0.00000000,0.00000000 /) ! Lv-Og

   !  pauling EN's 
      real(wp),parameter :: en(max_elem) = (/ &
      & 2.20,3.00, & ! H,He
      & 0.98,1.57,2.04,2.55,3.04,3.44,3.98,4.50, & ! Li-Ne
      & 0.93,1.31,1.61,1.90,2.19,2.58,3.16,3.50, & ! Na-Ar
      & 0.82,1.00, & ! K,Ca
      &           1.36,1.54,1.63,1.66,1.55,1.83,1.88,1.91,1.90,1.65, & ! Sc-Zn
      &           1.81,2.01,2.18,2.55,2.96,3.00, & ! Ga-Kr
      & 0.82,0.95, & ! Rb,Sr
      &           1.22,1.33,1.60,2.16,1.90,2.20,2.28,2.20,1.93,1.69, & ! Y-Cd
      &           1.78,1.96,2.05,2.10,2.66,2.60, & ! In-Xe
      & 0.79,0.89, & ! Cs,Ba
      &      1.10,1.12,1.13,1.14,1.15,1.17,1.18, & ! La-Eu
      &      1.20,1.21,1.22,1.23,1.24,1.25,1.26, & ! Gd-Yb
      &           1.27,1.30,1.50,2.36,1.90,2.20,2.20,2.28,2.54,2.00, & ! Lu-Hg
      &           1.62,2.33,2.02,2.00,2.20,2.20, & ! Tl-Rn
      ! only dummies below
      & 1.50,1.50, & ! Fr,Ra
      &      1.50,1.50,1.50,1.50,1.50,1.50,1.50, & ! Ac-Am
      &      1.50,1.50,1.50,1.50,1.50,1.50,1.50, & ! Cm-No
      &           1.50,1.50,1.50,1.50,1.50,1.50,1.50,1.50,1.50,1.50, & ! Rf-Cn
      &           1.50,1.50,1.50,1.50,1.50,1.50 /) ! Nh-Og
   ! same values in old arrangement
   ! &         2.200,3.000,0.980,1.570,2.040,2.550,3.040,3.440,3.980 &
   ! &        ,4.500,0.930,1.310,1.610,1.900,2.190,2.580,3.160,3.500 &
   ! &        ,0.820,1.000,1.360,1.540,1.630,1.660,1.550,1.830,1.880 &
   ! &        ,1.910,1.900,1.650,1.810,2.010,2.180,2.550,2.960,3.000 &
   ! &        ,0.820,0.950,1.220,1.330,1.600,2.160,1.900,2.200,2.280 &
   ! &        ,2.200,1.930,1.690,1.780,1.960,2.050,2.100,2.660,2.600 &
   ! &,0.79,0.89,1.10,1.12,1.13,1.14,1.15,1.17,1.18,1.20,1.21,1.22 &
   ! &,1.23,1.24,1.25,1.26,1.27,1.3,1.5,2.36,1.9,2.2,2.20,2.28,2.54 &
   ! &,2.00,1.62,2.33,2.02,2.0,2.2,2.2,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5/)

   !  D3 radii
   !   real(wp),parameter :: rcov(max_elem) = (/ &
   !  & 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865, &
   !  & 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527, &
   !  & 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820, &
   !  & 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730, &
   !  & 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923, &
   !  & 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188, &
   !  & 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349, &
   !  & 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216, &
   !  & 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717, &
   !  & 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967, &
   !  & 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625, &
   !  & 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657, &
   !  & 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833, &
   !  & 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098, &
   !  & 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878, &
   !  & 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790, &
   !  & 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584, &
   !  & 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289, &
   !  & 3.82984466, 3.85504098, 3.88023730, 3.90543362 /)

   !  covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
   !  188-197), values for metals decreased by 10 %
      real(wp),private,parameter :: rad(max_elem) = (/  &
      & 0.32,0.46, & ! H,He
      & 1.20,0.94,0.77,0.75,0.71,0.63,0.64,0.67, & ! Li-Ne
      & 1.40,1.25,1.13,1.04,1.10,1.02,0.99,0.96, & ! Na-Ar
      & 1.76,1.54, & ! K,Ca
      &           1.33,1.22,1.21,1.10,1.07,1.04,1.00,0.99,1.01,1.09, & ! Sc-Zn
      &           1.12,1.09,1.15,1.10,1.14,1.17, & ! Ga-Kr
      & 1.89,1.67, & ! Rb,Sr
      &           1.47,1.39,1.32,1.24,1.15,1.13,1.13,1.08,1.15,1.23, & ! Y-Cd
      &           1.28,1.26,1.26,1.23,1.32,1.31, & ! In-Xe
      & 2.09,1.76, & ! Cs,Ba
      &      1.62,1.47,1.58,1.57,1.56,1.55,1.51, & ! La-Eu
      &      1.52,1.51,1.50,1.49,1.49,1.48,1.53, & ! Gd-Yb
      &           1.46,1.37,1.31,1.23,1.18,1.16,1.11,1.12,1.13,1.32, & ! Lu-Hg
      &           1.30,1.30,1.36,1.31,1.38,1.42, & ! Tl-Rn
      & 2.01,1.81, & ! Fr,Ra
      &      1.67,1.58,1.52,1.53,1.54,1.55,1.49, & ! Ac-Am
      &      1.49,1.51,1.51,1.48,1.50,1.56,1.58, & ! Cm-No
      &           1.45,1.41,1.34,1.29,1.27,1.21,1.16,1.15,1.09,1.22, & ! Lr-Cn
      &           1.36,1.43,1.46,1.58,1.48,1.57 /) ! Nh-Og
      real(wp),parameter :: rcov(max_elem) = 4.0_wp/3.0_wp * rad / 0.52917726_wp


   !  r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
   !  the large number of digits is just to keep the results consistent
   !  with older versions. They should not imply any higher accuracy
   !  than the old values
   !   real(wp),parameter :: r4r2(max_elem) = (/ &
   !   &   2.00734898,  1.56637132,  5.01986934,  3.85379032,  3.64446594, &
   !   &   3.10492822,  2.71175247,  2.59361680,  2.38825250,  2.21522516, &
   !   &   6.58585536,  5.46295967,  5.65216669,  4.88284902,  4.29727576, &
   !   &   4.04108902,  3.72932356,  3.44677275,  7.97762753,  7.07623947, &
   !   &   6.60844053,  6.28791364,  6.07728703,  5.54643096,  5.80491167, &
   !   &   5.58415602,  5.41374528,  5.28497229,  5.22592821,  5.09817141, &
   !   &   6.12149689,  5.54083734,  5.06696878,  4.87005108,  4.59089647, &
   !   &   4.31176304,  9.55461698,  8.67396077,  7.97210197,  7.43439917, &
   !   &   6.58711862,  6.19536215,  6.01517290,  5.81623410,  5.65710424, &
   !   &   5.52640661,  5.44263305,  5.58285373,  7.02081898,  6.46815523, &
   !   &   5.98089120,  5.81686657,  5.53321815,  5.25477007, 11.02204549, &
   !   &  10.15679528,  9.35167836,  9.06926079,  8.97241155,  8.90092807, &
   !   &   8.85984840,  8.81736827,  8.79317710,  7.89969626,  8.80588454, &
   !   &   8.42439218,  8.54289262,  8.47583370,  8.45090888,  8.47339339, &
   !   &   7.83525634,  8.20702843,  7.70559063,  7.32755997,  7.03887381, &
   !   &   6.68978720,  6.05450052,  5.88752022,  5.70661499,  5.78450695, &
   !   &   7.79780729,  7.26443867,  6.78151984,  6.67883169,  6.39024318, &
   !   &   6.09527958, 11.79156076, 11.10997644,  9.51377795,  8.67197068, &
   !   &   8.77140725,  8.65402716,  8.53923501,  8.85024712 /)

   !  PBE0/def2-QZVP atomic values calculated by S. Grimme in Gaussian (2010)
   !  rare gases recalculated by J. Mewes with PBE0/aug-cc-pVQZ in Dirac (2018)
   !  He: 3.4698 -> 3.5544, Ne: 3.1036 -> 3.7943, Ar: 5.6004 -> 5.6638, 
   !  Kr: 6.1971 -> 6.2312, Xe: 7.5152 -> 8.8367
   !  not replaced but recalculated (PBE0/cc-pVQZ) were
   !   H: 8.0589 ->10.9359, Li:29.0974 ->39.7226, Be:14.8517 ->17.7460
   !  also new super heavies Cn,Nh,Fl,Lv,Og
      real(wp),private,parameter :: r2r4(max_elem) = (/  &
      &  8.0589, 3.4698, & ! H,He
      & 29.0974,14.8517,11.8799, 7.8715, 5.5588, 4.7566, 3.8025, 3.1036, & ! Li-Ne
      & 26.1552,17.2304,17.7210,12.7442, 9.5361, 8.1652, 6.7463, 5.6004, & ! Na-Ar
      & 29.2012,22.3934, & ! K,Ca
      &         19.0598,16.8590,15.4023,12.5589,13.4788, & ! Sc-
      &         12.2309,11.2809,10.5569,10.1428, 9.4907, & ! -Zn
      &                 13.4606,10.8544, 8.9386, 8.1350, 7.1251, 6.1971, & ! Ga-Kr
      & 30.0162,24.4103, & ! Rb,Sr
      &         20.3537,17.4780,13.5528,11.8451,11.0355, & ! Y-
      &         10.1997, 9.5414, 9.0061, 8.6417, 8.9975, & ! -Cd
      &                 14.0834,11.8333,10.0179, 9.3844, 8.4110, 7.5152, & ! In-Xe
      & 32.7622,27.5708, & ! Cs,Ba
      &         23.1671,21.6003,20.9615,20.4562,20.1010,19.7475,19.4828, & ! La-Eu
      &         15.6013,19.2362,17.4717,17.8321,17.4237,17.1954,17.1631, & ! Gd-Yb
      &         14.5716,15.8758,13.8989,12.4834,11.4421, & ! Lu-
      &         10.2671, 8.3549, 7.8496, 7.3278, 7.4820, & ! -Hg
      &                 13.5124,11.6554,10.0959, 9.7340, 8.8584, 8.0125, & ! Tl-Rn
      & 29.8135,26.3157, & ! Fr,Ra
      &         19.1885,15.8542,16.1305,15.6161,15.1226,16.1576, 0.0000, & ! Ac-Am
      &          0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, & ! Cm-No
      &          0.0000, 0.0000, 0.0000, 0.0000, 0.0000, & ! Lr-
      &          0.0000, 0.0000, 0.0000, 0.0000, 5.4929, & ! -Cn
      &                  6.7286, 6.5144,10.9169,10.3600, 9.4723, 8.6641 /) ! Nh-Og
      integer,private :: idum
      real(wp),parameter :: r4r2(max_elem) = &
      &  sqrt(0.5_wp*(r2r4*(/(sqrt(real(idum,wp)),idum=1,max_elem)/)))

      integer, dimension(max_elem)      :: refn ! for D4
      integer, dimension(max_elem)      :: refs ! for D3 (generated on-the-fly)
      integer, dimension(7,max_elem)    :: refc
      real(wp),dimension(7,max_elem)    :: refq
      real(wp),dimension(7,max_elem)    :: refh
      real(wp),dimension(7,max_elem)    :: dftq,pbcq,gffq,solq,clsq
      real(wp),dimension(7,max_elem)    :: dfth,pbch,gffh,solh,clsh
      real(wp),dimension(7,max_elem)    :: hcount 
      real(wp),dimension(7,max_elem)    :: ascale
      real(wp),dimension(7,max_elem)    :: refcovcn
      real(wp),dimension(7,max_elem)    :: refcn
      integer, dimension(7,max_elem)    :: refsys 
      real(wp),dimension(23,7,max_elem) :: alphaiw
      real(wp),dimension(23,7,max_elem) :: refal
      real(wp),dimension(8)       :: secq
      real(wp),dimension(8)       :: dfts,pbcs,gffs,sols,clss
      real(wp),dimension(8)       :: sscale
      real(wp),dimension(8)       :: seccn
      real(wp),dimension(8)       :: seccnd3
      real(wp),dimension(23,8)    :: secaiw

      include 'param_ref.inc'

contains

subroutine prmolc6(molc6,molc8,molpol,nat,at,  &
           &       cn,covcn,q,qlmom,c6ab,alpha,rvdw,hvol)
   use iso_fortran_env, only : id => output_unit, wp => real64
   use mctc_econv, only : autoaa
   implicit none
   real(wp),intent(in)  :: molc6
   real(wp),intent(in)  :: molc8
   real(wp),intent(in)  :: molpol
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in),optional :: cn(nat)
   real(wp),intent(in),optional :: covcn(nat)
   real(wp),intent(in),optional :: q(nat)
   real(wp),intent(in),optional :: qlmom(3,nat)
   real(wp),intent(in),optional :: c6ab(nat,nat)
   real(wp),intent(in),optional :: alpha(nat)
   real(wp),intent(in),optional :: rvdw(nat)
   real(wp),intent(in),optional :: hvol(nat)
   character(len=2),   external :: asym
   integer :: i
   if(present(cn).or.present(covcn).or.present(q).or.present(c6ab) &
   &   .or.present(alpha).or.present(rvdw).or.present(hvol)) then
   write(id,'(a)')
   write(id,'(''   #   Z   '')',advance='no')
   if(present(cn))   write(id,'(''        CN'')',advance='no')
   if(present(covcn))write(id,'(''     covCN'')',advance='no')
   if(present(q))    write(id,'(''         q'')',advance='no')
   if(present(qlmom))write(id,   '(''   n(s)'')',advance='no')
   if(present(qlmom))write(id,   '(''   n(p)'')',advance='no')
   if(present(qlmom))write(id,   '(''   n(d)'')',advance='no')
   if(present(c6ab)) write(id,'(''      C6AA'')',advance='no')
   if(present(alpha))write(id,'(''      α(0)'')',advance='no')
   if(present(rvdw)) write(id,'(''    RvdW/Å'')',advance='no')
   if(present(hvol)) write(id,'(''    relVol'')',advance='no')
   write(*,'(a)')
   do i=1,nat
      write(*,'(i4,1x,i3,1x,a2)',advance='no') &
      &     i,at(i),asym(at(i))
      if(present(cn))   write(id,'(f10.3)',advance='no')cn(i)
      if(present(covcn))write(id,'(f10.3)',advance='no')covcn(i)
      if(present(q))    write(id,'(f10.3)',advance='no')q(i)
      if(present(qlmom))write(id, '(f7.3)',advance='no')qlmom(1,i)
      if(present(qlmom))write(id, '(f7.3)',advance='no')qlmom(2,i)
      if(present(qlmom))write(id, '(f7.3)',advance='no')qlmom(3,i)
      if(present(c6ab)) write(id,'(f10.3)',advance='no')c6ab(i,i)
      if(present(alpha))write(id,'(f10.3)',advance='no')alpha(i)
      if(present(rvdw)) write(id,'(f10.3)',advance='no')rvdw(i)*autoaa
      if(present(hvol)) write(id,'(f10.3)',advance='no')hvol(i)
      write(*,'(a)')
   enddo
   endif
   write(id,'(/,1x,''Mol. C6AA /au·bohr⁶  :'',f18.6,'// &
   &         '/,1x,''Mol. C8AA /au·bohr⁸  :'',f18.6,'// &
   &         '/,1x,''Mol. α(0) /au        :'',f18.6,/)') &
   &          molc6,molc8,molpol
end subroutine prmolc6

subroutine mdisp_old(nat,ndim,at,xyz, &
           &         gw,c6abns,molc6,molc8,molpol,aout,cout,rout,vout)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat) 
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   real(wp),intent(out) :: molc6
   real(wp),intent(out) :: molc8
   real(wp),intent(out) :: molpol
   real(wp),intent(out),optional :: aout(23,nat)
   real(wp),intent(out),optional :: cout(nat,nat)
   real(wp),intent(out),optional :: rout(nat)
   real(wp),intent(out),optional :: vout(nat)

   integer  :: i,ii,ia,j,jj,ja,k,l
   integer, allocatable :: itbl(:,:)
   real(wp) :: qmod,oth,iz
   real(wp),allocatable :: rvdw(:)
   real(wp),allocatable :: phv(:)
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: aw(:,:)
   parameter (oth=1._wp/3._wp)
   
   allocate( rvdw(nat),phv(nat),c6ab(nat,nat),aw(23,nat), &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   molc6  = 0._wp
   molc8  = 0._wp
   molpol = 0._wp

   k = 0
   do i = 1, nat
      ia = at(i)
      iz = zeff(ia)
      do ii = 1, refs(ia)
         k = k+1
         itbl(ii,i) = k
         aw(:,i) = aw(:,i) + gw(k) * refal(:,ii,ia)
      enddo
!     van-der-Waals radius, alpha = 4/3 pi r**3 <=> r = (3/(4pi) alpha)**(1/3)
      rvdw(i) = (0.25_wp*thopi*aw(1,i))**oth
!     pseudo-hirshfeld volume
      phv(i) = aw(1,i)/refal(1,1,ia)
      c6ab(i,i) = thopi * trapzd(aw(:,i)**2)
      molpol = molpol + aw(1,i)
      molc6  = molc6  + c6ab(i,i)
      molc8 = molc8 + 3*r4r2(ia)**2*c6ab(i,i)
      do j = 1, i-1
         ja = at(j)
         c6ab(j,i) = thopi * trapzd(aw(:,i)*aw(:,j))
         c6ab(i,j) = c6ab(j,i)
         molc6 = molc6 + 2*c6ab(j,i)
         molc8 = molc8 + 6*r4r2(ia)*r4r2(ja)*c6ab(j,i)
      enddo
   enddo

   if (present(aout)) aout = aw
   if (present(vout)) vout = phv
   if (present(rout)) rout = rvdw
   if (present(cout)) cout = c6ab

end subroutine mdisp_old

subroutine mdisp(nat,ndim,at,q,xyz,g_a,g_c, &
           &     gw,c6abns,molc6,molc8,molpol,aout,cout,rout,vout)
   use iso_fortran_env, only : wp => real64
!  use dftd4, only : thopi,gam, &
!  &  trapzd,zeta,r4r2, &
!  &  refn,refq,refal
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat) 
   real(wp),intent(in)  :: xyz(3,nat) 
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   real(wp),intent(out) :: molc6
   real(wp),intent(out) :: molc8
   real(wp),intent(out) :: molpol
   real(wp),intent(out),optional :: aout(23,nat)
   real(wp),intent(out),optional :: cout(nat,nat)
   real(wp),intent(out),optional :: rout(nat)
   real(wp),intent(out),optional :: vout(nat)

   integer  :: i,ii,ia,j,jj,ja,k,l
   integer, allocatable :: itbl(:,:)
   real(wp) :: qmod,oth,iz
   real(wp),allocatable :: zetvec(:)
   real(wp),allocatable :: rvdw(:)
   real(wp),allocatable :: phv(:)
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: aw(:,:)
   parameter (oth=1._wp/3._wp)
   
   allocate( zetvec(ndim),rvdw(nat),phv(nat),c6ab(nat,nat),aw(23,nat), &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   molc6  = 0._wp
   molc8  = 0._wp
   molpol = 0._wp

   k = 0
   do i = 1, nat
      ia = at(i)
      iz = zeff(ia)
      do ii = 1, refn(ia)
         k = k+1
         itbl(ii,i) = k
         zetvec(k) = gw(k) * zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz)
         aw(:,i) = aw(:,i) + zetvec(k) * refal(:,ii,ia)
      enddo
!     van-der-Waals radius, alpha = 4/3 pi r**3 <=> r = (3/(4pi) alpha)**(1/3)
      rvdw(i) = (0.25_wp*thopi*aw(1,i))**oth
!     pseudo-hirshfeld volume
      phv(i) = aw(1,i)/refal(1,1,ia)
      c6ab(i,i) = thopi * trapzd(aw(:,i)**2)
      molpol = molpol + aw(1,i)
      molc6  = molc6  + c6ab(i,i)
      molc8 = molc8 + 3*r4r2(ia)**2*c6ab(i,i)
      do j = 1, i-1
         ja = at(j)
         c6ab(j,i) = thopi * trapzd(aw(:,i)*aw(:,j))
         c6ab(i,j) = c6ab(j,i)
         molc6 = molc6 + 2*c6ab(j,i)
         molc8 = molc8 + 6*r4r2(ia)*r4r2(ja)*c6ab(j,i)
      enddo
   enddo

   if (present(aout)) aout = aw
   if (present(vout)) vout = phv
   if (present(rout)) rout = rvdw
   if (present(cout)) cout = c6ab

end subroutine mdisp

!pure subroutine dist_r2(nat,xyz,r2)
!   use iso_fortran_env, only : wp => real64
!   implicit none
!   integer, intent(in)  :: nat
!   real(wp),intent(in)  :: xyz(3,nat)
!   real(wp),intent(out) :: r2(nat,nat)
!
!   integer :: i,j
!
!   r2 = 0._wp
!
!   do i=1,nat
!      do j=1,i-1
!         r2(i,j) = sum( (xyz(:,i)-xyz(:,j))**2 )
!         r2(j,i) = r2(i,j)
!      enddo
!   enddo
!
!end subroutine dist_r2
!
!pure subroutine dist_bj(oor6ab,oor8ab,r2,nat,at,par)
!   use iso_fortran_env, only : wp => real64
!!  use dftd4, only : r4r2
!   implicit none
!   integer, intent(in)  :: nat,at(nat)
!   type(dftd_parameter),intent(in) :: par
!   real(wp),intent(in)  :: r2(nat,nat)
!   real(wp),intent(out) :: oor6ab(nat,nat),oor8ab(nat,nat)
!
!   integer  :: i,j
!   real(wp) :: cutoff
!
!   oor6ab = 0._wp
!   oor8ab = 0._wp
!
!   do i=1,nat
!      do j=1,i-1
!         cutoff      = par%a1*sqrt(3._wp*r4r2(at(i))*r4r2(at(j)))+par%a2
!         oor6ab(i,j) = 1._wp/(r2(i,j)**3 + cutoff**6)
!         oor6ab(j,i) = oor6ab(i,j)
!         oor8ab(i,j) = 1._wp/(r2(i,j)**4 + cutoff**8)
!         oor8ab(j,i) = oor8ab(i,j)
!      enddo
!   enddo
!
!end subroutine dist_bj

pure elemental function zeta(a,c,qref,qmod)
   use iso_fortran_env, only : wp => real64
   implicit none
   real(wp),intent(in) :: qmod,qref
   real(wp),intent(in) :: a,c
   real(wp)            :: zeta

   intrinsic :: exp

   if (qmod.lt.0._wp) then
      zeta = exp( a )
   else
      zeta = exp( a * ( 1._wp - exp( c * ( 1._wp - qref/qmod ) ) ) )
   endif

end function zeta

pure elemental function dzeta(a,c,qref,qmod)
   use iso_fortran_env, only : wp => real64
!  use dftd4, only : zeta
   implicit none
   real(wp),intent(in) :: qmod,qref
   real(wp),intent(in) :: a,c
   real(wp)            :: dzeta

   intrinsic :: exp

   if (qmod.lt.0._wp) then
      dzeta = 0._wp
   else
      dzeta = - a * c * exp( c * ( 1._wp - qref/qmod ) ) &
      &           * zeta(a,c,qref,qmod) * qref / ( qmod**2 )
   endif

end function dzeta

pure function trapzd(pol)
   use iso_fortran_env, only : wp => real64
   implicit none
   real(wp),intent(in) :: pol(23)
   real(wp)            :: trapzd

   real(wp)            :: tmp1, tmp2
   real(wp),parameter  :: freq(23) = (/ &
&   0.000001_wp,0.050000_wp,0.100000_wp, &
&   0.200000_wp,0.300000_wp,0.400000_wp, &
&   0.500000_wp,0.600000_wp,0.700000_wp, &
&   0.800000_wp,0.900000_wp,1.000000_wp, &
&   1.200000_wp,1.400000_wp,1.600000_wp, &
&   1.800000_wp,2.000000_wp,2.500000_wp, &
&   3.000000_wp,4.000000_wp,5.000000_wp, &
&   7.500000_wp,10.00000_wp /)
!  just precalculate all weights and get the job done
   real(wp),parameter :: weights(23) = 0.5_wp * (/ &
&  ( freq (2) - freq (1) ),  &
&  ( freq (2) - freq (1) ) + ( freq (3) - freq (2) ),  &
&  ( freq (3) - freq (2) ) + ( freq (4) - freq (3) ),  &
&  ( freq (4) - freq (3) ) + ( freq (5) - freq (4) ),  &
&  ( freq (5) - freq (4) ) + ( freq (6) - freq (5) ),  &
&  ( freq (6) - freq (5) ) + ( freq (7) - freq (6) ),  &
&  ( freq (7) - freq (6) ) + ( freq (8) - freq (7) ),  &
&  ( freq (8) - freq (7) ) + ( freq (9) - freq (8) ),  &
&  ( freq (9) - freq (8) ) + ( freq(10) - freq (9) ),  &
&  ( freq(10) - freq (9) ) + ( freq(11) - freq(10) ),  &
&  ( freq(11) - freq(10) ) + ( freq(12) - freq(11) ),  &
&  ( freq(12) - freq(11) ) + ( freq(13) - freq(12) ),  &
&  ( freq(13) - freq(12) ) + ( freq(14) - freq(13) ),  &
&  ( freq(14) - freq(13) ) + ( freq(15) - freq(14) ),  &
&  ( freq(15) - freq(14) ) + ( freq(16) - freq(15) ),  &
&  ( freq(16) - freq(15) ) + ( freq(17) - freq(16) ),  &
&  ( freq(17) - freq(16) ) + ( freq(18) - freq(17) ),  &
&  ( freq(18) - freq(17) ) + ( freq(19) - freq(18) ),  &
&  ( freq(19) - freq(18) ) + ( freq(20) - freq(19) ),  &
&  ( freq(20) - freq(19) ) + ( freq(21) - freq(20) ),  &
&  ( freq(21) - freq(20) ) + ( freq(22) - freq(21) ),  &
&  ( freq(22) - freq(21) ) + ( freq(23) - freq(22) ),  &
&  ( freq(23) - freq(22) ) /)

!!  do average between trap(1)-trap(22) .and. trap(2)-trap(23)
!   tmp1 = 0.5_wp * ( &
!&  ( freq (2) - freq (1) ) * ( pol (2) + pol (1) )+ &
!&  ( freq (4) - freq (3) ) * ( pol (4) + pol (3) )+ &
!&  ( freq (6) - freq (5) ) * ( pol (6) + pol (5) )+ &
!&  ( freq (8) - freq (7) ) * ( pol (8) + pol (7) )+ &
!&  ( freq(10) - freq (9) ) * ( pol(10) + pol (9) )+ &
!&  ( freq(12) - freq(11) ) * ( pol(12) + pol(11) )+ &
!&  ( freq(14) - freq(13) ) * ( pol(14) + pol(13) )+ &
!&  ( freq(16) - freq(15) ) * ( pol(16) + pol(15) )+ &
!&  ( freq(18) - freq(17) ) * ( pol(18) + pol(17) )+ &
!&  ( freq(20) - freq(19) ) * ( pol(20) + pol(19) )+ &
!&  ( freq(22) - freq(21) ) * ( pol(22) + pol(21) ))
!   tmp2 = 0.5_wp * ( &
!&  ( freq (3) - freq (2) ) * ( pol (3) + pol (2) )+ &
!&  ( freq (5) - freq (4) ) * ( pol (5) + pol (4) )+ &
!&  ( freq (7) - freq (6) ) * ( pol (7) + pol (6) )+ &
!&  ( freq (9) - freq (8) ) * ( pol (9) + pol (8) )+ &
!&  ( freq(11) - freq(10) ) * ( pol(11) + pol(10) )+ &
!&  ( freq(13) - freq(12) ) * ( pol(13) + pol(12) )+ &
!&  ( freq(15) - freq(14) ) * ( pol(15) + pol(14) )+ &
!&  ( freq(17) - freq(16) ) * ( pol(17) + pol(16) )+ &
!&  ( freq(19) - freq(18) ) * ( pol(19) + pol(18) )+ &
!&  ( freq(21) - freq(20) ) * ( pol(21) + pol(20) )+ &
!&  ( freq(23) - freq(22) ) * ( pol(23) + pol(22) ))

   trapzd = sum(pol*weights)

end function trapzd

pure elemental function cngw(wf,cn,cnref)
   use iso_fortran_env, only : wp => real64
   implicit none
   real(wp),intent(in) :: wf,cn,cnref
   real(wp)            :: cngw ! CN-gaussian-weight

   intrinsic :: exp

   cngw = exp ( -wf * ( cn - cnref )**2 )

end function cngw

pure elemental function dcngw(wf,cn,cnref)
   use iso_fortran_env, only : wp => real64
!  use dftd4, only : cngw
   implicit none
   real(wp),intent(in) :: wf,cn,cnref
   real(wp) :: dcngw

   dcngw = 2*wf*(cnref-cn)*cngw(wf,cn,cnref)

end function dcngw

!* BJ damping function ala DFT-D3(BJ)
!  f(n,rab) = sn*rab**n/(rab**n + R0**n)  w/ R0 = a1*sqrt(C6/C8)+a2
!  see: https://doi.org/10.1002/jcc.21759
pure elemental function fdmpr_bj(n,r,c) result(fdmp)
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   real(wp) :: fdmp
   fdmp = 1.0_wp / ( r**n + c**n )
end function fdmpr_bj
pure elemental function fdmprdr_bj(n,r,c) result(dfdmp)
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   real(wp) :: dfdmp
   dfdmp = -n*r**(n-1) * fdmpr_bj(n,r,c)**2
end function fdmprdr_bj

!* original DFT-D3(0) damping
!  f(n,rab) = sn/(1+6*(4/3*R0/rab)**alp)  w/ R0 of unknown origin
pure elemental function fdmpr_zero(n,r,c,alp) result(fdmp)
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: fdmp
   fdmp = 1.0_wp / (r**n*(1 + six * (c/r)**(n+alp)))
end function fdmpr_zero
pure elemental function fdmprdr_zero(n,r,c,alp) result(dfdmp)
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: dfdmp
   dfdmp = -( n*r**(n-1)*(1+six*(c/r)**(alp)) &
             - alp*r**n/c*six*(c/r)**(alp-1) ) &
           * fdmpr_zero(n,r,c,alp)**2
!  fdmp = 1.0_wp / (r**n*(1 + 6.0_wp * (c/r)**(n+alp)))
end function fdmprdr_zero

!* fermi damping function from TS and MBD methods
!  f(n,rab) = sn/(1+exp[-alp*(rab/R0-1)]) w/ R0 as experimenal vdW-Radii
pure elemental function fdmpr_fermi(n,r,c,alp) result(fdmp)
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp) :: fdmp
   fdmp = 1.0_wp / (r**n*(1.0_wp+exp(-alp*(r/c - 1.0))))
end function fdmpr_fermi
pure elemental function fdmprdr_fermi(n,r,c,alp) result(dfdmp)
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp) :: dfdmp
   dfdmp = -(-alp/c*r**n*exp(-alp*(r/c - 1.0)) &
             + n*r**(n-1)*(1.0_wp+exp(-alp*(r/c - 1.0)))) &
             * fdmpr_fermi(n,r,c,alp)**2
end function fdmprdr_fermi

!* optimized power zero damping (M. Head-Gordon)
!  f(n,rab) = sn*rab**(n+alp)/(rab**(n+alp) + R0**(n+alp))
!  see: https://dx.doi.org/10.1021/acs.jpclett.7b00176
pure elemental function fdmpr_op(n,r,c,alp) result(fdmp)
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp) :: fdmp
   fdmp = r**alp / (r**(n+alp)*c**(n+alp))
end function fdmpr_op
pure elemental function fdmprdr_op(n,r,c,alp) result(dfdmp)
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   integer, intent(in)  :: alp
   real(wp) :: dfdmp
   dfdmp = (alp*r*(alp-1) - (n+alp)*r**alp*r**(n+alp-1)) &
           * fdmpr_op(n,r,c,alp)**2
!  fdmp = r**alp / (r**(n+alp)*c**(n+alp))
end function fdmprdr_op

!* Sherrill's M-zero damping function
!  f(n,rab) = sn/(1+6*(4/3*R0/rab+a2*R0)**(-alp))
!  see: https://dx.doi.org/10.1021/acs.jpclett.6b00780
pure elemental function fdmpr_zerom(n,r,c,rsn,alp) result(fdmp)
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   real(wp),intent(in)  :: rsn
   integer, intent(in)  :: alp
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: fdmp
   fdmp = 1.0_wp / (r**n*(1 + six * (r/c+rsn*c)**(-alp)))
end function fdmpr_zerom
pure elemental function fdmprdr_zerom(n,r,c,rsn,alp) result(dfdmp)
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: r
   real(wp),intent(in)  :: c
   real(wp),intent(in)  :: rsn
   integer, intent(in)  :: alp
   real(wp),parameter   :: six = 6.0_wp
   real(wp) :: dfdmp
   dfdmp = -( n*r**(n-1)*(1+six*(r/c+rsn*c)**(-alp)) &
              - alp*r**n/c*six*(r/c+rsn*c)**(-alp-1) ) &
           * fdmpr_zerom(n,r,c,rsn,alp)**2
end function fdmprdr_zerom

subroutine d3init(nat,at,ndim)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   integer, intent(out) :: ndim

   integer  :: i,ia,is,icn,j
   integer  :: cncount(0:15)
   real(wp) :: sec_al(23)

   intrinsic :: nint

   ndim = 0
   refal = 0.0_wp

!* set up refc und refal, also obtain the dimension of the dispmat
   do i = 1, nat
      cncount = 0
      ia = at(i)
      refs(ia) = 0
      do j = 1, refn(ia)
         is = refsys(j,ia)
         sec_al = sscale(is)*secaiw(:,is)
         icn =nint(refcn(j,ia))
!        print*,i,j,ia,refs(ia),icn,cncount(icn)
         if (cncount(icn).gt.0) exit
         refs(ia) = refs(ia) + 1
         cncount(icn) = cncount(icn) + 1
         refal(:,j,ia) = max(ascale(j,ia)*(alphaiw(:,j,ia)-hcount(j,ia)*sec_al),0.0_wp)
      enddo
      ndim = ndim + refs(ia)
   enddo

end subroutine d3init

subroutine d4init(nat,at,g_a,g_c,mode,ndim)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: g_a,g_c
   integer, intent(in)  :: mode
   integer, intent(out) :: ndim

   integer  :: i,ia,is,icn,j
   integer  :: cncount(0:15)
   real(wp) :: sec_al(23),iz

   intrinsic :: nint

   select case(mode)
   case(p_refq_hirshfeld,p_refq_periodic)
!     print'(1x,''* using PBE0/def2-TZVP Hirshfeld charges'')'
      refq = dftq
      refh = dfth
      secq = dfts
!  case(2)
!     refq = pbcq
!     refh = pbch
!     secq = pbcs
   case(p_refq_gasteiger)
!     print'(1x,''* using classical Gasteiger charges'')'
      refq = gffq
      refh = gffh
      secq = gffs
   case(p_refq_goedecker)
      refq = clsq
      refh = clsh
      secq = clss
   case(p_refq_gfn2xtb_gbsa_h2o)
!     print'(1x,''* using GFN2-xTB//GBSA(H2O) charges'')'
      refq = solq
      refh = solh
      secq = sols
   end select

   ndim = 0
   refal = 0.0_wp

!* set up refc und refal, also obtain the dimension of the dispmat
   do i = 1, nat
      cncount = 0
      cncount(0) = 1
      ia = at(i)
      do j = 1, refn(ia)
         is = refsys(j,ia)
         iz = zeff(is)
         sec_al = sscale(is)*secaiw(:,is) &
         &  * zeta(g_a,gam(is)*g_c,secq(is)+iz,refh(j,ia)+iz)
         icn =nint(refcn(j,ia))
         cncount(icn) = cncount(icn) + 1
         refal(:,j,ia) = max(ascale(j,ia)*(alphaiw(:,j,ia)-hcount(j,ia)*sec_al),0.0_wp)
      enddo
      do j = 1, refn(ia)
         icn = cncount(nint(refcn(j,ia)))
         refc(j,ia) = icn*(icn+1)/2
      enddo
      ndim = ndim + refn(ia)
   enddo

end subroutine d4init

subroutine d3dim(nat,at,ndim)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   integer, intent(out) :: ndim

   integer :: i

   ndim = 0

   do i = 1, nat
      ndim = ndim + refs(at(i))
   enddo

end subroutine d3dim

subroutine d4dim(nat,at,ndim)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   integer, intent(out) :: ndim

   integer :: i

   ndim = 0

   do i = 1, nat
      ndim = ndim + refn(at(i))
   enddo

end subroutine d4dim

subroutine d3(nat,ndim,at,wf,cn,gw,c6abns)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: wf
   real(wp),intent(in)  :: cn(nat)
   real(wp),intent(out) :: gw(ndim)
   real(wp),intent(out) :: c6abns(ndim,ndim)

   integer  :: i,ia,is,icn,ii,iii,j,jj,ja,k,l
   integer,allocatable :: itbl(:,:)
   real(wp) :: twf,norm,aiw(23)

   intrinsic :: maxval

   allocate( itbl(7,nat), source = 0 )

   gw = 0._wp
   c6abns = 0._wp

   k = 0
   do i = 1, nat
      do ii = 1, refs(at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      norm = 0.0_wp
      do ii = 1, refs(ia)
         norm = norm + cngw(wf,cn(i),refcn(ii,ia))
      enddo
      norm = 1._wp / norm
      do ii = 1, refs(ia)
         k = itbl(ii,i)
         gw(k) = cngw(wf,cn(i),refcn(ii,ia)) * norm
!    --- okay, if we run out of numerical iso_fortran_env, gw(k) will be NaN.
!        In case it is NaN, it will not match itself! So we can rescue
!        this exception. This can only happen for very high CNs.
         if (gw(k).ne.gw(k)) then
            if (maxval(refcn(:refs(ia),ia)).eq.refcn(ii,ia)) then
               gw(k) = 1.0_wp
            else
               gw(k) = 0.0_wp
            endif
         endif
         do j = 1, i-1
            ja = at(j)
            do jj = 1, refs(ja)
               l = itbl(jj,j)
               aiw = refal(:,ii,ia)*refal(:,jj,ja)
               c6abns(l,k) = thopi * trapzd(aiw)
               c6abns(k,l) = c6abns(l,k)
            enddo
         enddo
      enddo
   enddo

end subroutine d3

subroutine d4(nat,ndim,at,wf,g_a,g_c,covcn,gw,c6abns)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: wf,g_a,g_c
   real(wp),intent(in)  :: covcn(nat)
   real(wp),intent(out) :: gw(ndim)
   real(wp),intent(out) :: c6abns(ndim,ndim)

   integer  :: i,ia,is,icn,ii,iii,j,jj,ja,k,l
   integer,allocatable :: itbl(:,:)
   real(wp) :: twf,norm,aiw(23)

   intrinsic :: maxval

   allocate( itbl(7,nat), source = 0 )

   gw = 0._wp
   c6abns = 0._wp

   k = 0
   do i = 1, nat
      do ii = 1, refn(at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      norm = 0.0_wp
      do ii = 1, refn(ia)
         do iii = 1, refc(ii,ia)
            twf = iii*wf
            norm = norm + cngw(twf,covcn(i),refcovcn(ii,ia))
         enddo
      enddo
      norm = 1._wp / norm
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         do iii = 1, refc(ii,ia)
            twf = iii*wf
            gw(k) = gw(k) + cngw(twf,covcn(i),refcovcn(ii,ia)) * norm
         enddo
!    --- okay, if we run out of numerical iso_fortran_env, gw(k) will be NaN.
!        In case it is NaN, it will not match itself! So we can rescue
!        this exception. This can only happen for very high CNs.
         if (gw(k).ne.gw(k)) then
            if (maxval(refcovcn(:refn(ia),ia)).eq.refcovcn(ii,ia)) then
               gw(k) = 1.0_wp
            else
               gw(k) = 0.0_wp
            endif
         endif
         do j = 1, i-1
            ja = at(j)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               aiw = refal(:,ii,ia)*refal(:,jj,ja)
               c6abns(l,k) = thopi * trapzd(aiw)
               c6abns(k,l) = c6abns(l,k)
            enddo
         enddo
      enddo
   enddo

end subroutine d4

pure subroutine build_dispmat(nat,ndim,at,xyz,par,c6abns,dispmat)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   type(dftd_parameter),intent(in)  :: par
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   real(wp),intent(out) :: dispmat(ndim,ndim)

   integer  :: i,ii,ia,j,jj,ja,k,l
   integer, allocatable :: itbl(:,:)
   real(wp) :: c8abns,c10abns,r,r2,oor6,oor8,oor10,cutoff
   real(wp), parameter :: rthr = 72.0_wp ! slightly larger than in gradient

   allocate( itbl(7,nat), source = 0 )

   dispmat = 0.0_wp

   k = 0
   do i = 1, nat
      do ii = 1, refn(at(i))
         k = k + 1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         cutoff = par%a1*sqrt(3._wp*r4r2(at(i))*r4r2(at(j)))+par%a2
         r = norm2(xyz(:,j)-xyz(:,i))
         if (r.gt.rthr) cycle
!        oor6  = 1.0_wp/(r**6  + cutoff**6 )
!        oor8  = 1.0_wp/(r**8  + cutoff**8 )
!        oor10 = 1.0_wp/(r**10 + cutoff**10)
         oor6  = fdmpr_bj( 6,r,cutoff)
         oor8  = fdmpr_bj( 8,r,cutoff)
         oor10 = fdmpr_bj(10,r,cutoff)
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c8abns = 3.0_wp * r4r2(ia)*r4r2(ja) * c6abns(k,l)
               c10abns = 49.0_wp/40.0_wp * c8abns**2/c6abns(k,l)
               dispmat(k,l) = &
               &  - par%s6 * ( c6abns(k,l) * oor6 ) &
               &  - par%s8 * ( c8abns      * oor8 ) &
               &  - par%s8 * ( c10abns     * oor8 )
               dispmat(l,k) = dispmat(k,l)
            enddo
         enddo
      enddo
   enddo

end subroutine build_dispmat

subroutine build_wdispmat(nat,ndim,at,xyz,par,c6abns,gw,wdispmat)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   type(dftd_parameter),intent(in)  :: par
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(out) :: wdispmat(ndim,ndim)

   integer :: i,ii,ia,j,jj,ja,k,l
   integer, allocatable :: itbl(:,:)
   real(wp) :: c8abns,c10abns,r2,cutoff,oor6,oor8,oor10,r,gwgw,r4r2ij
   real(wp), parameter :: rthr = 72.0_wp ! slightly larger than in gradient
   real(wp), parameter :: gwcut = 1.0e-7_wp

   allocate( itbl(7,nat), source = 0 )
 
   wdispmat = 0.0_wp

   k = 0
   do i = 1, nat
      do ii = 1, refn(at(i))
         k = k + 1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         r4r2ij = 3.0_wp*r4r2(ia)*r4r2(ja)
         cutoff = par%a1*sqrt(r4r2ij)+par%a2
!        r2 = sum( (xyz(:,j)-xyz(:,i))**2 )
!        oor6  = 1.0_wp/(r2**3 + cutoff**6 )
!        oor8  = 1.0_wp/(r2**4 + cutoff**8 )
!        oor10 = 1.0_wp/(r2**5 + cutoff**10)
         r = norm2(xyz(:,j)-xyz(:,i))
         if (r.gt.rthr) cycle
         oor6  = fdmpr_bj( 6,r,cutoff)
         oor8  = fdmpr_bj( 8,r,cutoff)
         oor10 = fdmpr_bj(10,r,cutoff)
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               gwgw = gw(k)*gw(l)
               if (gwgw.lt.gwcut) cycle
               c8abns  = r4r2ij * c6abns(k,l)
               c10abns = 49.0_wp/40.0_wp * r4r2ij**2 * c6abns(k,l)
               wdispmat(k,l) = gw(k)*gw(l) * ( &
               &  - par%s6  * ( c6abns(k,l)  * oor6 ) &
               &  - par%s8  * ( c8abns       * oor8 ) &
               &  - par%s10 * ( c10abns      * oor10) )
               wdispmat(l,k) = wdispmat(k,l)
            enddo
         enddo
      enddo
   enddo

end subroutine build_wdispmat

subroutine disppot(nat,ndim,at,q,g_a,g_c,wdispmat,gw,hdisp)
   use iso_fortran_env, only : wp => real64
   use mctc_la, only : symv
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat)
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: wdispmat(ndim,ndim)
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(out) :: hdisp(nat)

   integer  :: i,ii,k,ia
   real(wp) :: qmod,iz
   real(wp),parameter   :: gw_cut = 1.0e-7_wp
   real(wp),allocatable :: zetavec(:)
   real(wp),allocatable :: zerovec(:)
   real(wp),allocatable :: dumvec(:)

   intrinsic :: sum,dble

   allocate( zetavec(ndim),zerovec(ndim),dumvec(ndim), source = 0._wp )

   zetavec = 0.0_wp
   zerovec = 0.0_wp
   dumvec  = 0.0_wp
   hdisp   = 0.0_wp

   k = 0
   do i = 1, nat
       ia = at(i)
       iz = zeff(ia)
       do ii = 1, refn(ia)
          k = k + 1
          if (gw(k).lt.gw_cut) cycle
          zerovec(k) = dzeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz)
          zetavec(k) =  zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz)
      enddo
   enddo
!  create vector -> dispmat(ndim,dnim) * zetavec(ndim) = dumvec(ndim) 
   call symv('U',ndim,1._wp,wdispmat,ndim,zetavec,1,0._wp,dumvec,1)
!  call gemv('N',ndim,ndim,1._wp,wdispmat,ndim,zetavec, &
!  &     1,0._wp,dumvec,1)
!  get atomic reference contributions
   k = 0
   do i = 1, nat
      ia = at(i)
      hdisp(i) = sum(dumvec(k+1:k+refn(ia))*zerovec(k+1:k+refn(ia)))
      k = k + refn(ia)
   enddo

   deallocate(zetavec,zerovec,dumvec)

end subroutine disppot

function edisp_scc(nat,ndim,at,q,g_a,g_c,wdispmat,gw) result(ed)
   use iso_fortran_env, only : wp => real64
   use mctc_la, only : symv
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat)
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: wdispmat(ndim,ndim)
   real(wp),intent(in)  :: gw(ndim)
   real(wp) :: ed

   integer  :: i,ii,k,ia
   real(wp) :: qmod,iz
   real(wp),parameter   :: gw_cut = 1.0e-7_wp
   real(wp),allocatable :: zetavec(:)
   real(wp),allocatable :: dumvec(:)

   intrinsic :: sum,dble

   allocate( zetavec(ndim),dumvec(ndim), source = 0._wp )

   ed = 0.0_wp

   k = 0
   do i = 1, nat
       ia = at(i)
       iz = zeff(ia)
       do ii = 1, refn(ia)
          k = k + 1
          if (gw(k).lt.gw_cut) cycle
          zetavec(k) =  zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz)
      enddo
   enddo
!  create vector -> dispmat(ndim,dnim) * zetavec(ndim) = dumvec(ndim) 
   call symv('U',ndim,0.5_wp,wdispmat,ndim,zetavec,1,0.0_wp,dumvec,1)
!  call gemv('N',ndim,ndim,0.5_wp,wdispmat,ndim,zetavec, &
!  &           1,0.0_wp,dumvec,1)
   ed = dot_product(dumvec,zetavec)

   deallocate(zetavec,dumvec)

end function edisp_scc

subroutine edisp_old(nat,ndim,at,xyz,par,gw,c6abns,mbd,E,aout,etwo,emany)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   integer, intent(in)  :: mbd
   real(wp),intent(out) :: E
   real(wp),intent(out),optional :: aout(23,nat)
   real(wp),intent(out),optional :: etwo
   real(wp),intent(out),optional :: emany

   integer  :: i,ii,ia,k,ij,l,j,jj,ja
   integer, allocatable :: itbl(:,:)
   real(wp) :: Embd,c6ij,c6ij_ns,oor6,oor8,oor10,r2,cutoff,r
   real(wp),allocatable :: dispmat(:,:)
   real(wp),allocatable :: c6ab(:)
   real(wp),allocatable :: aw(:,:)
   real(wp),allocatable :: oor6ab(:,:)

   intrinsic :: present,sqrt,sum
   
   allocate( aw(23,nat),oor6ab(nat,nat),c6ab(nat*(nat+1)/2), &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   e = 0.0_wp

   k = 0
   do i = 1, nat
      do ii = 1, refs(at(i))
         k = k + 1
         itbl(ii,i) = k
         aw(:,i) = aw(:,i) + gw(k) * refal(:,ii,at(i))
      enddo
   enddo

!$OMP parallel private(i,j,ia,ja,ij,k,l,r,oor6,oor8,oor10,cutoff,c6ij,c6ij_ns) &
!$omp&         shared(c6ab) reduction(+:E)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j
!        r = norm2(xyz(:,i)-xyz(:,j))
         r2 = sum( (xyz(:,i)-xyz(:,j))**2 )
         cutoff = par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+par%a2
         oor6 = 1._wp/(r2**3 + cutoff**6)
!        oor6 = fdmpr_bj( 6,r,cutoff)
         oor6ab(i,j) = oor6
         oor6ab(j,i) = oor6
         oor8  = 1._wp/(r2**4 + cutoff**8)
         oor10 = 1._wp/(r2**5 + cutoff**10)
!        oor8  = fdmpr_bj( 8,r,cutoff)
!        oor10 = fdmpr_bj(10,r,cutoff)
         c6ij = 0.0_wp
         do ii = 1, refs(ia)
            k = itbl(ii,i)
            do jj = 1, refs(ja)
               l = itbl(jj,j)
               c6ij = c6ij + gw(k)*gw(l)*c6abns(k,l)
            enddo
         enddo
         c6ab(ij) = c6ij
         E = E - c6ij*(par%s6*oor6 + par%s8*3._wp*r4r2(ia)*r4r2(ja)*oor8 &
         &      + par%s10*49.0_wp/40._wp*(3.0_wp*r4r2(ia)*r4r2(ja))**2*oor10 )
      enddo
   enddo
!$omp enddo
!$omp end parallel

   if (present(Etwo)) Etwo = E

   select case(mbd)
   case(p_mbd_rpalike) ! full RPA-like MBD
!     print'(1x,''* MBD effects calculated by RPA like scheme'')'
      call dispmb(Embd,aw,xyz,oor6ab,nat)
      Embd = par%s9*Embd
      E = E + Embd
   case(p_mbd_exact_atm) ! Axilrod-Teller-Muto three-body term
!     print'(1x,''* MBD effects calculated by ATM formula'')'
      call dispabc(nat,at,xyz,aw,par,Embd)
      E = E + Embd
   case(p_mbd_approx_atm) ! D3-like approximated ATM term
!     print'(1x,''* MBD effects approximated by ATM formula'')'
      call apprabc(nat,at,xyz,c6ab,par,Embd)
      E = E + Embd
   case default
      Embd = 0.0_wp
   end select

   if (present(Emany)) Emany = Embd

   if (present(aout)) then
      aout = 0._wp
      do i = 1, nat
         ia = at(i)
         do ii = 1, refs(ia)
            aout(:,i) = aout(:,i) + gw(k) * refal(:,ii,ia)
         enddo
      enddo
   endif

end subroutine edisp_old

subroutine edisp(nat,ndim,at,q,xyz,par,g_a,g_c, &
           &     gw,c6abns,mbd,E,aout,etwo,emany)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat) 
   real(wp),intent(in)  :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   integer, intent(in)  :: mbd
   real(wp),intent(out) :: E
   real(wp),intent(out),optional :: aout(23,nat)
   real(wp),intent(out),optional :: etwo
   real(wp),intent(out),optional :: emany

   integer  :: i,ii,ia,k,ij,l,j,jj,ja
   integer, allocatable :: itbl(:,:)
   real(wp) :: Embd,qmod,c6ij,c6ij_ns,oor6,oor8,oor10,r2,cutoff,iz,r
   real(wp),allocatable :: dispmat(:,:)
   real(wp),allocatable :: zetvec(:)
   real(wp),allocatable :: zerovec(:)
   real(wp),allocatable :: dumvec(:)
   real(wp),allocatable :: c6ab(:)
   real(wp),allocatable :: aw(:,:)
   real(wp),allocatable :: oor6ab(:,:)

   intrinsic :: present,sqrt,sum
   
   allocate( zetvec(ndim),aw(23,nat),oor6ab(nat,nat), &
   &         zerovec(ndim),c6ab(nat*(nat+1)/2), &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   e = 0.0_wp

   k = 0
   do i = 1, nat
      do ii = 1, refn(at(i))
         k = k + 1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      iz = zeff(ia)
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         zetvec(k) = gw(k) * zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz)
         zerovec(k) = gw(k) * zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,iz)
         aw(:,i) = aw(:,i) + zerovec(k) * refal(:,ii,ia)
      enddo
   enddo

!$OMP parallel private(i,j,ia,ja,ij,k,l,r,oor6,oor8,oor10,cutoff,c6ij,c6ij_ns) &
!$omp&         shared(c6ab) reduction(+:E)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j
         r = norm2(xyz(:,i)-xyz(:,j))
!        r2 = sum( (xyz(:,i)-xyz(:,j))**2 )
         cutoff = par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+par%a2
!        oor6 = 1._wp/(r2**3 + cutoff**6)
         oor6 = fdmpr_bj( 6,r,cutoff)
         oor6ab(i,j) = oor6
         oor6ab(j,i) = oor6
!        oor8  = 1._wp/(r2**4 + cutoff**8)
!        oor10 = 1._wp/(r2**5 + cutoff**10)
         oor8  = fdmpr_bj( 8,r,cutoff)
         oor10 = fdmpr_bj(10,r,cutoff)
         c6ij_ns = 0.0_wp
         c6ij = 0.0_wp
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c6ij_ns = c6ij_ns + zerovec(k)*zerovec(l)*c6abns(k,l)
               c6ij = c6ij + zetvec(k)*zetvec(l)*c6abns(k,l)
            enddo
         enddo
         c6ab(ij) = c6ij_ns
         E = E - c6ij*(par%s6*oor6 + par%s8*3._wp*r4r2(ia)*r4r2(ja)*oor8 &
         &      + par%s10*49.0_wp/40._wp*(3.0_wp*r4r2(ia)*r4r2(ja))**2*oor10 )
      enddo
   enddo
!$omp enddo
!$omp end parallel

   if (present(Etwo)) Etwo = E

   select case(mbd)
   case(p_mbd_rpalike) ! full RPA-like MBD
!     print'(1x,''* MBD effects calculated by RPA like scheme'')'
      call dispmb(Embd,aw,xyz,oor6ab,nat)
      Embd = par%s9*Embd
      E = E + Embd
   case(p_mbd_exact_atm) ! Axilrod-Teller-Muto three-body term
!     print'(1x,''* MBD effects calculated by ATM formula'')'
      call dispabc(nat,at,xyz,aw,par,Embd)
      E = E + Embd
   case(p_mbd_approx_atm) ! D3-like approximated ATM term
!     print'(1x,''* MBD effects approximated by ATM formula'')'
      call apprabc(nat,at,xyz,c6ab,par,Embd)
      E = E + Embd
   case default
      Embd = 0.0_wp
   end select

   if (present(Emany)) Emany = Embd

   if (present(aout)) then
      aout = 0._wp
      do i = 1, nat
         ia = at(i)
         do ii = 1, refn(ia)
            aout(:,i) = aout(:,i) + zetvec(k) * refal(:,ii,ia)
         enddo
      enddo
   endif

end subroutine edisp

!* compute D4 gradient
!  ∂E/∂rij = ∂/∂rij (W·D·W)
!          = ∂W/∂rij·D·W  + W·∂D/∂rij·W + W·D·∂W/∂rij
!  ∂W/∂rij = ∂(ζ·w)/∂rij = ζ·∂w/∂rij = ζ·∂w/∂CN·∂CN/∂rij
!  ∂ζ/∂rij = 0
subroutine dispgrad(nat,ndim,at,q,xyz, &
           &        par,wf,g_a,g_c, &
           &        c6abns,mbd, &
           &        g,eout,dq,aout)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat)
!  real(wp),intent(in)  :: cn(nat) ! calculate on-the-fly
   real(wp),intent(in)  :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: wf,g_a,g_c
!  real(wp),intent(in)  :: gw(ndim) ! calculate on-the-fly
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   integer, intent(in)  :: mbd
   real(wp),intent(inout)        :: g(3,nat)
   real(wp),intent(out),optional :: eout
   real(wp),intent(in), optional :: dq(3,nat,nat+1)
   real(wp),intent(out),optional :: aout(23,nat)


   integer  :: i,ii,iii,j,jj,k,l,ia,ja,ij
   integer, allocatable :: itbl(:,:)
   real(wp) :: iz
   real(wp) :: qmod,eabc,ed
   real(wp) :: norm,dnorm
   real(wp) :: dexpw,expw
   real(wp) :: twf,tgw,r4r2ij
   real(wp) :: rij(3),r,r2,r4,r6,r8,R0
   real(wp) :: oor6,oor8,oor10,door6,door8,door10
   real(wp) :: c8abns,disp,x1,x2,x3
   real(wp) :: c6ij,dic6ij,djc6ij,dizij,djzij
   real(wp) :: rcovij,expterm,den,dcndr
   real(wp) :: drdx(3),dtmp,gwk,dgwk
   real(wp),allocatable :: r2ab(:)
   real(wp),allocatable :: dc6dr(:)
   real(wp),allocatable :: dc6dcn(:)
   real(wp),allocatable :: zvec(:)
   real(wp),allocatable :: dzvec(:)
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: gw(:)
   real(wp),allocatable :: dgw(:)
   real(wp),allocatable :: dc6dq(:)
   real(wp),allocatable :: dzdq(:)
   real(wp) :: cn_thr,r_thr,gw_thr
   parameter(cn_thr = 1600.0_wp)
   parameter( r_thr=5000._wp)
   parameter(gw_thr=0.000001_wp)
   real(wp),parameter :: sqrtpi = 1.77245385091_wp
   real(wp),parameter :: hlfosqrtpi = 0.5_wp/1.77245385091_wp
!  timing
!  real(wp) :: time0,time1
!  real(wp) :: wall0,wall1

   intrinsic :: present,sqrt,sum,maxval,exp,abs

!  print'(" * Allocating local memory")'
   allocate( dc6dr(nat*(nat+1)/2),dc6dcn(nat),  &
   &         r2ab(nat*(nat+1)/2),cn(nat),  &
   &         zvec(ndim),dzvec(ndim),  &
   &         gw(ndim),dgw(ndim),dc6dq(nat),dzdq(ndim),  &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   ed = 0.0_wp

!  derivatives of the covalent coordination number
!  ∂CN/∂rij = ∂/∂rij k4·exp(-(∂EN+k5)²/(2·k6²))/(1+exp(-k1(-1)))
!  print'(" * Calculating CN")'
   call covncoord(nat,at,xyz,cn,cn_thr)

!  precalc
!  print'(" * Setting up index table")'
   k = 0
   do i = 1, nat
      do ii = 1, refn(at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

!  print'(" * Entering first OMP section")'
!$OMP parallel default(none) &
!$omp private(i,ii,iii,ia,iz,k,norm,dnorm,twf,tgw,dexpw,expw,gwk,dgwk)  &
!$omp shared (nat,at,refn,refc,refcovcn,itbl,refq,wf,cn,g_a,g_c,q) &
!$omp shared (gw,dgw,zvec,dzvec,dzdq)
!$omp do
   do i = 1, nat
      ia = at(i)
      iz = zeff(ia)
      norm  = 0.0_wp
      dnorm = 0.0_wp
      do ii=1,refn(ia)
         do iii = 1, refc(ii,ia)
            twf = iii*wf
            tgw = cngw(twf,cn(i),refcovcn(ii,ia))
            norm  =  norm + tgw
            dnorm = dnorm + 2*twf*(refcovcn(ii,ia)-cn(i))*tgw
         enddo
      enddo
      norm = 1._wp/norm
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         dexpw=0.0_wp
         expw=0.0_wp
         do iii = 1, refc(ii,ia)
            twf = wf*iii
            tgw = cngw(twf,cn(i),refcovcn(ii,ia))
            expw  =  expw + tgw
            dexpw = dexpw + 2*twf*(refcovcn(ii,ia)-cn(i))*tgw
         enddo

         ! save
         gwk = expw*norm
         if (gwk.ne.gwk) then
            if (maxval(refcovcn(:refn(ia),ia)).eq.refcovcn(ii,ia)) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         zvec(k) = zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz) * gwk
         ! NEW: q=0 for ATM
         gw(k) =  zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,iz) * gwk

         dgwk = dexpw*norm-expw*dnorm*norm**2
         if (dgwk.ne.dgwk) then
            dgwk = 0.0_wp
         endif
         dzvec(k) = zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz) * dgwk
         dzdq(k) = dzeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz) * gwk
         ! NEW: q=0 for ATM
         dgw(k) = zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,iz) * dgwk
      enddo
   enddo
!$omp end do
!$omp end parallel

!  print'(" * Entering second OMP section")'
!$OMP parallel default(none) &
!$omp private(i,j,ia,ja,ij,k,l,c6ij,dic6ij,djc6ij,disp,dizij,djzij,  &
!$omp         rij,r2,r,r4r2ij,r0,oor6,oor8,oor10,door6,door8,door10)  &
!$omp shared(nat,at,xyz,refn,itbl,zvec,dzvec,c6abns,par,dzdq) &
!$omp shared(r2ab) reduction(+:dc6dr,dc6dq,dc6dcn,ed)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j
         rij = xyz(:,j) - xyz(:,i)
         r2 = sum( rij**2 )
         r2ab(ij) = r2
         if (r2.gt.r_thr) cycle
         ! temps
         c6ij = 0.0_wp
         dic6ij = 0.0_wp
         djc6ij = 0.0_wp
         dizij = 0.0_wp
         djzij = 0.0_wp
         ! all refs
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*c6abns(k,l)
               dic6ij = dic6ij + dzvec(k)*zvec(l)*c6abns(k,l)
               djc6ij = djc6ij + zvec(k)*dzvec(l)*c6abns(k,l)
               dizij = dizij + dzdq(k)*zvec(l)*c6abns(k,l)
               djzij = djzij + zvec(k)*dzdq(l)*c6abns(k,l)
            enddo
         enddo

         r = sqrt(r2)

         r4r2ij = 3*r4r2(ia)*r4r2(ja)
         r0 = par%a1*sqrt(r4r2ij) + par%a2

         oor6 = 1._wp/(r2**3+r0**6)
         oor8 = 1._wp/(r2**4+r0**8)
         oor10 = 1._wp/(r2**5+r0**10)
         door6 = -6*r2**2*r*oor6**2
         door8 = -8*r2**3*r*oor8**2
         door10 = -10*r2**4*r*oor10**2
!        oor6   = fdmpr_bj( 6,r,r0)
!        oor8   = fdmpr_bj( 8,r,r0)
!        oor10  = fdmpr_bj(10,r,r0)
!        door6  = fdmprdr_bj( 6,r,r0)
!        door8  = fdmprdr_bj( 8,r,r0)
!        door10 = fdmprdr_bj(10,r,r0)

         disp = par%s6*oor6 + par%s8*r4r2ij*oor8 &
         &    + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10

         ! save
         dc6dq(i) = dc6dq(i) + dizij*disp
         dc6dq(j) = dc6dq(j) + djzij*disp
         dc6dcn(i) = dc6dcn(i) + dic6ij*disp
         dc6dcn(j) = dc6dcn(j) + djc6ij*disp
         dc6dr(ij) = dc6dr(ij) + c6ij*(par%s6*door6 + par%s8*r4r2ij*door8 &
         &                       + par%s10*49.0_wp/40.0_wp*r4r2ij**2*door10 )

         ed = ed - c6ij*disp
      enddo
   enddo
!$omp enddo
!$omp end parallel

!  select case(mbd)
!  case(1) ! full RPA-like MBD
!     print'(1x,''* MBD effects calculated by RPA like scheme'')'
!     call raise('W','MBD gradient not fully implemented yet')
!     call mbdgrad(nat,xyz,aw,daw,oor6ab,g,embd)
!  case(1,2) ! Axilrod-Teller-Muto three-body term
!     if(mbd.eq.1) then
!        call raise('W','MBD gradient not fully implemented yet')
!        print'(''MBD gradient not fully implemented yet'')'
!        print'(1x,''* calculate MBD effects with ATM formula instead'')'
!     else
!        print'(1x,''* MBD effects calculated by ATM formula'')'
!     endif
!     call dabcgrad(nat,ndim,at,xyz,par,dcn,gw,dgw,itbl,g,embd)
!  case(3) ! D3-like approximated ATM term
!     print'(1x,''* MBD effects approximated by ATM formula'')'
!  print'(" * Starting MBD gradient calculation")'
   if (mbd.ne.p_mbd_none) &
   &   call dabcappr(nat,ndim,at,xyz,par,  &
           &        r2ab,gw,dgw,c6abns,itbl,dc6dr,dc6dcn,eabc)
!  end select

!  print'(" * Entering third OMP section")'
!$OMP parallel default(none) &
!$omp private(i,j,ia,ja,ij,rij,r2,r,drdx,den,rcovij,  &
!$omp         expterm,dcndr,dtmp) reduction(+:g) &
!$omp shared(nat,at,xyz,dc6dr,dc6dcn)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j

         rij = xyz(:,j) - xyz(:,i)
         r2 = sum( rij**2 )
         r = sqrt(r2)
         drdx = rij/r

         if (r2.lt.cn_thr) then
            den=k4*exp(-(abs((en(ia)-en(ja)))+ k5)**2/k6 )
            rcovij=rcov(ia)+rcov(ja)
            !expterm=exp(-k1*(rcovij/r-1._wp))
            !dcndr=(-k1*rcovij*expterm*den)/(r2*((1._wp + expterm)**2))
            dcndr = -den*kn/sqrtpi/rcovij*exp(-kn**2*(r-rcovij)**2/rcovij**2)
         else
            dcndr=0.0_wp
         endif

         dtmp = dc6dr(ij) + (dc6dcn(i)+dc6dcn(j))*dcndr

         g(:,i) = g(:,i) + dtmp * drdx
         g(:,j) = g(:,j) - dtmp * drdx
      enddo
   enddo
!$omp enddo
!$omp end parallel

   if (present(dq)) then
      call dgemv('n',3*nat,nat,-1.0_wp,dq,3*nat,dc6dq,1,1.0_wp,g,1)
   endif

!  print*,ed,eabc

!  print'(" * Dispersion all done, saving variables")'
   if (present(eout)) eout = ed + eabc

   if (present(aout)) then
      aout = 0._wp
      do i = 1, nat
         ia = at(i)
         do ii = 1, refn(ia)
            aout(:,i) = aout(:,i) + zvec(k) * refal(:,ii,ia)
         enddo
      enddo
   endif


end subroutine dispgrad

subroutine apprabc(nat,at,xyz,c6ab,par,E)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: c6ab(nat*(nat+1)/2)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(out) :: E

   integer  :: i,j,k,ia,ja,ka,ij,ik,jk
   real(wp) :: rij(3),rjk(3),rik(3),r2ij,r2jk,r2ik,cij,cjk,cik,cijk
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk,fdmp
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
!  parameter(thf = 3._wp/4._wp)

   intrinsic :: sum,sqrt

   E = 0.0_wp

   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j
         rij  = xyz(:,j) - xyz(:,i)
         r2ij = sum(rij**2)
         cij  = par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+par%a2
         do k = 1, j-1
            ka = at(k)
            ik = i*(i-1)/2 + k
            jk = j*(j-1)/2 + k
            rik   = xyz(:,i) - xyz(:,k)
            r2ik  = sum(rik**2)
            cik   = par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ka))+par%a2
            rjk   = xyz(:,k) - xyz(:,j)
            r2jk  = sum(rjk**2)
            cjk   = par%a1*sqrt(3._wp*r4r2(ja)*r4r2(ka))+par%a2
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik
            c9ijk = par%s9*sqrt(c6ab(ij)*c6ab(jk)*c6ab(ik))
            atm = ( 0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp
            fdmp = one/(one+six*((cijk/rijk)**oth)**par%alp)
            oor9ijk = atm/rijk**3*fdmp
            E = E + c9ijk * oor9ijk
         enddo
      enddo
   enddo

end subroutine apprabc

subroutine dispabc(nat,at,xyz,aw,par,E)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: aw(23,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(out) :: E

   integer  :: i,j,k,ia,ja,ka
   real(wp) :: rij(3),rjk(3),rik(3),r2ij,r2jk,r2ik,cij,cjk,cik,cijk
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk,fdmp
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
   parameter(thf = 3._wp/4._wp)

   intrinsic :: sum,sqrt

   E = 0.0_wp

   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         rij  = xyz(:,j) - xyz(:,i)
         r2ij = sum(rij**2)
         cij  = (par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+par%a2)
         do k = 1, j-1
            ka = at(k)
            rik   = xyz(:,i) - xyz(:,k)
            r2ik  = sum(rik**2)
            cik   = (par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ka))+par%a2)
            rjk   = xyz(:,k) - xyz(:,j)
            r2jk  = sum(rjk**2)
            cjk   = (par%a1*sqrt(3._wp*r4r2(ja)*r4r2(ka))+par%a2)
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik
            c9ijk = par%s9*thopi*trapzd( aw(:,i)*aw(:,j)*aw(:,k) )
            atm = ( 0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp
            fdmp = one/(one+six*(thf*(cijk/rijk)**oth)**par%alp)
            oor9ijk = atm/rijk**3*fdmp
            E = E + c9ijk * oor9ijk
         enddo
      enddo
   enddo

end subroutine dispabc

subroutine abcappr(nat,ndim,at,xyz,g_a,g_c,par,gw,r2ab,c6abns,eabc)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: g_a,g_c
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: r2ab(nat*(nat+1)/2)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   real(wp),intent(out) :: eabc

   integer  :: i,ii,ia,j,jj,ja,k,kk,ka,l,m,n
   integer  :: ij,jk,ik
   integer, allocatable :: itbl(:,:)
   real(wp),allocatable :: c6ab(:),zvec(:),c(:)
   real(wp) :: r2ij,r2jk,r2ik,iz
   real(wp) :: cij,cjk,cik,cijk
   real(wp) :: fdmp,dtmp,oor9tmp,c9tmp
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk
   real(wp) :: drij(3),drjk(3),drik(3)
   real(wp) :: oorij,oorjk,oorik
   real(wp) :: dijfdmp,dikfdmp,djkfdmp
   real(wp) :: dijatm,dikatm,djkatm
   real(wp) :: dijoor9ijk,djkoor9ijk,dikoor9ijk
   real(wp) :: c6ij,dic6ij,djc6ij
   real(wp) :: dic9ijk,djc9ijk,dkc9ijk
   real(wp) :: x1,x2,x3,x4,x5,x6,x7,x8,x9
   real(wp) :: dum1,dum2,dum3
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
!  parameter(thf = 3._wp/4._wp)
   real(wp) :: r_thr,gw_thr
   parameter( r_thr=1600._wp)
   parameter(gw_thr=0.0001_wp)

   intrinsic :: sqrt

   allocate( c6ab(nat*(nat+1)/2),zvec(ndim),c(nat*(nat+1)/2),  &
   &         source = 0.0_wp )
   allocate( itbl(7,nat), source = 0 )

   eabc = 0.0_wp

!  precalc
   k = 0
   do i = 1, nat
      ia = at(i)
      iz = zeff(ia)
      do ii = 1, refn(ia)
         k = k+1
         itbl(ii,i) = k
         ! NEW: q=0 for ATM
         zvec(k) = zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,iz) * gw(k)
      enddo
   enddo

!$OMP parallel private(i,ia,j,ja,ij,r2ij,c6ij)  &
!$omp&         shared (c6ab,c)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
!        if(i.eq.j) cycle
         ja = at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         r2ij = r2ab(ij)
         c(ij) = par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+par%a2
         if(r2ij.gt.r_thr) cycle

         ! temps
         c6ij = 0.0_wp
         ! all refs
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*c6abns(k,l)
            enddo
         enddo
         ! save
         c6ab(ij) = c6ij
      enddo
   enddo
!$omp enddo
!$omp end parallel

!$OMP parallel private(i,j,ij,ia,ja,k,ka,ik,jk,atm,fdmp,  &
!$omp&                 r2ij,cij,r2ik,r2jk,cik,cjk,r2ijk,rijk,cijk, &
!$omp&                 c9ijk,oor9ijk) &
!$omp&         reduction(+:eabc)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         r2ij = r2ab(ij)
         if(r2ij.gt.r_thr) cycle
         cij  = c(ij)
         do k = 1, j-1
            ka = at(k)
            ik = i*(i-1)/2 + k
            jk = j*(j-1)/2 + k
            r2ik  = r2ab(ik)
            r2jk  = r2ab(jk)
            if((r2ik.gt.r_thr).or.(r2jk.gt.r_thr)) cycle
            cik   = c(ik)
            cjk   = c(jk)
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik

            atm = ((0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp)/(rijk**3)

            fdmp = one/(one+six*((cijk/rijk)**oth)**par%alp)

            c9ijk = par%s9*sqrt(c6ab(ij)*c6ab(jk)*c6ab(ik))

            oor9ijk = atm*fdmp
            eabc = eabc + c9ijk*oor9ijk

         enddo ! k/C
      enddo ! j/B
   enddo ! i/A
!$omp enddo
!$omp end parallel

   deallocate( c6ab,c,zvec )

end subroutine abcappr

subroutine dabcappr(nat,ndim,at,xyz,par,  &
                &        r2ab,zvec,dzvec,c6abns,itbl,dc6dr,dc6dcn,eout)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: r2ab(nat*(nat+1)/2)
   real(wp),intent(in)  :: zvec(ndim)
   real(wp),intent(in)  :: dzvec(ndim)
   real(wp),intent(in)  :: c6abns(ndim,ndim)
   integer, intent(in)  :: itbl(7,nat)
   real(wp),intent(inout)        :: dc6dr(nat*(nat+1)/2)
   real(wp),intent(inout)        :: dc6dcn(nat)
   real(wp),intent(out),optional :: eout

   integer  :: i,ii,ia,j,jj,ja,k,kk,ka,l,m,n
   integer  :: ij,jk,ik
   real(wp),allocatable :: c6ab(:),dc6ab(:,:)
   real(wp) :: r2ij,r2jk,r2ik
   real(wp) :: cij,cjk,cik,cijk
   real(wp) :: fdmp,dtmp,oor9tmp,c9tmp
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk
   real(wp) :: drij(3),drjk(3),drik(3)
   real(wp) :: oorij,oorjk,oorik
   real(wp) :: dijfdmp,dikfdmp,djkfdmp
   real(wp) :: dijatm,dikatm,djkatm
   real(wp) :: dijoor9ijk,djkoor9ijk,dikoor9ijk
   real(wp) :: c6ij,dic6ij,djc6ij
   real(wp) :: dic9ijk,djc9ijk,dkc9ijk
   real(wp) :: x1,x2,x3,x4,x5,x6,x7,x8,x9
   real(wp) :: eabc
   real(wp) :: dum1,dum2,dum3
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
!  parameter(thf = 3._wp/4._wp)
   real(wp) :: r_thr,gw_thr
   parameter( r_thr=1600._wp)
   parameter(gw_thr=0.0001_wp)

   intrinsic :: present,sqrt

   allocate( c6ab(nat*(nat+1)/2),dc6ab(nat,nat),  &
   &         source = 0.0_wp )

   eabc = 0.0_wp

!$OMP parallel default(none) &
!$omp private(i,ia,j,ja,ij,r2ij,c6ij,dic6ij,djc6ij,k,l)  &
!$omp shared (nat,at,r2ab,refn,itbl,c6abns,zvec,dzvec) &
!$omp shared (c6ab,dc6ab)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
!        if(i.eq.j) cycle
         ja = at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         r2ij = r2ab(ij)
         if(r2ij.gt.r_thr) cycle

         ! temps
         c6ij = 0.0_wp
         dic6ij = 0.0_wp
         djc6ij = 0.0_wp
         ! all refs
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*c6abns(k,l)
               dic6ij = dic6ij + dzvec(k)*zvec(l)*c6abns(k,l)
               djc6ij = djc6ij + zvec(k)*dzvec(l)*c6abns(k,l)
            enddo
         enddo
         ! save
         c6ab(ij) = c6ij
         dc6ab(i,j) = dic6ij
         dc6ab(j,i) = djc6ij
      enddo
   enddo
!$omp enddo
!$omp end parallel

!$OMP parallel default(none) &
!$omp private(i,j,ij,ia,ja,k,ka,ik,jk,oorjk,oorik,atm,fdmp,  &
!$omp         r2ij,cij,oorij,r2ik,r2jk,cik,cjk,r2ijk,rijk,cijk, &
!$omp         dijatm,djkatm,dikatm,dtmp,dijfdmp,djkfdmp,dikfdmp,  &
!$omp         c9ijk,oor9ijk,dic9ijk,djc9ijk,dkc9ijk) &
!$omp shared(nat,at,r2ab,par,c6ab,dc6ab) &
!$omp reduction(+:eabc,dc6dr) reduction(-:dc6dcn)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         r2ij = r2ab(ij)
         if(r2ij.gt.r_thr) cycle
         cij  = par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+par%a2
         oorij = 1._wp/sqrt(r2ij)
         do k = 1, j-1
            ka = at(k)
            ik = i*(i-1)/2 + k
            jk = j*(j-1)/2 + k
            r2ik  = r2ab(ik)
            r2jk  = r2ab(jk)
            if((r2ik.gt.r_thr).or.(r2jk.gt.r_thr)) cycle
            cik   = par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ka))+par%a2
            cjk   = par%a1*sqrt(3._wp*r4r2(ja)*r4r2(ka))+par%a2
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik
            oorjk = 1._wp/sqrt(r2jk)
            oorik = 1._wp/sqrt(r2ik)

            atm = ((0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp)/(rijk**3)
            dijatm=-0.375_wp*(r2ij**3+r2ij**2*(r2jk+r2ik) &
            &      +r2ij*(3._wp*r2jk**2+2._wp*r2jk*r2ik+3._wp*r2ik**2) &
            &      -5._wp*(r2jk-r2ik)**2*(r2jk+r2ik)) &
            &      /(r2ijk*rijk**3)*oorij
            djkatm=-0.375_wp*(r2jk**3+r2jk**2*(r2ik+r2ij) &
            &      +r2jk*(3._wp*r2ik**2+2._wp*r2ik*r2ij+3._wp*r2ij**2) &
            &      -5._wp*(r2ik-r2ij)**2*(r2ik+r2ij)) &
            &      /(r2ijk*rijk**3)*oorjk
            dikatm=-0.375_wp*(r2ik**3+r2ik**2*(r2jk+r2ij) &
            &      +r2ik*(3._wp*r2jk**2+2._wp*r2jk*r2ij+3._wp*r2ij**2) &
            &      -5._wp*(r2jk-r2ij)**2*(r2jk+r2ij)) &
            &      /(r2ijk*rijk**3)*oorik

            fdmp = one/(one+six*((cijk/rijk)**oth)**par%alp)
            dtmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2
            dijfdmp = dtmp*oorij
            djkfdmp = dtmp*oorjk
            dikfdmp = dtmp*oorik

            c9ijk = par%s9*sqrt(c6ab(ij)*c6ab(jk)*c6ab(ik))

            oor9ijk = atm*fdmp
            eabc = eabc + c9ijk*oor9ijk

            dc6dr(ij) = dc6dr(ij) + (atm*dijfdmp - dijatm*fdmp)*c9ijk
            dc6dr(ik) = dc6dr(ik) + (atm*dikfdmp - dikatm*fdmp)*c9ijk
            dc6dr(jk) = dc6dr(jk) + (atm*djkfdmp - djkatm*fdmp)*c9ijk
            dic9ijk = dc6ab(i,j)/c6ab(ij) + dc6ab(i,k)/c6ab(ik)
            djc9ijk = dc6ab(j,i)/c6ab(ij) + dc6ab(j,k)/c6ab(jk)
            dkc9ijk = dc6ab(k,j)/c6ab(jk) + dc6ab(k,i)/c6ab(ik)
            dc6dcn(i) = dc6dcn(i) - 0.5_wp*c9ijk*oor9ijk*dic9ijk
            dc6dcn(j) = dc6dcn(j) - 0.5_wp*c9ijk*oor9ijk*djc9ijk
            dc6dcn(k) = dc6dcn(k) - 0.5_wp*c9ijk*oor9ijk*dkc9ijk

         enddo ! k/C
      enddo ! j/B
   enddo ! i/A
!$omp enddo
!$omp end parallel

   if (present(eout)) eout=eabc

end subroutine dabcappr


!* here is the theory for the ATM-gradient (SAW, 180224)
! EABC = WA·WB·WC·DABC
! ∂EABC/∂X = ∂/∂X(WA·WB·WC·DABC)
!          = ∂WA/∂X·WB·WC·DABC + WA·∂WB/∂X·WC·DABC + WA·WB·∂WC/∂X·WC·DABC
!            + WA·WB·WC·∂DABC/∂X
! ∂/∂X =  ∂rAB/∂X·∂/∂rAB +  ∂rBC/∂X·∂/∂rBC +  ∂rCA/∂X·∂/∂rCA
!      = (δAX-δBX)∂/∂rAB + (δBX-δCX)∂/∂rBC + (δCX-δAX)∂/∂rCA
! ∂EABC/∂A = ∑A,ref ∑B,ref ∑C,ref
!            + (∂WA/∂rAB-∂WA/∂rCA)·WB·WC·DABC
!            + WA·∂WB/∂rAB·WC·DABC
!            - WA·WB·∂WC/∂rCA·DABC
!            + WA·WB·WC·(∂DABC/∂rAB-∂DABC/∂rCA)
! ∂EABC/∂B = ∑A,ref ∑B,ref ∑C,ref
!            - ∂WA/∂rAB·WB·WC·DABC
!            + WA·(∂WB/∂rBC-∂WB/∂rAB)·WC·DABC
!            + WA·WB·∂WC/∂rBC·DABC
!            + WA·WB·WC·(∂DABC/∂rBC-∂DABC/∂rAB)
! ∂EABC/∂C = ∑A,ref ∑B,ref ∑C,ref
!            + ∂WA/∂rCA·WB·WC·DABC
!            - WA·∂WB/∂rBC·WC·DABC
!            + WA·WB·(∂WC/∂rCA-∂WC/∂rBC)·DABC
!            + WA·WB·WC·(∂DABC/∂rCA-∂DABC/∂rBC)
! ∂WA/∂rAB = ∂CNA/∂rAB·∂WA/∂CNA w/ ζ=1 and WA=wA
! ATM = 3·cos(α)cos(β)cos(γ)+1
!     = 3/8(r²AB+r²BC-r²CA)(r²AB+r²CA-r²BC)(r²BC+r²CA-r²AB)/(r²BC·r²CA·r²AB)+1
! ∂ATM/∂rAB = 3/4(2r⁶AB-r⁶BC-r⁶CA-r⁴AB·r²BC-r⁴AB·r²CA+r⁴BC·r²CA+r²BC·r⁴CA)
!             /(r³AB·r²BC·r²CA)
! DABC = C9ABCns·f·ATM/(rAB·rBC·rCA)³
! f = 1/(1+6(¾·∛[RAB·RBC·RCA/(rAB·rBC·rCA)])¹⁶)
! ∂(f/(r³AB·r³BC·r³CA)/∂rAB = 
!   ⅓·((6·(16-9)·(¾·∛[RAB·RBC·RCA/(rAB·rBC·rCA)])¹⁶-9)·f²/(r⁴AB·r³BC·r³CA)
subroutine dabcgrad(nat,ndim,at,xyz,par,dcn,zvec,dzvec,itbl,g,eout)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: dcn(nat,nat)
   real(wp),intent(in)  :: zvec(ndim)
   real(wp),intent(in)  :: dzvec(ndim)
   integer, intent(in)  :: itbl(7,nat)
   real(wp),intent(inout)        :: g(3,nat)
   real(wp),intent(out),optional :: eout

   integer  :: i,ii,ia,j,jj,ja,k,kk,ka,l,m,n
   real(wp) :: rij(3),rjk(3),rik(3)
   real(wp) :: r2ij,r2jk,r2ik
   real(wp) :: cij,cjk,cik,cijk
   real(wp) :: fdmp,dtmp,oor9tmp,c9tmp
   real(wp) :: atm,r2ijk,c9ijk,oor9ijk,rijk
   real(wp) :: drij(3),drjk(3),drik(3)
   real(wp) :: oorij,oorjk,oorik
   real(wp) :: dijfdmp,dikfdmp,djkfdmp
   real(wp) :: dijatm,dikatm,djkatm
   real(wp) :: dijoor9ijk,djkoor9ijk,dikoor9ijk
   real(wp) :: x1,x2,x3,x4,x5,x6,x7,x8,x9
   real(wp) :: eabc
   real(wp) :: dum1,dum2,dum3
   real(wp) :: one,oth,six,thf
   parameter(one = 1._wp)
   parameter(oth = 1._wp/3._wp)
   parameter(six = 6._wp)
   parameter(thf = 3._wp/4._wp)
   real(wp) :: r_thr,gw_thr
   parameter( r_thr=1600._wp)
   parameter(gw_thr=0.0001_wp)

   intrinsic :: present,sqrt,sum

   eabc = 0._wp

!$omp parallel private(ia,ja,ka,l,m,n, &
!$omp&         rij,rjk,rik,r2ij,r2jk,r2ik, &
!$omp&         cij,cjk,cik,cijk, &
!$omp&         fdmp,dtmp,oor9tmp,c9tmp, &
!$omp&         atm,r2ijk,c9ijk,oor9ijk,rijk, &
!$omp&         drij,drjk,drik,oorij,oorjk,oorik, &
!$omp&         dijfdmp,dikfdmp,djkfdmp, &
!$omp&         dijatm,dikatm,djkatm, &
!$omp&         dijoor9ijk,djkoor9ijk,dikoor9ijk, &
!$omp&         x1,x2,x3,x4,x5,x6,x7,x8,x9) &
!$omp&         reduction(+:g,eabc)
!$omp do schedule(dynamic)
   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
!        if(i.eq.j) cycle
         ja = at(j)
!    --- all distances, cutoff radii ---
         rij  = xyz(:,j) - xyz(:,i)
         r2ij = sum(rij**2)
         if(r2ij.gt.r_thr) cycle
         cij  = (par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ja))+par%a2)
         oorij = 1._wp/sqrt(r2ij)
         do k = 1, j-1
!           if(k.eq.j) cycle
!           if(i.eq.k) cycle
            ka = at(k)
            rik   = xyz(:,i) - xyz(:,k)
            rjk   = xyz(:,k) - xyz(:,j)
            r2ik  = sum(rik**2)
            r2jk  = sum(rjk**2)
            if((r2ik.gt.r_thr).or.(r2jk.gt.r_thr)) cycle
            cik   = (par%a1*sqrt(3._wp*r4r2(ia)*r4r2(ka))+par%a2)
            cjk   = (par%a1*sqrt(3._wp*r4r2(ja)*r4r2(ka))+par%a2)
            r2ijk = r2ij*r2ik*r2jk
            rijk  = sqrt(r2ijk)
            cijk  = cij*cjk*cik
            oorjk = 1._wp/sqrt(r2jk)
            oorik = 1._wp/sqrt(r2ik)

            x2 = 0._wp
            x4 = 0._wp
            x6 = 0._wp
            c9ijk = 0._wp

!       --- sum up all references ---
            do ii = 1, refn(ia) ! refs of A
               l = itbl(ii,i)
               do jj = 1, refn(ja) ! refs of B
                  m = itbl(jj,j)
                  do kk = 1, refn(ka) ! refs of C
                     n = itbl(kk,k)
                     if ((zvec(l)*zvec(m)*zvec(n)).lt.gw_thr) cycle
                     c9tmp = par%s9*thopi*trapzd(refal(:,ii,ia)*refal(:,jj,ja) &
                     &                      *refal(:,kk,ka))

                     c9ijk = c9ijk + c9tmp*zvec(n)*zvec(m)*zvec(l)
!                --- intermediates ---
!                    ∂WA/∂CNA·WB·WC
                     x2 = x2 - dzvec(l)*zvec(m)*zvec(n)*c9tmp
!                    WA·∂WB/∂CNB·WC
                     x4 = x4 - dzvec(m)*zvec(l)*zvec(n)*c9tmp
!                    WA·WB·∂WC/∂CNC
                     x6 = x6 - dzvec(n)*zvec(m)*zvec(l)*c9tmp

                  enddo ! refs of k/C
               enddo ! refs of j/B
            enddo ! refs of i/A

!       --- geometrical term and r⁻³AB·r⁻³BC·r⁻³CA ---
!           ATM = 3·cos(α)cos(β)cos(γ)+1
!               = 3/8(r²AB+r²BC-r²CA)(r²AB+r²CA-r²BC)(r²BC+r²CA-r²AB)
!                 /(r²BC·r²CA·r²AB)+1
            atm = ((0.375_wp * (r2ij+r2jk-r2ik) &
            &                * (r2ij+r2ik-r2jk) &
            &                * (r2ik+r2jk-r2ij) / r2ijk ) + 1._wp)/(rijk**3)
            dijatm=-0.375_wp*(r2ij**3+r2ij**2*(r2jk+r2ik) &
            &      +r2ij*(3._wp*r2jk**2+2._wp*r2jk*r2ik+3._wp*r2ik**2) &
            &      -5._wp*(r2jk-r2ik)**2*(r2jk+r2ik)) &
            &      /(r2ijk*rijk**3)*oorij
            djkatm=-0.375_wp*(r2jk**3+r2jk**2*(r2ik+r2ij) &
            &      +r2jk*(3._wp*r2ik**2+2._wp*r2ik*r2ij+3._wp*r2ij**2) &
            &      -5._wp*(r2ik-r2ij)**2*(r2ik+r2ij)) &
            &      /(r2ijk*rijk**3)*oorjk
            dikatm=-0.375_wp*(r2ik**3+r2ik**2*(r2jk+r2ij) &
            &      +r2ik*(3._wp*r2jk**2+2._wp*r2jk*r2ij+3._wp*r2ij**2) &
            &      -5._wp*(r2jk-r2ij)**2*(r2jk+r2ij)) &
            &      /(r2ijk*rijk**3)*oorik

!       --- damping function ---
!           1/(1+6(¾·∛[RAB·RBC·RCA/(rAB·rBC·rCA)])¹⁶)
            fdmp = one/(one+six*(thf*(cijk/rijk)**oth)**par%alp)
            dtmp = -(oth*six*par%alp*(thf*(cijk/rijk)**oth)**par%alp)*fdmp**2
            dijfdmp = dtmp*oorij
            djkfdmp = dtmp*oorjk
            dikfdmp = dtmp*oorik

!       --- intermediates ---
!           ∂WA/∂rAB·WB·WC·DABC = ∂CNA/∂rAB·(∂WA/∂CNA·WB·WC)·DABC
            x1 = x2*dcn(i,j)*( atm*fdmp )
!           ∂WA/∂rCA·WB·WC·DABC = ∂CNA/∂rCA·(∂WA/∂CNA·WB·WC)·DABC
            x2 = x2*dcn(i,k)*( atm*fdmp )
!           WA·∂WB/∂rBC·WC·DABC = ∂CNB/∂rBC·(WA·∂WB/∂rBC·WC)·DABC
            x3 = x4*dcn(j,k)*( atm*fdmp )
!           WA·∂WB/∂rAB·WC·DABC = ∂CNB/∂rAB·(WA·∂WB/∂rAB·WC)·DABC
            x4 = x4*dcn(i,j)*( atm*fdmp )
!           WA·WB·∂WC/∂rCA·DABC = ∂CNC/∂rCA·(WA·WB·∂WC/∂rCA)·DABC
            x5 = x6*dcn(i,k)*( atm*fdmp )
!           WA·WB·∂WC/∂rBC·DABC = ∂CNC/∂rBC·(WA·WB·∂WC/∂rBC)·DABC
            x6 = x6*dcn(j,k)*( atm*fdmp )
!           WA·WB·WC·∂DABC/∂rAB
            x7 = c9ijk*( atm*dijfdmp-dijatm*fdmp )
!           WA·WB·WC·∂DABC/∂rBC
            x8 = c9ijk*( atm*djkfdmp-djkatm*fdmp )
!           WA·WB·WC·∂DABC/∂rCA
            x9 = c9ijk*( atm*dikfdmp-dikatm*fdmp )

!       --- build everything together ---
            eabc = eabc + c9ijk*atm*fdmp

!           ∂rAB/∂A = -∂rAB/∂B
            drij = rij*oorij
!           ∂rBC/∂B = -∂rBC/∂C
            drjk = rjk*oorjk
!           ∂rCA/∂C = -∂rCA/∂A
            drik = rik*oorik

!           ∂EABC/∂A =
!           + (∂WA/∂rAB-∂WA/∂rCA)·WB·WC·DABC
!           + WA·∂WB/∂rAB·WC·DABC
!           - WA·WB·∂WC/∂rCA·DABC
!           + WA·WB·WC·(∂DABC/∂rAB-∂DABC/∂rCA)
            g(:,i) = g(:,i) + ( &
            &        + (x1+x4+x7)*drij &
            &        - (x2+x5+x9)*drik )
!           ∂EABC/∂B =
!           - ∂WA/∂rAB·WB·WC·DABC
!           + WA·(∂WB/∂rBC-∂WB/∂rAB)·WC·DABC
!           + WA·WB·∂WC/∂rBC·DABC
!           + WA·WB·WC·(∂DABC/∂rBC-∂DABC/∂rAB)
            g(:,j) = g(:,j) + ( &
            &        - (x1+x4+x7)*drij &
            &        + (x3+x6+x8)*drjk )
!           ∂EABC/∂C =
!           + ∂WA/∂rCA·WB·WC·DABC
!           - WA·∂WB/∂rBC·WC·DABC
!           + WA·WB·(∂WC/∂rCA-∂WC/∂rBC)·DABC
!           + WA·WB·WC·(∂DABC/∂rCA-∂DABC/∂rBC)
            g(:,k) = g(:,k) + ( &
            &        + (x2+x5+x9)*drik &
            &        - (x3+x6+x8)*drjk )

         enddo ! k/C
      enddo ! j/B
   enddo ! i/A
!$omp enddo
!$omp endparallel

   if (present(eout)) eout=eabc

end subroutine dabcgrad


subroutine d4_numgrad(nat,ndim,at,q,xyz, &
           &          par,wf,g_a,g_c,mbd,numg)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: q(nat) 
   real(wp),intent(in)  :: xyz(3,nat) 
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: wf,g_a,g_c
   integer, intent(in)  :: mbd
   real(wp),intent(inout) :: numg(3,nat)

   integer  :: k,l
   real(wp) :: er,el,step,cn_thr
   real(wp),allocatable :: xyzl(:,:)
   real(wp),allocatable :: xyzr(:,:)
   real(wp),allocatable :: xyzdup(:,:)
   real(wp),allocatable :: c6abns(:,:)
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: gw(:)
   parameter(step=1e-4_wp)
   parameter(cn_thr=1600._wp)

   allocate( xyzr(3,nat),xyzl(3,nat), source=xyz )
   allocate( c6abns(ndim,ndim),cn(nat),gw(ndim), &
   &         source = 0._wp )

   do l = 1, 3
      do k = 1, nat
         er = 0._wp
         el = 0._wp
         cn=0._wp
         gw=0._wp
         c6abns=0._wp
         xyzr = xyz
         xyzr(l,k) = xyz(l,k) + step
         call covncoord(nat,at,xyzr,cn,cn_thr)
         call d4(nat,ndim,at,wf,g_a,g_c,cn,gw,c6abns)
         call edisp(nat,ndim,at,q,xyzr,par,g_a,g_c, &
         &          gw,c6abns,mbd,er)
         cn=0._wp
         gw=0._wp
         c6abns=0._wp
         xyzl = xyz
         xyzl(l,k) = xyz(l,k) - step
         call covncoord(nat,at,xyzl,cn,cn_thr)
         call d4(nat,ndim,at,wf,g_a,g_c,cn,gw,c6abns)
         call edisp(nat,ndim,at,q,xyzl,par,g_a,g_c, &
         &          gw,c6abns,mbd,el)
         numg(l,k) = numg(l,k) + (er-el)/(2*step)
      enddo
   enddo

end subroutine d4_numgrad

subroutine dispmb(E,aw,xyz,oor6ab,nat)
   use iso_fortran_env, only : wp => real64
   use mctc_la, only : syev,gemm
   implicit none
   integer, intent(in)  :: nat
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: aw(23,nat)
   real(wp),intent(in)  :: oor6ab(nat,nat)
   real(wp),intent(out) :: E

   integer  :: i,j,ii,jj,k
   integer  :: info
   real(wp) :: tau(3,3),spur(23),d_,d2,r(3),r2,alpha
   real(wp) :: two(23),atm(23),d3
   real(wp),allocatable :: T (:,:)
   real(wp),allocatable :: A (:,:)
   real(wp),allocatable :: AT(:,:)
   real(wp),allocatable :: F (:,:)
   real(wp),allocatable :: F_(:,:)
   real(wp),allocatable :: d (:)
   real(wp),allocatable :: w (:)

   intrinsic :: sum,sqrt,minval,log

   allocate( T(3*nat,3*nat),  A(3*nat,3*nat), AT(3*nat,3*nat), &
   &         F(3*nat,3*nat), F_(3*nat,3*nat),  d(3*nat), &
   &         w(12*nat), &
   &         source = 0.0_wp )

   spur = 0.0_wp

   do i = 1, 3*nat
      F(i,i) = 1.0_wp
   enddo

   do i = 1, nat
      do j  = 1, i-1
         r  = xyz(:,j) - xyz(:,i)
         r2 = sum(r**2)
         do ii = 1, 3
            tau(ii,ii) = (3*r(ii)*r(ii)-r2)/r2
            do jj = ii+1, 3
               tau(ii,jj) = (3*r(ii)*r(jj))/r2
               tau(jj,ii) = tau(ii,jj)
            enddo
         enddo
         tau = tau*sqrt(oor6ab(i,j))
         T(3*i-2:3*i,3*j-2:3*j) = tau
         T(3*j-2:3*j,3*i-2:3*i) = tau
      enddo
   enddo

   !call prmat(6,T,3*nat,3*nat,'T')

   do k = 1, 23
      A = 0.0_wp
      do i =  1, nat
         alpha = sqrt(aw(k,i))
         A(3*i-2,3*i-2) = alpha
         A(3*i-1,3*i-1) = alpha
         A(3*i  ,3*i  ) = alpha
      enddo

      AT  = 0.0d0 
      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,A,3*nat,T, &
  &             3*nat,0.0_wp,F_,3*nat)
      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,F_,3*nat,A, &
  &             3*nat,0.0_wp,AT,3*nat)

      F_ = F - AT

      d = 0.0d0
      call syev('N','U',3*nat,F_,3*nat,d,w,12*nat,info)
      if (info.ne.0) then
!        call raise('W','MBD eigenvalue not solvable')
         print'(1x,''* MBD eigenvalue not solvable'')'
         E = 0.0_wp
         return
      endif
      if (minval(d).le.0.0d0) then
!        call raise('W','Negative MBD eigenvalue occurred')
         print'(1x,''* Negative MBD eigenvalue occurred'')'
         E = 0.0_wp
         return
      endif

      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,AT,3*nat,AT, &
  &             3*nat,0.0_wp,F_,3*nat)
!     call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,F_,3*nat,AT, &
! &             3*nat,0.0_wp,A,3*nat)
       
      d_ = 1.0_wp; d2 = 0.0_wp!; d3 = 0.0_wp
      do i = 1, 3*nat
         d_ = d_ * d(i)
         d2 = d2 - F_(i,i)
!        d3 = d3 - A(i,i)
      enddo
      spur(k) = log(d_) - d2*0.5
!     two(k) = d2/2.0_wp
!     atm(k) = d3/3.0_wp
   enddo

   E = trapzd(spur)*ooTPI
   !print*,'     full contribution', trapzd(spur)*ooTPI
   !print*,' manybody contribution', trapzd(spur-two)*ooTPI
   !print*,'  twobody contribution', trapzd(two)*ootpi
   !print*,'threebody contribution', trapzd(atm)*ootpi

   deallocate(T,A,AT,F,F_,d)
end subroutine dispmb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING! !!
!! This implementation of the MBD gradient is incorrect, DO NOT USE! !!
!! !WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING! !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mbdgrad(nat,xyz,aw,daw,oor6ab,g,E)
   use iso_fortran_env, only : wp => real64
   use mctc_la, only : syev,gemm
   implicit none
   integer, intent(in)  :: nat
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: aw(23,nat)
   real(wp),intent(in)  :: daw(23,3,nat)
   real(wp),intent(in)  :: oor6ab(nat,nat)
   real(wp),intent(inout)        :: g(3,nat)
   real(wp),intent(out),optional :: E

   integer  :: i,ii,j,jj,k,kk,l,ll,m,n
   integer  :: info
   real(wp) :: spur(23)
   real(wp) :: fdmp,dfdmp
   real(wp) :: dum,dum2,ddum(3)
   real(wp) :: r(3),r2,drij(3)
   real(wp) :: alpha
   real(wp) :: tau(3,3),dtau(3,3,3)
   real(wp),allocatable :: dspur(:,:,:)
   real(wp),allocatable :: dA(:,:)
   real(wp),allocatable :: dT(:,:)
   real(wp),allocatable :: dAT(:,:)
   real(wp),allocatable :: nT(:,:,:) ! nabla T
   real(wp),allocatable :: T(:,:)
   real(wp),allocatable :: A(:,:)
   real(wp),allocatable :: AT(:,:)
   real(wp),allocatable :: dF(:,:)
   real(wp),allocatable :: F(:,:)
   real(wp),allocatable :: tmp1(:,:)
   real(wp),allocatable :: tmp2(:,:)
   real(wp),allocatable :: d(:)
   real(wp),allocatable :: w(:)

   intrinsic :: sum,sqrt,minval

   ! the MBD calculation in this RPA-like scheme needs quite a lot
   ! of memory, which might be problematic for big systems.
   ! We could use one or maybe to matrices less which would
   ! obfuscate the code
   allocate( T(3*nat,3*nat),A(3*nat,3*nat),AT(3*nat,3*nat), &
   &         F(3*nat,3*nat),dF(3*nat,3*nat), &
   &         tmp1(3*nat,3*nat),tmp2(3*nat,3*nat), &
   &         dT(3*nat,3*nat),dA(3*nat,3*nat), &
   &         dAT(3*nat,3*nat),nT(3*nat,3*nat,3), &
   &         d(3*nat),w(12*nat),dspur(23,3,nat), &
   &         source = 0.0_wp )

   spur = 0.0_wp

   do i = 1, 3*nat
      F(i,i) = 1.0_wp
   enddo

!-----------------------------------------------------------------------
!  interaction tensor setup                                    SAW 1708
!-----------------------------------------------------------------------
   do i = 1, nat
      do j  = 1, i-1
         r  = xyz(:,j) - xyz(:,i)
         r2 = sum(r**2)
         do ii = 1, 3
            tau(ii,ii) = (3*r(ii)**2-r2)/r2
            do jj = ii+1, 3
               tau(ii,jj) = (3*r(ii)*r(jj))/r2
               tau(jj,ii) = tau(ii,jj)
            enddo
         enddo

         fdmp = sqrt(oor6ab(i,j))
         !print*, fdmp
         tau = tau*fdmp
         T(3*i-2:3*i,3*j-2:3*j) = tau
         T(3*j-2:3*j,3*i-2:3*i) = tau
         
!-----------------------------------------------------------------------
!        derivative of interaction tensor                      SAW 1802
!-----------------------------------------------------------------------
         dfdmp = -3*oor6ab(i,j)**2/fdmp*r2**2 ! *sqrt(r2)

         do ii = 1, 3
            dtau(ii,ii,ii) = (15*r(ii)*(0.6_wp*r2-r(ii)**2))/r2**2
            do jj = ii+1, 3
               dtau(jj,ii,ii) = (15*r(jj)*(0.2_wp*r2-r(ii)**2))/r2**2
               dtau(ii,jj,ii) = dtau(jj,ii,ii)
               dtau(ii,ii,jj) = dtau(jj,ii,ii)
               dtau(jj,jj,ii) = (15*r(ii)*(0.2_wp*r2-r(jj)**2))/r2**2
               dtau(jj,ii,jj) = dtau(jj,jj,ii)
               dtau(ii,jj,jj) = dtau(jj,jj,ii)
               do kk = jj+1, 3
                  dtau(ii,jj,kk) = -(15*r(ii)*r(jj)*r(kk))/r2**2
                  dtau(jj,kk,ii) = dtau(ii,jj,kk)
                  dtau(kk,ii,jj) = dtau(ii,jj,kk)
                  dtau(kk,jj,ii) = dtau(ii,jj,kk)
                  dtau(ii,kk,jj) = dtau(ii,jj,kk)
                  dtau(jj,ii,kk) = dtau(ii,jj,kk)
               enddo
            enddo
         enddo

         !drij = r/sqrt(r2)

         dtau(:,:,1) = ( dtau(:,:,1)*fdmp + tau*dfdmp*r(1) )
         dtau(:,:,2) = ( dtau(:,:,2)*fdmp + tau*dfdmp*r(2) )
         dtau(:,:,3) = ( dtau(:,:,3)*fdmp + tau*dfdmp*r(3) )

         !print*,dtau
         nT(3*i-2:3*i,3*j-2:3*j,1:3) = dtau
         nT(3*j-2:3*j,3*i-2:3*i,1:3) = dtau

      enddo
   enddo
        !call prmat(6,T,3*nat,3*nat,'T')

!-----------------------------------------------------------------------
! RPA-like calculation of MBD energy                           SAW 1708
!-----------------------------------------------------------------------
! EMBD = 1/(2π)∫dω Tr{log[1-AT]} = 1/(2π)∫dω log[∏(i)Λii]
! E(2) = 1/(2π)∫dω ½ Tr{(AT)²} ! this two-body energy has to be removed
!-----------------------------------------------------------------------
   do k = 1, 23
      A = 0.0_wp
      do i =  1, nat
         alpha = sqrt(aw(k,i))
         A(3*i-2,3*i-2) = alpha
         A(3*i-1,3*i-1) = alpha
         A(3*i  ,3*i  ) = alpha
      enddo

      AT  = 0.0_wp
      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,A,3*nat,T, &
           &     3*nat,0.0_wp,tmp1,3*nat)
      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp1,3*nat,A, &
           &     3*nat,0.0_wp,AT,3*nat)
        !call prmat(6,AT,3*nat,3*nat,'AT')

      tmp1 = F - AT
        !call prmat(6,dF,3*nat,3*nat,'1-AT')

      d = 0.0d0
      call syev('V','U',3*nat,tmp1,3*nat,d,w,12*nat,info)
        !call prmat(6,tmp1,3*nat,3*nat,'F eigv.')
      if (info.ne.0) then
!        call raise('W','MBD eigenvalue not solvable')
         print'(1x,''* MBD eigenvalue not solvable'')'
         E = 0.0_wp
         return
      endif
      if (minval(d).le.0.0d0) then
!        call raise('W','Negative MBD eigenvalue occurred')
         print'(1x,''* Negative MBD eigenvalue occurred'')'
         E = 0.0_wp
         return
      endif

!     two-body contribution to energy
      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,AT,3*nat,AT, &
           &     3*nat,0.0_wp,tmp2,3*nat)
        !call prmat(6,dF,3*nat,3*nat,'dF')

       
!     get the trace
      dum = 1.0_wp; dum2 = 0.0_wp
      do i = 1, 3*nat
         dum = dum * d(i)
         dum2 = dum2 + tmp2(i,i)
      enddo
      spur(k) = log(dum) + 0.5_wp * dum2
      !print*, log(dum), 0.5_wp*dum2,spur(k)

!-----------------------------------------------------------------------
! MBD gradient calculation                                     SAW 1803
!-----------------------------------------------------------------------
! Some theory:
! EMBD = 1/(2π)∫dω Tr{log[1-AT]} = 1/(2π)∫dω log[∏(i)Λii]
! E(2) = 1/(2π)∫dω ½ Tr{(AT)²}
! ∇EMBD = 1/(2π)∫dω Tr{∇(log[1-AT])} = 1/(2π)∫dω Tr{(1-AT)⁻¹·(∇AT+A∇T)}
!       = 1/(2π)∫dω Tr{(1-A^½TA^½)⁻¹·((∇A^½)TA^½+A^½(∇T)A^½+A^½T(∇A^½))}
! (1-AT)⁻¹ = U·Λ⁻¹·U†  w/ U⁻¹ = U†  if (1-AT) is symmetrized
! Aijkl = δij·δkl·∑(ref)Wi·α
! ∂/∂X Aijkl = δij·δkl·∑(ref)∂Wi/∂CNi·∂CNi/∂rij·∂rij/∂X·α
!-----------------------------------------------------------------------
! eigenvectors are still saved on tmp1, eigenvalues are still on d
! we (over)use tmp2 to hold intermediate results
!-----------------------------------------------------------------------

!     get the inverse of 1-AT by using its eigenvalues
      dF = 0.0_wp
      do i = 1, 3*nat
         dF(i,i) = 1._wp/d(i)
      enddo
        !call prmat(6,dF,3*nat,3*nat,'dF')
!     (1-AT)⁻¹ = U·Λ⁻¹·U†
      call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp1,3*nat,dF, &
           &     3*nat,0.0_wp,tmp2,3*nat)
!     call gemm('N','T',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,tmp1, &
!          &     3*nat,0.0_wp,AT,3*nat)
        !call prmat(6,dF,3*nat,3*nat,'dF')
!     unfortunately this might not be enough, we have to substract the
!     twobody gradient for the dipole-dipole interaction, which is in
!     fact not easily accessable from here.
!     If this is correct,
!        E(2) = 1/(2π)∫dω ½ Tr{(AT)²}
!       ∇E(2) = 1/(2π)∫dω Tr{(AT)·(∇AT+A∇T)}
!     then the MBD gradient w/o two-body contrib. could be represented by
!     dF = dF - AT
!     or more easier by the use of dgemm's beta by replacing the last
!     dgemm by:
!     call gemm('N','T',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,tmp1, &
!          &     3*nat,1.0_wp,AT,3*nat)
!     But we would have to consider AT instead of dF in the following code
!     which would be at least confusing

!-----------------------------------------------------------------------
!     Efficient MBD-Gradient calculation, O(N³)                SAW 1803
!-----------------------------------------------------------------------
      do l = 1, 3
         dA = 0.0_wp
         do i = 1, nat
            dA(3*i-2,3*i-2) = 0.5*daw(k,l,i)/sqrt(aw(k,i))
            dA(3*i-1,3*i-1) = 0.5*daw(k,l,i)/sqrt(aw(k,i))
            dA(3*i  ,3*i  ) = 0.5*daw(k,l,i)/sqrt(aw(k,i))
         enddo

!        (∇A^½)TA^½
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,dA,3*nat,T, &
              &     3*nat,0.0_wp,tmp2,3*nat)
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,A, &
              &     3*nat,0.0_wp,dAT,3*nat)

!        (∇A^½)TA^½+A^½(∇T)A^½ (please note the use of beta=1.0!)
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,A,3*nat,nT(:,:,l), &
              &     3*nat,0.0_wp,tmp2,3*nat)
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,A, &
              &     3*nat,1.0_wp,dAT,3*nat)

!        last term and we have: (∇A^½)TA^½+A^½(∇T)A^½+A^½T(∇A^½)
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,A,3*nat,T, &
              &     3*nat,0.0_wp,tmp2,3*nat)
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,tmp2,3*nat,dA, &
              &     3*nat,1.0_wp,dAT,3*nat)

!        by multiplying the prefactor we get to:
!        ((1-A^½TA^½)⁻¹+A^½TA^½)·((∇A^½)TA^½+A^½(∇T)A^½+A^½T(∇A^½))
         call gemm('N','N',3*nat,3*nat,3*nat,1.0_wp,AT,3*nat,dAT, &
              &     3*nat,0.0_wp,tmp2,3*nat)
!        now the other term of
!        [(1-A^½TA^½)⁻¹+A^½TA^½),((∇A^½)TA^½+A^½(∇T)A^½+A^½T(∇A^½)]
         call gemm('N','N',3*nat,3*nat,3*nat,-1.0_wp,dAT,3*nat,AT, &
              &     3*nat,1.0_wp,tmp2,3*nat)

         do i = 1, nat
            dspur(k,l,i) = dspur(k,l,i) &
            &              + tmp2(3*i-2,3*i-2) &
            &              + tmp2(3*i-1,3*i-1) &
            &              + tmp2(3*i  ,3*i  )
         enddo ! all atoms
         
      enddo ! cartesian parts

   enddo ! k, imaginary frequencies

   E = trapzd(spur)*ootpi
   !print*, E
   do i = 1, nat
      do j = 1, 3
         g(j,i) = g(j,i) + ootpi*trapzd(dspur(:,j,i))
         !print*,ootpi*trapzd(dspur(:,j,i))
      enddo
   enddo

   deallocate( T,A,AT,F,tmp1,tmp2,dT,dA,dAT,dF,nT,d,w,dspur )

end subroutine mbdgrad

function refq2string(refq) result(string)
   implicit none
   integer,intent(in) :: refq
   character(len=:),allocatable :: string
   select case(refq)
   case default;                  string = 'unknown'
   case(p_refq_gfn2xtb);          string = 'GFN2'
   case(p_refq_gasteiger);        string = 'EEQ'
   case(p_refq_hirshfeld);        string = 'extern'
   case(p_refq_periodic);         string = 'EEQ'
   case(p_refq_gfn2xtb_gbsa_h2o); string = 'GFN2/GBSA'
   case(p_refq_goedecker);        string = 'EEQ'
   end select
end function refq2string

function lmbd2string(lmbd) result(string)
   implicit none
   integer,intent(in) :: lmbd
   character(len=:),allocatable :: string
   select case(lmbd)
   case default;           string = 'unknown'
   case(p_mbd_none);       string = 'none'
   case(p_mbd_rpalike);    string = 'RPA like'
   case(p_mbd_exact_atm);  string = 'ATM'
   case(p_mbd_approx_atm); string = 'ATM'
   end select
end function lmbd2string

! --- PBC
subroutine pbc_d4(nat,ndim,at,wf,g_a,g_c,covcn,gw,refc6)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: ndim
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: wf,g_a,g_c
   real(wp),intent(in)  :: covcn(nat)
   real(wp),intent(out) :: gw(ndim)
   real(wp),intent(out) :: refc6(ndim,ndim)

   integer  :: i,ia,is,icn,ii,iii,j,jj,ja,k,l
   integer,allocatable :: itbl(:,:)
   real(wp) :: twf,norm,aiw(23)

   intrinsic :: maxval

   allocate( itbl(7,nat), source = 0 )

   gw = 0._wp
   refc6 = 0._wp

   k = 0
   do i = 1, nat
      do ii = 1, refn(at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

   do i = 1, nat
      ia = at(i)
      norm = 0.0_wp
      do ii = 1, refn(ia)
         do iii = 1, refc(ii,ia)
            twf = iii*wf
            norm = norm + cngw(twf,covcn(i),refcovcn(ii,ia))
         enddo
      enddo
      norm = 1._wp / norm
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         do iii = 1, refc(ii,ia)
            twf = iii*wf
            gw(k) = gw(k) + cngw(twf,covcn(i),refcovcn(ii,ia)) * norm
         enddo
!    --- okay, if we run out of numerical iso_fortran_env, gw(k) will be NaN.
!        In case it is NaN, it will not match itself! So we can rescue
!        this exception. This can only happen for very high CNs.
         if (gw(k).ne.gw(k)) then
            if (maxval(refcovcn(:refn(ia),ia)).eq.refcovcn(ii,ia)) then
               gw(k) = 1.0_wp
            else
               gw(k) = 0.0_wp
            endif
         endif
         ! diagonal terms
         do jj = 1, ii
            l = itbl(jj,i)
            aiw = refal(:,ii,ia)*refal(:,jj,ia)
            refc6(l,k) = thopi * trapzd(aiw)
            refc6(k,l) = refc6(l,k)
         enddo
         ! offdiagonal terms
         do j = 1, i-1
            ja = at(j)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               aiw = refal(:,ii,ia)*refal(:,jj,ja)
               refc6(l,k) = thopi * trapzd(aiw)
               refc6(k,l) = refc6(l,k)
            enddo
         enddo
      enddo
   enddo

end subroutine pbc_d4

! compute D4 gradient under pbc
subroutine edisp_3d(mol,ndim,q,rep,atm_rep,r_thr,atm_thr,par,g_a,g_c,gw,refc6,mbd, &
      &             ed,etwo,embd)
   use iso_fortran_env, only : wp => real64
   use mctc_constants
   use tbdef_molecule
   implicit none
   type(tb_molecule),intent(in) :: mol
   integer, intent(in)  :: ndim
   real(wp),intent(in)  :: q(mol%n)
   integer, intent(in)  :: rep(3)
   integer, intent(in)  :: atm_rep(3)
   real(wp),intent(in)  :: r_thr
   real(wp),intent(in)  :: atm_thr
   integer              :: tx,ty,tz
   real(wp)             :: t(3)
   real(wp)             :: aiw(23)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: g_a,g_c
   real(wp),intent(in)  :: gw(ndim)
   real(wp),intent(in)  :: refc6(ndim,ndim)
   integer, intent(in)  :: mbd
   real(wp),intent(out) :: ed
   real(wp),intent(out),optional :: etwo
   real(wp),intent(out),optional :: embd

   integer  :: i,ii,iii,j,jj,k,l,ia,ja,ij
   integer, allocatable :: itbl(:,:)
   real(wp) :: iz
   real(wp) :: qmod,eabc
   real(wp) :: norm,dnorm
   real(wp) :: dexpw,expw
   real(wp) :: twf,tgw,r4r2ij
   real(wp) :: rij(3),r,r2,r4,r6,r8,R0
   real(wp) :: oor6,oor8,oor10,door6,door8,door10,cutoff
   real(wp) :: c8abns,disp,ddisp,x1,x2,x3
   real(wp) :: c6ii,c6ij,dic6ii,dic6ij,djc6ij,dizii,dizij,djzij
   real(wp) :: rcovij,expterm,den,dcndr

   real(wp) :: drdx(3),dtmp,gwk,dgwk
   real(wp),allocatable :: r2ab(:)
   real(wp),allocatable :: dc6dcn(:)
   real(wp),allocatable :: zetavec(:)
   real(wp),allocatable :: zerovec(:)
   real(wp) :: cn_thr,gw_thr

   parameter(cn_thr = 1600.0_wp)
   parameter(gw_thr=0.000001_wp)

   intrinsic :: present,sqrt,sum,maxval,exp,abs

   !  print'(" * Allocating local memory")'
   allocate( zetavec(ndim),zerovec(ndim), source = 0.0_wp )
   allocate( itbl(7,mol%n), source = 0 )

   ed = 0.0_wp
   eabc = 0.0_wp

   k = 0
   do i = 1, mol%n
      do ii = 1, refn(mol%at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

!  print'(" * Entering first OMP section")'
!$OMP parallel default(none) &
!$omp private(i,ii,iii,ia,iz,k)  &
!$omp shared (mol,refn,refc,refcovcn,itbl,refq,g_a,g_c,q) &
!$omp shared (gw,zetavec,zerovec,r_thr)
!$omp do
   do i = 1, mol%n
      ia = mol%at(i)
      iz = zeff(ia)
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         zetavec(k) = zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz) * gw(k)
         ! NEW: q=0 for ATM
         zerovec(k) =  zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,iz) * gw(k)

      enddo
   enddo
!$omp end do
!$omp end parallel

!$OMP parallel default(none) &
!$omp private(i,j,ia,ja,ij,k,l,c6ii,c6ij,disp,  &
!$omp         rij,r2,r,r4r2ij,r0,oor6,oor8,oor10, &
!$omp         t,tx,ty,tz)  &
!$omp shared(mol,refn,itbl,zetavec,refc6,par,rep,r_Thr) &
!$omp shared(r2ab) reduction(+:ed)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      ia = mol%at(i)
      ! temps
      c6ij   = 0.0_wp
      ! all refs
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         do jj = 1, refn(ia)
            l = itbl(jj,i)
            c6ij   = c6ij   +  zetavec(k) *  zetavec(l) * refc6(k,l)
         enddo
      enddo
      ! i in primitive cell with i in images
      r4r2ij = 3*r4r2(ia)*r4r2(ia)
      r0 = par%a1*sqrt(r4r2ij) + par%a2
      do concurrent(tx = -rep(1):rep(1), &
            &       ty = -rep(2):rep(2), &
            &       tz = -rep(3):rep(3))
         ! cycle i with i interaction in same cell
         if (tx.eq.0.and.ty.eq.0.and.tz.eq.0) cycle
         rij = tx*mol%lattice(:,1) + ty*mol%lattice(:,2) + tz*mol%lattice(:,3)
         r2  = sum( rij**2 )
         if (r2.gt.r_thr) cycle
         r   = sqrt(r2)
         oor6 = 1._wp/(r2**3+r0**6)
         oor8 = 1._wp/(r2**4+r0**8)
         oor10 = 1._wp/(r2**5+r0**10)
         disp = par%s6*oor6 + par%s8*r4r2ij*oor8   &
            & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
         ed = ed - c6ij*disp
      enddo ! tx
      ! over all j atoms
      do j = 1, i-1
         ja = mol%at(j)
         ! temps
         c6ij   = 0.0_wp
         ! all refs
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c6ij   = c6ij   +  zetavec(k) *  zetavec(l) * refc6(k,l)
            enddo
         enddo
         r4r2ij = 3*r4r2(ia)*r4r2(ja)
         r0 = par%a1*sqrt(r4r2ij) + par%a2
         do concurrent(tx = -rep(1):rep(1), &
               &       ty = -rep(2):rep(2), &
               &       tz = -rep(3):rep(3))
            t = tx*mol%lattice(:,1) + ty*mol%lattice(:,2) + tz*mol%lattice(:,3)
            rij = mol%xyz(:,i) - mol%xyz(:,j) + t
            r2 = sum(rij**2)
            if (r2.gt.r_thr) cycle
            r = sqrt(r2)
            oor6 = 1._wp/(r2**3+r0**6)
            oor8 = 1._wp/(r2**4+r0**8)
            oor10 = 1._wp/(r2**5+r0**10)
            disp = par%s6*oor6 + par%s8*r4r2ij*oor8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
            ed = ed - c6ij*disp
         enddo ! tx
      enddo
   enddo
!$omp enddo
!$omp end parallel

   eabc = 0.0_wp
   if (mbd > 0) then
      call abcappr_3d(mol,ndim,par,zerovec,refc6,itbl,atm_rep,atm_thr,eabc)
   endif

   !  print'(" * Dispersion all done, saving variables")'
   if (present(etwo)) etwo = ed
   if (present(embd)) embd = eabc
   ed = ed + eabc

end subroutine edisp_3d

!> @brief calculates threebody dispersion gradient from C6 coefficients
subroutine abcappr_3d(mol,ndim,par,zvec,refc6,itbl,rep,r_thr,eabc)
   use tbdef_molecule
   use pbc_tools
   implicit none
   type(tb_molecule),intent(in) :: mol !< molecular structure information
   integer, intent(in)  :: ndim
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: zvec(ndim)
   real(wp),intent(in)  :: refc6(ndim,ndim)
   integer, intent(in)  :: itbl(7,mol%n)
   integer, intent(in)  :: rep(3)
   real(wp),intent(in)  :: r_thr
   real(wp),intent(out) :: eabc

   integer  :: i,ii,ia,j,jj,ja,k,l
   integer  :: ij
   real(wp),allocatable :: c6ab(:)
   real(wp) :: r,r2ij
   real(wp) :: c6ij
   real(wp) :: gw_thr
   parameter(gw_thr=0.0001_wp)

   intrinsic :: present,sqrt

   allocate( c6ab(mol%n*(mol%n+1)/2), source = 0.0_wp )

   eabc = 0.0_wp

   !$OMP parallel default(none) &
   !$omp private(i,ia,j,ja,ij,r2ij,c6ij,k,l,r)  &
   !$omp shared (mol,refn,itbl,refc6,zvec) &
   !$omp shared (c6ab)
   !$omp do schedule(dynamic)
   do i = 1, mol%n
      ia = mol%at(i)
      ij = i*(i-1)/2+i

      ! temps
      c6ij = 0.0_wp
      ! all refs
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         do jj = 1, refn(ia)
            l = itbl(jj,i)
            c6ij = c6ij + zvec(k)*zvec(l)*refc6(k,l)
         enddo
      enddo
      ! save
      c6ab(ij) = c6ij

      do j = 1, i-1
         !        if(i.eq.j) cycle
         ja = mol%at(j)
         ij = i*(i-1)/2 + j

         !        first check if we want this contribution
         !if(mol%dist(j,i)**2.gt.r_thr) cycle

         ! temps
         c6ij = 0.0_wp
         ! all refs
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*refc6(k,l)
            enddo
         enddo
         ! save
         c6ab(ij) = c6ij
      enddo
   enddo
   !$omp enddo
   !$omp end parallel

   call abcappr_3d_dftd3_like_style(mol%n,mol%at,mol%xyz,par,r_thr,rep, &
      &                             mol%lattice,c6ab,eabc)
   !call abcappr_3d_wsc(mol,par,[2,2,2],r_thr,c6ab,eabc)
   !call abcappr_3d_bvk(mol,par,rep,r_thr,c6ab,eabc)

end subroutine abcappr_3d

subroutine abcappr_3d_dftd3_like_style(nat,at,xyz,par,thr,rep,dlat,c6ab,eabc)
   implicit none
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in) :: thr
   integer, intent(in) :: rep(3)
   real(wp),intent(in) :: c6ab(nat*(nat+1)/2)
   real(wp),intent(in) :: dlat(3,3)

   real(wp),intent(inout) :: eabc

   real(wp),parameter :: six = 6.0_wp, oth = 1.0_wp/3.0_wp

   integer iat,jat,kat
   real(wp) x1
   real(wp) r2,r
   real(wp) fdmp,tmp1

   real(wp) rij(3),rik(3),rjk(3)
   real(wp), allocatable,dimension(:,:,:,:) ::  dc6dr  !d(E)/d(r_ij) derivative wrt. dist. iat-jat
   !dCN(jat)/d(r_ij)
   real(wp) :: r9ijk
   real(wp) vec(3)
   integer ij,ik,jk

   real(wp),dimension(3) ::ijvec,ikvec,jkvec,t,s,dumvec
   integer tx,ty,tz,sx,sy,sz
   real(wp) rij2,rik2,rjk2,c9,c6ij,c6ik,c6jk,rijk,rijk3
   real(wp) :: cij,cjk,cik,cijk
   real(wp) time1,time2,rijk2,dc9,dfdmp,dang,ang
   integer,dimension(3) :: repmin,repmax

   allocate(dc6dr(-rep(3):rep(3),-rep(2):rep(2), &
      &           -rep(1):rep(1),nat*(nat+1)/2))
   dc6dr = 0.0_wp


   !        write(*,*)'!!!!!!!!!!    THREEBODY  GRADIENT  !!!!!!!!!!'
   eabc=0.0_wp
   !        write(*,*)'thr:',sqrt(thr)

   do iat=3,nat
      do jat=2,iat-1
         ij=iat*(iat-1)/2+jat
         ijvec=xyz(:,jat)-xyz(:,iat)

         c6ij=c6ab(ij)
         do kat=1,jat-1
            ik=iat*(iat-1)/2+kat
            jk=jat*(jat-1)/2+kat
            ikvec=xyz(:,kat)-xyz(:,iat)
            jkvec=xyz(:,kat)-xyz(:,jat)

            c6ik=c6ab(ik)
            c6jk=c6ab(jk)
            c9=-sqrt(c6ij*c6ik*c6jk)
            cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
            cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
            cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
            cijk = cij*cjk*cik

            do concurrent (tx=-rep(1):rep(1), &
                  &        ty=-rep(2):rep(2), &
                  &        tz=-rep(3):rep(3))
               repmin(1)=max(-rep(1),tx-rep(1))
               repmax(1)=min(+rep(1),tx+rep(1))
               repmin(2)=max(-rep(2),ty-rep(2))
               repmax(2)=min(+rep(2),ty+rep(2))
               repmin(3)=max(-rep(3),tz-rep(3))
               repmax(3)=min(+rep(3),tz+rep(3))
               t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
               rij2=SUM((ijvec+t)*(ijvec+t))
               if(rij2.gt.thr)cycle

               !rr0ij=sqrt(rij2)/r0ab(at(iat),at(jat))


               do concurrent (sx=repmin(1):repmax(1), &
                     &        sy=repmin(2):repmax(2), &
                     &        sz=repmin(3):repmax(3))
                  s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
                  rik2=SUM((ikvec+s)*(ikvec+s))
                  if(rik2.gt.thr)cycle

                  dumvec=jkvec+s-t
                  rjk2=SUM(dumvec*dumvec)
                  if(rjk2.gt.thr)cycle
                  !rr0ik=sqrt(rik2)/r0ab(at(iat),at(kat))
                  !rr0jk=sqrt(rjk2)/r0ab(at(jat),at(kat))
                  rijk2=(rij2*rjk2*rik2)
                  ! first calculate the three components for the energy calculation fdmp
                  ! and ang
                  !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
                  !damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

                  rijk=sqrt(rijk2)
                  fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
                  rijk3=rijk*rijk2
                  ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                     *(-rij2+rjk2+rik2)/(rijk3*rijk2) &
                     +1.0_wp/(rijk3)

                  r9ijk=ang*fdmp
                  eabc=eabc-r9ijk*c9

               enddo !sz
            enddo !tx
         enddo !kat
      enddo !jat
   enddo !iat

   ! Now the interaction with jat=iat of the triples iat,iat,kat
   DO iat=2,nat
      jat=iat
      ij=iat*(iat-1)/2+jat
      ijvec=0.0_wp

      c6ij=c6ab(ij)
      DO kat=1,iat-1
         jk=jat*(jat-1)/2+kat
         ik=jk

         c6ik=c6ab(ik)
         c6jk=c6ik
         ikvec=xyz(:,kat)-xyz(:,iat)
         jkvec=ikvec
         c9=-sqrt(c6ij*c6ik*c6jk)
         cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
         cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
         cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
         cijk = cij*cjk*cik
         do concurrent (tx=-rep(1):rep(1), &
               &        ty=-rep(2):rep(2), &
               &        tz=-rep(3):rep(3))
            repmin(1)=max(-rep(1),tx-rep(1))
            repmax(1)=min(+rep(1),tx+rep(1))
            repmin(2)=max(-rep(2),ty-rep(2))
            repmax(2)=min(+rep(2),ty+rep(2))
            repmin(3)=max(-rep(3),tz-rep(3))
            repmax(3)=min(+rep(3),tz+rep(3))
            IF (tx.eq.0 .and. ty.eq.0 .and. tz.eq.0) cycle
            t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
            dumvec=t
            rij2=SUM(dumvec*dumvec)
            if(rij2.gt.thr)cycle

            !rr0ij=sqrt(rij2)/r0ab(at(iat),at(jat))

            do concurrent (sx=repmin(1):repmax(1), &
                  &        sy=repmin(2):repmax(2), &
                  &        sz=repmin(3):repmax(3))
               ! every result * 0.5

               s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
               dumvec=ikvec+s
               dumvec=dumvec*dumvec
               rik2=SUM(dumvec)
               if(rik2.gt.thr)cycle

               dumvec=jkvec+s-t
               dumvec=dumvec*dumvec
               rjk2=SUM(dumvec)
               if(rjk2.gt.thr)cycle
               !rr0ik=sqrt(rik2)/r0ab(at(iat),at(kat))
               !rr0jk=sqrt(rjk2)/r0ab(at(jat),at(kat))


               rijk2=(rij2*rjk2*rik2)
               !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
               !damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

               rijk=sqrt(rijk2)
               fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
               rijk3=rijk*rijk2
               ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                  *(-rij2+rjk2+rik2)/(rijk3*rijk2) &
                  +1.0_wp/(rijk3)


               r9ijk=ang*fdmp/2.0_wp   !factor 1/2 for doublecounting
               eabc=eabc-r9ijk*c9

            enddo !sx

         enddo !tx
      enddo !kat
   enddo !iat
   ! And now kat=jat, but cycling throug all imagecells without t=s. and jat>iat going though all cells    (iat,jat,jat)
   ! But this counts only 1/2

   do iat=2,nat
      do jat=1,iat-1
         kat=jat
         ij=iat*(iat-1)/2+jat
         jk=jat*(jat-1)/2+kat
         ik=ij

         c6ij=c6ab(ij)
         c6ik=c6ij

         c6jk=c6ab(jk)
         ikvec=xyz(:,kat)-xyz(:,iat)
         ijvec=ikvec
         jkvec=0.0_wp
         cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
         cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
         cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
         cijk = cij*cjk*cik

         c9=-sqrt(c6ij*c6ik*c6jk)
         do concurrent(tx=-rep(1):rep(1), &
               &       ty=-rep(2):rep(2), &
               &       tz=-rep(3):rep(3))
            repmin(1)=max(-rep(1),tx-rep(1))
            repmax(1)=min(+rep(1),tx+rep(1))
            repmin(2)=max(-rep(2),ty-rep(2))
            repmax(2)=min(+rep(2),ty+rep(2))
            repmin(3)=max(-rep(3),tz-rep(3))
            repmax(3)=min(+rep(3),tz+rep(3))

            t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
            dumvec=ijvec+t
            dumvec=dumvec*dumvec
            rij2=SUM(dumvec)
            if(rij2.gt.thr)cycle

            !rr0ij=SQRT(rij2)/r0ab(at(iat),at(jat))

            do concurrent (sx=repmin(1):repmax(1), &
                  &        sy=repmin(2):repmax(2), &
                  &        sz=repmin(3):repmax(3))
               ! every result * 0.5
               IF (tx.eq.sx .and. ty.eq.sy  &
                  .and. tz.eq.sz) cycle
               s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
               dumvec=ikvec+s
               dumvec=dumvec*dumvec
               rik2=SUM(dumvec)
               if(rik2.gt.thr)cycle
               !rr0ik=SQRT(rik2)/r0ab(at(iat),at(kat))

               dumvec=jkvec+s-t
               dumvec=dumvec*dumvec
               rjk2=SUM(dumvec)
               if(rjk2.gt.thr)cycle
               !rr0jk=SQRT(rjk2)/r0ab(at(jat),at(kat))

               !              if (rij*rjk*rik.gt.thr)cycle

               rijk2=(rij2*rjk2*rik2)
               !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
               !damp9=1./(1.+6._wp*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

               rijk=sqrt(rijk2)
               fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
               rijk3=rijk*rijk2
               ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                  *(-rij2+rjk2+rik2)/(rijk2*rijk3) &
                  +1.0_wp/(rijk3)
               r9ijk=ang*fdmp/2.0_wp   !factor 1/2 for doublecounting
               eabc=eabc-r9ijk*c9

            enddo !sx

         enddo !tx
      enddo !kat
   enddo !iat


   ! and finally the self interaction iat=jat=kat all

   do iat=1,nat
      jat=iat
      kat=iat
      ijvec=0.0_wp
      ij=iat*(iat-1)/2+jat
      ik=iat*(iat-1)/2+kat
      jk=jat*(jat-1)/2+kat
      ikvec=ijvec
      jkvec=ikvec
      c6ij=c6ab(ij)
      c6ik=c6ij
      c6jk=c6ij
      c9=-sqrt(c6ij*c6ij*c6ij)
      cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
      cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
      cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
      cijk = cij*cjk*cik

      do concurrent ( tx=-rep(1):rep(1), &
            &         ty=-rep(2):rep(2), &
            &         tz=-rep(3):rep(3))
         repmin(1)=max(-rep(1),tx-rep(1))
         repmax(1)=min(+rep(1),tx+rep(1))
         repmin(2)=max(-rep(2),ty-rep(2))
         repmax(2)=min(+rep(2),ty+rep(2))
         repmin(3)=max(-rep(3),tz-rep(3))
         repmax(3)=min(+rep(3),tz+rep(3))
         if ((tx.eq.0) .and.(ty.eq.0) .and.(tz.eq.0))cycle !IF iat and jat are the same then cycle
         t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
         dumvec=t
         dumvec=dumvec*dumvec
         rij2=SUM(dumvec)
         if(rij2.gt.thr)cycle
         !rr0ij=SQRT(rij2)/r0ab(at(iat),at(jat))

         do concurrent (sx=repmin(1):repmax(1), &
               &        sy=repmin(2):repmax(2), &
               &        sz=repmin(3):repmax(3))
            if ((sx.eq.0) .and.( sy.eq.0) .and.( sz.eq.0))cycle !IF iat and kat are the same then cycle
            if ((sx.eq.tx) .and. (sy.eq.ty)  &
               .and. (sz.eq.tz)) cycle      !If kat and jat are the same then cycle

            ! every result * 1/6 becaues every triple is counted twice due to the two loops t and s going from -rep to rep -> *1/2
            !
            !plus 1/3 becaues every triple is three times in each unitcell
            s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
            dumvec=s
            dumvec=dumvec*dumvec
            rik2=SUM(dumvec)
            if(rik2.gt.thr)cycle
            !rr0ik=SQRT(rik2)/r0ab(at(iat),at(kat))

            dumvec=jkvec+s-t
            dumvec=dumvec*dumvec
            rjk2=SUM(dumvec)
            if(rjk2.gt.thr)cycle
            !rr0jk=SQRT(rjk2)/r0ab(at(jat),at(kat))

            rijk2=(rij2*rjk2*rik2)
            !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
            !damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

            rijk=sqrt(rijk2)
            fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
            rijk3=rijk*rijk2
            ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
               *(-rij2+rjk2+rik2)/(rijk2*rijk3) &
               +1.0_wp/(rijk3)
            r9ijk=ang*fdmp/6.0_wp
            eabc=eabc-c9*r9ijk

         enddo !sx
      enddo !tx


   enddo !iat


end subroutine abcappr_3d_dftd3_like_style

! compute D4 gradient under pbc
subroutine dispgrad_3d(mol,ndim,q,cn,dcndr,dcndL,rep,atm_rep,r_thr,atm_thr,par, &
      &                wf,g_a,g_c,refc6,mbd,g,sigma,eout,dqdr,dqdL,aout)
   use iso_fortran_env, only : wp => real64
   use mctc_constants
   use tbdef_molecule
   use pbc_tools
   implicit none
   type(tb_molecule),intent(in) :: mol
   integer, intent(in)  :: ndim
   real(wp),intent(in)  :: q(mol%n)
   real(wp),intent(in)  :: cn(mol%n)
   real(wp),intent(in)  :: dcndr(3,mol%n,mol%n)
   real(wp),intent(in)  :: dcndL(3,3,mol%n)
   integer, intent(in)  :: rep(3),atm_rep(3)
   real(wp),intent(in)  :: r_thr,atm_thr
   integer              :: tx,ty,tz
   real(wp)             :: t(3)
   real(wp)             :: aiw(23)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: wf,g_a,g_c
   real(wp),intent(in)  :: refc6(ndim,ndim)
   integer, intent(in)  :: mbd
   real(wp),intent(inout)        :: g(3,mol%n)
   real(wp),intent(inout)        :: sigma(3,3)
   real(wp),intent(out),optional :: eout
   real(wp),intent(in), optional :: dqdr(3,mol%n,mol%n+1)
   real(wp),intent(in), optional :: dqdL(3,3,mol%n+1)
   real(wp),intent(out),optional :: aout(23,mol%n)


   integer  :: i,ii,iii,j,jj,k,l,ia,ja,ij
   integer, allocatable :: itbl(:,:)
   real(wp) :: iz
   real(wp) :: qmod,eabc,ed
   real(wp) :: norm,dnorm
   real(wp) :: dexpw,expw
   real(wp) :: twf,tgw,r4r2ij
   real(wp) :: rij(3),r,r2,r4,r6,r8,R0
   real(wp) :: oor6,oor8,oor10,door6,door8,door10,cutoff
   real(wp) :: c8abns,disp,ddisp,x1,x2,x3
   real(wp) :: c6ii,c6ij,dic6ii,dic6ij,djc6ij,dizii,dizij,djzij
   real(wp) :: rcovij,expterm,den

   real(wp) :: drdx(3),dtmp,gwk,dgwk
   real(wp),allocatable :: r2ab(:)
   real(wp),allocatable :: dc6dcn(:)
   real(wp),allocatable :: zvec(:)
   real(wp),allocatable :: dzvec(:)
   real(wp),allocatable :: gw(:)
   real(wp),allocatable :: dgw(:)
   real(wp),allocatable :: dc6dq(:)
   real(wp),allocatable :: dzdq(:)
   real(wp) :: cn_thr,gw_thr

   parameter(cn_thr = 1600.0_wp)
   parameter(gw_thr=0.000001_wp)

   intrinsic :: present,sqrt,sum,maxval,exp,abs

   !  print'(" * Allocating local memory")'
   allocate( dc6dcn(mol%n),r2ab(mol%n*(mol%n+1)/2),dc6dq(mol%n),dzdq(ndim), &
      &      zvec(ndim),dzvec(ndim),gw(ndim),dgw(ndim), &
      &      source = 0.0_wp )
   allocate( itbl(7,mol%n), source = 0 )

   ed = 0.0_wp
   eabc = 0.0_wp

   k = 0
   do i = 1, mol%n
      do ii = 1, refn(mol%at(i))
         k = k+1
         itbl(ii,i) = k
      enddo
   enddo

!  print'(" * Entering first OMP section")'
!$OMP parallel default(none) &
!$omp private(i,ii,iii,ia,iz,k,norm,dnorm,twf,tgw,dexpw,expw,gwk,dgwk)  &
!$omp shared (mol,refn,refc,refcovcn,itbl,refq,wf,cn,g_a,g_c,q) &
!$omp shared (gw,dgw,zvec,dzvec,dzdq,r_thr)
!$omp do
   do i = 1, mol%n
      ia = mol%at(i)
      iz = zeff(ia)
      norm  = 0.0_wp
      dnorm = 0.0_wp
      do ii=1,refn(ia)
         do iii = 1, refc(ii,ia)
            twf = iii*wf
            tgw = cngw(twf,cn(i),refcovcn(ii,ia))
            norm  =  norm + tgw
            dnorm = dnorm + 2*twf*(refcovcn(ii,ia)-cn(i))*tgw
         enddo
      enddo
      norm = 1._wp/norm
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         dexpw=0.0_wp
         expw=0.0_wp
         do iii = 1, refc(ii,ia)
            twf = wf*iii
            tgw = cngw(twf,cn(i),refcovcn(ii,ia))
            expw  =  expw + tgw
            dexpw = dexpw + 2*twf*(refcovcn(ii,ia)-cn(i))*tgw
         enddo

         ! save
         gwk = expw*norm
         if (gwk.ne.gwk) then
            if (maxval(refcovcn(:refn(ia),ia)).eq.refcovcn(ii,ia)) then
               gwk = 1.0_wp
            else
               gwk = 0.0_wp
            endif
         endif
         zvec(k) = zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz) * gwk
         ! NEW: q=0 for ATM
         gw(k) =  zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,iz) * gwk

         dgwk = dexpw*norm-expw*dnorm*norm**2
         if (dgwk.ne.dgwk) then
            dgwk = 0.0_wp
         endif
         dzvec(k) = zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz) * dgwk
         dzdq(k) = dzeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,q(i)+iz) * gwk
         ! NEW: q=0 for ATM
         dgw(k) = zeta(g_a,gam(ia)*g_c,refq(ii,ia)+iz,iz) * dgwk
      enddo
   enddo
!$omp end do
!$omp end parallel

!$OMP parallel default(none) &
!$omp private(i,j,ia,ja,ij,k,l,c6ii,c6ij,dic6ii,dic6ij,djc6ij,disp,ddisp,dizii,dizij,djzij,  &
!$omp         rij,r2,r,r4r2ij,r0,oor6,oor8,oor10,door6,door8,door10, &
!$omp         t,tx,ty,tz,dtmp,drdx)  &
!$omp shared(mol,refn,itbl,zvec,dzvec,refc6,par,dzdq,rep,r_Thr) &
!$omp shared(r2ab) reduction(+:dc6dq,dc6dcn,ed,g,sigma)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      ia = mol%at(i)
      ! temps
      c6ij   = 0.0_wp
      dic6ij = 0.0_wp
      djc6ij = 0.0_wp
      dizij  = 0.0_wp
      djzij  = 0.0_wp
      ! all refs
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         do jj = 1, refn(ia)
            l = itbl(jj,i)
            c6ij   = c6ij   +  zvec(k) *  zvec(l) * refc6(k,l)
            dic6ij = dic6ij + dzvec(k) *  zvec(l) * refc6(k,l)
            djc6ij = djc6ij +  zvec(k) * dzvec(l) * refc6(k,l)
            dizij  = dizij  +  dzdq(k) *  zvec(l) * refc6(k,l)
            djzij  = djzij  +  zvec(k) *  dzdq(l) * refc6(k,l)
         enddo
      enddo
      ! i in primitive cell with i in images
      r4r2ij = 3*r4r2(ia)*r4r2(ia)
      r0 = par%a1*sqrt(r4r2ij) + par%a2
      do concurrent(tx = -rep(1):rep(1), &
            &       ty = -rep(2):rep(2), &
            &       tz = -rep(3):rep(3))
         ! cycle i with i interaction in same cell
         if (tx.eq.0.and.ty.eq.0.and.tz.eq.0) cycle
         t = [tx,ty,tz]
         rij = matmul(mol%lattice,t)
         r2  = sum( rij**2 )
         if (r2.gt.r_thr) cycle
         r   = sqrt(r2)
         oor6 = 1._wp/(r2**3+r0**6)
         oor8 = 1._wp/(r2**4+r0**8)
         oor10 = 1._wp/(r2**5+r0**10)
         door6 = -6*r2**2*r*oor6**2
         door8 = -8*r2**3*r*oor8**2
         door10 = -10*r2**4*r*oor10**2
         disp = par%s6*oor6 + par%s8*r4r2ij*oor8   &
            & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
         ddisp= par%s6*door6 + par%s8*r4r2ij*door8 &
            & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*door10
         ed = ed - c6ij*disp
         ! save this
         dtmp = c6ij*ddisp
         dc6dq(i)  = dc6dq(i)  + (dizij  + djzij )*disp
         dc6dcn(i) = dc6dcn(i) + (dic6ij + djc6ij)*disp
         drdx = rij/r
         sigma = sigma - dtmp * outer_prod_3x3(drdx,rij)
      enddo ! tx
      ! over all j atoms
      do j = 1, i-1
         ja = mol%at(j)
         ! temps
         c6ij   = 0.0_wp
         dic6ij = 0.0_wp
         djc6ij = 0.0_wp
         dizij  = 0.0_wp
         djzij  = 0.0_wp
         ! all refs
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c6ij   = c6ij   +  zvec(k) *  zvec(l) * refc6(k,l)
               dic6ij = dic6ij + dzvec(k) *  zvec(l) * refc6(k,l)
               djc6ij = djc6ij +  zvec(k) * dzvec(l) * refc6(k,l)
               dizij  = dizij  +  dzdq(k) *  zvec(l) * refc6(k,l)
               djzij  = djzij  +  zvec(k) *  dzdq(l) * refc6(k,l)
            enddo
         enddo
         r4r2ij = 3*r4r2(ia)*r4r2(ja)
         r0 = par%a1*sqrt(r4r2ij) + par%a2
         do concurrent(tx = -rep(1):rep(1), &
               &       ty = -rep(2):rep(2), &
               &       tz = -rep(3):rep(3))
            t = [tx,ty,tz]
            rij = mol%xyz(:,i) - mol%xyz(:,j) + matmul(mol%lattice,t)
            r2 = sum(rij**2)
            if (r2.gt.r_thr) cycle
            r = sqrt(r2)
            oor6 = 1._wp/(r2**3+r0**6)
            oor8 = 1._wp/(r2**4+r0**8)
            oor10 = 1._wp/(r2**5+r0**10)
            door6 = -6*r2**2*r*oor6**2
            door8 = -8*r2**3*r*oor8**2
            door10 = -10*r2**4*r*oor10**2
            disp = par%s6*oor6 + par%s8*r4r2ij*oor8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*oor10
            ddisp= par%s6*door6 + par%s8*r4r2ij*door8 &
               & + par%s10*49.0_wp/40.0_wp*r4r2ij**2*door10
            ed = ed - c6ij*disp
            ! save this
            dc6dq(i)  = dc6dq(i)  + dizij  *disp
            dc6dq(j)  = dc6dq(j)  + djzij  *disp
            dc6dcn(i) = dc6dcn(i) + dic6ij *disp
            dc6dcn(j) = dc6dcn(j) + djc6ij *disp
            dtmp = c6ij*ddisp
            drdx = rij/r
            g(:,i) = g(:,i) - dtmp * drdx
            g(:,j) = g(:,j) + dtmp * drdx

            sigma = sigma - dtmp * outer_prod_3x3(drdx,rij)
         enddo ! tx
      enddo
   enddo
!$omp enddo
!$omp end parallel

   eabc = 0.0_wp
   if (mbd > 0) then
      call dabcappr_3d(mol,ndim,par,gw,dgw,refc6,itbl,atm_rep,atm_thr, &
         &             g,sigma,dc6dcn,eabc)
      !dc6dcn = 0.0_wp
   endif

   if(present(dqdr)) then
      ! handle dqdr  :: gradient is exact e-11
      call dgemv('n',3*mol%n,mol%n,-1.0_wp,dqdr,3*mol%n,dc6dq,1,1.0_wp,g,1)
   endif
   if (mol%npbc > 0) then
      if(present(dqdL)) then
         ! handle dqdL
         call dgemv('n',3*3,mol%n,-1.0_wp,dqdL,3*3,dc6dq,1,1.0_wp,sigma,1)
      endif
   endif

   ! always handle dcndr :: gradient is exact e-11
   call dgemv('n',3*mol%n,mol%n,-1.0_wp,dcndr,3*mol%n,dc6dcn,1,1.0_wp,g,1)
   if (mol%npbc > 0) then
      call dgemv('n',3*3,mol%n,-1.0_wp,dcndL,3*3,dc6dcn,1,1.0_wp,sigma,1)
   endif

   !  print'(" * Dispersion all done, saving variables")'
   if (present(eout)) eout = ed + eabc

   if (present(aout)) then
      aout = 0._wp
      do i = 1, mol%n
         ia = mol%at(i)
         do ii = 1, refn(ia)
            aout(:,i) = aout(:,i) + zvec(k) * refal(:,ii,ia)
         enddo
      enddo
   endif
end subroutine dispgrad_3d

!> @brief calculates threebody dispersion gradient from C6 coefficients
subroutine dabcappr_3d(mol,ndim,par,zvec,dzvec,refc6,itbl,rep,r_thr, &
      &                g,sigma,dc6dcn,eabc)
   use tbdef_molecule
   use pbc_tools
   implicit none
   type(tb_molecule),intent(in) :: mol !< molecular structure information
   integer, intent(in)  :: ndim
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in)  :: zvec(ndim)
   real(wp),intent(in)  :: dzvec(ndim)
   real(wp),intent(in)  :: refc6(ndim,ndim)
   integer, intent(in)  :: itbl(7,mol%n)
   integer, intent(in)  :: rep(3)
   real(wp),intent(in)  :: r_thr
   real(wp),intent(inout) :: g(3,mol%n)
   real(wp),intent(inout) :: sigma(3,3)
   real(wp),intent(inout) :: dc6dcn(mol%n)
   real(wp),intent(out)   :: eabc

   integer  :: i,ii,ia,j,jj,ja,k,kk,ka,l,m,n
   integer  :: ij,jk,ik
   real(wp),allocatable :: c6ab(:),dc6ab(:,:)
   real(wp) :: r2ij,r2jk,r2ik,r
   real(wp) :: cii,cij,cjk,cik,ciii,ciij,cijk
   real(wp) :: c9iii,c9iij,c9ijk,oor9ijk,rijk
   real(wp) :: rij(3),rjk(3),rik(3)
   real(wp) :: dijfdmp,dikfdmp,djkfdmp
   real(wp) :: dijatm,dikatm,djkatm
   real(wp) :: dijoor9ijk,djkoor9ijk,dikoor9ijk
   real(wp) :: c6ij,dic6ij,djc6ij
   real(wp) :: dic9iii,dic9iij,djc9iij,dic9ijk,djc9ijk,dkc9ijk
   real(wp),parameter :: zero(3) = [0.0_wp,0.0_wp,0.0_wp]
   real(wp) :: gw_thr
   parameter(gw_thr=0.0001_wp)

   intrinsic :: present,sqrt

   allocate( c6ab(mol%n*(mol%n+1)/2),dc6ab(mol%n,mol%n), source = 0.0_wp )

   eabc = 0.0_wp

!$OMP parallel default(none) &
!$omp private(i,ia,j,ja,ij,r2ij,c6ij,dic6ij,djc6ij,k,l,r)  &
!$omp shared (mol,refn,itbl,refc6,zvec,dzvec) &
!$omp shared (c6ab,dc6ab)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      ia = mol%at(i)
      ij = i*(i-1)/2+i

      ! temps
      c6ij = 0.0_wp
      dic6ij = 0.0_wp
      djc6ij = 0.0_wp
      ! all refs
      do ii = 1, refn(ia)
         k = itbl(ii,i)
         do jj = 1, refn(ia)
            l = itbl(jj,i)
            c6ij = c6ij + zvec(k)*zvec(l)*refc6(k,l)
            dic6ij = dic6ij + dzvec(k)*zvec(l)*refc6(k,l)
            djc6ij = djc6ij + zvec(k)*dzvec(l)*refc6(k,l)
         enddo
      enddo
      ! save
      c6ab(ij) = c6ij
      dc6ab(i,i) = dic6ij! + djc6ij

      do j = 1, i-1
!        if(i.eq.j) cycle
         ja = mol%at(j)
         ij = i*(i-1)/2 + j

!        first check if we want this contribution
         !if(mol%dist(j,i)**2.gt.r_thr) cycle

         ! temps
         c6ij = 0.0_wp
         dic6ij = 0.0_wp
         djc6ij = 0.0_wp
         ! all refs
         do ii = 1, refn(ia)
            k = itbl(ii,i)
            do jj = 1, refn(ja)
               l = itbl(jj,j)
               c6ij = c6ij + zvec(k)*zvec(l)*refc6(k,l)
               dic6ij = dic6ij + dzvec(k)*zvec(l)*refc6(k,l)
               djc6ij = djc6ij + zvec(k)*dzvec(l)*refc6(k,l)
            enddo
         enddo
         ! save
         c6ab(ij) = c6ij
         dc6ab(i,j) = dic6ij
         dc6ab(j,i) = djc6ij
      enddo
   enddo
!$omp enddo
!$omp end parallel

   call dabcappr_3d_dftd3_like_style(mol%n,mol%at,mol%xyz,par,r_thr,rep,&
      &    mol%lattice,c6ab,dc6ab,eabc,dc6dcn,g,sigma)
   !call dabcappr_3d_wsc(mol,par,[2,2,2],r_thr,c6ab,dc6ab,g,dc6dcn,eabc)
   !call dabcappr_3d_bvk(mol,par,rep,r_thr,c6ab,dc6ab,g,dc6dcn,eabc)

end subroutine dabcappr_3d

subroutine dabcappr_3d_dftd3_like_style(nat,at,xyz,par,thr,rep,dlat,c6ab,dc6ab, &
      &    eabc,dc6dcn,g,sigma)
   use mctc_constants
   use pbc_tools
   implicit none
   integer, intent(in) :: nat
   integer, intent(in) :: at(nat)
   real(wp),intent(in) :: xyz(3,nat)
   type(dftd_parameter),intent(in) :: par
   real(wp),intent(in) :: thr
   integer, intent(in) :: rep(3)
   real(wp),intent(in) :: c6ab(nat*(nat+1)/2)
   real(wp),intent(in) :: dc6ab(nat,nat)    !dC6(iat,jat)/cCN(iat) in dc6ab(i,j) for ABC-grad
   real(wp),intent(in) :: dlat(3,3)

   real(wp),intent(inout) :: g(3,nat)
   real(wp),intent(inout) :: sigma(3,3)
   real(wp),intent(inout) :: eabc
   real(wp),intent(inout) :: dc6dcn(nat)

   real(wp),parameter :: six = 6.0_wp, oth = 1.0_wp/3.0_wp

   integer iat,jat,kat
   real(wp) x1
   real(wp) r2,r
   real(wp) fdmp,tmp1
   real(wp) :: rco,den,tmp,dtmp

   real(wp) rij(3),rik(3),rjk(3)
   real(wp), allocatable,dimension(:,:,:,:) ::  dc6dr  !d(E)/d(r_ij) derivative wrt. dist. iat-jat
   !dCN(jat)/d(r_ij)
   real(wp) :: r9ijk
   real(wp) vec(3)
   integer ij,ik,jk

   real(wp),dimension(3) ::ijvec,ikvec,jkvec,t,s,dumvec
   integer tx,ty,tz,sx,sy,sz
   real(wp) rij2,rik2,rjk2,c9,c6ij,c6ik,c6jk,rijk,rijk3
   real(wp) :: cij,cjk,cik,cijk
   real(wp) time1,time2,rijk2,dc9,dfdmp,dang,ang
   integer,dimension(3) :: repmin,repmax

   allocate(dc6dr(-rep(3):rep(3),-rep(2):rep(2), &
      &           -rep(1):rep(1),nat*(nat+1)/2))
   dc6dr = 0.0_wp


   !        write(*,*)'!!!!!!!!!!    THREEBODY  GRADIENT  !!!!!!!!!!'
   eabc=0.0_wp
   !        write(*,*)'thr:',sqrt(thr)

   iAt_ijk: do iat=3,nat
      jAt_ijk: do jat=2,iat-1
         ij=iat*(iat-1)/2+jat
         ijvec=xyz(:,jat)-xyz(:,iat)

         c6ij=c6ab(ij)
         kAt_ijk: do kat=1,jat-1
            ik=iat*(iat-1)/2+kat
            jk=jat*(jat-1)/2+kat
            ikvec=xyz(:,kat)-xyz(:,iat)
            jkvec=xyz(:,kat)-xyz(:,jat)

            c6ik=c6ab(ik)
            c6jk=c6ab(jk)
            c9=-sqrt(c6ij*c6ik*c6jk)
            cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
            cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
            cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
            cijk = cij*cjk*cik

            do concurrent (tx=-rep(1):rep(1), &
                  &        ty=-rep(2):rep(2), &
                  &        tz=-rep(3):rep(3))
               repmin(1)=max(-rep(1),tx-rep(1))
               repmax(1)=min(+rep(1),tx+rep(1))
               repmin(2)=max(-rep(2),ty-rep(2))
               repmax(2)=min(+rep(2),ty+rep(2))
               repmin(3)=max(-rep(3),tz-rep(3))
               repmax(3)=min(+rep(3),tz+rep(3))
               t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
               rij2=SUM((ijvec+t)*(ijvec+t))
               if(rij2.gt.thr)cycle

               !rr0ij=sqrt(rij2)/r0ab(at(iat),at(jat))


               do concurrent (sx=repmin(1):repmax(1), &
                     &        sy=repmin(2):repmax(2), &
                     &        sz=repmin(3):repmax(3))
                  s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
                  rik2=SUM((ikvec+s)*(ikvec+s))
                  if(rik2.gt.thr)cycle

                  dumvec=jkvec+s-t
                  rjk2=SUM(dumvec*dumvec)
                  if(rjk2.gt.thr)cycle
                  !rr0ik=sqrt(rik2)/r0ab(at(iat),at(kat))
                  !rr0jk=sqrt(rjk2)/r0ab(at(jat),at(kat))
                  rijk2=(rij2*rjk2*rik2)
                  ! first calculate the three components for the energy calculation fdmp
                  ! and ang
                  !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
                  !damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

                  rijk=sqrt(rijk2)
                  fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
                  rijk3=rijk*rijk2
                  ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                     *(-rij2+rjk2+rik2)/(rijk3*rijk2) &
                     +1.0_wp/(rijk3)

                  r9ijk=ang*fdmp
                  eabc=eabc-r9ijk*c9
                  !
                  !start calculating the gradient components dfdmp, dang and dc9

                  !dfdmp is the same for all three distances
                  !dfdmp=2._wp*alp9*(0.75_wp*r0av)**(alp9)*fdmp*fdmp
                  dfdmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2

                  !start calculating the derivatives of each part w.r.t. r_ij
                  r=sqrt(rij2)


                  dang=-0.375_wp*(rij2**3+rij2**2*(rjk2+rik2) &
                     +rij2*(3.0_wp*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2) &
                     -5.0*(rjk2-rik2)**2*(rjk2+rik2)) &
                     /(r*rijk3*rijk2)

                  tmp1=-dang*c9*fdmp+dfdmp/r*c9*ang
                  dc6dr(tz,ty,tx,ij) = dc6dr(tz,ty,tx,ij)-tmp1

                  !start calculating the derivatives of each part w.r.t. r_ik

                  r=sqrt(rik2)


                  dang=-0.375_wp*(rik2**3+rik2**2*(rjk2+rij2) &
                     +rik2*(3.0_wp*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2) &
                     -5.0*(rjk2-rij2)**2*(rjk2+rij2)) &
                     /(r*rijk3*rijk2)

                  tmp1=-dang*c9*fdmp+dfdmp/r*c9*ang
                  !                 tmp1=-dc9
                  dc6dr(sz,sy,sx,ik) = dc6dr(sz,sy,sx,ik)-tmp1

                  !
                  !start calculating the derivatives of each part w.r.t. r_jk

                  r=sqrt(rjk2)

                  dang=-0.375_wp*(rjk2**3+rjk2**2*(rik2+rij2) &
                     +rjk2*(3.0_wp*rik2**2+2.0*rik2*rij2+3.0*rij2**2) &
                     -5.0*(rik2-rij2)**2*(rik2+rij2)) &
                     /(r*rijk3*rijk2)

                  tmp1=-dang*c9*fdmp+dfdmp/r*c9*ang
                  dc6dr(sz-tz,sy-ty,sx-tx,jk) = dc6dr(sz-tz,sy-ty,sx-tx,jk)-tmp1

                  !calculating the CN derivative dE_disp(ijk)/dCN(i)

                  dc9=0.5_wp*c9*(dc6ab(iat,jat)/c6ij+dc6ab(iat,kat)/c6ik)
                  dc6dcn(iat)=dc6dcn(iat)+r9ijk*dc9

                  dc9=0.5_wp*c9*(dc6ab(jat,iat)/c6ij+dc6ab(jat,kat)/c6jk)
                  dc6dcn(jat)=dc6dcn(jat)+r9ijk*dc9

                  dc9=0.5_wp*c9*(dc6ab(kat,iat)/c6ik+dc6ab(kat,jat)/c6jk)
                  dc6dcn(kat)=dc6dcn(kat)+r9ijk*dc9


               enddo !sz
            enddo !tx
         enddo kAt_ijk
      enddo jAt_ijk
   enddo iAt_ijk

   ! Now the interaction with jat=iat of the triples iat,iat,kat
   iAt_iik: do iat=2,nat
      jat=iat
      ij=iat*(iat-1)/2+jat
      ijvec=0.0_wp

      c6ij=c6ab(ij)
      kAt_iik: do kat=1,iat-1
         jk=jat*(jat-1)/2+kat
         ik=jk

         c6ik=c6ab(ik)
         c6jk=c6ik
         ikvec=xyz(:,kat)-xyz(:,iat)
         jkvec=ikvec
         c9=-sqrt(c6ij*c6ik*c6jk)
         cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
         cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
         cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
         cijk = cij*cjk*cik
         do concurrent (tx=-rep(1):rep(1), &
               &        ty=-rep(2):rep(2), &
               &        tz=-rep(3):rep(3))
            repmin(1)=max(-rep(1),tx-rep(1))
            repmax(1)=min(+rep(1),tx+rep(1))
            repmin(2)=max(-rep(2),ty-rep(2))
            repmax(2)=min(+rep(2),ty+rep(2))
            repmin(3)=max(-rep(3),tz-rep(3))
            repmax(3)=min(+rep(3),tz+rep(3))
            IF (tx.eq.0 .and. ty.eq.0 .and. tz.eq.0) cycle
            t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
            dumvec=t
            rij2=SUM(dumvec*dumvec)
            if(rij2.gt.thr)cycle

            !rr0ij=sqrt(rij2)/r0ab(at(iat),at(jat))

            do concurrent (sx=repmin(1):repmax(1), &
                  &        sy=repmin(2):repmax(2), &
                  &        sz=repmin(3):repmax(3))
               ! every result * 0.5

               s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
               dumvec=ikvec+s
               dumvec=dumvec*dumvec
               rik2=SUM(dumvec)
               if(rik2.gt.thr)cycle

               dumvec=jkvec+s-t
               dumvec=dumvec*dumvec
               rjk2=SUM(dumvec)
               if(rjk2.gt.thr)cycle
               !rr0ik=sqrt(rik2)/r0ab(at(iat),at(kat))
               !rr0jk=sqrt(rjk2)/r0ab(at(jat),at(kat))


               rijk2=(rij2*rjk2*rik2)
               !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
               !damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

               rijk=sqrt(rijk2)
               fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
               rijk3=rijk*rijk2
               ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                  *(-rij2+rjk2+rik2)/(rijk3*rijk2) &
                  +1.0_wp/(rijk3)


               r9ijk=ang*fdmp/2.0_wp   !factor 1/2 for doublecounting
               eabc=eabc-r9ijk*c9

               !              iat=jat
               !dfdmp=2._wp*alp9*(0.75_wp*r0av)**(alp9)*fdmp*fdmp
               dfdmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2

               !start calculating the derivatives of each part w.r.t. r_ij
               r=sqrt(rij2)

               dang=-0.375_wp*(rij2**3+rij2**2*(rjk2+rik2) &
                  +rij2*(3.0_wp*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2) &
                  -5.0*(rjk2-rik2)**2*(rjk2+rik2)) &
                  /(r*rijk3*rijk2)

               tmp1=-dang*c9*fdmp+dfdmp/r*c9*ang
               dc6dr(tz,ty,tx,ij) = dc6dr(tz,ty,tx,ij)-tmp1/2.0

               !start calculating the derivatives of each part w.r.t. r_ik
               r=sqrt(rik2)


               dang=-0.375_wp*(rik2**3+rik2**2*(rjk2+rij2) &
                  +rik2*(3.0_wp*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2) &
                  -5.0*(rjk2-rij2)**2*(rjk2+rij2)) &
                  /(r*rijk3*rijk2)

               tmp1=-dang*c9*fdmp+dfdmp/r*c9*ang
               dc6dr(sz,sy,sx,ik) = dc6dr(sz,sy,sx,ik)-tmp1/2.0
               !
               !start calculating the derivatives of each part w.r.t. r_ik
               r=sqrt(rjk2)

               dang=-0.375_wp*(rjk2**3+rjk2**2*(rik2+rij2) &
                  +rjk2*(3.0_wp*rik2**2+2.0*rik2*rij2+3.0*rij2**2) &
                  -5.0*(rik2-rij2)**2*(rik2+rij2)) &
                  /(r*rijk3*rijk2)

               tmp1=-dang*c9*fdmp+dfdmp/r*c9*ang

               dc6dr(sz-tz,sy-ty,sx-tx,jk) = dc6dr(sz-tz,sy-ty,sx-tx,jk)-tmp1/2.0

               dc9=0.5_wp*c9*(dc6ab(iat,jat)/c6ij+dc6ab(iat,kat)/c6ik)
               dc6dcn(iat)=dc6dcn(iat)+r9ijk*dc9

               dc9=0.5_wp*c9*(dc6ab(jat,iat)/c6ij+dc6ab(jat,kat)/c6jk)
               dc6dcn(jat)=dc6dcn(jat)+r9ijk*dc9

               dc9=0.5_wp*c9*(dc6ab(kat,iat)/c6ik+dc6ab(kat,jat)/c6jk)
               dc6dcn(kat)=dc6dcn(kat)+r9ijk*dc9




            enddo !sx

         enddo !tx
      enddo kAt_iik
   enddo iAt_iik
   ! And now kat=jat, but cycling throug all imagecells without t=s. and jat>iat going though all cells    (iat,jat,jat)
   ! But this counts only 1/2

   iAt_ijj: do iat=2,nat
      jAt_ijj: do jat=1,iat-1
         kat=jat
         ij=iat*(iat-1)/2+jat
         jk=jat*(jat-1)/2+kat
         ik=ij

         c6ij=c6ab(ij)
         c6ik=c6ij

         c6jk=c6ab(jk)
         ikvec=xyz(:,kat)-xyz(:,iat)
         ijvec=ikvec
         jkvec=0.0_wp
         cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
         cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
         cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
         cijk = cij*cjk*cik

         c9=-sqrt(c6ij*c6ik*c6jk)
         do concurrent(tx=-rep(1):rep(1), &
               &       ty=-rep(2):rep(2), &
               &       tz=-rep(3):rep(3))
            repmin(1)=max(-rep(1),tx-rep(1))
            repmax(1)=min(+rep(1),tx+rep(1))
            repmin(2)=max(-rep(2),ty-rep(2))
            repmax(2)=min(+rep(2),ty+rep(2))
            repmin(3)=max(-rep(3),tz-rep(3))
            repmax(3)=min(+rep(3),tz+rep(3))

            t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
            dumvec=ijvec+t
            dumvec=dumvec*dumvec
            rij2=SUM(dumvec)
            if(rij2.gt.thr)cycle

            !rr0ij=SQRT(rij2)/r0ab(at(iat),at(jat))

            do concurrent (sx=repmin(1):repmax(1), &
                  &        sy=repmin(2):repmax(2), &
                  &        sz=repmin(3):repmax(3))
               ! every result * 0.5
               if (tx.eq.sx .and. ty.eq.sy .and. tz.eq.sz) cycle
               s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
               dumvec=ikvec+s
               dumvec=dumvec*dumvec
               rik2=SUM(dumvec)
               if(rik2.gt.thr)cycle
               !rr0ik=SQRT(rik2)/r0ab(at(iat),at(kat))

               dumvec=jkvec+s-t
               dumvec=dumvec*dumvec
               rjk2=SUM(dumvec)
               if(rjk2.gt.thr)cycle
               !rr0jk=SQRT(rjk2)/r0ab(at(jat),at(kat))

               !              if (rij*rjk*rik.gt.thr)cycle

               rijk2=(rij2*rjk2*rik2)
               !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
               !damp9=1./(1.+6._wp*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

               rijk=sqrt(rijk2)
               fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
               rijk3=rijk*rijk2
               ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
                  *(-rij2+rjk2+rik2)/(rijk2*rijk3) &
                  +1.0_wp/(rijk3)
               r9ijk=ang*fdmp/2.0_wp   !factor 1/2 for doublecounting
               eabc=eabc-r9ijk*c9


               !              jat=kat
               !dfdmp=2._wp*alp9*(0.75_wp*r0av)**(alp9)*fdmp*fdmp
               dfdmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2
               !start calculating the derivatives of each part w.r.t. r_ij
               r=sqrt(rij2)

               dang=-0.375_wp*(rij2**3+rij2**2*(rjk2+rik2) &
                  +rij2*(3.0_wp*rjk2**2+2.0_wp*rjk2*rik2+3.0_wp*rik2**2) &
                  -5.0_wp*(rjk2-rik2)**2*(rjk2+rik2)) &
                  /(r*rijk3*rijk2)

               tmp1=-dang*c9*fdmp+dfdmp/r*c9*ang
               dc6dr(tz,ty,tx,ij) = dc6dr(tz,ty,tx,ij)-tmp1/2.0_wp

               !start calculating the derivatives of each part w.r.t. r_ik
               r=sqrt(rik2)


               dang=-0.375_wp*(rik2**3+rik2**2*(rjk2+rij2) &
                  +rik2*(3.0_wp*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2) &
                  -5.0*(rjk2-rij2)**2*(rjk2+rij2)) &
                  /(r*rijk3*rijk2)

               tmp1=-dang*c9*fdmp+dfdmp/r*c9*ang
               !                 tmp1=-dc9
               dc6dr(sz,sy,sx,ik) = dc6dr(sz,sy,sx,ik)-tmp1/2.0_wp
               !
               !start calculating the derivatives of each part w.r.t. r_jk
               r=sqrt(rjk2)

               dang=-0.375_wp*(rjk2**3+rjk2**2*(rik2+rij2) &
                  +rjk2*(3.0_wp*rik2**2+2.0*rik2*rij2+3.0*rij2**2) &
                  -5.0_wp*(rik2-rij2)**2*(rik2+rij2)) &
                  /(r*rijk3*rijk2)

               tmp1=-dang*c9*fdmp+dfdmp/r*c9*ang
               dc6dr(sz-tz,sy-ty,sx-tx,jk) = dc6dr(sz-tz,sy-ty,sx-tx,jk)-tmp1/2.0_wp

               !calculating the CN derivative dE_disp(ijk)/dCN(i)

               dc9=0.5_wp*c9*(dc6ab(iat,jat)/c6ij+dc6ab(iat,kat)/c6ik)
               dc6dcn(iat)=dc6dcn(iat)+r9ijk*dc9

               dc9=0.5_wp*c9*(dc6ab(jat,iat)/c6ij+dc6ab(jat,kat)/c6jk)
               dc6dcn(jat)=dc6dcn(jat)+r9ijk*dc9

               dc9=0.5_wp*c9*(dc6ab(kat,iat)/c6ik+dc6ab(kat,jat)/c6jk)
               dc6dcn(kat)=dc6dcn(kat)+r9ijk*dc9




            enddo !sx

         enddo !tx
      enddo jAt_ijj
   enddo iAt_ijj

   ! And finally the self interaction iat=jat=kat all

   iAt_iii: do iat=1,nat
      jat=iat
      kat=iat
      ijvec=0.0_wp
      ij=iat*(iat-1)/2+jat
      ik=iat*(iat-1)/2+kat
      jk=jat*(jat-1)/2+kat
      ikvec=ijvec
      jkvec=ikvec
      c6ij=c6ab(ij)
      c6ik=c6ij
      c6jk=c6ij
      c9=-sqrt(c6ij*c6ij*c6ij)
      cij  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(jat)))+par%a2
      cik  = par%a1*sqrt(3._wp*r4r2(at(iat))*r4r2(at(kat)))+par%a2
      cjk  = par%a1*sqrt(3._wp*r4r2(at(jat))*r4r2(at(kat)))+par%a2
      cijk = cij*cjk*cik

      do concurrent ( tx=-rep(1):rep(1), &
            &         ty=-rep(2):rep(2), &
            &         tz=-rep(3):rep(3))
         repmin(1)=max(-rep(1),tx-rep(1))
         repmax(1)=min(+rep(1),tx+rep(1))
         repmin(2)=max(-rep(2),ty-rep(2))
         repmax(2)=min(+rep(2),ty+rep(2))
         repmin(3)=max(-rep(3),tz-rep(3))
         repmax(3)=min(+rep(3),tz+rep(3))
         ! IF iat and jat are the same then cycle
         if ((tx.eq.0) .and.(ty.eq.0) .and.(tz.eq.0))cycle
         t=tx*dlat(:,1)+ty*dlat(:,2)+tz*dlat(:,3)
         dumvec=t
         dumvec=dumvec*dumvec
         rij2=SUM(dumvec)
         if(rij2.gt.thr)cycle
         !rr0ij=SQRT(rij2)/r0ab(at(iat),at(jat))

         do concurrent (sx=repmin(1):repmax(1), &
               &        sy=repmin(2):repmax(2), &
               &        sz=repmin(3):repmax(3))
            ! if iat and kat are the same then cycle
            if ((sx.eq.0) .and.( sy.eq.0) .and.( sz.eq.0))cycle
            ! If kat and jat are the same then cycle
            if ((sx.eq.tx) .and. (sy.eq.ty) .and. (sz.eq.tz)) cycle

            ! every result * 1/6 becaues every triple is counted twice due
            ! to the two loops t and s going from -rep to rep -> *1/2
            !
            !plus 1/3 becaues every triple is three times in each unitcell
            s=sx*dlat(:,1)+sy*dlat(:,2)+sz*dlat(:,3)
            dumvec=s
            dumvec=dumvec*dumvec
            rik2=SUM(dumvec)
            if(rik2.gt.thr)cycle
            !rr0ik=SQRT(rik2)/r0ab(at(iat),at(kat))

            dumvec=jkvec+s-t
            dumvec=dumvec*dumvec
            rjk2=SUM(dumvec)
            if(rjk2.gt.thr)cycle
            !rr0jk=SQRT(rjk2)/r0ab(at(jat),at(kat))

            rijk2=(rij2*rjk2*rik2)
            !r0av=(rr0ij*rr0ik*rr0jk)**(1.0_wp/3.0_wp)
            !damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

            rijk=sqrt(rijk2)
            fdmp = 1.0_wp/(1.0_wp+six*((cijk/rijk)**oth)**par%alp)
            rijk3=rijk*rijk2
            ang=0.375_wp*(rij2+rjk2-rik2)*(rij2-rjk2+rik2) &
               *(-rij2+rjk2+rik2)/(rijk2*rijk3) &
               +1.0_wp/(rijk3)
            r9ijk=ang*fdmp/6.0_wp
            eabc=eabc-c9*r9ijk

            !                          iat=jat=kat
            !dfdmp=2._wp*alp9*(0.75_wp*r0av)**(alp9)*fdmp*fdmp
            dfdmp = -(oth*six*par%alp*((cijk/rijk)**oth)**par%alp)*fdmp**2
            !start calculating the derivatives of each part w.r.t. r_ij

            r=sqrt(rij2)
            dang=-0.375_wp*(rij2**3+rij2**2*(rjk2+rik2) &
               +rij2*(3.0_wp*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2) &
               -5.0*(rjk2-rik2)**2*(rjk2+rik2)) &
               /(r*rijk3*rijk2)


            tmp1=-dang*c9*fdmp+dfdmp/r*c9*ang
            dc6dr(tz,ty,tx,ij) = dc6dr(tz,ty,tx,ij)-tmp1/6.0_wp

            !start calculating the derivatives of each part w.r.t. r_ik

            r=sqrt(rik2)

            dang=-0.375_wp*(rik2**3+rik2**2*(rjk2+rij2) &
               +rik2*(3.0_wp*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2) &
               -5.0*(rjk2-rij2)**2*(rjk2+rij2)) &
               /(r*rijk3*rijk2)

            tmp1=-dang*c9*fdmp+dfdmp/r*c9*ang
            dc6dr(sz,sy,sx,ik) = dc6dr(sz,sy,sx,ik)-tmp1/6.0_wp
            !
            !start calculating the derivatives of each part w.r.t. r_jk

            r=sqrt(rjk2)
            dang=-0.375_wp*(rjk2**3+rjk2**2*(rik2+rij2) &
               +rjk2*(3.0_wp*rik2**2+2.0*rik2*rij2+3.0*rij2**2) &
               -5.0*(rik2-rij2)**2*(rik2+rij2)) &
               /(r*rijk3*rijk2)

            tmp1=-dang*c9*fdmp+dfdmp/r*c9*ang
            dc6dr(sz-tz,sy-ty,sx-tx,jk) = dc6dr(sz-tz,sy-ty,sx-tx,jk)-tmp1/6.0_wp


            !calculating the CN derivative dE_disp(ijk)/dCN(i)

            dc9=0.5_wp*c9*(dc6ab(iat,jat)/c6ij+dc6ab(iat,kat)/c6ik)
            dc6dcn(iat)=dc6dcn(iat)+r9ijk*dc9

            dc9=0.5_wp*c9*(dc6ab(jat,iat)/c6ij+dc6ab(jat,kat)/c6jk)
            dc6dcn(jat)=dc6dcn(jat)+r9ijk*dc9

            dc9=0.5_wp*c9*(dc6ab(kat,iat)/c6ik+dc6ab(kat,jat)/c6jk)
            dc6dcn(kat)=dc6dcn(kat)+r9ijk*dc9

         enddo !sx
      enddo !tx


   enddo iAt_iii

   ! After calculating all derivatives dE/dr_ij w.r.t. distances,
   ! the grad w.r.t. the coordinates is calculated dE/dr_ij * dr_ij/dxyz_i
!$omp parallel default(none) &
!$omp shared(nat,rep,dlat,thr,xyz,dc6dr) &
!$omp private(iat,jat,ij,t,rij,r2,r,x1,vec,tx,ty,tz) &
!$omp reduction(+:g,sigma)
!$omp do schedule(dynamic)
   iAt_ij: do iat=2,nat
      jAt_ij: do jat=1,iat-1
         ij=iat*(iat-1)/2+jat
         !       write(*,'(3E17.6,XX,2I2)'),dc6dr(0,0,-1:1,lin(iat,jat)),iat,jat
         do concurrent(tx=-rep(1):rep(1),&
               &       ty=-rep(2):rep(2),&
               &       tz=-rep(3):rep(3))
            t=[tx,ty,tz]

            rij=xyz(:,jat)-xyz(:,iat)+matmul(dlat,t)
            r2=sum(rij*rij)
            if (r2.gt.thr) cycle
            r=sqrt(r2)

            x1=dc6dr(tz,ty,tx,ij)!+dtmp*(dc6dcn(iat)+dc6dcn(jat))

            vec=x1*rij/r
            g(:,iat)=g(:,iat)+vec
            g(:,jat)=g(:,jat)-vec
            sigma = sigma - outer_prod_3x3(vec,rij)

         enddo !tx
      enddo jAt_ij
   enddo iAt_ij
!$omp enddo

!$omp do schedule(dynamic)
   iAt_ii: do iat=1,nat
      jat = iat
      ij=iat*(iat-1)/2+iat
      !       write(*,'(3E17.6,XX,2I2)'),dc6dr(0,0,-1:1,lin(iat,jat)),iat,jat
      do concurrent(tx=-rep(1):rep(1),&
            &       ty=-rep(2):rep(2),&
            &       tz=-rep(3):rep(3))
         if (tx==0 .and. ty==0 .and. tz==0) cycle
         t=[tx,ty,tz]

         rij=matmul(dlat,t)
         r2=sum(rij*rij)
         if (r2.gt.thr) cycle
         r=sqrt(r2)

         x1=dc6dr(tz,ty,tx,ij)!+dtmp*(dc6dcn(iat)+dc6dcn(jat))

         vec=x1*rij/r
         sigma = sigma - outer_prod_3x3(vec,rij)

      enddo !tx
   enddo iAt_ii
!$omp enddo
!$omp end parallel

end subroutine dabcappr_3d_dftd3_like_style


end module dftd4
