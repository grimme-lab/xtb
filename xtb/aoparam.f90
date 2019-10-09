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

module aoparam
   use iso_fortran_env, wp => real64
   implicit none
   private :: wp

   integer,private,parameter :: max_sh = 10
   type :: tb_parameter
      real(wp) :: en = 1.50_wp
      real(wp) :: mc = 0.0_wp
      real(wp) :: gam = 0.0_wp
      real(wp) :: gam3 = 0.0_wp
      real(wp) :: alp0 = 0.0_wp
      real(wp) :: rad = 0.0_wp
      real(wp) :: wll(max_sh) = 0.0_wp
      real(wp) :: rep(2) = 0.0_wp
      real(wp) :: polyr(4) = 0.0_wp
      real(wp) :: cxb = 0.0_wp
      real(wp) :: ao_exp(max_sh) = 0.0_wp
      real(wp) :: ao_lev(max_sh) = 0.0_wp
      real(wp) :: lpar(0:2) = 0.0_wp
      real(wp) :: kcnat(0:2) = 0.0_wp
      real(wp) :: kqat(3) = 0.0_wp
      real(wp) :: radaes = 5.0_wp
      real(wp) :: dpolc = 0.0_wp
      real(wp) :: qpolc = 0.0_wp
      integer  :: ao_pqn(max_sh) = 0
      integer  :: ao_l(max_sh) = 0
      integer  :: ao_n = 0
      integer  :: ao_typ(max_sh) = 0
      integer  :: metal = 0
      integer  :: cnval = 0
      character(len=30) :: timestp='------------------------------'
   end type tb_parameter

   character,private,parameter :: flag = '$'
   character,private,parameter :: space = ' '
   character,private,parameter :: equal = '='
   character,private,parameter :: hash = '#'
   character,private,parameter :: dot = '.'
   character(len=*),private,parameter :: flag_end = '$end'

!! ========================================================================
!  default parameter sets, included directly in the library/executable
   integer,private,parameter :: max_elem = 94
!! ========================================================================
!  non-selfconsistent vTB/xTB Hamiltonian
!!! ------------------------------------------------------------------------
!!  sTDA-xTB parametrisation for non-selfconsistent vTB calculation
!   real(wp),dimension(25),private :: stda1_globparameter
!   real(wp),dimension(max_elem,max_elem),private :: stda1_pairparameter = 1.0_wp
!   type(tb_parameter),dimension(max_elem),private :: stda1_atomparameter
!   include 'param_stda1.f90'
!!! ------------------------------------------------------------------------
!!  sTDA-xTB parametrisation for non-selfconsistent xTB calculation
!   real(wp),dimension(25),private :: stda2_globparameter
!   real(wp),dimension(max_elem,max_elem),private :: stda2_pairparameter = 1.0_wp
!   type(tb_parameter),dimension(max_elem),private :: stda2_atomparameter
!   include 'param_stda2.f90'
!!! ========================================================================
!!  GFN1 Hamiltonian
!!! ------------------------------------------------------------------------
!!  first geometry, frequency and non-covalent interactions parametrisation
!   real(wp),dimension(25),private :: gfn1_globparameter
!   real(wp),dimension(max_elem,max_elem),private :: gfn1_pairparameter = 1.0_wp
!   type(tb_parameter),dimension(max_elem),private :: gfn1_atomparameter
!   include 'param_gfn1.f90'
!!! ------------------------------------------------------------------------
!!  special purpose parametrisaton for ionisation potentials and electon
!!  affinities as extension on the GFN1 parametrisation
!   real(wp),dimension(25),private :: ipea_globparameter
!   real(wp),dimension(max_elem,max_elem),private :: ipea_pairparameter = 1.0_wp
!   type(tb_parameter),dimension(max_elem),private :: ipea_atomparameter
!   include 'param_ipea.f90'
!!! ========================================================================
!!  GFN2 Hamiltonian
!!! ------------------------------------------------------------------------
!!  second geometry, frequency and non-covalent interactions parametrisation
!!  with anisotropic electrostatics and selfconsistent dispersion
!   real(wp),dimension(25),private :: gfn2_globparameter
!   real(wp),dimension(max_elem,max_elem),private :: gfn2_pairparameter = 1.0_wp
!   type(tb_parameter),dimension(max_elem),private :: gfn2_atomparameter
!   include 'param_gfn2.f90'

   public

   real(wp) :: en(94)
   real(wp) :: mc(94)
   real(wp) :: gam(94)
   real(wp) :: gam3(94)
   real(wp) :: alp0(94)
   real(wp) :: rad(94)
   real(wp) :: wll(94,10)
   real(wp) :: rep(2,94) 
   real(wp) :: polyr(4,94)
   real(wp) :: cxb(94)
   real(wp) :: ao_exp(10,94)
   real(wp) :: ao_lev(10,94)
   real(wp) :: lpar(0:2,94)
   real(wp) :: kpair(94,94)
   !real(wp) :: gam3l(0:3)
   real(wp) :: kcnat(0:2,94)
   real(wp) :: kqat(3,94)
   real(wp) :: radaes(94)
   real(wp) :: dpolc(94)
   real(wp) :: qpolc(94)
   integer  :: ao_pqn(10,94)
   integer  :: ao_l(10,94)
   integer  :: ao_n(94)
   integer  :: ao_typ(10,94)
   integer  :: metal(94)
   integer  :: cnval(94)
   character(len=30) :: timestp(94)

   data metal / &
   &  0,                                                 0, &! H-He
   &  1, 1,                               1, 0, 0, 0, 0, 0, &! Li-Ne
   &  1, 1,                               1, 0, 0, 0, 0, 0, &! Na-Ar
   &  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, &! K-Kr
   &  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, &! Rb-Xe
   &  1, 1, &! Cs/Ba
   &        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &!La-Lu
   &        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, &! Lu-Rn
   &  1, 1, 1, 1, 1, 1, 1, 1  /! Fr-Pu

! pauling EN's 
      data en/ 2.200_wp,3.000_wp,0.980_wp,1.570_wp,2.040_wp,2.550_wp,3.040_wp,3.440_wp,3.980_wp &
     &          ,4.500_wp,0.930_wp,1.310_wp,1.610_wp,1.900_wp,2.190_wp,2.580_wp,3.160_wp,3.500_wp &
     &          ,0.820_wp,1.000_wp,1.360_wp,1.540_wp,1.630_wp,1.660_wp,1.550_wp,1.830_wp,1.880 &
     &          ,1.910_wp,1.900_wp,1.650_wp,1.810_wp,2.010_wp,2.180_wp,2.550_wp,2.960_wp,3.000_wp &
     &          ,0.820_wp,0.950_wp,1.220_wp,1.330_wp,1.600_wp,2.160_wp,1.900_wp,2.200_wp,2.280_wp &
     &          ,2.200_wp,1.930_wp,1.690_wp,1.780_wp,1.960_wp,2.050_wp,2.100_wp,2.660_wp,2.600_wp &
     &  ,0.79_wp,0.89_wp,1.10_wp,1.12_wp,1.13_wp,1.14_wp,1.15_wp,1.17_wp,1.18_wp,1.20_wp,1.21_wp,1.22_wp &
     &  ,1.23_wp,1.24_wp,1.25_wp,1.26_wp,1.27_wp,1.3_wp,1.5_wp,2.36_wp,1.9_wp,2.2_wp,2.20_wp,2.28_wp,2.54_wp &
     &  ,2.00_wp,1.62_wp,2.33_wp,2.02_wp,2.0_wp,2.2_wp,2.2_wp,8*1.5_wp/

! D3 radii      

      data radaes/ &
     &   0.80628308_wp, 1.15903197_wp, 3.02356173_wp, 2.36845659_wp, 1.94011865_wp, &
     &   1.88972601_wp, 1.78894056_wp, 1.58736983_wp, 1.61256616_wp, 1.68815527_wp, &
     &   3.52748848_wp, 3.14954334_wp, 2.84718717_wp, 2.62041997_wp, 2.77159820_wp, &
     &   2.57002732_wp, 2.49443835_wp, 2.41884923_wp, 4.43455700_wp, 3.88023730_wp, &
     &   3.35111422_wp, 3.07395437_wp, 3.04875805_wp, 2.77159820_wp, 2.69600923_wp, &
     &   2.62041997_wp, 2.51963467_wp, 2.49443835_wp, 2.54483100_wp, 2.74640188_wp, &
     &   2.82199085_wp, 2.74640188_wp, 2.89757982_wp, 2.77159820_wp, 2.87238349_wp, &
     &   2.94797246_wp, 4.76210950_wp, 4.20778980_wp, 3.70386304_wp, 3.50229216_wp, &
     &   3.32591790_wp, 3.12434702_wp, 2.89757982_wp, 2.84718717_wp, 2.84718717_wp, &
     &   2.72120556_wp, 2.89757982_wp, 3.09915070_wp, 3.22513231_wp, 3.17473967_wp, &
     &   3.17473967_wp, 3.09915070_wp, 3.32591790_wp, 3.30072128_wp, 5.26603625_wp, &
     &   4.43455700_wp, 4.08180818_wp, 3.70386304_wp, 3.98102289_wp, 3.95582657_wp, &
     &   3.93062995_wp, 3.90543362_wp, 3.80464833_wp, 3.82984466_wp, 3.80464833_wp, &
     &   3.77945201_wp, 3.75425569_wp, 3.75425569_wp, 3.72905937_wp, 3.85504098_wp, &
     &   3.67866672_wp, 3.45189952_wp, 3.30072128_wp, 3.09915070_wp, 2.97316878_wp, &
     &   2.92277614_wp, 2.79679452_wp, 2.82199085_wp, 2.84718717_wp, 3.32591790_wp, &
     &   3.27552496_wp, 3.27552496_wp, 3.42670319_wp, 3.30072128_wp, 3.47709584_wp, &
     &   3.57788113_wp, 5.06446567_wp, 4.56053862_wp, 4.20778980_wp, 3.98102289_wp, &
     &   3.82984466_wp, 3.85504098_wp, 3.88023730_wp, 3.90543362 /

!  Semiempirical Evaluation of the GlobalHardness of the Atoms of 103
!  Elementsof the Periodic Table Using the MostProbable Radii as their Size
!  Descriptors
!  DULAL C. GHOSH, NAZMUL ISLAM
!  2009 in Wiley InterScience (www.interscience.wiley.com).
!  DOI 10.1002/qua.22202
!  values in the paper multiplied by two because (ii:ii)=(IP-EA)=d^2 E/dN^2 
!  but the hardness definition they use is 1/2d^2 E/dN^2 (in Eh)
      data gam/ &
     &  0.47259288_wp,0.92203391_wp,0.17452888_wp,0.25700733_wp,0.33949086_wp,0.42195412_wp, &
     &  0.50438193_wp,0.58691863_wp,0.66931351_wp,0.75191607_wp,0.17964105_wp,0.22157276_wp, &
     &  0.26348578_wp,0.30539645_wp,0.34734014_wp,0.38924725_wp,0.43115670_wp,0.47308269_wp, &
     &  0.17105469_wp,0.20276244_wp,0.21007322_wp,0.21739647_wp,0.22471039_wp,0.23201501_wp, &
     &  0.23933969_wp,0.24665638_wp,0.25398255_wp,0.26128863_wp,0.26859476_wp,0.27592565_wp, &
     &  0.30762999_wp,0.33931580_wp,0.37235985_wp,0.40273549_wp,0.43445776_wp,0.46611708_wp, &
     &  0.15585079_wp,0.18649324_wp,0.19356210_wp,0.20063311_wp,0.20770522_wp,0.21477254_wp, &
     &  0.22184614_wp,0.22891872_wp,0.23598621_wp,0.24305612_wp,0.25013018_wp,0.25719937_wp, &
     &  0.28784780_wp,0.31848673_wp,0.34912431_wp,0.37976593_wp,0.41040808_wp,0.44105777_wp, &
     &  0.05019332_wp,0.06762570_wp,0.08504445_wp,0.10247736_wp,0.11991105_wp,0.13732772_wp, &
     &  0.15476297_wp,0.17218265_wp,0.18961288_wp,0.20704760_wp,0.22446752_wp,0.24189645_wp, &
     &  0.25932503_wp,0.27676094_wp,0.29418231_wp,0.31159587_wp,0.32902274_wp,0.34592298_wp, &
     &  0.36388048_wp,0.38130586_wp,0.39877476_wp,0.41614298_wp,0.43364510_wp,0.45104014_wp, &
     &  0.46848986_wp,0.48584550_wp,0.12526730_wp,0.14268677_wp,0.16011615_wp,0.17755889_wp, &
     &  0.19497557_wp,0.21240778_wp,0.07263525_wp,0.09422158_wp,0.09920295_wp,0.10418621_wp, &
     &  0.14235633_wp,0.16394294_wp,0.18551941_wp,0.22370139 /

! COVALENT RADII
! based on "Atomic Radii of the Elements," M. Mantina, R. Valero, 
! C. J. Cramer, and D. G. Truhlar,
! in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
! edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
! corrected Nov. 17, 2010 for the 92nd edition.

      data rad / &
     &  0.32_wp,0.37_wp,1.30_wp,0.99_wp,0.84_wp,0.75_wp,0.71_wp,0.64_wp,0.60_wp, &
     &  0.62_wp,1.60_wp,1.40_wp,1.24_wp,1.14_wp,1.09_wp,1.04_wp,1.00_wp,1.01_wp, &
     &  2.00_wp,1.74_wp,1.59_wp,1.48_wp,1.44_wp,1.30_wp,1.29_wp,1.24_wp,1.18_wp, &
     &  1.17_wp,1.22_wp,1.20_wp,1.23_wp,1.20_wp,1.20_wp,1.18_wp,1.17_wp,1.16_wp, &
     &  2.15_wp,1.90_wp,1.76_wp,1.64_wp,1.56_wp,1.46_wp,1.38_wp,1.36_wp,1.34_wp, &
     &  1.30_wp,1.36_wp,1.40_wp,1.42_wp,1.40_wp,1.40_wp,1.37_wp,1.36_wp,1.36_wp, &
     &  2.38_wp,2.06_wp,1.94_wp,1.84_wp,1.90_wp,1.88_wp,1.86_wp,1.85_wp,1.83_wp, &
     &  1.82_wp,1.81_wp,1.80_wp,1.79_wp,1.77_wp,1.77_wp,1.78_wp,1.74_wp,1.64_wp, &
     &  1.58_wp,1.50_wp,1.41_wp,1.36_wp,1.32_wp,1.30_wp,1.30_wp,1.32_wp,1.44_wp, &
     &  1.45_wp,1.50_wp,1.42_wp,1.48_wp,1.46_wp,2.42_wp,2.11_wp,2.01_wp,1.90_wp, &
     &  1.84_wp,1.83_wp,1.80_wp,1.80_wp /


!     data gam3/&
!    &-0.02448,0.178614,0.194034,0.154068,0.173892,0.16716,0.156306, &
!    &0.161466,0.163314,0.170862,0.256128,0.18906,0.146310,0.136686,&
!    &0.123558,0.122070,0.119424,0.115368,0.175938,0.152214,0.24030,&
!    &0.204654,0.186552,0.178824,0.174474,0.205038,0.188772,0.180462,&
!    &0.180354,0.176508,0.158250,0.139212,0.131868,0.119712,0.115476,&
!    &0.108444,0.091032,0.07698,0.102642,0.095346,0.088266,0.086364,&
!    &0.085254,0.088242,0.087774,0.088470,0.091314,0.090372,0.110862,&
!    &0.093588,0.079908,0.074082,0.069342,0.064638,0.077826,0.0590214,&
!    &0.073614,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
!    &0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
!    &0.000000,0.085236,0.07782,0.074112,0.07206,0.070188,0.069246,&
!    &0.072042,0.07038,0.072570,0.108096,0.084150,0.070350,0.06954,&
!    &0.064866,0.0596874,0.000000,0.000000,0.000000,0.000000,0.000000,&
!    &0.000000,0.000000,0.000000 /

!     data kcnat/&
!& 0.15_wp,0.15_wp,&
!& 0.0_wp,0.0_wp,0.0_wp,0.15_wp,0.15_wp,0.15_wp,0.05_wp,0.0_wp,             &! 2-10
!& 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp, 0.0_wp, 0.20_wp, 0.0_wp,             &! 11-18
!& 0.0_wp,0.0_wp,                                                     &! 19-20
!& 1.5_wp,1.5_wp,0.0_wp,0.0_wp,0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,-0.1_wp, 0.0_wp,&! 21-30
!& 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.2_wp, 0.0_wp,                            &! 31-36 
!& 0.0_wp,0.0_wp,                                                     &! 37-38
!& 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,&! 39-48
!& 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.1_wp, 0.0_wp,                            &! 49-54 
!& 0.0_wp,0.0_wp,                                                     &! 55-56
!& 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,&! 57-71
!& 0.0_wp,0.0_wp,0.0_wp,0.0_wp,                                         &! 57-71
!& 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,&! 72-80
!& 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp, 0.0_wp,                            &! 81-86 
!& 8*0.0_wp  /


      data    cnval   / & ! normal CN used in CN dep. AES damping
     & 1,                                                             1, &
     & 1, 2,                                           3, 3, 3, 2, 1, 1, &
     & 1, 2,                                           3, 3, 3, 3, 1, 1, &
     & 1, 2, 4,          4, 6, 6, 6, 6, 6, 4, 4,  2,   3, 3, 3, 3, 1, 1, &
     & 1, 2, 4,          4, 6, 6, 6, 6, 6, 4, 4,  2,   3, 3, 3, 3, 1, 1, &
     & 1, 2, 4,14*6,     4, 6, 6, 6, 6, 6, 4, 4,  2,   3, 3, 3, 3, 1, 1, &
     & 8*0 /

contains

subroutine use_parameterset(name,globpar,exist)
   implicit none
   character(len=*),intent(in) :: name
   logical,intent(out)  :: exist
   real(wp),intent(out) :: globpar(25)
   exist = .false.
   select case(name)
   case('.param_stda1.xtb')
!      call copy_parameterset(globpar,stda1_globparameter, &
!      &  stda1_pairparameter,stda1_atomparameter)
      call copy_stda1_parameterset(globpar)
   case('.param_stda2.xtb')
!      call copy_parameterset(globpar,stda2_globparameter, &
!      &  stda2_pairparameter,stda2_atomparameter)
      call copy_stda2_parameterset(globpar)
   case('.param_gfn.xtb')
!      call copy_parameterset(globpar,gfn1_globparameter, &
!      &  gfn1_pairparameter,gfn1_atomparameter)
      call copy_gfn1_parameterset(globpar)
   case('.param_ipea.xtb')
!      call copy_parameterset(globpar,ipea_globparameter, &
!      &  ipea_pairparameter,ipea_atomparameter)
      call copy_ipea_parameterset(globpar)
   case('.param_gfn2.xtb')
!      call copy_parameterset(globpar,gfn2_globparameter, &
!      &  gfn2_pairparameter,gfn2_atomparameter)
      call copy_gfn2_parameterset(globpar)
   case default
      return
   end select
   exist = .true.
contains
include 'paramcopy_stda1.inc'
include 'paramcopy_stda2.inc'
include 'paramcopy_gfn1.inc'
include 'paramcopy_ipea.inc'
include 'paramcopy_gfn2.inc'
end subroutine use_parameterset


subroutine copy_parameterset(globpar,globparameter,pairparameter,atomparameter)
   implicit none
   real(wp),intent(out) :: globpar(25)
   real(wp),intent(in)  :: globparameter(25)
   real(wp),intent(in)  :: pairparameter(max_elem,max_elem)
   type(tb_parameter),intent(in) :: atomparameter(max_elem)
   integer :: i

   globpar           = globparameter
   kpair             = pairparameter
   en(1:94)          = atomparameter % en
   mc(1:94)          = atomparameter % mc
   gam(1:94)         = atomparameter % gam
   gam3(1:94)        = atomparameter % gam3
   rad(1:94)         = atomparameter % rad
   do i = 1,10
      wll(1:94,i)    = atomparameter % wll(i)
   enddo
   rep(1,1:94)       = atomparameter % rep(1)
   rep(2,1:94)       = atomparameter % rep(2)
   do i = 1,4
      polyr(i,1:94)   = atomparameter % polyr(i)
   enddo
   cxb(1:94)         = atomparameter % cxb
   do i = 1,10
      ao_exp(i,1:94) = atomparameter % ao_exp(i)
      ao_lev(i,1:94) = atomparameter % ao_lev(i)
      ao_pqn(i,1:94) = atomparameter % ao_pqn(i)
      ao_l(i,1:94)   = atomparameter % ao_l(i)
      ao_typ(i,1:94) = atomparameter % ao_typ(i)
   enddo
   do i = 0,2
      lpar(i,1:94)   = atomparameter % lpar(i)
      kcnat(i,1:94)  = atomparameter % kcnat(i)
   enddo
   radaes(1:94)      = atomparameter % radaes
   dpolc(1:94)       = atomparameter % dpolc
   qpolc(1:94)       = atomparameter % qpolc
   ao_n(1:94)        = atomparameter % ao_n
   metal(1:94)       = atomparameter % metal
   cnval(1:94)       = atomparameter % cnval
   timestp(1:94)     = atomparameter % timestp

end subroutine copy_parameterset

subroutine rd_param(fname)
   use iso_fortran_env, only : output_unit,iostat_end
   use readin, only : mirror_line,get_value
   implicit none
   character(len=*),intent(in)  :: fname
   character(len=:),allocatable :: line
   character(len=:),allocatable :: key
   character(len=:),allocatable :: val
   character(len=:),allocatable :: newname
   integer :: id
   integer :: ic
   integer :: ie
   integer :: iz
   integer :: err
   logical :: exist

   inquire(file=fname,exist=exist)
   if (.not.exist) then
      call raise('S','could not find parameter file '''//fname//'''',1)
      return
   endif
   id = 42
   open(unit=id,file=fname,status='old')

   call mirror_line(id,-1,line,err)
   readparam: do
      if (index(line,flag).eq.1) then
         select case(line(2:))
         case('globpar'); call rdglobpar(line,id,err)
         case('pairpar'); call rdpairpar(line,id,err)
         case default
            if (index(line,'Z').eq.2) then
               ie = index(line,space)
               val = trim(adjustl(line(ie+1:6)))
               if (get_value(val,iz)) then
                  timestp(iz) = line(7:)
                  call rdelempar(line,id,iz,err)
               endif
            else
               call mirror_line(id,-1,line,err)
            endif
         end select
      else
         call mirror_line(id,-1,line,err)
      endif ! not a keyword -> ignore
      if (err.eq.iostat_end) exit readparam
   enddo readparam

   close(id)

end subroutine rd_param

subroutine rdglobpar(line,id,err)
   use readin, only : mirror_line,get_value
   use mctc_strings, only : lowercase
   implicit none
   character(len=:),allocatable,intent(out) :: line
   integer,intent(in)  :: id
   integer,intent(out) :: err
   character(len=:),allocatable :: key
   character(len=:),allocatable :: val
   integer :: ie
   real(wp) :: globpar(25)
   real(wp) :: ddum

   do
      call mirror_line(id,-1,line,err)
      if (err.eq.iostat_end) exit
      if (index(line,flag_end).ne.0) exit

      ! find the first space
      ie = index(line,space)
      if ((line.eq.'').or.(ie.eq.0)) cycle
      key = lowercase(trim(line(:ie-1)))
      val = trim(adjustl(line(ie+1:)))
      select case(key)
      case('ks');      if (get_value(val,ddum)) globpar( 1) = ddum
      case('kp');      if (get_value(val,ddum)) globpar( 2) = ddum
      case('kd');      if (get_value(val,ddum)) globpar( 3) = ddum
      case('kf');      if (get_value(val,ddum)) globpar( 4) = ddum
      case('kdiffa');  if (get_value(val,ddum)) globpar( 5) = ddum
      case('kdiffb');  if (get_value(val,ddum)) globpar( 6) = ddum
      case('wllscal'); if (get_value(val,ddum)) globpar( 7) = ddum
      case('gscal');   if (get_value(val,ddum)) globpar( 8) = ddum
      case('zcnf');    if (get_value(val,ddum)) globpar( 9) = ddum
      case('tscal');   if (get_value(val,ddum)) globpar(10) = ddum
      case('kcn');     if (get_value(val,ddum)) globpar(11) = ddum
      case('fpol');    if (get_value(val,ddum)) globpar(12) = ddum
      case('ken');     if (get_value(val,ddum)) globpar(13) = ddum
      case('lshift');  if (get_value(val,ddum)) globpar(14) = ddum
      case('lshifta'); if (get_value(val,ddum)) globpar(15) = ddum
      case('split');   if (get_value(val,ddum)) globpar(16) = ddum
      case('zpf');     if (get_value(val,ddum)) globpar(17) = ddum
      case('alphaj');  if (get_value(val,ddum)) globpar(18) = ddum
      case('kexpo');   if (get_value(val,ddum)) globpar(19) = ddum
      case('dispa');   if (get_value(val,ddum)) globpar(20) = ddum
      case('dispb');   if (get_value(val,ddum)) globpar(21) = ddum
      case('dispc');   if (get_value(val,ddum)) globpar(22) = ddum
      case('dispatm'); if (get_value(val,ddum)) globpar(23) = ddum
      case('xbdamp');  if (get_value(val,ddum)) globpar(24) = ddum
      case('xbrad');   if (get_value(val,ddum)) globpar(25) = ddum
      case default
      end select
   enddo

end subroutine rdglobpar

subroutine rdpairpar(line,id,err)
   use readin, only : mirror_line,get_value
   implicit none
   character(len=:),allocatable,intent(out) :: line
   integer,intent(in)  :: id
   integer,intent(out) :: err
   character(len=:),allocatable :: key
   character(len=:),allocatable :: val
   integer :: ie

   do
      call mirror_line(id,-1,line,err)
      if (err.eq.iostat_end) exit
      if (index(line,flag_end).ne.0) exit

      ! find the first space
      ie = index(line,space)
      if ((line.eq.'').or.(ie.eq.0)) cycle
      key = trim(line(:ie-1))
      val = trim(adjustl(line(ie+1:)))
   enddo

end subroutine rdpairpar

subroutine rdelempar(line,id,iz,err)
   use readin, only : mirror_line,get_value
   use mctc_strings, only : lowercase
   implicit none
   character(len=:),allocatable,intent(out) :: line
   integer,intent(in)  :: id
   integer,intent(in)  :: iz
   integer,intent(out) :: err
   character(len=:),allocatable :: key
   character(len=:),allocatable :: val
   integer :: ie
   real(wp) :: ddum

   do
      call mirror_line(id,-1,line,err)
      if (err.eq.iostat_end) return
      if (index(line,flag_end).ne.0) return

      ! find the equal sign
      ie = index(line,equal)
      if ((line.eq.'').or.(ie.eq.0)) cycle
      key = lowercase(trim(line(:ie-1)))
      val = trim(adjustl(line(ie+1:)))
      select case(key)
      case('ao')
      case('lev')
      case('exp')
      case('en');    if (get_value(val,ddum))    en  (iz) = ddum
      case('gam');   if (get_value(val,ddum))   gam  (iz) = ddum
      case('alpg');  if (get_value(val,ddum))  alp0  (iz) = ddum
      case('epr');   if (get_value(val,ddum))    mc  (iz) = ddum
      case('gam3');  if (get_value(val,ddum))  gam3  (iz) = ddum*0.1_wp
      case('cxb');   if (get_value(val,ddum))   cxb  (iz) = ddum*0.1_wp
      case('dpol');  if (get_value(val,ddum)) dpolc  (iz) = ddum*0.01_wp
      case('qpol');  if (get_value(val,ddum)) qpolc  (iz) = ddum*0.01_wp
      case('repa');  if (get_value(val,ddum))   rep(1,iz) = ddum
      case('repb');  if (get_value(val,ddum))   rep(2,iz) = ddum
      case('polys'); if (get_value(val,ddum)) polyr(1,iz) = ddum
      case('polyp'); if (get_value(val,ddum)) polyr(2,iz) = ddum
      case('polyd'); if (get_value(val,ddum)) polyr(3,iz) = ddum
      case('polyf'); if (get_value(val,ddum)) polyr(4,iz) = ddum
      case('lpars'); if (get_value(val,ddum))  lpar(0,iz) = ddum*0.1_wp
      case('lparp'); if (get_value(val,ddum))  lpar(1,iz) = ddum*0.1_wp
      case('lpard'); if (get_value(val,ddum))  lpar(2,iz) = ddum*0.1_wp
      case('kcns');  if (get_value(val,ddum)) kcnat(0,iz) = ddum*0.1_wp
      case('kcnp');  if (get_value(val,ddum)) kcnat(1,iz) = ddum*0.1_wp
      case('kcnd');  if (get_value(val,ddum)) kcnat(2,iz) = ddum*0.1_wp
      case('kqs');   if (get_value(val,ddum))  kqat(1,iz) = ddum
      case('kqp');   if (get_value(val,ddum))  kqat(2,iz) = ddum
      case('kqd');   if (get_value(val,ddum))  kqat(3,iz) = ddum
      end select
   enddo

end subroutine rdelempar

end module aoparam
