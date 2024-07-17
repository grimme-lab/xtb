! This file is part of xtb.
!
! Copyright (C) 2019-2020 Stefan Grimme
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

module xtb_gfnff_param
   use xtb_mctc_accuracy, only : wp, sp
   use xtb_gfnff_data, only : TGFFData, init
   use xtb_gfnff_generator, only : TGFFGenerator
   use xtb_gfnff_topology, only : TGFFTopology
   implicit none
   private
   public :: ini, gff_print, make_chrg
   public :: gfnff_set_param, gfnff_load_param, gfnff_read_param
   public :: gfnff_param_alloc, gfnff_param_dealloc, gfnff_thresholds
   public :: gffVersion


   !> Available versions of the force field implemented
   type :: TGFFVersionEnum

      !> Version published in Angew. Chem. 2020
      integer :: angewChem2020 = 1

      !> Adding improved amide description to Angew. Chem. 2020 version
      integer :: angewChem2020_1 = 2

      !> Fixing generation of topological charges
      integer :: angewChem2020_2 = 3

      !> Harmonic potential version of the Angew. Chem. 2020
      integer :: harmonic2020 = -1

      !> Periodic Boundary Conditions with Original Parameters
      integer :: mcgfnff2023 = 4

   end type TGFFVersionEnum

   !> Actual enumerator for force field versions
   type(TGFFVersionEnum), parameter :: gffVersion = TGFFVersionEnum()


   !========================================================================
   ! all fixed parameters/potetnial parameters/thresholds/cut-offs
   ! used in gfnff_ini and set in gfnff_set_param
   !------------------------------------------------------------------------

   !loicals
   logical  :: ini       !make a neghbour lsit or use existing one?
   logical  :: gff_print = .true. !shows timing of energy + gradient
   logical  :: make_chrg = .true. !generates new eeq chareges based on current geometry

   real(wp), parameter :: chi_angewChem2020(103) = [&
      & 1.227054_wp, 1.451412_wp, 0.813363_wp, 1.062841_wp, 1.186499_wp, &
      & 1.311555_wp, 1.528485_wp, 1.691201_wp, 1.456784_wp, 1.231037_wp, &
      & 0.772989_wp, 1.199092_wp, 1.221576_wp, 1.245964_wp, 1.248942_wp, &
      & 1.301708_wp, 1.312474_wp, 1.247701_wp, 0.781237_wp, 0.940834_wp, &
      & 0.950000_wp, 0.974455_wp, 0.998911_wp, 1.023366_wp, 1.047822_wp, &
      & 1.072277_wp, 1.096733_wp, 1.121188_wp, 1.145644_wp, 1.170099_wp, &
      & 1.205357_wp, 1.145447_wp, 1.169499_wp, 1.253293_wp, 1.329909_wp, &
      & 1.116527_wp, 0.950975_wp, 0.964592_wp, 0.897786_wp, 0.932824_wp, &
      & 0.967863_wp, 1.002901_wp, 1.037940_wp, 1.072978_wp, 1.108017_wp, &
      & 1.143055_wp, 1.178094_wp, 1.213132_wp, 1.205076_wp, 1.075529_wp, &
      & 1.206919_wp, 1.303658_wp, 1.332656_wp, 1.179317_wp, 0.789115_wp, &
      & 0.798704_wp, 0.993208_wp, 0.907847_wp, 0.836114_wp, 0.778008_wp, & ! 60
      & 0.733529_wp, 0.702678_wp, 0.685455_wp, 0.681858_wp, 0.691889_wp, &
      & 0.715548_wp, 0.752834_wp, 0.803747_wp, 0.868288_wp, 0.946457_wp, & ! 70
      & 1.038252_wp, 1.128780_wp, 1.129764_wp, 1.130747_wp, 1.131731_wp, &
      & 1.132714_wp, 1.133698_wp, 1.134681_wp, 1.135665_wp, 1.136648_wp, & ! 80
      & 1.061832_wp, 1.053084_wp, 1.207830_wp, 1.236314_wp, 1.310129_wp, & 
      & 1.157380_wp, 0.789115_wp, 0.798704_wp, 1.053384_wp, 1.056040_wp, & ! 90
      & 1.058772_wp, 1.061580_wp, 1.064463_wp, 1.067422_wp, 1.070456_wp, &
      & 1.073566_wp, 1.076751_wp, 1.080012_wp, 1.083349_wp, 1.086761_wp, & ! 100
      & 1.090249_wp, 1.093812_wp, 1.097451_wp & ! 103 Lr
      &]

   real(wp), parameter :: gam_angewChem2020(103) = [&
      &-0.448428_wp, 0.131022_wp, 0.571431_wp, 0.334622_wp,-0.089208_wp, &
      &-0.025895_wp,-0.027280_wp,-0.031236_wp,-0.159892_wp, 0.074198_wp, &
      & 0.316829_wp, 0.326072_wp, 0.069748_wp,-0.120184_wp,-0.193159_wp, &
      &-0.182428_wp,-0.064093_wp, 0.061914_wp, 0.318112_wp, 0.189248_wp, &
      &-0.104172_wp,-0.082038_wp,-0.059903_wp,-0.037769_wp,-0.015635_wp, &
      & 0.006500_wp, 0.028634_wp, 0.050768_wp, 0.072903_wp, 0.095037_wp, &
      & 0.131140_wp, 0.097006_wp,-0.065744_wp,-0.058394_wp, 0.063307_wp, &
      & 0.091652_wp, 0.386337_wp, 0.530677_wp,-0.030705_wp,-0.020787_wp, &
      &-0.010869_wp,-0.000951_wp, 0.008967_wp, 0.018884_wp, 0.028802_wp, &
      & 0.038720_wp, 0.048638_wp, 0.058556_wp, 0.036488_wp, 0.077711_wp, &
      & 0.077025_wp, 0.004547_wp, 0.039909_wp, 0.082630_wp, 0.485375_wp, &
      & 0.416264_wp,-0.007277_wp,-0.007862_wp,-0.008388_wp,-0.008855_wp, & ! 60
      &-0.009263_wp,-0.009611_wp,-0.009901_wp,-0.010132_wp,-0.010304_wp, &
      &-0.010417_wp,-0.010471_wp,-0.010465_wp,-0.010401_wp,-0.010278_wp, & ! 70
      &-0.010096_wp,-0.008716_wp,-0.006220_wp,-0.003724_wp,-0.001228_wp, &
      & 0.001267_wp, 0.003763_wp, 0.006259_wp, 0.008755_wp, 0.011251_wp, & ! 80
      & 0.020477_wp,-0.056566_wp, 0.051943_wp, 0.076708_wp, 0.000273_wp, &
      &-0.068929_wp, 0.485375_wp, 0.416264_wp,-0.009993_wp,-0.010075_wp, & ! 90
      &-0.010138_wp,-0.010179_wp,-0.010201_wp,-0.010201_wp,-0.010182_wp, &
      &-0.010141_wp,-0.010080_wp,-0.009999_wp,-0.009897_wp,-0.009775_wp, & ! 100
      &-0.009632_wp,-0.009468_wp,-0.009284_wp & ! 103
      &]

   real(wp), parameter :: cnf_angewChem2020(103) = [&
      & 0.008904_wp, 0.004641_wp, 0.048324_wp, 0.080316_wp,-0.051990_wp, &
      & 0.031779_wp, 0.132184_wp, 0.157353_wp, 0.064120_wp, 0.036540_wp, &
      &-0.000627_wp, 0.005412_wp, 0.018809_wp, 0.016329_wp, 0.012149_wp, &
      & 0.021484_wp, 0.014212_wp, 0.014939_wp, 0.003597_wp, 0.032921_wp, &
      &-0.021804_wp,-0.022797_wp,-0.023789_wp,-0.024782_wp,-0.025775_wp, &
      &-0.026767_wp,-0.027760_wp,-0.028753_wp,-0.029745_wp,-0.030738_wp, &
      &-0.004189_wp,-0.011113_wp,-0.021305_wp,-0.012311_wp, 0.049781_wp, &
      &-0.040533_wp, 0.012872_wp, 0.021056_wp,-0.003395_wp, 0.000799_wp, &
      & 0.004992_wp, 0.009186_wp, 0.013379_wp, 0.017573_wp, 0.021766_wp, &
      & 0.025960_wp, 0.030153_wp, 0.034347_wp,-0.000052_wp,-0.039776_wp, &
      & 0.006661_wp, 0.050424_wp, 0.068985_wp, 0.023470_wp,-0.024950_wp, &
      &-0.033006_wp, 0.062573_wp, 0.057117_wp, 0.052085_wp, 0.047480_wp, & ! 60
      & 0.043299_wp, 0.039544_wp, 0.036215_wp, 0.033311_wp, 0.030833_wp, &
      & 0.028780_wp, 0.027153_wp, 0.025951_wp, 0.025174_wp, 0.024823_wp, & ! 70
      & 0.024898_wp, 0.053300_wp, 0.047628_wp, 0.041955_wp, 0.036282_wp, &
      & 0.030610_wp, 0.024937_wp, 0.019264_wp, 0.013592_wp, 0.007919_wp, & ! 80
      & 0.006383_wp,-0.089155_wp,-0.001293_wp, 0.019269_wp, 0.074803_wp, &
      & 0.016657_wp,-0.024950_wp,-0.033006_wp, 0.065498_wp, 0.058290_wp, & ! 90
      & 0.052304_wp, 0.047540_wp, 0.043999_wp, 0.041680_wp, 0.040584_wp, &
      & 0.040710_wp, 0.042058_wp, 0.044628_wp, 0.048421_wp, 0.053436_wp, & ! 100
      & 0.059674_wp, 0.067134_wp, 0.075816_wp & ! 103
      &]

   real(wp), parameter :: alp_angewChem2020(103) = [&
      & 0.585069_wp, 0.432382_wp, 0.628636_wp, 0.743646_wp, 1.167323_wp, &
      & 0.903430_wp, 1.278388_wp, 0.905347_wp, 1.067014_wp, 2.941513_wp, &
      & 0.687680_wp, 0.792170_wp, 1.337040_wp, 1.251409_wp, 1.068295_wp, &
      & 1.186476_wp, 1.593532_wp, 2.056749_wp, 0.674196_wp, 0.868052_wp, &
      & 0.575052_wp, 0.613424_wp, 0.651796_wp, 0.690169_wp, 0.728541_wp, &
      & 0.766913_wp, 0.805285_wp, 0.843658_wp, 0.882030_wp, 0.920402_wp, &
      & 0.877178_wp, 1.422350_wp, 1.405901_wp, 1.646860_wp, 2.001970_wp, &
      & 2.301695_wp, 1.020617_wp, 0.634141_wp, 0.652752_wp, 0.668845_wp, &
      & 0.684938_wp, 0.701032_wp, 0.717125_wp, 0.733218_wp, 0.749311_wp, &
      & 0.765405_wp, 0.781498_wp, 0.797591_wp, 1.296844_wp, 1.534068_wp, &
      & 1.727781_wp, 1.926871_wp, 2.175548_wp, 2.177702_wp, 0.977079_wp, &
      & 0.770260_wp, 0.695335_wp, 0.653051_wp, 0.617284_wp, 0.588034_wp, & ! 60
      & 0.565301_wp, 0.549084_wp, 0.539383_wp, 0.536199_wp, 0.539532_wp, &
      & 0.549382_wp, 0.565747_wp, 0.588630_wp, 0.618029_wp, 0.653945_wp, & ! 70
      & 0.696377_wp, 0.757075_wp, 0.756778_wp, 0.756480_wp, 0.756183_wp, &
      & 0.755886_wp, 0.755589_wp, 0.755291_wp, 0.754994_wp, 0.754697_wp, & ! 80
      & 0.868029_wp, 1.684375_wp, 2.001040_wp, 2.067331_wp, 2.228923_wp, &
      & 1.874218_wp, 0.977079_wp, 0.770260_wp, 0.717658_wp, 0.732120_wp, & ! 90
      & 0.743405_wp, 0.751515_wp, 0.756447_wp, 0.758203_wp, 0.756783_wp, &
      & 0.752186_wp, 0.744413_wp, 0.733463_wp, 0.719337_wp, 0.702034_wp, & ! 100
      & 0.681554_wp, 0.657898_wp, 0.631066_wp & ! 103
      &]

   real(wp), parameter :: bond_angewChem2020(103) = [&
      & 0.417997_wp, 0.258490_wp, 0.113608_wp, 0.195935_wp, 0.231217_wp, &
      & 0.385248_wp, 0.379257_wp, 0.339249_wp, 0.330706_wp, 0.120319_wp, & ! 10
      & 0.127255_wp, 0.173647_wp, 0.183796_wp, 0.273055_wp, 0.249044_wp, &
      & 0.290653_wp, 0.218744_wp, 0.034706_wp, 0.136353_wp, 0.192467_wp, & ! 20
      & 0.335860_wp, 0.314452_wp, 0.293044_wp, 0.271636_wp, 0.250228_wp, &
      & 0.228819_wp, 0.207411_wp, 0.186003_wp, 0.164595_wp, 0.143187_wp, & ! 30
      & 0.212434_wp, 0.210451_wp, 0.219870_wp, 0.224618_wp, 0.272206_wp, &
      & 0.147864_wp, 0.150000_wp, 0.150000_wp, 0.329501_wp, 0.309632_wp, & ! 40
      & 0.289763_wp, 0.269894_wp, 0.250025_wp, 0.230155_wp, 0.210286_wp, &
      & 0.190417_wp, 0.170548_wp, 0.150679_wp, 0.192977_wp, 0.173411_wp, & ! 50
      & 0.186907_wp, 0.192891_wp, 0.223202_wp, 0.172577_wp, 0.150000_wp, &
      & 0.150000_wp, 0.370682_wp, 0.368511_wp, 0.366339_wp, 0.364168_wp, & ! 60
      & 0.361996_wp, 0.359825_wp, 0.357654_wp, 0.355482_wp, 0.353311_wp, &
      & 0.351139_wp, 0.348968_wp, 0.346797_wp, 0.344625_wp, 0.342454_wp, & ! 70
      & 0.340282_wp, 0.338111_wp, 0.305540_wp, 0.272969_wp, 0.240398_wp, &
      & 0.207828_wp, 0.175257_wp, 0.142686_wp, 0.110115_wp, 0.077544_wp, & ! 80
      & 0.108597_wp, 0.148422_wp, 0.183731_wp, 0.192274_wp, 0.127706_wp, &
      & 0.086756_wp, 0.150000_wp, 0.150000_wp, 0.370682_wp, 0.368511_wp, & ! 90
      & 0.366339_wp, 0.364168_wp, 0.361996_wp, 0.359825_wp, 0.357654_wp, &
      & 0.355482_wp, 0.353311_wp, 0.351139_wp, 0.348968_wp, 0.346797_wp, & ! 100
      & 0.344625_wp, 0.342454_wp, 0.340282_wp]

   real(wp), parameter :: repa_angewChem2020(103) = [&
      & 2.639785_wp, 3.575012_wp, 0.732142_wp, 1.159621_wp, 1.561585_wp, &
      & 1.762895_wp, 2.173015_wp, 2.262269_wp, 2.511112_wp, 3.577220_wp, &
      & 0.338845_wp, 0.693023_wp, 0.678792_wp, 0.804784_wp, 1.012178_wp, &
      & 1.103469_wp, 1.209798_wp, 1.167791_wp, 0.326946_wp, 0.595242_wp, &
      & 1.447860_wp, 1.414501_wp, 1.381142_wp, 1.347783_wp, 1.314424_wp, &
      & 1.281065_wp, 1.247706_wp, 1.214347_wp, 1.180988_wp, 1.147629_wp, &
      & 0.700620_wp, 0.721266_wp, 0.741789_wp, 0.857434_wp, 0.875583_wp, &
      & 0.835876_wp, 0.290625_wp, 0.554446_wp, 0.623980_wp, 0.696005_wp, &
      & 0.768030_wp, 0.840055_wp, 0.912081_wp, 0.984106_wp, 1.056131_wp, &
      & 1.128156_wp, 1.200181_wp, 1.272206_wp, 0.478807_wp, 0.479759_wp, &
      & 0.579840_wp, 0.595241_wp, 0.644458_wp, 0.655289_wp, 0.574626_wp, &
      & 0.560506_wp, 0.682723_wp, 0.684824_wp, 0.686925_wp, 0.689026_wp, & ! 60
      & 0.691127_wp, 0.693228_wp, 0.695329_wp, 0.697430_wp, 0.699531_wp, &
      & 0.701631_wp, 0.703732_wp, 0.705833_wp, 0.707934_wp, 0.710035_wp, & ! 70
      & 0.712136_wp, 0.714237_wp, 0.745751_wp, 0.777265_wp, 0.808779_wp, &
      & 0.840294_wp, 0.871808_wp, 0.903322_wp, 0.934836_wp, 0.966350_wp, & ! 80
      & 0.467729_wp, 0.486102_wp, 0.559176_wp, 0.557520_wp, 0.563373_wp, &
      &  0.484713_wp,  0.574626_wp,  0.560506_wp,  0.682723_wp,  0.684824_wp, & ! 90
      &  0.686925_wp,  0.689026_wp,  0.691127_wp,  0.693228_wp,  0.695329_wp, &
      &  0.697430_wp,  0.699531_wp,  0.701631_wp,  0.703732_wp,  0.705833_wp, &
      &  0.707934_wp,  0.710035_wp,  0.712136_wp &
      &]

   real(wp), parameter :: repan_angewChem2020(103) = [&
      & 1.071395_wp, 1.072699_wp, 1.416847_wp, 1.156187_wp, 0.682382_wp, &
      & 0.556380_wp, 0.746785_wp, 0.847242_wp, 0.997252_wp, 0.873051_wp, &
      & 0.322503_wp, 0.415554_wp, 0.423946_wp, 0.415776_wp, 0.486773_wp, &
      & 0.494532_wp, 0.705274_wp, 0.706778_wp, 0.311178_wp, 0.399439_wp, &
      & 0.440983_wp, 0.475582_wp, 0.510180_wp, 0.544779_wp, 0.579377_wp, &
      & 0.613976_wp, 0.648574_wp, 0.683173_wp, 0.717772_wp, 0.752370_wp, &
      & 0.429944_wp, 0.420053_wp, 0.384743_wp, 0.443762_wp, 0.538680_wp, &
      & 0.472196_wp, 0.423850_wp, 0.385815_wp, 0.249213_wp, 0.285604_wp, &
      & 0.321995_wp, 0.358387_wp, 0.394778_wp, 0.431169_wp, 0.467560_wp, &
      & 0.503952_wp, 0.540343_wp, 0.576734_wp, 0.333476_wp, 0.348734_wp, &
      & 0.358194_wp, 0.351053_wp, 0.404536_wp, 0.389847_wp, 0.302575_wp, &
      & 0.163290_wp, 0.187645_wp, 0.190821_wp, 0.193998_wp, 0.197174_wp, &
      & 0.200351_wp, 0.203527_wp, 0.206703_wp, 0.209880_wp, 0.213056_wp, &
      & 0.216233_wp, 0.219409_wp, 0.222585_wp, 0.225762_wp, 0.228938_wp, &
      & 0.232115_wp, 0.235291_wp, 0.282937_wp, 0.330583_wp, 0.378229_wp, &
      & 0.425876_wp, 0.473522_wp, 0.521168_wp, 0.568814_wp, 0.616460_wp, &
      & 0.242521_wp, 0.293680_wp, 0.320931_wp, 0.322666_wp, 0.333641_wp, &
      &  0.434163_wp,  0.302575_wp,  0.163290_wp,  0.187645_wp,  0.190821_wp, &
      &  0.193998_wp,  0.197174_wp,  0.200351_wp,  0.203527_wp,  0.206703_wp, &
      &  0.209880_wp,  0.213056_wp,  0.216233_wp,  0.219409_wp,  0.222585_wp, &
      &  0.225762_wp,  0.228938_wp,  0.232115_wp &
      &]

   real(wp), parameter :: angl_angewChem2020(103) = [&
      & 1.661808_wp, 0.300000_wp, 0.018158_wp, 0.029224_wp, 0.572683_wp, &
      & 0.771055_wp, 1.053577_wp, 2.159889_wp, 1.525582_wp, 0.400000_wp, &
      & 0.041070_wp, 0.028889_wp, 0.086910_wp, 0.494456_wp, 0.409204_wp, &
      & 0.864972_wp, 1.986025_wp, 0.491537_wp, 0.050168_wp, 0.072745_wp, &
      & 0.378334_wp, 0.346400_wp, 0.314466_wp, 0.282532_wp, 0.250598_wp, &
      & 0.218663_wp, 0.186729_wp, 0.154795_wp, 0.122861_wp, 0.090927_wp, &
      & 0.140458_wp, 0.653971_wp, 0.528465_wp, 0.420379_wp, 2.243492_wp, &
      & 0.400000_wp, 0.035341_wp, 0.022704_wp, 0.195060_wp, 0.188476_wp, &
      & 0.181892_wp, 0.175308_wp, 0.168724_wp, 0.162139_wp, 0.155555_wp, &
      & 0.148971_wp, 0.142387_wp, 0.135803_wp, 0.169779_wp, 0.265730_wp, &
      & 0.505495_wp, 0.398254_wp, 2.640752_wp, 0.568026_wp, 0.032198_wp, &
      & 0.036663_wp, 0.281449_wp, 0.280526_wp, 0.279603_wp, 0.278680_wp, &
      & 0.277757_wp, 0.276834_wp, 0.275911_wp, 0.274988_wp, 0.274065_wp, &
      & 0.273142_wp, 0.272219_wp, 0.271296_wp, 0.270373_wp, 0.269450_wp, &
      & 0.268528_wp, 0.267605_wp, 0.253760_wp, 0.239916_wp, 0.226071_wp, &
      & 0.212227_wp, 0.198382_wp, 0.184538_wp, 0.170693_wp, 0.156849_wp, &
      & 0.104547_wp, 0.313474_wp, 0.220185_wp, 0.415042_wp, 1.259822_wp, &
      &  0.400000_wp,  0.032198_wp,  0.036663_wp,  0.281449_wp,  0.280526_wp, &
      &  0.279603_wp,  0.278680_wp,  0.277757_wp,  0.276834_wp,  0.275911_wp, &
      &  0.274988_wp,  0.274065_wp,  0.273142_wp,  0.272219_wp,  0.271296_wp, &
      &  0.270373_wp,  0.269450_wp,  0.268528_wp &
      &]

   real(wp), parameter :: angl2_angewChem2020(103) = [&
      & 0.624197_wp, 0.600000_wp, 0.050000_wp, 0.101579_wp, 0.180347_wp, &
      & 0.755851_wp, 0.761551_wp, 0.813653_wp, 0.791274_wp, 0.400000_wp, &
      & 0.000000_wp, 0.022706_wp, 0.100000_wp, 0.338514_wp, 0.453023_wp, &
      & 0.603722_wp, 1.051121_wp, 0.547904_wp, 0.000000_wp, 0.059059_wp, &
      & 0.117040_wp, 0.118438_wp, 0.119836_wp, 0.121234_wp, 0.122632_wp, &
      & 0.124031_wp, 0.125429_wp, 0.126827_wp, 0.128225_wp, 0.129623_wp, &
      & 0.206779_wp, 0.466678_wp, 0.496442_wp, 0.617321_wp, 0.409933_wp, &
      & 0.400000_wp, 0.000000_wp, 0.000000_wp, 0.119120_wp, 0.118163_wp, &
      & 0.117206_wp, 0.116249_wp, 0.115292_wp, 0.114336_wp, 0.113379_wp, &
      & 0.112422_wp, 0.111465_wp, 0.110508_wp, 0.149917_wp, 0.308383_wp, &
      & 0.527398_wp, 0.577885_wp, 0.320371_wp, 0.568026_wp, 0.000000_wp, &
      & 0.000000_wp, 0.078710_wp, 0.079266_wp, 0.079822_wp, 0.080379_wp, &
      & 0.080935_wp, 0.081491_wp, 0.082047_wp, 0.082603_wp, 0.083159_wp, &
      & 0.083716_wp, 0.084272_wp, 0.084828_wp, 0.085384_wp, 0.085940_wp, &
      & 0.086496_wp, 0.087053_wp, 0.095395_wp, 0.103738_wp, 0.112081_wp, &
      & 0.120423_wp, 0.128766_wp, 0.137109_wp, 0.145451_wp, 0.153794_wp, &
      & 0.323570_wp, 0.233450_wp, 0.268137_wp, 0.307481_wp, 0.316447_wp, &
      &  0.400000_wp,  0.000000_wp,  0.000000_wp,  0.078710_wp,  0.079266_wp, &
      &  0.079822_wp,  0.080379_wp,  0.080935_wp,  0.081491_wp,  0.082047_wp, &
      &  0.082603_wp,  0.083159_wp,  0.083716_wp,  0.084272_wp,  0.084828_wp, &
      &  0.085384_wp,  0.085940_wp,  0.086496_wp &
      &]

   real(wp), parameter :: tors_angewChem2020(103) = [&
      & 0.100000_wp, 0.100000_wp, 0.100000_wp, 0.000000_wp, 0.121170_wp, &
      & 0.260028_wp, 0.222546_wp, 0.250620_wp, 0.256328_wp, 0.400000_wp, &
      & 0.115000_wp, 0.000000_wp, 0.103731_wp, 0.069103_wp, 0.104280_wp, &
      & 0.226131_wp, 0.300000_wp, 0.400000_wp, 0.124098_wp, 0.000000_wp, &
      & 0.105007_wp, 0.107267_wp, 0.109526_wp, 0.111786_wp, 0.114046_wp, &
      & 0.116305_wp, 0.118565_wp, 0.120825_wp, 0.123084_wp, 0.125344_wp, &
      & 0.395722_wp, 0.349100_wp, 0.147808_wp, 0.259811_wp, 0.400000_wp, &
      & 0.400000_wp, 0.112206_wp,-0.004549_wp, 0.198713_wp, 0.179472_wp, &
      & 0.160232_wp, 0.140991_wp, 0.121751_wp, 0.102510_wp, 0.083270_wp, &
      & 0.064029_wp, 0.044789_wp, 0.025548_wp, 0.202245_wp, 0.278223_wp, &
      & 0.280596_wp, 0.229057_wp, 0.300000_wp, 0.423199_wp, 0.090741_wp, &
      & 0.076783_wp, 0.310896_wp, 0.309131_wp, 0.307367_wp, 0.305602_wp, &
      & 0.303838_wp, 0.302073_wp, 0.300309_wp, 0.298544_wp, 0.296779_wp, &
      & 0.295015_wp, 0.293250_wp, 0.291486_wp, 0.289721_wp, 0.287957_wp, & ! 70
      & 0.286192_wp, 0.284427_wp, 0.257959_wp, 0.231490_wp, 0.205022_wp, &
      & 0.178553_wp, 0.152085_wp, 0.125616_wp, 0.099147_wp, 0.072679_wp, & ! 80
      & 0.203077_wp, 0.169346_wp, 0.090568_wp, 0.144762_wp, 0.231884_wp, &
      &  0.400000_wp,  0.090741_wp,  0.076783_wp,  0.310896_wp,  0.309131_wp, & ! 90
      &  0.307367_wp,  0.305602_wp,  0.303838_wp,  0.302073_wp,  0.300309_wp, &
      &  0.298544_wp,  0.296779_wp,  0.295015_wp,  0.293250_wp,  0.291486_wp, & ! 100
      &  0.289721_wp,  0.287957_wp,  0.286192_wp &
      &]

   real(wp), parameter :: tors2_angewChem2020(103) = [&
      & 1.618678_wp, 1.000000_wp, 0.064677_wp, 0.000000_wp, 0.965814_wp, &
      & 1.324709_wp, 1.079334_wp, 1.478599_wp, 0.304844_wp, 0.500000_wp, &
      & 0.029210_wp, 0.000000_wp, 0.417423_wp, 0.334275_wp, 0.817008_wp, &
      & 0.922181_wp, 0.356367_wp, 0.684881_wp, 0.029210_wp, 0.000000_wp, &
      & 0.035902_wp, 0.090952_wp, 0.146002_wp, 0.201052_wp, 0.256103_wp, &
      & 0.311153_wp, 0.366203_wp, 0.421253_wp, 0.476303_wp, 0.531353_wp, &
      & 0.482963_wp, 1.415893_wp, 1.146581_wp, 1.338448_wp, 0.376801_wp, &
      & 0.500000_wp, 0.027213_wp,-0.004549_wp, 0.003820_wp, 0.093011_wp, &
      & 0.182202_wp, 0.271393_wp, 0.360584_wp, 0.449775_wp, 0.538965_wp, &
      & 0.628156_wp, 0.717347_wp, 0.806538_wp, 0.077000_wp, 0.185110_wp, &
      & 0.432427_wp, 0.887811_wp, 0.267721_wp, 0.571662_wp, 0.000000_wp, &
      & 0.000000_wp, 0.122336_wp, 0.131176_wp, 0.140015_wp, 0.148855_wp, &
      & 0.157695_wp, 0.166534_wp, 0.175374_wp, 0.184214_wp, 0.193053_wp, &
      & 0.201893_wp, 0.210733_wp, 0.219572_wp, 0.228412_wp, 0.237252_wp, &
      & 0.246091_wp, 0.254931_wp, 0.387526_wp, 0.520121_wp, 0.652716_wp, &
      & 0.785311_wp, 0.917906_wp, 1.050500_wp, 1.183095_wp, 1.315690_wp, &
      & 0.219729_wp, 0.344830_wp, 0.331862_wp, 0.767979_wp, 0.536799_wp, &
      &  0.500000_wp,  0.000000_wp,  0.000000_wp,  0.122336_wp,  0.131176_wp, &
      &  0.140015_wp,  0.148855_wp,  0.157695_wp,  0.166534_wp,  0.175374_wp, &
      &  0.184214_wp,  0.193053_wp,  0.201893_wp,  0.210733_wp,  0.219572_wp, &
      &  0.228412_wp,  0.237252_wp,  0.246091_wp &
      &]

!----------------------------------------------------------------------------------------

   !========================================================================
   ! DATA
   !------------------------------------------------------------------------

   !Pauling EN
   ! Period 7: Pauling EN from https://en.wikipedia.org/wiki/Template:Periodic_table_(electronegativity_by_Pauling_scale)
   real(wp), parameter :: en(1:103) = [&
  &         2.200,3.000,0.980,1.570,2.040,2.550,3.040,3.440,3.980 &
  &        ,4.500,0.930,1.310,1.610,1.900,2.190,2.580,3.160,3.500 &
  &        ,0.820,1.000,1.360,1.540,1.630,1.660,1.550,1.830,1.880 &
  &        ,1.910,1.900,1.650,1.810,2.010,2.180,2.550,2.960,3.000 &
  &        ,0.820,0.950,1.220,1.330,1.600,2.160,1.900,2.200,2.280 &
  &        ,2.200,1.930,1.690,1.780,1.960,2.050,2.100,2.660,2.600 &
  &,0.79,0.89,1.10,1.12,1.13,1.14,1.15,1.17,1.18,1.20,1.21,1.22 &
  &,1.23,1.24,1.25,1.26,1.27,1.3,1.5,1.7,1.9,2.1,2.2,2.2,2.2 &   ! value of W-Au modified
  &,2.00,1.62,2.33,2.02,2.0,2.2,2.2 &
  &,0.79,0.9,1.1,1.3,1.5,1.38,1.36,1.28,1.13,1.28,1.3,1.3,1.3,1.3,1.3,1.3,1.3 & !Fr-Lr
  &]

   ! COVALENT RADII, used only in EHB, XB: abhgfnff_eg1,abhgfnff_eg2new,
   !   abhgfnff_eg3, and rbxgfnff_eg
   ! based on "Atomic Radii of the Elements," M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
   ! in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
   ! edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
   ! corrected Nov. 17, 2010 for the 92nd edition.
   ! Period 7: Same book, 95th edition
   real(wp), parameter :: rad(1:103) = [&
  &0.32D0,0.37D0,1.30D0,0.99D0,0.84D0,0.75D0,0.71D0,0.64D0,0.60D0,&
  &0.62D0,1.60D0,1.40D0,1.24D0,1.14D0,1.09D0,1.04D0,1.00D0,1.01D0,&
  &2.00D0,1.74D0,1.59D0,1.48D0,1.44D0,1.30D0,1.29D0,1.24D0,1.18D0,&
  &1.17D0,1.22D0,1.20D0,1.23D0,1.20D0,1.20D0,1.18D0,1.17D0,1.16D0,&
  &2.15D0,1.90D0,1.76D0,1.64D0,1.56D0,1.46D0,1.38D0,1.36D0,1.34D0,&
  &1.30D0,1.36D0,1.40D0,1.42D0,1.40D0,1.40D0,1.37D0,1.36D0,1.36D0,&
  &2.38D0,2.06D0,1.94D0,1.84D0,1.90D0,1.88D0,1.86D0,1.85D0,1.83D0,&
  &1.82D0,1.81D0,1.80D0,1.79D0,1.77D0,1.77D0,1.78D0,1.74D0,1.64D0,&
  &1.58D0,1.50D0,1.41D0,1.36D0,1.32D0,1.30D0,1.30D0,1.32D0,1.44D0,&
  &1.45D0,1.50D0,1.42D0,1.48D0,1.46D0, &
  &2.42D0,2.11D0,2.01D0,1.90D0,1.84D0,1.83D0,1.80D0,1.80D0,1.73D0,& !Fr-Am
  &1.68D0,1.68D0,1.68D0,1.65D0,1.67D0,1.73D0,1.76D0,1.61D0 & !Cm-Lr
  &]

   integer, parameter :: metal(103) = (/ &
  &0,                                                                0,&!He
  &1,1,                                               0, 0, 0, 0, 0, 0,&!Ne
  &1,1,                                               1, 0, 0, 0, 0, 0,&!Ar
  &1,1,2,                2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 0, 0, 0, 0, 0,&!Kr
  &1,2,2,                2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 1, 0, 0, 0, 0,&!Xe
  &1,2,2, 2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
  &                      2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 1, 1, 1, 0, 0,&!Rn
  &1,2,2, 2,2,2,2,2,2,2,2,2,2,2,2,2,2 & !Fr-Lr
  &/)
  ! At is NOT a metal, Po is borderline but slightly better as metal
   integer, private, parameter :: group(103) = (/ &
  &1,                                                                   8,&!He
  &1,2,                                                  3, 4, 5, 6, 7, 8,&!Ne
  &1,2,                                                  3, 4, 5, 6, 7, 8,&!Ar
  &1,2,-3,                 -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8,&!Kr
  &1,2,-3,                 -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8,&!Xe
  &1,2,-3,  -3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3, &
  &                        -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8,&!Rn
  &1,2,-13, -13,-13,-13,-13,-13,-13,-13,-13,-13,-13,-13,-13,-13,-13 & !Fr-Lr
  &/)
   integer, parameter :: normcn(103) = (/ &  ! only for non metals well defined
  &1,                                                                0,&!He
  &4,4,                                               4, 4, 4, 2, 1, 0,&!Ne
  &4,4,                                               4, 4, 4, 2, 1, 0,&!Ar
  &4,4,4,                4, 6, 6, 6, 6, 6, 6, 4, 4,   4, 4, 4, 4, 1, 0,&!Kr
  &4,4,4,                4, 6, 6, 6, 6, 6, 6, 4, 4,   4, 4, 4, 4, 1, 0,&!Xe
  &4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4, &
  &                      4, 6, 6, 6, 6, 6, 6, 6, 4,   4, 4, 4, 4, 1, 0,&!Rn
  &4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4 & !Fr-Lr
  &/)
   real(wp), parameter :: repz(103) =   (/ &
  &1.,                                                                    2.,&!He
  &1.,2.,                                                  3.,4.,5.,6.,7.,8.,&!Ne
  &1.,2.,                                                  3.,4.,5.,6.,7.,8.,&!Ar
  &1.,2.,3.,                 4.,5.,6.,7.,8.,9.,10.,11.,12.,3.,4.,5.,6.,7.,8.,&!Kr
  &1.,2.,3.,                 4.,5.,6.,7.,8.,9.,10.,11.,12.,3.,4.,5.,6.,7.,8.,&!Xe
  &1.,2.,3.,  3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3., &
  &                          4.,5.,6.,7.,8.,9.,10.,11.,12.,3.,4.,5.,6.,7.,8.,& !Rn
  &1.,2.,3.,  3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3. & !Fr-Lr
  &/)

   contains

   subroutine gfnff_set_param(n, gen, param)
     use xtb_mctc_accuracy, only : wp 
     use xtb_param_sqrtzr4r2, only : sqrtZr4r2
     implicit none
!    Dummy                      ,
     integer,intent(in)  :: n
     type(TGFFGenerator), intent(out) :: gen
     type(TGFFData), intent(inout) :: param
!    Stack
     integer   :: i,j,k
     real(wp)  :: dum

     call newGFNFFGenerator(gen)

     param%cnmax   = 4.4         ! max. CN considered ie all larger values smoothly set to this val
     param%atcuta  = 0.595d0     ! angle damping
     param%atcutt  = 0.505d0     ! torsion angle damping
     param%atcuta_nci  = 0.395d0 ! nci angle damping in HB term
     param%atcutt_nci  = 0.305d0 ! nci torsion angle damping in HB term
     param%repscalb= 1.7583      ! bonded rep. scaling
     param%repscaln= 0.4270_wp   ! non-bonded rep. scaling
     param%hbacut   =49.         ! HB angle cut-off
     param%hbscut   =22.         ! HB SR     "   "
     param%xbacut   =70.         ! same for XB
     param%xbscut   = 5.         !
     param%hbsf     = 1.         ! charge dep.
     param%hbst     =15.         ! 10 is better for S22, 20 better for HCN2 and S30L
     param%xbsf     =0.03        !
     param%xbst     =15.         !
     param%hbalp    = 6.         ! damp
     param%hblongcut=85.         ! values larger than 85 yield large RMSDs for P26
     param%hblongcut_xb=70.      ! values larger than 70 yield large MAD for HAL28
     param%hbabmix  =0.80        !
     param%hbnbcut  =11.20       !
     param%tors_hb   =0.94       ! torsion potential shift in HB term
     param%bend_hb   =0.20       ! bending potential shift in HB term
     param%vbond_scale=0.9       ! vbond(2) scaling for CN(H) = 1
     param%xhaci_globabh=0.268   ! A-H...B gen. scaling
     param%xhaci_coh=0.350       ! A-H...O=C gen. scaling
     param%xhaci_glob=1.50       ! acidity
     param%xhbas(:) = 0.0_wp
     param%xhbas( 6)=0.80d0      ! basicities (XB and HB), i.e., B...X-A or B...H..A
     param%xhbas( 7)=1.68d0
     param%xhbas( 8)=0.67d0
     param%xhbas( 9)=0.52d0
     param%xhbas(14)=4.0d0
     param%xhbas(15)=3.5d0
     param%xhbas(16)=2.0d0
     param%xhbas(17)=1.5d0
     param%xhbas(35)=1.5d0
     param%xhbas(53)=1.9d0
     param%xhbas(33)=param%xhbas(15)
     param%xhbas(34)=param%xhbas(16)
     param%xhbas(51)=param%xhbas(15)
     param%xhbas(52)=param%xhbas(16)
     param%xhaci(:) = 0.0_wp
     param%xhaci( 6)=0.75               ! HB acidities, a bit weaker for CH
     param%xhaci( 7)=param%xhaci_glob+0.1
     param%xhaci( 8)=param%xhaci_glob
     param%xhaci( 9)=param%xhaci_glob
     param%xhaci(15)=param%xhaci_glob
     param%xhaci(16)=param%xhaci_glob
     param%xhaci(17)=param%xhaci_glob+1.0
     param%xhaci(35)=param%xhaci_glob+1.0
     param%xhaci(53)=param%xhaci_glob+1.0
     param%xbaci(:) = 0.0_wp
     param%xbaci(15)=1.0d0              ! XB acidities
     param%xbaci(16)=1.0d0
     param%xbaci(17)=0.5d0
     param%xbaci(33)=1.2d0
     param%xbaci(34)=1.2d0
     param%xbaci(35)=0.9d0
     param%xbaci(51)=1.2d0
     param%xbaci(52)=1.2d0
     param%xbaci(53)=1.2d0

!    3B bond prefactors and D3 stuff
     k=0
     do i=1,103
        dum=dble(i)
        param%zb3atm(i)=-dum*gen%batmscal**(1.d0/3.d0)  ! inlcude pre-factor
        do j=1,i
           k=k+1
           dum=sqrtZr4r2(i)*sqrtZr4r2(j)*3.0d0
           param%d3r0(k)=(gen%d3a1*dsqrt(dum)+gen%d3a2)**2   ! save R0^2 for efficiency reasons
        enddo
     enddo
     param%zb3atm(1)=-0.25d0*gen%batmscal**(1.d0/3.d0) ! slightly better than 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! numerical precision settings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end subroutine gfnff_set_param


   subroutine gfnff_thresholds(accuracy, dispthr, cnthr, repthr, hbthr1, hbthr2)
      real(wp), intent(in) :: accuracy
      real(wp), intent(out) :: dispthr
      real(wp), intent(out) :: cnthr
      real(wp), intent(out) :: repthr
      real(wp), intent(out) :: hbthr1
      real(wp), intent(out) :: hbthr2
      dispthr=1500.0d0-log10(accuracy)*1000.d0
      cnthr  = 100.0d0-log10(accuracy)* 50.0d0
      repthr = 400.0d0-log10(accuracy)*100.0d0
      hbthr1 = 200.0d0-log10(accuracy)* 50.0d0
      hbthr2 = 400.0d0-log10(accuracy)* 50.0d0
   end subroutine gfnff_thresholds


   subroutine gfnff_load_param(version, param, exist)
     use xtb_mctc_accuracy, only : wp 
     use xtb_param_covalentRadD3, only : covalentRadD3
     implicit none
     integer, intent(in) :: version
     type(TGFFData), intent(out) :: param
     logical, intent(out) :: exist

     exist = .false.

     call init(param, 103)

     param%en(:) = en
     param%rad(:) = rad
     param%rcov(:) = covalentRadD3(1:103)
     param%metal(:) = metal
     param%group(:) = group
     param%normcn(:) = normcn
     param%repz(:) = repz

     select case(version)
     case(gffVersion%angewChem2020, gffVersion%angewChem2020_1, gffVersion%angewChem2020_2, &
                     & gffVersion%harmonic2020, gffVersion%mcgfnff2023)
        call loadGFNFFAngewChem2020(param)
        exist = .true.
     end select

   end subroutine gfnff_load_param


   subroutine loadGFNFFAngewChem2020(param)
     type(TGFFData), intent(inout) :: param
     param%chi(:) = chi_angewChem2020
     param%gam(:) = gam_angewChem2020
     param%cnf(:) = cnf_angewChem2020
     param%alp(:) = alp_angewChem2020
     param%bond(:) = bond_angewChem2020
     param%repa(:) = repa_angewChem2020
     param%repan(:) = repan_angewChem2020
     param%angl(:) = angl_angewChem2020
     param%angl2(:) = angl2_angewChem2020
     param%tors(:) = tors_angewChem2020
     param%tors2(:) = tors2_angewChem2020
   end subroutine loadGFNFFAngewChem2020


   subroutine gfnff_read_param(iunit, param)
     use xtb_mctc_accuracy, only : wp 
     use xtb_param_covalentRadD3, only : covalentRadD3
     implicit none
!    Dummy
     integer,intent(in)  :: iunit
     type(TGFFData), intent(out) :: param
!    Stack
     integer  :: i,nn
     real(wp) :: xx(20)
     character(len=256) :: atmp

     call init(param, 103)

     param%en(:) = en
     param%rad(:) = rad
     param%rcov(:) = covalentRadD3(1:103)
     param%metal(:) = metal
     param%group(:) = group
     param%normcn(:) = normcn
     param%repz(:) = repz

     do i=1,103
         read(iunit,'(a)')atmp
         call readl(atmp,xx,nn)
         param%chi(i)=xx(2)
         param%gam(i)=xx(3)
         param%cnf(i)=xx(4)
         param%alp(i)=xx(5)
         param%bond(i)=xx(6)
         param%repa(i)=xx(7)
         param%repan(i)=xx(8)
         param%angl(i)=xx(9)
         param%angl2(i)=xx(10)
         param%tors(i)=xx(11)
         param%tors2(i)=xx(12)
      enddo

   end subroutine gfnff_read_param

   subroutine gfnff_param_alloc(topo, neigh, n)
     use xtb_mctc_accuracy, only : wp 
     use xtb_gfnff_neighbor
     implicit none
!    Dummy
     type(TGFFTopology), intent(inout) :: topo
     type(TNeigh), intent(inout) :: neigh
     integer,intent(in) :: n

     if (.not.allocated(topo%alphanb)) allocate( topo%alphanb(n,n,neigh%numctr+1), source = 0.0d0 )
     if (.not.allocated(topo%chieeq)) allocate( topo%chieeq(n), source = 0.0d0 )
     if (.not.allocated(topo%gameeq)) allocate( topo%gameeq(n), source = 0.0d0 )
     if (.not.allocated(topo%alpeeq)) allocate( topo%alpeeq(n), source = 0.0d0 )
     if (.not.allocated(topo%qa)) allocate( topo%qa(n), source = 0.0d0 )
     if (.not.allocated(topo%zetac6)) allocate( topo%zetac6(n*(n+1)/2), source = 0.0d0 )
     if (.not.allocated(topo%xyze0)) allocate( topo%xyze0(3,n), source = 0.0d0 )
     if (.not.allocated(topo%b3list)) allocate( topo%b3list(5,1000*n), source = 0 )
     if (.not.allocated(topo%sTorsl)) allocate( topo%sTorsl(6,topo%nstors), source = 0 )
     if (.not.allocated(topo%fraglist)) allocate( topo%fraglist(n), source = 0 )
     if (.not.allocated(topo%qfrag)) allocate( topo%qfrag(n), source = 0.0d0 )
     if (.not.allocated(topo%hbatHl)) allocate( topo%hbatHl(2,n), source = 0 )
     if (.not.allocated(topo%hbbas)) allocate( topo%hbbas(n), source = 0.0d0 )
     if (.not.allocated(topo%hbaci)) allocate( topo%hbaci(n), source = 0.0d0 )
     if (.not.allocated(topo%hbatABl)) allocate( topo%hbatABl(2,n*(n+1)/2), source = 0 )
     if (.not.allocated(topo%xbatABl)) allocate( topo%xbatABl(5,topo%natxbAB), source = 0 )

     if (.not.allocated(neigh%blist)) allocate( neigh%blist(3,neigh%nbond), source = 0 )
     if (.not.allocated(neigh%nr_hb)) allocate( neigh%nr_hb(neigh%nbond), source = 0 )
     if (.not.allocated(topo%bond_hb_AH)) allocate( topo%bond_hb_AH(4,topo%bond_hb_nr), source = 0 )
     if (.not.allocated(topo%bond_hb_B)) allocate( topo%bond_hb_B(2,topo%b_max,topo%bond_hb_nr), source = 0 )
     if (.not.allocated(topo%bond_hb_Bn)) allocate( topo%bond_hb_Bn(topo%bond_hb_nr), source = 0 )
     if (.not.allocated(topo%alist)) allocate( topo%alist(5,topo%nangl_alloc), source = 0 )
     if (.not.allocated(topo%tlist)) allocate( topo%tlist(8,topo%ntors_alloc), source = 0 )
     if (.not.allocated(topo%vbond)) allocate( topo%vbond(3,topo%nbond_vbond), source = 0.0d0 )
     if (.not.allocated(topo%vangl)) allocate( topo%vangl(2,topo%nangl_alloc), source = 0.0d0 )
     if (.not.allocated(topo%vtors)) allocate( topo%vtors(2,topo%ntors_alloc), source = 0.0d0 )

   end subroutine gfnff_param_alloc

   subroutine gfnff_param_dealloc(topo)
     use xtb_mctc_accuracy, only : wp 
     implicit none
     type(TGFFTopology), intent(inout) :: topo
!    Dummy

     if (allocated(topo%alphanb)) deallocate( topo%alphanb )
     if (allocated(topo%chieeq)) deallocate( topo%chieeq )
     if (allocated(topo%gameeq)) deallocate( topo%gameeq )
     if (allocated(topo%alpeeq)) deallocate( topo%alpeeq )
     if (allocated(topo%qa)) deallocate( topo%qa )
     if (allocated(topo%zetac6)) deallocate( topo%zetac6 )
     if (allocated(topo%xyze0)) deallocate( topo%xyze0 )
     if (allocated(topo%b3list)) deallocate( topo%b3list )
     if (allocated(topo%fraglist)) deallocate( topo%fraglist )
     if (allocated(topo%qfrag)) deallocate( topo%qfrag )
     if (allocated(topo%hbatHl)) deallocate( topo%hbatHl )
     if (allocated(topo%hbbas)) deallocate( topo%hbbas )
     if (allocated(topo%hbaci)) deallocate( topo%hbaci )
     if (allocated(topo%hbatABl)) deallocate( topo%hbatABl )
     if (allocated(topo%xbatABl)) deallocate( topo%xbatABl )

     if (allocated(topo%blist)) deallocate( topo%blist )
     if (allocated(topo%nr_hb)) deallocate( topo%nr_hb )
     if (allocated(topo%bond_hb_AH)) deallocate( topo%bond_hb_AH )
     if (allocated(topo%bond_hb_B)) deallocate( topo%bond_hb_B )
     if (allocated(topo%bond_hb_Bn)) deallocate( topo%bond_hb_Bn )
     if (allocated(topo%alist)) deallocate( topo%alist )
     if (allocated(topo%tlist)) deallocate( topo%tlist )
     if (allocated(topo%vbond)) deallocate( topo%vbond )
     if (allocated(topo%vangl)) deallocate( topo%vangl )
     if (allocated(topo%vtors)) deallocate( topo%vtors )

   end subroutine gfnff_param_dealloc


   subroutine newGFNFFGenerator(gen)
     type(TGFFGenerator), intent(out) :: gen

     gen%cnmax   = 4.4         ! max. CN considered ie all larger values smoothly set to this val
     gen%linthr  = 160.        ! when is an angle close to linear ? (GEODEP) for metals values closer to 170 (than to 160) are better
                           ! but this occurs e.g. for Sc in unclear situations. So make it save (160)
     gen%fcthr   = 1.d-3       ! skip torsion and bending if potential is small
     gen%tdist_thr=12.         ! R threshold in Angstroem for cov distance estimated used in apprx EEQ
                           ! the following two parameters are critical for topo setup
     gen%rthr     =1.25        ! important bond determination threshold
                           ! large values yield more 1.23
     gen%rthr2    =1.00        ! decrease if a metal is present, larger values yield smaller CN
     gen%rqshrink =0.23        ! change of R0 for topo with charge qa, larger values yield smaller CN for metals in particular
     gen%hqabthr  =0.01        ! H charge (qa) threshold for H in HB list 18
     gen%qabthr  =0.10        ! AB charge (qa) threshold for AB in HB list, avoids HBs with positive atoms,
     ! larger val. better for S30L but worse in PubChem RMSD checks
     gen%srb1    = 0.3731      ! bond params
     gen%srb2    = 0.3171      !
     gen%srb3    = 0.2538      !
     gen%qrepscal= 0.3480      ! change of non-bonded rep. with q(topo)
     gen%nrepscal=-0.1270      !   "    "      "       "   CN
     gen%hhfac   = 0.6290      ! HH repulsion
     gen%hh13rep = 1.4580      !
     gen%hh14rep = 0.7080      !
     gen%bstren(1)=1.00d0      ! single bond
     gen%bstren(2)=1.24d0      ! double bond
     gen%bstren(3)=1.98d0      ! triple bond
     gen%bstren(4)=1.22d0      ! hyperval bond
     gen%bstren(5)=1.00d0      ! M-X
     gen%bstren(6)=0.78d0      ! M eta
     gen%bstren(7)=3.40d0      ! M-M
     gen%bstren(8)=3.40d0      ! M-M
     gen%qfacBEN =-0.54        ! bend FC change with polarity
     gen%qfacTOR =12.0d0       ! torsion FC change with polarity
     gen%fr3     =0.3          ! tors FC 3-ring
     gen%fr4     =1.0          ! tors FC 4-ring
     gen%fr5     =1.5          ! tors FC 5-ring
     gen%fr6     =5.7          ! tors FC 6-ring
     gen%torsf(1)=1.00         ! single bond
     gen%torsf(2)=1.18         ! pi bond
     gen%torsf(3)=1.05         ! improper
     gen%torsf(5)=0.50         ! pi part improper
     gen%torsf(6)=-0.90        ! extra sp3 C
     gen%torsf(7)= 0.70        ! extra sp3 N
     gen%torsf(8)=-2.00        ! extra sp3 O
     gen%fbs1    =0.50         ! small bend corr.
     gen%batmscal=0.30d0       ! bonded ATM scal
     gen%mchishift=-0.09d0
     gen%rabshift    =-0.110   ! gen shift
     gen%rabshifth   =-0.050   ! XH
     gen%hyper_shift = 0.03    ! hypervalent
     gen%hshift3     = -0.11   ! heavy
     gen%hshift4     = -0.11   !
     gen%hshift5     = -0.06   !
     gen%metal1_shift= 0.2     ! group 1+2 metals
     gen%metal2_shift= 0.15    ! TM
     gen%metal3_shift= 0.05    ! main group metals
     gen%eta_shift   = 0.040   ! eta bonded
     gen%qfacbm(0)   =1.0d0    ! bond charge dep.gff_srcs += 'gff/gfnff_input.f90'

     gen%qfacbm(1:2) =-0.2d0   !
     gen%qfacbm(  3) =0.70d0   !
     gen%qfacbm(  4) =0.50d0   !
     gen%qfacbm0     =0.047    !
     gen%rfgoed1  =1.175       ! topo dist scaling
     gen%htriple = 1.45d0      ! decrease Hueckel off-diag for triple bonds because they are less well conjugated 1.4
     gen%hueckelp2=1.00d0      ! increase pot depth depending on P
     gen%hueckelp3=-0.24d0     ! diagonal element change with qa
     gen%hdiag(5) =-0.5d0      ! diagonal element relative to C
     gen%hdiag(6) =0.00d0      !
     gen%hdiag(7) =0.14d0      !
     gen%hdiag(8) =-0.38d0     !
     gen%hdiag(9) =-0.29d0     !
     gen%hdiag(16)=-0.30d0     !
     gen%hdiag(17)=-0.30d0     !
     gen%hoffdiag(5)=0.5d0     ! Huckel off-diag constants
     gen%hoffdiag(6)=1.00d0    !
     gen%hoffdiag(7)=0.66d0    !
     gen%hoffdiag(8)=1.10d0    !
     gen%hoffdiag(9)=0.23d0    !
     gen%hoffdiag(16)=0.60d0   !
     gen%hoffdiag(17)=1.00d0   !
     gen%hiter   =0.700d0      ! iteration mixing
     gen%hueckelp=0.340d0      ! diagonal qa dep.
     gen%bzref   =0.370d0      ! ref P value R shift
     gen%bzref2  =0.315d0      !  "  "  "    k stretch
     gen%pilpf   =0.530d0      ! 2el diag shift
     gen%maxhiter=5            ! the Hückel iterations can diverge so take only a few steps
     gen%d3a1    = 0.58d0      ! D3, s8 fixed = 2
     gen%d3a2    = 4.80d0
     gen%split0  =0.670d0      ! mixing of sp^n with sp^n-1
     gen%fringbo =0.020d0      ! str ring size dep.
     gen%aheavy3 =89.          ! three coord. heavy eq. angle
     gen%aheavy4 =100.         ! four   "       "    "    "
     gen%split1=1.0d0-gen%split0
     gen%bsmat = -999.
     gen%bsmat(0,0)=gen%bstren(1)
     gen%bsmat(3,0)=gen%bstren(1)
     gen%bsmat(3,3)=gen%bstren(1)
     gen%bsmat(2,2)=gen%bstren(2)
     gen%bsmat(1,1)=gen%bstren(3)
     gen%bsmat(1,0)=gen%split0*gen%bstren(1)+gen%split1*gen%bstren(3)
     gen%bsmat(3,1)=gen%split0*gen%bstren(1)+gen%split1*gen%bstren(3)
     gen%bsmat(2,1)=gen%split0*gen%bstren(2)+gen%split1*gen%bstren(3)
     gen%bsmat(2,0)=gen%split0*gen%bstren(1)+gen%split1*gen%bstren(2)
     gen%bsmat(3,2)=gen%split0*gen%bstren(1)+gen%split1*gen%bstren(2)
     gen%bstren(9)=0.5*(gen%bstren(7)+gen%bstren(8))

   end subroutine newGFNFFGenerator


end module xtb_gfnff_param
