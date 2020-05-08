module gff_param
   use iso_fortran_env, only : wp => real64, sp => real32
   implicit none
   private :: wp
   public

   !========================================================================
   ! all fixed parameters/potetnial parameters/thresholds/cut-offs 
   ! used in gfnff_ini and set in gfnff_set_param
   !------------------------------------------------------------------------

   !loicals
   logical  :: ini       !make a neghbour lsit or use existing one?
   logical  :: gff_print = .true. !shows timing of energy + gradient
   logical  :: make_chrg = .true. !generates new eeq chareges based on current geometry

   real(wp), private, parameter :: p_gff_chi(86) = [&
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
      & 0.798704_wp, 1.127797_wp, 1.127863_wp, 1.127928_wp, 1.127994_wp, &
      & 1.128059_wp, 1.128125_wp, 1.128190_wp, 1.128256_wp, 1.128322_wp, &
      & 1.128387_wp, 1.128453_wp, 1.128518_wp, 1.128584_wp, 1.128649_wp, &
      & 1.128715_wp, 1.128780_wp, 1.129764_wp, 1.130747_wp, 1.131731_wp, &
      & 1.132714_wp, 1.133698_wp, 1.134681_wp, 1.135665_wp, 1.136648_wp, &
      & 1.061832_wp, 1.053084_wp, 1.207830_wp, 1.236314_wp, 1.310129_wp, &
      & 1.157380_wp]
   real(wp), private, parameter :: p_gff_gam(86) = [&
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
      & 0.416264_wp,-0.011212_wp,-0.011046_wp,-0.010879_wp,-0.010713_wp, &
      &-0.010546_wp,-0.010380_wp,-0.010214_wp,-0.010047_wp,-0.009881_wp, &
      &-0.009714_wp,-0.009548_wp,-0.009382_wp,-0.009215_wp,-0.009049_wp, &
      &-0.008883_wp,-0.008716_wp,-0.006220_wp,-0.003724_wp,-0.001228_wp, &
      & 0.001267_wp, 0.003763_wp, 0.006259_wp, 0.008755_wp, 0.011251_wp, &
      & 0.020477_wp,-0.056566_wp, 0.051943_wp, 0.076708_wp, 0.000273_wp, &
      &-0.068929_wp]
   real(wp), private, parameter :: p_gff_cnf(86) = [&
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
      &-0.033006_wp, 0.058973_wp, 0.058595_wp, 0.058217_wp, 0.057838_wp, &
      & 0.057460_wp, 0.057082_wp, 0.056704_wp, 0.056326_wp, 0.055948_wp, &
      & 0.055569_wp, 0.055191_wp, 0.054813_wp, 0.054435_wp, 0.054057_wp, &
      & 0.053679_wp, 0.053300_wp, 0.047628_wp, 0.041955_wp, 0.036282_wp, &
      & 0.030610_wp, 0.024937_wp, 0.019264_wp, 0.013592_wp, 0.007919_wp, &
      & 0.006383_wp,-0.089155_wp,-0.001293_wp, 0.019269_wp, 0.074803_wp, &
      & 0.016657_wp]
   real(wp), private, parameter :: p_gff_alp(86) = [&
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
      & 0.770260_wp, 0.757372_wp, 0.757352_wp, 0.757332_wp, 0.757313_wp, &
      & 0.757293_wp, 0.757273_wp, 0.757253_wp, 0.757233_wp, 0.757213_wp, &
      & 0.757194_wp, 0.757174_wp, 0.757154_wp, 0.757134_wp, 0.757114_wp, &
      & 0.757095_wp, 0.757075_wp, 0.756778_wp, 0.756480_wp, 0.756183_wp, &
      & 0.755886_wp, 0.755589_wp, 0.755291_wp, 0.754994_wp, 0.754697_wp, &
      & 0.868029_wp, 1.684375_wp, 2.001040_wp, 2.067331_wp, 2.228923_wp, &
      & 1.874218_wp]
   real(wp), private, parameter :: p_gff_bond(86) = [&
      & 0.417997_wp, 0.258490_wp, 0.113608_wp, 0.195935_wp, 0.231217_wp, &
      & 0.385248_wp, 0.379257_wp, 0.339249_wp, 0.330706_wp, 0.120319_wp, &
      & 0.127255_wp, 0.173647_wp, 0.183796_wp, 0.273055_wp, 0.249044_wp, &
      & 0.290653_wp, 0.218744_wp, 0.034706_wp, 0.136353_wp, 0.192467_wp, &
      & 0.335860_wp, 0.314452_wp, 0.293044_wp, 0.271636_wp, 0.250228_wp, &
      & 0.228819_wp, 0.207411_wp, 0.186003_wp, 0.164595_wp, 0.143187_wp, &
      & 0.212434_wp, 0.210451_wp, 0.219870_wp, 0.224618_wp, 0.272206_wp, &
      & 0.147864_wp, 0.150000_wp, 0.150000_wp, 0.329501_wp, 0.309632_wp, &
      & 0.289763_wp, 0.269894_wp, 0.250025_wp, 0.230155_wp, 0.210286_wp, &
      & 0.190417_wp, 0.170548_wp, 0.150679_wp, 0.192977_wp, 0.173411_wp, &
      & 0.186907_wp, 0.192891_wp, 0.223202_wp, 0.172577_wp, 0.150000_wp, &
      & 0.150000_wp, 0.370682_wp, 0.368511_wp, 0.366339_wp, 0.364168_wp, &
      & 0.361996_wp, 0.359825_wp, 0.357654_wp, 0.355482_wp, 0.353311_wp, &
      & 0.351139_wp, 0.348968_wp, 0.346797_wp, 0.344625_wp, 0.342454_wp, &
      & 0.340282_wp, 0.338111_wp, 0.305540_wp, 0.272969_wp, 0.240398_wp, &
      & 0.207828_wp, 0.175257_wp, 0.142686_wp, 0.110115_wp, 0.077544_wp, &
      & 0.108597_wp, 0.148422_wp, 0.183731_wp, 0.192274_wp, 0.127706_wp, &
      & 0.086756_wp]
   real(wp), private, parameter :: p_gff_repa(86) = [&
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
      & 0.560506_wp, 0.682723_wp, 0.684824_wp, 0.686925_wp, 0.689026_wp, &
      & 0.691127_wp, 0.693228_wp, 0.695329_wp, 0.697430_wp, 0.699531_wp, &
      & 0.701631_wp, 0.703732_wp, 0.705833_wp, 0.707934_wp, 0.710035_wp, &
      & 0.712136_wp, 0.714237_wp, 0.745751_wp, 0.777265_wp, 0.808779_wp, &
      & 0.840294_wp, 0.871808_wp, 0.903322_wp, 0.934836_wp, 0.966350_wp, &
      & 0.467729_wp, 0.486102_wp, 0.559176_wp, 0.557520_wp, 0.563373_wp, &
      & 0.484713_wp]
   real(wp), private, parameter :: p_gff_repan(86) = [&
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
      & 0.434163_wp]
   real(wp), private, parameter :: p_gff_angl(86) = [&
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
      & 0.400000_wp]
   real(wp), private, parameter :: p_gff_angl2(86) = [&
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
      & 0.400000_wp]
   real(wp), private, parameter :: p_gff_tors(86) = [&
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
      & 0.295015_wp, 0.293250_wp, 0.291486_wp, 0.289721_wp, 0.287957_wp, &
      & 0.286192_wp, 0.284427_wp, 0.257959_wp, 0.231490_wp, 0.205022_wp, &
      & 0.178553_wp, 0.152085_wp, 0.125616_wp, 0.099147_wp, 0.072679_wp, &
      & 0.203077_wp, 0.169346_wp, 0.090568_wp, 0.144762_wp, 0.231884_wp, &
      & 0.400000_wp]
   real(wp), private, parameter :: p_gff_tors2(86) = [&
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
      & 0.500000_wp]

   !common variables which are used in energy-gradient routines
   real(wp) :: repscaln,repscalb      ! repulsion scaling
   real(wp) :: atcuta,atcutt          ! bend/tors angle damping
   real(wp) :: atcuta_nci,atcutt_nci  ! bend/tors nci angle damping for HB term
   real(wp) :: hbacut,hbscut          ! damping HB
   real(wp) :: xbacut,xbscut          ! damping XB
   real(wp) :: hbalp                  ! damping HB/XB
   real(wp) :: hblongcut,hblongcut_xb ! damping HB/XB
   real(wp) :: hbst,hbsf,xbst,xbsf    ! charge scaling HB/XB
   real(wp) :: xhaci_globabh          ! HB AH-B
   real(wp) :: xhaci_coh              ! HB AH-O=C
   real(wp) :: xhaci_glob             ! acidity
   real(wp) :: hbabmix                ! HB AH-B
   real(wp) :: hbnbcut                ! new parameter for neighbour angle
   real(wp) :: tors_hb                ! new parameter for HB NCI angle term
   real(wp) :: bend_hb                ! new parameter for HB NCI torsion term
   real(wp) :: vbond_scale            ! new parameter for FC scaling of bonds in HB
   real(wp) :: cnmax                  ! max CN cut-off
!----------------------------------------------------------------------------------------   
   real(wp) :: efield(3)              ! electric field components
   
   !Thresholds
   real(wp) linthr        ! when is an angle close to linear ? (GEODEP) for metals values closer to 170 (than to 160) are better
   real(wp) fcthr         ! skip torsion and bending if potential is small
   real(sp) tdist_thr     ! R threshold in Angstroem for cov distance estimated used in apprx EEQ 
   real(wp) rthr          ! important bond determination threshold, large values yield more 1.23
   real(wp) rthr2         ! decrease if a metal is present, larger values yield smaller CN
   real(wp) rqshrink      ! change of R0 for topo with charge qa, larger values yield smaller CN for metals in particular
   real(wp) hqabthr       ! H charge (qa) threshold for H in HB list 18
   real(wp) qabthr        ! AB charge (qa) threshold for AB in HB list, avoids HBs with positive atoms,
                          ! larger val. better for S30L but worse in PubChem RMSD checks
   !Parameter
   real(wp) srb1,srb2,srb3
   real(wp) qrepscal      ! change of non-bonded rep. with q(topo)
   real(wp) nrepscal      !   "    "      "       "   CN
   real(wp) hhfac         ! HH repulsion
   real(wp) hh13rep
   real(wp) hh14rep
   real(wp) bstren(9)
   real(wp) qfacBEN       ! bend FC change with polarity 
   real(wp) qfacTOR       ! torsion FC change with polarity
   real(wp) fr3           ! tors FC 3-ring
   real(wp) fr4           ! tors FC 4-ring
   real(wp) fr5           ! tors FC 5-ring
   real(wp) fr6           ! tors FC 6-ring
   real(wp) torsf(8)      ! bonds     
   real(wp) fbs1          ! small bend corr.
   real(wp) batmscal      ! bonded ATM scal
   
   !Shifts
   real(wp) mchishift
   real(wp) rabshift     ! gen shift
   real(wp) rabshifth    ! XH
   real(wp) hyper_shift  ! hypervalent
   real(wp) hshift3      ! heavy
   real(wp) hshift4      !
   real(wp) hshift5      !
   real(wp) metal1_shift ! group 1+2 metals
   real(wp) metal2_shift ! TM
   real(wp) metal3_shift ! main group metals
   real(wp) eta_shift    ! eta bonded

   !Charge Param
   real(wp) qfacbm(0:4)  ! bond charge dependent
   real(wp) qfacbm0
   real(wp) rfgoed1      ! topo dist scaling

   !Hückel Param
   real(wp) htriple      ! decrease Hueckel off-diag for triple bonds because they are less well conjugated 1.4
   real(wp) hueckelp2    ! increase pot depth depending on P
   real(wp) hueckelp3    ! diagonal element change with qa  
   real(wp) hdiag(17)    ! diagonal element relative to C
   real(wp) hoffdiag(17) ! Huckel off-diag constants
   real(wp) hiter        ! iteration mixing
   real(wp) hueckelp     ! diagonal qa dep.
   real(wp) bzref        ! ref P value R shift
   real(wp) bzref2       !  "  "  "    k stretch
   real(wp) pilpf        ! 2el diag shift
   real(wp) maxhiter     ! the Hückel iterations can diverge so take only a few steps
   
   !D3 Param
   real(wp) d3a1         ! D3, s8 fixed = 2
   real(wp) d3a2    

   !Stuff
   real(wp) split0       ! mixing of sp^n with sp^n-1
   real(wp) split1       ! mixing of sp^n with sp^n-1
   real(wp) fringbo      ! str ring size dep.
   real(wp) aheavy3      ! three coord. heavy eq. angle
   real(wp) aheavy4      ! four   "       "    "    "
   real(wp) bsmat(0:3,0:3)

   !========================================================================
   ! parameters which are either rlisted below 
   ! or read in by gfnff_read_param
   !------------------------------------------------------------------------

   !general common stuff used in energy-gradient routines
   !> rep alpha bond
   real(wp) :: repa (86) = p_gff_repa,repan(86) = p_gff_repan
   real(wp) :: repz (86),zb3atm(86)             ! prefactor (Zval), 3atm bond term
   real(wp) :: xhaci(86),xhbas(86),xbaci(86)    ! HB/XB
   real(wp) :: rad(86)                          ! radii used for HB/XB damping and topology 
   real(wp) :: en (86)                          ! EN
   real(wp) :: chi(86) = p_gff_chi              ! EN dep. in EEQ.
   real(wp) :: gam(86) = p_gff_gam              ! EN dep. in EEQ.
   real(wp) :: cnf(86) = p_gff_cnf              ! EN dep. in EEQ.
   real(wp) :: alp(86) = p_gff_alp              ! EN dep. in EEQ.
   real(wp) :: bond(86) = p_gff_bond            ! Elem. bond param.
   real(wp) :: angl(86) = p_gff_angl            ! Elem. angular param.
   real(wp) :: angl2(86) = p_gff_angl2          ! Elem. angular param.
   real(wp) :: tors(86) = p_gff_tors            ! Elem. torsion param_alloc.
   real(wp) :: tors2(86) = p_gff_tors2          ! Elem. torsion param.
   real(wp) :: d3r0(86*87/2)                    ! BJ radii set in gnff_ini()
   
   !numerical precision cut-offs
   real(wp) :: cnthr,repthr,dispthr,hbthr1,hbthr2,accff

   integer  :: group(86),metal(86),normcn(86)   ! for assignment
   integer  :: ffmode
   

   !========================================================================
   ! parameters which are determined in gfnff_ini
   !------------------------------------------------------------------------

   !number of terms
   integer  :: nbond,nangl,ntors,nhb1,nhb2,nxb,nathbH,nathbAB,natxbAB,nbatm
   integer  :: nfrag
   integer  :: maxsystem   ! max. number of fragmentsfor hessian
   integer  :: bond_hb_nr  ! number of unique AH...B HB/bond terms
   integer  :: b_max      ! number of B atoms per unique AH bond

   !numbers that are rewritten, so must be stored for allocation
   integer  :: nbond_blist,nbond_vbond,nangl_alloc,ntors_alloc

   !file type read
   integer  :: read_file_type

   !lists
   integer,allocatable ::     nb(:,:)   ! neighbors nb(20,i) is the # neigbors
   integer,allocatable ::    bpair(:)   ! # of cov. between atoms
   integer,allocatable ::  blist(:,:)   ! bonded atoms
   integer,allocatable ::  alist(:,:)   ! angles
   integer,allocatable ::  tlist(:,:)   ! torsions
   integer,allocatable :: b3list(:,:)   ! bond atm   
   integer,allocatable :: hblist1(:,:)  ! HBs loose
   integer,allocatable :: hblist2(:,:)  ! HBs bonded
   integer,allocatable :: hblist3(:,:)  ! XBs
   !-----------------------------------------------
   integer,allocatable :: nr_hb(:)      ! Nr. of H bonds per O-H or N-H bond
   integer,allocatable :: bond_hb_AH(:,:) ! A, H atoms in bonds that are also part of HBs
   integer,allocatable :: bond_hb_B(:,:)  ! B atoms in bonds that are also part of HBs
   integer,allocatable :: bond_hb_Bn(:)   ! Nr. of B atoms for one AH bond pair
   !-----------------------------------------------
   integer,allocatable :: hbatABl(:,:)  ! AB atoms for HB
   integer,allocatable :: xbatABl(:,:)  ! AB atoms for XB
   integer,allocatable :: hbatHl (:)    ! H  atoms for HB
   integer,allocatable :: fraglist(:)   ! atoms in molecular fragments (for EEQ)
   integer,allocatable :: qpdb  (:)     ! atomic charge in residues from PDB file

   !potential parameters used in energy-gradient routine
   real(wp),allocatable:: vbond(:,:)    ! bonds
   real(wp),allocatable:: vangl(:,:)    ! angles
   real(wp),allocatable:: vtors(:,:)    ! torsions
   real(wp),allocatable:: chieeq(:)     ! atomic ENs for EEQ
   real(wp),allocatable:: gameeq(:)     ! atomic gamma for EEQ
   real(wp),allocatable:: alpeeq(:)     ! atomic alpha for EEQ, squared
   real(wp),allocatable:: alphanb(:)    ! non-bonded exponent for atom pairs
   real(wp),allocatable::    qa(:)      ! estimated atomic charges (fixed and obtained from topology EEQ)
   real(wp),allocatable::     q(:)      ! atomic charges (obtained from EEQ)
   real(wp),allocatable:: hbrefgeo(:,:) ! atom xyz, used to check for HB list update       
   real(wp),allocatable::    xyze0(:,:) ! atom xyz, starting geom. (for Efield energy)     
   real(wp),allocatable:: zetac6(:)     ! D4 scaling factor product 
   real(wp),allocatable:: qfrag (:)     ! fragment charge (for EEQ)
   real(wp),allocatable:: hbbas (:)     ! HB donor atom basicity


   !========================================================================
   ! DATA
   !------------------------------------------------------------------------

   data xhaci / 86 * 0 /
   data xhbas / 86 * 0 /
   data xbaci / 86 * 0 /

   !Pauling EN
   data en/ 2.200,3.000,0.980,1.570,2.040,2.550,3.040,3.440,3.980 &
  &        ,4.500,0.930,1.310,1.610,1.900,2.190,2.580,3.160,3.500 &
  &        ,0.820,1.000,1.360,1.540,1.630,1.660,1.550,1.830,1.880 &
  &        ,1.910,1.900,1.650,1.810,2.010,2.180,2.550,2.960,3.000 &
  &        ,0.820,0.950,1.220,1.330,1.600,2.160,1.900,2.200,2.280 &
  &        ,2.200,1.930,1.690,1.780,1.960,2.050,2.100,2.660,2.600 &
  &,0.79,0.89,1.10,1.12,1.13,1.14,1.15,1.17,1.18,1.20,1.21,1.22 &
  &,1.23,1.24,1.25,1.26,1.27,1.3,1.5,1.7,1.9,2.1,2.2,2.2,2.2 &   ! value of W-Au modified
  &,2.00,1.62,2.33,2.02,2.0,2.2,2.2/

   ! COVALENT RADII, used only in neighbor list determination
   ! based on "Atomic Radii of the Elements," M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
   ! in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
   ! edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
   ! corrected Nov. 17, 2010 for the 92nd edition.
   data rad /&
  &0.32D0,0.37D0,1.30D0,0.99D0,0.84D0,0.75D0,0.71D0,0.64D0,0.60D0,&
  &0.62D0,1.60D0,1.40D0,1.24D0,1.14D0,1.09D0,1.04D0,1.00D0,1.01D0,&
  &2.00D0,1.74D0,1.59D0,1.48D0,1.44D0,1.30D0,1.29D0,1.24D0,1.18D0,&
  &1.17D0,1.22D0,1.20D0,1.23D0,1.20D0,1.20D0,1.18D0,1.17D0,1.16D0,&
  &2.15D0,1.90D0,1.76D0,1.64D0,1.56D0,1.46D0,1.38D0,1.36D0,1.34D0,&
  &1.30D0,1.36D0,1.40D0,1.42D0,1.40D0,1.40D0,1.37D0,1.36D0,1.36D0,&
  &2.38D0,2.06D0,1.94D0,1.84D0,1.90D0,1.88D0,1.86D0,1.85D0,1.83D0,&
  &1.82D0,1.81D0,1.80D0,1.79D0,1.77D0,1.77D0,1.78D0,1.74D0,1.64D0,&
  &1.58D0,1.50D0,1.41D0,1.36D0,1.32D0,1.30D0,1.30D0,1.32D0,1.44D0,&
  &1.45D0,1.50D0,1.42D0,1.48D0,1.46D0/

   data metal / &
  &0,                                                          0,&!He
  &1,1,                                         0, 0, 0, 0, 0, 0,&!Ne
  &1,1,                                         1, 0, 0, 0, 0, 0,&!Ar
  &1,1,2,          2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 0, 0, 0, 0, 0,&!Kr
  &1,2,2,          2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 1, 0, 0, 0, 0,&!Xe
  &1,2,2,  14*2,   2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 1, 1, 1, 0, 0/ !Rn  ! At is NOT a metal, Po is borderline but slightly better as metal
   data group / &
  &1,                                                          8,&!He
  &1,2,                                         3, 4, 5, 6, 7, 8,&!Ne
  &1,2,                                         3, 4, 5, 6, 7, 8,&!Ar
  &1,2,-3,        -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8,&!Kr
  &1,2,-3,        -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8,&!Xe
  &1,2,-3, 14*-3, -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8/ !Rn
   data normcn/ &  ! only for non metals well defined
  &1,                                                          0,&!He
  &4,4,                                         4, 4, 4, 2, 1, 0,&!Ne
  &4,4,                                         4, 4, 4, 2, 1, 0,&!Ar
  &4,4,4,          4, 6, 6, 6, 6, 6, 6, 4, 4,   4, 4, 4, 4, 1, 0,&!Kr
  &4,4,4,          4, 6, 6, 6, 6, 6, 6, 4, 4,   4, 4, 4, 4, 1, 0,&!Xe
  &4,4,4,  14*4,   4, 6, 6, 6, 6, 6, 6, 6, 4,   4, 4, 4, 4, 1, 0/ !Rn
   data repz  / &
  &1.,                                                         2.,&!He
  &1.,2.,                                       3.,4.,5.,6.,7.,8.,&!Ne
  &1.,2.,                                       3.,4.,5.,6.,7.,8.,&!Ar
  &1.,2.,3.,      4.,5.,6.,7.,8.,9.,10.,11.,12.,3.,4.,5.,6.,7.,8.,&!Kr
  &1.,2.,3.,      4.,5.,6.,7.,8.,9.,10.,11.,12.,3.,4.,5.,6.,7.,8.,&!Xe
  &1.,2.,3.,14*3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,3.,4.,5.,6.,7.,8./ !Rn

   !Mass
   real(wp), private              :: ams(107)
   data  ams /  1.00790d0,  4.00260d0,  6.94000d0,  9.01218d0,&
  &10.81000d0, 12.01100d0, 14.00670d0, 15.99940d0, 18.99840d0,&
  &20.17900d0, 22.98977d0, 24.30500d0, 26.98154d0, 28.08550d0,&
  &30.97376d0, 32.06000d0, 35.45300d0, 39.94800d0, 39.09830d0,&
  &40.08000d0, 44.95590d0, 47.90000d0, 50.94150d0, 51.99600d0,&
  &54.93800d0, 55.84700d0, 58.93320d0, 58.71000d0, 63.54600d0,&
  &65.38000d0, 69.73500d0, 72.59000d0, 74.92160d0, 78.96000d0,&
  &79.90400d0, 83.80000d0, 85.46780d0, 87.62000d0, 88.90590d0,&
  &91.22000d0, 92.90640d0, 95.94000d0, 98.90620d0, 101.0700d0,&
  &102.9055d0, 106.4000d0, 107.8680d0, 112.4100d0, 114.8200d0,&
  &118.6900d0, 121.7500d0, 127.6000d0, 126.9045d0, 131.3000d0,&
  &132.9054d0, 137.3300d0,&
  &138.91,140.12,140.91,144.24,147.00,150.36,151.97,157.25,&
  &158.93,162.50,164.93,167.26,168.93,173.04,174.97,&
  &178.4900d0, 180.9479d0,&
  &183.8500d0, 186.2070d0, 190.2000d0, 192.2200d0, 195.0900d0,&
  &196.9665d0, 200.5900d0, 204.3700d0, 207.2000d0, 208.9804d0,&
  &209.,210.,222.,21*0.000d0/

   contains

   subroutine gfnff_set_param(n)  
     use iso_fortran_env, only :   wp => real64
     use gff_d3com                
     implicit none                
!    Dummy                      ,
     integer,intent(in)  :: n     
!    Stack                        
     integer   :: i,j,k           
     real(wp)  :: dum             
                                  
     hbnbcut = hbnbcut*10. ! new parameter for neighbour angle
     cnmax   = 4.4         ! max. CN considered ie all larger values smoothly set to this val
     linthr  = 160.        ! when is an angle close to linear ? (GEODEP) for metals values closer to 170 (than to 160) are better
                           ! but this occurs e.g. for Sc in unclear situations. So make it save (160)
     fcthr   = 1.d-3       ! skip torsion and bending if potential is small
     tdist_thr=12.         ! R threshold in Angstroem for cov distance estimated used in apprx EEQ 
                           ! the following two parameters are critical for topo setup
     rthr     =1.25        ! important bond determination threshold
                           ! large values yield more 1.23
     rthr2    =1.00        ! decrease if a metal is present, larger values yield smaller CN
     rqshrink =0.23        ! change of R0 for topo with charge qa, larger values yield smaller CN for metals in particular
     hqabthr  =0.01        ! H charge (qa) threshold for H in HB list 18
      qabthr  =0.10        ! AB charge (qa) threshold for AB in HB list, avoids HBs with positive atoms, 
                           ! larger val. better for S30L but worse in PubChem RMSD checks
     atcuta  = 0.595d0     ! angle damping
     atcutt  = 0.505d0     ! torsion angle damping 
     atcuta_nci  = 0.395d0 ! nci angle damping in HB term
     atcutt_nci  = 0.305d0 ! nci torsion angle damping in HB term
     srb1    = 0.3731      ! bond params
     srb2    = 0.3171      !      
     srb3    = 0.2538      !      
     repscalb= 1.7583      ! bonded rep. scaling
     repscaln= 0.4270      ! non-bonded rep. scaling
     qrepscal= 0.3480      ! change of non-bonded rep. with q(topo)
     nrepscal=-0.1270      !   "    "      "       "   CN
     hhfac   = 0.6290      ! HH repulsion
     hh13rep = 1.4580      !
     hh14rep = 0.7080      !
     bstren(1)=1.00d0      ! single bond
     bstren(2)=1.24d0      ! double bond
     bstren(3)=1.98d0      ! triple bond
     bstren(4)=1.22d0      ! hyperval bond
     bstren(5)=1.00d0      ! M-X
     bstren(6)=0.78d0      ! M eta
     bstren(7)=3.40d0      ! M-M
     bstren(8)=3.40d0      ! M-M
     qfacBEN =-0.54        ! bend FC change with polarity 
     qfacTOR =12.0d0       ! torsion FC change with polarity
     fr3     =0.3          ! tors FC 3-ring
     fr4     =1.0          ! tors FC 4-ring
     fr5     =1.5          ! tors FC 5-ring
     fr6     =5.7          ! tors FC 6-ring
     torsf(1)=1.00         ! single bond     
     torsf(2)=1.18         ! pi bond
     torsf(3)=1.05         ! improper
     torsf(5)=0.50         ! pi part improper
     torsf(6)=-0.90        ! extra sp3 C
     torsf(7)= 0.70        ! extra sp3 N
     torsf(8)=-2.00        ! extra sp3 O
     fbs1    =0.50         ! small bend corr.
     batmscal=-0.30d0      ! bonded ATM scal
     mchishift=-0.09d0
     rabshift    =-0.110   ! gen shift
     rabshifth   =-0.050   ! XH
     hyper_shift = 0.03    ! hypervalent
     hshift3     = -0.11   ! heavy
     hshift4     = -0.11   !
     hshift5     = -0.06   !
     metal1_shift= 0.2     ! group 1+2 metals
     metal2_shift= 0.15    ! TM
     metal3_shift= 0.05    ! main group metals
     eta_shift   = 0.040   ! eta bonded
     qfacbm(0)   =1.0d0    ! bond charge dep.gff_srcs += 'gff/gfnff_input.f90'

     qfacbm(1:2) =-0.2d0   !
     qfacbm(  3) =0.70d0   !
     qfacbm(  4) =0.50d0   !
     qfacbm0     =0.047    !
     rfgoed1  =1.175       ! topo dist scaling
     htriple = 1.45d0      ! decrease Hueckel off-diag for triple bonds because they are less well conjugated 1.4
     hueckelp2=1.00d0      ! increase pot depth depending on P
     hueckelp3=-0.24d0     ! diagonal element change with qa  
     hdiag(5) =-0.5d0      ! diagonal element relative to C
     hdiag(6) =0.00d0      ! 
     hdiag(7) =0.14d0      ! 
     hdiag(8) =-0.38d0     ! 
     hdiag(9) =-0.29d0     ! 
     hdiag(16)=-0.30d0     ! 
     hdiag(17)=-0.30d0     ! 
     hoffdiag(5)=0.5d0     ! Huckel off-diag constants
     hoffdiag(6)=1.00d0    ! 
     hoffdiag(7)=0.66d0    !
     hoffdiag(8)=1.10d0    ! 
     hoffdiag(9)=0.23d0    ! 
     hoffdiag(16)=0.60d0   ! 
     hoffdiag(17)=1.00d0   ! 
     hiter   =0.700d0      ! iteration mixing
     hueckelp=0.340d0      ! diagonal qa dep.
     bzref   =0.370d0      ! ref P value R shift
     bzref2  =0.315d0      !  "  "  "    k stretch
     pilpf   =0.530d0      ! 2el diag shift
     maxhiter=5            ! the Hückel iterations can diverge so take only a few steps
     d3a1    = 0.58d0      ! D3, s8 fixed = 2
     d3a2    = 4.80d0
     split0  =0.670d0      ! mixing of sp^n with sp^n-1
     fringbo =0.020d0      ! str ring size dep.
     aheavy3 =89.          ! three coord. heavy eq. angle
     aheavy4 =100.         ! four   "       "    "    "
     hbacut   =49.         ! HB angle cut-off
     hbscut   =22.         ! HB SR     "   "
     xbacut   =70.         ! same for XB 
     xbscut   = 5.         ! 
     hbsf     = 1.         ! charge dep.
     hbst     =15.         ! 10 is better for S22, 20 better for HCN2 and S30L
     xbsf     =0.03        !
     xbst     =15.         !
     hbalp    = 6.         ! damp
     hblongcut=85.         ! values larger than 85 yield large RMSDs for P26
     hblongcut_xb=70.      ! values larger than 70 yield large MAD for HAL28
     hbabmix  =0.80        ! 
     hbnbcut  =11.20       !
     tors_hb   =0.94       ! torsion potential shift in HB term
     bend_hb   =0.20       ! bending potential shift in HB term
     vbond_scale=0.9       ! vbond(2) scaling for CN(H) = 1
     xhaci_globabh=0.268   ! A-H...B gen. scaling
     xhaci_coh=0.350       ! A-H...O=C gen. scaling
     xhaci_glob=1.50       ! acidity 
     xhbas( 6)=0.80d0      ! basicities (XB and HB), i.e., B...X-A or B...H..A
     xhbas( 7)=1.68d0   
     xhbas( 8)=0.67d0  
     xhbas( 9)=0.52d0  
     xhbas(14)=4.0d0    
     xhbas(15)=3.5d0      
     xhbas(16)=2.0d0
     xhbas(17)=1.5d0
     xhbas(35)=1.5d0
     xhbas(53)=1.9d0
     xhbas(33)=xhbas(15)   
     xhbas(34)=xhbas(16)
     xhbas(51)=xhbas(15) 
     xhbas(52)=xhbas(16)
     xhaci( 6)=0.75               ! HB acidities, a bit weaker for CH
     xhaci( 7)=xhaci_glob+0.1  
     xhaci( 8)=xhaci_glob 
     xhaci( 9)=xhaci_glob 
     xhaci(15)=xhaci_glob 
     xhaci(16)=xhaci_glob 
     xhaci(17)=xhaci_glob+1.0 
     xhaci(35)=xhaci_glob+1.0 
     xhaci(53)=xhaci_glob+1.0 
     xbaci(15)=1.0d0              ! XB acidities
     xbaci(16)=1.0d0
     xbaci(17)=0.5d0
     xbaci(33)=1.2d0
     xbaci(34)=1.2d0
     xbaci(35)=0.9d0
     xbaci(51)=1.2d0
     xbaci(52)=1.2d0
     xbaci(53)=1.2d0
     split1=1.0d0-split0
     bsmat = -999.
     bsmat(0,0)=bstren(1)
     bsmat(3,0)=bstren(1)
     bsmat(3,3)=bstren(1)
     bsmat(2,2)=bstren(2)
     bsmat(1,1)=bstren(3)
     bsmat(1,0)=split0*bstren(1)+split1*bstren(3)
     bsmat(3,1)=split0*bstren(1)+split1*bstren(3)
     bsmat(2,1)=split0*bstren(2)+split1*bstren(3)
     bsmat(2,0)=split0*bstren(1)+split1*bstren(2)
     bsmat(3,2)=split0*bstren(1)+split1*bstren(2)
     bstren(9)=0.5*(bstren(7)+bstren(8))

!    3B bond prefactors and D3 stuff
     k=0
     do i=1,86
        dum=dble(i)
        zb3atm(i)=dum*batmscal**(1.d0/3.d0)  ! inlcude pre-factor
        do j=1,i
           k=k+1
           dum=r2r4(i)*r2r4(j)*3.0d0    
           d3r0(k)=(d3a1*dsqrt(dum)+d3a2)**2   ! save R0^2 for efficiency reasons
        enddo
     enddo
     zb3atm(1)=0.25d0*batmscal**(1.d0/3.d0) ! slightly better than 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! numerical precision settings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      accff  =0.1d0  ! global accuracy factor similar to acc in xtb used in SCF
      if(n.gt.10000) accff = 2.0d0

      dispthr=1500.0d0-log10(accff)*1000.d0
      cnthr  = 100.0d0-log10(accff)* 50.0d0
      repthr = 400.0d0-log10(accff)*100.0d0
      hbthr1 = 200.0d0-log10(accff)* 50.0d0
      hbthr2 = 400.0d0-log10(accff)* 50.0d0

   end subroutine gfnff_set_param

   subroutine gfnff_read_param(iunit)
     use iso_fortran_env, only : wp => real64
     implicit none
!    Dummy
     integer,intent(in)  :: iunit
!    Stack     
     integer  :: i,nn
     real(wp) :: xx(20)
     character(len=256) :: atmp
     
     do i=1,86
         read(iunit,'(a)')atmp
         call readl(atmp,xx,nn)
         chi (i)=xx(2)
         gam (i)=xx(3)
         cnf (i)=xx(4)
         alp (i)=xx(5)
         bond(i)=xx(6)
         repa(i)=xx(7)
        repan(i)=xx(8)
         angl(i)=xx(9)
        angl2(i)=xx(10)
         tors(i)=xx(11)
        tors2(i)=xx(12)
      enddo
     
   end subroutine gfnff_read_param

   subroutine gfnff_param_alloc(n)
     use iso_fortran_env, only : wp => real64
     implicit none
!    Dummy
     integer,intent(in) :: n
   
     if (.not.allocated(nb)) allocate( nb(20,n), source = 0 )      
     if (.not.allocated(bpair)) allocate( bpair(n*(n+1)/2), source = 0 )      
     if (.not.allocated(alphanb)) allocate( alphanb(n*(n+1)/2), source = 0.0d0 )
     if (.not.allocated(chieeq)) allocate( chieeq(n), source = 0.0d0 )
     if (.not.allocated(gameeq)) allocate( gameeq(n), source = 0.0d0 )
     if (.not.allocated(alpeeq)) allocate( alpeeq(n), source = 0.0d0 )
     if (.not.allocated(qa)) allocate( qa(n), source = 0.0d0 )
     if (.not.allocated(q)) allocate( q(n), source = 0.0d0 )
     if (.not.allocated(hbrefgeo)) allocate( hbrefgeo(3,n), source = 0.0d0 )
     if (.not.allocated(zetac6)) allocate( zetac6(n*(n+1)/2), source = 0.0d0 )
     if (.not.allocated(xyze0)) allocate( xyze0(3,n), source = 0.0d0 )
     if (.not.allocated(b3list)) allocate( b3list(3,1000*n), source = 0 )
     if (.not.allocated(fraglist)) allocate( fraglist(n), source = 0 )
     if (.not.allocated(qfrag)) allocate( qfrag(n), source = 0.0d0 )
     if (.not.allocated(hbatHl)) allocate( hbatHl(n), source = 0 )
     if (.not.allocated(hbbas)) allocate( hbbas(n), source = 0.0d0 )
     if (.not.allocated(hbatABl)) allocate( hbatABl(2,n*(n+1)/2), source = 0 )
     if (.not.allocated(xbatABl)) allocate( xbatABl(3,natxbAB), source = 0 )

     if (.not.allocated(blist)) allocate( blist(2,nbond_blist), source = 0 )
     if (.not.allocated(nr_hb)) allocate( nr_hb(nbond_blist), source = 0 )
     if (.not.allocated(bond_hb_AH)) allocate( bond_hb_AH(2,bond_hb_nr), source = 0 )
     if (.not.allocated(bond_hb_B)) allocate( bond_hb_B(b_max,bond_hb_nr), source = 0 )
     if (.not.allocated(bond_hb_Bn)) allocate( bond_hb_Bn(bond_hb_nr), source = 0 )
     if (.not.allocated(alist)) allocate( alist(3,nangl_alloc), source = 0 )
     if (.not.allocated(tlist)) allocate( tlist(5,ntors_alloc), source = 0 )
     if (.not.allocated(vbond)) allocate( vbond(3,nbond_vbond), source = 0.0d0 )
     if (.not.allocated(vangl)) allocate( vangl(2,nangl_alloc), source = 0.0d0 )
     if (.not.allocated(vtors)) allocate( vtors(2,ntors_alloc), source = 0.0d0 )
     if (.not.allocated(hblist1)) allocate( hblist1(3,nhb1), source = 0 )
     if (.not.allocated(hblist2)) allocate( hblist2(3,nhb2), source = 0 )
     if (.not.allocated(hblist3)) allocate( hblist3(3,nxb), source = 0 )

   end subroutine gfnff_param_alloc

   subroutine gfnff_param_dealloc()
     use iso_fortran_env, only : wp => real64
     implicit none
!    Dummy
   
     if (allocated(nb)) deallocate( nb )
     if (allocated(bpair)) deallocate( bpair )      
     if (allocated(alphanb)) deallocate( alphanb )
     if (allocated(chieeq)) deallocate( chieeq )
     if (allocated(gameeq)) deallocate( gameeq )
     if (allocated(alpeeq)) deallocate( alpeeq )
     if (allocated(qa)) deallocate( qa )
     if (allocated(q)) deallocate( q )
     if (allocated(hbrefgeo)) deallocate( hbrefgeo )
     if (allocated(zetac6)) deallocate( zetac6 )
     if (allocated(xyze0)) deallocate( xyze0 )
     if (allocated(b3list)) deallocate( b3list )
     if (allocated(fraglist)) deallocate( fraglist )
     if (allocated(qfrag)) deallocate( qfrag )
     if (allocated(hbatHl)) deallocate( hbatHl )
     if (allocated(hbbas)) deallocate( hbbas )
     if (allocated(hbatABl)) deallocate( hbatABl )
     if (allocated(xbatABl)) deallocate( xbatABl )

     if (allocated(blist)) deallocate( blist )
     if (allocated(nr_hb)) deallocate( nr_hb )
     if (allocated(bond_hb_AH)) deallocate( bond_hb_AH )
     if (allocated(bond_hb_B)) deallocate( bond_hb_B )
     if (allocated(bond_hb_Bn)) deallocate( bond_hb_Bn )
     if (allocated(alist)) deallocate( alist )
     if (allocated(tlist)) deallocate( tlist )
     if (allocated(vbond)) deallocate( vbond )
     if (allocated(vangl)) deallocate( vangl )
     if (allocated(vtors)) deallocate( vtors )
     if (allocated(hblist1)) deallocate( hblist1 )
     if (allocated(hblist2)) deallocate( hblist2 )
     if (allocated(hblist3)) deallocate( hblist3 )

   end subroutine gfnff_param_dealloc

end module gff_param
