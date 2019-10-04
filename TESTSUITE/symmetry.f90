subroutine test_symmetry_water
   use iso_fortran_env, wp => real64
   use assertion
   use thermo

   implicit none

   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp)           :: xyz(3,nat) = reshape(&
      [ 0.00000000000000_wp,  0.00000000000000_wp, -0.73578586109551_wp, &
      & 1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp, &
      &-1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp  &
      & ],shape(xyz))
   real(wp),parameter :: desy = 0.1_wp
   integer, parameter :: maxatdesy = 200
   character(len= 4)  :: pgroup = 'C1  '
   integer, parameter :: iunit = output_unit
   logical, parameter :: pr = .true.

   call rattle(nat,xyz,0.025_wp)

   call getsymmetry(pr,iunit,nat,at,xyz,desy,maxatdesy,pgroup)

   if (pgroup.ne."c2v") stop 1

   call terminate(afail)

end subroutine test_symmetry_water

subroutine test_symmetry_li8
   use iso_fortran_env, wp => real64
   use assertion
   use thermo

   implicit none

   integer, parameter :: nat = 8
   integer, parameter :: at(nat) = [3,3,3,3,3,3,3,3]
   real(wp)           :: xyz(3,nat) = reshape(&
      &[3.44376028745523_wp,-3.44376028745523_wp, 3.44376028745523_wp, &
      &-3.44376028745523_wp, 3.44376028745523_wp, 3.44376028745523_wp, &
      & 3.44376028745523_wp, 3.44376028745523_wp,-3.44376028745523_wp, &
      &-3.44376028745523_wp,-3.44376028745523_wp,-3.44376028745523_wp, &
      & 1.96829703889624_wp, 1.96829703889624_wp, 1.96829703889624_wp, &
      &-1.96829703889624_wp, 1.96829703889624_wp,-1.96829703889624_wp, &
      &-1.96829703889624_wp,-1.96829703889624_wp, 1.96829703889624_wp, &
      & 1.96829703889624_wp,-1.96829703889624_wp,-1.96829703889624_wp  &
      & ],shape(xyz))
   real(wp),parameter :: desy = 0.1_wp
   integer, parameter :: maxatdesy = 200
   character(len= 4)  :: pgroup = 'C1  '
   integer, parameter :: iunit = output_unit
   logical, parameter :: pr = .true.

   call rattle(nat,xyz,0.025_wp)

   call getsymmetry(pr,iunit,nat,at,xyz,desy,maxatdesy,pgroup)

   if (pgroup.ne."td") stop 1

   call terminate(afail)

end subroutine test_symmetry_li8

subroutine test_symmetry_pcl3
   use iso_fortran_env, wp => real64
   use assertion
   use thermo

   implicit none

   integer, parameter :: nat = 4
   integer, parameter :: at(nat) = [15,17,17,17]
   real(wp)           :: xyz(3,nat) = reshape([&
      &-0.3927247746_wp,    3.3399945638_wp,    0.0000000000_wp, &
      & 2.0622241137_wp,    3.7182629335_wp,    2.9906068785_wp, &
      & 2.0622241137_wp,    3.7182629335_wp,   -2.9906068785_wp, &
      &-2.0622241137_wp,    6.8478421617_wp,    0.0000000000_wp  &
      & ],shape(xyz))
   real(wp),parameter :: desy = 0.1_wp
   integer, parameter :: maxatdesy = 200
   character(len= 4)  :: pgroup = 'C1  '
   integer, parameter :: iunit = output_unit
   logical, parameter :: pr = .true.

   call rattle(nat,xyz,0.025_wp)

   call getsymmetry(pr,iunit,nat,at,xyz,desy,maxatdesy,pgroup)

   if (pgroup.ne."c3v") stop 1

   call terminate(afail)

end subroutine test_symmetry_pcl3

subroutine test_symmetry_c20
   use iso_fortran_env, wp => real64
   use assertion
   use thermo

   implicit none

   integer, parameter :: nat = 60
   integer, parameter :: at(nat) = [6,6,1,6,6,6,6,6,1,6,6,6,1,6,6,6,1,6, &
      & 6,6,1,6,1,6,6,1,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
      & 1,1,1,1,1,1,1,1,1,1]
   real(wp)           :: xyz(3,nat) = reshape(&
      &[2.41219194492747_wp,-4.91851303826010_wp,-0.22758538736182_wp, &
      & 4.82472026110659_wp,-5.14171048154136_wp, 1.32760630767135_wp, &
      & 4.77295499314845_wp,-3.75166762243947_wp, 2.85870210478409_wp, &
      & 7.25869507916783_wp,-4.83305017390035_wp,-0.19559043420593_wp, &
      & 7.47000803479937_wp,-2.36707851295709_wp,-1.68298442476956_wp, &
      & 7.18881867595212_wp, 0.01015538823813_wp,-0.08242607446421_wp, &
      & 7.47520667139093_wp, 2.42658976833392_wp,-1.62214091246904_wp, &
      & 7.23773612662799_wp, 4.85495021005436_wp,-0.07807403518017_wp, &
      & 7.43559045274112_wp, 6.47235924828665_wp,-1.35229557959750_wp, &
      & 4.77842386057702_wp, 5.12562513269823_wp, 1.41156683975546_wp, &
      & 2.39275611165074_wp, 4.92495322492098_wp,-0.18756854677183_wp, &
      &-0.02299796703722_wp, 5.26631335356544_wp, 1.34096856115698_wp, &
      &-0.03259777579228_wp, 7.13976327526854_wp, 2.22155637264196_wp, &
      &-2.41218627574907_wp, 4.91807273207114_wp,-0.22732082570321_wp, &
      &-4.82480340905644_wp, 5.14144591988275_wp, 1.32769890425186_wp, &
      &-7.25865917437130_wp, 4.83304639444808_wp,-0.19573972257043_wp, &
      &-8.86910079235165_wp, 4.98168658288247_wp, 1.09382261859367_wp, &
      &-7.46991354849272_wp, 2.36719567597733_wp,-1.68336236999613_wp, &
      &-7.18862025470816_wp,-0.01008357864508_wp,-0.08290606490196_wp, &
      &-7.47523312755680_wp,-2.42654819435900_wp,-1.62250751933882_wp, &
      &-6.07711735841522_wp,-2.42521404770918_wp,-3.14764689968070_wp, &
      &-7.23779281841198_wp,-4.85479147305920_wp,-0.07825355916279_wp, &
      &-8.82692399479177_wp,-4.97337745707617_wp, 1.24058441925196_wp, &
      &-4.77845409619515_wp,-5.12534923268283_wp, 1.41135896988084_wp, &
      &-2.39279957535179_wp,-4.92475669340316_wp,-0.18781421116911_wp, &
      &-2.32699175250018_wp,-3.07877394104155_wp,-1.12088538654272_wp, &
      & 0.02293560607483_wp,-5.26647398028673_wp, 1.34066431524959_wp, &
      & 0.03469348207365_wp,-3.90425355903932_wp, 2.90074284206239_wp, &
      & 0.03226518399289_wp,-7.13995602733409_wp, 2.22119543495058_wp, &
      &-2.45713530154589_wp,-6.32739646690686_wp,-1.70951995912754_wp, &
      &-4.70585081816968_wp,-3.70545814931202_wp, 2.91385376197235_wp, &
      &-4.79280467644828_wp,-6.95357611829985_wp, 2.37843009783733_wp, &
      &-7.43573029247495_wp,-6.47233090239466_wp,-1.35229935904976_wp, &
      &-9.31829625304309_wp,-2.41906487887277_wp,-2.56133857996547_wp, &
      &-5.33755035959683_wp, 0.00004535342719_wp, 0.83753040154718_wp, &
      &-8.58318711817218_wp, 0.01042183962286_wp, 1.44852232401017_wp, &
      &-6.05623399492070_wp, 2.33441459675016_wp,-3.19366362074260_wp, &
      &-9.30322001795493_wp, 2.32943327866387_wp,-2.64038771383021_wp, &
      &-7.43720049940634_wp, 6.41898960284169_wp,-1.51171476589386_wp, &
      &-4.86035482679441_wp, 6.98806739967728_wp, 2.25847784182628_wp, &
      &-4.77331782056597_wp, 3.75132747173555_wp, 2.85873045067608_wp, &
      &-2.44550214747184_wp, 6.30647530888968_wp,-1.76300676759273_wp, &
      &-2.33508355980120_wp, 3.06328196620415_wp,-1.14178764729857_wp, &
      &-0.03458765741021_wp, 3.90403813026017_wp, 2.90100551399486_wp, &
      & 2.45698223372912_wp, 6.32771394089719_wp,-1.70916469061455_wp, &
      & 2.32713915113855_wp, 3.07902716434336_wp,-1.12077767215315_wp, &
      & 4.70583570036062_wp, 3.70582664590793_wp, 2.91414289007068_wp, &
      & 4.79280845590055_wp, 6.95391059982537_wp, 2.37851324578718_wp, &
      & 8.82689186944751_wp, 4.97377807901634_wp, 1.24070536172447_wp, &
      & 9.31820365646258_wp, 2.41910456312156_wp,-2.56111748200792_wp, &
      & 6.07699263649045_wp, 2.42512334085480_wp,-3.14720281403947_wp, &
      & 8.58359340929075_wp,-0.01037837592181_wp, 1.44885113635729_wp, &
      & 5.33779224454184_wp,-0.00002078698746_wp, 0.83817101870623_wp, &
      & 6.05631525314441_wp,-2.33413302755636_wp,-3.19328000633762_wp, &
      & 9.30332395289224_wp,-2.32931422591750_wp,-2.64000031997296_wp, &
      & 7.43737435421056_wp,-6.41882330694200_wp,-1.51175822959492_wp, &
      & 8.86907244645965_wp,-4.98172626713126_wp, 1.09405505490801_wp, &
      & 4.86034915761601_wp,-6.98837353531080_wp, 2.25830209729592_wp, &
      & 2.33516859747718_wp,-3.06384888404402_wp,-1.14234889596004_wp, &
      & 2.44560041323075_wp,-6.30719907399858_wp,-1.76302566485406_wp  &
      & ],shape(xyz))
   real(wp),parameter :: desy = 0.25_wp
   integer, parameter :: maxatdesy = 200
   character(len= 4)  :: pgroup = 'C1  '
   integer, parameter :: iunit = output_unit
   logical, parameter :: pr = .true.

   call rattle(nat,xyz,0.025_wp)

   call getsymmetry(pr,iunit,nat,at,xyz,0.1_wp,maxatdesy,pgroup)

   if (pgroup.ne."c2") stop 1

   call getsymmetry(pr,iunit,nat,at,xyz,0.2_wp,maxatdesy,pgroup)

   if (pgroup.ne."c2v") stop 1

   call terminate(afail)

end subroutine test_symmetry_c20

subroutine rattle(nat,xyz,magnitude)
   use iso_fortran_env, wp => real64

   implicit none
   integer, intent(in) :: nat
   real(wp), intent(out) :: xyz(3,nat)
   real(wp), intent(in) :: magnitude
   real(wp) :: ran
   integer  :: i,j

   do i = 1, nat
      do j = 1, 3
         call random_number(ran)
         xyz(j,i) = xyz(j,i) + (-1)**(i+j)*magnitude*ran
      enddo
   enddo

end subroutine rattle
