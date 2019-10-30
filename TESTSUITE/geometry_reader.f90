!> @brief test disperion under 3D periodic boundary conditions
subroutine test_geometry_reader_file_poscar_sio2_3d
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-9_wp
   character(len=*),parameter :: file_poscar_sio2_3d = &
      & '("Si  O ",/,&
      & " 1.0000000000000000",/,&
      & "     4.6257000000000001    0.0000000000000000    0.0000000000000000",/,&
      & "     0.0000000000000000    4.6257000000000001    0.0000000000000000",/,&
      & "     0.0000000000000000    0.0000000000000000    4.6257000000000001",/,&
      & "   2   4",/,&
      & "Cartesian",/,&
      & "  0.44400000000000000  0.230000000000000000  0.470000000000000000",/,&
      & "  2.53128500000000001  2.445312850000000000  1.7884903000000000000",/,&
      & "  1.38049099700000002  1.83049099700000002   0.940000000000000000",/,&
      & "  1.08079400299999999  3.96177599700000003   1.474903000000000000",/,&
      & "  3.34207900299999999  3.03207900299999999   0.260000000000000000",/,&
      & "  3.26177599700000003  1.310079400299999999  1.684903000000000000")'
   integer :: iunit
   type(tb_molecule) :: mol

   open(newunit=iunit,status='scratch')
   write(iunit,file_poscar_sio2_3d)
   rewind(iunit)
   call mol%read(iunit, format=p_ftype%vasp)

   call assert_close(mol%volume,       667.92680030347_wp,thr)
   call assert_close(mol%cellpar(1),8.7413053236641_wp,thr)
   call assert_close(mol%cellpar(4),1.5707963267949_wp,thr)
   call assert_close(mol%rec_lat(1,1),0.71879256867621_wp,thr)

   call assert_close(mol%abc(1,3),0.29843937068984_wp,thr)
   call assert_close(mol%abc(2,6),0.28321754549582_wp,thr)
   call assert_close(mol%abc(3,1),0.10160624337938_wp,thr)

   call assert_close(mol%xyz(2,4),7.4866709068338_wp,thr)
   call assert_close(mol%xyz(1,5),6.3156134163815_wp,thr)
   call assert_close(mol%xyz(3,2),3.3797565299764_wp,thr)

   call mol%deallocate
   close(iunit,status='delete')

   call terminate(afail)

end subroutine test_geometry_reader_file_poscar_sio2_3d

subroutine test_geometry_reader_file_xmol_water_0d
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-10_wp
   character(len=*),parameter :: file_xmol_water_0d = &
      & '("9",/,&
      & "WATER27, (H2O)3",/,&
      & "O     1.1847029    1.1150792   -0.0344641 ",/,&
      & "H     0.4939088    0.9563767    0.6340089 ",/,&
      & "H     2.0242676    1.0811246    0.4301417 ",/,&
      & "O    -1.1469443    0.0697649    1.1470196 ",/,&
      & "H    -1.2798308   -0.5232169    1.8902833 ",/,&
      & "H    -1.0641398   -0.4956693    0.3569250 ",/,&
      & "O    -0.1633508   -1.0289346   -1.2401808 ",/,&
      & "H     0.4914771   -0.3248733   -1.0784838 ",/,&
      & "H    -0.5400907   -0.8496512   -2.1052499 ")'
   integer :: iunit
   type(tb_molecule) :: mol

   open(newunit=iunit,status='scratch')
   write(iunit,file_xmol_water_0d)
   rewind(iunit)
   call mol%read(iunit, format=p_ftype%xyz)

   call assert_eq(mol%n,9)

   call assert_close(mol%xyz(1,2), .93335227594625_wp,thr)
   call assert_close(mol%xyz(2,5),-.98873655304085_wp,thr)
   call assert_close(mol%xyz(1,9),-1.0206234107641_wp,thr)
   call assert_close(mol%xyz(3,3), .81284993236482_wp,thr)

   call mol%deallocate
   close(iunit,status='delete')

   call terminate(afail)

end subroutine test_geometry_reader_file_xmol_water_0d

subroutine test_geometry_reader_file_coord_general_0d
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-10_wp
   character(len=*),parameter :: file_coord_general_0d_p1 = &
      & '("$coord",/,&
      & "   -2.83449054141383   -0.81686825435058    9.21776298032051      c ",/,&
      & "   -4.14924213549981   -2.27671857946298    7.51444609965666      n ",/,&
      & "   -4.89221891191980   -4.38375126958911    8.73662159082462      c ",/,&
      & "   -4.08966629931927   -4.37844740936509   11.18148413192309      c ",/,&
      & "   -2.76412704180850   -2.06518595253845   11.47664619656747      c ",/,&
      & "   -5.27341047178455   -1.44694863069164    3.79024851342436      sb",/,&
      & "   -2.72022250443229   -4.11374003325320    2.02589332400802      c ",/,&
      & "   -1.68673671654791   -2.71845226880146   -1.23223099384719      ge",/,&
      & "   -0.33836616422483    0.68435923235649   -0.64698844235331      c ",/,&
      & "   -2.14498825611689    2.50532279558922    0.70167668533214      c ",/,&
      & "   -2.64687090621542    1.70512649775402    3.45551752220984      c ")'
   character(len=*),parameter :: file_coord_general_0d_p2 = &
      & '("   -1.15030356647053    5.16801694545513    0.70880247141464      c ",/,&
      & "   -2.66291175757117    7.37686607275785    0.92639065980257      c ",/,&
      & "   -1.04382559568649    9.50564213115659    0.99473552928660      c ",/,&
      & "    1.47621944151256    8.62625177628255    0.82443436691346      c ",/,&
      & "    1.41445008883828    5.95543737803193    0.65952274533898      c ",/,&
      & "   -0.54529064133928    7.51336038163097   -2.30620237507063      fe",/,&
      & "   -1.13319651552659    9.96175330724574   -5.27221419699432      c ",/,&
      & "   -2.96542350971078    8.01125651836652   -5.31023559376596      c ",/,&
      & "   -1.67962454162497    5.67404846994265   -5.54045935015240      c ",/,&
      & "    0.94883812131510    6.18164739055204   -5.64434905287388      c ",/,&
      & "    1.28552802255622    8.83144502472722   -5.47753977544342      c ",/,&
      & "   -4.98685670751547    8.26196598048072   -5.18209722842839      h ",/,&
      & "    3.07242445834905    9.81671524766024   -5.49709647604619      h ",/,&
      & "    2.43752344355029    4.79591597064188   -5.81797108098135      h ",/,&
      & "   -2.55593371358751    3.83162467300035   -5.63320144731592      h ",/,&
      & "   -1.51348788464750   11.95966574727045   -5.10636613657016      h ",/,&
      & "   -4.70456421906849    7.41962386476884    1.00495597944190      h ",/,&
      & "    3.15509261500725    9.78658966959972    0.80555833604958      h ",/,&
      & "    3.04787658110567    4.74120850762418    0.49603302321748      h ")'
   character(len=*),parameter :: file_coord_general_0d_p3 = &
      & '("   -1.63030870810950   11.45617324545936    1.12534439746916      h ",/,&
      & "   -3.94240915981432    2.55539341594866   -0.33639931534539      h ",/,&
      & "   -0.86691459289722    1.21603908477650    4.38240677647076      h ",/,&
      & "   -3.46915354865397    3.27986411922675    4.50330918125007      h ",/,&
      & "   -1.08871309678725   -4.32760875840763    3.27606487670684      h ",/,&
      & "   -3.63656038759555   -5.95179774391449    1.84321497866878      h ",/,&
      & "    0.88793735017662   -4.80885835927470   -2.83254792820997      c ",/,&
      & "   -3.95291614722229   -2.56392761112535   -3.02791467140108      h ",/,&
      & "    1.40971969416168    0.50275121457463    0.44670554030343      h ",/,&
      & "    0.21425359037146    1.48744908535502   -2.46754028898236      h ",/,&
      & "   -2.04433873645781    0.98960427068737    8.71776817187328      h ",/,&
      & "   -1.20957134640006   -0.97148940925269   14.44264243113319      br",/,&
      & "   -4.55398310845333   -6.66965646831917   13.38809200915616      cl",/,&
      & "   -6.21589279677564   -6.11621330224226    7.50219248533889      f ",/,&
      & "    0.96920569809100   -5.10734351153847   -5.44441038463856      c ",/,&
      & "    2.87016003827288   -6.50046794479532   -6.58587442261051      c ",/,&
      & "    4.70536813226671   -7.60983664254865   -5.09583320473503      c ",/,&
      & "    4.65728219279955   -7.40871401511349   -2.47746786874799      c ",/,&
      & "    2.74659902402950   -5.98087910118857   -1.38443645887459      c ",/,&
      & "   -0.45217183205771   -4.24042173747880   -6.63573003838007      h ",/,&
      & "    2.94369197343435   -6.69065625615941   -8.61819765502526      h ",/,&
      & "    6.81472990694335   -8.85338882847323   -6.37003101469018      n ",/,&
      & "    6.45935584379557   -8.70259871911758   -0.72993743886431      c ",/,&
      & "    2.72403966645085   -5.82586224146740    0.65750112692980      h ",/,&
      & "    8.86710352143466   -8.63574098892235   -5.40639790188841      o ",/,&
      & "    6.36615538013801   -9.93634291193171   -8.32176638301194      o ",/,&
      & "    6.73355918399491  -11.13036831778917   -1.30541276090452      o ",/,&
      & "    7.38825871001876   -7.68904799991807    1.06638373991088      o ",/,&
      & "    7.94931938543849  -11.87978475199051   -0.16550598493729      h ",/,&
      & "$end")'
   integer :: iunit
   type(tb_molecule) :: mol

   open(newunit=iunit,status='scratch')
   write(iunit,file_coord_general_0d_p1)
   write(iunit,file_coord_general_0d_p2)
   write(iunit,file_coord_general_0d_p3)
   rewind(iunit)

   call mol%read(iunit, format=p_ftype%tmol)

   call assert(.not.any(mol%pbc))
   call assert_eq(mol%n, 59)

   call mol%deallocate; close(iunit,status='delete')

   call terminate(afail)

end subroutine test_geometry_reader_file_coord_general_0d

subroutine test_geometry_reader_file_coord_CaF2_3d
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-10_wp
   character(len=*),parameter :: file_coord_CaF2_3d = &
      & '("$coord frac",/,&
      & "    0.2500000000      0.2500000000      0.2500000000      f",/,&
      & "    0.7500000000      0.7500000000      0.7500000000      f",/,&
      & "    0.0000000000      0.0000000000      0.0000000000      ca",/,&
      & "$user-defined bonds",/,&
      & "$lattice angs",/,&
      & "   3.153833580475253   1.115048555743951   1.931320751454818",/,&
      & "   0.000000000000000   3.345145667231851   1.931320751454818",/,&
      & "   0.000000000000000   0.000000000000000   3.862641502909638",/,&
      & "$periodic 3",/,&
      & "$end")'
   integer :: iunit
   type(tb_molecule) :: mol

   open(newunit=iunit,status='scratch')

   write(iunit,file_coord_CaF2_3d); rewind(iunit)
   ! reads in from cell parameters in bohr and coordinates in bohr
   call mol%read(iunit, format=p_ftype%tmol)

   call assert(all(mol%pbc))
   call assert_eq(mol%npbc,3)

   call assert_eq(mol%n, 3)

   call assert_close(mol%volume,       275.00126402469_wp,thr)

   call assert_close(mol%xyz(1,1),1.4899702891972_wp,thr)
   call assert_close(mol%xyz(2,1),2.1071361905157_wp,thr)
   call assert_close(mol%xyz(2,2),6.3214085715472_wp,thr)

   call assert_close(mol%rec_lat(1,2),-0.35141558150155_wp,thr)
   call assert_close(mol%rec_lat(3,3), 0.86078886234226_wp,thr)

   call assert_close(mol%cellpar(1),mol%cellpar(2),thr)
   call assert_close(mol%cellpar(4),mol%cellpar(6),thr)
   call assert_close(mol%cellpar(3),7.2993338808807_wp,thr)
   call assert_close(mol%cellpar(5),1.0471975511966_wp,thr)

   call mol%deallocate; close(iunit,status='delete')

   call terminate(afail)

end subroutine test_geometry_reader_file_coord_CaF2_3d

subroutine test_geometry_reader_file_coord_CaMgCO_3d
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-10_wp
   character(len=*),parameter :: file_coord_CaMgCO_3d = &
      & '("$cell",/,&
      & " 9.09903133 9.09903130512 30.4604956 90.0 90.0 120.000000127",/,&
      & "$coord",/,&
      & "   -0.579494558      0.0683589331    -7.5199348400     ca",/,&
      & "   -0.579494558      0.0683589331     7.7103129400     mg",/,&
      & "   -0.579494558      0.0683589331    -0.1028041720     c",/,&
      & "    1.738483670     -0.2050767990    -0.0875739247     o",/,&
      & "$periodic 3",/,&
      & "$end")'
   integer :: iunit
   type(tb_molecule) :: mol

   open(newunit=iunit,file='coord')

   write(iunit,file_coord_CaMgCO_3d); rewind(iunit)
   ! reads in from cell parameters in bohr and coordinates in bohr
   call mol%read(iunit, format=p_ftype%tmol)

   call assert(all(mol%pbc))
   call assert_eq(mol%npbc,3)

   call assert_eq(mol%n, 4)

   call assert_close(mol%volume,       2184.02656187534_wp,thr)
   call assert_close(mol%lattice(1,2),-4.54951567002654_wp,thr)
   call assert_close(mol%lattice(2,2), 7.87999224997948_wp,thr)
   call assert_close(mol%rec_lat(1,1),0.690533429252363_wp,thr)
   call assert_close(mol%rec_lat(2,1),0.398679663304106_wp,thr)
   call assert_close(mol%rec_lat(3,3),0.206273246164110_wp,thr)

   ! this value should get's wrapped back from -7.51993484
   call assert_close(mol%xyz(1,3),    8.5195367720000_wp,thr)

   call mol%deallocate; close(iunit,status='delete')

   call terminate(afail)

end subroutine test_geometry_reader_file_coord_CaMgCO_3d

subroutine test_geometry_reader_file_coord_C_2d
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-10_wp

   character(len=*),parameter :: file_coord_C_2d = &
      & '("$coord",/,&
      & "   -0.12918412100093      0.06210659750976     -2.13384498734326  c",/,&
      & "    0.12856915667443     -0.07403227791901      4.02358027265954  c",/,&
      & "   -0.12317720857511      2.75170732207802     -2.13345350602279  c",/,&
      & "    2.44816466162280      1.28612566399214      4.02317048854901  c",/,&
      & "$periodic 2",/,&
      & "$cell  angs",/,&
      & "    2.4809835980     2.4811430162   120.2612191150",/,&
      & "$end")'

   integer :: iunit
   type(tb_molecule) :: mol

   stop 77

end subroutine test_geometry_reader_file_coord_C_2d

subroutine test_geometry_reader_file_coord_C_1d
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-10_wp

   character(len=*),parameter :: file_coord_C_1d = &
      & '("$coord",/,&
      & "   1.36794785746435    13.45808943446053     8.83754983226359     c",/,&
      & "   3.69183290816438    13.13552229161569    10.16652201690950     c",/,&
      & "   1.36792668081267    10.38660504434782    13.04411926632965     c",/,&
      & "   3.69180534781206    11.55414582295511    12.33193380846742     c",/,&
      & "   1.36791549262702     3.53066844289674    10.38660588677206     c",/,&
      & "   1.36792046664920     7.73723910626293    13.45809224934817     c",/,&
      & "   3.69181279359489     6.40826717723392    13.13552570942280     c",/,&
      & "   1.36792009865062     3.11669712338516     7.73723850632628     c",/,&
      & "   3.69181515738094     3.43926499914873     6.40826580885474     c",/,&
      & "   3.69178443989294     4.24285720771059    11.55415026712869     c",/,&
      & "   1.36790824853106     6.18818490375705     3.53066863732142     c",/,&
      & "   3.69178194163078     5.02063901427657     4.24285736953327     c",/,&
      & "   1.36794124909207    13.04411858182861     6.18818324080182     c",/,&
      & "   1.36792249732236     8.83755133592807     3.11669686076913     c",/,&
      & "   3.69182456413952    10.16652118921143     3.43926084011816     c",/,&
      & "   3.69181444966104    12.33193631088573     5.02063847821044     c",/,&
      & "   6.01572566324028    13.45790756713123     8.83752222635545     c",/,&
      & "   8.33965926123256    13.13576644753615    10.16660228658307     c",/,&
      & "   6.01574747573805    10.38654070512969    13.04391961251944     c",/,&
      & "   8.33964066450677    11.55427002850905    12.33211653730939     c",/,&
      & "   6.01574728097580     3.53087013230607    10.38654217813321     c",/,&
      & "   6.01568913853645     7.73726406411719    13.45790864082374     c",/,&
      & "   8.33963586549168     6.40818371470975    13.13576911116618     c",/,&
      & "   6.01568179676984     3.11688332536281     7.73726611148835     c",/,&
      & "   8.33963704688671     3.43902559351770     6.40818390180453     c",/,&
      & "   8.33962496288127     4.24267007149867    11.55427031066552     c",/,&
      & "   6.01573464280675     6.18824653544318     3.53086861480278     c",/,&
      & "   8.33961857277245     5.02052001792996     4.24267413625204     c",/,&
      & "   6.01575677304189    13.04392044501564     6.18824448603611     c",/,&
      & "   6.01568344836224     8.83752193432504     3.11688171781516     c",/,&
      & "   8.33964228963694    10.16660428027860     3.43902155668011     c",/,&
      & "   8.33965118613331    12.33211762632282     5.02051902430387     c",/,&
      & "$periodic 1",/,&
      & "$cell",/,&
      & " 9.29556285275863798006 ",/,&
      & "$end")'

   integer :: iunit
   type(tb_molecule) :: mol
!

   stop 77

   !write(iunit,file_coord_C_1d); rewind(iunit)
   !call read_coord(iunit,mol)
   !call mol%write(istdout,'C_1d')
   !call mol%deallocate; rewind(iunit)

end subroutine test_geometry_reader_file_coord_C_1D

