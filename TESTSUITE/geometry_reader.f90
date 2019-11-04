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

subroutine test_geometry_reader_molfile_benzen_flat
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-9_wp
   character(len=*),parameter :: file_mol_benzen_2d = '(/,&
      & "  Mrv1823 10191918342D          ",/,/,&
      & " 12 12  0  0  0  0            999 V2000",/,&
      & "    1.0868    1.0574    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    0.3723    1.4699    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    1.0868    0.2324    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "   -0.3421    1.0574    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    0.3723   -0.1801    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "   -0.3421    0.2324    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    1.8013    1.4699    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    0.3723    2.2949    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    1.8013   -0.1801    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "   -1.0566    1.4699    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    0.3723   -1.0051    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "   -1.0566   -0.1801    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "  2  1  2  0  0  0  0",/,&
      & "  3  1  1  0  0  0  0",/,&
      & "  4  2  1  0  0  0  0",/,&
      & "  5  3  2  0  0  0  0",/,&
      & "  6  4  2  0  0  0  0",/,&
      & "  6  5  1  0  0  0  0",/,&
      & "  1  7  1  0  0  0  0",/,&
      & "  2  8  1  0  0  0  0",/,&
      & "  3  9  1  0  0  0  0",/,&
      & "  4 10  1  0  0  0  0",/,&
      & "  5 11  1  0  0  0  0",/,&
      & "  6 12  1  0  0  0  0",/,&
      & "M  END")'
   integer :: iunit
   type(tb_molecule) :: mol

   open(newunit=iunit,status='scratch')
   write(iunit, file_mol_benzen_2d)
   rewind(iunit)
   call mol%read(iunit, format=p_ftype%molfile)

   call terminate(0) ! should fail
end subroutine

subroutine test_geometry_reader_molfile_benzen
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-9_wp
   character(len=*),parameter :: file_mol_benzen = '(/,&
      & "  Mrv1823 10191918163D          ",/,/,&
      & " 12 12  0  0  0  0            999 V2000",/,&
      & "   -0.0090   -0.0157   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "   -0.7131    1.2038   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    1.3990   -0.0157   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "   -0.0090    2.4232   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    2.1031    1.2038   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    1.3990    2.4232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "   -0.5203   -0.9011   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "   -1.7355    1.2038    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    1.9103   -0.9011    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "   -0.5203    3.3087    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    3.1255    1.2038    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    1.9103    3.3087   -0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "  2  1  4  0  0  0  0",/,&
      & "  3  1  4  0  0  0  0",/,&
      & "  4  2  4  0  0  0  0",/,&
      & "  5  3  4  0  0  0  0",/,&
      & "  6  4  4  0  0  0  0",/,&
      & "  6  5  4  0  0  0  0",/,&
      & "  1  7  1  0  0  0  0",/,&
      & "  2  8  1  0  0  0  0",/,&
      & "  3  9  1  0  0  0  0",/,&
      & "  4 10  1  0  0  0  0",/,&
      & "  5 11  1  0  0  0  0",/,&
      & "  6 12  1  0  0  0  0",/,&
      & "M  END")'
   integer :: iunit
   type(tb_molecule) :: mol

   open(newunit=iunit,status='scratch')
   write(iunit, file_mol_benzen)
   rewind(iunit)
   call mol%read(iunit, format=p_ftype%molfile)

   call assert_eq(len(mol), 12)
   call assert_eq(len(mol%bonds), 12)
   call assert_eq(len(mol%info), 0)

   call assert_close(mol%xyz(1,1), -0.17007533543675E-01_wp, thr)
   call assert_close(mol%xyz(2,8), 2.2748520977640_wp, thr)

   call terminate(afail)
end subroutine

subroutine test_geometry_reader_file_sdf_h2o
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-9_wp
   character(len=*),parameter :: file_sdf_h2o = '(&
      & "962",/,&
      & "  Marvin  12300703363D          ",/,&
      & "",/,&
      & "  3  2  0  0  0  0            999 V2000",/,&
      & "   -0.2309   -0.3265    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    0.7484   -0.2843    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "   -0.5175    0.6108    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "  1  2  1  0  0  0  0",/,&
      & "  1  3  1  0  0  0  0",/,&
      & "M  END",/,&
      & "",/,&
      & "> <StdInChI>",/,&
      & "InChI=1S/H2O/h1H2",/,&
      & "",/,&
      & "> <StdInChIKey>",/,&
      & "XLYOFNOQVPJJNP-UHFFFAOYSA-N",/,&
      & "",/,&
      & "> <AuxInfo>",/,&
      & "1/0/N:1/rA:3nOHH/rB:s1;s1;/rC:-.2309,-.3265,0;.7484,-.2843,0;-.5175,.6108,0;",/,&
      & "",/,&
      & "> <Formula>",/,&
      & "H2 O",/,&
      & "",/,&
      & "> <Mw>",/,&
      & "18.01528",/,&
      & "",/,&
      & "> <SMILES>",/,&
      & "O([H])[H]",/,&
      & "",/,&
      & "> <CSID>",/,&
      & "937",/,&
      & "",/,&
      & "$$$$")'
   integer :: iunit
   type(tb_molecule) :: mol

   open(newunit=iunit,status='scratch')
   write(iunit, file_sdf_h2o)
   rewind(iunit)
   call mol%read(iunit, format=p_ftype%sdf)

   call assert_eq(len(mol), 3)
   call assert_eq(len(mol%bonds), 2)
   call assert_eq(len(mol%info), 22)

   call assert_close(mol%xyz(1,1), -0.43633772169273_wp, thr)

   call terminate(afail)
end subroutine

subroutine test_geometry_reader_file_sdf_benzen_hquery
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-9_wp
   character(len=*),parameter :: file_sdf_benzen = '(&
      & "benzen",/,&
      & "  xtb     10301913453D",/,&
      & "",/,&
      & "  6  6  0  0  0  0  0  0  0  0999 V2000",/,&
      & "   -0.0090   -0.0157    0.0000 C   0  0  0  2  0  0  0  0  0  0  0  0",/,&
      & "   -0.7131    1.2038    0.0000 C   0  0  0  2  0  0  0  0  0  0  0  0",/,&
      & "    1.3990   -0.0157    0.0000 C   0  0  0  2  0  0  0  0  0  0  0  0",/,&
      & "   -0.0090    2.4232    0.0000 C   0  0  0  2  0  0  0  0  0  0  0  0",/,&
      & "    2.1031    1.2038    0.0000 C   0  0  0  2  0  0  0  0  0  0  0  0",/,&
      & "    1.3990    2.4232    0.0000 C   0  0  0  2  0  0  0  0  0  0  0  0",/,&
      & "  1  2  4  0  0  0  0",/,&
      & "  1  3  4  0  0  0  0",/,&
      & "  2  4  4  0  0  0  0",/,&
      & "  3  5  4  0  0  0  0",/,&
      & "  4  6  4  0  0  0  0",/,&
      & "  5  6  4  0  0  0  0",/,&
      & "M  END",/,&
      & "$$$$")'
   integer :: iunit
   type(tb_molecule) :: mol

   open(newunit=iunit,status='scratch')
   write(iunit, file_sdf_benzen)
   rewind(iunit)
   call mol%read(iunit, format=p_ftype%sdf)

   call terminate(0) ! should fail
end subroutine

subroutine test_geometry_reader_file_sdf_h2o_flat
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-9_wp
   character(len=*),parameter :: file_sdf_h2o = '(&
      & "water",/,&
      & "  xtb     10301913452D          ",/,&
      & "",/,&
      & "  3  2  0  0  0  0            999 V2000",/,&
      & "   -0.2309   -0.3265    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "    0.7484   -0.2843    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "   -0.5175    0.6108    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0",/,&
      & "  1  2  1  0  0  0  0",/,&
      & "  1  3  1  0  0  0  0",/,&
      & "M  END",/,&
      & "> <SMILES>",/,&
      & "O([H])[H]",/,&
      & "",/,&
      & "$$$$")'
   integer :: iunit
   type(tb_molecule) :: mol

   open(newunit=iunit,status='scratch')
   write(iunit, file_sdf_h2o)
   rewind(iunit)
   call mol%read(iunit, format=p_ftype%sdf)

   call terminate(0) ! should fail
end subroutine

subroutine test_geometry_reader_file_pdb_4qxx_noh
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-9_wp
   character(len=*),parameter :: file_pdb_4qxx_p1 = '(&
      & "HEADER    PROTEIN FIBRIL                          22-JUL-14   4QXX",/,&
      & "TITLE     STRUCTURE OF THE AMYLOID FORMING PEPTIDE GNLVS (RESIDUES 26-30) FROM",/,&
      & "TITLE    2 THE EOSINOPHIL MAJOR BASIC PROTEIN (EMBP)",/,&
      & "DBREF  4QXX Z    1     5  UNP    P13727   PRG2_HUMAN     131    135",/,&
      & "SEQRES   1 Z    5  GLY ASN LEU VAL SER",/,&
      & "FORMUL   2  HOH   *2(H2 O)",/,&
      & "CRYST1    4.755   16.816   35.759  90.00  90.00  90.00 P 2 21 21     4",/,&
      & "ORIGX1      1.000000  0.000000  0.000000        0.00000",/,&
      & "ORIGX2      0.000000  1.000000  0.000000        0.00000",/,&
      & "ORIGX3      0.000000  0.000000  1.000000        0.00000",/,&
      & "SCALE1      0.210305  0.000000  0.000000        0.00000",/,&
      & "SCALE2      0.000000  0.059467  0.000000        0.00000",/,&
      & "SCALE3      0.000000  0.000000  0.027965        0.00000",/,&
      & "ATOM      1  N   GLY Z   1      -0.821  -2.072  16.609  1.00  9.93           N",/,&
      & "ATOM      2  CA  GLY Z   1      -1.705  -2.345  15.487  1.00  7.38           C",/,&
      & "ATOM      3  C   GLY Z   1      -0.968  -3.008  14.344  1.00  4.89           C",/,&
      & "ATOM      4  O   GLY Z   1       0.258  -2.982  14.292  1.00  5.05           O",/,&
      & "ATOM      5  N   ASN Z   2      -1.721  -3.603  13.425  1.00  3.53           N",/,&
      & "ATOM      6  CA  ASN Z   2      -1.141  -4.323  12.291  1.00  1.85           C",/,&
      & "ATOM      7  C   ASN Z   2      -1.748  -3.900  10.968  1.00  3.00           C",/,&
      & "ATOM      8  O   ASN Z   2      -2.955  -3.683  10.873  1.00  3.99           O",/,&
      & "ATOM      9  CB  ASN Z   2      -1.353  -5.827  12.446  1.00  5.03           C",/,&
      & "ATOM     10  CG  ASN Z   2      -0.679  -6.391  13.683  1.00  5.08           C",/,&
      & "ATOM     11  OD1 ASN Z   2       0.519  -6.202  13.896  1.00  6.10           O",/,&
      & "ATOM     12  ND2 ASN Z   2      -1.448  -7.087  14.506  1.00  8.41           N")'
   character(len=*),parameter :: file_pdb_4qxx_p2 = '(&
      & "ATOM     13  N   LEU Z   3      -0.907  -3.803   9.944  1.00  3.47           N",/,&
      & "ATOM     14  CA  LEU Z   3      -1.388  -3.576   8.586  1.00  3.48           C",/,&
      & "ATOM     15  C   LEU Z   3      -0.783  -4.660   7.709  1.00  3.29           C",/,&
      & "ATOM     16  O   LEU Z   3       0.437  -4.788   7.643  1.00  3.80           O",/,&
      & "ATOM     17  CB  LEU Z   3      -0.977  -2.185   8.081  1.00  3.88           C",/,&
      & "ATOM     18  CG  LEU Z   3      -1.524  -1.669   6.736  1.00  8.66           C",/,&
      & "ATOM     19  CD1 LEU Z   3      -1.225  -0.191   6.570  1.00  9.89           C")'
   character(len=*),parameter :: file_pdb_4qxx_p3 = '(&
      & "ATOM     20  CD2 LEU Z   3      -0.962  -2.409   5.541  1.00 13.56           C",/,&
      & "ATOM     21  N   VAL Z   4      -1.635  -5.424   7.029  1.00  3.17           N",/,&
      & "ATOM     22  CA  VAL Z   4      -1.165  -6.460   6.119  1.00  3.61           C",/,&
      & "ATOM     23  C   VAL Z   4      -1.791  -6.230   4.755  1.00  5.31           C",/,&
      & "ATOM     24  O   VAL Z   4      -3.014  -6.209   4.620  1.00  7.31           O",/,&
      & "ATOM     25  CB  VAL Z   4      -1.567  -7.872   6.593  1.00  5.31           C",/,&
      & "ATOM     26  CG1 VAL Z   4      -1.012  -8.934   5.633  1.00  6.73           C",/,&
      & "ATOM     27  CG2 VAL Z   4      -1.083  -8.120   8.018  1.00  5.48           C",/,&
      & "ATOM     28  N   SER Z   5      -0.966  -6.052   3.736  1.00  7.53           N",/,&
      & "ATOM     29  CA  SER Z   5      -1.526  -5.888   2.407  1.00 11.48           C",/,&
      & "ATOM     30  C   SER Z   5      -1.207  -7.085   1.529  1.00 16.35           C",/,&
      & "ATOM     31  O   SER Z   5      -0.437  -7.976   1.902  1.00 14.00           O",/,&
      & "ATOM     32  CB  SER Z   5      -1.031  -4.596   1.767  1.00 13.36           C",/,&
      & "ATOM     33  OG  SER Z   5       0.361  -4.652   1.540  1.00 15.80           O",/,&
      & "ATOM     34  OXT SER Z   5      -1.737  -7.178   0.429  1.00 17.09           O",/,&
      & "TER      35      SER Z   5",/,&
      & "HETATM   36  O   HOH Z 101       0.935  -5.175  16.502  1.00 18.83           O",/,&
      & "HETATM   37  O  AHOH Z 102       0.691  -8.408  17.879  0.91 56.55           O",/,&
      & "HETATM   38  O  BHOH Z 102      -0.788  -9.006  16.641  0.09 38.95           O",/,&
      & "END")'
   integer :: iunit
   type(tb_molecule) :: mol

   open(newunit=iunit,status='scratch')
   write(iunit, file_pdb_4qxx_p1)
   write(iunit, file_pdb_4qxx_p2)
   write(iunit, file_pdb_4qxx_p3)
   rewind(iunit)
   call mol%read(iunit, format=p_ftype%pdb)

   call terminate(0) ! should fail
end subroutine


subroutine test_geometry_reader_file_pdb_4qxx
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use assertion
   use tbmod_file_utils
   real(wp),parameter :: thr = 1.0e-9_wp
   character(len=*),parameter :: file_pdb_4qxx_p1 = '(&
      & "HEADER    PROTEIN FIBRIL                          22-JUL-14   4QXX",/,&
      & "TITLE     STRUCTURE OF THE AMYLOID FORMING PEPTIDE GNLVS (RESIDUES 26-30) FROM",/,&
      & "TITLE    2 THE EOSINOPHIL MAJOR BASIC PROTEIN (EMBP)",/,&
      & "DBREF  4QXX Z    1     5  UNP    P13727   PRG2_HUMAN     131    135",/,&
      & "SEQRES   1 Z    5  GLY ASN LEU VAL SER",/,&
      & "FORMUL   2  HOH   *2(H2 O)",/,&
      & "CRYST1    4.755   16.816   35.759  90.00  90.00  90.00 P 2 21 21     4",/,&
      & "ORIGX1      1.000000  0.000000  0.000000        0.00000",/,&
      & "ORIGX2      0.000000  1.000000  0.000000        0.00000",/,&
      & "ORIGX3      0.000000  0.000000  1.000000        0.00000",/,&
      & "SCALE1      0.210305  0.000000  0.000000        0.00000",/,&
      & "SCALE2      0.000000  0.059467  0.000000        0.00000",/,&
      & "SCALE3      0.000000  0.000000  0.027965        0.00000",/,&
      & "ATOM      1  N   GLY Z   1      -0.821  -2.072  16.609  1.00  9.93           N",/,&
      & "ATOM      2  CA  GLY Z   1      -1.705  -2.345  15.487  1.00  7.38           C",/,&
      & "ATOM      3  C   GLY Z   1      -0.968  -3.008  14.344  1.00  4.89           C",/,&
      & "ATOM      4  O   GLY Z   1       0.258  -2.982  14.292  1.00  5.05           O",/,&
      & "ATOM      5  HA2 GLY Z   1      -2.130  -1.405  15.135  1.00  0.00           H",/,&
      & "ATOM      6  HA3 GLY Z   1      -2.511  -2.999  15.819  1.00  0.00           H",/,&
      & "ATOM      7  H1  GLY Z   1      -1.364  -1.742  17.394  1.00  0.00           H",/,&
      & "ATOM      8  H2  GLY Z   1      -0.150  -1.365  16.344  1.00  0.00           H",/,&
      & "ATOM      9  H3  GLY Z   1      -0.334  -2.918  16.868  1.00  0.00           H")'
   character(len=*),parameter :: file_pdb_4qxx_p2 = '(&
      & "ATOM     10  N   ASN Z   2      -1.721  -3.603  13.425  1.00  3.53           N",/,&
      & "ATOM     11  CA  ASN Z   2      -1.141  -4.323  12.291  1.00  1.85           C",/,&
      & "ATOM     12  C   ASN Z   2      -1.748  -3.900  10.968  1.00  3.00           C",/,&
      & "ATOM     13  O   ASN Z   2      -2.955  -3.683  10.873  1.00  3.99           O",/,&
      & "ATOM     14  CB  ASN Z   2      -1.353  -5.827  12.446  1.00  5.03           C",/,&
      & "ATOM     15  CG  ASN Z   2      -0.679  -6.391  13.683  1.00  5.08           C",/,&
      & "ATOM     16  OD1 ASN Z   2       0.519  -6.202  13.896  1.00  6.10           O",/,&
      & "ATOM     17  ND2 ASN Z   2      -1.448  -7.087  14.506  1.00  8.41           N",/,&
      & "ATOM     18  H   ASN Z   2      -2.726  -3.557  13.512  1.00  0.00           H",/,&
      & "ATOM     19  HA  ASN Z   2      -0.070  -4.123  12.263  1.00  0.00           H",/,&
      & "ATOM     20  HB2 ASN Z   2      -0.945  -6.328  11.568  1.00  0.00           H",/,&
      & "ATOM     21  HB3 ASN Z   2      -2.423  -6.029  12.503  1.00  0.00           H",/,&
      & "ATOM     22 HD21 ASN Z   2      -2.427  -7.218  14.293  1.00  0.00           H",/,&
      & "ATOM     23 HD22 ASN Z   2      -1.056  -7.487  15.346  1.00  0.00           H",/,&
      & "ATOM     24  N   LEU Z   3      -0.907  -3.803   9.944  1.00  3.47           N",/,&
      & "ATOM     25  CA  LEU Z   3      -1.388  -3.576   8.586  1.00  3.48           C",/,&
      & "ATOM     26  C   LEU Z   3      -0.783  -4.660   7.709  1.00  3.29           C",/,&
      & "ATOM     27  O   LEU Z   3       0.437  -4.788   7.643  1.00  3.80           O",/,&
      & "ATOM     28  CB  LEU Z   3      -0.977  -2.185   8.081  1.00  3.88           C",/,&
      & "ATOM     29  CG  LEU Z   3      -1.524  -1.669   6.736  1.00  8.66           C",/,&
      & "ATOM     30  CD1 LEU Z   3      -1.225  -0.191   6.570  1.00  9.89           C",/,&
      & "ATOM     31  CD2 LEU Z   3      -0.962  -2.409   5.541  1.00 13.56           C",/,&
      & "ATOM     32  H   LEU Z   3       0.086  -3.888  10.109  1.00  0.00           H",/,&
      & "ATOM     33  HA  LEU Z   3      -2.475  -3.661   8.568  1.00  0.00           H",/,&
      & "ATOM     34  HB2 LEU Z   3      -1.284  -1.469   8.843  1.00  0.00           H",/,&
      & "ATOM     35  HB3 LEU Z   3       0.111  -2.162   8.026  1.00  0.00           H",/,&
      & "ATOM     36  HG  LEU Z   3      -2.606  -1.798   6.737  1.00  0.00           H",/,&
      & "ATOM     37 HD11 LEU Z   3      -1.623   0.359   7.423  1.00  0.00           H",/,&
      & "ATOM     38 HD12 LEU Z   3      -1.691   0.173   5.654  1.00  0.00           H",/,&
      & "ATOM     39 HD13 LEU Z   3      -0.147  -0.043   6.513  1.00  0.00           H",/,&
      & "ATOM     40 HD21 LEU Z   3      -1.168  -3.475   5.643  1.00  0.00           H",/,&
      & "ATOM     41 HD22 LEU Z   3      -1.429  -2.035   4.630  1.00  0.00           H",/,&
      & "ATOM     42 HD23 LEU Z   3       0.115  -2.250   5.489  1.00  0.00           H")'
   character(len=*),parameter :: file_pdb_4qxx_p3 = '(&
      & "ATOM     43  N   VAL Z   4      -1.635  -5.424   7.029  1.00  3.17           N",/,&
      & "ATOM     44  CA  VAL Z   4      -1.165  -6.460   6.119  1.00  3.61           C",/,&
      & "ATOM     45  C   VAL Z   4      -1.791  -6.230   4.755  1.00  5.31           C",/,&
      & "ATOM     46  O   VAL Z   4      -3.014  -6.209   4.620  1.00  7.31           O",/,&
      & "ATOM     47  CB  VAL Z   4      -1.567  -7.872   6.593  1.00  5.31           C",/,&
      & "ATOM     48  CG1 VAL Z   4      -1.012  -8.934   5.633  1.00  6.73           C",/,&
      & "ATOM     49  CG2 VAL Z   4      -1.083  -8.120   8.018  1.00  5.48           C",/,&
      & "ATOM     50  H   VAL Z   4      -2.628  -5.282   7.146  1.00  0.00           H",/,&
      & "ATOM     51  HA  VAL Z   4      -0.080  -6.402   6.034  1.00  0.00           H",/,&
      & "ATOM     52  HB  VAL Z   4      -2.655  -7.939   6.585  1.00  0.00           H",/,&
      & "ATOM     53 HG11 VAL Z   4      -1.303  -9.926   5.980  1.00  0.00           H",/,&
      & "ATOM     54 HG12 VAL Z   4      -1.414  -8.766   4.634  1.00  0.00           H",/,&
      & "ATOM     55 HG13 VAL Z   4       0.075  -8.864   5.603  1.00  0.00           H",/,&
      & "ATOM     56 HG21 VAL Z   4      -1.377  -9.121   8.333  1.00  0.00           H",/,&
      & "ATOM     57 HG22 VAL Z   4       0.003  -8.032   8.053  1.00  0.00           H",/,&
      & "ATOM     58 HG23 VAL Z   4      -1.529  -7.383   8.686  1.00  0.00           H",/,&
      & "ATOM     59  N   SER Z   5      -0.966  -6.052   3.736  1.00  7.53           N",/,&
      & "ATOM     60  CA  SER Z   5      -1.526  -5.888   2.407  1.00 11.48           C",/,&
      & "ATOM     61  C   SER Z   5      -1.207  -7.085   1.529  1.00 16.35           C",/,&
      & "ATOM     62  O   SER Z   5      -0.437  -7.976   1.902  1.00 14.00           O",/,&
      & "ATOM     63  CB  SER Z   5      -1.031  -4.596   1.767  1.00 13.36           C",/,&
      & "ATOM     64  OG  SER Z   5       0.361  -4.652   1.540  1.00 15.80           O",/,&
      & "ATOM     65  OXT SER Z   5      -1.737  -7.178   0.429  1.00 17.09           O",/,&
      & "ATOM     66  H   SER Z   5       0.033  -6.031   3.880  1.00  0.00           H",/,&
      & "ATOM     67  HA  SER Z   5      -2.610  -5.822   2.504  1.00  0.00           H",/,&
      & "ATOM     68  HB2 SER Z   5      -1.543  -4.449   0.816  1.00  0.00           H",/,&
      & "ATOM     69  HB3 SER Z   5      -1.254  -3.759   2.428  1.00  0.00           H",/,&
      & "ATOM     70  HG  SER Z   5       0.653  -3.831   1.137  1.00  0.00           H",/,&
      & "TER      71      SER Z   5")'
   character(len=*),parameter :: file_pdb_4qxx_p4 = '(&
      & "HETATM   72  O   HOH Z 101       0.935  -5.175  16.502  1.00 18.83           O",/,&
      & "HETATM   73  H1  HOH Z 101       0.794  -5.522  15.621  1.00  0.00           H",/,&
      & "HETATM   74  H2  HOH Z 101       1.669  -4.561  16.489  1.00  0.00           H",/,&
      & "HETATM   75  O  AHOH Z 102       0.691  -8.408  17.879  0.91 56.55           O",/,&
      & "HETATM   76  O  BHOH Z 102      -0.788  -9.006  16.641  0.09 38.95           O",/,&
      & "HETATM   77  H1 AHOH Z 102       1.392  -8.125  18.466  0.91  0.00           H",/,&
      & "HETATM   78  H1 BHOH Z 102      -1.351  -9.776  16.563  0.09  0.00           H",/,&
      & "HETATM   79  H2 AHOH Z 102       0.993  -8.356  16.972  0.91  0.00           H",/,&
      & "HETATM   80  H2 BHOH Z 102      -0.927  -8.594  17.494  0.09  0.00           H",/,&
      & "END")'
   integer :: iunit
   type(tb_molecule) :: mol

   open(newunit=iunit,status='scratch')
   write(iunit, file_pdb_4qxx_p1)
   write(iunit, file_pdb_4qxx_p2)
   write(iunit, file_pdb_4qxx_p3)
   write(iunit, file_pdb_4qxx_p4)
   rewind(iunit)
   call mol%read(iunit, format=p_ftype%pdb)

   call assert(allocated(mol%pdb))

   call assert_eq(len(mol), 79)

   call assert_eq(mol%at(14), 6)
   call assert_eq(mol%at(42), 1)
   call assert_eq(mol%at(71), 8)

   call assert_close(mol%xyz(1,3),  -1.8292547189197_wp, thr)
   call assert_close(mol%xyz(2,21), -11.393157748313_wp, thr)
   call assert_close(mol%xyz(3,35),  15.166940469059_wp, thr)
   call assert_close(mol%xyz(1,79), -1.7517759549985_wp, thr)

   call terminate(afail) ! should fail
end subroutine
