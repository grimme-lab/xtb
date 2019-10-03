!> test calculation of EEQ charges
subroutine test_eeq_model_water
   use iso_fortran_env, wp => real64
   use assertion
   use tbdef_molecule
   use tbdef_param
   use ncoord
   use eeq_model
   implicit none
   type(tb_molecule)    :: mol
   type(chrg_parameter) :: chrgeq

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.00000000000000_wp,  0.00000000000000_wp, -0.73578586109551_wp, &
      & 1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp, &
      &-1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp  &
      & ],shape(xyz))

   real(wp) :: es,sigma(3,3)
   real(wp),allocatable :: cn(:),dcndr(:,:,:),dcndL(:,:,:)
   real(wp),allocatable :: q(:),dqdr(:,:,:),dqdL(:,:,:),ges(:,:)

   allocate( cn(nat),dcndr(3,nat,nat),q(nat),dqdr(3,nat,nat+1),ges(3,nat) )
   es  = 0.0_wp
   ges = 0.0_wp

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp

   call dncoord_erf(mol%n,mol%at,mol%xyz,cn,dcndr)
   ! test for correct CN and correct symmetry in CN
   call assert_close(cn(1),1.9895544848535_wp,thr)
   call assert_close(cn(2),cn(3),             thr)
   call assert_close(cn(1),cn(2)+cn(3),       thr)

   ! test CN derivative, correct summation of diagonals
   call assert_close(dcndr(3,1,1),-0.80973198569003E-01_wp,thr)
   call assert_close(dcndr(1,2,1), 0.52891163425093E-01_wp,thr)
   call assert_close(dcndr(3,1,3),-0.40486599284501E-01_wp,thr)

   call new_charge_model_2019(chrgeq,mol%n,mol%at)
   call eeq_chrgeq(mol,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
      &            .false.,.true.,.true.)

   ! test electrostatic energy
   call assert_close(es,-0.64308088326667E-01_wp,thr)

   ! test molecular gradient of ES, also check symmetry
   call assert_close(ges(3,1),-0.44053032330503E-01_wp,thr)
   call assert_close(ges(3,2),ges(3,3),                thr)
   call assert_close(ges(1,2), 0.18102071270235E-01_wp,thr)

   ! test for charge constraint
   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1),-0.59582744684480_wp,thr)

   ! test derivative of partial charges
   call assert_close(dqdr(3,1,1),-0.41466764014389E+00_wp,thr)
   call assert_close(dqdr(1,2,1), 0.17196993288918E+00_wp,thr)
   call assert_close(dqdr(1,3,2), 0.18631993794394E-01_wp,thr)
   call assert_close(dqdr(2,1,3), 0.00000000000000E+00_wp,thr)

   call mol%deallocate

   ! done: everythings fine
   call terminate(afail)
end subroutine test_eeq_model_water

subroutine test_eeq_model_ewald
   use iso_fortran_env, wp => real64
   use assertion
   use tbdef_molecule
   use tbdef_param
   use eeq_model
   use ncoord
   use pbc_tools
   implicit none
   real(wp),parameter :: thr = 1.0e-9_wp
   integer, parameter :: nat = 6
   integer, parameter :: at(nat) = [14,14,8,8,8,8]
   real(wp),parameter :: abc(3,nat) = reshape(&
      &[.095985472469032_wp, .049722204206931_wp, 0.10160624337938_wp, &
      & 0.54722204206931_wp, 0.52863628207623_wp, 0.38664208660311_wp, &
      & 0.29843937068984_wp, 0.39572194413818_wp, 0.20321248675876_wp, &
      & 0.23364982659922_wp, 0.85647058758674_wp, 0.31884968761485_wp, &
      & 0.72250232459952_wp, 0.65548544066844_wp, .056207709103487_wp, &
      & 0.70514214000043_wp, 0.28321754549582_wp, 0.36424822189074_wp],&
      & shape(abc))
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[ 8.7413053236641_wp,  0.0000000000000_wp,  0.0000000000000_wp,   &
      &  0.0000000000000_wp,  8.7413053236641_wp,  0.0000000000000_wp,   &
      &  0.0000000000000_wp,  0.0000000000000_wp,  8.7413053236641_wp],  &
      & shape(lattice))
   integer, parameter :: wsc_rep(3) = [1,1,1]
   real(wp),parameter :: beta = 7.5_wp

   type(tb_molecule)    :: mol
   type(chrg_parameter) :: chrgeq
   real(wp)             :: energy
   real(wp)             :: sigma(3,3)
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: dcndr(:,:,:)
   real(wp),allocatable :: dcndL(:,:,:)
   real(wp),allocatable :: q(:)
   real(wp),allocatable :: dqdr(:,:,:)
   real(wp),allocatable :: dqdL(:,:,:)
   real(wp),allocatable :: gradient(:,:)
   integer              :: rep_cn(3)

   call mol%allocate(nat)
   mol%at   = at
   mol%abc  = abc
   mol%npbc = 3
   mol%pbc  = .true.
   mol%lattice = lattice
   mol%volume = dlat_to_dvol(lattice)
   call dlat_to_cell(lattice,mol%cellpar)
   call dlat_to_rlat(lattice,mol%rec_lat)
   call coord_trafo(nat,lattice,abc,mol%xyz)
   call mol%wrap_back

   allocate( cn(nat), dcndr(3,nat,nat), dcndL(3,3,nat), &
      &      q(nat), dqdr(3,nat,nat+1), dqdL(3,3,nat+1), &
      &      gradient(3,nat), source = 0.0_wp )
   cn    = 0.0_wp
   dcndr = 0.0_wp
   dcndL = 0.0_wp
   q     = 0.0_wp
   dqdr  = 0.0_wp
   dqdL  = 0.0_wp
   gradient = 0.0_wp
   energy = 0.0_wp
   sigma = 0.0_wp

   call generate_wsc(mol,mol%wsc,wsc_rep)

   call pbc_derfcoord(mol%n,mol%at,mol%xyz,mol%lattice,cn,dcndr,dcndL,900.0_wp)
   !call derfsum(mol%n,mol%at,mol%xyz,mol%lattice,cn,cnp,dcnpdr,dcnpdL, &
      !&         beta,rep_cn,900.0_wp)

   call assert_close(cn(2),3.6864725130236_wp,thr)
   call assert_close(cn(5),1.0523558225297_wp,thr)
   call assert_close(cn(6),1.1699488478421_wp,thr)

   call assert_close(dcndr(2,3,2),-0.24281192795725E-02_wp,thr)
   call assert_close(dcndr(1,6,6),-0.42240789876965E+00_wp,thr)
   call assert_close(dcndr(1,1,2),-0.71913132896636E-01_wp,thr)

   call assert_close(dcndL(1,3,4),-0.36068395425405_wp,thr)
   call assert_close(dcndL(3,3,1),-0.92808088688242_wp,thr)
   call assert_close(dcndL(2,1,3),-0.44871260782914_wp,thr)

   call new_charge_model_2019(chrgeq,mol%n,mol%at)

   call eeq_chrgeq(mol,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,energy,gradient,sigma,&
      &            .false.,.true.,.true.)

   call assert_close(energy,-0.90576568382295E-01_wp,thr)

   call assert_close(gradient(2,1), 0.59224535975576E-02_wp,thr)
   call assert_close(gradient(1,3),-0.17291053212435E-02_wp,thr)
   call assert_close(gradient(3,5), 0.13109339409759E-02_wp,thr)
   call assert_close(gradient(1,6), 0.10491055487530E-01_wp,thr)

   call assert_close(sigma(1,1),-0.69556402379101E-01_wp,thr)
   call assert_close(sigma(2,3),-0.11063392011917E-01_wp,thr)
   call assert_close(sigma(3,1), 0.15351359672016E-01_wp,thr)

   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1), 0.39808668429315_wp,thr)
   call assert_close(q(3),-0.16133275860565_wp,thr)
   call assert_close(q(4),-0.19939812633388_wp,thr)

   call assert_close(dqdr(1,4,2),-0.16619680140869E-01_wp,thr)
   call assert_close(dqdr(2,5,7),-0.71528076609028E-04_wp,thr)
   call assert_close(dqdr(1,2,2), 0.36177437391610E-01_wp,thr)
   call assert_close(dqdr(3,1,4),-0.39909176876716E-01_wp,thr)

   call assert_close(dqdL(2,3,2),-0.66131872796832E-02_wp,thr)
   call assert_close(dqdL(2,1,7), 0.66206886462644E-02_wp,thr)
   call assert_close(dqdL(3,2,1), 0.51752004657789E-01_wp,thr)
   call assert_close(dqdL(1,1,5),-0.29765138003366E-01_wp,thr)

   ! reset energy, gradient and stress
   energy = 0.0_wp
   gradient = 0.0_wp
   sigma = 0.0_wp

   call dncoord_logcn(mol%n,cn,dcndr,dcndL,cn_max=8.0_wp)

   call assert_close(cn(2),3.6735104771649_wp,thr)
   call assert_close(cn(5),1.0517307940776_wp,thr)
   call assert_close(cn(6),1.1692040350311_wp,thr)

   call assert_close(dcndr(2,3,2),-0.23960452277676E-02_wp,thr)
   call assert_close(dcndr(1,6,6),-0.42195185201461E+00_wp,thr)
   call assert_close(dcndr(1,1,2),-0.70963201989459E-01_wp,thr)

   call assert_close(dcndL(1,3,4),-0.36002844990181_wp,thr)
   call assert_close(dcndL(3,3,1),-0.92527392943064_wp,thr)
   call assert_close(dcndL(2,1,3),-0.44767072306065_wp,thr)


   call eeq_chrgeq(mol,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,energy,gradient,sigma,&
      &            .false.,.true.,.true.)

   call assert_close(energy,-0.90708295779596E-01_wp,thr)

   call assert_close(gradient(2,1), 0.59151782557369E-02_wp,thr)
   call assert_close(gradient(1,3),-0.17299335354888E-02_wp,thr)
   call assert_close(gradient(3,5), 0.12993726124624E-02_wp,thr)
   call assert_close(gradient(1,6), 0.10488641464264E-01_wp,thr)

   call assert_close(sigma(1,1),-0.69363778247591E-01_wp,thr)
   call assert_close(sigma(2,3),-0.11069891635844E-01_wp,thr)
   call assert_close(sigma(3,1), 0.15371056800913E-01_wp,thr)

   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1), 0.39826062215939_wp,thr)
   call assert_close(q(3),-0.16154676123047_wp,thr)
   call assert_close(q(4),-0.19955653474773_wp,thr)

   call assert_close(dqdr(1,4,2),-0.16414391687062E-01_wp,thr)
   call assert_close(dqdr(2,5,7),-0.72381798942512E-04_wp,thr)
   call assert_close(dqdr(1,2,2), 0.35974299411392E-01_wp,thr)
   call assert_close(dqdr(3,1,4),-0.39939521246564E-01_wp,thr)

   call assert_close(dqdL(2,3,2),-0.65259046811329E-02_wp,thr)
   call assert_close(dqdL(2,1,7), 0.66316093331740E-02_wp,thr)
   call assert_close(dqdL(3,2,1), 0.51675870143073E-01_wp,thr)
   call assert_close(dqdL(1,1,5),-0.29605046716685E-01_wp,thr)

   call mol%deallocate

   ! done
   call terminate(afail)

end subroutine test_eeq_model_ewald

subroutine test_eeq_model_gbsa
   use iso_fortran_env, wp => real64, istderr => error_unit
   use assertion
   use tbdef_molecule
   use tbdef_param
   use gbobc
   use ncoord
   use eeq_model
   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 8
   integer, parameter :: at(nat) = [82,1,1,1,1,52,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[0.00607452450000_wp,-4.06569893412500_wp, 0.00000000000000_wp, &
      &-0.10304014450000_wp,-0.81095583012500_wp, 0.00000000000000_wp, &
      &-1.48984957750000_wp,-5.20124851012500_wp, 2.65857423800000_wp, &
      &-1.48984957750000_wp,-5.20124851012500_wp,-2.65857423800000_wp, &
      & 3.11171824650000_wp,-5.04151186112500_wp, 0.00000000000000_wp, &
      & 0.00607452450000_wp, 5.31621802187500_wp, 0.00000000000000_wp, &
      & 2.16097631450000_wp, 7.52891897087500_wp, 0.00000000000000_wp, &
      &-2.20210431050000_wp, 7.47552665287500_wp, 0.00000000000000_wp], shape(xyz))
   integer, parameter :: iunit = output_unit
   real(wp),parameter :: temp = 298.15_wp

   type(tb_molecule)    :: mol
   type(chrg_parameter) :: chrgeq
   type(tb_solvent)     :: gbsa
   real(wp) :: es,gsolv,sigma(3,3)
   real(wp),allocatable :: cn(:),dcndr(:,:,:),dcndL(:,:,:)
   real(wp),allocatable :: q(:),dqdr(:,:,:),dqdL(:,:,:),ges(:,:)
   real(wp),allocatable :: fgb(:,:),fhb(:)

   allocate( cn(nat),dcndr(3,nat,nat),q(nat),dqdr(3,nat,nat+1),ges(3,nat), &
      &      fgb(nat,nat),fhb(nat) )
   es  = 0.0_wp
   ges = 0.0_wp

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp

   call dncoord_erf(mol%n,mol%at,mol%xyz,cn,dcndr)
   call assert_close(cn(1),3.9364408722599_wp,thr)
   call assert_close(cn(6),1.9734817313511_wp,thr)

   ! test CN derivative, correct summation of diagonals
   call assert_close(dcndr(2,1,5), 0.30966754386801E-01_wp,thr)
   call assert_close(dcndr(1,3,2), 0.00000000000000E+00_wp,thr)
   call assert_close(dcndr(2,1,6),-0.24290126624329E-05_wp,thr)

   call new_charge_model_2019(chrgeq,mol%n,mol%at)
   call eeq_chrgeq(mol,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
      &            .false.,.true.,.true.)

   ! test electrostatic energy
   call assert_close(es,-0.10559262715710E-01_wp,thr)

   ! test molecular gradient of ES, also check symmetry
   call assert_close(ges(2,5), 0.45485199335534E-04_wp,thr)
   call assert_close(ges(3,3),-0.62757274919919E-04_wp,thr)
   call assert_close(ges(1,2), 0.70062774412852E-05_wp,thr)

   ! test for charge constraint
   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1),0.21873929670215_wp,thr)

   ! test derivative of partial charges
   call assert_close(dqdr(2,1,6),-0.32458122507411E-02_wp,thr)
   call assert_close(dqdr(3,2,1), 0.00000000000000E+00_wp,thr)
   call assert_close(dqdr(1,3,2), 0.27025784902288E-03_wp,thr)
   call assert_close(dqdr(2,8,5), 0.36651303109315E-03_wp,thr)

   es = 0.0_wp
   ges = 0.0_wp
   q = 0.0_wp
   dqdr = 0.0_wp

   lgbsa = .true.
   call init_gbsa(iunit,'ch2cl2',0,temp,2,230)
   call new_gbsa(gbsa,mol%n,mol%at)
   call update_nnlist_gbsa(gbsa,mol%xyz,.false.)
   call compute_brad_sasa(gbsa,mol%xyz)

   es = gbsa%gsasa
   ges = gbsa%dsdr
   call eeq_chrgeq(mol,chrgeq,gbsa,cn,dcndr,q,dqdr,es,gsolv,ges, &
      &            .false.,.true.,.true.)

   ! test electrostatic energy
   call assert_close(es,-0.18723100944150E-01_wp,thr)

   ! test molecular gradient of ES, also check symmetry
   call assert_close(ges(2,5), 0.15098175186158E-03_wp,thr)
   call assert_close(ges(3,3),-0.42131892534669E-03_wp,thr)
   call assert_close(ges(1,2), 0.28893943115425E-04_wp,thr)

   ! test for charge constraint
   call assert_close(sum(q),0.0_wp,thr)
   call assert_close(q(1), 0.22653931001924E+00_wp,thr)
   call assert_close(q(2),-0.20322834228298E-01_wp,thr)
   call assert_close(q(3),-0.26032042423717E-01_wp,thr)
   call assert_close(q(4),-0.26032042423717E-01_wp,thr)
   call assert_close(q(5),-0.25976606121313E-01_wp,thr)
   call assert_close(q(6),-0.18580233359855E+00_wp,thr)
   call assert_close(q(7), 0.28815815576667E-01_wp,thr)
   call assert_close(q(8), 0.28810733199693E-01_wp,thr)

   !call assert_close(dqdr(2,1,6),-0.21098409244566E-02_wp,thr)
   !call assert_close(dqdr(3,2,1), 0.00000000000000E+00_wp,thr)
   !call assert_close(dqdr(1,3,2), 0.37653484040321E-03_wp,thr)
   !call assert_close(dqdr(2,8,5), 0.82628766939384E-03_wp,thr)

   call mol%deallocate

   ! done: everythings fine
   call terminate(afail)

end subroutine test_eeq_model_gbsa

subroutine test_eeq_model_hbond
   use iso_fortran_env, wp => real64
   use assertion
   stop 77
   call terminate(afail)
end subroutine test_eeq_model_hbond

subroutine test_eeq_model_salt
   use iso_fortran_env, wp => real64, istderr => error_unit
   use assertion
   use tbdef_molecule
   use tbdef_param
   use gbobc
   use ncoord
   use eeq_model
   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 12
   integer, parameter :: at(nat) = [9,1,6,6,6,6,6,6,9,1,1,9]
   real(wp),parameter :: xyz(3,nat) = reshape([&
      &-3.61707571665846_wp,-2.21326241566733_wp,-0.20914111039455_wp, &
      &-4.12085912923274_wp, 2.64361027892242_wp,-0.77360733295900_wp, &
      & 0.70571179547784_wp,-1.46804568562192_wp, 0.40260319330016_wp, &
      &-1.70284271199247_wp,-0.55881754958407_wp,-0.07245619686109_wp, &
      &-2.22721391221737_wp, 1.98624548249273_wp,-0.40689082126896_wp, &
      &-0.21016646823535_wp, 3.64828289401216_wp,-0.24136575031326_wp, &
      & 2.24393351740138_wp, 2.86731529508407_wp, 0.23334807696843_wp, &
      & 2.63467443489279_wp, 0.29590222227677_wp, 0.54871426109542_wp, &
      &-0.65510226134586_wp, 6.12711283170949_wp,-0.54503487549757_wp, &
      & 3.78023043231610_wp, 4.20058933929303_wp, 0.35763959771486_wp, &
      & 1.05737533649814_wp,-3.45504988937755_wp, 0.68084171455529_wp, &
      & 4.98699896397970_wp,-0.51782531551223_wp, 1.02289296328971_wp], shape(xyz))
   integer, parameter :: iunit = output_unit
   real(wp),parameter :: temp = 298.15_wp

   type(tb_molecule)    :: mol
   type(chrg_parameter) :: chrgeq
   type(tb_solvent)     :: gbsa
   real(wp) :: es,gsolv,sigma(3,3)
   real(wp),allocatable :: cn(:),dcndr(:,:,:),dcndL(:,:,:)
   real(wp),allocatable :: q(:),dqdr(:,:,:),dqdL(:,:,:),ges(:,:)
   real(wp),allocatable :: fgb(:,:),fhb(:)

   allocate( cn(nat),dcndr(3,nat,nat),q(nat),dqdr(3,nat,nat+1),ges(3,nat), &
      &      fgb(nat,nat),fhb(nat) )
   es  = 0.0_wp
   ges = 0.0_wp

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp

   call dncoord_erf(mol%n,mol%at,mol%xyz,cn,dcndr)

   call new_charge_model_2019(chrgeq,mol%n,mol%at)
   call eeq_chrgeq(mol,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
      &            .false.,.true.,.true.)

   ! test electrostatic energy
   call assert_close(es,-0.64149462344256E-01_wp,thr)

   ! test molecular gradient of ES, also check symmetry
   call assert_close(ges(2,7),-0.52237685114757E-03_wp,thr)
   call assert_close(ges(3,3),-0.10347464904462E-03_wp,thr)
   call assert_close(ges(1,9),-0.19202182317373E-03_wp,thr)

   ! test for charge constraint
   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1),-0.20644568417986E+00_wp,thr)
   call assert_close(q(8), 0.83875312119020E-01_wp,thr)

   ! test derivative of partial charges
   call assert_close(dqdr(2,2,6),-0.93834552520241E-02_wp,thr)
   call assert_close(dqdr(3,3,1),-0.24876158910325E-03_wp,thr)
   call assert_close(dqdr(1,5,2),-0.15306250365591E-01_wp,thr)
   call assert_close(dqdr(2,8,7),-0.92052318691156E-02_wp,thr)

   es = 0.0_wp
   ges = 0.0_wp
   q = 0.0_wp
   dqdr = 0.0_wp

   lgbsa = .true.
   lsalt = .true.
   ionst = 0.001_wp
   ion_rad = 1.0_wp
   call init_gbsa(iunit,'ch2cl2',0,temp,2,230)
   call new_gbsa(gbsa,mol%n,mol%at)
   call update_nnlist_gbsa(gbsa,mol%xyz,.false.)
   call compute_brad_sasa(gbsa,mol%xyz)

   es = gbsa%gsasa
   ges = gbsa%dsdr
   call eeq_chrgeq(mol,chrgeq,gbsa,cn,dcndr,q,dqdr,es,gsolv,ges, &
      &            .false.,.true.,.true.)

   ! test electrostatic energy
   call assert_close(es,-0.73041696453862E-01_wp,thr)

   ! test molecular gradient of ES, also check symmetry
   call assert_close(ges(2,7),-0.54144355306156E-03_wp,thr)
   call assert_close(ges(3,3),-0.13994375398115E-03_wp,thr)
   call assert_close(ges(1,9),-0.72309374268259E-04_wp,thr)

   ! test for charge constraint
   call assert_close(sum(q),0.0_wp,thr)
   call assert_close(q( 1),-0.22628341725089E+00_wp,thr)
   call assert_close(q( 2), 0.17181969742365E+00_wp,thr)
   call assert_close(q( 3),-0.31688852764536E-01_wp,thr)
   call assert_close(q( 4), 0.86131431150190E-01_wp,thr)
   call assert_close(q( 5),-0.31623215923196E-01_wp,thr)
   call assert_close(q( 6), 0.86034414210273E-01_wp,thr)
   call assert_close(q( 7),-0.31614086140679E-01_wp,thr)
   call assert_close(q( 8), 0.86132083253082E-01_wp,thr)
   call assert_close(q( 9),-0.22624788057207E+00_wp,thr)
   call assert_close(q(10), 0.17181804147276E+00_wp,thr)
   call assert_close(q(11), 0.17180510824691E+00_wp,thr)
   call assert_close(q(12),-0.22628332310549E+00_wp,thr)

   ! test derivative of partial charges (still incorrect)
   !call assert_close(dqdr(2,2,6),-0.88602533712134E-02_wp,thr)
   !call assert_close(dqdr(3,3,1),-0.38467542758235E-03_wp,thr)
   !call assert_close(dqdr(1,5,2),-0.13407102470387E-01_wp,thr)
   !call assert_close(dqdr(2,8,7),-0.84239153450647E-02_wp,thr)

   call mol%deallocate

   call terminate(afail)
end subroutine test_eeq_model_salt
