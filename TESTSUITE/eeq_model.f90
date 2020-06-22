subroutine test_eeq_model_gbsa
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use assertion
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_param
   use xtb_solv_gbsa
   use xtb_solv_input
   use xtb_solv_model
   use xtb_solv_state
   use xtb_disp_coordinationnumber, only : getCoordinationNumber, cnType
   use xtb_eeq
   use xtb_chargemodel
   use xtb_solv_kernel
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
   integer, parameter :: iunit = stdout
   real(wp),parameter :: temp = 298.15_wp

   type(TMolecule) :: mol
   type(chrg_parameter) :: chrgeq
   type(TBorn) :: gbsa
   type(TSolvInput) :: input
   type(TSolvModel) :: model
   type(TEnvironment) :: env
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp) :: es,gsolv,sigma(3,3)
   real(wp),allocatable :: cn(:),dcndr(:,:,:),dcndL(:,:,:)
   real(wp),allocatable :: q(:),dqdr(:,:,:),dqdL(:,:,:),ges(:,:)

   allocate(cn(nat),dcndr(3,nat,nat),dcndL(3,3,nat), &
      & q(nat),dqdr(3,nat,nat),ges(3,nat))
   es  = 0.0_wp
   ges = 0.0_wp

   call init(env)
   call init(mol, at, xyz)

   call getCoordinationNumber(mol, trans, 40.0_wp, cnType%erf, cn, dcndr, dcndL)
   call assert_close(cn(1),3.9364408722599_wp,thr)
   call assert_close(cn(6),1.9734817313511_wp,thr)

   ! test CN derivative, correct summation of diagonals
   call assert_close(dcndr(2,1,5),-0.30966754386801E-01_wp,thr)
   call assert_close(dcndr(1,3,2),-0.00000000000000E+00_wp,thr)
   call assert_close(dcndr(2,1,6), 0.24290126624329E-05_wp,thr)

   call new_charge_model_2019(chrgeq,mol%n,mol%at)
   call eeq_chrgeq(mol,env,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
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

   input = TSolvInput(solvent='ch2cl2', alpb=.false., kernel=gbKernel%still, &
      & state=solutionState%gsolv)
   call init(model, env, input, 2)
   call newBornModel(model, env, gbsa, mol%at)
   call gbsa%update(env, mol%at, mol%xyz)

   es = gbsa%gsasa + gbsa%gshift
   ges = gbsa%dsdr
   call eeq_chrgeq(mol,env,chrgeq,gbsa,cn,dcndr,q,dqdr,es,gsolv,ges, &
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
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use assertion
   stop 77
   call terminate(afail)
end subroutine test_eeq_model_hbond

subroutine test_eeq_model_salt
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_convert, only : aatoau
   use assertion
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_param
   use xtb_solv_gbsa
   use xtb_solv_input
   use xtb_solv_model
   use xtb_solv_state
   use xtb_disp_coordinationnumber, only : getCoordinationNumber, cnType
   use xtb_eeq
   use xtb_chargemodel
   use xtb_solv_kernel
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
   integer, parameter :: iunit = stdout
   real(wp),parameter :: temp = 298.15_wp

   type(TMolecule) :: mol
   type(chrg_parameter) :: chrgeq
   type(TBorn) :: gbsa
   type(TSolvInput) :: input
   type(TSolvModel) :: model
   type(TEnvironment) :: env
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp) :: es,gsolv,sigma(3,3)
   real(wp),allocatable :: cn(:),dcndr(:,:,:),dcndL(:,:,:)
   real(wp),allocatable :: q(:),dqdr(:,:,:),dqdL(:,:,:),ges(:,:)
   real(wp),allocatable :: fgb(:,:),fhb(:)

   allocate( cn(nat),dcndr(3,nat,nat),dcndL(3,3,nat), &
      & q(nat),dqdr(3,nat,nat),ges(3,nat))
   es  = 0.0_wp
   ges = 0.0_wp

   call init(env)
   call init(mol, at, xyz)

   call getCoordinationNumber(mol, trans, 40.0_wp, cnType%erf, cn, dcndr, dcndL)

   call new_charge_model_2019(chrgeq,mol%n,mol%at)
   call eeq_chrgeq(mol,env,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
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

   input = TSolvInput(solvent='ch2cl2', alpb=.false., kernel=gbKernel%still, &
      & state=solutionState%gsolv, ionStrength=0.001_wp, ionRad=1.0_wp*aatoau)
   call init(model, env, input, 2)
   call newBornModel(model, env, gbsa, mol%at)
   call gbsa%update(env, mol%at, mol%xyz)

   es = gbsa%gsasa + gbsa%gshift
   ges = gbsa%dsdr
   call eeq_chrgeq(mol,env,chrgeq,gbsa,cn,dcndr,q,dqdr,es,gsolv,ges, &
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
