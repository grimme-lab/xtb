!> test calculation of dispersion related properties
subroutine test_dftd4_properties
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use assertion
   use xtb_type_molecule
   use xtb_disp_dftd4
   implicit none
   type(TMolecule)       :: mol
   integer              :: ndim
   real(wp) :: molpol,molc6,molc8        ! molecular Polarizibility
   real(wp),allocatable :: gweights(:)   ! gaussian weights
   real(wp),allocatable :: refc6(:,:)    ! reference C6 coeffients
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: aw(:,:)

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.00000000000000_wp, 0.00000000000000_wp,-0.73578586109551_wp, &
      & 1.44183152868459_wp, 0.00000000000000_wp, 0.36789293054775_wp, &
      &-1.44183152868459_wp, 0.00000000000000_wp, 0.36789293054775_wp  &
      & ],shape(xyz))
   real(wp),parameter :: covcn(nat) = &
      [ 1.6105486019977_wp,  0.80527430099886_wp, 0.80527430099886_wp]
   real(wp),parameter :: q(nat) = &
      [-0.59582744708873_wp, 0.29791372354436_wp, 0.29791372354436_wp]
   real(wp),parameter :: g_a = 3.0_wp
   real(wp),parameter :: g_c = 2.0_wp
   real(wp),parameter :: wf  = 6.0_wp
   integer, parameter :: lmbd = p_mbd_approx_atm
   integer, parameter :: refqmode = p_refq_goedecker

   call init(mol, at, xyz)

   call d4init(g_a,g_c,refqmode)
   call d4dim(mol%n,mol%at,ndim)

   call assert_eq(ndim,8)

   allocate( gweights(ndim),refc6(ndim,ndim),&
             c6ab(mol%n,mol%n),aw(23,mol%n) )

   call d4(mol%n,ndim,mol%at,wf,g_a,g_c,covcn,gweights,refc6)
   call mdisp(mol%n,ndim,mol%at,q,mol%xyz,g_a,g_c,gweights,refc6, &
      &       molc6,molc8,molpol,aw,c6ab)

   call assert_close(molpol,9.4271529762816_wp,thr)
   call assert_close(molc6, 44.521546516541_wp,thr)
   call assert_close(molc8, 798.69642220617_wp,thr)

   call assert_close(aw(1,1),6.7482856791776_wp,thr)
   call assert_close(aw(4,2),1.1637689932984_wp,thr)
   call assert_close(aw(7,2),aw(7,3),           thr)

   call assert_close(c6ab(1,2),c6ab(2,1),         thr)
   call assert_close(c6ab(1,1),24.900853294042_wp,thr)
   call assert_close(c6ab(1,3),4.1779699081826_wp,thr)
   call assert_close(c6ab(2,2),c6ab(2,3),         thr)

   call assert_close(sum(gweights),3.0_wp,               thr)
   call assert_close(gweights(2),0.18388891750767E-01_wp,thr)
   call assert_close(gweights(7),0.21400765778468E-01_wp,thr)

   call assert_close(refc6(5,1),10.282421843343_wp,thr)
   call assert_close(refc6(8,6),3.0374149547985_wp,thr)

   call mol%deallocate

   ! done: everythings fine
   call terminate(afail)
end subroutine test_dftd4_properties

!> @brief test calculation of dispersion energies
subroutine test_dftd4_energies
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use assertion
   use xtb_type_molecule
   use xtb_type_param
   use xtb_disp_dftd4
   implicit none
   type(TMolecule)       :: mol
   integer  :: idum
   real(wp) :: energy

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.00000000000000_wp, 0.00000000000000_wp,-0.73578586109551_wp, &
      & 1.44183152868459_wp, 0.00000000000000_wp, 0.36789293054775_wp, &
      &-1.44183152868459_wp, 0.00000000000000_wp, 0.36789293054775_wp  &
      & ],shape(xyz))
   real(wp),parameter :: covcn(nat) = &
      [ 1.6105486019977_wp,  0.80527430099886_wp, 0.80527430099886_wp]
   real(wp),parameter :: q(nat) = &
      [-0.59582744708873_wp, 0.29791372354436_wp, 0.29791372354436_wp]
   real(wp),parameter :: g_a = 3.0_wp
   real(wp),parameter :: g_c = 2.0_wp
   real(wp),parameter :: wf  = 6.0_wp
   integer, parameter :: lmbd = p_mbd_approx_atm
   integer, parameter :: refqmode = p_refq_goedecker
   integer, parameter :: ndim = 8
   real(wp),parameter :: gweights(ndim) = &
      [ 0.15526686926080E-06_wp, 0.18388886232894E-01_wp, 0.89143504504233_wp, &
      & 0.90175913457907E-01_wp, 0.21400765336381E-01_wp, 0.97859923466362_wp, &
      & 0.21400765336381E-01_wp, 0.97859923466362_wp ]
   real(wp),parameter :: refc6(ndim,ndim) = reshape(&
      [ 0.0000000000000_wp,      0.0000000000000_wp,      0.0000000000000_wp,  &
      & 0.0000000000000_wp,      10.282422024500_wp,      6.7431228212696_wp,  &
      & 10.282422024500_wp,      6.7431228212696_wp,      0.0000000000000_wp,  &
      & 0.0000000000000_wp,      0.0000000000000_wp,      0.0000000000000_wp,  &
      & 12.052429296454_wp,      7.8894703511335_wp,      12.052429296454_wp,  &
      & 7.8894703511335_wp,      0.0000000000000_wp,      0.0000000000000_wp,  &
      & 0.0000000000000_wp,      0.0000000000000_wp,      13.246161891965_wp,  &
      & 8.6635841400632_wp,      13.246161891965_wp,      8.6635841400632_wp,  &
      & 0.0000000000000_wp,      0.0000000000000_wp,      0.0000000000000_wp,  &
      & 0.0000000000000_wp,      10.100325850238_wp,      6.6163452797181_wp,  &
      & 10.100325850238_wp,      6.6163452797181_wp,      10.282422024500_wp,  &
      & 12.052429296454_wp,      13.246161891965_wp,      10.100325850238_wp,  &
      & 0.0000000000000_wp,      0.0000000000000_wp,      7.6362416262742_wp,  &
      & 4.7593057612608_wp,      6.7431228212696_wp,      7.8894703511335_wp,  &
      & 8.6635841400632_wp,      6.6163452797181_wp,      0.0000000000000_wp,  &
      & 0.0000000000000_wp,      4.7593057612608_wp,      3.0374149102818_wp,  &
      & 10.282422024500_wp,      12.052429296454_wp,      13.246161891965_wp,  &
      & 10.100325850238_wp,      7.6362416262742_wp,      4.7593057612608_wp,  &
      & 0.0000000000000_wp,      0.0000000000000_wp,      6.7431228212696_wp,  &
      & 7.8894703511335_wp,      8.6635841400632_wp,      6.6163452797181_wp,  &
      & 4.7593057612608_wp,      3.0374149102818_wp,      0.0000000000000_wp,  &
      & 0.0000000000000_wp],     shape(refc6))
   type(dftd_parameter),parameter :: dparam_pwpb95 = dftd_parameter ( &
      &  s6=0.8200_wp, s8=-0.34639127_wp, a1=0.41080636_wp, a2=3.83878274_wp )
   type(dftd_parameter),parameter :: dparam_pbe    = dftd_parameter ( &
      &  s6=1.0000_wp, s8=0.95948085_wp, a1=0.38574991_wp, a2=4.80688534_wp )
   type(dftd_parameter),parameter :: dparam_random = dftd_parameter ( &
      &  s6=0.95_wp, s8=0.45_wp, s10=0.65_wp, s9=1.10_wp, a1=0.43_wp, a2=5.10_wp )

   call init(mol, at, xyz)

   call d4init(g_a,g_c,refqmode)
   call d4dim(mol%n,mol%at,idum)

   call assert_eq(idum,ndim)

   energy = +1.0_wp ! energy is intent(out)

   call edisp(mol%n,ndim,mol%at,q,mol%xyz,dparam_pwpb95,g_a,g_c, &
      &       gweights,refc6,lmbd,energy)
   call assert_close(energy,-0.22526819184723E-03_wp,thr)

   call edisp(mol%n,ndim,mol%at,q,mol%xyz,dparam_pbe,g_a,g_c, &
      &       gweights,refc6,lmbd,energy)
   call assert_close(energy,-0.19788865790096E-03_wp,thr)
   !-0.19558245089408E-03

   call edisp(mol%n,ndim,mol%at,q,mol%xyz,dparam_random,g_a,g_c, &
      &       gweights,refc6,lmbd,energy)
   call assert_close(energy,-0.11213581758666E-03_wp,thr)

   call mol%deallocate

   ! done: everythings fine
   call terminate(afail)
end subroutine test_dftd4_energies

subroutine test_dftd4_pbc_energies
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use assertion
   use xtb_type_molecule
   use xtb_type_param
   use xtb_disp_dftd4
   use xtb_eeq
   use xtb_chargemodel
   use xtb_disp_ncoord
   use xtb_pbc_tools
   implicit none
   real(wp),parameter :: thr = 1.0e-10_wp
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
   real(wp),parameter :: covcn(nat) = &
      &[1.7036310879369_wp, 2.8486795615265_wp, 1.4951928728082_wp, &
      & 1.3047953045826_wp, .81178816819290_wp, .90351747314633_wp]
   real(wp),parameter :: q(nat) = &
      &[.39808668429315_wp, .41552173661168_wp,-.16133275860565_wp, &
      &-.19939812633388_wp,-.20923526498153_wp,-.24364227098378_wp]
   integer, parameter :: wsc_rep(3) = [1,1,1]
   real(wp),parameter :: g_a = 3.0_wp
   real(wp),parameter :: g_c = 2.0_wp
   real(wp),parameter :: wf  = 6.0_wp
   integer, parameter :: lmbd = p_mbd_approx_atm
   integer, parameter :: refqmode = p_refq_goedecker
   real(wp),parameter :: rthr_atm = 1600.0_wp
   real(wp),parameter :: rthr_vdw = 4000.0_wp
   integer, parameter :: vdw_rep(3) = [8,8,8]
   integer, parameter :: atm_rep(3)  = [5,5,5]
   type(dftd_parameter),parameter :: dparam_pbe    = dftd_parameter ( &
      &  s6=1.0000_wp, s8=0.95948085_wp, a1=0.38574991_wp, a2=4.80688534_wp )
   type(dftd_parameter),parameter :: dparam_tpss   = dftd_parameter ( &
      &  s6=1.0000_wp, s8=1.76596355_wp, a1=0.42822303_wp, a2=4.54257102_wp )
   type(dftd_parameter),parameter :: dparam_random = dftd_parameter ( &
      &  s6=0.95_wp, s8=0.45_wp, s10=0.65_wp, s9=1.10_wp, a1=0.43_wp, a2=5.10_wp )
   real(wp),parameter :: step = 1.0e-4_wp, step2 = 0.5_wp/step

   real(wp),allocatable :: xyz(:,:)
   integer              :: i,j
   type(TMolecule)       :: mol
   integer              :: ndim
   real(wp),allocatable :: gweights(:)   ! gaussian weights
   real(wp),allocatable :: refc6(:,:)    ! reference C6 coeffients
   real(wp)             :: energy,e2,e3

   allocate(xyz(3, nat))
   call coord_trafo(nat,lattice,abc,xyz)
   call init(mol, at, xyz, lattice=lattice)

   call d4init(g_a,g_c,refqmode)
   call d4dim(mol%n,mol%at,ndim)
   call assert_eq(ndim,26)

   allocate( gweights(ndim),refc6(ndim,ndim) )

   call pbc_d4(mol%n,ndim,mol%at,wf,g_a,g_c,covcn,gweights,refc6)
   call assert_close(gweights( 4),0.78386266917391E-03_wp,thr)
   call assert_close(gweights(13),0.74644931926423E+00_wp,thr)
   call assert_close(gweights(22),0.41525004702946E+00_wp,thr)

   call assert_close(refc6( 5,23),46.767562973793_wp,thr)
   call assert_close(refc6(12,12),20.848889774951_wp,thr)
   call assert_close(refc6(25, 9),68.827000333505_wp,thr)

   ! energies are intent(out)
   energy = 1.0_wp
   e2     = 1.0_wp
   e3     = 1.0_wp

   call edisp_3d(mol,ndim,q,vdw_rep,atm_rep,rthr_vdw,rthr_atm,dparam_pbe, &
      &          g_a,g_c,gweights,refc6,lmbd,energy,e2,e3)

   call assert_close(e2,    -0.37164345120511E-01_wp,thr)
   call assert_close(e3,     0.19951867090604E-02_wp,thr)
   call assert_close(energy,-0.35169158411451E-01_wp,thr)

   call edisp_3d(mol,ndim,q,vdw_rep,atm_rep,rthr_vdw,rthr_atm,dparam_tpss, &
      &          g_a,g_c,gweights,refc6,lmbd,energy,e2,e3)

   call assert_close(e2,    -0.48642371573216E-01_wp,thr)
   call assert_close(e3,     0.19786934360753E-02_wp,thr)
   call assert_close(energy,-0.46663678137140E-01_wp,thr)

   call edisp_3d(mol,ndim,q,vdw_rep,atm_rep,rthr_vdw,rthr_atm,dparam_random, &
      &          g_a,g_c,gweights,refc6,lmbd,energy,e2,e3)

   call assert_close(e2,    -0.27778143076318E-01_wp,thr)
   call assert_close(e3,     0.16050327438669E-02_wp,thr)
   call assert_close(energy,-0.26173110332451E-01_wp,thr)

   call mol%deallocate

   call terminate(afail)

end subroutine test_dftd4_pbc_energies

subroutine test_dftd4_cell_gradient
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use assertion
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_param
   use xtb_disp_dftd4
   use xtb_eeq
   use xtb_chargemodel
   use xtb_disp_ncoord
   use xtb_pbc_tools
   implicit none
   real(wp),parameter :: thr = 1.0e-10_wp
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
   real(wp),parameter :: g_a = 3.0_wp
   real(wp),parameter :: g_c = 2.0_wp
   real(wp),parameter :: wf  = 6.0_wp
   integer, parameter :: lmbd = p_mbd_approx_atm
   integer, parameter :: refqmode = p_refq_goedecker
   real(wp),parameter :: rthr_cn  = 1600.0_wp
   real(wp),parameter :: rthr_vdw = 4000.0_wp
   integer, parameter :: vdw_rep(3) = [8,8,8]
   integer            :: cn_rep(3)  = [5,5,5]
   type(TEnvironment) :: env
   type(dftd_parameter),parameter :: dparam_pbe    = dftd_parameter ( &
   &  s6=1.0000_wp, s8=0.95948085_wp, a1=0.38574991_wp, a2=4.80688534_wp )
   real(wp),parameter :: step = 1.0e-4_wp, step2 = 0.5_wp/step

   integer              :: i,j
   type(TMolecule)    :: mol
   integer              :: ndim
   real(wp) :: molpol,molc6,molc8        ! molecular Polarizibility
   real(wp),allocatable :: gweights(:)   ! gaussian weights
   real(wp),allocatable :: refc6(:,:)    ! reference C6 coeffients
   real(wp),allocatable :: c6ab(:,:)
   real(wp),allocatable :: aw(:,:)
   type(chrg_parameter) :: chrgeq
   real(wp)             :: energy
   real(wp)             :: stmp(3,3)
   real(wp)             :: sigma(3,3)
   real(wp)             :: er,el,ees
   real(wp),allocatable :: xyz(:,:)
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: dcndr(:,:,:)
   real(wp),allocatable :: dcndL(:,:,:)
   real(wp),allocatable :: covcn(:)
   real(wp),allocatable :: dcovcndr(:,:,:)
   real(wp),allocatable :: dcovcndL(:,:,:)
   real(wp),allocatable :: q(:)
   real(wp),allocatable :: qr(:)
   real(wp),allocatable :: ql(:)
   real(wp),allocatable :: dqdr(:,:,:)
   real(wp),allocatable :: dqdL(:,:,:)
   real(wp),allocatable :: gradient(:,:)
   real(wp),allocatable :: ges(:,:)
   real(wp),allocatable :: numg(:,:)
   real(wp),allocatable :: numq(:,:,:)


   allocate( cn(nat), dcndr(3,nat,nat), q(nat), dqdr(3,nat,nat+1), &
      &      dcndL(3,3,nat), dqdL(3,3,nat+1), dcovcndL(3,3,nat), &
      &      gradient(3,nat), ges(3,nat), numg(3,nat), ql(nat), qr(nat), &
      &      numq(3,nat,nat), covcn(nat), dcovcndr(3,nat,nat), source = 0.0_wp )

   allocate(xyz(3, nat))
   call coord_trafo(nat,lattice,abc,xyz)
   call init(mol, at, xyz, lattice=lattice)

   call d4init(g_a,g_c,refqmode)
   call d4dim(mol%n,mol%at,ndim)

   allocate( gweights(ndim),refc6(ndim,ndim) )

   call print_pbcsum(stdout,mol)

   call pbc_derfcoord(mol%n,mol%at,mol%xyz,mol%lattice,cn,dcndr,dcndL,900.0_wp)

   call new_charge_model_2019(chrgeq,mol%n,mol%at)

   call eeq_chrgeq(mol,env,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,energy,gradient,sigma,&
      &            .false.,.false.,.true.)

   call pbc_dncoord_d4(mol%n,mol%at,mol%xyz,mol%lattice, &
      &                covcn,dcovcndr,dcovcndL,rthr_cn)

   call pbc_d4(mol%n,ndim,mol%at,wf,g_a,g_c,cn,gweights,refc6)

   energy = 0.0_wp
   gradient = 0.0_wp
   sigma = 0.0_wp

   call dispgrad_3d(mol,ndim,q,covcn,dcovcndr,dcovcndL,vdw_rep,cn_rep, &
      &             rthr_vdw,rthr_cn,dparam_pbe,wf,g_a,g_c,refc6,lmbd, &
      &             gradient,sigma,energy,dqdr,dqdL)

   call assert_close(energy,-0.40714027910528E-01_wp,thr)
   call assert_close(gradient(1,1),-0.99274382274332E-05_wp,thr)
   call assert_close(gradient(2,3),-0.13713496197771E-03_wp,thr)
   call assert_close(gradient(3,6),-0.24645103343381E-03_wp,thr)

   call assert_close(sigma(1,1), 0.43552159014448E-01_wp,thr)
   call assert_close(sigma(2,3),-0.26923442064386E-03_wp,thr)
   call assert_close(sigma(3,1),-0.75395807484337E-03_wp,thr)
   call assert_close(sigma(2,1),sigma(1,2),thr)

   call mol%deallocate

   call terminate(afail)

end subroutine test_dftd4_cell_gradient

!> @brief test the general wrapper for DFT-D4 calculations
subroutine test_dftd4_api
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use assertion
   use xtb_type_molecule
   use xtb_type_options
   use xtb_type_param
   use xtb_disp_dftd4
   use xtb_calculators
   implicit none
   type(TMolecule)  :: mol

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.00000000000000_wp,  0.00000000000000_wp, -0.73578586109551_wp, &
      & 1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp, &
      &-1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp  &
      & ],shape(xyz))
   type(dftd_parameter),parameter :: dparam_b2plyp = dftd_parameter ( &
      &  s6=0.6400_wp, s8=1.16888646_wp, a1=0.44154604_wp, a2=4.73114642_wp )
   type(dftd_parameter),parameter :: dparam_tpss   = dftd_parameter ( &
      &  s6=1.0000_wp, s8=1.76596355_wp, a1=0.42822303_wp, a2=4.54257102_wp )
   type(dftd_options),  parameter :: opt_1 = dftd_options ( &
      &  lmbd = p_mbd_approx_atm, refq = p_refq_goedecker, &
      &  wf = 6.0_wp, g_a = 3.0_wp, g_c = 2.0_wp, &
      &  lmolpol=.false., lenergy=.true., lgradient=.false., &
      &  verbose = .false., veryverbose = .false., silent = .false. )
   type(dftd_options),  parameter :: opt_2 = dftd_options ( &
      &  lmbd = p_mbd_approx_atm, refq = p_refq_goedecker, &
      &  wf = 6.0_wp, g_a = 3.0_wp, g_c = 2.0_wp, &
      &  lmolpol=.false., lenergy=.false., lgradient=.true., &
      &  verbose = .false., veryverbose = .false., silent = .true. )

   real(wp) :: energy
   real(wp),allocatable :: gradient(:,:)

   allocate( gradient(3,nat) )
   energy   = 0.0_wp
   gradient = 0.0_wp

   call init(mol, at, xyz)
   
   call d4_calculation(stdout,opt_1,mol,dparam_tpss,energy,gradient)
   call assert_close(energy,-0.26682682254336E-03_wp,thr)

   call d4_calculation(stdout,opt_2,mol,dparam_b2plyp,energy,gradient)
   call assert_close(energy,-0.13368190339570E-03_wp,thr)

   call assert_close(gradient(1,1), 0.00000000000000E+00_wp,thr)
   call assert_close(gradient(3,1), 0.39778648945254E-04_wp,thr)
   call assert_close(gradient(3,2),-0.19889324472627E-04_wp,thr)
   call assert_close(gradient(1,2),-gradient(1,3),          thr)

   call mol%deallocate

   ! done: everythings fine
   call terminate(afail)
end subroutine test_dftd4_api

!> @brief test the general wrapper for DFT-D4 calculations
subroutine test_dftd4_pbc_api
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use assertion
   use xtb_type_molecule
   use xtb_type_options
   use xtb_type_param
   use xtb_disp_dftd4
   use xtb_pbc_tools
   use xtb_calculators
   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
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
   type(dftd_parameter),parameter :: dparam_pbe    = dftd_parameter ( &
      &  s6=1.0000_wp, s8=0.95948085_wp, a1=0.38574991_wp, a2=4.80688534_wp )
   type(dftd_options),  parameter :: opt = dftd_options ( &
      &  lmbd = p_mbd_approx_atm, refq = p_refq_goedecker, &
      &  wf = 6.0_wp, g_a = 3.0_wp, g_c = 2.0_wp, &
      &  lmolpol=.false., lenergy=.false., lgradient=.true., &
      &  verbose = .false., veryverbose = .false., silent = .true. )

   integer              :: i,j
   type(TMolecule)    :: mol
   integer              :: ndim
   real(wp) :: energy
   real(wp) :: lattice_grad(3,3)
   real(wp),allocatable :: xyz(:,:)
   real(wp),allocatable :: gradient(:,:)

   allocate(xyz(3, nat))
   call coord_trafo(nat,lattice,abc,xyz)
   call init(mol, at, xyz, lattice=lattice)

   allocate( gradient(3,nat) )
   energy   = 0.0_wp
   gradient = 0.0_wp
   lattice_grad = 0.0_wp

   call d4_pbc_calculation(stdout,opt,mol,dparam_pbe,energy,gradient,lattice_grad)
   call assert_close(energy,-0.40714027910526E-01_wp,thr)

   call assert_close(norm2(gradient),0.16753166250887E-02_wp,thr)
   call assert_close(gradient(2,5),  0.39323973508111E-04_wp,thr)
   call assert_close(gradient(1,1), -0.99274382274340E-05_wp,thr)

   call assert_close(norm2(lattice_grad),0.85938071998117E-02_wp,thr)
   call assert_close(lattice_grad(2,1),  0.13258077458069E-03_wp,thr)
   call assert_close(lattice_grad(3,3),  0.49734876304364E-02_wp,thr)

   call mol%deallocate

   ! done: everythings fine
   call terminate(afail)
end subroutine test_dftd4_pbc_api


subroutine test_dftd4_pbc3d_neighbourlist
   use xtb_mctc_accuracy, only : wp
   use assertion

   use xtb_mctc_convert, only : aatoau

   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init, len
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_type_latticepoint, only : TLatticePoint, init
   use xtb_type_param, only : dftd_parameter

   use xtb_disp_coordinationnumber, only : getCoordinationNumber, cnType, &
      & cutCoordinationNumber
   use xtb_disp_dftd4, only : d4_gradient, p_refq_goedecker, d4init
   use xtb_disp_encharges, only : getENCharges

   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 32
   integer, parameter :: at(nat) = [8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[4.5853168464880421_wp,  4.9392326929575878_wp,  4.1894081210748118_wp,  &
      & 5.8862267152491423_wp,  1.1425258978245871_wp,  6.5015058204768126_wp,  &
      & 1.4284279616220412_wp,  4.6017511285540875_wp,  3.1465884436348119_wp,  &
      & 2.3323404704521411_wp,  1.4471801154820869_wp,  0.7121932185858125_wp,  &
      & 4.7333543155561415_wp,  2.7747291872305868_wp,  3.1951352976178122_wp,  &
      & 5.6617754419101418_wp,  2.2133191164485870_wp,  2.1235404838618126_wp,  &
      & 5.3107598618381422_wp,  2.6902988056185868_wp,  0.7384466319968125_wp,  &
      & 4.4947071071761426_wp,  3.9530790692635867_wp,  0.6801776747598124_wp,  &
      & 4.8005171923760424_wp,  4.9185102874975870_wp,  1.8186363449528122_wp,  &
      & 4.6951362687070421_wp,  4.2781752812835867_wp,  3.1816411821728123_wp,  &
      & 1.3838574419160412_wp,  5.2817805008910863_wp,  4.1482702947948136_wp,  &
      & 1.0268974195990415_wp,  4.7234752637800881_wp,  5.4989995400388123_wp,  &
      & 2.0852659694760409_wp,  5.0956317453800875_wp,  6.5351699846458127_wp,  &
      & 2.3344644666691412_wp, -0.0736561690909131_wp,  6.4245628001158135_wp,  &
      & 2.4894017448231409_wp,  0.6213510313930869_wp,  5.0967297417158131_wp,  &
      & 1.5745272273791413_wp,  0.1243470825760870_wp,  3.9731040773988129_wp,  &
      & 5.8221065925130420_wp,  5.3013563342055878_wp,  1.7264876737078123_wp,  &
      & 3.4487807319551416_wp,  3.6355832152975864_wp,  0.7429568016758125_wp,  &
      & 4.8499393376520423_wp,  3.4713855169305874_wp,  6.4691872586348129_wp,  &
      & 0.2495364434351412_wp,  2.4795455690160870_wp,  2.1043557230378123_wp,  &
      & 5.6691068338331423_wp,  1.1234174220755870_wp,  2.1414388326468128_wp,  &
      & 3.7072009289431418_wp,  2.4357632918535872_wp,  3.0094700999208119_wp,  &
      & 4.1414520030430415_wp,  5.7877262477775879_wp,  1.7803680119358125_wp,  &
      & 5.0142851411171421_wp,  2.4165926460955873_wp,  4.1857610486448129_wp,  &
      & 3.0280930003030413_wp,  4.6201081184690871_wp,  6.2533190952188136_wp,  &
      & 0.5863628696651412_wp,  0.5757236365910867_wp,  4.1021714214668128_wp,  &
      & 2.3776130524831411_wp,  1.6969724987740866_wp,  5.2327688986668139_wp,  &
      & 1.9486148363011413_wp,  0.4390675147070869_wp,  2.9999022491838123_wp,  &
      & 3.5312997625581413_wp,  0.4467415528495868_wp,  4.8114121395028135_wp,  &
      & 6.5089895990100421_wp,  5.2100409408535882_wp,  6.0066553789008132_wp,  &
      & 0.9001165013630412_wp,  3.6420787128610868_wp,  5.4413106648508132_wp,  &
      & 1.6012116650460413_wp,  5.6845471271780879_wp,  0.7675566847298124_wp], &
      & shape(xyz)) * aatoau
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[6.4411018522600_wp,    0.0492571261505_wp,    0.2192046129910_wp,  &
      & 0.0462076831749_wp,    6.6435057067500_wp,    0.1670513770770_wp,  &
      & 0.2262248220170_wp,   -0.9573234940220_wp,    6.7608039126200_wp], &
      & shape(lattice)) * aatoau
   type(dftd_parameter),parameter :: dparam_pbe = dftd_parameter ( &
      & s6=1.0_wp, s8=0.95948085_wp, a1=0.38574991_wp, a2=4.80688534_wp, s9=0.0_wp)
   real(wp),parameter :: g_a = 3.0_wp
   real(wp),parameter :: g_c = 2.0_wp
   real(wp),parameter :: wf  = 6.0_wp
   integer, parameter :: refqmode = p_refq_goedecker

   type(TMolecule) :: mol
   type(TEnvironment) :: env
   type(TLatticePoint) :: latp
   type(TNeighbourlist) :: neighList

   integer, allocatable :: neighs(:)
   real(wp), allocatable :: trans(:, :)

   integer :: ii, jj
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp) :: energy, sigma(3, 3), sdum(3, 3)
   real(wp), allocatable :: gradient(:, :), gdum(:, :)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: q(:), dqdr(:,:,:), dqdL(:,:,:)

   call init(env)

   call init(mol, at, xyz, lattice=lattice)

   call d4init(g_a,g_c,refqmode)
   allocate(neighs(nat))
   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat))
   allocate(q(nat), dqdr(3, nat, nat), dqdL(3, 3, nat))
   allocate(gradient(3, nat), gdum(3, nat))

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call init(latp, env, mol, 60.0_wp)  ! Fixed cutoff
   call latp%getLatticePoints(trans, 60.0_wp)
   call init(neighList, len(mol))
   call neighList%generate(env, mol%xyz, 60.0_wp, trans, .false.)

!   cn(:) = 0.0_wp
!   dcndr(:, :, :) = 0.0_wp
!   dcndL(:, :, :) = 0.0_wp

   q(:) = 0.0_wp
   dqdr(:, :, :) = 0.0_wp
   dqdL(:, :, :) = 0.0_wp

   call neighlist%getNeighs(neighs, 40.0_wp)
   call getCoordinationNumber(mol, neighs, neighList, cnType%erf, &
      & cn, dcndr, dcndL)
   call cutCoordinationNumber(mol%n, cn, dcndr, dcndL, maxCN=8.0_wp)
   call getENCharges(env, mol, cn, dcndr, dcndL, q, dqdr, dqdL)
   call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
      & cn, dcndr, dcndL)
   call neighlist%getNeighs(neighs, 60.0_wp)
   call d4_gradient(mol, neighs, neighList, dparam_pbe, g_a, g_c, wf, &
      & cn, dcndr, dcndL, q, dqdr, dqdL, energy, gradient, sigma)

   call assert_close(energy, -0.73670210332109E-01_wp, thr)
   call assert_close(norm2(gradient), 0.67876404541808E-03_wp, thr)

   call assert_close(gradient(1, 3), -0.48796348327688E-04_wp, thr)
   call assert_close(gradient(2, 7),  0.85632660892000E-04_wp, thr)
   call assert_close(gradient(3,12),  0.15057422955536E-04_wp, thr)

   call assert_close(sigma(1,1),  0.80650010267247E-01_wp, thr)
   call assert_close(sigma(2,1), -0.20505317558255E-03_wp, thr)
   call assert_close(sigma(3,2),  0.21880502391427E-02_wp, thr)

   ! check numerical gradient
   ! reduce the number of numerical gradient evaluations significantly
   do ii = 1, 2
      do jj = ii, 3, 2
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call neighlist%getNeighs(neighs, 40.0_wp)
         call getCoordinationNumber(mol, neighs, neighList, cnType%erf, &
            & cn, dcndr, dcndL)
         call cutCoordinationNumber(mol%n, cn, dcndr, dcndL, maxCN=8.0_wp)
         call getENCharges(env, mol, cn, dcndr, dcndL, q, dqdr, dqdL)
         call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
            & cn, dcndr, dcndL)
         call neighlist%getNeighs(neighs, 60.0_wp)
         call d4_gradient(mol, neighs, neighList, dparam_pbe, g_a, g_c, wf, &
            & cn, dcndr, dcndL, q, dqdr, dqdL, er, gdum, sdum)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step

         call mol%update
         call neighList%update(mol%xyz, trans)
         call neighlist%getNeighs(neighs, 40.0_wp)
         call getCoordinationNumber(mol, neighs, neighList, cnType%erf, &
            & cn, dcndr, dcndL)
         call cutCoordinationNumber(mol%n, cn, dcndr, dcndL, maxCN=8.0_wp)
         call getENCharges(env, mol, cn, dcndr, dcndL, q, dqdr, dqdL)
         call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
            & cn, dcndr, dcndL)
         call neighlist%getNeighs(neighs, 60.0_wp)
         call d4_gradient(mol, neighs, neighList, dparam_pbe, g_a, g_c, wf, &
            & cn, dcndr, dcndL, q, dqdr, dqdL, el, gdum, sdum)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call assert_close(gradient(jj, ii), (er - el)*step2, thr)
      end do
   end do

!   ! check numerical strain derivatives
!   trans2 = trans
!   eps = unity
!   do ii = 1, 3
!      do jj = 1, ii
!         eps(jj, ii) = eps(jj, ii) + step
!         mol%lattice(:, :) = matmul(eps, lattice)
!         mol%xyz(:, :) = matmul(eps, xyz)
!         trans(:, :) = matmul(eps, trans2)
!         call mol%update
!         !call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
!            !& cn, dcndr, dcndL)
!         call d4_atm_gradient(mol, neighs, neighlist, dparam_pbe, g_a, g_c, wf, &
!            & cn, dcndr, dcndL, er, gdum, sdum)
!
!         eps(jj, ii) = eps(jj, ii) - 2*step
!         mol%lattice(:, :) = matmul(eps, lattice)
!         mol%xyz(:, :) = matmul(eps, xyz)
!         trans(:, :) = matmul(eps, trans2)
!         call mol%update
!         call neighList%update(mol%xyz, trans)
!         !call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
!            !& cn, dcndr, dcndL)
!         call d4_atm_gradient(mol, neighs, neighlist, dparam_pbe, g_a, g_c, wf, &
!            & cn, dcndr, dcndL, el, gdum, sdum)
!
!         eps(jj, ii) = eps(jj, ii) + step
!
!         call assert_close(sigma(jj, ii), (er - el)*step2, thr)
!      end do
!   end do

   call terminate(afail)

end subroutine test_dftd4_pbc3d_neighbourlist


subroutine test_dftd4_pbc3d_latticepoints
   use xtb_mctc_accuracy, only : wp
   use assertion

   use xtb_mctc_convert, only : aatoau

   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init, len
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_type_latticepoint, only : TLatticePoint, init
   use xtb_type_wignerseitzcell, only : TWignerSeitzCell, init
   use xtb_type_param, only : dftd_parameter

   use xtb_disp_coordinationnumber, only : getCoordinationNumber, cnType, &
      & cutCoordinationNumber
   use xtb_disp_dftd4, only : d4_gradient, p_refq_goedecker, d4init
   use xtb_disp_encharges, only : getENCharges

   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
   real(wp),parameter :: thr2 = 1.0e-9_wp
   integer, parameter :: nat = 32
   integer, parameter :: at(nat) = [8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[4.5853168464880421_wp,  4.9392326929575878_wp,  4.1894081210748118_wp,  &
      & 5.8862267152491423_wp,  1.1425258978245871_wp,  6.5015058204768126_wp,  &
      & 1.4284279616220412_wp,  4.6017511285540875_wp,  3.1465884436348119_wp,  &
      & 2.3323404704521411_wp,  1.4471801154820869_wp,  0.7121932185858125_wp,  &
      & 4.7333543155561415_wp,  2.7747291872305868_wp,  3.1951352976178122_wp,  &
      & 5.6617754419101418_wp,  2.2133191164485870_wp,  2.1235404838618126_wp,  &
      & 5.3107598618381422_wp,  2.6902988056185868_wp,  0.7384466319968125_wp,  &
      & 4.4947071071761426_wp,  3.9530790692635867_wp,  0.6801776747598124_wp,  &
      & 4.8005171923760424_wp,  4.9185102874975870_wp,  1.8186363449528122_wp,  &
      & 4.6951362687070421_wp,  4.2781752812835867_wp,  3.1816411821728123_wp,  &
      & 1.3838574419160412_wp,  5.2817805008910863_wp,  4.1482702947948136_wp,  &
      & 1.0268974195990415_wp,  4.7234752637800881_wp,  5.4989995400388123_wp,  &
      & 2.0852659694760409_wp,  5.0956317453800875_wp,  6.5351699846458127_wp,  &
      & 2.3344644666691412_wp, -0.0736561690909131_wp,  6.4245628001158135_wp,  &
      & 2.4894017448231409_wp,  0.6213510313930869_wp,  5.0967297417158131_wp,  &
      & 1.5745272273791413_wp,  0.1243470825760870_wp,  3.9731040773988129_wp,  &
      & 5.8221065925130420_wp,  5.3013563342055878_wp,  1.7264876737078123_wp,  &
      & 3.4487807319551416_wp,  3.6355832152975864_wp,  0.7429568016758125_wp,  &
      & 4.8499393376520423_wp,  3.4713855169305874_wp,  6.4691872586348129_wp,  &
      & 0.2495364434351412_wp,  2.4795455690160870_wp,  2.1043557230378123_wp,  &
      & 5.6691068338331423_wp,  1.1234174220755870_wp,  2.1414388326468128_wp,  &
      & 3.7072009289431418_wp,  2.4357632918535872_wp,  3.0094700999208119_wp,  &
      & 4.1414520030430415_wp,  5.7877262477775879_wp,  1.7803680119358125_wp,  &
      & 5.0142851411171421_wp,  2.4165926460955873_wp,  4.1857610486448129_wp,  &
      & 3.0280930003030413_wp,  4.6201081184690871_wp,  6.2533190952188136_wp,  &
      & 0.5863628696651412_wp,  0.5757236365910867_wp,  4.1021714214668128_wp,  &
      & 2.3776130524831411_wp,  1.6969724987740866_wp,  5.2327688986668139_wp,  &
      & 1.9486148363011413_wp,  0.4390675147070869_wp,  2.9999022491838123_wp,  &
      & 3.5312997625581413_wp,  0.4467415528495868_wp,  4.8114121395028135_wp,  &
      & 6.5089895990100421_wp,  5.2100409408535882_wp,  6.0066553789008132_wp,  &
      & 0.9001165013630412_wp,  3.6420787128610868_wp,  5.4413106648508132_wp,  &
      & 1.6012116650460413_wp,  5.6845471271780879_wp,  0.7675566847298124_wp], &
      & shape(xyz)) * aatoau
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[6.4411018522600_wp,    0.0492571261505_wp,    0.2192046129910_wp,  &
      & 0.0462076831749_wp,    6.6435057067500_wp,    0.1670513770770_wp,  &
      & 0.2262248220170_wp,   -0.9573234940220_wp,    6.7608039126200_wp], &
      & shape(lattice)) * aatoau
   type(dftd_parameter),parameter :: dparam_pbe = dftd_parameter ( &
      & s6=1.0_wp, s8=0.95948085_wp, a1=0.38574991_wp, a2=4.80688534_wp, s9=0.0_wp)
   real(wp),parameter :: g_a = 3.0_wp
   real(wp),parameter :: g_c = 2.0_wp
   real(wp),parameter :: wf  = 6.0_wp
   integer, parameter :: refqmode = p_refq_goedecker

   type(TMolecule) :: mol
   type(TEnvironment) :: env
   type(TLatticePoint) :: latp
   type(TNeighbourlist) :: neighList

   integer, allocatable :: neighs(:)
   real(wp), allocatable :: trans(:, :), trans2(:, :)

   integer :: ii, jj
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp) :: energy, sigma(3, 3), sdum(3, 3)
   real(wp), allocatable :: gradient(:, :), gdum(:, :)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: q(:), dqdr(:,:,:), dqdL(:,:,:)

   call init(env)

   call init(mol, at, xyz, lattice=lattice)

   call d4init(g_a,g_c,refqmode)
   allocate(neighs(nat))
   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat))
   allocate(q(nat), dqdr(3, nat, nat), dqdL(3, 3, nat))
   allocate(gradient(3, nat), gdum(3, nat))

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

!   cn(:) = 0.0_wp
!   dcndr(:, :, :) = 0.0_wp
!   dcndL(:, :, :) = 0.0_wp

   q(:) = 0.0_wp
   dqdr(:, :, :) = 0.0_wp
   dqdL(:, :, :) = 0.0_wp

   call init(latp, env, mol, 60.0_wp)  ! Fixed cutoff
   call latp%getLatticePoints(trans, 40.0_wp)
   call init(neighList, len(mol))
   call neighList%generate(env, mol%xyz, 40.0_wp, trans, .false.)

   call neighlist%getNeighs(neighs, 40.0_wp)
   call getCoordinationNumber(mol, neighs, neighList, cnType%erf, &
      & cn, dcndr, dcndL)
   call cutCoordinationNumber(mol%n, cn, dcndr, dcndL, maxCN=8.0_wp)
   call getENCharges(env, mol, cn, dcndr, dcndL, q, dqdr, dqdL)
   call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
      & cn, dcndr, dcndL)
   call latp%getLatticePoints(trans, 60.0_wp)
   call d4_gradient(mol, trans, dparam_pbe, g_a, g_c, wf, 60.0_wp, &
      & cn, dcndr, dcndL, q, dqdr, dqdL, energy, gradient, sigma)

   call assert_close(energy, -0.73670210332109E-01_wp, thr)
   call assert_close(norm2(gradient), 0.67876404541808E-03_wp, thr)

   call assert_close(gradient(1, 3), -0.48796348327688E-04_wp, thr)
   call assert_close(gradient(2, 7),  0.85632660892001E-04_wp, thr)
   call assert_close(gradient(3,12),  0.15057422955536E-04_wp, thr)

   call assert_close(sigma(1,1),  0.80650010267246E-01_wp, thr)
   call assert_close(sigma(2,1), -0.20505317558255E-03_wp, thr)
   call assert_close(sigma(3,2),  0.21880502391427E-02_wp, thr)

   ! check numerical gradient
   do ii = 1, nat, 9
      do jj = mod(ii, 2) + 1, 3, 2
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call neighlist%getNeighs(neighs, 40.0_wp)
         call getCoordinationNumber(mol, neighs, neighList, cnType%erf, &
            & cn, dcndr, dcndL)
         call cutCoordinationNumber(mol%n, cn, dcndr, dcndL, maxCN=8.0_wp)
         call getENCharges(env, mol, cn, dcndr, dcndL, q, dqdr, dqdL)
         call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
            & cn, dcndr, dcndL)
         call d4_gradient(mol, trans, dparam_pbe, g_a, g_c, wf, 60.0_wp, &
            & cn, dcndr, dcndL, q, dqdr, dqdL, er, gdum, sdum)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step

         call mol%update
         call neighList%update(mol%xyz, trans)
         call neighlist%getNeighs(neighs, 40.0_wp)
         call getCoordinationNumber(mol, neighs, neighList, cnType%erf, &
            & cn, dcndr, dcndL)
         call cutCoordinationNumber(mol%n, cn, dcndr, dcndL, maxCN=8.0_wp)
         call getENCharges(env, mol, cn, dcndr, dcndL, q, dqdr, dqdL)
         call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
            & cn, dcndr, dcndL)
         call d4_gradient(mol, trans, dparam_pbe, g_a, g_c, wf, 60.0_wp, &
            & cn, dcndr, dcndL, q, dqdr, dqdL, el, gdum, sdum)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call assert_close(gradient(jj, ii), (er - el)*step2, thr2)
      end do
   end do

!   ! check numerical strain derivatives
!   trans2 = trans
!   eps = unity
!   do ii = 1, 3
!      do jj = 1, ii
!         er = 0.0_wp
!         el = 0.0_wp
!         eps(jj, ii) = eps(jj, ii) + step
!         mol%lattice(:, :) = matmul(eps, lattice)
!         mol%xyz(:, :) = matmul(eps, xyz)
!         trans(:, :) = matmul(eps, trans2)
!         call mol%update
!         call neighList%update(mol%xyz, trans)
!         call neighlist%getNeighs(neighs, 40.0_wp)
!         call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
!            & cn, dcndr, dcndL)
!         call d4_gradient(mol, trans, dparam_pbe, g_a, g_c, wf, 60.0_wp, &
!            & cn, dcndr, dcndL, er, gdum, sdum)
!
!         eps(jj, ii) = eps(jj, ii) - 2*step
!         mol%lattice(:, :) = matmul(eps, lattice)
!         mol%xyz(:, :) = matmul(eps, xyz)
!         trans(:, :) = matmul(eps, trans2)
!         call mol%update
!         call neighList%update(mol%xyz, trans)
!         call neighlist%getNeighs(neighs, 40.0_wp)
!         call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
!            & cn, dcndr, dcndL)
!         call d4_gradient(mol, trans, dparam_pbe, g_a, g_c, wf, 60.0_wp, &
!            & cn, dcndr, dcndL, el, gdum, sdum)
!
!         eps(jj, ii) = eps(jj, ii) + step
!
!         call assert_close(sigma(jj, ii), (er - el)*step2, thr2)
!      end do
!   end do

   call terminate(afail)

end subroutine test_dftd4_pbc3d_latticepoints


subroutine test_dftd4_pbc3d_threebody_neighs
   use xtb_mctc_accuracy, only : wp
   use assertion

   use xtb_mctc_convert, only : aatoau

   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init, len
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_type_latticepoint, only : TLatticePoint, init
   use xtb_type_wignerseitzcell, only : TWignerSeitzCell, init
   use xtb_type_param, only : dftd_parameter

   use xtb_disp_coordinationnumber, only : getCoordinationNumber, cnType
   use xtb_disp_dftd4, only : d4_atm_gradient, p_refq_goedecker, d4init
   use xtb_disp_encharges, only : getENCharges

   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 32
   integer, parameter :: at(nat) = [8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[4.5853168464880421_wp,  4.9392326929575878_wp,  4.1894081210748118_wp,  &
      & 5.8862267152491423_wp,  1.1425258978245871_wp,  6.5015058204768126_wp,  &
      & 1.4284279616220412_wp,  4.6017511285540875_wp,  3.1465884436348119_wp,  &
      & 2.3323404704521411_wp,  1.4471801154820869_wp,  0.7121932185858125_wp,  &
      & 4.7333543155561415_wp,  2.7747291872305868_wp,  3.1951352976178122_wp,  &
      & 5.6617754419101418_wp,  2.2133191164485870_wp,  2.1235404838618126_wp,  &
      & 5.3107598618381422_wp,  2.6902988056185868_wp,  0.7384466319968125_wp,  &
      & 4.4947071071761426_wp,  3.9530790692635867_wp,  0.6801776747598124_wp,  &
      & 4.8005171923760424_wp,  4.9185102874975870_wp,  1.8186363449528122_wp,  &
      & 4.6951362687070421_wp,  4.2781752812835867_wp,  3.1816411821728123_wp,  &
      & 1.3838574419160412_wp,  5.2817805008910863_wp,  4.1482702947948136_wp,  &
      & 1.0268974195990415_wp,  4.7234752637800881_wp,  5.4989995400388123_wp,  &
      & 2.0852659694760409_wp,  5.0956317453800875_wp,  6.5351699846458127_wp,  &
      & 2.3344644666691412_wp, -0.0736561690909131_wp,  6.4245628001158135_wp,  &
      & 2.4894017448231409_wp,  0.6213510313930869_wp,  5.0967297417158131_wp,  &
      & 1.5745272273791413_wp,  0.1243470825760870_wp,  3.9731040773988129_wp,  &
      & 5.8221065925130420_wp,  5.3013563342055878_wp,  1.7264876737078123_wp,  &
      & 3.4487807319551416_wp,  3.6355832152975864_wp,  0.7429568016758125_wp,  &
      & 4.8499393376520423_wp,  3.4713855169305874_wp,  6.4691872586348129_wp,  &
      & 0.2495364434351412_wp,  2.4795455690160870_wp,  2.1043557230378123_wp,  &
      & 5.6691068338331423_wp,  1.1234174220755870_wp,  2.1414388326468128_wp,  &
      & 3.7072009289431418_wp,  2.4357632918535872_wp,  3.0094700999208119_wp,  &
      & 4.1414520030430415_wp,  5.7877262477775879_wp,  1.7803680119358125_wp,  &
      & 5.0142851411171421_wp,  2.4165926460955873_wp,  4.1857610486448129_wp,  &
      & 3.0280930003030413_wp,  4.6201081184690871_wp,  6.2533190952188136_wp,  &
      & 0.5863628696651412_wp,  0.5757236365910867_wp,  4.1021714214668128_wp,  &
      & 2.3776130524831411_wp,  1.6969724987740866_wp,  5.2327688986668139_wp,  &
      & 1.9486148363011413_wp,  0.4390675147070869_wp,  2.9999022491838123_wp,  &
      & 3.5312997625581413_wp,  0.4467415528495868_wp,  4.8114121395028135_wp,  &
      & 6.5089895990100421_wp,  5.2100409408535882_wp,  6.0066553789008132_wp,  &
      & 0.9001165013630412_wp,  3.6420787128610868_wp,  5.4413106648508132_wp,  &
      & 1.6012116650460413_wp,  5.6845471271780879_wp,  0.7675566847298124_wp], &
      & shape(xyz)) * aatoau
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[6.4411018522600_wp,    0.0492571261505_wp,    0.2192046129910_wp,  &
      & 0.0462076831749_wp,    6.6435057067500_wp,    0.1670513770770_wp,  &
      & 0.2262248220170_wp,   -0.9573234940220_wp,    6.7608039126200_wp], &
      & shape(lattice)) * aatoau
   type(dftd_parameter),parameter :: dparam_pbe = dftd_parameter ( &
      & s6=1.0_wp, s8=0.95948085_wp, a1=0.38574991_wp, a2=4.80688534_wp, s9=1.0_wp)
   real(wp),parameter :: g_a = 3.0_wp
   real(wp),parameter :: g_c = 2.0_wp
   real(wp),parameter :: wf  = 6.0_wp
   integer, parameter :: refqmode = p_refq_goedecker

   type(TMolecule) :: mol
   type(TEnvironment) :: env
   type(TLatticePoint) :: latp
   type(TNeighbourlist) :: neighList

   integer, allocatable :: neighs(:)
   real(wp), allocatable :: trans(:, :), trans2(:, :)

   integer :: ii, jj
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp) :: energy, sigma(3, 3), sdum(3, 3)
   real(wp), allocatable :: gradient(:, :), gdum(:, :)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: q(:), dqdr(:,:,:), dqdL(:,:,:)

   call init(env)

   call init(mol, at, xyz)

   call d4init(g_a,g_c,refqmode)
   allocate(neighs(nat))
   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat))
   allocate(gradient(3, nat), gdum(3, nat))

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call init(latp, env, mol, 15.0_wp)  ! Fixed cutoff
   call latp%getLatticePoints(trans, 15.0_wp)
   call init(neighList, len(mol))
   call neighList%generate(env, mol%xyz, 15.0_wp, trans, .false.)

   call neighlist%getNeighs(neighs, 15.0_wp)
   call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
      & cn, dcndr, dcndL)
   call d4_atm_gradient(mol, neighs, neighlist, dparam_pbe, g_a, g_c, wf, &
      & cn, dcndr, dcndL, energy, gradient, sigma)

   call assert_close(energy, 0.44030151321408E-03_wp, thr)
   call assert_close(norm2(gradient), 0.54938536852646E-04_wp, thr)

   call assert_close(gradient(1, 3), -0.44803663456486E-05_wp, thr)
   call assert_close(gradient(2, 7), -0.11366396738396E-05_wp, thr)
   call assert_close(gradient(3,12), -0.76710649053134E-05_wp, thr)

   call assert_close(sigma(1,1), -0.18373605865520E-03_wp, thr)
   call assert_close(sigma(2,1), -0.33090837432067E-04_wp, thr)
   call assert_close(sigma(3,2),  0.34358178263607E-04_wp, thr)

   ! check numerical gradient
   do ii = 1, nat, 5
      do jj = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
            & cn, dcndr, dcndL)
         call d4_atm_gradient(mol, neighs, neighlist, dparam_pbe, g_a, g_c, wf, &
            & cn, dcndr, dcndL, er, gdum, sdum)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
            & cn, dcndr, dcndL)
         call d4_atm_gradient(mol, neighs, neighlist, dparam_pbe, g_a, g_c, wf, &
            & cn, dcndr, dcndL, el, gdum, sdum)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call assert_close(gradient(jj, ii), (er - el)*step2, thr)
      end do
   end do

!   ! check numerical strain derivatives
!   trans2 = trans
!   eps = unity
!   do ii = 1, 3
!      do jj = 1, ii
!         eps(jj, ii) = eps(jj, ii) + step
!         mol%lattice(:, :) = matmul(eps, lattice)
!         mol%xyz(:, :) = matmul(eps, xyz)
!         trans(:, :) = matmul(eps, trans2)
!         call mol%update
!         !call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
!            !& cn, dcndr, dcndL)
!         call d4_atm_gradient(mol, neighs, neighlist, dparam_pbe, g_a, g_c, wf, &
!            & cn, dcndr, dcndL, er, gdum, sdum)
!
!         eps(jj, ii) = eps(jj, ii) - 2*step
!         mol%lattice(:, :) = matmul(eps, lattice)
!         mol%xyz(:, :) = matmul(eps, xyz)
!         trans(:, :) = matmul(eps, trans2)
!         call mol%update
!         call neighList%update(mol%xyz, trans)
!         !call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
!            !& cn, dcndr, dcndL)
!         call d4_atm_gradient(mol, neighs, neighlist, dparam_pbe, g_a, g_c, wf, &
!            & cn, dcndr, dcndL, el, gdum, sdum)
!
!         eps(jj, ii) = eps(jj, ii) + step
!
!         call assert_close(sigma(jj, ii), (er - el)*step2, thr)
!      end do
!   end do

   call terminate(afail)
end subroutine test_dftd4_pbc3d_threebody_neighs


subroutine test_dftd4_pbc3d_threebody_latp
   use xtb_mctc_accuracy, only : wp
   use assertion

   use xtb_mctc_convert, only : aatoau

   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init, len
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_type_latticepoint, only : TLatticePoint, init
   use xtb_type_wignerseitzcell, only : TWignerSeitzCell, init
   use xtb_type_param, only : dftd_parameter

   use xtb_disp_coordinationnumber, only : getCoordinationNumber, cnType
   use xtb_disp_dftd4, only : d4_atm_gradient, p_refq_goedecker, d4init
   use xtb_disp_encharges, only : getENCharges

   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 32
   integer, parameter :: at(nat) = [8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
      & 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[4.5853168464880421_wp,  4.9392326929575878_wp,  4.1894081210748118_wp,  &
      & 5.8862267152491423_wp,  1.1425258978245871_wp,  6.5015058204768126_wp,  &
      & 1.4284279616220412_wp,  4.6017511285540875_wp,  3.1465884436348119_wp,  &
      & 2.3323404704521411_wp,  1.4471801154820869_wp,  0.7121932185858125_wp,  &
      & 4.7333543155561415_wp,  2.7747291872305868_wp,  3.1951352976178122_wp,  &
      & 5.6617754419101418_wp,  2.2133191164485870_wp,  2.1235404838618126_wp,  &
      & 5.3107598618381422_wp,  2.6902988056185868_wp,  0.7384466319968125_wp,  &
      & 4.4947071071761426_wp,  3.9530790692635867_wp,  0.6801776747598124_wp,  &
      & 4.8005171923760424_wp,  4.9185102874975870_wp,  1.8186363449528122_wp,  &
      & 4.6951362687070421_wp,  4.2781752812835867_wp,  3.1816411821728123_wp,  &
      & 1.3838574419160412_wp,  5.2817805008910863_wp,  4.1482702947948136_wp,  &
      & 1.0268974195990415_wp,  4.7234752637800881_wp,  5.4989995400388123_wp,  &
      & 2.0852659694760409_wp,  5.0956317453800875_wp,  6.5351699846458127_wp,  &
      & 2.3344644666691412_wp, -0.0736561690909131_wp,  6.4245628001158135_wp,  &
      & 2.4894017448231409_wp,  0.6213510313930869_wp,  5.0967297417158131_wp,  &
      & 1.5745272273791413_wp,  0.1243470825760870_wp,  3.9731040773988129_wp,  &
      & 5.8221065925130420_wp,  5.3013563342055878_wp,  1.7264876737078123_wp,  &
      & 3.4487807319551416_wp,  3.6355832152975864_wp,  0.7429568016758125_wp,  &
      & 4.8499393376520423_wp,  3.4713855169305874_wp,  6.4691872586348129_wp,  &
      & 0.2495364434351412_wp,  2.4795455690160870_wp,  2.1043557230378123_wp,  &
      & 5.6691068338331423_wp,  1.1234174220755870_wp,  2.1414388326468128_wp,  &
      & 3.7072009289431418_wp,  2.4357632918535872_wp,  3.0094700999208119_wp,  &
      & 4.1414520030430415_wp,  5.7877262477775879_wp,  1.7803680119358125_wp,  &
      & 5.0142851411171421_wp,  2.4165926460955873_wp,  4.1857610486448129_wp,  &
      & 3.0280930003030413_wp,  4.6201081184690871_wp,  6.2533190952188136_wp,  &
      & 0.5863628696651412_wp,  0.5757236365910867_wp,  4.1021714214668128_wp,  &
      & 2.3776130524831411_wp,  1.6969724987740866_wp,  5.2327688986668139_wp,  &
      & 1.9486148363011413_wp,  0.4390675147070869_wp,  2.9999022491838123_wp,  &
      & 3.5312997625581413_wp,  0.4467415528495868_wp,  4.8114121395028135_wp,  &
      & 6.5089895990100421_wp,  5.2100409408535882_wp,  6.0066553789008132_wp,  &
      & 0.9001165013630412_wp,  3.6420787128610868_wp,  5.4413106648508132_wp,  &
      & 1.6012116650460413_wp,  5.6845471271780879_wp,  0.7675566847298124_wp], &
      & shape(xyz)) * aatoau
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[6.4411018522600_wp,    0.0492571261505_wp,    0.2192046129910_wp,  &
      & 0.0462076831749_wp,    6.6435057067500_wp,    0.1670513770770_wp,  &
      & 0.2262248220170_wp,   -0.9573234940220_wp,    6.7608039126200_wp], &
      & shape(lattice)) * aatoau
   type(dftd_parameter),parameter :: dparam_pbe = dftd_parameter ( &
      & s6=1.0_wp, s8=0.95948085_wp, a1=0.38574991_wp, a2=4.80688534_wp, s9=1.0_wp)
   real(wp),parameter :: g_a = 3.0_wp
   real(wp),parameter :: g_c = 2.0_wp
   real(wp),parameter :: wf  = 6.0_wp
   integer, parameter :: refqmode = p_refq_goedecker

   type(TMolecule) :: mol
   type(TEnvironment) :: env
   type(TLatticePoint) :: latp
   type(TNeighbourlist) :: neighList

   integer, allocatable :: neighs(:)
   real(wp), allocatable :: trans(:, :), trans2(:, :)

   integer :: ii, jj
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp) :: energy, sigma(3, 3), sdum(3, 3)
   real(wp), allocatable :: gradient(:, :), gdum(:, :)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: q(:), dqdr(:,:,:), dqdL(:,:,:)

   call init(env)

   call init(mol, at, xyz)

   call d4init(g_a,g_c,refqmode)
   allocate(neighs(nat))
   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat))
   allocate(gradient(3, nat), gdum(3, nat))

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call init(latp, env, mol, 15.0_wp)  ! Fixed cutoff
   call latp%getLatticePoints(trans, 15.0_wp)
   call init(neighList, len(mol))
   call neighList%generate(env, mol%xyz, 15.0_wp, trans, .false.)

   call neighlist%getNeighs(neighs, 15.0_wp)
   call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
      & cn, dcndr, dcndL)
   call d4_atm_gradient(mol, trans, dparam_pbe, g_a, g_c, wf, 15.0_wp, &
      & cn, dcndr, dcndL, energy, gradient, sigma)

   call assert_close(energy, 0.44070137748975E-03_wp, thr)
   call assert_close(norm2(gradient), 0.54976067894467E-04_wp, thr)

   call assert_close(gradient(1, 3), -0.44417914354499E-05_wp, thr)
   call assert_close(gradient(2, 7), -0.11543946516301E-05_wp, thr)
   call assert_close(gradient(3,12), -0.76778342664578E-05_wp, thr)

   call assert_close(sigma(1,1), -0.18394332445576E-03_wp, thr)
   call assert_close(sigma(2,1), -0.33093630257770E-04_wp, thr)
   call assert_close(sigma(3,2),  0.34563088053137E-04_wp, thr)

   ! check numerical gradient
   do ii = 1, nat, 5
      do jj = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
            & cn, dcndr, dcndL)
         call d4_atm_gradient(mol, trans, dparam_pbe, g_a, g_c, wf, 15.0_wp, &
            & cn, dcndr, dcndL, er, gdum, sdum)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
            & cn, dcndr, dcndL)
         call d4_atm_gradient(mol, trans, dparam_pbe, g_a, g_c, wf, 15.0_wp, &
            & cn, dcndr, dcndL, el, gdum, sdum)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call assert_close(gradient(jj, ii), (er - el)*step2, thr)
      end do
   end do

!   ! check numerical strain derivatives
!   trans2 = trans
!   eps = unity
!   do ii = 1, 3
!      do jj = 1, ii
!         eps(jj, ii) = eps(jj, ii) + step
!         mol%lattice(:, :) = matmul(eps, lattice)
!         mol%xyz(:, :) = matmul(eps, xyz)
!         trans(:, :) = matmul(eps, trans2)
!         call mol%update
!         !call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
!            !& cn, dcndr, dcndL)
!         call d4_atm_gradient(mol, neighs, neighlist, dparam_pbe, g_a, g_c, wf, &
!            & cn, dcndr, dcndL, er, gdum, sdum)
!
!         eps(jj, ii) = eps(jj, ii) - 2*step
!         mol%lattice(:, :) = matmul(eps, lattice)
!         mol%xyz(:, :) = matmul(eps, xyz)
!         trans(:, :) = matmul(eps, trans2)
!         call mol%update
!         call neighList%update(mol%xyz, trans)
!         !call getCoordinationNumber(mol, neighs, neighList, cnType%cov, &
!            !& cn, dcndr, dcndL)
!         call d4_atm_gradient(mol, neighs, neighlist, dparam_pbe, g_a, g_c, wf, &
!            & cn, dcndr, dcndL, el, gdum, sdum)
!
!         eps(jj, ii) = eps(jj, ii) + step
!
!         call assert_close(sigma(jj, ii), (er - el)*step2, thr)
!      end do
!   end do

   call terminate(afail)
end subroutine test_dftd4_pbc3d_threebody_latp
