subroutine test_eeq_water
   use assertion
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_coulomb_gaussian, only : TGaussianSmeared, init
   use xtb_disp_coordinationnumber, only : getCoordinationNumber, cnType
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init, len
   use xtb_type_latticepoint, only : TLatticePoint, init
   use xtb_xtb_eeq
   implicit none
   real(wp), parameter :: thr = 1.0e-10_wp

   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp), parameter :: xyz(3,nat) = reshape([&
      & 0.00000000000000_wp,  0.00000000000000_wp, -0.73578586109551_wp, &
      & 1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp, &
      &-1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp],&
      & shape(xyz))
   real(wp), parameter :: chi(1, 2) = reshape(&
      & [1.56866440_wp, 1.23695041_wp], [1, 2])
   real(wp), parameter :: gam(1, 2) = reshape(&
      & [0.03151644_wp,-0.35015861_wp], [1, 2])
   real(wp), parameter :: kcn(1, 2) = reshape(&
      & [0.11689703_wp, 0.04916110_wp], [1, 2])
   real(wp), parameter :: alp(1, 2) = reshape(&
      & [1.23166285_wp, 0.55159092_wp], [1, 2])

   type(TMolecule) :: mol
   type(TLatticePoint) :: latp
   type(TEnvironment) :: env
   type(TGaussianSmeared) :: coulomb
   type(TENEquilibration) :: eeq

   logical :: exitRun
   real(wp) :: energy, sigma(3,3)
   real(wp), allocatable :: trans(:,:)
   real(wp), allocatable :: cn(:), dcndr(:,:,:), dcndL(:,:,:)
   real(wp), allocatable :: q(:), dqdr(:,:,:), dqdL(:,:,:)
   real(wp), allocatable :: gradient(:,:), qsh(:)

   call init(env)

   call init(mol, at, xyz)
   call init(latp, env, mol, 40.0_wp)  ! Fixed cutoff
   call env%checkpoint("Lattice point init failed")
   call init(coulomb, env, mol, alp)
   call env%checkpoint("Coulomb evaluator init failed")
   call init(eeq, env, chi, kcn, gam)
   call env%checkpoint("EN-Equilibration init failed")

   allocate(gradient(3, nat))
   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat))
   allocate(q(nat), dqdr(3, nat, nat), dqdL(3, 3, nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call latp%getLatticePoints(trans, 40.0_wp)

   call getCoordinationNumber(mol, trans, 40.0_wp, cnType%erf, &
      & cn, dcndr, dcndL)

   call assert_close(cn(1),1.9895544848535_wp,thr)
   call assert_close(cn(2),cn(3),             thr)
   call assert_close(cn(1),cn(2)+cn(3),       thr)

   call assert_close(dcndr(3,1,1), 0.80973198569003E-01_wp,thr)
   call assert_close(dcndr(1,2,1),-0.52891163425093E-01_wp,thr)
   call assert_close(dcndr(3,1,3), 0.40486599284501E-01_wp,thr)

   call eeq%chargeEquilibration(env, mol, coulomb, cn, dcndr, dcndL, &
      & energy, gradient, sigma, q, qsh, dqdr, dqdL)
   call env%checkpoint("Charge equilibration failed")

   ! test electrostatic energy
   call assert_close(energy,-0.64308088326667E-01_wp,thr)

   ! test molecular gradient of ES, also check symmetry
   call assert_close(gradient(3,1),-0.44053032330503E-01_wp,thr)
   call assert_close(gradient(3,2),gradient(3,3),           thr)
   call assert_close(gradient(1,2), 0.18102071270235E-01_wp,thr)

   ! test for charge constraint
   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1),-0.59582744684480_wp,thr)

   ! test derivative of partial charges
   call assert_close(dqdr(3,1,1),-0.41466764014389E+00_wp,thr)
   call assert_close(dqdr(1,2,1), 0.17196993288918E+00_wp,thr)
   call assert_close(dqdr(1,3,2), 0.18631993794394E-01_wp,thr)
   call assert_close(dqdr(2,1,3), 0.00000000000000E+00_wp,thr)

   call terminate(afail)
end subroutine test_eeq_water


subroutine test_eeq_ewald
   use assertion
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_coulomb_gaussian, only : TGaussianSmeared, init
   use xtb_disp_coordinationnumber, only : getCoordinationNumber, cnType, &
      & cutCoordinationNumber
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init, len
   use xtb_type_latticepoint, only : TLatticePoint, init
   use xtb_xtb_eeq
   use xtb_disp_ncoord
   use xtb_pbc_tools
   implicit none
   real(wp),parameter :: thr = 1.0e-9_wp
   real(wp),parameter :: thr2 = 1.0e-6_wp
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
   real(wp), allocatable :: xyz(:, :)
   real(wp), parameter :: chi(1, 2) = reshape(&
      & [ 1.12306840_wp, 1.56866440_wp], [1, 2])
   real(wp), parameter :: gam(1, 2) = reshape(&
      & [ 0.13155266_wp, 0.03151644_wp], [1, 2])
   real(wp), parameter :: kcn(1, 2) = reshape(&
      & [-0.04964052_wp, 0.11689703_wp], [1, 2])
   real(wp), parameter :: alp(1, 2) = reshape(&
      & [ 1.54075036_wp, 1.23166285_wp], [1, 2])

   type(TMolecule) :: mol
   type(TLatticePoint) :: latp
   type(TEnvironment) :: env
   type(TGaussianSmeared) :: coulomb
   type(TENEquilibration) :: eeq

   logical :: exitRun
   integer :: ii, jj, kk
   real(wp) :: energy, sigma(3,3)
   real(wp) :: er, el, eps(3, 3)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp), parameter :: step = 1.0e-6_wp, step2 = 0.5_wp/step
   real(wp), allocatable :: trans(:,:), trans2(:, :)
   real(wp), allocatable :: cn(:), dcndr(:,:,:), dcndL(:,:,:)
   real(wp), allocatable :: q(:), dqdr(:,:,:), dqdL(:,:,:)
   real(wp), allocatable :: gradient(:,:), qsh(:)
   real(wp), allocatable :: qr(:), ql(:)

   allocate(xyz(3, nat))
   call coord_trafo(nat, lattice, abc, xyz)

   call init(env)

   call init(mol, at, xyz, lattice=lattice)
   call init(latp, env, mol, 40.0_wp)  ! Fixed cutoff
   call env%checkpoint("Lattice point init failed")
   call init(coulomb, env, mol, alp)
   call env%checkpoint("Coulomb evaluator init failed")
   call init(eeq, env, chi, kcn, gam)
   call env%checkpoint("EN-Equilibration init failed")

   allocate(cn(nat), dcndr(3,nat,nat), dcndL(3,3,nat))
   allocate(q(nat), dqdr(3,nat,nat), dqdL(3,3,nat))
   allocate(gradient(3,nat))

   gradient = 0.0_wp
   energy = 0.0_wp
   sigma = 0.0_wp

   call latp%getLatticePoints(trans, 40.0_wp)

   call getCoordinationNumber(mol, trans, 40.0_wp, cnType%erf, &
      & cn, dcndr, dcndL)

   call assert_close(cn(2),3.6864725130236_wp,thr)
   call assert_close(cn(5),1.0523558225297_wp,thr)
   call assert_close(cn(6),1.1699488478421_wp,thr)

   call assert_close(dcndr(2,3,2), 0.24281192795725E-02_wp,thr)
   call assert_close(dcndr(1,6,6), 0.42240789876965E+00_wp,thr)
   call assert_close(dcndr(1,1,2), 0.71913132896636E-01_wp,thr)

   call assert_close(dcndL(1,3,4),-0.36068395425405_wp,thr)
   call assert_close(dcndL(3,3,1),-0.92808088688242_wp,thr)
   call assert_close(dcndL(2,1,3),-0.44871260782914_wp,thr)

   call eeq%chargeEquilibration(env, mol, coulomb, cn, dcndr, dcndL, &
      & energy, gradient, sigma, q, qsh, dqdr, dqdL)
   call env%checkpoint("Charge equilibration failed")

   call assert_close(energy,-0.13262207653027_wp,thr)

   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1), 0.56976165345597_wp,thr)
   call assert_close(q(3),-0.24549792737727_wp,thr)
   call assert_close(q(4),-0.30064628677070_wp,thr)

   if (afail > 0) call terminate(afail)

   ! check numerical gradient
   allocate(qr(nat), ql(nat))
   do ii = 1, nat
      do jj = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call coulomb%update(env, mol)
         call getCoordinationNumber(mol, trans, 40.0_wp, cnType%erf, &
            & cn, dcndr, dcndL)
         call eeq%chargeEquilibration(env, mol, coulomb, cn, dcndr, dcndL, &
            & er, qat=qr)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call coulomb%update(env, mol)
         call getCoordinationNumber(mol, trans, 40.0_wp, cnType%erf, &
            & cn, dcndr, dcndL)
         call eeq%chargeEquilibration(env, mol, coulomb, cn, dcndr, dcndL, &
            & el, qat=ql)
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call assert_close(gradient(jj, ii), (er - el)*step2, thr)
         do kk = 1, nat
            call assert_close(dqdr(jj, ii, kk), (qr(kk) - ql(kk))*step2, thr)
         end do
      end do
   end do

   if (afail > 0) call terminate(afail)

   ! check numerical strain derivatives
   trans2 = trans
   eps = unity
   do ii = 1, 3
      do jj = 1, ii
         er = 0.0_wp
         el = 0.0_wp
         eps(jj, ii) = eps(jj, ii) + step
         mol%lattice(:, :) = matmul(eps, lattice)
         mol%xyz(:, :) = matmul(eps, xyz)
         trans(:, :) = matmul(eps, trans2)
         call mol%update
         call coulomb%update(env, mol)
         call getCoordinationNumber(mol, trans, 40.0_wp, cnType%erf, &
            & cn, dcndr, dcndL)
         call eeq%chargeEquilibration(env, mol, coulomb, cn, dcndr, dcndL, &
            & er, qat=qr)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%lattice(:, :) = matmul(eps, lattice)
         mol%xyz(:, :) = matmul(eps, xyz)
         trans(:, :) = matmul(eps, trans2)
         call mol%update
         call coulomb%update(env, mol)
         call getCoordinationNumber(mol, trans, 40.0_wp, cnType%erf, &
            & cn, dcndr, dcndL)
         call eeq%chargeEquilibration(env, mol, coulomb, cn, dcndr, dcndL, &
            & el, qat=ql)

         eps(jj, ii) = eps(jj, ii) + step

         call assert_close(sigma(jj, ii), (er - el)*step2, thr2)
         do kk = 1, nat
            call assert_close(dqdL(jj, ii, kk), (qr(kk) - ql(kk))*step2, thr2)
         end do
      end do
   end do

   if (afail > 0) call terminate(afail)

   ! reset energy, gradient and stress
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   trans(:, :) = trans2

   mol%lattice(:, :) = lattice
   mol%xyz(:, :) = xyz
   call mol%update
   call coulomb%update(env, mol)
   call getCoordinationNumber(mol, trans, 40.0_wp, cnType%erf, &
      & cn, dcndr, dcndL)
   call cutCoordinationNumber(nat, cn, dcndr, dcndL, maxCN=8.0_wp)

   call assert_close(cn(2),3.6735104771649_wp,thr)
   call assert_close(cn(5),1.0517307940776_wp,thr)
   call assert_close(cn(6),1.1692040350311_wp,thr)

   call assert_close(dcndr(2,3,2), 0.23960452277676E-02_wp,thr)
   call assert_close(dcndr(1,6,6), 0.42195185201461E+00_wp,thr)
   call assert_close(dcndr(1,1,2), 0.70963201989459E-01_wp,thr)

   call assert_close(dcndL(1,3,4),-0.36002844990181_wp,thr)
   call assert_close(dcndL(3,3,1),-0.92527392943064_wp,thr)
   call assert_close(dcndL(2,1,3),-0.44767072306065_wp,thr)

   call eeq%chargeEquilibration(env, mol, coulomb, cn, dcndr, dcndL, &
      & energy, gradient, sigma, q, qsh, dqdr, dqdL)
   call env%checkpoint("Charge equilibration failed")

   call assert_close(energy,-0.13281818007229_wp,thr)

   call assert_close(gradient(2,1), 0.97364690950909E-02_wp,thr)
   call assert_close(gradient(1,3),-0.32536878267545E-02_wp,thr)
   call assert_close(gradient(3,5), 0.72645758604223E-03_wp,thr)
   call assert_close(gradient(1,6), 0.17141994546793E-01_wp,thr)

   call assert_close(sigma(1,1),-0.71215570232998E-01_wp,thr)
   call assert_close(sigma(2,3),-0.15833762501521E-01_wp,thr)
   call assert_close(sigma(3,1), 0.25423858284897E-01_wp,thr)

   call assert_close(sum(q),0.0_wp,            thr)
   call assert_close(q(1), 0.57003309593516_wp,thr)
   call assert_close(q(3),-0.24580661514550_wp,thr)
   call assert_close(q(4),-0.30087488142104_wp,thr)

   call assert_close(dqdr(1,4,2),-0.18544856267141E-01_wp,thr)
   call assert_close(dqdr(1,2,2), 0.58677578747262E-01_wp,thr)
   call assert_close(dqdr(3,1,4),-0.78746078692863E-01_wp,thr)

   call assert_close(dqdL(2,3,2),-0.21027678205985E-01_wp,thr)
   call assert_close(dqdL(3,2,1), 0.83475462871348E-01_wp,thr)
   call assert_close(dqdL(1,1,5), 0.81210808302498E-02_wp,thr)

   call terminate(afail)
end subroutine test_eeq_ewald
