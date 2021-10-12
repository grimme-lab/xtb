! This file is part of xtb.
! SPDX-Identifier: LGPL-3.0-or-later
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

module test_eeq
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed, &
      & skip_test
   implicit none
   private

   public :: collect_eeq

contains

!> Collect all exported unit tests
subroutine collect_eeq(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("water", test_eeq_water), &
      new_unittest("ewald", test_eeq_ewald), &
      new_unittest("gbsa", test_eeq_model_gbsa), &
      new_unittest("hbond", test_eeq_model_hbond), &
      new_unittest("salt", test_eeq_model_salt) &
      ]

end subroutine collect_eeq


subroutine test_eeq_water(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_coulomb_gaussian, only : TGaussianSmeared, init
   use xtb_disp_coordinationnumber, only : getCoordinationNumber, cnType
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init, len
   use xtb_type_latticepoint, only : TLatticePoint, init
   use xtb_xtb_eeq
   type(error_type), allocatable, intent(out) :: error
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

   call check(error, cn(1),1.9895544848535_wp,thr=thr)
   call check(error, cn(2),cn(3),             thr=thr)
   call check(error, cn(1),cn(2)+cn(3),       thr=thr)

   call check(error, dcndr(3,1,1), 0.80973198569003E-01_wp,thr=thr)
   call check(error, dcndr(1,2,1),-0.52891163425093E-01_wp,thr=thr)
   call check(error, dcndr(3,1,3), 0.40486599284501E-01_wp,thr=thr)

   call eeq%chargeEquilibration(env, mol, coulomb, cn, dcndr, dcndL, &
      & energy, gradient, sigma, q, qsh, dqdr, dqdL)
   call env%checkpoint("Charge equilibration failed")

   ! test electrostatic energy
   call check(error, energy,-0.64308088326667E-01_wp,thr=thr)

   ! test molecular gradient of ES, also check symmetry
   call check(error, gradient(3,1),-0.44053032330503E-01_wp,thr=thr)
   call check(error, gradient(3,2),gradient(3,3),           thr=thr)
   call check(error, gradient(1,2), 0.18102071270235E-01_wp,thr=thr)

   ! test for charge constraint
   call check(error, sum(q),0.0_wp,            thr=thr)
   call check(error, q(1),-0.59582744684480_wp,thr=thr)

   ! test derivative of partial charges
   call check(error, dqdr(3,1,1),-0.41466764014389E+00_wp,thr=thr)
   call check(error, dqdr(1,2,1), 0.17196993288918E+00_wp,thr=thr)
   call check(error, dqdr(1,3,2), 0.18631993794394E-01_wp,thr=thr)
   call check(error, dqdr(2,1,3), 0.00000000000000E+00_wp,thr=thr)

end subroutine test_eeq_water


subroutine test_eeq_ewald(error)
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
   type(error_type), allocatable, intent(out) :: error
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

   call check(error, cn(2),3.6864725130236_wp,thr=thr)
   call check(error, cn(5),1.0523558225297_wp,thr=thr)
   call check(error, cn(6),1.1699488478421_wp,thr=thr)

   call check(error, dcndr(2,3,2), 0.24281192795725E-02_wp,thr=thr)
   call check(error, dcndr(1,6,6), 0.42240789876965E+00_wp,thr=thr)
   call check(error, dcndr(1,1,2), 0.71913132896636E-01_wp,thr=thr)

   call check(error, dcndL(1,3,4),-0.36068395425405_wp,thr=thr)
   call check(error, dcndL(3,3,1),-0.92808088688242_wp,thr=thr)
   call check(error, dcndL(2,1,3),-0.44871260782914_wp,thr=thr)

   call eeq%chargeEquilibration(env, mol, coulomb, cn, dcndr, dcndL, &
      & energy, gradient, sigma, q, qsh, dqdr, dqdL)
   call env%checkpoint("Charge equilibration failed")

   call check(error, energy,-0.13262207653027_wp,thr=thr)

   call check(error, sum(q),0.0_wp,            thr=thr)
   call check(error, q(1), 0.56976165345597_wp,thr=thr)
   call check(error, q(3),-0.24549792737727_wp,thr=thr)
   call check(error, q(4),-0.30064628677070_wp,thr=thr)

   if (allocated(error)) return

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
         call check(error, gradient(jj, ii), (er - el)*step2, thr=thr)
         do kk = 1, nat
            call check(error, dqdr(jj, ii, kk), (qr(kk) - ql(kk))*step2, thr=thr)
         end do
      end do
   end do

   if (allocated(error)) return

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

         call check(error, sigma(jj, ii), (er - el)*step2, thr=thr2)
         do kk = 1, nat
            call check(error, dqdL(jj, ii, kk), (qr(kk) - ql(kk))*step2, thr=thr2)
         end do
      end do
   end do

   if (allocated(error)) return

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

   call check(error, cn(2),3.6735104771649_wp,thr=thr)
   call check(error, cn(5),1.0517307940776_wp,thr=thr)
   call check(error, cn(6),1.1692040350311_wp,thr=thr)

   call check(error, dcndr(2,3,2), 0.23960452277676E-02_wp,thr=thr)
   call check(error, dcndr(1,6,6), 0.42195185201461E+00_wp,thr=thr)
   call check(error, dcndr(1,1,2), 0.70963201989459E-01_wp,thr=thr)

   call check(error, dcndL(1,3,4),-0.36002844990181_wp,thr=thr)
   call check(error, dcndL(3,3,1),-0.92527392943064_wp,thr=thr)
   call check(error, dcndL(2,1,3),-0.44767072306065_wp,thr=thr)

   call eeq%chargeEquilibration(env, mol, coulomb, cn, dcndr, dcndL, &
      & energy, gradient, sigma, q, qsh, dqdr, dqdL)
   call env%checkpoint("Charge equilibration failed")

   call check(error, energy,-0.13281818007229_wp,thr=thr)

   call check(error, gradient(2,1), 0.97364690950909E-02_wp,thr=thr)
   call check(error, gradient(1,3),-0.32536878267545E-02_wp,thr=thr)
   call check(error, gradient(3,5), 0.72645758604223E-03_wp,thr=thr)
   call check(error, gradient(1,6), 0.17141994546793E-01_wp,thr=thr)

   call check(error, sigma(1,1),-0.71215570232998E-01_wp,thr=thr)
   call check(error, sigma(2,3),-0.15833762501521E-01_wp,thr=thr)
   call check(error, sigma(3,1), 0.25423858284897E-01_wp,thr=thr)

   call check(error, sum(q),0.0_wp,            thr=thr)
   call check(error, q(1), 0.57003309593516_wp,thr=thr)
   call check(error, q(3),-0.24580661514550_wp,thr=thr)
   call check(error, q(4),-0.30087488142104_wp,thr=thr)

   call check(error, dqdr(1,4,2),-0.18544856267141E-01_wp,thr=thr)
   call check(error, dqdr(1,2,2), 0.58677578747262E-01_wp,thr=thr)
   call check(error, dqdr(3,1,4),-0.78746078692863E-01_wp,thr=thr)

   call check(error, dqdL(2,3,2),-0.21027678205985E-01_wp,thr=thr)
   call check(error, dqdL(3,2,1), 0.83475462871348E-01_wp,thr=thr)
   call check(error, dqdL(1,1,5), 0.81210808302498E-02_wp,thr=thr)

end subroutine test_eeq_ewald

subroutine test_eeq_model_gbsa(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
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
   type(error_type), allocatable, intent(out) :: error

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
   call check(error, cn(1),3.9364408722599_wp,thr=thr)
   call check(error, cn(6),1.9734817313511_wp,thr=thr)

   ! test CN derivative, correct summation of diagonals
   call check(error, dcndr(2,1,5),-0.30966754386801E-01_wp,thr=thr)
   call check(error, dcndr(1,3,2),-0.00000000000000E+00_wp,thr=thr)
   call check(error, dcndr(2,1,6), 0.24290126624329E-05_wp,thr=thr)

   call new_charge_model_2019(chrgeq,mol%n,mol%at)
   call eeq_chrgeq(mol,env,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,es,ges,sigma, &
      &            .false.,.true.,.true.)

   ! test electrostatic energy
   call check(error, es,-0.10559262715710E-01_wp,thr=thr)

   ! test molecular gradient of ES, also check symmetry
   call check(error, ges(2,5), 0.45485199335534E-04_wp,thr=thr)
   call check(error, ges(3,3),-0.62757274919919E-04_wp,thr=thr)
   call check(error, ges(1,2), 0.70062774412852E-05_wp,thr=thr)

   ! test for charge constraint
   call check(error, sum(q),0.0_wp,            thr=thr)
   call check(error, q(1),0.21873929670215_wp,thr=thr)

   ! test derivative of partial charges
   call check(error, dqdr(2,1,6),-0.32458122507411E-02_wp,thr=thr)
   call check(error, dqdr(3,2,1), 0.00000000000000E+00_wp,thr=thr)
   call check(error, dqdr(1,3,2), 0.27025784902288E-03_wp,thr=thr)
   call check(error, dqdr(2,8,5), 0.36651303109315E-03_wp,thr=thr)

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
   call check(error, es,-0.18723100944150E-01_wp,thr=thr)

   ! test molecular gradient of ES, also check symmetry
   call check(error, ges(2,5), 0.15098175186158E-03_wp,thr=thr)
   call check(error, ges(3,3),-0.42131892534669E-03_wp,thr=thr)
   call check(error, ges(1,2), 0.28893943115425E-04_wp,thr=thr)

   ! test for charge constraint
   call check(error, sum(q),0.0_wp,thr=thr)
   call check(error, q(1), 0.22653931001924E+00_wp,thr=thr)
   call check(error, q(2),-0.20322834228298E-01_wp,thr=thr)
   call check(error, q(3),-0.26032042423717E-01_wp,thr=thr)
   call check(error, q(4),-0.26032042423717E-01_wp,thr=thr)
   call check(error, q(5),-0.25976606121313E-01_wp,thr=thr)
   call check(error, q(6),-0.18580233359855E+00_wp,thr=thr)
   call check(error, q(7), 0.28815815576667E-01_wp,thr=thr)
   call check(error, q(8), 0.28810733199693E-01_wp,thr=thr)

   !call check(error, dqdr(2,1,6),-0.21098409244566E-02_wp,thr=thr)
   !call check(error, dqdr(3,2,1), 0.00000000000000E+00_wp,thr=thr)
   !call check(error, dqdr(1,3,2), 0.37653484040321E-03_wp,thr=thr)
   !call check(error, dqdr(2,8,5), 0.82628766939384E-03_wp,thr=thr)

   call mol%deallocate

end subroutine test_eeq_model_gbsa

subroutine test_eeq_model_hbond(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   type(error_type), allocatable, intent(out) :: error
   call skip_test(error, "Not implemented")
end subroutine test_eeq_model_hbond

subroutine test_eeq_model_salt(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_convert, only : aatoau
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
   type(error_type), allocatable, intent(out) :: error

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
   call check(error, es,-0.64149462344256E-01_wp,thr=thr)

   ! test molecular gradient of ES, also check symmetry
   call check(error, ges(2,7),-0.52237685114757E-03_wp,thr=thr)
   call check(error, ges(3,3),-0.10347464904462E-03_wp,thr=thr)
   call check(error, ges(1,9),-0.19202182317373E-03_wp,thr=thr)

   ! test for charge constraint
   call check(error, sum(q),0.0_wp,            thr=thr)
   call check(error, q(1),-0.20644568417986E+00_wp,thr=thr)
   call check(error, q(8), 0.83875312119020E-01_wp,thr=thr)

   ! test derivative of partial charges
   call check(error, dqdr(2,2,6),-0.93834552520241E-02_wp,thr=thr)
   call check(error, dqdr(3,3,1),-0.24876158910325E-03_wp,thr=thr)
   call check(error, dqdr(1,5,2),-0.15306250365591E-01_wp,thr=thr)
   call check(error, dqdr(2,8,7),-0.92052318691156E-02_wp,thr=thr)

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
   call check(error, es,-0.73041696453862E-01_wp,thr=thr)

   ! test molecular gradient of ES, also check symmetry
   call check(error, ges(2,7),-0.54144355306156E-03_wp,thr=thr)
   call check(error, ges(3,3),-0.13994375398115E-03_wp,thr=thr)
   call check(error, ges(1,9),-0.72309374268259E-04_wp,thr=thr)

   ! test for charge constraint
   call check(error, sum(q),0.0_wp,thr=thr)
   call check(error, q( 1),-0.22628341725089E+00_wp,thr=thr)
   call check(error, q( 2), 0.17181969742365E+00_wp,thr=thr)
   call check(error, q( 3),-0.31688852764536E-01_wp,thr=thr)
   call check(error, q( 4), 0.86131431150190E-01_wp,thr=thr)
   call check(error, q( 5),-0.31623215923196E-01_wp,thr=thr)
   call check(error, q( 6), 0.86034414210273E-01_wp,thr=thr)
   call check(error, q( 7),-0.31614086140679E-01_wp,thr=thr)
   call check(error, q( 8), 0.86132083253082E-01_wp,thr=thr)
   call check(error, q( 9),-0.22624788057207E+00_wp,thr=thr)
   call check(error, q(10), 0.17181804147276E+00_wp,thr=thr)
   call check(error, q(11), 0.17180510824691E+00_wp,thr=thr)
   call check(error, q(12),-0.22628332310549E+00_wp,thr=thr)

   ! test derivative of partial charges (still incorrect)
   !call check(error, dqdr(2,2,6),-0.88602533712134E-02_wp,thr=thr)
   !call check(error, dqdr(3,3,1),-0.38467542758235E-03_wp,thr=thr)
   !call check(error, dqdr(1,5,2),-0.13407102470387E-01_wp,thr=thr)
   !call check(error, dqdr(2,8,7),-0.84239153450647E-02_wp,thr=thr)

   call mol%deallocate

end subroutine test_eeq_model_salt
end module test_eeq
