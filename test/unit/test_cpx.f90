
module test_cpx
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check, test_failed, skip_test
   use xtb_features, only : get_xtb_feature
   implicit none
   private

   public :: collect_cpx

contains

!> Collect all exported unit tests
subroutine collect_cpx(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("solvation", test_cpx_solv) &
      ]

end subroutine collect_cpx

subroutine test_cpx_solv(error)
   use xtb_solv_cpx, only : TCpcmx
   use mctc_env, only: err_type => error_type, wp

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment
   use xtb_solv_input
   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_kernel

   type(error_type), allocatable, intent(out) :: error
   type(err_type), allocatable :: err
   type(Tcpcmx) :: cpcmx
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   real(wp) :: energy, sigma(3,3)
   real(wp) :: hl_gap
   real(wp) :: energy_gas, total_energy
   real(wp),allocatable :: gradient(:,:)
   real(wp),parameter :: thr = 5.0e-5_wp
   integer, parameter :: nat = 11
   integer, parameter :: at(nat) = [7,1,1,1,15,6,6,6,7,7,7]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[4.84112686339230_wp,-1.15671648374837_wp, 0.00000000000000_wp, &
      & 5.87217277985151_wp,-0.70323701391366_wp, 1.54511189584035_wp, &
      & 5.87217277985151_wp,-0.70323701391366_wp,-1.54511189584035_wp, &
      & 4.61871932563459_wp,-3.05696083228965_wp, 0.00000000000000_wp, &
      & 0.23751589819013_wp, 1.28722286021158_wp, 0.00000000000000_wp, &
      & 0.00161949529588_wp,-0.89075264835377_wp, 2.57562301980395_wp, &
      &-3.06580664031769_wp, 2.16487781673832_wp, 0.00000000000000_wp, &
      & 0.00161949529588_wp,-0.89075264835377_wp,-2.57562301980395_wp, &
      &-5.06585309511506_wp, 3.00827014874972_wp, 0.00000000000000_wp, &
      & 0.00161949529588_wp,-2.15043086045241_wp, 4.34183098346579_wp, &
      & 0.00161949529588_wp,-2.15043086045241_wp,-4.34183098346579_wp], shape(xyz))
   type(scc_options),parameter :: opt = scc_options( &
      &  prlevel = 2, maxiter = 30, acc = 1.0_wp, etemp = 300.0_wp, grad = .true.,&
      &  solvent = "h2o")
   integer :: ich

   character(len=*), parameter :: cpxsolvent="water"

   call skip_test(error, 'Excluding CPCM-X test.')
   return
   if (.not.get_xtb_feature('cpcmx')) then
      call skip_test(error, 'CPCM-X libary not available.')
      return
   end if
   call init(env)

   call init(mol, at, xyz)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=2)
   call newWavefunction(env, mol, calc, chk)
   call addSolvationModel(env, calc, TSolvInput(solvent='infinity', alpb=.false., cosmo=.true.))

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)
   call cpcmx%setup(env,cpxsolvent)
   deallocate(calc%solvation)
   call calc%singlepoint(env, mol, chk, 0, .false., energy_gas, gradient, sigma, &
      & hl_gap, res)
   call cpcmx%calc_solv(env,cpxsolvent,energy_gas,0.4_wp,298.15_wp,500,0.0001_wp,total_energy)
   call cpcmx%print(.true.)
   call check_(error, total_energy-energy_gas, -0.1555665E-01_wp, thr=thr)
   call check_(error, total_energy, -22.00017093824971_wp, thr=thr)

end subroutine test_cpx_solv
end module test_cpx
