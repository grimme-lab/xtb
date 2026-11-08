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

module test_gfnff
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check, test_failed
   implicit none
   private

   public :: collect_gfnff

contains

!> Collect all exported unit tests
subroutine collect_gfnff(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("sp", test_gfnff_sp), &
      new_unittest("eeq-atom-order", test_eeq_atom_order), &
      new_unittest("harmonic", test_gfnff_harmonic), &
      new_unittest("hb", test_gfnff_hb), &
      new_unittest("gbsa", test_gfnff_gbsa), &
      new_unittest("mindless", test_gfnff_mindless_basic), &
      new_unittest("mindless-solvation", test_gfnff_mindless_solvation), &
      new_unittest("scaleup", test_gfnff_scaleup), &
      new_unittest("pdb", test_gfnff_pdb), &
      new_unittest("sdf", test_gfnff_sdf), &
      new_unittest("pbc", test_gfnff_pbc), &
      new_unittest("fragmentation", test_gfnff_fragmentation), &
      new_unittest("Ln_An", test_gfnff_LnAn_H) &
      ]

end subroutine collect_gfnff



!> N-methylacetamide: parameter generation must follow atoms, not list positions.
subroutine test_eeq_atom_order(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init
   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 12
   integer, parameter :: at(nat) = [6,6,8,7,6,1,1,1,1,1,1,1]
   ! Coordinates in bohr; the same molecule in each ordering.
   real(wp), parameter :: xyz(3,nat) = reshape([ &
     & 3.622312140781_wp, 0.584150709785_wp, 0.403823853883_wp, &
     & 0.841690132958_wp, 0.709612631882_wp, 0.875537937257_wp, &
     & 0.001690216911_wp, 2.103743880637_wp, 2.615851423658_wp, &
     & -0.857946784445_wp, -0.712426316667_wp, -0.678850260895_wp, &
     & -3.584295528389_wp, -0.636030709585_wp, -0.295866267463_wp, &
     & 4.066614256909_wp, 1.536668293196_wp, -1.411023959546_wp, &
     & 4.677338535329_wp, 1.531874376674_wp, 1.951448115885_wp, &
     & 4.231725963854_wp, -1.420334469809_wp, 0.299405229455_wp, &
     & -0.129625612682_wp, -1.837950193426_wp, -2.128976201128_wp, &
     & -4.051091096309_wp, -1.315129049540_wp, 1.633832596409_wp, &
     & -4.536218537284_wp, -1.872131986606_wp, -1.699163157357_wp, &
     & -4.282193687633_wp, 1.327952833458_wp, -0.537763562819_wp], [3,nat])
   integer, parameter :: orders(nat,3) = reshape([ &
      & 1,2,3,4,5,6,7,8,9,10,11,12, &
      & 12,11,10,9,8,7,6,5,4,3,2,1, &
      & 4,1,2,3,5,6,7,8,9,10,11,12], [nat,3])
   real(wp) :: reference(nat)
   integer :: pass, i, order(nat)
   logical :: failed

   do pass = 1,size(orders,2)
      order = orders(:,pass)
      block
         type(TEnvironment) :: env
         type(TMolecule) :: mol
         type(TGFFCalculator) :: calc
         call init(env)
         call init(mol,at(order),xyz(:,order))
         call newGFFCalculator(env,mol,calc,'---',.false.)
         call env%check(failed)
         call check_(error,.not.failed)
         if (allocated(error)) return
         if (pass == 1) then
            reference = calc%topo%gameeq
         else
            do i = 1,nat
               call check_(error,calc%topo%gameeq(i),reference(order(i)),thr=1e-8_wp)
               if (allocated(error)) return
            end do
         end if

      end block
   end do
end subroutine test_eeq_atom_order


subroutine test_gfnff_sp(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_systools
   use xtb_solv_gbsa
   use xtb_type_environment
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_data
   use xtb_gfnff_param
   use xtb_gfnff_setup
   use xtb_gfnff_eg
   use xtb_gfnff_ini
   use xtb_gfnff_neighbourlist
   use xtb_setparam
   use xtb_setmod
   use xtb_disp_dftd3param
   use xtb_disp_dftd4
   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 8
   integer, parameter :: at(nat) = [7,15,9,9,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[1.50040286526241_wp,-2.88219140061585_wp, 0.00000000000000_wp, &
      & 0.00000000000000_wp, 1.90142164792207_wp, 0.00000000000000_wp, &
      &-0.02649585010919_wp,-5.05651406856634_wp, 0.00000000000000_wp, &
      &-1.39762821979929_wp, 4.65236211997835_wp, 0.00000000000000_wp, &
      & 2.62205170116282_wp,-3.14316635901963_wp, 1.53958066744940_wp, &
      &-1.46489869067775_wp, 0.78885483581631_wp, 1.94964934855945_wp, &
      & 2.62205170116282_wp,-3.14316635901963_wp,-1.53958066744940_wp, &
      &-1.46489869067775_wp, 0.78885483581631_wp,-1.94964934855945_wp],&
      & shape(xyz))
   logical, parameter :: restart = .false.

   type(TMolecule)     :: mol
   type(TEnvironment)  :: env
   type(scc_results)   :: res_gff
   type(TGFFCalculator) :: calc
   type(TGFFNeighbourList) :: nlist
   type(TBorn), allocatable :: solvation

   real(wp) :: etot
   real(wp), allocatable :: g(:,:)
   real(wp) :: sigma(3, 3)
   character(len=:),allocatable :: fnv
   integer  :: ipar
   real(wp) :: efield(3) = [0.0_wp, 0.0_wp, 0.0_wp]

   logical  :: exist

   call init(env)
   call init(mol,at,xyz)

   call delete_file('charges')
   call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)

   call env%checkpoint("GFN-FF parameter setup failed")

   allocate( g(3,mol%n), source = 0.0_wp )
 
   call check_(error, calc%neigh%nbond,6)
   call check_(error, calc%topo%nangl,6)
   call check_(error, calc%topo%ntors,1)

   g = 0.0_wp
   gff_print=.true.

   call gfnff_eg(env,mol,gff_print,mol%n,nint(mol%chrg),mol%at,mol%xyz,sigma, &
      & g,etot,res_gff,calc%param,calc%topo,calc%neigh,nlist,efield,solvation,.true.,calc%version, &
      & calc%accuracy)

   call check_(error, res_gff%e_total,-0.76480130328909_wp, thr=thr)
   call check_(error, res_gff%gnorm,   0.06237476775485_wp, thr=thr)
   call check_(error, res_gff%e_bond, -0.74131050087160_wp, thr=thr)
   call check_(error, res_gff%e_angl,  0.00633910615752_wp, thr=thr)
   call check_(error, res_gff%e_tors,  0.00004724445432_wp, thr=thr)
   call check_(error, res_gff%e_es,   -0.05070333195447_wp, thr=thr*10)
   call check_(error, res_gff%e_disp, -0.00224146422313_wp, thr=thr)
   call check_(error, res_gff%e_rep,   0.03086605590295_wp, thr=thr)
   call check_(error, res_gff%e_hb,   -0.00003142616658_wp, thr=thr)
   call check_(error, res_gff%e_xb,   -0.00776698664545_wp, thr=thr)
   call check_(error, res_gff%e_batm, -0.00000000000000_wp, thr=thr)

   call mol%deallocate

end subroutine test_gfnff_sp


subroutine test_gfnff_harmonic(error)
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_gfnff_param, only : gffVersion

   type(error_type), allocatable, intent(out) :: error

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res
   real(wp), parameter :: thr = 1000*epsilon(1.0_wp)
   real(wp) :: energy, hl_gap, sigma(3, 3), gradient(3, 3)
   logical :: exitRun

   call init(env)
   call getMolecule(mol, "h2o")
   call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false., &
      & gffVersion%harmonic2020)

   call env%check(exitRun)
   call check_(error, .not.exitRun)
   if (exitRun) return

   ! Harmonic GFN-FF stores its bond list in the neighbor data.
   call check_(error, calc%neigh%nbond, 2)
   call check_(error, allocated(calc%neigh%blist))
   call check_(error, .not.allocated(calc%topo%blist))

   call calc%singlepoint(env, mol, chk, 0, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call env%check(exitRun)
   call check_(error, .not.exitRun)
   call check_(error, ieee_is_finite(energy))
   call check_(error, all(ieee_is_finite(gradient)))
   call check_(error, energy, 0.00476277749486647_wp, thr=thr)
   call check_(error, norm2(gradient), 0.0478775550094384_wp, thr=thr)

end subroutine test_gfnff_harmonic

subroutine test_gfnff_hb(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_systools
   use xtb_type_environment
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_data
   use xtb_gfnff_param
   use xtb_gfnff_setup
   use xtb_gfnff_eg
   use xtb_gfnff_ini
   use xtb_setparam
   use xtb_setmod
   use xtb_disp_dftd3param
   use xtb_disp_dftd4
   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 7
   integer, parameter :: at(nat) = [6,8,1,1,8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[-5.78520874132429_wp,-1.92328475821000_wp,-0.02944611115854_wp, &
      & -5.57801768832583_wp, 0.17912532844037_wp, 0.72444143178660_wp, &
      & -4.27822256938673_wp,-2.74845397256109_wp,-1.13038073598642_wp, &
      & -7.47879539136783_wp,-2.97570121473211_wp, 0.39488815557786_wp, &
      & -0.83005338399036_wp, 2.43458470560665_wp,-0.78566331969245_wp, &
      & -0.74201439536855_wp, 4.04199055249898_wp, 0.09144422329636_wp, &
      & -2.44679415487233_wp, 1.69392751177087_wp,-0.27417668699116_wp],&
      & shape(xyz))
   logical, parameter :: restart = .false.

   type(TMolecule)     :: mol
   type(TEnvironment)  :: env
   type(TRestart)      :: chk
   type(scc_results)   :: res_gff
   type(TGFFCalculator) :: calc

   real(wp) :: etot, sigma(3,3), gap
   real(wp), allocatable :: g(:,:)
   character(len=:),allocatable :: fnv
   integer  :: ipar

   logical  :: exist

   call init(env)
   call init(mol,at,xyz)

   call delete_file('charges')
   call newGFFCalculator(env, mol, calc, '---', .false.)

   call env%checkpoint("GFN-FF parameter setup failed")

   allocate( g(3,mol%n), source = 0.0_wp )
 
   call check_(error, calc%neigh%nbond,5)
   call check_(error, calc%topo%nangl,4)
   call check_(error, calc%topo%ntors,1)

   g = 0.0_wp
   gff_print=.true.

   call calc%singlepoint(env, mol, chk, 1, .false., etot, g, sigma, gap, res_gff)

   call check_(error, res_gff%e_total,-0.949706677118_wp, thr=thr)
   call check_(error, res_gff%gnorm,   0.001152722069_wp, thr=thr)
   call check_(error, res_gff%e_bond, -0.856707644691_wp, thr=thr)
   call check_(error, res_gff%e_angl,  0.000579711987_wp, thr=thr)
   call check_(error, res_gff%e_tors,  0.000000008811_wp, thr=thr)
   call check_(error, res_gff%e_es,   -0.152313816530_wp, thr=thr*10)
   call check_(error, res_gff%e_disp, -0.001251669186_wp, thr=thr)
   call check_(error, res_gff%e_rep,   0.066881023899_wp, thr=thr)
   call check_(error, res_gff%e_hb,   -0.006894292337_wp, thr=thr)
   call check_(error, res_gff%e_xb,   -0.000000000000_wp, thr=thr)
   call check_(error, res_gff%e_batm, -0.000000000000_wp, thr=thr)

   call mol%deallocate

end subroutine test_gfnff_hb

subroutine test_gfnff_gbsa(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_systools
   use xtb_solv_input
   use xtb_type_environment
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_data
   use xtb_gfnff_param
   use xtb_gfnff_setup
   use xtb_gfnff_eg
   use xtb_gfnff_ini
   use xtb_setparam
   use xtb_setmod
   use xtb_disp_dftd3param
   use xtb_disp_dftd4
   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_main_setup, only : addSolvationModel
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 7
   integer, parameter :: at(nat) = [6,8,1,1,8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[-5.78520874132429_wp,-1.92328475821000_wp,-0.02944611115854_wp, &
      & -5.57801768832583_wp, 0.17912532844037_wp, 0.72444143178660_wp, &
      & -4.27822256938673_wp,-2.74845397256109_wp,-1.13038073598642_wp, &
      & -7.47879539136783_wp,-2.97570121473211_wp, 0.39488815557786_wp, &
      & -0.83005338399036_wp, 2.43458470560665_wp,-0.78566331969245_wp, &
      & -0.74201439536855_wp, 4.04199055249898_wp, 0.09144422329636_wp, &
      & -2.44679415487233_wp, 1.69392751177087_wp,-0.27417668699116_wp],&
      & shape(xyz))
   type(scc_options),parameter :: opt = scc_options( &
      &  prlevel = 2, maxiter = 30, acc = 1.0_wp, etemp = 300.0_wp, grad = .true.,&
      &  solvent = "h2o")
   logical, parameter :: restart = .false.

   type(TMolecule)     :: mol
   type(TEnvironment)  :: env
   type(TRestart)      :: chk
   type(scc_results)   :: res_gff
   type(TGFFCalculator) :: calc

   real(wp) :: etot, sigma(3,3), gap
   real(wp), allocatable :: g(:,:)
   character(len=:),allocatable :: fnv
   integer  :: ipar

   logical  :: exist

   call init(env)
   call init(mol,at,xyz)

   call delete_file('charges')
   call newGFFCalculator(env, mol, calc, '---', .false.)
   call addSolvationModel(env, calc, TSolvInput(solvent=opt%solvent, alpb=.false., kernel=gbKernel%still))

   call env%checkpoint("GFN-FF parameter setup failed")

   allocate( g(3,mol%n), source = 0.0_wp )
 
   call check_(error, calc%neigh%nbond,5)
   call check_(error, calc%topo%nangl,4)
   call check_(error, calc%topo%ntors,1)

   g = 0.0_wp
   gff_print=.true.

   call calc%singlepoint(env, mol, chk, 1, .false., etot, g, sigma, gap, res_gff)

   call check_(error, res_gff%e_total,-0.964158673358_wp, thr=thr)
   call check_(error, res_gff%gnorm,   0.013624271716_wp, thr=thr)
   call check_(error, res_gff%e_bond, -0.856707644691_wp, thr=thr)
   call check_(error, res_gff%e_angl,  0.000579711987_wp, thr=thr)
   call check_(error, res_gff%e_tors,  0.000000008811_wp, thr=thr)
   call check_(error, res_gff%e_es,   -0.150043166563_wp, thr=thr*10)
   call check_(error, res_gff%e_disp, -0.001251669186_wp, thr=thr)
   call check_(error, res_gff%e_rep,   0.066881023899_wp, thr=thr)
   call check_(error, res_gff%e_hb,   -0.006894292337_wp, thr=thr)
   call check_(error, res_gff%e_xb,   -0.000000000000_wp, thr=thr)
   call check_(error, res_gff%e_batm, -0.000000000000_wp, thr=thr)
   call check_(error, res_gff%g_solv, -0.016722645252_wp, thr=thr)
   call check_(error, res_gff%g_sasa,  0.000126368690_wp, thr=thr)
   call check_(error, res_gff%g_hb,   -0.009238121227_wp, thr=thr)
   call check_(error, res_gff%g_born, -0.009468335867_wp, thr=thr)
   call check_(error, res_gff%g_shift, 0.001857443126_wp, thr=thr)

   call mol%deallocate

end subroutine test_gfnff_gbsa


subroutine test_gfnff_mindless_basic(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)

   character(len=*), parameter :: mindless(10) = [&
      & "mindless01", "mindless02", "mindless03", "mindless04", "mindless05", &
      & "mindless06", "mindless07", "mindless08", "mindless09", "mindless10"]
   real(wp), parameter :: ref_energies(10) = &
      &[-1.6235670473950_wp, -1.297406090454654_wp, -1.582507291161514_wp, &
      & -1.6086171137100_wp, -1.509052868925628_wp, -1.607621995045164_wp, &
      & -1.732814270402046_wp, -1.887533992836867_wp, -1.392403540539921_wp, &
      & -1.958370274153000_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.117147122151611_wp, 0.089331678173547_wp, 0.156871336157090_wp, &
      & 0.09820449598566_wp, 0.084673071047178_wp, 0.087877392753885_wp, &
      & 0.124633459584315_wp, 0.10062733351355_wp, 0.06347505376937_wp, &
      & 0.094455606244636_wp]

   call init(env)
   do iMol = 1, 10

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call delete_file('charges')
      call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)

   end do

end subroutine test_gfnff_mindless_basic


subroutine test_gfnff_mindless_solvation(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)

   character(len=*), parameter :: mindless(10) = [&
      & "mindless01", "mindless02", "mindless03", "mindless04", "mindless05", &
      & "mindless06", "mindless07", "mindless08", "mindless09", "mindless10"]
   character(len=*), parameter :: solvents(10) = [character(len=20) ::&
      & "h2o", "chcl3", "thf", "acetonitrile", "toluene", &
      & "ch2cl2", "ether", "methanol", "cs2", "dmso"]
   real(wp), parameter :: ref_energies(10) = &
      &[-1.6540916244790_wp, -1.314305650588927_wp, -1.614189945831727_wp, &
      & -1.6204572122425_wp, -1.528526520178969_wp, -1.629308985567609_wp, &
      & -1.756188857803479_wp, -1.903412312838537_wp, -1.406801065436834_wp, &
      & -1.972260175936260_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.115941224167542_wp, 0.088717961962545_wp, 0.151908311668963_wp, &
      & 0.09873899177968_wp, 0.084338387828248_wp, 0.086708754851036_wp, &
      & 0.119474890540696_wp, 0.10022282792406_wp, 0.06288311963180_wp, &
      & 0.095645684316527_wp]

   call init(env)
   do iMol = 1, 10

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call delete_file('charges')
      call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
      call addSolvationModel(env, calc, TSolvInput(solvent=trim(solvents(iMol)), &
         & alpb=mod(iMol, 2)==0, kernel=gbKernel%still))

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)

   end do

end subroutine test_gfnff_mindless_solvation


subroutine test_gfnff_scaleup(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)

   character(len=*), parameter :: molecules(5) = [character(len=20) ::&
      & "caffeine", "rivaroxaban", "grubbs", "remdesivir", "taxol"]
   character(len=*), parameter :: solvents(5) = [character(len=20) ::&
      & "h2o", "acetonitrile", "toluene", "ether", "dmso"]
   real(wp), parameter :: ref_energies(5) = &
      &[-4.691966834534834_wp, -8.880608838859230_wp, -13.407098044870112_wp, &
      & -13.822447571445171_wp, -20.577627233774599_wp]
   real(wp), parameter :: ref_gnorms(5) = &
      &[0.059472120579392_wp, 0.095224192845175_wp, 0.175362338835436_wp, &
      & 0.124794782905474_wp, 0.193682754414573_wp]

   call init(env)
   do iMol = 1, 5

      call getMolecule(mol, trim(molecules(iMol)))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call delete_file('charges')
      call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
      call addSolvationModel(env, calc, TSolvInput(solvent=trim(solvents(iMol)), &
         & alpb=mod(iMol, 2)==0, kernel=gbKernel%still))

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)

   end do

end subroutine test_gfnff_scaleup


subroutine test_gfnff_pdb(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)
   real(wp), parameter :: ref_energies(3) = &
      &[-12.358374073082047_wp, -12.442021492870769_wp, -12.441563890957545_wp]
   real(wp), parameter :: ref_gnorms(3) = &
      &[0.158543908308563_wp, 0.151470765364003_wp, 0.151492720120697_wp]

   call init(env)
   do iMol = 1, 3

      call getMolecule(mol, 'pdb-4qxx')

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call delete_file('charges')
      call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
      if (iMol > 1) then
         call addSolvationModel(env, calc, TSolvInput(solvent='h2o', &
            & alpb=iMol==3, kernel=gbKernel%still))
      end if

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)

   end do

end subroutine test_gfnff_pdb


subroutine test_gfnff_sdf(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel

   use xtb_setparam, only : set

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)
   real(wp), parameter :: ref_energies(3) = &
      &[-0.983213277503977_wp, -1.082452170542686_wp, -1.082483911234450_wp]
   real(wp), parameter :: ref_gnorms(3) = &
      &[0.001248107348587_wp, 0.005694085165979_wp, 0.005681761874571_wp]

   call init(env)
   do iMol = 1, 3

      call getMolecule(mol, 'bug332')
      set%ichrg = nint(mol%chrg)

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call delete_file('charges')
      call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
      if (iMol > 1) then
         call addSolvationModel(env, calc, TSolvInput(solvent='h2o', &
            & alpb=iMol==3, kernel=gbKernel%p16))
      end if

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)

   end do

end subroutine test_gfnff_sdf


subroutine test_gfnff_pbc(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp
   real(wp), parameter :: thr2 = 3.0e-3_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

    integer :: iMol, i
    logical :: exitRun
    real(wp) :: energy, energy2, hl_gap, sigma(3, 3)
    real(wp), allocatable :: gradient(:, :)
    real(wp), allocatable :: xyztmp(:,:), lattmp(:,:)
    integer, allocatable :: attmp(:)
    character(len=4), allocatable :: symtmp(:)
    logical, parameter :: pbc(3) = [.true., .true., .true. ]
    ! structures from X23, mcVOL22, and GFN-FF for Ln/An paper
    character(len=*), parameter :: pbc_strucs(3) = [&
       & "x06_b", "mcv15", "Th_15"]
    !> references for original GFN-FF
    real(wp), parameter :: ref_energies(3) = &
        & [-9.522300465997_wp, -11.059826760233_wp, -46.265672821_wp ]
    real(wp), parameter :: ref_gnorms(3) = &
        & [ 0.083513307537938_wp, 0.083870464024606_wp, 1.6127503967_wp ]
    real(wp), parameter :: ref_snorm(2) = &
        & [ 1.410826798_wp, 0.419545170_wp ]
    !> references for mcGFN-FF
    real(wp), parameter :: mcref_energies(2) = &
        & [-9.419561538766_wp, -10.850342427712_wp ]
    real(wp), parameter :: mcref_gnorms(2) = &
        & [ 0.071130542937_wp, 0.104285409436030_wp ]
    real(wp), parameter :: mcref_snorm(2) = &
        & [ 0.839555069_wp, 0.729174966856260_wp ]

    call init(env)
    do iMol = 1, 2

       !  Calucaltions with original parameterization  !
       ! load molecule from molstock
       call getMolecule(mol, pbc_strucs(iMol))
       ! reset energy, gradient, sigma
       energy=0.0_wp
       if (allocated(gradient)) deallocate(gradient)
       allocate(gradient(3, mol%n))
       sigma = 0.0_wp
       ! setup new calculator
       call delete_file('charges')
       ! original angewChem2020_2 version is default
       call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
       call env%check(exitRun)
       call check_(error, .not.exitRun)
       if (exitRun) exit
       ! run single point calculation
       call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
          & hl_gap, res)
       call env%check(exitRun)
       call check_(error, .not.exitRun)
       if (exitRun) exit
       ! check energy, gradient and sigma versus reference calculation
       call check_(error, energy, ref_energies(iMol), thr=thr)
       call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)
       call check_(error, norm2(sigma), ref_snorm(iMol), thr=thr)

       !  Calucaltions with mcGFN-FF parameterization  !
       ! reset energy, gradient, sigma
       energy=0.0_wp
       if (allocated(gradient)) deallocate(gradient)
       allocate(gradient(3, mol%n))
       sigma = 0.0_wp
       ! setup new calculator
       call delete_file('charges')
       ! mcGFN-FF is version=4
       call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false., 4)
       call env%check(exitRun)
       call check_(error, .not.exitRun)
       if (exitRun) exit
       ! run single point calculation
       call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
          & hl_gap, res)
       call env%check(exitRun)
       call check_(error, .not.exitRun)
       if (exitRun) exit
       ! check energy, gradient and sigma versus reference calculation
       call check_(error, energy, mcref_energies(iMol), thr=thr)
       call check_(error, norm2(gradient), mcref_gnorms(iMol), thr=thr)
       call check_(error, norm2(sigma), mcref_snorm(iMol), thr=thr)
    end do

    ! check GFN-FF calculation on periodic system with actinide (Th)
    ! load molecule "Th_15" from molstock
    call getMolecule(mol, pbc_strucs(3))
    ! reset energy, gradient, sigma
    energy=0.0_wp
    if (allocated(gradient)) deallocate(gradient)
    allocate(gradient(3, mol%n))
    sigma = 0.0_wp
    ! setup new calculator
    call delete_file('charges')
    ! original angewChem2020_2 version is default
    call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! run single point calculation
    call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
       & hl_gap, res)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! check energy and gradient versus reference calculation
    call check_(error, energy, ref_energies(3), thr=thr)
    call check_(error, norm2(gradient), ref_gnorms(3), thr=thr)


    !  check super cell scaling  !
    ! load molecule from molstock
    call getMolecule(mol, "mcv15")
    ! reset energy, gradient, sigma
    energy=0.0_wp
    if (allocated(gradient)) deallocate(gradient)
    allocate(gradient(3, mol%n))
    sigma = 0.0_wp
    ! setup new calculator
    call delete_file('charges')
    ! mcGFN-FF is version=4
    call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false., 4)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! run single point calculation
    call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
       & hl_gap, res)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! build 2x1x1 supercell
    ! coordinates
    if (allocated(xyztmp)) deallocate(xyztmp)
    allocate(xyztmp(3,mol%n))
    xyztmp = mol%xyz
    deallocate(mol%xyz)
    allocate(mol%xyz(3,2*mol%n), source=0.0_wp)
    mol%xyz(:,1:mol%n)=xyztmp
    do i=1, mol%n
       mol%xyz(:,i+mol%n) = xyztmp(:,i)+mol%lattice(:,1)
    enddo
    ! atom types
    allocate(attmp(mol%n), source=0)
    attmp=mol%at
    deallocate(mol%at)
    allocate(mol%at(2*mol%n))
    mol%at(1:mol%n)=attmp
    mol%at(mol%n+1:2*mol%n)=attmp
    ! symbols
    allocate(symtmp(2*mol%n))
    symtmp = mol%sym
    deallocate(mol%sym)
    allocate(character(len=len(symtmp)) :: mol%sym(2*mol%n))
    mol%sym(1:mol%n)= symtmp
    mol%sym(mol%n+1:2*mol%n) = symtmp
    ! lattice
    mol%lattice(:,1) = 2.0_wp*mol%lattice(:,1)
    ! number of atoms
    mol%n = 2*mol%n
    ! init mol
    deallocate(symtmp, xyztmp)
    allocate(symtmp(mol%n))
    allocate(xyztmp(3,mol%n), source=0.0_wp)
    allocate(lattmp(3,3), source=0.0_wp)
    symtmp=mol%sym
    xyztmp=mol%xyz
    lattmp=mol%lattice
    call init(mol, symtmp, xyztmp, 0.0_wp, 0, lattmp, pbc)

    energy2=0.0_wp ! energy of supercell
    if (allocated(gradient)) deallocate(gradient)
    allocate(gradient(3, mol%n))
    sigma = 0.0_wp
    ! setup new calculator
    call delete_file('charges')
    ! mcGFN-FF is version=4
    call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false., 4)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! run single point calculation
    call calc%singlepoint(env, mol, chk, 2, .false., energy2, gradient, sigma, &
       & hl_gap, res)
    call env%check(exitRun)
    call check_(error, .not.exitRun)

    ! scale down energy of supercell and compare
    energy2 = energy2/2.0_wp
    call check_(error, energy2, energy, thr=thr2)

end subroutine test_gfnff_pbc

subroutine test_gfnff_fragmentation(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment, init
   use xtb_gfnff_fraghess, only : fragmentize

   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 501
   integer, parameter :: maxsystem = 600
   integer, parameter :: maxmagnat = 500
   integer, parameter :: numnb = 2
   integer, parameter :: numctr = 1

   type(TEnvironment) :: env
   integer :: at(nat), neigh(numnb,nat,numctr), counts(nat)
   integer :: i, j, nsystem, shifted_nsystem
   integer, allocatable :: fragments(:,:), fragment_sizes(:)
   integer, allocatable :: shifted_fragments(:,:), shifted_sizes(:)
   real(wp) :: xyz(3,nat), shifted_xyz(3,nat)
   real(wp), allocatable :: jab(:)

   call init(env)
   at = 1
   neigh = 0
   allocate(jab(nat*(nat+1)/2), source=1.0_wp)
   do i = 1, nat
      jab(i*(i+1)/2) = 0.0_wp
      j = i - 1
      xyz(:,i) = [real(mod(j,11),wp), real(mod(j/11,11),wp), &
         & real(j/121,wp)]
   end do

   call fragmentize(nat, at, xyz, maxsystem, maxmagnat, jab, numnb, numctr, neigh, &
      & fragments, fragment_sizes, nsystem, env)
   shifted_xyz = xyz + 300.0_wp
   call fragmentize(nat, at, shifted_xyz, maxsystem, maxmagnat, jab, numnb, numctr, neigh, &
      & shifted_fragments, shifted_sizes, shifted_nsystem, env)

   call check_(error, shifted_nsystem, nsystem)
   if (allocated(error)) return
   call check_(error, all(shifted_sizes(1:nsystem) == fragment_sizes(1:nsystem)))
   if (allocated(error)) return
   call check_(error, all(shifted_fragments(:,1:nsystem) == fragments(:,1:nsystem)))
   if (allocated(error)) return

   counts = 0
   do i = 1, nsystem
      do j = 1, fragment_sizes(i)
         if (fragments(j,i) < 1 .or. fragments(j,i) > nat) then
            call check_(error, .false.)
            return
         end if
         counts(fragments(j,i)) = counts(fragments(j,i)) + 1
      end do
   end do
   call check_(error, all(counts == 1))
   if (allocated(error)) return
   call check_(error, maxval(fragment_sizes(1:nsystem)) <= maxmagnat)
   if (allocated(error)) return

   ! A crowded grid cell must be divided without losing disconnected atoms.
   xyz = 300.0_wp
   call fragmentize(nat, at, xyz, maxsystem, maxmagnat, jab, numnb, numctr, neigh, &
      & fragments, fragment_sizes, nsystem, env)
   call check_(error, nsystem, 2)
   if (allocated(error)) return
   call check_(error, sum(fragment_sizes(1:nsystem)), nat)
   if (allocated(error)) return
   call check_(error, maxval(fragment_sizes(1:nsystem)) <= maxmagnat)

end subroutine test_gfnff_fragmentation

subroutine test_gfnff_LnAn_H(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp ! threshold for energy
   real(wp), parameter :: qthr = 1.0e-4_wp ! threshold for charge

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

    integer :: i
    integer, parameter :: nat = 33
    logical :: exitRun
    real(wp) :: energy, charges(nat), hl_gap, sigma(3,3)
    real(wp), allocatable :: gradient(:, :)
    real(wp), allocatable :: xyztmp(:,:), lattmp(:,:)
    integer, allocatable :: attmp(:)
    character(len=4), allocatable :: symtmp(:)
    logical, parameter :: pbc(3) = [.true., .true., .true. ]
    ! structures from X23, mcVOL22, and GFN-FF for Ln/An paper
    character(len=*), parameter :: struc(1) = ["Ce_0a7745"]
    !> references for original GFN-FF
    real(wp), parameter :: ref_charges(nat) = [&
      &  0.93746177_wp, -0.40179759_wp, -0.00012766_wp, -0.01777159_wp, -0.26280860_wp, &
      & -0.01740112_wp, -0.00287629_wp, -0.36245214_wp,  0.19642020_wp,  0.20665288_wp, &
      &  0.09972845_wp,  0.08148308_wp,  0.10674954_wp,  0.09486206_wp,  0.17296671_wp, &
      &  0.11414769_wp,  0.09722964_wp,  0.10512679_wp,  0.08798849_wp,  0.22945410_wp, &
      &  0.15310280_wp,  0.88629154_wp, -0.48033324_wp, -0.48478265_wp, -0.49522170_wp, &
      & -0.14264511_wp, -0.12522326_wp, -0.01061975_wp, -0.40598965_wp,  0.13119126_wp, &
      &  0.11822500_wp,  0.13303034_wp,  0.25781424_wp &
      & ]
    real(wp), parameter :: ref_energy = -3.550569267815145_wp ! reference energy
    real(wp), parameter :: ref_gnorm = 0.213640754560263_wp ! reference gradient norm
    ! number of neighbors and sum of neighbor indices for atom 1 (Ce)
    integer, parameter :: ref_numnb_1 = 9
    integer, parameter :: ref_sumnb_1 = 177 ! 2+5+8+22+24+25+26+27+29 +9
    ! number of neighbors and sum of neighbor indices for atom 21 (H)
    integer, parameter :: ref_numnb_21 = 1
    integer, parameter :: ref_sumnb_21 = 9 ! 8 +1
    ! number of neighbors and sum of neighbor indices for atom 33 (H)
    integer, parameter :: ref_numnb_33 = 1
    integer, parameter :: ref_sumnb_33 = 30 ! 29 +1



    call init(env)

    ! check GFN-FF calculation on periodic system with actinide (Th)
    ! load molecule from molstock
    call getMolecule(mol, struc(1))
    ! reset energy and gradient
    energy=0.0_wp
    if (allocated(gradient)) deallocate(gradient)
    allocate(gradient(3, mol%n))
    ! setup new calculator
    call delete_file('charges')
    ! original angewChem2020_2 version is default
    call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! run single point calculation
    call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
       & hl_gap, res)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! check energy and gradient versus reference calculation
    call check_(error, energy, ref_energy, thr=thr)
    call check_(error, norm2(gradient), ref_gnorm, thr=thr)
    ! get charges from SP calculation
    charges = chk%nlist%q ! these charges are in the gfnff_charges file
    ! check charges of all atoms
    do i=1, nat
      ! compare to reference
      call check_(error, charges(i), ref_charges(i), thr=qthr)
    enddo
    ! The Ln/An-H bonds in GFN-FF are charge dependent
    ! But it is compared in the gfnff_ini where only topo%qa is available
    ! Therefore we check the neighbor list here additionally
    call check_(error, calc%neigh%nb(42,1,1), ref_numnb_1) ! check numnb of atom 1
    call check_(error, sum(calc%neigh%nb(:,1,1)), ref_sumnb_1) ! check sumnb of atom 1
    call check_(error, calc%neigh%nb(42,21,1), ref_numnb_21) ! check numnb of atom 21
    call check_(error, sum(calc%neigh%nb(:,21,1)), ref_sumnb_21) ! check sumnb of atom 21
    call check_(error, calc%neigh%nb(42,33,1), ref_numnb_33) ! check numnb of atom 33
    call check_(error, sum(calc%neigh%nb(:,33,1)), ref_sumnb_33) ! check sumnb of atom 33

end subroutine test_gfnff_LnAn_H

end module test_gfnff
