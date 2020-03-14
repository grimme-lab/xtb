subroutine test_wigner_seitz_0d
   use xtb_mctc_accuracy, only : wp
   use assertion
   use xtb_type_molecule
   use xtb_type_wsc
   implicit none
   type(TMolecule) :: mol
   type(tb_wsc) :: wsc

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.00000000000000_wp,  0.00000000000000_wp, -0.73578586109551_wp, &
      & 1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp, &
      &-1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp  &
      & ],shape(xyz))

   call init(mol, at, xyz)

   call generate_wsc(mol,wsc)

   call assert(allocated(wsc%lattr))
   call assert(allocated(wsc%itbl))

   call assert(all(wsc%itbl <= 1))
   call assert(all(wsc%lattr == 0))

   call terminate(afail)
end subroutine test_wigner_seitz_0D

subroutine test_wigner_seitz_3D
   use xtb_mctc_accuracy, only : wp
   use assertion

   use xtb_type_molecule
   use xtb_type_wsc

   use xtb_pbc_tools

   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
   ! CaF2
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [9,9,20]
   real(wp),parameter :: abc(3,nat) = reshape(&
      &[0.25_wp, 0.25_wp, 0.25_wp, &
      & 0.75_wp, 0.75_wp, 0.75_wp, &
      & 0.00_wp, 0.00_wp, 0.00_wp], shape(abc))
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[5.9598811567890_wp,      2.1071361905157_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      6.3214085715472_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      0.0000000000000_wp,      7.2993338808807_wp],   &
      & shape(lattice))

   type(TMolecule) :: mol
   type(tb_wsc) :: wsc
   real(wp), allocatable :: xyz(:, :)

   allocate(xyz(3, nat))
   call coord_trafo(nat,lattice,abc,xyz)
   call init(mol, at, xyz, lattice=lattice)

   call generate_wsc(mol,wsc)

   call assert(allocated(wsc%lattr))
   call assert(allocated(wsc%itbl))

   print*,wsc%itbl
   call assert(all(abs(wsc%lattr) <= 1))
   call assert_eq(wsc%itbl(1,1), 12)
   call assert_eq(wsc%itbl(2,2), 12)
   call assert_eq(wsc%itbl(3,3), 12)
   call assert_eq(wsc%itbl(1,2), 6)
   call assert_eq(wsc%itbl(1,3), 4)
   call assert_eq(wsc%itbl(2,3), 4)
   call assert_eq(wsc%itbl(1,2), wsc%itbl(2,1))
   call assert_eq(wsc%itbl(1,3), wsc%itbl(3,1))
   call assert_eq(wsc%itbl(2,3), wsc%itbl(3,2))

   call terminate(afail)
end subroutine test_wigner_seitz_3D
