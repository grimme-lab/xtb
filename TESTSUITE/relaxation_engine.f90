subroutine test_pbc_lancopt
   use iso_fortran_env, wp => real64
   use assertion

   use tbdef_molecule

   use model_phonons

   implicit none

   real(wp),parameter :: thr = 1.0e-9_wp
   integer, parameter :: nat = 6
   integer, parameter :: at(nat) = [22,22, 8, 8, 8, 8]
   real(wp),parameter :: xyz(3,nat) = reshape( &
      &[0.00000000_wp, 0.00000000_wp, 0.00000000_wp, &
      & 0.08499567_wp, 0.08499567_wp, 0.05476752_wp, &
      & 0.04795456_wp, 0.04795456_wp, 0.00000000_wp, &
      & 0.03704111_wp, 0.13295023_wp, 0.05476752_wp, &
      & 0.12203678_wp, 0.12203678_wp, 0.00000000_wp, &
      & 0.13295023_wp, 0.03704111_wp, 0.05476752_wp], shape(xyz))
   real(wp),parameter :: lattice(3,3) = reshape( &
      &[0.16999134_wp, 0.00000000_wp, 0.00000000_wp, &
      & 0.00000000_wp, 0.16999134_wp, 0.00000000_wp, &
      & 0.00000000_wp, 0.00000000_wp, 0.10953503_wp], shape(lattice))
   integer, parameter :: wsc_rep(3) = [1,1,1]

   type(tb_molecule) :: mol
   type(mp_options)  :: opt = mp_options(k_stretch = 0.04_wp, coupled = .true.)

   integer :: i,j,ij
   integer :: nvar,npvar
   real(wp), allocatable :: hess(:)
   real(wp), allocatable :: phon(:,:)
   real(wp), allocatable :: freq(:)

   integer  :: info
   integer  :: lwork
   integer  :: liwork
   integer, allocatable :: iwork(:)
   real(wp),allocatable :: aux(:)

   call mol%allocate(nat)
   mol%at   = at
   mol%xyz  = xyz
   mol%npbc = 3
   mol%pbc  = .true.
   mol%lattice = lattice
   call mol%update
   call generate_wsc(mol,mol%wsc,wsc_rep)
   call generate_wsl(mol,mol%wsl,wsc_rep)


   call terminate(afail+1)

end subroutine test_pbc_lancopt
