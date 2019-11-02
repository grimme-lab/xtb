module tbdef_hamiltonian
   use iso_fortran_env, only: wp => real64
   use mctc_param, only : pauling_en, chemical_hardness
   implicit none
   public :: tb_hamiltonian
   private


   integer, parameter :: max_elem = 94
   integer, parameter :: max_shell = 3
   integer, parameter :: max_prim = 12

   type :: tb_hamiltonian
      integer :: method = -1
      real(wp) :: en(max_elem) = pauling_en(:max_elem)
      real(wp) :: mc(max_elem) = 0.0_wp
      real(wp) :: gam(max_elem) = chemical_hardness(:max_elem)
      real(wp) :: gam3(max_elem) = 0.0_wp
      real(wp) :: alp0(max_elem) = 0.0_wp
      real(wp) :: wll(max_elem,10) = 0.0_wp
      real(wp) :: rep(2,max_elem)  = 0.0_wp
      real(wp) :: polyr(4,max_elem) = 0.0_wp
      real(wp) :: cxb(max_elem) = 0.0_wp
      real(wp) :: ao_exp(10,max_elem) = 0.0_wp
      real(wp) :: ao_lev(10,max_elem) = 0.0_wp
      real(wp) :: lpar(0:2,max_elem) = 0.0_wp
      real(wp) :: kpair(max_elem,max_elem) = 1.0_wp
      real(wp) :: kcnat(0:2,max_elem) = 0.0_wp
      real(wp) :: kqat(3,max_elem) = 0.0_wp
      real(wp) :: radaes(max_elem) = 5.0_wp
      real(wp) :: dpolc(max_elem) = 0.0_wp
      real(wp) :: qpolc(max_elem) = 0.0_wp
      integer  :: ao_pqn(10,max_elem) = 0
      integer  :: ao_l(10,max_elem) = 0
      integer  :: ao_n(max_elem) = 0
      integer  :: ao_typ(10,max_elem) = 0
      integer  :: cnval(max_elem) = [&
         & 1,                                                             1, &
         & 1, 2,                                           3, 3, 3, 2, 1, 1, &
         & 1, 2,                                           3, 3, 3, 3, 1, 1, &
         & 1, 2, 4,          4, 6, 6, 6, 6, 6, 4, 4,  2,   3, 3, 3, 3, 1, 1, &
         & 1, 2, 4,          4, 6, 6, 6, 6, 6, 4, 4,  2,   3, 3, 3, 3, 1, 1, &
         & 1, 2, 4, &
         &       6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
         &                   4, 6, 6, 6, 6, 6, 4, 4,  2,   3, 3, 3, 3, 1, 1, &
         & 0, 0, 0, 0, 0, 0, 0, 0 ]
      character(len=30) :: timestp(max_elem) = '------------------------------'
   contains
   end type tb_hamiltonian


contains


end module tbdef_hamiltonian
