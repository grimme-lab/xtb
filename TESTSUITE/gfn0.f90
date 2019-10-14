subroutine test_gfn0_sp
   use iso_fortran_env, wp => real64

   use assertion

   use mctc_systools
   use tbdef_options

   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data

   use setparam
   use aoparam
   use xbasis
   use peeq_module

   implicit none
   real(wp),parameter :: thr = 1.0e-7_wp
   real(wp),parameter :: thr2 = 1.0e-5_wp
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
   real(wp),parameter :: et = 300.0_wp
   integer, parameter :: maxiter = 50
   integer, parameter :: prlevel = 2
   logical, parameter :: lgrad = .true.
   logical, parameter :: restart = .false.
   real(wp),parameter :: acc = 1.0_wp

   type(gfn_parameter)   :: gfn
   type(tb_environment)  :: env
   type(tb_molecule)     :: mol
   type(scc_results)     :: res
   type(tb_basisset)     :: basis
   type(tb_wavefunction) :: wfn
   type(scc_parameter)   :: param

   real(wp) :: etot,egap,sigma(3,3)
   real(wp), allocatable :: g(:,:)
   character(len=:),allocatable :: fnv
   integer  :: ipar

   real(wp) :: globpar(25)
   logical  :: okpar,okbas,exist,diff

   call env%setup

   gfn_method = 0

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp
   call mol%set_nuclear_charge
   call mol%update

   wfn%nel = idint(sum(mol%z))
   wfn%nopen = 0

   allocate( g(3,mol%n), source = 0.0_wp )
 
   call use_parameterset('.param_gfn0.xtb',globpar,okpar)
   !call assert(okpar)
      call rdpath(env%xtbpath,'.param_gfn0.xtb',fnv,exist)
      ! maybe the user provides a local parameter file, this was always
      ! an option in `xtb', so we will give it a try
      if (.not.exist) fnv = '.param_gfn0.xtb'
      call open_file(ipar,fnv,'r')
      if (ipar.eq.-1) then
         ! at this point there is no chance to recover from this error
         ! THEREFORE, we have to kill the program
         call raise('E',"Parameter file '"//fnv//"' not found!")
      endif
      call read_gfn_param(ipar,globpar,.true.)
      call close_file(ipar)

   call set_gfn0_parameter(param,globpar,mol%n,mol%at)

   call xbasis0(mol%n,mol%at,basis)
   call xbasis_gfn0(mol%n,mol%at,basis,okbas,diff)
   call assert(okbas)

   call xbasis_cao2sao(mol%n,mol%at,basis)

   call assert_eq(basis%nshell,17)
   call assert_eq(basis%nbf,   30)
   call assert_eq(basis%nao,   29)

   call wfn%allocate(mol%n,basis%nshell,basis%nao)
   wfn%q = mol%chrg/real(mol%n,wp)

   g = 0.0_wp

   call peeq(output_unit,mol,wfn,basis,param, &
      &   egap,et,prlevel,lgrad,.false.,acc,etot,g,sigma,res)

   call assert(res%converged)

   call assert_close(res%e_total,-15.271927490420_wp,thr)
   call assert_close(res%gnorm,0.55600868499186E-01_wp,thr)
   ! value in electron volt
   call assert_close(res%hl_gap, 4.4910991783546_wp,1.0e-4_wp)
   call assert_close(res%e_elec,-15.238855850629_wp,thr)
   ! es and aes are not welldefined at this accuracy level
   call assert_close(res%e_es,   -0.091088841288_wp,thr*10)
   call assert_close(res%e_disp, -0.001331408911_wp,thr)
   call assert_close(res%e_rep,   0.071830283799_wp,thr)
   call assert_close(res%e_xb,   -0.012481673391_wp,thr)

   call mol%deallocate
   call wfn%deallocate
   call basis%deallocate

   call terminate(afail)
end subroutine test_gfn0_sp

subroutine test_gfn0_api
   use iso_fortran_env, wp => real64, istdout => output_unit
   use assertion

   use tbdef_options
   use tbdef_molecule
   use tbdef_param

   use pbc_tools

   use tb_calculators

   implicit none

   real(wp),parameter :: thr = 1.0e-7_wp
   integer, parameter :: nat = 7
   integer, parameter :: at(nat) = [51,1,1,1,1,16,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [-0.33839323285714_wp,-2.52035427985714_wp, 0.00000000000000_wp, &
      & 0.87126950614286_wp,-4.36564525285714_wp, 2.28275498100000_wp, &
      & 0.87126950614286_wp,-4.36564525285714_wp,-2.28275498100000_wp, &
      &-3.07540243785714_wp,-4.12867550385714_wp, 0.00000000000000_wp, &
      & 2.16321847414286_wp, 6.15526246914286_wp, 0.00000000000000_wp, &
      &-0.33839323285714_wp, 5.87053934414286_wp, 0.00000000000000_wp, &
      &-0.15356858285714_wp, 3.35451847614286_wp, 0.00000000000000_wp], shape(xyz))
   type(peeq_options),parameter :: opt = peeq_options( &
      &  prlevel = 2, ccm = .false., acc = 1.0_wp, etemp = 300.0_wp, grad = .true. )

   type(tb_molecule)    :: mol
   type(tb_environment) :: env
   type(gfn_parameter)  :: gfn

   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp) :: dum(3,3)
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call env%setup

   call mol%allocate(nat)
   mol%at   = at
   mol%xyz  = xyz
   call mol%set_nuclear_charge
   call mol%update

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call gfn0_calculation &
      (istdout,env,opt,mol,gfn,hl_gap,energy,gradient,dum,dum)

   call assert_close(hl_gap, 5.5384029314207_wp,thr)
   call assert_close(energy,-8.6908532561691_wp,thr)
   call assert_close(norm2(gradient),0.24709483673929E-01_wp,thr)

   call assert_close(gradient(1,3), 0.12460501676405E-02_wp,thr)
   call assert_close(gradient(2,1),-0.20307939298595E-01_wp,thr)
   call assert_close(gradient(1,6),-0.29075269403059E-02_wp,thr)
   call assert_close(gradient(3,5), 0.00000000000000E+00_wp,thr)

   call terminate(afail)

end subroutine test_gfn0_api
