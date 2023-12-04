! This file is part of ancopt.
! SPDX-Identifier: LGPL-3.0-or-later
!
! ancopt is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! ancopt is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with ancopt.  If not, see <https://www.gnu.org/licenses/>.

module xtb_pbc_optimizer_driver
   use xtb_mctc_accuracy, only : wp!, error_type, fatal_error
   use xtb_type_molecule, only : TMolecule
   use xtb_type_calculator
   use xtb_gfnff_calculator, only : TGFFCalculator
   use xtb_type_environment, only : TEnvironment
   use xtb_pbc_optimizer_filter_cart, only : cartesian_filter
   use xtb_type_restart, only : TRestart
   use xtb_type_data, only : scc_results
   use xtb_pbc_optimizer_lbfgs, only : optimizer_type, lbfgs_optimizer

   implicit none
   private

    public :: relax_pbc


 contains


!> Driver for performing geometry optimization
!subroutine    relax(self, ctx, mol, calc, filter, accuracy, verbosity)
subroutine relax_pbc(self, env, mol, chk, calc, filter, printlevel)
  use xtb_pbc, only: cross_product
  use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
  !> Instance of the optimization driver
  class(optimizer_type), intent(inout) :: self
  !> Calculation environment
  type(TEnvironment), intent(inout) :: env
  type(TEnvironment) :: env2
  !> Molecular structure data
  type(TMolecule), intent(inout) :: mol
  !> Wavefunction data
  type(TRestart), intent(inout) :: chk
  !> calculator type
  type(TGFFCalculator), allocatable :: gfnff
  class(TCalculator), intent(inout) :: calc
  !> Transformation filter for coordinates
  class(cartesian_filter), intent(in) :: filter
  !> Output verbosity
  integer, intent(in) :: printlevel
  ! Variables for/from singlepoint 
  !> Restart from previous results
  logical :: restart ! = .true.
  !> HOMO-LUMO gab
  real(wp) :: hlgap
  !> Detailed results
  type(scc_results) :: results
  !> optimization cycle failed
  logical :: fail
  !> singlepoint energy
  real(wp) :: energy
  !> gradient and stress tensor
  real(wp), allocatable :: gxyz(:, :), sigma(:, :)

  !> Convergence accuracy
  real(wp) :: accuracy
  integer :: nvar, step, this_step, i, stat
  real(wp) :: elast, ediff, gnorm, dnorm, ethr, gthr, dthr, gamax, damax, cvol
  real(wp), allocatable :: val, gcurr(:), glast(:), displ(:)
  logical :: econv, gconv, dconv, converged, f_exists
  !> Source of errors in the program unit
  character(len=*), parameter :: source = "relax_pbc"
  !> numerical gradient and stress tensor for debugging
  real(wp) :: numgrad(3,mol%n), numsig(3,3)
  logical, parameter :: debug = .false.

  real(wp) :: emin_global, xyz_emin(3,mol%n), latt_emin(3,3)
  character(6), parameter :: fname = "noName"
  
  nvar = filter%get_dimension()
  allocate(gcurr(nvar), glast(nvar), displ(nvar))
  allocate(gxyz(3, mol%n), sigma(3, 3), source=0.0_wp)
  val = huge(0.0_wp)
  emin_global = huge(0.0_wp)
  energy = 0.0_wp
  gcurr(:) = 0.0_wp
  glast(:) = 0.0_wp
  displ(:) = 0.0_wp

  accuracy = 1.0_wp
  call get_thr(accuracy, ethr, gthr, dthr)
  this_step = 20*nvar

  inquire( file="xtboptlog.cif", exist=f_exists )
  if(f_exists)then
    open(unit=133, iostat=stat, file="xtboptlog.cif", status='old')
    if (stat == 0) close(133, status='delete')
  endif
  open(133, file = "xtboptlog.cif", status = "new")

  ! start of optimization loop
  optLoop: do step = 1, this_step
    ! CYCLE header
    if (printlevel > 0) then
      write(env%unit,'(/,72("."),/,30(".")," CYCLE",i5,1x,30("."),/,72("."))')step
    end if

    ! 
    call filter%transform_structure(mol, displ) ! original

    dnorm = norm2(displ)
    damax = maxval(abs(gcurr))
    dconv = dnorm/nvar < dthr

    glast(:) = gcurr
    elast = val
    restart=.true.

    call calc%singlepoint(env, mol, chk, printlevel-1, restart, &
         & energy, gxyz, sigma, hlgap, results)

    call env%check(fail)
    if (fail) then
      call env%error('Single point calculation failed, aborting...', source)
      return
    endif

    ! calculate numerical gradient and stress tensor if requested
    if (debug) then
       numsig = 0.0_wp
       call num_sigma(env, mol, chk, calc, printlevel, sigma, numsig)
       numgrad = 0.0_wp
       call num_grad(env, mol, chk, calc, printlevel, gxyz, numgrad)
    endif


    ! calculate new gradients from gxyz and sigma from SP
    call filter%transform_derivative(mol, energy, gxyz, sigma, val, gcurr)

    ! set up values for convergence test
    ediff = val - elast
    econv = ediff <= epsilon(0.0_wp) .and. abs(ediff) < ethr
    gnorm = norm2(gcurr)
    gamax = maxval(abs(gcurr))
    gconv = gnorm/nvar < gthr
    converged = econv .and. gconv .and. dconv
    cvol=dot_product(mol%lattice(:,1), cross_product( &
                   & mol%lattice(:,2), mol%lattice(:,3)))/6.748334606_wp !in a0 now

    ! output key values of optimization step
    call step_summary(env, energy, ediff, gnorm, gamax, dnorm, damax, &
            cvol, mol%lattice, printlevel)

    ! save coordinates and lattice of  
    if (val.lt.emin_global) then
      emin_global = val
      xyz_emin  = mol%xyz
      latt_emin = mol%lattice
    endif
    ! write current structure to xtboptlog.cif
    call write_optlog(step, mol%n, mol%xyz, mol%lattice, mol%sym, val)

    ! check for convergence 
    if (converged) exit

    ! perform LBFGS step
    call self%step(val, gcurr, glast, displ, step)

  end do optLoop
  ! end of optimization loop

close(133) ! close xtboptlog.cif unit

! write structure of minimum energy to file
if (val.ne.emin_global) then
  inquire( file="xtbopt_emin.coord", exist=f_exists )
  if(f_exists)then
    open(unit=133, iostat=stat, file="xtbopt_emin.coord", status='old')
    if (stat == 0) then
      close(133, status='delete')
    else
      write(env%unit,*) 'A problem occured while handling xtbopt_emin.coord'
    endif
  endif
  open(133, file = "xtbopt_emin.coord", status = "new")
  write(133,*) '$coord'
  do i=1, mol%n
    write(133,'(3f20.10,a,a)') xyz_emin(:,i),'    ',mol%sym(i)
  enddo
  write(133,*) '$periodic 3'
  write(133,*) '$lattice'
  write(133,'(3f25.10)') latt_emin
  write(133,*) '$end'
  close(133)
endif

  ! Output message if not converged after maximum steps
  if (.not.converged) then
    write(*,*) ''
    write(*,*) 'Could not converge in',step,' steps.'
    call env%error("Could not converge geometry.",source)
  end if
  ! Output message if converged in time
  if (printlevel > 0 .and. converged) then
    write(*,*) ''
    write(env%unit,'(a,i6,a)') "Geometry optimization converged after"&
            &,step," steps."
  end if

end subroutine relax_pbc



subroutine step_summary(env, energy, ediff, gnorm, gamax, dnorm, damax, cvol, &
                & lattice, prlevel)
  !> Computational environment
  type(TEnvironment), intent(inout) :: env
  real(wp), intent(in) :: energy, ediff, gnorm, gamax, dnorm, damax, cvol
  real(wp), intent(in) :: lattice(3,3)
  integer, intent(in) :: prlevel

   if (prlevel <= 0) return

   if (prlevel > 1) then
     write(env%unit,'(a25,F25.15,a)')  " total energy   ",energy," Eh"
     write(env%unit,'(a25,E16.2,a)')    "   energy change",ediff," Eh"
     write(env%unit,*) ""
     write(env%unit,'(a25,F25.15,a)')  " gradient norm  ",gnorm," Eh/a0"
     write(env%unit,'(a25,F25.15,a)')    "   max. gradient",gamax," Eh/a0"
     write(env%unit,*) ""
     write(env%unit,'(a25,F25.15,a)')    " step length  ",dnorm," a0"
     write(env%unit,'(a25,F25.15,a)')  "   max. step    ",damax," a0"
     write(env%unit,*) ""
     write(env%unit,'(a25,F25.15,a)')    " unit cell volume  ",cvol," A^3"
     write(env%unit,*) ""
     write(env%unit,'(a25)')    " current lattice:  "
     write(env%unit,'(3F25.15)') lattice
   else
     write(env%unit,'(a25,F25.15,a)')  " total energy   ",energy," Eh"
   end if
end subroutine step_summary


subroutine get_thr(accuracy, ethr, gthr, dthr)
   real(wp), intent(in) :: accuracy
   real(wp), intent(out) :: ethr
   real(wp), intent(out) :: gthr
   real(wp), intent(out) :: dthr

   ethr = accuracy * 5.0e-7_wp
   gthr = accuracy * 1.0e-3_wp
   dthr = 5.0e-2_wp
end subroutine get_thr

!>
subroutine write_optlog(step, n, xyz, lattice, symbol, energy)
   integer, intent(in)  :: step
   integer, intent(in)  :: n
   real(wp), intent(in) :: xyz(3,n)
   real(wp), intent(in) :: lattice(3,3)
   character(*), intent(in) :: symbol(:)
   real(wp), intent(in) :: energy
   integer :: i
   real(wp) :: lat_len(3)
   real(wp) :: lat_ang(3)
   real(wp), parameter :: faua = 0.5291772109_wp

   call calculate_cell_parameters(lattice, lat_len, lat_ang)
   
   write(133,'(a,i0)') 'data_',step
   write(133,'(a,i6,a,f16.10)') '# Step:',step,'  Energy:',energy
   write(133,'(a)') "_symmetry_space_group_name_H-M   'P 1'"
   write(133,'(a,f12.4)') "_cell_length_a",lat_len(1)*faua
   write(133,'(a,f12.4)') "_cell_length_b",lat_len(2)*faua
   write(133,'(a,f12.4)') "_cell_length_c",lat_len(3)*faua
   write(133,'(a,f10.2)') "_cell_angle_alpha",lat_ang(1)
   write(133,'(a,f10.2)') "_cell_angle_beta",lat_ang(2)
   write(133,'(a,f10.2)') "_cell_angle_gamma",lat_ang(3)
   write(133,'(a)') "loop_"
   write(133,'(a)') "_atom_site_label"
   write(133,'(a)') "_atom_site_Cartn_x"
   write(133,'(a)') "_atom_site_Cartn_y"
   write(133,'(a)') "_atom_site_Cartn_z"
   do i=1, n
     write(133,'(A4,3f20.10)') symbol(i), xyz(1,i)*faua, xyz(2,i)*faua, xyz(3,i)*faua
   enddo
   write(133,*) ""
   
end subroutine write_optlog

!> Calculate the cell parameters for CIF format optimization trajectory
subroutine calculate_cell_parameters(lattice, lat_len, lat_ang)
  real(wp), intent(in) :: lattice(3,3)
  real(wp), intent(out) :: lat_len(3)
  real(wp), intent(out) :: lat_ang(3)

  ! Calculate lengths
  lat_len(1) = sqrt(lattice(1,1)**2 + lattice(2,1)**2 + lattice(3,1)**2) ! a
  lat_len(2) = sqrt(lattice(1,2)**2 + lattice(2,2)**2 + lattice(3,2)**2) ! b
  lat_len(3) = sqrt(lattice(1,3)**2 + lattice(2,3)**2 + lattice(3,3)**2) ! c

  ! Calculate angles (in degrees)
  lat_ang(1) = acosd(dot_product(lattice(:,2),lattice(:,3))/(lat_len(2)*lat_len(3))) ! angle "oposite" of a
  lat_ang(2) = acosd(dot_product(lattice(:,1),lattice(:,3))/(lat_len(1)*lat_len(3))) ! angle "oposite" of b
  lat_ang(3) = acosd(dot_product(lattice(:,1),lattice(:,2))/(lat_len(1)*lat_len(2))) ! angle "oposite" of c

end subroutine calculate_cell_parameters

! Function to convert radians to degrees degrees = rad*180/pi
pure function acosd(x) result(degrees)
  real(8), intent(in) :: x
  real(8) :: degrees
  degrees = acos(x) * 180.0_wp / acos(-1.0_wp)
end function acosd


!> numerical gradient for debugging purposes
subroutine num_grad(env, mol, chk, calc, printlevel, grad, numgrad)
  !> Calculation environment
  type(TEnvironment), intent(inout) :: env
  !> Molecular structure data
  type(TMolecule), intent(inout) :: mol
  !> Wavefunction data
  type(TRestart), intent(inout) :: chk
  !> calculator type
  class(TCalculator), intent(inout) :: calc
  !> Output verbosity
  integer, intent(in) :: printlevel
  !> gradient and sigma from previous SP calculation
  real(wp), intent(in) :: grad(3,mol%n)
  real(wp), intent(inout) :: numgrad(3,mol%n)
  ! Variables for/from singlepoint 
  !> Restart from previous results
  logical :: restart ! = .true.
  !> HOMO-LUMO gab
  real(wp) :: hlgap
  !> Detailed results
  type(scc_results) :: results
  !> singlepoint energy
  real(wp) :: energy
  !> gradient and stress tensor
  real(wp) :: sdum(3,3), gdum(3, mol%n), numg(3,mol%n)
  !> copy of original coordinates
  real(wp) :: coord(3,mol%n)
  !> copy of Wavefunction data
  type(TRestart) :: wf0
  !
  real(wp) :: nums(3,3), eps(3,3), er, el
  !
  integer :: i,j,k
  ! step size for numerical sigma
  real(wp),parameter    ::  step = 0.00001_wp, step2 = 0.5_wp/step ! for numerical gradient
  
  numgrad=0.0_wp
  gdum=0.0_wp
  sdum=0.0_wp
  coord = mol%xyz 
!  allocate( numg(3,mol%n),gdum(3,mol%n),sdum(3,3), source = 0.0_wp )
  wf0 = chk
  do i = 1, mol%n
     do j = 1, 3
        mol%xyz(j,i) = mol%xyz(j,i) + step
        chk = wf0
            call calc%singlepoint(env, mol, chk, printlevel-1, restart, &
               & er, gdum, sdum, hlgap, results)
!        call singlepoint &
!           &       (env,mol,chk,calc, &
!           &        egap,etemp,maxscciter,0,.true.,.true.,acc,er,gdum,sdum,res)
        mol%xyz(j,i) = mol%xyz(j,i) - 2*step
        chk = wf0
            call calc%singlepoint(env, mol, chk, printlevel-1, restart, &
               & el, gdum, sdum, hlgap, results)
!        call singlepoint &
!           &       (env,mol,chk,calc, &
!           &        egap,etemp,maxscciter,0,.true.,.true.,acc,el,gdum,sdum,res)
        mol%xyz(j,i) = mol%xyz(j,i) + step
        numg(j,i) = step2 * (er - el)
     enddo
  enddo
  print'(/,"analytical gradient")'
  print *, grad
  print'(/,"numerical gradient")'
  print *, numg
  print'(/,"difference gradient")'
  print*,grad-numg
  numgrad=numg

end subroutine num_grad

!> numerical stress tensor for debugging purposes
subroutine num_sigma(env, mol, chk, calc, printlevel, sigma, numsig)
  !> Calculation environment
  type(TEnvironment), intent(inout) :: env
  !> Molecular structure data
  type(TMolecule), intent(inout) :: mol
  !> Wavefunction data
  type(TRestart), intent(inout) :: chk
  !> calculator type
  class(TCalculator), intent(inout) :: calc
  !> Output verbosity
  integer, intent(in) :: printlevel
  !> gradient and sigma from previous SP calculation
  real(wp), intent(in) :: sigma(3,3)
  real(wp), intent(inout) :: numsig(3,3)
  ! Variables for/from singlepoint 
  !> Restart from previous results
  logical :: restart ! = .true.
  !> HOMO-LUMO gab
  real(wp) :: hlgap
  !> Detailed results
  type(scc_results) :: results
  !> singlepoint energy
  real(wp) :: energy
  !> gradient and stress tensor
  real(wp), allocatable :: gdum(:, :)
  ! dummies and copies
  real(wp) ::  sdum(3,3), latt(3,3)
  !> copy of original coordinates
  real(wp), allocatable :: coord(:,:)
  !> copy of Wavefunction data
  type(TRestart) :: wf0
  !
  real(wp) :: nums(3,3), eps(3,3), er, el
  !
  integer :: i,j,k
  ! step size for numerical sigma
  real(wp),parameter   :: sstep = 1.0_wp*10.0_wp**(-6), sstep2 = 0.5_wp/sstep
  

   ! ------------------------------------------------------------------------
   !> numerical sigma (=volume*stressTensor) for debugging purposes
      !  generate a warning to keep release versions from calculating numerical gradients
      print'(/,"analytical sigma")'
      print *, sigma
      sdum=0.0_wp
      if(.not.allocated(gdum)) allocate(gdum(3,mol%n), source = 0.0_wp )
      ! save original xyz and lattice
      allocate( coord(3,mol%n), source = mol%xyz )
      latt=mol%lattice
      wf0 = chk
      nums = 0.0_wp
      do j = 1, 3
        do i = 1, j
          ! Only eps_ij=step the rest equals zero: eps_(kl.ne.ij)=0
            eps=0.0_wp
            eps(i,j)=sstep
            ! adjust position vectors and lattice to get er
            do k=1,3
              mol%xyz(k,:)  =   kron(k,1)*mol%xyz(1,:) + eps(k,1)*mol%xyz(1,:) + &
                            &   kron(k,2)*mol%xyz(2,:) + eps(k,2)*mol%xyz(2,:) + &
                            &   kron(k,3)*mol%xyz(3,:) + eps(k,3)*mol%xyz(3,:) 
              mol%lattice(k,1) = (kron(k,1) + eps(k,1))*mol%lattice(1,1) + &
                               & (kron(k,2) + eps(k,2))*mol%lattice(2,1) + &
                               & (kron(k,3) + eps(k,3))*mol%lattice(3,1)
              mol%lattice(k,2) = (kron(k,1) + eps(k,1))*mol%lattice(1,2) + &
                               & (kron(k,2) + eps(k,2))*mol%lattice(2,2) + &
                               & (kron(k,3) + eps(k,3))*mol%lattice(3,2)
              mol%lattice(k,3) = (kron(k,1) + eps(k,1))*mol%lattice(1,3) + &
                               & (kron(k,2) + eps(k,2))*mol%lattice(2,3) + &
                               & (kron(k,3) + eps(k,3))*mol%lattice(3,3)
            enddo
            chk = wf0
            call calc%singlepoint(env, mol, chk, printlevel-1, restart, &
               & er, gdum, sdum, hlgap, results)

            ! reset coordinates and lattice
            mol%xyz=coord
            mol%lattice=latt
            ! adjust position vectors and lattice to get el
            do k=1,3
              mol%xyz(k,:)  =  kron(k,1)*mol%xyz(1,:) - eps(k,1)*mol%xyz(1,:) + &
                            &  kron(k,2)*mol%xyz(2,:) - eps(k,2)*mol%xyz(2,:) + &
                            &  kron(k,3)*mol%xyz(3,:) - eps(k,3)*mol%xyz(3,:) 
              mol%lattice(k,1) = (kron(k,1) - eps(k,1))*mol%lattice(1,1) + &
                               & (kron(k,2) - eps(k,2))*mol%lattice(2,1) + &
                               & (kron(k,3) - eps(k,3))*mol%lattice(3,1)
              mol%lattice(k,2) = (kron(k,1) - eps(k,1))*mol%lattice(1,2) + &
                               & (kron(k,2) - eps(k,2))*mol%lattice(2,2) + &
                               & (kron(k,3) - eps(k,3))*mol%lattice(3,2)
              mol%lattice(k,3) = (kron(k,1) - eps(k,1))*mol%lattice(1,3) + &
                               & (kron(k,2) - eps(k,2))*mol%lattice(2,3) + &
                               & (kron(k,3) - eps(k,3))*mol%lattice(3,3)
            enddo
            chk = wf0
            call calc%singlepoint(env, mol, chk, printlevel-1, restart, &
               & el, gdum, sdum, hlgap, results)

            ! numerical sigma (=volume*stressTensor)
            nums(i,j) = sstep2 * (er - el)  ! divide by 2 times step size
            nums(j,i) = nums(i,j)  ! stress tensor is symmetric
            ! reset coordinates and lattice
            mol%xyz=coord
            mol%lattice=latt
         enddo
      enddo

      print'(/,"numerical sigma")'
      print *, nums
      print'(/,"difference sigma")'
      print*,sigma-nums
      chk = wf0
      numsig=nums
end subroutine num_sigma


function kron(i,j) result(res_kronij)  ! kronecker delta
  integer, intent(in) :: i,j
  real(wp) :: res_kronij

  res_kronij = 0.0_wp
  if(i.eq.j) res_kronij=1.0_wp

end function


end module xtb_pbc_optimizer_driver
