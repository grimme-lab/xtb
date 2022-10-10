! This file is part of xtb.
!
! Copyright (C) 2022 Stefan Grimme, Christoph Plett
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

module xtb_docking_param
   use xtb_mctc_accuracy, only: wp
   use xtb_mctc_symbols, only: toSymbol
   use xtb_type_setvar
   use xtb_type_molecule, only: TMolecule, init
   use xtb_param_covalentRadD3, only: covalentRadD3
   use xtb_mctc_param, only: sqrt_z_r4_over_r2
   use xtb_splitparam, only: atmass
   implicit none

   private :: toSymbol

   !> Parameter
   real(wp) :: par_rep_scal,par_d3_a1,par_d3_a2,par_d3_s8,par_chrg_lp,par_chrg_sig,&
             & par_chrg_pi,par_pos1_lp,par_pos2_lp,par_pos_sig,par_es_damp,par_drude_fc,&
             & par_drude_damp, par_xh1, par_xh2
   real(wp) :: r0ab,r0ab6(94,94),r0ab8(94,94),rrab(94,94),r0scal(94),val_e(86),r0_atom(86),rcov(94),r2r4(94)
   integer :: lpatom(86), sigatom(86)

   !> # of generations
   integer :: maxgen = 10

   !> # of parents for genetic algo
   integer :: maxparent = 100 

   !> step size for RG grid
   real(wp) :: stepr = 2.5    

   !> Step size for angular grid
   real(wp) :: stepa = 45    

   !> Maximal number of points CMA search
   integer :: mxcma = 1000    

   !> Size for pocket clusteing
   integer :: mxcent_clust = 500  

   !> Include input in gene pool?
   logical :: incl_org = .false.  

   !> Search types
   logical :: stack_only = .false.
   logical :: pocket_only = .false.
   logical :: no_pocket = .true.

   !> Probe atom type
   integer :: probe_atom_type = 36

   !> # of final geo. opts
   integer :: n_opt = 15 

   !> Mode
   integer :: mode = 0            

   !> CS symmetric molecule A?
   logical :: cssym = .false.     

   character*80 :: solvent, XTBHOME
   character(len=:), allocatable :: natom_arg
   integer :: natom_molA = 0

   !> Optimization accuracy
   real(wp) :: acc = 1.0_wp
   real(wp) :: shift_geo !How much molB is shifted away from molA to determine Topo and D4 coefficients
   real(wp) :: pre_e_A = 0.0_wp
   real(wp) :: pre_e_B = 0.0_wp

   !>Settings
   logical :: fulle = .false.
   logical :: hess = .false.
   logical :: debug = .false.
   integer :: nfrag1 = 0
   logical :: samerand = .false.
   logical :: test = .false.
   character(len=:), allocatable :: optlvl

   !>docklmocommon
   integer :: maxlmo = 50000
   real(wp) :: lmoint(3, 50000, 2) = 0.0_wp
   integer :: lmoatom(4, 50000, 2) = 0

   !>dockmolcommon
   real(wp) :: dipol(2) = 0.0_wp
   real(wp) :: grotrav(2) = 0.0_wp
   real(wp) :: hrotrav(2) = 0.0_wp
   real(wp) :: zpve(2) = 0.0_wp
   real(wp) :: extb(2) = 0.0_wp

   real(wp) :: chrg(2) = 0.0_wp
   real(wp) :: uhf(2) = 0.0_wp
   real(wp) :: elumo(2) = 0.0_wp
   real(wp) :: ehomo(2) = 0.0_wp

   !>Stuff for constraints
   logical :: constraint_xyz = .false.
   logical :: auto_wall = .false.
   character(len=:), allocatable :: xcontrol

   !> Options
   logical :: pocket_grid = .false.
   logical :: angular_grid = .true.
   logical :: stack_grid = .true.  

   !> General Control
   integer :: gsolvstate_iff
   logical :: docking_ens = .false.

   !> Fixing for Directed Docking
   type(fix_setvar) :: directedset
   !> Scaling factor for repulstive potential
   real(wp),parameter :: pot_scal = 1.0_wp
   !> Exponent for repulsitve potential
   real(wp),parameter :: pot_expo = 2.0_wp
   !> 20 kcal/mol attractive potential
   real(wp) :: attractive_pot = -0.032_wp
   integer :: directed_type = 3
   !> Repulsive atom-centered potential
   integer, parameter :: p_atom_pot = 2
   !> Attractive atom-centered potential
   integer, parameter :: p_atom_att = 3 
   !Wall pot for directed docking
   integer, parameter :: p_wall_pot = 1 
   integer :: place_wall_pot
   !QCG mode (special treatment of wall potentials)
   logical :: qcg = .false.

   !> Drude
   real(wp) :: gam(94)
   data gam/&
  &0.47259288, 0.92203391, 0.17452888, 0.25700733, 0.33949086, 0.42195412,&
  &0.50438193, 0.58691863, 0.66931351, 0.75191607, 0.17964105, 0.22157276,&
  &0.26348578, 0.30539645, 0.34734014, 0.38924725, 0.43115670, 0.47308269,&
  &0.17105469, 0.20276244, 0.21007322, 0.21739647, 0.22471039, 0.23201501,&
  &0.23933969, 0.24665638, 0.25398255, 0.26128863, 0.26859476, 0.27592565,&
  &0.30762999, 0.33931580, 0.37235985, 0.40273549, 0.43445776, 0.46611708,&
  &0.15585079, 0.18649324, 0.19356210, 0.20063311, 0.20770522, 0.21477254,&
  &0.22184614, 0.22891872, 0.23598621, 0.24305612, 0.25013018, 0.25719937,&
  &0.28784780, 0.31848673, 0.34912431, 0.37976593, 0.41040808, 0.44105777,&
  &0.05019332, 0.06762570, 0.08504445, 0.10247736, 0.11991105, 0.13732772,&
  &0.15476297, 0.17218265, 0.18961288, 0.20704760, 0.22446752, 0.24189645,&
  &0.25932503, 0.27676094, 0.29418231, 0.31159587, 0.32902274, 0.34592298,&
  &0.36388048, 0.38130586, 0.39877476, 0.41614298, 0.43364510, 0.45104014,&
  &0.46848986, 0.48584550, 0.12526730, 0.14268677, 0.16011615, 0.17755889,&
  &0.19497557, 0.21240778, 0.07263525, 0.09422158, 0.09920295, 0.10418621,&
  &0.14235633, 0.16394294, 0.18551941, 0.22370139/

contains

   subroutine set_iff_param

      character(len=80) :: fname
      logical :: ex
      integer :: i, j

      rcov = covalentRadD3(1:94)
      r2r4 = sqrt_z_r4_over_r2(1:94)
      call valel(val_e)

      fname = trim(XTBHOME)//'.param.xtbiff'
      inquire (file=fname, exist=ex)

      if (ex) then
         write (*, *) 'reading parameter file', trim(fname)
         open (unit=2, file=fname)
         read (2, *) par_rep_scal
         read (2, *) par_d3_a1
         read (2, *) par_d3_a2
         read (2, *) par_d3_s8
         read (2, *) par_chrg_lp    ! / 1.8 for Z > 10 in setlmo
         read (2, *) par_chrg_sig
         read (2, *) par_chrg_pi    ! x 3.0 for deloc pi in setlmo
         read (2, *) par_pos1_lp
         read (2, *) par_es_damp
         read (2, *) par_xh1
         read (2, *) par_xh2
         read (2, *) par_drude_fc
         do i = 1, 86
            read (2, *) j, r0scal(i)
         end do
         close (2)
      else
         write (*, *) 'taking internal default parameters'
         par_rep_scal = 0.1232398
         par_d3_a1 = 0.4459020
         par_d3_a2 = 4.4216615
         par_d3_s8 = 2.0000000
         par_chrg_lp = 3.1064458
         par_chrg_sig = 2.7000000
         par_chrg_pi = 0.8892689
         par_pos1_lp = 3.7082486
         par_es_damp = 0.2298506
         par_xh1 = 0.3614582
         par_xh2 = 6.7354042
         par_drude_fc = 0.9811899
         r0scal(1) = 0.6500000
         r0scal(2) = 0.6500000
         r0scal(3) = 0.8500000
         r0scal(4) = 1.1000000
         r0scal(5) = 1.2000000
         r0scal(6) = 1.2606324
         r0scal(7) = 1.3744700
         r0scal(8) = 1.1377692
         r0scal(9) = 0.9300000
         r0scal(10) = 0.8500000
         r0scal(11) = 0.8500000
         r0scal(12) = 1.3000000
         r0scal(13) = 1.3500000
         r0scal(14) = 1.3500000
         r0scal(15) = 1.6000000
         r0scal(16) = 1.6000000
         r0scal(17) = 1.6000000
         r0scal(18) = 0.9000000
         r0scal(19) = 0.8500000
         r0scal(20) = 1.4000000
         r0scal(21) = 1.3500000
         r0scal(22) = 1.3500000
         r0scal(23) = 1.3500000
         r0scal(24) = 1.3500000
         r0scal(25) = 1.3500000
         r0scal(26) = 1.3500000
         r0scal(27) = 1.6000000
         r0scal(28) = 1.6000000
         r0scal(29) = 1.6000000
         r0scal(30) = 1.6000000
         r0scal(31) = 1.5500000
         r0scal(32) = 1.5500000
         r0scal(33) = 1.5500000
         r0scal(34) = 1.6000000
         r0scal(35) = 1.6000000
         r0scal(36) = 0.9000000
         r0scal(37) = 0.9500000
         r0scal(38) = 1.5000000
         r0scal(39) = 1.3500000
         r0scal(40) = 1.3500000
         r0scal(41) = 1.6000000
         r0scal(42) = 1.3500000
         r0scal(43) = 1.6500000
         r0scal(44) = 1.3500000
         r0scal(45) = 1.6000000
         r0scal(46) = 1.6000000
         r0scal(47) = 1.6000000
         r0scal(48) = 1.6000000
         r0scal(49) = 1.3500000
         r0scal(50) = 1.3500000
         r0scal(51) = 1.5500000
         r0scal(52) = 1.5500000
         r0scal(53) = 1.6000000
         r0scal(54) = 0.9500000
         r0scal(55) = 0.9500000
         r0scal(56) = 1.5000000
         r0scal(57) = 1.3500000
         r0scal(58) = 1.6000000
         r0scal(59) = 1.6000000
         r0scal(60) = 1.6000000
         r0scal(61) = 1.6000000
         r0scal(62) = 1.6000000
         r0scal(63) = 1.6000000
         r0scal(64) = 1.6000000
         r0scal(65) = 1.6000000
         r0scal(66) = 1.6000000
         r0scal(67) = 1.6000000
         r0scal(68) = 1.6000000
         r0scal(69) = 1.6000000
         r0scal(70) = 1.6000000
         r0scal(71) = 1.6000000
         r0scal(72) = 1.3500000
         r0scal(73) = 1.6000000
         r0scal(74) = 1.6000000
         r0scal(75) = 1.6000000
         r0scal(76) = 1.6000000
         r0scal(77) = 1.6000000
         r0scal(78) = 1.6000000
         r0scal(79) = 1.6000000
         r0scal(80) = 1.6000000
         r0scal(81) = 1.4500000
         r0scal(82) = 1.4500000
         r0scal(83) = 1.4500000
         r0scal(84) = 1.6000000
         r0scal(85) = 1.6000000
         r0scal(86) = 0.9500000

      end if

!     not fully fitted ("hand made")
!     four additional CT parameters in xtbdock_energy.f (chrgtransfer)
      par_drude_damp = 0.020! induction damp
      par_pos_sig = par_pos1_lp
      par_pos2_lp = 0.0

! which atom has LP ES interactions?
      lpatom = 0
      lpatom(6:9) = 1
      lpatom(14:17) = 1
      lpatom(32:35) = 1
      lpatom(50:53) = 1
      lpatom(82:85) = 1

! which atom has sigma quadrupole ES interactions?
      sigatom = 0
      sigatom(9) = 1
      sigatom(17) = 1
      sigatom(35) = 1
      sigatom(53) = 1
      sigatom(85) = 1

   end subroutine set_iff_param

   subroutine diptot(n1, n2, nl1, nl2, c1, c2, cl1, cl2, l1, l2, q1, q2, dipol)
      integer, intent(in) :: n1, n2
      integer, intent(in) :: nl1, nl2
      real(wp), intent(in) :: c1(3, n1)
      real(wp), intent(in) :: c2(3, n2)
      real(wp), intent(in) :: cl1(4, n1*10)
      real(wp), intent(in) :: cl2(4, n2*10)
      real(wp), intent(in) :: q1(n1)
      real(wp), intent(in) :: q2(n2)
      integer, intent(in) :: l1(n1*10)
      integer, intent(in) :: l2(n2*10)
      real(wp), intent(out) :: dipol(3)

      integer :: i

      dipol = 0
      do i = 1, nl1
         dipol(1:3) = dipol(1:3) + cl1(1:3, i)*cl1(4, i)    ! one off-charge at LMO
      end do
      do i = 1, n1
         dipol(1:3) = dipol(1:3) + c1(1:3, i)*q1(i)
      end do
      do i = 1, nl2
         dipol(1:3) = dipol(1:3) + cl2(1:3, i)*cl2(4, i)    ! one off-charge at LMO
      end do
      do i = 1, n2
         dipol(1:3) = dipol(1:3) + c2(1:3, i)*q2(i)
      end do

   end subroutine diptot

! compute center of mass(sum3) and moment of intertia and corresponding
! axis
! molw is the weigth, sum3 the CMA (all in a.u.)

   subroutine axis(numat, nat, coord, sum3, sumw, eig, evec)

      integer, intent(in) ::numat, nat(numat)
      real(wp), intent(in) :: coord(3, *)
      real(wp), intent(out) :: sum3(3), sumw, eig(3), evec(3, 3)

      real(wp) :: t(6)
      real(wp) :: x(numat), y(numat), z(numat)
      real(wp) :: sumwx, sumwy, sumwz
      real(wp) :: atmass
      integer :: i
      real(wp) ::  ams(107)
      data ams/1.00790d0, 4.00260d0, 6.94000d0, 9.01218d0,&
     & 10.81000d0, 12.01100d0, 14.00670d0, 15.99940d0, 18.99840d0,&
     & 20.17900d0, 22.98977d0, 24.30500d0, 26.98154d0, 28.08550d0,&
     & 30.97376d0, 32.06000d0, 35.45300d0, 39.94800d0, 39.09830d0,&
     & 40.08000d0, 44.95590d0, 47.90000d0, 50.94150d0, 51.99600d0,&
     & 54.93800d0, 55.84700d0, 58.93320d0, 58.71000d0, 63.54600d0,&
     & 65.38000d0, 69.73500d0, 72.59000d0, 74.92160d0, 78.96000d0,&
     & 79.90400d0, 83.80000d0, 85.46780d0, 87.62000d0, 88.90590d0,&
     & 91.22000d0, 92.90640d0, 95.94000d0, 98.90620d0, 101.0700d0,&
     & 102.9055d0, 106.4000d0, 107.8680d0, 112.4100d0, 114.8200d0,&
     & 118.6900d0, 121.7500d0, 127.6000d0, 126.9045d0, 131.3000d0,&
     & 132.9054d0, 137.3300d0, &
     & 138.91, 140.12, 140.91, 144.24, 147.00, 150.36, 151.97, 157.25,&
     & 158.93, 162.50, 164.93, 167.26, 168.93, 173.04, 174.97,&
     & 178.4900d0, 180.9479d0,&
     & 183.8500d0, 186.2070d0, 190.2000d0, 192.2200d0, 195.0900d0,&
     & 196.9665d0, 200.5900d0, 204.3700d0, 207.2000d0, 208.9804d0,&
     & 209., 210., 222., 21*0.000d0/

      sumw = 0
      sumwx = 0.d0
      sumwy = 0.d0
      sumwz = 0.d0

      do i = 1, numat
         atmass = ams(nat(i))
         sumw = sumw + atmass
         sumwx = sumwx + atmass*coord(1, i)
         sumwy = sumwy + atmass*coord(2, i)
         sumwz = sumwz + atmass*coord(3, i)
      end do

      sum3(1) = sumwx/sumw
      sum3(2) = sumwy/sumw
      sum3(3) = sumwz/sumw

      do i = 1, numat
         x(i) = coord(1, i) - sum3(1)
         y(i) = coord(2, i) - sum3(2)
         z(i) = coord(3, i) - sum3(3)
      end do

      do i = 1, 6
         t(i) = real(i, wp)*1.0e-10_wp
      end do

      do i = 1, numat
         atmass = ams(nat(i))
         t(1) = t(1) + atmass*(y(i)**2 + z(i)**2)
         t(2) = t(2) - atmass*x(i)*y(i)
         t(3) = t(3) + atmass*(z(i)**2 + x(i)**2)
         t(4) = t(4) - atmass*z(i)*x(i)
         t(5) = t(5) - atmass*y(i)*z(i)
         t(6) = t(6) + atmass*(x(i)**2 + y(i)**2)
      end do

      call rsp(t, 3, 3, eig, evec)
      eig = eig/sumw

   end subroutine axis

   subroutine cmadock(n, numat, nat, coord, sum3)

      integer, intent(in) :: n, numat, nat(numat)
      real(wp), intent(in) :: coord(3, numat)
      real(wp), intent(out) :: sum3(3)
      real(wp)  :: ams(107)
      real(wp) :: sumw, sumwx, sumwy, sumwz
      integer :: i
      data ams/1.00790d0, 4.00260d0, 6.94000d0, 9.01218d0,&
     & 10.81000d0, 12.01100d0, 14.00670d0, 15.99940d0, 18.99840d0,&
     & 20.17900d0, 22.98977d0, 24.30500d0, 26.98154d0, 28.08550d0,&
     & 30.97376d0, 32.06000d0, 35.45300d0, 39.94800d0, 39.09830d0,&
     & 40.08000d0, 44.95590d0, 47.90000d0, 50.94150d0, 51.99600d0,&
     & 54.93800d0, 55.84700d0, 58.93320d0, 58.71000d0, 63.54600d0,&
     & 65.38000d0, 69.73500d0, 72.59000d0, 74.92160d0, 78.96000d0,&
     & 79.90400d0, 83.80000d0, 85.46780d0, 87.62000d0, 88.90590d0,&
     & 91.22000d0, 92.90640d0, 95.94000d0, 98.90620d0, 101.0700d0,&
     & 102.9055d0, 106.4000d0, 107.8680d0, 112.4100d0, 114.8200d0,&
     & 118.6900d0, 121.7500d0, 127.6000d0, 126.9045d0, 131.3000d0,&
     & 132.9054d0, 137.3300d0, &
     & 138.91, 140.12, 140.91, 144.24, 147.00, 150.36, 151.97, 157.25,&
     & 158.93, 162.50, 164.93, 167.26, 168.93, 173.04, 174.97,&
     & 178.4900d0, 180.9479d0,&
     & 183.8500d0, 186.2070d0, 190.2000d0, 192.2200d0, 195.0900d0,&
     & 196.9665d0, 200.5900d0, 204.3700d0, 207.2000d0, 208.9804d0,&
     & 209., 210., 222., 21*0.000d0/
! atomic masses

      sumw = 1.e-20_wp
      sumwx = 0.0_wp
      sumwy = 0.0_wp
      sumwz = 0.0_wp

      do i = 1, n
         sumw = sumw + ams(nat(i))
         sumwx = sumwx + ams(nat(i))*coord(1, i)
         sumwy = sumwy + ams(nat(i))*coord(2, i)
         sumwz = sumwz + ams(nat(i))*coord(3, i)
      end do

      sum3(1) = sumwx/sumw
      sum3(2) = sumwy/sumw
      sum3(3) = sumwz/sumw

   end subroutine cmadock

   subroutine rcma(n1, xyz1, iz1, n2, xyz2, iz2, r, rmin)

      integer, intent(in) :: n1, n2
      real(wp), intent(in) :: xyz1(3, n1)
      real(wp), intent(in) :: xyz2(3, n2)
      integer, intent(in) :: iz1(n1)
      integer, intent(in) :: iz2(n2)
      real(wp), intent(out) :: r, rmin

      real(wp) :: x1(3)
      real(wp) :: x2(3)
      integer :: i, j
      real(wp) :: rr

      call cmadock(n1, n1, iz1, xyz1, x1)
      call cmadock(n2, n2, iz2, xyz2, x2)

      r = sqrt((x1(1) - x2(1))**2 + (x1(2) - x2(2))**2 + (x1(3) - x2(3))**2)

      rmin = 1.d+42
      do i = 1, n1
         do j = 1, n2
            rr = sqrt((xyz1(1, i) - xyz2(1, j))**2 +&
               &(xyz1(2, i) - xyz2(2, j))**2 +&
               &(xyz1(3, i) - xyz2(3, j))**2)
            if (rr .lt. rmin) rmin = rr
         end do
      end do

   end subroutine rcma

! xtbdock version
   subroutine rotmat(rxyz, rot)
      real(wp), intent(in) :: rxyz(6)
      real(wp), intent(out) :: rot(3, 3)

      real(wp) :: r1(3, 3)
      real(wp) :: r2(3, 3)
      real(wp) :: r3(3, 3)
      real(wp) :: tmp(3, 3)

      r1(1, 1) = 1.0
      r1(1, 2) = 0.0
      r1(1, 3) = 0.0
      r1(2, 1) = 0.0
      r1(2, 2) = cos(rxyz(4))
      r1(2, 3) = -sin(rxyz(4))
      r1(3, 1) = 0.0
      r1(3, 2) = sin(rxyz(4))
      r1(3, 3) = cos(rxyz(4))

      r2(2, 1) = 0.0
      r2(2, 2) = 1.0
      r2(2, 3) = 0.0
      r2(1, 1) = cos(rxyz(5))
      r2(1, 2) = 0.0
      r2(1, 3) = sin(rxyz(5))
      r2(3, 1) = -sin(rxyz(5))
      r2(3, 2) = 0.0
      r2(3, 3) = cos(rxyz(5))

      r3(3, 1) = 0.0
      r3(3, 2) = 0.0
      r3(3, 3) = 1.0
      r3(1, 1) = cos(rxyz(6))
      r3(1, 2) = -sin(rxyz(6))
      r3(1, 3) = 0.0
      r3(2, 1) = sin(rxyz(6))
      r3(2, 2) = cos(rxyz(6))
      r3(2, 3) = 0.0

      tmp = matmul(r2, r3)
      rot = matmul(r1, tmp)

   end subroutine rotmat

! rotate 3D-vector array around Euler angles ang
   subroutine rot3(n, xyz, rotm)

      integer, intent(in) :: n
      real(wp), intent(inout) :: xyz(3, n)
      real(wp), intent(in) :: rotm(3, 3)

      integer :: i

      xyz = matmul(rotm, xyz)

   end subroutine rot3

   subroutine rot4(n, xyz, rotm)

      integer, intent(in) :: n
      real(wp), intent(inout) :: xyz(4, 10*n)
      real(wp), intent(in) :: rotm(3, 3)

      real(wp) :: tmp(3, n)
      integer :: i

      tmp(1:3, 1:n) = xyz(1:3, 1:n)
      tmp = matmul(rotm, tmp)
      xyz(1:3, 1:n) = tmp(1:3, 1:n)

   end subroutine rot4

   subroutine wrc0(fname, n1, at1, xyz1)

      character(len=*), intent(in) :: fname
      integer, intent(in) :: n1, at1(n1)
      real(wp), intent(in) :: xyz1(3, n1)

      integer :: j

      open (unit=1, file=fname)
      write (1, '(''$coord'')')
      do j = 1, n1
         write (1, '(3F24.10,5x,a2)') xyz1(1:3, j), toSymbol(at1(j))
      end do
      write (1, '(''$end'')')
      close (1)

   end subroutine wrc0

   subroutine wrc(fname, n1, n2, at1, at2, xyz1, xyz2, icoord)

      character(len=*), intent(in) :: fname
      integer, intent(in) :: n1, at1(n1)
      integer, intent(in) :: n2, at2(n2)
      real(wp), intent(in) :: xyz1(3, n1)
      real(wp), intent(in) :: xyz2(3, n2)
      real(wp), intent(in) :: icoord(6)

      real(wp) ::   c2(3, n2)
      integer :: j

      call move2(n2, xyz2, c2, icoord)

      open (unit=1, file=fname)
      write (1, '(''$coord'')')
      do j = 1, n1
         write (1, '(3F24.10,5x,a2)') xyz1(1:3, j), toSymbol(at1(j))
      end do
      do j = 1, n2
         write (1, '(3F24.10,5x,a2)') c2(1:3, j), toSymbol(at2(j))
      end do
      write (1, '(''$end'')')
      close (1)

   end subroutine wrc

   subroutine wrc2(fname, mode, n1, n2, at1, at2, xyz1, xyz2, nf, found)

      character(len=*), intent(in) :: fname
      integer, intent(in) :: nf, mode
      integer, intent(in) :: n1, at1(n1)
      integer, intent(in) :: n2, at2(n2)
      real(wp), intent(in) :: found(7, *)
      real(wp), intent(in) :: xyz1(3, n1)
      real(wp), intent(in) :: xyz2(3, n2)

      real(wp) ::   c2(3, n2)
      real(wp) :: icoord(6)
      real(wp), parameter :: autoang = 0.52917726d0
      integer :: i, j

      open (unit=33, file=fname)
      do i = 1, nf
         write (33, *) n1 + n2
         write (33, '('' SCF done '',F16.8)') found(7, i)
         icoord(1:6) = found(1:6, i)
         call move2(n2, xyz2, c2, icoord)
         do j = 1, n1
            write (33, '(a2,3F24.10)') toSymbol(at1(j)), xyz1(1:3, j)*autoang
         end do
         do j = 1, n2
            if (mode .eq. 0) then
               write (33, '(a2,3F24.10)') toSymbol(at2(j)), c2(1:3, j)*autoang !  element symbol
            else
               write (33, '(a2,3F24.10)') toSymbol(5), c2(1:3, j)*autoang !  Bor for visualization
            end if
         end do
      end do
      close (33)

   end subroutine wrc2

   subroutine wrc3(fname, n1, n2, at1, at2, xyz1, xyz2, i, found)

      character(len=*), intent(in) :: fname
      integer, intent(in) :: i
      integer, intent(in) :: n1, at1(n1)
      integer, intent(in) :: n2, at2(n2)
      real(wp), intent(in) :: found(7, *)
      real(wp), intent(in) :: xyz1(3, n1)
      real(wp), intent(in) :: xyz2(3, n2)

      real(wp) ::   c2(3, n2)
      real(wp) :: icoord(6)
      real(wp), parameter :: autoang = 0.52917726d0
      real(wp), parameter :: au = 627.509541d0
      real(wp) :: f
      integer :: j

      f = 1.
      if (i .eq. 0) f = au
      open (unit=33, file=fname)
      write (33, *) n1 + n2
      write (33, '('' SCF done '',F16.8)') found(7, 1)*f
      icoord(1:6) = found(1:6, 1)
      call move2(n2, xyz2, c2, icoord)
      do j = 1, n1
         write (33, '(a2,3F24.10)') toSymbol(at1(j)), xyz1(1:3, j)*autoang
      end do
      do j = 1, n2
         write (33, '(a2,3F24.10)') toSymbol(at2(j)), c2(1:3, j)*autoang !  element symbol
      end do
      close (33)

   end subroutine wrc3

   subroutine wr_grid(fname, n1, n2, at1, element, xyz1, xyz2)

      character(len=*), intent(in) :: fname
      integer, intent(in) :: n1, at1(n1)
      integer, intent(in) :: n2
      integer, intent(in) :: element
      real(wp), intent(in) :: xyz1(3, n1)
      real(wp), intent(in) :: xyz2(3, n2)

      real(wp), parameter :: autoang = 0.52917726d0
      integer :: i, j

      open (unit=33, file=fname)
      write (33, *) n1 + n2
      write (33, *)
      do j = 1, n1
         write (33, '(a2,3F24.10)') toSymbol(at1(j)), xyz1(1:3, j)*autoang
      end do
      do j = 1, n2
         write (33, '(a2,3F24.10)') toSymbol(element), xyz2(1:3, j)*autoang !  Bor for visualization
      end do
      close (33)

   end subroutine wr_grid

! new cartesian cooordinates from internals
! move coordinates and LP centers
   subroutine move(n2, nl2, xyz, cl, xyznew, clnew, rotm, coord)

      integer, intent(in) :: n2, nl2
      real(wp), intent(in) :: xyz(3, n2)
      real(wp), intent(out) :: xyznew(3, n2)
      real(wp), intent(in) :: cl(4, n2*10)
      real(wp), intent(out) :: clnew(4, n2*10)
      real(wp), intent(in) :: coord(6)
      real(wp), intent(in) :: rotm(3, 3)

      integer :: k

!     rotate (its at the CMA)
      xyznew = xyz
      clnew = cl
      if (sum(abs(coord)) .lt. 1.d-12) return

      call rot3(n2, xyznew, rotm)
      call rot4(nl2, clnew, rotm)

!     shift
      do k = 1, 3
         xyznew(k, 1:n2) = xyznew(k, 1:n2) + coord(k)
         clnew(k, 1:nl2) = clnew(k, 1:nl2) + coord(k)
      end do

   end subroutine move

! move just coordinates
   subroutine move2(n2, xyz, xyznew, coord)

      integer, intent(in) :: n2
      real(wp), intent(in) :: xyz(3, n2)
      real(wp), intent(out) :: xyznew(3, n2)
      real(wp), intent(in) :: coord(6)

      real(wp) :: rotm(3, 3)
      integer :: k

!     rotate (its at the CMA)
      xyznew = xyz
      call rotmat(coord, rotm)
      call rot3(n2, xyznew, rotm)

!     shift
      do k = 1, 3
         xyznew(k, 1:n2) = xyznew(k, 1:n2) + coord(k)
      end do

   end subroutine move2

   subroutine rnorm(i, j, found, r)

      integer, intent(in) :: i, j
      real(wp), intent(in) :: found(7, *)
      real(wp), intent(out) :: r

      integer :: k
      real(wp) :: x

      r = 0
      do k = 1, 3
         r = r + (found(k, i) - found(k, j))**2
      end do
      do k = 4, 6
         x = found(k, i) - found(k, j)
         r = r + sin(0.5*x)**2
      end do

      r = sqrt(r)

   end subroutine rnorm

! random stuff
   subroutine rand6(mut, r1, r2, d) !Currently mut=0.5, r1=r2=1.0

      real(wp), intent(in) :: mut, r1, r2
      real(wp), intent(inout) :: d(6)

      real(wp), parameter :: pi2 = 2*3.14159265358979_wp
      real(wp) :: x, y, yy, f
      integer :: j, k

      call random_number(yy)
      if (yy .lt. mut) return    ! 50 % mutation rate

      do j = 1, 3                ! Shift position by either -x or +x (0<x<1)
         call random_number(y)
         f = 1.0
         if (y .lt. 0.5) f = -1.0
         call random_number(x)
         d(j) = d(j) + f*r1*x
      end do

      if (yy .gt. 0.95) then       !in 5 % of the cases turn around an random axis by 180 deg.
         call random_number(y)
         call irand(3, k)
         f = 0.5
         if (y .lt. 0.5) f = -0.5
         d(k + 3) = d(k + 3) + f*pi2
         return
      end if

      do j = 4, 6                  !Change angle by either -x or +x (0<x<1)
         call random_number(y)
         f = 1.0
         if (y .lt. 0.5) f = -1.0
         call random_number(x)
         d(j) = d(j) + f*r2*x
         if (d(j) .gt. pi2) d(j) = d(j) - pi2 !To not get lower than 0°
         if (d(j) .lt. 0) d(j) = d(j) + pi2   !To not get larger than 360°
      end do

   end subroutine rand6

! integer random number n=<irand<=1
   subroutine irand(n, k)
      integer, intent(in) :: n
      integer, intent(out) :: k

      real :: x, nx

      call random_number(x)
      nx = n*x + 1
      k = int(nx)
      if (k .gt. n) k = n

   end subroutine irand

   subroutine crossover(r, d)
      real(wp), intent(in) :: r
      real(wp), intent(inout) :: d(6)

      integer :: j
      real :: x

      do j = 1, 6
         call random_number(x)
         if (x .gt. r) then
            d(j) = x
         else
            cycle
         end if
      end do

   end subroutine crossover

! sorts
   subroutine sort6(n, e)

      integer, intent(in) :: n
      real(wp), intent(inout) :: e(7, n)

      integer :: ind(n), i
      real(wp) :: tmp(n), ftmp(7, n)

      ftmp = e
      do i = 1, n
         ind(i) = i
         tmp(i) = e(7, i)
      end do

      call Qsort(tmp, 1, n, ind) ! sort

      do i = 1, n
         e(1:7, i) = ftmp(1:7, ind(i))
      end do

   end subroutine sort6

! write complex coord file
   subroutine wrlmocoord(fname, n1, n2, xyz1, xyz2, at1, at2, nlmo1, nlmo2,&
  &                      lmo1, lmo2, rlmo1, rlmo2)

      character(len=*), intent(in) :: fname
      integer, intent(in) :: n1, nlmo1, at1(n1), lmo1(n1*10)
      integer, intent(in) :: n2, nlmo2, at2(n2), lmo2(n2*10)
      real(wp), intent(in) :: xyz1(3, n1), rlmo1(4, n1*10)
      real(wp), intent(in) :: xyz2(3, n2), rlmo2(4, n2*10)

      integer :: i

      open (unit=32, file=fname)
      write (32, '(''$coord'')')
      do i = 1, n1
         write (32, '(3F24.10,5x,a2)') xyz1(1:3, i), toSymbol(at1(i))
      end do
      do i = 1, n2
         write (32, '(3F24.10,5x,a2)') xyz2(1:3, i), toSymbol(at2(i))
      end do
      do i = 1, nlmo1
         write (32, '(3F24.10,5x,a2)') rlmo1(1:3, i), toSymbol(2)
      end do
      do i = 1, nlmo2
         write (32, '(3F24.10,5x,a2)') rlmo2(1:3, i), toSymbol(2)
      end do
      write (32, '(''$end'')')
      close (32)

   end subroutine wrlmocoord

   ! recursively split molecule xyz,at into non-cov. bound fragments
   subroutine mrec(molcount, xyz, nat, at, molvec)

      ! molcount: number of total fragments (increased during search)
      ! xyz: overall Cart. coordinates
      ! nat: overall number of atoms
      ! at: atomic number array
      ! molvec: assignment vector of atom to fragment

      real(wp), intent(in) :: xyz(3, nat)
      integer, intent(in) :: nat, at(nat)
      integer, intent(out) :: molvec(nat), molcount

      real(wp) :: cn(nat), bond(nat, nat)
      integer :: i
      logical :: taken(nat)

      molvec = 0
      molcount = 1
      taken = .false.

      call xcoord(nat, at, xyz, cn, bond)

      do i = 1, nat
         if (.not. taken(i)) then
            molvec(i) = molcount
            taken(i) = .true.
            call neighbours(i, xyz, at, taken, nat, cn, bond, molvec, molcount)
            molcount = molcount + 1
         end if
      end do
      molcount = molcount - 1
   end subroutine mrec

   recursive subroutine neighbours(i, xyz, iat, taken, nat, cn, bond,&
        &                                molvec, molcnt)

      real(wp), intent(in) :: xyz(3, nat)
      real(wp), intent(inout) :: cn(nat), bond(nat, nat)
      integer, intent(in) :: i, nat, molcnt
      integer, intent(inout) :: molvec(nat)
      logical :: taken(nat)

      real(wp) :: r, tr, xi(3), cntmp
      integer :: j, iat(nat), icn, k

      tr = 2.0d0

      xi(1:3) = xyz(1:3, i)
      icn = nint(cn(i))
      do k = 1, icn
         j = maxloc(bond(:, i), 1)
         bond(j, i) = 0.0d0
         if (i .eq. j) cycle
         if (.not. taken(j)) then
            molvec(j) = molcnt
            taken(j) = .true.
            call neighbours(j, xyz, iat, taken, nat, cn, bond, molvec, molcnt)
         end if
      end do
   end subroutine neighbours

   ! compute coordination numbers by adding an inverse damping function
   subroutine xcoord(natoms, iz, xyz, cn, bond)

      integer, intent(in) :: iz(natoms), natoms
      real(wp), intent(in) :: xyz(3, natoms)
      real(wp), intent(out) :: bond(natoms, natoms), cn(natoms)

      integer :: iat, i, j
      real(wp) :: xn

      bond = 0.0d0
      cn = 0.0d0

      do i = 1, natoms
         call xdamp(natoms, i, xyz, iz, rcov, xn, bond)
         cn(i) = xn
      end do

   end subroutine xcoord

   subroutine xdamp(natoms, i, xyz, iz, rcv, xn, bond)

      integer, intent(in) :: natoms, iz(natoms), i
      real(wp), intent(in) :: xyz(3, natoms), rcv(94)
      real(wp), intent(out) :: bond(natoms, natoms), xn

      real(wp) :: dx, dy, dz, r, damp, rr, rco, r2
      integer :: iat

      xn = 0.0d0
      do iat = 1, natoms
         if (iat .eq. i) cycle
         dx = xyz(1, iat) - xyz(1, i)
         dy = xyz(2, iat) - xyz(2, i)
         dz = xyz(3, iat) - xyz(3, i)
         r2 = dx*dx + dy*dy + dz*dz
         r = sqrt(r2)
         ! covalent distance in Bohr
         rco = (rcv(iz(i)) + rcv(iz(iat)))*0.75! this scaling reduces the size of the clusters
         ! counting function exponential has a better long-range behavior than MHGs inverse damping
         damp = 1.d0/(1.d0 + exp(-16.0d0*(rco/r - 1.0d0)))
         bond(iat, i) = damp
         xn = xn + damp
      end do

   end subroutine xdamp

   subroutine valel(z)

      real(wp), intent(out) :: z(86)
      z(1) = 1.0
      z(2) = 2.0
      z(3) = 1.0
      z(4) = 2.0
      z(5) = 3.0
      z(6) = 4.0
      z(7) = 5.0
      z(8) = 6.0
      z(9) = 7.0
      z(10) = 8.0
      z(11) = 1.0
      z(12) = 2.0
      z(13) = 3.0
      z(14) = 4.0
      z(15) = 5.0
      z(16) = 6.0
      z(17) = 7.0
      z(18) = 8.0
      z(19) = 1.0
      z(20) = 2.0
      z(21) = 3.0
      z(22) = 4.0
      z(23) = 5.0
      z(24) = 6.0
      z(25) = 7.0
      z(26) = 8.0
      z(27) = 9.0
      z(28) = 10.0
      z(29) = 11.0
      z(30) = 2.0
      z(31) = 3.0
      z(32) = 4.0
      z(33) = 5.0
      z(34) = 6.0
      z(35) = 7.0
      z(36) = 8.0
      z(37) = 1.0
      z(38) = 2.0
      z(39) = 3.0
      z(40) = 4.0
      z(41) = 5.0
      z(42) = 6.0
      z(43) = 7.0
      z(44) = 8.0
      z(45) = 9.0
      z(46) = 10.0
      z(47) = 11.0
      z(48) = 2.0
      z(49) = 3.0
      z(50) = 4.0
      z(51) = 5.0
      z(52) = 6.0
      z(53) = 7.0
      z(54) = 8.0
      z(55) = 1.0
      z(56) = 2.0
      z(57) = 3.0
      z(58) = 3.0
      z(59) = 3.0
      z(60) = 3.0
      z(61) = 3.0
      z(62) = 3.0
      z(63) = 3.0
      z(64) = 3.0
      z(65) = 3.0
      z(66) = 3.0
      z(67) = 3.0
      z(68) = 3.0
      z(69) = 3.0
      z(70) = 3.0
      z(71) = 3.0
      z(72) = 3.0
      z(73) = 3.0
      z(74) = 3.0
      z(75) = 3.0
      z(76) = 3.0
      z(77) = 3.0
      z(78) = 3.0
      z(79) = 3.0
      z(80) = 3.0
      z(81) = 3.0
      z(82) = 4.0
      z(83) = 5.0
      z(84) = 6.0
      z(85) = 7.0
      z(86) = 8.0

   end subroutine valel

   subroutine get_mol(istart, iend, comb, mol)

      integer, intent(in) :: istart, iend

      type(TMolecule), intent(in) :: comb

      type(TMolecule), intent(out) :: mol

      integer, allocatable :: at(:)
      real(wp), allocatable :: xyz(:, :)
      integer :: n

      write (*, *) 'istart,iend', istart, iend
      n = iend - istart + 1

      allocate (at(n), xyz(3, n))

      at(1:iend) = comb%at(istart:iend)
      xyz(1:3, 1:iend) = comb%xyz(1:3, istart:iend)

      write (*, *) 'all', istart, iend, at, xyz

      call init(mol, at, xyz)

      deallocate (at, xyz)
   end subroutine get_mol

   subroutine split_mol(molA, molB, size_list, list, comb)

      type(TMolecule), intent(in) :: comb
      integer, intent(in) :: size_list, list(size_list)

      type(TMolecule), intent(inout) :: molA, molB

      integer :: at1(molA%n), at2(molB%n)
      integer :: i
      real(wp) :: xyz1(3, molA%n), xyz2(3, molB%n)

      do i=1, comb%n
         if(any(list == i)) then
            at1 = comb%at(1:molA%n)
            xyz1(1:3, :) = comb%xyz(1:3, 1:molA%n)
         else
            at2 = comb%at(molA%n + 1:molA%n + molB%n)
            xyz2(1:3, :) = comb%xyz(1:3, molA%n + 1:molA%n + molB%n)
         end if
      end do

      call init(molA, at1, xyz1)
      call init(molB, at2, xyz2)

   end subroutine split_mol

   subroutine combine_mol(comb, molA, molB)

      !> Molecular structure data
      type(TMolecule), intent(in) :: molA, molB
      type(TMolecule), intent(inout) :: comb

      integer, allocatable :: at(:)
      real(wp), allocatable :: xyz(:, :)
      integer :: i, j

      comb%n = molA%n + molB%n
      allocate (at(comb%n), source=0)
      allocate (xyz(3, comb%n), source=0.0_wp)

      !> Adding coords
      do i = 1, molA%n
         xyz(1, i) = molA%xyz(1, i)
         xyz(2, i) = molA%xyz(2, i)
         xyz(3, i) = molA%xyz(3, i)
      end do

      j = 1
      do i = molA%n + 1, comb%n
         xyz(1, i) = molB%xyz(1, j) + shift_geo
         xyz(2, i) = molB%xyz(2, j) + shift_geo
         xyz(3, i) = molB%xyz(3, j) + shift_geo
         j = j + 1
      end do

      !> Adding Atom-Types
      do i = 1, molA%n
         at(i) = molA%at(i)
      end do

      j = 1
      do i = molA%n + 1, comb%n
         at(i) = molB%at(j)
         j = j + 1
      end do

      if (allocated(comb%sym)) deallocate (comb%sym)
      allocate (comb%sym(comb%n), source='    ')
      deallocate (comb%sym)
      call init(comb, at, xyz)
      comb%chrg = molA%chrg + molB%chrg
      comb%uhf = molA%uhf + molB%uhf
      deallocate (at)
      deallocate (xyz)

      end subroutine combine_mol

end module xtb_docking_param
