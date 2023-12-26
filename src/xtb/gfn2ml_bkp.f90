! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> GFN2-xTB parametrisation data
module xtb_xtb_gfn2ml
   use xtb_mctc_accuracy, only : wp
   use xtb_param_atomicrad, only : atomicRad
   use xtb_param_paulingen, only : paulingEN
   use xtb_type_param, only : TxTBParameter, dftd_parameter
   use xtb_xtb_data
   use xtb_xtb_gfn1, only : setGFN1ShellHardness
   use xtb_disp_dftd4, only : newD4Model, p_refq_gfn2xtb
   implicit none
   private

   public :: initGFN2ML

   interface initGFN2ML
      module procedure :: initData
      module procedure :: initRepulsion
      module procedure :: initDispersion
      module procedure :: initCoulomb
      module procedure :: initMultipole
      module procedure :: initHamiltonian
   end interface initGFN2ML

   interface loadArray
         module procedure :: loadArray1D
         module procedure :: loadArray2D
   end interface loadArray


   real(wp), parameter :: gam3shell(2, 0:3) = reshape(&
      &[1.0_wp, 1.0_wp, 0.5_wp, 0.5_wp, 0.25_wp, 0.25_wp, 0.25_wp, 0.25_wp], &
      & shape(gam3shell))

   real(wp), parameter :: kshell(4) = [1.85_wp, 2.23_wp, 2.23_wp, 2.23_wp]

   type(TxTBParameter), parameter :: gfn2Globals = TxTBParameter( &
      kshell = kshell, &
      enshell = 2.0_wp, &
      ksd = 2.0_wp, &
      kpd = 2.0_wp, &
      kdiff = 2.0_wp, &
      ipeashift = 1.78069_wp, &
      gam3shell = gam3shell, &
      aesshift =1.2_wp, &
      aesexp =  4.0_wp, &
      aesrmax = 5.0_wp, &
      alphaj = 2.0_wp, &
      aesdmp3 = 3.0_wp, &
      aesdmp5 = 4.0_wp )

   type(dftd_parameter) :: gfn2Disp = dftd_parameter(&
      s6=1.0_wp, s8=2.7_wp, a1=0.52_wp, a2=5.0_wp, s9=5.0_wp)


   !> Maximum number of elements supported by GFN2-xTB
   integer, parameter :: maxElem = 86


   integer, parameter :: gfn2Kinds(118) = [&
   &  1,                                                 1, &! H-He
   &  1, 1,                               1, 1, 1, 1, 1, 1, &! Li-Ne
   &  1, 1,                               1, 1, 1, 1, 1, 1, &! Na-Ar
   &  1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, &! K-Kr
   &  1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, &! Rb-Xe
   &  1, 1, &! Cs/Ba
   &        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &!La-Lu
   &        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, &! Lu-Rn
   &  1, 1, &
   &        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &!Fr-
   &        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1 ]! -Og

   ! ========================================================================
   ! REPULSION DATA
   !> Repulsion exponent for heavy elements
   real(wp), parameter :: kExp = 1.5_wp

   !> Repulsion exponent for light elements
   real(wp), parameter :: kExpLight = 1.0_wp

   !> Repulsion exponent
   real(wp), parameter :: rExp = 1.0_wp

   !> Exponents of repulsion term
   real(wp) :: repAlpha(1:maxElem)
   !> Effective nuclear charge
   real(wp) :: repZeff(1:maxElem)
   ! ========================================================================
   ! COULOMB DATA
   !> Atomic hardnesses used in second order electrostatics
   real(wp) :: chemicalHardness(1:maxElem)
   !> Third order electrostatics is shell resolved
   logical, parameter :: thirdOrderShellResolved = .true.

   !> Third order Hubbard derivatives
   real(wp) :: thirdOrderAtom(1:maxElem)   
   !> Scaling factors for shell electrostatics
   real(wp), parameter :: shellHardness(1:3, 1:maxElem) = reshape([&
      & 0.0_wp, 0.0000000_wp, 0.0000000_wp, 0.0_wp, 0.0000000_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1972612_wp, 0.0000000_wp, 0.0_wp, 0.9658467_wp, 0.0000000_wp, &
      & 0.0_wp, 0.3994080_wp, 0.0000000_wp, 0.0_wp, 0.1056358_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1164892_wp, 0.0000000_wp, 0.0_wp, 0.1497020_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1677376_wp, 0.0000000_wp, 0.0_wp, 0.1190576_wp,-0.3200000_wp, &
      & 0.0_wp, 0.1018894_wp, 0.0000000_wp, 0.0_wp, 1.4000000_wp,-0.0500000_wp, &
      & 0.0_wp,-0.0603699_wp, 0.2000000_wp, 0.0_wp,-0.5580042_wp,-0.2300000_wp, &
      & 0.0_wp,-0.1558060_wp,-0.3500000_wp, 0.0_wp,-0.1085866_wp,-0.2500000_wp, &
      & 0.0_wp, 0.4989400_wp, 0.5000000_wp, 0.0_wp,-0.0461133_wp,-0.0100000_wp, &
      & 0.0_wp, 0.3483655_wp, 0.0000000_wp, 0.0_wp, 1.5000000_wp,-0.2500000_wp, &
      & 0.0_wp,-0.0800000_wp,-0.2046716_wp, 0.0_wp,-0.3800000_wp,-0.4921114_wp, &
      & 0.0_wp,-0.4500000_wp,-0.0379088_wp, 0.0_wp,-0.4700000_wp, 0.7405872_wp, &
      & 0.0_wp,-0.6000000_wp, 0.0545811_wp, 0.0_wp,-0.6500000_wp, 0.4046615_wp, &
      & 0.0_wp,-0.6500000_wp,-0.2418493_wp, 0.0_wp,-0.6000000_wp,-0.0611188_wp, &
      & 0.0_wp, 0.0700000_wp, 1.3333066_wp, 0.0_wp, 0.0684343_wp, 0.0000000_wp, &
      & 0.0_wp,-0.5416555_wp,-0.3000000_wp, 0.0_wp,-0.3809089_wp,-0.1500000_wp, &
      & 0.0_wp,-0.4104743_wp,-0.5000000_wp, 0.0_wp, 0.1192113_wp,-0.2500000_wp, &
      & 0.0_wp, 0.5203002_wp, 0.4000000_wp, 0.0_wp,-0.2503223_wp,-0.0700000_wp, &
      & 0.0_wp, 0.9386493_wp, 0.0000000_wp, 0.0_wp, 1.5000000_wp,-0.2500000_wp, &
      & 0.0_wp,-0.4500000_wp,-0.3349288_wp, 0.0_wp,-0.1100000_wp,-0.4422630_wp, &
      & 0.0_wp,-0.0500000_wp,-0.3562950_wp, 0.0_wp,-0.3000000_wp,-0.4301371_wp, &
      & 0.0_wp,-0.6000000_wp, 0.3956819_wp, 0.0_wp,-0.6500000_wp,-0.3052305_wp, &
      & 0.0_wp,-0.6500000_wp,-0.1881774_wp, 0.0_wp,-0.6000000_wp, 0.0931707_wp, &
      & 0.0_wp,-0.0300000_wp, 0.8024848_wp, 0.0_wp, 0.2388669_wp, 0.0000000_wp, &
      & 0.0_wp,-0.5867460_wp,-0.2800000_wp, 0.0_wp,-0.5090746_wp,-0.0600000_wp, &
      & 0.0_wp,-0.6278501_wp,-0.5500000_wp, 0.0_wp,-0.1555334_wp, 0.0600000_wp, &
      & 0.0_wp,-0.0338735_wp, 0.3000000_wp, 0.0_wp,-0.2302667_wp,-0.2300000_wp, &
      & 0.0_wp, 0.2494305_wp, 0.0000000_wp, 0.0_wp, 2.2247532_wp,-0.2300000_wp, &
      & 0.0_wp,-0.3000000_wp,-0.4699666_wp, 0.0_wp,-0.3000000_wp,-0.5539659_wp, &
      & 0.0_wp,-0.2769230_wp,-0.5462784_wp, 0.0_wp,-0.2538460_wp,-0.5385909_wp, &
      & 0.0_wp,-0.2307691_wp,-0.5309034_wp, 0.0_wp,-0.2076921_wp,-0.5232158_wp, &
      & 0.0_wp,-0.1846151_wp,-0.5155283_wp, 0.0_wp,-0.1615381_wp,-0.5078408_wp, &
      & 0.0_wp,-0.1384612_wp,-0.5001533_wp, 0.0_wp,-0.1153842_wp,-0.4924658_wp, &
      & 0.0_wp,-0.0923072_wp,-0.4847782_wp, 0.0_wp,-0.0692302_wp,-0.4770907_wp, &
      & 0.0_wp,-0.0461533_wp,-0.4694032_wp, 0.0_wp,-0.0230763_wp,-0.4617157_wp, &
      & 0.0_wp, 0.0000007_wp,-0.4540282_wp, 0.0_wp, 0.1000000_wp,-0.4486165_wp, &
      & 0.0_wp, 0.0500000_wp,-0.3394380_wp, 0.0_wp, 0.3700000_wp,-0.3419199_wp, &
      & 0.0_wp,-0.6000000_wp, 0.6586864_wp, 0.0_wp,-0.6500000_wp, 0.1350223_wp, &
      & 0.0_wp,-0.6500000_wp,-0.0977957_wp, 0.0_wp,-0.6000000_wp,-0.0203212_wp, &
      & 0.0_wp,-0.6000000_wp, 0.0614126_wp, 0.0_wp,-0.5375121_wp, 0.0000000_wp, &
      & 0.0_wp,-0.7133401_wp, 0.0000000_wp, 0.0_wp, 0.7838251_wp, 0.0000000_wp, &
      & 0.0_wp,-0.6000000_wp, 0.0000000_wp, 0.0_wp,-0.8109155_wp, 0.0000000_wp, &
      & 0.0_wp,-0.2532073_wp, 0.2500000_wp, 0.0_wp,-0.0302388_wp,-0.2300000_wp],&
      & shape(shellHardness))

   ! ========================================================================
   ! MULTIPOLE DATA
   !> Damping parameter for charge-dipole interactions
   real(wp), parameter :: dipDamp = 3.0_wp

   !> Damping parameter for dipole-dipole, charge-quadrupole interactions
   real(wp), parameter :: quadDamp = 4.0_wp

   !> Coordination number shift
   real(wp), parameter :: cnShift = 1.2_wp

   !> Coordination number exponent for radii
   real(wp), parameter :: cnExp = 4.0_wp

   !> Maximum radius
   real(wp), parameter :: cnRMax = 5.0_wp

   !> Dipole exchange-correlation kernel
   real(wp) :: dipKernel(1:maxElem)
   !> Quadrupole exchange-correlation kernel
   real(wp) :: quadKernel(1:maxElem)
   ! ========================================================================
   ! HAMILTONIAN DATA
   !> Number of shells
   integer, parameter :: nShell(1:maxElem) = [&
      & 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, &
      & 2, 2, 2, 2, 3, 3]

   !> Angular momentum of each shell
   integer, parameter :: angShell(3, 1:maxElem) = reshape([&
      & 0, 0, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, &
      & 0, 1, 0,  0, 1, 0,  0, 1, 2,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 0,  0, 1, 2,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  0, 1, 0,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, &
      & 0, 1, 2,  0, 1, 2], shape(angShell))

   !> Principal quantum number of each shell
   integer, parameter :: principalQuantumNumber(3, 1:maxElem) = reshape([&
      & 1, 0, 0,  1, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0, &
      & 2, 2, 0,  2, 2, 0,  2, 2, 3,  3, 3, 0,  3, 3, 3,  3, 3, 3,  3, 3, 3, &
      & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  4, 4, 0,  4, 4, 3,  3, 4, 4, &
      & 3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4, &
      & 3, 4, 4,  4, 4, 0,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, &
      & 4, 4, 4,  5, 5, 0,  5, 5, 4,  4, 5, 5,  4, 5, 5,  4, 5, 5,  4, 5, 5, &
      & 4, 5, 5,  4, 5, 5,  4, 5, 5,  4, 5, 5,  4, 5, 5,  5, 5, 0,  5, 5, 5, &
      & 5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  6, 6, 0,  6, 6, 5, &
      & 5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6, &
      & 5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6, &
      & 5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6, &
      & 5, 6, 6,  5, 6, 6,  6, 6, 0,  6, 6, 0,  6, 6, 0,  6, 6, 0,  6, 6, 0, &
      & 6, 6, 5,  6, 6, 5], shape(principalQuantumNumber))

   !> Reference occupation of the atom
   real(wp), parameter :: referenceOcc(0:2, 1:maxElem) = reshape([&
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp,  1.0_wp, 3.0_wp, 0.0_wp, &
      & 1.5_wp, 3.5_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp, &
      & 2.0_wp, 6.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  1.5_wp, 2.5_wp, 0.0_wp,  1.5_wp, 3.5_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 2.0_wp,  1.0_wp, 1.0_wp, 3.0_wp,  1.0_wp, 1.0_wp, 4.0_wp, &
      & 1.0_wp, 1.0_wp, 5.0_wp,  1.0_wp, 1.0_wp, 6.0_wp,  1.0_wp, 1.0_wp, 7.0_wp, &
      & 1.0_wp, 1.0_wp, 8.0_wp,  1.0_wp, 0.0_wp,10.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  1.5_wp, 2.5_wp, 0.0_wp,  1.5_wp, 3.5_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 2.0_wp,  1.0_wp, 1.0_wp, 3.0_wp,  1.0_wp, 1.0_wp, 4.0_wp, &
      & 1.0_wp, 1.0_wp, 5.0_wp,  1.0_wp, 1.0_wp, 6.0_wp,  1.0_wp, 1.0_wp, 7.0_wp, &
      & 1.0_wp, 1.0_wp, 8.0_wp,  1.0_wp, 0.0_wp,10.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 2.0_wp, &
      & 1.0_wp, 1.0_wp, 3.0_wp,  1.0_wp, 1.0_wp, 4.0_wp,  1.0_wp, 1.0_wp, 5.0_wp, &
      & 1.0_wp, 1.0_wp, 6.0_wp,  1.0_wp, 1.0_wp, 7.0_wp,  1.0_wp, 1.0_wp, 8.0_wp, &
      & 1.0_wp, 0.0_wp,10.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp, &
      & 2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp, &
      & 2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp], shape(referenceOcc))

   !> Shell polynomials to scale Hamiltonian elements
   real(wp) :: shellPoly(1:4, 1:maxElem)   
   !> Coordination number dependence of the atomic levels
   real(wp) :: kCN(1:4, 1:maxElem)
   !> Atomic level information
   real(wp) :: selfEnergy(3, 1:maxElem)
   !> Exponent of the Slater function
   real(wp) :: slaterExponent(3, 1:maxElem)

   character(len=*), parameter :: quadKernel_path='gfn2ml/quadKernel.dat'
   character(len=*), parameter :: chemicalHardness_path='gfn2ml/chemicalHardness.dat'
   character(len=*), parameter :: dipKernel_path='gfn2ml/dipKernel.dat'
   character(len=*), parameter :: repAlpha_path='gfn2ml/repAlpha.dat'
   character(len=*), parameter :: repZeff_path='gfn2ml/repZeff.dat'
   character(len=*), parameter :: shellPoly_path='gfn2ml/shellPoly.dat'
   character(len=*), parameter :: slaterExponent_path='gfn2ml/slaterExponent.dat'
   character(len=*), parameter :: thirdOrderAtom_path='gfn2ml/thirdOrderAtom.dat'
   character(len=*), parameter :: selfEnergy_path='gfn2ml/selfEnergy.dat'
   character(len=*), parameter :: kCN_path='gfn2ml/KCN.dat'

contains

subroutine loadArray1D(filepath, array)
      implicit none
      character(len=*), intent(in) :: filepath
      real(wp), intent(inout) :: array(:)
      integer :: iunit, j, r
      r = size(array)
      open(newunit=iunit, file=filepath, action='read')
      do j = 1, r
         read(iunit, *) array(j)
      enddo
end subroutine loadArray1D

subroutine loadArray2D(filepath, array)
      implicit none
      character(len=*), intent(in) :: filepath
      real(wp), intent(inout) :: array(:, :)
      integer :: iunit, j, s(2)
      s = shape(array)
      open(newunit=iunit, file=filepath, action='read')
      do j = 1, s(2)
         read(iunit, *) array(:, j)
      enddo
end subroutine loadArray2D


subroutine initData(self)

   !> Data instance
   type(TxTBData), intent(out) :: self
   print*,'ALLAHOAKBAR KHOMEINI RAHBAR'
   call loadArray('gfn2ml/chemicalHardness.dat', chemicalHardness)
   call loadArray('gfn2ml/dipKernel.dat', dipKernel) 
   call loadArray('gfn2ml/quadKernel.dat', quadKernel) 
   call loadArray('gfn2ml/repAlpha.dat', repAlpha) 
   call loadArray('gfn2ml/repZeff.dat', repZeff) 
   call loadArray('gfn2ml/shellPoly.dat', shellPoly) 
   call loadArray('gfn2ml/slaterExponent.dat', slaterExponent) 
   call loadArray('gfn2ml/thirdOrderAtom.dat', thirdOrderAtom) 
   call loadArray('gfn2ml/selfEnergy.dat', selfEnergy) 
   call loadArray('gfn2ml/KCN.dat', kCN) 

   print*,   chemicalHardness(1)
   print*,   chemicalHardness(8)
   print*,   dipKernel(1)
   print*,   dipKernel(8)
   print*,   KCN(:, 1)
   print*,   KCN(:, 8)
   print*,   quadKernel(1)
   print*,   quadKernel(8)
   print*,   repAlpha(1)
   print*,   repAlpha(8)
   print*,   repZeff(1)
   print*,   repZeff(8)
   print*,   selfEnergy(:, 1)
   print*,   selfEnergy(:, 8)
   print*,   shellPoly(:, 1)
   print*,   shellPoly(:, 8)
   print*,   slaterExponent(:, 1)
   print*,   slaterExponent(:, 8)
   print*,   thirdOrderAtom(1)
   print*,   thirdOrderAtom(8)

   self%name = 'GFN2-xTB_ML'
   self%doi = 'NULL'
   self%level = 2
   self%nShell = nShell(:maxElem)
   self%ipeashift = gfn2Globals%ipeashift * 0.1_wp

   call initGFN2ML(self%repulsion)
   call initGFN2ML(self%dispersion)
   call initGFN2ML(self%coulomb, self%nShell)
   allocate(self%multipole)
   call initGFN2ML(self%multipole)
   call initGFN2ML(self%hamiltonian, self%nShell)

end subroutine initData


subroutine initRepulsion(self)

   !> Data instance
   type(TRepulsionData), intent(out) :: self

   call init(self, kExp, kExpLight, rExp, 0.0_wp, repAlpha, repZeff)

end subroutine initRepulsion


subroutine initDispersion(self)

   !> Data instance
   type(TDispersionData), intent(out) :: self

   self%dpar = gfn2Disp
   self%g_a = 3.0_wp
   self%g_c = 2.0_wp
   self%wf  = 6.0_wp

   call newD4Model(self%dispm,self%g_a, self%g_c, p_refq_gfn2xtb)

end subroutine initDispersion


subroutine initCoulomb(self, nShell)

   !> Data instance
   type(TCoulombData), intent(out) :: self

   !> Number of shells
   integer, intent(in) :: nShell(:)

   self%gExp = gfn2Globals%alphaj
   self%chemicalHardness = chemicalHardness
   self%thirdOrderAtom = thirdOrderAtom
   allocate(self%thirdOrderShell(maxval(nShell), size(nShell)))
   call setGFN2ThirdOrderShell(self%thirdOrderShell, nShell, angShell, &
      & thirdOrderAtom, gfn2Globals%gam3shell)
   allocate(self%shellHardness(maxval(nShell), maxElem))
   call setGFN1ShellHardness(self%shellHardness, nShell, angShell, &
      & chemicalHardness, shellHardness)

end subroutine initCoulomb


subroutine initMultipole(self)

   !> Data instance
   type(TMultipoleData), intent(out) :: self

   call init(self, cnShift, cnExp, cnRMax, dipDamp, quadDamp, &
      & dipKernel, quadKernel)

end subroutine initMultipole


subroutine initHamiltonian(self, nShell)

   !> Data instance
   type(THamiltonianData), intent(out) :: self

   !> Number of shells
   integer, intent(in) :: nShell(:)

   integer :: mShell, nPrim, lAng
   integer :: iZp, iSh, jSh
   logical :: valShell(0:3)

   mShell = maxval(nShell)
   self%angShell = angShell(:mShell, :maxElem)

   do iSh = 0, 3
      do jSh = 0, 3
         self%kScale(jSh, iSh) = 0.5_wp * (gfn2Globals%kShell(iSh) &
            & + gfn2Globals%kShell(jSh))
      end do
   end do
   self%kScale(0,2) = gfn2Globals%ksd
   self%kScale(2,0) = gfn2Globals%ksd
   self%kScale(1,2) = gfn2Globals%kpd
   self%kScale(2,1) = gfn2Globals%kpd
   self%kDiff = gfn2Globals%kDiff
   do iSh = 0, 3
      do jSh = 0, 3
         self%enScale(jSh, iSh) = 0.005_wp * (gfn2Globals%enshell(iSh) &
            & + gfn2Globals%enshell(jSh))
      end do
   end do
   self%enScale4 = gfn2Globals%enscale4
   self%wExp = 0.5_wp

   self%electronegativity = paulingEN(:maxElem)
   self%atomicRad = atomicRad(:maxElem)
   self%shellPoly = shellPoly(:, :maxElem)
   self%selfEnergy = selfEnergy(:mShell, :maxElem)
   self%slaterExponent = slaterExponent(:mShell, :maxElem)
   self%principalQuantumNumber = principalQuantumNumber(:mShell, :maxElem)

   allocate(self%kCN(mShell, maxElem))
   call angToShellData(self%kCN, nShell, self%angShell, kCN)

   allocate(self%pairParam(maxElem, maxElem))
   self%pairParam = 1.0_wp

   allocate(self%valenceShell(mShell, maxElem))
   call generateValenceShellData(self%valenceShell, nShell, self%angShell)

   allocate(self%referenceOcc(mShell, maxElem))
   call setGFN2ReferenceOcc(self, nShell)

   allocate(self%numberOfPrimitives(mShell, maxElem))
   call setGFN2NumberOfPrimitives(self, nShell)

end subroutine initHamiltonian


subroutine setGFN2ReferenceOcc(self, nShell)

   !> Data instance
   type(THamiltonianData), intent(inout) :: self

   !> Number of shells
   integer, intent(in) :: nShell(:)

   integer :: lAng, iZp, iSh
   logical :: valShell(0:3)

   self%referenceOcc(:, :) = 0.0_wp
   do iZp = 1, maxElem
      do iSh = 1, nShell(iZp)
         lAng = self%angShell(iSh, iZp)
         if (self%valenceShell(iSh, iZp) /= 0) then
            self%referenceOcc(iSh, iZp) = referenceOcc(lAng, iZp)
         end if
      end do
   end do

end subroutine setGFN2ReferenceOcc


subroutine setGFN2NumberOfPrimitives(self, nShell)

   !> Data instance
   type(THamiltonianData), intent(inout) :: self

   !> Number of shells
   integer, intent(in) :: nShell(:)

   integer :: nPrim, iZp, iSh

   do iZp = 1, maxElem
      do iSh = 1, nShell(iZp)
         nPrim = 0
         if (iZp <= 2) then
            select case(self%angShell(iSh, iZp))
            case(0)
               nPrim = 3
            case(1)
               nPrim = 4
            end select
         else
            select case(self%angShell(iSh, iZp))
            case(0)
               if (self%principalQuantumNumber(iSh, iZp) > 5) then
                  nPrim = 6
               else
                  nPrim = 4
               end if
            case(1)
               if (self%principalQuantumNumber(iSh, iZp) > 5) then
                  nPrim = 6
               else
                  nPrim = 4
               end if
            case(2)
               nPrim = 3
            case(3)
               nPrim = 4
            end select
         end if
         self%numberOfPrimitives(iSh, iZp) = nPrim
      end do
   end do

end subroutine setGFN2NumberOfPrimitives


subroutine setGFN2ThirdOrderShell(thirdOrderShell, nShell, angShell, &
      & thirdOrderAtom, gam3Shell)

   real(wp), intent(out) :: thirdOrderShell(:, :)

   integer, intent(in) :: nShell(:)

   integer, intent(in) :: angShell(:, :)

   real(wp), intent(in) :: thirdOrderAtom(:)

   real(wp), intent(in) :: gam3Shell(:, 0:)

   integer :: nElem, iZp, iSh, lAng, iKind

   nElem = min(size(thirdOrderShell, dim=2), size(nShell), size(angShell, dim=2), &
      & size(thirdOrderAtom))

   thirdOrderShell(:, :) = 0.0_wp
   do iZp = 1, maxElem
      iKind = gfn2Kinds(iZp)
      do iSh = 1, nShell(iZp)
         lAng = angShell(iSh, iZp)
         thirdOrderShell(iSh, iZp) = thirdOrderAtom(iZp) * gam3Shell(iKind, lAng)
      end do
   end do

end subroutine setGFN2ThirdOrderShell


end module xtb_xtb_gfn2ml
