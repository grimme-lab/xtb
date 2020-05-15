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

!> Data for the generation of a force field topology
module xtb_gfnff_generator
   use xtb_mctc_accuracy, only : wp, sp
   implicit none
   private

   public :: TGFFGenerator


   !> Generator for the force field topology
   type TGFFGenerator

      !> when is an angle close to linear ? (GEODEP)
      !  for metals values closer to 170 (than to 160) are better
      real(wp) linthr

      !> skip torsion and bending if potential is small
      real(wp) fcthr

      !> R threshold in Angstroem for cov distance estimated used in apprx EEQ
      real(sp) tdist_thr

      !> important bond determination threshold, large values yield more 1.23
      real(wp) rthr

      !> decrease if a metal is present, larger values yield smaller CN
      real(wp) rthr2

      !> change of R0 for topo with charge qa
      !  larger values yield smaller CN for metals in particular
      real(wp) rqshrink

      !> H charge (qa) threshold for H in HB list 18
      real(wp) hqabthr

      !> AB charge (qa) threshold for AB in HB list
      !  - avoids HBs with positive atoms,
      !  - larger val. better for S30L but worse in PubChem RMSD checks
      real(wp) qabthr

      !> Parameter
      real(wp) srb1
      real(wp) srb2
      real(wp) srb3

      !> change of non-bonded rep. with q(topo)
      real(wp) qrepscal

      !> change of non-bonded rep. with CN
      real(wp) nrepscal

      !> HH repulsion
      real(wp) hhfac
      real(wp) hh13rep
      real(wp) hh14rep
      real(wp) bstren(9)

      !> bend FC change with polarity
      real(wp) qfacBEN

      !> torsion FC change with polarity
      real(wp) qfacTOR

      !> tors FC 3-ring
      real(wp) fr3

      !> tors FC 4-ring
      real(wp) fr4

      !> tors FC 5-ring
      real(wp) fr5

      !> tors FC 6-ring
      real(wp) fr6

      !> bonds
      real(wp) torsf(8)

      !> small bend corr.
      real(wp) fbs1

      !> bonded ATM scal
      real(wp) batmscal

      !> Shifts
      real(wp) mchishift

      !> gen shift
      real(wp) rabshift

      !> XH
      real(wp) rabshifth

      !> hypervalent
      real(wp) hyper_shift

      !> heavy
      real(wp) hshift3
      real(wp) hshift4
      real(wp) hshift5

      !> group 1+2 metals
      real(wp) metal1_shift

      !> TM
      real(wp) metal2_shift

      !> main group metals
      real(wp) metal3_shift

      !> eta bonded
      real(wp) eta_shift

      !> Charge Param
      real(wp) qfacbm(0:4)

      !> bond charge dependent
      real(wp) qfacbm0

      !> topo dist scaling
      real(wp) rfgoed1

      !> Hückel Param
      !  decrease Hueckel off-diag for triple bonds because they are less well conjugated 1.4
      real(wp) htriple

      !> increase pot depth depending on P
      real(wp) hueckelp2

      !> diagonal element change with qa
      real(wp) hueckelp3

      !> diagonal element relative to C
      real(wp) hdiag(17)

      !> Huckel off-diag constants
      real(wp) hoffdiag(17)

      !> iteration mixing
      real(wp) hiter

      !> diagonal qa dep.
      real(wp) hueckelp

      !> ref P value R shift
      real(wp) bzref

      !> ref P value k stretch
      real(wp) bzref2

      !> 2el diag shift
      real(wp) pilpf

      !> the Hückel iterations can diverge so take only a few steps
      real(wp) maxhiter

      !> D3 Param
      real(wp) d3a1

      !> D3
      real(wp) d3a2

      !> mixing of sp^n with sp^n-1
      real(wp) split0

      !> mixing of sp^n with sp^n-1
      real(wp) split1

      !> str ring size dep.
      real(wp) fringbo

      !> three coord. heavy eq. angle
      real(wp) aheavy3

      !> four coord. heavy eq. angle
      real(wp) aheavy4
      real(wp) bsmat(0:3,0:3)

      !> max CN cut-off
      real(wp) :: cnmax

   end type TGFFGenerator


contains


end module xtb_gfnff_generator
