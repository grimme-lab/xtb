! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

module mctcpar_chemical_hardnesses
   use xtb_mctc_accuracy, only : wp
   implicit none

   integer, private, parameter :: max_elem = 118
!  Semiempirical Evaluation of the GlobalHardness of the Atoms of 103
!  Elements of the Periodic Table Using the Most Probable Radii as
!  their Size Descriptors DULAL C. GHOSH, NAZMUL ISLAM 2009 in
!  Wiley InterScience (www.interscience.wiley.com).
!  DOI 10.1002/qua.22202
!  values in the paper multiplied by two because
!  (ii:ii)=(IP-EA)=d^2 E/dN^2 but the hardness
!  definition they use is 1/2d^2 E/dN^2 (in Eh)
   real(wp),public, parameter :: chemical_hardness(max_elem) = [ &
  &0.47259288_Wp,0.92203391_Wp,0.17452888_Wp,0.25700733_Wp,0.33949086_wp,0.42195412_wp, & ! H-C
  &0.50438193_Wp,0.58691863_Wp,0.66931351_Wp,0.75191607_Wp,0.17964105_wp,0.22157276_wp, & ! N-Mg
  &0.26348578_Wp,0.30539645_Wp,0.34734014_Wp,0.38924725_Wp,0.43115670_wp,0.47308269_wp, & ! Al-Ar
  &0.17105469_Wp,0.20276244_Wp,0.21007322_Wp,0.21739647_Wp,0.22471039_wp,0.23201501_wp, & ! Ca-Cr
  &0.23933969_Wp,0.24665638_Wp,0.25398255_Wp,0.26128863_Wp,0.26859476_wp,0.27592565_wp, & ! Mn-Zn
  &0.30762999_Wp,0.33931580_Wp,0.37235985_Wp,0.40273549_Wp,0.43445776_wp,0.46611708_wp, & ! Ga-Kr
  &0.15585079_Wp,0.18649324_Wp,0.19356210_Wp,0.20063311_Wp,0.20770522_wp,0.21477254_wp, & ! Rb-Mo
  &0.22184614_Wp,0.22891872_Wp,0.23598621_Wp,0.24305612_Wp,0.25013018_wp,0.25719937_wp, & ! Tc-Cd
  &0.28784780_Wp,0.31848673_Wp,0.34912431_Wp,0.37976593_Wp,0.41040808_wp,0.44105777_wp, & ! In-Xe
  &0.05019332_Wp,0.06762570_Wp,0.08504445_Wp,0.10247736_Wp,0.11991105_wp,0.13732772_wp, & ! Cs-Nd
  &0.15476297_Wp,0.17218265_Wp,0.18961288_Wp,0.20704760_Wp,0.22446752_wp,0.24189645_wp, & ! Pm-Dy
  &0.25932503_Wp,0.27676094_Wp,0.29418231_Wp,0.31159587_Wp,0.32902274_wp,0.34592298_wp, & ! Ho-Hf
  &0.36388048_Wp,0.38130586_Wp,0.39877476_Wp,0.41614298_Wp,0.43364510_wp,0.45104014_wp, & ! Ta-Pt
  &0.46848986_Wp,0.48584550_Wp,0.12526730_Wp,0.14268677_Wp,0.16011615_wp,0.17755889_wp, & ! Au-Po
  &0.19497557_Wp,0.21240778_Wp,0.07263525_Wp,0.09422158_Wp,0.09920295_wp,0.10418621_wp, & ! At-Th
  &0.14235633_Wp,0.16394294_Wp,0.18551941_Wp,0.22370139_Wp,0.00000000_wp,0.00000000_wp, & ! Pa-Cm
  &0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_wp,0.00000000_wp, & ! Bk-No
  &0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_wp,0.00000000_wp, & ! Rf-Mt
  &0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_wp,0.00000000_wp, & ! Ds-Mc
  &0.00000000_Wp,0.00000000_Wp,0.00000000_Wp,0.00000000_Wp ] ! Lv-Og
end module mctcpar_chemical_hardnesses
