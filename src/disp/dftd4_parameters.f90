! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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

module xtb_disp_dftd4param
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   implicit none
   public

   real(wp) :: thopi,ootpi
   parameter ( thopi = 3._wp/pi )
   parameter ( ootpi = 0.5_wp/pi )

   integer, parameter :: p_refq_gfn2xtb          = 0
   integer, parameter :: p_refq_gasteiger        = 1
   integer, parameter :: p_refq_hirshfeld        = 2
   integer, parameter :: p_refq_periodic         = 3
   integer, parameter :: p_refq_gfn2xtb_gbsa_h2o = 4
   integer, parameter :: p_refq_goedecker        = 5

   integer, parameter :: p_mbd_none       = 0
   integer, parameter :: p_mbd_rpalike    = 1
   integer, parameter :: p_mbd_exact_atm  = 2
   integer, parameter :: p_mbd_approx_atm = 3

   integer, private, parameter :: max_elem = 118
   real(wp), parameter :: zeff(max_elem) = [ &
      &   1,                                                 2,  & ! H-He
      &   3, 4,                               5, 6, 7, 8, 9,10,  & ! Li-Ne
      &  11,12,                              13,14,15,16,17,18,  & ! Na-Ar
      &  19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,  & ! K-Kr
      &   9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,  & ! Rb-Xe
      &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Cs-Lu
      &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26, & ! Hf-Rn
      !  just copy & paste from above
      &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Fr-Lr
      &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26  ] ! Rf-Og

   integer, dimension(max_elem)      :: refn ! for D4
   real(wp),dimension(7,max_elem)    :: refq
   real(wp),dimension(7,max_elem)    :: refh
   real(wp),dimension(7,max_elem)    :: dftq,pbcq,gffq,solq,clsq
   real(wp),dimension(7,max_elem)    :: dfth,pbch,gffh,solh,clsh
   real(wp),dimension(7,max_elem)    :: hcount 
   real(wp),dimension(7,max_elem)    :: ascale
   real(wp),dimension(7,max_elem)    :: refcovcn
   real(wp),dimension(7,max_elem)    :: refcn
   integer, dimension(7,max_elem)    :: refsys 
   real(wp),dimension(23,7,max_elem) :: alphaiw
   real(wp),dimension(17)       :: secq
   real(wp),dimension(17)       :: dfts,pbcs,gffs,sols,clss
   real(wp),dimension(17)       :: sscale
   real(wp),dimension(17)       :: seccn
   real(wp),dimension(17)       :: seccnd3
   real(wp),dimension(23,17)    :: secaiw

   include 'param_ref.fh'

end module xtb_disp_dftd4param
