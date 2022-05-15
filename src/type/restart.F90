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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

!> This data types wraps common restart data for calculator, the main purpose
!  of this data type is to keep the calculator free from data that is changed
!  in the course of the calculation.
module xtb_type_restart
#if WITH_TBLITE
   use tblite_wavefunction, only : wavefunction_type
#endif
   use xtb_type_wavefunction, only : TWavefunction
   use xtb_gfnff_neighbourlist, only : TGFFNeighbourList
   implicit none
   private

   public :: TRestart


   !> Restart wrapper type
   type :: TRestart
      !> Tight binding wavefunction
      type(TWavefunction) :: wfn
      !> Force field topology
      type(TGFFNeighbourList) :: nlist
#if WITH_TBLITE
      !> TBLite cache
      type(wavefunction_type) :: tblite
#endif
   end type TRestart

contains


end module xtb_type_restart
