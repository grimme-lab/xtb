!> abstract calculator that hides implementation details
!  from calling codes
module tbdef_calculator
   use iso_fortran_env, only: wp => real64
   use tbdef_basisset
   use tbdef_param
   implicit none

   public :: tb_calculator
   private

   type :: tb_calculator
      type(tb_basisset), allocatable :: basis
      type(scc_parameter), allocatable :: param
   end type tb_calculator

end module tbdef_calculator
