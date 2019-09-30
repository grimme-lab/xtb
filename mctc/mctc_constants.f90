module mctc_constants
   use iso_fortran_env, only : wp => real64
   implicit none
   private
   real(wp),public,parameter :: pi = 3.1415926535897932384626433832795029_wp
   real(wp),public,parameter :: pi4 = 3.1415926535897932384626433832795029_wp*4._wp
!  √π
   real(wp),public,parameter ::  sqrtpi  = sqrt(pi)
!  2×π
   real(wp),public,parameter :: twopi = 2.0_wp * pi
!  4×π
   real(wp),public,parameter :: fourpi = 4.0_wp * pi
!  π/2
   real(wp),public,parameter :: pihalf = 0.5_wp * pi
!  Boltzmann constant in Eh/K
   real(wp),public,parameter :: kB = 3.166808578545117e-06_wp
!  speed of light c in vacuum in a.u.
   real(wp),public,parameter :: lightspeed = 137.0359990740_wp
end module mctc_constants
