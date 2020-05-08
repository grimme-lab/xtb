module gff_mdcom 
      use iso_fortran_env, only : wp => real64
      implicit none
      private :: wp
      public

      real(wp) :: tmax = 100.
      real(wp) :: tsoll= 298.
      real(wp) :: tstep0=2.5 
      real(wp) :: dumpmd=100.
      real(wp) :: hmass =4. 
      logical  :: thermostat = .true. 

end
