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

module gff_gffcom
      use iso_fortran_env, only : wp => real64
      implicit none
      private :: wp
      public

!common variables which are used in energy-gradient routines
      real(wp) :: repscaln,repscalb  ! repulsion scaling
      real(wp) :: atcuta,atcutt      ! bend/tors angle damping
      real(wp) :: hbacut,hbscut      ! damping HB
      real(wp) :: xbacut,xbscut      ! damping XB
      real(wp) :: hbalp,hblongcut    ! damping HB/XB
      real(wp) :: hbst,hbsf,xbst,xbsf! charge scaling HB/XB
      real(wp) :: xhaci_globabh      ! HB AH-B
      real(wp) :: hbabmix            ! HB AH-B
      real(wp) :: hbnbcut            ! new parameter for neighbour angle
      real(wp) :: cnmax              ! max CN cut-off              
      real(wp) :: efield(3)          ! electric field components   
!general common stuff used in energy-gradient routines
      real(wp) :: repa (86),repz (86),zb3atm(86)   ! rep alpha bond, prefactor (Zval), 3atm bond term
      real(wp) :: xhaci(86),xhbas(86),xbaci(86)    ! HB/XB
      real(wp) :: rad(86)                          ! radii used for HB/XB damping and topology 
      real(wp) :: en (86)                          ! EN
      real(wp) :: cnf(86)                          ! EN dep. in EEQ. xi,gam,alp set for atoms in gnff_ini
      real(wp) :: d3r0(86*87/2)                    ! BJ radii set in gnff_ini()
!numerical precision cut-offs
      real(wp) :: cnthr,repthr,dispthr,hbthr1,hbthr2,accff

      integer  :: group(86),metal(86),normcn(86)   ! for assignment
      integer  :: ffmode
!number of terms
      integer  :: nbond,nangl,ntors,nhb1,nhb2,nxb,nathbH,nathbAB,natxbAB,nbatm
      integer  :: nfrag
!file type read
      integer  :: read_file_type
!lists
      integer,allocatable ::     nb(:,:)   ! neighbors nb(20,i) is the # neigbors
      integer,allocatable ::    bpair(:)   ! # of cov. between atoms
      integer,allocatable ::  blist(:,:)   ! bonded atoms
      integer,allocatable ::  alist(:,:)   ! angles
      integer,allocatable ::  tlist(:,:)   ! torsions
      integer,allocatable :: b3list(:,:)   ! bond atm   
      integer,allocatable :: hblist1(:,:)  ! HBs loose
      integer,allocatable :: hblist2(:,:)  ! HBs bonded
      integer,allocatable :: hblist3(:,:)  ! XBs
      integer,allocatable ::hbatABl(:,:)   ! AB atoms for HB
      integer,allocatable ::xbatABl(:,:)   ! AB atoms for XB
      integer,allocatable :: hbatHl (:)    ! H  atoms for HB
      integer,allocatable :: fraglist(:)   ! atoms in molecular fragments (for EEQ)
      integer,allocatable :: qpdb  (:)     ! atomic charge in residues from PDB file
!potential parameters used in energy-gradient routine
      real(wp),allocatable:: vbond(:,:)   ! bonds
      real(wp),allocatable:: vangl(:,:)   ! angles
      real(wp),allocatable:: vtors(:,:)   ! torsions
      real(wp),allocatable:: chieeq(:)    ! atomic ENs for EEQ
      real(wp),allocatable:: gameeq(:)    ! atomic gamma for EEQ
      real(wp),allocatable:: alpeeq(:)    ! atomic alpha for EEQ, squared
      real(wp),allocatable:: alphanb(:)   ! non-bonded exponent for atom pairs
      real(wp),allocatable::    qa(:)     ! estimated atomic charges (fixed and obtained from topology EEQ)
      real(wp),allocatable::     q(:)     ! atomic charges (obtained from EEQ)
      real(wp),allocatable:: hbrefgeo(:,:)! atom xyz, used to check for HB list update       
      real(wp),allocatable::    xyze0(:,:)! atom xyz, starting geom. (for Efield energy)     
      real(wp),allocatable:: zetac6(:)    ! D4 scaling factor product 
      real(wp),allocatable:: qfrag (:)    ! fragment charge (for EEQ)

     data xhaci / 86 * 0 /
     data xhbas / 86 * 0 /
     data xbaci / 86 * 0 /

! Pauling EN
      data en/ 2.200,3.000,0.980,1.570,2.040,2.550,3.040,3.440,3.980 &
     &        ,4.500,0.930,1.310,1.610,1.900,2.190,2.580,3.160,3.500 &
     &        ,0.820,1.000,1.360,1.540,1.630,1.660,1.550,1.830,1.880 &
     &        ,1.910,1.900,1.650,1.810,2.010,2.180,2.550,2.960,3.000 &
     &        ,0.820,0.950,1.220,1.330,1.600,2.160,1.900,2.200,2.280 &
     &        ,2.200,1.930,1.690,1.780,1.960,2.050,2.100,2.660,2.600 &
     &,0.79,0.89,1.10,1.12,1.13,1.14,1.15,1.17,1.18,1.20,1.21,1.22 &
     &,1.23,1.24,1.25,1.26,1.27,1.3,1.5,1.7,1.9,2.1,2.2,2.2,2.2 &   ! value of W-Au modified
!org &,1.23,1.24,1.25,1.26,1.27,1.3,1.5,2.36,1.9,2.2,2.20,2.28,2.54 &
     &,2.00,1.62,2.33,2.02,2.0,2.2,2.2/
! COVALENT RADII, used only in neighbor list determination
! based on "Atomic Radii of the Elements," M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
! in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
! edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
! corrected Nov. 17, 2010 for the 92nd edition.
      data rad /&
     &0.32D0,0.37D0,1.30D0,0.99D0,0.84D0,0.75D0,0.71D0,0.64D0,0.60D0,&
     &0.62D0,1.60D0,1.40D0,1.24D0,1.14D0,1.09D0,1.04D0,1.00D0,1.01D0,&
     &2.00D0,1.74D0,1.59D0,1.48D0,1.44D0,1.30D0,1.29D0,1.24D0,1.18D0,&
     &1.17D0,1.22D0,1.20D0,1.23D0,1.20D0,1.20D0,1.18D0,1.17D0,1.16D0,&
     &2.15D0,1.90D0,1.76D0,1.64D0,1.56D0,1.46D0,1.38D0,1.36D0,1.34D0,&
     &1.30D0,1.36D0,1.40D0,1.42D0,1.40D0,1.40D0,1.37D0,1.36D0,1.36D0,&
     &2.38D0,2.06D0,1.94D0,1.84D0,1.90D0,1.88D0,1.86D0,1.85D0,1.83D0,&
     &1.82D0,1.81D0,1.80D0,1.79D0,1.77D0,1.77D0,1.78D0,1.74D0,1.64D0,&
     &1.58D0,1.50D0,1.41D0,1.36D0,1.32D0,1.30D0,1.30D0,1.32D0,1.44D0,&
     &1.45D0,1.50D0,1.42D0,1.48D0,1.46D0/

      data metal / &
     &0,                                                          0,&!He
     &1,1,                                         0, 0, 0, 0, 0, 0,&!Ne
     &1,1,                                         1, 0, 0, 0, 0, 0,&!Ar
     &1,1,2,          2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 0, 0, 0, 0, 0,&!Kr
     &1,2,2,          2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 1, 0, 0, 0, 0,&!Xe
     &1,2,2,  14*2,   2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 1, 1, 1, 0, 0/ !Rn  ! At is NOT a metal, Po is borderline but slightly better as metal
      data group / &
     &1,                                                          8,&!He
     &1,2,                                         3, 4, 5, 6, 7, 8,&!Ne
     &1,2,                                         3, 4, 5, 6, 7, 8,&!Ar
     &1,2,-3,        -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8,&!Kr
     &1,2,-3,        -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8,&!Xe
     &1,2,-3, 14*-3, -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8/ !Rn
      data normcn/ &  ! only for non metals well defined
     &1,                                                          0,&!He
     &4,4,                                         4, 4, 4, 2, 1, 0,&!Ne
     &4,4,                                         4, 4, 4, 2, 1, 0,&!Ar
     &4,4,4,          4, 6, 6, 6, 6, 6, 6, 4, 4,   4, 4, 4, 4, 1, 0,&!Kr
     &4,4,4,          4, 6, 6, 6, 6, 6, 6, 4, 4,   4, 4, 4, 4, 1, 0,&!Xe
     &4,4,4,  14*4,   4, 6, 6, 6, 6, 6, 6, 6, 4,   4, 4, 4, 4, 1, 0/ !Rn
      data repz  / &
     &1.,                                                         2.,&!He
     &1.,2.,                                       3.,4.,5.,6.,7.,8.,&!Ne
     &1.,2.,                                       3.,4.,5.,6.,7.,8.,&!Ar
     &1.,2.,3.,      4.,5.,6.,7.,8.,9.,10.,11.,12.,3.,4.,5.,6.,7.,8.,&!Kr
     &1.,2.,3.,      4.,5.,6.,7.,8.,9.,10.,11.,12.,3.,4.,5.,6.,7.,8.,&!Xe
     &1.,2.,3.,14*3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,3.,4.,5.,6.,7.,8./ !Rn

     !Mass
     real(wp)              :: ams(107)
     data  ams /  1.00790d0,  4.00260d0,  6.94000d0,  9.01218d0,&
     &10.81000d0, 12.01100d0, 14.00670d0, 15.99940d0, 18.99840d0,&
     &20.17900d0, 22.98977d0, 24.30500d0, 26.98154d0, 28.08550d0,&
     &30.97376d0, 32.06000d0, 35.45300d0, 39.94800d0, 39.09830d0,&
     &40.08000d0, 44.95590d0, 47.90000d0, 50.94150d0, 51.99600d0,&
     &54.93800d0, 55.84700d0, 58.93320d0, 58.71000d0, 63.54600d0,&
     &65.38000d0, 69.73500d0, 72.59000d0, 74.92160d0, 78.96000d0,&
     &79.90400d0, 83.80000d0, 85.46780d0, 87.62000d0, 88.90590d0,&
     &91.22000d0, 92.90640d0, 95.94000d0, 98.90620d0, 101.0700d0,&
     &102.9055d0, 106.4000d0, 107.8680d0, 112.4100d0, 114.8200d0,&
     &118.6900d0, 121.7500d0, 127.6000d0, 126.9045d0, 131.3000d0,&
     &132.9054d0, 137.3300d0,&
     &138.91,140.12,140.91,144.24,147.00,150.36,151.97,157.25,&
     &158.93,162.50,164.93,167.26,168.93,173.04,174.97,&
     &178.4900d0, 180.9479d0,&
     &183.8500d0, 186.2070d0, 190.2000d0, 192.2200d0, 195.0900d0,&
     &196.9665d0, 200.5900d0, 204.3700d0, 207.2000d0, 208.9804d0,&
     &209.,210.,222.,21*0.000d0/

end module gff_gffcom
