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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! computes CM5 charges in ccm5(nat) from Hirshfeld charges
! in array chir(nat)
! nat      : # of atoms
! iz(nat)  : ordinal numbers
! echo     : logical for printout
! q(3,nat) : Cartesian coord. in Bohr
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine docm5(nat,iz,echo,q,chir,ccm5)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in)    :: nat
      integer, intent(in)    :: iz(nat)
      real(wp),intent(inout) :: q(3,nat)
      real(wp),intent(in)    :: chir(nat)
      real(wp),allocatable   :: ccm5a(:)
      logical, intent(in)    :: echo
      real(wp),allocatable   :: dccm5(:,:,:)
      real(wp) :: ccm5,dhirx,dhiry,dhirz,dcm5x,dcm5y,dcm5z,dcm5,dhir
      real(wp),parameter :: bohr=0.52917726_wp
      integer  :: i

      allocate( dccm5(3,nat,nat),ccm5a(nat), source = 0.0_wp )

      q = q * Bohr
!
      call cm5mod(nat,iz,chir,q,ccm5,dhirx,dhiry,dhirz,
     & dcm5x,dcm5y,dcm5z,ccm5a,dccm5)
      dcm5=dsqrt(dcm5x**2+dcm5y**2+dcm5z**2)
      dhir=dsqrt(dhirx**2+dhiry**2+dhirz**2)
!
      if(echo) then
      write (*,'(/,a,/,a)')
     & ' Charges (in a.u.) from CM5PAC (June 22, 2013)',
     & ' -----------------------------------------------'
      write (*,'(a,/,a)')
     & ' Center     Atomic      CM5         Hirshfeld',            
     & ' Number     Number      Charge      Charge'                 
      write (*,'(a)')
     & ' -----------------------------------------------'
        do i=1,nat    
        write (*,'(i5,6x,i5,5x,f11.6,x,f11.6)') i,iz(i),Ccm5(i),chir(i)
        enddo 
      write (*,'(a)')
     & ' -----------------------------------------------'
!
      write (*,'(/,a,/,a)')
     & ' Dipole moment (in Debye)',
     & ' -----------------------------------------------'
      write (*,'(a,/,a)') 
     & '                 X        Y        Z     Total', 
     & ' -----------------------------------------------'
      write (*,'(a,4f9.4)') ' cm5       ',dcm5x,dcm5y,dcm5z,dcm5
      write (*,'(a,4f9.4)') ' hirshfeld ',dhirx,dhiry,dhirz,dhir
      write (*,'(a,/)')
     & ' -----------------------------------------------'
      endif
 
      q = q / bohr
  
      dccm5 = dccm5 * bohr
      
      return
      end subroutine docm5

      subroutine prepcm5(nat,iz,q,chir,ccm5,ccm5a,dccm5)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(3,nat),IZ(nat),CHIR(nat)
      DIMENSION CCM5(nat),CCM5A(nat),DCCM5(3,nat,nat)
      logical echo
      parameter (bohr=0.52917726)

      q = q * Bohr
!
      CALL CM5MOD(NAT,IZ,CHIR,Q,CCM5,DHIRX,DHIRY,DHIRZ,
     $ DCM5X,DCM5Y,DCM5Z,CCM5A,DCCM5)
!
      q = q / Bohr
  
      DCCM5 = DCCM5 * Bohr
      
      RETURN
      END

      SUBROUTINE CM5MOD(NAT,IZ,CHIR,Q,CCM5,DHIRX,DHIRY,DHIRZ,
     & DCM5X,DCM5Y,DCM5Z,CCM5A,DCCM5)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in)  :: nat
      integer,PARAMETER    :: MZ=118
      real(wp),intent(in)  :: Q(3,nat)
      integer, intent(in)  :: IZ(nat)
      real(wp),intent(in)  :: CHIR(nat)
      real(wp),intent(out) :: CCM5(nat)
      real(wp),intent(out) :: CCM5A(nat)
      real(wp),intent(out) :: DCCM5(3,NAT,nat)
      real(wp),intent(out) :: dhirx,dhiry,dhirz
      real(wp),intent(out) :: dcm5x,dcm5y,dcm5z
      real(wp) :: RAD(MZ),A0(MZ),D(MZ,MZ),alp,q1,q2,q3
      real(wp) :: dis,bkk,bkkd,bkkda
      integer  :: i,j,k,k1,k2
!
! COVALENT RADII
!
! based on "Atomic Radii of the Elements," M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
! in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
! edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
! corrected Nov. 17, 2010 for the 92nd edition.
!
      RAD(1)=0.32D0
      RAD(2)=0.37D0
      RAD(3)=1.30D0
      RAD(4)=0.99D0
      RAD(5)=0.84D0
      RAD(6)=0.75D0
      RAD(7)=0.71D0
      RAD(8)=0.64D0
      RAD(9)=0.60D0
      RAD(10)=0.62D0
      RAD(11)=1.60D0
      RAD(12)=1.40D0
      RAD(13)=1.24D0
      RAD(14)=1.14D0
      RAD(15)=1.09D0
      RAD(16)=1.04D0
      RAD(17)=1.00D0
      RAD(18)=1.01D0
      RAD(19)=2.00D0
      RAD(20)=1.74D0
      RAD(21)=1.59D0
      RAD(22)=1.48D0
      RAD(23)=1.44D0
      RAD(24)=1.30D0
      RAD(25)=1.29D0
      RAD(26)=1.24D0
      RAD(27)=1.18D0
      RAD(28)=1.17D0
      RAD(29)=1.22D0
      RAD(30)=1.20D0
      RAD(31)=1.23D0
      RAD(32)=1.20D0
      RAD(33)=1.20D0
      RAD(34)=1.18D0
      RAD(35)=1.17D0
      RAD(36)=1.16D0
      RAD(37)=2.15D0
      RAD(38)=1.90D0
      RAD(39)=1.76D0
      RAD(40)=1.64D0
      RAD(41)=1.56D0
      RAD(42)=1.46D0
      RAD(43)=1.38D0
      RAD(44)=1.36D0
      RAD(45)=1.34D0
      RAD(46)=1.30D0
      RAD(47)=1.36D0
      RAD(48)=1.40D0
      RAD(49)=1.42D0
      RAD(50)=1.40D0
      RAD(51)=1.40D0
      RAD(52)=1.37D0
      RAD(53)=1.36D0
      RAD(54)=1.36D0
      RAD(55)=2.38D0
      RAD(56)=2.06D0
      RAD(57)=1.94D0
      RAD(58)=1.84D0
      RAD(59)=1.90D0
      RAD(60)=1.88D0
      RAD(61)=1.86D0
      RAD(62)=1.85D0
      RAD(63)=1.83D0
      RAD(64)=1.82D0
      RAD(65)=1.81D0
      RAD(66)=1.80D0
      RAD(67)=1.79D0
      RAD(68)=1.77D0
      RAD(69)=1.77D0
      RAD(70)=1.78D0
      RAD(71)=1.74D0
      RAD(72)=1.64D0
      RAD(73)=1.58D0
      RAD(74)=1.50D0
      RAD(75)=1.41D0
      RAD(76)=1.36D0
      RAD(77)=1.32D0
      RAD(78)=1.30D0
      RAD(79)=1.30D0
      RAD(80)=1.32D0
      RAD(81)=1.44D0
      RAD(82)=1.45D0
      RAD(83)=1.50D0
      RAD(84)=1.42D0
      RAD(85)=1.48D0
      RAD(86)=1.46D0
      RAD(87)=2.42D0
      RAD(88)=2.11D0
      RAD(89)=2.01D0
      RAD(90)=1.90D0
      RAD(91)=1.84D0
      RAD(92)=1.83D0
      RAD(93)=1.80D0
      RAD(94)=1.80D0
      RAD(95)=1.73D0
      RAD(96)=1.68D0
      RAD(97)=1.68D0
      RAD(98)=1.68D0
      RAD(99)=1.65D0
      RAD(100)=1.67D0
      RAD(101)=1.73D0
      RAD(102)=1.76D0
      RAD(103)=1.61D0
      RAD(104)=1.57D0
      RAD(105)=1.49D0
      RAD(106)=1.43D0
      RAD(107)=1.41D0
      RAD(108)=1.34D0
      RAD(109)=1.29D0
      RAD(110)=1.28D0
      RAD(111)=1.21D0
      RAD(112)=1.22D0
      RAD(113)=1.36D0
      RAD(114)=1.43D0
      RAD(115)=1.62D0
      RAD(116)=1.75D0
      RAD(117)=1.65D0
      RAD(118)=1.57D0
!
! CM5 MODEL PARAMETERS
!
      DO I=1,MZ
       A0(I)=0.D0
       DO J=1,MZ
        D(I,J)=0.D0
       ENDDO
      ENDDO
! ATOMWISE PARAMETERS
      A0(  1)= 0.0056
      A0(  2)=-0.1543
      A0(  3)= 0.0000
      A0(  4)= 0.0333
      A0(  5)=-0.1030
      A0(  6)=-0.0446
      A0(  7)=-0.1072
      A0(  8)=-0.0802
      A0(  9)=-0.0629
      A0( 10)=-0.1088
      A0( 11)= 0.0184
      A0( 12)= 0.0000
      A0( 13)=-0.0726
      A0( 14)=-0.0790
      A0( 15)=-0.0756
      A0( 16)=-0.0565
      A0( 17)=-0.0444
      A0( 18)=-0.0767
      A0( 19)= 0.0130
      A0( 20)= 0.0000
      A0( 21)= 0.0000
      A0( 22)= 0.0000
      A0( 23)= 0.0000
      A0( 24)= 0.0000
      A0( 25)= 0.0000
      A0( 26)= 0.0000
      A0( 27)= 0.0000
      A0( 28)= 0.0000
      A0( 29)= 0.0000
      A0( 30)= 0.0000
      A0( 31)=-0.0512
      A0( 32)=-0.0557
      A0( 33)=-0.0533
      A0( 34)=-0.0399
      A0( 35)=-0.0313
      A0( 36)=-0.0541
      A0( 37)= 0.0092
      A0( 38)= 0.0000
      A0( 39)= 0.0000
      A0( 40)= 0.0000
      A0( 41)= 0.0000
      A0( 42)= 0.0000
      A0( 43)= 0.0000
      A0( 44)= 0.0000
      A0( 45)= 0.0000
      A0( 46)= 0.0000
      A0( 47)= 0.0000
      A0( 48)= 0.0000
      A0( 49)=-0.0361
      A0( 50)=-0.0393
      A0( 51)=-0.0376
      A0( 52)=-0.0281
      A0( 53)=-0.0220
      A0( 54)=-0.0381
      A0( 55)= 0.0065
      A0( 56)= 0.0000
      A0( 57)= 0.0000
      A0( 58)= 0.0000
      A0( 59)= 0.0000
      A0( 60)= 0.0000
      A0( 61)= 0.0000
      A0( 62)= 0.0000
      A0( 63)= 0.0000
      A0( 64)= 0.0000
      A0( 65)= 0.0000
      A0( 66)= 0.0000
      A0( 67)= 0.0000
      A0( 68)= 0.0000
      A0( 69)= 0.0000
      A0( 70)= 0.0000
      A0( 71)= 0.0000
      A0( 72)= 0.0000
      A0( 73)= 0.0000
      A0( 74)= 0.0000
      A0( 75)= 0.0000
      A0( 76)= 0.0000
      A0( 77)= 0.0000
      A0( 78)= 0.0000
      A0( 79)= 0.0000
      A0( 80)= 0.0000
      A0( 81)=-0.0255
      A0( 82)=-0.0277
      A0( 83)=-0.0265
      A0( 84)=-0.0198
      A0( 85)=-0.0155
      A0( 86)=-0.0269
      A0( 87)= 0.0046
      A0( 88)= 0.0000
      A0( 89)= 0.0000
      A0( 90)= 0.0000
      A0( 91)= 0.0000
      A0( 92)= 0.0000
      A0( 93)= 0.0000
      A0( 94)= 0.0000
      A0( 95)= 0.0000
      A0( 96)= 0.0000
      A0( 97)= 0.0000
      A0( 98)= 0.0000
      A0( 99)= 0.0000
      A0(100)= 0.0000
      A0(101)= 0.0000
      A0(102)= 0.0000
      A0(103)= 0.0000
      A0(104)= 0.0000
      A0(105)= 0.0000
      A0(106)= 0.0000
      A0(107)= 0.0000
      A0(108)= 0.0000
      A0(109)= 0.0000
      A0(110)= 0.0000
      A0(111)= 0.0000
      A0(112)= 0.0000
      A0(113)=-0.0179
      A0(114)=-0.0195
      A0(115)=-0.0187
      A0(116)=-0.0140
      A0(117)=-0.0110
      A0(118)=-0.0189
!
      DO K1=1,MZ
       DO K2=K1+1,MZ
        D(K1,K2)=A0(K1)-A0(K2)
       ENDDO
      ENDDO
! PAIRWISE PARAMETERS
      D( 1, 6)= 0.0502                   
      D( 1, 7)= 0.1747              
      D( 1, 8)= 0.1671             
      D( 6, 7)= 0.0556             
      D( 6, 8)= 0.0234             
      D( 7, 8)=-0.0346              
!
      DO I=1,MZ
       DO J=I+1,MZ
        D(J,I)=-D(I,J)
       ENDDO
      ENDDO
! ALPHA
      ALP=2.4740
! C-COEFFICIENT: 0.7050   ! ALREADY INCLUDED IN A0
!
!
      DO K1=1,NAT
       DO K=1,NAT
         DCCM5(1,K,K1)=0.D0
         DCCM5(2,K,K1)=0.D0
         DCCM5(3,K,K1)=0.D0
       ENDDO
      ENDDO

      DO K=1,NAT   
       CCM5A(K)=0.D0
       DO K1=1,NAT   
        IF (IZ(K).NE.IZ(K1)) THEN
            Q1=(Q(1,K)-Q(1,K1))
            Q2=(Q(2,K)-Q(2,K1))
            Q3=(Q(3,K)-Q(3,K1))
            DIS=DSQRT(Q1*Q1+Q2*Q2+Q3*Q3)
            BKK=DEXP(-ALP*(DIS-RAD(IZ(K))-RAD(IZ(K1))))
            BKKD=BKK*D(IZ(K),IZ(K1))
            CCM5A(K)=CCM5A(K)+BKKD
            BKKDA=BKKD*ALP/DIS
            DCCM5(1,K1,K )=BKKDA*Q1
            DCCM5(2,K1,K )=BKKDA*Q2
            DCCM5(3,K1,K )=BKKDA*Q3
           ENDIF
          ENDDO
         ENDDO
         DO K = 1, NAT
          CCM5(K)=CHIR(K)+CCM5A(K)
      ENDDO
      DHIRX=0.D0
      DHIRY=0.D0
      DHIRZ=0.D0
      DCM5X=0.D0
      DCM5Y=0.D0
      DCM5Z=0.D0
      DO J=1,NAT   
       DHIRX=DHIRX+Q(1,J)*CHIR(J)*4.803242D0
       DHIRY=DHIRY+Q(2,J)*CHIR(J)*4.803242D0
       DHIRZ=DHIRZ+Q(3,J)*CHIR(J)*4.803242D0
       DCM5X=DCM5X+Q(1,J)*CCM5(J)*4.803242D0
       DCM5Y=DCM5Y+Q(2,J)*CCM5(J)*4.803242D0
       DCM5Z=DCM5Z+Q(3,J)*CCM5(J)*4.803242D0
      ENDDO

      RETURN
      END
