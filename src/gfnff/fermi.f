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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c fermi smearing      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE FERMISMEAR(prt,NORBS,NEL,T,eig,occ,fod,e_fermi,s)
      IMPLICIT NONE                        
      integer norbs 
      integer nel   
      real*8  eig(norbs)
      real*8  occ(norbs)
      real*8  t
      real*8  fod
      real*8  e_fermi
      LOGICAL PRT

      real*8 boltz,bkt,occt,total_number,thr
      real*8 total_dfermi,dfermifunct,fermifunct,s,change_fermi

      PARAMETER (BOLTZ = 3.166808578545117E-06*27.2113957)
      PARAMETER (thr   = 1.D-9)
      integer ncycle,i,j,m,k,i1,i2

      if(nel+1.gt.norbs) return

      BKT = BOLTZ*T
 
      E_FERMI = 0.5*(EIG(NEL)+EIG(NEL+1))
      OCCT=NEL

      NCYCLE = 0
 10     TOTAL_NUMBER = 0.0
        TOTAL_DFERMI = 0.0
        NCYCLE = NCYCLE+1
        DO I = 1, NORBS
          FERMIFUNCT = 0.0
          if((EIG(I)-E_FERMI)/BKT.lt.50) then 
            FERMIFUNCT = 1.0/(EXP((EIG(I)-E_FERMI)/BKT)+1.0)
            DFERMIFUNCT = EXP((EIG(I)-E_FERMI)/BKT) /                
     .      (BKT*(EXP((EIG(I)-E_FERMI)/BKT)+1.0)**2)
          ELSE
            DFERMIFUNCT = 0.0
          END IF
          OCC(I) = FERMIFUNCT
          TOTAL_NUMBER = TOTAL_NUMBER + FERMIFUNCT
          TOTAL_DFERMI = TOTAL_DFERMI + DFERMIFUNCT
        END DO
        CHANGE_FERMI = (OCCT-TOTAL_NUMBER)/TOTAL_DFERMI
        E_FERMI = E_FERMI+CHANGE_FERMI
      IF(ABS(OCCT-TOTAL_NUMBER).GT.thr.AND.NCYCLE.LT.200) GOTO 10

      fod=0
      s  =0
      do i=1,norbs
      if(occ(i).gt.thr.and.1.0D00-OCC(I).gt.thr)
     .   S=S+OCC(I)*LOG(OCC(I))+(1.0d0-OCC(I))*LOG(1.0D00-OCC(I))
         if(eig(i).lt.e_fermi)then
            fod=fod+1.0d0-occ(i)
         else
            fod=fod+      occ(i)
         endif
      enddo
      s=s*3.166808578545117E-06*t

      IF(PRT)THEN
      WRITE(*,'('' T,E(Fermi),NFOD : '',2F10.3,F10.6)') T,E_FERMI,fod
      ENDIF

      RETURN
      END

