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
! You should have received a copy of the GNU Lesser General Public Licen
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

subroutine getvdwxy(rx,ry,rz, c66, s6,r0, vdw)
   !cc Ableitung nach rx und ry
   Implicit Real*8 (a-h, o-z)
   integer k, l
   !    write(*,*) 's6:', s6
   avdw=20.0
   t1 = s6 * C66
   t2 = rx ** 2
   t3 = ry ** 2
   t4 = rz ** 2
   t5 = t2 + t3 + t4
   t6 = t5 ** 2
   t7 = t6 ** 2
   t11 = sqrt(t5)
   t12 = 0.1D1 / r0
   t16 = exp(-avdw * (t11 * t12 - 0.1D1))
   t17 = 0.1D1 + t16
   t25 = t17 ** 2
   t26 = 0.1D1 / t25
   t35 = 0.1D1 / t7
   t40 = avdw ** 2
   t41 = r0 ** 2
   t43 = t40 / t41
   t44 = t16 ** 2
   t56 = -0.48D2 * t1 / t7 / t5 / t17 * rx * ry + 0.13D2 * t1 / t11 /&
      & t7 * t26 * rx * avdw * t12 * ry * t16 - 0.2D1 * t1 * t35 / t25 /&
      &t17 * t43 * rx * t44 * ry + t1 * t35 * t26 * t43 * rx * ry * t16
   vdw=t56
end subroutine

subroutine getvdwxx(rx, ry, rz, c66, s6, r0, vdw)
   !cc Ableitung nach rx und rx
   Implicit Real*8 (a-h, o-z)
   avdw=20.0
   !      write(*,*) 's6:', s6
   t1 = s6 * C66
   t2 = rx ** 2
   t3 = ry ** 2
   t4 = rz ** 2
   t5 = t2 + t3 + t4
   t6 = t5 ** 2
   t7 = t6 ** 2
   t10 = sqrt(t5)
   t11 = 0.1D1 / r0
   t15 = exp(-avdw * (t10 * t11 - 0.1D1))
   t16 = 0.1D1 + t15
   t17 = 0.1D1 / t16
   t24 = t16 ** 2
   t25 = 0.1D1 / t24
   t29 = t11 * t15
   t33 = 0.1D1 / t7
   t41 = avdw ** 2
   t42 = r0 ** 2
   t44 = t41 / t42
   t45 = t15 ** 2
   t62 = -0.48D2 * t1 / t7 / t5 * t17 * t2 + 0.13D2 * t1 / t10 / t7 *&
      & t25 * t2 * avdw * t29 + 0.6D1 * t1 * t33 * t17 - 0.2D1 * t1 * t33&
      & / t24 / t16 * t44 * t2 * t45 - t1 / t10 / t6 / t5 * t25 * avdw *&
      &t29 + t1 * t33 * t25 * t44 * t2 * t15
   vdw=t62
end subroutine

! used in ancopt Lindh
logical function rcutoff(cart,katom,latom)
   implicit none
   Real*8 Cart(3,*),xkl,ykl,zkl,rkl2
   integer katom,latom
   rcutoff=.false.
   xkl=Cart(1,kAtom)-Cart(1,lAtom)
   ykl=Cart(2,kAtom)-Cart(2,lAtom)
   zkl=Cart(3,kAtom)-Cart(3,lAtom)
   rkl2 = xkl**2 + ykl**2 + zkl**2
   if(rkl2.gt.70.0)rcutoff=.true.
end function

! used in thermodynamic Lindh
logical function rcutoff1(cart,katom,latom)
   implicit none
   Real*8 Cart(3,*),xkl,ykl,zkl,rkl2
   integer katom,latom
   rcutoff1=.false.
   xkl=Cart(1,kAtom)-Cart(1,lAtom)
   ykl=Cart(2,kAtom)-Cart(2,lAtom)
   zkl=Cart(3,kAtom)-Cart(3,lAtom)
   rkl2 = xkl**2 + ykl**2 + zkl**2
   if(rkl2.gt.150.0)rcutoff1=.true.
end function

! used in thermodynamic Lindh
logical function rcutoff2(cart,katom,latom)
   implicit none
   Real*8 Cart(3,*),xkl,ykl,zkl,rkl2
   integer katom,latom
   rcutoff2=.false.
   xkl=Cart(1,kAtom)-Cart(1,lAtom)
   ykl=Cart(2,kAtom)-Cart(2,lAtom)
   zkl=Cart(3,kAtom)-Cart(3,lAtom)
   rkl2 = xkl**2 + ykl**2 + zkl**2
   if(rkl2.gt.40.0)rcutoff2=.true.
end function

INTEGER FUNCTION iTabRow(i)
   INTEGER :: i
   !
   iTabRow=0
   If (i.gt. 0 .and. i.le. 2) Then
      iTabRow=1
   Else If (i.gt. 2 .and. i.le.10) Then
      iTabRow=2
   Else If (i.gt.10 .and. i.le.18) Then
      iTabRow=3
   Else If (i.gt.18 .and. i.le.36) Then
      iTabRow=3
   Else If (i.gt.36 .and. i.le.54) Then
      iTabRow=3
   Else If (i.gt.54 .and. i.le.86) Then
      iTabRow=3
   Else If (i.gt.86) Then
      iTabRow=3
   End If
   !
   Return
End function

Subroutine Trsn(xyz,nCent,Tau,Bt,lWrite,lWarn,Label,dBt,ldB)
   !***********************************************************************
   !                                                                      *
   ! Reference: Molecular Vibrations, E. Bright Wilson, Jr, J. C. Decicius*
   !             nd Paul C. Cross, Sec. 4-1, Eq. 20-24                    *
   !                                                                      *
   ! R.Lindh May-June '96                                                 *
   !***********************************************************************
   Implicit Integer (i-n)
   Implicit Real*8 (a-h,o-z)
   !     include "common/real.inc"
   !comdeck real.inc $Revision: 2002.3 $
   Real*8 Zero, One, Two, Three, Four, Five, Six, Seven,&
      &       Eight, RNine, Ten, Half, Pi, SqrtP2, TwoP34,&
      &       TwoP54, One2C2
   Parameter(Zero =0.0D0, One  =1.0D0, Two=2.0D0, Three=3.0D0,&
      &          Four =4.0D0, Five =5.0D0, Six=6.0D0, Seven=7.0D0,&
      &          Eight=8.0D0, rNine=9.0D0, Ten=1.0D1, Half=0.5D0,&
      &          Pi    =3.141592653589793D0,&
      &          SqrtP2=0.8862269254527579D0,&
      &          TwoP34=0.2519794355383808D0,&
      &          TwoP54=5.914967172795612D0,&
      &          One2C2=0.2662567690426443D-04)

   Real*8 Bt(3,nCent), xyz(3,nCent), Rij(3), Eij(3), Rjk(3), Ejk(3),&
      &       Rkl(3), Ekl(3), Rijk(3), Eijk(3), dBt(3,nCent,3,nCent),&
      &       BRij(3,2), dBRij(3,2,3,2), BRjk(3,2), dBRjk(3,2,3,2),&
      &       BRkl(3,2), dBRkl(3,2,3,2), Bf2(3,3), dum(3,4,3,4),&
      &       Bf3(3,3)
   Logical lWrite, lWarn, ldB
   Character*8 Label
   !
   !     Call qEnter('Trsn')
   mCent=2
   Call Strtch(xyz(1,1),mCent,Rij1,BRij,.False.,Label,dBRij,ldB)
   Call Strtch(xyz(1,2),mCent,Rjk1,BRjk,.False.,Label,dBRjk,ldB)
   Call Strtch(xyz(1,3),mCent,Rkl1,BRkl,.False.,Label,dBRkl,ldB)
   mCent=3
   Call Bend(xyz(1,1),mCent,Fi2,Bf2,.False.,.False.,Label,Dum,&
      &          .False.)
   SinFi2=Sin(Fi2)
   CosFi2=Cos(Fi2)
   Call Bend(xyz(1,2),mCent,Fi3,Bf3,.False.,.False.,Label,Dum,&
      &          .False.)
   SinFi3=Sin(Fi3)
   CosFi3=Cos(Fi3)
   !
   !     Get the angle between the two planes, i.e. the
   !     angle between the normal vectors.
   !
   !     r123 * r234 = CosTau
   !
   CosTau = ( ( BRij(2,1)*BRjk(3,2) - BRij(3,1)*BRjk(2,2) ) *&
      &           ( BRjk(2,1)*BRkl(3,2) - BRjk(3,1)*BRkl(2,2) ) +&
      &           ( BRij(3,1)*BRjk(1,2) - BRij(1,1)*BRjk(3,2) ) *&
      &           ( BRjk(3,1)*BRkl(1,2) - BRjk(1,1)*BRkl(3,2) ) +&
      &           ( BRij(1,1)*BRjk(2,2) - BRij(2,1)*BRjk(1,2) ) *&
      &           ( BRjk(1,1)*BRkl(2,2) - BRjk(2,1)*BRkl(1,2) ) )&
      &         / (SinFi2*SinFi3)
   !
   !     For the vector product of the two vectors. This
   !     will give a vector parallell to e23. The direction
   !     relative to e23 defines the sign.
   !
   !     e123 X e234 = SinTau * e23
   !
   SinTau = ( BRij(1,2) * (BRjk(2,1)*BRkl(3,2)-BRjk(3,1)*BRkl(2,2))&
      &         + BRij(2,2) * (BRjk(3,1)*BRkl(1,2)-BRjk(1,1)*BRkl(3,2))&
      &         + BRij(3,2) * (BRjk(1,1)*BRkl(2,2)-BRjk(2,1)*BRkl(1,2)) )&
      &         / (SinFi2*SinFi3)
   !
   !     (-Pi < Tau <= Pi)
   !
   Tau = ATan2(SinTau,CosTau)
   If (Abs(Tau).eq.Pi) Tau=Pi
   !
   dTau = 180.0D+00*Tau/Pi
   dFi2 = 180.0D+00*Fi2/Pi
   dFi3 = 180.0D+00*Fi3/Pi
   If (lWarn) Then
      If (dTau.gt.177.5 .or. dTau.lt.-177.5) Then
         Write (*,*) ' Warning: dihedral angle close to'&
            &         //' end of range'
      End If
      If (dFi2.gt.177.5 .or. dFi2.lt.2.5) Then
         Write (*,*) ' Warning: bond angle close to'&
            &         //' end of range'
      End If
      If (dFi3.gt.177.5 .or. dFi3.lt.2.5) Then
         Write (*,*) ' Warning: bond angle close to'&
            &         //' end of range'
      End If
   End If
   If (LWRITE) Write (*,1) Label,dTau,Tau
   1     FORMAT(1X,A,' : Dihedral Angle=',F10.4,&
      & '/degree,',F10.4,'/rad')
   !
   !---- Compute the WDC matrix.
   !
   Do ix = 1, 3
      iy=ix+1
      If (iy.gt.3) iy=iy-3
      iz=iy+1
      If (iz.gt.3) iz=iz-3
      Bt(ix,1) = (BRij(iy,2)*BRjk(iz,2)-BRij(iz,2)*BRjk(iy,2))&
         &           / (Rij1*SinFi2**2)
      Bt(ix,4) = (BRkl(iy,1)*BRjk(iz,1)-BRkl(iz,1)*BRjk(iy,1))&
         &           / (Rkl1*SinFi3**2)
      Bt(ix,2) = -( (Rjk1-Rij1*CosFi2) * Bt(ix,1)&
         &             +         Rkl1*CosFi3  * Bt(ix,4))/Rjk1
      Bt(ix,3) = - ( Bt(ix,1)+Bt(ix,2)+Bt(ix,4))
   End Do
   !
   If (ldB) Then
      !
      !------- Compute the derivative of the WDC matrix.
      !
      Do ix = 1, 3
         iy=ix+1
         If (iy.gt.3) iy=iy-3
         iz=iy+1
         If (iz.gt.3) iz=iz-3
         Do jx = 1, ix
            jy=jx+1
            If (jy.gt.3) jy=jy-3
            jz=jy+1
            If (jz.gt.3) jz=jz-3
            !
            dBt(ix,1,jx,1) =(  dBRij(ix,1,jy,2)*BRjk(jz,2)&
               &                       - dBRij(ix,1,jz,2)*BRjk(jy,2)&
               &                       - Bt(jx,1)*(BRij(ix,1)*SinFi2**2&
               &                       + Rij1*Two*SinFi2*CosFi2*Bf2(ix,1)) )&
               &                       / (Rij1*SinFi2**2)
            dBt(ix,1,jx,2) =-( (-BRij(ix,1)*CosFi2&
               &                         + Rij1*SinFi2*Bf2(ix,1) ) * Bt(jx,1)&
               &                         + (Rjk1-Rij1*CosFi2) * dBt(ix,1,jx,1) )&
               &                       / Rjk1
            dBt(jx,2,ix,1) = dBt(ix,1,jx,2)
            dBt(ix,1,jx,4) = Zero
            dBt(jx,4,ix,1) = dBt(ix,1,jx,4)
            dBt(ix,1,jx,3) = - (dBt(ix,1,jx,1) + dBt(ix,1,jx,2) )
            dBt(jx,3,ix,1) = dBt(ix,1,jx,3)
            dBt(ix,4,jx,4) =(  dBRkl(ix,2,jy,1)*BRjk(jz,1)&
               &                       - dBRkl(ix,2,jz,1)*BRjk(jy,1)&
               &                       - Bt(jx,4)*(BRkl(ix,2)*SinFi3**2&
               &                       + Rkl1*Two*SinFi3*CosFi3*Bf3(ix,3)) )&
               &                       / (Rkl1*SinFi3**2)
            dBt(ix,4,jx,3) =-( (-BRkl(ix,2)*CosFi3&
               &                         + Rkl1*SinFi3*Bf3(ix,3) ) * Bt(jx,4)&
               &                         + (Rjk1-Rkl1*CosFi3) * dBt(ix,4,jx,4) )&
               &                       / Rjk1
            dBt(jx,3,ix,4) = dBt(ix,4,jx,3)
            dBt(ix,4,jx,2) = - ( dBt(ix,4,jx,4) + dBt(ix,4,jx,3) )
            dBt(jx,2,ix,4) = dBt(ix,4,jx,2)
            If (ix.ne.jx) Then
               dBt(jx,1,ix,1) = dBt(ix,1,jx,1)
               dBt(ix,4,jx,1) = Zero
               dBt(jx,4,ix,4) = dBt(ix,4,jx,4)
               dBt(jx,1,ix,4) = dBt(ix,4,jx,1)
               dBt(jx,1,ix,2) =-( (-BRij(jx,1)*CosFi2&
                  &                            + Rij1*SinFi2*Bf2(jx,1) ) * Bt(ix,1)&
                  &                            + (Rjk1-Rij1*CosFi2) * dBt(jx,1,ix,1))&
                  &                          / Rjk1
               dBt(ix,2,jx,1) = dBt(jx,1,ix,2)
               dBt(ix,3,jx,1) = - ( dBt(ix,1,jx,1) + dBt(ix,2,jx,1)&
                  &                          + dBt(ix,4,jx,1) )
               dBt(jx,1,ix,3) = dBt(ix,3,jx,1)
               dBt(jx,4,ix,3) =-( (-BRkl(jx,2)*CosFi3&
                  &                            + Rkl1*SinFi3*Bf3(jx,3) ) * Bt(ix,4)&
                  &                            + (Rjk1-Rkl1*CosFi3) * dBt(jx,4,ix,4))&
                  &                          / Rjk1
               dBt(ix,3,jx,4) = dBt(jx,4,ix,3)
               dBt(ix,2,jx,4) = - ( dBt(ix,4,jx,4) + dBt(ix,3,jx,4) )
               dBt(jx,4,ix,2) = dBt(ix,2,jx,4)
            End If
            dBt(ix,2,jx,3) = - ( ( BRjk(ix,1)&
               &                           + Rkl1*SinFi3*Bf3(ix,1) ) * Bt(jx,4)&
               &                         + ( Rjk1 - Rkl1*CosFi3 ) * dBt(ix,2,jx,4)&
               &                         + ( BRij(ix,2)*CosFi2&
               &                           - Rij1*SinFi2*Bf2(ix,2) ) * Bt(jx,1)&
               &                         +  Rij1*CosFi2 * dBt(ix,2,jx,1)&
               &                         + Bt(jx,3) * BRjk(ix,1) ) / Rjk1
            dBt(jx,3,ix,2) = dBt(ix,2,jx,3)
            dBt(ix,2,jx,2) = - ( dBt(ix,2,jx,1) + dBt(ix,2,jx,4)&
               &                         + dBt(ix,2,jx,3) )
            dBt(ix,3,jx,3) = - ( dBt(ix,2,jx,3) + dBt(ix,1,jx,3)&
               &                         + dBt(ix,4,jx,3) )
            If (ix.ne.jx) Then
               dBt(ix,3,jx,2) = - ( dBt(ix,2,jx,2) + dBt(ix,1,jx,2)&
                  &                            + dBt(ix,4,jx,2) )
               dBt(jx,2,ix,3) = dBt(ix,3,jx,2)
               dBt(jx,2,ix,2) = dBt(ix,2,jx,2)
               dBt(jx,3,ix,3) = dBt(ix,3,jx,3)
            End If
            !
         End Do
      End Do
      !
   End If
   !     Call qExit('Trsn')
   Return
End      subroutine

Subroutine Strtch(xyz,nCent,Avst,B,lWrite,Label,dB,ldB)
   Implicit Integer (i-n)
   Implicit Real*8 (a-h,o-z)
   !      include "common/real.inc"
   !comdeck real.inc $Revision: 2002.3 $
   Real*8 Zero, One, Two, Three, Four, Five, Six, Seven,&
      &       Eight, RNine, Ten, Half, Pi, SqrtP2, TwoP34,&
      &       TwoP54, One2C2
   Parameter(Zero =0.0D0, One  =1.0D0, Two=2.0D0, Three=3.0D0,&
      &          Four =4.0D0, Five =5.0D0, Six=6.0D0, Seven=7.0D0,&
      &          Eight=8.0D0, rNine=9.0D0, Ten=1.0D1, Half=0.5D0,&
      &          Pi    =3.141592653589793D0,&
      &          SqrtP2=0.8862269254527579D0,&
      &          TwoP34=0.2519794355383808D0,&
      &          TwoP54=5.914967172795612D0,&
      &          One2C2=0.2662567690426443D-04)

   Real*8 B(3,nCent), xyz(3,nCent), dB(3,nCent,3,nCent), R(3)
   Logical lWrite, ldB
   Character*8 Label
   !      include "common/angstr.inc"
   !comdeck angstr.inc $Revision: 2002.3 $
   !
   !     Conversion factor angstrom to bohr from the IUPAC
   !     publication
   !     .529177249(24) angstrom / bohr
   !     "Quantities, Units and Symbols in Physical Chemistry"
   !     I. Mills, T. Cvitas, K. Homann, N. Kallay and
   !     K. Kuchitsu, Blackwell Scientific Publications,
   !     Oxford, 1988.
   !
   Data Angstr/0.529177249D+00/
   !
   R(1)=xyz(1,2)-xyz(1,1)
   R(2)=xyz(2,2)-xyz(2,1)
   R(3)=xyz(3,2)-xyz(3,1)
   R2=R(1)**2+R(2)**2+R(3)**2
   RR=Sqrt(R2)
   Avst=RR
   !
   aRR=RR*Angstr
   If (lWrite) Write (*,'(1X,A,A,2(F10.6,A))') Label,&
      &      ' : Bond Length=',aRR,' / Angstrom',RR,' / bohr'
   !
   !---- Compute the WDC B-matrix.
   !
   B(1,1)=-R(1)/RR
   B(2,1)=-R(2)/RR
   B(3,1)=-R(3)/RR
   !.... Utilize translational invariance.
   B(1,2)=-B(1,1)
   B(2,2)=-B(2,1)
   B(3,2)=-B(3,1)
   !
   !---- Compute the cartesian derivative of the B-matrix.
   !
   If (ldB) Then
      !
      Do i=1,3
         Do j=1,i
            If (i.eq.j) Then
               dB(i,1,j,1)= (One-B(j,1)*B(i,1))/RR
            Else
               dB(i,1,j,1)= (-B(j,1)*B(i,1))/RR
            End If
            dB(j,1,i,1)= dB(i,1,j,1)
            !
            dB(i,2,j,1)=-dB(i,1,j,1)
            dB(j,1,i,2)= dB(i,2,j,1)
            !
            dB(i,1,j,2)=-dB(i,1,j,1)
            dB(j,2,i,1)= dB(i,1,j,2)
            !
            dB(i,2,j,2)=-dB(i,2,j,1)
            dB(j,2,i,2)=dB(i,2,j,2)
         End Do
      End Do
      !
   End If
   !     Call qExit('Strtch')
   !     Call GetMem('Exit Strtch','Chec','Real',ipMass,2*msAtom)
   Return
End subroutine

!fordeck bend $Revision: 5.3 $
Subroutine Bend(xyz,nCent,Fir,Bf,lWrite,lWarn,Label,dBf,ldB)
   Implicit Integer (i-n)
   Implicit Real*8  (a-h,o-z)
   !      include "common/real.inc"
   !comdeck real.inc $Revision: 2002.3 $
   Real*8 Zero, One, Two, Three, Four, Five, Six, Seven,&
      &       Eight, RNine, Ten, Half, Pi, SqrtP2, TwoP34,&
      &       TwoP54, One2C2
   Parameter(Zero =0.0D0, One  =1.0D0, Two=2.0D0, Three=3.0D0,&
      &          Four =4.0D0, Five =5.0D0, Six=6.0D0, Seven=7.0D0,&
      &          Eight=8.0D0, rNine=9.0D0, Ten=1.0D1, Half=0.5D0,&
      &          Pi    =3.141592653589793D0,&
      &          SqrtP2=0.8862269254527579D0,&
      &          TwoP34=0.2519794355383808D0,&
      &          TwoP54=5.914967172795612D0,&
      &          One2C2=0.2662567690426443D-04)

   Real*8   Bf(3,nCent), xyz(3,nCent), dBf(3,nCent,3,nCent),&
      &        BRij(3,2), dBRij(3,2,3,2),&
      &        BRjk(3,2), dBRjk(3,2,3,2)
   Logical lWrite, ldB, lWarn
   Character*8 Label
   !
   !     Call QEnter('Bend')
   !
   mCent=2
   Call Strtch(xyz(1,1),mCent,Rij1,BRij,.False.,Label,dBRij,ldB)
   Call Strtch(xyz(1,2),mCent,Rjk1,BRjk,.False.,Label,dBRjk,ldB)
   Co=Zero
   Crap=Zero
   Do i = 1, 3
      Co=Co+BRij(i,1)*BRjk(i,2)
      Crap=Crap+(BRjk(i,2)+BRij(i,1))**2
   End Do
   !
   !.... Special care for cases close to linearity
   !
   If (Sqrt(Crap).lt.1.0D-6) Then
      Fir=Pi-ArSin(Sqrt(Crap))
      Si=Sqrt(Crap)
   Else
      Fir=ArCos(Co)
      Si=Sqrt(One-Co**2)
   End If
   !
   If (Abs(Fir-Pi).lt.1.0d-13) Then
      Fir=Pi
      Return
   End If
   dFir=180.0D0*Fir/Pi
   If ((Abs(dFir).gt.177.5 .or. Abs(dFir).lt.2.5).and.lWarn)&
      &   Write (*,*) ' Valence angle close to end in '//&
      &               'range of definition'
   If (lWrite) Write (*,'(1X,A,A,F10.4,A,F10.6,A)') Label,&
      &            ' : Angle=', dFir,'/degree, ',Fir,'/rad'
   !
   !---- Compute the WDC B-matrix
   !
   !     Bf=-11.1111
   Do i = 1, 3
      Bf(i,1)= (Co*BRij(i,1)-BRjk(i,2))/(Si*Rij1)
      Bf(i,3)= (Co*BRjk(i,2)-BRij(i,1))/(Si*Rjk1)
      !....... Utilize translational invariance.
      Bf(i,2)=-(Bf(i,1)+Bf(i,3))
   End Do
   !     Call RecPrt('Bf',' ',Bf,9,1)
   !
   !---- Compute the cartesian derivative of the B-Matrix.
   !
   If (ldB) Then
      !
      !        dBf=-11.11111
      Do i = 1, 3
         Do j = 1, i
            dBf(i,1,j,1)=( -Si*Bf(i,1)*BRij(j,1)&
               &                        +Co*dBRij(i,1,j,1)&
               &                        -Bf(j,1)*(Co*Bf(i,1)*Rij1&
               &                        +Si*BRij(i,1)) ) / (Si*Rij1)
            dBf(i,1,j,3)=(-Si*Bf(i,1)*BRjk(j,2)&
               &                       +dBRij(i,1,j,2)&
               &                       -Bf(j,3)*Co*Bf(i,1)*Rjk1)&
               &                       / (Si*Rjk1)
            !              Write (*,*) '13',dBf(i,1,j,3), i, j
            dBf(i,3,j,1)=(-Si*Bf(i,3)*BRij(j,1)&
               &                       +dBRjk(i,2,j,1)&
               &                       -Bf(j,1)*Co*Bf(i,3)*Rij1)&
               &                       / (Si*Rij1)
            dBf(i,3,j,3)=( -Si*Bf(i,3)*BRjk(j,2)&
               &                        +Co*dBRjk(i,2,j,2)&
               &                        -Bf(j,3)*(Co*Bf(i,3)*Rjk1&
               &                        +Si*BRjk(i,2)) ) / (Si*Rjk1)
            !
            dBf(j,1,i,1)=dBf(i,1,j,1)
            dBf(j,3,i,1)=dBf(i,1,j,3)
            dBf(j,1,i,3)=dBf(i,3,j,1)
            dBf(j,3,i,3)=dBf(i,3,j,3)
            !
            dBf(i,1,j,2)=-(dBf(i,1,j,1)+dBf(i,1,j,3))
            dBf(j,2,i,1)=dBf(i,1,j,2)
            dBf(j,1,i,2)=-(dBf(j,1,i,1)+dBf(j,1,i,3))
            dBf(i,2,j,1)=dBf(j,1,i,2)
            dBf(i,3,j,2)=-(dBf(i,3,j,1)+dBf(i,3,j,3))
            dBf(j,2,i,3)=dBf(i,3,j,2)
            dBf(j,3,i,2)=-(dBf(j,3,i,1)+dBf(j,3,i,3))
            dBf(i,2,j,3)=dBf(j,3,i,2)
            !
            dBf(i,2,j,2)=-(dBf(i,2,j,1)+dBf(i,2,j,3))
            dBf(j,2,i,2)=dBf(i,2,j,2)
            !
         End Do
      End Do
      !        Call RecPrt('dBf','(9F9.1)',dBf,9,9)
      !
   End If
   !
   !     Call QExit('Bend')
   Return
End subroutine

!fordeck arsin.f $Revision: 5.3 $
Function arSin(Arg)
   Implicit Real*8 (a-h,o-z)
   Real*8 ArSin
   !      include "common/real.inc"
   !comdeck real.inc $Revision: 2002.3 $
   Real*8 Zero, One, Two, Three, Four, Five, Six, Seven,&
      &       Eight, RNine, Ten, Half, Pi, SqrtP2, TwoP34,&
      &       TwoP54, One2C2
   Parameter(Zero =0.0D0, One  =1.0D0, Two=2.0D0, Three=3.0D0,&
      &          Four =4.0D0, Five =5.0D0, Six=6.0D0, Seven=7.0D0,&
      &          Eight=8.0D0, rNine=9.0D0, Ten=1.0D1, Half=0.5D0,&
      &          Pi    =3.141592653589793D0,&
      &          SqrtP2=0.8862269254527579D0,&
      &          TwoP34=0.2519794355383808D0,&
      &          TwoP54=5.914967172795612D0,&
      &          One2C2=0.2662567690426443D-04)

   A=Arg
   IF(ABS(A).GT.One) Then
      PRINT 3,A
      3        FORMAT(1X,'Warning argument of aSin= ',1F21.18)
      A=Sign(One,A)
   End If
   !
   ArSin=ASin(A)
   Return
End function

!fordeck arcos.f $Revision: 5.3 $
Function arCos(Arg)
   Implicit Real*8 (a-h,o-z)
   Real*8 ArCos
   !      include "common/real.inc"
   !comdeck real.inc $Revision: 2002.3 $
   Real*8 Zero, One, Two, Three, Four, Five, Six, Seven,&
      &       Eight, RNine, Ten, Half, Pi, SqrtP2, TwoP34,&
      &       TwoP54, One2C2
   Parameter(Zero =0.0D0, One  =1.0D0, Two=2.0D0, Three=3.0D0,&
      &          Four =4.0D0, Five =5.0D0, Six=6.0D0, Seven=7.0D0,&
      &          Eight=8.0D0, rNine=9.0D0, Ten=1.0D1, Half=0.5D0,&
      &          Pi    =3.141592653589793D0,&
      &          SqrtP2=0.8862269254527579D0,&
      &          TwoP34=0.2519794355383808D0,&
      &          TwoP54=5.914967172795612D0,&
      &          One2C2=0.2662567690426443D-04)

   A=Arg
   IF (ABS(A).GT.One) Then
      !         PRINT 3,A
      !3        FORMAT(1X,'Warning argument of aCos= ',1F21.18)
      A=Sign(One,A)
   End If
   !
   ArCos=ACos(A)
   Return
End function
