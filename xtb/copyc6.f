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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c copy from machine generated data statements inside pars.f
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine copyc6
      use iso_fortran_env, wp => real64
      use d3param
      implicit none

      integer nlines
      integer iat,jat,iadr,jadr,nn,kk
      include 'pars.inc'

!     define r2r4 and rcov at this point

      call setr2r4(r2r4)
      call setrcov(rcov)

      c6ab=-1
      mxc=0
!     process all entries in pars.f 
      kk=1
      do nn=1,nlines
       iat=int(pars(kk+1))
       jat=int(pars(kk+2))
       call limit(iat,jat,iadr,jadr)
       mxc(iat)=max(mxc(iat),iadr)
       mxc(jat)=max(mxc(jat),jadr)

       c6ab(iat,jat,iadr,jadr,1)=pars(kk)  
       c6ab(iat,jat,iadr,jadr,2)=pars(kk+3)
       c6ab(iat,jat,iadr,jadr,3)=pars(kk+4)

       c6ab(jat,iat,jadr,iadr,1)=pars(kk) 
       c6ab(jat,iat,jadr,iadr,2)=pars(kk+4)
       c6ab(jat,iat,jadr,iadr,3)=pars(kk+3)
       kk=(nn*5)+1
      enddo

      end subroutine copyc6

      subroutine limit(iat,jat,iadr,jadr)
      implicit none
      integer iat,jat,iadr,jadr,i
      iadr=1
      jadr=1
      i=100
 10   if(iat.gt.100) then
         iat=iat-100
         iadr=iadr+1
         goto 10
      endif

      i=100
 20   if(jat.gt.100) then
         jat=jat-100
         jadr=jadr+1
         goto 20
      endif

      end subroutine limit

      subroutine setrcov(rcov)
      use iso_fortran_env, wp => real64
      real(wp) rcov(94)

      rcov( 1 )= 0.80628308
      rcov( 2 )= 1.15903197
      rcov( 3 )= 3.02356173
      rcov( 4 )= 2.36845659
      rcov( 5 )= 1.94011865
      rcov( 6 )= 1.88972601
      rcov( 7 )= 1.78894056
      rcov( 8 )= 1.58736983
      rcov( 9 )= 1.61256616
      rcov( 10 )= 1.68815527
      rcov( 11 )= 3.52748848
      rcov( 12 )= 3.14954334
      rcov( 13 )= 2.84718717
      rcov( 14 )= 2.62041997
      rcov( 15 )= 2.77159820
      rcov( 16 )= 2.57002732
      rcov( 17 )= 2.49443835
      rcov( 18 )= 2.41884923
      rcov( 19 )= 4.43455700
      rcov( 20 )= 3.88023730
      rcov( 21 )= 3.35111422
      rcov( 22 )= 3.07395437
      rcov( 23 )= 3.04875805
      rcov( 24 )= 2.77159820
      rcov( 25 )= 2.69600923
      rcov( 26 )= 2.62041997
      rcov( 27 )= 2.51963467
      rcov( 28 )= 2.49443835
      rcov( 29 )= 2.54483100
      rcov( 30 )= 2.74640188
      rcov( 31 )= 2.82199085
      rcov( 32 )= 2.74640188
      rcov( 33 )= 2.89757982
      rcov( 34 )= 2.77159820
      rcov( 35 )= 2.87238349
      rcov( 36 )= 2.94797246
      rcov( 37 )= 4.76210950
      rcov( 38 )= 4.20778980
      rcov( 39 )= 3.70386304
      rcov( 40 )= 3.50229216
      rcov( 41 )= 3.32591790
      rcov( 42 )= 3.12434702
      rcov( 43 )= 2.89757982
      rcov( 44 )= 2.84718717
      rcov( 45 )= 2.84718717
      rcov( 46 )= 2.72120556
      rcov( 47 )= 2.89757982
      rcov( 48 )= 3.09915070
      rcov( 49 )= 3.22513231
      rcov( 50 )= 3.17473967
      rcov( 51 )= 3.17473967
      rcov( 52 )= 3.09915070
      rcov( 53 )= 3.32591790
      rcov( 54 )= 3.30072128
      rcov( 55 )= 5.26603625
      rcov( 56 )= 4.43455700
      rcov( 57 )= 4.08180818
      rcov( 58 )= 3.70386304
      rcov( 59 )= 3.98102289
      rcov( 60 )= 3.95582657
      rcov( 61 )= 3.93062995
      rcov( 62 )= 3.90543362
      rcov( 63 )= 3.80464833
      rcov( 64 )= 3.82984466
      rcov( 65 )= 3.80464833
      rcov( 66 )= 3.77945201
      rcov( 67 )= 3.75425569
      rcov( 68 )= 3.75425569
      rcov( 69 )= 3.72905937
      rcov( 70 )= 3.85504098
      rcov( 71 )= 3.67866672
      rcov( 72 )= 3.45189952
      rcov( 73 )= 3.30072128
      rcov( 74 )= 3.09915070
      rcov( 75 )= 2.97316878
      rcov( 76 )= 2.92277614
      rcov( 77 )= 2.79679452
      rcov( 78 )= 2.82199085
      rcov( 79 )= 2.84718717
      rcov( 80 )= 3.32591790
      rcov( 81 )= 3.27552496
      rcov( 82 )= 3.27552496
      rcov( 83 )= 3.42670319
      rcov( 84 )= 3.30072128
      rcov( 85 )= 3.47709584
      rcov( 86 )= 3.57788113
      rcov( 87 )= 5.06446567
      rcov( 88 )= 4.56053862
      rcov( 89 )= 4.20778980
      rcov( 90 )= 3.98102289
      rcov( 91 )= 3.82984466
      rcov( 92 )= 3.85504098
      rcov( 93 )= 3.88023730
      rcov( 94 )= 3.90543362

      return 

      end subroutine setrcov
 
      subroutine setr2r4(r2r4)
      use iso_fortran_env, wp => real64
      real(wp) r2r4(94)

      r2r4( 1 )= 2.00734898
      r2r4( 2 )= 1.56637132
      r2r4( 3 )= 5.01986934
      r2r4( 4 )= 3.85379032
      r2r4( 5 )= 3.64446594
      r2r4( 6 )= 3.10492822
      r2r4( 7 )= 2.71175247
      r2r4( 8 )= 2.59361680
      r2r4( 9 )= 2.38825250
      r2r4( 10 )= 2.21522516
      r2r4( 11 )= 6.58585536
      r2r4( 12 )= 5.46295967
      r2r4( 13 )= 5.65216669
      r2r4( 14 )= 4.88284902
      r2r4( 15 )= 4.29727576
      r2r4( 16 )= 4.04108902
      r2r4( 17 )= 3.72932356
      r2r4( 18 )= 3.44677275
      r2r4( 19 )= 7.97762753
      r2r4( 20 )= 7.07623947
      r2r4( 21 )= 6.60844053
      r2r4( 22 )= 6.28791364
      r2r4( 23 )= 6.07728703
      r2r4( 24 )= 5.54643096
      r2r4( 25 )= 5.80491167
      r2r4( 26 )= 5.58415602
      r2r4( 27 )= 5.41374528
      r2r4( 28 )= 5.28497229
      r2r4( 29 )= 5.22592821
      r2r4( 30 )= 5.09817141
      r2r4( 31 )= 6.12149689
      r2r4( 32 )= 5.54083734
      r2r4( 33 )= 5.06696878
      r2r4( 34 )= 4.87005108
      r2r4( 35 )= 4.59089647
      r2r4( 36 )= 4.31176304
      r2r4( 37 )= 9.55461698
      r2r4( 38 )= 8.67396077
      r2r4( 39 )= 7.97210197
      r2r4( 40 )= 7.43439917
      r2r4( 41 )= 6.58711862
      r2r4( 42 )= 6.19536215
      r2r4( 43 )= 6.01517290
      r2r4( 44 )= 5.81623410
      r2r4( 45 )= 5.65710424
      r2r4( 46 )= 5.52640661
      r2r4( 47 )= 5.44263305
      r2r4( 48 )= 5.58285373
      r2r4( 49 )= 7.02081898
      r2r4( 50 )= 6.46815523
      r2r4( 51 )= 5.98089120
      r2r4( 52 )= 5.81686657
      r2r4( 53 )= 5.53321815
      r2r4( 54 )= 5.25477007
      r2r4( 55 )= 11.02204549
      r2r4( 56 )= 10.15679528
      r2r4( 57 )= 9.35167836
      r2r4( 58 )= 9.06926079
      r2r4( 59 )= 8.97241155
      r2r4( 60 )= 8.90092807
      r2r4( 61 )= 8.85984840
      r2r4( 62 )= 8.81736827
      r2r4( 63 )= 8.79317710
      r2r4( 64 )= 7.89969626
      r2r4( 65 )= 8.80588454
      r2r4( 66 )= 8.42439218
      r2r4( 67 )= 8.54289262
      r2r4( 68 )= 8.47583370
      r2r4( 69 )= 8.45090888
      r2r4( 70 )= 8.47339339
      r2r4( 71 )= 7.83525634
      r2r4( 72 )= 8.20702843
      r2r4( 73 )= 7.70559063
      r2r4( 74 )= 7.32755997
      r2r4( 75 )= 7.03887381
      r2r4( 76 )= 6.68978720
      r2r4( 77 )= 6.05450052
      r2r4( 78 )= 5.88752022
      r2r4( 79 )= 5.70661499
      r2r4( 80 )= 5.78450695
      r2r4( 81 )= 7.79780729
      r2r4( 82 )= 7.26443867
      r2r4( 83 )= 6.78151984
      r2r4( 84 )= 6.67883169
      r2r4( 85 )= 6.39024318
      r2r4( 86 )= 6.09527958
      r2r4( 87 )= 11.79156076
      r2r4( 88 )= 11.10997644
      r2r4( 89 )= 9.51377795
      r2r4( 90 )= 8.67197068
      r2r4( 91 )= 8.77140725
      r2r4( 92 )= 8.65402716
      r2r4( 93 )= 8.53923501
      r2r4( 94 )= 8.85024712

      return 

      end subroutine setr2r4
         
