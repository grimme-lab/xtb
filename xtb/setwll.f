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

cccccccccccccccccccccccccccccccccccccc
c set atomic spin constants in Eh 
c lower triangle of s,p,d * s,p,d
cccccccccccccccccccccccccccccccccccccc
      subroutine setwll_pbe(f)    ! f=1.9
      use mctc_econv, only : autoev
      use aoparam
      implicit none

      real*8 f   
      integer lin,ss,sp,sd,pp,pd,dd,ff

      wll(1:94,1:10)=0
c index
      ss=1
      sp=2
      sd=3
      pp=lin(2,2)
      pd=lin(2,3)
      dd=6
      ff=10
        
ccccccccccccccccccccccc
c data      
ccccccccccccccccccccccc

      wll(1,ss)=-0.12427000
                 
      wll(2,ss)=-0.08653250
                 
      wll(3,ss)=-0.02672500
                 
      wll(4,ss)=-0.02292750
                 
      wll(5,ss)=-0.02767000
      wll(5,sp)=-0.0223200000
      wll(5,pp)=-0.02008500
                 
      wll(6,ss)=-0.03189750
      wll(6,sp)=-0.0258375000
      wll(6,pp)=-0.02334500
                 
      wll(7,ss)=-0.03551000
      wll(7,sp)=-0.0286600000
      wll(7,pp)=-0.02598000
                 
      wll(8,ss)=-0.03646750
      wll(8,sp)=-0.0302525000
      wll(8,pp)=-0.02802750
                 
      wll(9,ss)=-0.03716000
      wll(9,sp)=-0.0313825000
      wll(9,pp)=-0.02996750
       
      wll(10,ss)=-0.03835000
      wll(10,sp)=-0.0326350000
      wll(10,pp)=-0.03166500
       
      wll(11,ss)=-0.02015000
       
      wll(12,ss)=-0.01656500
       
      wll(13,ss)=-0.01847500
      wll(13,sp)=-0.0142500000
      wll(13,pp)=-0.01485000
       
      wll(14,ss)=-0.02031500
      wll(14,sp)=-0.0160375000
      wll(14,pp)=-0.01592000
       
      wll(15,ss)=-0.02173250
      wll(15,sp)=-0.0175750000
      wll(15,pp)=-0.01704500
       
      wll(16,ss)=-0.02189750
      wll(16,sp)=-0.0175987500
      wll(16,pp)=-0.01615750
       
      wll(17,ss)=-0.02196250
      wll(17,sp)=-0.0178575000
      wll(17,pp)=-0.01619750
       
      wll(18,ss)=-0.02221000
      wll(18,sp)=-0.0183425000
      wll(18,pp)=-0.01659000
       
      wll(19,ss)=-0.01433750
       
      wll(20,ss)=-0.01181250
       
      wll(31,ss)=-0.01731250
      wll(31,sp)=-0.0129837500
      wll(31,pp)=-0.01434500
       
      wll(32,ss)=-0.01780000
      wll(32,sp)=-0.0138462500
      wll(32,pp)=-0.01463250
       
      wll(33,ss)=-0.01821500
      wll(33,sp)=-0.0145262500
      wll(33,pp)=-0.01500000
                  
      wll(34,ss)=-0.01831000
      wll(34,sp)=-0.0145500000
      wll(34,pp)=-0.01416000
                  
      wll(35,ss)=-0.01816750
      wll(35,sp)=-0.0144787500
      wll(35,pp)=-0.01381000
                  
      wll(36,ss)=-0.01812750
      wll(36,sp)=-0.0146025000
      wll(36,pp)=-0.01377000
       
      wll(37,ss)=-0.01383000
       
      wll(38,ss)=-0.01065000
       
      wll(49,ss)=-0.01414250
      wll(49,sp)=-0.0106275000
      wll(49,pp)=-0.01257750
                  
      wll(50,ss)=-0.01429750
      wll(50,sp)=-0.0110262500
      wll(50,pp)=-0.01247500
                  
      wll(51,ss)=-0.01447000
      wll(51,sp)=-0.0113687500
      wll(51,pp)=-0.01247250
                  
      wll(52,ss)=-0.01471500
      wll(52,sp)=-0.0115250000
      wll(52,pp)=-0.01189250
                  
      wll(53,ss)=-0.01466250
      wll(53,sp)=-0.0114325000
      wll(53,pp)=-0.01156250
                  
      wll(54,ss)=-0.01463500
      wll(54,sp)=-0.0114462500
      wll(54,pp)=-0.01141500
       
      wll(55,ss)=-0.01214500
       
      wll(56,ss)=-0.00928750
       
      wll(81,ss)=-0.01328500
      wll(81,sp)=-0.0091625000
      wll(81,pp)=-0.01248000
                  
      wll(82,ss)=-0.01347250
      wll(82,sp)=-0.0094275000
      wll(82,pp)=-0.01219250
                  
      wll(83,ss)=-0.01359750
      wll(83,sp)=-0.0096312500
      wll(83,pp)=-0.01200500
                  
      wll(84,ss)=-0.01383000
      wll(84,sp)=-0.0098412500
      wll(84,pp)=-0.01104250
                  
      wll(85,ss)=-0.01384250
      wll(85,sp)=-0.0097612500
      wll(85,pp)=-0.01079500
                  
      wll(86,ss)=-0.01387250
      wll(86,sp)=-0.0097450000
      wll(86,pp)=-0.01062000

ccccccccccccccccccccccccccccccccccccccccc
       
      wll(21,ss)=-0.01275750
      wll(21,sp)=0
      wll(21,sd)=-0.0048125000
      wll(21,pp)=0
      wll(21,pd)=0
      wll(21,dd)=-0.01238500
       
      wll(22,ss)=-0.01342500
      wll(22,sp)=0
      wll(22,sd)=-0.0044387500
      wll(22,pp)=0
      wll(22,pd)=0
      wll(22,dd)=-0.01354250
       
      wll(23,ss)=-.01396250
      wll(23,sp)=0
      wll(23,sd)=-.0041950000
      wll(23,pp)=0
      wll(23,pd)=0
      wll(23,dd)=-.01436000
       
      wll(24,ss)=-0.01442250
      wll(24,sp)=0
      wll(24,sd)=-0.0040387500
      wll(24,pp)=0
      wll(24,pd)=0
      wll(24,dd)=-0.01499750
       
      wll(25,ss)=-0.01485750
      wll(25,sp)=0
      wll(25,sd)=-0.0039337500
      wll(25,pp)=0
      wll(25,pd)=0
      wll(25,dd)=-0.01552250
       
      wll(26,ss)=-0.01540000
      wll(26,sp)=0
      wll(26,sd)=-0.0036000000
      wll(26,pp)=0
      wll(26,pd)=0
      wll(26,dd)=-0.01679000
       
      wll(27,ss)=-0.01580000
      wll(27,sp)=0
      wll(27,sd)=-0.0032912500
      wll(27,pp)=0
      wll(27,pd)=0
      wll(27,dd)=-0.01762500
       
      wll(28,ss)=-0.01616000
      wll(28,sp)=0
      wll(28,sd)=-0.0030625000
      wll(28,pp)=0
      wll(28,pd)=0
      wll(28,dd)=-0.01825500
       
      wll(29,ss)=-0.01656000
      wll(29,sp)=0
      wll(29,sd)=-0.0028575000
      wll(29,pp)=0
      wll(29,pd)=0
      wll(29,dd)=-0.01880000
       
      wll(30,ss)=-0.01684000
      wll(30,sp)=0
      wll(30,sd)=-0.0027650000
      wll(30,pp)=0
      wll(30,pd)=0
      wll(30,dd)=-0.01926000
       
      wll(39,ss)=-0.01140750
      wll(39,sp)=0
      wll(39,sd)=-0.0067850000
      wll(39,pp)=0
      wll(39,pd)=0
      wll(39,dd)=-0.00966500
       
      wll(40,ss)=-0.01194000
      wll(40,sp)=0
      wll(40,sd)=-0.0061900000
      wll(40,pp)=0
      wll(40,pd)=0
      wll(40,dd)=-0.01054750
       
      wll(41,ss)=-0.01233750
      wll(41,sp)=0
      wll(41,sd)=-0.0057312500
      wll(41,pp)=0
      wll(41,pd)=0
      wll(41,dd)=-0.01115750
       
      wll(42,ss)=-0.01265250
      wll(42,sp)=0
      wll(42,sd)=-0.0053500000
      wll(42,pp)=0
      wll(42,pd)=0
      wll(42,dd)=-0.01167250
       
      wll(43,ss)=-0.01291750
      wll(43,sp)=0
      wll(43,sd)=-0.0050150000
      wll(43,pp)=0
      wll(43,pd)=0
      wll(43,dd)=-0.01217250
       
      wll(44,ss)=-0.01322250
      wll(44,sp)=0
      wll(44,sd)=-0.0047550000
      wll(44,pp)=0
      wll(44,pd)=0
      wll(44,dd)=-0.01274000
       
      wll(45,ss)=-0.01332750
      wll(45,sp)=0
      wll(45,sd)=-0.0043787500
      wll(45,pp)=0
      wll(45,pd)=0
      wll(45,dd)=-0.01286750
       
      wll(46,ss)=-0.01346750
      wll(46,sp)=0
      wll(46,sd)=-0.0041000000
      wll(46,pp)=0
      wll(46,pd)=0
      wll(46,dd)=-0.01302500
       
      wll(47,ss)=-0.01370500
      wll(47,sp)=0
      wll(47,sd)=-0.0038637500
      wll(47,pp)=0
      wll(47,pd)=0
      wll(47,dd)=-0.01319750
       
      wll(48,ss)=-0.01388500
      wll(48,sp)=0
      wll(48,sd)=-0.0036875000
      wll(48,pp)=0
      wll(48,pd)=0
      wll(48,dd)=-0.01340000
       
      wll(57,ss)=-0.00992750
      wll(57,sp)=0
      wll(57,sd)=-0.0059287500
      wll(57,pp)=0
      wll(57,pd)=0
      wll(57,dd)=-0.00896250
       
      wll(72,ss)=-0.01228500
      wll(72,sp)=0
      wll(72,sd)=-0.0077550000
      wll(72,pp)=0
      wll(72,pd)=0
      wll(72,dd)=-0.01025500
       
      wll(73,ss)=-0.01252500
      wll(73,sp)=0
      wll(73,sd)=-0.0073250000
      wll(73,pp)=0
      wll(73,pd)=0
      wll(73,dd)=-0.01075000
       
      wll(74,ss)=-0.01271750
      wll(74,sp)=0
      wll(74,sd)=-0.0069475000
      wll(74,pp)=0
      wll(74,pd)=0
      wll(74,dd)=-0.01110750
       
      wll(75,ss)=-0.01287250
      wll(75,sp)=0
      wll(75,sd)=-0.0065937500
      wll(75,pp)=0
      wll(75,pd)=0
      wll(75,dd)=-0.01142250
       
      wll(76,ss)=-0.01299000
      wll(76,sp)=0
      wll(76,sd)=-0.0062725000
      wll(76,pp)=0
      wll(76,pd)=0
      wll(76,dd)=-0.01184750
       
      wll(77,ss)=-0.01292000
      wll(77,sp)=0
      wll(77,sd)=-0.0058575000
      wll(77,pp)=0
      wll(77,pd)=0
      wll(77,dd)=-0.01185500
       
      wll(78,ss)=-0.01285750
      wll(78,sp)=0
      wll(78,sd)=-0.0055112500
      wll(78,pp)=0
      wll(78,pd)=0
      wll(78,dd)=-0.01180000
       
      wll(79,ss)=-0.01291000
      wll(79,sp)=0
      wll(79,sd)=-0.0052475000
      wll(79,pp)=0
      wll(79,pd)=0
      wll(79,dd)=-0.01179250
       
      wll(80,ss)=-0.01292500
      wll(80,sp)=0
      wll(80,sd)=-0.0050612500
      wll(80,pp)=0
      wll(80,pd)=0
      wll(80,dd)=-0.01184000

c f-elements
      wll(58:71,ss)=-0.01
      wll(90:94,ss)=-0.01
      wll(58:71,ff)=-0.01
      wll(90:94,ff)=-0.01

c convert to eV      
      wll = wll * autoev * f

      end

