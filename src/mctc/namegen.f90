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

!> Generate a new random name
module xtb_mctc_namegen
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: generateName


   character(len=1), parameter :: delimiter = '-'

   character(len=*), parameter :: tWord(61) = [ &
      &"talkative  ", "tall       ", "tame       ", "tan        ", "tangible   ", &
      &"tart       ", "tasty      ", "tattered   ", "taut       ", "tedious    ", &
      &"teeming    ", "tempting   ", "tender     ", "tense      ", "tepid      ", &
      &"terrible   ", "terrific   ", "testy      ", "thankful   ", "that       ", &
      &"these      ", "thick      ", "thin       ", "third      ", "thirsty    ", &
      &"this       ", "thorny     ", "thorough   ", "those      ", "thoughtful ", &
      &"threadbare ", "thrifty    ", "thunderous ", "tidy       ", "tight      ", &
      &"timely     ", "tinted     ", "tiny       ", "tired      ", "torn       ", &
      &"total      ", "tough      ", "tragic     ", "trained    ", "traumatic  ", &
      &"treasured  ", "tremendous ", "triangular ", "tricky     ", "trifling   ", &
      &"trim       ", "trivial    ", "troubled   ", "true       ", "trusting   ", &
      &"trustworthy", "trusty     ", "truthful   ", "tubby      ", "turbulent  ", &
      &"twin       "]

   character(len=*), parameter :: bWord(53) = [ &
      & "back      ", "backs     ", "bail      ", "balance   ", "balances  ", &
      & "balloon   ", "balloons  ", "ban       ", "bans      ", "bandage   ", &
      & "bandages  ", "bank      ", "bare      ", "bargain   ", "battle    ", &
      & "battles   ", "beam      ", "bear      ", "beat      ", "beats     ", &
      & "bend      ", "bends     ", "benefit   ", "benefits  ", "blame     ", &
      & "blast     ", "blasts    ", "bleach    ", "block     ", "bloom     ", &
      & "blow      ", "board     ", "boards    ", "bother    ", "bounce    ", &
      & "bounces   ", "bow       ", "bows      ", "box       ", "boxes     ", &
      & "bread     ", "break     ", "breaks    ", "breed     ", "broadcast ", &
      & "broadcasts", "brush     ", "brushes   ", "bump      ", "bumps     ", &
      & "burn      ", "burns     ", "buy       "]


contains


!> Generate a new random name
subroutine generateName(name)

   !> Generated name
   character(len=:), allocatable, intent(out) :: name

   real :: rNum(2)
   integer :: idx(2)

   call random_number(rNum)
   idx = int([size(tWord), size(bWord)] * rNum) + 1
   name = trim(tWord(idx(1))) // delimiter // trim(bWord(idx(2)))

end subroutine generateName


end module xtb_mctc_namegen
