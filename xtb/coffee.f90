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

subroutine get_COFFEE(mol)
   use iso_fortran_env, only : wp => real64
   use tbdef_molecule
   type(tb_molecule), intent(inout) :: mol
   integer, parameter :: nat = 24
   integer, parameter :: at(nat) = [6,7,6,7,6,6,6,8,7,6,8,7,6,6, &
      &                             1,1,1,1,1,1,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[ 2.02799738646442_wp,  0.09231312124713_wp, -0.14310895950963_wp, &
      &  4.75011007621000_wp,  0.02373496014051_wp, -0.14324124033844_wp, &
      &  6.33434307654413_wp,  2.07098865582721_wp, -0.14235306905930_wp, &
      &  8.72860718071825_wp,  1.38002919517619_wp, -0.14265542523943_wp, &
      &  8.65318821103610_wp, -1.19324866489847_wp, -0.14231527453678_wp, &
      &  6.23857175648671_wp, -2.08353643730276_wp, -0.14218299370797_wp, &
      &  5.63266886875962_wp, -4.69950321056008_wp, -0.13940509630299_wp, &
      &  3.44931709749015_wp, -5.48092386085491_wp, -0.14318454855466_wp, &
      &  7.77508917214346_wp, -6.24427872938674_wp, -0.13107140408805_wp, &
      & 10.30229550927022_wp, -5.39739796609292_wp, -0.13672168520430_wp, &
      & 12.07410272485492_wp, -6.91573621641911_wp, -0.13666499342053_wp, &
      & 10.70038521493902_wp, -2.79078533715849_wp, -0.14148379504141_wp, &
      & 13.24597858727017_wp, -1.76969072232377_wp, -0.14218299370797_wp, &
      &  7.40891694074004_wp, -8.95905928176407_wp, -0.11636933482904_wp, &
      &  1.38702118184179_wp,  2.05575746325296_wp, -0.14178615122154_wp, &
      &  1.34622199478497_wp, -0.86356704498496_wp,  1.55590600570783_wp, &
      &  1.34624089204623_wp, -0.86133716815647_wp, -1.84340893849267_wp, &
      &  5.65596919189118_wp,  4.00172183859480_wp, -0.14131371969009_wp, &
      & 14.67430918222276_wp, -3.26230980007732_wp, -0.14344911021228_wp, &
      & 13.50897177220290_wp, -0.60815166181684_wp,  1.54898960808727_wp, &
      & 13.50780014200488_wp, -0.60614855212345_wp, -1.83214617078268_wp, &
      &  5.41408424778406_wp, -9.49239668625902_wp, -0.11022772492007_wp, &
      &  8.31919801555568_wp, -9.74947502841788_wp,  1.56539243085954_wp, &
      &  8.31511620712388_wp, -9.76854236502758_wp, -1.79108242206824_wp],&
      &  shape(xyz))

   call mol%allocate(nat)
   mol%at   = at
   mol%xyz  = xyz
   mol%chrg = -1

   call mol%update

end subroutine get_COFFEE

subroutine get_CUP(iunit,fmt)
   implicit none
   integer, intent(in) :: iunit
   character(len=*),intent(in) :: fmt
   write(iunit,fmt) &
      & "            )", &
      & "         ) (", &
      & "        (   ) (", &
      & "     _..,)----,)._", &
      & "  .-;'-.,_____,.-';", &
      & " (( |             |", &
      & "  `))             ;", &
      & "   ` \           /", &
      & "  .-' `,._____.,' '-.", &
      & " (     '-------'     )", &
      & "  `-=.._________..--'"
end subroutine get_CUP
