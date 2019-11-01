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

module iniq
use iso_fortran_env, only : wp => real64
implicit none

contains

!! ========================================================================
!  driver for calculation of Grimme-Gasteiger charges, which come in 
!  three flavors: GG(vTB), GG(GFN1) and GG(GFN2). They differ in the way
!  the coordination number is taken into account.
!  INPUT
!  nat  :: number of atoms
!  nel  :: number of electrons
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  z    :: nuclear charges
!  chrg :: system charge
!  OUTPUT
!  q    :: partial charges
!  cn   :: coordination number
!  PARAMETER
!  kchrg1  :: for GFN1-type Grimme-Gasteiger charges
!  OPTIONS
!  version :: vTB(=0),GFN1(=1),GFN2(>1)
!  verbose :: prints some header and makes a bit of analysis in the end
subroutine iniqcn(nat,nel,at,z,xyz,chrg,kchrg1,q,cn,version,verbose)
   use iso_fortran_env, only : id => output_unit

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: nel
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: z(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   integer, intent(in)  :: chrg
   real(wp),intent(in)  :: kchrg1
   real(wp),intent(out) :: q(nat)
   real(wp),intent(out) :: cn(nat)
   integer, intent(in)  :: version
   logical, intent(in),optional :: verbose

   integer  :: i
   real(wp) :: dipole(3)

   logical :: pr
   if (present(verbose)) then
      pr = verbose
   else
      pr = .false.
   endif

!  print some header if requested
   if (pr) then
      write(id,'(a)')
      write(id,'(''doing Grimme-Gasteiger EN-charges'')')
   endif

!  select the variant of Grimme-Gasteiger charges based on the version number
   if (version.le.0) then
      call iniqcn_vtb(nat,nel,at,z,xyz,chrg,q,cn)
   else if (version.eq.1) then
      if (pr) write(id,'(''using D3 CN'')')
      call iniqcn_gfn1(nat,nel,at,z,xyz,chrg,kchrg1,q,cn)
   else
      if (pr) write(id,'(''using GFN CN'')')
      call iniqcn_gfn2(nat,nel,at,z,xyz,chrg,q,cn)
   endif

 
!  calculate the dipole moment of the guess charges if requested
   if (pr) then
      dipole = 0.0_wp
      do i = 1, nat
         dipole(1)=dipole(1)+xyz(1,i)*q(i)
         dipole(2)=dipole(2)+xyz(2,i)*q(i)
         dipole(3)=dipole(3)+xyz(3,i)*q(i)
      enddo
      write(id,'(1x,''sum q :'',1x,d14.7)') sum(q)
      write(id,'(1x,''point charge moment (au)'')')
      write(id,'(5x,''X'',7x,''Y'',7x,''Z'')')
      write(id,'(3f9.4,''  total (Debye): '',f8.3)') &
      &  dipole(1), dipole(2), dipole(3), &
      &  sqrt(sum(dipole**2))*2.5418_wp
   endif

end subroutine iniqcn

!! ========================================================================
!  non-iterative, distance dependent Gasteiger-Grimme type charges
!  for the vTB part in sTDA-xTB calculations
!  INPUT
!  nat  :: number of atoms
!  nel  :: number of electrons
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  z    :: nuclear charges
!  en   :: atomic EN (taken from aoparam)
!  chrg :: system charge
!  OUTPUT
!  q    :: partial charges
!  cn   :: coordination number
pure subroutine iniqcn_vtb(nat,nel,at,z,xyz,chrg,q,cn)

!  get data from parameter modules
   use aoparam, only : en,metal

!  get interface to ncoord
   use ncoord, only : ncoord_d3

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: nel
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: z(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   integer, intent(in)  :: chrg
   real(wp),intent(out) :: q(nat)
   real(wp),intent(out) :: cn(nat)

   integer :: i
   real(wp),allocatable :: ena(:)

!  in case it is a single atom, quick return
   if (nat.eq.1) then
      q(1) = real(chrg,wp)
      cn(1) = 0.0_wp
      return
   endif

   allocate( ena(nat), source = 0.0_wp )

!  get the coordination number
   call ncoord_d3(nat,at,xyz,cn)

!  uncorrected electronegativites
   do i = 1, nat
      ena(i) = en(at(i))
      if (metal(at(i)).gt.0) cn(i) = 0.0_wp
   enddo

!  the neutral part, start from nuclear charges
   q = z

!  now distribute the electrons ena dependent
   call gasteiger_partition(nat,at,xyz,ena,q)

!  normalize to Nel
   q = q * (sum(z)-real(chrg,wp))/sum(z)

!  make partial charges
   q = z - q

end subroutine iniqcn_vtb

!! ========================================================================
!  non-iterative, distance dependent Gasteiger-Grimme type charges
!  for GFN1-xTB calculations
!  INPUT
!  nat  :: number of atoms
!  nel  :: number of electrons
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  z    :: nuclear charges
!  en   :: atomic EN (taken from aoparam)
!  chrg :: system charge
!  OUTPUT
!  q    :: partial charges
!  cn   :: coordination number
!  PARAMETER:
!  kchrg1 (via INPUT)
pure subroutine iniqcn_gfn1(nat,nel,at,z,xyz,chrg,kchrg1,q,cn)

!  get data from parameter modules
   use aoparam, only : en,metal

!  get interface to ncoord
   use ncoord, only : ncoord_d3

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: nel
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: z(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   integer, intent(in)  :: chrg
   real(wp),intent(in)  :: kchrg1
   real(wp),intent(out) :: q(nat)
   real(wp),intent(out) :: cn(nat)

   integer :: i
   real(wp),allocatable :: ena(:)

!  in case it is a single atom, quick return
   if (nat.eq.1) then
      q(1) = real(chrg,wp)
      cn(1) = 0.0_wp
      return
   endif

   allocate( ena(nat), source = 0.0_wp )

!  get the coordination number
   call ncoord_d3(nat,at,xyz,cn)

!  correct the CN with effective coordination numbers
   do i = 1, nat
      if (metal(at(i)).gt.0) then
         ena(i) = 0.0_wp
      else
         ena(i) = en(at(i)) - kchrg1*sqrt(cn(i))
      endif
   enddo

!  the neutral part, start from nuclear charges
   q = z

!  now distribute the electrons ena dependent
   call gasteiger_partition(nat,at,xyz,ena,q)

!  normalize to Nel
   q = q * (sum(z)-real(chrg,wp))/sum(z)

!  make partial charges
   q = z - q

end subroutine iniqcn_gfn1

!! ========================================================================
!  non-iterative, distance dependent Gasteiger-Grimme type charges
!  for GFN2-xTB calculations
!  INPUT
!  nat  :: number of atoms
!  nel  :: number of electrons
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  z    :: nuclear charges
!  en   :: atomic EN (taken from aoparam)
!  chrg :: system charge
!  OUTPUT
!  q    :: partial charges
!  cn   :: coordination number
pure subroutine iniqcn_gfn2(nat,nel,at,z,xyz,chrg,q,cn)

!  get data from parameter modules
   use aoparam, only : en,metal

!  get interface to ncoord
   use ncoord, only : ncoord_gfn

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: nel
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: z(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   integer, intent(in)  :: chrg
   real(wp),intent(out) :: q(nat)
   real(wp),intent(out) :: cn(nat)

   real(wp),parameter   :: kchrg2 = 0.2_wp

   integer :: i
   real(wp),allocatable :: ena(:)

!  in case it is a single atom, quick return
   if (nat.eq.1) then
      q(1) = real(chrg,wp)
      cn(1) = 0.0_wp
      return
   endif

   allocate( ena(nat), source = 0.0_wp )

!  get the coordination number
   call ncoord_gfn(nat,at,xyz,cn)

!  correct the CN with effective coordination numbers
   do i = 1, nat
      if (metal(at(i)).gt.0) then
         ena(i) = kchrg2*sqrt(cn(i)) ! avoid to big q for metals
      else
         ena(i) = en(at(i)) + kchrg2*sqrt(cn(i))
      endif
   enddo

!  the neutral part, start from nuclear charges
   q = z

!  now distribute the electrons ena dependent
   call gasteiger_partition(nat,at,xyz,ena,q)

!  normalize to Nel
   q = q * (sum(z)-real(chrg,wp))/sum(z)

!  make partial charges
   q = z - q

end subroutine iniqcn_gfn2

!! ========================================================================
!  non-iterative, distance dependent Gasteiger-Grimme type charges
!  INPUT
!  nat  :: number of atoms
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  en   :: systemspecific EN
!  OUTPUT
!  q    :: EN based population
pure subroutine gasteiger_partition(nat,at,xyz,en,q)

!  get data from parameter modules
   use dftd4, only : rcov

   implicit none

   integer, intent(in)    :: nat
   integer, intent(in)    :: at(nat)
   real(wp),intent(in)    :: xyz(3,nat)
   real(wp),intent(in)    :: en(nat)
   real(wp),intent(inout) :: q(nat)

   integer  :: i,j
   real(wp) :: tmp,r2,rav,den,rexp

   do i = 1, nat
      tmp = 0.0_wp
      do j = 1, nat
         if (i.eq.j) cycle
         r2 = sqrt(sum( (xyz(:,j) - xyz(:,i))**2 ))
         rav = 0.5_wp * ( rcov(at(i))+rcov(at(j)) )
         den = en(i) - en(j)
         rexp = (rav/r2)**6
         tmp = tmp + den*rexp
      enddo
      q(i) = q(i) + tmp
   enddo

end subroutine gasteiger_partition
end module iniq
