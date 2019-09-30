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

module scanparam
   use iso_fortran_env, only : wp => real64
   use tbdef_setvar
   implicit none
   private :: wp
   public

   integer,parameter :: p_scan_sequential = 1
   integer,parameter :: p_scan_concerted = 2
   integer  :: scan_mode = p_scan_sequential

   integer  :: maxscan = 0
   integer  :: maxconstr = 0
   real(wp),allocatable :: valscan(:)
   real(wp),allocatable :: valconstr(:)
   real(wp) :: fcconstr = 0.05_wp
   real(wp) :: pmf ! whatever
   real(wp) :: springexpo = 2.0_wp
   integer  :: nscan = 0
   integer  :: nconstr = 0
   integer, allocatable :: atconstr(:,:)
   integer  :: iconstr = 0
   integer  :: zconstr = 0
   logical  :: lconstr_all_torsions = .false.
   logical  :: lconstr_all_bonds    = .false.
   logical  :: lconstr_all_angles   = .false.

   type(constr_setvar) :: potset
   integer, parameter  :: p_typeid_dist     = 1
   integer, parameter  :: p_typeid_angle    = 2
   integer, parameter  :: p_typeid_dihedral = 3
   integer, parameter  :: p_typeid_cma      = 4
   integer, parameter  :: p_typeid_zaxis    = 5

   type tb_scan
      integer  :: nscan = 0
      integer  :: iconstr = 0
      real(wp) :: valconstr = 0.0_wp
      real(wp),allocatable :: valscan(:)
   end type tb_scan

   type(tb_scan),allocatable :: scan_list(:)

contains

subroutine init_constr(n,at)
   implicit none
   integer,intent(in) :: n
   integer,intent(in) :: at(n)
   if (lconstr_all_bonds)    maxconstr = maxconstr+n*(1+n)/2
   if (lconstr_all_angles)   maxconstr = maxconstr+n*(1+n)*(2+n)/6
   if (lconstr_all_torsions) maxconstr = maxconstr+n*(1+n)*(2+n)*(3+n)/24
   call clear_constr
   allocate( valconstr(maxconstr),  source = 0.0_wp )
   valconstr = 0.0_wp
   allocate( atconstr(4,maxconstr), source = 0 )
   atconstr = 0

   call potset%allocate(n,maxconstr,fc=fcconstr,expo=springexpo)
!  new model
   if (allocated(potset%fname)) then
      allocate( potset%xyz(3,n), source = 0.0_wp )
      call read_reference(potset%fname,n,at,potset%xyz)
   endif
end subroutine init_constr

subroutine init_scan
   call clear_scan
   allocate( scan_list(maxscan) )
end subroutine init_scan

subroutine clear_scan
   if (allocated(scan_list)) deallocate(scan_list)
end subroutine clear_scan

subroutine clear_constr
   if (allocated(valconstr)) deallocate(valconstr)
   if (allocated(atconstr))  deallocate(atconstr)
end subroutine clear_constr

subroutine read_reference(fname,n,at,xyz)
   implicit none
   character(len=*),intent(in) :: fname
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(out) :: xyz(3,n)

   integer  :: n1
   integer, allocatable :: at1(:)

   n1 = n
   allocate( at1(n), source = 0 )

   call rdcoord(fname,n1,xyz,at1)

   if (n1.ne.n) &
      call raise('E',"Atom number missmatch in constraint reference geometry!",1)
   if (any(at.ne.at1)) &
      call raise('E',"Atom type missmatch in constraint reference geometry!",1)

end subroutine read_reference

subroutine setup_constrain_pot(n,at,xyz)
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: xyz(3,n)

   integer  :: i,ii,j,jj,ij
   real(wp) :: rij

   if (potset%pos%n.gt.0) then
      ij = 0
      do i = 1, potset%pos%n
         ii = potset%pos%atoms(i)
         do j = 1, i-1
            jj = potset%pos%atoms(j)
            ij = ij+1
            rij = norm2(xyz(1:3,jj)-xyz(1:3,ii))
            potset%pos%val(ij) = rij
         enddo
      enddo
      potset%pos%fc = potset%pos%fc / real(potset%pos%n-1,wp)
   endif

end subroutine setup_constrain_pot

subroutine pot_info(iunit,n,at,xyz)
   use mctc_constants
   use mctc_econv
   implicit none
   integer, intent(in)  :: iunit
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: xyz(3,n)

   character(len=2),external :: asym
   real(wp),external :: valijkl

   integer :: i,ii,j,k,l,m,mm
   real(wp) :: val0,val

   real(wp) :: e,efix,gfix
   real(wp),allocatable :: g(:,:)

   allocate( g(3,n), source = 0.0_wp )


   efix = 0.0_wp
   gfix = 0.0_wp

   if (potset%n.gt.0 .or. potset%pos%n.gt.0) then
      call generic_header(iunit,"Constraints",49,10)
      write(iunit,'(a)')
   endif
   if (potset%pos%n.gt.0) then
      write(iunit,'(1x,"*",1x,i0,1x,a)') potset%pos%n,"constrained atom positions"
      if (allocated(potset%xyz).and.allocated(potset%fname)) then
         write(iunit,'(3x,a,1x,"''",a,"''")') "positions refering to",&
            potset%fname
      else
         write(iunit,'(3x,a)') "positions refering to input geometry"
      endif
      write(iunit,'(a)')
      write(iunit,'(5x,"#",3x,"Z",3x,32x,"position/Å",6x,"displ./Å")')
      do i = 1, potset%pos%n
         ii = potset%pos%atoms(i)
         write(iunit,'(i6,1x,i3,1x,a2)',advance='no') ii,at(ii),asym(at(ii))
         if (allocated(potset%xyz)) then
            write(iunit,'(4f14.7)') &
               potset%xyz(:,ii)*autoaa, norm2(potset%xyz(:,ii)-xyz(:,ii))*autoaa
         else
            write(iunit,'(4f14.7)') &
               xyz(:,ii)*autoaa, 0.0_wp
         endif
      enddo
      write(iunit,'(a)')
      write(iunit,'(5x,a,1x,i0,1x,a)') "applying",potset%pos%n*(potset%pos%n-1)/2,&
         "atom pairwise harmonic potentials"
      write(iunit,'(5x,a,f14.7)') "  applied force constant per pair:", &
         potset%pos%fc
      write(iunit,'(5x,a,f14.7)') "effective force constant per atom:", &
         potset%pos%fc*(potset%pos%n-1)

      g = 0.0_wp
      e = 0.0_wp
      call constrain_pos(potset%pos,n,at,xyz,g,e)
      efix = efix+e
      gfix = gfix+norm2(g)
      write(iunit,'(5x,a,2f14.7)')"    constraining energy/grad norm:", &
         e,norm2(g)

      write(iunit,'(a)')
   endif

   if (potset%dist%n.gt.0) then
      write(iunit,'(1x,"*",1x,i0,1x,a)') potset%dist%n,"constrained distances"
      write(iunit,'(a)')
      write(iunit,'(2(5x,"#",3x,"Z",3x),26x,8x,"value/Å",6x,"actual/Å")')
      do m = 1, potset%dist%n
         mm = 2*m-1
         i = potset%dist%atoms(mm)
         j = potset%dist%atoms(mm+1)
         val0 = potset%dist%val(m)
         val = norm2(xyz(:,i)-xyz(:,j))
         write(iunit,'(2(i6,1x,i3,1x,a2),26x,1x,2f14.7)') &
            i,at(i),asym(at(i)), &
            j,at(j),asym(at(j)), &
            val0*autoaa, val*autoaa
      enddo
      write(iunit,'(a)')
      write(iunit,'(5x,a,f14.7)') "  constraining potential exponent:", &
         potset%dist%expo
      write(iunit,'(5x,a,f14.7)') " applied force constant per dist.:", &
         potset%dist%fc
      write(iunit,'(5x,a,f14.7)') "effective force constant per atom:", &
         potset%dist%fc/2.0_wp

      g = 0.0_wp
      e = 0.0_wp
      call constrain_dist(potset%dist,n,at,xyz,g,e)
      efix = efix+e
      gfix = gfix+norm2(g)
      write(iunit,'(5x,a,2f14.7)')"    constraining energy/grad norm:", &
         e,norm2(g)

      write(iunit,'(a)')
   endif

   if (potset%angle%n.gt.0) then
      write(iunit,'(1x,"*",1x,i0,1x,a)') potset%angle%n,"constrained angles"
      write(iunit,'(a)')
      write(iunit,'(3(5x,"#",3x,"Z",3x),13x,8x,"value/°",6x,"actual/°")')
      do m = 1, potset%angle%n
         mm = 3*m-2
         i = potset%angle%atoms(mm)
         j = potset%angle%atoms(mm+1)
         k = potset%angle%atoms(mm+2)
         val0 = potset%angle%val(m)
         call bangl(xyz,i,j,k,val)
         write(iunit,'(3(i6,1x,i3,1x,a2),13x,1x,2f14.7)') &
            i,at(i),asym(at(i)), &
            j,at(j),asym(at(j)), &
            k,at(k),asym(at(k)), &
            val0*180.0_wp/pi, val*180.0_wp/pi
      enddo
      write(iunit,'(a)')
      write(iunit,'(5x,a,f14.7)') " applied force constant per angle:", &
         potset%angle%fc
      write(iunit,'(5x,a,f14.7)') "effective force constant per atom:", &
         potset%angle%fc/3.0_wp

      g = 0.0_wp
      e = 0.0_wp
      call constrain_angle(potset%angle,n,at,xyz,g,e)
      efix = efix+e
      gfix = gfix+norm2(g)
      write(iunit,'(5x,a,2f14.7)')"    constraining energy/grad norm:", &
         e,norm2(g)

      write(iunit,'(a)')
   endif

   if (potset%dihedral%n.gt.0) then
      write(iunit,'(1x,"*",1x,i0,1x,a)') potset%dihedral%n, &
         "constrained dihedral angles"
      write(iunit,'(a)')
      write(iunit,'(4(5x,"#",3x,"Z",3x),8x,"value/°",6x,"actual/°")')
      do m = 1, potset%dihedral%n
         mm = 4*m-3
         i = potset%dihedral%atoms(mm)
         j = potset%dihedral%atoms(mm+1)
         k = potset%dihedral%atoms(mm+2)
         l = potset%dihedral%atoms(mm+3)
         val0 = potset%dihedral%val(m)
         val = valijkl(n,xyz,i,j,k,l)
         write(iunit,'(4(i6,1x,i3,1x,a2),1x,2f14.7)') &
            i,at(i),asym(at(i)), &
            j,at(j),asym(at(j)), &
            k,at(k),asym(at(k)), &
            l,at(l),asym(at(l)), &
            val0*180.0_wp/pi, val*180.0_wp/pi
      enddo
      write(iunit,'(a)')
      write(iunit,'(5x,a,f14.7)') " applied force constant per angle:", &
         potset%dihedral%fc
      write(iunit,'(5x,a,f14.7)') "effective force constant per atom:", &
         potset%dihedral%fc/4.0_wp

      g = 0.0_wp
      e = 0.0_wp
      call constrain_dihedral(potset%dihedral,n,at,xyz,g,e)
      efix = efix+e
      gfix = gfix+norm2(g)
      write(iunit,'(5x,a,2f14.7)')"    constraining energy/grad norm:", &
         e,norm2(g)

      write(iunit,'(a)')
   endif
   if (potset%n.gt.0 .or. potset%pos%n.gt.0) then
      write(iunit,'(5x,a,2f14.7)')"total constraint energy/grad norm:", &
         efix, gfix
      write(iunit,'(a)')
   endif

end subroutine pot_info

subroutine constrain_all_bonds(n,at,xyz)
   use iso_fortran_env, wp => real64
   use mctc_param
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: xyz(3,n)

   integer  :: i,j
   integer  :: ioffset
   real(wp) :: r,r0
   real(wp),parameter :: f = 1.2_wp

   do i = 1, n
      do j = 1, i-1
         r  = norm2(xyz(:,i)-xyz(:,j))
         r0 = rad(at(j))+rad(at(i))
         if (r.lt.f*r0) then
            ioffset = 2*potset%dist%n
            potset%dist%n = potset%dist%n+1
            potset%dist%atoms(ioffset+1) = j
            potset%dist%atoms(ioffset+2) = i
            potset%dist%val(potset%dist%n) = r
         endif
      enddo
   enddo

end subroutine constrain_all_bonds

subroutine constrain_all_angles(n,at,xyz)
   use iso_fortran_env, wp => real64
   use mctc_constants
   use mctc_param
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: xyz(3,n)

   integer  :: i,j,k
   integer  :: ioffset
   real(wp) :: phi
   integer, allocatable :: bond(:,:)
   real(wp),parameter   :: thr = 0.2_wp

   allocate( bond(n,n), source = 0 )
   call get_bonds(n,at,xyz,bond)

   do i = 1, n
      if (bond(i,i).lt.2) cycle
      do j = 1, i-1
         if (i.eq.j) cycle
         if (bond(i,i).lt.2) cycle
         if (bond(j,i).lt.1) cycle
         do k = 1, j-1
            if (i.eq.k .or. j.eq.k) cycle
            if (bond(k,j).lt.1) cycle
            call bangl(xyz,k,j,i,phi)
            if (abs(pi-phi).lt.thr) cycle
            ioffset = 3*potset%angle%n
            potset%angle%n = potset%angle%n+1
            potset%angle%atoms(ioffset+1) = k
            potset%angle%atoms(ioffset+2) = j
            potset%angle%atoms(ioffset+3) = i
            potset%angle%val(potset%angle%n) = phi
         enddo
      enddo
   enddo

end subroutine constrain_all_angles

subroutine constrain_all_torsions(n,at,xyz)
   use iso_fortran_env, wp => real64
   use mctc_constants
   use mctc_param
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: xyz(3,n)

   integer  :: i,j,k,l
   integer  :: ioffset
   real(wp) :: phi,thijk,thjkl
   integer, allocatable :: bond(:,:)
   real(wp),parameter   :: thr = 0.2_wp
   real(wp),external    :: valijkl

   allocate( bond(n,n), source = 0 )
   call get_bonds(n,at,xyz,bond)

   do i = 1, n
      if (bond(i,i).lt.2) cycle
      do j = 1, n
         if (i.eq.j) cycle
         if (bond(j,i).lt.1 .or. bond(j,j).lt.2) cycle
         do k = 1, n
            if (i.eq.k .or. j.eq.k) cycle
            if (bond(k,j).lt.1 .or. bond(k,k).lt.2) cycle
            do l = 1, n
               if (bond(l,k).lt.1) cycle
               if (i.eq.l .or. j.eq.l .or. k.eq.l) cycle
               call bangl(xyz,k,j,i,thijk)
               call bangl(xyz,l,k,j,thjkl)
               if (abs(pi-thijk).lt.thr .or. abs(pi-thjkl).lt.thr) cycle
               phi = valijkl(n,xyz,l,k,j,i)
               ioffset = 4*potset%dihedral%n
               potset%dihedral%n = potset%dihedral%n+1
               potset%dihedral%atoms(ioffset+1) = l
               potset%dihedral%atoms(ioffset+2) = k
               potset%dihedral%atoms(ioffset+3) = j
               potset%dihedral%atoms(ioffset+4) = i
               potset%dihedral%val(potset%dihedral%n) = phi
            enddo
         enddo
      enddo
   enddo

end subroutine constrain_all_torsions

end module scanparam
