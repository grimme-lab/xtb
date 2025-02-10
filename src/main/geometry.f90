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
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

subroutine main_geometry(iunit,mol)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule
   implicit none
   integer, intent(in)  :: iunit
   type(TMolecule), intent(in) :: mol

   if (mol%npbc == 0) then
      call print_geosum(iunit,mol%n,mol%at,mol%sym,mol%xyz)
   else
      call print_pbcsum(iunit,mol)
   endif

end subroutine main_geometry

subroutine print_pbcsum(iunit,mol)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_mctc_convert
   use xtb_type_molecule
   implicit none
   integer, intent(in)  :: iunit
   type(TMolecule),intent(in) :: mol

   integer  :: i
   real(wp) :: conv

   conv = autoaa

   call generic_header(iunit,"Geometry Summary",49,10)
   write(iunit,'(a)')

   ! atomic coordinates
   write(iunit,'(1x,"*",1x,i0,1x,a)') mol%n,"atoms in unit cell"
   write(iunit,'(a)')
   write(iunit,'(5x,"#",3x,"Z",5x,32x,"position/Å",8x,"charge")')
   do i = 1, mol%n
      write(iunit,'(i6,1x,i3,1x,a4)',advance='no') i,mol%at(i),mol%sym(i)
      write(iunit,'(3f14.7)',advance='no') mol%xyz(:,i)*conv
      write(iunit,'(f14.7)') mol%z(i)
   enddo
   write(iunit,'(a)')

   ! periodicity
   write(iunit,'(1x,"*",1x,i0,a)') mol%npbc,"D periodic system"
   write(iunit,'(a)')

   if (mol%npbc > 0) then
      ! cell parameters
      write(iunit,'(1x,"*",1x,a)') "cell parameter"
      write(iunit,'(a)')
      write(iunit,'(a12,2a15,2x,3a11)') &
         "|a|/Å", "|b|/Å", "|c|/Å", "α/°", "β/°", "γ/°"
      write(iunit,'(f13.7,2f14.7,1x,3f9.3)') &
         mol%cellpar(:3)*conv,mol%cellpar(4:)*180.0_wp/pi
      write(iunit,'(a)')

      ! direct lattice (transformation abc -> xyz)
      write(iunit,'(1x,"*",1x,a)') "direct lattice/Å"
      write(iunit,'(a)')
      write(iunit,'(12x,a,3f14.7)') "a",mol%lattice(:,1)*conv
      write(iunit,'(12x,a,3f14.7)') "b",mol%lattice(:,2)*conv
      write(iunit,'(12x,a,3f14.7)') "c",mol%lattice(:,3)*conv
      write(iunit,'(a)')

      ! reciprocal lattice
      write(iunit,'(1x,"*",1x,a)') "reciprocal lattice/Å⁻¹"
      write(iunit,'(a)')
      write(iunit,'(11x,a,3f14.7)') "a*",mol%rec_lat(:,1)/conv
      write(iunit,'(11x,a,3f14.7)') "b*",mol%rec_lat(:,2)/conv
      write(iunit,'(11x,a,3f14.7)') "c*",mol%rec_lat(:,3)/conv
      write(iunit,'(a)')

      ! geometry in fractional coordinates
      write(iunit,'(1x,"*",1x,a)') "geometry in fractional coordinates"
      write(iunit,'(a)')
      write(iunit,'(5x,"#",3x,"Z",5x,20x,"fractional coordinates")')
      do i = 1, mol%n
         write(iunit,'(i6,1x,i3,1x,a4)',advance='no') i,mol%at(i),mol%sym(i)
         write(iunit,'(3f14.7)',advance='no') mol%abc(:,i)
         write(iunit,'(a)')
      enddo
      write(iunit,'(a)')

      ! volume of unit cell
      write(iunit,'(1x,"*",1x,a,1x,"=",f14.7)') "volume of direct unit cell/Å³", &
         mol%volume*conv**3
      write(iunit,'(a)')

      ! mass density of unit cell
      write(iunit,'(1x,"*",1x,a,1x,"=",f14.7)') "mass density of direct unit cell/g·cm⁻³", &
         sum(mol%atmass)*metokg*1.0e27_wp/(mol%volume*conv**3)
      write(iunit,'(a)')
   endif

end subroutine print_pbcsum


subroutine print_geosum(iunit,n,at,sym,xyz)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_mctc_convert
   use xtb_setparam
   use xtb_splitparam
   use xtb_disp_ncoord
   use xtb_approxrab
   implicit none
   integer, intent(in)  :: iunit
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   character(len=*), intent(in) :: sym(n)
   real(wp),intent(in)  :: xyz(3,n)

   integer  :: i,j,k,l,m
   integer  :: iat,jat
   integer  :: imax,imin

   integer(8) :: maxdist,maxbend,maxtrsn
   integer    :: maxelem
   integer    :: ndist,nbend,ntrsn,nelem
   integer, allocatable :: bond(:,:)
   integer, allocatable :: id(:,:)
   real(wp),allocatable :: dist(:)
   integer, allocatable :: ndel(:,:)
   real(wp),allocatable :: distel(:,:)
   real(wp),allocatable :: maxdel(:,:)
   real(wp),allocatable :: mindel(:,:)
   integer, allocatable :: ib(:,:)
   real(wp),allocatable :: bend(:)
   integer, allocatable :: it(:,:)
   real(wp),allocatable :: trsn(:)
   real(wp),allocatable :: cn(:)
   real(wp) :: moments(3)

   maxdist = n*(n-1)/2
   maxbend = n*(n-1)*(n-2)/6
   maxtrsn = n*(n-1)*(n-2)*(n-3)/24
   maxelem = min(maxval(at),118)

   allocate( cn(n), source = 0.0_wp )
   allocate( bond(n,n), source = 0 )
   call ncoord_erf(n,at,xyz,cn,900.0_wp)
   call approx_bonds(n,at,xyz,cn,bond,0.0537_wp)
   !call get_bonds(n,at,xyz,bond)

   if (set%pr_moments.or.set%pr_distances.or.set%pr_angles.or.set%pr_torsions) then
      call generic_header(iunit,"Geometry Summary",49,10)
      write(iunit,'(a)')
   endif

   if (set%pr_moments) then
      call print_moments(iunit,n,atmass,xyz)
   endif

   if (set%pr_distances) then
      allocate( dist(maxdist),distel(maxelem,maxelem), &
         maxdel(maxelem,maxelem),mindel(maxelem,maxelem), source = 0.0_wp )
      allocate( id(2,maxdist),ndel(maxelem,maxelem), source = 0 )
      call calc_distances(n,at,xyz,bond,maxdist,ndist,dist,id, &
         maxelem,nelem,ndel,distel,maxdel,mindel)
      if (ndist.gt.0) then
         call print_distances(iunit,n,at,sym,ndist,dist,id)
         call print_elem_dist(iunit,maxelem,nelem,ndel,distel,maxdel,mindel)
      endif
      deallocate( dist,distel,maxdel,mindel,id,ndel)
   endif

   if (set%pr_angles) then
      allocate( bend(maxbend), source = 0.0_wp )
      allocate( ib(3,maxbend), source = 0 )
      call calc_angles(n,at,xyz,bond,maxbend,nbend,bend,ib)
      if (nbend.gt.0) then
         call print_angles(iunit,n,at,sym,nbend,bend,ib)
      endif
      deallocate( bend,ib )
   endif

   if (set%pr_torsions) then
      allocate( trsn(maxtrsn), source = 0.0_wp )
      allocate( it(4,maxtrsn), source = 0 )
      call calc_torsions(n,at,xyz,bond,maxtrsn,ntrsn,trsn,it)
      if (ntrsn.gt.0) then
         call print_torsions(iunit,n,at,sym,ntrsn,trsn,it)
      endif
      deallocate( trsn,it )
   endif

end subroutine print_geosum

subroutine print_elem_dist(iunit,n,nbond,bond,dist,maxd,mind)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   use xtb_mctc_symbols, only : toSymbol
   implicit none
   integer, intent(in)  :: iunit
   integer, intent(in)  :: n
   integer, intent(in)  :: nbond
   integer, intent(in)  :: bond(n,n)
   real(wp),intent(in)  :: dist(n,n)
   real(wp),intent(in)  :: maxd(n,n)
   real(wp),intent(in)  :: mind(n,n)

   integer :: iat,jat

   write(iunit,'(1x,"*",1x,i0,1x,a)') nbond,"distinct bonds (by element types)"
   write(iunit,'(a)')
   write(iunit,'(2(3x,"Z",3x),10x,"#",3x,"av. dist./Å",8x,"max./Å",8x,"min./Å")')
   do iat = 1, n
      do jat = 1, iat
         if (bond(jat,iat).lt.1) cycle
         write(iunit,'(2(1x,i3,1x,a2),i11,3f14.7)') &
            jat,toSymbol(jat),iat,toSymbol(iat),bond(jat,iat),dist(jat,iat)*autoaa, &
            maxd(jat,iat)*autoaa,mind(jat,iat)*autoaa
      enddo
   enddo
   write(iunit,'(a)')

end subroutine print_elem_dist

subroutine print_distances(iunit,n,at,sym,ndist,dist,id)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   use xtb_mctc_symbols, only : toSymbol
   implicit none
   integer, intent(in)  :: iunit
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   character(len=*), intent(in) :: sym(n)
   integer, intent(in)  :: ndist
   real(wp),intent(in)  :: dist(ndist)
   integer, intent(in)  :: id(2,ndist)

   integer :: i,j,m
   integer :: imax,imin

   write(iunit,'(1x,"*",1x,i0,1x,a)') ndist,"selected distances"
   write(iunit,'(a)')
   write(iunit,'(2(5x,"#",3x,"Z",5x),30x,8x,"value/Å")')
   imax = maxloc(dist,1)
   imin = minloc(dist,1)
   do m = 1, ndist
      i = id(1,m)
      j = id(2,m)
      write(iunit,'(2(i6,1x,i3,1x,a4),30x,1x,f14.7)',advance='no') &
         i,at(i),sym(i), &
         j,at(j),sym(j), &
         dist(m)*autoaa
      if (imax.eq.m) then
         write(iunit,'(1x,"(max)")')
      elseif (imin.eq.m) then
         write(iunit,'(1x,"(min)")')
      else
         write(iunit,'(a)')
      endif
   enddo
   write(iunit,'(a)')

end subroutine print_distances

subroutine print_angles(iunit,n,at,sym,nbend,bend,ib)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_mctc_symbols, only : toSymbol
   implicit none
   integer, intent(in)  :: iunit
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   character(len=*), intent(in) :: sym(n)
   integer, intent(in)  :: nbend
   real(wp),intent(in)  :: bend(nbend)
   integer, intent(in)  :: ib(3,nbend)

   integer :: i,j,k,m
   integer :: imax,imin

   write(iunit,'(1x,"*",1x,i0,1x,a)') nbend,"selected angles"
   write(iunit,'(a)')
   write(iunit,'(3(5x,"#",3x,"Z",5x),13x,10x,"value/°")')
   do m = 1, nbend
      i = ib(1,m)
      j = ib(2,m)
      k = ib(3,m)
      write(iunit,'(3(i6,1x,i3,1x,a4),13x,1x,f14.7)') &
         i,at(i),sym(i), &
         j,at(j),sym(j), &
         k,at(k),sym(k), &
         bend(m)*180.0_wp/pi
   enddo
   write(iunit,'(a)')

end subroutine print_angles

subroutine print_torsions(iunit,n,at,sym,ntrsn,trsn,it)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_mctc_symbols, only : toSymbol
   implicit none
   integer, intent(in)  :: iunit
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   character(len=*), intent(in) :: sym(n)
   integer, intent(in)  :: ntrsn
   real(wp),intent(in)  :: trsn(ntrsn)
   integer, intent(in)  :: it(4,ntrsn)

   integer :: i,j,k,l,m
   integer :: imax,imin

   write(iunit,'(1x,"*",1x,i0,1x,a)') ntrsn, "selected dihedral angles"
   write(iunit,'(a)')
   write(iunit,'(4(5x,"#",3x,"Z",5x),8x,"value/°")')
   do m = 1, ntrsn
      i = it(1,m)
      j = it(2,m)
      k = it(3,m)
      l = it(4,m)
      write(iunit,'(4(i6,1x,i3,1x,a4),1x,f14.7)') &
         i,at(i),sym(i), &
         j,at(j),sym(j), &
         k,at(k),sym(k), &
         l,at(l),sym(l), &
         trsn(m)*180.0_wp/pi
   enddo
   write(iunit,'(a)')

end subroutine print_torsions

subroutine calc_distances(n,at,xyz,bond,maxdist,ndist,dist,id, &
      maxelem,nelem,ndel,distel,maxdel,mindel)
   use xtb_mctc_accuracy, only : wp
   implicit none
   integer,    intent(in)  :: n
   integer,    intent(in)  :: at(n)
   real(wp),   intent(in)  :: xyz(3,n)
   integer,    intent(in)  :: bond(n,n)
   integer(8), intent(in)  :: maxdist
   integer,    intent(out) :: ndist
   real(wp),   intent(out) :: dist(maxdist)
   integer,    intent(out) :: id(2,maxdist)
   integer,    intent(in)  :: maxelem
   integer,    intent(out) :: nelem
   integer,    intent(out) :: ndel(maxelem,maxelem)
   real(wp),   intent(out) :: distel(maxelem,maxelem)
   real(wp),   intent(out) :: maxdel(maxelem,maxelem)
   real(wp),   intent(out) :: mindel(maxelem,maxelem)

   integer  :: i,j,m
   integer  :: iat,jat
   real(wp) :: r

   ndist = 0
   nelem = 0
   id = 0
   dist = 0.0_wp
   ndel = 0
   distel = 0.0_wp
   maxdel = -1.0e42_wp
   mindel = +1.0e42_wp

   m = 0
   get_dist: do i = 1, n
      iat = at(i)
      do j = 1, i-1
         jat = at(j)
         r  = sqrt(sum((xyz(:,i)-xyz(:,j))**2))
         if (bond(j,i).gt.0) then
            if (m.ge.maxdist) exit get_dist
            m = m+1
            ndel(iat,jat) = ndel(iat,jat)+1
            maxdel(iat,jat) = max(maxdel(iat,jat),r)
            mindel(iat,jat) = min(mindel(iat,jat),r)
            distel(iat,jat) = distel(iat,jat)+r
            id(1,m) = j
            id(2,m) = i
            dist(m) = r
         endif
      enddo
   enddo get_dist
   ndist = m
   m = 0
   do iat = 1, maxelem
      do jat = 1, iat
         if (iat.ne.jat) then
            ndel(jat,iat) = ndel(jat,iat)+ndel(iat,jat)
            ndel(iat,jat) = ndel(jat,iat)
         endif
         if (ndel(jat,iat).lt.1) cycle
         m = m+1
         distel(jat,iat) = (distel(iat,jat)+distel(jat,iat)) &
            /real(ndel(jat,iat),wp)
         if (iat.eq.jat) distel(jat,iat) = distel(jat,iat)*0.5_wp
         distel(iat,jat) = distel(jat,iat)
         maxdel(jat,iat) = max(maxdel(iat,jat),maxdel(jat,iat))
         maxdel(iat,jat) = maxdel(jat,iat)
         mindel(jat,iat) = min(mindel(iat,jat),mindel(jat,iat))
         mindel(iat,jat) = mindel(jat,iat)
      enddo
   enddo
   nelem = m

end subroutine calc_distances

subroutine calc_angles(n,at,xyz,bond,maxbend,nbend,bend,ib)
   use xtb_mctc_accuracy, only : wp
   implicit none
   integer,    intent(in)  :: n
   integer,    intent(in)  :: at(n)
   real(wp),   intent(in)  :: xyz(3,n)
   integer,    intent(in)  :: bond(n,n)
   integer(8), intent(in)  :: maxbend
   integer,    intent(out) :: nbend
   real(wp),   intent(out) :: bend(maxbend)
   integer,    intent(out) :: ib(3,maxbend)

   integer  :: i,j,k,m
   real(wp) :: phi

   nbend = 0
   ib = 0
   bend = 0.0_wp

   m = 0
   get_bend: do i = 1, n
      if (bond(i,i).lt.2) cycle
      do j = 1, i-1
         if (i.eq.j) cycle
         if (bond(i,i).lt.2) cycle
         if (bond(j,i).lt.1) cycle
         do k = 1, j-1
            if (i.eq.k .or. j.eq.k) cycle
            if (bond(k,j).lt.1) cycle
            call bangl(xyz,k,j,i,phi)
            if (m.ge.maxbend) exit get_bend
            m = m+1
            ib(1,m) = k
            ib(2,m) = j
            ib(3,m) = i
            bend(m) = phi
         enddo
      enddo
   enddo get_bend
   nbend = m

end subroutine calc_angles

subroutine calc_torsions(n,at,xyz,bond,maxtrsn,ntrsn,trsn,it)
   use xtb_mctc_accuracy, only : wp
   implicit none
   integer,    intent(in)  :: n
   integer,    intent(in)  :: at(n)
   real(wp),   intent(in)  :: xyz(3,n)
   integer,    intent(in)  :: bond(n,n)
   integer(8), intent(in)  :: maxtrsn
   integer,    intent(out) :: ntrsn
   real(wp),   intent(out) :: trsn(maxtrsn)
   integer,    intent(out) :: it(4,maxtrsn)

   integer  :: i,j,k,l,m
   real(wp) :: phi
   real(wp),external :: valijkl

   m = 0
   get_trsn: do i = 1, n
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
               phi = valijkl(n,xyz,l,k,j,i)
               if (m.ge.maxtrsn) exit get_trsn
               m = m+1
               it(1,m) = l
               it(2,m) = k
               it(3,m) = j
               it(4,m) = i
               trsn(m) = phi
            enddo
         enddo
      enddo
   enddo get_trsn
   ntrsn = m

end subroutine calc_torsions

subroutine get_bonds(n,at,xyz,bond)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_param, only: rad => covalent_radius_2009
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: xyz(3,n)
   integer, intent(out) :: bond(n,n)

   integer  :: i,j,k
   real(wp) :: r,f,r0

   bond = 0

   do i = 1, n
      f = 1.3_wp
      k = 0
      do while(k.eq.0 .and. f.lt.1.5_wp)
         do j = 1, i-1
            r = sqrt(sum((xyz(:,j)-xyz(:,i))**2))
            r0 = rad(at(i))+rad(at(j))
            if (r.lt.f*r0) then
               k = k+1
               bond(j,i) = 1
               bond(i,j) = 1
            endif
         enddo
         f = f*1.1_wp
      enddo
      bond(i,i) = k
   enddo

end subroutine get_bonds

subroutine print_moments(iunit,n,atmass,xyz)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   implicit none
   integer, intent(in)  :: iunit
   integer, intent(in)  :: n
   real(wp),intent(in)  :: atmass(n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp) :: molmass
   real(wp) :: cma(3)
   real(wp) :: moments(3)

   molmass = molecular_mass(n,atmass)
   cma = center_of_mass(n,atmass,xyz)
   moments = moments_of_inertia(n,atmass,xyz)

   write(iunit,'(a," : ", f16.7)') "      molecular mass/u   ", molmass*autoamu
   write(iunit,'(a," : ",3f16.7)') "   center of mass at/Å   ", cma*autoaa
   write(iunit,'(a," : ",4x,3e16.7)') "  moments of inertia/u·Å²", moments*autoaa**2*autoamu
   write(iunit,'(a," : ",4x,3e16.7)') "rotational constants/cm⁻¹", 0.5_wp*autorcm/moments
   write(iunit,'(a)')

contains

pure function moments_of_inertia(n,atmass,xyz) result(moments)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_lapack, only : lapack_spev
   use xtb_mctc_convert
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: atmass(n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp) :: moments(3)
   real(wp) :: center(3),mass,t(6),work(9),tmp(3,3)
   real(wp) :: x,x2,y,y2,z,z2
   integer  :: iat,info

   center = center_of_mass(n,atmass,xyz)

   t = 0.0_wp

   do iat = 1, n
      mass = atmass(iat)*amutoau
      x = xyz(1,iat)-center(1); x2 = x**2
      y = xyz(2,iat)-center(2); y2 = y**2
      z = xyz(3,iat)-center(3); z2 = z**2
      t(1) = t(1) + mass * (y2+z2)
      t(2) = t(2) - mass * x*y
      t(3) = t(3) + mass * (x2+z2)
      t(4) = t(4) - mass * x*z
      t(5) = t(5) - mass * y*z
      t(6) = t(6) + mass * (x2+y2)
   enddo

   call lapack_spev('N','U',3,t,moments,tmp,3,work,info)
   if (info.ne.0) moments = -1.0_wp

end function moments_of_inertia

pure function center_of_mass(n,atmass,xyz) result(center)
   use xtb_mctc_accuracy, only : wp
   implicit none
   integer, intent(in)  :: n
   real(wp),intent(in)  :: atmass(n)
   real(wp),intent(in)  :: xyz(3,n)
   real(wp) :: center(3)
   integer  :: idir
   center = 0.0_wp
   do idir = 1, 3
      center(idir) = sum(atmass*xyz(idir,:))
   enddo
   center = center/sum(atmass)
end function center_of_mass

pure function molecular_mass(n,atmass) result(molmass)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   implicit none
   integer, intent(in) :: n
   real(wp),intent(in) :: atmass(n)
   real(wp) :: molmass
   molmass = sum(atmass)*amutoau
end function molecular_mass

end subroutine print_moments


subroutine check_cold_fusion(env, mol, cold_fusion)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment
   use xtb_type_molecule
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol
   logical, intent(out) :: cold_fusion
   integer :: iat, jat
   character(len=10) :: a10
   character(len=20) :: a20
   cold_fusion = .false.
   do iat = 1, len(mol)
      do jat = 1, iat-1
         if (mol%dist(jat, iat) < 1.0e-9_wp) then
            cold_fusion = .true.
            write(a20, '(a,i0,"-",a,i0)') &
               &  trim(mol%sym(jat)), jat, trim(mol%sym(iat)), iat
            write(a10, '(es10.3)') mol%dist(jat, iat)
            call env%error("Found *very* short distance of "//a10//" for "//&
               &           trim(a20))
         endif
      enddo
   enddo
end subroutine check_cold_fusion
