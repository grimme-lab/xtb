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

module xtb_sphereparam
   use xtb_mctc_accuracy, only : wp
   implicit none
   private :: wp
   public

   integer,parameter :: p_type_polynomial = 1
   integer,parameter :: p_type_logfermi   = 2
   integer  :: spherepot_type = p_type_polynomial

   integer,parameter :: p_cent_zero = 1
   integer,parameter :: p_cent_com  = 2
   integer,parameter :: p_cent_cog  = 3
   integer  :: spherecent = p_cent_zero

   integer,parameter :: p_auto_dist = 1
   integer,parameter :: p_auto_dens = 2
   integer,parameter :: p_auto_force  = 3
   integer  :: sphereauto = p_auto_dist

   integer  :: maxwalls = 0

!  old, only needed for QCG mode of docking
   real(wp) :: boxr = -1.0_wp
   real(wp) :: rabc(3) = (/-1.0_wp,-1.0_wp,-1.0_wp/)
   integer  :: sphere = -1

   integer  :: sphere_alpha = 30
   real(wp) :: sphere_beta  = 6.0_wp
   real(wp) :: sphere_temp  = 300.0_wp
   real(wp) :: sphere_autoscale = 1.0_wp
   real(wp) :: sphere_shift = 3.5_wp

   type :: tb_wall
      integer  :: fragment = 0
      integer,allocatable :: list(:)
      real(wp) :: radius(3) = 0.0_wp
      real(wp) :: center(3) = 0.0_wp
   end type tb_wall

   integer :: number_walls
   type(tb_wall),allocatable :: wpot(:)

   interface set_sphere_radius
      module procedure set_sphere_radius_iso
      module procedure set_sphere_radius_aniso
   end interface

   interface get_sphere_radius
      module procedure get_sphere_radius_list
      module procedure get_sphere_radius_fragment
      module procedure get_sphere_radius_all
   end interface

!  I need to implement this first...
!   interface get_ellipsoid_radius
!      module procedure get_ellipsoid_radius_list
!      module procedure get_ellipsoid_radius_fragment
!      module procedure get_ellipsoid_radius_all
!   end interface

!  original wall potential
   interface polynomial_cavity
      module procedure polynomial_cavity_list
      module procedure polynomial_cavity_frag
      module procedure polynomial_cavity_all
   end interface polynomial_cavity

! `bias function' by M. Shiga and M. Masia, J. Chem. Phys. 139, 044120 (2013).
!  http://dx.doi.org/10.1063/1.4816629
!  E = -Σi kT·log{1/(1+exp[-α(R0-Ri)])} = Σi kT·log{1+exp[-α(R0-Ri)]}
!  NOTE: since this is defined only for isotropic spheres, therefore we cannot
!        use anisotropic ellipsoid axis, but have to change the metric of the
!        cartesian coordinates (this is done by redefining the scalarproduct)
   interface logfermi_cavity
      module procedure logfermi_cavity_list
      module procedure logfermi_cavity_frag
      module procedure logfermi_cavity_all
   end interface logfermi_cavity

contains

!! ========================================================================
!  get some space for wall potentials, maxwalls is counted in rdcontrol
!  in setparam.f90 and is supposed to be the maximal number of possible
!  walls
subroutine init_walls
   call clear_walls
   allocate ( wpot(maxwalls) )
end subroutine init_walls

!! ========================================================================
!  clean up the heap space for reallocating or termination of the program
subroutine clear_walls
   if (allocated(wpot)) deallocate(wpot)
end subroutine clear_walls

!! ========================================================================
!  argument setter for interaction with this module in setparam.f90,
!  isotropic and anisotropic spheres are overloaded in one function
subroutine set_sphere_radius_iso(radius,center,nlist,list,fragment)
   implicit none
   real(wp),intent(in) :: radius
   real(wp),intent(in),optional :: center(3)
   integer, intent(in),optional :: nlist,list(:)
   integer, intent(in),optional :: fragment
   number_walls = number_walls + 1
   if (number_walls.gt.maxwalls) & ! This should never happen
   &  call raise('E','Number of wall potentials exceeded provided array size')
   wpot(number_walls)%radius = radius
   if (present(center)) wpot(number_walls)%center = center
   if (present(list).and.present(nlist)) &
   &  allocate( wpot(number_walls)%list(nlist), source = list(1:nlist) )
   if (present(fragment)) wpot(number_walls)%fragment = fragment
end subroutine set_sphere_radius_iso
subroutine set_sphere_radius_aniso(radius,center,nlist,list,fragment)
   implicit none
   real(wp), intent(in) :: radius(3)
   real(wp),intent(in),optional :: center(3)
   integer, intent(in),optional :: nlist,list(:)
   integer, intent(in),optional :: fragment
   number_walls = number_walls + 1
   if (number_walls.gt.maxwalls) & ! If this happens, you used it wrong!
   &  call raise('E','Number of wall potentials exceeded provided array size')
   wpot(number_walls)%radius = radius
   if (present(center)) wpot(number_walls)%center = center
   if (present(list).and.present(nlist)) &
   &  allocate( wpot(number_walls)%list(nlist), source = list(1:nlist) )
   if (present(fragment)) wpot(number_walls)%fragment = fragment
end subroutine set_sphere_radius_aniso

!! ========================================================================
!  determine sphere radius such that all atoms fit in
!  requires cema trafo (not done)
!! ------------------------------------------------------------------------
!  this is the old subroutine, not used anywhere anymore, I think...
subroutine getsphererad(n,at,xyz)
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   real(wp),intent(in)  :: xyz(3,n)

   real(wp) :: dum(3),rx,ry,rz
   real(wp) :: x,y,z,f,rr,r
   integer  :: i,j

!  call axis3(0,n,at,xyz,xyz2,dum)

!  xyz = xyz2

   r=0.0_wp
   do i=1,n-1
      do j=i+1,n
         rx=xyz(1,i)-xyz(1,j)
         ry=xyz(2,i)-xyz(2,j)
         rz=xyz(3,i)-xyz(3,j)
         rr=sqrt(rx**2+ry**2+rz**2)
         if(rr.gt.r) r=rr
      enddo
   enddo

   boxr=0.5*r+3.5d0 ! assume some vdW radius
!  this is wrong, because the cavity uses a relative measure for the constrain,
!  but here it is treated as an absolute measure, using this routine will
!  artifically compress larger molecules. I should fix this, dunno how.

end subroutine getsphererad

!! --------------------------------------------------------------[SAW1809]-
!  new version
subroutine get_sphere_radius_list(nat,at,xyz,nlist,list,center,radius,do_trafo)
   use xtb_axis, only : axis3
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   integer, intent(in)  :: nlist
   integer, intent(in)  :: list(nlist)
   real(wp),intent(out) :: center(3)
   real(wp),intent(out) :: radius
   logical, intent(in)  :: do_trafo

   integer  :: i,j
   real(wp) :: max_distance,distance,dum(3)
   real(wp),allocatable :: coord(:,:)
   logical  :: trafo
   optional :: do_trafo

   center = 0.0_wp

   if (present(do_trafo)) then
      trafo = do_trafo
   else
      trafo = .false.
   endif

   if (trafo) then
      allocate( coord(3,nat), source = 0.0_wp )
      call axis3(0,nat,at,xyz,coord,dum)
   else
      allocate( coord(3,nat), source = xyz )
   endif

   max_distance = 0.0_wp ! ~equals diameter of sphere
   do i = 1, nlist
      do j = 1, i-1
         distance = sqrt(sum((coord(:,list(j))-coord(:,list(i)))**2))
         max_distance = max(max_distance,distance)
      enddo
   enddo

   select case(spherepot_type)
   case(p_type_polynomial)
      radius = (0.5_wp/0.7_wp * max_distance) * sphere_autoscale
   case(p_type_logfermi)
      radius = (0.5_wp * max_distance + sphere_shift) * sphere_autoscale
   end select

end subroutine get_sphere_radius_list

subroutine get_sphere_radius_fragment(nat,at,xyz,center,fragment,radius,do_trafo)
   use xtb_axis, only : axis3
   use xtb_splitparam
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   integer, intent(in)  :: fragment
   real(wp),intent(out) :: center(3)
   real(wp),intent(out) :: radius
   logical, intent(in)  :: do_trafo

   integer  :: i,j
   real(wp) :: max_distance,distance,dum(3)
   real(wp),allocatable :: coord(:,:)
   logical  :: trafo
   optional :: do_trafo

   center = 0.0_wp

   if (present(do_trafo)) then
      trafo = do_trafo
   else
      trafo = .false.
   endif

   if (trafo) then
      allocate( coord(3,nat), source = 0.0_wp )
      call axis3(0,nat,at,xyz,coord,dum)
   else
      allocate( coord(3,nat), source = xyz )
   endif

   max_distance = 0.0_wp ! ~equals diameter of sphere
   do i = 1, nat
      do j = 1, i-1
         if (splitlist(i).ne.fragment.or.fragment.ne.splitlist(j)) cycle
         distance = sqrt(sum((coord(:,j)-coord(:,i))**2))
         max_distance = max(max_distance,distance)
      enddo
   enddo

   select case(spherepot_type)
   case(p_type_polynomial)
      radius = (0.5_wp/0.7_wp * max_distance) * sphere_autoscale
   case(p_type_logfermi)
      radius = (0.5_wp * max_distance + sphere_shift) * sphere_autoscale
   end select

end subroutine get_sphere_radius_fragment

subroutine get_sphere_radius_all(nat,at,xyz,center,radius,do_trafo)
   use xtb_axis, only : axis3
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: center(3)
   real(wp),intent(out) :: radius
   logical, intent(in)  :: do_trafo

   integer  :: i,j
   real(wp) :: max_distance,distance,dum(3)
   real(wp),allocatable :: coord(:,:)
   logical  :: trafo
   optional :: do_trafo

   center = 0.0_wp

   if (present(do_trafo)) then
      trafo = do_trafo
   else
      trafo = .false.
   endif

   if (trafo) then
      allocate( coord(3,nat), source = 0.0_wp )
      call axis3(0,nat,at,xyz,coord,dum)
   else
      allocate( coord(3,nat), source = xyz )
   endif

   max_distance = 0.0_wp ! ~equals diameter of sphere
   do i = 1, nat
      do j = 1, i-1
         distance = sqrt(sum((coord(:,j)-coord(:,i))**2))
         max_distance = max(max_distance,distance)
      enddo
   enddo

   select case(spherepot_type)
   case(p_type_polynomial)
      radius = (0.5_wp/0.7_wp * max_distance) * sphere_autoscale
   case(p_type_logfermi)
      radius = (0.5_wp * max_distance + sphere_shift) * sphere_autoscale
   end select

end subroutine get_sphere_radius_all

!! ========================================================================
subroutine logfermi_cavity_list(nat,at,xyz,nlist,list,temp,alpha,center,radius,&
   &                            efix,gfix)
   use xtb_mctc_constants, only : kB
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   integer, intent(in)  :: nlist
   integer, intent(in)  :: list(nlist)

   real(wp),intent(in)  :: temp      ! temperature
   real(wp),intent(in)  :: alpha     ! potential steepness
   real(wp),intent(in)  :: center(3) ! aufpunkt of wall potential
   real(wp),intent(in)  :: radius(3) ! radius of the sphere (anisotropic)

   integer  :: i,iat
   real(wp) :: r(3),w(3),dist,R0,expterm,fermi

   real(wp),intent(inout) :: efix
   real(wp),intent(inout) :: gfix(3,nat)

   R0 = maxval(radius)
   w  = R0/radius ! for anisotropy

   do i = 1, nlist
      iat = list(i)
      r = w*(xyz(:,iat) - center)
      dist = sqrt(sum(r**2))
      expterm = exp(alpha*(dist-R0))
      fermi = 1.0_wp/(1.0_wp+expterm)
      efix = efix + kB*temp * log( 1.0_wp+expterm )
      gfix(:,iat) = gfix(:,iat) + kB*temp * alpha*expterm*fermi * (r*w)/(dist+1.0e-14_wp)
   enddo

end subroutine logfermi_cavity_list

subroutine logfermi_cavity_frag(nat,at,xyz,fragment,temp,alpha,center,radius,&
   &                            efix,gfix)
   use xtb_mctc_constants, only : kB
   use xtb_splitparam
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   integer, intent(in)  :: fragment

   real(wp),intent(in)  :: temp      ! temperature
   real(wp),intent(in)  :: alpha     ! potential steepness
   real(wp),intent(in)  :: center(3) ! aufpunkt of wall potential
   real(wp),intent(in)  :: radius(3) ! radius of the sphere (anisotropic)

   integer  :: i
   real(wp) :: r(3),w(3),dist,R0,expterm,fermi

   real(wp),intent(inout) :: efix
   real(wp),intent(inout) :: gfix(3,nat)

   R0 = maxval(radius)
   w  = R0/radius ! for anisotropy

   do i = 1, nat
      if (splitlist(i).ne.fragment) cycle
      r = w*(xyz(:,i) - center)
      dist = sqrt(sum(r**2))
      expterm = exp(alpha*(dist-R0))
      fermi = 1.0_wp/(1.0_wp+expterm)
      efix = efix + kB*temp * log( 1.0_wp+expterm )
      gfix(:,i) = gfix(:,i) + kB*temp * alpha*expterm*fermi * (r*w)/(dist+1.0e-14_wp)
   enddo

end subroutine logfermi_cavity_frag

subroutine logfermi_cavity_all(nat,at,xyz,temp,alpha,center,radius,&
   &                           efix,gfix)
   use xtb_mctc_constants, only : kB
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)

   real(wp),intent(in)  :: temp      ! temperature
   real(wp),intent(in)  :: alpha     ! potential steepness
   real(wp),intent(in)  :: center(3) ! aufpunkt of wall potential
   real(wp),intent(in)  :: radius(3) ! radius of the sphere (anisotropic)

   integer  :: i
   real(wp) :: r(3),w(3),dist,R0,expterm,fermi

   real(wp),intent(inout) :: efix
   real(wp),intent(inout) :: gfix(3,nat)

   R0 = maxval(radius)
   w  = R0/radius ! for anisotropy

   do i = 1, nat
      r = w*(xyz(:,i) - center)
      dist = sqrt(sum(r**2))
      expterm = exp(alpha*(dist-R0))
      fermi = 1.0_wp/(1.0_wp+expterm)
      efix = efix + kB*temp * log( 1.0_wp+expterm )
      gfix(:,i) = gfix(:,i) + kB*temp * alpha*expterm*fermi * (r*w)/(dist+1.0e-14_wp)
   enddo

end subroutine logfermi_cavity_all

!! ========================================================================
subroutine polynomial_cavity_list(nat,at,xyz,nlist,list,alpha,center,radius,&
   &                              efix,gfix)
   use xtb_mctc_constants, only : kB
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   integer, intent(in)  :: nlist
   integer, intent(in)  :: list(nlist)

   integer, intent(in)  :: alpha     ! potential steepness
   real(wp),intent(in)  :: center(3) ! aufpunkt of wall potential
   real(wp),intent(in)  :: radius(3) ! radius of the sphere (anisotropic)

   integer  :: i,iat
   real(wp) :: r(3),w(3),dist,R0,polyterm

   real(wp),intent(inout) :: efix
   real(wp),intent(inout) :: gfix(3,nat)

   R0 = maxval(radius)
   w  = R0/radius ! for anisotropy

   do i = 1, nlist
      iat = list(i)
      r = w*(xyz(:,iat) - center)
      dist = sqrt(sum(r**2))
      polyterm = (dist/R0)**alpha
      efix = efix + polyterm
      gfix(:,iat) = gfix(:,iat) + alpha*polyterm * (r*w)/(dist**2+1.0e-14_wp)
   enddo

end subroutine polynomial_cavity_list

subroutine polynomial_cavity_frag(nat,at,xyz,fragment,alpha,center,radius,&
   &                              efix,gfix)
   use xtb_mctc_constants, only : kB
   use xtb_splitparam
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   integer, intent(in)  :: fragment

   integer, intent(in)  :: alpha     ! potential steepness
   real(wp),intent(in)  :: center(3) ! aufpunkt of wall potential
   real(wp),intent(in)  :: radius(3) ! radius of the sphere (anisotropic)

   integer  :: i
   real(wp) :: r(3),w(3),dist,R0,polyterm

   real(wp),intent(inout) :: efix
   real(wp),intent(inout) :: gfix(3,nat)

   R0 = maxval(radius)
   w  = R0/radius ! for anisotropy

   do i = 1, nat
      if (splitlist(i).ne.fragment) cycle
      r = w*(xyz(:,i) - center)
      dist = sqrt(sum(r**2))
      polyterm = (dist/R0)**alpha
      efix = efix + polyterm
      gfix(:,i) = gfix(:,i) + alpha*polyterm * (r*w)/(dist**2+1.0e-14_wp)
   enddo

end subroutine polynomial_cavity_frag

subroutine polynomial_cavity_all(nat,at,xyz,alpha,center,radius,&
   &                             efix,gfix)
   use xtb_mctc_constants, only : kB
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)

   integer, intent(in)  :: alpha     ! potential steepness
   real(wp),intent(in)  :: center(3) ! aufpunkt of wall potential
   real(wp),intent(in)  :: radius(3) ! radius of the sphere (anisotropic)

   integer  :: i
   real(wp) :: r(3),w(3),dist,R0,polyterm

   real(wp),intent(inout) :: efix
   real(wp),intent(inout) :: gfix(3,nat)

   R0 = maxval(radius)
   w  = R0/radius ! for anisotropy

   do i = 1, nat
      r = w*(xyz(:,i) - center)
      dist = sqrt(sum(r**2))
      polyterm = (dist/R0)**alpha
      efix = efix + polyterm
      gfix(:,i) = gfix(:,i) + alpha*polyterm * (r*w)/(dist**2+1.0e-14_wp)
   enddo

end subroutine polynomial_cavity_all

!! ========================================================================
!  energy and gradient of the wall potentials
subroutine cavity_egrad(nat,at,xyz,efix,gfix)
   use xtb_mctc_constants, only : kB
   implicit none
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)

   real(wp),intent(inout) :: efix
   real(wp),intent(inout) :: gfix(3,nat)

   integer :: i,nlist

   if (.not.allocated(wpot)) return

   do i = 1, number_walls
      select case(spherepot_type)
      case default ! make sure that this never happens, REALLY!
         call raise('E','Internal error in sphereparam.f90, please report this.')
      case(p_type_polynomial)
         if (wpot(i)%fragment.gt.0) then
            call polynomial_cavity(nat,at,xyz,wpot(i)%fragment,sphere_alpha, &
                 wpot(i)%center,wpot(i)%radius,efix,gfix)
         else if (allocated(wpot(i)%list)) then
            nlist = size(wpot(i)%list,1)
            call polynomial_cavity(nat,at,xyz,nlist,wpot(i)%list,sphere_alpha, &
                 wpot(i)%center,wpot(i)%radius,efix,gfix)
         else
            call polynomial_cavity(nat,at,xyz,sphere_alpha, &
                 wpot(i)%center,wpot(i)%radius,efix,gfix)
         endif
      case(p_type_logfermi)
         if (wpot(i)%fragment.gt.0) then
            call logfermi_cavity(nat,at,xyz,wpot(i)%fragment,sphere_temp, &
                 sphere_beta,wpot(i)%center,wpot(i)%radius,efix,gfix)
         else if (allocated(wpot(i)%list)) then
            nlist = size(wpot(i)%list,1)
            call logfermi_cavity(nat,at,xyz,nlist,wpot(i)%list,sphere_temp, &
                 sphere_beta,wpot(i)%center,wpot(i)%radius,efix,gfix)
         else
            call logfermi_cavity(nat,at,xyz,sphere_temp,sphere_beta, &
                 wpot(i)%center,wpot(i)%radius,efix,gfix)
         endif
      end select
   end do

end subroutine cavity_egrad

!! ========================================================================
!  old subroutines from spherepot.f
!! ========================================================================

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! energy and gradient of wall
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine cavityeg(n,xyz,e,g)
   implicit none
   real(wp) :: xyz(3,n)
   real(wp) :: g  (3,n)
   real(wp) :: e
   integer n,i,a
   real(wp) :: rij(3),dr,r,term,dum,r0,r2,rx,ry,rz,term1,term3,a2,am1

   if(sphere.lt.1) return

   if(sphere.eq.2) then
      a  =15
      am1=14
      a2=2.0d0*a
      do i=1,n
         rx=(xyz(1,i)/rabc(1))**2
         ry=(xyz(2,i)/rabc(2))**2
         rz=(xyz(3,i)/rabc(3))**2
         term1=rx+ry+rz
         term =term1**a
         term3=a2*term1**am1
         e=e+term
         g(1,i)=g(1,i)+xyz(1,i)*term3/rabc(1)**2
         g(2,i)=g(2,i)+xyz(2,i)*term3/rabc(2)**2
         g(3,i)=g(3,i)+xyz(3,i)*term3/rabc(3)**2
      enddo
      return
   endif

   if(sphere.eq.1) then
      a=30
      do i=1,n
         rij=xyz(:,i)
         r2=sum(rij*rij)
         r=sqrt(r2)
         term=(r/boxr)**a
         e=e+term
         dum=term*a/r2
         g(1,i)=g(1,i)+rij(1)*dum
         g(2,i)=g(2,i)+rij(2)*dum
         g(3,i)=g(3,i)+rij(3)*dum
      enddo
   endif

end subroutine cavityeg

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! energy of wall
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine cavitye(n,xyz,e)
   implicit none
   real(wp) :: xyz(3,n)
   real(wp) :: e
   integer n,i,a
   real(wp) :: rij(3),dr,r,term,dum,r0,r2,rx,ry,rz

   if(sphere.lt.1) return

   if(sphere.eq.2) then
      a=15
      do i=1,n
         rx=(xyz(1,i)/rabc(1))**2
         ry=(xyz(2,i)/rabc(2))**2
         rz=(xyz(3,i)/rabc(3))**2
         term=(rx+ry+rz)**a
         e=e+term
      enddo
      return
   endif

   if(sphere.eq.1) then
      a=30
      do i=1,n
         rij=xyz(:,i)
         r2=sum(rij*rij)
         r=sqrt(r2)
         term=(r/boxr)**a
         e=e+term
      enddo
      return
   endif

end subroutine cavitye
end module xtb_sphereparam
