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

!! ========================================================================
!  call as:
!  interface
!     pure subroutine ncoord_name(nat,at,xyz,cn,thr)
!     use iso_fortran_env, only : wp => real64
!     implicit none
!     integer, intent(in)  :: nat
!     integer, intent(in)  :: at(nat)
!     real(wp),intent(in)  :: xyz(3,nat)
!     real(wp),intent(out) :: cn(nat)
!     real(wp),intent(in),optional :: thr
!  end interface
!  interface
!     pure subroutine dncoord_name(nat,at,xyz,cn,dcn,thr)
!     use iso_fortran_env, only : wp => real64
!     implicit none
!     integer, intent(in)  :: nat
!     integer, intent(in)  :: at(nat)
!     real(wp),intent(in)  :: xyz(3,nat)
!     real(wp),intent(out) :: cn(nat)
!     real(wp),intent(out) :: dcn(3,nat,nat)
!     real(wp),intent(in),optional :: thr
!  end interface
!! ========================================================================
module ncoord
   use iso_fortran_env, only : wp => real64
   implicit none

   real(wp),private,parameter :: cnthr = 1600.0_wp

   real(wp),parameter :: k1 = 16.0_wp
   real(wp),parameter :: k2 = 4.0_wp/3.0_wp

   real(wp),parameter :: ka=10.0_wp
   real(wp),parameter :: kb=20.0_wp
   real(wp),parameter :: r_shift=2.0_wp

   real(wp),parameter :: k4=4.10451_wp
   real(wp),parameter :: k5=19.08857_wp
   real(wp),parameter :: k6=2*11.28174_wp**2

   real(wp),parameter :: kn=7.50_wp
   real(wp),parameter :: kr=0.25_wp
   real(wp),parameter :: ke=0.05_wp

   integer,private,parameter :: max_elem = 118

!  covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
!  188-197), values for metals decreased by 10 %
   real(wp),private,parameter :: rad(max_elem) = (/  &
   & 0.32,0.46, & ! H,He
   & 1.20,0.94,0.77,0.75,0.71,0.63,0.64,0.67, & ! Li-Ne
   & 1.40,1.25,1.13,1.04,1.10,1.02,0.99,0.96, & ! Na-Ar
   & 1.76,1.54, & ! K,Ca
   &           1.33,1.22,1.21,1.10,1.07,1.04,1.00,0.99,1.01,1.09, & ! Sc-Zn
   &           1.12,1.09,1.15,1.10,1.14,1.17, & ! Ga-Kr
   & 1.89,1.67, & ! Rb,Sr
   &           1.47,1.39,1.32,1.24,1.15,1.13,1.13,1.08,1.15,1.23, & ! Y-Cd
   &           1.28,1.26,1.26,1.23,1.32,1.31, & ! In-Xe
   & 2.09,1.76, & ! Cs,Ba
   &      1.62,1.47,1.58,1.57,1.56,1.55,1.51, & ! La-Eu
   &      1.52,1.51,1.50,1.49,1.49,1.48,1.53, & ! Gd-Yb
   &           1.46,1.37,1.31,1.23,1.18,1.16,1.11,1.12,1.13,1.32, & ! Lu-Hg
   &           1.30,1.30,1.36,1.31,1.38,1.42, & ! Tl-Rn
   & 2.01,1.81, & ! Fr,Ra
   &      1.67,1.58,1.52,1.53,1.54,1.55,1.49, & ! Ac-Am
   &      1.49,1.51,1.51,1.48,1.50,1.56,1.58, & ! Cm-No
   &           1.45,1.41,1.34,1.29,1.27,1.21,1.16,1.15,1.09,1.22, & ! Lr-Cn
   &           1.22,1.29,1.46,1.58,1.48,1.41 /) ! Nh-Og
   real(wp),parameter :: rcov(max_elem) = 4.0_wp/3.0_wp*rad/0.52917726_wp

!  pauling EN's
   real(wp),parameter :: en(max_elem) = (/ &
   & 2.20,3.00, & ! H,He
   & 0.98,1.57,2.04,2.55,3.04,3.44,3.98,4.50, & ! Li-Ne
   & 0.93,1.31,1.61,1.90,2.19,2.58,3.16,3.50, & ! Na-Ar
   & 0.82,1.00, & ! K,Ca
   &           1.36,1.54,1.63,1.66,1.55,1.83,1.88,1.91,1.90,1.65, & ! Sc-Zn
   &           1.81,2.01,2.18,2.55,2.96,3.00, & ! Ga-Kr
   & 0.82,0.95, & ! Rb,Sr
   &           1.22,1.33,1.60,2.16,1.90,2.20,2.28,2.20,1.93,1.69, & ! Y-Cd
   &           1.78,1.96,2.05,2.10,2.66,2.60, & ! In-Xe
   & 0.79,0.89, & ! Cs,Ba
   &      1.10,1.12,1.13,1.14,1.15,1.17,1.18, & ! La-Eu
   &      1.20,1.21,1.22,1.23,1.24,1.25,1.26, & ! Gd-Yb
   &           1.27,1.30,1.50,2.36,1.90,2.20,2.20,2.28,2.54,2.00, & ! Lu-Hg
   &           1.62,2.33,2.02,2.00,2.20,2.20, & ! Tl-Rn
   ! only dummies below
   & 1.50,1.50, & ! Fr,Ra
   &      1.50,1.50,1.50,1.50,1.50,1.50,1.50, & ! Ac-Am
   &      1.50,1.50,1.50,1.50,1.50,1.50,1.50, & ! Cm-No
   &           1.50,1.50,1.50,1.50,1.50,1.50,1.50,1.50,1.50,1.50, & ! Rf-Cn
   &           1.50,1.50,1.50,1.50,1.50,1.50 /) ! Nh-Og


   interface get_d3_cn
      module procedure ncoord_d3_driver
      module procedure ncoord_d3
      module procedure dncoord_d3
   end interface
   interface get_erf_cn
      module procedure ncoord_erf_driver
      module procedure ncoord_erf
      module procedure dncoord_erf
   end interface
   interface get_d4_cn
      module procedure ncoord_d4_driver
      module procedure ncoord_d4
      module procedure dncoord_d4
      module procedure pbc_ncoord_d4
      module procedure pbc_dncoord_d4
   end interface
   interface get_gfn_cn
      module procedure ncoord_gfn_driver
      module procedure ncoord_gfn
      module procedure dncoord_gfn
   end interface

contains

pure subroutine ncoord_d3_driver(mol,cn,dcndr,thr)
   use tbdef_molecule
   implicit none
   !> molecular structure information
   type(tb_molecule), intent(in) :: mol
   !> coordination number
   real(wp), intent(out) :: cn(:)
   !> derivative of coordination number w.r.t. nuclear positions
   real(wp), intent(out), optional :: dcndr(:,:,:)
   !> quadratic distance threshold for sum termination
   real(wp), intent(in), optional :: thr
   real(wp) :: cn_thr

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

!   if (mol%npbc > 0) then
!      if (present(dcndr).and.present(dcndL)) then
!      else
!      endif
!   else
      if (present(dcndr)) then
         call dncoord_d3(mol%n,mol%at,mol%xyz,cn,dcndr,cn_thr)
      else
         call ncoord_d3(mol%n,mol%at,mol%xyz,cn,cn_thr)
      endif
!   endif

end subroutine ncoord_d3_driver

pure subroutine ncoord_gfn_driver(mol,cn,dcndr,thr)
   use tbdef_molecule
   implicit none
   !> molecular structure information
   type(tb_molecule), intent(in) :: mol
   !> coordination number
   real(wp), intent(out) :: cn(:)
   !> derivative of coordination number w.r.t. nuclear positions
   real(wp), intent(out), optional :: dcndr(:,:,:)
   !> quadratic distance threshold for sum termination
   real(wp), intent(in), optional :: thr
   real(wp) :: cn_thr

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

!   if (mol%npbc > 0) then
!      if (present(dcndr).and.present(dcndL)) then
!      else
!      endif
!   else
      if (present(dcndr)) then
         call dncoord_gfn(mol%n,mol%at,mol%xyz,cn,dcndr,cn_thr)
      else
         call ncoord_gfn(mol%n,mol%at,mol%xyz,cn,cn_thr)
      endif
!   endif

end subroutine ncoord_gfn_driver

pure subroutine ncoord_erf_driver(mol,cn,dcndr,dcndL,thr)
   use tbdef_molecule
   implicit none
   !> molecular structure information
   type(tb_molecule), intent(in) :: mol
   !> coordination number
   real(wp), intent(out) :: cn(:)
   !> derivative of coordination number w.r.t. nuclear positions
   real(wp), intent(out), optional :: dcndr(:,:,:)
   !> derivative of coordination number w.r.t. lattice vectors
   real(wp), intent(out), optional :: dcndL(:,:,:)
   !> quadratic distance threshold for sum termination
   real(wp), intent(in), optional :: thr
   real(wp) :: cn_thr

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   if (mol%npbc > 0) then
      if (present(dcndr).and.present(dcndL)) then
         call pbc_derfcoord(mol%n,mol%at,mol%xyz,mol%lattice,cn,dcndr,dcndL,cn_thr)
      else
         call pbc_erfcoord(mol%n,mol%at,mol%xyz,mol%lattice,cn,cn_thr)
      endif
   else
      if (present(dcndr)) then
         call dncoord_erf(mol%n,mol%at,mol%xyz,cn,dcndr,cn_thr)
      else
         call ncoord_erf(mol%n,mol%at,mol%xyz,cn,cn_thr)
      endif
   endif

end subroutine ncoord_erf_driver

pure subroutine ncoord_d4_driver(mol,cn,dcndr,dcndL,thr)
   use tbdef_molecule
   implicit none
   !> molecular structure information
   type(tb_molecule), intent(in) :: mol
   !> coordination number
   real(wp), intent(out) :: cn(:)
   !> derivative of coordination number w.r.t. nuclear positions
   real(wp), intent(out), optional :: dcndr(:,:,:)
   !> derivative of coordination number w.r.t. lattice vectors
   real(wp), intent(out), optional :: dcndL(:,:,:)
   !> quadratic distance threshold for sum termination
   real(wp), intent(in), optional :: thr
   real(wp) :: cn_thr

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   if (mol%npbc > 0) then
      if (present(dcndr).and.present(dcndL)) then
         call pbc_dncoord_d4(mol%n,mol%at,mol%xyz,mol%lattice,cn,dcndr,dcndL,cn_thr)
      else
         call pbc_ncoord_d4(mol%n,mol%at,mol%xyz,mol%lattice,cn,cn_thr)
      endif
   else
      if (present(dcndr)) then
         call dncoord_d4(mol%n,mol%at,mol%xyz,cn,dcndr,cn_thr)
      else
         call ncoord_d4(mol%n,mol%at,mol%xyz,cn,cn_thr)
      endif
   endif

end subroutine ncoord_d4_driver

! ========================================================================
!> original D3 type coordination number from 2010
pure subroutine ncoord_d3(nat,at,xyz,cn,thr)

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j
   real(wp) :: rij(3), r, rco, den, tmp, r2

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn = 0.0_wp

   do i = 1, nat
      do j = 1, i-1
         rij = xyz(:,j) - xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
!        covalent distance in bohr
         rco=rcov(at(j)) + rcov(at(i))
!        counting function exponential has a better long-range
!        behavior than MHGs inverse damping
         tmp=exp_count(k1,r,rco)
         cn(i)=cn(i)+tmp
         cn(j)=cn(j)+tmp
      enddo
   enddo

end subroutine ncoord_d3

! ========================================================================
!  original D3 type coordination number from 2010
!  INPUT
!  nat  :: number of atoms
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  OUTPUT
!  cn   :: coordination number
!  dcndr  :: derivative of coordination number w.r.t. atom position
!  PARAMETERS: k1,k2
!  NOTE: k2 is already included in rcov
pure subroutine dncoord_d3(nat,at,xyz,cn,dcndr,thr)

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(out) :: dcndr(3,nat,nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i, j
   real(wp) :: r, r2, rij(3)
   real(wp) :: rcovij
   real(wp) :: expterm
   real(wp) :: dtmp, tmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn  = 0._wp
   dcndr = 0._wp

   do i = 1, nat
      do j = 1, i-1
         rij = xyz(:,j) - xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r = sqrt(r2)
         rcovij=(rcov(at(i))+rcov(at(j)))
         tmp = exp_count(k1,r,rcovij)
         dtmp = dexp_count(k1,r,rcovij)
         cn(i) = cn(i) + tmp
         cn(j) = cn(j) + tmp
         dcndr(:,i,i)=-dtmp*rij/r + dcndr(:,i,i) ! FIXME
         dcndr(:,j,j)= dtmp*rij/r + dcndr(:,j,j) ! FIXME
         dcndr(:,i,j)= dtmp*rij/r
         dcndr(:,j,i)=-dtmp*rij/r
      enddo
   enddo

end subroutine dncoord_d3

! gradients for pbc coordination number with exp function
pure subroutine pbc_dncoord_d3(nat,at,xyz,lat,cn,dcndr,dcndL,thr)
   use mctc_constants
   use pbc_tools
   use pbc, only : get_realspace_cutoff

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: lat(3,3)
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(out) :: dcndr(3,nat,nat)
   real(wp),intent(out) :: dcndL(3,3,nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j,tx,ty,tz
   real(wp) :: rij(3), r, rco, den, rr, r2, t(3)
   integer  :: rep_cn(3)
   real(wp) :: tmp,dtmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   call get_realspace_cutoff(lat,cn_thr,rep_cn)

   cn = 0.0_wp
   dcndr = 0.0_wp
   dcndL = 0.0_wp

   do i = 1, nat
      do j = 1, i-1 ! loop over all atoms for PBC case
         do concurrent(tx = -rep_cn(1):rep_cn(1), &
               &       ty = -rep_cn(2):rep_cn(2), &
               &       tz = -rep_cn(3):rep_cn(3))
            t = tx*lat(:,1) + ty*lat(:,2) + tz*lat(:,3)
            rij = xyz(:,j) - xyz(:,i) + t
            r2  = sum(rij**2)
            if (r2.gt.cn_thr) cycle
            r=sqrt(r2)
!           covalent distance in bohr
            rco=rcov(at(j)) + rcov(at(i))
            tmp=exp_count(k1,r,rco)
            dtmp = dexp_count(k1,r,rco)
            cn(i)=cn(i) + tmp
            cn(j)=cn(j) + tmp
            dcndr(:,i,i)= dtmp*rij/r + dcndr(:,i,i)
            dcndr(:,j,j)=-dtmp*rij/r + dcndr(:,j,j)
            dcndr(:,i,j)= dtmp*rij/r + dcndr(:,i,j)
            dcndr(:,j,i)=-dtmp*rij/r + dcndr(:,j,i)
            dcndL(:,:,j)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,j)
            dcndL(:,:,i)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,i)
         enddo
      enddo
      do concurrent(tx = -rep_cn(1):rep_cn(1), &
            &       ty = -rep_cn(2):rep_cn(2), &
            &       tz = -rep_cn(3):rep_cn(3), &
            &       tx.ne.0 .or. ty.ne.0 .or. tz.ne.0)
         t = tx*lat(:,1) + ty*lat(:,2) + tz*lat(:,3)
         rij = t
         r2  = sum(rij**2)
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
!        covalent distance in bohr
         rco=2*rcov(at(i))
         tmp=exp_count(k1,r,rco)
         dtmp = dexp_count(k1,r,rco)
         cn(i)=cn(i) + tmp
         dcndL(:,:,i)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,i)
      enddo
   enddo

end subroutine pbc_dncoord_d3



! ========================================================================
!  modified D3 type coordination number from 2018
!  INPUT
!  nat  :: number of atoms
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  OUTPUT
!  cn   :: coordination number
!  PARAMETERS: kn,k2
!  NOTE: k2 is already included in rcov
pure subroutine ncoord_erf(nat,at,xyz,cn,thr)

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j
   real(wp) :: rij(3), r, rco, den, rr, r2, tmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn = 0.0_wp

   do i = 1, nat
      do j = 1, i-1
         rij = xyz(:,j) - xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
!        covalent distance in bohr
         rco=rcov(at(j)) + rcov(at(i))
!        error function has an even better long range behavior
         tmp = erf_count(kn,r,rco)
         cn(i)=cn(i)+tmp
         cn(j)=cn(j)+tmp
      enddo
   enddo

end subroutine ncoord_erf

! ========================================================================
!  modified D3 type coordination number from 2018
!  INPUT
!  nat  :: number of atoms
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  OUTPUT
!  cn   :: coordination number
!  dcndr  :: derivative of coordination number w.r.t. atom position
!  PARAMETERS: kn,k2
!  NOTE: k2 is already included in rcov
pure subroutine dncoord_erf(nat,at,xyz,cn,dcndr,thr)

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(out) :: dcndr(3,nat,nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j
   real(wp) :: r,r2,rij(3)
   real(wp) :: rcovij
   real(wp) :: dtmp,tmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn  = 0._wp
   dcndr = 0._wp

   do i = 1, nat
      do j = 1, i-1
         rij = xyz(:,j) - xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r = sqrt(r2)
         rcovij=(rcov(at(i))+rcov(at(j)))
         tmp = erf_count(kn,r,rcovij)
         dtmp = derf_count(kn,r,rcovij)
         cn(i) = cn(i) + tmp
         cn(j) = cn(j) + tmp
         dcndr(:,i,i)= dtmp*rij/r + dcndr(:,i,i)
         dcndr(:,j,j)=-dtmp*rij/r + dcndr(:,j,j)
         dcndr(:,i,j)= dtmp*rij/r
         dcndr(:,j,i)=-dtmp*rij/r
      enddo
   enddo

end subroutine dncoord_erf

! ========================================================================
!  covalent coordination number of the DFT-D4 method
!  INPUT
!  nat  :: number of atoms
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  OUTPUT
!  cn   :: coordination number
!  PARAMETERS: k1,k2,k4,k5,k6
!  NOTE: k2 is already included in rcov
pure subroutine ncoord_d4(nat,at,xyz,cn,thr)

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j
   real(wp) :: rij(3), r, rco, den, rr, r2, xn, tmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn = 0._wp

   do i=1,nat
      do j=1,i-1
         rij = xyz(:,j) - xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
!        covalent distance in bohr
         rco=rcov(at(j)) + rcov(at(i))
         rr=rco/r
         den=k4*exp(-(abs((en(at(i))-en(at(j))))+ k5)**2/k6 )
!        counting function exponential has a better long-range
!        behavior than MHGs inverse damping
         !tmp = den/(1.d0+exp(-k1*(rr-1.0d0)))
!        error function has an even better long range behavior
         tmp = den * erf_count(kn,r,rco)
         cn(i)=cn(i)+tmp
         cn(j)=cn(j)+tmp
      enddo
   enddo

end subroutine ncoord_d4

! ========================================================================
!  derivative of the covalent coordination number of the DFT-D4 method
!  NOTE: the derivative is inlined in the dispgrad in dftd4 by hand
!  INPUT
!  nat  :: number of atoms
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  OUTPUT
!  cn   :: coordination number
!  dcndr  :: derivative of coordination number w.r.t. atom position
!  PARAMETERS: k1,k2,k4,k5,k6
!  NOTE: k2 is already included in rcov
pure subroutine dncoord_d4(nat,at,xyz,cn,dcndr,thr)
   use mctc_constants

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(out) :: dcndr(3,nat,nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i, j, ia, ja
   real(wp) :: r, r2, rij(3)
   real(wp) :: rcovij,den
   real(wp) :: dtmp, tmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn  = 0._wp
   dcndr = 0._wp

   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         rij = xyz(:,j) - xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r = sqrt(r2)
         rcovij=(rcov(ia)+rcov(ja))
         den=k4*exp(-(abs((en(ia)-en(ja)))+ k5)**2/k6 )
         tmp = den * erf_count(kn,r,rcovij)
         dtmp = den * derf_count(kn,r,rcovij)
         cn(i) = cn(i) + tmp
         cn(j) = cn(j) + tmp
         dcndr(:,i,i)=-dtmp*rij/r + dcndr(:,i,i)
         dcndr(:,j,j)= dtmp*rij/r + dcndr(:,j,j)
         dcndr(:,i,j)=-dtmp*rij/r
         dcndr(:,j,i)= dtmp*rij/r
      enddo
   enddo

end subroutine dncoord_d4

! same for the pbc case
pure subroutine pbc_ncoord_d4(nat,at,xyz,lat,cn,thr)
   use mctc_constants
   use pbc_tools
   use pbc, only : get_realspace_cutoff
   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   integer              :: rep_cn(3)
   integer              :: tx,ty,tz
   real(wp)             :: t(3)
   real(wp),intent(in)  :: lat(3,3)
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j
   real(wp) :: rij(3), r, rco, den, r2, xn, tmp, dtmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn  = 0.0_wp

   call get_realspace_cutoff(lat,cn_thr,rep_cn)

   do i=1,nat
      do j=1,i-1
         do concurrent(tx = -rep_cn(1):rep_cn(1), &
               &       ty = -rep_cn(2):rep_cn(2), &
               &       tz = -rep_cn(3):rep_cn(3))
            t = [tx,ty,tz]
            rij = xyz(:,j) - xyz(:,i) + matmul(lat,t)
            r2  = sum(rij**2)
            if (r2.gt.cn_thr) cycle
            r=sqrt(r2)
            ! covalent distance in bohr
            rco=rcov(at(j)) + rcov(at(i))
            den=k4*exp(-(abs((en(at(i))-en(at(j))))+ k5)**2/k6 )
            tmp = den * erf_count(kn,r,rco)
            cn(i)=cn(i)+tmp
            cn(j)=cn(j)+tmp
         enddo
      enddo
      do concurrent(tx = -rep_cn(1):rep_cn(1), &
            &       ty = -rep_cn(2):rep_cn(2), &
            &       tz = -rep_cn(3):rep_cn(3), &
            &       tx.ne.0 .or. ty.ne.0 .or. tz.ne.0)
         t = [tx,ty,tz]
         rij = matmul(lat,t)
         r2  = sum(rij**2)
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
         ! covalent distance in bohr
         rco=2*rcov(at(i))
         den=k4*exp(-k5**2/k6 )
         tmp = den * erf_count(kn,r,rco)
         cn(i)=cn(i)+tmp
      enddo
   enddo

end subroutine pbc_ncoord_d4

! same for the pbc case
pure subroutine pbc_dncoord_d4(nat,at,xyz,lat,cn,dcndr,dcndL,thr)
   use mctc_constants
   use pbc_tools
   use pbc, only : get_realspace_cutoff
   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   integer              :: rep_cn(3)
   integer              :: tx,ty,tz
   real(wp)             :: t(3)
   real(wp),intent(in)  :: lat(3,3)
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(out) :: dcndr(3,nat,nat)
   real(wp),intent(out) :: dcndL(3,3,nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j
   real(wp) :: rij(3), r, rco, den, r2, xn, tmp, dtmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn  = 0.0_wp
   dcndr = 0.0_wp
   dcndL = 0.0_wp

   call get_realspace_cutoff(lat,cn_thr,rep_cn)

   do i=1,nat
      do j=1,i-1
         do concurrent(tx = -rep_cn(1):rep_cn(1), &
               &       ty = -rep_cn(2):rep_cn(2), &
               &       tz = -rep_cn(3):rep_cn(3))
            t = [tx,ty,tz]
            rij = xyz(:,j) - xyz(:,i) + matmul(lat,t)
            r2  = sum(rij**2)
            if (r2.gt.cn_thr) cycle
            r=sqrt(r2)
            ! covalent distance in bohr
            rco=rcov(at(j)) + rcov(at(i))
            den=k4*exp(-(abs((en(at(i))-en(at(j))))+ k5)**2/k6 )
            tmp = den * erf_count(kn,r,rco)
            dtmp = den * derf_count(kn,r,rco)
            cn(i)=cn(i)+tmp
            cn(j)=cn(j)+tmp
            dcndr(:,i,i)=-dtmp*rij/r + dcndr(:,i,i)
            dcndr(:,j,j)= dtmp*rij/r + dcndr(:,j,j)
            dcndr(:,i,j)=-dtmp*rij/r + dcndr(:,i,j)
            dcndr(:,j,i)= dtmp*rij/r + dcndr(:,j,i)
            dcndL(:,:,j)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,j)
            dcndL(:,:,i)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,i)
         enddo
      enddo
      do concurrent(tx = -rep_cn(1):rep_cn(1), &
            &       ty = -rep_cn(2):rep_cn(2), &
            &       tz = -rep_cn(3):rep_cn(3))
         ! avoid self interaction
         if ((tx.eq.0).and.(ty.eq.0).and.(tz.eq.0)) cycle
         t = [tx,ty,tz]
         rij = matmul(lat,t)
         r2  = sum(rij**2)
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
         ! covalent distance in bohr
         rco=2*rcov(at(i))
         den=k4*exp(-k5**2/k6)
         tmp = den * erf_count(kn,r,rco)
         dtmp = den * derf_count(kn,r,rco)
         cn(i)=cn(i)+tmp
         dcndL(:,:,i)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,i)
      enddo
   enddo

end subroutine pbc_dncoord_d4

! ========================================================================
!  GFN2-xTB coordination number for CN dependent parts of the Hamiltonian
!  it's similar to the D3 coordination number but doubly damped, to avoid
!  large tail contributions in dense packed systems
!  INPUT
!  nat  :: number of atoms
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  OUTPUT
!  cn   :: coordination number
!  PARAMETERS: ka,kb,k2,r_shift
!  NOTE: k2 is already included in rcov
pure subroutine ncoord_gfn(nat,at,xyz,cn,thr)

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i, j, ia, ja
   real(wp) :: r, r2, rij(3)
   real(wp) :: rcovij
   real(wp) :: expterm1, expterm2
   real(wp) :: tmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn  = 0.0_wp

   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         rij = xyz(:,j) - xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r = sqrt(r2)
         rcovij=(rcov(ia)+rcov(ja))
         expterm1=exp_count(ka,r,rcovij)
         expterm2=exp_count(kb,r,rcovij+r_shift)
         tmp = expterm1*expterm2
         cn(i) = cn(i) + tmp
         cn(j) = cn(j) + tmp
      enddo
   enddo

end subroutine ncoord_gfn

! ========================================================================
!  GFN2-xTB coordination number for CN dependent parts of the Hamiltonian
!  it's similar to the D3 coordination number but doubly damped, to avoid
!  large tail contributions in dense packed systems
!  INPUT
!  nat  :: number of atoms
!  at   :: ordinal number of atoms
!  xyz  :: coordinates in Bohr
!  OUTPUT
!  cn   :: coordination number
!  dcndr  :: derivative of coordination number w.r.t. atom position
!  PARAMETERS: ka,kb,k2,r_shift
!  NOTE: k2 is already included in rcov
pure subroutine dncoord_gfn(nat,at,xyz,cn,dcndr,thr)

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(out) :: dcndr(3,nat,nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i, j, ia, ja
   real(wp) :: r, r2, rij(3)
   real(wp) :: rcovij
   real(wp) :: expterm1, expterm2
   real(wp) :: dtmp, tmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   cn  = 0._wp
   dcndr = 0._wp

   do i = 1, nat
      ia = at(i)
      do j = 1, i-1
         ja = at(j)
         rij = xyz(:,j) - xyz(:,i)
         r2  = sum( rij**2 )
         if (r2.gt.cn_thr) cycle
         r = sqrt(r2)
         rcovij=(rcov(ia)+rcov(ja))
         expterm1=exp(-ka*(rcovij/r-1._wp))
         expterm2=exp(-kb*((rcovij+r_shift)/r-1._wp))
         tmp = 1._wp/(1._wp+expterm1)/(1._wp+expterm2)
         dtmp = (-ka*rcovij*expterm1) &
         &  /(r2*((expterm1+1._wp)**2))/(1._wp+expterm2) &
         &  + (-kb*(rcovij+r_shift)*expterm2) &
         &  /(r2*((expterm2+1._wp)**2))/(1._wp+expterm1)
         cn(i) = cn(i) + tmp
         cn(j) = cn(j) + tmp
         dcndr(:,i,i)=-dtmp*rij/r + dcndr(:,i,i)
         dcndr(:,j,j)= dtmp*rij/r + dcndr(:,j,j)
         dcndr(:,i,j)= dtmp*rij/r
         dcndr(:,j,i)=-dtmp*rij/r
      enddo
   enddo

end subroutine dncoord_gfn

! pbc coordination number with error function
pure subroutine pbc_erfcoord(nat,at,xyz,lat,cn,thr)
   use pbc, only : get_realspace_cutoff

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: lat(3,3)
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j,tx,ty,tz
   real(wp) :: rij(3), r, rco, den, rr, r2, t(3)
   integer  :: rep_cn(3)
   real(wp) :: tmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   call get_realspace_cutoff(lat,cn_thr,rep_cn)

   cn = 0.0_wp

   do i = 1, nat
      do j = 1, i-1 ! loop over all atoms for PBC case
         do concurrent(tx = -rep_cn(1):rep_cn(1), &
               &       ty = -rep_cn(2):rep_cn(2), &
               &       tz = -rep_cn(3):rep_cn(3))
            ! avoid self interaction
            if ((j.eq.i).and.(tx.eq.0).and.(ty.eq.0).and.(tz.eq.0)) cycle
            t = tx*lat(:,1) + ty*lat(:,2) + tz*lat(:,3)
            rij = xyz(:,j) - xyz(:,i) + t
            r2  = sum(rij**2)
            if (r2.gt.cn_thr) cycle
            r=sqrt(r2)
!           covalent distance in bohr
            rco=rcov(at(j)) + rcov(at(i))

            tmp=erf_count(kn,r,rco)
            cn(i)=cn(i) + tmp
            cn(j)=cn(j) + tmp
         enddo
      enddo
      do concurrent(tx = -rep_cn(1):rep_cn(1), &
            &       ty = -rep_cn(2):rep_cn(2), &
            &       tz = -rep_cn(3):rep_cn(3), &
            &       tx.ne.0 .or. ty.ne.0 .or. tz.ne.0)
         t = tx*lat(:,1) + ty*lat(:,2) + tz*lat(:,3)
         rij = t
         r2  = sum(rij**2)
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
!        covalent distance in bohr
         rco=2*rcov(at(i))

         tmp=erf_count(kn,r,rco)
         cn(i)=cn(i) + tmp
      enddo
   enddo

end subroutine pbc_erfcoord

! gradients for pbc coordination number with error function
pure subroutine pbc_derfcoord(nat,at,xyz,lat,cn,dcndr,dcndL,thr)
   use mctc_constants
   use pbc_tools
   use pbc, only : get_realspace_cutoff

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: lat(3,3)
   real(wp),intent(out) :: cn(nat)
   real(wp),intent(out) :: dcndr(3,nat,nat)
   real(wp),intent(out) :: dcndL(3,3,nat)
   real(wp),intent(in),optional :: thr
   real(wp) :: cn_thr

   integer  :: i,j,tx,ty,tz
   real(wp) :: rij(3), r, rco, den, rr, r2, t(3)
   integer  :: rep_cn(3)
   real(wp) :: tmp,dtmp

   if (present(thr)) then
      cn_thr = thr
   else
      cn_thr = cnthr
   endif

   call get_realspace_cutoff(lat,cn_thr,rep_cn)

   cn = 0.0_wp
   dcndr = 0.0_wp
   dcndL = 0.0_wp

   do i = 1, nat
      do j = 1, i-1 ! loop over all atoms for PBC case
         do concurrent(tx = -rep_cn(1):rep_cn(1), &
               &       ty = -rep_cn(2):rep_cn(2), &
               &       tz = -rep_cn(3):rep_cn(3))
            t = tx*lat(:,1) + ty*lat(:,2) + tz*lat(:,3)
            rij = xyz(:,j) - xyz(:,i) + t
            r2  = sum(rij**2)
            if (r2.gt.cn_thr) cycle
            r=sqrt(r2)
!           covalent distance in bohr
            rco=rcov(at(j)) + rcov(at(i))
            tmp=erf_count(kn,r,rco)
            dtmp = derf_count(kn,r,rco)
            cn(i)=cn(i) + tmp
            cn(j)=cn(j) + tmp
            dcndr(:,i,i)= dtmp*rij/r + dcndr(:,i,i)
            dcndr(:,j,j)=-dtmp*rij/r + dcndr(:,j,j)
            dcndr(:,i,j)= dtmp*rij/r + dcndr(:,i,j)
            dcndr(:,j,i)=-dtmp*rij/r + dcndr(:,j,i)
            dcndL(:,:,j)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,j)
            dcndL(:,:,i)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,i)
         enddo
      enddo
      do concurrent(tx = -rep_cn(1):rep_cn(1), &
            &       ty = -rep_cn(2):rep_cn(2), &
            &       tz = -rep_cn(3):rep_cn(3), &
            &       tx.ne.0 .or. ty.ne.0 .or. tz.ne.0)
         t = tx*lat(:,1) + ty*lat(:,2) + tz*lat(:,3)
         rij = t
         r2  = sum(rij**2)
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
!        covalent distance in bohr
         rco=2*rcov(at(i))
         tmp=erf_count(kn,r,rco)
         dtmp = derf_count(kn,r,rco)
         cn(i)=cn(i) + tmp
         dcndL(:,:,i)= dtmp*outer_prod_3x3(rij,rij)/r + dcndL(:,:,i)
      enddo
   enddo

end subroutine pbc_derfcoord

!> cutoff function for large coordination numbers
subroutine dncoord_logcn(n,cn,dcndr,dcndL,cn_max)
   implicit none
   !> number of atoms
   integer, intent(in) :: n
   !> on input coordination number, on output modified CN
   real(wp), intent(inout) :: cn(n)
   !> on input derivative of CN w.r.t. cartesian coordinates,
   !  on output derivative of modified CN
   real(wp), intent(inout), optional :: dcndr(3,n,n)
   !> on input derivative of CN w.r.t. strain deformation,
   !  on output derivative of modified CN
   real(wp), intent(inout), optional :: dcndL(3,3,n)
   !> maximum CN (not strictly obeyed)
   real(wp), intent(in), optional :: cn_max
   !  local
   real(wp) :: cnmax
   real(wp) :: dcnpdcn
   integer  :: i

   if (present(cn_max)) then
      cnmax = max(cn_max,0.0_wp)
   else
      cnmax = 4.5_wp
   endif

   if (present(dcndL)) then
      do i = 1, n
         dcnpdcn = dlog_cn_cut(cn(i),cnmax)
         dcndL(:,:,i) = dcnpdcn*dcndL(:,:,i)
      enddo
   endif

   if (present(dcndr)) then
      do i = 1, n
         dcnpdcn = dlog_cn_cut(cn(i),cnmax)
         dcndr(:,:,i) = dcnpdcn*dcndr(:,:,i)
      enddo
   endif

   do i = 1, n
      cn(i) = log_cn_cut(cn(i),cnmax)
   enddo

end subroutine dncoord_logcn

pure elemental function log_cn_cut(cn,cnmax) result(cnp)
   real(wp), intent(in) :: cn
   real(wp), intent(in) :: cnmax
   real(wp) :: cnp
   cnp = log(1.0_wp + exp(cnmax)) - log(1.0_wp + exp(cnmax - cn))
end function log_cn_cut

pure elemental function dlog_cn_cut(cn,cnmax) result(dcnpdcn)
   real(wp), intent(in) :: cn
   real(wp), intent(in) :: cnmax
   real(wp) :: dcnpdcn
   dcnpdcn = exp(cnmax)/(exp(cnmax) + exp(cn))
end function dlog_cn_cut

pure elemental function erf_count(k,r,r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   count = 0.5_wp * (1.0_wp + erf(-k*(r-r0)/r0))
end function erf_count

pure elemental function derf_count(k,r,r0) result(count)
   use mctc_constants
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   count = -k/sqrtpi/r0*exp(-k**2*(r-r0)**2/r0**2)
end function derf_count

pure elemental function exp_count(k,r,r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   count =1.0_wp/(1.0_wp+exp(-k*(r0/r-1.0_wp)))
end function exp_count

pure elemental function dexp_count(k,r,r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   real(wp) :: expterm
   expterm=exp(-k*(r0/r-1._wp))
   count = (-k*r0*expterm)/(r**2*((expterm+1._wp)**2))
end function dexp_count

end module ncoord
