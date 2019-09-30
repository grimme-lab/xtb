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

module grid_module
   use iso_fortran_env, wp => real64

   private

   real(wp),private,parameter :: pi = 3.141592653589793_wp
   real(wp),private,parameter :: fourpi = 4.0_wp*pi
   real(wp),private,parameter :: third = 1.0_wp/3.0_wp
   real(wp),private,parameter :: third2 = 2.0_wp/3.0_wp
   real(wp),private,parameter :: small = 1e-10_wp

   ! grid type
   public :: tb_grid
   ! spline
   public :: spline
   ! lebedev grids
   public :: LD0006, LD0014, LD0026, LD0038, LD0050, LD0074, LD0086,&
              LD0110, LD0146, LD0170, LD0194, LD0230, LD0266, LD0302,&
              LD0350, LD0434, LD0590, LD0770, LD0974, LD1202, LD1454,&
              LD1730, LD2030, LD2354, LD2702, LD3074, LD3470, LD3890,&
              LD4334, LD4802, LD5294, LD5810, GEN_OH

   type tb_grid
      integer :: n
      integer, allocatable :: at(:)
      ! weights and coordinates of the gridpoints
      real(wp),allocatable :: w(:), x(:,:)
      real(wp),allocatable :: rho(:,:), b(:,:)
   contains
      procedure :: new => gengrid
      procedure :: delete => deallocate_grid
   end type tb_grid

contains

subroutine deallocate_grid(grid)
   implicit none
   class(tb_grid) :: grid
   grid%n = 0
   if (allocated(grid%w))   deallocate(grid%w)
   if (allocated(grid%x))   deallocate(grid%x)
   if (allocated(grid%rho)) deallocate(grid%rho)
   if (allocated(grid%b))   deallocate(grid%b)
end subroutine deallocate_grid

subroutine gengrid(grid,nat,at,xyz)
   implicit none
   ! GRID POINTS, NUCLEAR WEIGHTS, AND INTEGRATION WEIGHTS.

   class(tb_grid) :: grid
   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)

   real(wp) :: rr(nat*(1+nat)/2),rmid,r,r1,r2,hypr,vp0,vpsum,vpi
   integer  :: i,j,k,kk,err
   real(wp),allocatable :: rads(:),wrads(:),xang(:),yang(:),zang(:),wang(:)
   integer  :: nr,nang,ir,il,mang,mr
   real(wp) :: cutoff(nat,nat),x(3)
   real(wp),allocatable :: gridrl(:,:,:),gridx(:,:,:,:)

   call grid%delete

   ! interatomic distances
   rr = 0.0_wp
   do i = 1, nat
      rr(i*(i-1)/2) = 0.0_wp
      do j = 1, i-1
         rr(i*(i-1)/2+j) = norm2(xyz(:,i)-xyz(:,j))
      enddo
   enddo

   ! allocate space for the grid
   grid%n = 0
   do i = 1, nat
      grid%n = grid%n + at2nr(at(i)) * at2nang(at(i))
   enddo
   allocate( grid%w(grid%n),grid%x(3,grid%n), stat=err )

   ! allocate work arrays
   mr = -1
   mang = -1
   do i = 1, nat
      mang = max(mang,at2nang(at(i)))
      mr = max(mr,at2nr(at(i)))
   end do
   allocate( gridrl(mang,mr,nat),gridx(3,mang,mr,nat) )
   allocate( rads(mr),wrads(mr), stat=err )
   if (err.ne.0) &
   &  call raise('E','(gengrid) could not allocate memory for radial grid',1)
   allocate( xang(mang),yang(mang),zang(mang),wang(mang), stat=err )
   if (err.ne.0) &
   &  call raise('E','(gengrid) could not allocate memory for angular grid',1)

   ! Precompute the grid weights with multiple threads. The job has to be
   ! split in two because the nodes have to be positioned in the array in 
   ! the correct order 
   !$omp parallel default(none) shared(nat,at,xyz,rr,gridrl,gridx) &
   !$omp private(nr,nang,rmid,ir,r,il,x,j,k,r1,r2,hypr,cutoff,vp0,vpsum,vpi) &
   !$omp firstprivate(rads,wrads,xang,yang,zang,wang)
   !$omp do
   do i = 1, nat

      ! radial grid
      nr = at2nr(at(i))
      nang = at2nang(at(i))
      rmid = 1.0_wp/real(at(i),wp)**third
      call rgrid(nr,rmid,rads,wrads)

      ! angular grid

      select case(nang)
      case(   6); call ld0006(xang,yang,zang,wang,nang)
      case(  14); call ld0014(xang,yang,zang,wang,nang)
      case(  26); call ld0026(xang,yang,zang,wang,nang)
      case(  38); call ld0038(xang,yang,zang,wang,nang)
      case(  50); call ld0050(xang,yang,zang,wang,nang)
      case(  74); call ld0074(xang,yang,zang,wang,nang)
      case(  86); call ld0086(xang,yang,zang,wang,nang)
      case( 110); call ld0110(xang,yang,zang,wang,nang)
      case( 146); call ld0146(xang,yang,zang,wang,nang)
      case( 170); call ld0170(xang,yang,zang,wang,nang)
      case( 194); call ld0194(xang,yang,zang,wang,nang)
      case( 230); call ld0230(xang,yang,zang,wang,nang)
      case( 266); call ld0266(xang,yang,zang,wang,nang)
      case( 302); call ld0302(xang,yang,zang,wang,nang)
      case( 350); call ld0350(xang,yang,zang,wang,nang)
      case( 434); call ld0434(xang,yang,zang,wang,nang)
      case( 590); call ld0590(xang,yang,zang,wang,nang)
      case( 770); call ld0770(xang,yang,zang,wang,nang)
      case( 974); call ld0974(xang,yang,zang,wang,nang)
      case(1202); call ld1202(xang,yang,zang,wang,nang)
      case(1454); call ld1454(xang,yang,zang,wang,nang)
      case(1730); call ld1730(xang,yang,zang,wang,nang)
      case(2030); call ld2030(xang,yang,zang,wang,nang)
      case(2354); call ld2354(xang,yang,zang,wang,nang)
      case(2702); call ld2702(xang,yang,zang,wang,nang)
      case(3074); call ld3074(xang,yang,zang,wang,nang)
      case(3470); call ld3470(xang,yang,zang,wang,nang)
      case(3890); call ld3890(xang,yang,zang,wang,nang)
      case(4334); call ld4334(xang,yang,zang,wang,nang)
      case(4802); call ld4802(xang,yang,zang,wang,nang)
      case(5294); call ld5294(xang,yang,zang,wang,nang)
      case(5810); call ld5810(xang,yang,zang,wang,nang)
      case default
         call raise('E',"(gengrid) unknown value of nang",1)
      end select

      ! 3d grid, do not parallelize to get the nodes in order
      do ir = 1, nr
         r = rads(ir)
         do il = 1, nang
            x = xyz(:,i) + r * (/xang(il),yang(il),zang(il)/)
            do j = 1, nat
               do k = 1, j-1
                  r1 = norm2(x-xyz(:,j))
                  r2 = norm2(x-xyz(:,k))
                  hypr = (r1-r2) / rr(j*(j-1)/2+k)
                  hypr = 1.5_wp*hypr-0.5_wp*hypr**3
                  hypr = 1.5_wp*hypr-0.5_wp*hypr**3
                  hypr = 1.5_wp*hypr-0.5_wp*hypr**3
                  hypr = 1.5_wp*hypr-0.5_wp*hypr**3
                  cutoff(j,k) = 0.5_wp * (1.0_wp-hypr)
                  cutoff(k,j) = 0.5_wp * (1.0_wp+hypr)
               enddo
               cutoff(j,j) = 1.0_wp
            enddo
            cutoff(1,1) = 1.0_wp
            vp0 = 1.0_wp
            vpsum = 0.0_wp
            do j = 1, nat
               vp0=vp0*cutoff(i,j)
               vpi=1.0_wp
               do k = 1, nat
                  vpi = vpi * cutoff(j,k)
               enddo
               vpsum = vpsum + vpi
            enddo
            gridrl(il,ir,i) = vp0/vpsum * wrads(ir) * wang(il)
            gridx(:,il,ir,i) = x
         enddo
      enddo
   enddo
   !$omp enddo
   !$omp end parallel

   ! clean up
   deallocate(rads,wrads,xang,yang,zang,wang)

   ! fill the 3d grid
   kk = 0
   do i = 1, nat
      nr = at2nr(at(i))
      nang = at2nang(at(i))
      do ir = 1, nr
         do il = 1, nang
            kk = kk + 1
            grid%w(kk) = gridrl(il,ir,i)
            grid%x(:,kk) = gridx(:,il,ir,i)
         enddo
      enddo
   enddo

end subroutine gengrid

pure function get_xl(dd,l) result(xl)
   implicit none
   integer, intent(in)  :: l(3)
   real(wp),intent(in)  :: dd(3)
   real(wp) :: xl(3,0:2)
   integer  :: ix
   real(wp) :: ddi,xl2
   do ix = 1, 3
      ddi = dd(ix)
      select case(l(ix))
      case(0) ! s
         xl(ix,0) = 1.0_wp
         xl(ix,1) = 0.0_wp
         xl(ix,2) = 0.0_wp
      case(1) ! p
         xl(ix,0) = ddi
         xl(ix,1) = 1.0_wp
         xl(ix,2) = 0.0_wp
      case(2) ! d
         xl(ix,0) = ddi * ddi
         xl(ix,1) = 2.0_wp * ddi
         xl(ix,2) = 2.0_wp
      case(3) ! f
         xl(ix,0) = ddi * ddi * ddi
         xl(ix,1) = 3.0_wp * ddi * ddi
         xl(ix,2) = 6.0_wp * ddi
      case(4) ! g
         xl2 = ddi * ddi
         xl(ix,0) = xl2 * xl2
         xl(ix,1) = 4.0_wp * xl2 * ddi
         xl(ix,2) = 12.0_wp * xl2
      case(5) ! h
         xl2 = ddi * ddi
         xl(ix,0) = xl2 * xl2 * ddi
         xl(ix,1) = 5.0_wp * xl2 * xl2
         xl(ix,2) = 20.0_wp * xl2 * ddi
      end select
   enddo
   return
end function get_xl

function bhole(rho,quad,hnorm) result(b)
   implicit none

   real(wp),intent(in)  :: rho,quad,hnorm
   real(wp) :: b 

   real(wp) :: rhs,x0,shift,x,x1,expo,prefac,alf,f,df
   integer  :: i

   rhs=third2*(pi*rho/hnorm)**third2*rho/quad
   x0=2.0_wp
   shift=1.0_wp
   if(rhs.le.0.0_wp) then
      do i=1,16
         x=x0-shift
         call xfuncs(x,rhs,f,df)
         if(f.lt.0.0_wp) goto 11
         shift=0.1_wp*shift
      enddo
   else
      do i=1,16
         x=x0+shift
         call xfuncs(x,rhs,f,df)
         if(f.gt.0.0_wp) goto 11
         shift=0.1_wp*shift
      enddo
   endif
   call raise('E','(bhole) newton algorithm fails to initialize!',1)
   return
11 continue
   do i=1,100
      call xfuncs(x,rhs,f,df)
      x1=x-f/df
      if(abs(x1-x).lt.1.d-10) goto 22
      x=x1
   enddo
   call raise('E','(bhole) newton algorithm fails to converge!',1)
   return
22 x=x1
   expo=exp(-x)
   prefac=rho/expo
   alf=(8.0_wp*pi*prefac/hnorm)**third
   b=x/alf
end function bhole

pure elemental subroutine xfuncs(x,rhs,f,df)
   real(wp),intent(in)  :: x,rhs
   real(wp),intent(out) :: f,df

   real(wp) :: expo23

   expo23=exp(-2.0_wp/3.0_wp*x)
   f = x*expo23/(x-2.0_wp) - rhs
   df=2.0_wp/3.0_wp*(2.0_wp*x-x**2-3.0_wp)/(x-2.0_wp)**2*expo23
end subroutine xfuncs



pure elemental function at2nr(iat) result(nr)
   integer, intent(in) :: iat
   integer :: nr

   if (iat.gt.54) then
      nr = 140
   else if (iat.gt.36) then
      nr = 120
   else if (iat.gt.18) then
      nr = 100
   else if (iat.gt.10) then
      nr = 80
   else if (iat.gt.2) then
      nr = 60
   else
      nr = 40
   endif

end function at2nr

pure elemental function at2nang(iat) result(nang)
   integer, intent(in) :: iat
   integer :: nang

   nang = 194

endfunction at2nang

subroutine rgrid(n,rmid,r,wintr)
   implicit none
   !
   !     radial grid and integration weights,
   !     derivatives of variable r with respect to variable q.
   !
   !     the q-grid is uniform on the interval (0,+1). transformation is
   !
   !                        q
   !         r  =  rmid ---------
   !                    ( 1 - q )
   !
   !     also,
   !     3 and 7-point finite difference matrices for d2(wrt)r on q-grid.
   !

   integer, intent(in) :: n
   real(wp),intent(in) :: rmid
   real(wp),intent(out) :: r(n), wintr(n)

   real(wp) :: h, q
   integer  :: i

   h = 1.0_wp/real(n+1,wp)
   do i=1,n
      q=h*i
      r(i) = rmid * q / (1.0_wp-q)
      wintr(i) = fourpi * h * r(i)**2 * rmid / (1.0_wp-q)**2
   enddo

end subroutine rgrid

pure subroutine spline (h,y,a,b,c,n,bcr)
   !
   !     natural cubic spline through data points on the uniform r-grid.
   !     returns coefficients a,b,c defining the cubic polynomial to the
   !     right of each knot (including a(0),b(0),c(0)).
   !
   !     input:
   !         h    - interval between knots
   !         y    - array of ordinate values
   !         n    - the number of data points
   !         bcr  - value of function at r=infinity
   !
   !     output:
   !         a    - array of coefficients of x
   !         b    - array of coefficients of x**2
   !         c    - array of coefficients of x**3
   !
   implicit real(wp) (a-h,o-z)
   dimension y(0:n),a(0:n),b(0:n),c(0:n)
   dimension d(n),e(n),bb(n)
   intent(in) :: h, y, n, bcr
   intent(out) :: a, b, c

   thrd=1._wp/3._wp
   thrd2=2._wp*thrd
   hinv=1._wp/h
   h2inv=hinv*hinv

   do i=1,n-1
      d(i)=4._wp
      e(i)=1._wp
      bb(i)=3._wp*h2inv*(y(i-1)-2._wp*y(i)+y(i+1))
   enddo
   d(n)=4._wp
   e(n)=0._wp
   bb(n)=3._wp*h2inv*(y(n-1)-2._wp*y(n)+bcr)

   call dptsl(n,d,e,bb)

   a(0)=hinv*(y(1)-y(0))-thrd*h*bb(1)
   b(0)=0._wp
   c(0)=thrd*hinv*bb(1)
   do i=1,n-1
      a(i)=hinv*(y(i+1)-y(i))-thrd2*h*bb(i)-thrd*h*bb(i+1)
      b(i)=bb(i)
      c(i)=thrd*hinv*(bb(i+1)-bb(i))
   enddo
   a(n)=hinv*(bcr-y(n))-thrd2*h*bb(n)
   b(n)=bb(n)
   c(n)=-thrd*hinv*bb(n)

end subroutine spline

pure subroutine dptsl(n,d,e,b)
   integer, intent(in)    :: n
   real(wp),intent(inout) :: d(:)
   real(wp),intent(in)    :: e(:)
   real(wp),intent(inout) :: b(:)
!
!  dptsl given a positive definite tridiagonal matrix and a right
!  hand side will find the solution.
!
!  on entry
!
!     n        integer
!              is the order of the tridiagonal matrix.
!
!     d        real(wp)(n)
!              is the diagonal of the tridiagonal matrix.
!              on output d is destroyed.
!
!     e        real(wp)(n)
!              is the offdiagonal of the tridiagonal matrix.
!              e(1) through e(n-1) should contain the
!              offdiagonal.
!
!     b        real(wp)(n)
!              is the right hand side vector.
!
!  on return
!
!     b        contains the solution.
!
!  linpack. this version dated 08/14/78 .
!  jack dongarra, argonne national laboratory.
!
!  no externals
!  fortran mod
!
!  internal variables
!
   integer k,kbm1,ke,kf,kp1,nm1,nm1d2
   real(wp) t1,t2
!
!  check for 1 x 1 case
!
   if (n .eq. 1) then
      b(1) = b(1)/d(1)
      return
   endif
   nm1 = n - 1
   nm1d2 = nm1/2
   if (n .ne. 2) then
      kbm1 = n - 1
!
!     zero top half of subdiagonal and bottom half of
!     superdiagonal
!
      do k = 1, nm1d2
         t1 = e(k)/d(k)
         d(k+1) = d(k+1) - t1*e(k)
         b(k+1) = b(k+1) - t1*b(k)
         t2 = e(kbm1)/d(kbm1+1)
         d(kbm1) = d(kbm1) - t2*e(kbm1)
         b(kbm1) = b(kbm1) - t2*b(kbm1+1)
         kbm1 = kbm1 - 1
      enddo
   endif
   kp1 = nm1d2 + 1
!
!  clean up for possible 2 x 2 block at center
!
   if (mod(n,2) .eq. 0) then
      t1 = e(kp1)/d(kp1)
      d(kp1+1) = d(kp1+1) - t1*e(kp1)
      b(kp1+1) = b(kp1+1) - t1*b(kp1)
      kp1 = kp1 + 1
   endif
!
!  back solve starting at the center, going towards the top
!  and bottom
!
   b(kp1) = b(kp1)/d(kp1)
   if (n .ne. 2) then
      k = kp1 - 1
      ke = kp1 + nm1d2 - 1
      do kf = kp1, ke
         b(k) = (b(k) - e(k)*b(k+1))/d(k)
         b(kf+1) = (b(kf+1) - e(kf)*b(kf))/d(kf+1)
         k = k - 1
      enddo
   endif
   if (mod(n,2) .eq. 0) b(1) = (b(1) - e(1)*b(2))/d(1)
end subroutine

pure subroutine gen_oh(code, num, x, y, z, w, a, b, v)
   implicit logical(a-z)
   real(wp),intent(inout) :: x(*),y(*),z(*),w(*)
   real(wp),intent(inout) :: a,b,v
   integer, intent(in)    :: code
   integer, intent(inout) :: num
   real(wp) :: c
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated from C to fortran77 by hand.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   !vw
   !vw    Given a point on a sphere (specified by a and b), generate all
   !vw    the equivalent points under Oh symmetry, making grid points with
   !vw    weight v.
   !vw    The variable num is increased by the number of different points
   !vw    generated.
   !vw
   !vw    Depending on code, there are 6...48 different but equivalent
   !vw    points.
   !vw
   !vw    code=1:   (0,0,1) etc                                (  6 points)
   !vw    code=2:   (0,a,a) etc, a=1/sqrt(2)                   ( 12 points)
   !vw    code=3:   (a,a,a) etc, a=1/sqrt(3)                   (  8 points)
   !vw    code=4:   (a,a,b) etc, b=sqrt(1-2 a^2)               ( 24 points)
   !vw    code=5:   (a,b,0) etc, b=sqrt(1-a^2), a input        ( 24 points)
   !vw    code=6:   (a,b,c) etc, c=sqrt(1-a^2-b^2), a/b input  ( 48 points)
   !vw
   goto (1,2,3,4,5,6) code
   return
!  call raise('E','(Gen_Oh) Invalid Code',1)
1  continue
   a=1.0_wp
   x(1) =  a
   y(1) =  0.0_wp
   z(1) =  0.0_wp
   w(1) =  v
   x(2) = -a
   y(2) =  0.0_wp
   z(2) =  0.0_wp
   w(2) =  v
   x(3) =  0.0_wp
   y(3) =  a
   z(3) =  0.0_wp
   w(3) =  v
   x(4) =  0.0_wp
   y(4) = -a
   z(4) =  0.0_wp
   w(4) =  v
   x(5) =  0.0_wp
   y(5) =  0.0_wp
   z(5) =  a
   w(5) =  v
   x(6) =  0.0_wp
   y(6) =  0.0_wp
   z(6) = -a
   w(6) =  v
   num=num+6
   return
   !vw
2  continue
   a=sqrt(0.5_wp)
   x( 1) =  0_wp
   y( 1) =  a
   z( 1) =  a
   w( 1) =  v
   x( 2) =  0_wp
   y( 2) = -a
   z( 2) =  a
   w( 2) =  v
   x( 3) =  0_wp
   y( 3) =  a
   z( 3) = -a
   w( 3) =  v
   x( 4) =  0_wp
   y( 4) = -a
   z( 4) = -a
   w( 4) =  v
   x( 5) =  a
   y( 5) =  0_wp
   z( 5) =  a
   w( 5) =  v
   x( 6) = -a
   y( 6) =  0_wp
   z( 6) =  a
   w( 6) =  v
   x( 7) =  a
   y( 7) =  0_wp
   z( 7) = -a
   w( 7) =  v
   x( 8) = -a
   y( 8) =  0_wp
   z( 8) = -a
   w( 8) =  v
   x( 9) =  a
   y( 9) =  a
   z( 9) =  0_wp
   w( 9) =  v
   x(10) = -a
   y(10) =  a
   z(10) =  0_wp
   w(10) =  v
   x(11) =  a
   y(11) = -a
   z(11) =  0_wp
   w(11) =  v
   x(12) = -a
   y(12) = -a
   z(12) =  0_wp
   w(12) =  v
   num=num+12
   return
   !vw
3  continue
   a = sqrt(1._wp/3._wp)
   x(1) =  a
   y(1) =  a
   z(1) =  a
   w(1) =  v
   x(2) = -a
   y(2) =  a
   z(2) =  a
   w(2) =  v
   x(3) =  a
   y(3) = -a
   z(3) =  a
   w(3) =  v
   x(4) = -a
   y(4) = -a
   z(4) =  a
   w(4) =  v
   x(5) =  a
   y(5) =  a
   z(5) = -a
   w(5) =  v
   x(6) = -a
   y(6) =  a
   z(6) = -a
   w(6) =  v
   x(7) =  a
   y(7) = -a
   z(7) = -a
   w(7) =  v
   x(8) = -a
   y(8) = -a
   z(8) = -a
   w(8) =  v
   num=num+8
   return
   !vw
4  continue
   b = sqrt(1._wp - 2._wp*a*a)
   x( 1) =  a
   y( 1) =  a
   z( 1) =  b
   w( 1) =  v
   x( 2) = -a
   y( 2) =  a
   z( 2) =  b
   w( 2) =  v
   x( 3) =  a
   y( 3) = -a
   z( 3) =  b
   w( 3) =  v
   x( 4) = -a
   y( 4) = -a
   z( 4) =  b
   w( 4) =  v
   x( 5) =  a
   y( 5) =  a
   z( 5) = -b
   w( 5) =  v
   x( 6) = -a
   y( 6) =  a
   z( 6) = -b
   w( 6) =  v
   x( 7) =  a
   y( 7) = -a
   z( 7) = -b
   w( 7) =  v
   x( 8) = -a
   y( 8) = -a
   z( 8) = -b
   w( 8) =  v
   x( 9) =  a
   y( 9) =  b
   z( 9) =  a
   w( 9) =  v
   x(10) = -a
   y(10) =  b
   z(10) =  a
   w(10) =  v
   x(11) =  a
   y(11) = -b
   z(11) =  a
   w(11) =  v
   x(12) = -a
   y(12) = -b
   z(12) =  a
   w(12) =  v
   x(13) =  a
   y(13) =  b
   z(13) = -a
   w(13) =  v
   x(14) = -a
   y(14) =  b
   z(14) = -a
   w(14) =  v
   x(15) =  a
   y(15) = -b
   z(15) = -a
   w(15) =  v
   x(16) = -a
   y(16) = -b
   z(16) = -a
   w(16) =  v
   x(17) =  b
   y(17) =  a
   z(17) =  a
   w(17) =  v
   x(18) = -b
   y(18) =  a
   z(18) =  a
   w(18) =  v
   x(19) =  b
   y(19) = -a
   z(19) =  a
   w(19) =  v
   x(20) = -b
   y(20) = -a
   z(20) =  a
   w(20) =  v
   x(21) =  b
   y(21) =  a
   z(21) = -a
   w(21) =  v
   x(22) = -b
   y(22) =  a
   z(22) = -a
   w(22) =  v
   x(23) =  b
   y(23) = -a
   z(23) = -a
   w(23) =  v
   x(24) = -b
   y(24) = -a
   z(24) = -a
   w(24) =  v
   num=num+24
   return
   !vw
5  continue
   b=sqrt(1._wp-a*a)
   x( 1) =  a
   y( 1) =  b
   z( 1) =  0._wp
   w( 1) =  v
   x( 2) = -a
   y( 2) =  b
   z( 2) =  0._wp
   w( 2) =  v
   x( 3) =  a
   y( 3) = -b
   z( 3) =  0._wp
   w( 3) =  v
   x( 4) = -a
   y( 4) = -b
   z( 4) =  0._wp
   w( 4) =  v
   x( 5) =  b
   y( 5) =  a
   z( 5) =  0._wp
   w( 5) =  v
   x( 6) = -b
   y( 6) =  a
   z( 6) =  0._wp
   w( 6) =  v
   x( 7) =  b
   y( 7) = -a
   z( 7) =  0._wp
   w( 7) =  v
   x( 8) = -b
   y( 8) = -a
   z( 8) =  0._wp
   w( 8) =  v
   x( 9) =  a
   y( 9) =  0._wp
   z( 9) =  b
   w( 9) =  v
   x(10) = -a
   y(10) =  0._wp
   z(10) =  b
   w(10) =  v
   x(11) =  a
   y(11) =  0._wp
   z(11) = -b
   w(11) =  v
   x(12) = -a
   y(12) =  0._wp
   z(12) = -b
   w(12) =  v
   x(13) =  b
   y(13) =  0._wp
   z(13) =  a
   w(13) =  v
   x(14) = -b
   y(14) =  0._wp
   z(14) =  a
   w(14) =  v
   x(15) =  b
   y(15) =  0._wp
   z(15) = -a
   w(15) =  v
   x(16) = -b
   y(16) =  0._wp
   z(16) = -a
   w(16) =  v
   x(17) =  0._wp
   y(17) =  a
   z(17) =  b
   w(17) =  v
   x(18) =  0._wp
   y(18) = -a
   z(18) =  b
   w(18) =  v
   x(19) =  0._wp
   y(19) =  a
   z(19) = -b
   w(19) =  v
   x(20) =  0._wp
   y(20) = -a
   z(20) = -b
   w(20) =  v
   x(21) =  0._wp
   y(21) =  b
   z(21) =  a
   w(21) =  v
   x(22) =  0._wp
   y(22) = -b
   z(22) =  a
   w(22) =  v
   x(23) =  0._wp
   y(23) =  b
   z(23) = -a
   w(23) =  v
   x(24) =  0._wp
   y(24) = -b
   z(24) = -a
   w(24) =  v
   num=num+24
   return
   !vw
6  continue
   c=sqrt(1._wp - a*a - b*b)
   x( 1) =  a
   y( 1) =  b
   z( 1) =  c
   w( 1) =  v
   x( 2) = -a
   y( 2) =  b
   z( 2) =  c
   w( 2) =  v
   x( 3) =  a
   y( 3) = -b
   z( 3) =  c
   w( 3) =  v
   x( 4) = -a
   y( 4) = -b
   z( 4) =  c
   w( 4) =  v
   x( 5) =  a
   y( 5) =  b
   z( 5) = -c
   w( 5) =  v
   x( 6) = -a
   y( 6) =  b
   z( 6) = -c
   w( 6) =  v
   x( 7) =  a
   y( 7) = -b
   z( 7) = -c
   w( 7) =  v
   x( 8) = -a
   y( 8) = -b
   z( 8) = -c
   w( 8) =  v
   x( 9) =  a
   y( 9) =  c
   z( 9) =  b
   w( 9) =  v
   x(10) = -a
   y(10) =  c
   z(10) =  b
   w(10) =  v
   x(11) =  a
   y(11) = -c
   z(11) =  b
   w(11) =  v
   x(12) = -a
   y(12) = -c
   z(12) =  b
   w(12) =  v
   x(13) =  a
   y(13) =  c
   z(13) = -b
   w(13) =  v
   x(14) = -a
   y(14) =  c
   z(14) = -b
   w(14) =  v
   x(15) =  a
   y(15) = -c
   z(15) = -b
   w(15) =  v
   x(16) = -a
   y(16) = -c
   z(16) = -b
   w(16) =  v
   x(17) =  b
   y(17) =  a
   z(17) =  c
   w(17) =  v
   x(18) = -b
   y(18) =  a
   z(18) =  c
   w(18) =  v
   x(19) =  b
   y(19) = -a
   z(19) =  c
   w(19) =  v
   x(20) = -b
   y(20) = -a
   z(20) =  c
   w(20) =  v
   x(21) =  b
   y(21) =  a
   z(21) = -c
   w(21) =  v
   x(22) = -b
   y(22) =  a
   z(22) = -c
   w(22) =  v
   x(23) =  b
   y(23) = -a
   z(23) = -c
   w(23) =  v
   x(24) = -b
   y(24) = -a
   z(24) = -c
   w(24) =  v
   x(25) =  b
   y(25) =  c
   z(25) =  a
   w(25) =  v
   x(26) = -b
   y(26) =  c
   z(26) =  a
   w(26) =  v
   x(27) =  b
   y(27) = -c
   z(27) =  a
   w(27) =  v
   x(28) = -b
   y(28) = -c
   z(28) =  a
   w(28) =  v
   x(29) =  b
   y(29) =  c
   z(29) = -a
   w(29) =  v
   x(30) = -b
   y(30) =  c
   z(30) = -a
   w(30) =  v
   x(31) =  b
   y(31) = -c
   z(31) = -a
   w(31) =  v
   x(32) = -b
   y(32) = -c
   z(32) = -a
   w(32) =  v
   x(33) =  c
   y(33) =  a
   z(33) =  b
   w(33) =  v
   x(34) = -c
   y(34) =  a
   z(34) =  b
   w(34) =  v
   x(35) =  c
   y(35) = -a
   z(35) =  b
   w(35) =  v
   x(36) = -c
   y(36) = -a
   z(36) =  b
   w(36) =  v
   x(37) =  c
   y(37) =  a
   z(37) = -b
   w(37) =  v
   x(38) = -c
   y(38) =  a
   z(38) = -b
   w(38) =  v
   x(39) =  c
   y(39) = -a
   z(39) = -b
   w(39) =  v
   x(40) = -c
   y(40) = -a
   z(40) = -b
   w(40) =  v
   x(41) =  c
   y(41) =  b
   z(41) =  a
   w(41) =  v
   x(42) = -c
   y(42) =  b
   z(42) =  a
   w(42) =  v
   x(43) =  c
   y(43) = -b
   z(43) =  a
   w(43) =  v
   x(44) = -c
   y(44) = -b
   z(44) =  a
   w(44) =  v
   x(45) =  c
   y(45) =  b
   z(45) = -a
   w(45) =  v
   x(46) = -c
   y(46) =  b
   z(46) = -a
   w(46) =  v
   x(47) =  c
   y(47) = -b
   z(47) = -a
   w(47) =  v
   x(48) = -c
   y(48) = -b
   z(48) = -a
   w(48) =  v
   num=num+48
   return
end subroutine

SUBROUTINE LD0006(X,Y,Z,W,N)
   real(wp),intent(out) :: X(   6)
   real(wp),intent(out) :: Y(   6)
   real(wp),intent(out) :: Z(   6)
   real(wp),intent(out) :: W(   6)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV    6-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.1666666666666667D+0
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0014(X,Y,Z,W,N)
   real(wp),intent(out) :: X(  14)
   real(wp),intent(out) :: Y(  14)
   real(wp),intent(out) :: Z(  14)
   real(wp),intent(out) :: W(  14)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV   14-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.6666666666666667D-1
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.7500000000000000D-1
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0026(X,Y,Z,W,N)
   real(wp),intent(out) :: X(  26)
   real(wp),intent(out) :: Y(  26)
   real(wp),intent(out) :: Z(  26)
   real(wp),intent(out) :: W(  26)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV   26-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.4761904761904762D-1
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.3809523809523810D-1
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.3214285714285714D-1
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0038(X,Y,Z,W,N)
   real(wp),intent(out) :: X(  38)
   real(wp),intent(out) :: Y(  38)
   real(wp),intent(out) :: Z(  38)
   real(wp),intent(out) :: W(  38)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV   38-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.9523809523809524D-2
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.3214285714285714D-1
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4597008433809831D+0
   V=0.2857142857142857D-1
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0050(X,Y,Z,W,N)
   real(wp),intent(out) :: X(  50)
   real(wp),intent(out) :: Y(  50)
   real(wp),intent(out) :: Z(  50)
   real(wp),intent(out) :: W(  50)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV   50-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.1269841269841270D-1
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.2257495590828924D-1
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.2109375000000000D-1
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3015113445777636D+0
   V=0.2017333553791887D-1
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0074(X,Y,Z,W,N)
   real(wp),intent(out) :: X(  74)
   real(wp),intent(out) :: Y(  74)
   real(wp),intent(out) :: Z(  74)
   real(wp),intent(out) :: W(  74)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV   74-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.5130671797338464D-3
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.1660406956574204D-1
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=-0.2958603896103896D-1
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4803844614152614D+0
   V=0.2657620708215946D-1
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3207726489807764D+0
   V=0.1652217099371571D-1
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0086(X,Y,Z,W,N)
   real(wp),intent(out) :: X(  86)
   real(wp),intent(out) :: Y(  86)
   real(wp),intent(out) :: Z(  86)
   real(wp),intent(out) :: W(  86)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV   86-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.1154401154401154D-1
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.1194390908585628D-1
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3696028464541502D+0
   V=0.1111055571060340D-1
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6943540066026664D+0
   V=0.1187650129453714D-1
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3742430390903412D+0
   V=0.1181230374690448D-1
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0110(X,Y,Z,W,N)
   real(wp),intent(out) :: X( 110)
   real(wp),intent(out) :: Y( 110)
   real(wp),intent(out) :: Z( 110)
   real(wp),intent(out) :: W( 110)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV  110-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.3828270494937162D-2
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.9793737512487512D-2
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1851156353447362D+0
   V=0.8211737283191111D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6904210483822922D+0
   V=0.9942814891178103D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3956894730559419D+0
   V=0.9595471336070963D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4783690288121502D+0
   V=0.9694996361663028D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0146(X,Y,Z,W,N)
   real(wp),intent(out) :: X( 146)
   real(wp),intent(out) :: Y( 146)
   real(wp),intent(out) :: Z( 146)
   real(wp),intent(out) :: W( 146)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV  146-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.5996313688621381D-3
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.7372999718620756D-2
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.7210515360144488D-2
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6764410400114264D+0
   V=0.7116355493117555D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4174961227965453D+0
   V=0.6753829486314477D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1574676672039082D+0
   V=0.7574394159054034D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1403553811713183D+0
   B=0.4493328323269557D+0
   V=0.6991087353303262D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0170(X,Y,Z,W,N)
   real(wp),intent(out) :: X( 170)
   real(wp),intent(out) :: Y( 170)
   real(wp),intent(out) :: Z( 170)
   real(wp),intent(out) :: W( 170)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV  170-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.5544842902037365D-2
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.6071332770670752D-2
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.6383674773515093D-2
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2551252621114134D+0
   V=0.5183387587747790D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6743601460362766D+0
   V=0.6317929009813725D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4318910696719410D+0
   V=0.6201670006589077D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2613931360335988D+0
   V=0.5477143385137348D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4990453161796037D+0
   B=0.1446630744325115D+0
   V=0.5968383987681156D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0194(X,Y,Z,W,N)
   real(wp),intent(out) :: X( 194)
   real(wp),intent(out) :: Y( 194)
   real(wp),intent(out) :: Z( 194)
   real(wp),intent(out) :: W( 194)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV  194-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.1782340447244611D-2
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.5716905949977102D-2
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.5573383178848738D-2
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6712973442695226D+0
   V=0.5608704082587997D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2892465627575439D+0
   V=0.5158237711805383D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4446933178717437D+0
   V=0.5518771467273614D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1299335447650067D+0
   V=0.4106777028169394D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3457702197611283D+0
   V=0.5051846064614808D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1590417105383530D+0
   B=0.8360360154824589D+0
   V=0.5530248916233094D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0230(X,Y,Z,W,N)
   real(wp),intent(out) :: X( 230)
   real(wp),intent(out) :: Y( 230)
   real(wp),intent(out) :: Z( 230)
   real(wp),intent(out) :: W( 230)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV  230-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=-0.5522639919727325D-1
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.4450274607445226D-2
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4492044687397611D+0
   V=0.4496841067921404D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2520419490210201D+0
   V=0.5049153450478750D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6981906658447242D+0
   V=0.3976408018051883D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6587405243460960D+0
   V=0.4401400650381014D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4038544050097660D-1
   V=0.1724544350544401D-1
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5823842309715585D+0
   V=0.4231083095357343D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3545877390518688D+0
   V=0.5198069864064399D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2272181808998187D+0
   B=0.4864661535886647D+0
   V=0.4695720972568883D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0266(X,Y,Z,W,N)
   real(wp),intent(out) :: X( 266)
   real(wp),intent(out) :: Y( 266)
   real(wp),intent(out) :: Z( 266)
   real(wp),intent(out) :: W( 266)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV  266-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=-0.1313769127326952D-2
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=-0.2522728704859336D-2
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.4186853881700583D-2
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7039373391585475D+0
   V=0.5315167977810885D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1012526248572414D+0
   V=0.4047142377086219D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4647448726420539D+0
   V=0.4112482394406990D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3277420654971629D+0
   V=0.3595584899758782D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6620338663699974D+0
   V=0.4256131351428158D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8506508083520399D+0
   V=0.4229582700647240D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3233484542692899D+0
   B=0.1153112011009701D+0
   V=0.4080914225780505D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2314790158712601D+0
   B=0.5244939240922365D+0
   V=0.4071467593830964D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0302(X,Y,Z,W,N)
   real(wp),intent(out) :: X( 302)
   real(wp),intent(out) :: Y( 302)
   real(wp),intent(out) :: Z( 302)
   real(wp),intent(out) :: W( 302)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV  302-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.8545911725128148D-3
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.3599119285025571D-2
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3515640345570105D+0
   V=0.3449788424305883D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6566329410219612D+0
   V=0.3604822601419882D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4729054132581005D+0
   V=0.3576729661743367D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9618308522614784D-1
   V=0.2352101413689164D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2219645236294178D+0
   V=0.3108953122413675D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7011766416089545D+0
   V=0.3650045807677255D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2644152887060663D+0
   V=0.2982344963171804D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5718955891878961D+0
   V=0.3600820932216460D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2510034751770465D+0
   B=0.8000727494073952D+0
   V=0.3571540554273387D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1233548532583327D+0
   B=0.4127724083168531D+0
   V=0.3392312205006170D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0350(X,Y,Z,W,N)
   real(wp),intent(out) :: X( 350)
   real(wp),intent(out) :: Y( 350)
   real(wp),intent(out) :: Z( 350)
   real(wp),intent(out) :: W( 350)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV  350-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.3006796749453936D-2
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.3050627745650771D-2
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7068965463912316D+0
   V=0.1621104600288991D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4794682625712025D+0
   V=0.3005701484901752D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1927533154878019D+0
   V=0.2990992529653774D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6930357961327123D+0
   V=0.2982170644107595D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3608302115520091D+0
   V=0.2721564237310992D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6498486161496169D+0
   V=0.3033513795811141D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1932945013230339D+0
   V=0.3007949555218533D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3800494919899303D+0
   V=0.2881964603055307D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2899558825499574D+0
   B=0.7934537856582316D+0
   V=0.2958357626535696D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9684121455103957D-1
   B=0.8280801506686862D+0
   V=0.3036020026407088D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1833434647041659D+0
   B=0.9074658265305127D+0
   V=0.2832187403926303D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0434(X,Y,Z,W,N)
   real(wp),intent(out) :: X( 434)
   real(wp),intent(out) :: Y( 434)
   real(wp),intent(out) :: Z( 434)
   real(wp),intent(out) :: W( 434)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV  434-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.5265897968224436D-3
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.2548219972002607D-2
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.2512317418927307D-2
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6909346307509111D+0
   V=0.2530403801186355D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1774836054609158D+0
   V=0.2014279020918528D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4914342637784746D+0
   V=0.2501725168402936D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6456664707424256D+0
   V=0.2513267174597564D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2861289010307638D+0
   V=0.2302694782227416D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7568084367178018D-1
   V=0.1462495621594614D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3927259763368002D+0
   V=0.2445373437312980D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8818132877794288D+0
   V=0.2417442375638981D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9776428111182649D+0
   V=0.1910951282179532D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2054823696403044D+0
   B=0.8689460322872412D+0
   V=0.2416930044324775D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5905157048925271D+0
   B=0.7999278543857286D+0
   V=0.2512236854563495D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5550152361076807D+0
   B=0.7717462626915901D+0
   V=0.2496644054553086D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9371809858553722D+0
   B=0.3344363145343455D+0
   V=0.2236607760437849D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0590(X,Y,Z,W,N)
   real(wp),intent(out) :: X( 590)
   real(wp),intent(out) :: Y( 590)
   real(wp),intent(out) :: Z( 590)
   real(wp),intent(out) :: W( 590)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV  590-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.3095121295306187D-3
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.1852379698597489D-2
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7040954938227469D+0
   V=0.1871790639277744D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6807744066455243D+0
   V=0.1858812585438317D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6372546939258752D+0
   V=0.1852028828296213D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5044419707800358D+0
   V=0.1846715956151242D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4215761784010967D+0
   V=0.1818471778162769D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3317920736472123D+0
   V=0.1749564657281154D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2384736701421887D+0
   V=0.1617210647254411D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1459036449157763D+0
   V=0.1384737234851692D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6095034115507196D-1
   V=0.9764331165051050D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6116843442009876D+0
   V=0.1857161196774078D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3964755348199858D+0
   V=0.1705153996395864D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1724782009907724D+0
   V=0.1300321685886048D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5610263808622060D+0
   B=0.3518280927733519D+0
   V=0.1842866472905286D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4742392842551980D+0
   B=0.2634716655937950D+0
   V=0.1802658934377451D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5984126497885380D+0
   B=0.1816640840360209D+0
   V=0.1849830560443660D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3791035407695563D+0
   B=0.1720795225656878D+0
   V=0.1713904507106709D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2778673190586244D+0
   B=0.8213021581932511D-1
   V=0.1555213603396808D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5033564271075117D+0
   B=0.8999205842074875D-1
   V=0.1802239128008525D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0770(X,Y,Z,W,N)
   real(wp),intent(out) :: X( 770)
   real(wp),intent(out) :: Y( 770)
   real(wp),intent(out) :: Z( 770)
   real(wp),intent(out) :: W( 770)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV  770-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.2192942088181184D-3
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.1436433617319080D-2
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.1421940344335877D-2
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5087204410502360D-1
   V=0.6798123511050502D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1228198790178831D+0
   V=0.9913184235294912D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2026890814408786D+0
   V=0.1180207833238949D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2847745156464294D+0
   V=0.1296599602080921D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3656719078978026D+0
   V=0.1365871427428316D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4428264886713469D+0
   V=0.1402988604775325D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5140619627249735D+0
   V=0.1418645563595609D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6306401219166803D+0
   V=0.1421376741851662D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6716883332022612D+0
   V=0.1423996475490962D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6979792685336881D+0
   V=0.1431554042178567D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1446865674195309D+0
   V=0.9254401499865368D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3390263475411216D+0
   V=0.1250239995053509D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5335804651263506D+0
   V=0.1394365843329230D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6944024393349413D-1
   B=0.2355187894242326D+0
   V=0.1127089094671749D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2269004109529460D+0
   B=0.4102182474045730D+0
   V=0.1345753760910670D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8025574607775339D-1
   B=0.6214302417481605D+0
   V=0.1424957283316783D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1467999527896572D+0
   B=0.3245284345717394D+0
   V=0.1261523341237750D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1571507769824727D+0
   B=0.5224482189696630D+0
   V=0.1392547106052696D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2365702993157246D+0
   B=0.6017546634089558D+0
   V=0.1418761677877656D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7714815866765732D-1
   B=0.4346575516141163D+0
   V=0.1338366684479554D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3062936666210730D+0
   B=0.4908826589037616D+0
   V=0.1393700862676131D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3822477379524787D+0
   B=0.5648768149099500D+0
   V=0.1415914757466932D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD0974(X,Y,Z,W,N)
   real(wp),intent(out) :: X( 974)
   real(wp),intent(out) :: Y( 974)
   real(wp),intent(out) :: Z( 974)
   real(wp),intent(out) :: W( 974)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV  974-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.1438294190527431D-3
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.1125772288287004D-2
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4292963545341347D-1
   V=0.4948029341949241D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1051426854086404D+0
   V=0.7357990109125470D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1750024867623087D+0
   V=0.8889132771304384D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2477653379650257D+0
   V=0.9888347838921435D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3206567123955957D+0
   V=0.1053299681709471D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3916520749849983D+0
   V=0.1092778807014578D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4590825874187624D+0
   V=0.1114389394063227D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5214563888415861D+0
   V=0.1123724788051555D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6253170244654199D+0
   V=0.1125239325243814D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6637926744523170D+0
   V=0.1126153271815905D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6910410398498301D+0
   V=0.1130286931123841D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7052907007457760D+0
   V=0.1134986534363955D-2
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1236686762657990D+0
   V=0.6823367927109931D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2940777114468387D+0
   V=0.9454158160447096D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4697753849207649D+0
   V=0.1074429975385679D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6334563241139567D+0
   V=0.1129300086569132D-2
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5974048614181342D-1
   B=0.2029128752777523D+0
   V=0.8436884500901954D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1375760408473636D+0
   B=0.4602621942484054D+0
   V=0.1075255720448885D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3391016526336286D+0
   B=0.5030673999662036D+0
   V=0.1108577236864462D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1271675191439820D+0
   B=0.2817606422442134D+0
   V=0.9566475323783357D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2693120740413512D+0
   B=0.4331561291720157D+0
   V=0.1080663250717391D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1419786452601918D+0
   B=0.6256167358580814D+0
   V=0.1126797131196295D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6709284600738255D-1
   B=0.3798395216859157D+0
   V=0.1022568715358061D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7057738183256172D-1
   B=0.5517505421423520D+0
   V=0.1108960267713108D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2783888477882155D+0
   B=0.6029619156159187D+0
   V=0.1122790653435766D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1979578938917407D+0
   B=0.3589606329589096D+0
   V=0.1032401847117460D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2087307061103274D+0
   B=0.5348666438135476D+0
   V=0.1107249382283854D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4055122137872836D+0
   B=0.5674997546074373D+0
   V=0.1121780048519972D-2
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD1202(X,Y,Z,W,N)
   real(wp),intent(out) :: X(1202)
   real(wp),intent(out) :: Y(1202)
   real(wp),intent(out) :: Z(1202)
   real(wp),intent(out) :: W(1202)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV 1202-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.1105189233267572D-3
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.9205232738090741D-3
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.9133159786443561D-3
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3712636449657089D-1
   V=0.3690421898017899D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9140060412262223D-1
   V=0.5603990928680660D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1531077852469906D+0
   V=0.6865297629282609D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2180928891660612D+0
   V=0.7720338551145630D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2839874532200175D+0
   V=0.8301545958894795D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3491177600963764D+0
   V=0.8686692550179628D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4121431461444309D+0
   V=0.8927076285846890D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4718993627149127D+0
   V=0.9060820238568219D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5273145452842337D+0
   V=0.9119777254940867D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6209475332444019D+0
   V=0.9128720138604181D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6569722711857291D+0
   V=0.9130714935691735D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6841788309070143D+0
   V=0.9152873784554116D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7012604330123631D+0
   V=0.9187436274321654D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1072382215478166D+0
   V=0.5176977312965694D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2582068959496968D+0
   V=0.7331143682101417D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4172752955306717D+0
   V=0.8463232836379928D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5700366911792503D+0
   V=0.9031122694253992D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9827986018263947D+0
   B=0.1771774022615325D+0
   V=0.6485778453163257D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9624249230326228D+0
   B=0.2475716463426288D+0
   V=0.7435030910982369D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9402007994128811D+0
   B=0.3354616289066489D+0
   V=0.7998527891839054D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9320822040143202D+0
   B=0.3173615246611977D+0
   V=0.8101731497468018D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9043674199393299D+0
   B=0.4090268427085357D+0
   V=0.8483389574594331D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8912407560074747D+0
   B=0.3854291150669224D+0
   V=0.8556299257311812D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8676435628462708D+0
   B=0.4932221184851285D+0
   V=0.8803208679738260D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8581979986041619D+0
   B=0.4785320675922435D+0
   V=0.8811048182425720D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8396753624049856D+0
   B=0.4507422593157064D+0
   V=0.8850282341265444D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8165288564022188D+0
   B=0.5632123020762100D+0
   V=0.9021342299040653D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8015469370783529D+0
   B=0.5434303569693900D+0
   V=0.9010091677105086D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7773563069070351D+0
   B=0.5123518486419871D+0
   V=0.9022692938426915D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7661621213900394D+0
   B=0.6394279634749102D+0
   V=0.9158016174693465D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7553584143533510D+0
   B=0.6269805509024392D+0
   V=0.9131578003189435D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7344305757559503D+0
   B=0.6031161693096310D+0
   V=0.9107813579482705D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7043837184021765D+0
   B=0.5693702498468441D+0
   V=0.9105760258970126D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD1454(X,Y,Z,W,N)
   real(wp),intent(out) :: X(1454)
   real(wp),intent(out) :: Y(1454)
   real(wp),intent(out) :: Z(1454)
   real(wp),intent(out) :: W(1454)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV 1454-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.7777160743261247D-4
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.7557646413004701D-3
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3229290663413854D-1
   V=0.2841633806090617D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8036733271462222D-1
   V=0.4374419127053555D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1354289960531653D+0
   V=0.5417174740872172D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1938963861114426D+0
   V=0.6148000891358593D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2537343715011275D+0
   V=0.6664394485800705D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3135251434752570D+0
   V=0.7025039356923220D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3721558339375338D+0
   V=0.7268511789249627D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4286809575195696D+0
   V=0.7422637534208629D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4822510128282994D+0
   V=0.7509545035841214D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5320679333566263D+0
   V=0.7548535057718401D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6172998195394274D+0
   V=0.7554088969774001D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6510679849127481D+0
   V=0.7553147174442808D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6777315251687360D+0
   V=0.7564767653292297D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6963109410648741D+0
   V=0.7587991808518730D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7058935009831749D+0
   V=0.7608261832033027D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9955546194091857D+0
   V=0.4021680447874916D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9734115901794209D+0
   V=0.5804871793945964D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9275693732388626D+0
   V=0.6792151955945159D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8568022422795103D+0
   V=0.7336741211286294D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7623495553719372D+0
   V=0.7581866300989608D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5707522908892223D+0
   B=0.4387028039889501D+0
   V=0.7538257859800743D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5196463388403083D+0
   B=0.3858908414762617D+0
   V=0.7483517247053123D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4646337531215351D+0
   B=0.3301937372343854D+0
   V=0.7371763661112059D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4063901697557691D+0
   B=0.2725423573563777D+0
   V=0.7183448895756934D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3456329466643087D+0
   B=0.2139510237495250D+0
   V=0.6895815529822191D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2831395121050332D+0
   B=0.1555922309786647D+0
   V=0.6480105801792886D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2197682022925330D+0
   B=0.9892878979686097D-1
   V=0.5897558896594636D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1564696098650355D+0
   B=0.4598642910675510D-1
   V=0.5095708849247346D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6027356673721295D+0
   B=0.3376625140173426D+0
   V=0.7536906428909755D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5496032320255096D+0
   B=0.2822301309727988D+0
   V=0.7472505965575118D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4921707755234567D+0
   B=0.2248632342592540D+0
   V=0.7343017132279698D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4309422998598483D+0
   B=0.1666224723456479D+0
   V=0.7130871582177445D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3664108182313672D+0
   B=0.1086964901822169D+0
   V=0.6817022032112776D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2990189057758436D+0
   B=0.5251989784120085D-1
   V=0.6380941145604121D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6268724013144998D+0
   B=0.2297523657550023D+0
   V=0.7550381377920310D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5707324144834607D+0
   B=0.1723080607093800D+0
   V=0.7478646640144802D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5096360901960365D+0
   B=0.1140238465390513D+0
   V=0.7335918720601220D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4438729938312456D+0
   B=0.5611522095882537D-1
   V=0.7110120527658118D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6419978471082389D+0
   B=0.1164174423140873D+0
   V=0.7571363978689501D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5817218061802611D+0
   B=0.5797589531445219D-1
   V=0.7489908329079234D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD1730(X,Y,Z,W,N)
   real(wp),intent(out) :: X(1730)
   real(wp),intent(out) :: Y(1730)
   real(wp),intent(out) :: Z(1730)
   real(wp),intent(out) :: W(1730)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV 1730-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.6309049437420976D-4
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.6398287705571748D-3
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.6357185073530720D-3
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2860923126194662D-1
   V=0.2221207162188168D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7142556767711522D-1
   V=0.3475784022286848D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1209199540995559D+0
   V=0.4350742443589804D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1738673106594379D+0
   V=0.4978569136522127D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2284645438467734D+0
   V=0.5435036221998053D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2834807671701512D+0
   V=0.5765913388219542D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3379680145467339D+0
   V=0.6001200359226003D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3911355454819537D+0
   V=0.6162178172717512D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4422860353001403D+0
   V=0.6265218152438485D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4907781568726057D+0
   V=0.6323987160974212D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5360006153211468D+0
   V=0.6350767851540569D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6142105973596603D+0
   V=0.6354362775297107D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6459300387977504D+0
   V=0.6352302462706235D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6718056125089225D+0
   V=0.6358117881417972D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6910888533186254D+0
   V=0.6373101590310117D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7030467416823252D+0
   V=0.6390428961368665D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8354951166354646D-1
   V=0.3186913449946576D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2050143009099486D+0
   V=0.4678028558591711D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3370208290706637D+0
   V=0.5538829697598626D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4689051484233963D+0
   V=0.6044475907190476D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5939400424557334D+0
   V=0.6313575103509012D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1394983311832261D+0
   B=0.4097581162050343D-1
   V=0.4078626431855630D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1967999180485014D+0
   B=0.8851987391293348D-1
   V=0.4759933057812725D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2546183732548967D+0
   B=0.1397680182969819D+0
   V=0.5268151186413440D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3121281074713875D+0
   B=0.1929452542226526D+0
   V=0.5643048560507316D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3685981078502492D+0
   B=0.2467898337061562D+0
   V=0.5914501076613073D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4233760321547856D+0
   B=0.3003104124785409D+0
   V=0.6104561257874195D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4758671236059246D+0
   B=0.3526684328175033D+0
   V=0.6230252860707806D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5255178579796463D+0
   B=0.4031134861145713D+0
   V=0.6305618761760796D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5718025633734589D+0
   B=0.4509426448342351D+0
   V=0.6343092767597889D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2686927772723415D+0
   B=0.4711322502423248D-1
   V=0.5176268945737826D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3306006819904809D+0
   B=0.9784487303942695D-1
   V=0.5564840313313692D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3904906850594983D+0
   B=0.1505395810025273D+0
   V=0.5856426671038980D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4479957951904390D+0
   B=0.2039728156296050D+0
   V=0.6066386925777091D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5027076848919780D+0
   B=0.2571529941121107D+0
   V=0.6208824962234458D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5542087392260217D+0
   B=0.3092191375815670D+0
   V=0.6296314297822907D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6020850887375187D+0
   B=0.3593807506130276D+0
   V=0.6340423756791859D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4019851409179594D+0
   B=0.5063389934378671D-1
   V=0.5829627677107342D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4635614567449800D+0
   B=0.1032422269160612D+0
   V=0.6048693376081110D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5215860931591575D+0
   B=0.1566322094006254D+0
   V=0.6202362317732461D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5758202499099271D+0
   B=0.2098082827491099D+0
   V=0.6299005328403779D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6259893683876795D+0
   B=0.2618824114553391D+0
   V=0.6347722390609353D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5313795124811891D+0
   B=0.5263245019338556D-1
   V=0.6203778981238834D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5893317955931995D+0
   B=0.1061059730982005D+0
   V=0.6308414671239979D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6426246321215801D+0
   B=0.1594171564034221D+0
   V=0.6362706466959498D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6511904367376113D+0
   B=0.5354789536565540D-1
   V=0.6375414170333233D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD2030(X,Y,Z,W,N)
   real(wp),intent(out) :: X(2030)
   real(wp),intent(out) :: Y(2030)
   real(wp),intent(out) :: Z(2030)
   real(wp),intent(out) :: W(2030)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV 2030-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.4656031899197431D-4
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.5421549195295507D-3
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2540835336814348D-1
   V=0.1778522133346553D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6399322800504915D-1
   V=0.2811325405682796D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1088269469804125D+0
   V=0.3548896312631459D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1570670798818287D+0
   V=0.4090310897173364D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2071163932282514D+0
   V=0.4493286134169965D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2578914044450844D+0
   V=0.4793728447962723D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3085687558169623D+0
   V=0.5015415319164265D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3584719706267024D+0
   V=0.5175127372677937D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4070135594428709D+0
   V=0.5285522262081019D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4536618626222638D+0
   V=0.5356832703713962D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4979195686463577D+0
   V=0.5397914736175170D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5393075111126999D+0
   V=0.5416899441599930D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6115617676843916D+0
   V=0.5419308476889938D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6414308435160159D+0
   V=0.5416936902030596D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6664099412721607D+0
   V=0.5419544338703164D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6859161771214913D+0
   V=0.5428983656630975D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6993625593503890D+0
   V=0.5442286500098193D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7062393387719380D+0
   V=0.5452250345057301D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7479028168349763D-1
   V=0.2568002497728530D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1848951153969366D+0
   V=0.3827211700292145D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3059529066581305D+0
   V=0.4579491561917824D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4285556101021362D+0
   V=0.5042003969083574D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5468758653496526D+0
   V=0.5312708889976025D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6565821978343439D+0
   V=0.5438401790747117D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1253901572367117D+0
   B=0.3681917226439641D-1
   V=0.3316041873197344D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1775721510383941D+0
   B=0.7982487607213301D-1
   V=0.3899113567153771D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2305693358216114D+0
   B=0.1264640966592335D+0
   V=0.4343343327201309D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2836502845992063D+0
   B=0.1751585683418957D+0
   V=0.4679415262318919D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3361794746232590D+0
   B=0.2247995907632670D+0
   V=0.4930847981631031D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3875979172264824D+0
   B=0.2745299257422246D+0
   V=0.5115031867540091D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4374019316999074D+0
   B=0.3236373482441118D+0
   V=0.5245217148457367D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4851275843340022D+0
   B=0.3714967859436741D+0
   V=0.5332041499895321D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5303391803806868D+0
   B=0.4175353646321745D+0
   V=0.5384583126021542D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5726197380596287D+0
   B=0.4612084406355461D+0
   V=0.5411067210798852D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2431520732564863D+0
   B=0.4258040133043952D-1
   V=0.4259797391468714D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3002096800895869D+0
   B=0.8869424306722721D-1
   V=0.4604931368460021D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3558554457457432D+0
   B=0.1368811706510655D+0
   V=0.4871814878255202D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4097782537048887D+0
   B=0.1860739985015033D+0
   V=0.5072242910074885D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4616337666067458D+0
   B=0.2354235077395853D+0
   V=0.5217069845235350D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5110707008417874D+0
   B=0.2842074921347011D+0
   V=0.5315785966280310D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5577415286163795D+0
   B=0.3317784414984102D+0
   V=0.5376833708758905D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6013060431366950D+0
   B=0.3775299002040700D+0
   V=0.5408032092069521D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3661596767261781D+0
   B=0.4599367887164592D-1
   V=0.4842744917904866D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4237633153506581D+0
   B=0.9404893773654421D-1
   V=0.5048926076188130D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4786328454658452D+0
   B=0.1431377109091971D+0
   V=0.5202607980478373D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5305702076789774D+0
   B=0.1924186388843570D+0
   V=0.5309932388325743D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5793436224231788D+0
   B=0.2411590944775190D+0
   V=0.5377419770895208D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6247069017094747D+0
   B=0.2886871491583605D+0
   V=0.5411696331677717D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4874315552535204D+0
   B=0.4804978774953206D-1
   V=0.5197996293282420D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5427337322059053D+0
   B=0.9716857199366665D-1
   V=0.5311120836622945D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5943493747246700D+0
   B=0.1465205839795055D+0
   V=0.5384309319956951D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6421314033564943D+0
   B=0.1953579449803574D+0
   V=0.5421859504051886D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6020628374713980D+0
   B=0.4916375015738108D-1
   V=0.5390948355046314D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6529222529856881D+0
   B=0.9861621540127005D-1
   V=0.5433312705027845D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD2354(X,Y,Z,W,N)
   real(wp),intent(out) :: X(2354)
   real(wp),intent(out) :: Y(2354)
   real(wp),intent(out) :: Z(2354)
   real(wp),intent(out) :: W(2354)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV 2354-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.3922616270665292D-4
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.4703831750854424D-3
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.4678202801282136D-3
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2290024646530589D-1
   V=0.1437832228979900D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5779086652271284D-1
   V=0.2303572493577644D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9863103576375984D-1
   V=0.2933110752447454D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1428155792982185D+0
   V=0.3402905998359838D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1888978116601463D+0
   V=0.3759138466870372D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2359091682970210D+0
   V=0.4030638447899798D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2831228833706171D+0
   V=0.4236591432242211D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3299495857966693D+0
   V=0.4390522656946746D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3758840802660796D+0
   V=0.4502523466626247D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4204751831009480D+0
   V=0.4580577727783541D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4633068518751051D+0
   V=0.4631391616615899D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5039849474507313D+0
   V=0.4660928953698676D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5421265793440747D+0
   V=0.4674751807936953D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6092660230557310D+0
   V=0.4676414903932920D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6374654204984869D+0
   V=0.4674086492347870D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6615136472609892D+0
   V=0.4674928539483207D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6809487285958127D+0
   V=0.4680748979686447D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6952980021665196D+0
   V=0.4690449806389040D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7041245497695400D+0
   V=0.4699877075860818D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6744033088306065D-1
   V=0.2099942281069176D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1678684485334166D+0
   V=0.3172269150712804D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2793559049539613D+0
   V=0.3832051358546523D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3935264218057639D+0
   V=0.4252193818146985D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5052629268232558D+0
   V=0.4513807963755000D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6107905315437531D+0
   V=0.4657797469114178D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1135081039843524D+0
   B=0.3331954884662588D-1
   V=0.2733362800522836D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1612866626099378D+0
   B=0.7247167465436538D-1
   V=0.3235485368463559D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2100786550168205D+0
   B=0.1151539110849745D+0
   V=0.3624908726013453D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2592282009459942D+0
   B=0.1599491097143677D+0
   V=0.3925540070712828D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3081740561320203D+0
   B=0.2058699956028027D+0
   V=0.4156129781116235D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3564289781578164D+0
   B=0.2521624953502911D+0
   V=0.4330644984623263D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4035587288240703D+0
   B=0.2982090785797674D+0
   V=0.4459677725921312D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4491671196373903D+0
   B=0.3434762087235733D+0
   V=0.4551593004456795D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4928854782917489D+0
   B=0.3874831357203437D+0
   V=0.4613341462749918D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5343646791958988D+0
   B=0.4297814821746926D+0
   V=0.4651019618269806D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5732683216530990D+0
   B=0.4699402260943537D+0
   V=0.4670249536100625D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2214131583218986D+0
   B=0.3873602040643895D-1
   V=0.3549555576441708D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2741796504750071D+0
   B=0.8089496256902013D-1
   V=0.3856108245249010D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3259797439149485D+0
   B=0.1251732177620872D+0
   V=0.4098622845756882D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3765441148826891D+0
   B=0.1706260286403185D+0
   V=0.4286328604268950D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4255773574530558D+0
   B=0.2165115147300408D+0
   V=0.4427802198993945D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4727795117058430D+0
   B=0.2622089812225259D+0
   V=0.4530473511488561D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5178546895819012D+0
   B=0.3071721431296201D+0
   V=0.4600805475703138D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5605141192097460D+0
   B=0.3508998998801138D+0
   V=0.4644599059958017D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6004763319352512D+0
   B=0.3929160876166931D+0
   V=0.4667274455712508D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3352842634946949D+0
   B=0.4202563457288019D-1
   V=0.4069360518020356D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3891971629814670D+0
   B=0.8614309758870850D-1
   V=0.4260442819919195D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4409875565542281D+0
   B=0.1314500879380001D+0
   V=0.4408678508029063D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4904893058592484D+0
   B=0.1772189657383859D+0
   V=0.4518748115548597D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5375056138769549D+0
   B=0.2228277110050294D+0
   V=0.4595564875375116D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5818255708669969D+0
   B=0.2677179935014386D+0
   V=0.4643988774315846D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6232334858144959D+0
   B=0.3113675035544165D+0
   V=0.4668827491646946D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4489485354492058D+0
   B=0.4409162378368174D-1
   V=0.4400541823741973D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5015136875933150D+0
   B=0.8939009917748489D-1
   V=0.4514512890193797D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5511300550512623D+0
   B=0.1351806029383365D+0
   V=0.4596198627347549D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5976720409858000D+0
   B=0.1808370355053196D+0
   V=0.4648659016801781D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6409956378989354D+0
   B=0.2257852192301602D+0
   V=0.4675502017157673D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5581222330827514D+0
   B=0.4532173421637160D-1
   V=0.4598494476455523D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6074705984161695D+0
   B=0.9117488031840314D-1
   V=0.4654916955152048D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6532272537379033D+0
   B=0.1369294213140155D+0
   V=0.4684709779505137D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6594761494500487D+0
   B=0.4589901487275583D-1
   V=0.4691445539106986D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD2702(X,Y,Z,W,N)
   real(wp),intent(out) :: X(2702)
   real(wp),intent(out) :: Y(2702)
   real(wp),intent(out) :: Z(2702)
   real(wp),intent(out) :: W(2702)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV 2702-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.2998675149888161D-4
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.4077860529495355D-3
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2065562538818703D-1
   V=0.1185349192520667D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5250918173022379D-1
   V=0.1913408643425751D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8993480082038376D-1
   V=0.2452886577209897D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1306023924436019D+0
   V=0.2862408183288702D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1732060388531418D+0
   V=0.3178032258257357D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2168727084820249D+0
   V=0.3422945667633690D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2609528309173586D+0
   V=0.3612790520235922D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3049252927938952D+0
   V=0.3758638229818521D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3483484138084404D+0
   V=0.3868711798859953D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3908321549106406D+0
   V=0.3949429933189938D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4320210071894814D+0
   V=0.4006068107541156D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4715824795890053D+0
   V=0.4043192149672723D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5091984794078453D+0
   V=0.4064947495808078D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5445580145650803D+0
   V=0.4075245619813152D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6072575796841768D+0
   V=0.4076423540893566D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6339484505755803D+0
   V=0.4074280862251555D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6570718257486958D+0
   V=0.4074163756012244D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6762557330090709D+0
   V=0.4077647795071246D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6911161696923790D+0
   V=0.4084517552782530D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7012841911659961D+0
   V=0.4092468459224052D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7064559272410020D+0
   V=0.4097872687240906D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6123554989894765D-1
   V=0.1738986811745028D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1533070348312393D+0
   V=0.2659616045280191D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2563902605244206D+0
   V=0.3240596008171533D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3629346991663361D+0
   V=0.3621195964432943D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4683949968987538D+0
   V=0.3868838330760539D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5694479240657952D+0
   V=0.4018911532693111D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6634465430993955D+0
   V=0.4089929432983252D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1033958573552305D+0
   B=0.3034544009063584D-1
   V=0.2279907527706409D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1473521412414395D+0
   B=0.6618803044247135D-1
   V=0.2715205490578897D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1924552158705967D+0
   B=0.1054431128987715D+0
   V=0.3057917896703976D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2381094362890328D+0
   B=0.1468263551238858D+0
   V=0.3326913052452555D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2838121707936760D+0
   B=0.1894486108187886D+0
   V=0.3537334711890037D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3291323133373415D+0
   B=0.2326374238761579D+0
   V=0.3700567500783129D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3736896978741460D+0
   B=0.2758485808485768D+0
   V=0.3825245372589122D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4171406040760013D+0
   B=0.3186179331996921D+0
   V=0.3918125171518296D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4591677985256915D+0
   B=0.3605329796303794D+0
   V=0.3984720419937579D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4994733831718418D+0
   B=0.4012147253586509D+0
   V=0.4029746003338211D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5377731830445096D+0
   B=0.4403050025570692D+0
   V=0.4057428632156627D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5737917830001331D+0
   B=0.4774565904277483D+0
   V=0.4071719274114857D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2027323586271389D+0
   B=0.3544122504976147D-1
   V=0.2990236950664119D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2516942375187273D+0
   B=0.7418304388646328D-1
   V=0.3262951734212878D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3000227995257181D+0
   B=0.1150502745727186D+0
   V=0.3482634608242413D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3474806691046342D+0
   B=0.1571963371209364D+0
   V=0.3656596681700892D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3938103180359209D+0
   B=0.1999631877247100D+0
   V=0.3791740467794218D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4387519590455703D+0
   B=0.2428073457846535D+0
   V=0.3894034450156905D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4820503960077787D+0
   B=0.2852575132906155D+0
   V=0.3968600245508371D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5234573778475101D+0
   B=0.3268884208674639D+0
   V=0.4019931351420050D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5627318647235282D+0
   B=0.3673033321675939D+0
   V=0.4052108801278599D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5996390607156954D+0
   B=0.4061211551830290D+0
   V=0.4068978613940934D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3084780753791947D+0
   B=0.3860125523100059D-1
   V=0.3454275351319704D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3589988275920223D+0
   B=0.7928938987104867D-1
   V=0.3629963537007920D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4078628415881973D+0
   B=0.1212614643030087D+0
   V=0.3770187233889873D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4549287258889735D+0
   B=0.1638770827382693D+0
   V=0.3878608613694378D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5000278512957279D+0
   B=0.2065965798260176D+0
   V=0.3959065270221274D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5429785044928199D+0
   B=0.2489436378852235D+0
   V=0.4015286975463570D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5835939850491711D+0
   B=0.2904811368946891D+0
   V=0.4050866785614717D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6216870353444856D+0
   B=0.3307941957666609D+0
   V=0.4069320185051913D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4151104662709091D+0
   B=0.4064829146052554D-1
   V=0.3760120964062763D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4649804275009218D+0
   B=0.8258424547294755D-1
   V=0.3870969564418064D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5124695757009662D+0
   B=0.1251841962027289D+0
   V=0.3955287790534055D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5574711100606224D+0
   B=0.1679107505976331D+0
   V=0.4015361911302668D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5998597333287227D+0
   B=0.2102805057358715D+0
   V=0.4053836986719548D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6395007148516600D+0
   B=0.2518418087774107D+0
   V=0.4073578673299117D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5188456224746252D+0
   B=0.4194321676077518D-1
   V=0.3954628379231406D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5664190707942778D+0
   B=0.8457661551921499D-1
   V=0.4017645508847530D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6110464353283153D+0
   B=0.1273652932519396D+0
   V=0.4059030348651293D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6526430302051563D+0
   B=0.1698173239076354D+0
   V=0.4080565809484880D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6167551880377548D+0
   B=0.4266398851548864D-1
   V=0.4063018753664651D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6607195418355383D+0
   B=0.8551925814238349D-1
   V=0.4087191292799671D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD3074(X,Y,Z,W,N)
   real(wp),intent(out) :: X(3074)
   real(wp),intent(out) :: Y(3074)
   real(wp),intent(out) :: Z(3074)
   real(wp),intent(out) :: W(3074)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV 3074-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.2599095953754734D-4
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.3603134089687541D-3
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.3586067974412447D-3
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1886108518723392D-1
   V=0.9831528474385880D-4
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4800217244625303D-1
   V=0.1605023107954450D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8244922058397242D-1
   V=0.2072200131464099D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1200408362484023D+0
   V=0.2431297618814187D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1595773530809965D+0
   V=0.2711819064496707D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2002635973434064D+0
   V=0.2932762038321116D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2415127590139982D+0
   V=0.3107032514197368D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2828584158458477D+0
   V=0.3243808058921213D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3239091015338138D+0
   V=0.3349899091374030D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3643225097962194D+0
   V=0.3430580688505218D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4037897083691802D+0
   V=0.3490124109290343D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4420247515194127D+0
   V=0.3532148948561955D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4787572538464938D+0
   V=0.3559862669062833D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5137265251275234D+0
   V=0.3576224317551411D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5466764056654611D+0
   V=0.3584050533086076D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6054859420813535D+0
   V=0.3584903581373224D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6308106701764562D+0
   V=0.3582991879040586D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6530369230179584D+0
   V=0.3582371187963125D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6718609524611158D+0
   V=0.3584353631122350D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6869676499894013D+0
   V=0.3589120166517785D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6980467077240748D+0
   V=0.3595445704531601D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7048241721250522D+0
   V=0.3600943557111074D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5591105222058232D-1
   V=0.1456447096742039D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1407384078513916D+0
   V=0.2252370188283782D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2364035438976309D+0
   V=0.2766135443474897D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3360602737818170D+0
   V=0.3110729491500851D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4356292630054665D+0
   V=0.3342506712303391D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5321569415256174D+0
   V=0.3491981834026860D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6232956305040554D+0
   V=0.3576003604348932D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9469870086838469D-1
   B=0.2778748387309470D-1
   V=0.1921921305788564D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1353170300568141D+0
   B=0.6076569878628364D-1
   V=0.2301458216495632D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1771679481726077D+0
   B=0.9703072762711040D-1
   V=0.2604248549522893D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2197066664231751D+0
   B=0.1354112458524762D+0
   V=0.2845275425870697D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2624783557374927D+0
   B=0.1750996479744100D+0
   V=0.3036870897974840D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3050969521214442D+0
   B=0.2154896907449802D+0
   V=0.3188414832298066D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3472252637196021D+0
   B=0.2560954625740152D+0
   V=0.3307046414722089D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3885610219026360D+0
   B=0.2965070050624096D+0
   V=0.3398330969031360D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4288273776062765D+0
   B=0.3363641488734497D+0
   V=0.3466757899705373D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4677662471302948D+0
   B=0.3753400029836788D+0
   V=0.3516095923230054D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5051333589553359D+0
   B=0.4131297522144286D+0
   V=0.3549645184048486D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5406942145810492D+0
   B=0.4494423776081795D+0
   V=0.3570415969441392D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5742204122576457D+0
   B=0.4839938958841502D+0
   V=0.3581251798496118D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1865407027225188D+0
   B=0.3259144851070796D-1
   V=0.2543491329913348D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2321186453689432D+0
   B=0.6835679505297343D-1
   V=0.2786711051330776D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2773159142523882D+0
   B=0.1062284864451989D+0
   V=0.2985552361083679D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3219200192237254D+0
   B=0.1454404409323047D+0
   V=0.3145867929154039D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3657032593944029D+0
   B=0.1854018282582510D+0
   V=0.3273290662067609D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4084376778363622D+0
   B=0.2256297412014750D+0
   V=0.3372705511943501D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4499004945751427D+0
   B=0.2657104425000896D+0
   V=0.3448274437851510D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4898758141326335D+0
   B=0.3052755487631557D+0
   V=0.3503592783048583D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5281547442266309D+0
   B=0.3439863920645423D+0
   V=0.3541854792663162D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5645346989813992D+0
   B=0.3815229456121914D+0
   V=0.3565995517909428D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5988181252159848D+0
   B=0.4175752420966734D+0
   V=0.3578802078302898D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2850425424471603D+0
   B=0.3562149509862536D-1
   V=0.2958644592860982D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3324619433027876D+0
   B=0.7330318886871096D-1
   V=0.3119548129116835D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3785848333076282D+0
   B=0.1123226296008472D+0
   V=0.3250745225005984D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4232891028562115D+0
   B=0.1521084193337708D+0
   V=0.3355153415935208D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4664287050829722D+0
   B=0.1921844459223610D+0
   V=0.3435847568549328D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5078458493735726D+0
   B=0.2321360989678303D+0
   V=0.3495786831622488D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5473779816204180D+0
   B=0.2715886486360520D+0
   V=0.3537767805534621D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5848617133811376D+0
   B=0.3101924707571355D+0
   V=0.3564459815421428D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6201348281584888D+0
   B=0.3476121052890973D+0
   V=0.3578464061225468D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3852191185387871D+0
   B=0.3763224880035108D-1
   V=0.3239748762836212D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4325025061073423D+0
   B=0.7659581935637135D-1
   V=0.3345491784174287D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4778486229734490D+0
   B=0.1163381306083900D+0
   V=0.3429126177301782D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5211663693009000D+0
   B=0.1563890598752899D+0
   V=0.3492420343097421D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5623469504853703D+0
   B=0.1963320810149200D+0
   V=0.3537399050235257D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6012718188659246D+0
   B=0.2357847407258738D+0
   V=0.3566209152659172D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6378179206390117D+0
   B=0.2743846121244060D+0
   V=0.3581084321919782D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4836936460214534D+0
   B=0.3895902610739024D-1
   V=0.3426522117591512D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5293792562683797D+0
   B=0.7871246819312640D-1
   V=0.3491848770121379D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5726281253100033D+0
   B=0.1187963808202981D+0
   V=0.3539318235231476D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6133658776169068D+0
   B=0.1587914708061787D+0
   V=0.3570231438458694D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6515085491865307D+0
   B=0.1983058575227646D+0
   V=0.3586207335051714D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5778692716064976D+0
   B=0.3977209689791542D-1
   V=0.3541196205164025D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6207904288086192D+0
   B=0.7990157592981152D-1
   V=0.3574296911573953D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6608688171046802D+0
   B=0.1199671308754309D+0
   V=0.3591993279818963D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6656263089489130D+0
   B=0.4015955957805969D-1
   V=0.3595855034661997D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD3470(X,Y,Z,W,N)
   real(wp),intent(out) :: X(3470)
   real(wp),intent(out) :: Y(3470)
   real(wp),intent(out) :: Z(3470)
   real(wp),intent(out) :: W(3470)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV 3470-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.2040382730826330D-4
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.3178149703889544D-3
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1721420832906233D-1
   V=0.8288115128076110D-4
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4408875374981770D-1
   V=0.1360883192522954D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7594680813878681D-1
   V=0.1766854454542662D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1108335359204799D+0
   V=0.2083153161230153D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1476517054388567D+0
   V=0.2333279544657158D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1856731870860615D+0
   V=0.2532809539930247D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2243634099428821D+0
   V=0.2692472184211158D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2633006881662727D+0
   V=0.2819949946811885D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3021340904916283D+0
   V=0.2920953593973030D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3405594048030089D+0
   V=0.2999889782948352D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3783044434007372D+0
   V=0.3060292120496902D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4151194767407910D+0
   V=0.3105109167522192D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4507705766443257D+0
   V=0.3136902387550312D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4850346056573187D+0
   V=0.3157984652454632D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5176950817792470D+0
   V=0.3170516518425422D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5485384240820989D+0
   V=0.3176568425633755D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6039117238943308D+0
   V=0.3177198411207062D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6279956655573113D+0
   V=0.3175519492394733D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6493636169568952D+0
   V=0.3174654952634756D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6677644117704504D+0
   V=0.3175676415467654D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6829368572115624D+0
   V=0.3178923417835410D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6946195818184121D+0
   V=0.3183788287531909D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7025711542057026D+0
   V=0.3188755151918807D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7066004767140119D+0
   V=0.3191916889313849D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5132537689946062D-1
   V=0.1231779611744508D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1297994661331225D+0
   V=0.1924661373839880D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2188852049401307D+0
   V=0.2380881867403424D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3123174824903457D+0
   V=0.2693100663037885D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4064037620738195D+0
   V=0.2908673382834366D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4984958396944782D+0
   V=0.3053914619381535D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5864975046021365D+0
   V=0.3143916684147777D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6686711634580175D+0
   V=0.3187042244055363D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8715738780835950D-1
   B=0.2557175233367578D-1
   V=0.1635219535869790D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1248383123134007D+0
   B=0.5604823383376681D-1
   V=0.1968109917696070D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1638062693383378D+0
   B=0.8968568601900765D-1
   V=0.2236754342249974D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2035586203373176D+0
   B=0.1254086651976279D+0
   V=0.2453186687017181D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2436798975293774D+0
   B=0.1624780150162012D+0
   V=0.2627551791580541D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2838207507773806D+0
   B=0.2003422342683208D+0
   V=0.2767654860152220D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3236787502217692D+0
   B=0.2385628026255263D+0
   V=0.2879467027765895D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3629849554840691D+0
   B=0.2767731148783578D+0
   V=0.2967639918918702D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4014948081992087D+0
   B=0.3146542308245309D+0
   V=0.3035900684660351D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4389818379260225D+0
   B=0.3519196415895088D+0
   V=0.3087338237298308D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4752331143674377D+0
   B=0.3883050984023654D+0
   V=0.3124608838860167D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5100457318374018D+0
   B=0.4235613423908649D+0
   V=0.3150084294226743D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5432238388954868D+0
   B=0.4574484717196220D+0
   V=0.3165958398598402D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5745758685072442D+0
   B=0.4897311639255524D+0
   V=0.3174320440957372D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1723981437592809D+0
   B=0.3010630597881105D-1
   V=0.2182188909812599D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2149553257844597D+0
   B=0.6326031554204694D-1
   V=0.2399727933921445D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2573256081247422D+0
   B=0.9848566980258631D-1
   V=0.2579796133514652D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2993163751238106D+0
   B=0.1350835952384266D+0
   V=0.2727114052623535D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3407238005148000D+0
   B=0.1725184055442181D+0
   V=0.2846327656281355D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3813454978483264D+0
   B=0.2103559279730725D+0
   V=0.2941491102051334D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4209848104423343D+0
   B=0.2482278774554860D+0
   V=0.3016049492136107D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4594519699996300D+0
   B=0.2858099509982883D+0
   V=0.3072949726175648D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4965640166185930D+0
   B=0.3228075659915428D+0
   V=0.3114768142886460D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5321441655571562D+0
   B=0.3589459907204151D+0
   V=0.3143823673666223D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5660208438582166D+0
   B=0.3939630088864310D+0
   V=0.3162269764661535D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5980264315964364D+0
   B=0.4276029922949089D+0
   V=0.3172164663759821D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2644215852350733D+0
   B=0.3300939429072552D-1
   V=0.2554575398967435D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3090113743443063D+0
   B=0.6803887650078501D-1
   V=0.2701704069135677D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3525871079197808D+0
   B=0.1044326136206709D+0
   V=0.2823693413468940D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3950418005354029D+0
   B=0.1416751597517679D+0
   V=0.2922898463214289D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4362475663430163D+0
   B=0.1793408610504821D+0
   V=0.3001829062162428D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4760661812145854D+0
   B=0.2170630750175722D+0
   V=0.3062890864542953D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5143551042512103D+0
   B=0.2545145157815807D+0
   V=0.3108328279264746D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5509709026935597D+0
   B=0.2913940101706601D+0
   V=0.3140243146201245D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5857711030329428D+0
   B=0.3274169910910705D+0
   V=0.3160638030977130D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6186149917404392D+0
   B=0.3623081329317265D+0
   V=0.3171462882206275D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3586894569557064D+0
   B=0.3497354386450040D-1
   V=0.2812388416031796D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4035266610019441D+0
   B=0.7129736739757095D-1
   V=0.2912137500288045D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4467775312332510D+0
   B=0.1084758620193165D+0
   V=0.2993241256502206D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4883638346608543D+0
   B=0.1460915689241772D+0
   V=0.3057101738983822D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5281908348434601D+0
   B=0.1837790832369980D+0
   V=0.3105319326251432D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5661542687149311D+0
   B=0.2212075390874021D+0
   V=0.3139565514428167D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6021450102031452D+0
   B=0.2580682841160985D+0
   V=0.3161543006806366D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6360520783610050D+0
   B=0.2940656362094121D+0
   V=0.3172985960613294D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4521611065087196D+0
   B=0.3631055365867002D-1
   V=0.2989400336901431D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4959365651560963D+0
   B=0.7348318468484350D-1
   V=0.3054555883947677D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5376815804038283D+0
   B=0.1111087643812648D+0
   V=0.3104764960807702D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5773314480243768D+0
   B=0.1488226085145408D+0
   V=0.3141015825977616D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6148113245575056D+0
   B=0.1862892274135151D+0
   V=0.3164520621159896D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6500407462842380D+0
   B=0.2231909701714456D+0
   V=0.3176652305912204D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5425151448707213D+0
   B=0.3718201306118944D-1
   V=0.3105097161023939D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5841860556907931D+0
   B=0.7483616335067346D-1
   V=0.3143014117890550D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6234632186851500D+0
   B=0.1125990834266120D+0
   V=0.3168172866287200D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6602934551848843D+0
   B=0.1501303813157619D+0
   V=0.3181401865570968D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6278573968375105D+0
   B=0.3767559930245720D-1
   V=0.3170663659156037D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6665611711264577D+0
   B=0.7548443301360158D-1
   V=0.3185447944625510D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD3890(X,Y,Z,W,N)
   real(wp),intent(out) :: X(3890)
   real(wp),intent(out) :: Y(3890)
   real(wp),intent(out) :: Z(3890)
   real(wp),intent(out) :: W(3890)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV 3890-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.1807395252196920D-4
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.2848008782238827D-3
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.2836065837530581D-3
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1587876419858352D-1
   V=0.7013149266673816D-4
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4069193593751206D-1
   V=0.1162798021956766D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7025888115257997D-1
   V=0.1518728583972105D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1027495450028704D+0
   V=0.1798796108216934D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1371457730893426D+0
   V=0.2022593385972785D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1727758532671953D+0
   V=0.2203093105575464D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2091492038929037D+0
   V=0.2349294234299855D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2458813281751915D+0
   V=0.2467682058747003D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2826545859450066D+0
   V=0.2563092683572224D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3191957291799622D+0
   V=0.2639253896763318D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3552621469299578D+0
   V=0.2699137479265108D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3906329503406230D+0
   V=0.2745196420166739D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4251028614093031D+0
   V=0.2779529197397593D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4584777520111870D+0
   V=0.2803996086684265D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4905711358710193D+0
   V=0.2820302356715842D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5212011669847385D+0
   V=0.2830056747491068D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5501878488737995D+0
   V=0.2834808950776839D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6025037877479342D+0
   V=0.2835282339078929D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6254572689549016D+0
   V=0.2833819267065800D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6460107179528248D+0
   V=0.2832858336906784D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6639541138154251D+0
   V=0.2833268235451244D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6790688515667495D+0
   V=0.2835432677029253D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6911338580371512D+0
   V=0.2839091722743049D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6999385956126490D+0
   V=0.2843308178875841D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7053037748656896D+0
   V=0.2846703550533846D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4732224387180115D-1
   V=0.1051193406971900D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1202100529326803D+0
   V=0.1657871838796974D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2034304820664855D+0
   V=0.2064648113714232D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2912285643573002D+0
   V=0.2347942745819741D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3802361792726768D+0
   V=0.2547775326597726D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4680598511056146D+0
   V=0.2686876684847025D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5528151052155599D+0
   V=0.2778665755515867D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6329386307803041D+0
   V=0.2830996616782929D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8056516651369069D-1
   B=0.2363454684003124D-1
   V=0.1403063340168372D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1156476077139389D+0
   B=0.5191291632545936D-1
   V=0.1696504125939477D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1520473382760421D+0
   B=0.8322715736994519D-1
   V=0.1935787242745390D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1892986699745931D+0
   B=0.1165855667993712D+0
   V=0.2130614510521968D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2270194446777792D+0
   B=0.1513077167409504D+0
   V=0.2289381265931048D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2648908185093273D+0
   B=0.1868882025807859D+0
   V=0.2418630292816186D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3026389259574136D+0
   B=0.2229277629776224D+0
   V=0.2523400495631193D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3400220296151384D+0
   B=0.2590951840746235D+0
   V=0.2607623973449605D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3768217953335510D+0
   B=0.2951047291750847D+0
   V=0.2674441032689209D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4128372900921884D+0
   B=0.3307019714169930D+0
   V=0.2726432360343356D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4478807131815630D+0
   B=0.3656544101087634D+0
   V=0.2765787685924545D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4817742034089257D+0
   B=0.3997448951939695D+0
   V=0.2794428690642224D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5143472814653344D+0
   B=0.4327667110812024D+0
   V=0.2814099002062895D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5454346213905650D+0
   B=0.4645196123532293D+0
   V=0.2826429531578994D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5748739313170252D+0
   B=0.4948063555703345D+0
   V=0.2832983542550884D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1599598738286342D+0
   B=0.2792357590048985D-1
   V=0.1886695565284976D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1998097412500951D+0
   B=0.5877141038139065D-1
   V=0.2081867882748234D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2396228952566202D+0
   B=0.9164573914691377D-1
   V=0.2245148680600796D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2792228341097746D+0
   B=0.1259049641962687D+0
   V=0.2380370491511872D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3184251107546741D+0
   B=0.1610594823400863D+0
   V=0.2491398041852455D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3570481164426244D+0
   B=0.1967151653460898D+0
   V=0.2581632405881230D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3949164710492144D+0
   B=0.2325404606175168D+0
   V=0.2653965506227417D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4318617293970503D+0
   B=0.2682461141151439D+0
   V=0.2710857216747087D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4677221009931678D+0
   B=0.3035720116011973D+0
   V=0.2754434093903659D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5023417939270955D+0
   B=0.3382781859197439D+0
   V=0.2786579932519380D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5355701836636128D+0
   B=0.3721383065625942D+0
   V=0.2809011080679474D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5672608451328771D+0
   B=0.4049346360466055D+0
   V=0.2823336184560987D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5972704202540162D+0
   B=0.4364538098633802D+0
   V=0.2831101175806309D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2461687022333596D+0
   B=0.3070423166833368D-1
   V=0.2221679970354546D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2881774566286831D+0
   B=0.6338034669281885D-1
   V=0.2356185734270703D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3293963604116978D+0
   B=0.9742862487067941D-1
   V=0.2469228344805590D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3697303822241377D+0
   B=0.1323799532282290D+0
   V=0.2562726348642046D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4090663023135127D+0
   B=0.1678497018129336D+0
   V=0.2638756726753028D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4472819355411712D+0
   B=0.2035095105326114D+0
   V=0.2699311157390862D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4842513377231437D+0
   B=0.2390692566672091D+0
   V=0.2746233268403837D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5198477629962928D+0
   B=0.2742649818076149D+0
   V=0.2781225674454771D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5539453011883145D+0
   B=0.3088503806580094D+0
   V=0.2805881254045684D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5864196762401251D+0
   B=0.3425904245906614D+0
   V=0.2821719877004913D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6171484466668390D+0
   B=0.3752562294789468D+0
   V=0.2830222502333124D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3350337830565727D+0
   B=0.3261589934634747D-1
   V=0.2457995956744870D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3775773224758284D+0
   B=0.6658438928081572D-1
   V=0.2551474407503706D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4188155229848973D+0
   B=0.1014565797157954D+0
   V=0.2629065335195311D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4586805892009344D+0
   B=0.1368573320843822D+0
   V=0.2691900449925075D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4970895714224235D+0
   B=0.1724614851951608D+0
   V=0.2741275485754276D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5339505133960747D+0
   B=0.2079779381416412D+0
   V=0.2778530970122595D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5691665792531440D+0
   B=0.2431385788322288D+0
   V=0.2805010567646741D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6026387682680377D+0
   B=0.2776901883049853D+0
   V=0.2822055834031040D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6342676150163307D+0
   B=0.3113881356386632D+0
   V=0.2831016901243473D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4237951119537067D+0
   B=0.3394877848664351D-1
   V=0.2624474901131803D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4656918683234929D+0
   B=0.6880219556291447D-1
   V=0.2688034163039377D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5058857069185980D+0
   B=0.1041946859721635D+0
   V=0.2738932751287636D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5443204666713996D+0
   B=0.1398039738736393D+0
   V=0.2777944791242523D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5809298813759742D+0
   B=0.1753373381196155D+0
   V=0.2806011661660987D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6156416039447128D+0
   B=0.2105215793514010D+0
   V=0.2824181456597460D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6483801351066604D+0
   B=0.2450953312157051D+0
   V=0.2833585216577828D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5103616577251688D+0
   B=0.3485560643800719D-1
   V=0.2738165236962878D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5506738792580681D+0
   B=0.7026308631512033D-1
   V=0.2778365208203180D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5889573040995292D+0
   B=0.1059035061296403D+0
   V=0.2807852940418966D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6251641589516930D+0
   B=0.1414823925236026D+0
   V=0.2827245949674705D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6592414921570178D+0
   B=0.1767207908214530D+0
   V=0.2837342344829828D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5930314017533384D+0
   B=0.3542189339561672D-1
   V=0.2809233907610981D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6309812253390175D+0
   B=0.7109574040369549D-1
   V=0.2829930809742694D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6666296011353230D+0
   B=0.1067259792282730D+0
   V=0.2841097874111479D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6703715271049922D+0
   B=0.3569455268820809D-1
   V=0.2843455206008783D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD4334(X,Y,Z,W,N)
   real(wp),intent(out) :: X(4334)
   real(wp),intent(out) :: Y(4334)
   real(wp),intent(out) :: Z(4334)
   real(wp),intent(out) :: W(4334)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV 4334-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.1449063022537883D-4
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.2546377329828424D-3
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1462896151831013D-1
   V=0.6018432961087496D-4
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3769840812493139D-1
   V=0.1002286583263673D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6524701904096891D-1
   V=0.1315222931028093D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9560543416134648D-1
   V=0.1564213746876724D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1278335898929198D+0
   V=0.1765118841507736D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1613096104466031D+0
   V=0.1928737099311080D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1955806225745371D+0
   V=0.2062658534263270D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2302935218498028D+0
   V=0.2172395445953787D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2651584344113027D+0
   V=0.2262076188876047D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2999276825183209D+0
   V=0.2334885699462397D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3343828669718798D+0
   V=0.2393355273179203D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3683265013750518D+0
   V=0.2439559200468863D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4015763206518108D+0
   V=0.2475251866060002D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4339612026399770D+0
   V=0.2501965558158773D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4653180651114582D+0
   V=0.2521081407925925D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4954893331080803D+0
   V=0.2533881002388081D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5243207068924930D+0
   V=0.2541582900848261D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5516590479041704D+0
   V=0.2545365737525860D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6012371927804176D+0
   V=0.2545726993066799D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6231574466449819D+0
   V=0.2544456197465555D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6429416514181271D+0
   V=0.2543481596881064D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6604124272943595D+0
   V=0.2543506451429194D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6753851470408250D+0
   V=0.2544905675493763D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6876717970626160D+0
   V=0.2547611407344429D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6970895061319234D+0
   V=0.2551060375448869D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7034746912553310D+0
   V=0.2554291933816039D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7067017217542295D+0
   V=0.2556255710686343D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4382223501131123D-1
   V=0.9041339695118195D-4
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1117474077400006D+0
   V=0.1438426330079022D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1897153252911440D+0
   V=0.1802523089820518D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2724023009910331D+0
   V=0.2060052290565496D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3567163308709902D+0
   V=0.2245002248967466D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4404784483028087D+0
   V=0.2377059847731150D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5219833154161411D+0
   V=0.2468118955882525D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5998179868977553D+0
   V=0.2525410872966528D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6727803154548222D+0
   V=0.2553101409933397D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7476563943166086D-1
   B=0.2193168509461185D-1
   V=0.1212879733668632D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1075341482001416D+0
   B=0.4826419281533887D-1
   V=0.1472872881270931D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1416344885203259D+0
   B=0.7751191883575742D-1
   V=0.1686846601010828D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1766325315388586D+0
   B=0.1087558139247680D+0
   V=0.1862698414660208D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2121744174481514D+0
   B=0.1413661374253096D+0
   V=0.2007430956991861D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2479669443408145D+0
   B=0.1748768214258880D+0
   V=0.2126568125394796D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2837600452294113D+0
   B=0.2089216406612073D+0
   V=0.2224394603372113D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3193344933193984D+0
   B=0.2431987685545972D+0
   V=0.2304264522673135D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3544935442438745D+0
   B=0.2774497054377770D+0
   V=0.2368854288424087D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3890571932288154D+0
   B=0.3114460356156915D+0
   V=0.2420352089461772D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4228581214259090D+0
   B=0.3449806851913012D+0
   V=0.2460597113081295D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4557387211304052D+0
   B=0.3778618641248256D+0
   V=0.2491181912257687D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4875487950541643D+0
   B=0.4099086391698978D+0
   V=0.2513528194205857D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5181436529962997D+0
   B=0.4409474925853973D+0
   V=0.2528943096693220D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5473824095600661D+0
   B=0.4708094517711291D+0
   V=0.2538660368488136D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5751263398976174D+0
   B=0.4993275140354637D+0
   V=0.2543868648299022D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1489515746840028D+0
   B=0.2599381993267017D-1
   V=0.1642595537825183D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1863656444351767D+0
   B=0.5479286532462190D-1
   V=0.1818246659849308D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2238602880356348D+0
   B=0.8556763251425254D-1
   V=0.1966565649492420D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2612723375728160D+0
   B=0.1177257802267011D+0
   V=0.2090677905657991D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2984332990206190D+0
   B=0.1508168456192700D+0
   V=0.2193820409510504D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3351786584663333D+0
   B=0.1844801892177727D+0
   V=0.2278870827661928D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3713505522209120D+0
   B=0.2184145236087598D+0
   V=0.2348283192282090D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4067981098954663D+0
   B=0.2523590641486229D+0
   V=0.2404139755581477D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4413769993687534D+0
   B=0.2860812976901373D+0
   V=0.2448227407760734D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4749487182516394D+0
   B=0.3193686757808996D+0
   V=0.2482110455592573D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5073798105075426D+0
   B=0.3520226949547602D+0
   V=0.2507192397774103D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5385410448878654D+0
   B=0.3838544395667890D+0
   V=0.2524765968534880D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5683065353670530D+0
   B=0.4146810037640963D+0
   V=0.2536052388539425D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5965527620663510D+0
   B=0.4443224094681121D+0
   V=0.2542230588033068D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2299227700856157D+0
   B=0.2865757664057584D-1
   V=0.1944817013047896D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2695752998553267D+0
   B=0.5923421684485993D-1
   V=0.2067862362746635D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3086178716611389D+0
   B=0.9117817776057715D-1
   V=0.2172440734649114D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3469649871659077D+0
   B=0.1240593814082605D+0
   V=0.2260125991723423D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3845153566319655D+0
   B=0.1575272058259175D+0
   V=0.2332655008689523D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4211600033403215D+0
   B=0.1912845163525413D+0
   V=0.2391699681532458D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4567867834329882D+0
   B=0.2250710177858171D+0
   V=0.2438801528273928D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4912829319232061D+0
   B=0.2586521303440910D+0
   V=0.2475370504260665D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5245364793303812D+0
   B=0.2918112242865407D+0
   V=0.2502707235640574D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5564369788915756D+0
   B=0.3243439239067890D+0
   V=0.2522031701054241D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5868757697775287D+0
   B=0.3560536787835351D+0
   V=0.2534511269978784D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6157458853519617D+0
   B=0.3867480821242581D+0
   V=0.2541284914955151D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3138461110672113D+0
   B=0.3051374637507278D-1
   V=0.2161509250688394D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3542495872050569D+0
   B=0.6237111233730755D-1
   V=0.2248778513437852D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3935751553120181D+0
   B=0.9516223952401907D-1
   V=0.2322388803404617D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4317634668111147D+0
   B=0.1285467341508517D+0
   V=0.2383265471001355D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4687413842250821D+0
   B=0.1622318931656033D+0
   V=0.2432476675019525D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5044274237060283D+0
   B=0.1959581153836453D+0
   V=0.2471122223750674D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5387354077925727D+0
   B=0.2294888081183837D+0
   V=0.2500291752486870D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5715768898356105D+0
   B=0.2626031152713945D+0
   V=0.2521055942764682D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6028627200136111D+0
   B=0.2950904075286713D+0
   V=0.2534472785575503D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6325039812653463D+0
   B=0.3267458451113286D+0
   V=0.2541599713080121D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3981986708423407D+0
   B=0.3183291458749821D-1
   V=0.2317380975862936D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4382791182133300D+0
   B=0.6459548193880908D-1
   V=0.2378550733719775D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4769233057218166D+0
   B=0.9795757037087952D-1
   V=0.2428884456739118D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5140823911194238D+0
   B=0.1316307235126655D+0
   V=0.2469002655757292D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5496977833862983D+0
   B=0.1653556486358704D+0
   V=0.2499657574265851D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5837047306512727D+0
   B=0.1988931724126510D+0
   V=0.2521676168486082D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6160349566926879D+0
   B=0.2320174581438950D+0
   V=0.2535935662645334D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6466185353209440D+0
   B=0.2645106562168662D+0
   V=0.2543356743363214D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4810835158795404D+0
   B=0.3275917807743992D-1
   V=0.2427353285201535D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5199925041324341D+0
   B=0.6612546183967181D-1
   V=0.2468258039744386D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5571717692207494D+0
   B=0.9981498331474143D-1
   V=0.2500060956440310D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5925789250836378D+0
   B=0.1335687001410374D+0
   V=0.2523238365420979D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6261658523859670D+0
   B=0.1671444402896463D+0
   V=0.2538399260252846D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6578811126669331D+0
   B=0.2003106382156076D+0
   V=0.2546255927268069D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5609624612998100D+0
   B=0.3337500940231335D-1
   V=0.2500583360048449D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5979959659984670D+0
   B=0.6708750335901803D-1
   V=0.2524777638260203D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6330523711054002D+0
   B=0.1008792126424850D+0
   V=0.2540951193860656D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6660960998103972D+0
   B=0.1345050343171794D+0
   V=0.2549524085027472D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6365384364585819D+0
   B=0.3372799460737052D-1
   V=0.2542569507009158D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6710994302899275D+0
   B=0.6755249309678028D-1
   V=0.2552114127580376D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD4802(X,Y,Z,W,N)
   real(wp),intent(out) :: X(4802)
   real(wp),intent(out) :: Y(4802)
   real(wp),intent(out) :: Z(4802)
   real(wp),intent(out) :: W(4802)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV 4802-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.9687521879420705D-4
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.2307897895367918D-3
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.2297310852498558D-3
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2335728608887064D-1
   V=0.7386265944001919D-4
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4352987836550653D-1
   V=0.8257977698542210D-4
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6439200521088801D-1
   V=0.9706044762057630D-4
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9003943631993181D-1
   V=0.1302393847117003D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1196706615548473D+0
   V=0.1541957004600968D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1511715412838134D+0
   V=0.1704459770092199D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1835982828503801D+0
   V=0.1827374890942906D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2165081259155405D+0
   V=0.1926360817436107D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2496208720417563D+0
   V=0.2008010239494833D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2827200673567900D+0
   V=0.2075635983209175D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3156190823994346D+0
   V=0.2131306638690909D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3481476793749115D+0
   V=0.2176562329937335D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3801466086947226D+0
   V=0.2212682262991018D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4114652119634011D+0
   V=0.2240799515668565D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4419598786519751D+0
   V=0.2261959816187525D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4714925949329543D+0
   V=0.2277156368808855D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4999293972879466D+0
   V=0.2287351772128336D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5271387221431248D+0
   V=0.2293490814084085D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5529896780837761D+0
   V=0.2296505312376273D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6000856099481712D+0
   V=0.2296793832318756D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6210562192785175D+0
   V=0.2295785443842974D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6401165879934240D+0
   V=0.2295017931529102D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6571144029244334D+0
   V=0.2295059638184868D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6718910821718863D+0
   V=0.2296232343237362D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6842845591099010D+0
   V=0.2298530178740771D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6941353476269816D+0
   V=0.2301579790280501D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7012965242212991D+0
   V=0.2304690404996513D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7056471428242644D+0
   V=0.2307027995907102D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4595557643585895D-1
   V=0.9312274696671092D-4
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1049316742435023D+0
   V=0.1199919385876926D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1773548879549274D+0
   V=0.1598039138877690D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2559071411236127D+0
   V=0.1822253763574900D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3358156837985898D+0
   V=0.1988579593655040D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4155835743763893D+0
   V=0.2112620102533307D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4937894296167472D+0
   V=0.2201594887699007D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5691569694793316D+0
   V=0.2261622590895036D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6405840854894251D+0
   V=0.2296458453435705D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7345133894143348D-1
   B=0.2177844081486067D-1
   V=0.1006006990267000D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1009859834044931D+0
   B=0.4590362185775188D-1
   V=0.1227676689635876D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1324289619748758D+0
   B=0.7255063095690877D-1
   V=0.1467864280270117D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1654272109607127D+0
   B=0.1017825451960684D+0
   V=0.1644178912101232D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1990767186776461D+0
   B=0.1325652320980364D+0
   V=0.1777664890718961D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2330125945523278D+0
   B=0.1642765374496765D+0
   V=0.1884825664516690D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2670080611108287D+0
   B=0.1965360374337889D+0
   V=0.1973269246453848D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3008753376294316D+0
   B=0.2290726770542238D+0
   V=0.2046767775855328D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3344475596167860D+0
   B=0.2616645495370823D+0
   V=0.2107600125918040D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3675709724070786D+0
   B=0.2941150728843141D+0
   V=0.2157416362266829D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4001000887587812D+0
   B=0.3262440400919066D+0
   V=0.2197557816920721D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4318956350436028D+0
   B=0.3578835350611916D+0
   V=0.2229192611835437D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4628239056795531D+0
   B=0.3888751854043678D+0
   V=0.2253385110212775D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4927563229773636D+0
   B=0.4190678003222840D+0
   V=0.2271137107548774D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5215687136707969D+0
   B=0.4483151836883852D+0
   V=0.2283414092917525D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5491402346984905D+0
   B=0.4764740676087880D+0
   V=0.2291161673130077D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5753520160126075D+0
   B=0.5034021310998277D+0
   V=0.2295313908576598D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1388326356417754D+0
   B=0.2435436510372806D-1
   V=0.1438204721359031D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1743686900537244D+0
   B=0.5118897057342652D-1
   V=0.1607738025495257D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2099737037950268D+0
   B=0.8014695048539634D-1
   V=0.1741483853528379D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2454492590908548D+0
   B=0.1105117874155699D+0
   V=0.1851918467519151D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2807219257864278D+0
   B=0.1417950531570966D+0
   V=0.1944628638070613D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3156842271975842D+0
   B=0.1736604945719597D+0
   V=0.2022495446275152D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3502090945177752D+0
   B=0.2058466324693981D+0
   V=0.2087462382438514D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3841684849519686D+0
   B=0.2381284261195919D+0
   V=0.2141074754818308D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4174372367906016D+0
   B=0.2703031270422569D+0
   V=0.2184640913748162D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4498926465011892D+0
   B=0.3021845683091309D+0
   V=0.2219309165220329D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4814146229807701D+0
   B=0.3335993355165720D+0
   V=0.2246123118340624D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5118863625734701D+0
   B=0.3643833735518232D+0
   V=0.2266062766915125D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5411947455119144D+0
   B=0.3943789541958179D+0
   V=0.2280072952230796D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5692301500357246D+0
   B=0.4234320144403542D+0
   V=0.2289082025202583D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5958857204139576D+0
   B=0.4513897947419260D+0
   V=0.2294012695120025D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2156270284785766D+0
   B=0.2681225755444491D-1
   V=0.1722434488736947D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2532385054909710D+0
   B=0.5557495747805614D-1
   V=0.1830237421455091D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2902564617771537D+0
   B=0.8569368062950249D-1
   V=0.1923855349997633D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3266979823143256D+0
   B=0.1167367450324135D+0
   V=0.2004067861936271D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3625039627493614D+0
   B=0.1483861994003304D+0
   V=0.2071817297354263D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3975838937548699D+0
   B=0.1803821503011405D+0
   V=0.2128250834102103D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4318396099009774D+0
   B=0.2124962965666424D+0
   V=0.2174513719440102D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4651706555732742D+0
   B=0.2445221837805913D+0
   V=0.2211661839150214D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4974752649620969D+0
   B=0.2762701224322987D+0
   V=0.2240665257813102D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5286517579627517D+0
   B=0.3075627775211328D+0
   V=0.2262439516632620D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5586001195731895D+0
   B=0.3382311089826877D+0
   V=0.2277874557231869D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5872229902021319D+0
   B=0.3681108834741399D+0
   V=0.2287854314454994D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6144258616235123D+0
   B=0.3970397446872839D+0
   V=0.2293268499615575D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2951676508064861D+0
   B=0.2867499538750441D-1
   V=0.1912628201529828D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3335085485472725D+0
   B=0.5867879341903510D-1
   V=0.1992499672238701D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3709561760636381D+0
   B=0.8961099205022284D-1
   V=0.2061275533454027D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4074722861667498D+0
   B=0.1211627927626297D+0
   V=0.2119318215968572D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4429923648839117D+0
   B=0.1530748903554898D+0
   V=0.2167416581882652D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4774428052721736D+0
   B=0.1851176436721877D+0
   V=0.2206430730516600D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5107446539535904D+0
   B=0.2170829107658179D+0
   V=0.2237186938699523D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5428151370542935D+0
   B=0.2487786689026271D+0
   V=0.2260480075032884D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5735699292556964D+0
   B=0.2800239952795016D+0
   V=0.2277098884558542D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6029253794562866D+0
   B=0.3106445702878119D+0
   V=0.2287845715109671D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6307998987073145D+0
   B=0.3404689500841194D+0
   V=0.2293547268236294D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3752652273692719D+0
   B=0.2997145098184479D-1
   V=0.2056073839852528D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4135383879344028D+0
   B=0.6086725898678011D-1
   V=0.2114235865831876D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4506113885153907D+0
   B=0.9238849548435643D-1
   V=0.2163175629770551D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4864401554606072D+0
   B=0.1242786603851851D+0
   V=0.2203392158111650D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5209708076611709D+0
   B=0.1563086731483386D+0
   V=0.2235473176847839D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5541422135830122D+0
   B=0.1882696509388506D+0
   V=0.2260024141501235D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5858880915113817D+0
   B=0.2199672979126059D+0
   V=0.2277675929329182D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6161399390603444D+0
   B=0.2512165482924867D+0
   V=0.2289102112284834D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6448296482255090D+0
   B=0.2818368701871888D+0
   V=0.2295027954625118D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4544796274917948D+0
   B=0.3088970405060312D-1
   V=0.2161281589879992D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4919389072146628D+0
   B=0.6240947677636835D-1
   V=0.2201980477395102D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5279313026985183D+0
   B=0.9430706144280313D-1
   V=0.2234952066593166D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5624169925571135D+0
   B=0.1263547818770374D+0
   V=0.2260540098520838D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5953484627093287D+0
   B=0.1583430788822594D+0
   V=0.2279157981899988D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6266730715339185D+0
   B=0.1900748462555988D+0
   V=0.2291296918565571D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6563363204278871D+0
   B=0.2213599519592567D+0
   V=0.2297533752536649D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5314574716585696D+0
   B=0.3152508811515374D-1
   V=0.2234927356465995D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5674614932298185D+0
   B=0.6343865291465561D-1
   V=0.2261288012985219D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6017706004970264D+0
   B=0.9551503504223951D-1
   V=0.2280818160923688D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6343471270264178D+0
   B=0.1275440099801196D+0
   V=0.2293773295180159D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6651494599127802D+0
   B=0.1593252037671960D+0
   V=0.2300528767338634D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6050184986005704D+0
   B=0.3192538338496105D-1
   V=0.2281893855065666D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6390163550880400D+0
   B=0.6402824353962306D-1
   V=0.2295720444840727D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6711199107088448D+0
   B=0.9609805077002909D-1
   V=0.2303227649026753D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6741354429572275D+0
   B=0.3211853196273233D-1
   V=0.2304831913227114D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD5294(X,Y,Z,W,N)
   real(wp),intent(out) :: X(5294)
   real(wp),intent(out) :: Y(5294)
   real(wp),intent(out) :: Z(5294)
   real(wp),intent(out) :: W(5294)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV 5294-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.9080510764308163D-4
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.2084824361987793D-3
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2303261686261450D-1
   V=0.5011105657239616D-4
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3757208620162394D-1
   V=0.5942520409683854D-4
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5821912033821852D-1
   V=0.9564394826109721D-4
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8403127529194872D-1
   V=0.1185530657126338D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1122927798060578D+0
   V=0.1364510114230331D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1420125319192987D+0
   V=0.1505828825605415D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1726396437341978D+0
   V=0.1619298749867023D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2038170058115696D+0
   V=0.1712450504267789D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2352849892876508D+0
   V=0.1789891098164999D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2668363354312461D+0
   V=0.1854474955629795D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2982941279900452D+0
   V=0.1908148636673661D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3295002922087076D+0
   V=0.1952377405281833D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3603094918363593D+0
   V=0.1988349254282232D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3905857895173920D+0
   V=0.2017079807160050D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4202005758160837D+0
   V=0.2039473082709094D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4490310061597227D+0
   V=0.2056360279288953D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4769586160311491D+0
   V=0.2068525823066865D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5038679887049750D+0
   V=0.2076724877534488D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5296454286519961D+0
   V=0.2081694278237885D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5541776207164850D+0
   V=0.2084157631219326D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5990467321921213D+0
   V=0.2084381531128593D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6191467096294587D+0
   V=0.2083476277129307D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6375251212901849D+0
   V=0.2082686194459732D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6540514381131168D+0
   V=0.2082475686112415D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6685899064391510D+0
   V=0.2083139860289915D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6810013009681648D+0
   V=0.2084745561831237D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6911469578730340D+0
   V=0.2087091313375890D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6988956915141736D+0
   V=0.2089718413297697D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7041335794868720D+0
   V=0.2092003303479793D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7067754398018567D+0
   V=0.2093336148263241D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3840368707853623D-1
   V=0.7591708117365267D-4
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9835485954117399D-1
   V=0.1083383968169186D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1665774947612998D+0
   V=0.1403019395292510D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2405702335362910D+0
   V=0.1615970179286436D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3165270770189046D+0
   V=0.1771144187504911D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3927386145645443D+0
   V=0.1887760022988168D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4678825918374656D+0
   V=0.1973474670768214D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5408022024266935D+0
   V=0.2033787661234659D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6104967445752438D+0
   V=0.2072343626517331D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6760910702685738D+0
   V=0.2091177834226918D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6655644120217392D-1
   B=0.1936508874588424D-1
   V=0.9316684484675566D-4
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9446246161270182D-1
   B=0.4252442002115869D-1
   V=0.1116193688682976D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1242651925452509D+0
   B=0.6806529315354374D-1
   V=0.1298623551559414D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1553438064846751D+0
   B=0.9560957491205369D-1
   V=0.1450236832456426D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1871137110542670D+0
   B=0.1245931657452888D+0
   V=0.1572719958149914D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2192612628836257D+0
   B=0.1545385828778978D+0
   V=0.1673234785867195D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2515682807206955D+0
   B=0.1851004249723368D+0
   V=0.1756860118725188D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2838535866287290D+0
   B=0.2160182608272384D+0
   V=0.1826776290439367D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3159578817528521D+0
   B=0.2470799012277111D+0
   V=0.1885116347992865D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3477370882791392D+0
   B=0.2781014208986402D+0
   V=0.1933457860170574D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3790576960890540D+0
   B=0.3089172523515731D+0
   V=0.1973060671902064D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4097938317810200D+0
   B=0.3393750055472244D+0
   V=0.2004987099616311D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4398256572859637D+0
   B=0.3693322470987730D+0
   V=0.2030170909281499D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4690384114718480D+0
   B=0.3986541005609877D+0
   V=0.2049461460119080D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4973216048301053D+0
   B=0.4272112491408562D+0
   V=0.2063653565200186D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5245681526132446D+0
   B=0.4548781735309936D+0
   V=0.2073507927381027D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5506733911803888D+0
   B=0.4815315355023251D+0
   V=0.2079764593256122D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5755339829522475D+0
   B=0.5070486445801855D+0
   V=0.2083150534968778D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1305472386056362D+0
   B=0.2284970375722366D-1
   V=0.1262715121590664D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1637327908216477D+0
   B=0.4812254338288384D-1
   V=0.1414386128545972D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1972734634149637D+0
   B=0.7531734457511935D-1
   V=0.1538740401313898D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2308694653110130D+0
   B=0.1039043639882017D+0
   V=0.1642434942331432D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2643899218338160D+0
   B=0.1334526587117626D+0
   V=0.1729790609237496D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2977171599622171D+0
   B=0.1636414868936382D+0
   V=0.1803505190260828D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3307293903032310D+0
   B=0.1942195406166568D+0
   V=0.1865475350079657D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3633069198219073D+0
   B=0.2249752879943753D+0
   V=0.1917182669679069D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3953346955922727D+0
   B=0.2557218821820032D+0
   V=0.1959851709034382D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4267018394184914D+0
   B=0.2862897925213193D+0
   V=0.1994529548117882D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4573009622571704D+0
   B=0.3165224536636518D+0
   V=0.2022138911146548D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4870279559856109D+0
   B=0.3462730221636496D+0
   V=0.2043518024208592D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5157819581450322D+0
   B=0.3754016870282835D+0
   V=0.2059450313018110D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5434651666465393D+0
   B=0.4037733784993613D+0
   V=0.2070685715318472D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5699823887764627D+0
   B=0.4312557784139123D+0
   V=0.2077955310694373D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5952403350947741D+0
   B=0.4577175367122110D+0
   V=0.2081980387824712D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2025152599210369D+0
   B=0.2520253617719557D-1
   V=0.1521318610377956D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2381066653274425D+0
   B=0.5223254506119000D-1
   V=0.1622772720185755D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2732823383651612D+0
   B=0.8060669688588620D-1
   V=0.1710498139420709D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3080137692611118D+0
   B=0.1099335754081255D+0
   V=0.1785911149448736D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3422405614587601D+0
   B=0.1399120955959857D+0
   V=0.1850125313687736D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3758808773890420D+0
   B=0.1702977801651705D+0
   V=0.1904229703933298D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4088458383438932D+0
   B=0.2008799256601680D+0
   V=0.1949259956121987D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4410450550841152D+0
   B=0.2314703052180836D+0
   V=0.1986161545363960D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4723879420561312D+0
   B=0.2618972111375892D+0
   V=0.2015790585641370D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5027843561874343D+0
   B=0.2920013195600270D+0
   V=0.2038934198707418D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5321453674452458D+0
   B=0.3216322555190551D+0
   V=0.2056334060538251D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5603839113834030D+0
   B=0.3506456615934198D+0
   V=0.2068705959462289D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5874150706875146D+0
   B=0.3789007181306267D+0
   V=0.2076753906106002D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6131559381660038D+0
   B=0.4062580170572782D+0
   V=0.2081179391734803D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2778497016394506D+0
   B=0.2696271276876226D-1
   V=0.1700345216228943D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3143733562261912D+0
   B=0.5523469316960465D-1
   V=0.1774906779990410D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3501485810261827D+0
   B=0.8445193201626464D-1
   V=0.1839659377002642D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3851430322303653D+0
   B=0.1143263119336083D+0
   V=0.1894987462975169D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4193013979470415D+0
   B=0.1446177898344475D+0
   V=0.1941548809452595D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4525585960458567D+0
   B=0.1751165438438091D+0
   V=0.1980078427252384D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4848447779622947D+0
   B=0.2056338306745660D+0
   V=0.2011296284744488D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5160871208276894D+0
   B=0.2359965487229226D+0
   V=0.2035888456966776D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5462112185696926D+0
   B=0.2660430223139146D+0
   V=0.2054516325352142D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5751425068101757D+0
   B=0.2956193664498032D+0
   V=0.2067831033092635D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6028073872853596D+0
   B=0.3245763905312779D+0
   V=0.2076485320284876D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6291338275278409D+0
   B=0.3527670026206972D+0
   V=0.2081141439525255D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3541797528439391D+0
   B=0.2823853479435550D-1
   V=0.1834383015469222D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3908234972074657D+0
   B=0.5741296374713106D-1
   V=0.1889540591777677D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4264408450107590D+0
   B=0.8724646633650199D-1
   V=0.1936677023597375D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4609949666553286D+0
   B=0.1175034422915616D+0
   V=0.1976176495066504D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4944389496536006D+0
   B=0.1479755652628428D+0
   V=0.2008536004560983D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5267194884346086D+0
   B=0.1784740659484352D+0
   V=0.2034280351712291D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5577787810220990D+0
   B=0.2088245700431244D+0
   V=0.2053944466027758D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5875563763536670D+0
   B=0.2388628136570763D+0
   V=0.2068077642882360D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6159910016391269D+0
   B=0.2684308928769185D+0
   V=0.2077250949661599D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6430219602956268D+0
   B=0.2973740761960252D+0
   V=0.2082062440705320D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4300647036213646D+0
   B=0.2916399920493977D-1
   V=0.1934374486546626D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4661486308935531D+0
   B=0.5898803024755659D-1
   V=0.1974107010484300D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5009658555287261D+0
   B=0.8924162698525409D-1
   V=0.2007129290388658D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5344824270447704D+0
   B=0.1197185199637321D+0
   V=0.2033736947471293D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5666575997416371D+0
   B=0.1502300756161382D+0
   V=0.2054287125902493D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5974457471404752D+0
   B=0.1806004191913564D+0
   V=0.2069184936818894D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6267984444116886D+0
   B=0.2106621764786252D+0
   V=0.2078883689808782D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6546664713575417D+0
   B=0.2402526932671914D+0
   V=0.2083886366116359D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5042711004437253D+0
   B=0.2982529203607657D-1
   V=0.2006593275470817D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5392127456774380D+0
   B=0.6008728062339922D-1
   V=0.2033728426135397D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5726819437668618D+0
   B=0.9058227674571398D-1
   V=0.2055008781377608D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6046469254207278D+0
   B=0.1211219235803400D+0
   V=0.2070651783518502D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6350716157434952D+0
   B=0.1515286404791580D+0
   V=0.2080953335094320D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6639177679185454D+0
   B=0.1816314681255552D+0
   V=0.2086284998988521D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5757276040972253D+0
   B=0.3026991752575440D-1
   V=0.2055549387644668D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6090265823139755D+0
   B=0.6078402297870770D-1
   V=0.2071871850267654D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6406735344387661D+0
   B=0.9135459984176636D-1
   V=0.2082856600431965D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6706397927793709D+0
   B=0.1218024155966590D+0
   V=0.2088705858819358D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6435019674426665D+0
   B=0.3052608357660639D-1
   V=0.2083995867536322D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6747218676375681D+0
   B=0.6112185773983089D-1
   V=0.2090509712889637D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

pure SUBROUTINE LD5810(X,Y,Z,W,N)
   real(wp),intent(out) :: X(5810)
   real(wp),intent(out) :: Y(5810)
   real(wp),intent(out) :: Z(5810)
   real(wp),intent(out) :: W(5810)
   integer, intent(out) :: N
   real(wp) :: A,B,V
   !VW
   !VW    LEBEDEV 5810-POINT ANGULAR GRID
   !VW
   !hvd
   !hvd   This subroutine is part of a set of subroutines that generate
   !hvd   Lebedev grids [1-6] for integration on a sphere. The original 
   !hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
   !hvd   translated into fortran by Dr. Christoph van Wuellen.
   !hvd   This subroutine was translated using a C to fortran77 conversion
   !hvd   tool written by Dr. Christoph van Wuellen.
   !hvd
   !hvd   Users of this code are asked to include reference [1] in their
   !hvd   publications, and in the user- and programmers-manuals 
   !hvd   describing their codes.
   !hvd
   !hvd   This code was distributed through CCL (http://www.ccl.net/).
   !hvd
   !hvd   [1] V.I. Lebedev, and D.N. Laikov
   !hvd       "A quadrature formula for the sphere of the 131st
   !hvd        algebraic order of accuracy"
   !hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
   !hvd
   !hvd   [2] V.I. Lebedev
   !hvd       "A quadrature formula for the sphere of 59th algebraic
   !hvd        order of accuracy"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
   !hvd
   !hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
   !hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
   !hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
   !hvd
   !hvd   [4] V.I. Lebedev
   !hvd       "Spherical quadrature formulas exact to orders 25-29"
   !hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
   !hvd
   !hvd   [5] V.I. Lebedev
   !hvd       "Quadratures on a sphere"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
   !hvd       1976, pp. 10-24. 
   !hvd
   !hvd   [6] V.I. Lebedev
   !hvd       "Values of the nodes and weights of ninth to seventeenth 
   !hvd        order Gauss-Markov quadrature formulae invariant under the
   !hvd        octahedron group with inversion"
   !hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
   !hvd       1975, pp. 44-51.
   !hvd
   N=1
   V=0.9735347946175486D-5
   Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.1907581241803167D-3
   Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
   V=0.1901059546737578D-3
   Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1182361662400277D-1
   V=0.3926424538919212D-4
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3062145009138958D-1
   V=0.6667905467294382D-4
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5329794036834243D-1
   V=0.8868891315019135D-4
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7848165532862220D-1
   V=0.1066306000958872D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1054038157636201D+0
   V=0.1214506743336128D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1335577797766211D+0
   V=0.1338054681640871D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1625769955502252D+0
   V=0.1441677023628504D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1921787193412792D+0
   V=0.1528880200826557D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2221340534690548D+0
   V=0.1602330623773609D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2522504912791132D+0
   V=0.1664102653445244D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2823610860679697D+0
   V=0.1715845854011323D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3123173966267560D+0
   V=0.1758901000133069D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3419847036953789D+0
   V=0.1794382485256736D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3712386456999758D+0
   V=0.1823238106757407D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3999627649876828D+0
   V=0.1846293252959976D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4280466458648093D+0
   V=0.1864284079323098D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4553844360185711D+0
   V=0.1877882694626914D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4818736094437834D+0
   V=0.1887716321852025D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5074138709260629D+0
   V=0.1894381638175673D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5319061304570707D+0
   V=0.1898454899533629D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5552514978677286D+0
   V=0.1900497929577815D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5981009025246183D+0
   V=0.1900671501924092D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6173990192228116D+0
   V=0.1899837555533510D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6351365239411131D+0
   V=0.1899014113156229D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6512010228227200D+0
   V=0.1898581257705106D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6654758363948120D+0
   V=0.1898804756095753D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6778410414853370D+0
   V=0.1899793610426402D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6881760887484110D+0
   V=0.1901464554844117D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6963645267094598D+0
   V=0.1903533246259542D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7023010617153579D+0
   V=0.1905556158463228D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.7059004636628753D+0
   V=0.1907037155663528D-3
   Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3552470312472575D-1
   V=0.5992997844249967D-4
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.9151176620841283D-1
   V=0.9749059382456978D-4
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1566197930068980D+0
   V=0.1241680804599158D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2265467599271907D+0
   V=0.1437626154299360D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2988242318581361D+0
   V=0.1584200054793902D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3717482419703886D+0
   V=0.1694436550982744D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4440094491758889D+0
   V=0.1776617014018108D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5145337096756642D+0
   V=0.1836132434440077D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5824053672860230D+0
   V=0.1876494727075983D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6468283961043370D+0
   V=0.1899906535336482D-3
   Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6095964259104373D-1
   B=0.1787828275342931D-1
   V=0.8143252820767350D-4
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.8811962270959388D-1
   B=0.3953888740792096D-1
   V=0.9998859890887728D-4
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1165936722428831D+0
   B=0.6378121797722990D-1
   V=0.1156199403068359D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1460232857031785D+0
   B=0.8985890813745037D-1
   V=0.1287632092635513D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1761197110181755D+0
   B=0.1172606510576162D+0
   V=0.1398378643365139D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2066471190463718D+0
   B=0.1456102876970995D+0
   V=0.1491876468417391D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2374076026328152D+0
   B=0.1746153823011775D+0
   V=0.1570855679175456D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2682305474337051D+0
   B=0.2040383070295584D+0
   V=0.1637483948103775D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2989653312142369D+0
   B=0.2336788634003698D+0
   V=0.1693500566632843D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3294762752772209D+0
   B=0.2633632752654219D+0
   V=0.1740322769393633D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3596390887276086D+0
   B=0.2929369098051601D+0
   V=0.1779126637278296D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3893383046398812D+0
   B=0.3222592785275512D+0
   V=0.1810908108835412D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4184653789358347D+0
   B=0.3512004791195743D+0
   V=0.1836529132600190D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4469172319076166D+0
   B=0.3796385677684537D+0
   V=0.1856752841777379D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4745950813276976D+0
   B=0.4074575378263879D+0
   V=0.1872270566606832D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5014034601410262D+0
   B=0.4345456906027828D+0
   V=0.1883722645591307D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5272493404551239D+0
   B=0.4607942515205134D+0
   V=0.1891714324525297D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5520413051846366D+0
   B=0.4860961284181720D+0
   V=0.1896827480450146D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5756887237503077D+0
   B=0.5103447395342790D+0
   V=0.1899628417059528D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1225039430588352D+0
   B=0.2136455922655793D-1
   V=0.1123301829001669D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1539113217321372D+0
   B=0.4520926166137188D-1
   V=0.1253698826711277D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1856213098637712D+0
   B=0.7086468177864818D-1
   V=0.1366266117678531D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2174998728035131D+0
   B=0.9785239488772918D-1
   V=0.1462736856106918D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2494128336938330D+0
   B=0.1258106396267210D+0
   V=0.1545076466685412D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2812321562143480D+0
   B=0.1544529125047001D+0
   V=0.1615096280814007D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3128372276456111D+0
   B=0.1835433512202753D+0
   V=0.1674366639741759D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3441145160177973D+0
   B=0.2128813258619585D+0
   V=0.1724225002437900D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3749567714853510D+0
   B=0.2422913734880829D+0
   V=0.1765810822987288D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4052621732015610D+0
   B=0.2716163748391453D+0
   V=0.1800104126010751D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4349335453522385D+0
   B=0.3007127671240280D+0
   V=0.1827960437331284D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4638776641524965D+0
   B=0.3294470677216479D+0
   V=0.1850140300716308D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4920046410462687D+0
   B=0.3576932543699155D+0
   V=0.1867333507394938D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5192273554861704D+0
   B=0.3853307059757764D+0
   V=0.1880178688638289D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5454609081136522D+0
   B=0.4122425044452694D+0
   V=0.1889278925654758D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5706220661424140D+0
   B=0.4383139587781027D+0
   V=0.1895213832507346D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5946286755181518D+0
   B=0.4634312536300553D+0
   V=0.1898548277397420D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.1905370790924295D+0
   B=0.2371311537781979D-1
   V=0.1349105935937341D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2242518717748009D+0
   B=0.4917878059254806D-1
   V=0.1444060068369326D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2577190808025936D+0
   B=0.7595498960495142D-1
   V=0.1526797390930008D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2908724534927187D+0
   B=0.1036991083191100D+0
   V=0.1598208771406474D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3236354020056219D+0
   B=0.1321348584450234D+0
   V=0.1659354368615331D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3559267359304543D+0
   B=0.1610316571314789D+0
   V=0.1711279910946440D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3876637123676956D+0
   B=0.1901912080395707D+0
   V=0.1754952725601440D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4187636705218842D+0
   B=0.2194384950137950D+0
   V=0.1791247850802529D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4491449019883107D+0
   B=0.2486155334763858D+0
   V=0.1820954300877716D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4787270932425445D+0
   B=0.2775768931812335D+0
   V=0.1844788524548449D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5074315153055574D+0
   B=0.3061863786591120D+0
   V=0.1863409481706220D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5351810507738336D+0
   B=0.3343144718152556D+0
   V=0.1877433008795068D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5619001025975381D+0
   B=0.3618362729028427D+0
   V=0.1887444543705232D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5875144035268046D+0
   B=0.3886297583620408D+0
   V=0.1894009829375006D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6119507308734495D+0
   B=0.4145742277792031D+0
   V=0.1897683345035198D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2619733870119463D+0
   B=0.2540047186389353D-1
   V=0.1517327037467653D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.2968149743237949D+0
   B=0.5208107018543989D-1
   V=0.1587740557483543D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3310451504860488D+0
   B=0.7971828470885599D-1
   V=0.1649093382274097D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3646215567376676D+0
   B=0.1080465999177927D+0
   V=0.1701915216193265D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3974916785279360D+0
   B=0.1368413849366629D+0
   V=0.1746847753144065D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4295967403772029D+0
   B=0.1659073184763559D+0
   V=0.1784555512007570D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4608742854473447D+0
   B=0.1950703730454614D+0
   V=0.1815687562112174D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4912598858949903D+0
   B=0.2241721144376724D+0
   V=0.1840864370663302D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5206882758945558D+0
   B=0.2530655255406489D+0
   V=0.1860676785390006D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5490940914019819D+0
   B=0.2816118409731066D+0
   V=0.1875690583743703D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5764123302025542D+0
   B=0.3096780504593238D+0
   V=0.1886453236347225D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6025786004213506D+0
   B=0.3371348366394987D+0
   V=0.1893501123329645D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6275291964794956D+0
   B=0.3638547827694396D+0
   V=0.1897366184519868D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3348189479861771D+0
   B=0.2664841935537443D-1
   V=0.1643908815152736D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.3699515545855295D+0
   B=0.5424000066843495D-1
   V=0.1696300350907768D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4042003071474669D+0
   B=0.8251992715430854D-1
   V=0.1741553103844483D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4375320100182624D+0
   B=0.1112695182483710D+0
   V=0.1780015282386092D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4699054490335947D+0
   B=0.1402964116467816D+0
   V=0.1812116787077125D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5012739879431952D+0
   B=0.1694275117584291D+0
   V=0.1838323158085421D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5315874883754966D+0
   B=0.1985038235312689D+0
   V=0.1859113119837737D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5607937109622117D+0
   B=0.2273765660020893D+0
   V=0.1874969220221698D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5888393223495521D+0
   B=0.2559041492849764D+0
   V=0.1886375612681076D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6156705979160163D+0
   B=0.2839497251976899D+0
   V=0.1893819575809276D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6412338809078123D+0
   B=0.3113791060500690D+0
   V=0.1897794748256767D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4076051259257167D+0
   B=0.2757792290858463D-1
   V=0.1738963926584846D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4423788125791520D+0
   B=0.5584136834984293D-1
   V=0.1777442359873466D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4760480917328258D+0
   B=0.8457772087727143D-1
   V=0.1810010815068719D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5085838725946297D+0
   B=0.1135975846359248D+0
   V=0.1836920318248129D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5399513637391218D+0
   B=0.1427286904765053D+0
   V=0.1858489473214328D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5701118433636380D+0
   B=0.1718112740057635D+0
   V=0.1875079342496592D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5990240530606021D+0
   B=0.2006944855985351D+0
   V=0.1887080239102310D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6266452685139695D+0
   B=0.2292335090598907D+0
   V=0.1894905752176822D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6529320971415942D+0
   B=0.2572871512353714D+0
   V=0.1898991061200695D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.4791583834610126D+0
   B=0.2826094197735932D-1
   V=0.1809065016458791D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5130373952796940D+0
   B=0.5699871359683649D-1
   V=0.1836297121596799D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5456252429628476D+0
   B=0.8602712528554394D-1
   V=0.1858426916241869D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5768956329682385D+0
   B=0.1151748137221281D+0
   V=0.1875654101134641D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6068186944699046D+0
   B=0.1442811654136362D+0
   V=0.1888240751833503D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6353622248024907D+0
   B=0.1731930321657680D+0
   V=0.1896497383866979D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6624927035731797D+0
   B=0.2017619958756061D+0
   V=0.1900775530219121D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5484933508028488D+0
   B=0.2874219755907391D-1
   V=0.1858525041478814D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.5810207682142106D+0
   B=0.5778312123713695D-1
   V=0.1876248690077947D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6120955197181352D+0
   B=0.8695262371439526D-1
   V=0.1889404439064607D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6416944284294319D+0
   B=0.1160893767057166D+0
   V=0.1898168539265290D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6697926391731260D+0
   B=0.1450378826743251D+0
   V=0.1902779940661772D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6147594390585488D+0
   B=0.2904957622341456D-1
   V=0.1890125641731815D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6455390026356783D+0
   B=0.5823809152617197D-1
   V=0.1899434637795751D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6747258588365477D+0
   B=0.8740384899884715D-1
   V=0.1904520856831751D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   A=0.6772135750395347D+0
   B=0.2919946135808105D-1
   V=0.1905534498734563D-3
   Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
   N=N-1
   RETURN
END subroutine

end module grid_module
