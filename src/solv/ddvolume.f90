! This file is part of xtb.
! SPDX-Identifier: LGPL-3.0-or-later
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
!
! MS, 09/2023: Initial implementation


!> domain-decomposed volume calculation for molecule
module xtb_solv_ddvolume
    use xtb_mctc_accuracy, only: wp
    implicit none
    private

    public :: ddvolume

    type :: sphere_type
        real(wp) :: x,y,z,r
        integer, allocatable :: neighbours(:)
        integer :: num_neighbours
        real(wp) :: volume
    contains
        procedure :: init => sphere_init
    end type sphere_type

    interface distance2
        module procedure :: distance2_spheres
        module procedure :: distance2_spherepoint
    end interface

contains 

    !> Calculates the volumes of an array of xyz-coordinates with a given radii and resolution based on a custom dd algorithm
    subroutine ddvolume(xyz,r,resolution,volume)
        !> xyz-coordinates (3,nat)
        real(wp), dimension(:,:), intent(in) :: xyz
        !> Radius (nat)
        real(wp), dimension(:), intent(in) :: r
        !> Resolution
        real(wp), intent(in) :: resolution
        !> Volume (nat)
        real(wp), dimension(:), allocatable, intent(out) :: volume

        !> Transfer input to sphere type
        type(sphere_type), allocatable :: spheres(:)
        integer :: i,s,x,y,z
        real(wp) :: x_min,x_max,y_min,y_max,z_min,z_max
    
        allocate(spheres(size(xyz,2)))
        allocate(volume(size(xyz,2)))
        volume = 0.0_wp
        !! Initialize spheres
        do i=1,size(r)
            call spheres(i)%init(xyz(1,i),xyz(2,i),xyz(3,i),r(i),size(xyz,2))
        end do

        call find_neighbours(spheres)
        !$OMP PARALLEL DO DEFAULT(NONE) SHARED(spheres, volume,resolution) PRIVATE(s,x,y,z,x_min,x_max,y_min,y_max,z_min,z_max)
        do s=1,size(spheres)
            x_min=(spheres(s)%x-spheres(s)%r)
            x_max=(spheres(s)%x+spheres(s)%r)
            y_min=(spheres(s)%y-spheres(s)%r)
            y_max=(spheres(s)%y+spheres(s)%r)
            z_min=(spheres(s)%z-spheres(s)%r)
            z_max=(spheres(s)%z+spheres(s)%r)
            do x = floor(x_min/resolution), ceiling(x_max/resolution)
                do y = floor(y_min/resolution), ceiling(y_max/resolution)
                    do z = floor(z_min/resolution), ceiling(z_max/resolution)
                        if (is_inside(x*resolution,y*resolution,z*resolution,spheres,s)) then
                            spheres(s)%volume = spheres(s)%volume + resolution**3
                        end if
                    end do
                end do
            end do
            volume(s) = spheres(s)%volume
      end do
      !$OMP END PARALLEL DO
    end subroutine ddvolume

    !> Initializes a sphere
    subroutine sphere_init(self,x,y,z,r,count)
        class(sphere_type), intent(inout) :: self
        !> Coordinates and radii
        real(wp), intent(in) :: x,y,z,r
        !> Number of maximum possible neighbours (amount of spheres)
        integer, intent(in) :: count

        self%x = x
        self%y = y
        self%z = z
        self%r = r
        self%num_neighbours = 0
        self%volume = 0.0_wp
        allocate(self%neighbours(count))
        self%neighbours = -1

    end subroutine sphere_init

    !> Set up neighbour list for a given array of spheres
    subroutine find_neighbours(spheres)
        !> Array of Sphere Type
        type(sphere_type), intent(inout) :: spheres(:)

        integer :: i,j,k
        real(wp) :: dist2

        do i=1,size(spheres)-1
            do j=i+1,size(spheres)
                dist2 = distance2(spheres(i),spheres(j))
                if (dist2 < (spheres(i)%r+spheres(j)%r)**2) then
                    do k=1,size(spheres(i)%neighbours)
                        if (spheres(i)%neighbours(k) == -1) then
                            spheres(i)%neighbours(k) = j
                            spheres(i)%num_neighbours = spheres(i)%num_neighbours + 1
                            exit
                        end if
                    end do
                    do k=1,size(spheres(j)%neighbours)
                        if (spheres(j)%neighbours(k) == -1) then
                            spheres(j)%neighbours(k) = i
                            spheres(j)%num_neighbours = spheres(j)%num_neighbours + 1
                            exit
                        end if
                    end do
                end if
            end do
        end do
    end subroutine

    !> Quadratic distance between two spheres
    function distance2_spheres(sphere1,sphere2) result(distance2)
        type(sphere_type), intent(in) :: sphere1,sphere2
        real(wp) :: distance2

        distance2 = (sphere1%x-sphere2%x)**2 + (sphere1%y-sphere2%y)**2 + (sphere1%z-sphere2%z)**2
    end function distance2_spheres

    !> Quadratic distance between a sphere and a point
    function distance2_spherepoint(sphere1,x,y,z) result(distance2)
        type(sphere_type), intent(in) :: sphere1
        real(wp), intent(in) :: x,y,z
        real(wp) :: distance2

        distance2 = (sphere1%x-x)**2 + (sphere1%y-y)**2 + (sphere1%z-z)**2
    end function distance2_spherepoint

    !> Checks if a point is inside a sphere or not
    function is_inside(x,y,z,spheres,this) result(inside)
        implicit none
        !> Array of Sphere Type (all spheres)
        type(sphere_type), intent(in) :: spheres(:)
        !> Index of the sphere to check
        integer, intent(in) :: this
        !> Coordinates of the point to check
        real(wp), intent(in) :: x,y,z
        !> Is it inside?
        logical :: inside


        integer :: i,n
        real(wp) :: dist2, distn2


        !! Check if point is even inside the sphere
        dist2 = distance2(spheres(this),x,y,z)
        if (dist2 <= spheres(this)%r**2) then
            inside = .true.
        else
            inside = .false.
            return
        end if
    
        !! Check if point is nearer to one of the neighbours than to the sphere itself
        do i=1,spheres(this)%num_neighbours
            n = spheres(this)%neighbours(i)
            distn2 = distance2(spheres(n),x,y,z)
            if (distn2 < dist2 .and. distn2 <= spheres(n)%r**2) then
                inside = .false.
                return
            end if
        end do
      end function is_inside


end module xtb_solv_ddvolume




