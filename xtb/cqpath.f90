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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Repositionierung entlang des Reaktionspfads
!
! Methoden:
!   m_repos
!       1 - symmetrische Verteilung
!       2 - unsymmetrische Verteilung, höhere Dichte am pot. TS punkten
!   m_spline
!       1 - cubischer spline
!       2 - gedämpfter cubischer spline (constrained cubic spline)
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cqpath_spline_repos(nat,sn,m_repos,m_spline,xyz,e,g)
    implicit none
    integer nat,sn,m_repos,m_spline
    real*8 xyz(3,nat,sn)
    real*8 e(sn)
    real*8 g(sn)

    integer :: i,j,k,l
    real(8), allocatable, dimension(:) :: strecke_X
    real(8), allocatable, dimension(:) :: strecke_N
    real(8), allocatable, dimension(:,:) :: sp_atom
    real(8), allocatable, dimension(:) :: xyz_0,xyz_1,xyz_diff
    real(8), allocatable, dimension(:) :: x,y,x2,y2
    real(8) :: strecke_X_inc, temp_strecke_X

    allocate (strecke_X(sn))
    allocate (strecke_N(sn))
    allocate (sp_atom(3,sn))
    allocate (xyz_0(3))
    allocate (xyz_1(3))
    allocate (xyz_diff(3))
    allocate (x(sn))
    allocate (y(sn))
    allocate (x2(sn))
    allocate (y2(sn))

!     print *,'cqpath_spline_repos'
!     print *,'repos  : ',m_repos
!     print *,'spline : ',m_spline

!       loop über alle atome um ein array der punte eines atoms entlang des pfads zu bekommen
    do i=1,nat
!           loop über den pfad
    do j=1,sn
        sp_atom(1,j) = xyz(1,i,j)
        sp_atom(2,j) = xyz(2,i,j)
        sp_atom(3,j) = xyz(3,i,j)
    end do
!           X-Achse für den Spline generieren = Strecke
    xyz_0(1) = sp_atom(1,1)
    xyz_0(2) = sp_atom(2,1)
    xyz_0(3) = sp_atom(3,1)
    strecke_X(1) = 0.0
    do j=2,sn
        xyz_1(1) = sp_atom(1,j)
        xyz_1(2) = sp_atom(2,j)
        xyz_1(3) = sp_atom(3,j)
        xyz_diff = xyz_1 - xyz_0
        xyz_0 = xyz_1
        strecke_X(j) = strecke_X(j-1) + sqrt(sum(xyz_diff**2))
    end do

!           Für x,y,z in abh. von strecke_X einen cube Spline machen
    do k=1,3
        ! X = strecke_X , Y = [x|y|z]
        do l=1,sn
            x(l) = strecke_X(l)
            y(l) = sp_atom(k,l)
        end do
        
! Repositionieren auf x
        call cqpath_repos(sn,m_repos,x,e,g,x2)
        
! y2-Werte auf dem spline an den repositionierten Punkten x2 bestimmen.
        if (m_spline == 1) then
            call cqpath_cube_spline(sn,x,y,x2,y2)
        end if
        if (m_spline == 2) then
            call cqpath_damped_cube_spline(sn,x,y,x2,y2)
        end if
        
        do l=1,sn
            sp_atom(k,l) = y2(l)
        end do
    end do
!           Neuen Punkte übertragen
    do j=2,sn-1
        xyz(1,i,j) = sp_atom(1,j)
        xyz(2,i,j) = sp_atom(2,j)
        xyz(3,i,j) = sp_atom(3,j)
    end do
    end do

    deallocate (strecke_X)
    deallocate (strecke_N)
    deallocate (sp_atom)
    deallocate (xyz_0)
    deallocate (xyz_1)
    deallocate (xyz_diff)
    deallocate (x)
    deallocate (y)
    deallocate (x2)
    deallocate (y2)
end

!!!!!!!!!!!!!! REPOSITIONIERUNG AUF DEM PFAD !!!!!!!!!!!!!!!!
subroutine cqpath_repos(n,m_repos,x,e,g,x2)
! in    n  - dimension von array x,y
!       x  - array x
!       e  - energien entlang des pfades
!       g  - gradienten entlang des pfades
! out   x2 - symmetrisch Punkte aus array x
!
! Methoden:
!   m_repos
!       1 - symmetrische Verteilung
!       2 - unsymmetrische Verteilung, höhere Dichte am pot. TS punkten
!
!
!
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m_repos
    real(8), intent(in) :: x(n)
    real(8), intent(in) :: e(n)
    real(8), intent(in) :: g(n)
    real(8), intent(out) :: x2(n)
    
    real(8), allocatable :: e2(:)
    real(8) x_wert, x_anfang, x_ende, x_schritt
    real(8) emin,e2sum,temp_wert
    integer :: i
    
    ! Symmetrisch verteilen
    ! dx(i) = strecke/punkte
    if (m_repos == 1) then
        x_anfang = x(1)
        x_ende = x(n)
        x_schritt = (x_ende-x_anfang)/(REAL(n,8)-1.0)
        
        x_wert = x_anfang
        x2(1) = x_anfang
        do i=2,n-1
            x_wert = x_wert + x_schritt
            x2(i) = x_wert
        end do
        x2(n) = x_ende
    end if
    
    ! Asymmetrisch verteilen, in Abh. von e
    ! 1. e -> e' = -1.0*e                       ! umkehren
    ! 2. e' -> e'' = e' - e_min * 0.99          ! verschieben damit e_ts nicht = 0 ergibt
    ! 3. dx(i) = e''(i)/e''_summe               ! normieren
    if (m_repos == 2) then
        allocate(e2(n))
    
        x_anfang = x(1)
        x_ende = x(n)
        x2(1) = x_anfang
        x2(n) = x_ende
        x_schritt = x_ende-x_anfang

        emin = minval(e)*0.99
        do i=1,n
            e2(i) = (-1.0)*e(i)-emin
        end do
        e2sum = sum(e2)
        do i=1,n
            e2(i) = e2(i)/e2sum
        end do
        
        
        do i=2,n-1
!             x2(i) = x_anfang+x_schritt*e2(i)
            x_wert = x_wert + x_schritt*e2(i)
            x2(i) = x_wert
        end do
        
        deallocate(e2)
    end if
    
    ! Asymmetrisch verteilen, in Abh. von e
    ! 1. e -> e' = e*e                          ! umkehren
    ! 2. e' -> e'' = e' - e_min * 0.99          ! verschieben damit e_ts nicht = 0 ergibt
    ! 3. dx(i) = e''(i)/e''_summe               ! normieren
    if (m_repos == 3) then
        allocate(e2(n))
    
        x_anfang = x(1)
        x_ende = x(n)
        x2(1) = x_anfang
        x2(n) = x_ende
        x_schritt = x_ende-x_anfang

        emin = minval(e)*0.99
        do i=1,n
            e2(i) = e(i)*e(i)-emin
        end do
        e2sum = sum(e2)
        do i=1,n
            e2(i) = e2(i)/e2sum
        end do
        
        
        do i=2,n-1
!             x2(i) = x_anfang+x_schritt*e2(i)
            x_wert = x_wert + x_schritt*e2(i)
            x2(i) = x_wert
        end do
        
        deallocate(e2)
    end if
    
    ! Asymmetrisch verteilen, in Abh. von e
    ! 1. e -> e' = -1.0*e                       ! umkehren
    ! 2. e' -> e'' = (e' - e_min * 0.99)**2     ! verschieben damit e_ts nicht = 0 ergibt
    ! 3. dx(i) = e''(i)/e''_summe               ! normieren
    if (m_repos == 4) then
        allocate(e2(n))
    
        x_anfang = x(1)
        x_ende = x(n)
        x2(1) = x_anfang
        x2(n) = x_ende
        x_schritt = x_ende-x_anfang

        emin = minval(e)*0.99
        do i=1,n
            temp_wert = (-1.0)*e(i)-emin
            e2(i) = temp_wert*temp_wert
        end do
        e2sum = sum(e2)
        do i=1,n
            e2(i) = e2(i)/e2sum
        end do
        
        
        do i=2,n-1
!             x2(i) = x_anfang+x_schritt*e2(i)
            x_wert = x_wert + x_schritt*e2(i)
            x2(i) = x_wert
        end do
        
        deallocate(e2)
    end if
end

!!!!!!!!!!!!!! Y-Werte entlang der Repositionierung !!!!!!!!!!!!!!!!!!!!
subroutine cqpath_cube_spline(n,x,y,x2,y2)
! in    n  - dimension von array x,y
!       x  - array x
!       y  - array y
! out   x2 - symmetrisch Punkte aus array xx (eigentlich überflüssig)
!       y2 - symmetrische Punkte aus array yy
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: x(n),y(n),x2(n)
    real(8), intent(out) :: y2(n)
    
    real(8), dimension(:), allocatable :: ypp
    real(8) y_wert, y_wert_ypval, y_wert_yppval
    integer :: i

    allocate (ypp(n))
    call spline_cubic_set ( n, x, y, 0, 0.0d0, 0, 0.0d0, ypp )
    y2(1) = y(1)
    do i=2,n-1
        call spline_cubic_val ( n, x, y, ypp, x2(i), y_wert, y_wert_ypval, y_wert_yppval )
        y2(i) = y_wert
    end do
    y2(n) = y(n)
    deallocate (ypp)
    return
end

subroutine cqpath_damped_cube_spline(n,x,y,x2,y2)
! in    n  - dimension von array x,y
!       x  - array x
!       y  - array y
! out   x2 - symmetrisch Punkte aus array xx (eigentlich überflüssig)
!       y2 - symmetrische Punkte aus array yy
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: x(n),y(n),x2(n)
    real(8), intent(out) :: y2(n)
    
    real(8), dimension(:,:), allocatable :: abcd
    real(8) y_wert
    integer :: i
    
    allocate (abcd(n,4))    
    call cqpath_damped_cube_spline_set(n,x,y,abcd)
    
    y2(1) = y(1)
    do i=2,n-1
        call cqpath_damped_cube_spline_val(n,abcd,x,x2(i),y_wert)
        y2(i) = y_wert
    end do
    y2(n) = y(n)
    
    deallocate (abcd)
    return
end

!!!!!!!!!!!!!!!!!!! SPLINE FUNKTIONEN !!!!!!!!!!!!!!!!!!!!!!!!

subroutine cqpath_damped_cube_spline_set(n,x,y,abcd)
! in    n  - dimension von array x,y
!       x  - array x
!       y  - array y
! out   abcd - array of polynom coefficients
!
! Based on Constrained Cubic Spline Interpolation for Chemical Engineering Application
! by CJC Kruger
!
! Rewritten and optimized in fortran
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: x(n)
    real(8), intent(in) :: y(n)
    real(8), intent(out) :: abcd(n,4)

    real(8) fs_x0,fs_x1,fss_x0,fss_x1
    real(8) x0,x1,x2,y0,y1,y2
    real(8) x1_x0,y1_y0,x0_2,x0_3,x1_x0_2,x1_2
    real(8) steigung
    real(8), dimension(:), allocatable :: dx,dy,f1
    
    integer :: i,j,nf
    
    nf = n-1
    allocate (dx(nf))
    allocate (dy(nf))
    allocate (f1(n))
    
    do i=1,nf
        dx(i) = x(i+1)-x(i)
        dy(i) = y(i+1)-y(i)
    end do
    
    do i=2,nf
        steigung = dy(i-1)*dy(i)
        if (steigung > 0.0) then
            f1(i) = 2.0/(dx(i)/dy(i)+dx(i-1)/dy(i-1))
        else if (steigung <= 0.0) then
            f1(i) = 0.0
        end if
    end do
    
    f1(1) = 3.0*dy(1)/(2.0*dx(1)) - f1(2)/2.0
    f1(n) = 3.0*dy(n-1)/(2.0*dx(n-1)) - f1(n-1)/2.0
    
    do i=2,n
        x1_x0 = dx(i-1)
        y1_y0 = dy(i-1)
        x0_2 = x(i-1)*x(i-1)
        x0_3 = x0_2*x(i-1)
        x1_x0_2 = x1_x0*x1_x0
        x1_2 = x(i)*x(i)
        
        fss_x0 = -2.0*(f1(i) + 2.0*f1(i-1))/(x1_x0)+6.0*(y1_y0)/x1_x0_2
        fss_x1 = 2.0*(2.0*f1(i) + f1(i-1))/(x1_x0)-6.0*(y1_y0)/x1_x0_2
        
        abcd(i,4) = (fss_x1 - fss_x0)/(6.0*(x1_x0))
        abcd(i,3) = (x(i)*fss_x0-x(i-1)*fss_x1)/(2.0*x1_x0)
        abcd(i,2) = ((y1_y0)-abcd(i,3)*(x1_2-x0_2)-abcd(i,4)*(x1_2*x(i)-x0_3))/(x1_x0)
        abcd(i,1) = y(i-1) - abcd(i,2)*x(i-1) - abcd(i,3) * x0_2 - abcd(i,4) * x0_3
    end do

    deallocate(dx)
    deallocate(dy)
    deallocate(f1)
    return
end

subroutine cqpath_damped_cube_spline_val(n,abcd,x,x1,y1)
! in    n    - dimension von array x,y
!       abcd - array of polynom coefficients
!       x    - x
! out   y    - f(x,abcd)
!
! Based on Constrained Cubic Spline Interpolation for Chemical Engineering Application
! by CJC Kruger
!
! Rewritten and optimized in fortran
    integer, intent(in) :: n
    real(8), intent(in) :: abcd(n,4)
    real(8), intent(in) :: x(n)
    real(8), intent(in) :: x1
    real(8), intent(out) :: y1

    integer :: j
    real(8) :: xx
    
    do j=2,n
        if (x1 < x(j)) then
            if (x1 >= x(j-1)) then
                xx = x1*x1
                y1 = abcd(j,1)+abcd(j,2)*x1+abcd(j,3)*xx+abcd(j,4)*xx*x1
            end if
        end if
    end do
    return
end

subroutine cqpath_interpolate_nm(k,von,bis,nat,nstr,file_xyz,pn,xyz)
! Generierung der Strukturen zwischen Start und Zielstruktur
! Linerare interpolation
    implicit none
    integer, intent(in) :: k,von,bis
    integer, intent(in) :: nat,nstr,pn
    real(8), dimension(:,:,:), intent(in) :: file_xyz(3,nat,nstr)
    real(8), dimension(:,:,:), intent(inout) :: xyz(3,nat,pn)

    real(8), allocatable, dimension(:,:) :: xyz_diff,xyz_temp,xyz_start,xyz_ziel
    integer :: i,j

    allocate (xyz_diff(3,nat))
    allocate (xyz_temp(3,nat))
    allocate (xyz_start(3,nat))
    allocate (xyz_ziel(3,nat))
    
    do j=1,nat
        xyz_start(1,j) = file_xyz(1,j,k-1)
        xyz_start(2,j) = file_xyz(2,j,k-1)
        xyz_start(3,j) = file_xyz(3,j,k-1)
        xyz_ziel(1,j) = file_xyz(1,j,k)
        xyz_ziel(2,j) = file_xyz(2,j,k)
        xyz_ziel(3,j) = file_xyz(3,j,k)
    end do
    
    xyz_diff = xyz_ziel - xyz_start
    xyz_diff = xyz_diff / (bis-von)
    
    xyz_temp = xyz_start
    do j=1,nat
        xyz(1,j,von) = xyz_start(1,j)
        xyz(2,j,von) = xyz_start(2,j)
        xyz(3,j,von) = xyz_start(3,j)
    end do
    do i=von+1,bis-1
      xyz_temp = xyz_temp + xyz_diff
      do j=1,nat
        xyz(1,j,i) = xyz_temp(1,j)
        xyz(2,j,i) = xyz_temp(2,j)
        xyz(3,j,i) = xyz_temp(3,j)
      end do
    end do
    do j=1,nat
        xyz(1,j,bis) = xyz_ziel(1,j)
        xyz(2,j,bis) = xyz_ziel(2,j)
        xyz(3,j,bis) = xyz_ziel(3,j)
    end do
    
    deallocate (xyz_diff)
    deallocate (xyz_temp)
    deallocate (xyz_start)
    deallocate (xyz_ziel)
end

subroutine cqpath_read_pathfile_parameter(fn,nl,nat,nstr)
    implicit none
    character*20, intent(in) :: fn
    integer, intent(out) :: nl,nat,nstr
    integer :: u
    integer :: i,so,sr
    character*100 :: a
    character*20 :: fn2
    real(8) :: b,c,d
    
    fn2 = TRIM(adjustl(fn))
!     print *,'#',TRIM(fn2),'#'
    call open_file(u,trim(fn2),'r')
    if (u.eq.-1) return
    read (u,*,iostat=sr) nat
    i = 1
    do
      read (u,fmt='(A)',iostat=sr) a
      if (sr<0) exit
      if (sr>0) stop 'error in cqpath_read_pathfile_parameter'   
      i = i + 1
    end do
    
    nl = i
    nstr = nl/(nat+2)
    
    call close_file (u)
    return
end

subroutine cqpath_read_pathfile(fn,nl,nat,nstr,file_xyz,iat,energy)
    implicit none
    character*20, intent(in) :: fn
    integer, intent(in) :: nat,nl,nstr
    real(8), dimension(:,:,:), intent(inout) :: file_xyz(3,nat,nstr)
    real(8), dimension(:,:,:), intent(inout) :: energy(nstr)
    integer, dimension(:), intent(out) :: iat(nat)
    
    integer :: u
    integer :: i,j,k,sr,so,x,nn
    integer :: cqpathe2i
    real*8  :: xx(10)
    character*100 :: a
    character*2 :: b
    character*20 :: fn2

    fn2 = TRIM(adjustl(fn))
!   print *,'*',TRIM(fn2),'*',nl,nat,nstr
    call open_file(u,TRIM(fn2),'r')
    if (u.eq.-1) return
    read (u,*,iostat=sr) x
    read (u,fmt='(A)',iostat=sr) a
    call readl(a,xx,nn)
    k=1
    energy(k)=xx(1)
    do j=1,nat
        read (u,*) b,file_xyz(1,j,k),file_xyz(2,j,k),file_xyz(3,j,k)
        iat(j) = cqpathe2i( b )
!         print *, iat(j),file_xyz(1,j,k),file_xyz(2,j,k),file_xyz(3,j,k)
    end do
    do k=2,nstr
        read (u,*,iostat=sr) x
        read (u,fmt='(A)',iostat=sr) a
        call readl(a,xx,nn)
        energy(k)=xx(1)
        do j=1,nat
            read (u,*) b,file_xyz(1,j,k),file_xyz(2,j,k),file_xyz(3,j,k)
!             print *, iat(j),file_xyz(1,j,k),file_xyz(2,j,k),file_xyz(3,j,k)
        end do
    end do
    call close_file (u)
    return
end

character(2) function cqpathtohigher( s )
! wandelt Kleinbuchstaben in Großbuchstaben um
    implicit none
    character(2), intent(in) :: s
    character(2)             :: sout
    integer :: ic, i
    character(26), Parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    
    sout = s
    do i = 1, LEN_TRIM(s)
        ic = INDEX(low, s(i:i))
        if (ic > 0) sout(i:i) = high(ic:ic)
    end do
    cqpathtohigher = sout
end function cqpathtohigher

integer function cqpathe2i( cin )
! wandelt eine Zeichenkette in iat um
    implicit none
    character(2), intent(in) :: cin
    character(2) :: cqpathtohigher,c
    integer :: iout
    c = cqpathtohigher(cin)
    if (c == 'H') iout = 1
    if (c == 'HE') iout = 2
    if (c == 'LI') iout = 3
    if (c == 'BE') iout = 4
    if (c == 'B') iout = 5
    if (c == 'C') iout = 6
    if (c == 'N') iout = 7
    if (c == 'O') iout = 8
    if (c == 'F') iout = 9
    if (c == 'NE') iout = 10  
    if (c == 'NA') iout = 11
    if (c == 'MG') iout = 12
    if (c == 'AL') iout = 13
    if (c == 'SI') iout = 14
    if (c == 'P') iout = 15
    if (c == 'S') iout = 16
    if (c == 'CL') iout = 17
    if (c == 'AR') iout = 18
    if (c == 'K') iout = 19
    if (c == 'CA') iout = 20
    if (c == 'SC') iout = 21
    if (c == 'TI') iout = 22
    if (c == 'V') iout = 23
    if (c == 'CR') iout = 24
    if (c == 'MN') iout = 25
    if (c == 'FE') iout = 26
    if (c == 'CO') iout = 27
    if (c == 'NI') iout = 28
    if (c == 'CU') iout = 29
    if (c == 'ZN') iout = 30
    if (c == 'GA') iout = 31
    if (c == 'GE') iout = 32
    if (c == 'AS') iout = 33
    if (c == 'SE') iout = 34
    if (c == 'BR') iout = 35
    if (c == 'KR') iout = 36
    if (c == 'RB') iout = 37
    if (c == 'SR') iout = 38
    if (c == 'Y') iout = 39
    if (c == 'ZR') iout = 40
    if (c == 'NB') iout = 41
    if (c == 'MO') iout = 42
    if (c == 'TC') iout = 43
    if (c == 'RU') iout = 44
    if (c == 'RH') iout = 45
    if (c == 'PD') iout = 46
    if (c == 'AG') iout = 47
    if (c == 'CD') iout = 48
    if (c == 'IN') iout = 49
    if (c == 'SN') iout = 50
    if (c == 'SB') iout = 51
    if (c == 'TE') iout = 52
    if (c == 'I') iout = 53
    if (c == 'XE') iout = 54
    if (c == 'CS') iout = 55
    if (c == 'BA') iout = 56
    if (c == 'LA') iout = 57
    if (c == 'CE') iout = 58
    if (c == 'PR') iout = 59
    if (c == 'ND') iout = 60
    if (c == 'PM') iout = 61
    if (c == 'SM') iout = 62
    if (c == 'EU') iout = 63
    if (c == 'GD') iout = 64
    if (c == 'TB') iout = 65
    if (c == 'DY') iout = 66
    if (c == 'HO') iout = 67
    if (c == 'ER') iout = 68
    if (c == 'TM') iout = 69
    if (c == 'YB') iout = 70
    if (c == 'LU') iout = 71
    if (c == 'HF') iout = 72
    if (c == 'TA') iout = 73
    if (c == 'W') iout = 74
    if (c == 'RE') iout = 75
    if (c == 'OS') iout = 76
    if (c == 'IR') iout = 77
    if (c == 'PT') iout = 78
    if (c == 'AU') iout = 79
    if (c == 'HG') iout = 80
    if (c == 'TL') iout = 81
    if (c == 'PB') iout = 82
    if (c == 'BI') iout = 83
    if (c == 'PO') iout = 84
    if (c == 'AT') iout = 85
    if (c == 'RN') iout = 86
    cqpathe2i = iout
    return
end function cqpathe2i

