!-----------------------------------------------------------------------------
!
! compliance matrix
!
! BFS-based connectivity using covalent radii (Alvarez 2008)
! Automatic NB/NC selection with angle quality control (>20 deg)
! LINBEND detection based on true chemical bond angles
! Fixed (e1,e2) frame for Decius linear bend coordinates
! Dihedral wrapping to (-pi,pi] within the finite-difference calculations
!
! compliance.f90  -- C = F^{-1} via SVD pseudoinverse of G and F,
!                    handles redundant coordinate sets (symmetric tops etc.)
!
! run_compl.f90   -- Driver with ordered output
!                    (Bonds -> Angles -> Dihedrals -> Linear bend coordinates)
!                    Full compliance matrix written to compliance.dat:
!                    diagonal C_ii, 1/C_ii, and top-10 off-diagonal couplings
!                    per coordinate, sorted by |C_ij| descending.
!
! Ref.: K. Brandhorst, J. Grunenberg, Chem. Soc. Rev. 37 (2008), 1558.
!       J. Grunenberg, Chem. Sci. 6 (2015), 4086.
!
! SG (with Claude), 05/26
!-----------------------------------------------------------------------------

   subroutine compliance_driver(n, at, xyz, hess)
      use compliance_matrix
      use zmat_bmatrix
      integer,  intent(in) :: n
      integer,  intent(in) :: at(n)
      real(8),  intent(in) :: xyz(3,n)
      real(8),  intent(in) :: hess(3*n,3*n)
      integer :: nint, istat
      integer, allocatable :: na(:), nb(:), nc(:), ctype(:), qoff(:)
      real(8), allocatable :: e1f(:,:), e2f(:,:), bmat(:,:), compl(:,:)

      write(*,*)
      write(*,*) '             ======================================='
      write(*,*) '             |                                     |'
      write(*,*) '             |       compliance constants          |'
      write(*,*) '             |                                     |'
      write(*,*) '             ======================================='
      write(*,*)
      write(*,*) 'Ref.: K. Brandhorst, J. Grunenberg, Chem. Soc. Rev. 37 (2008), 1558.'

      allocate(na(n), nb(n), nc(n), ctype(n), qoff(n), e1f(3,n), e2f(3,n))
      call setup_zmat(xyz, n, at, na, nb, nc, ctype, e1f, e2f, qoff, nint)
      allocate(bmat(nint, 3*n), compl(nint, nint))
      call compute_bmatrix(xyz, n, na, nb, nc, ctype, e1f, e2f, qoff, nint, bmat, 5d-3) ! numerical deriv. step
      call compute_compliance(hess, bmat, n, nint, compl, istat)
      call print_compl(n, at, xyz, na, nb, nc, ctype, e1f, e2f, qoff, nint, compl)

   end subroutine compliance_driver


!-----------------------------------------------------------------------------
! PRINT_COMPL  -- stdout: diagonal elements, ordered bonds/angles/dihedrals/bends
!                file:   full matrix via write_compliance_dat
!-----------------------------------------------------------------------------

   subroutine print_compl(n, at, xyz, na, nb, nc, ctype, e1f, e2f, qoff, nint, C)
      use zmat_bmatrix
      integer,  intent(in) :: n, nint
      integer,  intent(in) :: at(n), na(n), nb(n), nc(n), ctype(n), qoff(n)
      real(8),  intent(in) :: xyz(3,n), e1f(3,n), e2f(3,n), C(nint,nint)
      integer       :: i, o, k
      integer       :: idx_ord(nint)
      real(8)       :: q(nint), cc, mu, freq
      character(20) :: s
      character(2)  :: elsym
      external        elsym
      real(8)       :: atmass
      external        atmass
      ! local mode frequency: nu_a = fac * sqrt(k_a[a.u.] / mu[a.u.])
      !   k_a = 1/C_ii  in Eh/a0^2
      !   mu  = m_A*m_B/(m_A+m_B)  in amu  (NOT converted to me)
      !   fac = 1/(2*pi*c) * sqrt(Eh/(a0^2 * amu))
      !       = 219474.631 * sqrt(me/amu)
      !       = 219474.631 / sqrt(1822.888)  =  5140.487  -- WRONG if mu in me
      !   Correct: keep mu in amu, use fac below (CODATA 2018):
      !     1/(2*pi*c) * sqrt(Eh/(a0^2*amu)) = 5140.4869 cm^-1
      !   Derivation:
      !     Eh/a0^2 = 1556.89 N/m
      !     1 amu   = 1.66054e-27 kg
      !     sqrt(1556.89/1.66054e-27) / (2*pi*2.99792e10) = 5140.49 cm^-1
      real(8), parameter :: fac = 5140.4869d0  ! cm^-1, mu must be in amu

      call cartesian_to_int(xyz, n, na, nb, nc, ctype, e1f, e2f, qoff, nint, q)

      write(*,*)
      write(*,*) 'units: Hartree, Bohr, radian'
      write(*,*) '1 Eh/a0^2 (1/C = relaxed force constant) = 15.570 N/cm'
      write(*,*) 'local mode frequency nu_loc = 5140.487*sqrt(1/(mu[amu]*C[a.u.])) cm^-1'
      write(*,*) 'Ref.: Cremer, Kraka, Zou, J. Chem. Theory Comput. 8 (2012) 2864.'

      k = 0

      ! --- 1) bond stretches ---
      ! local mode frequency (Kraka/Cremer = 1/C_ii route):
      ! nu_a = fac * sqrt(1/(mu_AB * C_ii))   [cm^-1]
      ! mu_AB = m_A*m_B/(m_A+m_B)  in atomic mass units -> convert to me
      write(*,'(a)') &
         '     type                atoms                   coord value' // &
         '       C       1/C    nu_loc/cm-1'
      do i = 2, n
         k = k + 1;  o = qoff(i);  idx_ord(k) = o+1
         cc = C(o+1,o+1);  s = 'bond stretch'
         mu = atmass(at(i)) * atmass(at(na(i))) &
                   / (atmass(at(i)) + atmass(at(na(i))))   ! reduced mass in amu
         freq = fac * sqrt(1d0 / (mu * cc))
         write(*,'(i4,1x,a14,2(a2,i3,3x),16x,4f10.2)') &
            k, s, elsym(at(i)),i, elsym(at(na(i))),na(i), &
            q(o+1), cc, 1d0/cc, freq
      end do

      ! --- 2) valence angles ---
      do i = 2, n
         if (ctype(i)/=COORD_ANGLE .and. ctype(i)/=COORD_DIHEDRAL) cycle
         k = k + 1;  o = qoff(i);  idx_ord(k) = o+2
         cc = C(o+2,o+2);  s = 'angle'
         write(*,'(i4,1x,a14,3(a2,i3,3x),8x,3f10.4)') &
            k, s, elsym(at(i)),i, elsym(at(na(i))),na(i), elsym(at(nb(i))),nb(i), &
            q(o+2), cc, 1d0/cc
      end do

      ! --- 3) dihedrals ---
      do i = 2, n
         if (ctype(i)/=COORD_DIHEDRAL) cycle
         k = k + 1;  o = qoff(i);  idx_ord(k) = o+3
         cc = C(o+3,o+3);  s = 'dihedral'
         write(*,'(i4,1x,a14,4(a2,i3,3x),3f10.4)') &
            k, s, elsym(at(i)),i, elsym(at(na(i))),na(i), &
            elsym(at(nb(i))),nb(i), elsym(at(nc(i))),nc(i), &
            q(o+3), cc, 1d0/cc
      end do

      ! --- 4) linear bending coordinates ---
      do i = 2, n
         if (ctype(i)/=COORD_LINBEND) cycle
         o = qoff(i)
         k = k + 1;  idx_ord(k) = o+2
         cc = C(o+2,o+2);  s = 'lin.bend(e1)'
         write(*,'(i4,1x,a14,3(a2,i3,3x),8x,3f10.4)') &
            k, s, elsym(at(nb(i))),nb(i), elsym(at(na(i))),na(i), elsym(at(i)),i, &
            q(o+2), cc, 1d0/cc
         k = k + 1;  idx_ord(k) = o+3
         cc = C(o+3,o+3);  s = 'lin.bend(e2)'
         write(*,'(i4,1x,a14,3(a2,i3,3x),8x,3f10.4)') &
            k, s, elsym(at(nb(i))),nb(i), elsym(at(na(i))),na(i), elsym(at(i)),i, &
            q(o+3), cc, 1d0/cc
      end do

      call write_compliance_dat(n, at, na, nb, nc, ctype, qoff, nint, C, idx_ord, k)

   end subroutine print_compl


!-----------------------------------------------------------------------------
! WRITE_COMPLIANCE_DAT
!   Writes compliance.dat.  For each coordinate (in output order):
!     - header line with type/atoms, C_ii, 1/C_ii
!     - top NCOUP off-diagonal |C_ij| couplings, sorted descending
!-----------------------------------------------------------------------------

   subroutine write_compliance_dat(n, at, na, nb, nc, ctype, qoff, &
                                   nint, C, idx_ord, ncoord)
      use zmat_bmatrix
      integer,  intent(in) :: n, nint, ncoord
      integer,  intent(in) :: at(n), na(n), nb(n), nc(n), ctype(n), qoff(n)
      real(8),  intent(in) :: C(nint,nint)
      integer,  intent(in) :: idx_ord(ncoord)
      !--- local ---
      integer, parameter :: NCOUP = 20
      integer, parameter :: iunit = 42
      integer  :: i, j, p, q, ii, jj, nc_act
      integer  :: jsort(ncoord)
      real(8)  :: aval(ncoord), tmp_r
      integer  :: tmp_i
      character(20) :: lbl(ncoord)

      ! Build labels (uses internal elsym, no external dependency)
      call build_labels(n, at, na, nb, nc, ctype, qoff, ncoord, idx_ord, lbl)

      open(iunit, file='compliance.dat', status='replace')

      write(iunit,'(a)') '#'
      write(iunit,'(a)') '# compliance.dat'
      write(iunit,'(a)') '#'
      write(iunit,'(a)') '# units: C   in Bohr^2/Hartree (a0^2/Eh)'
      write(iunit,'(a)') '#        1/C in Eh/a0^2  (relaxed force constant)'
      write(iunit,'(a)') '#        conversion: 1 Eh/a0^2 = 15.570 N/cm'
      write(iunit,'(a)') '#'
      write(iunit,'(a)') '# Ref.: K. Brandhorst, J. Grunenberg,'
      write(iunit,'(a)') '#       Chem. Soc. Rev. 37 (2008) 1558.'
      write(iunit,'(a)') '#'
      write(iunit,'(a,i6)') '# number of internal coordinates :', ncoord
      write(iunit,'(a,i4)')  '# top couplings shown per coord  :', NCOUP
      write(iunit,'(a)') '#'
      write(iunit,'(a)') '# columns: coord_j  label_j  C_ij  [<-> coord_i label_i]'
      write(iunit,'(a)') '#'

      do i = 1, ncoord
         ii = idx_ord(i)

         write(iunit,'(a)')  ''
         write(iunit,'(a,i4,2x,a20,a,f12.6,a,f12.6)') &
            '# coord ', i, lbl(i), &
            '   C_ii=', C(ii,ii), '   1/C_ii=', 1d0/C(ii,ii)

         ! diagonal entry
         write(iunit,'(2x,i4,2x,a20,2f14.6,a)') &
            i, lbl(i), C(ii,ii), 1d0/C(ii,ii), '  (diagonal)'

         ! collect off-diagonal |C_ij|
         nc_act = 0
         do j = 1, ncoord
            if (j == i) cycle
            nc_act = nc_act + 1
            jsort(nc_act) = j
            aval(nc_act)  = abs(C(ii, idx_ord(j)))
         end do

         ! insertion sort descending
         do p = 2, nc_act
            tmp_r = aval(p);  tmp_i = jsort(p);  q = p-1
            do while (q >= 1 .and. aval(q) < tmp_r)
               aval(q+1) = aval(q);  jsort(q+1) = jsort(q);  q = q-1
            end do
            aval(q+1) = tmp_r;  jsort(q+1) = tmp_i
         end do

         ! write top NCOUP couplings
         do p = 1, min(NCOUP, nc_act)
            j  = jsort(p);  jj = idx_ord(j)
            if (aval(p) < 1d-12) exit
            write(iunit,'(2x,i4,2x,a20,f14.6,a,i4,2x,a20)') &
               j, lbl(j), C(ii,jj), '  <-> ', i, lbl(i)
         end do

      end do

      write(iunit,'(a)') ''
      write(iunit,'(a)') '# end of compliance.dat'
      close(iunit)

      write(*,*)
      write(*,'(a,i4,a)') &
         ' compliance matrix written to compliance.dat  (', ncoord, ' coordinates)'

   end subroutine write_compliance_dat


!-----------------------------------------------------------------------------
! BUILD_LABELS  -- human-readable coordinate labels for compliance.dat
!-----------------------------------------------------------------------------

   subroutine build_labels(n, at, na, nb, nc, ctype, qoff, ncoord, idx_ord, lbl)
      use zmat_bmatrix
      integer,      intent(in)  :: n, ncoord
      integer,      intent(in)  :: at(n), na(n), nb(n), nc(n), ctype(n), qoff(n)
      integer,      intent(in)  :: idx_ord(ncoord)
      character(20),intent(out) :: lbl(ncoord)
      integer          :: i, k
      character(80)    :: buf
      character(2)     :: elsym
      external           elsym

      k = 0;  lbl = '??'

      do i = 2, n                                      ! bonds
         k = k + 1
         write(buf,'(a,a2,i0,a,a2,i0)') &
            'bond ', elsym(at(i)),i, '-', elsym(at(na(i))),na(i)
         lbl(k) = buf(1:20)
      end do
      do i = 2, n                                      ! angles
         if (ctype(i)/=COORD_ANGLE .and. ctype(i)/=COORD_DIHEDRAL) cycle
         k = k + 1
         write(buf,'(a,a2,i0,a,a2,i0,a,a2,i0)') &
            'ang ', elsym(at(i)),i,'-',elsym(at(na(i))),na(i),'-',elsym(at(nb(i))),nb(i)
         lbl(k) = buf(1:20)
      end do
      do i = 2, n                                      ! dihedrals
         if (ctype(i)/=COORD_DIHEDRAL) cycle
         k = k + 1
         write(buf,'(a,a2,i0,a,a2,i0,a,a2,i0,a,a2,i0)') &
            'dih ', elsym(at(i)),i,'-',elsym(at(na(i))),na(i),'-', &
            elsym(at(nb(i))),nb(i),'-',elsym(at(nc(i))),nc(i)
         lbl(k) = buf(1:20)
      end do
      do i = 2, n                                      ! linear bends
         if (ctype(i)/=COORD_LINBEND) cycle
         k = k + 1
         write(buf,'(a,a2,i0,a,a2,i0,a,a2,i0)') &
            'lb1 ',elsym(at(nb(i))),nb(i),'-',elsym(at(na(i))),na(i),'-',elsym(at(i)),i
         lbl(k) = buf(1:20)
         k = k + 1
         write(buf,'(a,a2,i0,a,a2,i0,a,a2,i0)') &
            'lb2 ',elsym(at(nb(i))),nb(i),'-',elsym(at(na(i))),na(i),'-',elsym(at(i)),i
         lbl(k) = buf(1:20)
      end do

   end subroutine build_labels


!-----------------------------------------------------------------------------
! ATMASS  -- standard atomic weight in amu for Z=1..86
!   Values: IUPAC 2021 abridged standard atomic weights.
!-----------------------------------------------------------------------------

   pure function atmass(iz) result(m)
      integer, intent(in) :: iz
      real(8) :: m
      real(8), parameter :: mass(86) = [ &
         1.008d0,   4.003d0,   6.941d0,   9.012d0,  10.811d0, &  !  1  H-B
        12.011d0,  14.007d0,  15.999d0,  18.998d0,  20.180d0, &  !  6  C-Ne
        22.990d0,  24.305d0,  26.982d0,  28.086d0,  30.974d0, &  ! 11  Na-P
        32.065d0,  35.453d0,  39.948d0,  39.098d0,  40.078d0, &  ! 16  S-Ca
        44.956d0,  47.867d0,  50.942d0,  51.996d0,  54.938d0, &  ! 21  Sc-Mn
        55.845d0,  58.933d0,  58.693d0,  63.546d0,  65.380d0, &  ! 26  Fe-Zn
        69.723d0,  72.630d0,  74.922d0,  78.971d0,  79.904d0, &  ! 31  Ga-Br
        83.798d0,  85.468d0,  87.620d0,  88.906d0,  91.224d0, &  ! 36  Kr-Zr
        92.906d0,  95.960d0,  98.000d0, 101.070d0, 102.906d0, &  ! 41  Nb-Rh
       106.420d0, 107.868d0, 112.411d0, 114.818d0, 118.710d0, &  ! 46  Pd-Sn
       121.760d0, 127.600d0, 126.904d0, 131.293d0, 132.905d0, &  ! 51  Sb-Cs
       137.327d0, 138.905d0, 140.116d0, 140.908d0, 144.242d0, &  ! 56  Ba-Nd
       145.000d0, 150.360d0, 151.964d0, 157.250d0, 158.925d0, &  ! 61  Pm-Tb
       162.500d0, 164.930d0, 167.259d0, 168.934d0, 173.054d0, &  ! 66  Dy-Yb
       174.967d0, 178.490d0, 180.948d0, 183.840d0, 186.207d0, &  ! 71  Lu-Re
       190.230d0, 192.217d0, 195.084d0, 196.967d0, 200.590d0, &  ! 76  Os-Hg
       204.383d0, 207.200d0, 208.980d0, 209.000d0, 210.000d0, &  ! 81  Tl-At
       222.000d0 ]                                                  ! 86  Rn
      if (iz >= 1 .and. iz <= 86) then
         m = mass(iz)
      else
         m = dble(iz) * 2.0d0   ! crude fallback for heavy elements
      end if
   end function atmass


!-----------------------------------------------------------------------------
! ELSYM  -- 2-char element symbol for atomic number Z=1..86
!-----------------------------------------------------------------------------

   function elsym(iz) result(s)
      integer, intent(in) :: iz
      character(2) :: s
      character(2) :: sym(86)
      data sym / &
         'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
         'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
         'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
         'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
         'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
         'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
         'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
         'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
         'Tl','Pb','Bi','Po','At','Rn' /
      if (iz >= 1 .and. iz <= 86) then
         s = sym(iz)
      else
         write(s,'(i2)') iz
      end if
   end function elsym
