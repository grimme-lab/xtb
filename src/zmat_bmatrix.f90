module xtb_zmat_bmatrix
   use xtb_mctc_accuracy, only : wp
   !-------------------------------------------------------------------------
   ! Z-matrix internal coordinates and numerical Wilson B-matrix.
   !
   ! Handles all molecular topologies:
   !   - diatomics
   !   - globally linear molecules        (e.g. CO2, HCN, HCCH)
   !   - molecules with linear functional groups (e.g. Ni(CO)4, metal carbonyls)
   !   - standard non-linear molecules
   !
   ! Key idea: per-atom coordinate type stored in ctype(i):
   !   BOND     (1): bond i--NA(i)                          contributes 1 DOF
   !   ANGLE    (2): valence angle i-NA(i)-NB(i)            contributes 1 DOF
   !   DIHEDRAL (3): dihedral i-NA-NB-NC                    contributes 1 DOF
   !   LINBEND  (4): 2 Decius linear bends at NA(i)         contributes 2 DOF
   !                 replaces ANGLE+DIHEDRAL for linear centres
   !
   ! The total number of internal coordinates nint is therefore not simply
   ! 3N-6 or 3N-5 but depends on the actual molecular topology:
   !   nint = (N-1) bonds + N_angle angles + N_dihedral dihedrals
   !          + 2*N_linbend linear bending pairs
   ! and is determined automatically by setup_zmat.
   !
   ! Connectivity is determined from covalent radii (Alvarez, Dalton Trans.
   ! 2008) with a tolerance factor of 1.3.  A BFS traversal of the bond
   ! graph defines the Z-matrix ordering and the NA parent array, ensuring
   ! chemically correct reference chains even for metal complexes where a
   ! nearest-neighbour heuristic fails.  NB and NC are chosen to maximise
   ! the reference angle (> 20 deg threshold) to avoid near-singular
   ! B-matrix rows.
   !
   ! LINBEND detection scans all covalent neighbours of NA(i): if any
   ! i--NA(i)--j angle exceeds 175 deg, atom i is assigned LINBEND.
   ! This correctly identifies linear centres in functional groups
   ! (e.g. M-C=O) without requiring a globally linear molecule.
   ! For each LINBEND centre a fixed orthonormal frame (e1,e2) perpendicular
   ! to the bond axis is built ONCE at the reference geometry and stored in
   ! e1f/e2f.  This frame MUST NOT be recomputed at perturbed geometries or
   ! the finite-difference B-matrix becomes meaningless.
   !
   ! Coordinate vector layout (flat vector q, length nint):
   !   For each atom i=2..n, offset qoff(i) (0-based):
   !     always:          q(qoff(i)+1) = bond length i--NA(i)   [Angstrom]
   !     if ANGLE:        q(qoff(i)+2) = valence angle           [rad]
   !     if DIHEDRAL:     q(qoff(i)+2) = valence angle           [rad]
   !                      q(qoff(i)+3) = dihedral                [rad, 0..2pi)
   !     if LINBEND:      q(qoff(i)+2) = beta_1 (in-plane  Decius bend) [rad]
   !                      q(qoff(i)+3) = beta_2 (out-of-plane Decius bend) [rad]
   !
   ! The B-matrix is computed by central finite differences with a default
   ! step of 1e-3 Angstrom (adjustable via optional argument).  Dihedral
   ! differences are wrapped to (-pi, pi] to handle the 0/2pi discontinuity.
   !
   ! Units throughout: bond lengths in Angstrom, angles in radians.
   ! Compliance constants from the companion xtb_compliance module are
   ! in Bohr^2/Hartree (= a0^2/Eh).  The relaxed force constant 1/C_ii
   ! is in Eh/a0^2; conversion to N/cm: 1 Eh/a0^2 = 15.5693 N/cm.
   !
   ! Public API:
   !   COORD_BOND, COORD_ANGLE, COORD_DIHEDRAL, COORD_LINBEND  (parameters)
   !   setup_zmat(xyz, n, at, na, nb, nc, ctype, e1f, e2f, qoff, nint)
   !   cartesian_to_int(xyz, n, na, nb, nc, ctype, e1f, e2f, qoff, nint, q)
   !   compute_bmatrix(xyz, n, na, nb, nc, ctype, e1f, e2f, qoff, nint, B, step)
   !   is_linear(xyz, n)   -- global linearity check (utility)
   !-------------------------------------------------------------------------
   implicit none
   private

   real(wp), parameter :: pi       = 4.0d0*atan(1.0d0)
   real(wp), parameter :: two_pi   = 2.0d0*pi
   real(wp), parameter :: DEGREE   = 1.0d0     ! keep radians
   ! Threshold for detecting a linear centre: cos(theta) > cos(lin_tol)
   ! i.e. theta > 180 - lin_tol_deg.  Use ~5 degrees margin.
   real(wp), parameter :: lin_tol  = 0.08726646d0  ! 5 degrees in radians

   integer, parameter, public :: COORD_BOND     = 1
   integer, parameter, public :: COORD_ANGLE    = 2
   integer, parameter, public :: COORD_DIHEDRAL = 3
   integer, parameter, public :: COORD_LINBEND  = 4

   public :: is_linear
   public :: setup_zmat
   public :: cartesian_to_int
   public :: compute_bmatrix

contains

   !==========================================================================
   ! IS_LINEAR: global linearity check (utility; used by run_compl for info)
   !==========================================================================
   function is_linear(xyz, n) result(lin)
      integer, intent(in) :: n
      real(wp), intent(in) :: xyz(3,n)
      logical :: lin
      real(wp) :: cx,cy,cz,dx,dy,dz
      real(wp) :: Ixx,Iyy,Izz,Ixy,Ixz,Iyz,Itens(3,3),evals(3),work(20)
      integer :: i, info
      if (n <= 2) then; lin = .true.; return; end if
      cx=sum(xyz(1,1:n))/n; cy=sum(xyz(2,1:n))/n; cz=sum(xyz(3,1:n))/n
      Ixx=0d0;Iyy=0d0;Izz=0d0;Ixy=0d0;Ixz=0d0;Iyz=0d0
      do i=1,n
         dx=xyz(1,i)-cx; dy=xyz(2,i)-cy; dz=xyz(3,i)-cz
         Ixx=Ixx+dy*dy+dz*dz; Iyy=Iyy+dx*dx+dz*dz; Izz=Izz+dx*dx+dy*dy
         Ixy=Ixy-dx*dy; Ixz=Ixz-dx*dz; Iyz=Iyz-dy*dz
      end do
      Itens(1,:)=[Ixx,Ixy,Ixz]; Itens(2,:)=[Ixy,Iyy,Iyz]; Itens(3,:)=[Ixz,Iyz,Izz]
      call dsyev('N','U',3,Itens,3,evals,work,20,info)
      lin = (evals(1) < 0.01d0*(abs(evals(3))+1d-30))
   end function is_linear


   !==========================================================================
   ! SETUP_ZMAT
   !   Determine connectivity (NA,NB,NC), per-atom coordinate type (ctype),
   !   the fixed perpendicular frames (e1f,e2f) for LINBEND centres, the
   !   offset array (qoff) into the flat coordinate vector, and nint.
   !
   !   at(n): atomic numbers (Z).  Used for covalent-radius bond detection.
   !
   !   Connectivity algorithm:
   !   Step 1a: build bond table using covalent radii (Alvarez 2008).
   !            Bond i-j if dist(i,j) < 1.3*(rcov(i)+rcov(j))  [Bohr].
   !   Step 1b: build Z-matrix chain from bond table.
   !            NA(i) = nearest bonded atom already placed (1..i-1).
   !            Fallback to geometrically nearest if no bonded partner placed.
   !            NB(i) = NA(NA(i)),  NC(i) = NB(NA(i)).
   !==========================================================================
   subroutine setup_zmat(xyz, n, at, na, nb, nc, ctype, e1f, e2f, qoff, nint)
      integer, intent(in)  :: n
      real(wp), intent(in)  :: xyz(3,n)
      integer, intent(in)  :: at(n)
      integer, intent(out) :: na(n), nb(n), nc(n)
      integer, intent(out) :: ctype(n)
      real(wp), intent(out) :: e1f(3,n), e2f(3,n)
      integer, intent(out) :: qoff(n)
      integer, intent(out) :: nint

      logical :: bonded(n,n), inqueue(n), is_linbend
      integer :: i, j, k, q, order(n), parent_arr(n), head, tail, placed
      real(wp) :: r, rmin, angl, u(3), thr
      real(wp), parameter :: bond_f    = 1.3d0
      real(wp), parameter :: min_angle = 0.34906585d0  ! 20 deg

      ! ------------------------------------------------------------------
      ! Step 1a: bond table from covalent radii
      ! ------------------------------------------------------------------
      bonded = .false.
      do i = 1, n
         do j = i+1, n
            thr = bond_f * (rcov(at(i)) + rcov(at(j)))
            r = sqrt( (xyz(1,i)-xyz(1,j))**2 &
                    + (xyz(2,i)-xyz(2,j))**2 &
                    + (xyz(3,i)-xyz(3,j))**2 )
            if (r < thr) then
               bonded(i,j) = .true.
               bonded(j,i) = .true.
            end if
         end do
      end do

      ! ------------------------------------------------------------------
      ! Step 1b: BFS on bond graph to build Z-matrix ordering + NA array.
      !   order(p) = atom at Z-matrix position p; parent(i) = NA(i).
      !   BFS guarantees every atom connects to an already-listed atom.
      !   Disconnected atoms (metals with no covalent bonds) appended last.
      ! ------------------------------------------------------------------
      order=0; parent_arr=0; inqueue=.false.
      head=1; tail=1; placed=0
      order(1)=1; inqueue(1)=.true.; placed=1; tail=2

      bfs_loop: do while (head < tail)
         i = order(head); head = head+1
         do j = 1, n
            if (inqueue(j)) cycle
            if (.not. bonded(i,j)) cycle
            order(tail)=j; parent_arr(j)=i
            inqueue(j)=.true.; tail=tail+1; placed=placed+1
         end do
      end do bfs_loop

      ! Append unreached atoms (fragments / isolated metals)
      do i = 1, n
         if (inqueue(i)) cycle
         rmin=1d30; k=order(1)
         do j=1,placed
            r=(xyz(1,i)-xyz(1,order(j)))**2+(xyz(2,i)-xyz(2,order(j)))**2 &
             +(xyz(3,i)-xyz(3,order(j)))**2
            if (r<rmin) then; rmin=r; k=order(j); end if
         end do
         order(tail)=i; parent_arr(i)=k
         inqueue(i)=.true.; tail=tail+1; placed=placed+1
      end do

      na=0
      do i=2,n; na(i)=parent_arr(i); end do

      ! ------------------------------------------------------------------
      ! NB(i): best angle-reference atom (not NA(i), not i itself).
      !   Start with NA(NA(i)); validate angle > min_angle.
      !   Fall back to best_nb() search if needed.
      ! NC(i): similarly from NA(NB(i)).
      ! ------------------------------------------------------------------
      nb=0; nc=0

      do i=3,n
         k=0
         if (na(i)>0) k=na(na(i))
         if (k>0 .and. k/=i) then
            angl=bangle(xyz,i,na(i),k)
            if (angl<min_angle .or. pi-angl<lin_tol*0.5d0) k=0
         else
            k=0
         end if
         if (k==0) k=best_nb(xyz,i,na(i),n,min_angle)
         nb(i)=k
      end do

      do i=4,n
         k=0
         if (nb(i)>0) k=na(nb(i))
         if (k>0 .and. k/=i .and. k/=na(i) .and. k/=nb(i)) then
            ! Validate NC: the dihedral i-NA-NB-NC is only well-defined if
            ! NC is not collinear with the NA-NB axis, i.e. the angle
            ! NC--NB--NA must not be near 0 or 180 degrees.
            angl=bangle(xyz,k,nb(i),na(i))
            if (angl<min_angle .or. pi-angl<lin_tol) k=0
         else
            k=0
         end if
         if (k==0) then
            ! Search ALL atoms (not just j<i) for a valid NC: j not collinear
            ! with the NA-NB axis.  Prefer nearest atom with good angle.
            rmin=1d30
            do j=1,n
               if (j==i.or.j==na(i).or.j==nb(i)) cycle
               angl=bangle(xyz,j,nb(i),na(i))
               if (angl<min_angle .or. pi-angl<lin_tol) cycle
               r=(xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2
               if (r<rmin) then; rmin=r; k=j; end if
            end do
            ! Last-resort: relax angle filter, take nearest non-degenerate atom
            if (k==0) then
               rmin=1d30
               do j=1,n
                  if (j==i.or.j==na(i).or.j==nb(i)) cycle
                  angl=bangle(xyz,j,nb(i),na(i))
                  if (pi-angl<lin_tol*0.5d0) cycle   ! exclude exactly linear
                  r=(xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2
                  if (r<rmin) then; rmin=r; k=j; end if
               end do
            end if
         end if
         nc(i)=k
      end do

      na(1)=0; nb(1)=0; nc(1)=0
                nb(2)=0; nc(2)=0
      if (n>=3)          nc(3)=0
      ! ------------------------------------------------------------------
      ! Step 2: coordinate types and perpendicular frames.
      !
      ! LINBEND detection: atom i gets LINBEND if the chemical angle
      ! i -- na(i) -- X is near 180° for X = any bonded neighbour of
      ! na(i) other than i itself.  This is independent of nb(i) (which
      ! was chosen to avoid near-linear triples).
      ! For globally linear centres all bonded pairs of na(i) are ~180°.
      ! ------------------------------------------------------------------
      e1f=0d0; e2f=0d0
      ctype(1) = 0
      do i = 2, n
         ctype(i) = COORD_BOND
         if (i < 3) cycle

         ! Check if na(i) is a linear centre w.r.t. atom i
         ! by scanning bonded neighbours of na(i)
         is_linbend = .false.
         do j = 1, n
            if (j == i) cycle
            if (.not. bonded(na(i), j)) cycle
            angl = bangle(xyz, i, na(i), j)
            if (pi - angl < lin_tol) then
               is_linbend = .true.
               exit
            end if
         end do

         if (is_linbend) then
            ctype(i) = COORD_LINBEND
            u = xyz(:,i) - xyz(:,na(i))
            u = u / sqrt(dot_product(u,u))
            call perp_frame(u, e1f(:,i), e2f(:,i))
         else if (i < 4) then
            ctype(i) = COORD_ANGLE
         else
            ctype(i) = COORD_DIHEDRAL
         end if
      end do

      ! ------------------------------------------------------------------
      ! Step 3: offset array and nint
      ! ------------------------------------------------------------------
      q = 0; qoff(1) = 0
      do i = 2, n
         qoff(i) = q
         select case (ctype(i))
            case (COORD_BOND)     ; q = q + 1
            case (COORD_ANGLE)    ; q = q + 2
            case (COORD_DIHEDRAL) ; q = q + 3
            case (COORD_LINBEND)  ; q = q + 3
         end select
      end do
      nint = q

   end subroutine setup_zmat


   !==========================================================================
   ! RCOV  -- covalent radius in Bohr for atomic number Z (Alvarez 2008)
   !==========================================================================
   pure function rcov(iz) result(r)
      integer, intent(in) :: iz
      real(wp) :: r
      real(wp), parameter :: ang2bohr = 1.0d0/0.529177210903d0
      ! Covalent radii in Angstrom, Z=1..86
      real(wp), parameter :: rc(86) = [ &
         0.31d0,0.28d0,1.28d0,0.96d0,0.84d0,0.73d0,0.71d0,0.66d0, &  !  1  H-O
         0.57d0,0.58d0,1.66d0,1.41d0,1.21d0,1.11d0,1.07d0,1.05d0, &  !  9  F-S
         1.02d0,1.06d0,2.03d0,1.76d0,1.70d0,1.60d0,1.53d0,1.39d0, &  ! 17  Cl-Cr
         1.61d0,1.52d0,1.50d0,1.24d0,1.32d0,1.22d0,1.22d0,1.20d0, &  ! 25  Mn-Ge
         1.19d0,1.20d0,1.20d0,1.16d0,2.20d0,1.95d0,1.90d0,1.75d0, &  ! 33  As-Zr
         1.64d0,1.54d0,1.47d0,1.46d0,1.42d0,1.39d0,1.45d0,1.44d0, &  ! 41  Nb-Cd
         1.42d0,1.39d0,1.39d0,1.38d0,1.39d0,1.40d0,2.44d0,2.15d0, &  ! 49  In-Ba
         2.07d0,2.04d0,2.03d0,2.01d0,1.99d0,1.98d0,1.98d0,1.96d0, &  ! 57  La-Gd
         1.94d0,1.92d0,1.92d0,1.89d0,1.90d0,1.87d0,1.87d0,1.75d0, &  ! 65  Tb-Hf
         1.70d0,1.62d0,1.51d0,1.44d0,1.41d0,1.36d0,1.36d0,1.32d0, &  ! 73  Ta-Hg
         1.45d0,1.46d0,1.48d0,1.40d0,1.50d0,1.50d0 ]                  ! 81  Tl-Rn
      if (iz >= 1 .and. iz <= 86) then
         r = rc(iz) * ang2bohr
      else
         r = 2.0d0
      end if
   end function rcov


   !==========================================================================
   ! CARTESIAN_TO_INT
   !   Evaluate the flat internal coordinate vector q(1:nint).
   !   The perpendicular frames e1f/e2f must be the ones from setup_zmat
   !   (i.e. fixed at the reference geometry).
   !==========================================================================
   subroutine cartesian_to_int(xyz, n, na, nb, nc, ctype, e1f, e2f, qoff, nint, q)
      integer, intent(in)  :: n, nint
      real(wp), intent(in)  :: xyz(3,n)
      integer, intent(in)  :: na(n), nb(n), nc(n), ctype(n), qoff(n)
      real(wp), intent(in)  :: e1f(3,n), e2f(3,n)
      real(wp), intent(out) :: q(nint)

      integer :: i, j, k, l, l_used, o
      real(wp) :: angl, tol, vi(3), vk(3), ri, rk

      q = 0d0

      do i = 2, n
         j = na(i); k = nb(i); l = nc(i)
         o = qoff(i)   ! 0-based offset => q(o+1) is first coord of atom i

         ! Bond (always)
         q(o+1) = dist(xyz, i, j)

         select case (ctype(i))

         case (COORD_BOND)
            ! nothing more

         case (COORD_ANGLE)
            q(o+2) = bangle(xyz, i, j, k) * DEGREE

         case (COORD_DIHEDRAL)
            q(o+2) = bangle(xyz, i, j, k) * DEGREE
            ! Check J-K-L angle for dihedral reference
            tol = 0.2617994d0
            call bangle_val(xyz, j, k, l, angl)
            l_used = l
            if (angl > pi-tol .or. angl < tol) &
               l_used = find_dihedral_atom(xyz, i, j, k, l, tol)
            q(o+3) = dihed(xyz, i, j, k, l_used) * DEGREE

         case (COORD_LINBEND)
            ! Two Decius linear bending coordinates at centre j=NA(i)
            ! vi = unit vector j->i,  vk = unit vector j->k=NB(i)
            vi = xyz(:,i) - xyz(:,j)
            ri = sqrt(dot_product(vi,vi)); vi = vi/ri
            vk = xyz(:,k) - xyz(:,j)
            rk = sqrt(dot_product(vk,vk)); vk = vk/rk
            q(o+2) = dot_product(e1f(:,i), vi) + dot_product(e1f(:,i), vk)
            q(o+3) = dot_product(e2f(:,i), vi) + dot_product(e2f(:,i), vk)

         end select
      end do

   end subroutine cartesian_to_int


   !==========================================================================
   ! COMPUTE_BMATRIX
   !   Numerical central-difference Wilson B-matrix, shape (nint, 3n).
   !   The perpendicular frames are held fixed (from setup_zmat) so that
   !   LINBEND finite differences are well-defined.
   !==========================================================================
   subroutine compute_bmatrix(xyz, n, na, nb, nc, ctype, e1f, e2f, qoff, nint, B, step_in)
      integer, intent(in)           :: n, nint
      real(wp), intent(in)           :: xyz(3,n)
      integer, intent(in)           :: na(n), nb(n), nc(n), ctype(n), qoff(n)
      real(wp), intent(in)           :: e1f(3,n), e2f(3,n)
      real(wp), intent(out)          :: B(:,:)
      real(wp), intent(in), optional :: step_in

      integer :: i, ii, iii
      real(wp) :: step
      real(wp) :: xyz_w(3,n), qp(nint), qm(nint)

      step  = 1d-3
      if (present(step_in)) step = step_in

      B     = 0d0
      xyz_w = xyz
      iii   = 0

      do i = 1, n
         do ii = 1, 3
            iii = iii + 1

            xyz_w(ii,i) = xyz_w(ii,i) + step
            call cartesian_to_int(xyz_w,n,na,nb,nc,ctype,e1f,e2f,qoff,nint,qp)

            xyz_w(ii,i) = xyz_w(ii,i) - 2d0*step
            call cartesian_to_int(xyz_w,n,na,nb,nc,ctype,e1f,e2f,qoff,nint,qm)

            xyz_w(ii,i) = xyz_w(ii,i) + step

            ! Raw finite differences
            B(:,iii) = (qp - qm) / (2d0*step)

            ! (dihedral wrapping done in wrap_dihedral_rows after loop)

         end do
      end do

      ! Redo dihedral wrapping cleanly (the inline above has a units error)
      ! Recompute properly: dihedral column wrapping
      call wrap_dihedral_rows(B, n, nint, ctype, qoff, step)

   end subroutine compute_bmatrix


   !==========================================================================
   ! PRIVATE helpers
   !==========================================================================

   !--------------------------------------------------------------------------
   ! Wrap dihedral rows of B to correct for 2pi discontinuity.
   ! The raw finite difference (qp-qm)/(2h) is wrong if qp-qm ~ 2pi.
   ! We recompute: dq/(2h) where dq is wrapped to (-pi,pi].
   ! This needs the original qp,qm — but we only have B=dq/(2h).
   ! Equivalently: if |B(row,col)| * 2h > pi, wrap B by +-pi/h.
   !--------------------------------------------------------------------------
   subroutine wrap_dihedral_rows(B, n, nint, ctype, qoff, step)
      integer, intent(in)    :: n, nint
      real(wp), intent(inout) :: B(nint,3*n)
      integer, intent(in)    :: ctype(n), qoff(n)
      real(wp), intent(in)    :: step
      integer :: i, col, row
      real(wp) :: dq

      do i = 2, n
         if (ctype(i) /= COORD_DIHEDRAL) cycle
         row = qoff(i) + 3   ! dihedral is 3rd entry for this atom
         do col = 1, 3*n
            dq = B(row,col) * 2d0*step   ! recover original dq
            if (dq >  pi) dq = dq - two_pi
            if (dq < -pi) dq = dq + two_pi
            B(row,col) = dq / (2d0*step)
         end do
      end do

   end subroutine wrap_dihedral_rows


   pure function dist(xyz, i, j) result(r)
      real(wp), intent(in) :: xyz(:,:)
      integer, intent(in) :: i, j
      real(wp) :: r
      r = sqrt((xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2)
   end function dist

   pure function bangle(xyz, i, j, k) result(angle)
      real(wp), intent(in) :: xyz(:,:)
      integer, intent(in) :: i, j, k
      real(wp) :: angle, d2ij, d2jk, d2ik, xy, temp
      d2ij=(xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2
      d2jk=(xyz(1,j)-xyz(1,k))**2+(xyz(2,j)-xyz(2,k))**2+(xyz(3,j)-xyz(3,k))**2
      d2ik=(xyz(1,i)-xyz(1,k))**2+(xyz(2,i)-xyz(2,k))**2+(xyz(3,i)-xyz(3,k))**2
      xy=sqrt(d2ij*d2jk)
      temp=max(-1d0,min(1d0,0.5d0*(d2ij+d2jk-d2ik)/xy))
      angle=acos(temp)
   end function bangle

   pure subroutine bangle_val(xyz, i, j, k, angle)
      real(wp), intent(in)  :: xyz(:,:)
      integer, intent(in)  :: i, j, k
      real(wp), intent(out) :: angle
      angle = bangle(xyz,i,j,k)
   end subroutine bangle_val

   pure function dihed(xyz, i, j, k, l) result(angle)
      real(wp), intent(in) :: xyz(:,:)
      integer, intent(in) :: i, j, k, l
      real(wp) :: angle
      real(wp) :: xi1,xj1,xl1,yi1,yj1,yl1,zi1,zj1,zl1
      real(wp) :: djk,cosa,ddd,yxd,cph,sph,cth,sth
      real(wp) :: xi2,xl2,yi2,yj2,yl2,yi3,yl3
      xi1=xyz(1,i)-xyz(1,k); yi1=xyz(2,i)-xyz(2,k); zi1=xyz(3,i)-xyz(3,k)
      xj1=xyz(1,j)-xyz(1,k); yj1=xyz(2,j)-xyz(2,k); zj1=xyz(3,j)-xyz(3,k)
      xl1=xyz(1,l)-xyz(1,k); yl1=xyz(2,l)-xyz(2,k); zl1=xyz(3,l)-xyz(3,k)
      djk=sqrt(xj1**2+yj1**2+zj1**2)
      cosa=max(-1d0,min(1d0,zj1/djk)); ddd=1d0-cosa**2
      if (ddd<=0d0) then
         xi2=xi1;xl2=xl1;yi2=yi1;yl2=yl1;cth=cosa;sth=0d0
      else
         yxd=djk*sqrt(ddd)
         if (yxd>1d-6) then
            cph=yj1/yxd; sph=xj1/yxd
            xi2=xi1*cph-yi1*sph; xl2=xl1*cph-yl1*sph
            yi2=xi1*sph+yi1*cph; yj2=xj1*sph+yj1*cph
            yl2=xl1*sph+yl1*cph; cth=cosa; sth=yj2/djk
         else
            xi2=xi1;xl2=xl1;yi2=yi1;yl2=yl1;cth=cosa;sth=0d0
         end if
      end if
      yi3=yi2*cth-zi1*sth; yl3=yl2*cth-zl1*sth
      angle=dang_pure(xl2,yl3,xi2,yi3)
      if (angle<0d0)    angle=two_pi+angle
      if (angle>=two_pi) angle=0d0
   end function dihed

   pure function dang_pure(a1,a2,b1,b2) result(r)
      real(wp), intent(in) :: a1,a2,b1,b2
      real(wp) :: r,ra1,ra2,rb1,rb2,an,bn,sn,cn
      real(wp), parameter :: z=1d-6
      r=0d0
      if (abs(a1)<z.and.abs(a2)<z) return
      if (abs(b1)<z.and.abs(b2)<z) return
      an=1d0/sqrt(a1**2+a2**2); bn=1d0/sqrt(b1**2+b2**2)
      ra1=a1*an; ra2=a2*an; rb1=b1*bn; rb2=b2*bn
      sn=ra1*rb2-ra2*rb1; cn=max(-1d0,min(1d0,ra1*rb1+ra2*rb2))
      r=acos(cn)
      if (abs(r)<4d-4) then; r=0d0; return; end if
      if (sn>0d0) r=two_pi-r
      r=-r
   end function dang_pure

   pure function find_dihedral_atom(xyz, ii, j, k, l_def, tol_in) result(lb)
      real(wp), intent(in) :: xyz(:,:)
      integer, intent(in) :: ii, j, k, l_def
      real(wp), intent(in) :: tol_in
      integer :: lb
      integer :: i1
      real(wp) :: r, rmin, angl, tol
      tol=tol_in; rmin=100d0; lb=l_def
      do i1=1,ii-1
         if (i1==j.or.i1==k) cycle
         r=(xyz(1,i1)-xyz(1,k))**2+(xyz(2,i1)-xyz(2,k))**2+(xyz(3,i1)-xyz(3,k))**2
         if (r<rmin) then
            angl=bangle(xyz,j,k,i1)
            if (angl<pi-tol.and.angl>tol) then; rmin=r; lb=i1; end if
         end if
      end do
      if (rmin>99d0) then
         tol=0.087266d0; rmin=100d0
         do i1=1,ii-1
            if (i1==j.or.i1==k) cycle
            r=(xyz(1,i1)-xyz(1,k))**2+(xyz(2,i1)-xyz(2,k))**2+(xyz(3,i1)-xyz(3,k))**2
            if (r<rmin) then
               angl=bangle(xyz,j,k,i1)
               if (angl<pi-tol.and.angl>tol) then; rmin=r; lb=i1; end if
            end if
         end do
      end if
   end function find_dihedral_atom

   !--------------------------------------------------------------------------
   ! BEST_NB: find the atom j in 1..i-1 (j /= na_i, j /= i) that gives the
   ! largest angle j-na_i-i, subject to angle > min_angle and < pi-lin_tol.
   ! Returns 0 if no valid atom found (caller should handle).
   !--------------------------------------------------------------------------
   pure function best_nb(xyz, i, na_i, n, min_angle) result(k)
      real(wp), intent(in) :: xyz(:,:)
      integer, intent(in) :: i, na_i, n
      real(wp), intent(in) :: min_angle
      integer :: k, j
      real(wp) :: angl, best, rmin, r
      k = 0; best = 0d0
      ! First pass: largest angle in (min_angle, pi-lin_tol/2)
      do j = 1, i-1
         if (j == i .or. j == na_i) cycle
         angl = bangle(xyz, i, na_i, j)
         if (angl < min_angle .or. pi - angl < lin_tol*0.5d0) cycle
         if (angl > best) then; best = angl; k = j; end if
      end do
      ! Fallback: geometrically nearest (no angle filter)
      if (k == 0) then
         rmin = 1d30
         do j = 1, i-1
            if (j == i .or. j == na_i) cycle
            r = (xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2
            if (r < rmin) then; rmin=r; k=j; end if
         end do
      end if
   end function best_nb

   pure subroutine perp_frame(u, e1, e2)
      real(wp), intent(in)  :: u(3)
      real(wp), intent(out) :: e1(3), e2(3)
      real(wp) :: ref(3), norm
      if (abs(u(1)) < 0.9d0) then
         ref = [1d0,0d0,0d0]
      else
         ref = [0d0,0d0,1d0]
      end if
      e1 = ref - dot_product(ref,u)*u
      norm = sqrt(dot_product(e1,e1)); e1 = e1/norm
      e2(1)=u(2)*e1(3)-u(3)*e1(2)
      e2(2)=u(3)*e1(1)-u(1)*e1(3)
      e2(3)=u(1)*e1(2)-u(2)*e1(1)
   end subroutine perp_frame

end module xtb_zmat_bmatrix
