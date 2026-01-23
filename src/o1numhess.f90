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

!> O1 numerical Hessian utilities
module xtb_o1numhess
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autoaa
   use xtb_mctc_blas, only : mctc_gemm, mctc_nrm2, mctc_dot
   implicit none
   private

   public :: adj_list
   public :: gen_local_hessian, lr_loop, get_neighbor_list, gen_displdir, swart

   type :: adj_list
      integer, allocatable :: neighbors(:)
   end type adj_list

contains

!> Main routine to recover local Hessian
subroutine gen_local_hessian(ndispl_final, distmat, displdir, g, dmax, hess_out)
   !> Number of displacements
   integer, intent(in) :: ndispl_final
   !> Distance matrix between atoms
   real(wp), intent(in) :: distmat(:, :)
   !> Displacement directions
   real(wp), intent(in) :: displdir(:, :)
   !> Gradient derivatives
   real(wp), intent(in) :: g(:, :)
   !> Maximum distance threshold
   real(wp), intent(in) :: dmax
   !> Output Hessian matrix
   real(wp), intent(out) :: hess_out(:, :)

   real(wp), parameter :: lam = 1.0e-2_wp, bet = 1.5_wp, ddmax = 5.0_wp
   
   ! Local work arrays
   real(wp), allocatable :: W2(:, :), rhs(:, :), rhsv(:), A(:, :), unit_vec(:), f1(:, :), f(:, :), tmp(:, :), tmp2(:, :), tmp3(:, :)
   logical, allocatable :: mask(:, :)
   integer, allocatable :: ipiv(:)
   integer :: i, j, k, l, ndim, N, info

   N = size(distmat, 1)

   ! Calculate Regularization Term W2
   allocate(W2(N, N))
   W2 = lam * max(0.0_wp, distmat(:, :) - dmax)**(2.0_wp * bet)

   ! Calculate rhs
   allocate(rhs(N, N))
   rhs = 0.0_wp
   call mctc_gemm(g(:, :ndispl_final), displdir(:, :ndispl_final), rhs, transb="t")

   ! Masks and Packing
   allocate(mask(N, N))
   mask = (distmat < (dmax + ddmax))
   do i = 2, N
      mask(i, 1:i-1) = .false.
   end do
   
   ! RHS Vector (b in Ax=b)
   rhsv = pack_sym(rhs, mask)
   ndim = size(rhsv)

   ! Compute A
   ! TODO: this needs to become more (RAM) efficient (use iterative solver)
   allocate(A(ndim, ndim), unit_vec(ndim), tmp(N, N), tmp2(N, N), f1(N, N))
   A = 0.0_wp
   do i = 1, ndim
      unit_vec = 0.0_wp
      unit_vec(i) = 1.0_wp
      tmp = unpack_sym(unit_vec, mask, N)
      call mctc_gemm(tmp, displdir(:, :ndispl_final), tmp2)
      call mctc_gemm(tmp2, displdir(:, :ndispl_final), f1, transb="t")
      f1 = (f1 + transpose(f1)) / 2.0_wp
      f = W2 * tmp
      A(:, i) = pack_sym(f1 + f, mask)
   end do

   ! Solve
   allocate(ipiv(ndim))
   call dgesv(ndim, 1, A, ndim, ipiv, rhsv, ndim, info)

   ! Recover Hessian from vector
   hess_out = unpack_sym(rhsv, mask, N)
end subroutine gen_local_hessian

!> Corrects Hessian hnum using a symmetric, low-rank update
!> so that g approx hnum * displdir
subroutine lr_loop(ndispl, g, hess_out, displdir, final_err)
   integer, intent(in) :: ndispl
   real(wp), intent(in) :: g(:, :)    ! Input Gradients
   real(wp), intent(inout) :: hess_out(:, :) ! Hessian to correct
   real(wp), intent(in) :: displdir(:, :)    ! Displacement directions
   real(wp), intent(out) :: final_err ! Final residual error

   real(wp), parameter :: mingrad_LR = 1.0e-3_wp
   real(wp), parameter :: thresh_LR = 1.0e-8_wp
   integer, parameter :: maxiter_LR = 100

   ! Local variables
   real(wp), allocatable :: resid(:, :), hcorr(:, :), tmp(:, :)
   real(wp) :: dampfac, err0, err, norm_g
   integer :: i, j, it, N
   
   ! BLAS Helper
   real(wp), external :: dnrm2

   N = size(g, 1)

   dampfac = 1.0_wp
   err0 = huge(1.0_wp)
   
   norm_g = mctc_nrm2(g(:, :ndispl))

   allocate(hcorr(N, N), tmp(N, N))
   ! 2. Iterative Correction Loop
   loop_lr: do it = 1, maxiter_LR
      
      call mctc_gemm(hess_out, displdir(:, :ndispl), tmp)
      resid = g(:, :ndispl) - tmp

      err = mctc_nrm2(resid)
      
      if (err < thresh_LR) then
            ! Converged successfully
            exit loop_lr
            
      else if (abs(err - err0) < thresh_LR * err0) then
            ! print *, 'Warning: Gradients cannot be reproduced by symmetric Hessian (Stagnation).'
            exit loop_lr
            
      else if (err > err0 .and. err > norm_g) then
            ! Divergence detected
            dampfac = dampfac * 0.5_wp
      end if
      
      ! hcorr = matmul(resid, transpose(displdir(:, :ndispl)))
      call mctc_gemm(resid, displdir(:, :ndispl), hcorr, transb="t")
      hcorr = 0.5_wp * (hcorr + transpose(hcorr))
      hess_out = hess_out + dampfac * hcorr
      
      err0 = err
      
   end do loop_lr

   final_err = err

end subroutine lr_loop

!> Helper to unpack vector to Symmetric Matrix
function unpack_sym(v, mask, n) result(H)
   real(wp), intent(in) :: v(:)
   logical, intent(in) :: mask(n, n)
   integer, intent(in) :: n
   real(wp) :: H(n, n)
   integer :: i, j
   
   H = 0.0_wp
   H = unpack(v, mask, field=0.0_wp)
   
   ! Symmetrize
   do i = 2, n
      H(i, 1:i-1) = H(1:i-1, i)
   end do
end function unpack_sym

function pack_sym(m, mask) result(v)
   real(wp), intent(in) :: m(:, :)
   logical, intent(in) :: mask(:, :)
   real(wp), allocatable :: v(:)
   integer :: i, j, n
   
   ! symmetrize, then pack
   v = pack((m + transpose(m)) * 0.5_wp, mask)
end function pack_sym

subroutine get_neighbor_list(distmat, dmax, nblist)
   real(wp), intent(in) :: distmat(:, :)
   real(wp), intent(in) :: dmax
   type(adj_list), allocatable, intent(out) :: nblist(:)

   integer, allocatable :: labels(:)
   real(wp), allocatable :: comp_dist(:, :)
   integer, allocatable :: mst_matrix(:, :)
   real(wp) :: d, min_d
   real(wp), parameter :: eps = 1.0e-8_wp
   integer :: i, j, ncomp, N

   N = size(distmat, 1)
   allocate(nblist(N))

   ! 1. Calculate Initial Neighbors
   do i = 1, N
      do j = 1, N
            d = distmat(i, j)
            if (d < dmax) then
               call add_neighbor(nblist(i), j)
            end if
      end do
   end do

   ! 2. Identify Connected Components (DFS)
   allocate(labels(N))
   labels = 0
   ncomp = 0
   
   do i = 1, N
      if (labels(i) == 0) then
            ncomp = ncomp + 1
            call dfs_label(i, N, nblist, labels, ncomp)
      end if
   end do

   if (ncomp == 1) return ! Graph is already connected

   ! 3. Distance between components
   !    For every pair of components, find the minimum distance
   allocate(comp_dist(ncomp, ncomp))
   comp_dist = huge(1.0_wp)

   do i = 1, N - 1
      do j = i + 1, N
            if (labels(i) /= labels(j)) then
               d = distmat(i, j)
               if (d < comp_dist(labels(i), labels(j))) then
                  comp_dist(labels(i), labels(j)) = d
                  comp_dist(labels(j), labels(i)) = d
               end if
            end if
      end do
   end do

   ! 4. Minimum Spanning Tree (Prim's Algorithm) on Components
   !    Returns symmetric matrix: 1 if connected in MST, 0 otherwise
   call prim_mst(ncomp, comp_dist, mst_matrix)

   ! 5. Stitching: Add necessary links closer than MST distance + eps
   do i = 1, N - 1
      do j = i + 1, N
            if (labels(i) /= labels(j)) then
               ! If these components are connected in the MST
               if (mst_matrix(labels(i), labels(j)) == 1) then
                  ! Get the min distance required to bridge them
                  min_d = comp_dist(labels(i), labels(j))
                  
                  ! If this pair provides that bridge (handling degeneracy)
                  if (distmat(i, j) <= min_d + eps) then
                        call add_neighbor_unique(nblist(i), j)
                        call add_neighbor_unique(nblist(j), i)
                  end if
               end if
            end if
      end do
   end do
end subroutine get_neighbor_list

! --- Helper: Add to dynamic array ---
subroutine add_neighbor(list, val)
   type(adj_list), intent(inout) :: list
   integer, intent(in) :: val
   integer, allocatable :: tmp(:)
   integer :: sz

   if (.not. allocated(list%neighbors)) then
      allocate(list%neighbors(1))
      list%neighbors(1) = val
   else
      sz = size(list%neighbors)
      allocate(tmp(sz + 1))
      tmp(1:sz) = list%neighbors
      tmp(sz + 1) = val
      call move_alloc(tmp, list%neighbors)
   end if
end subroutine add_neighbor

! --- Helper: Add only if not present ---
subroutine add_neighbor_unique(list, val)
   type(adj_list), intent(inout) :: list
   integer, intent(in) :: val
   integer :: k
   
   if (allocated(list%neighbors)) then
      do k = 1, size(list%neighbors)
            if (list%neighbors(k) == val) return
      end do
   end if
   call add_neighbor(list, val)
end subroutine add_neighbor_unique

! --- Helper: Recursive DFS for labeling ---
recursive subroutine dfs_label(u, n, nblist, labels, comp_id)
   integer, intent(in) :: u, n, comp_id
   type(adj_list), intent(in) :: nblist(:)
   integer, intent(inout) :: labels(:)
   integer :: k, v

   labels(u) = comp_id
   if (.not. allocated(nblist(u)%neighbors)) return

   do k = 1, size(nblist(u)%neighbors)
      v = nblist(u)%neighbors(k)
      if (labels(v) == 0) then
            call dfs_label(v, n, nblist, labels, comp_id)
      end if
   end do
end subroutine dfs_label

! --- Helper: Prim's Algorithm for MST ---
subroutine prim_mst(nc, dists, adj_mst)
   integer, intent(in) :: nc
   real(wp), intent(in) :: dists(nc, nc)
   integer, allocatable, intent(out) :: adj_mst(:, :)
   
   real(wp) :: min_val, key(nc)
   integer :: parent(nc)
   logical :: mst_set(nc)
   integer :: i, count, u, v

   allocate(adj_mst(nc, nc))
   adj_mst = 0
   
   key = huge(1.0_wp)
   parent = 0
   mst_set = .false.
   
   key(1) = 0.0_wp
   parent(1) = -1

   do count = 1, nc - 1
      ! Pick minimum key vertex not yet including in MST
      min_val = huge(1.0_wp)
      u = -1
      do i = 1, nc
            if (.not. mst_set(i) .and. key(i) < min_val) then
               min_val = key(i)
               u = i
            end if
      end do
      
      if (u == -1) exit
      mst_set(u) = .true.

      ! Update adjacent vertices
      do v = 1, nc
            if (dists(u, v) > 0.0_wp .and. .not. mst_set(v) .and. dists(u, v) < key(v)) then
               parent(v) = u
               key(v) = dists(u, v)
            end if
      end do
   end do

   ! Convert parent array to adjacency matrix
   do i = 2, nc
      if (parent(i) /= -1) then
            adj_mst(i, parent(i)) = 1
            adj_mst(parent(i), i) = 1
      end if
   end do
end subroutine prim_mst

subroutine gen_displdir(n, ndispl0, h0, max_nb, nblist, nbcounts, &
                        eps, eps2, displdir, ndispl_final)
   integer, intent(in) :: n, ndispl0, max_nb
   real(wp), intent(in) :: h0(n,n)
   type(adj_list), intent(in) :: nblist(:)
   integer, intent(in) :: nbcounts(n)           ! Actual number of neighbors per atom
   real(wp), intent(in) :: eps, eps2
   
   real(wp), intent(inout) :: displdir(n, n)
   integer, intent(out) :: ndispl_final

   ! Local variables
   integer :: i, j, k, p, q, nnb, info, n_curr, idx, local_max_ind, locind
   integer :: nb_idx(max_nb)
   real(wp) :: ev(n), coverage(n), locev(max_nb), submat(max_nb, max_nb)
   real(wp) :: projmat(max_nb, n), eye(max_nb, max_nb)
   real(wp) :: vec_subset(max_nb, n)
   real(wp) :: loceigs(max_nb)
   real(wp) :: norm_ev1, norm_ev2, v_norm, d_dot
   integer, allocatable :: iwork(:)
   real(wp), allocatable :: work(:), tmp(:, :)
   logical :: early_break

   ! Initialize
   eye = 0.0_wp
   do i = 1, max_nb
      eye(i, i) = 1.0_wp
   end do

   ! Workspace for LAPACK (allocate generously)
   allocate(work(10*max_nb + 10*n)) 
   allocate(iwork(8*max_nb))

   early_break = .true.
   ndispl_final = ndispl0

   allocate(tmp(max_nb, max_nb))
   ! --- Outer Loop: Generate new directions ---
   do n_curr = ndispl0, n - 1
      ! n_curr: number of existing displacements at this point
      
      ev = 0.0_wp
      coverage = 0.0_wp

      ! --- Inner Loop: Iterate over atoms/DoFs ---
      do j = 1, n
            nnb = nbcounts(j)
            nb_idx(:) = nblist(j)%neighbors

            ! Skip if subspace saturated
            if (nnb <= n_curr) cycle

            ! 1. Extract submatrix H0 (submat)
            do p = 1, nnb
               do q = 1, nnb
                  submat(q, p) = h0(nb_idx(q), nb_idx(p))
               end do
            end do

            ! 2. Local Projection
            ! Form matrix A = displdir[neighbors, 0:n_curr]
            do p = 1, n_curr
               do q = 1, nnb
                  vec_subset(q, p) = displdir(nb_idx(q), p)
               end do
            end do

            if (n_curr > 0) then
               projmat(:nnb, :n_curr) = -orth(vec_subset(:nnb, :n_curr)) ! TODO: why minus?
               call mctc_gemm(projmat(:nnb, :n_curr), projmat(:nnb, :n_curr), tmp(:nnb, :nnb), transb="t")
               projmat(:nnb, :nnb) = eye(:nnb, :nnb) - tmp(:nnb, :nnb)
            else
               projmat(:nnb, :nnb) = eye(:nnb, :nnb)
            end if
            
            ! submat = P * submat * P.T
            call mctc_gemm(submat(:nnb, :nnb), projmat(:nnb, :nnb), tmp(:nnb, :nnb), transb="t")
            call mctc_gemm(projmat(:nnb, :nnb), tmp(:nnb, :nnb), submat(:nnb, :nnb))
            
            ! Symmetrize
            submat = 0.5_wp * (submat + transpose(submat))

            ! 3. Diagonalization
            ! dsyev: computes eigenvalues and eigenvectors in ascending order
            call dsyev('V', 'U', nnb, submat, max_nb, loceigs, work, 10*max_nb, info)
            
            ! Find the index of maximum eigenvalue (first occurrence for ties, like Python's argmax)
            ! dsyev sorts ascending, so normally we'd take the last, but for ties we need first max
            locind = maxloc(loceigs(1:nnb), dim=1)
            locev(1:nnb) = submat(1:nnb, locind)

            ! 4. Patching / Phase fixing
            ! Calculate norms for sign decision
            norm_ev1 = 0.0_wp
            norm_ev2 = 0.0_wp
            do p = 1, nnb
               idx = nb_idx(p)
               norm_ev1 = norm_ev1 + ((coverage(idx)*ev(idx) + locev(p))/(coverage(idx)+1.0_wp))**2
               norm_ev2 = norm_ev2 + ((coverage(idx)*ev(idx) - locev(p))/(coverage(idx)+1.0_wp))**2
            end do
            norm_ev1 = sqrt(norm_ev1)
            norm_ev2 = sqrt(norm_ev2)
            
            ! Apply update
            if (norm_ev1 > norm_ev2 + eps) then
               do p = 1, nnb
                     idx = nb_idx(p)
                     ev(idx) = (coverage(idx)*ev(idx) + locev(p))/(coverage(idx)+1.0_wp)
               end do
            else if (norm_ev1 < norm_ev2 - eps) then
               do p = 1, nnb
                     idx = nb_idx(p)
                     ev(idx) = (coverage(idx)*ev(idx) - locev(p))/(coverage(idx)+1.0_wp)
               end do
            else
               ! Deterministic sign fix based on max element
               locind = maxloc(abs(locev(1:nnb)), dim=1)
               if (locev(locind) > 0.0_wp) then
                  do p = 1, nnb
                        idx = nb_idx(p)
                        ev(idx) = (coverage(idx)*ev(idx) + locev(p))/(coverage(idx)+1.0_wp)
                  end do
               else
                  do p = 1, nnb
                        idx = nb_idx(p)
                        ev(idx) = (coverage(idx)*ev(idx) - locev(p))/(coverage(idx)+1.0_wp)
                  end do
               end if
            end if

            ! Update coverage
            do p = 1, nnb
               coverage(nb_idx(p)) = coverage(nb_idx(p)) + 1.0_wp
            end do
      end do ! End J loop

      ! Project out previous columns from global ev
      do k = 1, n_curr
            d_dot = mctc_dot(ev, displdir(:, k))
            ev = ev - d_dot * displdir(:, k)
      end do
      ! --- Check Norm ---
      v_norm = norm2(ev)
      
      if (v_norm < eps2) exit

      ! Normalize and store
      ev = ev / v_norm
      displdir(:, n_curr + 1) = ev

      ndispl_final = n_curr + 1
   end do
end subroutine gen_displdir

function orth(A, tol_in) result(Q)
   real(wp), intent(in) :: A(:,:)
   real(wp), intent(in), optional :: tol_in
   real(wp), allocatable :: Q(:,:)
   
   real(wp), allocatable :: Acopy(:,:), U(:,:), S(:), VT(:,:), work(:)
   real(wp) :: tol
   integer :: m, n, lwork, info, i, k, rk
   
   m = size(A, 1)
   n = size(A, 2)
   k = min(m, n)
   
   allocate(Acopy(m, n), U(m, k), S(k), VT(k, n), work(1))
   
   Acopy = A
   
   call dgesvd('S', 'N', m, n, Acopy, m, S, U, m, VT, k, work, -1, info)
   lwork = int(work(1))
   deallocate(work)
   allocate(work(lwork))
   
   call dgesvd('S', 'N', m, n, Acopy, m, S, U, m, VT, k, work, lwork, info)
   
   if (present(tol_in)) then
      tol = tol_in
   else
      tol = max(m, n) * S(1) * epsilon(1.0_wp)
   end if
   
   rk = 0
   do i = 1, k
      if (S(i) > tol) rk = rk + 1
   end do
   
   allocate(Q(m, rk))
   Q = U(:, 1:rk)
   
   deallocate(Acopy, U, S, VT, work)
end function orth

!> Calculates a modified Swart model Hessian
subroutine swart(xyz, at, hess_out)
   !> coords
   real(wp), intent(in) :: xyz(:, :)
   !> ordinal numbers
   integer, intent(in) :: at(:)
   !> the full model hessian
   real(wp), intent(inout) :: hess_out(:, :)

   ! covalent radii
   real(wp), parameter :: cov(103) = [&
      & 0.32_wp, 0.46_wp, 1.33_wp, 1.02_wp, 0.85_wp, 0.75_wp, 0.71_wp, 0.63_wp, 0.64_wp, 0.67_wp, &
      & 1.55_wp, 1.39_wp, 1.26_wp, 1.16_wp, 1.11_wp, 1.03_wp, 0.99_wp, 0.96_wp, 1.96_wp, 1.71_wp, 1.48_wp, &
      & 1.36_wp, 1.34_wp, 1.22_wp, 1.19_wp, 1.16_wp, 1.10_wp, 1.11_wp, 1.12_wp, 1.18_wp, 1.24_wp, 1.21_wp, &
      & 1.21_wp, 1.16_wp, 1.14_wp, 1.17_wp, 2.10_wp, 1.85_wp, 1.63_wp, 1.54_wp, 1.47_wp, 1.38_wp, 1.28_wp, &
      & 1.25_wp, 1.25_wp, 1.20_wp, 1.28_wp, 1.36_wp, 1.42_wp, 1.40_wp, 1.40_wp, 1.36_wp, 1.33_wp, 1.31_wp, &
      & 2.32_wp, 1.96_wp, 1.80_wp, 1.63_wp, 1.76_wp, 1.74_wp, 1.73_wp, 1.72_wp, 1.68_wp, 1.69_wp, 1.68_wp, &
      & 1.67_wp, 1.66_wp, 1.65_wp, 1.64_wp, 1.70_wp, 1.62_wp, 1.52_wp, 1.46_wp, 1.37_wp, 1.31_wp, 1.29_wp, &
      & 1.22_wp, 1.23_wp, 1.24_wp, 1.33_wp, 1.44_wp, 1.44_wp, 1.51_wp, 1.45_wp, 1.47_wp, 1.42_wp, 2.23_wp, &
      & 2.01_wp, 1.86_wp, 1.75_wp, 1.69_wp, 1.70_wp, 1.71_wp, 1.72_wp, 1.66_wp, 1.66_wp, 1.68_wp, 1.68_wp, &
      & 1.65_wp, 1.67_wp, 1.73_wp, 1.76_wp, 1.61_wp] / 0.529177249_wp
   real(wp), parameter :: wthr = 0.3_wp, f = 0.12_wp, tolth = 0.2_wp, eps1 = wthr**2, eps2 = wthr**2 / exp(1.0_wp)

   real(wp) :: equildist, Hint, bmat6(6), bmat9(9), bmat29(2, 9), outer6(6, 6), outer9(9, 9), s_ijjk, costh, sinth, th1, scalelin
   real(wp), allocatable :: screenfunc(:, :)
   integer :: i, j, k, nat, N, i1, i2, j1, j2, k1, k2
   
   nat = size(xyz, 2)
   N = 3*nat

   hess_out = 0.0_wp

   allocate(screenfunc(nat, nat))
   do i = 1, nat
      do j = i + 1, nat
         equildist = cov(at(i)) + cov(at(j))
         screenfunc(i, j) = exp(1.0_wp - norm2(xyz(:, i) - xyz(:, j)) / equildist)
         screenfunc(j, i) = screenfunc(i, j)
      end do
   end do

   do i = 1, nat
      do j = i + 1, nat
         Hint = 0.35_wp * screenfunc(i, j)**3
         bmat6 = bmat_bond(xyz(:, i) - xyz(:, j))
         outer6 = (spread(bmat6, dim=2, ncopies=6) * spread(bmat6, dim=1, ncopies=6))
         i1 = 3 * i - 2
         i2 = 3 * i
         j1 = 3 * j - 2
         j2 = 3 * j
         hess_out(i1:i2, i1:i2) = hess_out(i1:i2, i1:i2) + Hint * outer6(1:3, 1:3)
         hess_out(i1:i2, j1:j2) = hess_out(i1:i2, j1:j2) + Hint * outer6(1:3, 4:6)
         hess_out(j1:j2, i1:i2) = hess_out(j1:j2, i1:i2) + Hint * outer6(4:6, 1:3)
         hess_out(j1:j2, j1:j2) = hess_out(j1:j2, j1:j2) + Hint * outer6(4:6, 4:6)
      end do
   end do

   do i = 1, nat
      do j = 1, nat
         if (i == j) cycle
         if (screenfunc(i, j) < eps2) cycle
         do k = i + 1, nat
            if (k == j) cycle
            s_ijjk = screenfunc(i, j) * screenfunc(j, k)
            if (s_ijjk < eps1) cycle

            costh = cosangle(xyz(:, i) - xyz(:, j), xyz(:, k) - xyz(:, j))
            sinth = sqrt(max(0.0_wp, 1.0_wp - costh**2))
            Hint = 0.075_wp * s_ijjk**2 * (f + (1 - f) * sinth)**2
            bmat9 = bmat_angle(xyz(:, i) - xyz(:, j), xyz(:, k) - xyz(:, j))

            if (costh > 1.0_wp - tolth) then
               th1 = 1.0_wp - costh
            else
               th1 = 1.0_wp + costh
            end if

            i1 = 3 * i - 2
            i2 = 3 * i
            j1 = 3 * j - 2
            j2 = 3 * j
            k1 = 3 * k - 2
            k2 = 3 * k
            if (th1 < tolth) then
               scalelin = (1.0_wp - (th1 / tolth)**2)**2
               if (costh > 1.0_wp - tolth) then
                  bmat29 = bmat_linangle(xyz(:, i) - xyz(:, j), xyz(:, k) - xyz(:, j))
                  bmat9 = scalelin * bmat29(1, :) + (1.0_wp - scalelin) * bmat9
                  outer9 = Hint * spread(bmat29(2, :), dim=2, ncopies=9) * spread(bmat29(2, :), dim=1, ncopies=9)
                  hess_out(i1:i2, i1:i2) = hess_out(i1:i2, i1:i2) + outer9(1:3, 1:3)
                  hess_out(i1:i2, j1:j2) = hess_out(i1:i2, j1:j2) + outer9(1:3, 4:6)
                  hess_out(i1:i2, k1:k2) = hess_out(i1:i2, k1:k2) + outer9(1:3, 7:9)
                  hess_out(j1:j2, i1:i2) = hess_out(j1:j2, i1:i2) + outer9(4:6, 1:3)
                  hess_out(j1:j2, j1:j2) = hess_out(j1:j2, j1:j2) + outer9(4:6, 4:6)
                  hess_out(j1:j2, k1:k2) = hess_out(j1:j2, k1:k2) + outer9(4:6, 7:9)
                  hess_out(k1:k2, i1:i2) = hess_out(k1:k2, i1:i2) + outer9(7:9, 1:3)
                  hess_out(k1:k2, j1:j2) = hess_out(k1:k2, j1:j2) + outer9(7:9, 4:6)
                  hess_out(k1:k2, k1:k2) = hess_out(k1:k2, k1:k2) + outer9(7:9, 7:9)
               else
                  bmat9 = (1.0_wp - scalelin) * bmat9
               end if
            end if

            outer9 = Hint * spread(bmat9, dim=2, ncopies=9) * spread(bmat9, dim=1, ncopies=9)
            hess_out(i1:i2, i1:i2) = hess_out(i1:i2, i1:i2) + outer9(1:3, 1:3)
            hess_out(i1:i2, j1:j2) = hess_out(i1:i2, j1:j2) + outer9(1:3, 4:6)
            hess_out(i1:i2, k1:k2) = hess_out(i1:i2, k1:k2) + outer9(1:3, 7:9)
            hess_out(j1:j2, i1:i2) = hess_out(j1:j2, i1:i2) + outer9(4:6, 1:3)
            hess_out(j1:j2, j1:j2) = hess_out(j1:j2, j1:j2) + outer9(4:6, 4:6)
            hess_out(j1:j2, k1:k2) = hess_out(j1:j2, k1:k2) + outer9(4:6, 7:9)
            hess_out(k1:k2, i1:i2) = hess_out(k1:k2, i1:i2) + outer9(7:9, 1:3)
            hess_out(k1:k2, j1:j2) = hess_out(k1:k2, j1:j2) + outer9(7:9, 4:6)
            hess_out(k1:k2, k1:k2) = hess_out(k1:k2, k1:k2) + outer9(7:9, 7:9)
         end do
      end do
   end do
end subroutine swart

function bmat_bond(vec) result(bmat)
   real(wp), intent(in) :: vec(3)

   real(wp) :: l, bmat(6)
   
   bmat = 0.0_wp
   l = norm2(vec)

   bmat(1:3) = vec(:) / l
   bmat(4:6) = -vec(:) / l
end function bmat_bond

function bmat_angle(vec1, vec2) result(bmat)
   real(wp), intent(in) :: vec1(3), vec2(3)
   real(wp) :: bmat(9)
   real(wp) :: l1, l2, nvec1(3), nvec2(3)
   real(wp) :: dl(2, 6), dnvec(2, 3, 6), dinprod(9)
   real(wp) :: dot_n1n2
   integer :: ii

   l1 = norm2(vec1)
   l2 = norm2(vec2)
   nvec1 = vec1 / l1
   nvec2 = vec2 / l2

   dl = 0.0_wp
   dl(1, 1:3) = nvec1
   dl(1, 4:6) = -nvec1
   dl(2, 1:3) = nvec2
   dl(2, 4:6) = -nvec2

   dnvec = 0.0_wp
   do ii = 1, 6
      dnvec(1, 1:3, ii) = -nvec1 * dl(1, ii) / l1
      dnvec(2, 1:3, ii) = -nvec2 * dl(2, ii) / l2
   end do
   do ii = 1, 3
      dnvec(1, ii, ii) = dnvec(1, ii, ii) + 1.0_wp/l1
      dnvec(2, ii, ii) = dnvec(2, ii, ii) + 1.0_wp/l2
      dnvec(1, ii, ii+3) = dnvec(1, ii, ii+3) - 1.0_wp/l1
      dnvec(2, ii, ii+3) = dnvec(2, ii, ii+3) - 1.0_wp/l2
   end do

   dinprod = 0.0_wp
   do ii = 1, 3
      dinprod(ii) = mctc_dot(dnvec(1, :, ii), nvec2)
      dinprod(ii+3) = mctc_dot(dnvec(1, :, ii+3), nvec2) + mctc_dot(dnvec(2, :, ii+3), nvec1)
      dinprod(ii+6) = mctc_dot(dnvec(2, :, ii), nvec1)
   end do

   dot_n1n2 = mctc_dot(nvec1, nvec2)
   bmat = -dinprod / sqrt(max(1.0e-15_wp, 1.0_wp - dot_n1n2**2))
end function bmat_angle

function bmat_linangle(vec1, vec2) result(bmat)
   real(wp), intent(in) :: vec1(3), vec2(3)
   real(wp) :: bmat(2,9)
   real(wp) :: l1, l2, nvec1(3), nvec2(3)
   real(wp) :: vn(3), vn2(3), nvn
   real(wp), parameter :: xaxis(3) = [1.0_wp, 0.0_wp, 0.0_wp], yaxis(3) = [0.0_wp, 1.0_wp, 0.0_wp]

   l1 = norm2(vec1)
   l2 = norm2(vec2)
   nvec1 = vec1 / l1
   nvec2 = vec2 / l2

   vn(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
   vn(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
   vn(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
   nvn = norm2(vn)

   if (nvn < 1.0e-15_wp) then
      vn = xaxis - mctc_dot(xaxis, vec1) / l1**2 * vec1
      nvn = norm2(vn)
      if (nvn < 1.0e-15_wp) then
         vn = yaxis - mctc_dot(yaxis, vec1) / l1**2 * vec1
         nvn = norm2(vn)
      end if
   end if
   vn = vn / nvn

   vn2(1) = (vec1(2) - vec2(2)) * vn(3) - (vec1(3) - vec2(3)) * vn(2)
   vn2(2) = (vec1(3) - vec2(3)) * vn(1) - (vec1(1) - vec2(1)) * vn(3)
   vn2(3) = (vec1(1) - vec2(1)) * vn(2) - (vec1(2) - vec2(2)) * vn(1)
   vn2 = vn2 / norm2(vn2)

   bmat = 0.0_wp
   bmat(2, 1:3) = vn / l1
   bmat(2, 7:9) = vn / l2
   bmat(2, 4:6) = -bmat(2, 1:3) - bmat(2, 7:9)
   bmat(1, 1:3) = vn2 / l1
   bmat(1, 7:9) = vn2 / l2
   bmat(1, 4:6) = -bmat(1, 1:3) - bmat(1, 7:9)
end function bmat_linangle

function cosangle(vec1, vec2) result(cos_theta)
    implicit none
    real(wp), intent(in) :: vec1(3), vec2(3)
    real(wp) :: cos_theta
    
    cos_theta = mctc_dot(vec1, vec2) / (norm2(vec1) * norm2(vec2))
end function cosangle
end module xtb_o1numhess
