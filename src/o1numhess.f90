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
   implicit none
   private

   public :: adj_list
   public :: gen_local_hessian, lr_loop, get_neighbor_list, gen_displdir

   type :: adj_list
      integer, allocatable :: neighbors(:)
   end type adj_list

contains

!> Main routine to recover local Hessian
subroutine gen_local_hessian(distmat, displdir, g, dmax, hess_out)
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
   real(wp), allocatable :: W2(:, :), rhs(:, :), rhsv(:), A(:, :), D(:, :), Dprime(:, :), unit_vec(:), f1(:, :), f(:, :), tmp(:, :)
   logical, allocatable :: mask(:, :)
   integer :: i, j, k, l, ndim, ndispl, N, info

   ndispl = size(displdir, 1)
   N = size(distmat, 1)

   ! Calculate Regularization Term W2
   allocate(W2(N, N))
   W2 = lam * max(0.0_wp, distmat(:, :) - dmax)**(2.0_wp * bet)

   ! Calculate rhs
   allocate(rhs(N, N))
   rhs = 0.0_wp
   ! call dgemm('N', 'T', N, N, ndispl, 1.0_wp, g, N, displdir, N, 0.0_wp, rhs, N)
   rhs = matmul(g, transpose(displdir))
   ! rhs = 0.5_wp * (rhs + transpose(rhs))

   ! Masks and Packing
   allocate(mask(N, N))
   mask = (distmat < (dmax + ddmax))
   do i = 1, N
      do j = 1, i - 1
         mask(i, j) = .false.
      end do
   end do
   
   ! RHS Vector (b in Ax=b)
   ! rhsv = pack(transpose(rhs), mask)
   rhsv = pack_sym(rhs, mask)
   ndim = size(rhsv)

   ! Compute A
   allocate(A(ndim, ndim), unit_vec(ndim), tmp(N, N))
   A = 0.0_wp
   do i = 1, ndim
      unit_vec = 0.0_wp
      unit_vec(i) = 1.0_wp
      tmp = unpack_sym(unit_vec, mask, N)
      f1 = matmul(matmul(tmp, displdir), transpose(displdir))
      f1 = (f1 + transpose(f1)) / 2.0_wp
      f = W2 * tmp
      A(:, i) = pack_sym(f1 + f, mask)
   end do

   ! Solve
   call dposv('U', ndim, 1, A, ndim, rhsv, ndim, info)

   ! Recover Hessian from vector
   hess_out = unpack_sym(rhsv, mask, N)
end subroutine gen_local_hessian

!> Corrects Hessian hnum using a symmetric, low-rank update
!> so that g approx hnum * displdir
subroutine lr_loop(ndispl, g_in, hess_out, d_in, final_err)
   integer, intent(in) :: ndispl
   real(wp), intent(in) :: g_in(:, :)    ! Input Gradients
   real(wp), intent(inout) :: hess_out(:, :) ! Hessian to correct
   real(wp), intent(in) :: d_in(:, :)    ! Displacement directions
   real(wp), intent(out) :: final_err ! Final residual error

   real(wp), parameter :: mingrad_LR = 1.0e-3_wp
   real(wp), parameter :: thresh_LR = 1.0e-8_wp
   integer, parameter :: maxiter_LR = 100

   ! Local variables
   real(wp), allocatable :: g_local(:, :), displdir_local(:, :)
   real(wp), allocatable :: resid(:, :), hcorr(:, :)
   real(wp) :: dampfac, err0, err, norm_g, norm_gi, val
   integer :: i, j, it, N
   
   ! BLAS Helper
   real(wp), external :: dnrm2

   N = size(g_in, 1)

   allocate(g_local(N, N), displdir_local(N, N))
   allocate(resid(N, N), hcorr(N, N))

   ! 1. Weighting Step
   g_local = 0.0_wp
   displdir_local = 0.0_wp
   
   do i = 1, ndispl
      ! Python: norm_gi = max(norm(g[:,i]), mingrad)
      norm_gi = max(dnrm2(n, g_in(:, i), 1), mingrad_LR)
      
      ! Scale columns
      g_local(:, i) = (mingrad_LR / norm_gi) * g_in(:, i)
      displdir_local(:, i) = (mingrad_LR / norm_gi) * d_in(:, i)
   end do

   dampfac = 1.0_wp
   err0 = huge(1.0_wp)
   
   ! Calculate Frobenius norm of entire matrix G (treated as vector of size n*ndispl)
   norm_g = dnrm2(N * ndispl, g_local, 1)

   ! 2. Iterative Correction Loop
   loop_lr: do it = 1, maxiter_LR
      
      resid = g_local ! Copy g to resid first
      
      ! Call DGEMM: C = alpha*A*B + beta*C
      ! resid = (-1.0) * hnum * d + (1.0) * resid
      call dgemm('N', 'N', N, ndispl, N, -1.0_wp, hess_out, N, displdir_local, N, 1.0_wp, resid, N)
      
      err = dnrm2(N * N, resid, 1)
      
      if (err < thresh_LR) then
            ! Converged successfully
            exit loop_lr
            
      else if (abs(err - err0) < thresh_LR * err0) then
            print *, 'Warning: Gradients cannot be reproduced by symmetric Hessian (Stagnation).'
            exit loop_lr
            
      else if (err > err0 .and. err > norm_g) then
            ! Divergence detected
            dampfac = dampfac * 0.5_wp
            ! (Optional: Print warning if verbose)
            ! print *, 'Damping factor reduced to', dampfac
      end if
      
      ! (Optional: Print iteration status)
      ! print *, 'Iter', it, 'Error', err

      ! hcorr = resid * d^T
      ! DGEMM: C = alpha*A*B^T + beta*C
      call dgemm('N', 'T', N, N, ndispl, 1.0_wp, resid, N, displdir_local, N, 0.0_wp, hcorr, N)
      
      ! Symmetrize hcorr and Apply Update to hnum
      ! hnum = hnum + dampfac * 0.5 * (hcorr + hcorr^T)
      do j = 1, N
            do i = 1, N
               ! Average off-diagonals
               val = 0.5_wp * (hcorr(i, j) + hcorr(j, i))
               hess_out(i, j) = hess_out(i, j) + (dampfac * val)
            end do
      end do
      
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
   do i = 1, n
      do j = 1, i - 1
         H(i, j) = H(j, i)
      end do
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
   real(8), intent(in) :: dists(nc, nc)
   integer, allocatable, intent(out) :: adj_mst(:, :)
   
   real(8) :: min_val, key(nc)
   integer :: parent(nc)
   logical :: mst_set(nc)
   integer :: i, count, u, v

   allocate(adj_mst(nc, nc))
   adj_mst = 0
   
   key = huge(1.0d0)
   parent = 0
   mst_set = .false.
   
   key(1) = 0.0d0
   parent(1) = -1

   do count = 1, nc - 1
      ! Pick minimum key vertex not yet including in MST
      min_val = huge(1.0d0)
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
            if (dists(u, v) > 0.0d0 .and. .not. mst_set(v) .and. dists(u, v) < key(v)) then
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
                        eps, eps2, displdir, ndispl_final, displdir0)
   integer, intent(in) :: n, ndispl0, max_nb
   real(wp), intent(in) :: h0(n,n)
   type(adj_list), intent(in) :: nblist(:)
   integer, intent(in) :: nbcounts(n)           ! Actual number of neighbors per atom
   real(wp), intent(in) :: eps, eps2
   
   real(wp), allocatable, intent(out) :: displdir(:, :)
   integer, intent(out) :: ndispl_final
   real(wp), intent(in), optional :: displdir0(n, ndispl0)

   ! Local variables
   integer :: i, j, k, p, q, nnb, info, n_curr, idx, local_max_ind, locind
   integer :: nb_idx(max_nb)
   real(wp) :: ev(n), coverage(n), locev(max_nb), submat(max_nb, max_nb)
   real(wp) :: projmat(max_nb, max_nb), eye(max_nb, max_nb)
   real(wp) :: vec_subset(max_nb, n), U(max_nb, max_nb), VT(max_nb, max_nb)
   real(wp) :: S(max_nb), loceigs(max_nb)
   real(wp) :: norm1, norm2, v_norm, d_dot
   integer, allocatable :: iwork(:)
   real(wp), allocatable :: work(:), displdir_tmp(:, :)
   real(wp) :: norm_locev_max
   logical :: early_break

   ! Initialize
   allocate(displdir_tmp(n, n))
   displdir_tmp = 0.0_wp
   if (present(displdir0)) then
      displdir_tmp(:, 1:ndispl0) = displdir0
   end if

   ! Workspace for LAPACK (allocate generously)
   allocate(work(10*max_nb + 10*n)) 
   allocate(iwork(8*max_nb))

   early_break = .true.
   ndispl_final = ndispl0

   ! --- Outer Loop: Generate new directions ---
   ! Corresponds to Python: for i in range(N-Ndispl0)
   do i = 1, n - ndispl0
      
      n_curr = ndispl0 + i - 1 ! Number of existing vectors
      ev = 0.0_wp
      coverage = 0.0_wp

      ! --- Inner Loop: Iterate over atoms/DoFs ---
      do j = 1, n
            nnb = nbcounts(j)
            nb_idx(1:nnb) = nblist(j)%neighbors(1:nnb)

            ! Skip if subspace saturated (heuristic from Python code)
            if (nnb <= n_curr) cycle

            ! 1. Extract submatrix H0 (submat)
            do p = 1, nnb
               do q = 1, nnb
                  submat(q, p) = h0(nb_idx(q), nb_idx(p))
               end do
            end do

            ! 2. Local Projection (orth replacement)
            ! Form matrix A = displdir[neighbors, 0:n_curr]
            do p = 1, n_curr
               do q = 1, nnb
                  vec_subset(q, p) = displdir_tmp(nb_idx(q), p)
               end do
            end do

            ! Perform SVD on vec_subset to find basis: A = U * S * VT
            ! We define rank based on S > eps
            if (n_curr > 0) then
               call dgesdd('S', nnb, n_curr, vec_subset, max_nb, S, U, max_nb, &
                           VT, max_nb, work, -1, iwork, info) ! Query size
               call dgesdd('S', nnb, n_curr, vec_subset, max_nb, S, U, max_nb, &
                           VT, max_nb, work, int(work(1)), iwork, info)
            else
               U = 0.0_wp ! No existing vectors
               S = 0.0_wp
            end if

            ! Construct Projector: P = I - sum(u*u.T) for significant singular values
            eye = 0.0_wp
            do p = 1, nnb; eye(p,p) = 1.0_wp; end do
            
            if (n_curr > 0) then
               do k = 1, min(nnb, n_curr)
                  if (S(k) > 1.0d-13) then ! Numerical threshold for rank
                        do p = 1, nnb
                           do q = 1, nnb
                              eye(q, p) = eye(q, p) - U(q, k) * U(p, k)
                           end do
                        end do
                  end if
               end do
            end if
            
            ! submat = P * submat * P.T
            ! Step A: temp = P * submat
            call dgemm('N', 'N', nnb, nnb, nnb, 1.0_wp, eye, max_nb, submat, max_nb, 0.0_wp, projmat, max_nb)
            ! Step B: submat = temp * P.T
            call dgemm('N', 'T', nnb, nnb, nnb, 1.0_wp, projmat, max_nb, eye, max_nb, 0.0_wp, submat, max_nb)
            
            ! Symmetrize
            do p = 1, nnb
               do q = 1, nnb
                     submat(q,p) = 0.5_wp * (submat(q,p) + submat(p,q))
               end do
            end do

            ! 3. Diagonalization (Eigen decomposition)
            ! dsyev: computes eigenvalues and eigenvectors in ascending order
            call dsyev('V', 'U', nnb, submat, max_nb, loceigs, work, 10*max_nb, info)
            
            ! Find the index of maximum eigenvalue (first occurrence for ties, like Python's argmax)
            ! dsyev sorts ascending, so normally we'd take the last, but for ties we need first max
            locind = maxloc(loceigs(1:nnb), dim=1)
            locev(1:nnb) = submat(1:nnb, locind)

            ! 4. Patching / Phase fixing
            ! First, find index of maximum absolute value in locev (for deterministic sign fix)
            local_max_ind = 1
            norm_locev_max = abs(locev(1))
            do p = 2, nnb
               if (abs(locev(p)) > norm_locev_max) then
                  norm_locev_max = abs(locev(p))
                  local_max_ind = p
               end if
            end do
            
            ! Calculate norms for sign decision
            norm1 = 0.0_wp
            norm2 = 0.0_wp
            do p = 1, nnb
               idx = nb_idx(p)
               norm1 = norm1 + ((coverage(idx)*ev(idx) + locev(p))/(coverage(idx)+1.0_wp))**2
               norm2 = norm2 + ((coverage(idx)*ev(idx) - locev(p))/(coverage(idx)+1.0_wp))**2
            end do
            norm1 = sqrt(norm1)
            norm2 = sqrt(norm2)
            
            ! Apply update
            if (norm1 > norm2 + eps) then
               do p = 1, nnb
                     idx = nb_idx(p)
                     ev(idx) = (coverage(idx)*ev(idx) + locev(p))/(coverage(idx)+1.0_wp)
               end do
            else if (norm1 < norm2 - eps) then
               do p = 1, nnb
                     idx = nb_idx(p)
                     ev(idx) = (coverage(idx)*ev(idx) - locev(p))/(coverage(idx)+1.0_wp)
               end do
            else
               ! Deterministic sign fix based on max element
               if (locev(local_max_ind) > 0.0_wp) then
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

      ! --- Gram-Schmidt Orthogonalization ---
      ! Project out previous columns from global ev
      do k = 1, n_curr
            ! d = dot(ev, displdir(:,k))
            d_dot = dot_product(ev, displdir_tmp(:, k))
            ! ev = ev - d * displdir(:,k)
            ev = ev - d_dot * displdir_tmp(:, k)
      end do

      ! --- Check Norm ---
      v_norm = sqrt(dot_product(ev, ev))
      
      if (v_norm < eps2) then
            early_break = .true.
            exit ! Break out of i loop
      else
            early_break = .false.
      end if

      ! Normalize and store
      ev = ev / v_norm
      displdir_tmp(:, n_curr + 1) = ev

      ndispl_final = n_curr + 1
   end do

   allocate(displdir(n, ndispl_final))
   displdir(:, :) = displdir_tmp(:, 1:ndispl_final)

end subroutine gen_displdir
end module xtb_o1numhess
