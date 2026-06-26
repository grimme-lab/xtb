! This file is part of xtb.
!
! Copyright (C) 2026 xtb GPU contributors
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

!> Low-level GPU runtime flag + the CUDA-C eigensolver interface, in a *leaf*
!> module so even core routines (e.g. xtb_scc_core) can route their
!> diagonalization to the GPU without creating a dependency cycle.
!>
!> `gpu_use` is set by the `--gpu` CLI flag. `gpu_sygvd_batch` (src/gpu/gpu_eig.cu,
!> compiled by nvcc and linked in only when WITH_GPU_SHIM is defined) solves one
!> or a batch of generalized symmetric-definite eigenproblems on the GPU.
module xtb_gpu_runtime
   use xtb_mctc_accuracy, only : wp
   use iso_c_binding, only : c_long_long
#ifdef WITH_GPU_SHIM
   use iso_c_binding, only : c_int, c_double
#endif
   implicit none
   private

   public :: gpu_use, gpu_solve, gpu_density, gpu_build_h1, gpu_mpopsh
   public :: gpu_shim_available, gpu_solve_count
   public :: gpu_scf_open, gpu_scf_solve, gpu_scf_finish, gpu_scf_close
   public :: gpu_scf_get_vectors
   public :: gpu_grad_dsdqh0
#ifdef WITH_GPU_SHIM
   public :: gpu_sygvd_batch
#endif

   !> Route diagonalizations to the GPU (cuSolver shim). Set by `--gpu`.
   logical :: gpu_use = .false.

#ifdef WITH_GPU_SHIM
   !> H_k C_k = S_k C_k diag(W_k), k=1..nbatch. H in: Hamiltonian blocks (n*n
   !> column-major); out: eigenvectors. S: overlap blocks. W out: eigenvalues.
   interface
      function gpu_sygvd_batch(n, nbatch, H, S, W) result(rc) &
            & bind(C, name="gpu_sygvd_batch")
         import :: c_int, c_double
         integer(c_int), value :: n, nbatch
         real(c_double), intent(inout) :: H(*)
         real(c_double), intent(in)    :: S(*)
         real(c_double), intent(out)   :: W(*)
         integer(c_int) :: rc
      end function gpu_sygvd_batch

      function gpu_sygvd_solve_count() result(count) &
            & bind(C, name="gpu_sygvd_solve_count")
         import :: c_long_long
         integer(c_long_long) :: count
      end function gpu_sygvd_solve_count

      function gpu_density_matrix(n, C, f, P) result(rc) &
            & bind(C, name="gpu_density_matrix")
         import :: c_int, c_double
         integer(c_int), value :: n
         real(c_double), intent(in) :: C(*)
         real(c_double), intent(in) :: f(*)
         real(c_double), intent(out) :: P(*)
         integer(c_int) :: rc
      end function gpu_density_matrix

      function gpu_build_isotropic_h1(n, nmat, nshell, matlist, H0, S, &
            & shift, ao2sh, autoev, H) result(rc) &
            & bind(C, name="gpu_build_isotropic_h1")
         import :: c_int, c_double
         integer(c_int), value :: n, nmat, nshell
         integer(c_int), intent(in) :: matlist(*)
         real(c_double), intent(in) :: H0(*), S(*), shift(*)
         integer(c_int), intent(in) :: ao2sh(*)
         real(c_double), value :: autoev
         real(c_double), intent(out) :: H(*)
         integer(c_int) :: rc
      end function gpu_build_isotropic_h1

      function gpu_mulliken_shell(n, nshell, ao2sh, S, P, qsh) result(rc) &
            & bind(C, name="gpu_mulliken_shell")
         import :: c_int, c_double
         integer(c_int), value :: n, nshell
         integer(c_int), intent(in) :: ao2sh(*)
         real(c_double), intent(in) :: S(*), P(*)
         real(c_double), intent(out) :: qsh(*)
         integer(c_int) :: rc
      end function gpu_mulliken_shell

      function gpu_scf_open_i(n, nshell, nmat, gfn1, H0, S, matlist, ao2sh) &
            & result(rc) bind(C, name="gpu_scf_open")
         import :: c_int, c_double
         integer(c_int), value :: n, nshell, nmat, gfn1
         real(c_double), intent(in) :: H0(*), S(*)
         integer(c_int), intent(in) :: matlist(*), ao2sh(*)
         integer(c_int) :: rc
      end function gpu_scf_open_i

      function gpu_scf_solve_i(n, nshell, H, shift, autoev, emo) result(rc) &
            & bind(C, name="gpu_scf_solve")
         import :: c_int, c_double
         integer(c_int), value :: n, nshell
         real(c_double), intent(in) :: H(*), shift(*)
         real(c_double), value :: autoev
         real(c_double), intent(out) :: emo(*)
         integer(c_int) :: rc
      end function gpu_scf_solve_i

      function gpu_scf_finish_i(n, nshell, focc, P, qsh) result(rc) &
            & bind(C, name="gpu_scf_finish")
         import :: c_int, c_double
         integer(c_int), value :: n, nshell
         real(c_double), intent(in) :: focc(*)
         real(c_double), intent(out) :: P(*), qsh(*)
         integer(c_int) :: rc
      end function gpu_scf_finish_i

      function gpu_scf_get_vectors_i(n, C) result(rc) &
            & bind(C, name="gpu_scf_get_vectors")
         import :: c_int, c_double
         integer(c_int), value :: n
         real(c_double), intent(out) :: C(*)
         integer(c_int) :: rc
      end function gpu_scf_get_vectors_i

      function gpu_scf_close_i() result(rc) bind(C, name="gpu_scf_close")
         import :: c_int
         integer(c_int) :: rc
      end function gpu_scf_close_i

      function gpu_build_dsdqh0_i(nat, nao) result(rc) &
            & bind(C, name="gpu_build_dsdqh0")
         import :: c_int
         integer(c_int), value :: nat, nao
         integer(c_int) :: rc
      end function gpu_build_dsdqh0_i
   end interface
#endif

contains

   !> Whether this build was compiled with the GPU shim linked in.
   pure logical function gpu_shim_available()
#ifdef WITH_GPU_SHIM
      gpu_shim_available = .true.
#else
      gpu_shim_available = .false.
#endif
   end function gpu_shim_available

   !> Number of generalized eigenproblems completed by cuSolver in this process.
   integer(c_long_long) function gpu_solve_count()
#ifdef WITH_GPU_SHIM
      gpu_solve_count = gpu_sygvd_solve_count()
#else
      gpu_solve_count = 0
#endif
   end function gpu_solve_count

   !> Solve one generalized symmetric-definite eigenproblem H C = S C diag(W) on
   !> the GPU.  H in: Hamiltonian (n x n); out: eigenvectors.  S: overlap.
   !> W out: ascending eigenvalues.  ok = .false. if no shim / GPU failure.
   !> Plain Fortran interface so non-preprocessed (.f90) core code can call it.
   subroutine gpu_solve(n, H, S, W, ok)
      integer, intent(in) :: n
      real(wp), intent(inout) :: H(n, n)
      real(wp), intent(in)    :: S(n, n)
      real(wp), intent(out)   :: W(n)
      logical, intent(out)    :: ok
#ifdef WITH_GPU_SHIM
      integer(c_int) :: rc
      logical, save :: announced = .false.
      rc = gpu_sygvd_batch(int(n, c_int), 1_c_int, H, S, W)
      ok = (rc == 0)
      if (ok .and. .not.announced) then
         write(*, '(a)') " GPU cuSolver SCC enabled"
         announced = .true.
      end if
#else
      ok = .false.
#endif
   end subroutine gpu_solve

   !> Form P = C diag(f) C^T on the GPU. Returns ok=.false. when the CUDA shim
   !> is unavailable or reports an error, allowing a transparent CPU fallback.
   subroutine gpu_density(n, C, f, P, ok)
      integer, intent(in) :: n
      real(wp), intent(in) :: C(n, n)
      real(wp), intent(in) :: f(n)
      real(wp), intent(out) :: P(n, n)
      logical, intent(out) :: ok
#ifdef WITH_GPU_SHIM
      integer(c_int) :: rc
      rc = gpu_density_matrix(int(n, c_int), C, f, P)
      ok = (rc == 0)
#else
      ok = .false.
#endif
   end subroutine gpu_density

   !> Build the isotropic GFN1 charge-dependent Hamiltonian on the GPU.
   subroutine gpu_build_h1(n, nmat, nshell, matlist, H0, S, shift, ao2sh, &
         & autoev, H, ok)
      integer, intent(in) :: n, nmat, nshell
      integer, intent(in) :: matlist(2, nmat), ao2sh(n)
      real(wp), intent(in) :: H0(nmat), S(n, n), shift(nshell), autoev
      real(wp), intent(out) :: H(n, n)
      logical, intent(out) :: ok
#ifdef WITH_GPU_SHIM
      integer(c_int) :: rc
      rc = gpu_build_isotropic_h1(int(n, c_int), int(nmat, c_int), &
         & int(nshell, c_int), matlist, H0, S, shift, ao2sh, autoev, H)
      ok = (rc == 0)
#else
      ok = .false.
#endif
   end subroutine gpu_build_h1

   !> Compute shell-resolved Mulliken populations on the GPU.
   subroutine gpu_mpopsh(n, nshell, ao2sh, S, P, qsh, ok)
      integer, intent(in) :: n, nshell
      integer, intent(in) :: ao2sh(n)
      real(wp), intent(in) :: S(n, n), P(n, n)
      real(wp), intent(out) :: qsh(nshell)
      logical, intent(out) :: ok
#ifdef WITH_GPU_SHIM
      integer(c_int) :: rc
      rc = gpu_mulliken_shell(int(n, c_int), int(nshell, c_int), &
         & ao2sh, S, P, qsh)
      ok = (rc == 0)
#else
      ok = .false.
#endif
   end subroutine gpu_mpopsh

   !> Open a resident SCF session: keep S (and, for GFN1, H0/matlist/ao2sh)
   !> resident on the device for the whole SCF.  For GFN1 set gfn1=.true. so the
   !> charge-dependent Hamiltonian is rebuilt on-device each cycle; for GFN2 set
   !> gfn1=.false. and pass the host-built H into gpu_scf_solve.
   subroutine gpu_scf_open(n, nshell, nmat, gfn1, H0, S, matlist, ao2sh, ok)
      integer, intent(in) :: n, nshell, nmat
      logical, intent(in) :: gfn1
      real(wp), intent(in) :: H0(*), S(n, n)
      integer, intent(in) :: matlist(*), ao2sh(n)
      logical, intent(out) :: ok
#ifdef WITH_GPU_SHIM
      integer(c_int) :: rc, gflag
      gflag = 0_c_int
      if (gfn1) gflag = 1_c_int
      rc = gpu_scf_open_i(int(n, c_int), int(nshell, c_int), int(nmat, c_int), &
         & gflag, H0, S, matlist, ao2sh)
      ok = (rc == 0)
#else
      ok = .false.
#endif
   end subroutine gpu_scf_open

   !> One SCF diagonalization with resident data.  GFN1: H is unused (built on
   !> the device from the resident H0/S and `shift`).  GFN2: H is the host-built
   !> Hamiltonian; `shift` is unused.  Eigenvectors stay resident; emo returns
   !> the eigenvalues.  ok=.false. on any failure (caller falls back to CPU).
   subroutine gpu_scf_solve(n, nshell, H, shift, autoev, emo, ok)
      integer, intent(in) :: n, nshell
      real(wp), intent(in) :: H(n, n), shift(nshell), autoev
      real(wp), intent(out) :: emo(n)
      logical, intent(out) :: ok
#ifdef WITH_GPU_SHIM
      integer(c_int) :: rc
      logical, save :: announced = .false.
      rc = gpu_scf_solve_i(int(n, c_int), int(nshell, c_int), H, shift, &
         & real(autoev, c_double), emo)
      ok = (rc == 0)
      if (ok .and. .not.announced) then
         write(*, '(a)') " GPU cuSolver SCC enabled"
         announced = .true.
      end if
#else
      ok = .false.
#endif
   end subroutine gpu_scf_solve

   !> Density (from the resident eigenvectors + occupations) and shell-Mulliken
   !> charges.  P is copied back (the host electro() energy needs it); qsh out.
   subroutine gpu_scf_finish(n, nshell, focc, P, qsh, ok)
      integer, intent(in) :: n, nshell
      real(wp), intent(in) :: focc(n)
      real(wp), intent(out) :: P(n, n), qsh(nshell)
      logical, intent(out) :: ok
#ifdef WITH_GPU_SHIM
      integer(c_int) :: rc
      rc = gpu_scf_finish_i(int(n, c_int), int(nshell, c_int), focc, P, qsh)
      ok = (rc == 0)
#else
      ok = .false.
#endif
   end subroutine gpu_scf_finish

   !> Fetch the resident eigenvectors (last solve) back to the host C, once,
   !> after convergence -- needed for the gradient's energy-weighted density.
   subroutine gpu_scf_get_vectors(n, C, ok)
      integer, intent(in) :: n
      real(wp), intent(out) :: C(n, n)
      logical, intent(out) :: ok
#ifdef WITH_GPU_SHIM
      integer(c_int) :: rc
      rc = gpu_scf_get_vectors_i(int(n, c_int), C)
      ok = (rc == 0)
#else
      ok = .false.
#endif
   end subroutine gpu_scf_get_vectors

   !> GPU analytical gradient (build_dSDQH0).  ok=.true. means the GPU filled
   !> g/sigma/dhdcn and the caller must skip the CPU build_dSDQH0; ok=.false.
   !> (current scaffold / no shim) means fall back to the CPU routine.
   subroutine gpu_grad_dsdqh0(nat, nao, ok)
      integer, intent(in) :: nat, nao
      logical, intent(out) :: ok
#ifdef WITH_GPU_SHIM
      ok = (gpu_build_dsdqh0_i(int(nat, c_int), int(nao, c_int)) == 0)
#else
      ok = .false.
#endif
   end subroutine gpu_grad_dsdqh0

   !> Release the resident SCF session.
   subroutine gpu_scf_close()
#ifdef WITH_GPU_SHIM
      integer(c_int) :: rc
      rc = gpu_scf_close_i()
#endif
   end subroutine gpu_scf_close

end module xtb_gpu_runtime
