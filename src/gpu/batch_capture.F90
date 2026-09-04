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

!> Capture of per-molecule generalized eigenproblems (H C = S C eps) at the
!> point where xtb diagonalizes them, so the batched eigensolver can be exercised
!> and numerically validated against the proven per-system path on REAL matrices.
!>
!> This is the validatable first half of the `--gpu-batch` INTEGRATION SEAM (see
!> xtb_gpu_batch): routing the *production* energy through a cross-molecule batched
!> solve requires splitting the monolithic single-point into build -> solve ->
!> finish phases. Before doing that invasive refactor, we prove the batched kernel
!> + the padding scheme reproduce the per-system spectrum bit-for-bit.
!>
!> The capture hook in the core single-point path (e.g. xtb_peeq) is INERT by
!> default: `gpu_capture_store` returns immediately unless `gpu_capture_enable`
!> has been called (only the `--gpu-batch` driver does so). The plain CPU/GPU
!> single-point thus pays at most one branch and behaves exactly as before.
!>
!> Memory is bounded: at most `cap` systems are retained (others are counted as
!> dropped and reported by the driver -- no silent truncation), so enabling
!> capture over a large screening set cannot blow up host memory.
module xtb_gpu_batch_capture
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: gpu_capture_enable, gpu_capture_disable, gpu_capture_active
   public :: gpu_capture_store, gpu_capture_clear
   public :: gpu_capture_count, gpu_capture_seen, gpu_capture_get
   public :: TCapturedSystem

   !> One captured generalized symmetric-definite eigenproblem.
   type :: TCapturedSystem
      !> Order of this system (number of AOs).
      integer :: nao = 0
      !> Symmetric Hamiltonian (full nao x nao, as handed to the solver).
      real(wp), allocatable :: H(:,:)
      !> Symmetric-positive-definite overlap (full nao x nao).
      real(wp), allocatable :: S(:,:)
   end type TCapturedSystem

   !> Capture is off until the batch driver enables it.
   logical :: enabled = .false.
   !> Maximum number of systems retained (host-memory bound).
   integer :: cap = 64
   !> Number of systems actually stored (<= cap).
   integer :: stored = 0
   !> Number of store calls observed (for dropped-count reporting).
   integer :: seen = 0
   !> Retained systems.
   type(TCapturedSystem), allocatable :: systems(:)

contains

!> Turn capture on and reset the buffer. `maxsys` caps how many systems are kept.
subroutine gpu_capture_enable(maxsys)
   integer, intent(in), optional :: maxsys
   if (present(maxsys)) cap = max(1, maxsys)
   if (allocated(systems)) deallocate(systems)
   allocate(systems(cap))
   stored = 0
   seen = 0
   enabled = .true.
end subroutine gpu_capture_enable

!> Turn capture off (retained systems remain available until cleared).
subroutine gpu_capture_disable()
   enabled = .false.
end subroutine gpu_capture_disable

!> Whether capture is currently armed.
pure logical function gpu_capture_active()
   gpu_capture_active = enabled
end function gpu_capture_active

!> Store one (H, S) pair. INERT unless capture is enabled; drops silently past
!> the cap (the drop is counted; `gpu_capture_seen` vs `gpu_capture_count`
!> exposes it to the caller).
subroutine gpu_capture_store(nao, H, S)
   integer, intent(in) :: nao
   real(wp), intent(in) :: H(nao, nao)
   real(wp), intent(in) :: S(nao, nao)
   if (.not. enabled) return
   if (nao <= 0) return
   seen = seen + 1
   if (stored >= cap) return
   stored = stored + 1
   systems(stored)%nao = nao
   systems(stored)%H = H(1:nao, 1:nao)
   systems(stored)%S = S(1:nao, 1:nao)
end subroutine gpu_capture_store

!> Number of systems retained.
pure integer function gpu_capture_count()
   gpu_capture_count = stored
end function gpu_capture_count

!> Number of systems seen at the hook (>= count when some were dropped).
pure integer function gpu_capture_seen()
   gpu_capture_seen = seen
end function gpu_capture_seen

!> Read back a retained system (1 <= i <= gpu_capture_count()).
subroutine gpu_capture_get(i, sys)
   integer, intent(in) :: i
   type(TCapturedSystem), intent(out) :: sys
   if (i < 1 .or. i > stored) return
   sys%nao = systems(i)%nao
   sys%H = systems(i)%H
   sys%S = systems(i)%S
end subroutine gpu_capture_get

!> Free all retained systems.
subroutine gpu_capture_clear()
   if (allocated(systems)) deallocate(systems)
   stored = 0
   seen = 0
end subroutine gpu_capture_clear

end module xtb_gpu_batch_capture
