! This file is part of xtb.
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

!> Runtime control for the CUDA acceleration paths.
!>
!> One switch (--gpu [dev]) turns on every GPU path in xtb:
!>   * generalised eigensolver in SCC (cuSOLVER Dsygvd/Ssygvd, precision ramp)
!>   * standard eigensolver for the ANC Hessian (cuSOLVER Dsyevd, FP64)
!>
!> Every path falls back to the CPU implementation on any failure (no device,
!> out of memory, solver error), so results stay correct without a GPU.
!>
!> The device state is process global, so the GPU is only used from serial
!> code: inside an OpenMP parallel region (the numerical Hessian runs one SCF
!> per thread) every query below answers "use the CPU".
!>
!> The C entry points below come from src/gpu/xtb_gpu_lin.c when the project
!> is configured with -Dgpu=true, and from src/gpu/xtb_gpu_stub.c otherwise.
!> The stub reports "not compiled" and every query returns false, so no source
!> outside this module needs preprocessor guards.
module xtb_gpu_ctrl
   use, intrinsic :: iso_c_binding
   use xtb_mctc_accuracy, only : wp
   !$ use omp_lib, only : omp_in_parallel
   implicit none
   private

   public :: TGPUConfig, gpuConfig
   public :: gpuSetup, gpuBanner, gpuReport
   public :: gpuCompiledIn, gpuEigActive
   public :: gpuSygvd, gpuSyevd, gpuRampReset, gpuAcceptDensity

   !> Everything the user can tune.  Filled from the command line, then
   !> overridden by XTB_GPU_* environment variables (debugging escape hatch).
   type :: TGPUConfig
      !> master switch, set by --gpu
      logical :: enabled = .false.
      !> CUDA device ordinal (--gpu <dev>).  Negative means auto: pick the
      !> largest visible device by SM count.  An explicit <dev> is interpreted
      !> in CUDA's own ordering, which is "fastest first" rather than PCI
      !> order; export CUDA_DEVICE_ORDER=PCI_BUS_ID if you mean the latter.
      integer :: device = -1
      !> use the GPU eigensolvers
      logical :: eig = .true.
      !> smallest matrix dimension worth sending to the GPU.  Below this the
      !> transfer plus launch latency costs more than the solve saves.
      integer :: nmin = 700
      !> precision ramp: run the SCC eigensolve in FP32 while the Hamiltonian
      !> is still moving, switch permanently to FP64 once max|dH_ii| drops
      !> below this (eV).  Set to zero to stay in FP64 throughout.
      real(wp) :: fp32_tol = 1.0e-2_wp
      !> print a timing breakdown at the end of the run
      logical :: stats = .false.
   end type TGPUConfig

   type(TGPUConfig), save :: gpuConfig

   interface
      function xtb_gpu_compiled() bind(C, name="xtb_gpu_compiled") result(stat)
         import :: c_int
         integer(c_int) :: stat
      end function xtb_gpu_compiled

      function xtb_gpu_open(dev, verbose) bind(C, name="xtb_gpu_open") result(stat)
         import :: c_int
         integer(c_int), value :: dev, verbose
         integer(c_int) :: stat
      end function xtb_gpu_open

      function xtb_gpu_sygvd_c(n, amat, bmat, eval, tol) &
            & bind(C, name="xtb_gpu_sygvd") result(stat)
         import :: c_int, c_double
         integer(c_int), value :: n
         real(c_double), intent(inout) :: amat(*)
         real(c_double), intent(in) :: bmat(*)
         real(c_double), intent(out) :: eval(*)
         real(c_double), value :: tol
         integer(c_int) :: stat
      end function xtb_gpu_sygvd_c

      function xtb_gpu_syevd_c(n, amat, eval) &
            & bind(C, name="xtb_gpu_syevd") result(stat)
         import :: c_int, c_double
         integer(c_int), value :: n
         real(c_double), intent(inout) :: amat(*)
         real(c_double), intent(out) :: eval(*)
         integer(c_int) :: stat
      end function xtb_gpu_syevd_c

      !> device actually in use, resolved from --gpu (or auto), -1 if closed
      function xtb_gpu_device() bind(C, name="xtb_gpu_device") result(dev)
         import :: c_int
         integer(c_int) :: dev
      end function xtb_gpu_device

      subroutine xtb_gpu_ramp_reset() bind(C, name="xtb_gpu_ramp_reset")
      end subroutine xtb_gpu_ramp_reset

      function xtb_gpu_accept_density() bind(C, name="xtb_gpu_accept_density") result(ok)
         import :: c_int
         integer(c_int) :: ok
      end function xtb_gpu_accept_density

      subroutine xtb_gpu_report() bind(C, name="xtb_gpu_report")
      end subroutine xtb_gpu_report

      subroutine xtb_gpu_name(buf, buflen) bind(C, name="xtb_gpu_name")
         import :: c_char, c_int
         character(kind=c_char), intent(inout) :: buf(*)
         integer(c_int), value :: buflen
      end subroutine xtb_gpu_name
   end interface

contains

!> Was the binary built with CUDA support?
logical function gpuCompiledIn()
   gpuCompiledIn = xtb_gpu_compiled() /= 0
end function gpuCompiledIn

!> Apply the command line request, then the environment overrides, and open
!> the device.  Called once from the driver after argument parsing.
!> stat is non-zero if --gpu was asked for but the device cannot be used:
!>   1 = binary built without -Dgpu=true, 2 = no usable CUDA device.
subroutine gpuSetup(stat)
   integer, intent(out) :: stat

   stat = 0
   call envLogical("XTB_GPU", gpuConfig%enabled)
   call envInteger("XTB_GPU_DEV", gpuConfig%device)
   call envLogical("XTB_GPU_EIG", gpuConfig%eig)
   call envInteger("XTB_GPU_NMIN", gpuConfig%nmin)
   call envReal("XTB_GPU_FP32_TOL", gpuConfig%fp32_tol)
   call envLogical("XTB_GPU_STATS", gpuConfig%stats)

   if (.not. gpuConfig%enabled) return

   if (xtb_gpu_compiled() == 0) then
      gpuConfig%enabled = .false.
      stat = 1
      return
   end if

   if (xtb_gpu_open(int(gpuConfig%device, c_int), &
         & merge(1_c_int, 0_c_int, gpuConfig%stats)) /= 0) then
      gpuConfig%enabled = .false.
      stat = 2
      return
   end if

end subroutine gpuSetup

!> One line for the program header, e.g.
!>    device 0 (NVIDIA RTX 5000 Ada Generation, sm89, 31.5 GB), eigensolver
subroutine gpuBanner(line)
   character(len=:), allocatable, intent(out) :: line
   character(kind=c_char) :: buf(128)
   character(len=128) :: name
   character(len=160) :: head
   character(len=64) :: paths
   integer :: i

   if (.not. gpuConfig%enabled) then
      line = ""
      return
   end if

   buf = c_null_char
   call xtb_gpu_name(buf, 128_c_int)
   name = ""
   do i = 1, 128
      if (buf(i) == c_null_char) exit
      name(i:i) = buf(i)
   end do

   paths = "none active"
   if (gpuConfig%eig) paths = "eigensolver"

   write(head, '(a,i0,3a)') "device ", int(xtb_gpu_device()), " (", trim(name), ")"
   line = trim(head)//", "//trim(paths)

end subroutine gpuBanner

!> Timing breakdown, only with --gpu and XTB_GPU_STATS=1.
subroutine gpuReport()
   if (gpuConfig%enabled .and. gpuConfig%stats) call xtb_gpu_report()
end subroutine gpuReport

!> Should a matrix of dimension n go to the GPU eigensolver?
logical function gpuEigActive(n)
   integer, intent(in) :: n
   !> cuSOLVER's [S|D]sygvd take the workspace size as an int; it is about 3 n^2
   !> elements and would overflow above n ~ 26,700.  Larger matrices stay on the CPU.
   integer, parameter :: nmax = 26000
   gpuEigActive = gpuConfig%enabled .and. gpuConfig%eig .and. n >= gpuConfig%nmin &
      & .and. n <= nmax .and. .not.inParallel()
end function gpuEigActive

!> Generalised symmetric eigenproblem A x = lambda B x, itype=1, jobz='v',
!> uplo='u'.  On success A holds the eigenvectors, as with LAPACK dsygvd, but
!> B is left untouched, so one copy of the overlap serves a whole SCF.
!> Returns 0 on success.  A negative value means "nothing was touched, use the
!> CPU path"; a positive one means A was partly overwritten and is lost.
function gpuSygvd(n, amat, bmat, eval) result(stat)
   integer, intent(in) :: n
   real(wp), intent(inout), contiguous :: amat(:, :)
   real(wp), intent(in), contiguous :: bmat(:, :)
   real(wp), intent(out), contiguous :: eval(:)
   integer :: stat
   stat = int(xtb_gpu_sygvd_c(int(n, c_int), amat, bmat, eval, &
      & real(gpuConfig%fp32_tol, c_double)))
end function gpuSygvd

!> Standard symmetric eigenproblem, jobz='v', uplo='u', always FP64.
!> Same return convention as gpuSygvd.
function gpuSyevd(n, amat, eval) result(stat)
   integer, intent(in) :: n
   real(wp), intent(inout), contiguous :: amat(:, :)
   real(wp), intent(out), contiguous :: eval(:)
   integer :: stat
   stat = int(xtb_gpu_syevd_c(int(n, c_int), amat, eval))
end function gpuSyevd

!> Tell the precision ramp that a new SCF starts here.  Without this a
!> geometry optimisation would stay in FP64 from the second cycle onwards.
subroutine gpuRampReset()
   if (gpuConfig%enabled .and. gpuConfig%eig .and. .not.inParallel()) &
      & call xtb_gpu_ramp_reset()
end subroutine gpuRampReset

!> Asked by the SCC once it has converged: may it stop on this density?  False
!> if the latest solve was FP32; the ramp then switches to FP64 and the SCC has
!> to iterate once more, so a converged SCF never ends on an FP32 solve.
logical function gpuAcceptDensity()
   gpuAcceptDensity = .true.
   if (gpuConfig%enabled .and. gpuConfig%eig .and. .not.inParallel()) &
      & gpuAcceptDensity = xtb_gpu_accept_density() /= 0
end function gpuAcceptDensity

!> The C backend keeps one cuSOLVER handle, one set of device buffers and one
!> precision ramp for the whole process, so it must never be entered from
!> several threads at once.
logical function inParallel()
   inParallel = .false.
   !$ inParallel = omp_in_parallel()
end function inParallel

! --------------------------------------------------------------- env helpers

subroutine envLogical(key, val)
   character(len=*), intent(in) :: key
   logical, intent(inout) :: val
   character(len=8) :: buf
   integer :: stat
   buf = ""
   call get_environment_variable(key, buf, status=stat)
   if (stat /= 0) return
   select case(trim(buf))
   case("1", "true", "TRUE", "yes", "on"); val = .true.
   case("0", "false", "FALSE", "no", "off"); val = .false.
   end select
end subroutine envLogical

subroutine envInteger(key, val)
   character(len=*), intent(in) :: key
   integer, intent(inout) :: val
   character(len=32) :: buf
   integer :: stat, tmp
   buf = ""
   call get_environment_variable(key, buf, status=stat)
   if (stat /= 0 .or. len_trim(buf) == 0) return
   read(buf, *, iostat=stat) tmp
   if (stat == 0) val = tmp
end subroutine envInteger

subroutine envReal(key, val)
   character(len=*), intent(in) :: key
   real(wp), intent(inout) :: val
   character(len=32) :: buf
   integer :: stat
   real(wp) :: tmp
   buf = ""
   call get_environment_variable(key, buf, status=stat)
   if (stat /= 0 .or. len_trim(buf) == 0) return
   read(buf, *, iostat=stat) tmp
   if (stat == 0) val = tmp
end subroutine envReal

end module xtb_gpu_ctrl
