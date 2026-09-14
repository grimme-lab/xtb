! Validate + benchmark the CUDA-C batched eigensolver (called from gfortran via
! iso_c_binding) against CPU LAPACK, on many small generalized eigenproblems --
! the "screen lots of small molecules" regime. Answers: does GPU beat CPU here?
program gpu_bench
   use iso_c_binding
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   interface
      integer(c_int) function gpu_sygvd_batch(n, nbatch, H, S, W) &
            bind(C, name="gpu_sygvd_batch")
         import :: c_int, c_double
         integer(c_int), value :: n, nbatch
         real(c_double), intent(in) :: H(*), S(*)
         real(c_double), intent(out) :: W(*)
      end function
   end interface

   integer :: n, nbatch, k, i, j, info, rc
   real(dp), allocatable :: H(:,:,:), S(:,:,:), Wgpu(:,:), Wcpu(:)
   real(dp), allocatable :: A(:,:), B(:,:), work(:)
   integer,  allocatable :: iwork(:)
   real(dp) :: qw(1), r, worst, tgpu, tcpu
   integer  :: qiw(1), lwork, liwork, seed
   integer(8) :: c0, c1, cr

   n = 48; nbatch = 200          ! ~typical small molecule, screening batch size
   call get_command_argument_count_compat(n, nbatch)

   allocate(H(n,n,nbatch), S(n,n,nbatch), Wgpu(n,nbatch), Wcpu(n))
   ! Build nbatch random symmetric H and SPD S (S = M^T M + n I, well conditioned)
   seed = 12345
   allocate(A(n,n), B(n,n))
   do k = 1, nbatch
      do j = 1, n
         do i = 1, n
            A(i,j) = ran(seed); B(i,j) = ran(seed)
         end do
      end do
      H(:,:,k) = 0.5_dp*(A + transpose(A))
      S(:,:,k) = matmul(transpose(B), B)
      do i = 1, n
         S(i,i,k) = S(i,i,k) + real(n, dp)
      end do
   end do

   ! ---- GPU (batched shim) ---- wall-clock
   call system_clock(c0, cr)
   rc = gpu_sygvd_batch(n, nbatch, H, S, Wgpu)
   call system_clock(c1); tgpu = real(c1 - c0, dp) / real(cr, dp)
   if (rc /= 0) print '(a,i0)', "WARNING: gpu_sygvd_batch returned ", rc

   ! ---- CPU (LAPACK loop) + parity check ---- wall-clock
   worst = 0.0_dp
   call dsygvd(1,'V','U',n,H(:,:,1),n,S(:,:,1),n,Wcpu,qw,-1,qiw,-1,info)
   lwork = max(1,int(qw(1))); liwork = max(1,qiw(1))
   allocate(work(lwork), iwork(liwork))
   call system_clock(c0, cr)
   do k = 1, nbatch
      A = H(:,:,k); B = S(:,:,k)            ! work on copies (GPU got its own copies)
      call dsygvd(1,'V','U',n,A,n,B,n,Wcpu,work,lwork,iwork,liwork,info)
      r = maxval(abs(Wcpu - Wgpu(:,k)))
      if (r > worst) worst = r
   end do
   call system_clock(c1); tcpu = real(c1 - c0, dp) / real(cr, dp)

   print '(a,i0,a,i0,a)', "batch: ", nbatch, " systems of size ", n, " (generalized sym-def eig)"
   print '(a,es10.3,a)',  "max |W_gpu - W_cpu| = ", worst, " (parity vs LAPACK)"
   print '(a,f8.3,a)',    "GPU (cuSolver loop) : ", tgpu, " s"
   print '(a,f8.3,a)',    "CPU (LAPACK loop)   : ", tcpu, " s"
   if (tgpu > 0.0_dp) print '(a,f6.2,a)', "speedup CPU/GPU     : ", tcpu/tgpu, "x"

contains
   ! tiny LCG so the test is deterministic and self-contained
   real(dp) function ran(s)
      integer, intent(inout) :: s
      s = mod(1103515245*s + 12345, 2147483647)
      ran = real(s, dp) / 2147483647.0_dp - 0.5_dp
   end function ran

   subroutine get_command_argument_count_compat(n, nbatch)
      integer, intent(inout) :: n, nbatch
      character(len=32) :: a
      if (command_argument_count() >= 1) then
         call get_command_argument(1, a); read(a,*) n
      end if
      if (command_argument_count() >= 2) then
         call get_command_argument(2, a); read(a,*) nbatch
      end if
   end subroutine
end program gpu_bench
