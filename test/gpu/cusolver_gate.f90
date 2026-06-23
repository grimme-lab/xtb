! Standalone GPU numerical gate: solve real GFN0 generalized eigenproblems
!   H C = S C eps   (itype=1, vectors, upper)
! via cuSolver (cusolverDnDsygvd, GPU) and LAPACK (dsygvd, CPU), and compare the
! eigenvalues. Mirrors the cuSolver path in xtb/src/gpu/batched_eig.F90. Also
! repeats the comparison with the bucket-padding xtb's --gpu-batch uses (large H
! diagonal + identity S on the padded block) to prove padding decouples on GPU.
!
! Input file (XTB_DUMP_HS from the CPU build): count, then per system nao, H, S.
program cusolver_gate
   use cusolverDn
   use cublas
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   type tsys
      integer :: n = 0
      real(dp), allocatable :: H(:,:), S(:,:)
   end type tsys

   type(tsys), allocatable :: sys(:)
   type(cusolverDnHandle) :: h
   character(len=512) :: path
   integer :: nsys, k, i, j, n, m, maxn, istat, length, stat, u
   real(dp), allocatable :: eref(:), egpu(:), A(:,:), B(:,:)
   real(dp), allocatable :: Hp(:,:), Sp(:,:)
   real(dp) :: dev, worst_direct, worst_padded
   real(dp), parameter :: padDiag = 1.0d6

   ! --- read the dumped matrices ---
   call get_environment_variable("XTB_DUMP_HS", path, length, stat)
   if (stat /= 0 .or. length == 0) path = "/tmp/gfn0_hs.txt"
   open(newunit=u, file=trim(path), status='old', action='read', iostat=stat)
   if (stat /= 0) then
      print *, "ERROR: cannot open ", trim(path); stop 1
   end if
   read(u,*) nsys
   allocate(sys(nsys))
   maxn = 0
   do k = 1, nsys
      read(u,*) n
      sys(k)%n = n
      allocate(sys(k)%H(n,n), sys(k)%S(n,n))
      read(u,*) sys(k)%H
      read(u,*) sys(k)%S
      maxn = max(maxn, n)
   end do
   close(u)
   print '(a,i0,a,i0)', "loaded ", nsys, " real GFN0 eigenproblems; max nao = ", maxn

   istat = cusolverDnCreate(h)
   if (istat /= 0) then; print *, "cusolverDnCreate failed", istat; stop 1; end if

   ! ================================================================
   ! Test 1: per-system, bare m x m   (cuSolver vs LAPACK)
   ! ================================================================
   print *, ""
   print *, "Test 1: cuSolver vs LAPACK, per system (bare)"
   print *, "-------------------------------------------------------"
   worst_direct = 0.0d0
   do k = 1, nsys
      n = sys(k)%n
      allocate(eref(n), egpu(n))
      ! LAPACK reference
      A = sys(k)%H; B = sys(k)%S
      call cpu_sygvd(n, A, B, eref, stat)
      ! cuSolver
      A = sys(k)%H; B = sys(k)%S
      call gpu_sygvd(h, n, A, B, egpu, istat)
      dev = maxval(abs(egpu - eref))
      worst_direct = max(worst_direct, dev)
      print '(a,i3,a,i6,a,es10.3)', "  sys ", k, "  n=", n, "  max|dEPS|=", dev
      deallocate(eref, egpu)
   end do
   print '(a,es10.3)', "  worst |dEPS| (direct) = ", worst_direct

   ! ================================================================
   ! Test 2: padded to bucket size maxn (cuSolver) vs bare LAPACK
   !         (exactly the xtb --gpu-batch padding scheme)
   ! ================================================================
   print *, ""
   print *, "Test 2: cuSolver on padded bucket vs LAPACK (bare)"
   print *, "-------------------------------------------------------"
   worst_padded = 0.0d0
   do k = 1, nsys
      n = sys(k)%n
      allocate(eref(n))
      A = sys(k)%H; B = sys(k)%S
      call cpu_sygvd(n, A, B, eref, stat)

      allocate(Hp(maxn,maxn), Sp(maxn,maxn), egpu(maxn))
      Hp = 0.0d0; Sp = 0.0d0
      Hp(1:n,1:n) = sys(k)%H
      Sp(1:n,1:n) = sys(k)%S
      do i = n+1, maxn
         Hp(i,i) = padDiag
         Sp(i,i) = 1.0d0
      end do
      call gpu_sygvd(h, maxn, Hp, Sp, egpu, istat)
      dev = maxval(abs(egpu(1:n) - eref(1:n)))   ! lowest n must match bare spectrum
      worst_padded = max(worst_padded, dev)
      print '(a,i3,a,i6,a,i6,a,es10.3)', "  sys ", k, "  n=", n, " -> padded ", maxn, &
         &  "  max|dEPS|=", dev
      deallocate(eref, egpu, Hp, Sp)
   end do
   print '(a,es10.3)', "  worst |dEPS| (padded) = ", worst_padded

   istat = cusolverDnDestroy(h)

   print *, ""
   print *, "======================================================="
   if (worst_direct <= 1.0d-6 .and. worst_padded <= 1.0d-6) then
      print '(a)', " GPU GATE: PASS  (cuSolver == LAPACK within 1e-6 eV,"
      print '(a)', "                  bare and padded, on real GFN0 matrices)"
   else
      print '(a)', " GPU GATE: FAIL  (deviation exceeds 1e-6 eV)"
   end if
   print *, "======================================================="

contains

   !> CPU reference: LAPACK dsygvd (itype=1, vectors, upper).
   subroutine cpu_sygvd(n, A, B, w, info)
      integer, intent(in) :: n
      real(dp), intent(inout) :: A(n,n), B(n,n)
      real(dp), intent(out) :: w(n)
      integer, intent(out) :: info
      real(dp), allocatable :: work(:)
      integer, allocatable :: iwork(:)
      real(dp) :: qw(1)
      integer :: qiw(1), lwork, liwork
      call dsygvd(1,'V','U',n,A,n,B,n,w,qw,-1,qiw,-1,info)
      lwork = max(1,int(qw(1))); liwork = max(1,qiw(1))
      allocate(work(lwork), iwork(liwork))
      call dsygvd(1,'V','U',n,A,n,B,n,w,work,lwork,iwork,liwork,info)
      deallocate(work, iwork)
   end subroutine cpu_sygvd

   !> GPU: cusolverDnDsygvd (itype=1, vectors, upper). A,B staged on device via
   !> OpenACC; device pointers handed to cuSolver with host_data use_device.
   subroutine gpu_sygvd(handle, n, A, B, w, istat)
      type(cusolverDnHandle), intent(in) :: handle
      integer, intent(in) :: n
      real(dp), intent(inout) :: A(n,n), B(n,n)
      real(dp), intent(out) :: w(n)
      integer, intent(out) :: istat
      real(dp), allocatable :: work(:)
      real(dp) :: dummy(1)
      integer :: lwork, devInfo
      istat = cusolverDnDsygvd_bufferSize(handle, CUSOLVER_EIG_TYPE_1, &
         & CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, dummy, n, &
         & dummy, n, dummy, lwork)
      allocate(work(lwork))
      devInfo = 0
      !$acc data copy(A, B, w) create(work, devInfo)
      !$acc host_data use_device(A, B, w, work, devInfo)
      istat = cusolverDnDsygvd(handle, CUSOLVER_EIG_TYPE_1, &
         & CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, n, A, n, B, n, &
         & w, work, lwork, devInfo)
      !$acc end host_data
      !$acc end data
      deallocate(work)
   end subroutine gpu_sygvd

end program cusolver_gate
