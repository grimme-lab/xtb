#!/usr/bin/env bash
export NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
export PATH="$NVHPC/compilers/bin:/usr/local/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="$NVHPC/compilers/lib:$NVHPC/math_libs/lib64:$LD_LIBRARY_PATH"
cd /tmp
cat > acctest.f90 <<'F90'
program acctest
  implicit none
  integer, parameter :: n = 1000000
  real(8), allocatable :: a(:)
  real(8) :: s
  integer :: i
  allocate(a(n))
  !$acc parallel loop
  do i = 1, n
     a(i) = sqrt(real(i,8)) * 2.0d0
  end do
  s = 0.0d0
  !$acc parallel loop reduction(+:s)
  do i = 1, n
     s = s + a(i)
  end do
  print '(a,f0.3)', "sum = ", s
  print *, "OpenACC GPU kernel ran."
end program
F90

echo "=== compile: -acc -gpu=cc86 (default LLVM) ==="
nvfortran -acc -gpu=cc86 -Minfo=accel -o acctest_llvm acctest.f90 2>&1 | grep -iE "acctest|Generating|error" | head
echo "--- run (LLVM) ---"; ./acctest_llvm 2>&1 | head

echo
echo "=== compile: -acc -gpu=cc86 -Mnollvm (classic) ==="
nvfortran -acc -gpu=cc86 -Mnollvm -Minfo=accel -o acctest_classic acctest.f90 2>&1 | grep -iE "acctest|Generating|error|nollvm|not supported" | head
echo "--- run (classic) ---"; ./acctest_classic 2>&1 | head
