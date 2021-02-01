# This file is part of xtb.
#
# Copyright (C) 2019-2020 Sebastian Ehlert
#
# xtb is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with xtb.  If not, see <https://www.gnu.org/licenses/>.

fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')
dependencies = []

if fc.get_id() != cc.get_id()
  warning('FC and CC are not from the same vendor')
endif

if fc.get_id() == 'gcc'
  add_project_arguments('-fdefault-real-8', language: 'fortran')
  add_project_arguments('-fdefault-double-8', language: 'fortran')
  add_project_arguments('-ffree-line-length-none', language: 'fortran')
  add_project_arguments('-fbacktrace', language: 'fortran')
elif fc.get_id() == 'intel'
  add_project_arguments('-axAVX2',    language: 'fortran')
  add_project_arguments('-r8',        language: 'fortran')
  add_project_arguments('-traceback', language: 'fortran')
elif fc.get_id() == 'pgi' or fc.get_id() == 'nvidia_hpc'
  add_project_arguments(
    '-Mpreprocess', '-Mbackslash', '-Mallocatable=03', '-traceback', '-r8',
    language: 'fortran'
  )

  if get_option('gpu')
    add_project_arguments('-acc', '-Minfo=accel', '-DXTB_GPU', language: 'fortran')
    add_project_link_arguments('-acc', '-Minfo=accel', language: 'fortran')

    gpu_arch = get_option('gpu_arch') 
    add_project_arguments('-ta=tesla:cc@0@'.format(gpu_arch), language: 'fortran')
    add_project_link_arguments('-ta=tesla:cc@0@'.format(gpu_arch), language: 'fortran')

    if get_option('cusolver')
      add_project_arguments('-Mcudalib=cusolver,cublas', '-DUSE_CUSOLVER', '-DUSE_CUBLAS', language: 'fortran')
      add_project_link_arguments('-Mcudalib=cusolver,cublas', language: 'fortran')
    endif
  endif
endif

# fix compiliation problems with of symmetry/symmetry_i.c
add_project_arguments('-D_Float128=__float128', language: 'c')

## ========================================== ##
## LIBRARIES
## ========================================== ##

la_backend = get_option('la_backend')
if la_backend == 'mkl' or la_backend == 'mkl-static'
  add_project_arguments('-DWITH_MKL', language: 'fortran')
  if la_backend == 'mkl-static'
    add_project_link_arguments('-static', language: 'fortran')
    add_project_link_arguments('-static', language: 'c') # icc will do linking
  endif

  mkl_rt_dep = fc.find_library('mkl_rt', required: true)

  if fc.get_id() == 'gcc'
    libmkl_exe = [fc.find_library('mkl_gf_lp64')]
    if get_option('openmp')
      libmkl_exe += fc.find_library('mkl_gnu_thread')
    endif
  elif fc.get_id() == 'intel'
    libmkl_exe = [fc.find_library('mkl_intel_lp64')]
    if get_option('openmp')
      libmkl_exe += fc.find_library('mkl_intel_thread')
    endif
  elif fc.get_id() == 'pgi' or fc.get_id() == 'nvidia_hpc'
    libmkl_exe = [fc.find_library('mkl_intel_lp64')]
    if get_option('openmp')
      libmkl_exe += fc.find_library('mkl_pgi_thread')
    endif
  endif
  if not get_option('openmp')
    libmkl_exe += fc.find_library('mkl_sequential')
  endif
  libmkl_exe += fc.find_library('mkl_core')
  lib_deps += mkl_rt_dep
  exe_deps += libmkl_exe

elif la_backend == 'mkl-rt'
  add_project_arguments('-DWITH_MKL', language: 'fortran')

  mkl_rt_dep = fc.find_library('mkl_rt', required: true)
  if fc.get_id() == 'intel'
    exe_deps += fc.find_library('ifcore')
  endif

  dependencies += mkl_rt_dep

elif la_backend == 'openblas'
  # search for OpenBLAS
  blas_dep = dependency('openblas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('openblas', required: true)
  endif
  dependencies += blas_dep
  # some OpenBLAS versions can provide lapack, check if we can find dsygvd
  openblas_provides_lapack = fc.links(
    files('openblas_lapack_test.f90'),
    dependencies: blas_dep,
  )
  # otherwise we fall back to LAPACK
  if not openblas_provides_lapack
    lapack_dep = dependency('lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack', required: true)
    endif
    dependencies += lapack_dep
  endif

elif la_backend == 'custom'
  foreach lib: get_option('custom_libraries')
    dependencies += fc.find_library(lib)
  endforeach

else
  # Find BLAS (usually netlib, but in conda also OpenBLAS/MKL)
  blas_dep = dependency('blas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('blas', required: true)
  endif
  dependencies += blas_dep
  # Find LAPACK (usually netlib, but in conda also MKL)
  lapack_dep = dependency('lapack', required: false)
  if not lapack_dep.found()
    lapack_dep = fc.find_library('lapack', required: true)
  endif
  dependencies += lapack_dep

endif

if get_option('openmp')
  omp_dep = dependency('openmp', required: fc.get_id() != 'intel' and fc.get_id() != 'nvidia_hpc')
  if not omp_dep.found()
    if fc.get_id() == 'intel'
      message('Using -qopenmp to use OpenMP with Intel compilers')
      omp_dep = declare_dependency(
        compile_args: '-qopenmp',
        link_args: '-qopenmp',
      )
    else
      message('Using -mp to use OpenMP with NVHPC compilers')
      omp_dep = declare_dependency(
        compile_args: '-mp',
        link_args: '-mp',
      )
    endif
  endif
  dependencies += omp_dep
endif

dependencies += dependency('threads')

if get_option('nvtx')
  dependencies += fc.find_library('nvToolsExt', required: true)
endif

# distribute dependencies for shared object and static executable
lib_deps += dependencies
exe_deps += dependencies
