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
fc_id = fc.get_id()
cc_id = cc.get_id()

if fc_id != cc_id
  warning('FC and CC are not from the same vendor')
endif

fopts = []
if fc_id == 'gcc'
  fopts = [
    '-fdefault-real-8',
    '-fdefault-double-8',
    '-ffree-line-length-none',
    '-fbacktrace',
  ]
elif fc_id == 'intel' or fc_id == 'intel-llvm'
  # ifort and ifx largely share the same options
  if get_option('buildtype') == 'release'
    opt_level = [
      '-Ofast',
      '-axAVX2',
      '-mtune=core-avx2',
      '-fma',
    ]
    # Not available with ifx
    if fc_id == 'intel'
      opt_level += ['-ip']
    endif  
  else
    opt_level = [
      '-axAVX2',
    ]
  endif
  fopts = opt_level + [
    '-r8',
    '-traceback',
  ]
elif fc_id == 'intel-cl'
  if get_option('buildtype') == 'release'
    opt_level = [
      '/O3',
      '/Qip',
      '/QaxAVX2',
      '/tune:core-avx2',
      '/Qfma',
    ]
  else
    opt_level = [
      '/QaxAVX2',
    ]
  endif
  fopts = opt_level + [
    '/fpp',
    '/4R8',
    '/traceback',
  ]
elif fc_id == 'pgi' or fc_id == 'nvidia_hpc'
  fopts = [
    '-Mpreprocess',
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    '-r8',
  ]

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
add_project_arguments(fopts, language: 'fortran')

# fix compiliation problems with of symmetry/symmetry_i.c
add_project_arguments('-D_Float128=__float128', language: 'c')

## ========================================== ##
## LIBRARIES
## ========================================== ##

if get_option('openmp')
  omp_dep = dependency('openmp', required: fc_id != 'intel' and fc_id != 'intel-llvm' and fc_id != 'nvidia_hpc')
  if not omp_dep.found()
    if fc_id == 'intel' or fc_id == 'intel-llvm'
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
  lib_deps += omp_dep
endif

lib_deps += dependency('threads')


lapack_vendor = get_option('lapack')
if lapack_vendor == 'auto'
  if fc_id == 'intel' or fc_id == 'intel-llvm'
    lapack_vendor = 'mkl'
  endif
endif

if lapack_vendor == 'mkl'
  mkl_dep = []
  if fc_id == 'intel' or fc_id == 'intel-llvm'
    mkl_dep += cc.find_library('mkl_intel_lp64')
    if get_option('openmp')
      mkl_dep += cc.find_library('mkl_intel_thread')
    endif
  elif fc_id == 'gcc'
    mkl_dep += cc.find_library('mkl_gf_lp64')
    if get_option('openmp')
      mkl_dep += cc.find_library('mkl_gnu_thread')
    endif
  else
    error('MKL not supported for this compiler')
  endif
  if not get_option('openmp')
    mkl_dep += cc.find_library('mkl_tbb_thread')
  endif
  mkl_dep += cc.find_library('mkl_core')
  lib_deps += mkl_dep

elif lapack_vendor == 'mkl-rt'
  mkl_dep = fc.find_library('mkl_rt')
  lib_deps += mkl_dep

elif lapack_vendor == 'openblas'
  openblas_dep = dependency('openblas', required: false)
  if not openblas_dep.found()
    openblas_dep = fc.find_library('openblas')
  endif
  lib_deps += openblas_dep
  if not fc.links('external dsytrs; call dsytrs(); end', dependencies: openblas_dep)
    lapack_dep = dependency('lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack')
    endif
    lib_deps += lapack_dep
  endif

elif lapack_vendor == 'custom'
  custom_deps = []
  libs = get_option('custom_libraries')
  if libs[0].startswith('-L')
    foreach lib: libs
      if lib != libs[0]
        custom_deps += cc.find_library(lib, dirs: libs[0].substring(2))
      endif
    endforeach
  else
    foreach lib: libs
      custom_deps += cc.find_library(lib)
    endforeach
  endif
  if (not fc.links('external dsytrs; call dsytrs(); end', dependencies:  (get_option('openmp') ? [custom_deps, omp_dep] : [custom_deps])))
    error('Custom LAPACK libraries do not link')
  elif (not fc.links('external dsytrs; call dgemm(); end', dependencies: (get_option('openmp') ? [custom_deps, omp_dep] : [custom_deps])))
    error('Custom BLAS libraries do not link')
  endif
  lib_deps += custom_deps

else
  lapack_dep = dependency('lapack', required: false)
  if not lapack_dep.found()
    lapack_dep = fc.find_library('lapack')
  endif
  lib_deps += lapack_dep
  blas_dep = dependency('blas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('blas')
  endif
  lib_deps += blas_dep
endif

if get_option('nvtx')
  lib_deps += fc.find_library('nvToolsExt', required: true)
endif

# Create the tool chain library as subproject
mctc_dep = dependency(
  'mctc-lib',
  fallback: ['mctc-lib', 'mctc_dep'],
  default_options: ['default_library=static'],
)
lib_deps += mctc_dep

# Get light-weight tight-binding framework dependency
tblite_dep = dependency(
  'tblite',
  fallback: ['tblite', 'tblite_dep'],
  default_options: ['default_library=static', 'api=false'],
  required: get_option('tblite')
)
lib_deps += tblite_dep

# Get multicharge dependency
multicharge_dep = dependency(
  'multicharge',
  fallback: ['multicharge', 'multicharge_dep'],
  required: get_option('tblite')
)
lib_deps += multicharge_dep

# Get DFTD dependency
dftd4_dep = dependency(
  'dftd4',
  fallback: ['dftd4', 'dftd4_dep'],
  required: get_option('tblite')
)
lib_deps += dftd4_dep

# Get CPCM-X
cpx_dep = dependency(
  'cpx',
  fallback: ['cpx', 'cpx_dep'],
  default_options: ['default_library=static'],
  required: get_option('cpcmx')
)
lib_deps += cpx_dep
