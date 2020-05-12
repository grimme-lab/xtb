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

if get_option('static') and fc.get_id() == 'intel'
  if get_option('la_backend') != 'mkl'
    error('Static linked binaries are only supported with Intel MKL')
  endif
endif

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
  if get_option('static')
    add_project_link_arguments('-static', language: 'fortran')
    add_project_link_arguments('-static', language: 'c') # icc will do linking
  endif
elif fc.get_id() == 'pgi'
  add_project_arguments(
    '-Mbackslash', '-Mallocatable=03', '-traceback', '-r8',
    language: 'fortran'
  )
endif

if cc.get_id() == 'gcc'
  add_project_arguments('-std=c11', language: 'c')
elif cc.get_id() == 'intel'
  add_project_arguments('-std=c11', language: 'c')
endif

# fix compiliation problems with of symmetry/symmetry_i.c
add_project_arguments('-D_Float128=__float128', language: 'c')

## ========================================== ##
## LIBRARIES
## ========================================== ##

if get_option('la_backend') == 'mkl'
  add_project_arguments('-DWITH_MKL', language: 'fortran')

  mkl_rt_dep = fc.find_library('mkl_rt', required: true)

  if get_option('static')
    if fc.get_id() == 'gcc'
      libmkl_exe = [fc.find_library('mkl_gf_lp64')]
      libmkl_exe += fc.find_library('mkl_gnu_thread')
    elif fc.get_id() == 'intel'
      libmkl_exe = [fc.find_library('mkl_intel_lp64')]
      libmkl_exe += fc.find_library('mkl_intel_thread')
    elif fc.get_id() == 'pgi'
      libmkl_exe = [fc.find_library('mkl_intel_lp64')]
      libmkl_exe += fc.find_library('mkl_pgi_thread')
    endif
    libmkl_exe += fc.find_library('mkl_core')
    libmkl_sha = [fc.find_library('mkl_rt')]
    dependencies_sha += mkl_rt_dep
    dependencies_exe += libmkl_exe
  else
    dependencies += mkl_rt_dep
  endif

elif get_option('la_backend') == 'openblas'
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

elif get_option('la_backend') == 'custom'
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
  if fc.get_id() == 'gcc'
    add_project_arguments('-fopenmp', language: 'fortran')
    add_project_link_arguments('-fopenmp', language: 'fortran')
    add_project_link_arguments('-fopenmp', language: 'c')
  elif fc.get_id() == 'intel'
    add_project_arguments('-qopenmp', language: 'fortran')
    add_project_link_arguments('-qopenmp', language: 'fortran')
    add_project_link_arguments('-qopenmp', language: 'c')
  endif
endif

dependencies += dependency('threads')

# distribute dependencies for shared object and static executable
dependencies_sha += dependencies
dependencies_exe += dependencies
