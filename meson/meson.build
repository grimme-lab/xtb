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
  if fc.get_id() == 'gcc'
    libmkl_exe = [fc.find_library('mkl_gf_lp64')]
    libmkl_exe += fc.find_library('mkl_gnu_thread')
  elif fc.get_id() == 'intel'
    libmkl_exe = [fc.find_library('mkl_intel_lp64')]
    libmkl_exe += fc.find_library('mkl_intel_thread')
  endif
  libmkl_exe += fc.find_library('mkl_core')
  libmkl_sha = [fc.find_library('mkl_rt')]
  dependencies_sha += libmkl_sha
  dependencies_exe += libmkl_exe
elif get_option('la_backend') == 'openblas'
  dependencies += fc.find_library('openblas', required: true)
  dependencies += fc.find_library('lapack', required: true)
elif get_option('la_backend') == 'custom'
  foreach lib: get_option('custom_libraries')
    dependencies += fc.find_library(lib)
  endforeach
else
  dependencies += fc.find_library('blas', required: true)
  dependencies += fc.find_library('lapack', required: true)
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
