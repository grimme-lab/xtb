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

add_subdirectory("blas")
add_subdirectory("lapack")

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

# MCTC library
list(APPEND srcs
  "${dir}/accuracy.f90"
  "${dir}/boundaryconditions.f90"
  "${dir}/chartools.f90"
  "${dir}/convert.f90"
  "${dir}/filetypes.f90"
  "${dir}/io.f90"
  "${dir}/math.f90"
  "${dir}/namegen.f90"
  "${dir}/resize.f90"
  "${dir}/runtypes.f90"
  "${dir}/search.f90"
  "${dir}/sort.f90"
  "${dir}/symbols.f90"
  "${dir}/systools.F90"
  "${dir}/thresholds.f90"
  "${dir}/version.f90"
  "${dir}/mctc_global.F90"
  "${dir}/mctc_strings.f90"
  "${dir}/mctc_constants.f90"
  "${dir}/mctc_param.f90"
  "${dir}/param/atomic_masses.f90"
  "${dir}/param/chemical_hardnesses.f90"
  "${dir}/param/covalent_radii.f90"
  "${dir}/param/electronegativities.f90"
  "${dir}/param/pse.f90"
  "${dir}/param/r4r2_expectation_values.f90"
  "${dir}/mctc_timings.f90"
  "${dir}/mctc_filetools.f90"
  "${dir}/linalg.f90"
  "${dir}/lapack.f90"
  "${dir}/blas.f90"
  "${dir}/mctc_init.F90"
  "${dir}/error.f90"
  "${dir}/signal.c"
)

set(srcs ${srcs} PARENT_SCOPE)
