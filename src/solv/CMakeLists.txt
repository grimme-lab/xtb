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

add_subdirectory("ddcosmo")

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/born.f90"
  "${dir}/cm5.f90"
  "${dir}/cosmo.f90"
  "${dir}/gbsa.f90"
  "${dir}/input.F90"
  "${dir}/kernel.f90"
  "${dir}/lebedev.f90"
  "${dir}/model.f90"
  "${dir}/sasa.f90"
  "${dir}/state.f90"
  "${dir}/cpx.F90"
  "${dir}/ddvolume.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
