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

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/atomicmass.f90"
  "${dir}/atomicrad.f90"
  "${dir}/covalentradd3.f90"
  "${dir}/paulingen.f90"
  "${dir}/sqrtzr4r2.f90"
  "${dir}/vdwradd3.f90"
  "${dir}/vdwradcosmo.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
