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

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND prog
  "${dir}/argparser.f90"
  "${dir}/dock.f90"
  "${dir}/info.f90"
  "${dir}/irmod.f90"
  "${dir}/main.F90"
  "${dir}/primary.f90"
  "${dir}/submodules.f90"
  "${dir}/thermo.f90"
  "${dir}/topology.f90"
)

set(prog ${prog} PARENT_SCOPE)
