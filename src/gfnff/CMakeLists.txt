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

list(APPEND srcs
  "${dir}/calculator.f90"
  "${dir}/data.f90"
  "${dir}/frag_hess.f90"
  "${dir}/gdisp0.f90"
  "${dir}/generator.f90"
  "${dir}/gfnff_eg.f90"
  "${dir}/gfnff_ini.f90"
  "${dir}/gfnff_ini2.f90"
  "${dir}/gfnff_param.f90"
  "${dir}/gfnff_qm.f90"
  "${dir}/gfnff_rab.f"
  "${dir}/gfnff_setup.f90"
  "${dir}/math.f"
  "${dir}/mrec.f90"
  "${dir}/neighbourlist.f90"
  "${dir}/shake_module.f90"
  "${dir}/struc_converter.f90"
  "${dir}/topology.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
