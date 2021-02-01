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
  "${dir}/basisset.f90"
  "${dir}/calculator.f90"
  "${dir}/coulomb.f90"
  "${dir}/data.f90"
  "${dir}/dispersion.f90"
  "${dir}/eeq.f90"
  "${dir}/gfn0.f90"
  "${dir}/gfn1.f90"
  "${dir}/gfn2.f90"
  "${dir}/halogen.f90"
  "${dir}/hamiltonian.f90"
  "${dir}/hamiltonian_gpu.f90"
  "${dir}/multipole.f90"
  "${dir}/repulsion.F90"
  "${dir}/thirdorder.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
