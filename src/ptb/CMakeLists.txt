# This file is part of xtb.
#
# Copyright (C) 2022 Marcel Müller 
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
  "${dir}/calculator.F90"
  "${dir}/param.F90"
  "${dir}/data.f90"
  "${dir}/vdzp.F90"
  "${dir}/twostepscf.F90"
  "${dir}/mmlpopanalysis.F90"
  "${dir}/ncoord.F90"
  "${dir}/corebasis.F90"
  "${dir}/hamiltonian.F90"
  "${dir}/guess.F90"
  "${dir}/paulixc.F90"
  "${dir}/integral_types.F90"
  "${dir}/integrals.F90"
  "${dir}/coulomb.F90"
  "${dir}/plusu.F90"
  "${dir}/property.F90"
  "${dir}/response.F90"
  "${dir}/solver.F90"
)

set(srcs ${srcs} PARENT_SCOPE)
