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
  "${dir}/anc.f90"
  "${dir}/atomlist.f90"
  "${dir}/basisset.f90"
  "${dir}/buffer.f90"
  "${dir}/calculator.f90"
  "${dir}/coulomb.f90"
  "${dir}/data.f90"
  "${dir}/dispersion_model.f90"
  "${dir}/environment.f90"
  "${dir}/fragments.f90"
  "${dir}/identitymap.f90"
  "${dir}/iohandler.f90"
  "${dir}/molecule.f90"
  "${dir}/neighbourlist.f90"
  "${dir}/latticepoint.f90"
  "${dir}/options.f90"
  "${dir}/param.f90"
  "${dir}/pcem.f90"
  "${dir}/reader.F90"
  "${dir}/restart.F90"
  "${dir}/setvar.f90"
  "${dir}/solvation.f90"
  "${dir}/timer.f90"
  "${dir}/topology.f90"
  "${dir}/wavefunction.f90"
  "${dir}/wignerseitzcell.f90"
  "${dir}/wsc.f90"
  "${dir}/vendordata.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
