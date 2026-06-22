# This file is part of xtb.
#
# Copyright (C) 2026 xtb GPU contributors
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

# GPU / high-throughput batch driver. These sources compile in the CPU build
# too: the batched eigensolver has a LAPACK fallback and the batch driver uses
# the standard single-point path; the cuSolver/OpenACC code is guarded by the
# USE_CUSOLVER / _OPENACC macros set only by the GPU build (see top-level
# CMakeLists.txt WITH_GPU/WITH_CUSOLVER options).
set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/batch_capture.F90"
  "${dir}/batched_eig.F90"
  "${dir}/batch_driver.F90"
)

set(srcs ${srcs} PARENT_SCOPE)
