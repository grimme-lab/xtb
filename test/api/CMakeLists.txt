# This file is part of xtb.
# SPDX-Identifier: LGPL-3.0-or-later
#
# xtb is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with xtb.  If not, see <https://www.gnu.org/licenses/>.

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(xtb_tests_capi "${dir}/c_api_example.c")
target_link_libraries(xtb_tests_capi
  PUBLIC
  lib-xtb-shared
)
set_target_properties(xtb_tests_capi PROPERTIES
  Fortran_MODULE_DIRECTORY ${xtb-mod}
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)

set(XTB-CTESTER ${PROJECT_BINARY_DIR}/xtb_tests_capi)

add_test("${PROJECT_NAME}/CAPI" ${XTB-CTESTER})
set_tests_properties(
   "${PROJECT_NAME}/CAPI"
  PROPERTIES
  ENVIRONMENT XTBPATH=${PROJECT_SOURCE_DIR}
)
