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

add_executable(
  "${PROJECT_NAME}_tests_capi"
  "c_api_example.c"
)
target_link_libraries(
  "${PROJECT_NAME}_tests_capi"
  PUBLIC
  "lib-${PROJECT_NAME}-shared"
)

add_test(
  "${PROJECT_NAME}/CAPI"
  "${PROJECT_NAME}_tests_capi"
)
set_tests_properties(
  "${PROJECT_NAME}/CAPI"
  PROPERTIES
  ENVIRONMENT XTBPATH=${PROJECT_SOURCE_DIR}
)
