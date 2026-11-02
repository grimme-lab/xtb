# This file is part of xtb.
# SPDX-Identifier: LGPL-3.0-or-later
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

# ddX is a transitive dependency of tblite. tblite v0.8.0 still pins ddX v1.0.0,
# which produces NaN values in solvation calculations (fixed in ddX v1.0.2).
# Providing the ddx::ddx target here prevents tblite from retrieving its own copy.
set(_lib "ddx")
set(_pkg "DDX")
set(_url "https://github.com/ddsolvation/ddX")
set(_rev "v1.0.2")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
   set("${_pkg}_FIND_METHOD" "cmake" "subproject" "fetch")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/xtb-utils.cmake")

# The installed config file is named ddXConfig.cmake, which a lookup for the
# lowercase package name cannot find; check it explicitly to enforce the version.
if(NOT TARGET "ddx::ddx" AND "cmake" IN_LIST "${_pkg}_FIND_METHOD")
   find_package("ddX" 1.0.2 CONFIG QUIET)
endif()

set(DDX_BUILD_EXAMPLES OFF)
if(NOT WITH_TESTS)
   set(BUILD_TESTING OFF)
endif()

xtb_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}" "${_rev}")

unset(DDX_BUILD_EXAMPLES)
unset(BUILD_TESTING)
unset(_lib)
unset(_pkg)
unset(_url)
unset(_rev)
