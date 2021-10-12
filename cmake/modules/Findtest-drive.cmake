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

#[[.rst:
Find test-drive
---------------

Makes the test-drive project available.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``test-drive::test-drive``
  The test-drive library


Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variables:

``TEST_DRIVE_FOUND``
  True if the test-drive library is available

``TEST_DRIVE_SOURCE_DIR``
  Path to the source directory of the test-drive project,
  only set if the project is included as source.

``TEST_DRIVE_BINARY_DIR``
  Path to the binary directory of the test-drive project,
  only set if the project is included as source.

Cache variables
^^^^^^^^^^^^^^^

The following cache variables may be set to influence the library detection:

``TEST_DRIVE_FIND_METHOD``
  Methods to find or make the project available. Available methods are
  - ``cmake``: Try to find via CMake config file
  - ``pkgconf``: Try to find via pkg-config file
  - ``subproject``: Use source in subprojects directory
  - ``fetch``: Fetch the source from upstream

``TEST_DRIVE_DIR``
  Used for searching the CMake config file

``TEST_DRIVE_SUBPROJECT``
  Directory to find the test-drive subproject, relative to the project root

#]]

set(_lib "test-drive")
set(_pkg "TEST_DRIVE")
set(_url "https://github.com/fortran-lang/test-drive")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  if(DEFINED "${PROJECT_NAME}-dependency-method")
    set("${_pkg}_FIND_METHOD" "${${PROJECT_NAME}-dependency-method}")
  else()
    set("${_pkg}_FIND_METHOD" "cmake" "pkgconf" "subproject" "fetch")
  endif()
  set("_${_pkg}_FIND_METHOD")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/xtb-utils.cmake")

xtb_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}")

if(DEFINED "_${_pkg}_FIND_METHOD")
  unset("${_pkg}_FIND_METHOD")
  unset("_${_pkg}_FIND_METHOD")
endif()
unset(_lib)
unset(_pkg)
unset(_url)
