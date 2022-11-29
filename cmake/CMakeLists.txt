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

option(WITH_OpenMP "Enable support for shared memory parallelisation with OpenMP" TRUE)
option(WITH_TBLITE "Use tblite library as backend for xTB" TRUE)
if(NOT DEFINED "${PROJECT_NAME}-dependeny-method")
  set(
    "${PROJECT_NAME}-dependency-method"
    "subproject" "cmake" "pkgconf" "fetch"
  )
endif()

set(
  module-dir
  "${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}"
)
set(module-dir "${module-dir}" PARENT_SCOPE)

# Set build type as CMake does not provide defaults
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(
    CMAKE_BUILD_TYPE "RelWithDebInfo"
    CACHE STRING "Build type to be used."
    FORCE
  )
  message(
    STATUS
    "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified."
  )
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/modules/"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fbacktrace")
  set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialect "-axAVX2 -r8 -traceback")
  set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  set(dialect "-Mbackslash -Mallocatable=03 -r8 -traceback")
endif()
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}" PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}" PARENT_SCOPE)

# Populate xtb_version.fh
set(version ${PROJECT_VERSION})
execute_process(COMMAND git rev-parse HEAD
  RESULT_VARIABLE git_return
  OUTPUT_VARIABLE commit)
if(git_return)
  set(commit "unknown-commit")
else()
  string(REGEX REPLACE "\n$" "" commit ${commit})
endif()
string(TIMESTAMP date "%m/%d/%Y")
set(author $ENV{USERNAME})
set(origin ${CMAKE_HOST_SYSTEM_NAME})
configure_file(
  "${PROJECT_SOURCE_DIR}/assets/templates/version.f90"
  "${PROJECT_BINARY_DIR}/xtb_version.fh"
  @ONLY
)
