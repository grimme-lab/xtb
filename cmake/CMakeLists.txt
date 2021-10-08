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

cmake_minimum_required(VERSION 3.9)

if(NOT DEFINED "${PROJECT_NAME}-dependeny-method")
  set(
    "${PROJECT_NAME}-dependency-method"
    "subproject" "cmake" "pkgconf" "fetch"
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
  if(NOT CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8.0.0)
    set(DISABLE_OMP ON PARENT_SCOPE)
  endif()
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


# Program version, split to get major version for SOVERSION
string(REGEX MATCHALL "[0-9]+" version-split "${PROJECT_VERSION}")
list(GET version-split 0 XTB_VERSION_MAJOR)
set(XTB_VERSION_MAJOR "${XTB_VERSION_MAJOR}" PARENT_SCOPE)


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
  "${xtb-dir}/assets/templates/version.f90"
  "${CMAKE_CURRENT_BINARY_DIR}/xtb_version.fh"
  @ONLY
)

set(xtb-config-dir "${CMAKE_CURRENT_BINARY_DIR}" PARENT_SCOPE)
