# Thi file is part of xtb.
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

# Some user-configurable features
option(WITH_OpenMP "Enable support for shared memory parallelisation with OpenMP" TRUE)
option(WITH_TBLITE "Use tblite library as backend for xTB" TRUE)
option(WITH_CPCMX "Use CPCM-X solvation library for xTB" TRUE)


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

# Add modules to the CMake build 
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)

# specify module installation directory
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/modules/"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# Compiler-specific configurations 
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialects "-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fbacktrace")
  set(bounds "-fbounds-check -ffpe-trap=invalid,zero,overflow")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialects "-axAVX2 -r8 -traceback")
  set(bounds "-check all -fpe0")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  set(dialects "-Mbackslash -Mallocatable=03 -r8 -traceback")
endif()

# Customize compiler flags
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}" PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialects}" PARENT_SCOPE)


# Populate xtb_version.fh with metadata
set(version ${PROJECT_VERSION})
execute_process(COMMAND git rev-parse HEAD
  RESULT_VARIABLE git_return
  OUTPUT_VARIABLE commit)
if(git_return)
  set(commit "unknown-commit")
else()
  string(REGEX REPLACE "\n$" "" commit ${commit})
endif()
string(TIMESTAMP date "%Y/%m/%d")
set(author $ENV{USERNAME})
set(origin ${CMAKE_HOST_SYSTEM_NAME})
configure_file(
  "${PROJECT_SOURCE_DIR}/assets/templates/version.f90"
  "${PROJECT_BINARY_DIR}/xtb_version.fh"
  @ONLY
)
