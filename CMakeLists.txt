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

cmake_minimum_required(VERSION 3.17)

# Setup the xtb Project
project(
  "xtb"
  VERSION "6.7.1"
  LANGUAGES "C" "Fortran")

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(WITH_OpenMP
       "Enable support for shared memory parallelisation with OpenMP" ON)
option(WITH_TBLITE "Use tblite library as backend for xTB" ON)
option(WITH_CPCMX "Use CPCM-X solvation library for xTB" ON)
option(WITH_TESTS "Enable compilation of unit tests" ON)

# Follow GNU conventions for installing directories
include(GNUInstallDirs)
set(CMAKE_INSTALL_Fortran_MODULES
    "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}"
    CACHE STRING "Installation directory for xTB Fortran module files")

# Set a default for single-configuration generators
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Build type to be used." FORCE)
  message(
    STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified.")
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
include("${PROJECT_SOURCE_DIR}/cmake/modules/xtb-utils.cmake")

# Compiler-specific options, applied only to xTB targets
set(xtb-fortran-options)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(xtb-fortran-options
      -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fbacktrace
      "$<$<CONFIG:Debug>:-fbounds-check>")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(xtb-fortran-options -axAVX2 -r8 -traceback
                          "$<$<CONFIG:Debug>:-check;all>")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  set(xtb-fortran-options -Mbackslash -Mallocatable=03 -r8 -traceback)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
  set(xtb-fortran-options -fdefault-real-8 -fdefault-double-8)
endif()

# Populate xtb_version.fh with metadata
set(version "${PROJECT_VERSION}")
execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  RESULT_VARIABLE git_return
  OUTPUT_VARIABLE commit
  OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(git_return)
  set(commit "unknown-commit")
endif()
string(TIMESTAMP date "%Y/%m/%d")
set(author "$ENV{USERNAME}")
set(origin "${CMAKE_HOST_SYSTEM_NAME}")
configure_file("${PROJECT_SOURCE_DIR}/assets/templates/version.f90"
               "${PROJECT_BINARY_DIR}/include/xtb_version.fh" @ONLY)

# Check a specific CMake targets for xtb build & execute corresponding Find
# scripts
if(NOT TARGET "mctc-lib::mctc-lib")
  find_package("mctc-lib" REQUIRED)
endif()

if(NOT TARGET "tblite::tblite" AND WITH_TBLITE)
  find_package("tblite" REQUIRED)
endif()

foreach(
  _xtb_dep IN
  ITEMS toml-f
        test-drive
        mstore
        multicharge
        dftd4
        s-dftd3
        tblite)
  xtb_promote_target("${_xtb_dep}")
endforeach()
unset(_xtb_dep)

if(NOT TARGET "cpcmx::cpcmx" AND WITH_CPCMX)
  find_package("cpcmx" REQUIRED)
endif()

foreach(_xtb_dep IN ITEMS cpcmx numsa toml-f test-drive mctc-lib)
  xtb_promote_target("${_xtb_dep}")
endforeach()
unset(_xtb_dep)

if(WITH_TESTS AND NOT TARGET "test-drive::test-drive")
  find_package("test-drive" REQUIRED)
endif()

# Sources: initialize program sources (prog) and library sources (srcs) empty
set(prog)
set(srcs)

# add filenames to the list variables
add_subdirectory("src")
add_subdirectory("symmetry")

# CMake modules for third-party software
if(NOT TARGET "OpenMP::OpenMP_Fortran" AND WITH_OpenMP)
  find_package("OpenMP" REQUIRED)
endif()
find_package("LAPACK" REQUIRED)

# Library type is selected by BUILD_SHARED_LIBS
add_library("${PROJECT_NAME}-lib" ${srcs})
add_library("${PROJECT_NAME}::${PROJECT_NAME}" ALIAS "${PROJECT_NAME}-lib")

target_link_libraries(
  "${PROJECT_NAME}-lib"
  PUBLIC ${LAPACK_LIBRARIES} "$<$<BOOL:${WITH_OpenMP}>:OpenMP::OpenMP_Fortran>"
         "mctc-lib::mctc-lib" "$<$<BOOL:${WITH_CPCMX}>:cpcmx::cpcmx>"
         "$<$<BOOL:${WITH_TBLITE}>:tblite::tblite>")

target_compile_options(
  "${PROJECT_NAME}-lib"
  PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:${xtb-fortran-options}>")
target_compile_definitions(
  "${PROJECT_NAME}-lib" PUBLIC "$<$<BOOL:${WITH_TBLITE}>:WITH_TBLITE>"
                               "$<$<BOOL:${WITH_CPCMX}>:WITH_CPCMX>")

set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
             POSITION_INDEPENDENT_CODE ON
             OUTPUT_NAME "${PROJECT_NAME}"
             EXPORT_NAME "${PROJECT_NAME}"
             VERSION "${PROJECT_VERSION}"
             SOVERSION "${PROJECT_VERSION_MAJOR}")

target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_Fortran_MODULES}>)

# ##############################################################################
# Executables #
# ##############################################################################

add_executable(${PROJECT_NAME}-exe ${prog})

target_link_libraries(${PROJECT_NAME}-exe
                      PRIVATE "${PROJECT_NAME}::${PROJECT_NAME}")

target_compile_options(
  "${PROJECT_NAME}-exe"
  PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:${xtb-fortran-options}>")

set_target_properties(
  ${PROJECT_NAME}-exe
  PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include
             OUTPUT_NAME "${PROJECT_NAME}")

# ##############################################################################
# Installation rules #
# ##############################################################################
# API for C/C++
install(FILES "${PROJECT_SOURCE_DIR}/include/xtb.h" TYPE INCLUDE)

install(
  DIRECTORY "${PROJECT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_Fortran_MODULES}"
  FILES_MATCHING
  PATTERN "*.mod"
  PATTERN "*.smod")

# xtb-parameters
install(
  FILES "${PROJECT_SOURCE_DIR}/param_gfn0-xtb.txt"
        "${PROJECT_SOURCE_DIR}/param_gfn1-xtb.txt"
        "${PROJECT_SOURCE_DIR}/param_gfn1-si-xtb.txt"
        "${PROJECT_SOURCE_DIR}/param_gfn2-xtb.txt"
        "${PROJECT_SOURCE_DIR}/param_ipea-xtb.txt"
        "${PROJECT_SOURCE_DIR}/.param_gfnff.xtb"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}")

# Build output artifacts
install(TARGETS "${PROJECT_NAME}-lib" "${PROJECT_NAME}-exe"
        EXPORT "${PROJECT_NAME}-targets")

# Find modules used by the installed package configuration
install(DIRECTORY "${PROJECT_SOURCE_DIR}/cmake/modules/"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

# CMake package files
include(CMakePackageConfigHelpers)

# -- Config version file
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY AnyNewerVersion)

# -- Config file
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

# -- Install config and configVersion
install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

# -- Targets file
install(
  EXPORT "${PROJECT_NAME}-targets"
  NAMESPACE "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

# Tests
if(WITH_TESTS)
  enable_testing()
  add_subdirectory("test")
endif()
