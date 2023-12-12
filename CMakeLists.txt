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
option(WITH_OBJECT "To build using object library" TRUE)
option(INSTALL_MODULES "Install Fortran module files to include directory." FALSE)

# Buggy CMake versions
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.27.0)
   set(WITH_OBJECT FALSE)
endif()

# Setup the xtb Project
project(
  "xtb"
  VERSION "6.6.1" 
  LANGUAGES "C" "Fortran"
)

enable_testing()

# Follow GNU conventions for installing directories
include(GNUInstallDirs)


# Include CMake specific configurations
add_subdirectory("cmake")

# Check a specific CMake targets for xtb build & execute corresponding Find scripts
if(NOT TARGET "mctc-lib::mctc-lib")
  find_package("mctc-lib" REQUIRED)
endif()

if(NOT TARGET "tblite::tblite" AND WITH_TBLITE)
  find_package("tblite" REQUIRED)
endif()

if(NOT TARGET "cpcmx::cpcmx" AND WITH_CPCMX)
  find_package("cpcmx" REQUIRED)
endif()

if(NOT TARGET "test-drive::test-drive")
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
find_package("BLAS" REQUIRED)


if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
   file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
endif()

##################
# OBJECT LIBRARY #
##################
if(WITH_OBJECT)
   
   add_library(
      "${PROJECT_NAME}-object"
      OBJECT
      ${srcs}
   )


   # customize object library
   set_target_properties(
      "${PROJECT_NAME}-object"
      PROPERTIES
      Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
      POSITION_INDEPENDENT_CODE ON
   )

   # link object library against mctc-lib &
   # & conditionally against OpenMP, tblite, cpcmx
   target_link_libraries(
      "${PROJECT_NAME}-object"
      PUBLIC
      "mctc-lib::mctc-lib"
      "$<$<BOOL:${WITH_CPCMX}>:cpcmx::cpcmx>"
      "$<$<BOOL:${WITH_TBLITE}>:tblite::tblite>"
      "$<$<BOOL:${WITH_OpenMP}>:OpenMP::OpenMP_Fortran>"
   )

   # include directories
   target_include_directories(
      "${PROJECT_NAME}-object"
      PUBLIC
      ${xtb-config-dir}
      ${PROJECT_BINARY_DIR}
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
   )

endif()

##################
# Static Library #
##################
if(WITH_OBJECT)
   add_library(
      "lib-${PROJECT_NAME}-static"
      STATIC
      $<TARGET_OBJECTS:${PROJECT_NAME}-object>
   )
else()
   add_library(
      "lib-${PROJECT_NAME}-static"
      STATIC
      ${srcs}
   )
endif()


target_link_libraries(
   "lib-${PROJECT_NAME}-static"
   PUBLIC
   ${BLAS_LIBRARIES}
   ${LAPACK_LIBRARIES}
   "$<$<BOOL:${WITH_OpenMP}>:OpenMP::OpenMP_Fortran>"
   "mctc-lib::mctc-lib"
   "$<$<BOOL:${WITH_CPCMX}>:cpcmx::cpcmx>"
   "$<$<BOOL:${WITH_TBLITE}>:tblite::tblite>"
)


 set_target_properties(
   "lib-${PROJECT_NAME}-static"
   PROPERTIES
   Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
   ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
   POSITION_INDEPENDENT_CODE ON
   OUTPUT_NAME "${PROJECT_NAME}"
 )


target_include_directories(
  "lib-${PROJECT_NAME}-static"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)


##################
# Shared Library #
##################
if (WITH_OBJECT) 
   add_library(
      "lib-${PROJECT_NAME}-shared"
      SHARED
      $<TARGET_OBJECTS:${PROJECT_NAME}-object>
   )

   target_link_libraries(
      "lib-${PROJECT_NAME}-shared"
      PUBLIC
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES}
      "$<$<BOOL:${WITH_OpenMP}>:OpenMP::OpenMP_Fortran>"
      "mctc-lib::mctc-lib"
      "$<$<BOOL:${WITH_CPCMX}>:cpcmx::cpcmx>"
      "$<$<BOOL:${WITH_TBLITE}>:tblite::tblite>"
   )

   set_target_properties(
      "lib-${PROJECT_NAME}-shared"
      PROPERTIES
      Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
      LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
      OUTPUT_NAME "${PROJECT_NAME}"
      VERSION "${PROJECT_VERSION}"
      SOVERSION "${PROJECT_VERSION_MAJOR}"
   )

   target_include_directories(
      "lib-${PROJECT_NAME}-shared"
      PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
   )
endif()

###############
# Executables #
###############

add_executable(
  ${PROJECT_NAME}-exe
  ${prog}
)

target_link_libraries(
   ${PROJECT_NAME}-exe
   PRIVATE
   "lib-${PROJECT_NAME}-static"
)
set_target_properties(
  ${PROJECT_NAME}-exe
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  OUTPUT_NAME "${PROJECT_NAME}"
)

target_include_directories(
   ${PROJECT_NAME}-exe 
   PRIVATE 
   ${PROJECT_SOURCE_DIR}/include
   $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

######################
# Installation rules #
######################
# API for C/C++
install(
   FILES
   "${PROJECT_SOURCE_DIR}/include/xtb.h"
   DESTINATION
   "${CMAKE_INSTALL_INCLUDEDIR}"
)

if (INSTALL_MODULES)
   install(
      DIRECTORY
      "${PROJECT_BINARY_DIR}/include/"
      DESTINATION
      "${CMAKE_INSTALL_INCLUDEDIR}"
   )
endif()

# xtb-parameters 
install(
   FILES
   "${PROJECT_SOURCE_DIR}/param_gfn0-xtb.txt"
   "${PROJECT_SOURCE_DIR}/param_gfn1-xtb.txt"
   "${PROJECT_SOURCE_DIR}/param_gfn1-si-xtb.txt"
   "${PROJECT_SOURCE_DIR}/param_gfn2-xtb.txt"
   "${PROJECT_SOURCE_DIR}/param_ipea-xtb.txt"
   "${PROJECT_SOURCE_DIR}/.param_gfnff.xtb"
   DESTINATION
   "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}"
)

# Build output artifacts
if (WITH_OBJECT)
   install(
      TARGETS
      "lib-${PROJECT_NAME}-static"
      "lib-${PROJECT_NAME}-shared"
      "${PROJECT_NAME}-exe"
      EXPORT "${PROJECT_NAME}-targets"
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
      INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"   
   )
else()
   install(
      TARGETS
      "lib-${PROJECT_NAME}-static"
      "${PROJECT_NAME}-exe"
      EXPORT "${PROJECT_NAME}-targets"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
      INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
   )
endif()

# CMake package files
include(CMakePackageConfigHelpers)

# -- Config version file
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY AnyNewerVersion
)

# -- Config file
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# -- Install config and configVersion
install(
  FILES
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# -- Targets file
# -- This makes the project importable from the build directory
export(
   EXPORT "${PROJECT_NAME}-targets"
   FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-targets.cmake"
)

# -- This makes the project importable from the install directory
install(
   EXPORT "${PROJECT_NAME}-targets"
   FILE "${PROJECT_NAME}-targets.cmake"
   DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
   
# Tests
add_subdirectory("test")
