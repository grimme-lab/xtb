cmake_minimum_required(VERSION 3.9)

# Setup the XTB Project
project(xtb
  VERSION 6.4.1
  LANGUAGES C Fortran
)
enable_testing()

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# Include CMake specific configurations
set(xtb-dir "${CMAKE_CURRENT_SOURCE_DIR}")
add_subdirectory("cmake")

# Sources: initialize program sources (prog) and library sources (srcs) empty
set(prog)
set(srcs)

add_subdirectory("src")
add_subdirectory("symmetry")

# Find dependencies
if (NOT DISABLE_OMP)
  find_package(OpenMP REQUIRED)
endif()
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

# Object library
set(xtb-mod "${CMAKE_CURRENT_BINARY_DIR}/xtb-mod")
add_library(xtb-object OBJECT ${srcs})
set_target_properties(xtb-object PROPERTIES
  Fortran_MODULE_DIRECTORY ${xtb-mod}
  POSITION_INDEPENDENT_CODE ON
)
target_include_directories(xtb-object
  PUBLIC
  ${xtb-config-dir}
  ${CMAKE_CURRENT_BINARY_DIR}
  $<BUILD_INTERFACE:${xtb-dir}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)

# Static Library
add_library(lib-xtb-static STATIC $<TARGET_OBJECTS:xtb-object> )
target_link_libraries(lib-xtb-static
  PUBLIC
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>
)
set_target_properties(lib-xtb-static PROPERTIES
  Fortran_MODULE_DIRECTORY ${xtb-mod}
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  POSITION_INDEPENDENT_CODE ON
  OUTPUT_NAME "${PROJECT_NAME}"
)
target_include_directories(lib-xtb-static
  PUBLIC
  $<BUILD_INTERFACE:${xtb-dir}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)

# Shared Library
add_library(lib-xtb-shared SHARED $<TARGET_OBJECTS:xtb-object> )
target_link_libraries(lib-xtb-shared
  PUBLIC
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>
)
set_target_properties(lib-xtb-shared PROPERTIES
  Fortran_MODULE_DIRECTORY ${xtb-mod}
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${XTB_VERSION_MAJOR}"
)
target_include_directories(lib-xtb-shared
  PUBLIC
  $<BUILD_INTERFACE:${xtb-dir}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)

# Executables
add_executable(xtb-exe ${prog})
target_link_libraries(xtb-exe
  PRIVATE
  lib-xtb-static
)
set_target_properties(xtb-exe PROPERTIES
  Fortran_MODULE_DIRECTORY ${xtb-mod}
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  OUTPUT_NAME "${PROJECT_NAME}"
)
target_include_directories(xtb-exe PRIVATE "${xtb-dir}/include")

# Install
install(FILES "${xtb-dir}/include/xtb.h"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(FILES
  "${xtb-dir}/param_gfn0-xtb.txt"
  "${xtb-dir}/param_gfn1-xtb.txt"
  "${xtb-dir}/param_gfn2-xtb.txt"
  "${xtb-dir}/param_ipea-xtb.txt"
  "${xtb-dir}/.param_gfnff.xtb"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}"
)
install(TARGETS lib-xtb-static lib-xtb-shared xtb-exe
  EXPORT xtbTargets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
# CMake package files
include(CMakePackageConfigHelpers)
# -- Config version file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/xtbConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)
# -- Config file
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/xtbConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
# -- Install config and configVersion
install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/xtbConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/xtbConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
# -- Targets file
# -- This makes the project importable from the build directory
export(EXPORT xtbTargets FILE "${CMAKE_CURRENT_BINARY_DIR}/xtbTargets.cmake")
# -- This makes the project importable from the install directory
install(EXPORT xtbTargets FILE "xtbTargets.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

add_subdirectory("TESTSUITE")
