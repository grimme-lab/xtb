cmake_minimum_required(VERSION 3.9)
enable_language(Fortran)

# Load lists of source files
include(Files.cmake)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fbacktrace")
    set(bounds "-fbounds-check")
    if(NOT CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8.0.0)
      set(DISABLE_OMP ON)
    endif()
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-axAVX2 -r8 -traceback")
    set(bounds "-check bounds")
endif()
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

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
configure_file("${XTB_ROOT}/xtb/version.f90" "${CMAKE_CURRENT_BINARY_DIR}/xtb_version.fh" @ONLY)

# Find dependencies
if (NOT DISABLE_OMP)
  find_package(OpenMP REQUIRED)
endif()
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

# Object library
add_library(xtb_object OBJECT ${XTB_SOURCES})
target_link_libraries(xtb_object
 PUBLIC
  $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>
)
set_target_properties(xtb_object PROPERTIES
  Fortran_MODULE_DIRECTORY xtb-mod
  POSITION_INDEPENDENT_CODE ON
)
target_include_directories(xtb_object
 PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}
  $<BUILD_INTERFACE:${XTB_ROOT}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/xtb>
)

# Static Library
add_library(xtb_static STATIC $<TARGET_OBJECTS:xtb_object> )
target_link_libraries(xtb_static
 PUBLIC
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>
)
set_target_properties(xtb_static PROPERTIES
  Fortran_MODULE_DIRECTORY xtb-mod
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  POSITION_INDEPENDENT_CODE ON
)
target_include_directories(xtb_static
 PUBLIC
  $<BUILD_INTERFACE:${XTB_ROOT}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/xtb>
)

# Shared Library
add_library(xtb_shared SHARED $<TARGET_OBJECTS:xtb_object> )
target_link_libraries(xtb_shared
 PUBLIC
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  $<$<BOOL:${OpenMP_Fortran_FOUND}>:OpenMP::OpenMP_Fortran>
)
set_target_properties(xtb_shared PROPERTIES
  Fortran_MODULE_DIRECTORY xtb-mod
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
target_include_directories(xtb_shared
 PUBLIC
  $<BUILD_INTERFACE:${XTB_ROOT}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/xtb>
)

# Executables
add_executable(xtb "${XTB_ROOT}/xtb/program_main.f90")
target_link_libraries(xtb
 PRIVATE
  xtb_static
)
set_target_properties(xtb PROPERTIES
  Fortran_MODULE_DIRECTORY xtb-mod
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
target_include_directories(xtb PRIVATE "${XTB_ROOT}/include")

# Install
install(FILES "${XTB_ROOT}/include/xtb.h"
        DESTINATION include/xtb/
)
install(FILES
 "${XTB_ROOT}/.param_gfn0.xtb"
 "${XTB_ROOT}/.param_gfn.xtb"
 "${XTB_ROOT}/.param_gfn2.xtb"
 DESTINATION share/xtb/
)
install(TARGETS xtb_static xtb_shared xtb
  EXPORT xtbTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
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
  "${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/xtbConfig.cmake"
  INSTALL_DESTINATION share/cmake/xtb
)
# -- Install config and configVersion
install(
 FILES
  "${CMAKE_CURRENT_BINARY_DIR}/xtbConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/xtbConfigVersion.cmake"
 DESTINATION share/cmake/xtb
)
# -- Targets file
# -- This makes the project importable from the build directory
export(EXPORT xtbTargets FILE "${CMAKE_CURRENT_BINARY_DIR}/xtbTargets.cmake")
# -- This makes the project importable from the install directory
install(EXPORT xtbTargets FILE "xtbTargets.cmake"
  DESTINATION share/cmake/xtb
)


# Build Tests
add_executable(xtb_tests_static ${XTB_F_TEST_SOURCES})
target_link_libraries(xtb_tests_static
 PUBLIC
  xtb_static
)
set_target_properties(xtb_tests_static PROPERTIES
  Fortran_MODULE_DIRECTORY xtb-mod
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
add_executable(xtb_tests_shared ${XTB_F_TEST_SOURCES})
target_link_libraries(xtb_tests_shared
 PUBLIC
  xtb_shared
)
set_target_properties(xtb_tests_shared PROPERTIES
  Fortran_MODULE_DIRECTORY xtb-mod
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
add_executable(xtb_tests_capi ${XTB_C_TEST_SOURCES})
target_link_libraries(xtb_tests_capi
 PUBLIC
  xtb_shared
)
set_target_properties(xtb_tests_capi PROPERTIES
  Fortran_MODULE_DIRECTORY xtb-mod
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
include(Tests.cmake)
