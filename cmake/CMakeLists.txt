cmake_minimum_required(VERSION 3.9)
enable_language(Fortran)

# Load lists of source files
include(Files.cmake)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#    set(dialect "-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fbacktrace")
    set(dialect "-ffree-line-length-none -fbacktrace")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-axAVX2 -r8 -traceback")
    set(bounds "-check bounds")
endif()
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

# TODO
set(version "6.2.2")
set(commit "dummy")
set(date "dummy")
set(author "author@origin")
configure_file("${XTB_ROOT}/xtb/version.f90" "${XTB_ROOT}/xtb/xtb_version.fh" @ONLY)

# Find dependencies
find_package(OpenMP REQUIRED)
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
  $<BUILD_INTERFACE:${XTB_ROOT}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/xtb>
)

# Static Library
add_library(xtb_static STATIC $<TARGET_OBJECTS:xtb_object> )
target_link_libraries(xtb_static
 PUBLIC
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
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
)
set_target_properties(xtb PROPERTIES
  Fortran_MODULE_DIRECTORY xtb-mod
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
target_include_directories(xtb PRIVATE "${XTB_ROOT}/include")

# Install
install(FILES "${XTB_ROOT}/include/xtb.h"
        DESTINATION include/xtb/xtb.h
)
install(TARGETS xtb_static xtb_shared xtb
  EXPORT xtbTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
)

# Build Tests
# TODO
