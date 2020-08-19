cmake_minimum_required(VERSION 3.9)

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
