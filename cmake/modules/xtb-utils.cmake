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

# Handling of subproject dependencies
macro(
   "xtb_find_package"
   package
   methods
   url
   rev
)
string(TOLOWER "${package}" _pkg_lc)
string(TOUPPER "${package}" _pkg_uc)

# iterate through all methods
foreach(method ${methods})

   if(TARGET "${package}::${package}")
      break()
   endif()

   # cmake case
   if("${method}" STREQUAL "cmake")
      if(DEFINED "${_pkg_uc}_DIR")
         set("_${_pkg_uc}_DIR")
         set("${package}_DIR" "${_pkg_uc}_DIR")
      endif()
      find_package("${package}" CONFIG QUIET)
      if("${package}_FOUND")
         message(STATUS "Found ${package} via CMake config")
         break()
      endif()

   # pkgconf case
   elseif("${method}" STREQUAL "pkgconf")
      find_package("PkgConfig" QUIET) # built-in Find script
      pkg_check_modules("${_pkg_uc}" QUIET "${package}") # check if it is a pkg-config module
      if("${_pkg_uc}_FOUND")
         message(STATUS "Found ${package} via pkg-config")
         add_library("${package}::${package}" INTERFACE IMPORTED) # interface library
         target_link_libraries(
            "${package}::${package}"
            INTERFACE
            "${${_pkg_uc}_LINK_LIBRARIES}"
         )
         target_include_directories(
            "${package}::${package}"
            INTERFACE
            "${${_pkg_uc}_INCLUDE_DIRS}"
         )
         break()
      endif()

   # subproject case
   elseif("${method}" STREQUAL "subproject")
      if(NOT DEFINED "${_pkg_uc}_SUBPROJECT")
         set("${_pkg_uc}_SUBPROJECT" "subprojects/${package}")
      endif()
      set("${_pkg_uc}_SOURCE_DIR" "${PROJECT_SOURCE_DIR}/${${_pkg_uc}_SUBPROJECT}")
      set("${_pkg_uc}_BINARY_DIR" "${PROJECT_BINARY_DIR}/${${_pkg_uc}_SUBPROJECT}")

      # if can be configured from the subprojects dir
      if(EXISTS "${${_pkg_uc}_SOURCE_DIR}/CMakeLists.txt")
         message(STATUS "Include ${package} from ${${_pkg_uc}_SUBPROJECT}")

         add_subdirectory(
            "${${_pkg_uc}_SOURCE_DIR}"
            "${${_pkg_uc}_BINARY_DIR}"
         )

         # create interface directory and manage it's dependencies
         add_library("${package}::${package}" INTERFACE IMPORTED)
         target_link_libraries("${package}::${package}" INTERFACE "${package}")

         # We need the module directory in the subproject before we finish the configure stage
         if(NOT EXISTS "${${_pkg_uc}_BINARY_DIR}/include")
            file(MAKE_DIRECTORY "${${_pkg_uc}_BINARY_DIR}/include")
         endif()

         break()
      endif()

   # fetch from url case
   elseif("${method}" STREQUAL "fetch")
      message(STATUS "Retrieving ${package} from ${url}")
      include(FetchContent) # module for fetching from repo
      if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(FETCHCONTENT_QUIET FALSE)
      endif()
      FetchContent_Declare(
         "${_pkg_lc}"
         GIT_REPOSITORY "${url}"
         GIT_TAG "${rev}"
      )
      FetchContent_MakeAvailable("${_pkg_lc}")

      add_library("${package}::${package}" INTERFACE IMPORTED)
      target_link_libraries("${package}::${package}" INTERFACE "${package}")

      if(NOT EXISTS "${${_pkg_lc}_BINARY_DIR}/include")
         file(MAKE_DIRECTORY "${${_pkg_lc}_BINARY_DIR}/include")
      endif()

      break()
   endif()

endforeach()

unset(_pkg_lc)
unset(_pkg_uc)

# sanity check
if(NOT TARGET "${package}::${package}")
   message(FATAL_ERROR "Could not find dependency ${package}")
endif()
endmacro()
