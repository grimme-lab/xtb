# This file is part of xtb.
#
# Copyright (C) 2019-2021 Sebastian Ehlert
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

if(NOT TARGET "test-drive::test-drive")
  find_package("test-drive" REQUIRED)
endif()

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

set(
  tests
  "atomlist"
  "coordinationnumber"
  "coulomb"
  "dftd3"
  "dftd4"
  "eeq"
  "geometry-reader"
  "gfn0"
  "gfn1"
  "gfn2"
  "gfnff"
  "hessian"
  "latticepoint"
  "molecule"
  "pbc-tools"
  "peeq"
  "repulsion"
  "symmetry"
  "thermo"
  "wsc"
)
set(
  test-srcs
  "main.f90"
  "molstock.f90"
)
foreach(t IN LISTS tests)
  string(MAKE_C_IDENTIFIER ${t} t) 
  list(APPEND test-srcs "test_${t}.f90")
endforeach()

add_executable(
  "${PROJECT_NAME}-tester"
  "${test-srcs}"
)
target_link_libraries(
  "${PROJECT_NAME}-tester"
  PRIVATE
  lib-xtb-static
  "test-drive::test-drive"
)
set_target_properties(
  "${PROJECT_NAME}-tester"
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${xtb-mod}
)

foreach(t IN LISTS tests)
  add_test("${PROJECT_NAME}/${t}" "${PROJECT_NAME}-tester" "${t}")
  set_tests_properties(
    "${PROJECT_NAME}/${t}"
    PROPERTIES
    ENVIRONMENT XTBPATH=${PROJECT_SOURCE_DIR}
  )
endforeach()


set(XTB-EXE ${PROJECT_BINARY_DIR}/xtb)

add_test("EXE_Argparser_print_version" ${XTB-EXE} --version)
add_test("EXE_Argparser_print_help" ${XTB-EXE} --help)
add_test("EXE_Argparser_print_license" ${XTB-EXE} --license)
add_test("EXE_Info" ${XTB-EXE} info
  "${xtb-dir}/assets/inputs/coord/caffeine.coord"
  "${xtb-dir}/assets/inputs/coord/quartz.3d.coord"
  "${xtb-dir}/assets/inputs/vasp/ammonia.vasp"
  "${xtb-dir}/assets/inputs/xyz/taxol.xyz"
)
add_test("EXE_Singlepoint" ${XTB-EXE} --coffee --strict --norestart --namespace test1)
add_test("EXE_IP/EA" ${XTB-EXE} --coffee --gfn 2 --vipea --strict --norestart --namespace test2)
add_test("EXE_GFN0-xTB" ${XTB-EXE} --coffee --gfn 0 --strict --norestart --namespace test3)
add_test("EXE_GFN1-xTB" ${XTB-EXE} --coffee --gfn 1 --strict --norestart --namespace test4)
add_test("EXE_GFN2-xTB/GBSA" ${XTB-EXE} --coffee --gfn 2 --strict --gbsa h2o --norestart --namespace test5)
add_test("EXE_GFN2-FF" ${XTB-EXE} --coffee --gfnff --strict --norestart --namespace test6)

set_tests_properties(
  "EXE_Singlepoint"
  "EXE_GFN0-xTB"
  "EXE_IP/EA"
  "EXE_GFN1-xTB"
  "EXE_GFN2-xTB/GBSA"
  "EXE_GFN2-FF"
  PROPERTIES
  ENVIRONMENT XTBPATH=${PROJECT_SOURCE_DIR}
)
