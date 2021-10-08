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
  "hessian"
  "latticepoint"
  "molecule"
  "peeq"
  "repulsion"
  "symmetry"
  "thermo"
  "wsc"
)
set(
  test-srcs
  "main.f90"
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


# Build Tests
set(test-srcs2
  "${dir}/tests_peeq.f90"
  "${dir}/assertion.f90"
  "${dir}/molstock.f90"
  "${dir}/pbc_tools.f90"
  "${dir}/gfn2.f90"
  "${dir}/gfn1.f90"
  "${dir}/gfn0.f90"
  "${dir}/gfnff.f90"
)

add_executable(xtb_tests_static ${test-srcs2})
target_link_libraries(xtb_tests_static
  PUBLIC
  lib-xtb-static
)
set_target_properties(xtb_tests_static PROPERTIES
  Fortran_MODULE_DIRECTORY ${xtb-mod}
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)

set(XTB-EXE ${PROJECT_BINARY_DIR}/xtb)
set(XTB-TESTER ${PROJECT_BINARY_DIR}/xtb_tests_static)

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

add_test("PBC_tools_convert" ${XTB-TESTER} pbc_tools convert)
add_test("PBC_tools_cutoff" ${XTB-TESTER} pbc_tools cutoff)
add_test("GFN2-xTB_Basic" ${XTB-TESTER} gfn2 basic)
add_test("GFN2-xTB_Solvation" ${XTB-TESTER} gfn2 solvation)
add_test("GFN2-xTB_SCC" ${XTB-TESTER} gfn2 scc)
add_test("GFN2-xTB_API" ${XTB-TESTER} gfn2 api)
add_test("GFN2-xTB_API_GBSA" ${XTB-TESTER} gfn2 gbsa)
add_test("GFN2-xTB_API_GBSA_salt" ${XTB-TESTER} gfn2 salt)
add_test("GFN2-xTB_d_Metal" ${XTB-TESTER} gfn2 d-metal)
add_test("GFN2-xTB_COSMO" ${XTB-TESTER} gfn2 cosmo)
#add_test("GFN2-xTB_API_PCEM" ${XTB-TESTER} gfn2 pcem)
add_test("GFN1-xTB_Basic" ${XTB-TESTER} gfn1 basic)
add_test("GFN1-xTB_Solvation" ${XTB-TESTER} gfn1 solvation)
add_test("GFN1-xTB_SCC" ${XTB-TESTER} gfn1 scc)
add_test("GFN1-xTB_API" ${XTB-TESTER} gfn1 api)
add_test("GFN1-xTB_API_XB" ${XTB-TESTER} gfn1 xb)
add_test("GFN1-xTB_API_PBC" ${XTB-TESTER} gfn1 pbc3d)
add_test("GFN1-xTB_API_GBSA" ${XTB-TESTER} gfn1 gbsa)
add_test("GFN1-xTB_API_PCEM" ${XTB-TESTER} gfn1 pcem)
add_test("GFN1-xTB_IPEA" ${XTB-TESTER} gfn1 ipea)
add_test("GFN0-xTB_Basic" ${XTB-TESTER} gfn0 basic)
add_test("GFN0-xTB_Solvation" ${XTB-TESTER} gfn0 solvation)
add_test("GFN0-xTB_SP" ${XTB-TESTER} gfn0 sp)
add_test("GFN0-xTB_API" ${XTB-TESTER} gfn0 api)
add_test("GFN0-xTB_SRB" ${XTB-TESTER} gfn0 srb)

add_test("GFN-FF_Basic" ${XTB-TESTER} gfnff basic)
add_test("GFN-FF_Solvation" ${XTB-TESTER} gfnff solvation)
add_test("GFN-FF_SP" ${XTB-TESTER} gfnff sp)
add_test("GFN-FF_HB" ${XTB-TESTER} gfnff hb)
add_test("GFN-FF_GBSA" ${XTB-TESTER} gfnff gbsa)
add_test("GFN-FF_Scale_up" ${XTB-TESTER} gfnff scaleup)
add_test("GFN-FF_PDB" ${XTB-TESTER} gfnff pdb)
add_test("GFN-FF_SDF" ${XTB-TESTER} gfnff sdf)

set_tests_properties(
  "EXE_Singlepoint"
  "EXE_GFN0-xTB"
  "EXE_IP/EA"
  "EXE_GFN1-xTB"
  "EXE_GFN2-xTB/GBSA"
  "EXE_GFN2-FF"
  "GFN2-xTB_SCC"
  "GFN2-xTB_API"
  "GFN2-xTB_API_GBSA"
  "GFN2-xTB_API_GBSA_salt"
  #"GFN2-xTB_API_PCEM"
  "GFN1-xTB_SCC"
  "GFN1-xTB_API"
  "GFN1-xTB_API_XB"
  "GFN1-xTB_API_PBC"
  "GFN1-xTB_API_GBSA"
  "GFN1-xTB_API_PCEM"
  "GFN1-xTB_IPEA"
  "GFN0-xTB_Basic"
  "GFN0-xTB_Solvation"
  "GFN0-xTB_SP"
  "GFN0-xTB_API"
  "GFN0-xTB_SRB"
  PROPERTIES
  ENVIRONMENT XTBPATH=${PROJECT_SOURCE_DIR}
)
