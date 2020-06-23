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

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

set(tests
  "${dir}/tests_peeq.f90"
  "${dir}/assertion.f90"
  "${dir}/molstock.f90"
  "${dir}/tbdef_molecule.f90"
  "${dir}/latticepoint.f90"
  "${dir}/tbdef_wsc.f90"
  "${dir}/geometry_reader.f90"
  "${dir}/coordinationnumber.f90"
  "${dir}/coulomb.f90"
  "${dir}/eeq.f90"
  "${dir}/eeq_model.f90"
  "${dir}/pbc_tools.f90"
  "${dir}/dftd3.f90"
  "${dir}/dftd4.f90"
  "${dir}/gfn2.f90"
  "${dir}/gfn1.f90"
  "${dir}/gfn0.f90"
  "${dir}/gfnff.f90"
  "${dir}/hessian.f90"
  "${dir}/peeq.f90"
  "${dir}/repulsion.f90"
  "${dir}/symmetry.f90"
  "${dir}/thermo.f90"
  "${dir}/tbdef_atomlist.f90"
)

# Build Tests
add_executable(xtb_tests_static ${tests})
target_link_libraries(xtb_tests_static
  PUBLIC
  lib-xtb-static
)
set_target_properties(xtb_tests_static PROPERTIES
  Fortran_MODULE_DIRECTORY ${xtb-mod}
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_executable(xtb_tests_capi "${xtb-dir}/TESTSUITE/c_api_example.c")
target_link_libraries(xtb_tests_capi
  PUBLIC
  lib-xtb-shared
)
set_target_properties(xtb_tests_capi PROPERTIES
  Fortran_MODULE_DIRECTORY ${xtb-mod}
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

set(XTB-EXE ${CMAKE_BINARY_DIR}/xtb)
set(XTB-TESTER ${CMAKE_BINARY_DIR}/xtb_tests_static)
set(XTB-CTESTER ${CMAKE_BINARY_DIR}/xtb_tests_capi)

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

add_test("Molecule_axis" ${XTB-TESTER} xtb_type_molecule axis)
add_test("Molecule_MIC" ${XTB-TESTER} xtb_type_molecule mic)
add_test("LatticePoint" ${XTB-TESTER} latticepoint pbc3d)
add_test("coord_3D" ${XTB-TESTER} geometry_reader coord_3d_a)
add_test("coord_3D" ${XTB-TESTER} geometry_reader coord_3d_b)
add_test("coord_0D" ${XTB-TESTER} geometry_reader coord_0d)
add_test("Xmol__0D" ${XTB-TESTER} geometry_reader xmol_0d)
add_test("POSCAR" ${XTB-TESTER} geometry_reader poscar_3d)
add_test("molfile" ${XTB-TESTER} geometry_reader molfile)
add_test("SDF" ${XTB-TESTER} geometry_reader sdfile)
add_test("PDB" ${XTB-TESTER} geometry_reader pdb)
add_test("genFormat" ${XTB-TESTER} geometry_reader gen)
add_test("PBC_tools_convert" ${XTB-TESTER} pbc_tools convert)
add_test("PBC_tools_cutoff" ${XTB-TESTER} pbc_tools cutoff)
add_test("Symmetry_Water" ${XTB-TESTER} symmetry water)
add_test("Symmetry_Li8" ${XTB-TESTER} symmetry li8)
add_test("Symmetry_PCl3" ${XTB-TESTER} symmetry pcl3)
add_test("Symmetry_large" ${XTB-TESTER} symmetry c20)
add_test("Thermo_axis" ${XTB-TESTER} thermo axis)
add_test("Thermo_calculation" ${XTB-TESTER} thermo calc)
add_test("Thermo_print" ${XTB-TESTER} thermo print)
add_test("Repulsion_0d" ${XTB-TESTER} repulsion cluster)
add_test("Repulsion_3d" ${XTB-TESTER} repulsion pbc3d)
add_test("Coulomb_point_0d" ${XTB-TESTER} coulomb point_0d)
add_test("Coulomb_point_3d" ${XTB-TESTER} coulomb point_3d)
add_test("Coulomb_gfn1_0d" ${XTB-TESTER} coulomb gfn1_0d)
add_test("Coulomb_gfn1_3d" ${XTB-TESTER} coulomb gfn1_3d)
add_test("Coulomb_gfn2_0d" ${XTB-TESTER} coulomb gfn2_0d)
add_test("Coulomb_gfn2_3d" ${XTB-TESTER} coulomb gfn2_3d)
add_test("Coulomb_gauss_0d" ${XTB-TESTER} coulomb gaussian_0d)
add_test("Coulomb_gauss_3d" ${XTB-TESTER} coulomb gaussian_3d)
add_test("EEQ_model_water" ${XTB-TESTER} eeq_model water)
add_test("EEQ_model_3D_Ewald" ${XTB-TESTER} eeq_model ewald)
add_test("EEQ_model_GBSA" ${XTB-TESTER} eeq_model gbsa)
add_test("EEQ_model_GBSA_salt" ${XTB-TESTER} eeq_model salt)
add_test("CN_latticepoints" ${XTB-TESTER} ncoord pbc3dlatp)
add_test("CN_neighbourlist" ${XTB-TESTER} ncoord pbc3dneighs)
add_test("DFTD3_latticepoints" ${XTB-TESTER} dftd3 pbc3dlatp)
add_test("DFTD3_neighbourlist" ${XTB-TESTER} dftd3 pbc3dneighs)
add_test("DFTD3_threebody_lp" ${XTB-TESTER} dftd3 pbc3datmlatp)
add_test("DFTD3_threebody_nl" ${XTB-TESTER} dftd3 pbc3datmneighs)
add_test("DFTD4_latticepoints" ${XTB-TESTER} dftd4 pbc3dlatp)
add_test("DFTD4_neighbourlist" ${XTB-TESTER} dftd4 pbc3dneighs)
add_test("DFTD4_threebody_lp" ${XTB-TESTER} dftd4 pbc3datmlatp)
add_test("DFTD4_threebody_nl" ${XTB-TESTER} dftd4 pbc3datmneighs)
add_test("GFN2-xTB_Basic" ${XTB-TESTER} gfn2 basic)
add_test("GFN2-xTB_Solvation" ${XTB-TESTER} gfn2 solvation)
add_test("GFN2-xTB_SCC" ${XTB-TESTER} gfn2 scc)
add_test("GFN2-xTB_API" ${XTB-TESTER} gfn2 api)
add_test("GFN2-xTB_API_GBSA" ${XTB-TESTER} gfn2 gbsa)
add_test("GFN2-xTB_API_GBSA_salt" ${XTB-TESTER} gfn2 salt)
add_test("GFN2-xTB_d_Metal" ${XTB-TESTER} gfn2 d-metal)
add_test("GFN2-xTB_freq" ${XTB-TESTER} hessian gfn2)
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
add_test("GFN1-xTB_freq" ${XTB-TESTER} hessian gfn1)
add_test("GFN1-xTB_IPEA" ${XTB-TESTER} gfn1 ipea)
add_test("GFN0-xTB_Basic" ${XTB-TESTER} gfn0 basic)
add_test("GFN0-xTB_Solvation" ${XTB-TESTER} gfn0 solvation)
add_test("GFN0-xTB_SP" ${XTB-TESTER} gfn0 sp)
add_test("GFN0-xTB_API" ${XTB-TESTER} gfn0 api)
add_test("GFN0-xTB_SRB" ${XTB-TESTER} gfn0 srb)
add_test("GFN0-xTB_SP_PBC" ${XTB-TESTER} peeq sp)
add_test("GFN0-xTB_API_PBC" ${XTB-TESTER} peeq api)
add_test("GFN0-xTB_SRB_PBC" ${XTB-TESTER} peeq srb)

add_test("GFN-FF_Basic" ${XTB-TESTER} gfnff basic)
add_test("GFN-FF_Solvation" ${XTB-TESTER} gfnff solvation)
add_test("GFN-FF_SP" ${XTB-TESTER} gfnff sp)
add_test("GFN-FF_HB" ${XTB-TESTER} gfnff hb)
add_test("GFN-FF_GBSA" ${XTB-TESTER} gfnff gbsa)
add_test("GFN-FF_Scale_up" ${XTB-TESTER} gfnff scaleup)
add_test("GFN-FF_PDB" ${XTB-TESTER} gfnff pdb)
add_test("GFN-FF_SDF" ${XTB-TESTER} gfnff sdf)

add_test("CAPI" ${XTB-CTESTER})

set_tests_properties(EXE_Singlepoint
  "EXE_GFN0-xTB"
  "EXE_IP/EA"
  "EXE_GFN1-xTB"
  "EXE_GFN2-xTB/GBSA"
  "EXE_GFN2-FF"
  "CAPI"
  "GFN2-xTB_SCC"
  "GFN2-xTB_API"
  "GFN2-xTB_API_GBSA"
  "GFN2-xTB_API_GBSA_salt"
  "GFN2-xTB_freq"
  #"GFN2-xTB_API_PCEM"
  "GFN1-xTB_SCC"
  "GFN1-xTB_API"
  "GFN1-xTB_API_XB"
  "GFN1-xTB_API_PBC"
  "GFN1-xTB_API_GBSA"
  "GFN1-xTB_API_PCEM"
  "GFN1-xTB_IPEA"
  "GFN1-xTB_freq"
  "GFN0-xTB_Basic"
  "GFN0-xTB_Solvation"
  "GFN0-xTB_SP"
  "GFN0-xTB_API"
  "GFN0-xTB_SRB"
  "GFN0-xTB_SP_PBC"
  "GFN0-xTB_API_PBC"
  "GFN0-xTB_SRB_PBC"
  PROPERTIES
  ENVIRONMENT XTBPATH=${CMAKE_CURRENT_SOURCE_DIR}/..
)
