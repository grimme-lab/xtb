set(TEST_EXE_LOCATION ../xtb)
set(TEST_STATIC_LOCATION ../xtb_tests_static)
set(TEST_SHARED_LOCATION ../xtb_tests_shared)
set(TEST_CAPI_LOCATION ../xtb_tests_capi)

add_test(EXE_Argparser_print_version ${TEST_EXE_LOCATION} --version)
add_test(EXE_Argparser_print_help ${TEST_EXE_LOCATION} --help)
add_test(EXE_Argparser_print_license ${TEST_EXE_LOCATION} --license)
add_test(EXE_Info ${TEST_EXE_LOCATION} info
   "${XTB_ROOT}/assets/inputs/coord/caffeine.coord"
   "${XTB_ROOT}/assets/inputs/coord/quartz.3d.coord"
   "${XTB_ROOT}/assets/inputs/vasp/ammonia.vasp"
   "${XTB_ROOT}/assets/inputs/xyz/taxol.xyz"
)
add_test(EXE_Singlepoint ${TEST_EXE_LOCATION} --coffee --strict --norestart)
add_test(EXE_IP/EA ${TEST_EXE_LOCATION} --coffee --gfn 2 --vipea --strict --norestart)
add_test(EXE_GFN0-xTB ${TEST_EXE_LOCATION} --coffee --gfn 0 --strict --norestart)
add_test(EXE_GFN1-xTB ${TEST_EXE_LOCATION} --coffee --gfn 1 --strict --norestart)
add_test(EXE_GFN2-xTB/GBSA ${TEST_EXE_LOCATION} --coffee --gfn 2 --strict --gbsa h2o --norestart)

add_test(STATIC_Molecule_axis ${TEST_STATIC_LOCATION} xtb_type_molecule axis)
add_test(STATIC_Molecule_MIC ${TEST_STATIC_LOCATION} xtb_type_molecule mic)
add_test(STATIC_LatticePoint ${TEST_STATIC_LOCATION} latticepoint pbc3d)
add_test(STATIC_coord_3D ${TEST_STATIC_LOCATION} geometry_reader coord_3d_a)
add_test(STATIC_coord_3D ${TEST_STATIC_LOCATION} geometry_reader coord_3d_b)
#add_test(STATIC_coord_2D ${TEST_STATIC_LOCATION} geometry_reader coord_2d)
#add_test(STATIC_coord_1D ${TEST_STATIC_LOCATION} geometry_reader coord_1d)
add_test(STATIC_coord_0D ${TEST_STATIC_LOCATION} geometry_reader coord_0d)
add_test(STATIC_Xmol__0D ${TEST_STATIC_LOCATION} geometry_reader xmol_0d)
add_test(STATIC_POSCAR ${TEST_STATIC_LOCATION} geometry_reader poscar_3d)
add_test(STATIC_molfile ${TEST_STATIC_LOCATION} geometry_reader molfile)
add_test(STATIC_SDF ${TEST_STATIC_LOCATION} geometry_reader sdfile)
add_test(STATIC_PDB ${TEST_STATIC_LOCATION} geometry_reader pdb)
add_test(STATIC_genFormat ${TEST_STATIC_LOCATION} geometry_reader gen)
add_test(STATIC_PBC_tools_convert ${TEST_STATIC_LOCATION} pbc_tools convert)
add_test(STATIC_PBC_tools_cutoff ${TEST_STATIC_LOCATION} pbc_tools cutoff)
add_test(STATIC_Symmetry_Water ${TEST_STATIC_LOCATION} symmetry water)
add_test(STATIC_Symmetry_Li8 ${TEST_STATIC_LOCATION} symmetry li8)
add_test(STATIC_Symmetry_PCl3 ${TEST_STATIC_LOCATION} symmetry pcl3)
add_test(STATIC_Symmetry_large ${TEST_STATIC_LOCATION} symmetry c20)
add_test(STATIC_Thermo_axis ${TEST_STATIC_LOCATION} thermo axis)
add_test(STATIC_Thermo_calculation ${TEST_STATIC_LOCATION} thermo calc)
add_test(STATIC_Thermo_print ${TEST_STATIC_LOCATION} thermo print)
add_test(STATIC_Repulsion_0d ${TEST_STATIC_LOCATION} repulsion cluster)
add_test(STATIC_Repulsion_3d ${TEST_STATIC_LOCATION} repulsion pbc3d)
add_test(STATIC_Coulomb_point_0d ${TEST_STATIC_LOCATION} coulomb point_0d)
add_test(STATIC_Coulomb_point_3d ${TEST_STATIC_LOCATION} coulomb point_3d)
add_test(STATIC_Coulomb_gfn1_0d ${TEST_STATIC_LOCATION} coulomb gfn1_0d)
add_test(STATIC_Coulomb_gfn1_3d ${TEST_STATIC_LOCATION} coulomb gfn1_3d)
add_test(STATIC_Coulomb_gfn2_0d ${TEST_STATIC_LOCATION} coulomb gfn2_0d)
add_test(STATIC_Coulomb_gfn2_3d ${TEST_STATIC_LOCATION} coulomb gfn2_3d)
add_test(STATIC_Coulomb_gauss_0d ${TEST_STATIC_LOCATION} coulomb gaussian_0d)
add_test(STATIC_Coulomb_gauss_3d ${TEST_STATIC_LOCATION} coulomb gaussian_3d)
add_test(STATIC_EEQ_model_water ${TEST_STATIC_LOCATION} eeq_model water)
add_test(STATIC_EEQ_model_3D_Ewald ${TEST_STATIC_LOCATION} eeq_model ewald)
add_test(STATIC_EEQ_model_GBSA ${TEST_STATIC_LOCATION} eeq_model gbsa)
add_test(STATIC_EEQ_model_GBSA_salt ${TEST_STATIC_LOCATION} eeq_model salt)
#add_test(STATIC_EEQ_model_GBSA_H-bond ${TEST_STATIC_LOCATION} eeq_model hbond)
add_test(STATIC_CN_latticepoints ${TEST_STATIC_LOCATION} ncoord pbc3dlatp)
add_test(STATIC_CN_neighbourlist ${TEST_STATIC_LOCATION} ncoord pbc3dneighs)
add_test(STATIC_DFTD3_latticepoints ${TEST_STATIC_LOCATION} dftd3 pbc3dlatp)
add_test(STATIC_DFTD3_neighbourlist ${TEST_STATIC_LOCATION} dftd3 pbc3dneighs)
add_test(STATIC_DFTD3_threebody_lp ${TEST_STATIC_LOCATION} dftd3 pbc3datmlatp)
add_test(STATIC_DFTD3_threebody_nl ${TEST_STATIC_LOCATION} dftd3 pbc3datmneighs)
add_test(STATIC_DFTD4_latticepoints ${TEST_STATIC_LOCATION} dftd4 pbc3dlatp)
add_test(STATIC_DFTD4_neighbourlist ${TEST_STATIC_LOCATION} dftd4 pbc3dneighs)
add_test(STATIC_DFTD4_threebody_lp ${TEST_STATIC_LOCATION} dftd4 pbc3datmlatp)
add_test(STATIC_DFTD4_threebody_nl ${TEST_STATIC_LOCATION} dftd4 pbc3datmneighs)
add_test(STATIC_Dispersion_energies ${TEST_STATIC_LOCATION} dftd4 energies)
add_test(STATIC_Dispersion_energies_PBC ${TEST_STATIC_LOCATION} dftd4 pbc_disp)
add_test(STATIC_Dispersion_API ${TEST_STATIC_LOCATION} dftd4 api)
add_test(STATIC_GFN2-xTB_SCC ${TEST_STATIC_LOCATION} gfn2 scc)
add_test(STATIC_GFN2-xTB_API ${TEST_STATIC_LOCATION} gfn2 api)
add_test(STATIC_GFN2-xTB_API_GBSA ${TEST_STATIC_LOCATION} gfn2 gbsa)
add_test(STATIC_GFN2-xTB_API_GBSA_salt ${TEST_STATIC_LOCATION} gfn2 salt)
add_test(STATIC_GFN2-xTB_API_PCEM ${TEST_STATIC_LOCATION} gfn2 pcem)
add_test(STATIC_GFN1-xTB_SCC ${TEST_STATIC_LOCATION} gfn1 scc)
add_test(STATIC_GFN1-xTB_API ${TEST_STATIC_LOCATION} gfn1 api)
add_test(STATIC_GFN1-xTB_API_XB ${TEST_STATIC_LOCATION} gfn1 xb)
add_test(STATIC_GFN1-xTB_API_PBC ${TEST_STATIC_LOCATION} gfn1 pbc3d)
add_test(STATIC_GFN1-xTB_API_GBSA ${TEST_STATIC_LOCATION} gfn1 gbsa)
add_test(STATIC_GFN1-xTB_API_PCEM ${TEST_STATIC_LOCATION} gfn1 pcem)
add_test(STATIC_GFN0-xTB_SP ${TEST_STATIC_LOCATION} gfn0 sp)
add_test(STATIC_GFN0-xTB_API ${TEST_STATIC_LOCATION} gfn0 api)
add_test(STATIC_GFN0-xTB_SRB ${TEST_STATIC_LOCATION} gfn0 srb)
add_test(STATIC_GFN0-xTB_SP_PBC ${TEST_STATIC_LOCATION} peeq sp)
add_test(STATIC_GFN0-xTB_API_PBC ${TEST_STATIC_LOCATION} peeq api)
add_test(STATIC_GFN0-xTB_SRB_PBC ${TEST_STATIC_LOCATION} peeq srb)

add_test(SHARED_Molecule_axis ${TEST_SHARED_LOCATION} xtb_type_molecule axis)
add_test(SHARED_Molecule_MIC ${TEST_SHARED_LOCATION} xtb_type_molecule mic)
add_test(SHARED_LatticePoint ${TEST_SHARED_LOCATION} latticepoint pbc3d)
add_test(SHARED_coord_3D ${TEST_SHARED_LOCATION} geometry_reader coord_3d_a)
add_test(SHARED_coord_3D ${TEST_SHARED_LOCATION} geometry_reader coord_3d_b)
#add_test(SHARED_coord_2D ${TEST_SHARED_LOCATION} geometry_reader coord_2d)
#add_test(SHARED_coord_1D ${TEST_SHARED_LOCATION} geometry_reader coord_1d)
add_test(SHARED_coord_0D ${TEST_SHARED_LOCATION} geometry_reader coord_0d)
add_test(SHARED_Xmol__0D ${TEST_SHARED_LOCATION} geometry_reader xmol_0d)
add_test(SHARED_POSCAR ${TEST_SHARED_LOCATION} geometry_reader poscar_3d)
add_test(SHARED_molfile ${TEST_SHARED_LOCATION} geometry_reader molfile)
add_test(SHARED_SDF ${TEST_SHARED_LOCATION} geometry_reader sdfile)
add_test(SHARED_PDB ${TEST_SHARED_LOCATION} geometry_reader pdb)
add_test(SHARED_genFormat ${TEST_SHARED_LOCATION} geometry_reader gen)
add_test(SHARED_PBC_tools_convert ${TEST_SHARED_LOCATION} pbc_tools convert)
add_test(SHARED_PBC_tools_cutoff ${TEST_SHARED_LOCATION} pbc_tools cutoff)
add_test(SHARED_Symmetry_Water ${TEST_SHARED_LOCATION} symmetry water)
add_test(SHARED_Symmetry_Li8 ${TEST_SHARED_LOCATION} symmetry li8)
add_test(SHARED_Symmetry_PCl3 ${TEST_SHARED_LOCATION} symmetry pcl3)
add_test(SHARED_Symmetry_large ${TEST_SHARED_LOCATION} symmetry c20)
add_test(SHARED_Thermo_axis ${TEST_SHARED_LOCATION} thermo axis)
add_test(SHARED_Thermo_calculation ${TEST_SHARED_LOCATION} thermo calc)
add_test(SHARED_Thermo_print ${TEST_SHARED_LOCATION} thermo print)
add_test(SHARED_repulsion_0d ${TEST_SHARED_LOCATION} repulsion cluster)
add_test(SHARED_repulsion_3d ${TEST_SHARED_LOCATION} repulsion pbc3d)
add_test(SHARED_Coulomb_point_0d ${TEST_SHARED_LOCATION} coulomb point_0d)
add_test(SHARED_Coulomb_point_3d ${TEST_SHARED_LOCATION} coulomb point_3d)
add_test(SHARED_Coulomb_gfn1_0d ${TEST_SHARED_LOCATION} coulomb gfn1_0d)
add_test(SHARED_Coulomb_gfn1_3d ${TEST_SHARED_LOCATION} coulomb gfn1_3d)
add_test(SHARED_Coulomb_gfn2_0d ${TEST_SHARED_LOCATION} coulomb gfn2_0d)
add_test(SHARED_Coulomb_gfn2_3d ${TEST_SHARED_LOCATION} coulomb gfn2_3d)
add_test(SHARED_Coulomb_gauss_0d ${TEST_SHARED_LOCATION} coulomb gaussian_0d)
add_test(SHARED_Coulomb_gauss_3d ${TEST_SHARED_LOCATION} coulomb gaussian_3d)
add_test(SHARED_EEQ_model_water ${TEST_SHARED_LOCATION} eeq_model water)
add_test(SHARED_EEQ_model_3D_Ewald ${TEST_SHARED_LOCATION} eeq_model ewald)
add_test(SHARED_EEQ_model_GBSA ${TEST_SHARED_LOCATION} eeq_model gbsa)
add_test(SHARED_EEQ_model_GBSA_salt ${TEST_SHARED_LOCATION} eeq_model salt)
#add_test(SHARED_EEQ_model_GBSA_H-bond ${TEST_SHARED_LOCATION} eeq_model hbond)
add_test(SHARED_CN_latticepoints ${TEST_SHARED_LOCATION} ncoord pbc3dlatp)
add_test(SHARED_CN_neighbourlist ${TEST_SHARED_LOCATION} ncoord pbc3dneighs)
add_test(SHARED_DFTD3_latticepoints ${TEST_SHARED_LOCATION} dftd3 pbc3dlatp)
add_test(SHARED_DFTD3_neighbourlist ${TEST_SHARED_LOCATION} dftd3 pbc3dneighs)
add_test(SHARED_DFTD3_threebody_lp ${TEST_SHARED_LOCATION} dftd3 pbc3datmlatp)
add_test(SHARED_DFTD3_threebody_nl ${TEST_SHARED_LOCATION} dftd3 pbc3datmneighs)
add_test(SHARED_DFTD4_latticepoints ${TEST_SHARED_LOCATION} dftd4 pbc3dlatp)
add_test(SHARED_DFTD4_neighbourlist ${TEST_SHARED_LOCATION} dftd4 pbc3dneighs)
add_test(SHARED_DFTD4_threebody_lp ${TEST_SHARED_LOCATION} dftd4 pbc3datmlatp)
add_test(SHARED_DFTD4_threebody_nl ${TEST_SHARED_LOCATION} dftd4 pbc3datmneighs)
add_test(SHARED_Dispersion_energies ${TEST_SHARED_LOCATION} dftd4 energies)
add_test(SHARED_Dispersion_energies_PBC ${TEST_SHARED_LOCATION} dftd4 pbc_disp)
add_test(SHARED_Dispersion_API ${TEST_SHARED_LOCATION} dftd4 api)
add_test(SHARED_GFN2-xTB_SCC ${TEST_SHARED_LOCATION} gfn2 scc)
add_test(SHARED_GFN2-xTB_API ${TEST_SHARED_LOCATION} gfn2 api)
add_test(SHARED_GFN2-xTB_API_GBSA ${TEST_SHARED_LOCATION} gfn2 gbsa)
add_test(SHARED_GFN2-xTB_API_GBSA_salt ${TEST_SHARED_LOCATION} gfn2 salt)
add_test(SHARED_GFN2-xTB_API_PCEM ${TEST_SHARED_LOCATION} gfn2 pcem)
add_test(SHARED_GFN1-xTB_SCC ${TEST_SHARED_LOCATION} gfn1 scc)
add_test(SHARED_GFN1-xTB_API ${TEST_SHARED_LOCATION} gfn1 api)
add_test(SHARED_GFN1-xTB_API_XB ${TEST_SHARED_LOCATION} gfn1 xb)
add_test(SHARED_GFN1-xTB_API_PBC ${TEST_SHARED_LOCATION} gfn1 pbc3d)
add_test(SHARED_GFN1-xTB_API_GBSA ${TEST_SHARED_LOCATION} gfn1 gbsa)
add_test(SHARED_GFN1-xTB_API_PCEM ${TEST_SHARED_LOCATION} gfn1 pcem)
add_test(SHARED_GFN0-xTB_SP ${TEST_SHARED_LOCATION} gfn0 sp)
add_test(SHARED_GFN0-xTB_API ${TEST_SHARED_LOCATION} gfn0 api)
add_test(SHARED_GFN0-xTB_SRB ${TEST_SHARED_LOCATION} gfn0 srb)
add_test(SHARED_GFN0-xTB_SP_PBC ${TEST_SHARED_LOCATION} peeq sp)
add_test(SHARED_GFN0-xTB_API_PBC ${TEST_SHARED_LOCATION} peeq api)
add_test(SHARED_GFN0-xTB_SRB_PBC ${TEST_SHARED_LOCATION} peeq srb)

add_test(CAPI ${TEST_CAPI_LOCATION})

set_tests_properties(EXE_Singlepoint
   EXE_GFN0-xTB
   EXE_IP/EA
   EXE_GFN1-xTB
   EXE_GFN2-xTB/GBSA
   CAPI
   STATIC_GFN2-xTB_SCC
   STATIC_GFN2-xTB_API
   STATIC_GFN2-xTB_API_GBSA
   STATIC_GFN2-xTB_API_GBSA_salt
   STATIC_GFN2-xTB_API_PCEM
   STATIC_GFN1-xTB_SCC
   STATIC_GFN1-xTB_API
   STATIC_GFN1-xTB_API_XB
   STATIC_GFN1-xTB_API_PBC
   STATIC_GFN1-xTB_API_GBSA
   STATIC_GFN1-xTB_API_PCEM
   STATIC_GFN0-xTB_SP
   STATIC_GFN0-xTB_API
   STATIC_GFN0-xTB_SRB
   STATIC_GFN0-xTB_SP_PBC
   STATIC_GFN0-xTB_API_PBC
   STATIC_GFN0-xTB_SRB_PBC
   SHARED_GFN2-xTB_SCC
   SHARED_GFN2-xTB_API
   SHARED_GFN2-xTB_API_GBSA
   SHARED_GFN2-xTB_API_GBSA_salt
   SHARED_GFN2-xTB_API_PCEM
   SHARED_GFN1-xTB_SCC
   SHARED_GFN1-xTB_API
   SHARED_GFN1-xTB_API_XB
   SHARED_GFN1-xTB_API_PBC
   SHARED_GFN1-xTB_API_GBSA
   SHARED_GFN1-xTB_API_PCEM
   SHARED_GFN0-xTB_SP
   SHARED_GFN0-xTB_API
   SHARED_GFN0-xTB_SRB
   SHARED_GFN0-xTB_SP_PBC
   SHARED_GFN0-xTB_API_PBC
   SHARED_GFN0-xTB_SRB_PBC
   PROPERTIES
   ENVIRONMENT XTBPATH=${CMAKE_CURRENT_SOURCE_DIR}/..)
