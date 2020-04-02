cmake_minimum_required(VERSION 3.9)
set(XTB_ROOT ${PROJECT_SOURCE_DIR})
set(XTB_SOURCES
  "${CMAKE_CURRENT_BINARY_DIR}/xtb_version.fh"

  "${XTB_ROOT}/src/prog/argparser.f90"
  "${XTB_ROOT}/src/prog/info.f90"
  "${XTB_ROOT}/src/prog/main.F90"
  "${XTB_ROOT}/src/prog/submodules.f90"

  # MCTC library
  "${XTB_ROOT}/src/mctc/accuracy.f90"
  "${XTB_ROOT}/src/mctc/boundaryconditions.f90"
  "${XTB_ROOT}/src/mctc/chartools.f90"
  "${XTB_ROOT}/src/mctc/convert.f90"
  "${XTB_ROOT}/src/mctc/filetypes.f90"
  "${XTB_ROOT}/src/mctc/io.f90"
  "${XTB_ROOT}/src/mctc/math.f90"
  "${XTB_ROOT}/src/mctc/namegen.f90"
  "${XTB_ROOT}/src/mctc/resize.f90"
  "${XTB_ROOT}/src/mctc/runtypes.f90"
  "${XTB_ROOT}/src/mctc/search.f90"
  "${XTB_ROOT}/src/mctc/sort.f90"
  "${XTB_ROOT}/src/mctc/symbols.f90"
  "${XTB_ROOT}/src/mctc/thresholds.f90"
  "${XTB_ROOT}/src/mctc/version.f90"
  "${XTB_ROOT}/src/mctc/mctc_global.f90"
  "${XTB_ROOT}/src/mctc/mctc_systools.f90"
  "${XTB_ROOT}/src/mctc/mctc_strings.f90"
  "${XTB_ROOT}/src/mctc/mctc_constants.f90"
  "${XTB_ROOT}/src/mctc/mctc_param.f90"
  "${XTB_ROOT}/src/mctc/param/atomic_masses.f90"
  "${XTB_ROOT}/src/mctc/param/chemical_hardnesses.f90"
  "${XTB_ROOT}/src/mctc/param/covalent_radii.f90"
  "${XTB_ROOT}/src/mctc/param/electronegativities.f90"
  "${XTB_ROOT}/src/mctc/param/pse.f90"
  "${XTB_ROOT}/src/mctc/param/r4r2_expectation_values.f90"
  "${XTB_ROOT}/src/mctc/mctc_timings.f90"
  "${XTB_ROOT}/src/mctc/mctc_filetools.f90"
  "${XTB_ROOT}/src/mctc/mctc_la.f90"
  "${XTB_ROOT}/src/mctc/mctc_init.f90"
  "${XTB_ROOT}/src/mctc/error.f90"
  "${XTB_ROOT}/src/mctc/signal.c"

  # Class definitions
  "${XTB_ROOT}/src/type/setvar.f90"
  "${XTB_ROOT}/src/type/wavefunction.f90"
  "${XTB_ROOT}/src/type/basisset.f90"
  "${XTB_ROOT}/src/type/molecule.f90"
  "${XTB_ROOT}/src/type/solvent.f90"
  "${XTB_ROOT}/src/type/param.f90"
  "${XTB_ROOT}/src/type/data.f90"
  "${XTB_ROOT}/src/type/anc.f90"
  "${XTB_ROOT}/src/type/timer.f90"
  "${XTB_ROOT}/src/type/pcem.f90"
  "${XTB_ROOT}/src/type/wsc.f90"
  "${XTB_ROOT}/src/type/options.f90"
  "${XTB_ROOT}/src/type/calculator.f90"
  "${XTB_ROOT}/src/type/atomlist.f90"
  "${XTB_ROOT}/src/type/topology.f90"
  "${XTB_ROOT}/src/type/fragments.f90"
  "${XTB_ROOT}/src/type/dispersion_model.f90"
  "${XTB_ROOT}/src/type/environment.f90"
  "${XTB_ROOT}/src/type/vendordata.f90"
  "${XTB_ROOT}/src/type/reader.f90"
  "${XTB_ROOT}/src/type/identitymap.f90"
  "${XTB_ROOT}/src/type/buffer.f90"
  "${XTB_ROOT}/src/type/neighbourlist.f90"
  "${XTB_ROOT}/src/type/wignerseitzcell.f90"
  "${XTB_ROOT}/src/type/latticepoint.f90"
  "${XTB_ROOT}/src/type/coulomb.f90"

  # Global data
  "${XTB_ROOT}/src/gfn0param.f90"
  "${XTB_ROOT}/src/sphereparam.f90"
  "${XTB_ROOT}/src/scanparam.f90"
  "${XTB_ROOT}/src/splitparam.f90"
  "${XTB_ROOT}/src/symparam.f90"
  "${XTB_ROOT}/src/fixparam.f90"
  "${XTB_ROOT}/src/aoparam.f90"
  "${XTB_ROOT}/src/setparam.f90"

  # Header and I/O
  "${XTB_ROOT}/src/readin.f90"
  "${XTB_ROOT}/src/filetools.f90"
  "${XTB_ROOT}/src/set_module.f90"
  "${XTB_ROOT}/src/constrain_param.f90"
  "${XTB_ROOT}/src/header.f90"
  "${XTB_ROOT}/src/printout.f90"
  "${XTB_ROOT}/src/xhelp.f90"
  "${XTB_ROOT}/src/expire.f90"
  "${XTB_ROOT}/src/define.f90"
  "${XTB_ROOT}/src/readl.f90"
  "${XTB_ROOT}/src/readl2.f90"
  "${XTB_ROOT}/src/rdcoord2.f90"
  "${XTB_ROOT}/src/wrbas.f90"
  "${XTB_ROOT}/src/wrmo.f90"
  "${XTB_ROOT}/src/main/property.f90"
  "${XTB_ROOT}/src/main/json.f90"
  "${XTB_ROOT}/src/main/geometry.f90"
  "${XTB_ROOT}/src/main/setup.f90"
  "${XTB_ROOT}/src/gfn_prparam.f90"
  "${XTB_ROOT}/src/wrgbw.f90"
  "${XTB_ROOT}/src/enso_printout.f90"
  "${XTB_ROOT}/src/read_gfn_param.f90"

  "${XTB_ROOT}/src/io/reader.f90"
  "${XTB_ROOT}/src/io/reader/ctfile.f90"
  "${XTB_ROOT}/src/io/reader/gaussian.f90"
  "${XTB_ROOT}/src/io/reader/genformat.f90"
  "${XTB_ROOT}/src/io/reader/pdb.f90"
  "${XTB_ROOT}/src/io/reader/turbomole.f90"
  "${XTB_ROOT}/src/io/reader/vasp.f90"
  "${XTB_ROOT}/src/io/reader/xyz.f90"
  "${XTB_ROOT}/src/io/writer.f90"
  "${XTB_ROOT}/src/io/writer/ctfile.f90"
  "${XTB_ROOT}/src/io/writer/gaussian.f90"
  "${XTB_ROOT}/src/io/writer/genformat.f90"
  "${XTB_ROOT}/src/io/writer/orca.f90"
  "${XTB_ROOT}/src/io/writer/pdb.f90"
  "${XTB_ROOT}/src/io/writer/turbomole.f90"
  "${XTB_ROOT}/src/io/writer/vasp.f90"
  "${XTB_ROOT}/src/io/writer/xyz.f90"

  # Parameters
  "${XTB_ROOT}/src/gfn_paramset.f90"
  "${XTB_ROOT}/src/disp/dftd3_parameters.f90"
  "${XTB_ROOT}/src/disp/dftd4_parameters.f90"
  "${XTB_ROOT}/src/charge_model.f90"
  "${XTB_ROOT}/src/param/atomicmass.f90"
  "${XTB_ROOT}/src/param/atomicrad.f90"
  "${XTB_ROOT}/src/param/covalentradd3.f90"
  "${XTB_ROOT}/src/param/paulingen.f90"
  "${XTB_ROOT}/src/param/sqrtzr4r2.f90"
  "${XTB_ROOT}/src/param/vdwradd3.f90"

  # Initial guess
  "${XTB_ROOT}/src/iniq.f90"
  "${XTB_ROOT}/src/eeq_model.f90"

  # Basis set
  "${XTB_ROOT}/src/xtb/basisset.f90"
  "${XTB_ROOT}/src/xtb/calculator.f90"
  "${XTB_ROOT}/src/xtb/data.f90"
  "${XTB_ROOT}/src/xtb/gfn0.f90"
  "${XTB_ROOT}/src/xtb/gfn1.f90"
  "${XTB_ROOT}/src/xtb/gfn2.f90"
  "${XTB_ROOT}/src/xtb/halogen.f90"
  "${XTB_ROOT}/src/xtb/repulsion.f90"
  "${XTB_ROOT}/src/slater.f90"

  # Dispersion
  "${XTB_ROOT}/src/disp/ncoord.f90"
  "${XTB_ROOT}/src/disp/coordinationnumber.f90"
  "${XTB_ROOT}/src/disp/dftd3.f"
  "${XTB_ROOT}/src/disp/dftd3.f90"
  "${XTB_ROOT}/src/disp/dftd4.f90"

  # Coulomb
  "${XTB_ROOT}/src/coulomb/ewald.f90"
  "${XTB_ROOT}/src/coulomb/gaussian.f90"
  "${XTB_ROOT}/src/coulomb/klopmanohno.f90"

  # Continuum solvation
  "${XTB_ROOT}/src/solv/gbobc.f90"
  "${XTB_ROOT}/src/solv/cm5.f90"

  # SCC
  # integrals and electrostatics
  "${XTB_ROOT}/src/intpack.f90"
  "${XTB_ROOT}/src/intgrad.f90"
  "${XTB_ROOT}/src/aespot.f90"

  # iterator + gradient
  "${XTB_ROOT}/src/scc_core.f90"
  "${XTB_ROOT}/src/grad_core.f90"
  "${XTB_ROOT}/src/scf_module.f90"

  # driver
  "${XTB_ROOT}/src/single.f90"

  # convergence accelerator
  "${XTB_ROOT}/src/broyden.f90"

  # misc
  "${XTB_ROOT}/src/pseudodiag.f90"
  "${XTB_ROOT}/src/embedding.f90"

  # vTB + EHT/GFN0/PEEQ
  "${XTB_ROOT}/src/peeq_module.f90"

  # periodic boundary conditions
  "${XTB_ROOT}/src/pbc_tools.f90"
  "${XTB_ROOT}/src/generate_wsc.f90"
  "${XTB_ROOT}/src/lidep.f90"

  # Misc/unsorted
  "${XTB_ROOT}/src/elem.f90"
  "${XTB_ROOT}/src/pqn.f90"
  "${XTB_ROOT}/src/lin_mod.f90"
  "${XTB_ROOT}/src/lin.f90"
  "${XTB_ROOT}/src/printmold.f90"
  "${XTB_ROOT}/src/surfac.f"
  "${XTB_ROOT}/src/esp.f"
  "${XTB_ROOT}/src/shiftlp.f90"
  "${XTB_ROOT}/src/blowsy.f90"
  "${XTB_ROOT}/src/onetri.f90"
  "${XTB_ROOT}/src/local.f90"
  "${XTB_ROOT}/src/makel.f90"
  "${XTB_ROOT}/src/foden.f90"
  "${XTB_ROOT}/src/hlex.f90"
  "${XTB_ROOT}/src/cube.f90"
  "${XTB_ROOT}/src/pbc.f90"
  "${XTB_ROOT}/src/qpot.f90"
  "${XTB_ROOT}/src/metadynamic.f90"
  "${XTB_ROOT}/src/ls_rmsd.f90"
  "${XTB_ROOT}/src/drsp.f"
  "${XTB_ROOT}/src/stm.f"
  "${XTB_ROOT}/src/axis_trafo.f90"
  "${XTB_ROOT}/src/coffee.f90"
  "${XTB_ROOT}/src/tmgrad.f90"
  "${XTB_ROOT}/src/thermo.f90"
  "${XTB_ROOT}/src/constr.f90"
  "${XTB_ROOT}/src/getsymnum.f90"
  "${XTB_ROOT}/src/solv/lebedev.f90"
  "${XTB_ROOT}/src/qmdff.f90"
  "${XTB_ROOT}/src/neighbor.f90"
  "${XTB_ROOT}/src/zmatpr.f90"
  "${XTB_ROOT}/src/fragment.f90"
  "${XTB_ROOT}/src/restart.f90"
  "${XTB_ROOT}/src/lopt.f90"
  "${XTB_ROOT}/src/qsort.f90"
  "${XTB_ROOT}/src/locmode.f90"
  "${XTB_ROOT}/src/intmodes.f90"
  "${XTB_ROOT}/src/wrmodef.f90"
  "${XTB_ROOT}/src/hessian.F90"
  "${XTB_ROOT}/src/timing.f90"
  "${XTB_ROOT}/src/prmat.f90"
  "${XTB_ROOT}/src/dipole.f90"
  "${XTB_ROOT}/src/dtrafo.f90"
  "${XTB_ROOT}/src/constrain_pot.f90"
  "${XTB_ROOT}/src/basic_geo.f90"
  "${XTB_ROOT}/src/approxrab.f90"

  # Ancopt
  "${XTB_ROOT}/src/geoopt_driver.f90"
  "${XTB_ROOT}/src/optimizer.f90"
  "${XTB_ROOT}/src/relaxation_engine.f90"
  "${XTB_ROOT}/src/bfgs.f90"
  "${XTB_ROOT}/src/david.f"
  "${XTB_ROOT}/src/david2.f90"
  "${XTB_ROOT}/src/bias_path.f90"
  "${XTB_ROOT}/src/geosum.f90"
  "${XTB_ROOT}/src/model_hessian.f90"
  "${XTB_ROOT}/src/lindh.f90"
  "${XTB_ROOT}/src/scan_driver.f90"
  "${XTB_ROOT}/src/symtranslib.f90"

  # Dynamic
  "${XTB_ROOT}/src/getname.f90"
  "${XTB_ROOT}/src/rmrottr.f90"
  "${XTB_ROOT}/src/matinv.f90"
  "${XTB_ROOT}/src/eqrot.f90"
  "${XTB_ROOT}/src/rmsd.f90"
  "${XTB_ROOT}/src/dynamic.f90"
  "${XTB_ROOT}/src/shake_module.f90"

  # Scan
  "${XTB_ROOT}/src/ifind.f90"
  "${XTB_ROOT}/src/spline2.f90"
  "${XTB_ROOT}/src/spline3.f90"
  "${XTB_ROOT}/src/anharmlib.f90"
  "${XTB_ROOT}/src/pocketscan.f90"
  "${XTB_ROOT}/src/screening.f90"
  "${XTB_ROOT}/src/modef.f90"
  "${XTB_ROOT}/src/mdoptim.f90"
  "${XTB_ROOT}/src/cqpath.f90"

  "${XTB_ROOT}/symmetry/symmetry.f90"
  "${XTB_ROOT}/symmetry/symmetry_i.c"

  # API
  "${XTB_ROOT}/src/calculator.f90"
  "${XTB_ROOT}/src/gfn0_calculator.f90"
  "${XTB_ROOT}/src/gfn1_calculator.f90"
  "${XTB_ROOT}/src/gfn2_calculator.f90"
  "${XTB_ROOT}/src/api/interface.f90"
  "${XTB_ROOT}/src/api/structs.f90"
  "${XTB_ROOT}/src/api/preload.f90"
  "${XTB_ROOT}/src/api/utils.f90"

  "${XTB_ROOT}/src/extern/mopac.f90"
  "${XTB_ROOT}/src/extern/orca.f90"
  "${XTB_ROOT}/src/extern/turbomole.f90"
)

set(XTB_F_TEST_SOURCES
  "${XTB_ROOT}/TESTSUITE/assertion.f90"
  "${XTB_ROOT}/TESTSUITE/coulomb.f90"
  "${XTB_ROOT}/TESTSUITE/dftd3.f90"
  "${XTB_ROOT}/TESTSUITE/dftd4.f90"
  "${XTB_ROOT}/TESTSUITE/eeq_model.f90"
  "${XTB_ROOT}/TESTSUITE/geometry_reader.f90"
  "${XTB_ROOT}/TESTSUITE/gfn0.f90"
  "${XTB_ROOT}/TESTSUITE/gfn1.f90"
  "${XTB_ROOT}/TESTSUITE/gfn2.f90"
  "${XTB_ROOT}/TESTSUITE/latticepoint.f90"
  "${XTB_ROOT}/TESTSUITE/pbc_tools.f90"
  "${XTB_ROOT}/TESTSUITE/peeq.f90"
  "${XTB_ROOT}/TESTSUITE/symmetry.f90"
  "${XTB_ROOT}/TESTSUITE/tbdef_atomlist.f90"
  "${XTB_ROOT}/TESTSUITE/tbdef_molecule.f90"
  "${XTB_ROOT}/TESTSUITE/tbdef_wsc.f90"
  "${XTB_ROOT}/TESTSUITE/tests_peeq.f90"
  "${XTB_ROOT}/TESTSUITE/thermo.f90"
)

set(XTB_C_TEST_SOURCES
  "${XTB_ROOT}/TESTSUITE/c_api_example.c"
)
