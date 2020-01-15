cmake_minimum_required(VERSION 3.9)
set(XTB_ROOT ${PROJECT_SOURCE_DIR})
set(XTB_SOURCES
  # MCTC library
  "${XTB_ROOT}/mctc/mctc_global.f90"
  "${XTB_ROOT}/mctc/mctc_systools.f90"
  "${XTB_ROOT}/mctc/mctc_strings.f90"
  "${XTB_ROOT}/mctc/mctc_constants.f90"
  "${XTB_ROOT}/mctc/mctc_econv.f90"
  "${XTB_ROOT}/mctc/mctc_param.f90"
  "${XTB_ROOT}/mctc/param/atomic_masses.f90"
  "${XTB_ROOT}/mctc/param/chemical_hardnesses.f90"
  "${XTB_ROOT}/mctc/param/covalent_radii.f90"
  "${XTB_ROOT}/mctc/param/electronegativities.f90"
  "${XTB_ROOT}/mctc/param/pse.f90"
  "${XTB_ROOT}/mctc/param/r4r2_expectation_values.f90"
  "${XTB_ROOT}/mctc/mctc_timings.f90"
  "${XTB_ROOT}/mctc/mctc_filetools.f90"
  "${XTB_ROOT}/mctc/mctc_la.f90"
  "${XTB_ROOT}/mctc/mctc_init.f90"
  "${XTB_ROOT}/mctc/mctc_resize_arrays.f90"
  "${XTB_ROOT}/mctc/error.f90"
  "${XTB_ROOT}/mctc/signal.c"

  # Class definitions
  "${XTB_ROOT}/xtb/tbdef_setvar.f90"
  "${XTB_ROOT}/xtb/tbdef_wavefunction.f90"
  "${XTB_ROOT}/xtb/tbdef_basisset.f90"
  "${XTB_ROOT}/xtb/tbdef_molecule.f90"
  "${XTB_ROOT}/xtb/tbdef_solvent.f90"
  "${XTB_ROOT}/xtb/tbdef_param.f90"
  "${XTB_ROOT}/xtb/tbdef_data.f90"
  "${XTB_ROOT}/xtb/tbdef_anc.f90"
  "${XTB_ROOT}/xtb/tbdef_timer.f90"
  "${XTB_ROOT}/xtb/tbdef_pcem.f90"
  "${XTB_ROOT}/xtb/tbdef_wsc.f90"
  "${XTB_ROOT}/xtb/tbdef_options.f90"
  "${XTB_ROOT}/xtb/tbdef_calculator.f90"
  "${XTB_ROOT}/xtb/tbdef_atomlist.f90"
  "${XTB_ROOT}/xtb/tbdef_topology.f90"
  "${XTB_ROOT}/xtb/tbdef_fragments.f90"
  "${XTB_ROOT}/xtb/tbdef_buffer.f90"

  # Global data
  "${XTB_ROOT}/xtb/gfn0param.f90"
  "${XTB_ROOT}/xtb/sphereparam.f90"
  "${XTB_ROOT}/xtb/scanparam.f90"
  "${XTB_ROOT}/xtb/splitparam.f90"
  "${XTB_ROOT}/xtb/symparam.f90"
  "${XTB_ROOT}/xtb/fixparam.f90"
  "${XTB_ROOT}/xtb/aoparam.f90"
  "${XTB_ROOT}/xtb/setparam.f90"

  # Header and I/O
  "${XTB_ROOT}/xtb/readin.f90"
  "${XTB_ROOT}/xtb/filetools.f90"
  "${XTB_ROOT}/xtb/set_module.f90"
  "${XTB_ROOT}/xtb/constrain_param.f90"
  "${XTB_ROOT}/xtb/argparser.F90"
  "${XTB_ROOT}/xtb/header.f90"
  "${XTB_ROOT}/xtb/printout.f90"
  "${XTB_ROOT}/xtb/xhelp.f90"
  "${XTB_ROOT}/xtb/expire.f90"
  "${XTB_ROOT}/xtb/define.f90"
  "${XTB_ROOT}/xtb/readl.f90"
  "${XTB_ROOT}/xtb/readl2.f90"
  "${XTB_ROOT}/xtb/rdcoord2.f90"
  "${XTB_ROOT}/xtb/geometry_reader.f90"
  "${XTB_ROOT}/xtb/wrbas.f90"
  "${XTB_ROOT}/xtb/wrmo.f90"
  "${XTB_ROOT}/xtb/main_property.f90"
  "${XTB_ROOT}/xtb/main_json.f90"
  "${XTB_ROOT}/xtb/main_geometry.f90"
  "${XTB_ROOT}/xtb/gfn_prparam.f90"
  "${XTB_ROOT}/xtb/wrgbw.f90"
  "${XTB_ROOT}/xtb/enso_printout.f90"
  "${XTB_ROOT}/xtb/read_gfn_param.f90"
  "${XTB_ROOT}/xtb/molecule_reader.f90"
  "${XTB_ROOT}/xtb/molecule_writer.f90"
  "${XTB_ROOT}/xtb/tbmod_file_utils.f90"

  # Parameters
  "${XTB_ROOT}/xtb/gfn_paramset.f90"
  "${XTB_ROOT}/xtb/dftd3_parameters.f90"
  "${XTB_ROOT}/xtb/setwll.f"
  "${XTB_ROOT}/xtb/charge_model.f90"

  # Initial guess
  "${XTB_ROOT}/xtb/iniq.f90"
  "${XTB_ROOT}/xtb/eeq_model.f90"

  # Basis set
  "${XTB_ROOT}/xtb/xbasis.f90"
  "${XTB_ROOT}/xtb/printbas.f90"
  "${XTB_ROOT}/xtb/gauss.f90"

  # Dispersion
  "${XTB_ROOT}/xtb/ncoord.f90"
  "${XTB_ROOT}/xtb/dftd3.f"
  "${XTB_ROOT}/xtb/dftd4.f90"
  "${XTB_ROOT}/xtb/dispplot.f90"

  # Continuum solvation
  "${XTB_ROOT}/xtb/gbobc.f90"
  "${XTB_ROOT}/xtb/cm5.f90"

  # SCC
  # integrals and electrostatics
  "${XTB_ROOT}/xtb/intpack.f90"
  "${XTB_ROOT}/xtb/intgrad.f90"
  "${XTB_ROOT}/xtb/aespot.f90"

  # iterator + gradient
  "${XTB_ROOT}/xtb/scc_core.f90"
  "${XTB_ROOT}/xtb/grad_core.f90"
  "${XTB_ROOT}/xtb/scf_module.f90"

  # driver
  "${XTB_ROOT}/xtb/single.f90"

  # convergence accelerator
  "${XTB_ROOT}/xtb/broyden.f"

  # misc
  "${XTB_ROOT}/xtb/pseudodiag.f90"
  "${XTB_ROOT}/xtb/embedding.f90"

  # vTB + EHT/GFN0/PEEQ
  "${XTB_ROOT}/xtb/peeq_module.f90"

  # periodic boundary conditions
  "${XTB_ROOT}/xtb/pbc_tools.f90"
  "${XTB_ROOT}/xtb/generate_wsc.f90"
  "${XTB_ROOT}/xtb/lidep.f90"

  # Misc/unsorted
  "${XTB_ROOT}/xtb/elem.f90"
  "${XTB_ROOT}/xtb/asym.f90"
  "${XTB_ROOT}/xtb/pqn.f90"
  "${XTB_ROOT}/xtb/lin_mod.f90"
  "${XTB_ROOT}/xtb/lin.f90"
  "${XTB_ROOT}/xtb/printmold.f"
  "${XTB_ROOT}/xtb/surfac.f"
  "${XTB_ROOT}/xtb/esp.f"
  "${XTB_ROOT}/xtb/shiftlp.f"
  "${XTB_ROOT}/xtb/blowsy.f90"
  "${XTB_ROOT}/xtb/onetri.f"
  "${XTB_ROOT}/xtb/exchange.f"
  "${XTB_ROOT}/xtb/local.f90"
  "${XTB_ROOT}/xtb/makel.f90"
  "${XTB_ROOT}/xtb/foden.f90"
  "${XTB_ROOT}/xtb/hlex.f90"
  "${XTB_ROOT}/xtb/cube.f90"
  "${XTB_ROOT}/xtb/pbc.f90"
  "${XTB_ROOT}/xtb/qpot.f90"
  "${XTB_ROOT}/xtb/metadynamic.f90"
  "${XTB_ROOT}/xtb/ls_rmsd.f90"
  "${XTB_ROOT}/xtb/drsp.f"
  "${XTB_ROOT}/xtb/stm.f"
  "${XTB_ROOT}/xtb/axis_trafo.f90"
  "${XTB_ROOT}/xtb/coffee.f90"
  "${XTB_ROOT}/xtb/write_geometry.f90"
  "${XTB_ROOT}/xtb/tmgrad.f90"
  "${XTB_ROOT}/xtb/thermo.f90"
  "${XTB_ROOT}/xtb/constr.f"
  "${XTB_ROOT}/xtb/getsymnum.f90"
  "${XTB_ROOT}/xtb/grid_module.f90"
  "${XTB_ROOT}/xtb/qmdff.f90"
  "${XTB_ROOT}/xtb/qmexternal.f90"
  "${XTB_ROOT}/xtb/qcextern.f90"
  "${XTB_ROOT}/xtb/neighbor.f"
  "${XTB_ROOT}/xtb/zmatpr.f"
  "${XTB_ROOT}/xtb/fragment.f"
  "${XTB_ROOT}/xtb/restart.f90"
  "${XTB_ROOT}/xtb/lopt.f90"
  "${XTB_ROOT}/xtb/qsort.f90"
  "${XTB_ROOT}/xtb/locmode.f"
  "${XTB_ROOT}/xtb/intmodes.f"
  "${XTB_ROOT}/xtb/wrmodef.f"
  "${XTB_ROOT}/xtb/hessian.F90"
  "${XTB_ROOT}/xtb/xbond.f90"
  "${XTB_ROOT}/xtb/timing.f90"
  "${XTB_ROOT}/xtb/prmat.f"
  "${XTB_ROOT}/xtb/dipole.f90"
  "${XTB_ROOT}/xtb/dtrafo.f90"
  "${XTB_ROOT}/xtb/constrain_pot.f90"
  "${XTB_ROOT}/xtb/basic_geo.f90"
  "${XTB_ROOT}/xtb/approxrab.f90"

  # Ancopt
  "${XTB_ROOT}/xtb/geoopt_driver.f90"
  "${XTB_ROOT}/xtb/optimizer.f90"
  "${XTB_ROOT}/xtb/relaxation_engine.f90"
  "${XTB_ROOT}/xtb/bfgs.f90"
  "${XTB_ROOT}/xtb/david.f"
  "${XTB_ROOT}/xtb/david2.f90"
  "${XTB_ROOT}/xtb/bias_path.f90"
  "${XTB_ROOT}/xtb/geosum.f"
  "${XTB_ROOT}/xtb/model_hessian.f90"
  "${XTB_ROOT}/xtb/lindh.f"
  "${XTB_ROOT}/xtb/scan_driver.f90"
  "${XTB_ROOT}/xtb/symtranslib.f"

  # Dynamic
  "${XTB_ROOT}/xtb/getname.f"
  "${XTB_ROOT}/xtb/rmrottr.f"
  "${XTB_ROOT}/xtb/matinv.f"
  "${XTB_ROOT}/xtb/eqrot.f"
  "${XTB_ROOT}/xtb/rmsd.f"
  "${XTB_ROOT}/xtb/dynamic.f90"
  "${XTB_ROOT}/xtb/shake_module.f90"

  # Scan
  "${XTB_ROOT}/xtb/ifind.f"
  "${XTB_ROOT}/xtb/spline2.f90"
  "${XTB_ROOT}/xtb/spline3.f90"
  "${XTB_ROOT}/xtb/anharmlib.f"
  "${XTB_ROOT}/xtb/pocketscan.f"
  "${XTB_ROOT}/xtb/screening.f90"
  "${XTB_ROOT}/xtb/modef.f90"
  "${XTB_ROOT}/xtb/mdoptim.f90"
  "${XTB_ROOT}/xtb/cqpath.f90"

  "${XTB_ROOT}/symmetry/symmetry.f90"
  "${XTB_ROOT}/symmetry/symmetry_i.c"

  # API
  "${XTB_ROOT}/xtb/calculator.f90"
  "${XTB_ROOT}/xtb/gfn0_calculator.f90"
  "${XTB_ROOT}/xtb/gfn1_calculator.f90"
  "${XTB_ROOT}/xtb/gfn2_calculator.f90"
  "${XTB_ROOT}/xtb/c_api.f90"
)

set(XTB_F_TEST_SOURCES
  "${XTB_ROOT}/TESTSUITE/assertion.f90"
  "${XTB_ROOT}/TESTSUITE/dftd4.f90"
  "${XTB_ROOT}/TESTSUITE/eeq_model.f90"
  "${XTB_ROOT}/TESTSUITE/geometry_reader.f90"
  "${XTB_ROOT}/TESTSUITE/gfn0.f90"
  "${XTB_ROOT}/TESTSUITE/gfn1.f90"
  "${XTB_ROOT}/TESTSUITE/gfn2.f90"
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
