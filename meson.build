# This file is part of xtb.
#
# Copyright (C) 2019-2020 Sebastian Ehlert
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

# extended tight binding program package
project('xtb',
        'fortran', 'c',
        version: '6.2.2',
        license: 'LGPL3',
        meson_version: '>=0.49',
        default_options: ['includedir=include/xtb'])

fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')

if fc.get_id() != cc.get_id()
  warning('FC and CC are not from the same vendor')
endif

# build a configuration data containing all the important data to propagate
# it to the automatically generated files
config = configuration_data({
  'name': meson.project_name(),
  'description': 'Semiempirical Extended Tight-Binding Program Package',
  'version': meson.project_version(),
  'commit': run_command(find_program('git'),'show','-s','--format=%h').stdout().strip(),
  'date': run_command(find_program('date'),'-I').stdout().strip(),
  'author': run_command(find_program('whoami')).stdout().strip(),
  'origin': run_command(find_program('hostname')).stdout().strip(),
  'prefix': get_option('prefix'),
  'mandir': get_option('mandir'),
  'bindir': get_option('bindir'),
  'libdir': get_option('libdir'),
  'includedir': get_option('includedir'),
  'datadir': get_option('datadir')/meson.project_name(),
})
configure_file(input : 'xtb'/'version.f90', output : 'xtb_version.fh',
  configuration : config)

if fc.get_id() == 'gcc'
  add_project_arguments('-fdefault-real-8', language: 'fortran')
  add_project_arguments('-fdefault-double-8', language: 'fortran')
  add_project_arguments('-ffree-line-length-none', language: 'fortran')
  add_project_arguments('-fbacktrace', language: 'fortran')
elif fc.get_id() == 'intel'
  add_project_arguments('-axAVX2',    language: 'fortran')
  add_project_arguments('-r8',        language: 'fortran')
  add_project_arguments('-traceback', language: 'fortran')
  add_project_link_arguments('-static', language: 'fortran')
  add_project_link_arguments('-static', language: 'c') # icc will do linking
endif

if cc.get_id() == 'gcc'
  add_project_arguments('-std=c11', language: 'c')
elif cc.get_id() == 'intel'
  add_project_arguments('-std=c11', language: 'c')
endif

# fix compiliation problems with of symmetry/symmetry_i.c
add_project_arguments('-D_Float128=__float128', language: 'c')

## ========================================== ##
## LIBRARIES
## ========================================== ##
dependencies = []
dependencies_exe = []
dependencies_sha = []

if get_option('la_backend') == 'mkl'
  add_project_arguments('-DWITH_MKL', language: 'fortran')
  if fc.get_id() == 'gcc'
    libmkl_exe = [fc.find_library('mkl_gf_lp64')]
    libmkl_exe += fc.find_library('mkl_gnu_thread')
  elif fc.get_id() == 'intel'
    libmkl_exe = [fc.find_library('mkl_intel_lp64')]
    libmkl_exe += fc.find_library('mkl_intel_thread')
  endif
  libmkl_exe += fc.find_library('mkl_core')
  libmkl_exe += fc.find_library('iomp5')
  libmkl_sha = [fc.find_library('mkl_rt')]
  dependencies_sha += libmkl_sha
  dependencies_exe += libmkl_exe
elif get_option('la_backend') == 'openblas'
  dependencies += fc.find_library('openblas', required: true)
  dependencies += fc.find_library('lapack', required: true)
elif get_option('la_backend') == 'custom'
  foreach lib: get_option('custom_libraries')
    dependencies += fc.find_library(lib)
  endforeach
else
  dependencies += fc.find_library('blas', required: true)
  dependencies += fc.find_library('lapack', required: true)
endif

if get_option('openmp')
  if fc.get_id() == 'gcc'
    add_project_arguments('-fopenmp', language: 'fortran')
    add_project_link_arguments('-fopenmp', language: 'fortran')
    add_project_link_arguments('-fopenmp', language: 'c')
  elif fc.get_id() == 'intel'
    add_project_arguments('-qopenmp', language: 'fortran')
    add_project_link_arguments('-qopenmp', language: 'fortran')
    add_project_link_arguments('-qopenmp', language: 'c')
  endif
endif
dependencies += dependency('threads')

# distribute dependencies for shared object and static executable
dependencies_sha += dependencies
dependencies_exe += dependencies

## ========================================== ##
## SOURCE
## ========================================== ##
mctc_srcs  = []
xtb_srcs   = []
xtb_main   = 'xtb/program_main.f90'

# MCTC library
mctc_srcs += 'mctc/mctc_global.f90'
mctc_srcs += 'mctc/mctc_systools.f90'
mctc_srcs += 'mctc/mctc_strings.f90'
mctc_srcs += 'mctc/mctc_constants.f90'
mctc_srcs += 'mctc/mctc_econv.f90'
mctc_srcs += 'mctc/mctc_param.f90'
mctc_srcs += 'mctc/param/atomic_masses.f90'
mctc_srcs += 'mctc/param/chemical_hardnesses.f90'
mctc_srcs += 'mctc/param/covalent_radii.f90'
mctc_srcs += 'mctc/param/electronegativities.f90'
mctc_srcs += 'mctc/param/pse.f90'
mctc_srcs += 'mctc/param/r4r2_expectation_values.f90'
mctc_srcs += 'mctc/mctc_timings.f90'
mctc_srcs += 'mctc/mctc_filetools.f90'
mctc_srcs += 'mctc/mctc_la.f90'
mctc_srcs += 'mctc/mctc_init.f90'
mctc_srcs += 'mctc/mctc_resize_arrays.f90'
mctc_srcs += 'mctc/error.f90'
mctc_srcs += 'mctc/signal.c'

# class definitions
xtb_srcs += 'xtb/tbdef_setvar.f90'
xtb_srcs += 'xtb/tbdef_wavefunction.f90'
xtb_srcs += 'xtb/tbdef_basisset.f90'
xtb_srcs += 'xtb/tbdef_molecule.f90'
xtb_srcs += 'xtb/tbdef_solvent.f90'
xtb_srcs += 'xtb/tbdef_param.f90'
xtb_srcs += 'xtb/tbdef_data.f90'
xtb_srcs += 'xtb/tbdef_anc.f90'
xtb_srcs += 'xtb/tbdef_timer.f90'
xtb_srcs += 'xtb/tbdef_pcem.f90'
xtb_srcs += 'xtb/tbdef_wsc.f90'
xtb_srcs += 'xtb/tbdef_options.f90'
xtb_srcs += 'xtb/tbdef_calculator.f90'
xtb_srcs += 'xtb/tbdef_atomlist.f90'
xtb_srcs += 'xtb/tbdef_topology.f90'
xtb_srcs += 'xtb/tbdef_fragments.f90'
xtb_srcs += 'xtb/tbdef_buffer.f90'

# global data
xtb_srcs += 'xtb/gfn0param.f90'
xtb_srcs += 'xtb/sphereparam.f90'
xtb_srcs += 'xtb/scanparam.f90'
xtb_srcs += 'xtb/splitparam.f90'
xtb_srcs += 'xtb/symparam.f90'
xtb_srcs += 'xtb/fixparam.f90'
xtb_srcs += 'xtb/aoparam.f90'
xtb_srcs += 'xtb/setparam.f90'

# header and I/O
xtb_srcs += 'xtb/readin.f90'
xtb_srcs += 'xtb/filetools.f90'
xtb_srcs += 'xtb/set_module.f90'
xtb_srcs += 'xtb/constrain_param.f90'
xtb_srcs += 'xtb/argparser.F90'
xtb_srcs += 'xtb/header.f90'
xtb_srcs += 'xtb/printout.f90'
xtb_srcs += 'xtb/xhelp.f90'
xtb_srcs += 'xtb/expire.f90'
xtb_srcs += 'xtb/define.f90'
xtb_srcs += 'xtb/readl.f90'
xtb_srcs += 'xtb/readl2.f90'
xtb_srcs += 'xtb/rdcoord2.f90'
xtb_srcs += 'xtb/geometry_reader.f90'
xtb_srcs += 'xtb/wrbas.f90'
xtb_srcs += 'xtb/wrmo.f90'
xtb_srcs += 'xtb/main_property.f90'
xtb_srcs += 'xtb/main_json.f90'
xtb_srcs += 'xtb/main_geometry.f90'
xtb_srcs += 'xtb/gfn_prparam.f90'
xtb_srcs += 'xtb/wrgbw.f90'
xtb_srcs += 'xtb/enso_printout.f90'
xtb_srcs += 'xtb/read_gfn_param.f90'
xtb_srcs += 'xtb/molecule_reader.f90'
xtb_srcs += 'xtb/molecule_writer.f90'
xtb_srcs += 'xtb/tbmod_file_utils.f90'

# parameters
xtb_srcs += 'xtb/gfn_paramset.f90'
xtb_srcs += 'xtb/dftd3_parameters.f90'
xtb_srcs += 'xtb/setwll.f'
xtb_srcs += 'xtb/charge_model.f90'

# initial guess
xtb_srcs += 'xtb/iniq.f90'
xtb_srcs += 'xtb/eeq_model.f90'

# basis set
xtb_srcs += 'xtb/xbasis.f90'
xtb_srcs += 'xtb/printbas.f90'
xtb_srcs += 'xtb/gauss.f90'

# dispersion
xtb_srcs += 'xtb/ncoord.f90'
xtb_srcs += 'xtb/dftd3.f'
xtb_srcs += 'xtb/dftd4.f90'
xtb_srcs += 'xtb/dispplot.f90'

# continuum solvation
xtb_srcs += 'xtb/gbobc.f90'
xtb_srcs += 'xtb/cm5.f90'

# SCC
# integrals and electrostatics
xtb_srcs += 'xtb/intpack.f90'
xtb_srcs += 'xtb/intgrad.f90'
xtb_srcs += 'xtb/aespot.f90'
# iterator + gradient
xtb_srcs += 'xtb/scc_core.f90'
xtb_srcs += 'xtb/grad_core.f90'
xtb_srcs += 'xtb/scf_module.f90'
# driver
xtb_srcs += 'xtb/single.f90'
# convergence accelerator
xtb_srcs += 'xtb/broyden.f'
# misc
xtb_srcs += 'xtb/pseudodiag.f90'
xtb_srcs += 'xtb/embedding.f90'

# vTB + EHT/GFN0/PEEQ
xtb_srcs += 'xtb/peeq_module.f90'

# periodic boundary conditions
xtb_srcs += 'xtb/pbc_tools.f90'
xtb_srcs += 'xtb/generate_wsc.f90'
xtb_srcs += 'xtb/lidep.f90'

# misc/unsorted
xtb_srcs += 'xtb/elem.f90'
xtb_srcs += 'xtb/asym.f90'
xtb_srcs += 'xtb/pqn.f90'
xtb_srcs += 'xtb/lin_mod.f90'
xtb_srcs += 'xtb/lin.f90'
xtb_srcs += 'xtb/printmold.f'
xtb_srcs += 'xtb/surfac.f'
xtb_srcs += 'xtb/esp.f'
xtb_srcs += 'xtb/shiftlp.f'
xtb_srcs += 'xtb/blowsy.f90'
xtb_srcs += 'xtb/onetri.f'
xtb_srcs += 'xtb/exchange.f'
xtb_srcs += 'xtb/local.f90'
xtb_srcs += 'xtb/makel.f90'
xtb_srcs += 'xtb/foden.f90'
xtb_srcs += 'xtb/hlex.f90'
xtb_srcs += 'xtb/cube.f90'
xtb_srcs += 'xtb/pbc.f90'
xtb_srcs += 'xtb/qpot.f90'
xtb_srcs += 'xtb/metadynamic.f90'
xtb_srcs += 'xtb/ls_rmsd.f90'
xtb_srcs += 'xtb/drsp.f'
xtb_srcs += 'xtb/stm.f'
xtb_srcs += 'xtb/axis_trafo.f90'
xtb_srcs += 'xtb/coffee.f90'
xtb_srcs += 'xtb/write_geometry.f90'
xtb_srcs += 'xtb/tmgrad.f90'
xtb_srcs += 'xtb/thermo.f90'
xtb_srcs += 'xtb/constr.f'
xtb_srcs += 'xtb/getsymnum.f90'
xtb_srcs += 'xtb/grid_module.f90'
xtb_srcs += 'xtb/qmdff.f90'
xtb_srcs += 'xtb/qmexternal.f90'
xtb_srcs += 'xtb/qcextern.f90'
xtb_srcs += 'xtb/neighbor.f'
xtb_srcs += 'xtb/zmatpr.f'
xtb_srcs += 'xtb/fragment.f'
xtb_srcs += 'xtb/restart.f90'
xtb_srcs += 'xtb/lopt.f90'
xtb_srcs += 'xtb/qsort.f90'
xtb_srcs += 'xtb/locmode.f'
xtb_srcs += 'xtb/intmodes.f'
xtb_srcs += 'xtb/wrmodef.f'
xtb_srcs += 'xtb/hessian.F90'
xtb_srcs += 'xtb/xbond.f90'
xtb_srcs += 'xtb/timing.f90'
xtb_srcs += 'xtb/prmat.f'
xtb_srcs += 'xtb/dipole.f90'
xtb_srcs += 'xtb/dtrafo.f90'
xtb_srcs += 'xtb/constrain_pot.f90'
xtb_srcs += 'xtb/basic_geo.f90'
xtb_srcs += 'xtb/approxrab.f90'

# ancopt
xtb_srcs += 'xtb/geoopt_driver.f90'
xtb_srcs += 'xtb/optimizer.f90'
xtb_srcs += 'xtb/relaxation_engine.f90'
xtb_srcs += 'xtb/bfgs.f90'
xtb_srcs += 'xtb/david.f'
xtb_srcs += 'xtb/david2.f90'
xtb_srcs += 'xtb/bias_path.f90'
xtb_srcs += 'xtb/geosum.f'
xtb_srcs += 'xtb/model_hessian.f90'
xtb_srcs += 'xtb/lindh.f'
xtb_srcs += 'xtb/scan_driver.f90'
xtb_srcs += 'xtb/symtranslib.f'

# dynamic
xtb_srcs += 'xtb/getname.f'
xtb_srcs += 'xtb/rmrottr.f'
xtb_srcs += 'xtb/matinv.f'
xtb_srcs += 'xtb/eqrot.f'
xtb_srcs += 'xtb/rmsd.f'
xtb_srcs += 'xtb/dynamic.f90'
xtb_srcs += 'xtb/shake_module.f90'

# scan
xtb_srcs += 'xtb/ifind.f'
xtb_srcs += 'xtb/spline2.f90'
xtb_srcs += 'xtb/spline3.f90'
xtb_srcs += 'xtb/anharmlib.f'
xtb_srcs += 'xtb/pocketscan.f'
xtb_srcs += 'xtb/screening.f90'
xtb_srcs += 'xtb/modef.f90'
xtb_srcs += 'xtb/mdoptim.f90'
xtb_srcs += 'xtb/cqpath.f90'

xtb_srcs += 'symmetry/symmetry.f90'
xtb_srcs += 'symmetry/symmetry_i.c'

# API
xtb_srcs += 'xtb/calculator.f90'
xtb_srcs += 'xtb/gfn0_calculator.f90'
xtb_srcs += 'xtb/gfn1_calculator.f90'
xtb_srcs += 'xtb/gfn2_calculator.f90'
xtb_srcs += 'xtb/c_api.f90'

incdir = include_directories('include')

xtb_srcs += mctc_srcs

## ========================================== ##
## TARGETS
## ========================================== ##
# create a static library from all sources
xtb_lib_static = static_library(meson.project_name(),
                                xtb_srcs,
                                include_directories: incdir,
                                install: true,
                                pic: true)

xtb_exe = executable(meson.project_name(),
                     xtb_main,
                     dependencies: dependencies_exe,
                     include_directories: incdir,
                     link_with: xtb_lib_static,
                     install: true)

xtb_lib_shared = shared_library(meson.project_name(),
                                version: meson.project_version(),
                                dependencies: dependencies_sha,
                                include_directories: incdir,
                                link_whole: xtb_lib_static,
                                link_args: '-fopenmp',
                                install: true)

xtb_dep_static = [declare_dependency(link_whole: xtb_lib_static), dependencies_exe]
xtb_dep = [declare_dependency(link_with: xtb_lib_shared), dependencies_sha]

## ========================================== ##
## INSTALL
## ========================================== ##
install_headers('include/xtb.h')

asciidoc = find_program('a2x', required: false)
if asciidoc.found()
  run_command(asciidoc, 'man/man1/xtb.1.txt', '--format=manpage')
  run_command(asciidoc, 'man/man7/xcontrol.7.txt', '--format=manpage')
  install_man(['man/man1/xtb.1', 'man/man7/xcontrol.7'])
endif

xtb_parameter_files = ['.param_gfn0.xtb',
                       '.param_gfn2.xtb',
                       '.param_gfn.xtb']

install_data(xtb_parameter_files)

if get_option('install_modules')
  xtb_objs = []
  foreach file: xtb_srcs
    xtb_objs += '@0@.o'.format('_'.join(file.split('/')))
  endforeach
  install_subdir(meson.current_build_dir()/'@0@@sta'.format(meson.project_name()),
                 install_dir: 'include',
                 strip_directory: true,
                 exclude_files: xtb_objs)
endif

# additional misc. files
subdir('assets')

## ========================================== ##
## TESTSUITE
## ========================================== ##
# make sure the correct library is loaded
xtbenv = environment()
xtbenv.prepend('LD_LIBRARY_PATH', meson.current_build_dir())
xtbenv.prepend('PYTHONPATH', meson.current_source_dir() / 'python')
xtbenv.set('XTBPATH', meson.current_source_dir())

# all tests are defined in a separate meson.build
subdir('TESTSUITE')

py_xtb = meson.current_source_dir() / 'python' / 'xtb'

pytest = find_program('pytest', required: false)
if pytest.found()
  test('pytest: xtb.py', pytest, args: ['--pyargs', 'xtb'], env: xtbenv)
endif

mypy = find_program('mypy', required: false)
if mypy.found()
  test('mypy: xtb.py', mypy, args: [py_xtb, '--ignore-missing-imports'], env: xtbenv)
endif

pylint = find_program('pylint', required: false)
if pylint.found()
  test('pylint: xtb.py', pylint, args: [py_xtb, '--disable=invalid-names,duplicate-code'], env: xtbenv)
endif

flake8 = find_program('flake8', required: false)
if flake8.found()
  test('flake8: xtb.py', flake8, args: [py_xtb, '--max-line-length=83'], env: xtbenv)
endif
