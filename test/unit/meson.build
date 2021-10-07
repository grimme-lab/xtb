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

# Create test-drive as subproject for testing
testdrive_prj = subproject(
  'test-drive',
  version: '>=0.3',
  required: not meson.is_subproject(),
  default_options: [
    'default_library=static',
    'openmp=false',
  ],
)
# If we do not find test-drive and are a subproject, we just skip testing
if not testdrive_prj.found()
  subdir_done()
endif
testdrive_dep = testdrive_prj.get_variable('testdrive_dep')

tests = [
  'coulomb',
  'dftd3',
  'dftd4',
  'eeq',
  'geometry-reader',
  'hessian',
  'molecule',
  'thermo',
  'wsc',
]

test_srcs = files(
  'main.f90',
)
foreach t : tests
  test_srcs += files('test_@0@.f90'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [xtb_dep_static, testdrive_dep],
  link_language: 'fortran',
)

foreach t : tests
  test(t, tester, args: t, suite: 'unit', env: xtbenv)
endforeach

xtb_test = files(
  'assertion.f90',
  'molstock.f90',
  'coordinationnumber.f90',
  'gfn0.f90',
  'gfn1.f90',
  'gfn2.f90',
  'gfnff.f90',
  'latticepoint.f90',
  'pbc_tools.f90',
  'peeq.f90',
  'repulsion.f90',
  'symmetry.f90',
  'tests_peeq.f90',
  'tbdef_atomlist.f90',
)

xtb_test = executable(
  'xtb_test',
  sources: xtb_test,
  dependencies: [xtb_dep_static, testdrive_dep],
  link_language: 'fortran',
)

# some very basic checks to see if the executable reacts at all
test('Argparser print version', xtb_exe, args: '--version', env: xtbenv)
test('Argparser print help', xtb_exe, args: '--help', env: xtbenv)
test('Argparser print license', xtb_exe, args: '--license', env: xtbenv)
test('Argparser no arguments', xtb_exe, should_fail: true, env: xtbenv)
test('Info', xtb_exe, env: xtbenv, args: [
  'info',
  files('../../assets/inputs/coord/caffeine.coord'),
  files('../../assets/inputs/coord/quartz.3d.coord'),
  files('../../assets/inputs/vasp/ammonia.vasp'),
  files('../../assets/inputs/xyz/taxol.xyz'),
])

test_timeout = get_option('test_timeout')

test('Singlepoint', xtb_exe, args: ['--coffee', '--strict', '--norestart', '--namespace', 'test1'], env: xtbenv, timeout: test_timeout)
if get_option('la_backend') != 'netlib'
  test('Geometry opt.', xtb_exe, args: ['--coffee', '--opt', '--strict', '--norestart', '--namespace', 'test2'], env: xtbenv, timeout: test_timeout)
endif
test('IP/EA', xtb_exe, args: ['--coffee', '--gfn', '2', '--vipea', '--strict', '--norestart', '--namespace', 'test3'], env: xtbenv, timeout: test_timeout)
test('GFN0-xTB', xtb_exe, args: ['--coffee', '--gfn', '0', '--strict', '--norestart', '--namespace', 'test4'], env: xtbenv, timeout: test_timeout)
test('GFN1-xTB', xtb_exe, args: ['--coffee', '--gfn', '1', '--strict', '--norestart', '--namespace', 'test5'], env: xtbenv, timeout: test_timeout)
test('GFN2-xTB/GBSA', xtb_exe, args: ['--coffee', '--gfn', '2', '--strict', '--gbsa', 'h2o', '--norestart', '--namespace', 'test6'], env: xtbenv, timeout: test_timeout)
test('GFN-FF', xtb_exe, args: ['--coffee', '--gfnff', '--strict', '--norestart', '--namespace', 'test7'], env: xtbenv, timeout: test_timeout)

# more specific tests are implemented by the tester binary
test('Lattice points PBC-3D', xtb_test, args: ['latticepoint', 'pbc3d'], env: xtbenv, timeout: test_timeout)

if fc.get_id() == 'intel'
  test('IO atom list', xtb_test, args: ['xtb_type_atomlist', 'list'], env: xtbenv, timeout: test_timeout)
endif

test('PBC tools convert',xtb_test, args: ['pbc_tools', 'convert'], env: xtbenv, timeout: test_timeout)
test('PBC tools cutoff', xtb_test, args: ['pbc_tools', 'cutoff'], env: xtbenv, timeout: test_timeout)

test('Symmetry Water', xtb_test, args: ['symmetry', 'water'], env: xtbenv, timeout: test_timeout)
test('Symmetry Li8', xtb_test, args: ['symmetry', 'li8'], env: xtbenv, timeout: test_timeout)
test('Symmetry PCl3', xtb_test, args: ['symmetry', 'pcl3'], env: xtbenv, timeout: test_timeout)
test('Symmetry large', xtb_test, args: ['symmetry', 'c20'], env: xtbenv, timeout: test_timeout)

test('Repulsion 0D', xtb_test, args: ['repulsion', 'cluster'], env: xtbenv, timeout: test_timeout)
test('Repulsion 3D', xtb_test, args: ['repulsion', 'pbc3d'], env: xtbenv, timeout: test_timeout)

test('CN latticepoints', xtb_test, args: ['ncoord', 'pbc3dlatp'], env: xtbenv, timeout: test_timeout)
test('CN neighbourlist', xtb_test, args: ['ncoord', 'pbc3dneighs'], env: xtbenv, timeout: test_timeout)

test('GFN2-xTB Basic', xtb_test, args: ['gfn2', 'basic'], env: xtbenv, timeout: test_timeout)
test('GFN2-xTB Solvation', xtb_test, args: ['gfn2', 'solvation'], env: xtbenv, timeout: test_timeout)
test('GFN2-xTB SCC', xtb_test, args: ['gfn2', 'scc'], env: xtbenv, timeout: test_timeout)
test('GFN2-xTB API', xtb_test, args: ['gfn2', 'api'], env: xtbenv, timeout: test_timeout)
test('GFN2-xTB API (GBSA)', xtb_test, args: ['gfn2', 'gbsa'], env: xtbenv, timeout: test_timeout)
test('GFN2-xTB API (GBSA+salt)', xtb_test, args: ['gfn2', 'salt'], env: xtbenv, timeout: test_timeout)
test('GFN2-xTB d-Metal', xtb_test, args: ['gfn2', 'd-metal'], env: xtbenv, timeout: test_timeout)
test('GFN2-xTB COSMO', xtb_test, args: ['gfn2', 'cosmo'], env: xtbenv, timeout: test_timeout)
#test('GFN2-xTB API (PCEM)', xtb_test, args: ['gfn2', 'pcem'], env: xtbenv, timeout: test_timeout)

test('GFN1-xTB Basic', xtb_test, args: ['gfn1', 'basic'], env: xtbenv, timeout: test_timeout)
test('GFN1-xTB Solvation', xtb_test, args: ['gfn1', 'solvation'], env: xtbenv, timeout: test_timeout)
test('GFN1-xTB SCC', xtb_test, args: ['gfn1', 'scc'], env: xtbenv, timeout: test_timeout)
test('GFN1-xTB API', xtb_test, args: ['gfn1', 'api'], env: xtbenv, timeout: test_timeout)
test('GFN1-xTB XB', xtb_test, args: ['gfn1', 'xb'], env: xtbenv, timeout: test_timeout)
test('GFN1-xTB PBC', xtb_test, args: ['gfn1', 'pbc3d'], env: xtbenv, timeout: test_timeout)
test('GFN1-xTB API (GBSA)', xtb_test, args: ['gfn1', 'gbsa'], env: xtbenv, timeout: test_timeout)
test('GFN1-xTB API (PCEM)', xtb_test, args: ['gfn1', 'pcem'], env: xtbenv, timeout: test_timeout)
test('GFN1-xTB IPEA', xtb_test, args: ['gfn1', 'ipea'], env: xtbenv, timeout: test_timeout)
test('GFN1-xTB COSMO', xtb_test, args: ['gfn1', 'cosmo'], env: xtbenv, timeout: test_timeout)

test('GFN0-xTB Basic', xtb_test, args: ['gfn0', 'basic'], env: xtbenv, timeout: test_timeout)
test('GFN0-xTB Solvation', xtb_test, args: ['gfn0', 'solvation'], env: xtbenv, timeout: test_timeout)
test('GFN0-xTB SP', xtb_test, args: ['gfn0', 'sp'], env: xtbenv, timeout: test_timeout)
test('GFN0-xTB API', xtb_test, args: ['gfn0', 'api'], env: xtbenv, timeout: test_timeout)
test('GFN0-xTB SRB', xtb_test, args: ['gfn0', 'srb'], env: xtbenv, timeout: test_timeout)
test('GFN0-xTB SP  (PBC)', xtb_test, args: ['peeq', 'sp'], env: xtbenv, timeout: test_timeout)
test('GFN0-xTB API (PBC)', xtb_test, args: ['peeq', 'api'], env: xtbenv, timeout: test_timeout)
test('GFN0-xTB SRB (PBC)', xtb_test, args: ['peeq', 'srb'], env: xtbenv, timeout: test_timeout)

test('GFN-FF Basic', xtb_test, args: ['gfnff', 'basic'], env: xtbenv, timeout: test_timeout)
test('GFN-FF Solvation', xtb_test, args: ['gfnff', 'solvation'], env: xtbenv, timeout: test_timeout)
test('GFN-FF SP', xtb_test, args: ['gfnff', 'sp'], env: xtbenv, timeout: test_timeout)
test('GFN-FF HB', xtb_test, args: ['gfnff', 'hb'], env: xtbenv, timeout: test_timeout)
test('GFN-FF GBSA', xtb_test, args: ['gfnff', 'gbsa'], env: xtbenv, timeout: test_timeout)
test('GFN-FF Scale-up', xtb_test, args: ['gfnff', 'scaleup'], env: xtbenv, timeout: test_timeout)
test('GFN-FF PDB', xtb_test, args: ['gfnff', 'pdb'], env: xtbenv, timeout: test_timeout)
test('GFN-FF SDF', xtb_test, args: ['gfnff', 'sdf'], env: xtbenv, timeout: test_timeout)
