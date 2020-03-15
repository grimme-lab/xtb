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

xtb_test = [files('tests_peeq.f90')]
xtb_test += files('assertion.f90')
xtb_test += files('tbdef_molecule.f90')
xtb_test += files('latticepoint.f90')
xtb_test += files('tbdef_wsc.f90')
xtb_test += files('geometry_reader.f90')
xtb_test += files('eeq_model.f90')
xtb_test += files('pbc_tools.f90')
xtb_test += files('dftd4.f90')
xtb_test += files('gfn2.f90')
xtb_test += files('gfn1.f90')
xtb_test += files('gfn0.f90')
xtb_test += files('peeq.f90')
xtb_test += files('symmetry.f90')
xtb_test += files('thermo.f90')
xtb_test += files('tbdef_atomlist.f90')

xtb_test = executable('xtb_test',
                      xtb_test,
                      dependencies: dependencies_exe,
                      include_directories: incdir,
                      link_with: xtb_lib_static)

# some very basic checks to see if the executable reacts at all
test('Argparser: print version', xtb_exe, args: '--version', env: xtbenv)
test('Argparser: print help', xtb_exe, args: '--help', env: xtbenv)
test('Argparser: print license', xtb_exe, args: '--license', env: xtbenv)
test('Argparser: no arguments', xtb_exe, should_fail: true, env: xtbenv)

test('Singlepoint', xtb_exe, args: ['--coffee', '--strict', '--norestart'], env: xtbenv)
if get_option('la_backend') != 'netlib'
  test('Geometry opt.', xtb_exe, args: ['--coffee', '--opt', '--strict', '--norestart'], env: xtbenv)
endif
test('IP/EA', xtb_exe, args: ['--coffee', '--gfn', '2', '--vipea', '--strict', '--norestart'], env: xtbenv)
test('GFN0-xTB', xtb_exe, args: ['--coffee', '--gfn', '0', '--strict', '--norestart'], env: xtbenv)
test('GFN1-xTB', xtb_exe, args: ['--coffee', '--gfn', '1', '--strict', '--norestart'], env: xtbenv)
test('GFN2-xTB/GBSA', xtb_exe, args: ['--coffee', '--gfn', '2', '--strict', '--gbsa', 'h2o', '--norestart'], env: xtbenv)

# more specific tests are implemented by the tester binary
test('Molecule: axis', xtb_test, args: ['xtb_type_molecule', 'axis'], env: xtbenv)
test('Molecule: MIC', xtb_test, args: ['xtb_type_molecule', 'mic'], env: xtbenv)
test('Lattice points: PBC-3D', xtb_test, args: ['latticepoint', 'pbc3d'], env: xtbenv)

if fc.get_id() == 'intel'
  test('Wigner-Seitz Cell (0D)', xtb_test, args: ['xtb_type_wsc', '0d'], env: xtbenv)
endif
test('Wigner-Seitz Cell (3D)', xtb_test, args: ['xtb_type_wsc', '3d'], env: xtbenv)

if fc.get_id() == 'intel'
  test('IO: atom list', xtb_test, args: ['xtb_type_atomlist', 'list'], env: xtbenv)
endif

test('coord 3D', xtb_test, args: ['geometry_reader', 'coord_3d_a'], env: xtbenv)
test('coord 3D', xtb_test, args: ['geometry_reader', 'coord_3d_b'], env: xtbenv)
test('coord 2D', xtb_test, args: ['geometry_reader', 'coord_2d'], env: xtbenv)
test('coord 1D', xtb_test, args: ['geometry_reader', 'coord_1d'], env: xtbenv)
test('coord 0D', xtb_test, args: ['geometry_reader', 'coord_0d'], env: xtbenv)
test('Xmol  0D', xtb_test, args: ['geometry_reader', 'xmol_0d'], env: xtbenv)
test('POSCAR', xtb_test, args: ['geometry_reader', 'poscar_3d'], env: xtbenv)
test('molfile', xtb_test, args: ['geometry_reader', 'molfile'], env: xtbenv)
test('molfile flat', xtb_test, should_fail: true, args: ['geometry_reader', 'molfile_flat'], env: xtbenv)
test('SDF', xtb_test, args: ['geometry_reader', 'sdfile'], env: xtbenv)
test('SDF flat', xtb_test, should_fail: true, args: ['geometry_reader', 'sdfile_flat'], env: xtbenv)
test('SDF no H', xtb_test, should_fail: true, args: ['geometry_reader', 'sdfile_noh'], env: xtbenv)
test('PDB', xtb_test, args: ['geometry_reader', 'pdb'], env: xtbenv)
test('PDB no H', xtb_test, should_fail: true, args: ['geometry_reader', 'pdb_noh'], env: xtbenv)
test('genFormat', xtb_test, args: ['geometry_reader', 'gen'], env: xtbenv)

test('PBC tools: convert',xtb_test, args: ['pbc_tools', 'convert'], env: xtbenv)
test('PBC tools: cutoff', xtb_test, args: ['pbc_tools', 'cutoff'], env: xtbenv)

test('Symmetry: Water', xtb_test, args: ['symmetry', 'water'], env: xtbenv)
test('Symmetry: Li8', xtb_test, args: ['symmetry', 'li8'], env: xtbenv)
test('Symmetry: PCl3', xtb_test, args: ['symmetry', 'pcl3'], env: xtbenv)
test('Symmetry: large', xtb_test, args: ['symmetry', 'c20'], env: xtbenv)

test('Thermo: axis', xtb_test, args: ['thermo', 'axis'], env: xtbenv)
test('Thermo: calculation', xtb_test, args: ['thermo', 'calc'], env: xtbenv)
test('Thermo: print', xtb_test, args: ['thermo', 'print'], env: xtbenv)

test('EEQ model: water', xtb_test,args: ['eeq_model', 'water'], env: xtbenv)
test('EEQ model: 3D Ewald', xtb_test,args: ['eeq_model', 'ewald'], env: xtbenv)
test('EEQ model: GBSA', xtb_test,args: ['eeq_model', 'gbsa'], env: xtbenv)
test('EEQ model: GBSA (salt)', xtb_test,args: ['eeq_model', 'salt'], env: xtbenv)
test('EEQ model: GBSA (H-bond)', xtb_test,args: ['eeq_model', 'hbond'], env: xtbenv)

test('Dispersion: properties', xtb_test,args: ['dftd4', 'properties'], env: xtbenv)
test('Dispersion: energies', xtb_test,args: ['dftd4', 'energies'], env: xtbenv)
test('Dispersion: energies (PBC)', xtb_test,args: ['dftd4', 'pbc_disp'], env: xtbenv)
#test('Dispersion: cell gradients', xtb_test,args: ['dftd4', 'cell_gradient'] env: xtbenv)
test('Dispersion: API', xtb_test,args: ['dftd4', 'api'], env: xtbenv)
#test('Dispersion: API (PBC)', xtb_test,args: ['dftd4', 'pbc_api'] env: xtbenv)

test('GFN2-xTB: SCC', xtb_test, args: ['gfn2', 'scc'], env: xtbenv)
test('GFN2-xTB: API', xtb_test, args: ['gfn2', 'api'], env: xtbenv)
test('GFN2-xTB: API (GBSA)', xtb_test, args: ['gfn2', 'gbsa'], env: xtbenv)
test('GFN2-xTB: API (GBSA+salt)', xtb_test, args: ['gfn2', 'salt'], env: xtbenv)
test('GFN2-xTB: API (PCEM)', xtb_test, args: ['gfn2', 'pcem'], env: xtbenv)
if cc.get_id() == 'intel'
   test('GFN2-xTB: C', executable('xtb_c_test', files('c_api_example.c'),
                                  include_directories: incdir,
                                  link_with: xtb_lib_static,
                                  dependencies: dependencies_exe),
                                  env: xtbenv)
endif

test('GFN1-xTB: SCC', xtb_test, args: ['gfn1', 'scc'], env: xtbenv)
test('GFN1-xTB: API', xtb_test, args: ['gfn1', 'api'], env: xtbenv)
test('GFN1-xTB: XB', xtb_test, args: ['gfn1', 'xb'], env: xtbenv)
test('GFN1-xTB: API (GBSA)', xtb_test, args: ['gfn1', 'gbsa'], env: xtbenv)
test('GFN1-xTB: API (PCEM)', xtb_test, args: ['gfn1', 'pcem'], env: xtbenv)

test('GFN0-xTB: SP', xtb_test, args: ['gfn0', 'sp'], env: xtbenv)
test('GFN0-xTB: API', xtb_test, args: ['gfn0', 'api'], env: xtbenv)
test('GFN0-xTB: SRB', xtb_test, args: ['gfn0', 'srb'], env: xtbenv)
test('GFN0-xTB: SP  (PBC)', xtb_test, args: ['peeq', 'sp'], env: xtbenv)
test('GFN0-xTB: API (PBC)', xtb_test, args: ['peeq', 'api'], env: xtbenv)
test('GFN0-xTB: SRB (PBC)', xtb_test, args: ['peeq', 'srb'], env: xtbenv)
