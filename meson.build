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
project(
  'xtb',
  'fortran', 'c',
  version: '6.2.3',
  license: 'LGPL3',
  meson_version: '>=0.51',
  default_options: [
    'includedir=include/xtb',
  ],
)

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

configure_file(
  input: files('assets/templates/version.f90'),
  output: 'xtb_version.fh',
  configuration : config,
)

## ========================================== ##
## Configuration and dependencies
## ========================================== ##
dependencies = []
dependencies_exe = []
dependencies_sha = []

subdir('meson')

## ========================================== ##
## Sources
## ========================================== ##
srcs = []
prog = []
main = []

subdir('src')

srcs += 'symmetry/symmetry.f90'
srcs += 'symmetry/symmetry_i.c'

incdir = include_directories('include')

# Build target

# create a static library from all sources
xtb_lib_static = static_library(
  meson.project_name(),
  srcs,
  include_directories: incdir,
  install: true,
  pic: true,
)

xtb_exe = executable(
  meson.project_name(),
  prog,
  dependencies: dependencies_exe,
  include_directories: incdir,
  link_with: xtb_lib_static,
  link_language: 'fortran',
  install: true,
)

xtb_lib_shared = shared_library(
  meson.project_name(),
  version: meson.project_version(),
  dependencies: dependencies_sha,
  include_directories: incdir,
  link_whole: xtb_lib_static,
  install: true
)

xtb_dep_static = [
  declare_dependency(
    link_whole: xtb_lib_static,
    ),
  dependencies_exe,
]

xtb_dep = [
  declare_dependency(
    link_with: xtb_lib_shared,
  ),
  dependencies_sha,
]

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

xtb_parameter_files = [
  '.param_gfn0.xtb',
  '.param_gfn2.xtb',
  '.param_gfn.xtb',
]

install_data(xtb_parameter_files)

if get_option('install_modules')
  xtb_objs = []
  foreach file: srcs
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

if get_option('python')
  py_xtb = meson.current_source_dir() / 'python' / 'xtb'

  has_numpy = run_command('python3', '-c', '"import numpy"')
  has_ase = run_command('python3', '-c', '"import ase"')
  pytest = find_program('pytest', required: false)
  if pytest.found() and has_numpy and has_ase
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
endif
