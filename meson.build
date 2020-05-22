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
  version: '6.3.0',
  license: 'LGPL3',
  meson_version: '>=0.51',
  default_options: [
    'includedir=include/xtb',
  ],
)

commit = get_option('build_name')
git = find_program('git', required: false)
if git.found()
  git_commit = run_command(git, 'show', '-s', '--format=%h')
  if git_commit.returncode() == 0
    commit = git_commit.stdout().strip()
  endif
endif

# build a configuration data containing all the important data to propagate
# it to the automatically generated files
config = configuration_data({
  'name': meson.project_name(),
  'description': 'Semiempirical Extended Tight-Binding Program Package',
  'version': meson.project_version(),
  'commit': commit,
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

xtb_inc = meson.current_source_dir() / 'include'
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

xtb_dep_static = declare_dependency(
  link_with: xtb_lib_static,
  include_directories: [incdir, xtb_lib_static.private_dir_include()],
  dependencies: dependencies_exe,
)

xtb_dep = declare_dependency(
  link_with: xtb_lib_shared,
  include_directories: incdir,
  dependencies: dependencies_sha,
)

xtb_header = files('include/xtb.h')

## ========================================== ##
## INSTALL
## ========================================== ##
install_headers(xtb_header)

asciidoc = find_program('a2x', required: false)
if asciidoc.found()
  run_command(asciidoc, 'man/man1/xtb.1.txt', '--format=manpage')
  run_command(asciidoc, 'man/man7/xcontrol.7.txt', '--format=manpage')
  install_man(['man/man1/xtb.1', 'man/man7/xcontrol.7'])
endif

xtb_parameter_files = [
  'param_gfn0-xtb.txt',
  'param_gfn2-xtb.txt',
  'param_gfn1-xtb.txt',
  'param_ipea-xtb.txt',
  '.param_gfnff.xtb',
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
  pymod = import('python')
  python = pymod.find_installation(
    'python3',
    required: false,
  )
  if python.found()
    has_numpy = run_command(python, '-c', '"import numpy"').returncode() == 0
    has_pytest = run_command(python, '-c', '"import pytest"').returncode() == 0
    if has_numpy and has_pytest
      test('pytest: xtb', python, args: ['-m', 'pytest', '--pyargs', 'xtb'], env: xtbenv)
    endif
  endif
endif
