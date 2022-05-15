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
  version: '6.5.1',
  license: 'LGPL3-3.0-or-later',
  meson_version: '>=0.57.2',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
    'optimization=2',
    'c_std=c11',
  ],
)

install = not meson.is_subproject()
install_modules = install and get_option('install_modules')

commit = get_option('build_name')
git = find_program('git', required: false)
if git.found()
  git_commit = run_command(git, 'show', '-s', '--format=%h')
  if git_commit.returncode() == 0
    commit = git_commit.stdout().strip()
  endif
endif

# We rely on the existence of Python since meson will need a Python interpreter to run,
# this way we have a cross-platform way to obtain relevant data for our build.
pymod = import('python')
python = pymod.find_installation('python3', required: false)
if not python.found()
  python = find_program('python3', 'python')
endif

# build a configuration data containing all the important data to propagate
# it to the automatically generated files
config = configuration_data({
  'name': meson.project_name(),
  'description': 'Semiempirical Extended Tight-Binding Program Package',
  'version': meson.project_version(),
  'commit': commit,
  'date': run_command(python, '-c', 'import datetime; print(datetime.date.today().strftime("%Y-%m-%d"))').stdout().strip(),
  'author': run_command(python, '-c', 'import getpass; print(getpass.getuser())').stdout().strip(),
  'origin': run_command(python, '-c', 'import socket; print(socket.gethostname())').stdout().strip(),
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
lib_deps = []

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

if get_option('nvtx')
  srcs += 'src/nvtx.f90'
endif

xtb_inc = meson.current_source_dir() / 'include'
incdir = include_directories('include')

# Build target
xtb_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  include_directories: incdir,
  install: install,
  fortran_args: [
    '-DWITH_TBLITE=@0@'.format(tblite_dep.found() ? 1 : 0),
  ],
  link_language: 'fortran',
)

xtb_dep = declare_dependency(
  link_with: xtb_lib,
  include_directories: [incdir, xtb_lib.private_dir_include()],
  dependencies: lib_deps,
  version: meson.project_version(),
)

xtb_common = get_option('default_library') == 'both' ? xtb_lib.get_static_lib() : xtb_lib

xtb_dep_static = declare_dependency(
  link_with: xtb_common,
  include_directories: [incdir, xtb_lib.private_dir_include()],
  dependencies: lib_deps,
  version: meson.project_version(),
)

xtb_exe = executable(
  meson.project_name(),
  sources: prog,
  dependencies: xtb_dep_static,
  link_language: 'fortran',
  install: install,
)

xtb_header = files('include/xtb.h')

## ========================================== ##
## INSTALL
## ========================================== ##
if install
  install_headers(xtb_header)
endif

asciidoc = find_program('asciidoctor', required: false)
if asciidoc.found()
  xtb_man = []
  xtb_man += configure_file(
    command: [asciidoc, '-b', 'manpage', '@INPUT@', '-o', '@OUTPUT@'],
    input: files('man/xtb.1.adoc'),
    output: '@BASENAME@',
  )
  xtb_man += configure_file(
    command: [asciidoc, '-b', 'manpage', '@INPUT@', '-o', '@OUTPUT@'],
    input: files('man/xcontrol.7.adoc'),
    output: '@BASENAME@',
  )
  if install
    install_man(xtb_man)
  endif
endif

xtb_parameter_files = [
  'param_gfn0-xtb.txt',
  'param_gfn2-xtb.txt',
  'param_gfn1-xtb.txt',
  'param_gfn1-si-xtb.txt',
  'param_ipea-xtb.txt',
  '.param_gfnff.xtb',
]

if install
  install_data(xtb_parameter_files)
endif

if install_modules
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
if install
  subdir('assets')
endif

## ========================================== ##
## TESTSUITE
## ========================================== ##
# make sure the correct library is loaded
xtbenv = environment()
xtbenv.prepend('LD_LIBRARY_PATH', meson.current_build_dir())
xtbenv.prepend('PYTHONPATH', meson.current_source_dir() / 'python')
xtbenv.set('XTBPATH', meson.current_source_dir())

# all tests are defined in a separate meson.build
subdir('test')
