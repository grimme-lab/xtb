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
  version: '6.4.1',
  license: 'LGPL3',
  meson_version: '>=0.51',
  default_options: [
    'includedir=include/xtb',
    'buildtype=debugoptimized',
    'default_library=both',
    'optimization=2',
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
exe_deps = []
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

# create a static library from all sources,
# this is the common base of all targets created for xtb
xtb_common = static_library(
  '_@0@'.format(meson.project_name()),
  sources: srcs,
  include_directories: incdir,
  dependencies: exe_deps,
  install: false,
  pic: true,
)

xtb_dep_static = declare_dependency(
  link_with: xtb_common,
  include_directories: [incdir, xtb_common.private_dir_include()],
  dependencies: exe_deps,
  version: meson.project_version(),
)

if fc.get_id() == 'pgi' or fc.get_id() == 'nvidia_hpc'
   warning('Shared library extension cannot be compiled with PGI')
else
   xtb_lib = library(
     meson.project_name(),
     version: meson.project_version(),
     dependencies: lib_deps,
     include_directories: incdir,
     link_whole: xtb_common,
     install: install,
   )

   xtb_dep = declare_dependency(
     link_with: xtb_lib,
     include_directories: incdir,
     dependencies: lib_deps,
     version: meson.project_version(),
   )
endif

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
subdir('TESTSUITE')
