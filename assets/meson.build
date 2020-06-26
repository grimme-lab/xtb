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

# misc. files for xtb installation
setup_files = [
  'config_env.bash',
  'config_env.csh',
]
foreach file: setup_files
  install_data(configure_file(input: files('templates'/file),
               output: file, configuration: config))
endforeach

# generate a pkg-config file for xtb
pkg_config = '@0@.pc'.format(meson.project_name())
install_data(configure_file(input: files('templates'/pkg_config),
             output: pkg_config, configuration: config),
             install_dir: get_option('libdir')/'pkgconfig')

# generate a tcl environment module
# The convention is to name the file after the version and put it in a directory
# with the program name, but we will name the intermediate tcl file differently.
# On installation the file will be copied to the data_dir (share/xtb) and renamed
# to the version number.
env_module = '@0@.tcl'.format(meson.project_name())
install_data(
  configure_file(
    input: files('templates/env-module.tcl'),
    output: env_module,
    configuration: config,
  ),
  rename: meson.project_version(),
  install_dir: get_option('datadir')/'modules'/'modulefiles'/meson.project_name()
)
