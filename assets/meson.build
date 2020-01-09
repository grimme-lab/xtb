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
