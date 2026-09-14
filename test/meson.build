# This file is part of xtb.
# SPDX-Identifier: LGPL-3.0-or-later
#
# xtb is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with xtb.  If not, see <https://www.gnu.org/licenses/>.

subdir('api')
subdir('unit')

if get_option('gpu_shim')
  python3 = find_program('python3', required: true)
  test(
    'gfn1-gfn2-gpu-sp-grad-opt',
    python3,
    args: [
      files('gpu/gfn12_gpu_gate.py'),
      '--cpu', xtb_exe,
      '--gpu', xtb_exe,
      '--fixture', files('gpu/water.xyz'),
    ],
    env: xtbenv,
    suite: 'gpu',
    timeout: 300,
  )
endif
