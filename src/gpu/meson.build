# This file is part of xtb.
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

# CUDA acceleration, enabled with -Dgpu=true and switched on at runtime by
# --gpu.  With -Dgpu=false a stub is compiled in whose entry points all report
# failure, so the binary neither links against nor needs the CUDA runtime.

srcs += files('gpu_ctrl.f90')

gpu_objects = []

if get_option('gpu')
  cuda_prefix = get_option('gpu_prefix')
  cuda_inc = include_directories(cuda_prefix / 'include')
  # lib64 is a system CUDA toolkit, lib is a conda installed one
  cuda_libdir = [cuda_prefix / 'lib64', cuda_prefix / 'lib',
                 cuda_prefix / 'targets/x86_64-linux/lib']

  foreach cudalib : ['cudart', 'cusolver', 'cublas']
    lib_deps += cc.find_library(cudalib, dirs: cuda_libdir, required: true)
  endforeach

  # The backend is plain C against the cuSOLVER API, so no nvcc is needed.
  # omp_dep matters: the FP64 -> FP32 staging cast touches 2 x n^2 doubles per
  # SCC iteration and its `#pragma omp parallel for` is silently ignored
  # without it.
  gpu_objects += static_library(
    'xtb_gpu_lin',
    sources: files('xtb_gpu_lin.c'),
    include_directories: cuda_inc,
    dependencies: omp_dep,
    pic: true,
  ).extract_all_objects(recursive: false)
else
  srcs += files('xtb_gpu_stub.c')
endif
