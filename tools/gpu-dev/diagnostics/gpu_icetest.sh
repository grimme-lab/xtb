#!/usr/bin/env bash
export NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
export PATH="$NVHPC/compilers/bin:/usr/local/bin:/usr/bin:/bin"
cd /mnt/e/Prasanna/xTB/xtb/build-gpu
SRC=../subprojects/toml-f/src/tomlf/ser.f90
INC="-Isubprojects/toml-f/libtoml-f.a.p -Isubprojects/toml-f -I../subprojects/toml-f -module subprojects/toml-f/libtoml-f.a.p"
for opt in "-O0" "-O2" "-O1 -Mnollvm" "-O0 -Mnollvm"; do
  echo "=== trying: $opt ==="
  if nvfortran $INC -Mbackslash -Mallocatable=03 $opt -o /tmp/ser_test.o -c "$SRC" 2>&1 | grep -qi "Internal compiler error"; then
    echo "   ICE"
  else
    echo "   OK (no ICE) with: $opt"
  fi
done
