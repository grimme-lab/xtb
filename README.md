# Semiempirical Extended Tight-Binding Program Package

[![License](https://img.shields.io/github/license/grimme-lab/xtb)](https://github.com/grimme-lab/xtb/blob/master/COPYING)
[![Latest Version](https://img.shields.io/github/v/release/grimme-lab/xtb)](https://github.com/grimme-lab/xtb/releases/latest)
[![DOI](https://zenodo.org/badge/211856832.svg)](https://zenodo.org/badge/latestdoi/211856832)
[![Github Downloads All Releases](https://img.shields.io/github/downloads/grimme-lab/xtb/total)](https://github.com/grimme-lab/xtb/releases)
[![Gitter](https://badges.gitter.im/xtb-dev/community.svg)](https://gitter.im/xtb-dev/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

This is the offical repository of the `xtb` program package developed by the Grimme group in Bonn.

<div align="center">
<img src="./assets/logo/xtb.svg" alt="Extended Tight Binding" width="220">
</div>

## Installation

[![Build Status](https://img.shields.io/travis/com/grimme-lab/xtb?logo=linux&logoColor=white)](https://travis-ci.com/grimme-lab/xtb)
[![Build Status](https://img.shields.io/github/workflow/status/grimme-lab/xtb/CI?logo=apple&logoColor=white)](https://github.com/grimme-lab/xtb/actions)

Statically linked binaries (Intel Compiler 17.0.7) can be found at the [latest release page](https://github.com/grimme-lab/xtb/releases/latest).
There is also a version of the shared library, which requires the Math Kernel Library and additional Intel specific libraries to be installed.

`xtb` is routinely compiled with Intel Parallel Studio 17 on our clusters in Bonn, successful builds on OSX have been performed as well.
We have not tried to build `xtb` on Windows so far.
It is also possible to compile `xtb` with GCC (version 8), but we recommend to use binaries compiled with Intel.

This projects supports two build systems, meson and CMake.
A short guide on the usage of each is given here, follow the linked instructions for a more detailed guide ([meson guide](./meson/README.adoc), [CMake guide](./cmake/README.adoc)).

### Meson

Using [meson](https://mesonbuild.com/) as build system requires you to install a fairly new version like 0.51 or newer.
To use the default backend of meson you have to install [ninja](https://ninja-build.org/) version 1.7 or newer.

```bash
export FC=ifort CC=icc
meson setup build --buildtype release --optimization 2
ninja -C build test
```

Make sure the testsuite is running without errors.

To install the `xtb` binaries to `/usr/local` use (might require `sudo`)

```bash
ninja -C build_intel install
```

For more information on the build with meson see the instructions [here](./meson/README.adoc).

### CMake

The CMake build system requires both make and CMake to be installed, the latter has to be version 3.9 or newer.

Building `xtb` with CMake works with the following chain of commands:

```bash
mkdir build
pushd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
ctest
popd
```

To install the `xtb` binaries to `/usr/local` use (might require `sudo`)

```bash
make -C build_intel install
```

For more detailed information on the build with CMake see the instructions [here](./cmake/README.adoc).

### Conda

[![Conda Version](https://img.shields.io/conda/vn/conda-forge/xtb.svg)](https://anaconda.org/conda-forge/xtb)

Installing `xtb` from the `conda-forge` channel can be achieved by adding `conda-forge` to your channels with:

```
conda config --add channels conda-forge
```

Once the `conda-forge` channel has been enabled, `xtb` can be installed with:

```
conda install xtb
```

It is possible to list all of the versions of `xtb` available on your platform with:

```
conda search xtb --channel conda-forge
```

## Documentation

[![Documentation Status](https://readthedocs.org/projects/xtb-docs/badge/?version=latest)](https://xtb-docs.readthedocs.io/en/latest/?badge=latest)

The `xtb` documentation is hosted at [read-the-docs](https://xtb-docs.readthedocs.io/en/latest/contents.html).

## Contributing

See our [contributing guidelines](CONTRIBUTING.md).

## Citations

for GFN-xTB:
- S. Grimme, C. Bannwarth, P. Shushkov, *J. Chem. Theory Comput.*, **2017**, 13, 1989-2009.
  DOI: [10.1021/acs.jctc.7b00118](https://dx.doi.org/10.1021/acs.jctc.7b00118)
- C. Bannwarth, S. Ehlert and S. Grimme., *J. Chem. Theory Comput.*, **2019**, 15, 1652-1671.
  DOI: [10.1021/acs.jctc.8b01176](https://dx.doi.org/10.1021/acs.jctc.8b01176)
- P. Pracht, E. Caldeweyher, S. Ehlert, S. Grimme, *ChemRxiv*, **2019**, preprint.
  DOI: [10.26434/chemrxiv.8326202.v1](https://dx.doi.org/10.26434/chemrxiv.8326202.v1)

for DFT-D4:
- E. Caldeweyher, C. Bannwarth and S. Grimme, *J. Chem. Phys.*, **2017**, 147, 034112.
  DOI: [10.1063/1.4993215](https://dx.doi.org/10.1063/1.4993215)
- E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth and S. Grimme, *J. Chem. Phys.*,
  **2019**, 150, 154122. DOI: [10.1063/1.5090222](https://dx.doi.org/10.1063/1.5090222)
- E. Caldeweyher, J.-M. Mewes, S. Ehlert and S. Grimme, *Phys. Chem. Chem. Phys.*, **2020**, just accepted. 
  DOI: [10.1039/D0CP00502A](https://dx.doi.org/10.1039/D0CP00502A) 

for sTDA-xTB:
- S. Grimme and C. Bannwarth, *J. Chem. Phys.*, **2016**, 145, 054103.
  DOI: [10.1063/1.4959605](https://dx.doi.org/10.1063/1.4959605)

in the mass-spec context:
- V. Asgeirsson, C. Bauer and S. Grimme, *Chem. Sci.*, **2017**, 8, 4879.
  DOI: [10.1039/c7sc00601b](https://dx.doi.org/10.1039/c7sc00601b)
- J. Koopman and S. Grimme, *ACS Omega*, **2019**, 4, 12, 15120-15133.",&
  DOI: [10.1021/acsomega.9b02011](https://dx.doi.org/10.1021/acsomega.9b02011)

for metadynamics refer to:
- S. Grimme, *J. Chem. Theory Comput.*, **2019**, 155, 2847-2862.
  DOI: [10.1021/acs.jctc.9b00143](https://dx.doi.org/10.1021/acs.jctc.9b00143)

## License

`xtb` is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

`xtb` is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
GNU Lesser General Public License for more details.
