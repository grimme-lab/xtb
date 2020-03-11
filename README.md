# Semiempirical Extended Tight-Binding Program Package

[![Build Status](https://img.shields.io/travis/com/grimme-lab/xtb?logo=linux&logoColor=white)](https://travis-ci.com/grimme-lab/xtb)
[![Build Status](https://img.shields.io/github/workflow/status/grimme-lab/xtb/CI?logo=apple&logoColor=white)](https://github.com/grimme-lab/xtb/actions)
[![License](https://img.shields.io/github/license/grimme-lab/xtb)](https://github.com/grimme-lab/xtb/blob/master/COPYING)
[![Latest Version](https://img.shields.io/github/v/release/grimme-lab/xtb)](https://github.com/grimme-lab/xtb/releases/latest)
[![Github Downloads All Releases](https://img.shields.io/github/downloads/grimme-lab/xtb/total)](https://github.com/grimme-lab/xtb/releases)

This is the offical repository of the `xtb` program package developed by the Grimme group in Bonn.

<img src="./assets/logo/xtb.svg" alt="Extended Tight Binding" width="120">

## Installation

Statically linked binaries (Intel Compiler 17.0.7) can be found at the [latest release page](https://github.com/grimme-lab/xtb/releases/latest).
There is also a version of the shared library, which requires the Math Kernel Library and additional Intel specific libraries to be installed.

To compile `xtb` from source install Intel Parallel Studio 17 or later, or GCC version 8 or later.

There are two ways of building xtb.

### Meson
Using [`meson`](https://mesonbuild.com/) as build system requires you to install a fairly new version like 0.49 or newer.
To use the default backend of `meson` you have to install [`ninja`](https://ninja-build.org/) version 1.5 or newer.

```bash
export FC=ifort CC=icc
meson setup build_intel --buildtype release --optimization=2
ninja -C build_intel test
```

Make sure the testsuite is running without errors.
`xtb` is routinely compiled with Intel Parallel Studio 17 on our clusters in Bonn,
but we have not tried to compile it on either OSX or Windows so far.
It is also possible to compile `xtb` with GCC (version 8), but we recommend to use binaries compiled with Intel.

To install the `xtb` binaries to `/usr/local` use (might require `sudo`)

```bash
ninja -C build_intel install
```

For a local installation (or if you want to pack a release), modify the
configuration by using

```bash
meson configure build_intel --prefix=/
DESTDIR=$HOME/.local ninja -C build_intel install
```

The build system will generate configuration files in `$DESTDIR/share/xtb` to
be sourced in your `.bashrc` or `.cshrc` which will make `xtb` and the
parameter files available. Also a pkg-config file is generated to make
`xtb` available in other projects as dependency.

### CMake

The CMake build system requires both `make` and `cmake` to be installed, the
latter has to be version 3.9 or newer.

Building xtb with CMake works with the following chain of commands:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local ..
make
make install
```

The compiled tests can be executed using the following commands:

```bash
cd build && XTBPATH=$PWD/.. OMP_NUM_THREADS=1 make test
```

(This also works when not installing the program, but just building it.)

## Documentation

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

for sTDA-xTB:
- S. Grimme and C. Bannwarth, *J. Chem. Phys.*, **2016**, 145, 054103.
  DOI: [10.1063/1.4959605](https://dx.doi.org/10.1063/1.4959605)

in the mass-spec context:
- V. Asgeirsson, C. Bauer and S. Grimme, *Chem. Sci.*, **2017**, 8, 4879.
  DOI: [10.1039/c7sc00601b](https://dx.doi.org/10.1039/c7sc00601b)

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
