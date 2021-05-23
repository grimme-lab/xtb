# Semiempirical Extended Tight-Binding Program Package

[![License](https://img.shields.io/github/license/grimme-lab/xtb)](https://github.com/grimme-lab/xtb/blob/master/COPYING)
[![Latest Version](https://img.shields.io/github/v/release/grimme-lab/xtb)](https://github.com/grimme-lab/xtb/releases/latest)
[![DOI](https://img.shields.io/badge/DOI-10.1002%2Fwcms.1493-blue)](https://doi.org/10.1002/wcms.1493)
[![Github Downloads All Releases](https://img.shields.io/github/downloads/grimme-lab/xtb/total)](https://github.com/grimme-lab/xtb/releases)

This is the offical repository of the `xtb` program package developed by the Grimme group in Bonn.

<div align="center">
<img src="./assets/logo/xtb.svg" alt="Extended Tight Binding" width="220">
</div>


## Installation

[![Build Status](https://img.shields.io/github/workflow/status/grimme-lab/xtb/CI)](https://github.com/grimme-lab/xtb/actions)

Statically linked binaries (Intel Compiler 18.0.2) can be found at the [latest release page](https://github.com/grimme-lab/xtb/releases/latest).
There is also a version of the shared library, which requires the Math Kernel Library and additional Intel specific libraries to be installed.
Bleeding edge releases of the latest source from this repository are available on the [continuous release tag](https://github.com/grimme-lab/xtb/releases/tag/bleed).

`xtb` is routinely compiled with Intel Parallel Studio 18 and newer on our clusters in Bonn.
It is also possible to compile `xtb` with GCC (version 7.5 or newer), but we recommend to use binaries compiled with Intel.
Successful builds on OSX with GCC via homebrew have been performed as well.
The NVHPC compilers (version 20.9, 20.11) can be used to compile for CPU and GPU.
To compile on Windows we recommend the MinGW toolchain (installable with [MSYS2](https://msys2.org)) or, if a POSIX environment is preferred, [Cygwin](https://cygwin.com).
It has been reported that `xtb` can be compiled using Intel Fortran on Windows as well, but official support is not yet established.

This projects supports two build systems, meson and CMake.
A short guide on the usage of each is given here, follow the linked instructions for a more detailed information ([meson guide](./meson/README.adoc), [CMake guide](./cmake/README.adoc)).


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
ninja -C build install
```

For more information on the build with meson see the instructions [here](./meson/README.adoc).


### CMake

The CMake build system requires both make and CMake to be installed, the latter has to be version 3.9 or newer.

Building `xtb` with CMake works with the following chain of commands:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
make -C build
make -C build test
```

To install the `xtb` binaries to `/usr/local` use (might require `sudo`)

```bash
make -C build install
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

Please read our [contributing guidelines](CONTRIBUTING.md)
before contributing to this project.


### Contributors

We are developing this program to make our research possible.
Many of the features that `xtb` has today have been added because there
was a dire need for them and we had many contributors who made these
features reality:

- P. Atkinson ([@patrickatkinson](https://github.com/patrickatkinson))
- C. Bannwarth ([@cbannwarth](https://github.com/cbannwarth))
- F. Bohle ([@fabothch](https://github.com/fabothch))
- [G. Brandenburg](http://www.gerit-brandenburg.de/) ([@gbrandenburg](https://github.com/gbrandenburg))
- [E. Caldeweyher](https://eikecaldeweyher.de/) ([@f3rmion](https://github.com/f3rmion))
- M. Checinski
- S. Dohm ([@thch-dohm](https://github.com/thch-dohm))
- S. Ehlert ([@awvwgk](https://github.com/awvwgk))
- S. Ehrlich
- I. Gerasimov ([@FulgurIgor](https://github.com/fulgurigor))
- [S. Grimme](https://www.chemie.uni-bonn.de/pctc/mulliken-center/grimme/) ([@stefangrimme](https://github.com/stefangrimme))
- J. Koopman ([@JayTheDog](https://github.com/jaythedog))
- C. Lavinge ([@clavigne](https://github.com/clavigne))
- S. Lehtola ([@susilehtola](https://github.com/susilehtola))
- F. März
- M. Müller ([@marcelmbn](https://github.com/marcelmbn))
- F. Musil ([@felixmusil](https://github.com/felixmusil))
- H. Neugebauer ([@haneug](https://github.com/haneug))
- J. Pisarek
- C. Plett ([@cplett](https://github.com/cplett))
- P. Pracht ([@pprcht](https://github.com/pprcht))
- J. Seibert ([@liljay42](https://github.com/liljay42))
- P. Shushkov
- S. Spicher ([@sespic](https://github.com/sespic))
- M. Stahn ([@MtoLStoN](https://github.com/mtolston))
- M. Steiner ([@steinmig](https://github.com/steinmig))
- T. Strunk ([@timostrunk](https://github.com/timostrunk))
- J. Stückrath ([@jbstueckrath](https://github.com/jbstueckrath))
- T. Rose ([@Thomas3R](https://github.com/thomas3r))
- J. Unsleber ([@nabbelbabbel](https://github.com/nabbelbabbel))

Contributors are listed in alphabetical order.
Some contributions predate the GitHub release of this project and are not visible in the repository commit history.
For the contributor data from the commit history since then look [here](https://github.com/grimme-lab/xtb/graphs/contributors).


## Citations

General Reference to `xtb` and the implemented GFN methods:
- C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme
  *WIREs Comput. Mol. Sci.*, **2020**, 11, e01493.
  DOI: [10.1002/wcms.1493](https://doi.org/10.1002/wcms.1493)

for GFN-xTB:
- S. Grimme, C. Bannwarth, P. Shushkov, *J. Chem. Theory Comput.*, **2017**, 13, 1989-2009.
  DOI: [10.1021/acs.jctc.7b00118](https://dx.doi.org/10.1021/acs.jctc.7b00118)
- C. Bannwarth, S. Ehlert and S. Grimme., *J. Chem. Theory Comput.*, **2019**, 15, 1652-1671.
  DOI: [10.1021/acs.jctc.8b01176](https://dx.doi.org/10.1021/acs.jctc.8b01176)
- P. Pracht, E. Caldeweyher, S. Ehlert, S. Grimme, *ChemRxiv*, **2019**, preprint.
  DOI: [10.26434/chemrxiv.8326202.v1](https://dx.doi.org/10.26434/chemrxiv.8326202.v1)

for GFN-FF:
- S. Spicher and S. Grimme, *Angew. Chem. Int. Ed.*, **2020**, 59, 15665–15673
  DOI: [10.1002/anie.202004239](https://doi.org/10.1002/anie.202004239)

for DFT-D4:
- E. Caldeweyher, C. Bannwarth and S. Grimme, *J. Chem. Phys.*, **2017**, 147, 034112.
  DOI: [10.1063/1.4993215](https://dx.doi.org/10.1063/1.4993215)
- E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth and S. Grimme, *J. Chem. Phys.*,
  **2019**, 150, 154122. DOI: [10.1063/1.5090222](https://dx.doi.org/10.1063/1.5090222)
- E. Caldeweyher, J.-M. Mewes, S. Ehlert and S. Grimme, *Phys. Chem. Chem. Phys.*, **2020**, 22, 8499-8512. 
  DOI: [10.1039/D0CP00502A](https://dx.doi.org/10.1039/D0CP00502A) 

for sTDA-xTB:
- S. Grimme and C. Bannwarth, *J. Chem. Phys.*, **2016**, 145, 054103.
  DOI: [10.1063/1.4959605](https://dx.doi.org/10.1063/1.4959605)

in the mass-spec context:
- V. Asgeirsson, C. Bauer and S. Grimme, *Chem. Sci.*, **2017**, 8, 4879.
  DOI: [10.1039/c7sc00601b](https://dx.doi.org/10.1039/c7sc00601b)
- J. Koopman and S. Grimme, *ACS Omega*, **2019**, 4, 12, 15120-15133.
  DOI: [10.1021/acsomega.9b02011](https://dx.doi.org/10.1021/acsomega.9b02011)

for metadynamics refer to:
- S. Grimme, *J. Chem. Theory Comput.*, **2019**, 155, 2847-2862.
  DOI: [10.1021/acs.jctc.9b00143](https://dx.doi.org/10.1021/acs.jctc.9b00143)
  
for SPH calculations refer to:
- S. Spicher and S. Grimme, *J. Chem. Theory Comput.*, **2021**, 17, 1701–1714.
  DOI: [10.1021/acs.jctc.0c01306](https://doi.org/10.1021/acs.jctc.0c01306)  

All references are available in [bibtex format](./assets/references.bib).


## License

`xtb` is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

`xtb` is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
GNU Lesser General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in `xtb` by you, as defined in the
GNU Lesser General Public license, shall be licensed as above, without any
additional terms or conditions.
