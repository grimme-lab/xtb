# Semiempirical Extended Tight-Binding Program Package

[![License](https://img.shields.io/github/license/grimme-lab/xtb)](https://github.com/grimme-lab/xtb/blob/main/COPYING)
[![Latest Version](https://img.shields.io/github/v/release/grimme-lab/xtb)](https://github.com/grimme-lab/xtb/releases/latest)
[![DOI](https://img.shields.io/badge/DOI-10.1002%2Fwcms.1493-blue)](https://doi.org/10.1002/wcms.1493)
[![Github Downloads All Releases](https://img.shields.io/github/downloads/grimme-lab/xtb/total)](https://github.com/grimme-lab/xtb/releases)

This is the official repository of the `xtb` program package developed by the Grimme group in Bonn.

<div align="center">
<img src="./assets/logo/xtb.svg" alt="Extended Tight Binding" width="220">
</div>

## Table of Contents
- [Installation](#installation)
  - [Meson](#meson)
  - [CMake](#cmake)
  - [Conda](#conda)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [Citations](#citations)
- [License](#license)

## Installation

[![Build Status](https://img.shields.io/github/actions/workflow/status/grimme-lab/xtb/fortran-build.yml?branch=main)](https://github.com/grimme-lab/xtb/actions)

Statically linked binaries are provided on the [latest release page](https://github.com/grimme-lab/xtb/releases/latest), and bleeding-edge static binaries (as well as the source distribution archive) built from the latest `main` branch are available on the [continuous release tag (`bleed`)](https://github.com/grimme-lab/xtb/releases/tag/bleed) for five platforms:
- Linux (`x86_64`, compiled with Intel oneAPI `ifx`/`icx` and MKL; `aarch64`, compiled with GCC and OpenBLAS)
- macOS (`arm64` and `x86_64`, compiled with GCC and OpenBLAS)
- Windows (`x86_64`, compiled with MinGW-w64 UCRT64 GCC and OpenBLAS)

The `xtb` program and library are also packaged on conda-forge for Linux (`x86_64`, `aarch64`, `ppc64le`) and macOS (`x86_64`, `arm64`).
For Homebrew users a custom tap is available at [grimme-lab/homebrew-qc](https://github.com/grimme-lab/homebrew-qc) providing prebuilt macOS/`x86_64` binaries; for macOS/`arm64`, binaries will be compiled on installation automatically.

This project supports two build systems, meson and CMake.
A short guide on the usage of each is given here, follow the linked instructions for more detailed information ([meson guide](./meson/README.adoc), [CMake guide](./cmake/README.adoc)).

**Compilers**:
  1. `ifx`, `icx`, `icpx`
  2. `gfortran`, `gcc`, `g++`
  3. `ifort` (<=2021.10.0), `icc` (<=2021.10.0)

Known issues with combinations:
  1. GCC-13 + MacOS + AArch64


### Meson

Using [meson](https://mesonbuild.com/) as build system requires you to install a fairly new version like 0.62 or newer (except 1.8.0 and 1.8.1).
To use the default backend of meson you have to install [ninja](https://ninja-build.org/) version 1.7 or newer.

```bash
export FC=ifx CC=icx
meson setup build --buildtype release --optimization 2
ninja -C build test
```

> [!IMPORTANT]
> Compilation with `meson` on macOS differs slightly from the protocol for Linux-based systems. Different BLAS libraries can lead to deviating results in rare cases – please stick to the following instructions.

<details>
  <summary><b>Setting up meson on macOS</b></summary>

#### Compiling with meson on macOS

1. **Use Homebrew for Package Management**: Install dependencies like `gcc`, `gfortran`, and `openblas` using Homebrew.
[Further information](https://brew.sh/) on how to setup `brew`.
Example:
   ```bash
   brew install gcc gfortran openblas
   ```
2. **meson setup call with appropriate environment variables**: Use the following adapted `meson setup` call to compile `xtb` on macOS (adjusting the OpenBLAS path and the installed Homebrew compiler version, e.g., `gfortran-14`/`gcc-14`, for your system).
   ```bash
   LDFLAGS="-L/opt/homebrew/opt/openblas/lib" CPPFLAGS="-I/opt/homebrew/opt/openblas/include" FC=gfortran-14 CC=gcc-14 meson setup _build --buildtype release -Dlapack=openblas
   ```
</details>
<br>

Make sure the testsuite is running without errors.

To install the `xtb` binaries to `/usr/local` use (might require `sudo`)

```bash
ninja -C build install
```

For more information on the build with meson see the instructions [here](./meson/README.adoc).


### CMake

The CMake build system requires both make and CMake to be installed, the latter has to be version 3.18 or newer.

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

The `xtb` documentation is hosted at [read-the-docs](https://xtb-docs.readthedocs.io/en/latest/).


## Contributing

Please read our [contributing guidelines](CONTRIBUTING.md)
before contributing to this project.


### Contributors

We are developing this program to make our research possible.
Many of the features that `xtb` has today have been added because there
was a dire need for them and we had many contributors who made these
features reality:

- P. Atkinson ([@patrickatkinson](https://github.com/patrickatkinson))
- T. Balduf ([@TyBalduf](https://github.com/TyBalduf))
- [C. Bannwarth](https://www.ipc.rwth-aachen.de/cms/IPC/Das-Institut/IPC-Arbeitsgruppen/~onnkh) ([@cbannwarth](https://github.com/cbannwarth))
- F. Bohle ([@fabothch](https://github.com/fabothch))
- [G. Brandenburg](http://www.gerit-brandenburg.de/) ([@gbrandenburg](https://github.com/gbrandenburg))
- [E. Caldeweyher](https://eikecaldeweyher.de/) ([@f3rmion](https://github.com/f3rmion))
- M. Checinski
- S. Dohm ([@thch-dohm](https://github.com/thch-dohm))
- S. Ehlert ([@awvwgk](https://github.com/awvwgk))
- S. Ehrlich
- M. Friede ([@marvinfriede](https://github.com/marvinfriede))
- T. Froitzheim ([@thfroitzheim](https://github.com/thfroitzheim))
- I. Gerasimov ([@foxtran](https://github.com/foxtran))
- J. Gorges ([@gorges97](https://github.com/gorges97))
- [S. Grimme](https://www.chemie.uni-bonn.de/pctc/mulliken-center/grimme/) ([@stefangrimme](https://github.com/stefangrimme))
- C. Hölzer ([@hoelzerC](https://github.com/hoelzerc))
- A. Katbashev ([@Albkat](https://github.com/albkat))
- J. Kohn ([@demonic-daisy](https://github.com/demonic-daisy))
- J. Koopman ([@JayTheDog](https://github.com/jaythedog))
- C. Lavigne ([@clavigne](https://github.com/clavigne))
- S. Lehtola ([@susilehtola](https://github.com/susilehtola))
- F. März
- M. Müller ([@marcelmbn](https://github.com/marcelmbn))
- F. Musil ([@felixmusil](https://github.com/felixmusil))
- H. Neugebauer ([@haneug](https://github.com/haneug))
- J. Pisarek
- C. Plett ([@cplett](https://github.com/cplett))
- P. Pracht ([@pprcht](https://github.com/pprcht))
- F. Pultar ([@pultar](https://github.com/pultar))
- T. Rose ([@Thomas3R](https://github.com/thomas3r))
- J. Seibert ([@liljay42](https://github.com/liljay42))
- L. M. Seidler ([@lmseidler](https://github.com/lmseidler))
- C. E. Selzer ([@chselz](https://github.com/chselz))
- P. Shushkov
- S. Spicher ([@sespic](https://github.com/sespic))
- M. Stahn ([@MtoLStoN](https://github.com/mtolston))
- M. Steiner ([@steinmig](https://github.com/steinmig))
- T. Strunk ([@timostrunk](https://github.com/timostrunk))
- J. Stückrath ([@jbstueckrath](https://github.com/jbstueckrath))
- J. Unsleber ([@nabbelbabbel](https://github.com/nabbelbabbel))

Contributors are listed in alphabetical order.
Some contributions predate the GitHub release of this project and are not visible in the repository commit history.
For the contributor data from the commit history since then look [here](https://github.com/grimme-lab/xtb/graphs/contributors).


## Citations

General Reference to `xtb` and the implemented GFN methods:
- C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme
  *WIREs Comput. Mol. Sci.*, **2021**, 11, e1493.
  DOI: [10.1002/wcms.1493](https://doi.org/10.1002/wcms.1493)

for the `xtb` software:
- A. Katbashev, M. Stahn, T. Rose, V. Alizadeh, M. Friede, C. Plett, P. Steinbach, S. Ehlert, *J. Phys. Chem. A*, **2025**, 129, 2667–2682.
  DOI: [10.1021/acs.jpca.4c08263](https://doi.org/10.1021/acs.jpca.4c08263)

for GFN-xTB:
- S. Grimme, C. Bannwarth, P. Shushkov, *J. Chem. Theory Comput.*, **2017**, 13, 1989-2009.
  DOI: [10.1021/acs.jctc.7b00118](https://dx.doi.org/10.1021/acs.jctc.7b00118)
- C. Bannwarth, S. Ehlert and S. Grimme., *J. Chem. Theory Comput.*, **2019**, 15, 1652-1671.
  DOI: [10.1021/acs.jctc.8b01176](https://dx.doi.org/10.1021/acs.jctc.8b01176)
- P. Pracht, E. Caldeweyher, S. Ehlert, S. Grimme, *ChemRxiv*, **2019**, preprint.
  DOI: [10.26434/chemrxiv.8326202.v1](https://dx.doi.org/10.26434/chemrxiv.8326202.v1)

for spGFNn-xTB (spin-polarization):
- H. Neugebauer, B. Bädorf, S. Ehlert, A. Hansen, S. Grimme, *J. Comput. Chem.*, **2023**, 44, 2120–2129.
  DOI: [10.1002/jcc.27185](https://doi.org/10.1002/jcc.27185)

for GFN-FF:
- S. Spicher and S. Grimme, *Angew. Chem. Int. Ed.*, **2020**, 59, 15665–15673
  DOI: [10.1002/anie.202004239](https://doi.org/10.1002/anie.202004239)
- T. Rose, M. Bursch, J.-M. Mewes, S. Grimme, *Inorg. Chem.*, **2024**, 63, 19364–19374.
  DOI: [10.1021/acs.inorgchem.4c03215](https://doi.org/10.1021/acs.inorgchem.4c03215)

for periodic GFN-FF (mcGFN-FF):
- S. Grimme and T. Rose, *Z. Naturforsch. B*, **2024**, 79, 191–200.
  DOI: [10.1515/znb-2023-0088](https://doi.org/10.1515/znb-2023-0088)

for PTB:
- S. Grimme, M. Müller, A. Hansen, *J. Chem. Phys.* **2023**, 158, 124111.
  DOI: [10.1063/5.0137838](https://doi.org/10.1063/5.0137838)

for CEH:
- M. Müller, T. Froitzheim, A. Hansen, S. Grimme, *J. Phys. Chem. A*, **2024**, 128, 10723–10736.
  DOI: [10.1021/acs.jpca.4c06989](https://doi.org/10.1021/acs.jpca.4c06989)

for GBSA and ALPB implicit solvation:
- S. Ehlert, M. Stahn, S. Spicher, S. Grimme,
  *J. Chem. Theory Comput.*, **2021**, 17, 4250-4261
  DOI: [10.1021/acs.jctc.1c00471](https://doi.org/10.1021/acs.jctc.1c00471)

for CPCM-X implicit solvation:
- M. Stahn, S. Ehlert, S. Grimme, *J. Phys. Chem. A*, **2023**, 127, 7036–7043.
  DOI: [10.1021/acs.jpca.3c04382](https://doi.org/10.1021/acs.jpca.3c04382)

for ddX implicit solvation:
- M. Nottoli, M. F. Herbst, A. Mikhalev, A. Jha, F. Lipparini, B. Stamm,
  *WIREs Comput. Mol. Sci.*, **2024**, 14, e1726.
  DOI: [10.1002/wcms.1726](https://doi.org/10.1002/wcms.1726)

for DFT-D4:
- E. Caldeweyher, C. Bannwarth and S. Grimme, *J. Chem. Phys.*, **2017**, 147, 034112.
  DOI: [10.1063/1.4993215](https://dx.doi.org/10.1063/1.4993215)
- E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth and S. Grimme, *J. Chem. Phys.*,
  **2019**, 150, 154122. DOI: [10.1063/1.5090222](https://dx.doi.org/10.1063/1.5090222)
- E. Caldeweyher, J.-M. Mewes, S. Ehlert and S. Grimme, *Phys. Chem. Chem. Phys.*, **2020**, 22, 8499-8512. 
  DOI: [10.1039/D0CP00502A](https://dx.doi.org/10.1039/D0CP00502A)
- L. Wittmann, I. Gordiy, M. Friede, B. Helmich-Paris, S. Grimme, A. Hansen and M. Bursch, *Phys. Chem. Chem. Phys.*, **2024**, 26, 21379-21394.
  DOI: [10.1039/D4CP01514B](https://doi.org/10.1039/D4CP01514B)

for sTDA-xTB:
- S. Grimme and C. Bannwarth, *J. Chem. Phys.*, **2016**, 145, 054103.
  DOI: [10.1063/1.4959605](https://dx.doi.org/10.1063/1.4959605)

for FOD:
- S. Grimme and A. Hansen, *Angew. Chem. Int. Ed.*, **2015**, 54, 12308–12313.
  DOI: [10.1002/anie.201501887](https://doi.org/10.1002/anie.201501887)
- C. A. Bauer, A. Hansen, S. Grimme, *Chem. Eur. J.*, **2017**, 23, 6150–6164.
  DOI: [10.1002/chem.201604682](https://doi.org/10.1002/chem.201604682)

in the mass-spec context:
- V. Asgeirsson, C. Bauer and S. Grimme, *Chem. Sci.*, **2017**, 8, 4879.
  DOI: [10.1039/c7sc00601b](https://dx.doi.org/10.1039/c7sc00601b)
- J. Koopman and S. Grimme, *ACS Omega*, **2019**, 4, 12, 15120-15133.
  DOI: [10.1021/acsomega.9b02011](https://dx.doi.org/10.1021/acsomega.9b02011)
- J. Koopman and S. Grimme, *J. Am. Soc. Mass Spectrom.*, **2021**, 32, 7, 1735-1751.
  DOI: [10.1021/jasms.1c00098](https://dx.doi.org/10.1021/jasms.1c00098)
- J. Gorges and S. Grimme, *Phys. Chem. Chem. Phys.*, **2025**, 27, 6899-6911.
  DOI: [10.1039/D5CP00316D](https://doi.org/10.1039/D5CP00316D)

for metadynamics refer to:
- S. Grimme, *J. Chem. Theory Comput.*, **2019**, 15, 2847-2862.
  DOI: [10.1021/acs.jctc.9b00143](https://dx.doi.org/10.1021/acs.jctc.9b00143)
  
for SPH calculations refer to:
- S. Spicher and S. Grimme, *J. Chem. Theory Comput.*, **2021**, 17, 1701–1714.
  DOI: [10.1021/acs.jctc.0c01306](https://doi.org/10.1021/acs.jctc.0c01306) 

for ONIOM refer to:
- C. Plett, A. Katbashev, S. Ehlert, S. Grimme, M. Bursch, *Phys. Chem. Chem. Phys.*, **2023**, 25, 17860-17868.
  DOI: [10.1039/D3CP02178E](https://doi.org/10.1039/D3CP02178E)

for DIPRO refer to:
- J. T. Kohn, N. Gildemeister, S. Grimme, D. Fazzi, A. Hansen, *J. Chem. Phys.*, **2023**, 159, 144106.
  DOI: [10.1063/5.0167484](https://doi.org/10.1063/5.0167484)

for O1NumHess refer to:
- B. Wang, S. Luo, Z. Wang, W. Liu, *J. Chem. Theory Comput.*, **2025**, 21, 10893-10909.
  DOI: [10.1021/acs.jctc.5c01354](https://doi.org/10.1021/acs.jctc.5c01354)

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
