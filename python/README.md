# Python Wrapper for the Extended Tight Binding Program

This is the Python-side of the wrapper around the `xtb` library.

This wrapper provides access to the three Hamiltonians (GFN0-xTB, GFN1-xTB
and GFN2-xTB) implemented in the `xtb` program, making accessable
energies, gradients and related properties like dipole moments, partial charges
and bond orders.

It is not yet a full replacement for the `xtb` binary since it is mainly
intended to supplement its functionality.

## Usage

For the everyday usa we recommend to use `xtb` together with the
Atomic Simulation Model (ASE) like

```python
from ase.io import read
from xtb import GFN2
mol = read('coord')
mol.set_calculator(GFN2())
energy = mol.get_potential_energy()
forces = mol.get_forces()
charges = mol.get_charges()
```

The three calculators implement the corresponding Hamiltonians with the
following properties

|                    |`GFN0`|`GFN1`|`GFN2`|
|--------------------|------|------|------|
| energy             | x    | x    | x    |
| forces             | x    | x    | x    |
| stress tensor      | x    |      |      |
| dipole moment      |      | x    | x    |
| partial charges    |      | x    | x    |
| atomic dipoles     |      |      | x    |
| atomic quadrupoles |      |      | x    |
| bond order         |      | x    | x    |

The atomic dipoles and quadrupoles as well as the bond orders are currently
not accessable from the atoms object and must be taken directely from the
calculator object.

Additionally there is a solvation model calculator based on a
generalized Born model with a solvent accessable surface area, `GBSA`,
which can calculate Born radii and the surface area by itself
and together with an extended tight-binding calculator also a
solvation free energy. Use it as

```python
from ase.io import read
from xtb import GFN2
from xtb.solvation import GBSA
mol = read('coord')
mol.set_calculator(GBSA(solvent='ch2cl2', calc=GFN2()))
gsolv = mol.get_potential_energy()
```

which evaluates the GBSA model and performs two GFN2-xTB calculations.

## Technical Details

We provide a multilayered API for `xtb` starting from the Fortran
side with a set of standalone calculator subroutines in `xtb/calculator.f90`,
these calculators are exposed by using the `iso_c_binding` module as C-API
(see `xtb/c_api.f90` and `include/xtb.h`).
All non-Fortran compatible languages (that makes all languages except for Fortran),
should use this C-API to interface with `xtb`.

To use the C-API in Python we decided to define the interface on the Python
side rather than on the Fortran side using the `ctypes` module.
This interface (or Python-API) hides some of the implementation dependent details
from the user and provides access to all main calculators implemented in `xtb`.

Using the Python-API I wrapped everything up for the Atomic Simulation Environment
(ASE) by implementing a Calculator that handles the conversion between the
Atoms objects and the bundle of ndarrays required for Python-API.

As a Python user I recommend using the ASE Calculators gives access to a
feature-rich environment for computational chemistry.
For more information on ASE I refer to their [detailed documentation](https://wiki.fysik.dtu.dk/ase/).

I you want to interface your Python program with `xtb` and are afraid to add
a such heavy module as ASE as your dependency, you can interface directely
to the Python-API which lives in `xtb.interface` and does only require
`numpy` and `ctypes`.

Coming from every other C-compatible language I would recommend to wrap
the C-API in a similar way like I did for Python.
