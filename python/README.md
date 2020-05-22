# Python Wrapper for the Extended Tight Binding Program

This is a basic `ctypes` Python API around the `xtb` library.

This wrapper provides access to the three Hamiltonians (GFN0-xTB, GFN1-xTB
and GFN2-xTB and experimental access to the GFN-FF) implemented in the `xtb`
program, making accessible energies, gradients and related properties like
dipole moments, partial charges and bond orders.

It is not yet a full replacement for the `xtb` binary since it is mainly
intended to supplement its functionality.

## Usage

To use the C-API import the wrapper classes around the `xtb` objects:

```python
from xtb.interface import Calculator, Results, Param
import numpy as np
# atomic numbers range from 1 to 86
numbers = np.array(
    [6, 7, 6, 7, 6, 6, 6, 8, 7, 6, 8, 7, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
)
# cartesian coordinates in atomic units (Bohr)
positions = np.array(
    [
        [ 2.02799738646442,  0.09231312124713, -0.14310895950963],
        [ 4.75011007621000,  0.02373496014051, -0.14324124033844],
        [ 6.33434307654413,  2.07098865582721, -0.14235306905930],
        [ 8.72860718071825,  1.38002919517619, -0.14265542523943],
        [ 8.65318821103610, -1.19324866489847, -0.14231527453678],
        [ 6.23857175648671, -2.08353643730276, -0.14218299370797],
        [ 5.63266886875962, -4.69950321056008, -0.13940509630299],
        [ 3.44931709749015, -5.48092386085491, -0.14318454855466],
        [ 7.77508917214346, -6.24427872938674, -0.13107140408805],
        [10.30229550927022, -5.39739796609292, -0.13672168520430],
        [12.07410272485492, -6.91573621641911, -0.13666499342053],
        [10.70038521493902, -2.79078533715849, -0.14148379504141],
        [13.24597858727017, -1.76969072232377, -0.14218299370797],
        [ 7.40891694074004, -8.95905928176407, -0.11636933482904],
        [ 1.38702118184179,  2.05575746325296, -0.14178615122154],
        [ 1.34622199478497, -0.86356704498496,  1.55590600570783],
        [ 1.34624089204623, -0.86133716815647, -1.84340893849267],
        [ 5.65596919189118,  4.00172183859480, -0.14131371969009],
        [14.67430918222276, -3.26230980007732, -0.14344911021228],
        [13.50897177220290, -0.60815166181684,  1.54898960808727],
        [13.50780014200488, -0.60614855212345, -1.83214617078268],
        [ 5.41408424778406, -9.49239668625902, -0.11022772492007],
        [ 8.31919801555568, -9.74947502841788,  1.56539243085954],
        [ 8.31511620712388, -9.76854236502758, -1.79108242206824],
    ],
)

calc = Calculator(Param.GFN2xTB, numbers, positions)
results = calc.singlepoint()

print(results.get_energy())
```

The Python classes provide pretty much a one to one translation from the
C-API to Python, with the exception that the calculation environment and
molecular structure data are absorbed into the single point calculator
for convenience.

To update the `xtb` objects with new coordinates or lattice parameters use

```python
calc.update(positions)
results = calc.singlepoint(results)
```

Previously calculated results objects can be used to restart single point
calculations.
To change the molecular structure data the calculator has to be reconstructed,
which means the number of atoms and the atomic types are immutable for an
already existing calculator.

## Troubleshooting

Importing the shared library via `ctypes` is somewhat errorprone, in case the
default setup of the library breaks use

```python
from xtb.interface import _libxtb, XTBLibrary, load_library
_libxtb = XTBLibrary(libarary=load_library('path/to/libxtb')
```

to point the interface to the correct library.
