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
"""Wrapper around the C-API of the xtb shared library."""

from ctypes import (
    c_void_p,
    POINTER,
    c_int,
    c_double,
    c_bool,
    c_char_p,
    cdll,
    CDLL,
    byref,
    pointer,
)
from distutils.sysconfig import get_config_vars
from enum import Enum, auto
from typing import List, Optional
import sys
import os.path as op
import numpy as np


_libxtb = None


def get_shared_lib_extension():
    """Try to figure out which extension a shared library should have on the
    given platform. This code is borrowed from numpy and slightly modified."""
    if sys.platform.startswith("linux") or sys.platform.startswith("gnukfreebsd"):
        return ".so"
    if sys.platform.startswith("darwin"):
        return ".dylib"
    if sys.platform.startswith("win"):
        return ".dll"
    confvars = get_config_vars()
    # SO is deprecated in 3.3.1, use EXT_SUFFIX instead
    so_ext = confvars.get("EXT_SUFFIX", None)
    if so_ext is None:
        so_ext = confvars.get("SO", "")
    if "SOABI" in confvars:
        # Does nothing unless SOABI config var exists
        so_ext = so_ext.replace("." + confvars.get("SOABI"), "", 1)
    return so_ext


def load_library(libname: str) -> CDLL:
    """load library cross-platform compatible."""

    if not op.splitext(libname)[1]:
        # Try to load library with platform-specific name
        so_ext = get_shared_lib_extension()
        libname_ext = libname + so_ext
    else:
        libname_ext = libname

    return cdll.LoadLibrary(libname_ext)


class Param(Enum):
    """Possible parametrisations for the Calculator class"""

    GFN2xTB = auto()
    GFN1xTB = auto()
    GFN0xTB = auto()
    GFNFF = auto()


VEnvironment = c_void_p
VMolecule = c_void_p
VCalculator = c_void_p
VResults = c_void_p

_XTB_API_6_3 = {
    "xtb_newEnvironment": (VEnvironment, []),
    #"xtb_delEnvironment": (None, [POINTER(VEnvironment)]),
    "xtb_checkEnvironment": (c_int, [VEnvironment]),
    "xtb_showEnvironment": (None, [VEnvironment, c_char_p]),
    "xtb_setOutput": (None, [VEnvironment, c_char_p]),
    "xtb_releaseOutput": (None, [VEnvironment]),
    "xtb_setVerbosity": (None, [VEnvironment, c_int]),
    "xtb_newMolecule": (
        VMolecule,
        [
            VEnvironment,
            POINTER(c_int),
            POINTER(c_int),
            POINTER(c_double),
            POINTER(c_double),
            POINTER(c_int),
            POINTER(c_double),
            POINTER(c_bool),
        ],
    ),
    #"xtb_delMolecule": (None, [POINTER(VMolecule)]),
    "xtb_updateMolecule": (
        None,
        [VEnvironment, VMolecule, POINTER(c_double), POINTER(c_double),],
    ),
    "xtb_newCalculator": (VCalculator, []),
    #"xtb_delCalculator": (None, [POINTER(VCalculator)]),
    "xtb_loadGFN0xTB": (None, [VEnvironment, VMolecule, VCalculator, c_char_p]),
    "xtb_loadGFN1xTB": (None, [VEnvironment, VMolecule, VCalculator, c_char_p]),
    "xtb_loadGFN2xTB": (None, [VEnvironment, VMolecule, VCalculator, c_char_p]),
    "xtb_loadGFNFF": (None, [VEnvironment, VMolecule, VCalculator, c_char_p]),
    "xtb_setSolvent": (
        None,
        [
            VEnvironment,
            VCalculator,
            c_char_p,
            POINTER(c_int),
            POINTER(c_double),
            POINTER(c_int),
        ],
    ),
    "xtb_releaseSolvent": (None, [VEnvironment, VCalculator]),
    "xtb_setExternalCharges": (
        None,
        [
            VEnvironment,
            VCalculator,
            POINTER(c_int),
            POINTER(c_int),
            POINTER(c_double),
            POINTER(c_double),
        ],
    ),
    "xtb_releaseExternalCharges": (None, [VEnvironment, VCalculator]),
    "xtb_singlepoint": (None, [VEnvironment, VMolecule, VCalculator, VResults]),
    "xtb_newResults": (VResults, []),
    #"xtb_delResults": (None, [POINTER(VResults)]),
    "xtb_getEnergy": (None, [VEnvironment, VResults, POINTER(c_double)]),
    "xtb_getGradient": (None, [VEnvironment, VResults, POINTER(c_double)]),
    "xtb_getVirial": (None, [VEnvironment, VResults, POINTER(c_double)]),
    "xtb_getDipole": (None, [VEnvironment, VResults, POINTER(c_double)]),
    "xtb_getBondOrders": (None, [VEnvironment, VResults, POINTER(c_double)]),
}


class XTBLibrary:
    """Shared library instance"""

    def __init__(self, library: Optional[CDLL] = None):
        if library is not None:
            if isinstance(library, XTBLibrary):
                self._lib = library._lib
            else:
                self._lib = library
        else:
            self._lib = load_library("libxtb")

        for name, signature in _XTB_API_6_3.items():
            restype, argtypes = signature
            _set(self._lib, name, restype, argtypes)


class Environment:
    """Calculation environment"""

    _env = None

    def __init__(self, library: Optional[CDLL] = None):
        """Create new xtb calculation environment object"""
        global _libxtb
        if _libxtb is None:
            _libxtb = XTBLibrary(library)
        self._env = c_void_p(_libxtb._lib.xtb_newEnvironment())

    def __del__(self):
        """Delete a xtb calculation environment object"""
        global _libxtb
        if self._env is not None:
            _libxtb._lib.xtb_delEnvironment(byref(self._env))

    def check(self) -> int:
        """Check current status of calculation environment"""
        global _libxtb
        return _libxtb._lib.xtb_checkEnvironment(self._env)

    def show(self, message: str) -> None:
        """Show and empty error stack"""
        global _libxtb
        _message = message.encode()
        _libxtb._lib.xtb_showEnvironment(self._env, _message)

    def set_output(self, filename: str) -> None:
        """Bind output from this environment"""
        global _libxtb
        _filename = filename.encode()
        _libxtb._lib.xtb_setOutput(self._env, _filename)

    def release_output(self) -> None:
        """Release output unit from this environment"""
        global _libxtb
        _libxtb._lib.xtb_releaseOutput(self._env)

    def set_verbosity(self, verbosity: int) -> None:
        """Set verbosity of calculation output"""
        global _libxtb
        _libxtb._lib.xtb_setVerbosity(self._env, verbosity)


class Molecule(Environment):
    """Molecular structure data"""

    _mol = None

    def __init__(
        self,
        numbers: List[int],
        positions: List[float],
        charge: Optional[float] = None,
        uhf: Optional[int] = None,
        lattice: Optional[List[float]] = None,
        periodic: Optional[List[bool]] = None,
        library: Optional[CDLL] = None,
    ):
        """Create new molecular structure data"""
        global _libxtb
        Environment.__init__(self, library)
        if isinstance(positions, np.ndarray):
            if positions.size % 3 != 0:
                raise ValueError("Expected tripels of cartesian coordinates")
        else:
            if len(positions) % 3 != 0:
                raise ValueError("Expected tripels of cartesian coordinates")

        if isinstance(positions, np.ndarray):
            if 3 * len(numbers) != positions.size:
                raise ValueError("Dimension missmatch between numbers and postions")
        else:
            if 3 * len(numbers) != len(positions):
                raise ValueError("Dimension missmatch between numbers and postions")

        self._natoms = len(numbers)
        _numbers = np.array(numbers, dtype="i4")
        _positions = np.array(positions, dtype="float")

        if lattice is not None:
            if len(lattice) != 9:
                raise ValueError("Invalid lattice provided")
            _lattice = np.array(lattice, dtype="float")
        else:
            _lattice = None

        if periodic is not None:
            if len(periodic) != 3:
                raise ValueError("Invalid periodicity provided")
            _periodic = np.array(periodic, dtype="bool")
        else:
            _periodic = None

        self._mol = c_void_p(_libxtb._lib.xtb_newMolecule(
            self._env,
            c_int(self._natoms),
            _cast(POINTER(c_int), _numbers),
            _cast(POINTER(c_double), _positions),
            None if charge is None else c_double(charge),
            None if uhf is None else c_int(uhf),
            _cast(POINTER(c_double), _lattice),
            _cast(POINTER(c_bool), _periodic),
        ))

        if self.check() != 0:
            raise ValueError("Could not initialize molecular structure data")

    def __del__(self):
        """Delete molecular structure data"""
        global _libxtb
        Environment.__del__(self)
        if self._mol is not None:
            _libxtb._lib.xtb_delMolecule(byref(self._mol))

    def __len__(self):
        return self._natoms

    def update(
        self, positions: List[float], lattice: Optional[List[float]] = None,
    ):
        """Update coordinates and lattice parameters"""
        global _libxtb

        if 3 * len(self) != len(positions):
            raise ValueError("Dimension missmatch for postions")
        _positions = np.array(positions, dtype="float")

        if lattice is not None:
            if len(lattice) != 9:
                raise ValueError("Invalid lattice provided")
            _lattice = np.array(lattice, dtype="float")
        else:
            _lattice = None

        _libxtb._lib.xtb_updateMolecule(
            self._env,
            self._mol,
            _cast(POINTER(c_double), _positions),
            _cast(POINTER(c_double), _lattice),
        )

        if self.check() != 0:
            raise ValueError("Could not update molecular structure data")


class Results(Environment):
    """Calculation results"""

    _res = None

    def __init__(self, mol: Molecule, library: Optional[CDLL] = None):
        """Create new singlepoint results object"""
        global _libxtb
        Environment.__init__(self, library)
        self._res = c_void_p(_libxtb._lib.xtb_newResults())
        self._natoms = len(mol)

    def __del__(self):
        """Delete singlepoint results object"""
        global _libxtb
        Environment.__del__(self)
        if self._res is not None:
            _libxtb._lib.xtb_delResults(byref(self._res))

    def __len__(self):
        return self._natoms

    def get_energy(self):
        """Query singlepoint results object for energy"""
        global _libxtb
        _energy = c_double(0.0)
        _libxtb._lib.xtb_getEnergy(self._env, self._res, _energy)
        if self.check() != 0:
            raise ValueError("Energy is not available")
        return _energy.value

    def get_gradient(self):
        """Query singlepoint results object for gradient"""
        global _libxtb
        _gradient = np.zeros((len(self), 3))
        _libxtb._lib.xtb_getGradient(
            self._env, self._res, _cast(POINTER(c_double), _gradient)
        )
        if self.check() != 0:
            raise ValueError("Gradient is not available")
        return _gradient

    def get_virial(self):
        """Query singlepoint results object for virial"""
        global _libxtb
        _virial = np.zeros((3, 3))
        _libxtb._lib.xtb_getVirial(self._env, self._res, _cast(POINTER(c_double), _virial))
        if self.check() != 0:
            raise ValueError("Virial is not available")
        return _virial

    def get_dipole(self):
        """Query singlepoint results object for dipole"""
        global _libxtb
        _dipole = np.zeros(3)
        _libxtb._lib.xtb_getDipole(self._env, self._res, _cast(POINTER(c_double), _dipole))
        if self.check() != 0:
            raise ValueError("Dipole is not available")
        return _dipole

    def get_charges(self):
        """Query singlepoint results object for partial charges"""
        global _libxtb
        _charges = np.zeros(len(self))
        _libxtb._lib.xtb_getCharges(
            self._env, self._res, _cast(POINTER(c_double), _charges)
        )
        if self.check() != 0:
            raise ValueError("Charges are not available")
        return _charges

    def get_bond_orders(self):
        """Query singlepoint results object for bond orders"""
        global _libxtb
        _bond_orders = np.zeros((len(self), len(self)))
        _libxtb._lib.xtb_getBondOrders(
            self._env, self._res, _cast(POINTER(c_double), _bond_orders)
        )
        if self.check() != 0:
            raise ValueError("Bond orders are not available")
        return _bond_orders


class Calculator(Molecule):
    """Singlepoint calculator"""

    _calc = None

    def __init__(
        self,
        param: Param,
        numbers: List[int],
        positions: List[float],
        charge: Optional[float] = None,
        uhf: Optional[int] = None,
        lattice: Optional[List[float]] = None,
        periodic: Optional[List[bool]] = None,
        solvent: Optional[str] = None,
        library: Optional[CDLL] = None,
    ):
        """Create new calculator object"""
        global _libxtb
        Molecule.__init__(
            self, numbers, positions, charge, uhf, lattice, periodic, library,
        )

        self._loader = {
            Param.GFN2xTB: _libxtb._lib.xtb_loadGFN2xTB,
            Param.GFN1xTB: _libxtb._lib.xtb_loadGFN1xTB,
            Param.GFN0xTB: _libxtb._lib.xtb_loadGFN0xTB,
            Param.GFNFF: _libxtb._lib.xtb_loadGFNFF,
        }

        self._calc = c_void_p(_libxtb._lib.xtb_newCalculator())
        self._load(param)
        if solvent is not None:
            self.set_solvent(solvent)

    def __del__(self):
        """Delete calculator object"""
        global _libxtb
        Molecule.__del__(self)
        if self._calc is not None:
            _libxtb._lib.xtb_delCalculator(byref(self._calc))

    def _load(self, param: Param):
        """Load parametrisation data into calculator"""

        self._loader[param](
            self._env, self._mol, self._calc, None,
        )

        if self.check() != 0:
            raise ValueError("Could not load parametrisation data")

    def set_solvent(self, solvent: Optional[str]):
        """Add/Remove solvation model to/from calculator"""
        global _libxtb

        if solvent is None:
            _libxtb._lib.xtb_releaseSolvent(
                self._env, self._calc,
            )
        else:
            _solvent = solvent.encode()

            _libxtb._lib.xtb_setSolvent(
                self._env, self._calc, _solvent, None, None, None,
            )

        if self.check() != 0:
            raise ValueError("Could not load parametrisation data")

    def singlepoint(self, res: Optional[Results] = None) -> Results:
        """Perform singlepoint calculation"""
        global _libxtb

        _res = Results(self, _libxtb._lib) if res is None else res
        _libxtb._lib.xtb_singlepoint(
            self._env, self._mol, self._calc, _res._res,
        )

        if self.check() != 0:
            raise ValueError("Single point calculation failed")

        return _res


def _cast(ctype, array):
    """Cast a numpy array to a FFI pointer"""
    return None if array is None else array.ctypes.data_as(ctype)


def _set(lib, name, restype, argtypes):
    """Setup the ctypes signature for a function in a shared library"""
    func = lib.__getattr__(name)
    func.restype = restype
    func.argtypes = argtypes
