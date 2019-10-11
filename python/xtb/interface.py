# This file is part of xtb.
#
# Copyright (C) 2019 Sebastian Ehlert
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

from typing import Optional

from ctypes import Structure, c_int, c_double, c_bool, c_char_p, c_char, \
                   POINTER, cdll, CDLL

import os.path as op
import numpy as np

# seems like ctypeslib is not always available
try:
    as_ctype = np.ctypeslib.as_ctypes_type  # pylint:disable=invalid-name
except AttributeError:
    as_ctype = None  # pylint:disable=invalid-name

__all__ = ['SCCOptions', 'PEEQOptions', 'XTBLibrary', 'load_library']


def load_library(libname: str) -> CDLL:
    """load library cross-platform compatible."""

    if not op.splitext(libname)[1]:
        # Try to load library with platform-specific name
        from numpy.distutils.misc_util import get_shared_lib_extension
        so_ext = get_shared_lib_extension()
        libname_ext = libname + so_ext
    else:
        libname_ext = libname

    return cdll.LoadLibrary(libname_ext)


class _Structure_(Structure):  # pylint: disable=invalid-name,protected-access
    """patch Structure class to allow returning it arguments as dict."""
    def to_dict(self) -> dict:
        """return structure as dictionary."""
        return {key: getattr(self, key) for key, _ in self._fields_}

    def to_list(self) -> list:
        """return structure as list. Order is the same as in structure."""
        return [getattr(self, key) for key, _ in self._fields_]


class SCCOptions(_Structure_):
    """Options for evaluating a SCC-Hamiltonian."""
    _fields_ = [
        ('print_level', c_int),
        ('parallel', c_int),
        ('accuracy', c_double),
        ('electronic_temperature', c_double),
        ('gradient', c_bool),
        ('restart', c_bool),
        ('ccm', c_bool),
        ('max_iterations', c_int),
        ('solvent', c_char*20),
    ]


class PEEQOptions(_Structure_):
    """Options for evaluating a EEQ-Hamiltonian."""
    _fields_ = [
        ('print_level', c_int),
        ('parallel', c_int),
        ('accuracy', c_double),
        ('electronic_temperature', c_double),
        ('gradient', c_bool),
        ('ccm', c_bool),
        ('solvent', c_char*20),
    ]


def check_ndarray(array: np.ndarray, ctype, size: int, name="array") -> None:
    """check if we got the correct array data"""
    if not isinstance(array, np.ndarray):
        raise ValueError("{} must be of type ndarray".format(name))
    if array.size != size:
        raise ValueError("{} does not have the correct size of {}"
                         .format(name, size))
    if as_ctype is not None:
        if as_ctype(array.dtype) != ctype:
            raise ValueError("{} must be of {} compatible type"
                             .format(name, ctype))


class XTBLibrary:
    """wrapper for the xtb shared library"""

    # define periodic GFN0-xTB interface
    _GFN0_PBC_calculation_ = (
        POINTER(c_int),  # number of atoms
        POINTER(c_int),  # atomic numbers, dimension(number of atoms)
        POINTER(c_double),  # molecular charge
        POINTER(c_int),  # number of unpaired electrons
        POINTER(c_double),  # cartesian coordinates, dimension(3*number of atoms)
        POINTER(c_double),  # lattice parameters, dimension(9)
        POINTER(c_bool),  # periodicity of the system
        POINTER(PEEQOptions),
        c_char_p,  # output file name
        POINTER(c_double),  # energy
        POINTER(c_double),  # gradient, dimension(3*number of atoms)
        POINTER(c_double),  # stress tensor, dimension(9)
        POINTER(c_double),  # lattice gradient, dimension(9)
    )

    # define GFN0-xTB interface
    _GFN0_calculation_ = (
        POINTER(c_int),  # number of atoms
        POINTER(c_int),  # atomic numbers, dimension(number of atoms)
        POINTER(c_double),  # molecular charge
        POINTER(c_int),  # number of unpaired electrons
        POINTER(c_double),  # cartesian coordinates, dimension(3*number of atoms)
        POINTER(PEEQOptions),
        c_char_p,  # output file name
        POINTER(c_double),  # energy
        POINTER(c_double),  # gradient, dimension(3*number of atoms)
    )

    # define GFN1-xTB interface
    _GFN1_calculation_ = (
        POINTER(c_int),  # number of atoms
        POINTER(c_int),  # atomic numbers, dimension(number of atoms)
        POINTER(c_double),  # molecular charge
        POINTER(c_int),  # number of unpaired electrons
        POINTER(c_double),  # cartesian coordinates, dimension(3*number of atoms)
        POINTER(SCCOptions),
        c_char_p,  # output file name
        POINTER(c_double),  # energy
        POINTER(c_double),  # gradient, dimension(3*number of atoms)
    )

    # define GFN2-xTB interface
    _GFN2_calculation_ = (
        POINTER(c_int),  # number of atoms
        POINTER(c_int),  # atomic numbers, dimension(number of atoms)
        POINTER(c_double),  # molecular charge
        POINTER(c_int),  # number of unpaired electrons
        POINTER(c_double),  # cartesian coordinates, dimension(3*number of atoms)
        POINTER(SCCOptions),
        c_char_p,  # output file name
        POINTER(c_double),  # energy
        POINTER(c_double),  # gradient, dimension(3*number of atoms)
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
    )

    _GBSA_model_preload_ = (
        POINTER(c_double),  # dielectric data
        POINTER(c_double),  # molar mass (g/mol)
        POINTER(c_double),  # solvent density (g/cm^3)
        POINTER(c_double),  # Born radii
        POINTER(c_double),  # Atomic surfaces
        POINTER(c_double),  # Gshift (gsolv=reference vs. gsolv)
        POINTER(c_double),  # offset parameter (fitted)
        POINTER(c_double),
        POINTER(c_double),  # Surface tension (mN/m=dyn/cm), dimension(94)
        POINTER(c_double),  # dielectric descreening parameters, dimension(94)
        POINTER(c_double),  # hydrogen bond strength, dimension(94)
    )

    _GBSA_calculation_ = (
        POINTER(c_int),  # number of atoms
        POINTER(c_int),  # atomic numbers, dimension(number of atoms)
        POINTER(c_double),  # cartesian coordinates, dimension(3*number of atoms)
        c_char_p,  # solvent name
        POINTER(c_int),  # reference state (0: gsolv=reference, 1: gsolv)
        POINTER(c_double),  # temperature (only for reference=0 or 2
        POINTER(c_int),  # method (1 or 2)
        POINTER(c_int),  # angular grid size (available Lebedev-grids)
        c_char_p,  # output file name
        POINTER(c_double),  # Born radii, dimension(number of atoms)
        POINTER(c_double),  # SASA, dimension(number of atoms)
    )

    _lebedev_grids_ = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266,
                       302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354,
                       2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810]

    def __init__(self, library: Optional[CDLL] = None):
        """construct library from CDLL object."""
        if library is not None:
            self.library = library
        else:
            self.library = load_library('libxtb')

        self._set_argtypes_()

    def _set_argtypes_(self) -> None:
        """define all interfaces."""
        self.library.GFN0_calculation.argtypes = self._GFN0_calculation_
        self.library.GFN1_calculation.argtypes = self._GFN1_calculation_
        self.library.GFN2_calculation.argtypes = self._GFN2_calculation_
        self.library.GFN0_PBC_calculation.argtypes = self._GFN0_PBC_calculation_
        self.library.GBSA_model_preload.argtypes = self._GBSA_model_preload_
        self.library.GBSA_calculation.argtypes = self._GBSA_calculation_

    # pylint: disable=invalid-name, too-many-arguments, too-many-locals
    def GFN0Calculation(self, natoms: int, numbers, positions, options: dict,
                        charge: float = 0.0, magnetic_moment: int = 0,
                        output: str = "-", cell=None, pbc=None) -> dict:
        """wrapper for calling the GFN0 Calculator from the library."""
        periodic = cell is not None and pbc is not None

        check_ndarray(numbers, c_int, natoms, "numbers")
        check_ndarray(positions, c_double, 3*natoms, "positions")
        if periodic:
            check_ndarray(cell, c_double, 9, "cell")
            check_ndarray(pbc, c_bool, 3, "pbc")

        energy = c_double(0.0)
        gradient = np.zeros((natoms, 3), dtype=c_double)

        # turn all strings to binary data, such that ctypes will not complain
        l_options = {key: val.encode('utf-8') if isinstance(val, str) else val
                     for key, val in options.items()}

        if periodic:
            cell_gradient = np.zeros((3, 3), dtype=c_double)
            stress_tensor = np.zeros((3, 3), dtype=c_double)
            args = [
                c_int(natoms),
                numbers.ctypes.data_as(POINTER(c_int)),
                c_double(charge),
                c_int(magnetic_moment),
                positions.ctypes.data_as(POINTER(c_double)),
                cell.ctypes.data_as(POINTER(c_double)),
                pbc.ctypes.data_as(POINTER(c_bool)),
                PEEQOptions(**l_options),
                output.encode('utf-8'),
                energy,
                gradient.ctypes.data_as(POINTER(c_double)),
                stress_tensor.ctypes.data_as(POINTER(c_double)),
                cell_gradient.ctypes.data_as(POINTER(c_double)),
            ]
            stat = self.library.GFN0_PBC_calculation(*args)
        else:
            args = [
                c_int(natoms),
                numbers.ctypes.data_as(POINTER(c_int)),
                c_double(charge),
                c_int(magnetic_moment),
                positions.ctypes.data_as(POINTER(c_double)),
                PEEQOptions(**l_options),
                output.encode('utf-8'),
                energy,
                gradient.ctypes.data_as(POINTER(c_double)),
            ]
            stat = self.library.GFN0_calculation(*args)

        if stat != 0:
            raise RuntimeError("GFN0 calculation failed in xtb.")

        results = {
            'energy': energy.value,
            'gradient': gradient,
        }
        if periodic:
            results['cell gradient'] = cell_gradient
            results['stress tensor'] = stress_tensor
        return results

    # pylint: disable=invalid-name, too-many-arguments, too-many-locals
    def GFN1Calculation(self, natoms: int, numbers, positions, options: dict,
                        charge: float = 0.0, magnetic_moment: int = 0,
                        output: str = "-") -> dict:
        """wrapper for calling the GFN1 Calculator from the library."""

        check_ndarray(numbers, c_int, natoms, "numbers")
        check_ndarray(positions, c_double, 3*natoms, "positions")

        energy = c_double(0.0)
        gradient = np.zeros((natoms, 3), dtype=c_double)
        dipole = np.zeros(3, dtype=c_double)
        charges = np.zeros(natoms, dtype=c_double)
        wiberg = np.zeros((natoms, natoms), dtype=c_double)

        # turn all strings to binary data, such that ctypes will not complain
        l_options = {key: val.encode('utf-8') if isinstance(val, str) else val
                     for key, val in options.items()}

        args = [
            c_int(natoms),
            numbers.ctypes.data_as(POINTER(c_int)),
            c_double(charge),
            c_int(magnetic_moment),
            positions.ctypes.data_as(POINTER(c_double)),
            SCCOptions(**l_options),
            output.encode('utf-8'),
            energy,
            gradient.ctypes.data_as(POINTER(c_double)),
            dipole.ctypes.data_as(POINTER(c_double)),
            charges.ctypes.data_as(POINTER(c_double)),
            wiberg.ctypes.data_as(POINTER(c_double)),
        ]
        stat = self.library.GFN1_calculation(*args)

        if stat != 0:
            raise RuntimeError("GFN1 calculation failed in xtb.")

        return {
            'energy': energy.value,
            'gradient': gradient,
            'dipole moment': dipole,
            'charges': charges,
            'wiberg': wiberg,
        }

    # pylint: disable=invalid-name, too-many-arguments, too-many-locals
    def GFN2Calculation(self, natoms: int, numbers, positions, options: dict,
                        charge: float = 0.0, magnetic_moment: int = 0,
                        output: str = "-") -> dict:
        """wrapper for calling the GFN2 Calculator from the library."""

        check_ndarray(numbers, c_int, natoms, "numbers")
        check_ndarray(positions, c_double, 3*natoms, "positions")

        energy = c_double(0.0)
        gradient = np.zeros((natoms, 3), dtype=c_double)
        dipole = np.zeros(3, dtype=c_double)
        charges = np.zeros(natoms, dtype=c_double)
        dipoles = np.zeros((natoms, 3), dtype=c_double)
        quadrupoles = np.zeros((natoms, 6), dtype=c_double)
        wiberg = np.zeros((natoms, natoms), dtype=c_double)

        # turn all strings to binary data, such that ctypes will not complain
        l_options = {key: val.encode('utf-8') if isinstance(val, str) else val
                     for key, val in options.items()}

        args = [
            c_int(natoms),
            numbers.ctypes.data_as(POINTER(c_int)),
            c_double(charge),
            c_int(magnetic_moment),
            positions.ctypes.data_as(POINTER(c_double)),
            SCCOptions(**l_options),
            output.encode('utf-8'),
            energy,
            gradient.ctypes.data_as(POINTER(c_double)),
            dipole.ctypes.data_as(POINTER(c_double)),
            charges.ctypes.data_as(POINTER(c_double)),
            dipoles.ctypes.data_as(POINTER(c_double)),
            quadrupoles.ctypes.data_as(POINTER(c_double)),
            wiberg.ctypes.data_as(POINTER(c_double)),
        ]
        stat = self.library.GFN2_calculation(*args)

        if stat != 0:
            raise RuntimeError("GFN2 calculation failed in xtb.")

        return {
            'energy': energy.value,
            'gradient': gradient,
            'dipole moment': dipole,
            'charges': charges,
            'dipoles': dipoles,
            'quadrupoles': quadrupoles,
            'wiberg': wiberg,
        }

    # pylint: disable=invalid-name, too-many-arguments, too-many-locals
    def GBSACalculation(self, natoms: int, numbers, positions,
                        options: dict, output: str) -> dict:
        """wrapper for calling the GBSA Calculator from the library."""

        check_ndarray(numbers, c_int, natoms, "numbers")
        check_ndarray(positions, c_double, 3*natoms, "positions")

        if 'solvent' not in options:
            raise ValueError("solvent has to present in options")
        ref_state = options['reference'] if 'reference' in options else 0
        if ref_state not in [0, 1, 2]:
            raise ValueError("reference state {} is no defined".format(ref_state))
        temp = options['temperature'] if 'temperature' in options else 298.15
        if temp <= 0.0:
            raise ValueError("negative temperature does not make sense")
        grid = options['grid size'] if 'grid size' in options else 230
        if grid not in self._lebedev_grids_:
            raise ValueError("angular grid {} is no Lebedev grid".format(grid))
        method = options['method'] if 'method' in options else 2
        if method not in [1, 2]:
            raise ValueError("method {} is no available".format(method))

        born = np.zeros(natoms, dtype=c_double)
        sasa = np.zeros(natoms, dtype=c_double)

        args = (
            c_int(natoms),
            numbers.ctypes.data_as(POINTER(c_int)),
            positions.ctypes.data_as(POINTER(c_double)),
            options['solvent'].encode('utf-8'),
            c_int(ref_state),
            c_double(temp),
            c_int(method),
            c_int(grid),
            output.encode('utf-8'),
            born.ctypes.data_as(POINTER(c_double)),
            sasa.ctypes.data_as(POINTER(c_double)),
        )
        stat = self.library.GBSA_calculation(*args)

        if stat != 0:
            raise RuntimeError("GBSA calculation failed in xtb.")
        return {
            'born': born,
            'sasa': sasa,
        }

    def GBSA_model_preload(self, epsv: float, smass: float, rhos: float, c1: float,
                           rprobe: float, gshift: float, soset: float,
                           gamscale, sx, tmp) -> None:
        """preload parameters into GBSA"""

        check_ndarray(gamscale, c_double, 94, "gamscale")
        check_ndarray(sx, c_double, 94, "sx")
        check_ndarray(tmp, c_double, 94, "tmp")

        args = (
            c_double(epsv),
            c_double(smass),
            c_double(rhos),
            c_double(c1),
            c_double(rprobe),
            c_double(gshift),
            c_double(soset),
            None,
            gamscale.ctypes.data_as(POINTER(c_double)),
            sx.ctypes.data_as(POINTER(c_double)),
            tmp.ctypes.data_as(POINTER(c_double)),
        )

        stat = self.library.GBSA_model_preload(*args)

        if stat != 0:
            raise RuntimeError("GBSA parameters could not be loaded by xtb.")
