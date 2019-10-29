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
"""ASE Calculator implementation for the xtb program."""

from __future__ import print_function
from typing import List

from ctypes import c_int, c_double, c_bool

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Hartree, Bohr

import numpy as np

from xtb.interface import XTBLibrary
__all__ = ["GFN2", "GFN1", "GFN0", "XTB"]


class XTB(Calculator):
    """Base calculator for xtb related methods."""
    implemented_properties = []  # type: List[str]
    default_parameters = {}  # type: dict

    _debug = False

    # pylint: disable=too-many-arguments
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label=None, atoms=None, library=None, **kwargs):
        """Construct the XTB base calculator object."""

        Calculator.__init__(self, restart, ignore_bad_restart_file, label, atoms,
                            **kwargs)

        if library is not None:
            self.library = library
        else:
            self.library = XTBLibrary()

        def missing_library_call(**kwargs) -> None:
            """Raise an error if not set by child class."""
            raise NotImplementedError("Child class must replace function call!")

        self.library_call = missing_library_call

        # loads the default parameters and updates with actual values
        self.parameters = self.get_default_parameters()
        # now set all parameters
        self.set(**kwargs)

    def output_file_name(self) -> str:
        """create output file name from label as ctype"""
        if self.label is not None:
            return self.label + ".out"
        return "-"  # standard output

    def create_arguments(self) -> dict:
        """create a list of arguments."""
        raise NotImplementedError("Child class must implement argument list!")

    def store_results(self, results: dict) -> None:
        """store results after successful calculation."""
        raise NotImplementedError("Child class must implement result processing!")

    # pylint: disable=dangerous-default-value
    def calculate(self, atoms=None, properties: List[str] = None,
                  system_changes: List[str] = all_changes) -> None:
        """calculation interface to libxtb"""

        if not properties:
            properties = ['energy']
        Calculator.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

        kwargs = self.create_arguments()
        results = self.library_call(**kwargs)
        self.store_results(results)


class GFN0(XTB):
    """actual implementation of the GFN0-xTB calculator."""
    implemented_properties = [
        'energy',
        'free_energy',
        'forces',
        'stress',
    ]
    default_parameters = {
        'print_level': 2,
        'parallel': 0,
        'accuracy': 1.0,
        'electronic_temperature': 300.0,
        'gradient': True,
        'ccm': True,
        'solvent': 'none',
    }

    # pylint: disable=too-many-arguments
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='gfn0', atoms=None, library=None, **kwargs):
        """Construct the GFN0-xTB calculator object."""

        XTB.__init__(self, restart, ignore_bad_restart_file, label, atoms, library,
                     **kwargs)

        self.library_call = self.library.GFN0Calculation

    def create_arguments(self) -> dict:
        """create a list of arguments."""
        kwargs = {
            'natoms': len(self.atoms),
            'numbers': np.array(self.atoms.get_atomic_numbers(), dtype=c_int),
            'charge': self.atoms.get_initial_charges().sum().round(),
            'magnetic_moment': int(self.atoms.get_initial_magnetic_moments()
                                   .sum().round()),
            'positions': np.array(self.atoms.get_positions()/Bohr, dtype=c_double),
            'cell': np.array(self.atoms.get_cell()/Bohr, dtype=c_double),
            'pbc': np.array(self.atoms.get_pbc(), dtype=c_bool),
            'options': self.parameters,
            'output': self.output_file_name(),
        }
        return kwargs

    def store_results(self, results: dict) -> None:
        """store results after successful calculation."""
        self.results['energy'] = results['energy']*Hartree
        self.results['free_energy'] = self.results['energy']
        self.results['forces'] = -results['gradient']*Hartree/Bohr
        if 'stress tensor' in results:
            self.results['stress'] = results['stress tensor']*Hartree/Bohr**3


class GFN1(XTB):
    """actual implementation of the GFN1-xTB calculator."""
    implemented_properties = [
        'energy',
        'free_energy',
        'forces',
        'dipole',
        'charges',
        'wiberg',
    ]
    default_parameters = {
        'print_level': 2,
        'parallel': 0,
        'accuracy': 1.0,
        'electronic_temperature': 300.0,
        'gradient': True,
        'restart': True,
        'max_iterations': 250,
        'solvent': 'none',
        'ccm': True,
    }

    # pylint: disable=too-many-arguments
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='gfn1', atoms=None, library=None, **kwargs):
        """Construct the GFN1-xTB calculator object."""

        XTB.__init__(self, restart, ignore_bad_restart_file, label, atoms, library,
                     **kwargs)

        self.library_call = self.library.GFN1Calculation

    def create_arguments(self) -> dict:
        """create a list of arguments."""
        if any(self.atoms.pbc):
            raise NotImplementedError("GFN1-xTB is not available with PBC!")
        kwargs = {
            'natoms': len(self.atoms),
            'numbers': np.array(self.atoms.get_atomic_numbers(), dtype=c_int),
            'charge': self.atoms.get_initial_charges().sum().round(),
            'magnetic_moment': int(self.atoms.get_initial_magnetic_moments()
                                   .sum().round()),
            'positions': np.array(self.atoms.get_positions()/Bohr, dtype=c_double),
            'options': self.parameters,
            'output': self.output_file_name(),
        }
        return kwargs

    def store_results(self, results: dict) -> None:
        """store results after successful calculation."""
        self.results['energy'] = results['energy']*Hartree
        self.results['free_energy'] = self.results['energy']
        self.results['forces'] = -results['gradient']*Hartree/Bohr
        self.results['dipole'] = results['dipole moment']*Bohr
        self.results['charges'] = results['charges']
        self.results['wiberg'] = results['wiberg']


class GFN2(XTB):
    """actual implementation of the GFN2-xTB calculator."""
    implemented_properties = [
        'energy',
        'free_energy',
        'forces',
        'dipole',
        'charges',
        'dipoles',
        'quadrupoles',
        'wiberg',
    ]
    default_parameters = {
        'print_level': 2,
        'parallel': 0,
        'accuracy': 1.0,
        'electronic_temperature': 300.0,
        'gradient': True,
        'restart': True,
        'max_iterations': 250,
        'solvent': 'none',
        'ccm': True,
    }

    # pylint: disable=too-many-arguments
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='gfn2', atoms=None, library=None, **kwargs):
        """Construct the GFN2-xTB calculator object."""

        XTB.__init__(self, restart, ignore_bad_restart_file, label, atoms, library,
                     **kwargs)

        self.library_call = self.library.GFN2Calculation

    def create_arguments(self) -> dict:
        """create a list of arguments."""
        if any(self.atoms.pbc):
            raise NotImplementedError("GFN2-xTB is not available with PBC!")
        kwargs = {
            'natoms': len(self.atoms),
            'numbers': np.array(self.atoms.get_atomic_numbers(), dtype=c_int),
            'charge': self.atoms.get_initial_charges().sum().round(),
            'magnetic_moment': int(self.atoms.get_initial_magnetic_moments()
                                   .sum().round()),
            'positions': np.array(self.atoms.get_positions()/Bohr, dtype=c_double),
            'options': self.parameters,
            'output': self.output_file_name(),
        }
        return kwargs

    def store_results(self, results: dict) -> None:
        """store results after successful calculation."""
        self.results['energy'] = results['energy']*Hartree
        self.results['free_energy'] = self.results['energy']
        self.results['forces'] = -results['gradient']*Hartree/Bohr
        self.results['dipole'] = results['dipole moment']*Bohr
        self.results['charges'] = results['charges']
        self.results['dipoles'] = results['dipoles']*Bohr
        self.results['quadrupoles'] = results['quadrupoles']*Bohr**2
        self.results['wiberg'] = results['wiberg']
