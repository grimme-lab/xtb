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
"""All solvation related implementations for xtb."""

from __future__ import print_function
from typing import List

from ctypes import c_int, c_double

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Bohr

import numpy as np

from xtb.interface import XTBLibrary
from xtb.calculators import XTB, GFN0, GFN1, GFN2

__all__ = ['GBSA', 'GSOLV', 'GSOLV_REFERENCE', 'MOL1BAR']

GSOLV = 0
GSOLV_REFERENCE = 1
MOL1BAR = 2


class GBSA(Calculator):
    """Interface to the GBSA method"""
    implemented_properties = [
        'energy',
        'born radii',
        'sasa',
    ]
    default_parameters = {
        'solvent': None,
        'temperature': 298.15,
        'reference': GSOLV,
        'grid size': 230,
        'method': 2,
    }
    calc = None
    _debug = False

    # pylint: disable=too-many-arguments
    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='gbsa', atoms=None, library=None, calc=None, **kwargs):
        """Construct the XTB base calculator object."""

        Calculator.__init__(self, restart, ignore_bad_restart_file, label, atoms,
                            **kwargs)

        if calc is not None:
            if not isinstance(calc, XTB):
                raise ValueError("Only xTB-calculators are supported in GBSA.")
            self.calc = calc
            # borrow library from xTB calculator
            self.library = self.calc.library
        else:
            if library is not None:
                self.library = library
            else:
                self.library = XTBLibrary()

        # loads the default parameters and updates with actual values
        self.parameters = self.get_default_parameters()
        # if we have a calculator we adjust the method keyword
        if isinstance(calc, (GFN2, GFN0)):
            self.set(method=2)  # GFN0 is treated as GFN2 (no parameters yet)
        if isinstance(calc, GFN1):
            self.set(method=1)
        # now set all parameters
        self.set(**kwargs)

        if self.parameters['solvent'] is None:
            raise ValueError("Solvent keyword is missing, cannot setup GBSA.")

    def output_file_name(self) -> str:
        """create output file name from label as ctype"""
        if self.label is not None:
            return self.label + ".out"
        return "-"  # standard output

    def create_arguments(self) -> dict:
        """create a list of arguments."""
        if any(self.atoms.pbc):
            raise NotImplementedError("GBSA is not available with PBC!")
        kwargs = {
            'natoms': len(self.atoms),
            'numbers': np.array(self.atoms.get_atomic_numbers(), dtype=c_int),
            'positions': np.array(self.atoms.get_positions()/Bohr, dtype=c_double),
            'options': self.parameters,
            'output': self.output_file_name(),
        }
        return kwargs

    # pylint: disable=dangerous-default-value
    def calculate(self, atoms=None, properties: List[str] = None,
                  system_changes: List[str] = all_changes) -> None:
        """perform calculation with libxtb"""

        if not properties:
            properties = ['energy']
        Calculator.calculate(self, atoms, properties, system_changes)

        kwargs = self.create_arguments()
        results = self.library.GBSACalculation(**kwargs)
        self.results['born radii'] = results['born']*Bohr
        self.results['sasa'] = results['sasa']*Bohr**2

        if self.calc is not None:
            # evaluate system in gas phase
            self.calc.set(solvent="none")
            self.atoms.set_calculator(self.calc)
            vacuum = self.atoms.get_potential_energy()
            # evaluate system in solution
            self.calc.reset()
            self.calc.set(solvent=self.parameters['solvent'])
            self.atoms.set_calculator(self.calc)
            solution = self.atoms.get_potential_energy()
            # derive properties like solvation free energy
            self.results['energy'] = solution - vacuum
