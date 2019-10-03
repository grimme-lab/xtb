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
"""Tests for the GBSA API."""


class Testing:
    from ase import Atoms

    octane = Atoms(
        "C8H18",
        [
            [0.00000000, 0.00000000, -0.78216596],
            [0.00000000, 0.00000000, 0.78216596],
            [1.13685895, -0.86573057, -1.33679656],
            [0.18131519, 1.41741401, -1.33679656],
            [-1.31817414, -0.55168345, -1.33679656],
            [1.13685895, 0.86573057, 1.33679656],
            [0.18131519, -1.41741401, 1.33679656],
            [-1.31817414, 0.55168345, 1.33679656],
            [2.10137753, -0.60632401, -0.89602505],
            [0.95640314, -1.92812772, -1.16563611],
            [1.21642453, -0.71505916, -2.41638078],
            [1.19160602, 1.79233328, -1.16563611],
            [0.01104713, 1.41098413, -2.41638078],
            [-0.52559677, 2.12300833, -0.89602505],
            [-1.57578076, -1.51668432, -0.89602505],
            [-2.14800916, 0.13579445, -1.16563611],
            [-1.22747167, -0.69592497, -2.41638078],
            [2.10137753, 0.60632401, 0.89602505],
            [0.95640314, 1.92812772, 1.16563611],
            [1.21642453, 0.71505916, 2.41638078],
            [1.19160602, -1.79233328, 1.16563611],
            [0.01104713, -1.41098413, 2.41638078],
            [-0.52559677, -2.12300833, 0.89602505],
            [-1.57578076, 1.51668432, 0.89602505],
            [-2.14800916, -0.13579445, 1.16563611],
            [-1.22747167, 0.69592497, 2.41638078],
        ],
    )

    def test_gfn1_gbsa_water(self):
        from pytest import approx
        from xtb.solvation import GBSA
        from xtb.calculators import GFN1

        born = [
            3.14547888,
            3.14547888,
            2.64032576,
            2.64032576,
            2.64032576,
            2.64032576,
            2.64032576,
            2.64032576,
            2.0475772,
            2.04527513,
            2.02622599,
            2.04527513,
            2.02622599,
            2.0475772,
            2.0475772,
            2.04527513,
            2.02622599,
            2.0475772,
            2.04527513,
            2.02622599,
            2.04527513,
            2.02622599,
            2.0475772,
            2.0475772,
            2.04527513,
            2.02622599,
        ]
        sasa = [
            -0.0,
            -0.0,
            1.81990386,
            1.83838909,
            1.81669523,
            1.81990386,
            1.83838909,
            1.81669523,
            11.40692384,
            12.57434505,
            14.97110513,
            12.56651178,
            14.98457022,
            11.40554208,
            11.41537855,
            12.52852953,
            14.95384903,
            11.40692384,
            12.57434505,
            14.97110513,
            12.56651178,
            14.98457022,
            11.40554208,
            11.41537855,
            12.52852953,
            14.95384903,
        ]

        mol = self.octane.copy()
        gbsa = GBSA(solvent="water", calc=GFN1())

        mol.set_calculator(gbsa)
        gsolv = mol.get_potential_energy()

        assert approx(gsolv) == 0.048946043429737074
        assert approx(gbsa.get_property("born radii")) == born
        assert approx(gbsa.get_property("sasa")) == sasa
