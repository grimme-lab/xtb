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

"""Tests for the ctypes interface to xtb."""


def test_library():
    """check if we can find the library and it looks okay"""
    from xtb.interface import load_library

    # check if we can load this one
    lib = load_library("libxtb")
    # check if we find some functions
    assert lib.GFN0_calculation is not None
    assert lib.GFN1_calculation is not None
    assert lib.GFN2_calculation is not None

    from xtb.interface import XTBLibrary

    libxtb = XTBLibrary(library=lib)
    assert libxtb.library is lib
    # check if we find some functions
    assert hasattr(libxtb, "GFN0Calculation")
    assert hasattr(libxtb, "GFN1Calculation")
    assert hasattr(libxtb, "GFN2Calculation")

    libxtb = XTBLibrary()
    # check if we find some functions
    assert hasattr(libxtb, "GFN0Calculation")
    assert hasattr(libxtb, "GFN1Calculation")
    assert hasattr(libxtb, "GFN2Calculation")


class Testing:
    """test all defined interfaces"""
    from ctypes import c_int, c_double
    from numpy import array
    from xtb.interface import XTBLibrary

    natoms = 24
    numbers = array(
        [6, 7, 6, 7, 6, 6, 6, 8, 7, 6, 8, 7, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        dtype=c_int,
    )
    positions = array(
        [
            [2.02799738646442, 0.09231312124713, -0.14310895950963],
            [4.75011007621000, 0.02373496014051, -0.14324124033844],
            [6.33434307654413, 2.07098865582721, -0.14235306905930],
            [8.72860718071825, 1.38002919517619, -0.14265542523943],
            [8.65318821103610, -1.19324866489847, -0.14231527453678],
            [6.23857175648671, -2.08353643730276, -0.14218299370797],
            [5.63266886875962, -4.69950321056008, -0.13940509630299],
            [3.44931709749015, -5.48092386085491, -0.14318454855466],
            [7.77508917214346, -6.24427872938674, -0.13107140408805],
            [10.30229550927022, -5.39739796609292, -0.13672168520430],
            [12.07410272485492, -6.91573621641911, -0.13666499342053],
            [10.70038521493902, -2.79078533715849, -0.14148379504141],
            [13.24597858727017, -1.76969072232377, -0.14218299370797],
            [7.40891694074004, -8.95905928176407, -0.11636933482904],
            [1.38702118184179, 2.05575746325296, -0.14178615122154],
            [1.34622199478497, -0.86356704498496, 1.55590600570783],
            [1.34624089204623, -0.86133716815647, -1.84340893849267],
            [5.65596919189118, 4.00172183859480, -0.14131371969009],
            [14.67430918222276, -3.26230980007732, -0.14344911021228],
            [13.50897177220290, -0.60815166181684, 1.54898960808727],
            [13.50780014200488, -0.60614855212345, -1.83214617078268],
            [5.41408424778406, -9.49239668625902, -0.11022772492007],
            [8.31919801555568, -9.74947502841788, 1.56539243085954],
            [8.31511620712388, -9.76854236502758, -1.79108242206824],
        ],
        dtype=c_double,
    )
    lib = XTBLibrary()

    def test_gfn0_interface(self):
        """check if the GFN0-xTB interface is working correctly."""
        thr = 1.0e-8
        from pytest import approx

        options = {
            "print_level": 1,
            "parallel": 0,
            "accuracy": 1.0,
            "electronic_temperature": 300.0,
            "gradient": True,
            "ccm": True,
            "solvent": "none",
        }
        kwargs = {
            "natoms": self.natoms,
            "numbers": self.numbers,
            "charge": 0.0,
            "positions": self.positions,
            "options": options,
            "output": "-",
            "magnetic_moment": 0,
        }
        results = self.lib.GFN0Calculation(**kwargs)
        assert results
        gnorm = abs(results["gradient"].flatten()).mean()
        assert approx(results["energy"], thr) == -40.90885036015837
        assert approx(gnorm, thr) == 0.00510542946913145

    def test_gfn1_interface(self):
        """check if the GFN1-xTB interface is working correctly."""
        thr = 1.0e-8
        from pytest import approx

        options = {
            "print_level": 1,
            "parallel": 0,
            "accuracy": 1.0,
            "electronic_temperature": 300.0,
            "gradient": True,
            "restart": False,
            "ccm": True,
            "max_iterations": 30,
            "solvent": "none",
        }
        args = (self.natoms, self.numbers, self.positions, options, 0.0, 0, "-")
        results = self.lib.GFN1Calculation(*args)
        assert results
        gnorm = abs(results["gradient"].flatten()).mean()
        assert approx(results["energy"], thr) == -44.50970242016477
        assert approx(gnorm, thr) == 0.004842110219908473

    def test_gfn2_gbsa_interface(self):
        """check if the GFN2-xTB interface is working correctly."""
        thr = 1.0e-8
        from pytest import approx

        options = {
            "print_level": 1,
            "parallel": 0,
            "accuracy": 1.0,
            "electronic_temperature": 300.0,
            "gradient": True,
            "restart": False,
            "ccm": True,
            "max_iterations": 30,
            "solvent": "ch2cl2",
        }
        results = self.lib.GFN2Calculation(
            natoms=self.natoms,
            numbers=self.numbers,
            charge=0.0,
            positions=self.positions,
            options=options,
            output="-",
            magnetic_moment=0,
        )
        assert results
        gnorm = abs(results["gradient"].flatten()).mean()
        assert approx(results["energy"], thr) == -42.170584528028506
        assert approx(gnorm, thr) == 0.004685246019376688
        assert approx(results["charges"].sum(), thr) == 0.0
