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

"""Tests for the ASE interface to xtb."""


class Testing:
    """Unit tests for xTB methods."""

    def test_g2_gfn0(self):  # pylint: disable=no-self-use
        """Short test for GFN0-xTB on a subset of the G2."""
        from pytest import approx
        from ase.collections import g2
        from xtb.calculators import GFN0 as GFN

        thr = 1.0e-7
        subset = {
            "cyclobutene": (-318.7866120321342, 0.18172370289855405),
            "CH3ONO": (-358.73445196501973, 0.409421901652166),
            "SiH3": (-93.9024809777723, 0.27588327621925385),
            "C3H6_D3h": (-262.5099874403974, 0.10397633188045725),
            "CO2": (-236.48715959452625, 0.7513923402713432),
            "NO": (-165.74409502219441, 0.06375042027444998),
            "SiO": (-136.30493006096225, 0.41855575434941694),
            "C3H4_D2d": (-231.20208222258086, 0.11557009795103503),
            "COF2": (-357.1115761281822, 0.45554166948678626),
            "2-butyne": (-319.609294498653, 0.1654237936464811),
            "C2H5": (-188.90399485898985, 0.08291942020794801),
            "BF3": (-356.10418929536525, 0.15840229870826247),
            "N2O": (-246.88770058584996, 1.5739076279140762),
        }

        calc = GFN(accuracy=0.01)
        for name, data in subset.items():
            atoms = g2[name]
            energy, gnorm = data
            atoms.set_calculator(calc)
            assert energy == approx(atoms.get_potential_energy(), thr)
            assert gnorm == approx(abs(atoms.get_forces().flatten()).mean(), thr)

    def test_g2_gfn1(self):  # pylint: disable=no-self-use
        """Short test for GFN1-xTB on a subset of the G2."""
        from pytest import approx
        from ase.collections import g2
        from xtb.calculators import GFN1 as GFN

        thr = 1.0e-7
        subset = {
            "C4H4NH": (-389.83886433843026, 0.11845899422716655),
            "CH3SCH3": (-303.42619191556383, 0.06645628082687234),
            "SiH2_s3B1d": (-74.46528605312987, 0.08310436985665531),
            "CH3CO": (-284.3582641984081, 0.20638671880968007),
            "CO": (-183.19840869756732, 0.9977431424545925),
            "H2CO": (-213.4823225876121, 0.39657080790031),
            "CH3COOH": (-429.92204679760215, 0.24893030739923092),
            "HCF3": (-492.1472682009942, 0.2744291603498096),
            "CS2": (-257.45587723444686, 0.1020065624570768),
            "SiH2_s1A1d": (-77.00356692788125, 0.08145702394953913),
            "C4H4S": (-388.5560517704511, 0.16524901027806613),
        }

        calc = GFN(accuracy=0.01)
        for name, data in subset.items():
            atoms = g2[name]
            energy, gnorm = data
            atoms.set_calculator(calc)
            assert energy == approx(atoms.get_potential_energy(), thr)
            assert gnorm == approx(abs(atoms.get_forces().flatten()).mean(), thr)

    def test_g2_gfn2(self):  # pylint: disable=no-self-use
        """Short test for GFN2-xTB on a subset of the G2."""
        from pytest import approx
        from ase.collections import g2
        from xtb.calculators import GFN2 as GFN

        thr = 1.0e-7
        subset = {
            "F2O": (-362.11924747368295, 0.7143257863045354),
            "SO2": (-309.1999361327123, 0.8740554743818797),
            "H2CCl2": (-333.2265933353902, 0.046195744939883536),
            "CF3CN": (-580.4856945730149, 0.6213782830693683),
            "HCN": (-149.6789224980101, 1.015081156428575),
            "C2H6NH": (-292.6330765170243, 0.09848149585026728),
            "OCS": (-257.1846815478052, 0.549785128276834),
            "ClO": (-230.53675733271504, 0.6103015273199656),
            "C3H8": (-285.74053465625167, 0.08676710697381718),
            "HF": (-142.12158161038408, 0.024290012553866525),
            "O2": (-215.0347545913049, 0.8388557044315267),
            "SO": (-196.79890058895995, 0.8295962355319549),
            "NH": (-87.11950843824161, 0.18471418108510287),
            "C2F4": (-630.1300665999896, 0.11452334981560626),
            "NF3": (-461.2701116859826, 0.4165894734601184),
            "CH2_s3B1d": (-79.94618154681717, 0.1971528027462761),
            "CH3CH2Cl": (-309.61292017128613, 0.05533781132943324),
        }

        calc = GFN(accuracy=0.01)
        for name, data in subset.items():
            atoms = g2[name]
            energy, gnorm = data
            atoms.set_calculator(calc)
            assert energy == approx(atoms.get_potential_energy(), thr)
            assert gnorm == approx(abs(atoms.get_forces().flatten()).mean(), thr)
