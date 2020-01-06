#!/usr/bin/env python
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

"""Wrapper script to perform geometry optimizations with xtb Python wrapper."""

if __name__ == "__main__":
    import argparse

    def get_cmd_args():
        """Commandline arguments for xtb wrapper in script mode"""
        parser = argparse.ArgumentParser(
            prog="xtb", description="Wrapper for xTB calculation."
        )

        parser.add_argument(
            "-f",
            "--format",
            dest="filetype",
            action="store",
            default=None,
            help="Format of input geometry (automatic)",
        )

        parser.add_argument(dest="filename", action="store", help="Input geometry")

        parser.add_argument(
            "-x",
            "--method",
            dest="method",
            action="store",
            default="gfn2",
            help="SQM method for calculation (gfn2)",
        )

        parser.add_argument(
            "-c",
            "--charge",
            dest="charge",
            action="store",
            type=int,
            default=0,
            help="Total charge (0)",
        )

        parser.add_argument(
            "--etemp",
            dest="etemp",
            action="store",
            type=float,
            default=300.0,
            help="Electronic temperature (300K)",
        )

        parser.add_argument(
            "--accuracy",
            dest="acc",
            action="store",
            type=float,
            default=1.0,
            help="Calculation accuracy (1.0)",
        )

        parser.add_argument(
            "--maxiter",
            dest="maxiter",
            action="store",
            type=int,
            default=250,
            help="Maximum number of SCC iterations (250)",
        )

        parser.add_argument(
            "--lcovb",
            dest="rnn",
            action="store",
            type=float,
            default=5.0,
            help="Longest covalent bond for preconditioning (5.0)",
        )

        parser.add_argument(
            "--econv",
            dest="econv",
            action="store",
            type=float,
            default=5.4423e-4,
            help="Convergence threshold for energy",
        )

        parser.add_argument(
            "--gconv",
            dest="fconv",
            action="store",
            type=float,
            default=3.88938e-05,
            help="Convergence threshold for forces",
        )

        parser.add_argument(
            "--linesearch",
            dest="armijo",
            action="store_true",
            help="Use line search"
        )

        parser.add_argument(
            "--precon",
            dest="precon",
            action="store_true",
            help="Use preconditioning"
        )

        parser.add_argument(
            "--optcell",
            dest="optcell",
            action="store_true",
            help="Relax cellparameters",
        )

        parser.add_argument(
            "--logfile",
            dest="logfile",
            action="store",
            default="-",
            help="File for logging information (STDOUT)",
        )

        parser.add_argument(
            "--trajectory",
            dest="trajectory",
            action="store",
            default="xtbopt.traj",
            help="File for trajectory output (xtbopt.traj)",
        )

        parser.add_argument(
            "--solvent",
            dest="solvent",
            action="store",
            default="none",
            help="Solvent for GBSA (none)",
        )

        return parser.parse_args()

    from math import sqrt
    from ase.io import read, write
    from ase.optimize.precon import Exp, PreconFIRE
    from ase.constraints import ExpCellFilter
    from ase.units import Hartree
    from xtb import GFN1, GFN2, GFN0

    # overwrite convergence thresholds of PreconFIRE optimizer
    class PatchedOptimizer(PreconFIRE):
        """Patches the ASE optimizer to use a different convergence threshold"""

        econv = None
        fconv = None
        elast = None

        def initialize(self):
            PreconFIRE.initialize(self)
            self.elast = None

        def run(
            self, steps=100000000, econv=0.00054423, fconv=3.88938e-05
        ):  # pylint: disable=arguments-differ
            self.econv = econv
            self.fconv = fconv
            PreconFIRE.run(self, 0.05, steps)

        def converged(self, forces=None):

            # get current energy and check for convergence in energy change
            ecurr = self.atoms.get_potential_energy()
            if self.elast is not None:
                ediff = ecurr - self.elast
                econverged = ediff < 0.0 and abs(ediff) < self.econv
            else:
                econverged = False
            # save potential energy
            self.elast = ecurr

            # check for convergence in forces
            if forces is None:
                forces = self.atoms.get_forces()
            fnorm = sqrt((forces ** 2).sum())
            fconverged = fnorm < self.fconv

            # print("--> norm(F):", fnorm, "converged?", fconverged)
            # print("--> energy :", ecurr, "converged?", econverged)

            return econverged and fconverged

    # pylint: disable=invalid-name
    args = get_cmd_args()

    if args.filetype is not None:
        mol = read(args.filename, format=str(args.filetype))
    else:
        mol = read(args.filename)

    parameters = {
        "accuracy": args.acc,
        "electronic_temperature": args.etemp,
        "max_iterations": args.maxiter,
        "solvent": args.solvent,
        "print_level": 2,
    }

    if args.method == "gfn0":
        calc = GFN0(**parameters)
    elif args.method == "gfn1":
        if mol.pbc.any():
            raise Exception("GFN1-xTB is not available with PBC")
        calc = GFN1(**parameters)
    elif args.method == "gfn2":
        if mol.pbc.any():
            raise Exception("GFN2-xTB is not available with PBC")
        calc = GFN2(**parameters)
    else:
        raise Exception("Method not implemented.")
    mol.set_calculator(calc)

    e = mol.get_potential_energy()

    print("Initial energy: eV, Eh", e, e / Hartree)

    if args.precon:
        precon = Exp(A=3, r_NN=args.rnn + 0.1, r_cut=args.rnn + 0.6)
    else:
        precon = None

    if args.optcell:
        sf = ExpCellFilter(mol)
        relax = PatchedOptimizer(
            sf,
            precon=precon,
            trajectory=args.trajectory,
            logfile=args.logfile,
            use_armijo=args.armijo,
        )
    else:
        relax = PatchedOptimizer(
            mol,
            precon=precon,
            trajectory=args.trajectory,
            logfile=args.logfile,
            use_armijo=args.armijo,
        )

    try:
        relax.run(econv=args.econv, fconv=args.fconv)
    except KeyboardInterrupt:
        print("User got impatient")
    except RuntimeError:
        print("Optimization terminated due to internal error")

    e = mol.get_potential_energy()
    print("Final energy: eV, Eh", e, e / Hartree)

    if mol.pbc.any():
        write("xtbopt.POSCAR", mol, format="vasp")
    else:
        write("xtbopt.xyz", mol, format="xyz")
