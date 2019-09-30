#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import ase
import ctypes

from ctypes import cdll, Structure, c_double, c_int, c_bool, POINTER, c_char_p

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.units import Hartree, Bohr

class SCC_options(Structure):
    _fields_ = [('prlevel', c_int),
                ('parallel',c_int),
                ('acc',     c_double),
                ('etemp',   c_double),
                ('grad',    c_bool),
                ('restart', c_bool),
                ('maxiter', c_int),
                ('solvent', ctypes.c_char*20)]

class PEEQ_options(Structure):
    _fields_ = [('prlevel', c_int),
                ('parallel',c_int),
                ('acc',     c_double),
                ('etemp',   c_double),
                ('grad',    c_bool),
                ('ccm',     c_bool)]

class XTB(Calculator):
    implemented_properties = ['energy', 'free_energy', 'forces']
    command = None

    def __init__(self, charge = None, accuracy = None, temperature = None,
            print_level = None,
            restart=None, ignore_bad_restart_file=False,
            label="saw", atoms=None, command=None, debug=False, **kwargs):
        """Construct the xtb-calculator object."""

        self._debug = debug
        self._lib = None
        self.label = None
        self.parameters = None
        self.results = None
        self.atoms = None

        if charge is not None:
            self.charge = charge
        else:
            self.charge = 0.0

        if accuracy is not None:
            self.accuracy = accuracy
        else:
            self.accuracy = 1.0

        if temperature is not None:
            self.temperature = temperature
        else:
            self.temperature = 300.0

        if print_level is not None:
            self.print_level = print_level
        else:
            self.print_level = 2

        if command is not None:
            self.command = command
        else:
            self.command = 'xtb' # default

        Calculator.__init__(self, restart, ignore_bad_restart_file,
                label, atoms, **kwargs)

        if restart is not None:
            try:
                self.read(restart)
            except:
                if ignore_bad_restart_file:
                    self.reset()
                else:
                    raise

    def set(self, **kwargs):
        """Set parameters like set(key1=value1, key2=value2, ...)."""
        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write(self, label):
        """Write atoms, parameters and calculated results into restart files."""
        if self._debug:
            print("Writting restart to: ", label)
        self.atoms.write(label + '_restart.traj')
        self.parameters.write(label + '_params.ase')
        open(label + '_results.ase', 'w').write(repr(self.results))

    def read(self, label):
        """Read atoms, parameters and calculated results from restart files."""
        self.atoms = ase.io.read(label + '_restart.traj')
        self.parameters = Parameters.read(label + '_params.ase')
        results_txt = open(label + '_results.ase').read()
        self.results = eval(results_txt, {'array': np.array})

    def _open_lib(self):
        # load xtb-library
        self._lib = cdll.LoadLibrary('libxtb.so')

        c_int_p = POINTER(c_int)
        c_bool_p = POINTER(c_bool)
        c_double_p = POINTER(c_double)

        # define periodic GFN0-xTB interface
        self._lib.GFN0_PBC_calculation.argtypes = [c_int_p, c_int_p,
            c_double_p, c_double_p, c_double_p,
            c_bool_p, POINTER(PEEQ_options), c_char_p,
            c_double_p, c_double_p, c_double_p]

        # define GFN0-xTB interface
        self._lib.GFN0_calculation.argtypes = [c_int_p, c_int_p,
            c_double_p, c_double_p, POINTER(PEEQ_options), c_char_p,
            c_double_p, c_double_p]

        # define GFN1-xTB interface
        self._lib.GFN1_calculation.argtypes = [c_int_p, c_int_p,
            c_double_p, c_double_p, POINTER(SCC_options), c_char_p,
            c_double_p, c_double_p]

        # define GFN2-xTB interface
        self._lib.GFN2_calculation.argtypes = [c_int_p, c_int_p,
            c_double_p, c_double_p, POINTER(SCC_options), c_char_p,
            c_double_p, c_double_p, c_double_p, c_double_p, c_double_p,
            c_double_p, c_double_p]

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """Do the calculation."""

        if not properties:
            properties =  ['energy']
        Calculator.calculate(self, atoms, properties, system_changes)

class GFN0(XTB):
    implemented_properties = ['energy', 'free_energy', 'forces', 'stress']
    command = None

    def __init__(self, charge = None, accuracy = None, temperature = None,
            print_level = None,
            restart=None, ignore_bad_restart_file=False,
            label="saw", atoms=None, command=None, debug=False, **kwargs):
        """Construct the xtb-calculator object."""

        XTB.__init__(self, charge, accuracy, temperature, print_level,
                restart, ignore_bad_restart_file, label, atoms, **kwargs)


    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """Do the calculation."""

        if not properties:
            properties =  ['energy']
        XTB.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

        if self._lib is None:
            self._open_lib()

        natoms  = c_int(atoms.get_number_of_atoms())
        attyp   = np.array(atoms.get_atomic_numbers(), dtype=np.int32)
        attyp_p = attyp.ctypes.data_as(POINTER(c_int))
        coord   = atoms.get_positions()/Bohr # from Angstrom to Bohr
        coord_p = coord.ctypes.data_as(POINTER(c_double))
        charge  = c_double(self.charge)

        energy   = c_double(0.0)
        gradient = np.zeros((3,natoms.value),order='F')
        grad_p   = gradient.ctypes.data_as(POINTER(c_double))

        opt = PEEQ_options(c_int(self.print_level),
                           c_int(0), # parallel -> automatic
                           c_double(self.accuracy),
                           c_double(self.temperature),
                           c_bool(True), # gradient -> always
                           c_bool(True)) # CCM -> yes

        outfile = "gfn0.out".encode('utf-8')

        stat = self._lib.GFN0_calculation(natoms,attyp_p,charge,
                coord_p,opt,outfile,energy,grad_p)
        # in case Fortran behaves we find a useful return value
        if (stat != 0):
            raise Exception("xtb terminated in error.")

        self.results['energy'] = energy.value * Hartree
        self.results['free_energy'] = energy.value * Hartree
        self.results['forces'] = - gradient.T * Hartree / Bohr

class GFN1(XTB):
    implemented_properties = ['energy', 'free_energy', 'forces']
    command = None

    def __init__(self, charge = None, accuracy = None, temperature = None,
            max_iterations = None, print_level = None, solvent = None,
            restart=None, ignore_bad_restart_file=False,
            label="saw", atoms=None, command=None, debug=False, **kwargs):
        """Construct the GFN2-calculator object."""

        if max_iterations is not None:
            self.max_iterations = max_iterations
        else:
            self.max_iterations = 250

        if solvent is not None:
            self.solvent = solvent
        else:
            self.solvent = "none"

        XTB.__init__(self, charge, accuracy, temperature, print_level,
                restart, ignore_bad_restart_file, label, atoms, **kwargs)

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """Do the calculation."""

        if not properties:
            properties =  ['energy']
        XTB.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

        if self._lib is None:
            self._open_lib()

        natoms  = c_int(atoms.get_number_of_atoms())
        attyp   = np.array(atoms.get_atomic_numbers(), dtype=np.int32)
        attyp_p = attyp.ctypes.data_as(POINTER(c_int))
        coord   = atoms.get_positions()/Bohr # from Angstrom to Bohr
        coord_p = coord.ctypes.data_as(POINTER(c_double))
        charge  = c_double(self.charge)

        energy   = c_double(0.0)
        gradient = np.zeros((3,natoms.value),order='F')
        grad_p   = gradient.ctypes.data_as(POINTER(c_double))

        opt = SCC_options(c_int(self.print_level),
                          c_int(0), # parallel -> automatic
                          c_double(self.accuracy),
                          c_double(self.temperature),
                          c_bool(True),  # gradient -> always
                          c_bool(False), # restart -> no
                          c_int(self.max_iterations),
                          self.solvent.encode('utf-8'))

        outfile = "gfn1.out".encode('utf-8')

        stat = self._lib.GFN1_calculation(natoms,attyp_p,charge,
                coord_p,opt,outfile,energy,grad_p)
        # in case Fortran behaves we find a useful return value
        if (stat != 0):
            raise Exception("xtb terminated in error.")


        self.results['energy'] = energy.value * Hartree
        self.results['free_energy'] = energy.value * Hartree
        self.results['forces'] = - gradient.T * Hartree / Bohr

class GFN2(XTB):
    implemented_properties = ['energy', 'free_energy', 'forces', 'dipole',
            'charges', 'dipoles', 'quadrupoles', 'wiberg']
    command = None

    def __init__(self, charge = None, accuracy = None, temperature = None,
            max_iterations = None, print_level = None, solvent = None,
            restart=None, ignore_bad_restart_file=False,
            label="saw", atoms=None, command=None, debug=False, **kwargs):
        """Construct the GFN2-calculator object."""

        if max_iterations is not None:
            self.max_iterations = max_iterations
        else:
            self.max_iterations = 250

        if solvent is not None:
            self.solvent = solvent
        else:
            self.solvent = "none"

        XTB.__init__(self, charge, accuracy, temperature, print_level,
                restart, ignore_bad_restart_file, label, atoms, **kwargs)

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """Do the calculation."""

        if not properties:
            properties =  ['energy']
        XTB.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

        if self._lib is None:
            self._open_lib()

        natoms  = c_int(atoms.get_number_of_atoms())
        attyp   = np.array(atoms.get_atomic_numbers(), dtype=np.int32)
        attyp_p = attyp.ctypes.data_as(POINTER(c_int))
        coord   = atoms.get_positions()/Bohr # from Angstrom to Bohr
        coord_p = coord.ctypes.data_as(POINTER(c_double))
        charge  = c_double(self.charge)
        self.results['charges']     = np.zeros(natoms.value)
        self.results['dipoles']     = np.zeros((3,natoms.value),order='F')
        self.results['quadrupoles'] = np.zeros((6,natoms.value),order='F')
        self.results['dipole'] = np.zeros(3)
        self.results['wiberg'] = np.zeros((natoms.value,natoms.value),order='F')
        q_p      = self.results['charges'].ctypes.data_as(POINTER(c_double))
        dipm_p   = self.results['dipoles'].ctypes.data_as(POINTER(c_double))
        qp_p     = self.results['quadrupoles'].ctypes.data_as(POINTER(c_double))
        dipole_p = self.results['dipole'].ctypes.data_as(POINTER(c_double))
        wbo_p    = self.results['wiberg'].ctypes.data_as(POINTER(c_double))

        energy   = c_double(0.0)
        gradient = np.zeros((3,natoms.value),order='F')
        grad_p   = gradient.ctypes.data_as(POINTER(c_double))

        opt = SCC_options(c_int(self.print_level),
                          c_int(0), # parallel -> automatic
                          c_double(self.accuracy),
                          c_double(self.temperature),
                          c_bool(True), # gradient -> always
                          c_bool(True), # restart  -> if possible
                          c_int(self.max_iterations),
                          self.solvent.encode('utf-8'))

        outfile = "gfn2.out".encode('utf-8')

        stat = self._lib.GFN2_calculation(natoms,attyp_p,charge,
            coord_p,opt,outfile,energy,grad_p,dipole_p,q_p,dipm_p,qp_p,wbo_p)
        # in case Fortran behaves we find a useful return value
        if (stat != 0):
            raise Exception("xtb terminated in error.")

        self.results['energy'] = energy.value * Hartree
        #print("-> GFN2:", energy.value * Hartree)
        self.results['free_energy'] = energy.value * Hartree
        self.results['forces'] = - gradient.T * Hartree / Bohr

        # some unit conversion for the sake of consistency
        self.results['dipole']  = self.results['dipole']  * Bohr
        self.results['dipoles'] = self.results['dipoles'] * Bohr
        self.results['quadrupoles'] = self.results['quadrupoles'] * Bohr*Bohr

class GFN0_PBC(XTB):
    implemented_properties = ['energy', 'free_energy', 'forces', 'stress']
    command = None

    def __init__(self, charge = None, accuracy = None, temperature = None,
            print_level = None,
            restart=None, ignore_bad_restart_file=False,
            label="saw", atoms=None, command=None, debug=False, **kwargs):
        """Construct the xtb-calculator object."""

        XTB.__init__(self, charge, accuracy, temperature, print_level,
                restart, ignore_bad_restart_file, label, atoms, **kwargs)


    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        """Do the calculation."""

        if not properties:
            properties =  ['energy']
        XTB.calculate(self, atoms, properties, system_changes)

        if self._debug:
            print("system_changes:", system_changes)

        if self._lib is None:
            self._open_lib()

        natoms  = c_int(atoms.get_number_of_atoms())
        attyp   = np.array(atoms.get_atomic_numbers(), dtype=np.int32)
        attyp_p = attyp.ctypes.data_as(POINTER(c_int))
        coord   = atoms.get_positions()/Bohr # from Angstrom to Bohr
        coord_p = coord.ctypes.data_as(POINTER(c_double))
        lattice = atoms.get_cell()/Bohr # from Angstrom to Bohr
        dlat_p  = lattice.ctypes.data_as(POINTER(c_double))
        charge  = c_double(self.charge)

        energy   = c_double(0.0)
        gradient = np.zeros((3,natoms.value),order='F')
        grad_p   = gradient.ctypes.data_as(POINTER(c_double))
        gradlatt = np.zeros((3,3),order='F')
        glat_p   = gradlatt.ctypes.data_as(POINTER(c_double))

        opt = PEEQ_options(c_int(self.print_level),
                           c_int(0), # parallel -> automatic
                           c_double(self.accuracy),
                           c_double(self.temperature),
                           c_bool(True), # gradient -> always
                           c_bool(True)) # CCM -> yes

        #pbc = np.array([True,True,True], dtype=np.bool)
        pbc = atoms.get_pbc()
        pbc_p = pbc.ctypes.data_as(POINTER(c_bool))

        outfile = "gfn0.out".encode('utf-8')

        # try/except is not working on Fortran
        stat = self._lib.GFN0_PBC_calculation(natoms,attyp_p,charge,
                coord_p,dlat_p,pbc_p,opt,outfile,energy,grad_p,glat_p)
        # in case Fortran behaves we find a useful return value
        if (stat != 0):
            raise Exception("xtb terminated in error.")

        self.results['energy'] = energy.value * Hartree
        self.results['free_energy'] = energy.value * Hartree
        self.results['forces'] = - gradient.T * Hartree / Bohr
        # xtb returns the cell gradient, so we have to backtransform
        # taking into account that the stress tensor should be symmetric,
        # we interface F and C ordered arrays and have to transpose we
        # can simply drop all reshape and transpose steps by writing
        stress = np.dot(self.atoms.cell,
                gradlatt * Hartree / Bohr / self.atoms.get_volume())
        self.results['stress'] = stress.flat[[0, 4, 8, 5, 2, 1]]

if __name__ == "__main__":

    import argparse

    def get_cmd_args():
        parser = argparse.ArgumentParser(prog='xtb',
                description='Wrapper for xTB calculation.')

        parser.add_argument("-f","--format",dest="filetype",action="store",
                default=None,help="Format of input geometry (automatic)")

        parser.add_argument(dest="filename",action="store",
                help="Input geometry")

        parser.add_argument("-x","--method",dest="method",action="store",
                default='gfn2',help="SQM method for calculation (gfn2)")

        parser.add_argument("-c","--charge",dest="charge",action="store",
                type=int,default=0,help="Total charge (0)")

        parser.add_argument("--etemp",dest="etemp",action="store",
                type=float,default=300.0,help="Electronic temperature (300K)")

        parser.add_argument("--accuracy",dest="acc",action="store",
                type=float,default=1.0,help="Calculation accuracy (1.0)")

        parser.add_argument("--maxiter",dest="maxiter",action="store",
                type=int,default=250,help="Maximum number of SCC iterations (250)")

        parser.add_argument("--lcovb",dest="rnn",action="store",
                type=float,default=5.0,help="Longest covalent bond for preconditioning (5.0)")

        parser.add_argument("--econv",dest="econv",action="store",
                type=float,default=5.4423e-4,help="Convergence threshold for energy")

        parser.add_argument("--gconv",dest="fconv",action="store",
                type=float,default=3.88938e-05,help="Convergence threshold for forces")

        parser.add_argument("--linesearch",dest="armijo",action="store_true",
                help="Use line search")

        parser.add_argument("--precon",dest="precon",action="store_true",
                help="Use preconditioning")

        parser.add_argument("--optcell",dest="optcell",action="store_true",
                help="Relax cellparameters")

        parser.add_argument("--logfile",dest="logfile",action="store",
                default='-',help="File for logging information (STDOUT)")

        parser.add_argument("--trajectory",dest="trajectory",action="store",
                default='xtbopt.traj',help="File for trajectory output (xtbopt.traj)")

        parser.add_argument("--solvent",dest="solvent",action="store",
                default="none",help="Solvent for GBSA (none)")

        return parser.parse_args()

    import numpy
    import ase
    from ase.io import read, write
    from ase.units import *
    from ase.optimize.precon import Exp, PreconFIRE
    from ase.constraints import ExpCellFilter

    # overwrite convergence thresholds of PreconFIRE optimizer
    class PatchedOptimizer(PreconFIRE):

        def initialize(self):
            PreconFIRE.initialize(self)
            self.elast = None

        def run(self, steps=100000000, econv = 0.00054423,fconv = 3.88938e-05):
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
            fnorm = sqrt((forces**2).sum())
            fconverged = fnorm < self.fconv

            #print("--> norm(F):", fnorm, "converged?", fconverged)
            #print("--> energy :", ecurr, "converged?", econverged)

            return econverged and fconverged

    args = get_cmd_args()
    
    if args.filetype is not None:
        mol = ase.io.read(args.filename, format = str(args.filetype))
    else:
        mol = ase.io.read(args.filename)

    if args.method == 'gfn0':
        if mol.pbc.any():
            calc = GFN0_PBC(charge=args.charge,accuracy=args.acc,
                    temperature=args.etemp,print_level=2)
        else:
            calc = GFN0(charge=args.charge,accuracy=args.acc,temperature=args.etemp,
                    print_level=2)
    elif args.method == 'gfn1':
        if mol.pbc.any():
            raise Exception("GFN1-xTB is not available with PBC")
        else:
            calc = GFN1(charge=args.charge,accuracy=args.acc,temperature=args.etemp,
                    max_iterations=args.maxiter,print_level=2,solvent=args.solvent)
    elif args.method == 'gfn2':
        if mol.pbc.any():
            raise Exception("GFN2-xTB is not available with PBC")
        else:
            calc = GFN2(charge=args.charge,accuracy=args.acc,temperature=args.etemp,
                    max_iterations=args.maxiter,print_level=2,solvent=args.solvent)
    else:
        raise Exception("Method not implemented.")
    mol.set_calculator(calc)

    e = mol.get_potential_energy()

    print('Initial energy: eV, Eh',e, e/Hartree)

    if args.precon:
        precon = Exp(A=3, r_NN=args.rnn + 0.1, r_cut=args.rnn + 0.6)
    else:
        precon = None

    if args.optcell:
        sf = ExpCellFilter(mol)
        relax = PatchedOptimizer(sf, precon=precon,
                trajectory=args.trajectory, logfile=args.logfile, use_armijo=args.armijo)
    else:
        relax = PatchedOptimizer(mol, precon=precon,
                trajectory=args.trajectory, logfile=args.logfile, use_armijo=args.armijo)

    try:
        relax.run(econv=args.econv,fconv=args.fconv)
    except KeyboardInterrupt:
        print('User got impatient')
    except RuntimeError:
        print('Optimization terminated due to internal error')

    e = mol.get_potential_energy()
    print('Final energy: eV, Eh',e, e/Hartree)

    if (mol.pbc.any()):
        ase.io.write("xtbopt.POSCAR",mol,format='vasp')
    else:
        ase.io.write("xtbopt.xyz",mol,format='xyz')
