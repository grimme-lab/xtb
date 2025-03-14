! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

subroutine citation(iunit)
integer,intent(in) :: iunit
write(iunit,'(3x,a)') &
   "Cite this work as:", &
   "* C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht,",&
   "  J. Seibert, S. Spicher, S. Grimme, WIREs Comput. Mol. Sci., 2020, 11,",&
   "  e01493. DOI: 10.1002/wcms.1493",&
   "",&
   "for GFN2-xTB:",&
   "* C. Bannwarth, S. Ehlert and S. Grimme., J. Chem. Theory Comput., 2019,",&
   "  15, 1652-1671. DOI: 10.1021/acs.jctc.8b01176",&
   "for GFN1-xTB:",&
   "* S. Grimme, C. Bannwarth, P. Shushkov, J. Chem. Theory Comput., 2017,",&
   "  13, 1989-2009. DOI: 10.1021/acs.jctc.7b00118", &
   "for GFN0-xTB:",&
   "* P. Pracht, E. Caldeweyher, S. Ehlert, S. Grimme, ChemRxiv, 2019, preprint.",&
   "  DOI: 10.26434/chemrxiv.8326202.v1",&
   "for GFN-FF:",&
   "* S. Spicher and S. Grimme, Angew. Chem. Int. Ed., 2020, 59, 15665-15673.",&
   "  DOI: 10.1002/anie.202004239", &
   "",&
   "for ALPB and GBSA implicit solvation:",&
   "* S. Ehlert, M. Stahn, S. Spicher, S. Grimme, J. Chem. Theory Comput.,", &
   "  2021, 17, 4250-4261. DOI: 10.1021/acs.jctc.1c00471", &
   "",&
   "for ddCOSMO and CPCM-X implicit solvation:",&
   "* M. Stahn, S. Ehlert, S. Grimme, J. Phys. Chem. A,", &
   "  2023, XX, XXXX-XXXX. DOI: 10.1021/acs.jpca.3c04382", &
   "",&
   "for DFT-D4:",&
   "* E. Caldeweyher, C. Bannwarth and S. Grimme, J. Chem. Phys., 2017,",&
   "  147, 034112. DOI: 10.1063/1.4993215", &
   "* E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher,", &
   "  C. Bannwarth and S. Grimme, J. Chem. Phys., 2019, 150, 154122.", &
   "  DOI: 10.1063/1.5090222", &
   "* E. Caldeweyher, J.-M. Mewes, S. Ehlert and S. Grimme, Phys. Chem. Chem. Phys.",&
   "  2020, 22, 8499-8512. DOI: 10.1039/D0CP00502A",&
   "",&
   "for sTDA-xTB:",&
   "* S. Grimme and C. Bannwarth, J. Chem. Phys., 2016, 145, 054103.",&
   "  DOI: 10.1063/1.4959605",&
   "",&
   "in the mass-spec context:",&
   "* V. Asgeirsson, C. Bauer and S. Grimme, Chem. Sci., 2017, 8, 4879.",&
   "  DOI: 10.1039/c7sc00601b",&
   "* J. Koopman and S. Grimme, ACS Omega 2019, 4, 12, 15120-15133.",&
   "  DOI: 10.1021/acsomega.9b02011",&
   "",&
   "for metadynamics refer to:",&
   "* S. Grimme, J. Chem. Theory Comput., 2019, 155, 2847-2862", &
   "  DOI: 10.1021/acs.jctc.9b00143", &
   "",&
   "for SPH calculations refer to:",&
   "* S. Spicher and S. Grimme, J. Chem. Theory Comput., 2021, 17, 1701-1714", &
   "  DOI: 10.1021/acs.jctc.0c01306", &
   "",&
   "for ONIOM refer to:",&
   "* C. Plett, A. Katbashev, S. Ehlert, S. Grimme, M. Bursch,",&
   "  Phys. Chem. Chem. Phys., 2023, 25, 17860-17868. DOI: 10.1039/D3CP02178E",&
   "",&
   "for DIPRO refer to:",&
   "* J. Kohn, N. Gildemeister, S. Grimme, D. Fazzi, A. Hansen,",&
   "  J. Chem. Phys., 2023, just accepted.",&
   "",&
   "for PTB refer to:",&
   "* S. Grimme, M. Mueller, A. Hansen, J. Chem. Phys., 2023, 158, 124111.",&
   "  DOI: 10.1063/5.0137838",&
   "",&
   "with help from (in alphabetical order)",&
   "P. Atkinson, C. Bannwarth, F. Bohle, G. Brandenburg, E. Caldeweyher", &
   "M. Checinski, S. Dohm, S. Ehlert, S. Ehrlich, I. Gerasimov, C. Hölzer", &
   "A. Katbashev, J. Kohn, J. Koopman, C. Lavigne, S. Lehtola, F. März, M. Müller,", &
   "F. Musil, H. Neugebauer, J. Pisarek, C. Plett, P. Pracht, F. Pultar,", &
   "J. Seibert, P. Shushkov, S. Spicher, M. Stahn, M. Steiner, T. Strunk,", &
   "J. Stückrath, T. Rose, and J. Unsleber", &
   ""
end subroutine citation

subroutine help(iunit)
   implicit none
   integer, intent(in) :: iunit
   write(iunit,'(a)') &
   "Usage: xtb [options] <geometry> [options]", &
   "",&
   "<geometry> may be provided as valid TM coordinate file (*coord in Bohr),",&
   "in xmol format (*.xyz in Ångström), sdf or mol file format, PDB format",&
   "genFormat input or Vasp's POSCAR format.",&
   "",&
   "For a full list visit: https://grimme-lab.github.io/mctc-lib/page/index.html", &
   "",&
   "Options:",&
   "",&
   "-c, --chrg INT",&
   "    specify molecular charge as INT, overrides .CHRG file and xcontrol option",&
   "",&
   "--ceh <grad> [REAL]",&
   "    calculate CEH (Charge-Extended Hückel model) charges and write them to the ceh.charges file",&
   "    optionally, calculate numerical gradients and write them to the ceh.charges.numgrad file",&
   "    with an adjustable step size for the numerical gradients.",&
   "",&
   "-u, --uhf INT",&
   "    specify number of unpaired electrons as INT, overrides .UHF file and xcontrol option",&
   "",&
   "--gfn INT",&
   "    specify parametrisation of GFN-xTB (default = 2)",&
   "",&
   "--gfnff, --gff",&
   "    specify parametrisation of GFN-FF",&
   "",&
   "--tblite,",&
   "    use tblite library as implementation for xTB, please note that solvation is not yet implemented within tblite",&
   "",&
   "--ptb,",&
   "    performs single-point calculation with the density tight-binding method PTB.", &
   "    Provides electronic structure and properties, such as, e.g., atomic charges, bond orders, and dipole moments,", & 
   "    but does not provide any energy-related properties, such as, e.g., total energy, nuclear gradients, or vibrational frequencies.", &
   "",&
   "--spinpol,",&
   "    enables spin-polarization for xTB methods (tblite required)",&
   "",&
   "--oniom METHOD LIST", &
   "    use subtractive embedding via ONIOM method. 'METHOD' is given as 'inner:outer'", &
   "    where 'inner' can be 'orca', 'turbomole', 'gfn2', 'gfn1', or 'gfnff' and", &
   "    'outer' can be 'gfn2', 'gfn1', or 'gfnff'.", &
   "    The inner region is given as a comma separated indices directly in the commandline", &
   "    or in a file with each index on a separate line.", &
   "",&
   "--etemp REAL",&
   "    electronic temperature (default = 300K)",&
   "",&
   "--esp",&
   "    calculate electrostatic potential on VdW-grid",&
   "",&
   "--stm",&
   "    calculate STM image",&
   "",&
   "-a, --acc REAL",&
   "    accuracy for SCC calculation, lower is better (default = 1.0)",&
   "",&
   "--vparam FILE",&
   "    Parameter file for vTB calculation",&
   "",&
   "--xparam FILE",&
   "    Parameter file for xTB calculation (not used)",&
   "",&
   "--alpb SOLVENT [STATE]",&
   "    analytical linearized Poisson-Boltzmann (ALPB) model,",&
   "    available solvents are acetone, acetonitrile, aniline, benzaldehyde,",&
   "    benzene, ch2cl2, chcl3, cs2, dioxane, dmf, dmso, ether, ethylacetate, furane,", &
   "    hexandecane, hexane, methanol, nitromethane, octanol, woctanol, phenol, toluene,",&
   "    thf, water.",&
   "    The solvent input is not case-sensitive. The Gsolv",&
   "    reference state can be chosen as reference, bar1M, or gsolv (default).",&
   "",&
   "-g, --gbsa SOLVENT [STATE]",&
   "    generalized born (GB) model with solvent accessable surface (SASA) model,",&
   "    available solvents are acetone, acetonitrile, benzene (only GFN1-xTB), CH2Cl2,",&
   "    CHCl3, CS2, DMF (only GFN2-xTB), DMSO, ether, H2O, methanol,",&
   "    n-hexane (only GFN2-xTB), THF and toluene.",&
   "    The solvent input is not case-sensitive.", &
   "    The Gsolv reference state can be chosen as reference, bar1M, or gsolv (default).",&
   "",&
   "--cosmo SOLVENT/EPSILON",&
   "    domain decomposition conductor-like screening model (ddCOSMO),",&
   "    available solvents are all solvents that are available for alpb.",&
   "    Additionally, the dielectric constant can be set manually or an ideal conductor", &
   "    can be chosen by setting epsilon to infinity.",&
   "",&
   "--tmcosmo SOLVENT/EPSILON",&
   "    same as --cosmo, but uses TM convention for writing the .cosmo files.",&
   "",&
   "--cpcmx SOLVENT",&
   "    extended conduction-like polarizable continuum solvation model (CPCM-X),",&
   "    available solvents are all solvents included in the Minnesota Solvation Database.",&
   "",&
   "--cma",&
   "    shifts molecule to center of mass and transforms cartesian coordinates into the",&
   "    coordinate system of the principle axis (not affected by ‘isotopes’-file).",&
   "",&
   "--pop",&
   "    requests printout of Mulliken population analysis",&
   "",&
   "--molden",&
   "    requests printout of molden file",&
   "",&
   "--dipole",&
   "    requests dipole printout",&
   "",&
   "--raman", &
   "    requests Raman spectrum calculation via combination of GFN2-xTB and PTB",&
   "    using the temperature 'REAL' (default 298.15 K) and the wavelength of the",& 
   "    incident laser which must be given in nm 'REAL' (default 514 nm).",&
   "",&
   "--alpha",&
   "    requests the extension of electrical properties to static molecular dipole polarizabilities. ",&
   "",&
   "--wbo",&
   "    requests Wiberg bond order printout",&
   "",&
   "--lmo",&
   "    requests localization of orbitals",&
   "",&
   "--fod",&
   "    requests FOD calculation",&
   "",&
   "--scc, --sp",&
   "    performs a single point calculation",&
   "",&
   "--vip",&
   "    performs calculation of ionisation potential. This needs the param_ipea-xtb.txt",&
   "    parameters and a GFN1 Hamiltonian.",&
   "",&
   "--vea",&
   "    performs calculation of electron affinity. This needs the param_ipea-xtb.txt",&
   "    parameters and a GFN1 Hamiltonian.",&
   "",&
   "--vipea",&
   "    performs calculation of electron affinity and ionisation potential.",&
   "    This needs the param_ipea-xtb.txt parameters and a GFN1 Hamiltonian.",&
   "",&
   "--vfukui",&
   "    performs calculation of Fukui indices.",&
   "",&
   "--vomega",&
   "    performs calculation of electrophilicity index. This needs the param_ipea-xtb.txt",&
   "    parameters and a GFN1 Hamiltonian.",&
   "",&
   "--grad",&
   "    performs a gradient calculation",&
   "",&
   "-o, --opt [LEVEL]",&
   "    call ancopt(3) to perform a geometry optimization, levels from crude, sloppy,", &
   "    loose, normal (default), tight, verytight to extreme can be chosen",&
   "",&
   "--hess",&
   "    perform a numerical hessian calculation on input geometry",&
   "",&
   "--ohess [LEVEL]",&
   "    perform a numerical hessian calculation on an ancopt(3) optimized geometry",&
   "",&
   "--bhess [LEVEL]",&
   "    perform a biased numerical hessian calculation on an ancopt(3) optimized geometry",&
   "",&
   "--md",&
   "    molecular dynamics simulation on start geometry",&
   "",&
   "--metadyn [int]",&
   "    meta dynamics simulation on start geometry, saving int snapshots of the", &
   "    trajectory to bias the simulation",&
   "",&
   "--omd",&
   "    molecular dynamics simulation on ancopt(3) optimized geometry,", &
   "    a loose optimization level will be chosen",&
   "",&
   "--metaopt [LEVEL]",&
   "    call ancopt(3) to perform a geometry optimization, then try to find other", &
   "    minimas by meta dynamics",&
   "",&
   "--path [FILE]",&
   "    use meta dynamics to calculate a path from the input geometry to the given", &
   "    product structure",&
   "",&
   "--modef INT",&
   "    modefollowing algorithm. INT specifies the mode that should be used for the", &
   "    modefollowing.",&
   "",&
   "-I, --input FILE",&
   "    use FILE as input source for xcontrol(7) instructions",&
   "",&
   "--namespace STRING",&
   "    give this xtb(1) run a namespace. All files, even temporary ones, will be named",&
   "    according to STRING (might not work everywhere).",&
   "",&
   "--[no]copy",&
   "    copies the xcontrol file at startup (default = true)",&
   "",&
   "--[no]restart",&
   "    restarts calculation from xtbrestart (default = true)",&
   "",&
   "-P, --parallel INT",&
   "    number of parallel processes",&
   "",&
   "--define",&
   "    performs automatic check of input and terminate",&
   "",&
   "--json",&
   "    write xtbout.json file",&
   "",&
   "--citation",&
   "    print citation and terminate",&
   "",&
   "--license",&
   "    print license and terminate",&
   "",&
   "-v, --verbose",&
   "    be more verbose (not supported in every unit)",&
   "",&
   "-s, --silent",&
   "    clutter the screen less (not supported in every unit)",&
   "",&
   "--ceasefiles",&
   "    reduce the amount of output and files written",&
   "",&
   "--strict",&
   "    turns all warnings into hard errors",&
   "",&
   "-h, --help",&
   "    show help page",&
   "",&
   "Useful Maschine Settings:",&
   "",&
   "export MKL_NUM_THREADS=<NCORE>",&
   "export OMP_NUM_THREADS=<NCORE>,1",&
   "export OMP_STACKSIZE=4G",&
   "ulimit -s unlimited",&
   "",&
   "Output Conventions:",&
   "",&
   "total energies are given in atomic units (Eh)",&
   "gaps/HL energies are given in eV",&
   "",&
   "More information can be obtained by `man 1 xtb` and `man 7 xcontrol`",&
   "or at https://xtb-docs.readthedocs.io/en/latest/contents.html",&
   ""
end subroutine help

subroutine help_legacy
   use, intrinsic :: iso_fortran_env, only : id => output_unit
   write(id,'(''Usage: xtb <geometry> [options]'',/)')

   write(id,'(''<geometry> may be provided as'','//&
   &    'x,''valid TM coordinate file (*coord in Bohr) or'','//&
   &    '/,''in xmol format (*xyz in Ångström).'',/)')

   write(id,'(''Options:'',/)')

   write(id,'(3x,''-c, --chrg <int>  '','// &
   &          'x,''Set charge to molecule'','//&
   &      '/,22x,''overrides charge in .CHRG file'')')

   write(id,'(3x,''-u, --uhf <int>   '','// &
   &          'x,''Number pf unpaired electrons or 2S'','//&
   &      '/,22x,''overrides setting in .UHF file'')')

   write(id,'(3x,''    --nox/--nodiff'','// &
   &          'x,''skip second, extended part in sTDA-xTB'')')

   write(id,'(3x,''    --pop         '','// &
   &          'x,''do population analysis'')')

   write(id,'(3x,''    --mowr <real> '','// &
   &          'x,''cut MO write above (default 3 Eh)'')')

   write(id,'(3x,''    --mold[en]    '','// &
   &          'x,''write formatted molden file'')')

   write(id,'(3x,''    --dip[ole]    '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --wbo         '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --lmo         '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --rdchrg      '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --md          '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --hessl       '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --optts       '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --ohessl      '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --mdopt       '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --path        '','// &
   &          'x,''deactivated feature'')')

   write(id,'(3x,''    --vip         '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --vea         '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --tm          '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --orca        '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --LS          '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --mopac       '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --qmdff       '','// &
   &          'x,''hidden feature'')')

   write(id,'(3x,''    --parx <file> '','// &
   &          'x,''read parameters for sTDA-xTB calculation'',' // &
   &      '/,22x,''(default: $XTBHOME/.param_stda2.xtb)'')')

   write(id,'(3x,''    --parv <file> '','// &
   &          'x,''read parameters for vTB part in sTDA-xTB'',' // &
   &      '/,22x,''(default: $XTBHOME/.param_stda1.xtb)'')')

   write(id,'(3x,''    --xtemp <real>'','// &
   &          'x,''electronic temperature for xTB (default 0 K)'')')

   write(id,'(3x,''    --etemp <real>'','// &
   &          'x,''electronic temperature for GFN (default 300 K)'')')

   write(id,'(3x,''    --fod         '','// &
   &          'x,''calculate the FOD and write molden.input file'','// &
   &      '/,22x,''with appropriate occupations for plotting'','// &
   &      '/,22x,''(defualt T = 12500 K)'')')

   write(id,'(3x,''    --gfn1        '','// &
   &          'x,''use GFN1-xTB parametrisation'')')

   write(id,'(3x,''    --gfn2        '','// &
   &          'x,''use GFN2-xTB parametrisation'')')

   write(id,'(3x,''    --gfn2d3      '','// &
   &          'x,''GFN2-xTB with D3 instead of D4 dispersion'')')

   write(id,'(3x,''    --scc         '','// &
   &          'x,''do single point only'')')

   write(id,'(3x,''    --grad        '','// &
   &          'x,''calculate GFNn-xTB gradient'')')

   write(id,'(3x,''    --acc <real>  '','// &
   &          'x,''GFNn-xTB accuracy (default 1.0)'')')

   write(id,'(3x,''    --opt [level] '','// &
   &          'x,''optimize at GFNn-xTB level, level can be one of'',' // &
   &      '/,22x,''crude, vloose, loose, tight, vtight, extreme'')')

   write(id,'(3x,''    --hess        '','// &
   &          'x,''compute Hessian at GFNn-xTB level'')')

   write(id,'(3x,''    --ohess       '','// &
   &          'x,''optimize and compute Hessian'')')

   write(id,'(3x,''    --omd         '','// &
   &          'x,''optimize and do MD'')')

   write(id,'(3x,''    --dipro <thr> '','// &
   &          'x,''calculate dipro coupling integrals'')')

   write(id,'(3x,''    --siman       '','// &
   &          'x,''conformational search'')')

   write(id,'(3x,''    --screen      '','// &
   &          'x,''opt. loop over ensemble'')')

   write(id,'(3x,''    --gmd         '','// &
   &          'x,''annealed MD for GMD procedure'')')

   write(id,'(3x,''    --modef       '','// &
   &          'x,''vibrational mode following'')')

   write(id,'(3x,''-g, --gbsa [solv] [state]'','// &
   &      '/,22x,''use GBSA implicit solvent for solvent [solv]'',' // &
   &      '/,22x,''and solvation state [state] = reference, bar1M'','// &
   &      '/,22x,''(default: bar1M)'')')

   write(id,'(3x,''    --help        '','// &
   &          'x,''show this message'')')

   write(id,'(/,''Useful Maschine Settings:'',/)')

   write(id,'(3x,''setenv MKL_NUM_THREADS <NCORE>'')')
   write(id,'(3x,''setenv OMP_STACKSIZE 500m'')')
   write(id,'(3x,''limit stacksize unlimited'')')

   write(id,'(/,''Output Conventions:'',/)')

   write(id,'(3x,''total energies are given in atomic units (Eh)'')')
   write(id,'(3x,''gaps/HL energies are given in eV'')')

   write(id,'(/,''Please read the manual carefully'',' // &
   &         'x,''(just kidding... there is no manual)'',/)')
end subroutine help_legacy

subroutine definebanner
   use, intrinsic :: iso_fortran_env, only : id => output_unit
   write(id,"(""           _"")")
   write(id,"(""          | |"")")
   write(id,"(""  _ __ ___| |_ _   _ _ __ _ __        __"")")
   write(id,"("" | '__/ _ \ __| | | | '__| '_ \      / _|"")")
   write(id,"("" | | |  __/ |_| |_| | |  | | | |___ | |_"")")
   write(id,"("" |_|  \___|\__|\__,_|_|  |_| |_/ _ \|  _|"")")
   write(id,"(""         _       __ _         | (_) | |"")")
   write(id,"(""        | |     / _(_)         \___/|_|"")")
   write(id,"(""      __| | ___| |_ _ _ __   ___"")")
   write(id,"(""     / _` |/ _ \  _| | '_ \ / _ \"")")
   write(id,"(""    | (_| |  __/ | | | | | |  __/"")")
   write(id,"(""     \__,_|\___|_| |_|_| |_|\___|"")")
end subroutine definebanner
